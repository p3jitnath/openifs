#include "oifs_cuda_radiation.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

namespace {

std::size_t layer_index(int g, int lev, int col, int ng, int nlev) {
  return g + static_cast<std::size_t>(ng) *
      (lev + static_cast<std::size_t>(nlev) * col);
}

std::size_t interface_index(int g, int lev, int col, int ng, int nlev) {
  return g + static_cast<std::size_t>(ng) *
      (lev + static_cast<std::size_t>(nlev + 1) * col);
}

struct Inputs {
  int ng = 14, nbands = 6, nlev = 17, ncol = 257;
  std::vector<double> mu0, od, ssa, asymmetry, albedo_direct;
  std::vector<double> albedo_diffuse, incoming, cloud_fraction;
  std::vector<double> total_cloud_cover, od_scaling, od_cloud;
  std::vector<double> ssa_cloud, asymmetry_cloud;
  std::vector<int> band_from_g;

  Inputs() {
    const std::size_t layer = static_cast<std::size_t>(ng) * nlev * ncol;
    const std::size_t gpcol = static_cast<std::size_t>(ng) * ncol;
    const std::size_t cloud_layer = static_cast<std::size_t>(nbands) * nlev * ncol;
    mu0.resize(ncol);
    od.resize(layer); ssa.resize(layer); asymmetry.resize(layer);
    albedo_direct.resize(gpcol); albedo_diffuse.resize(gpcol); incoming.resize(gpcol);
    cloud_fraction.resize(static_cast<std::size_t>(ncol) * nlev);
    total_cloud_cover.resize(ncol); od_scaling.resize(layer);
    od_cloud.resize(cloud_layer); ssa_cloud.resize(cloud_layer);
    asymmetry_cloud.resize(cloud_layer); band_from_g.resize(ng);

    for (int g = 0; g < ng; ++g) band_from_g[g] = 1 + (g * nbands) / ng;
    for (int col = 0; col < ncol; ++col) {
      mu0[col] = col % 11 == 0 ? -0.05 : 0.12 + 0.83 * ((col % 37) / 36.0);
      total_cloud_cover[col] = mu0[col] > 0.0 ? (col % 9) / 10.0 : 0.0;
      for (int g = 0; g < ng; ++g) {
        const auto gp = g + static_cast<std::size_t>(ng) * col;
        albedo_direct[gp] = 0.03 + 0.25 * ((g + 2 * col) % 13) / 12.0;
        albedo_diffuse[gp] = 0.05 + 0.30 * ((2 * g + col) % 11) / 10.0;
        incoming[gp] = 0.4 + 1.1 * ((3 * g + col) % 17) / 16.0;
      }
      for (int lev = 0; lev < nlev; ++lev) {
        cloud_fraction[col + static_cast<std::size_t>(ncol) * lev] =
            (col + 3 * lev) % 7 == 0 ? 0.0 : 0.15 + 0.1 * ((col + lev) % 6);
        for (int g = 0; g < ng; ++g) {
          const auto at = layer_index(g, lev, col, ng, nlev);
          od[at] = 0.002 + 0.12 * ((g + 2 * lev + col) % 23) / 22.0;
          ssa[at] = 0.08 + 0.84 * ((2 * g + lev + col) % 19) / 18.0;
          asymmetry[at] = -0.05 + 0.77 * ((g + lev + 2 * col) % 17) / 16.0;
          od_scaling[at] = 0.15 + 1.45 * ((3 * g + lev + col) % 21) / 20.0;
        }
        for (int band = 0; band < nbands; ++band) {
          const auto at = band + static_cast<std::size_t>(nbands) *
              (lev + static_cast<std::size_t>(nlev) * col);
          od_cloud[at] = 0.01 + 0.55 * ((band + lev + col) % 15) / 14.0;
          ssa_cloud[at] = 0.72 + 0.26 * ((2 * band + lev + col) % 13) / 12.0;
          asymmetry_cloud[at] = 0.45 + 0.42 * ((band + 2 * lev + col) % 11) / 10.0;
        }
      }
    }
  }
};

struct Outputs {
  std::array<std::vector<double>, 6> profile;
  std::array<std::vector<double>, 4> surface;
  Outputs(int ng, int nlev, int ncol) {
    for (auto& values : profile) values.resize(static_cast<std::size_t>(ncol) * (nlev + 1));
    for (auto& values : surface) values.resize(static_cast<std::size_t>(ng) * ncol);
  }
};

void optics(const Inputs& in, bool cloudy, bool delta, double threshold,
            std::array<std::vector<double>, 5>& coefficient) {
  const std::size_t count = static_cast<std::size_t>(in.ng) * in.nlev * in.ncol;
  for (auto& values : coefficient) values.resize(count);
  for (int col = 0; col < in.ncol; ++col) {
    for (int lev = 0; lev < in.nlev; ++lev) {
      for (int gp = 0; gp < in.ng; ++gp) {
        const auto at = layer_index(gp, lev, col, in.ng, in.nlev);
        if (in.mu0[col] <= 0.0) {
          for (auto& values : coefficient) values[at] = 0.0;
          continue;
        }
        double od = in.od[at], ssa = in.ssa[at], asym = in.asymmetry[at];
        if (cloudy && in.cloud_fraction[col + static_cast<std::size_t>(in.ncol) * lev] >= threshold) {
          const int band = in.band_from_g[gp] - 1;
          const auto cloud = band + static_cast<std::size_t>(in.nbands) *
              (lev + static_cast<std::size_t>(in.nlev) * col);
          const double cloud_od = in.od_scaling[at] * in.od_cloud[cloud];
          od += cloud_od;
          ssa = od > 0.0 ? (in.ssa[at] * in.od[at] + in.ssa_cloud[cloud] * cloud_od) / od : 0.0;
          const double scattering_od = ssa * od;
          asym = scattering_od > 0.0
              ? (in.asymmetry[at] * in.ssa[at] * in.od[at]
                 + in.asymmetry_cloud[cloud] * in.ssa_cloud[cloud] * cloud_od) / scattering_od
              : 0.0;
        }
        if (delta) {
          const double forward = asym * asym;
          od *= 1.0 - ssa * forward;
          ssa = ssa * (1.0 - forward) / (1.0 - ssa * forward);
          asym /= 1.0 + asym;
        }
        const double factor = 0.75 * asym;
        const double gamma1 = 2.0 - ssa * (1.25 + factor);
        const double gamma2 = ssa * (0.75 - factor);
        const double gamma3 = 0.5 - in.mu0[col] * factor;
        const double gamma4 = 1.0 - gamma3;
        const double alpha1 = gamma1 * gamma4 + gamma2 * gamma3;
        const double alpha2 = gamma1 * gamma3 + gamma2 * gamma4;
        const double k = std::sqrt(std::max((gamma1 - gamma2) * (gamma1 + gamma2), 1.0e-12));
        double mu = in.mu0[col];
        if (std::abs(1.0 - k * mu) < 1000.0 * 2.2204460492503131e-16)
          mu *= 1.0 - 1000.0 * 2.2204460492503131e-16;
        const double kmu = k * mu;
        const double exp0 = std::exp(-std::max(od / mu, 0.0));
        const double exp1 = std::exp(-k * od), exp2 = exp1 * exp1;
        const double two_k_exp = 2.0 * k * exp1;
        double denominator = 1.0 / (k + gamma1 + (k - gamma1) * exp2);
        coefficient[0][at] = gamma2 * (1.0 - exp2) * denominator;
        coefficient[1][at] = two_k_exp * denominator;
        coefficient[4][at] = exp0;
        denominator *= mu * ssa / (1.0 - kmu * kmu);
        const double direct_ref = denominator * ((1.0 - kmu) * (alpha2 + k * gamma3)
            - (1.0 + kmu) * (alpha2 - k * gamma3) * exp2
            - two_k_exp * (gamma3 - alpha2 * mu) * exp0);
        coefficient[2][at] = std::max(0.0, std::min(direct_ref, 1.0));
        const double direct_diffuse = denominator * (two_k_exp * (gamma4 + alpha1 * mu)
            - exp0 * ((1.0 + kmu) * (alpha1 + k * gamma4)
            - (1.0 - kmu) * (alpha1 - k * gamma4) * exp2));
        coefficient[3][at] = std::max(0.0, std::min(direct_diffuse, 1.0 - coefficient[2][at]));
      }
    }
  }
}

void adding(const Inputs& in, const std::array<std::vector<double>, 5>& c,
            std::array<std::vector<double>, 3>& flux) {
  const std::size_t count = static_cast<std::size_t>(in.ng) * (in.nlev + 1) * in.ncol;
  for (auto& values : flux) values.assign(count, 0.0);
  std::vector<double> albedo(in.nlev + 1), source(in.nlev + 1), inverse(in.nlev);
  for (int col = 0; col < in.ncol; ++col) {
    if (in.mu0[col] <= 0.0) continue;
    for (int gp = 0; gp < in.ng; ++gp) {
      const auto gpcol = gp + static_cast<std::size_t>(in.ng) * col;
      flux[2][interface_index(gp, 0, col, in.ng, in.nlev)] = in.incoming[gpcol];
      for (int lev = 0; lev < in.nlev; ++lev) {
        const auto layer = layer_index(gp, lev, col, in.ng, in.nlev);
        flux[2][interface_index(gp, lev + 1, col, in.ng, in.nlev)] =
            flux[2][interface_index(gp, lev, col, in.ng, in.nlev)] * c[4][layer];
      }
      albedo[in.nlev] = in.albedo_diffuse[gpcol];
      source[in.nlev] = in.albedo_direct[gpcol]
          * flux[2][interface_index(gp, in.nlev, col, in.ng, in.nlev)] * in.mu0[col];
      for (int lev = in.nlev - 1; lev >= 0; --lev) {
        const auto layer = layer_index(gp, lev, col, in.ng, in.nlev);
        inverse[lev] = 1.0 / (1.0 - albedo[lev + 1] * c[0][layer]);
        albedo[lev] = c[0][layer] + c[1][layer] * c[1][layer]
            * albedo[lev + 1] * inverse[lev];
        source[lev] = c[2][layer] * flux[2][interface_index(gp, lev, col, in.ng, in.nlev)]
            + c[1][layer] * (source[lev + 1] + albedo[lev + 1] * c[3][layer]
            * flux[2][interface_index(gp, lev, col, in.ng, in.nlev)]) * inverse[lev];
      }
      flux[0][interface_index(gp, 0, col, in.ng, in.nlev)] = source[0];
      for (int lev = 0; lev < in.nlev; ++lev) {
        const auto layer = layer_index(gp, lev, col, in.ng, in.nlev);
        const auto above = interface_index(gp, lev, col, in.ng, in.nlev);
        const auto below = interface_index(gp, lev + 1, col, in.ng, in.nlev);
        flux[1][below] = (c[1][layer] * flux[1][above] + c[0][layer] * source[lev + 1]
            + c[3][layer] * flux[2][above]) * inverse[lev];
        flux[0][below] = albedo[lev + 1] * flux[1][below] + source[lev + 1];
        flux[2][above] *= in.mu0[col];
      }
      flux[2][interface_index(gp, in.nlev, col, in.ng, in.nlev)] *= in.mu0[col];
    }
  }
}

Outputs reference(const Inputs& in, bool delta) {
  std::array<std::vector<double>, 5> clear_coefficients, cloud_coefficients;
  std::array<std::vector<double>, 3> clear_flux, cloud_flux;
  optics(in, false, delta, 1.0e-6, clear_coefficients);
  adding(in, clear_coefficients, clear_flux);
  optics(in, true, delta, 1.0e-6, cloud_coefficients);
  adding(in, cloud_coefficients, cloud_flux);
  Outputs out(in.ng, in.nlev, in.ncol);
  for (int col = 0; col < in.ncol; ++col) {
    for (int lev = 0; lev <= in.nlev; ++lev) {
      const auto profile = col + static_cast<std::size_t>(in.ncol) * lev;
      double up0 = 0.0, diffuse0 = 0.0, direct0 = 0.0;
      double up1 = 0.0, diffuse1 = 0.0, direct1 = 0.0;
      for (int gp = 0; gp < in.ng; ++gp) {
        const auto at = interface_index(gp, lev, col, in.ng, in.nlev);
        up0 += clear_flux[0][at]; diffuse0 += clear_flux[1][at]; direct0 += clear_flux[2][at];
        up1 += cloud_flux[0][at]; diffuse1 += cloud_flux[1][at]; direct1 += cloud_flux[2][at];
      }
      const double cloud = in.total_cloud_cover[col], clear = 1.0 - cloud;
      out.profile[0][profile] = up0;
      out.profile[1][profile] = diffuse0 + direct0;
      out.profile[2][profile] = direct0;
      out.profile[3][profile] = cloud * up1 + clear * up0;
      out.profile[4][profile] = cloud * (diffuse1 + direct1) + clear * (diffuse0 + direct0);
      out.profile[5][profile] = cloud * direct1 + clear * direct0;
    }
    for (int gp = 0; gp < in.ng; ++gp) {
      const auto surface = gp + static_cast<std::size_t>(in.ng) * col;
      const auto at = interface_index(gp, in.nlev, col, in.ng, in.nlev);
      const double cloud = in.total_cloud_cover[col], clear = 1.0 - cloud;
      out.surface[0][surface] = clear_flux[1][at];
      out.surface[1][surface] = clear_flux[2][at];
      out.surface[2][surface] = cloud * cloud_flux[1][at] + clear * clear_flux[1][at];
      out.surface[3][surface] = cloud * cloud_flux[2][at] + clear * clear_flux[2][at];
    }
  }
  return out;
}

double compare(const std::vector<double>& expected, const std::vector<double>& actual,
               const char* name) {
  double worst = 0.0;
  std::size_t worst_at = 0;
  for (std::size_t i = 0; i < expected.size(); ++i) {
    const double ratio = std::abs(expected[i] - actual[i]) /
        (5.0e-10 + 2.0e-12 * std::max(std::abs(expected[i]), std::abs(actual[i])));
    if (ratio > worst) { worst = ratio; worst_at = i; }
  }
  if (worst > 1.0)
    std::fprintf(stderr, "%s differs at %zu: cpu=%.17g gpu=%.17g normalized_error=%g\n",
                 name, worst_at, expected[worst_at], actual[worst_at], worst);
  return worst;
}

int run_case(const Inputs& in, bool delta) {
  Outputs gpu(in.ng, in.nlev, in.ncol);
  const int status = oifs_cuda_sw_compute_dp(
      in.ng, in.nbands, in.nlev, in.ncol, delta, in.mu0.data(), in.od.data(),
      in.ssa.data(), in.asymmetry.data(), in.albedo_direct.data(),
      in.albedo_diffuse.data(), in.incoming.data(), in.band_from_g.data(), 1.0e-6,
      in.cloud_fraction.data(), in.total_cloud_cover.data(), in.od_scaling.data(),
      in.od_cloud.data(), in.ssa_cloud.data(), in.asymmetry_cloud.data(),
      gpu.profile[0].data(), gpu.profile[1].data(), gpu.profile[2].data(),
      gpu.profile[3].data(), gpu.profile[4].data(), gpu.profile[5].data(),
      gpu.surface[0].data(), gpu.surface[1].data(), gpu.surface[2].data(),
      gpu.surface[3].data());
  if (status != 0) {
    std::fprintf(stderr, "CUDA shortwave failed (%d): %s\n", status,
                 oifs_cuda_radiation_last_error());
    return 1;
  }
  const Outputs cpu = reference(in, delta);
  static const char* profile_names[] = {"sw_up_clear", "sw_dn_clear", "sw_dn_direct_clear",
      "sw_up", "sw_dn", "sw_dn_direct"};
  static const char* surface_names[] = {"sw_dn_diffuse_surf_clear_g",
      "sw_dn_direct_surf_clear_g", "sw_dn_diffuse_surf_g", "sw_dn_direct_surf_g"};
  double worst = 0.0;
  for (int i = 0; i < 6; ++i) worst = std::max(worst, compare(cpu.profile[i], gpu.profile[i], profile_names[i]));
  for (int i = 0; i < 4; ++i) worst = std::max(worst, compare(cpu.surface[i], gpu.surface[i], surface_names[i]));
  if (worst > 1.0) return 2;
  std::printf("CUDA parity passed (delta=%d, %d columns): normalized max error=%g\n",
              delta, in.ncol, worst);
  return 0;
}

}  // namespace

int main() {
  if (!oifs_cuda_radiation_available()) {
    std::fprintf(stderr, "No CUDA device available: %s\n", oifs_cuda_radiation_last_error());
    return 77;
  }
  const Inputs inputs;
  int status = run_case(inputs, false);
  if (status == 0) status = run_case(inputs, true);
  oifs_cuda_radiation_finalize();
  return status;
}
