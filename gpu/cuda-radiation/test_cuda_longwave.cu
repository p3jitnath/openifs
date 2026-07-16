#include "oifs_cuda_radiation.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <vector>

namespace {

std::size_t layer(int g, int lev, int col, int ng, int nlev) {
  return g + static_cast<std::size_t>(ng) *
      (lev + static_cast<std::size_t>(nlev) * col);
}

std::size_t interface_at(int g, int lev, int col, int ng, int nlev) {
  return g + static_cast<std::size_t>(ng) *
      (lev + static_cast<std::size_t>(nlev + 1) * col);
}

struct Inputs {
  int ng = 16, nbands = 5, nlev = 13, ncol = 259;
  double threshold = 1.0e-6;
  std::vector<double> od, ssa, asymmetry, planck, emission, albedo;
  std::vector<double> cloud_fraction, total_cloud_cover, od_scaling;
  std::vector<double> od_cloud, ssa_cloud, asymmetry_cloud;
  std::vector<int> band_from_g;

  Inputs() {
    const auto layers = static_cast<std::size_t>(ng) * nlev * ncol;
    const auto interfaces = static_cast<std::size_t>(ng) * (nlev + 1) * ncol;
    const auto gpcols = static_cast<std::size_t>(ng) * ncol;
    const auto cloud_layers = static_cast<std::size_t>(nbands) * nlev * ncol;
    od.resize(layers); ssa.resize(layers); asymmetry.resize(layers);
    planck.resize(interfaces); emission.resize(gpcols); albedo.resize(gpcols);
    cloud_fraction.resize(static_cast<std::size_t>(ncol) * nlev);
    total_cloud_cover.resize(ncol); od_scaling.resize(layers);
    od_cloud.resize(cloud_layers); ssa_cloud.resize(cloud_layers);
    asymmetry_cloud.resize(cloud_layers); band_from_g.resize(ng);
    for (int g = 0; g < ng; ++g) band_from_g[g] = 1 + (g * nbands) / ng;
    for (int col = 0; col < ncol; ++col) {
      total_cloud_cover[col] = col % 17 == 0 ? 0.0 :
          (col % 19 == 0 ? 5.0e-7 : 0.05 + 0.9 * ((col % 23) / 22.0));
      for (int g = 0; g < ng; ++g) {
        const auto gp = g + static_cast<std::size_t>(ng) * col;
        emission[gp] = 0.2 + 1.5 * ((g + 3 * col) % 29) / 28.0;
        albedo[gp] = 0.01 + 0.08 * ((2 * g + col) % 13) / 12.0;
        for (int lev = 0; lev <= nlev; ++lev) {
          planck[interface_at(g, lev, col, ng, nlev)] =
              0.1 + 1.9 * ((g + 2 * lev + col) % 37) / 36.0;
        }
      }
      for (int lev = 0; lev < nlev; ++lev) {
        cloud_fraction[col + static_cast<std::size_t>(ncol) * lev] =
            (col + 2 * lev) % 7 == 0 ? 0.0 : 0.1 + 0.12 * ((col + lev) % 7);
        for (int g = 0; g < ng; ++g) {
          const auto at = layer(g, lev, col, ng, nlev);
          od[at] = (g + lev + col) % 31 == 0 ? 4.0e-4 :
              0.004 + 0.45 * ((g + 2 * lev + col) % 27) / 26.0;
          ssa[at] = 0.01 + 0.65 * ((2 * g + lev + col) % 21) / 20.0;
          asymmetry[at] = -0.1 + 0.75 * ((g + lev + 3 * col) % 19) / 18.0;
          od_scaling[at] = 0.05 + 1.5 * ((g + 3 * lev + col) % 25) / 24.0;
        }
        for (int band = 0; band < nbands; ++band) {
          const auto at = band + static_cast<std::size_t>(nbands) *
              (lev + static_cast<std::size_t>(nlev) * col);
          od_cloud[at] = 0.002 + 0.7 * ((band + lev + col) % 17) / 16.0;
          ssa_cloud[at] = 0.15 + 0.72 * ((2 * band + lev + col) % 15) / 14.0;
          asymmetry_cloud[at] = 0.05 + 0.78 * ((band + 2 * lev + col) % 13) / 12.0;
        }
      }
    }
  }
};

struct Outputs {
  std::array<std::vector<double>, 4> profile;
  std::array<std::vector<double>, 2> surface;
  std::vector<double> derivative;
  Outputs(int ng, int nlev, int ncol) {
    for (auto& v : profile) v.resize(static_cast<std::size_t>(ncol) * (nlev + 1));
    for (auto& v : surface) v.resize(static_cast<std::size_t>(ng) * ncol);
    derivative.resize(static_cast<std::size_t>(ncol) * (nlev + 1));
  }
};

void coefficients(double od, double ssa, double asym, double top, double bottom,
                  bool scattering, double& ref, double& trans,
                  double& source_up, double& source_dn) {
  constexpr double diffusivity = 1.66;
  if (!scattering) {
    ref = 0.0;
    const double scaled = diffusivity * od;
    if (od > 1.0e-3) {
      trans = std::exp(-scaled);
      const double c = (bottom - top) / scaled;
      source_up = c + top - trans * (c + bottom);
      source_dn = -c + bottom - trans * (-c + top);
    } else {
      trans = 1.0 - scaled;
      source_up = scaled * 0.5 * (top + bottom);
      source_dn = source_up;
    }
    return;
  }
  const double f = 0.5 * diffusivity * ssa;
  const double gamma1 = diffusivity - f * (1.0 + asym);
  const double gamma2 = f * (1.0 - asym);
  const double k = std::sqrt(std::max((gamma1 - gamma2) * (gamma1 + gamma2), 1.0e-12));
  if (od > 1.0e-3) {
    const double e = std::exp(-k * od), e2 = e * e;
    const double inv = 1.0 / (k + gamma1 + (k - gamma1) * e2);
    ref = gamma2 * (1.0 - e2) * inv;
    trans = 2.0 * k * e * inv;
    const double c = (bottom - top) / (od * (gamma1 + gamma2));
    source_up = c + top - ref * (-c + top) - trans * (c + bottom);
    source_dn = -c + bottom - ref * (c + bottom) - trans * (-c + top);
  } else {
    ref = gamma2 * od;
    trans = (1.0 - k * od) / (1.0 + od * (gamma1 - k));
    source_up = (1.0 - ref - trans) * 0.5 * (top + bottom);
    source_dn = source_up;
  }
}

void make_optics(const Inputs& in, bool cloudy, bool aerosol_scattering,
                 bool cloud_scattering, std::array<std::vector<double>, 4>& out) {
  const auto count = static_cast<std::size_t>(in.ng) * in.nlev * in.ncol;
  for (auto& v : out) v.resize(count);
  for (int col = 0; col < in.ncol; ++col) {
    for (int lev = 0; lev < in.nlev; ++lev) {
      for (int g = 0; g < in.ng; ++g) {
        const auto at = layer(g, lev, col, in.ng, in.nlev);
        double od = in.od[at];
        double ssa = aerosol_scattering ? in.ssa[at] : 0.0;
        double asym = aerosol_scattering ? in.asymmetry[at] : 0.0;
        const bool cloud_layer = cloudy &&
            in.cloud_fraction[col + static_cast<std::size_t>(in.ncol) * lev] >= in.threshold;
        if (cloud_layer) {
          const int band = in.band_from_g[g] - 1;
          const auto ca = band + static_cast<std::size_t>(in.nbands) *
              (lev + static_cast<std::size_t>(in.nlev) * col);
          const double cod = in.od_scaling[at] * in.od_cloud[ca];
          const double total_od = od + cod;
          if (cloud_scattering) {
            const double aerosol_sod = aerosol_scattering ? ssa * od : 0.0;
            const double cloud_sod = in.ssa_cloud[ca] * cod;
            const double total_sod = aerosol_sod + cloud_sod;
            ssa = total_od > 0.0 ? total_sod / total_od : 0.0;
            asym = total_sod > 0.0
                ? ((aerosol_scattering ? in.asymmetry[at] * aerosol_sod : 0.0)
                   + in.asymmetry_cloud[ca] * cloud_sod) / total_sod : 0.0;
          } else {
            ssa = 0.0; asym = 0.0;
          }
          od = total_od;
        }
        coefficients(od, ssa, asym, in.planck[interface_at(g, lev, col, in.ng, in.nlev)],
                     in.planck[interface_at(g, lev + 1, col, in.ng, in.nlev)],
                     aerosol_scattering || (cloud_layer && cloud_scattering),
                     out[0][at], out[1][at], out[2][at], out[3][at]);
      }
    }
  }
}

void adding(const Inputs& in, bool scattering,
            const std::array<std::vector<double>, 4>& c,
            std::array<std::vector<double>, 2>& flux) {
  const auto count = static_cast<std::size_t>(in.ng) * (in.nlev + 1) * in.ncol;
  flux[0].assign(count, 0.0); flux[1].assign(count, 0.0);
  std::vector<double> alb(in.nlev + 1), source(in.nlev + 1), inverse(in.nlev);
  for (int col = 0; col < in.ncol; ++col) for (int g = 0; g < in.ng; ++g) {
    const auto gp = g + static_cast<std::size_t>(in.ng) * col;
    if (!scattering) {
      for (int lev = 0; lev < in.nlev; ++lev) {
        const auto above = interface_at(g, lev, col, in.ng, in.nlev);
        flux[1][above + in.ng] = c[1][layer(g, lev, col, in.ng, in.nlev)] * flux[1][above]
            + c[3][layer(g, lev, col, in.ng, in.nlev)];
      }
      const auto surface = interface_at(g, in.nlev, col, in.ng, in.nlev);
      flux[0][surface] = in.emission[gp] + in.albedo[gp] * flux[1][surface];
      for (int lev = in.nlev - 1; lev >= 0; --lev) {
        const auto above = interface_at(g, lev, col, in.ng, in.nlev);
        flux[0][above] = c[1][layer(g, lev, col, in.ng, in.nlev)] * flux[0][above + in.ng]
            + c[2][layer(g, lev, col, in.ng, in.nlev)];
      }
      continue;
    }
    alb[in.nlev] = in.albedo[gp]; source[in.nlev] = in.emission[gp];
    for (int lev = in.nlev - 1; lev >= 0; --lev) {
      const auto at = layer(g, lev, col, in.ng, in.nlev);
      inverse[lev] = 1.0 / (1.0 - alb[lev + 1] * c[0][at]);
      alb[lev] = c[0][at] + c[1][at] * c[1][at] * alb[lev + 1] * inverse[lev];
      source[lev] = c[2][at] + c[1][at] *
          (source[lev + 1] + alb[lev + 1] * c[3][at]) * inverse[lev];
    }
    flux[0][interface_at(g, 0, col, in.ng, in.nlev)] = source[0];
    for (int lev = 0; lev < in.nlev; ++lev) {
      const auto at = layer(g, lev, col, in.ng, in.nlev);
      const auto above = interface_at(g, lev, col, in.ng, in.nlev);
      const auto below = above + in.ng;
      flux[1][below] = (c[1][at] * flux[1][above] + c[0][at] * source[lev + 1]
          + c[3][at]) * inverse[lev];
      flux[0][below] = alb[lev + 1] * flux[1][below] + source[lev + 1];
    }
  }
}

Outputs reference(const Inputs& in, bool aerosol_scattering, bool cloud_scattering) {
  std::array<std::vector<double>, 4> clear_c, cloud_c;
  std::array<std::vector<double>, 2> clear_f, cloud_f;
  make_optics(in, false, aerosol_scattering, cloud_scattering, clear_c);
  make_optics(in, true, aerosol_scattering, cloud_scattering, cloud_c);
  adding(in, aerosol_scattering, clear_c, clear_f);
  adding(in, aerosol_scattering || cloud_scattering, cloud_c, cloud_f);
  Outputs out(in.ng, in.nlev, in.ncol);
  std::vector<double> derivative_g(in.ng);
  for (int col = 0; col < in.ncol; ++col) {
    const double cloud = in.total_cloud_cover[col];
    const bool has_cloud = cloud >= in.threshold;
    for (int lev = 0; lev <= in.nlev; ++lev) {
      double up0 = 0.0, dn0 = 0.0, up1 = 0.0, dn1 = 0.0;
      for (int g = 0; g < in.ng; ++g) {
        const auto at = interface_at(g, lev, col, in.ng, in.nlev);
        up0 += clear_f[0][at]; dn0 += clear_f[1][at];
        up1 += cloud_f[0][at]; dn1 += cloud_f[1][at];
      }
      const auto p = col + static_cast<std::size_t>(in.ncol) * lev;
      out.profile[0][p] = up0; out.profile[1][p] = dn0;
      out.profile[2][p] = has_cloud ? cloud * up1 + (1.0 - cloud) * up0 : up0;
      out.profile[3][p] = has_cloud ? cloud * dn1 + (1.0 - cloud) * dn0 : dn0;
    }
    for (int g = 0; g < in.ng; ++g) {
      const auto gp = g + static_cast<std::size_t>(in.ng) * col;
      const auto at = interface_at(g, in.nlev, col, in.ng, in.nlev);
      out.surface[0][gp] = clear_f[1][at];
      out.surface[1][gp] = has_cloud
          ? cloud * cloud_f[1][at] + (1.0 - cloud) * clear_f[1][at] : clear_f[1][at];
    }
    const auto& trans = has_cloud ? cloud_c[1] : clear_c[1];
    const auto& up = has_cloud ? cloud_f[0] : clear_f[0];
    double total = 0.0;
    for (int g = 0; g < in.ng; ++g) total += up[interface_at(g, in.nlev, col, in.ng, in.nlev)];
    for (int g = 0; g < in.ng; ++g)
      derivative_g[g] = up[interface_at(g, in.nlev, col, in.ng, in.nlev)] / total;
    out.derivative[col + static_cast<std::size_t>(in.ncol) * in.nlev] = 1.0;
    for (int lev = in.nlev - 1; lev >= 0; --lev) {
      double sum = 0.0;
      for (int g = 0; g < in.ng; ++g) {
        derivative_g[g] *= trans[layer(g, lev, col, in.ng, in.nlev)];
        sum += derivative_g[g];
      }
      out.derivative[col + static_cast<std::size_t>(in.ncol) * lev] = sum;
    }
    if (has_cloud && cloud < 1.0 - in.threshold) {
      total = 0.0;
      for (int g = 0; g < in.ng; ++g) total += clear_f[0][interface_at(g, in.nlev, col, in.ng, in.nlev)];
      for (int g = 0; g < in.ng; ++g)
        derivative_g[g] = clear_f[0][interface_at(g, in.nlev, col, in.ng, in.nlev)] / total;
      for (int lev = in.nlev - 1; lev >= 0; --lev) {
        double sum = 0.0;
        for (int g = 0; g < in.ng; ++g) {
          derivative_g[g] *= clear_c[1][layer(g, lev, col, in.ng, in.nlev)];
          sum += derivative_g[g];
        }
        const auto p = col + static_cast<std::size_t>(in.ncol) * lev;
        out.derivative[p] = cloud * out.derivative[p] + (1.0 - cloud) * sum;
      }
    }
  }
  return out;
}

double compare(const std::vector<double>& expected, const std::vector<double>& actual,
               const char* name) {
  double worst = 0.0; std::size_t where = 0;
  for (std::size_t i = 0; i < expected.size(); ++i) {
    const double ratio = std::abs(expected[i] - actual[i]) /
        (1.0e-8 + 2.0e-12 * std::max(std::abs(expected[i]), std::abs(actual[i])));
    if (ratio > worst) { worst = ratio; where = i; }
  }
  if (worst > 1.0) std::fprintf(stderr,
      "%s differs at %zu: cpu=%.17g gpu=%.17g normalized_error=%g\n",
      name, where, expected[where], actual[where], worst);
  return worst;
}

int run_case(const Inputs& in, bool aerosol, bool cloud) {
  Outputs gpu(in.ng, in.nlev, in.ncol);
  const int status = oifs_cuda_lw_compute_dp(
      in.ng, in.nbands, in.nlev, in.ncol, aerosol, cloud, true,
      in.od.data(), in.ssa.data(), in.asymmetry.data(), in.planck.data(),
      in.emission.data(), in.albedo.data(), in.band_from_g.data(), in.threshold,
      in.cloud_fraction.data(), in.total_cloud_cover.data(), in.od_scaling.data(),
      in.od_cloud.data(), in.ssa_cloud.data(), in.asymmetry_cloud.data(),
      gpu.profile[0].data(), gpu.profile[1].data(), gpu.profile[2].data(),
      gpu.profile[3].data(), gpu.surface[0].data(), gpu.surface[1].data(),
      gpu.derivative.data());
  if (status != 0) {
    std::fprintf(stderr, "CUDA longwave failed (%d): %s\n", status,
                 oifs_cuda_radiation_last_error());
    return 1;
  }
  const Outputs cpu = reference(in, aerosol, cloud);
  const char* profiles[] = {"lw_up_clear", "lw_dn_clear", "lw_up", "lw_dn"};
  double worst = 0.0;
  for (int i = 0; i < 4; ++i)
    worst = std::max(worst, compare(cpu.profile[i], gpu.profile[i], profiles[i]));
  worst = std::max(worst, compare(cpu.surface[0], gpu.surface[0], "lw_dn_surf_clear_g"));
  worst = std::max(worst, compare(cpu.surface[1], gpu.surface[1], "lw_dn_surf_g"));
  worst = std::max(worst, compare(cpu.derivative, gpu.derivative, "lw_derivatives"));
  if (worst > 1.0) return 2;
  std::printf("CUDA longwave parity passed (aerosol=%d cloud=%d): normalized max error=%g\n",
              aerosol, cloud, worst);
  return 0;
}

}  // namespace

int main() {
  if (!oifs_cuda_radiation_available()) {
    std::fprintf(stderr, "No CUDA device available: %s\n",
                 oifs_cuda_radiation_last_error());
    return 77;
  }
  const Inputs inputs;
  int status = 0;
  for (int aerosol = 0; aerosol <= 1 && status == 0; ++aerosol)
    for (int cloud = 0; cloud <= 1 && status == 0; ++cloud)
      status = run_case(inputs, aerosol != 0, cloud != 0);
  oifs_cuda_radiation_finalize();
  return status;
}
