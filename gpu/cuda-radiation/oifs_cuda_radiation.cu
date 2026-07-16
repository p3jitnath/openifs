#include "oifs_cuda_radiation.h"

#include <cuda_runtime.h>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <mutex>
#include <sstream>
#include <string>

namespace {

struct Workspace {
  int ng = 0;
  int nbands = 0;
  int nlev = 0;
  int ncol = 0;
  cudaStream_t stream = nullptr;

  double *mu0 = nullptr, *od = nullptr, *ssa = nullptr, *asymmetry = nullptr;
  double *albedo_direct = nullptr, *albedo_diffuse = nullptr, *incoming_sw = nullptr;
  int* band_from_g = nullptr;
  double *cloud_fraction = nullptr, *total_cloud_cover = nullptr;
  double* od_scaling = nullptr;
  double *od_cloud = nullptr, *ssa_cloud = nullptr, *asymmetry_cloud = nullptr;
  double *ref_diff = nullptr, *trans_diff = nullptr, *ref_dir = nullptr;
  double *trans_dir_diff = nullptr, *trans_dir_dir = nullptr;
  double *albedo = nullptr, *source = nullptr, *inv_denominator = nullptr;
  double *flux_up = nullptr, *flux_dn_diffuse = nullptr, *flux_dn_direct = nullptr;
  double *flux_up_cloud = nullptr, *flux_dn_diffuse_cloud = nullptr;
  double* flux_dn_direct_cloud = nullptr;
  double *sw_up_clear = nullptr, *sw_dn_clear = nullptr;
  double *sw_dn_direct_clear = nullptr, *sw_up = nullptr, *sw_dn = nullptr;
  double* sw_dn_direct = nullptr;
  double *sw_dn_diffuse_surf_clear_g = nullptr;
  double *sw_dn_direct_surf_clear_g = nullptr;
  double *sw_dn_diffuse_surf_g = nullptr, *sw_dn_direct_surf_g = nullptr;
};

Workspace workspace;

struct LongwaveWorkspace {
  int ng = 0;
  int nbands = 0;
  int nlev = 0;
  int ncol = 0;
  cudaStream_t stream = nullptr;

  double *od = nullptr, *ssa = nullptr, *asymmetry = nullptr;
  double *planck_hl = nullptr, *emission = nullptr, *albedo_surface = nullptr;
  int* band_from_g = nullptr;
  double *cloud_fraction = nullptr, *total_cloud_cover = nullptr;
  double *od_scaling = nullptr, *od_cloud = nullptr;
  double *ssa_cloud = nullptr, *asymmetry_cloud = nullptr;
  double *ref_clear = nullptr, *trans_clear = nullptr;
  double *source_up_clear = nullptr, *source_dn_clear = nullptr;
  double *ref_cloud = nullptr, *trans_cloud = nullptr;
  double *source_up_cloud = nullptr, *source_dn_cloud = nullptr;
  double *albedo = nullptr, *source = nullptr, *inv_denominator = nullptr;
  double *flux_up_clear = nullptr, *flux_dn_clear = nullptr;
  double *flux_up_cloud = nullptr, *flux_dn_cloud = nullptr;
  double *lw_up_clear = nullptr, *lw_dn_clear = nullptr;
  double *lw_up = nullptr, *lw_dn = nullptr;
  double *lw_dn_surf_clear_g = nullptr, *lw_dn_surf_g = nullptr;
  double *derivative_g = nullptr, *lw_derivatives = nullptr;
};

LongwaveWorkspace longwave_workspace;

struct CloudWorkspace {
  int ng = 0;
  int nlev = 0;
  int ncol = 0;
  int pdf_ncdf = 0;
  int pdf_nfsd = 0;
  cudaStream_t stream = nullptr;

  int* seeds = nullptr;
  double *active = nullptr, *cloud_fraction = nullptr;
  double *overlap_parameter = nullptr, *fractional_std = nullptr;
  double* pdf_values = nullptr;
  double *od_scaling = nullptr, *total_cloud_cover = nullptr;
  double *cumulative_cover = nullptr, *pair_cover = nullptr;
  double *overhang = nullptr, *overlap_inhom = nullptr;
  double* random_top = nullptr;
  double *random_cloud = nullptr, *random_inhom1 = nullptr;
  double* random_inhom2 = nullptr;
  std::uint32_t* random_state = nullptr;
};

CloudWorkspace cloud_workspaces[2];
CloudWorkspace* ready_cloud_workspace = nullptr;
std::mutex workspace_mutex;
std::string last_error;
int selected_device = -1;

void free_ptr(void* ptr) {
  if (ptr) cudaFree(ptr);
}

void release_workspace() {
  free_ptr(workspace.mu0);
  free_ptr(workspace.od);
  free_ptr(workspace.ssa);
  free_ptr(workspace.asymmetry);
  free_ptr(workspace.albedo_direct);
  free_ptr(workspace.albedo_diffuse);
  free_ptr(workspace.incoming_sw);
  free_ptr(workspace.band_from_g);
  free_ptr(workspace.cloud_fraction);
  free_ptr(workspace.total_cloud_cover);
  free_ptr(workspace.od_scaling);
  free_ptr(workspace.od_cloud);
  free_ptr(workspace.ssa_cloud);
  free_ptr(workspace.asymmetry_cloud);
  free_ptr(workspace.ref_diff);
  free_ptr(workspace.trans_diff);
  free_ptr(workspace.ref_dir);
  free_ptr(workspace.trans_dir_diff);
  free_ptr(workspace.trans_dir_dir);
  free_ptr(workspace.albedo);
  free_ptr(workspace.source);
  free_ptr(workspace.inv_denominator);
  free_ptr(workspace.flux_up);
  free_ptr(workspace.flux_dn_diffuse);
  free_ptr(workspace.flux_dn_direct);
  free_ptr(workspace.flux_up_cloud);
  free_ptr(workspace.flux_dn_diffuse_cloud);
  free_ptr(workspace.flux_dn_direct_cloud);
  free_ptr(workspace.sw_up_clear);
  free_ptr(workspace.sw_dn_clear);
  free_ptr(workspace.sw_dn_direct_clear);
  free_ptr(workspace.sw_up);
  free_ptr(workspace.sw_dn);
  free_ptr(workspace.sw_dn_direct);
  free_ptr(workspace.sw_dn_diffuse_surf_clear_g);
  free_ptr(workspace.sw_dn_direct_surf_clear_g);
  free_ptr(workspace.sw_dn_diffuse_surf_g);
  free_ptr(workspace.sw_dn_direct_surf_g);
  if (workspace.stream) cudaStreamDestroy(workspace.stream);
  workspace = Workspace{};
}

void release_longwave_workspace() {
#define FREE_LW(member) free_ptr(longwave_workspace.member)
  FREE_LW(od); FREE_LW(ssa); FREE_LW(asymmetry); FREE_LW(planck_hl);
  FREE_LW(emission); FREE_LW(albedo_surface); FREE_LW(band_from_g);
  FREE_LW(cloud_fraction); FREE_LW(total_cloud_cover); FREE_LW(od_scaling);
  FREE_LW(od_cloud); FREE_LW(ssa_cloud); FREE_LW(asymmetry_cloud);
  FREE_LW(ref_clear); FREE_LW(trans_clear); FREE_LW(source_up_clear);
  FREE_LW(source_dn_clear); FREE_LW(ref_cloud); FREE_LW(trans_cloud);
  FREE_LW(source_up_cloud); FREE_LW(source_dn_cloud); FREE_LW(albedo);
  FREE_LW(source); FREE_LW(inv_denominator); FREE_LW(flux_up_clear);
  FREE_LW(flux_dn_clear); FREE_LW(flux_up_cloud); FREE_LW(flux_dn_cloud);
  FREE_LW(lw_up_clear); FREE_LW(lw_dn_clear); FREE_LW(lw_up); FREE_LW(lw_dn);
  FREE_LW(lw_dn_surf_clear_g); FREE_LW(lw_dn_surf_g);
  FREE_LW(derivative_g); FREE_LW(lw_derivatives);
#undef FREE_LW
  if (longwave_workspace.stream) cudaStreamDestroy(longwave_workspace.stream);
  longwave_workspace = LongwaveWorkspace{};
}

void release_cloud_workspace(CloudWorkspace& cloud_workspace) {
  if (ready_cloud_workspace == &cloud_workspace) ready_cloud_workspace = nullptr;
#define FREE_CLOUD(member) free_ptr(cloud_workspace.member)
  FREE_CLOUD(seeds); FREE_CLOUD(active); FREE_CLOUD(cloud_fraction);
  FREE_CLOUD(overlap_parameter); FREE_CLOUD(fractional_std);
  FREE_CLOUD(pdf_values); FREE_CLOUD(od_scaling);
  FREE_CLOUD(total_cloud_cover); FREE_CLOUD(cumulative_cover);
  FREE_CLOUD(pair_cover); FREE_CLOUD(overhang); FREE_CLOUD(overlap_inhom);
  FREE_CLOUD(random_top); FREE_CLOUD(random_cloud); FREE_CLOUD(random_inhom1);
  FREE_CLOUD(random_inhom2); FREE_CLOUD(random_state);
#undef FREE_CLOUD
  if (cloud_workspace.stream) cudaStreamDestroy(cloud_workspace.stream);
  cloud_workspace = CloudWorkspace{};
}

void release_cloud_workspaces() {
  release_cloud_workspace(cloud_workspaces[0]);
  release_cloud_workspace(cloud_workspaces[1]);
}

bool cuda_ok(cudaError_t status, const char* operation) {
  if (status == cudaSuccess) return true;
  std::ostringstream message;
  message << operation << ": " << cudaGetErrorString(status);
  last_error = message.str();
  return false;
}

bool select_device() {
  if (selected_device >= 0) return true;

  int count = 0;
  if (!cuda_ok(cudaGetDeviceCount(&count), "cudaGetDeviceCount") || count == 0) {
    if (count == 0) last_error = "no CUDA device is visible";
    return false;
  }

  int device = 0;
  if (const char* value = std::getenv("OIFS_CUDA_DEVICE")) {
    char* end = nullptr;
    const long parsed = std::strtol(value, &end, 10);
    if (end == value || *end != '\0' || parsed < 0 || parsed >= count) {
      last_error = "OIFS_CUDA_DEVICE is not a visible CUDA device index";
      return false;
    }
    device = static_cast<int>(parsed);
  }
  if (!cuda_ok(cudaSetDevice(device), "cudaSetDevice")) return false;
  selected_device = device;
  return true;
}

template <typename T>
bool allocate(T*& ptr, std::size_t count, const char* name) {
  if (cuda_ok(cudaMalloc(reinterpret_cast<void**>(&ptr), count * sizeof(T)), name)) return true;
  return false;
}

bool ensure_workspace(int ng, int nbands, int nlev, int ncol) {
  if (workspace.ng == ng && workspace.nbands == nbands &&
      workspace.nlev == nlev && workspace.ncol >= ncol) return true;

  release_workspace();
  workspace.ng = ng;
  workspace.nbands = nbands;
  workspace.nlev = nlev;
  workspace.ncol = ncol;
  const std::size_t layer = static_cast<std::size_t>(ng) * nlev * ncol;
  const std::size_t interface = static_cast<std::size_t>(ng) * (nlev + 1) * ncol;
  const std::size_t gpcol = static_cast<std::size_t>(ng) * ncol;
  const std::size_t cloud_layer = static_cast<std::size_t>(nbands) * nlev * ncol;
  const std::size_t fraction = static_cast<std::size_t>(ncol) * nlev;
  const std::size_t profile = static_cast<std::size_t>(ncol) * (nlev + 1);

  if (!cuda_ok(cudaStreamCreateWithFlags(&workspace.stream, cudaStreamNonBlocking),
               "cudaStreamCreate")) {
    release_workspace();
    return false;
  }
#define ALLOCATE(member, count)                                                   \
  do {                                                                            \
    if (!allocate(workspace.member, count, "cudaMalloc(" #member ")")) {         \
      release_workspace();                                                        \
      return false;                                                               \
    }                                                                             \
  } while (false)
  ALLOCATE(mu0, ncol);
  ALLOCATE(od, layer); ALLOCATE(ssa, layer); ALLOCATE(asymmetry, layer);
  ALLOCATE(albedo_direct, gpcol); ALLOCATE(albedo_diffuse, gpcol); ALLOCATE(incoming_sw, gpcol);
  ALLOCATE(band_from_g, ng);
  ALLOCATE(cloud_fraction, fraction); ALLOCATE(total_cloud_cover, ncol);
  ALLOCATE(od_scaling, layer);
  ALLOCATE(od_cloud, cloud_layer); ALLOCATE(ssa_cloud, cloud_layer); ALLOCATE(asymmetry_cloud, cloud_layer);
  ALLOCATE(ref_diff, layer); ALLOCATE(trans_diff, layer); ALLOCATE(ref_dir, layer);
  ALLOCATE(trans_dir_diff, layer); ALLOCATE(trans_dir_dir, layer);
  ALLOCATE(albedo, interface); ALLOCATE(source, interface); ALLOCATE(inv_denominator, layer);
  ALLOCATE(flux_up, interface); ALLOCATE(flux_dn_diffuse, interface); ALLOCATE(flux_dn_direct, interface);
  ALLOCATE(flux_up_cloud, interface); ALLOCATE(flux_dn_diffuse_cloud, interface);
  ALLOCATE(flux_dn_direct_cloud, interface);
  ALLOCATE(sw_up_clear, profile); ALLOCATE(sw_dn_clear, profile);
  ALLOCATE(sw_dn_direct_clear, profile); ALLOCATE(sw_up, profile);
  ALLOCATE(sw_dn, profile); ALLOCATE(sw_dn_direct, profile);
  ALLOCATE(sw_dn_diffuse_surf_clear_g, gpcol);
  ALLOCATE(sw_dn_direct_surf_clear_g, gpcol);
  ALLOCATE(sw_dn_diffuse_surf_g, gpcol); ALLOCATE(sw_dn_direct_surf_g, gpcol);
#undef ALLOCATE
  return true;
}

bool ensure_longwave_workspace(int ng, int nbands, int nlev, int ncol) {
  if (longwave_workspace.ng == ng && longwave_workspace.nbands == nbands &&
      longwave_workspace.nlev == nlev && longwave_workspace.ncol >= ncol) return true;

  release_longwave_workspace();
  longwave_workspace.ng = ng;
  longwave_workspace.nbands = nbands;
  longwave_workspace.nlev = nlev;
  longwave_workspace.ncol = ncol;
  const std::size_t layer = static_cast<std::size_t>(ng) * nlev * ncol;
  const std::size_t interface = static_cast<std::size_t>(ng) * (nlev + 1) * ncol;
  const std::size_t gpcol = static_cast<std::size_t>(ng) * ncol;
  const std::size_t cloud_layer = static_cast<std::size_t>(nbands) * nlev * ncol;
  const std::size_t fraction = static_cast<std::size_t>(ncol) * nlev;
  const std::size_t profile = static_cast<std::size_t>(ncol) * (nlev + 1);

  if (!cuda_ok(cudaStreamCreateWithFlags(&longwave_workspace.stream, cudaStreamNonBlocking),
               "cudaStreamCreate(longwave)")) {
    release_longwave_workspace();
    return false;
  }
#define ALLOCATE_LW(member, count)                                                \
  do {                                                                            \
    if (!allocate(longwave_workspace.member, count, "cudaMalloc(lw " #member ")")) { \
      release_longwave_workspace();                                               \
      return false;                                                               \
    }                                                                             \
  } while (false)
  ALLOCATE_LW(od, layer); ALLOCATE_LW(ssa, layer); ALLOCATE_LW(asymmetry, layer);
  ALLOCATE_LW(planck_hl, interface); ALLOCATE_LW(emission, gpcol);
  ALLOCATE_LW(albedo_surface, gpcol); ALLOCATE_LW(band_from_g, ng);
  ALLOCATE_LW(cloud_fraction, fraction); ALLOCATE_LW(total_cloud_cover, ncol);
  ALLOCATE_LW(od_scaling, layer); ALLOCATE_LW(od_cloud, cloud_layer);
  ALLOCATE_LW(ssa_cloud, cloud_layer); ALLOCATE_LW(asymmetry_cloud, cloud_layer);
  ALLOCATE_LW(ref_clear, layer); ALLOCATE_LW(trans_clear, layer);
  ALLOCATE_LW(source_up_clear, layer); ALLOCATE_LW(source_dn_clear, layer);
  ALLOCATE_LW(ref_cloud, layer); ALLOCATE_LW(trans_cloud, layer);
  ALLOCATE_LW(source_up_cloud, layer); ALLOCATE_LW(source_dn_cloud, layer);
  ALLOCATE_LW(albedo, interface); ALLOCATE_LW(source, interface);
  ALLOCATE_LW(inv_denominator, layer); ALLOCATE_LW(flux_up_clear, interface);
  ALLOCATE_LW(flux_dn_clear, interface); ALLOCATE_LW(flux_up_cloud, interface);
  ALLOCATE_LW(flux_dn_cloud, interface); ALLOCATE_LW(lw_up_clear, profile);
  ALLOCATE_LW(lw_dn_clear, profile); ALLOCATE_LW(lw_up, profile);
  ALLOCATE_LW(lw_dn, profile); ALLOCATE_LW(lw_dn_surf_clear_g, gpcol);
  ALLOCATE_LW(lw_dn_surf_g, gpcol); ALLOCATE_LW(derivative_g, gpcol);
  ALLOCATE_LW(lw_derivatives, profile);
#undef ALLOCATE_LW
  return true;
}

bool ensure_cloud_workspace(
    int ng, int nlev, int ncol, int pdf_ncdf, int pdf_nfsd,
    CloudWorkspace*& selected) {
  for (auto& candidate : cloud_workspaces) {
    if (candidate.ng == ng && candidate.nlev == nlev &&
        candidate.ncol >= ncol && candidate.pdf_ncdf == pdf_ncdf &&
        candidate.pdf_nfsd == pdf_nfsd) {
      selected = &candidate;
      return true;
    }
  }

  selected = nullptr;
  for (auto& candidate : cloud_workspaces) {
    if (!candidate.stream) {
      selected = &candidate;
      break;
    }
  }
  if (!selected) selected = &cloud_workspaces[0];
  release_cloud_workspace(*selected);
  auto& cloud_workspace = *selected;
  cloud_workspace.ng = ng;
  cloud_workspace.nlev = nlev;
  cloud_workspace.ncol = ncol;
  cloud_workspace.pdf_ncdf = pdf_ncdf;
  cloud_workspace.pdf_nfsd = pdf_nfsd;
  const std::size_t profile = static_cast<std::size_t>(ncol) * nlev;
  const std::size_t interfaces = static_cast<std::size_t>(ncol) * (nlev - 1);
  const std::size_t scaling = static_cast<std::size_t>(ng) * nlev * ncol;
  const std::size_t subcolumns = static_cast<std::size_t>(ng) * ncol;
  const std::size_t random_state = static_cast<std::size_t>(607) * ncol;
  const std::size_t pdf_size = static_cast<std::size_t>(pdf_ncdf) * pdf_nfsd;

  if (!cuda_ok(cudaStreamCreateWithFlags(&cloud_workspace.stream, cudaStreamNonBlocking),
               "cudaStreamCreate(cloud)")) {
    release_cloud_workspace(cloud_workspace);
    return false;
  }
#define ALLOCATE_CLOUD(member, count)                                             \
  do {                                                                            \
    if (!allocate(cloud_workspace.member, count, "cudaMalloc(cloud " #member ")")) { \
      release_cloud_workspace(cloud_workspace);                                   \
      return false;                                                               \
    }                                                                             \
  } while (false)
  ALLOCATE_CLOUD(seeds, ncol); ALLOCATE_CLOUD(active, ncol);
  ALLOCATE_CLOUD(cloud_fraction, profile);
  ALLOCATE_CLOUD(overlap_parameter, interfaces);
  ALLOCATE_CLOUD(fractional_std, profile); ALLOCATE_CLOUD(pdf_values, pdf_size);
  ALLOCATE_CLOUD(od_scaling, scaling); ALLOCATE_CLOUD(total_cloud_cover, ncol);
  ALLOCATE_CLOUD(cumulative_cover, profile); ALLOCATE_CLOUD(pair_cover, profile);
  ALLOCATE_CLOUD(overhang, profile); ALLOCATE_CLOUD(overlap_inhom, profile);
  ALLOCATE_CLOUD(random_top, subcolumns); ALLOCATE_CLOUD(random_cloud, profile);
  ALLOCATE_CLOUD(random_inhom1, profile);
  ALLOCATE_CLOUD(random_inhom2, profile); ALLOCATE_CLOUD(random_state, random_state);
#undef ALLOCATE_CLOUD
  return true;
}

__device__ inline std::size_t layer_index(int g, int lev, int col, int ng, int nlev) {
  return g + static_cast<std::size_t>(ng) * (lev + static_cast<std::size_t>(nlev) * col);
}

__device__ inline std::size_t interface_index(int g, int lev, int col, int ng, int nlev) {
  return g + static_cast<std::size_t>(ng) * (lev + static_cast<std::size_t>(nlev + 1) * col);
}

__global__ void optics_kernel(
    int ng, int nbands, int nlev, int ncol, bool cloudy, bool do_delta_scaling,
    double cloud_fraction_threshold, const double* mu0, const double* od,
    const double* ssa, const double* asymmetry, const int* band_from_g,
    const double* cloud_fraction, const double* od_scaling,
    const double* od_cloud, const double* ssa_cloud, const double* asymmetry_cloud,
    double* ref_diff, double* trans_diff, double* ref_dir,
    double* trans_dir_diff, double* trans_dir_dir) {
  const std::size_t linear = blockIdx.x * static_cast<std::size_t>(blockDim.x) + threadIdx.x;
  const std::size_t count = static_cast<std::size_t>(ng) * nlev * ncol;
  if (linear >= count) return;
  const int g = linear % ng;
  const int lev = (linear / ng) % nlev;
  const int col = linear / (static_cast<std::size_t>(ng) * nlev);
  const double cosine = mu0[col];
  if (cosine <= 0.0) {
    ref_diff[linear] = trans_diff[linear] = ref_dir[linear] = 0.0;
    trans_dir_diff[linear] = trans_dir_dir[linear] = 0.0;
    return;
  }

  double optical_depth = od[linear];
  double single_scattering_albedo = ssa[linear];
  double asymmetry_factor = asymmetry[linear];
  if (cloudy && cloud_fraction[col + static_cast<std::size_t>(ncol) * lev] >= cloud_fraction_threshold) {
    const int band = band_from_g[g] - 1;
    const std::size_t cloud = band + static_cast<std::size_t>(nbands) *
        (lev + static_cast<std::size_t>(nlev) * col);
    const double cloud_od = od_scaling[linear] * od_cloud[cloud];
    optical_depth += cloud_od;
    single_scattering_albedo = optical_depth > 0.0
        ? (ssa[linear] * od[linear] + ssa_cloud[cloud] * cloud_od) / optical_depth : 0.0;
    const double scattering_depth = single_scattering_albedo * optical_depth;
    asymmetry_factor = scattering_depth > 0.0
        ? (asymmetry[linear] * ssa[linear] * od[linear]
           + asymmetry_cloud[cloud] * ssa_cloud[cloud] * cloud_od) / scattering_depth : 0.0;
  }
  if (do_delta_scaling) {
    const double forward_lobe = asymmetry_factor * asymmetry_factor;
    optical_depth *= 1.0 - single_scattering_albedo * forward_lobe;
    single_scattering_albedo = single_scattering_albedo * (1.0 - forward_lobe)
        / (1.0 - single_scattering_albedo * forward_lobe);
    asymmetry_factor /= 1.0 + asymmetry_factor;
  }

  const double factor = 0.75 * asymmetry_factor;
  const double gamma1 = 2.0 - single_scattering_albedo * (1.25 + factor);
  const double gamma2 = single_scattering_albedo * (0.75 - factor);
  const double gamma3 = 0.5 - cosine * factor;
  const double gamma4 = 1.0 - gamma3;
  const double alpha1 = gamma1 * gamma4 + gamma2 * gamma3;
  const double alpha2 = gamma1 * gamma3 + gamma2 * gamma4;
  const double k = sqrt(fmax((gamma1 - gamma2) * (gamma1 + gamma2), 1.0e-12));
  double mu = cosine;
  if (fabs(1.0 - k * mu) < 1000.0 * 2.2204460492503131e-16)
    mu *= 1.0 - 1000.0 * 2.2204460492503131e-16;
  const double k_mu = k * mu;
  const double exponential0 = exp(-fmax(optical_depth / mu, 0.0));
  const double exponential = exp(-k * optical_depth);
  const double exponential2 = exponential * exponential;
  const double k2exponential = 2.0 * k * exponential;
  double denominator = 1.0 / (k + gamma1 + (k - gamma1) * exponential2);
  ref_diff[linear] = gamma2 * (1.0 - exponential2) * denominator;
  trans_diff[linear] = k2exponential * denominator;
  trans_dir_dir[linear] = exponential0;
  denominator *= mu * single_scattering_albedo / (1.0 - k_mu * k_mu);
  const double direct_reflectance = denominator *
      ((1.0 - k_mu) * (alpha2 + k * gamma3)
       - (1.0 + k_mu) * (alpha2 - k * gamma3) * exponential2
       - k2exponential * (gamma3 - alpha2 * mu) * exponential0);
  ref_dir[linear] = fmax(0.0, fmin(direct_reflectance, 1.0));
  const double direct_diffuse = denominator *
      (k2exponential * (gamma4 + alpha1 * mu)
       - exponential0 * ((1.0 + k_mu) * (alpha1 + k * gamma4)
       - (1.0 - k_mu) * (alpha1 - k * gamma4) * exponential2));
  trans_dir_diff[linear] = fmax(0.0, fmin(direct_diffuse, 1.0 - ref_dir[linear]));
}

__global__ void adding_kernel(
    int ng, int nlev, int ncol, const double* mu0, const double* incoming_sw,
    const double* albedo_direct_surface, const double* albedo_diffuse_surface,
    const double* ref_diff, const double* trans_diff, const double* ref_dir,
    const double* trans_dir_diff, const double* trans_dir_dir,
    double* albedo, double* source, double* inv_denominator,
    double* flux_up, double* flux_dn_diffuse, double* flux_dn_direct) {
  const std::size_t linear = blockIdx.x * static_cast<std::size_t>(blockDim.x) + threadIdx.x;
  const std::size_t count = static_cast<std::size_t>(ng) * ncol;
  if (linear >= count) return;
  const int g = linear % ng;
  const int col = linear / ng;
  const double cosine = mu0[col];
  if (cosine <= 0.0) {
    for (int lev = 0; lev <= nlev; ++lev) {
      const auto at = interface_index(g, lev, col, ng, nlev);
      flux_up[at] = flux_dn_diffuse[at] = flux_dn_direct[at] = 0.0;
    }
    return;
  }

  const std::size_t gpcol = g + static_cast<std::size_t>(ng) * col;
  auto at0 = interface_index(g, 0, col, ng, nlev);
  flux_dn_direct[at0] = incoming_sw[gpcol];
  for (int lev = 0; lev < nlev; ++lev) {
    const auto layer = layer_index(g, lev, col, ng, nlev);
    const auto above = interface_index(g, lev, col, ng, nlev);
    const auto below = interface_index(g, lev + 1, col, ng, nlev);
    flux_dn_direct[below] = flux_dn_direct[above] * trans_dir_dir[layer];
  }
  auto surface = interface_index(g, nlev, col, ng, nlev);
  albedo[surface] = albedo_diffuse_surface[gpcol];
  source[surface] = albedo_direct_surface[gpcol] * flux_dn_direct[surface] * cosine;
  for (int lev = nlev - 1; lev >= 0; --lev) {
    const auto layer = layer_index(g, lev, col, ng, nlev);
    const auto above = interface_index(g, lev, col, ng, nlev);
    const auto below = interface_index(g, lev + 1, col, ng, nlev);
    inv_denominator[layer] = 1.0 / (1.0 - albedo[below] * ref_diff[layer]);
    albedo[above] = ref_diff[layer] + trans_diff[layer] * trans_diff[layer]
        * albedo[below] * inv_denominator[layer];
    source[above] = ref_dir[layer] * flux_dn_direct[above]
        + trans_diff[layer] * (source[below] + albedo[below] * trans_dir_diff[layer]
        * flux_dn_direct[above]) * inv_denominator[layer];
  }
  flux_dn_diffuse[at0] = 0.0;
  flux_up[at0] = source[at0];
  for (int lev = 0; lev < nlev; ++lev) {
    const auto layer = layer_index(g, lev, col, ng, nlev);
    const auto above = interface_index(g, lev, col, ng, nlev);
    const auto below = interface_index(g, lev + 1, col, ng, nlev);
    flux_dn_diffuse[below] = (trans_diff[layer] * flux_dn_diffuse[above]
        + ref_diff[layer] * source[below] + trans_dir_diff[layer]
        * flux_dn_direct[above]) * inv_denominator[layer];
    flux_up[below] = albedo[below] * flux_dn_diffuse[below] + source[below];
    flux_dn_direct[above] *= cosine;
  }
  flux_dn_direct[surface] *= cosine;
}

__global__ void reduce_profiles_kernel(
    int ng, int nlev, int ncol, const double* total_cloud_cover,
    const double* flux_up_clear_g, const double* flux_dn_diffuse_clear_g,
    const double* flux_dn_direct_clear_g, const double* flux_up_cloud_g,
    const double* flux_dn_diffuse_cloud_g, const double* flux_dn_direct_cloud_g,
    double* sw_up_clear, double* sw_dn_clear, double* sw_dn_direct_clear,
    double* sw_up, double* sw_dn, double* sw_dn_direct) {
  const std::size_t linear = blockIdx.x * static_cast<std::size_t>(blockDim.x) + threadIdx.x;
  const std::size_t count = static_cast<std::size_t>(ncol) * (nlev + 1);
  if (linear >= count) return;

  const int col = linear % ncol;
  const int lev = linear / ncol;
  double up_clear = 0.0;
  double diffuse_clear = 0.0;
  double direct_clear = 0.0;
  double up_cloud = 0.0;
  double diffuse_cloud = 0.0;
  double direct_cloud = 0.0;
  for (int g = 0; g < ng; ++g) {
    const auto at = interface_index(g, lev, col, ng, nlev);
    up_clear += flux_up_clear_g[at];
    diffuse_clear += flux_dn_diffuse_clear_g[at];
    direct_clear += flux_dn_direct_clear_g[at];
    up_cloud += flux_up_cloud_g[at];
    diffuse_cloud += flux_dn_diffuse_cloud_g[at];
    direct_cloud += flux_dn_direct_cloud_g[at];
  }

  const double cloud = total_cloud_cover[col];
  const double clear = 1.0 - cloud;
  sw_up_clear[linear] = up_clear;
  sw_dn_clear[linear] = diffuse_clear + direct_clear;
  sw_dn_direct_clear[linear] = direct_clear;
  sw_up[linear] = cloud * up_cloud + clear * up_clear;
  sw_dn[linear] = cloud * (diffuse_cloud + direct_cloud)
      + clear * (diffuse_clear + direct_clear);
  sw_dn_direct[linear] = cloud * direct_cloud + clear * direct_clear;
}

__global__ void reduce_surface_kernel(
    int ng, int nlev, int ncol, const double* total_cloud_cover,
    const double* flux_dn_diffuse_clear, const double* flux_dn_direct_clear,
    const double* flux_dn_diffuse_cloud, const double* flux_dn_direct_cloud,
    double* sw_dn_diffuse_surf_clear_g, double* sw_dn_direct_surf_clear_g,
    double* sw_dn_diffuse_surf_g, double* sw_dn_direct_surf_g) {
  const std::size_t linear = blockIdx.x * static_cast<std::size_t>(blockDim.x) + threadIdx.x;
  const std::size_t count = static_cast<std::size_t>(ng) * ncol;
  if (linear >= count) return;

  const int g = linear % ng;
  const int col = linear / ng;
  const auto at = interface_index(g, nlev, col, ng, nlev);
  const double diffuse_clear = flux_dn_diffuse_clear[at];
  const double direct_clear = flux_dn_direct_clear[at];
  const double cloud = total_cloud_cover[col];
  const double clear = 1.0 - cloud;
  sw_dn_diffuse_surf_clear_g[linear] = diffuse_clear;
  sw_dn_direct_surf_clear_g[linear] = direct_clear;
  sw_dn_diffuse_surf_g[linear] = cloud * flux_dn_diffuse_cloud[at]
      + clear * diffuse_clear;
  sw_dn_direct_surf_g[linear] = cloud * flux_dn_direct_cloud[at]
      + clear * direct_clear;
}

__device__ inline void longwave_no_scattering(
    double optical_depth, double planck_top, double planck_bottom,
    double& reflectance, double& transmittance,
    double& source_up, double& source_dn) {
  constexpr double diffusivity = 1.66;
  reflectance = 0.0;
  const double scaled_depth = diffusivity * optical_depth;
  if (optical_depth > 1.0e-3) {
    transmittance = exp(-scaled_depth);
    const double coefficient = (planck_bottom - planck_top) / scaled_depth;
    const double up_top = coefficient + planck_top;
    const double up_bottom = coefficient + planck_bottom;
    const double dn_top = -coefficient + planck_top;
    const double dn_bottom = -coefficient + planck_bottom;
    source_up = up_top - transmittance * up_bottom;
    source_dn = dn_bottom - transmittance * dn_top;
  } else {
    transmittance = 1.0 - scaled_depth;
    source_up = scaled_depth * 0.5 * (planck_top + planck_bottom);
    source_dn = source_up;
  }
}

__device__ inline void longwave_scattering(
    double optical_depth, double single_scattering_albedo,
    double asymmetry_factor, double planck_top, double planck_bottom,
    double& reflectance, double& transmittance,
    double& source_up, double& source_dn) {
  constexpr double diffusivity = 1.66;
  const double factor = (diffusivity * 0.5) * single_scattering_albedo;
  const double gamma1 = diffusivity - factor * (1.0 + asymmetry_factor);
  const double gamma2 = factor * (1.0 - asymmetry_factor);
  const double exponent = sqrt(fmax((gamma1 - gamma2) * (gamma1 + gamma2), 1.0e-12));
  if (optical_depth > 1.0e-3) {
    const double exponential = exp(-exponent * optical_depth);
    const double exponential2 = exponential * exponential;
    const double inverse = 1.0 /
        (exponent + gamma1 + (exponent - gamma1) * exponential2);
    reflectance = gamma2 * (1.0 - exponential2) * inverse;
    transmittance = 2.0 * exponent * exponential * inverse;
    const double coefficient = (planck_bottom - planck_top) /
        (optical_depth * (gamma1 + gamma2));
    const double up_top = coefficient + planck_top;
    const double up_bottom = coefficient + planck_bottom;
    const double dn_top = -coefficient + planck_top;
    const double dn_bottom = -coefficient + planck_bottom;
    source_up = up_top - reflectance * dn_top - transmittance * up_bottom;
    source_dn = dn_bottom - reflectance * up_bottom - transmittance * dn_top;
  } else {
    reflectance = gamma2 * optical_depth;
    transmittance = (1.0 - exponent * optical_depth) /
        (1.0 + optical_depth * (gamma1 - exponent));
    source_up = (1.0 - reflectance - transmittance) *
        0.5 * (planck_top + planck_bottom);
    source_dn = source_up;
  }
}

__global__ void longwave_optics_kernel(
    int ng, int nbands, int nlev, int ncol, bool cloudy,
    bool aerosol_scattering, bool cloud_scattering,
    double cloud_fraction_threshold, const double* od, const double* ssa,
    const double* asymmetry, const double* planck_hl,
    const int* band_from_g, const double* cloud_fraction,
    const double* od_scaling, const double* od_cloud,
    const double* ssa_cloud, const double* asymmetry_cloud,
    double* reflectance, double* transmittance,
    double* source_up, double* source_dn) {
  const std::size_t linear = blockIdx.x * static_cast<std::size_t>(blockDim.x) + threadIdx.x;
  const std::size_t count = static_cast<std::size_t>(ng) * nlev * ncol;
  if (linear >= count) return;
  const int g = linear % ng;
  const int lev = (linear / ng) % nlev;
  const int col = linear / (static_cast<std::size_t>(ng) * nlev);
  const auto top = interface_index(g, lev, col, ng, nlev);
  const auto bottom = interface_index(g, lev + 1, col, ng, nlev);

  double optical_depth = od[linear];
  double single_scattering_albedo = aerosol_scattering ? ssa[linear] : 0.0;
  double asymmetry_factor = aerosol_scattering ? asymmetry[linear] : 0.0;
  const bool cloud_layer = cloudy &&
      cloud_fraction[col + static_cast<std::size_t>(ncol) * lev] >= cloud_fraction_threshold;
  if (cloud_layer) {
    const int band = band_from_g[g] - 1;
    const auto cloud = band + static_cast<std::size_t>(nbands) *
        (lev + static_cast<std::size_t>(nlev) * col);
    const double cloud_depth = od_scaling[linear] * od_cloud[cloud];
    const double total_depth = optical_depth + cloud_depth;
    if (cloud_scattering) {
      const double aerosol_scattering_depth = aerosol_scattering
          ? single_scattering_albedo * optical_depth : 0.0;
      const double cloud_scattering_depth = ssa_cloud[cloud] * cloud_depth;
      const double total_scattering_depth =
          aerosol_scattering_depth + cloud_scattering_depth;
      single_scattering_albedo = total_depth > 0.0
          ? total_scattering_depth / total_depth : 0.0;
      asymmetry_factor = total_scattering_depth > 0.0
          ? ((aerosol_scattering ? asymmetry[linear] * aerosol_scattering_depth : 0.0)
             + asymmetry_cloud[cloud] * cloud_scattering_depth)
                / total_scattering_depth
          : 0.0;
    } else {
      single_scattering_albedo = 0.0;
      asymmetry_factor = 0.0;
    }
    optical_depth = total_depth;
  }

  const bool scattering = aerosol_scattering || (cloud_layer && cloud_scattering);
  if (scattering) {
    longwave_scattering(optical_depth, single_scattering_albedo, asymmetry_factor,
                        planck_hl[top], planck_hl[bottom], reflectance[linear],
                        transmittance[linear], source_up[linear], source_dn[linear]);
  } else {
    longwave_no_scattering(optical_depth, planck_hl[top], planck_hl[bottom],
                           reflectance[linear], transmittance[linear],
                           source_up[linear], source_dn[linear]);
  }
}

__global__ void longwave_adding_kernel(
    int ng, int nlev, int ncol, bool scattering,
    const double* reflectance, const double* transmittance,
    const double* source_up, const double* source_dn,
    const double* emission_surface, const double* albedo_surface,
    double* albedo, double* source, double* inv_denominator,
    double* flux_up, double* flux_dn) {
  const std::size_t linear = blockIdx.x * static_cast<std::size_t>(blockDim.x) + threadIdx.x;
  const std::size_t count = static_cast<std::size_t>(ng) * ncol;
  if (linear >= count) return;
  const int g = linear % ng;
  const int col = linear / ng;
  const auto top = interface_index(g, 0, col, ng, nlev);
  const auto surface = interface_index(g, nlev, col, ng, nlev);

  flux_dn[top] = 0.0;
  if (!scattering) {
    for (int lev = 0; lev < nlev; ++lev) {
      const auto layer = layer_index(g, lev, col, ng, nlev);
      const auto above = interface_index(g, lev, col, ng, nlev);
      const auto below = interface_index(g, lev + 1, col, ng, nlev);
      flux_dn[below] = transmittance[layer] * flux_dn[above] + source_dn[layer];
    }
    flux_up[surface] = emission_surface[linear] + albedo_surface[linear] * flux_dn[surface];
    for (int lev = nlev - 1; lev >= 0; --lev) {
      const auto layer = layer_index(g, lev, col, ng, nlev);
      const auto above = interface_index(g, lev, col, ng, nlev);
      const auto below = interface_index(g, lev + 1, col, ng, nlev);
      flux_up[above] = transmittance[layer] * flux_up[below] + source_up[layer];
    }
    return;
  }

  albedo[surface] = albedo_surface[linear];
  source[surface] = emission_surface[linear];
  for (int lev = nlev - 1; lev >= 0; --lev) {
    const auto layer = layer_index(g, lev, col, ng, nlev);
    const auto above = interface_index(g, lev, col, ng, nlev);
    const auto below = interface_index(g, lev + 1, col, ng, nlev);
    inv_denominator[layer] = 1.0 / (1.0 - albedo[below] * reflectance[layer]);
    albedo[above] = reflectance[layer] + transmittance[layer] * transmittance[layer]
        * albedo[below] * inv_denominator[layer];
    source[above] = source_up[layer] + transmittance[layer]
        * (source[below] + albedo[below] * source_dn[layer]) * inv_denominator[layer];
  }
  flux_up[top] = source[top];
  for (int lev = 0; lev < nlev; ++lev) {
    const auto layer = layer_index(g, lev, col, ng, nlev);
    const auto above = interface_index(g, lev, col, ng, nlev);
    const auto below = interface_index(g, lev + 1, col, ng, nlev);
    flux_dn[below] = (transmittance[layer] * flux_dn[above]
        + reflectance[layer] * source[below] + source_dn[layer])
        * inv_denominator[layer];
    flux_up[below] = albedo[below] * flux_dn[below] + source[below];
  }
}

__global__ void reduce_longwave_kernel(
    int ng, int nlev, int ncol, double cloud_fraction_threshold,
    const double* total_cloud_cover, const double* flux_up_clear_g,
    const double* flux_dn_clear_g, const double* flux_up_cloud_g,
    const double* flux_dn_cloud_g, double* lw_up_clear, double* lw_dn_clear,
    double* lw_up, double* lw_dn) {
  const std::size_t linear = blockIdx.x * static_cast<std::size_t>(blockDim.x) + threadIdx.x;
  const std::size_t count = static_cast<std::size_t>(ncol) * (nlev + 1);
  if (linear >= count) return;
  const int col = linear % ncol;
  const int lev = linear / ncol;
  double up_clear = 0.0, dn_clear = 0.0, up_cloud = 0.0, dn_cloud = 0.0;
  for (int g = 0; g < ng; ++g) {
    const auto at = interface_index(g, lev, col, ng, nlev);
    up_clear += flux_up_clear_g[at];
    dn_clear += flux_dn_clear_g[at];
    up_cloud += flux_up_cloud_g[at];
    dn_cloud += flux_dn_cloud_g[at];
  }
  const double cloud = total_cloud_cover[col];
  lw_up_clear[linear] = up_clear;
  lw_dn_clear[linear] = dn_clear;
  if (cloud >= cloud_fraction_threshold) {
    lw_up[linear] = cloud * up_cloud + (1.0 - cloud) * up_clear;
    lw_dn[linear] = cloud * dn_cloud + (1.0 - cloud) * dn_clear;
  } else {
    lw_up[linear] = up_clear;
    lw_dn[linear] = dn_clear;
  }
}

__global__ void reduce_longwave_surface_kernel(
    int ng, int nlev, int ncol, double cloud_fraction_threshold,
    const double* total_cloud_cover, const double* flux_dn_clear,
    const double* flux_dn_cloud, double* lw_dn_surf_clear_g,
    double* lw_dn_surf_g) {
  const std::size_t linear = blockIdx.x * static_cast<std::size_t>(blockDim.x) + threadIdx.x;
  const std::size_t count = static_cast<std::size_t>(ng) * ncol;
  if (linear >= count) return;
  const int g = linear % ng;
  const int col = linear / ng;
  const auto surface = interface_index(g, nlev, col, ng, nlev);
  const double clear_value = flux_dn_clear[surface];
  const double cloud = total_cloud_cover[col];
  lw_dn_surf_clear_g[linear] = clear_value;
  lw_dn_surf_g[linear] = cloud >= cloud_fraction_threshold
      ? cloud * flux_dn_cloud[surface] + (1.0 - cloud) * clear_value
      : clear_value;
}

__global__ void longwave_derivatives_kernel(
    int ng, int nlev, int ncol, double cloud_fraction_threshold,
    const double* total_cloud_cover, const double* trans_clear,
    const double* trans_cloud, const double* flux_up_clear,
    const double* flux_up_cloud, double* derivative_g,
    double* lw_derivatives) {
  const int col = blockIdx.x * blockDim.x + threadIdx.x;
  if (col >= ncol) return;
  const double cloud = total_cloud_cover[col];
  const bool has_cloud = cloud >= cloud_fraction_threshold;
  const double* trans = has_cloud ? trans_cloud : trans_clear;
  const double* up = has_cloud ? flux_up_cloud : flux_up_clear;
  double total = 0.0;
  for (int g = 0; g < ng; ++g)
    total += up[interface_index(g, nlev, col, ng, nlev)];
  for (int g = 0; g < ng; ++g) {
    const auto gp = g + static_cast<std::size_t>(ng) * col;
    derivative_g[gp] = up[interface_index(g, nlev, col, ng, nlev)] / total;
  }
  lw_derivatives[col + static_cast<std::size_t>(ncol) * nlev] = 1.0;
  for (int lev = nlev - 1; lev >= 0; --lev) {
    double sum = 0.0;
    for (int g = 0; g < ng; ++g) {
      const auto gp = g + static_cast<std::size_t>(ng) * col;
      derivative_g[gp] *= trans[layer_index(g, lev, col, ng, nlev)];
      sum += derivative_g[gp];
    }
    lw_derivatives[col + static_cast<std::size_t>(ncol) * lev] = sum;
  }
  if (has_cloud && cloud < 1.0 - cloud_fraction_threshold) {
    total = 0.0;
    for (int g = 0; g < ng; ++g)
      total += flux_up_clear[interface_index(g, nlev, col, ng, nlev)];
    for (int g = 0; g < ng; ++g) {
      const auto gp = g + static_cast<std::size_t>(ng) * col;
      derivative_g[gp] = flux_up_clear[interface_index(g, nlev, col, ng, nlev)] / total;
    }
    for (int lev = nlev - 1; lev >= 0; --lev) {
      double sum = 0.0;
      for (int g = 0; g < ng; ++g) {
        const auto gp = g + static_cast<std::size_t>(ng) * col;
        derivative_g[gp] *= trans_clear[layer_index(g, lev, col, ng, nlev)];
        sum += derivative_g[gp];
      }
      const auto at = col + static_cast<std::size_t>(ncol) * lev;
      lw_derivatives[at] = cloud * lw_derivatives[at] + (1.0 - cloud) * sum;
    }
  }
}

__device__ inline void cloud_rng_generate(std::uint32_t* state) {
  constexpr std::uint32_t mask = 0x3fffffffU;
  for (int j = 0; j < 273; ++j)
    state[j] = mask & (state[j] + state[j + 334]);
  for (int j = 273; j < 607; ++j)
    state[j] = mask & (state[j] + state[j - 273]);
}

__device__ inline double cloud_rng_next(std::uint32_t* state, int& used) {
  if (used >= 607) {
    cloud_rng_generate(state);
    used = 0;
  }
  return static_cast<double>(state[used++]) * (1.0 / 1073741824.0);
}

__device__ __constant__ std::uint32_t cloud_rng_jump_604[32] = {
  0x29c7fc9aU, 0x538ff934U, 0xa71ff268U, 0x4e3fe47fU,
  0x9c7fc8feU, 0x38ff9153U, 0x71ff22a6U, 0xe3fe454cU,
  0xc7fc8a37U, 0x8ff914c1U, 0x1ff2292dU, 0x3fe4525aU,
  0x7fc8a4b4U, 0xff914968U, 0xff22927fU, 0xfe452451U,
  0xfc8a480dU, 0xf91490b5U, 0xf22921c5U, 0xe4524325U,
  0xc8a486e5U, 0x91490d65U, 0x22921a65U, 0x452434caU,
  0x8a486994U, 0x1490d387U, 0x2921a70eU, 0x52434e1cU,
  0xa4869c38U, 0x490d38dfU, 0x921a71beU, 0x2434e3d3U
};

__device__ inline std::uint32_t cloud_rng_step(std::uint32_t value) {
  return (value << 1U) ^ ((value & 0x80000000U) ? 0xafU : 0U);
}

__device__ inline std::uint32_t cloud_rng_advance_604(std::uint32_t value) {
  std::uint32_t advanced = 0;
#pragma unroll
  for (int bit = 0; bit < 32; ++bit)
    if (value & (1U << bit)) advanced ^= cloud_rng_jump_604[bit];
  return advanced;
}

__device__ void cloud_rng_generate_warp(std::uint32_t* state, int lane) {
  constexpr std::uint32_t mask = 0x3fffffffU;
  for (int j = lane; j < 273; j += 32)
    state[j] = mask & (state[j] + state[j + 334]);
  __syncwarp();
  for (int j = 273 + lane; j < 546; j += 32)
    state[j] = mask & (state[j] + state[j - 273]);
  __syncwarp();
  for (int j = 546 + lane; j < 607; j += 32)
    state[j] = mask & (state[j] + state[j - 273]);
  __syncwarp();
}

__device__ void cloud_rng_initialize_warp(
    int seed, std::uint32_t* state, int lane) {
  constexpr std::uint32_t seed_mask = 123459876U;
  std::int32_t signed_value = static_cast<std::int32_t>(
      static_cast<std::uint32_t>(seed) ^ seed_mask);
  if (signed_value < 0) signed_value = -signed_value;
  std::uint32_t initial_value = static_cast<std::uint32_t>(signed_value);
  if (initial_value == 0) initial_value = seed_mask;
  for (int spin = 0; spin < 64; ++spin)
    initial_value = cloud_rng_step(initial_value);

  std::uint32_t value = initial_value;
  for (int plane = 0; plane < lane && lane < 29; ++plane)
    value = cloud_rng_advance_604(value);
  for (int j = 2; j <= 605; ++j) {
    const unsigned int plane_bits = __ballot_sync(
        0xffffffffU, lane < 29 && (value & 0x80000000U));
    if (lane == 0) state[j] = (plane_bits & 0x1fffffffU) << 1U;
    if (lane < 29) value = cloud_rng_step(value);
  }
  if (lane == 0) {
    state[0] = 0;
    state[1] = (initial_value & 0x1fffffffU) << 1U;
    state[606] = (initial_value >> 29U) & 0x7U;
    state[501] |= 1U;
  }
  __syncwarp();

  // The 999-value warm-up triggers two state generations and leaves the
  // second generated state positioned at element 392.
  cloud_rng_generate_warp(state, lane);
  cloud_rng_generate_warp(state, lane);
}

__device__ inline double cloud_beta_to_alpha(
    double beta, double frac1, double frac2) {
  if (beta >= 1.0) return 1.0;
  const double difference = fabs(frac1 - frac2);
  return beta + (1.0 - beta) * difference /
      (difference + 1.0 / beta - 1.0);
}

__device__ inline double cloud_pdf_sample(
    double fsd, double cdf, int ncdf, int nfsd, double fsd1,
    double inv_fsd_interval, const double* values) {
  double weighted_cdf = cdf * (ncdf - 1) + 1.0;
  int icdf = static_cast<int>(weighted_cdf);
  icdf = icdf < 1 ? 1 : (icdf > ncdf - 1 ? ncdf - 1 : icdf);
  weighted_cdf = fmax(0.0, fmin(weighted_cdf - icdf, 1.0));
  double weighted_fsd = (fsd - fsd1) * inv_fsd_interval + 1.0;
  int ifsd = static_cast<int>(weighted_fsd);
  ifsd = ifsd < 1 ? 1 : (ifsd > nfsd - 1 ? nfsd - 1 : ifsd);
  weighted_fsd = fmax(0.0, fmin(weighted_fsd - ifsd, 1.0));
  const int c0 = icdf - 1;
  const int f0 = ifsd - 1;
  return (1.0 - weighted_cdf) * (1.0 - weighted_fsd) * values[c0 + ncdf * f0]
      + (1.0 - weighted_cdf) * weighted_fsd * values[c0 + ncdf * (f0 + 1)]
      + weighted_cdf * (1.0 - weighted_fsd) * values[c0 + 1 + ncdf * f0]
      + weighted_cdf * weighted_fsd * values[c0 + 1 + ncdf * (f0 + 1)];
}

__global__ void cloud_generator_kernel(
    int ng, int nlev, int ncol, int overlap_scheme, bool beta_overlap,
    const int* seeds, const double* active, double fraction_threshold,
    const double* fraction, const double* overlap_parameter,
    double decorrelation_scaling, const double* fractional_std,
    int pdf_ncdf, int pdf_nfsd, double pdf_fsd1,
    double pdf_inv_fsd_interval, const double* pdf_values,
    double* od_scaling, double* total_cloud_cover,
    double* cumulative_cover, double* pair_cover, double* overhang,
    double* overlap_inhom, double* random_top, double* random_cloud,
    double* random_inhom1, double* random_inhom2,
    std::uint32_t* random_state) {
  const int lane = threadIdx.x & 31;
  const int warps_per_block = blockDim.x / 32;
  const int col = blockIdx.x * warps_per_block + threadIdx.x / 32;
  if (col >= ncol) return;
  const auto profile = [=](int lev) {
    return col + static_cast<std::size_t>(ncol) * lev;
  };
  const auto scaling = [=](int g, int lev) {
    return g + static_cast<std::size_t>(ng) *
        (lev + static_cast<std::size_t>(nlev) * col);
  };
  constexpr double max_cloud_fraction = 1.0 - 2.2204460492503131e-15;
  double cover = 0.0;
  int begin = 0;
  int end = 0;
  int generate_randoms = 0;
  if (lane == 0) {
    for (int lev = 0; lev < nlev; ++lev)
      for (int g = 0; g < ng; ++g) od_scaling[scaling(g, lev)] = 0.0;
    total_cloud_cover[col] = 0.0;
    if (active[col] > 0.0) {
      double cumulative_product = 1.0 - fraction[profile(0)];
      cumulative_cover[profile(0)] = fraction[profile(0)];
      for (int lev = 0; lev < nlev - 1; ++lev) {
        const double upper = fraction[profile(lev)];
        const double lower = fraction[profile(lev + 1)];
        double pair;
        if (overlap_scheme == 0) {
          pair = upper > lower ? upper : lower;
        } else {
          double alpha = overlap_parameter[profile(lev)];
          if (beta_overlap) alpha = cloud_beta_to_alpha(alpha, upper, lower);
          const double maximum = upper > lower ? upper : lower;
          pair = alpha * maximum + (1.0 - alpha) *
              (upper + lower - upper * lower);
        }
        pair_cover[profile(lev)] = pair;
        if (upper >= max_cloud_fraction) cumulative_product = 0.0;
        else cumulative_product *= (1.0 - pair) / (1.0 - upper);
        cumulative_cover[profile(lev + 1)] = 1.0 - cumulative_product;
        overhang[profile(lev)] = cumulative_cover[profile(lev + 1)]
            - cumulative_cover[profile(lev)];
      }
      cover = cumulative_cover[profile(nlev - 1)];
      if (cover >= fraction_threshold) {
        total_cloud_cover[col] = cover;
        while (begin < nlev && fraction[profile(begin)] <= 0.0) ++begin;
        if (begin < nlev) {
          end = begin;
          for (int lev = begin + 1; lev < nlev; ++lev)
            if (fraction[profile(lev)] > 0.0) end = lev;
          for (int lev = 0; lev < nlev - 1; ++lev)
            overlap_inhom[profile(lev)] = overlap_parameter[profile(lev)];
          for (int lev = begin; lev < end; ++lev) {
            const double overlap = overlap_parameter[profile(lev)];
            if (overlap > 0.0)
              overlap_inhom[profile(lev)] = pow(
                  overlap, 1.0 / decorrelation_scaling);
          }
          generate_randoms = 1;
        }
      }
    }
  }
  generate_randoms = __shfl_sync(0xffffffffU, generate_randoms, 0);
  if (!generate_randoms) return;

  std::uint32_t* state = random_state + static_cast<std::size_t>(607) * col;
  cloud_rng_initialize_warp(seeds[col], state, lane);
  if (lane != 0) return;
  int used = 392;
  for (int g = 0; g < ng; ++g)
    random_top[g + static_cast<std::size_t>(ng) * col] = cloud_rng_next(state, used);
  for (int g = 0; g < ng; ++g) {
    const double trigger = random_top[g + static_cast<std::size_t>(ng) * col] * cover;
    int trigger_level = begin;
    while (trigger > cumulative_cover[profile(trigger_level)] &&
           trigger_level < end) ++trigger_level;

    const int cloud_random_count = end + 1 - trigger_level;
    for (int j = 0; j < cloud_random_count; ++j)
      random_cloud[profile(j)] = cloud_rng_next(state, used);
    int random_index = 0;
    int layers_to_scale = 1;
    for (int level = trigger_level + 1; level <= end + 1; ++level) {
      bool fill_scaling = false;
      if (level <= end) {
        const double random = random_cloud[profile(random_index++)];
        const int above = level - 1;
        if (layers_to_scale > 0) {
          if (random * fraction[profile(above)] <
              fraction[profile(level)] + fraction[profile(above)]
                  - pair_cover[profile(above)]) {
            ++layers_to_scale;
          } else {
            fill_scaling = true;
          }
        } else if (random * (cumulative_cover[profile(above)]
                       - fraction[profile(above)]) <
                   pair_cover[profile(above)] - overhang[profile(above)]
                       - fraction[profile(above)]) {
          layers_to_scale = 1;
        }
      } else {
        fill_scaling = true;
      }

      if (fill_scaling) {
        for (int j = 0; j < layers_to_scale; ++j)
          random_inhom1[profile(j)] = cloud_rng_next(state, used);
        for (int j = 0; j < layers_to_scale; ++j)
          random_inhom2[profile(j)] = cloud_rng_next(state, used);
        const int first_level = level - layers_to_scale;
        for (int j = 1; j < layers_to_scale; ++j) {
          if (random_inhom2[profile(j)] <
              overlap_inhom[profile(first_level + j - 1)])
            random_inhom1[profile(j)] = random_inhom1[profile(j - 1)];
        }
        for (int j = 0; j < layers_to_scale; ++j) {
          const int cloud_level = first_level + j;
          od_scaling[scaling(g, cloud_level)] = cloud_pdf_sample(
              fractional_std[profile(cloud_level)], random_inhom1[profile(j)],
              pdf_ncdf, pdf_nfsd, pdf_fsd1, pdf_inv_fsd_interval, pdf_values);
        }
        layers_to_scale = 0;
      }
    }
  }
}

bool copy_to_device(void* dst, const void* src, std::size_t bytes, const char* name) {
  return cuda_ok(cudaMemcpyAsync(dst, src, bytes, cudaMemcpyHostToDevice, workspace.stream), name);
}

bool copy_to_host(void* dst, const void* src, std::size_t bytes, const char* name) {
  return cuda_ok(cudaMemcpyAsync(dst, src, bytes, cudaMemcpyDeviceToHost, workspace.stream), name);
}

bool copy_to_longwave_device(void* dst, const void* src, std::size_t bytes, const char* name) {
  return cuda_ok(cudaMemcpyAsync(dst, src, bytes, cudaMemcpyHostToDevice,
                                 longwave_workspace.stream), name);
}

bool copy_from_longwave_device(void* dst, const void* src, std::size_t bytes, const char* name) {
  return cuda_ok(cudaMemcpyAsync(dst, src, bytes, cudaMemcpyDeviceToHost,
                                 longwave_workspace.stream), name);
}

bool copy_to_cloud_device(void* dst, const void* src, std::size_t bytes,
                          cudaStream_t stream, const char* name) {
  return cuda_ok(cudaMemcpyAsync(dst, src, bytes, cudaMemcpyHostToDevice,
                                 stream), name);
}

bool copy_from_cloud_device(void* dst, const void* src, std::size_t bytes,
                            cudaStream_t stream, const char* name) {
  return cuda_ok(cudaMemcpyAsync(dst, src, bytes, cudaMemcpyDeviceToHost,
                                 stream), name);
}

}  // namespace

extern "C" int oifs_cuda_radiation_available(void) {
  std::lock_guard<std::mutex> lock(workspace_mutex);
  last_error.clear();
  return select_device() ? 1 : 0;
}

extern "C" int oifs_cuda_sw_compute_dp(
    int ng, int nbands, int nlev, int ncol, int do_delta_scaling,
    const double* mu0, const double* od, const double* ssa, const double* asymmetry,
    const double* albedo_direct, const double* albedo_diffuse, const double* incoming_sw,
    const int* band_from_g, double cloud_fraction_threshold,
    const double* cloud_fraction, const double* total_cloud_cover,
    const double* od_scaling,
    const double* od_cloud, const double* ssa_cloud, const double* asymmetry_cloud,
    double* sw_up_clear, double* sw_dn_clear, double* sw_dn_direct_clear,
    double* sw_up, double* sw_dn, double* sw_dn_direct,
    double* sw_dn_diffuse_surf_clear_g, double* sw_dn_direct_surf_clear_g,
    double* sw_dn_diffuse_surf_g, double* sw_dn_direct_surf_g) {
  std::lock_guard<std::mutex> lock(workspace_mutex);
  last_error.clear();
  if (ng <= 0 || nbands <= 0 || nlev <= 0 || ncol <= 0) {
    last_error = "invalid shortwave dimensions";
    return 1;
  }
  if (!select_device()) return 2;
  if (!ensure_workspace(ng, nbands, nlev, ncol)) return 2;
  const std::size_t layer = static_cast<std::size_t>(ng) * nlev * ncol;
  const std::size_t gpcol = static_cast<std::size_t>(ng) * ncol;
  const std::size_t cloud_layer = static_cast<std::size_t>(nbands) * nlev * ncol;
  const std::size_t fraction = static_cast<std::size_t>(ncol) * nlev;
  const std::size_t profile = static_cast<std::size_t>(ncol) * (nlev + 1);
  CloudWorkspace* cloud = ready_cloud_workspace;
  const bool use_device_cloud = cloud && cloud->ng == ng &&
      cloud->nlev == nlev && cloud->ncol >= ncol;
  const double* device_od_scaling = use_device_cloud ? cloud->od_scaling
                                                     : workspace.od_scaling;
  const double* device_cloud_cover = use_device_cloud ? cloud->total_cloud_cover
                                                      : workspace.total_cloud_cover;
  ready_cloud_workspace = nullptr;
#define COPY_IN(member, src, count) if (!copy_to_device(workspace.member, src, (count) * sizeof(*(src)), "copy " #member)) return 3
  COPY_IN(mu0, mu0, ncol);
  COPY_IN(od, od, layer); COPY_IN(ssa, ssa, layer); COPY_IN(asymmetry, asymmetry, layer);
  COPY_IN(albedo_direct, albedo_direct, gpcol); COPY_IN(albedo_diffuse, albedo_diffuse, gpcol);
  COPY_IN(incoming_sw, incoming_sw, gpcol); COPY_IN(band_from_g, band_from_g, ng);
  COPY_IN(cloud_fraction, cloud_fraction, fraction);
  if (!use_device_cloud) {
    COPY_IN(total_cloud_cover, total_cloud_cover, ncol);
    COPY_IN(od_scaling, od_scaling, layer);
  }
  COPY_IN(od_cloud, od_cloud, cloud_layer); COPY_IN(ssa_cloud, ssa_cloud, cloud_layer);
  COPY_IN(asymmetry_cloud, asymmetry_cloud, cloud_layer);
#undef COPY_IN
  constexpr int block_size = 256;
  const int optics_blocks = static_cast<int>((layer + block_size - 1) / block_size);
  const int adding_blocks = static_cast<int>((gpcol + block_size - 1) / block_size);
  const int profile_blocks = static_cast<int>((profile + block_size - 1) / block_size);
  optics_kernel<<<optics_blocks, block_size, 0, workspace.stream>>>(
      ng, nbands, nlev, ncol, false, do_delta_scaling != 0, cloud_fraction_threshold,
      workspace.mu0, workspace.od, workspace.ssa, workspace.asymmetry,
      workspace.band_from_g, workspace.cloud_fraction, device_od_scaling,
      workspace.od_cloud, workspace.ssa_cloud, workspace.asymmetry_cloud,
      workspace.ref_diff, workspace.trans_diff, workspace.ref_dir,
      workspace.trans_dir_diff, workspace.trans_dir_dir);
  adding_kernel<<<adding_blocks, block_size, 0, workspace.stream>>>(
      ng, nlev, ncol, workspace.mu0, workspace.incoming_sw,
      workspace.albedo_direct, workspace.albedo_diffuse, workspace.ref_diff,
      workspace.trans_diff, workspace.ref_dir, workspace.trans_dir_diff,
      workspace.trans_dir_dir, workspace.albedo, workspace.source,
      workspace.inv_denominator, workspace.flux_up, workspace.flux_dn_diffuse,
      workspace.flux_dn_direct);
  optics_kernel<<<optics_blocks, block_size, 0, workspace.stream>>>(
      ng, nbands, nlev, ncol, true, do_delta_scaling != 0, cloud_fraction_threshold,
      workspace.mu0, workspace.od, workspace.ssa, workspace.asymmetry,
      workspace.band_from_g, workspace.cloud_fraction, device_od_scaling,
      workspace.od_cloud, workspace.ssa_cloud, workspace.asymmetry_cloud,
      workspace.ref_diff, workspace.trans_diff, workspace.ref_dir,
      workspace.trans_dir_diff, workspace.trans_dir_dir);
  adding_kernel<<<adding_blocks, block_size, 0, workspace.stream>>>(
      ng, nlev, ncol, workspace.mu0, workspace.incoming_sw,
      workspace.albedo_direct, workspace.albedo_diffuse, workspace.ref_diff,
      workspace.trans_diff, workspace.ref_dir, workspace.trans_dir_diff,
      workspace.trans_dir_dir, workspace.albedo, workspace.source,
      workspace.inv_denominator, workspace.flux_up_cloud,
      workspace.flux_dn_diffuse_cloud, workspace.flux_dn_direct_cloud);
  reduce_profiles_kernel<<<profile_blocks, block_size, 0, workspace.stream>>>(
      ng, nlev, ncol, device_cloud_cover, workspace.flux_up,
      workspace.flux_dn_diffuse, workspace.flux_dn_direct,
      workspace.flux_up_cloud, workspace.flux_dn_diffuse_cloud,
      workspace.flux_dn_direct_cloud, workspace.sw_up_clear,
      workspace.sw_dn_clear, workspace.sw_dn_direct_clear, workspace.sw_up,
      workspace.sw_dn, workspace.sw_dn_direct);
  reduce_surface_kernel<<<adding_blocks, block_size, 0, workspace.stream>>>(
      ng, nlev, ncol, device_cloud_cover, workspace.flux_dn_diffuse,
      workspace.flux_dn_direct, workspace.flux_dn_diffuse_cloud,
      workspace.flux_dn_direct_cloud, workspace.sw_dn_diffuse_surf_clear_g,
      workspace.sw_dn_direct_surf_clear_g, workspace.sw_dn_diffuse_surf_g,
      workspace.sw_dn_direct_surf_g);
  if (!cuda_ok(cudaGetLastError(), "shortwave kernel launch")) return 4;
#define COPY_PROFILE(dst, member) if (!copy_to_host(dst, workspace.member, profile * sizeof(double), "copy " #member)) return 5
  COPY_PROFILE(sw_up_clear, sw_up_clear); COPY_PROFILE(sw_dn_clear, sw_dn_clear);
  COPY_PROFILE(sw_dn_direct_clear, sw_dn_direct_clear);
  COPY_PROFILE(sw_up, sw_up); COPY_PROFILE(sw_dn, sw_dn);
  COPY_PROFILE(sw_dn_direct, sw_dn_direct);
#undef COPY_PROFILE
#define COPY_SURFACE(dst, member) if (!copy_to_host(dst, workspace.member, gpcol * sizeof(double), "copy " #member)) return 5
  COPY_SURFACE(sw_dn_diffuse_surf_clear_g, sw_dn_diffuse_surf_clear_g);
  COPY_SURFACE(sw_dn_direct_surf_clear_g, sw_dn_direct_surf_clear_g);
  COPY_SURFACE(sw_dn_diffuse_surf_g, sw_dn_diffuse_surf_g);
  COPY_SURFACE(sw_dn_direct_surf_g, sw_dn_direct_surf_g);
#undef COPY_SURFACE
  if (!cuda_ok(cudaStreamSynchronize(workspace.stream), "shortwave synchronize")) return 6;
  return 0;
}

extern "C" int oifs_cuda_lw_compute_dp(
    int ng, int nbands, int nlev, int ncol,
    int do_aerosol_scattering, int do_cloud_scattering, int do_derivatives,
    const double* od, const double* ssa, const double* asymmetry,
    const double* planck_hl, const double* emission, const double* albedo,
    const int* band_from_g, double cloud_fraction_threshold,
    const double* cloud_fraction, const double* total_cloud_cover,
    const double* od_scaling, const double* od_cloud,
    const double* ssa_cloud, const double* asymmetry_cloud,
    double* lw_up_clear, double* lw_dn_clear, double* lw_up, double* lw_dn,
    double* lw_dn_surf_clear_g, double* lw_dn_surf_g,
    double* lw_derivatives) {
  std::lock_guard<std::mutex> lock(workspace_mutex);
  last_error.clear();
  if (ng <= 0 || nbands <= 0 || nlev <= 0 || ncol <= 0) {
    last_error = "invalid longwave dimensions";
    return 1;
  }
  if (!select_device()) return 2;
  if (!ensure_longwave_workspace(ng, nbands, nlev, ncol)) return 2;
  auto& lw = longwave_workspace;
  const std::size_t layer = static_cast<std::size_t>(ng) * nlev * ncol;
  const std::size_t interface = static_cast<std::size_t>(ng) * (nlev + 1) * ncol;
  const std::size_t gpcol = static_cast<std::size_t>(ng) * ncol;
  const std::size_t cloud_layer = static_cast<std::size_t>(nbands) * nlev * ncol;
  const std::size_t fraction = static_cast<std::size_t>(ncol) * nlev;
  const std::size_t profile = static_cast<std::size_t>(ncol) * (nlev + 1);
  CloudWorkspace* cloud = ready_cloud_workspace;
  const bool use_device_cloud = cloud && cloud->ng == ng &&
      cloud->nlev == nlev && cloud->ncol >= ncol;
  const double* device_od_scaling = use_device_cloud ? cloud->od_scaling
                                                     : lw.od_scaling;
  const double* device_cloud_cover = use_device_cloud ? cloud->total_cloud_cover
                                                      : lw.total_cloud_cover;
  ready_cloud_workspace = nullptr;
#define COPY_LW_IN(member, src, count)                                            \
  if (!copy_to_longwave_device(lw.member, src, (count) * sizeof(*(src)),          \
                               "copy lw " #member)) return 3
  COPY_LW_IN(od, od, layer);
  if (do_aerosol_scattering) {
    COPY_LW_IN(ssa, ssa, layer);
    COPY_LW_IN(asymmetry, asymmetry, layer);
  }
  COPY_LW_IN(planck_hl, planck_hl, interface);
  COPY_LW_IN(emission, emission, gpcol);
  COPY_LW_IN(albedo_surface, albedo, gpcol);
  COPY_LW_IN(band_from_g, band_from_g, ng);
  COPY_LW_IN(cloud_fraction, cloud_fraction, fraction);
  if (!use_device_cloud) {
    COPY_LW_IN(total_cloud_cover, total_cloud_cover, ncol);
    COPY_LW_IN(od_scaling, od_scaling, layer);
  }
  COPY_LW_IN(od_cloud, od_cloud, cloud_layer);
  if (do_cloud_scattering) {
    COPY_LW_IN(ssa_cloud, ssa_cloud, cloud_layer);
    COPY_LW_IN(asymmetry_cloud, asymmetry_cloud, cloud_layer);
  }
#undef COPY_LW_IN

  constexpr int block_size = 256;
  const int layer_blocks = static_cast<int>((layer + block_size - 1) / block_size);
  const int column_blocks = static_cast<int>((gpcol + block_size - 1) / block_size);
  const int profile_blocks = static_cast<int>((profile + block_size - 1) / block_size);
  const int derivative_blocks = (ncol + block_size - 1) / block_size;
  longwave_optics_kernel<<<layer_blocks, block_size, 0, lw.stream>>>(
      ng, nbands, nlev, ncol, false, do_aerosol_scattering != 0,
      do_cloud_scattering != 0, cloud_fraction_threshold, lw.od, lw.ssa,
      lw.asymmetry, lw.planck_hl, lw.band_from_g, lw.cloud_fraction,
      device_od_scaling, lw.od_cloud, lw.ssa_cloud, lw.asymmetry_cloud,
      lw.ref_clear, lw.trans_clear, lw.source_up_clear, lw.source_dn_clear);
  longwave_adding_kernel<<<column_blocks, block_size, 0, lw.stream>>>(
      ng, nlev, ncol, do_aerosol_scattering != 0, lw.ref_clear, lw.trans_clear,
      lw.source_up_clear, lw.source_dn_clear, lw.emission, lw.albedo_surface,
      lw.albedo, lw.source, lw.inv_denominator, lw.flux_up_clear, lw.flux_dn_clear);
  longwave_optics_kernel<<<layer_blocks, block_size, 0, lw.stream>>>(
      ng, nbands, nlev, ncol, true, do_aerosol_scattering != 0,
      do_cloud_scattering != 0, cloud_fraction_threshold, lw.od, lw.ssa,
      lw.asymmetry, lw.planck_hl, lw.band_from_g, lw.cloud_fraction,
      device_od_scaling, lw.od_cloud, lw.ssa_cloud, lw.asymmetry_cloud,
      lw.ref_cloud, lw.trans_cloud, lw.source_up_cloud, lw.source_dn_cloud);
  longwave_adding_kernel<<<column_blocks, block_size, 0, lw.stream>>>(
      ng, nlev, ncol, (do_aerosol_scattering || do_cloud_scattering) != 0,
      lw.ref_cloud, lw.trans_cloud, lw.source_up_cloud, lw.source_dn_cloud,
      lw.emission, lw.albedo_surface, lw.albedo, lw.source, lw.inv_denominator,
      lw.flux_up_cloud, lw.flux_dn_cloud);
  reduce_longwave_kernel<<<profile_blocks, block_size, 0, lw.stream>>>(
      ng, nlev, ncol, cloud_fraction_threshold, device_cloud_cover,
      lw.flux_up_clear, lw.flux_dn_clear, lw.flux_up_cloud, lw.flux_dn_cloud,
      lw.lw_up_clear, lw.lw_dn_clear, lw.lw_up, lw.lw_dn);
  reduce_longwave_surface_kernel<<<column_blocks, block_size, 0, lw.stream>>>(
      ng, nlev, ncol, cloud_fraction_threshold, device_cloud_cover,
      lw.flux_dn_clear, lw.flux_dn_cloud, lw.lw_dn_surf_clear_g,
      lw.lw_dn_surf_g);
  if (do_derivatives) {
    longwave_derivatives_kernel<<<derivative_blocks, block_size, 0, lw.stream>>>(
        ng, nlev, ncol, cloud_fraction_threshold, device_cloud_cover,
        lw.trans_clear, lw.trans_cloud, lw.flux_up_clear, lw.flux_up_cloud,
        lw.derivative_g, lw.lw_derivatives);
  } else {
    cudaMemsetAsync(lw.lw_derivatives, 0, profile * sizeof(double), lw.stream);
  }
  if (!cuda_ok(cudaGetLastError(), "longwave kernel launch")) return 4;
#define COPY_LW_OUT(dst, member, count)                                           \
  if (!copy_from_longwave_device(dst, lw.member, (count) * sizeof(*(dst)),        \
                                  "copy lw " #member)) return 5
  COPY_LW_OUT(lw_up_clear, lw_up_clear, profile);
  COPY_LW_OUT(lw_dn_clear, lw_dn_clear, profile);
  COPY_LW_OUT(lw_up, lw_up, profile);
  COPY_LW_OUT(lw_dn, lw_dn, profile);
  COPY_LW_OUT(lw_dn_surf_clear_g, lw_dn_surf_clear_g, gpcol);
  COPY_LW_OUT(lw_dn_surf_g, lw_dn_surf_g, gpcol);
  COPY_LW_OUT(lw_derivatives, lw_derivatives, profile);
#undef COPY_LW_OUT
  if (!cuda_ok(cudaStreamSynchronize(lw.stream), "longwave synchronize")) return 6;
  return 0;
}

extern "C" int oifs_cuda_cloud_compute_dp(
    int ng, int nlev, int ncol, int overlap_scheme, int is_beta_overlap,
    int device_handoff,
    const int* seeds, const double* active, double frac_threshold,
    const double* cloud_fraction, const double* overlap_parameter,
    double decorrelation_scaling, const double* fractional_std,
    int pdf_ncdf, int pdf_nfsd, double pdf_fsd1,
    double pdf_inv_fsd_interval, const double* pdf_values,
    double* od_scaling, double* total_cloud_cover) {
  std::lock_guard<std::mutex> lock(workspace_mutex);
  last_error.clear();
  ready_cloud_workspace = nullptr;
  if (ng <= 0 || nlev <= 1 || ncol <= 0 || pdf_ncdf < 2 || pdf_nfsd < 2) {
    last_error = "invalid cloud-generator dimensions";
    return 1;
  }
  if (overlap_scheme != 0 && overlap_scheme != 1) {
    last_error = "cloud overlap scheme is not implemented by CUDA";
    return 7;
  }
  if (!select_device()) return 2;
  CloudWorkspace* selected_cloud = nullptr;
  if (!ensure_cloud_workspace(ng, nlev, ncol, pdf_ncdf, pdf_nfsd,
                              selected_cloud)) return 2;
  auto& cloud = *selected_cloud;
  const std::size_t profile = static_cast<std::size_t>(ncol) * nlev;
  const std::size_t interfaces = static_cast<std::size_t>(ncol) * (nlev - 1);
  const std::size_t scaling = static_cast<std::size_t>(ng) * nlev * ncol;
  const std::size_t pdf_size = static_cast<std::size_t>(pdf_ncdf) * pdf_nfsd;
#define COPY_CLOUD_IN(member, src, count)                                        \
  if (!copy_to_cloud_device(cloud.member, src, (count) * sizeof(*(src)),          \
                            cloud.stream, "copy cloud " #member)) return 3
  COPY_CLOUD_IN(seeds, seeds, ncol);
  COPY_CLOUD_IN(active, active, ncol);
  COPY_CLOUD_IN(cloud_fraction, cloud_fraction, profile);
  COPY_CLOUD_IN(overlap_parameter, overlap_parameter, interfaces);
  COPY_CLOUD_IN(fractional_std, fractional_std, profile);
  COPY_CLOUD_IN(pdf_values, pdf_values, pdf_size);
#undef COPY_CLOUD_IN
  constexpr int block_size = 64;
  constexpr int columns_per_block = block_size / 32;
  const int blocks = (ncol + columns_per_block - 1) / columns_per_block;
  cloud_generator_kernel<<<blocks, block_size, 0, cloud.stream>>>(
      ng, nlev, ncol, overlap_scheme, is_beta_overlap != 0, cloud.seeds,
      cloud.active, frac_threshold, cloud.cloud_fraction,
      cloud.overlap_parameter, decorrelation_scaling, cloud.fractional_std,
      pdf_ncdf, pdf_nfsd, pdf_fsd1, pdf_inv_fsd_interval, cloud.pdf_values,
      cloud.od_scaling, cloud.total_cloud_cover, cloud.cumulative_cover,
      cloud.pair_cover, cloud.overhang, cloud.overlap_inhom, cloud.random_top,
      cloud.random_cloud, cloud.random_inhom1, cloud.random_inhom2,
      cloud.random_state);
  if (!cuda_ok(cudaGetLastError(), "cloud-generator kernel launch")) return 4;
  if (!device_handoff && !copy_from_cloud_device(
          od_scaling, cloud.od_scaling, scaling * sizeof(double),
          cloud.stream, "copy cloud od_scaling")) return 5;
  if (!copy_from_cloud_device(total_cloud_cover, cloud.total_cloud_cover,
                              ncol * sizeof(double), cloud.stream,
                              "copy cloud cover")) return 5;
  if (!cuda_ok(cudaStreamSynchronize(cloud.stream),
               "cloud-generator synchronize")) return 6;
  ready_cloud_workspace = selected_cloud;
  return 0;
}

extern "C" void oifs_cuda_radiation_finalize(void) {
  std::lock_guard<std::mutex> lock(workspace_mutex);
  release_workspace();
  release_longwave_workspace();
  release_cloud_workspaces();
}

extern "C" const char* oifs_cuda_radiation_last_error(void) {
  return last_error.c_str();
}
