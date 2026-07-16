#include "oifs_cuda_radiation.h"

#include <cuda_runtime.h>

#include <algorithm>
#include <cmath>
#include <cstddef>
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
#define COPY_IN(member, src, count) if (!copy_to_device(workspace.member, src, (count) * sizeof(*(src)), "copy " #member)) return 3
  COPY_IN(mu0, mu0, ncol);
  COPY_IN(od, od, layer); COPY_IN(ssa, ssa, layer); COPY_IN(asymmetry, asymmetry, layer);
  COPY_IN(albedo_direct, albedo_direct, gpcol); COPY_IN(albedo_diffuse, albedo_diffuse, gpcol);
  COPY_IN(incoming_sw, incoming_sw, gpcol); COPY_IN(band_from_g, band_from_g, ng);
  COPY_IN(cloud_fraction, cloud_fraction, fraction);
  COPY_IN(total_cloud_cover, total_cloud_cover, ncol);
  COPY_IN(od_scaling, od_scaling, layer);
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
      workspace.band_from_g, workspace.cloud_fraction, workspace.od_scaling,
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
      workspace.band_from_g, workspace.cloud_fraction, workspace.od_scaling,
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
      ng, nlev, ncol, workspace.total_cloud_cover, workspace.flux_up,
      workspace.flux_dn_diffuse, workspace.flux_dn_direct,
      workspace.flux_up_cloud, workspace.flux_dn_diffuse_cloud,
      workspace.flux_dn_direct_cloud, workspace.sw_up_clear,
      workspace.sw_dn_clear, workspace.sw_dn_direct_clear, workspace.sw_up,
      workspace.sw_dn, workspace.sw_dn_direct);
  reduce_surface_kernel<<<adding_blocks, block_size, 0, workspace.stream>>>(
      ng, nlev, ncol, workspace.total_cloud_cover, workspace.flux_dn_diffuse,
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
  COPY_LW_IN(total_cloud_cover, total_cloud_cover, ncol);
  COPY_LW_IN(od_scaling, od_scaling, layer);
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
      lw.od_scaling, lw.od_cloud, lw.ssa_cloud, lw.asymmetry_cloud,
      lw.ref_clear, lw.trans_clear, lw.source_up_clear, lw.source_dn_clear);
  longwave_adding_kernel<<<column_blocks, block_size, 0, lw.stream>>>(
      ng, nlev, ncol, do_aerosol_scattering != 0, lw.ref_clear, lw.trans_clear,
      lw.source_up_clear, lw.source_dn_clear, lw.emission, lw.albedo_surface,
      lw.albedo, lw.source, lw.inv_denominator, lw.flux_up_clear, lw.flux_dn_clear);
  longwave_optics_kernel<<<layer_blocks, block_size, 0, lw.stream>>>(
      ng, nbands, nlev, ncol, true, do_aerosol_scattering != 0,
      do_cloud_scattering != 0, cloud_fraction_threshold, lw.od, lw.ssa,
      lw.asymmetry, lw.planck_hl, lw.band_from_g, lw.cloud_fraction,
      lw.od_scaling, lw.od_cloud, lw.ssa_cloud, lw.asymmetry_cloud,
      lw.ref_cloud, lw.trans_cloud, lw.source_up_cloud, lw.source_dn_cloud);
  longwave_adding_kernel<<<column_blocks, block_size, 0, lw.stream>>>(
      ng, nlev, ncol, (do_aerosol_scattering || do_cloud_scattering) != 0,
      lw.ref_cloud, lw.trans_cloud, lw.source_up_cloud, lw.source_dn_cloud,
      lw.emission, lw.albedo_surface, lw.albedo, lw.source, lw.inv_denominator,
      lw.flux_up_cloud, lw.flux_dn_cloud);
  reduce_longwave_kernel<<<profile_blocks, block_size, 0, lw.stream>>>(
      ng, nlev, ncol, cloud_fraction_threshold, lw.total_cloud_cover,
      lw.flux_up_clear, lw.flux_dn_clear, lw.flux_up_cloud, lw.flux_dn_cloud,
      lw.lw_up_clear, lw.lw_dn_clear, lw.lw_up, lw.lw_dn);
  reduce_longwave_surface_kernel<<<column_blocks, block_size, 0, lw.stream>>>(
      ng, nlev, ncol, cloud_fraction_threshold, lw.total_cloud_cover,
      lw.flux_dn_clear, lw.flux_dn_cloud, lw.lw_dn_surf_clear_g,
      lw.lw_dn_surf_g);
  if (do_derivatives) {
    longwave_derivatives_kernel<<<derivative_blocks, block_size, 0, lw.stream>>>(
        ng, nlev, ncol, cloud_fraction_threshold, lw.total_cloud_cover,
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

extern "C" void oifs_cuda_radiation_finalize(void) {
  std::lock_guard<std::mutex> lock(workspace_mutex);
  release_workspace();
  release_longwave_workspace();
}

extern "C" const char* oifs_cuda_radiation_last_error(void) {
  return last_error.c_str();
}
