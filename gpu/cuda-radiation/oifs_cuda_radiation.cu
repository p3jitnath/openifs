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

bool copy_to_device(void* dst, const void* src, std::size_t bytes, const char* name) {
  return cuda_ok(cudaMemcpyAsync(dst, src, bytes, cudaMemcpyHostToDevice, workspace.stream), name);
}

bool copy_to_host(void* dst, const void* src, std::size_t bytes, const char* name) {
  return cuda_ok(cudaMemcpyAsync(dst, src, bytes, cudaMemcpyDeviceToHost, workspace.stream), name);
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

extern "C" void oifs_cuda_radiation_finalize(void) {
  std::lock_guard<std::mutex> lock(workspace_mutex);
  release_workspace();
}

extern "C" const char* oifs_cuda_radiation_last_error(void) {
  return last_error.c_str();
}
