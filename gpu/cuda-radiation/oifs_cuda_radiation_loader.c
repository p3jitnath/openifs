#include <dlfcn.h>
#include <pthread.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>

typedef int (*available_fn)(void);
typedef int (*compute_fn)(
    int, int, int, int, int, const double *, const double *, const double *,
    const double *, const double *, const double *, const double *, const int *,
    double, const double *, const double *, const double *, const double *,
    const double *, const double *, double *, double *, double *, double *,
    double *, double *, double *, double *, double *, double *);
typedef int (*compute_lw_fn)(
    int, int, int, int, int, int, int, const double *, const double *,
    const double *, const double *, const double *, const double *, const int *,
    double, const double *, const double *, const double *, const double *,
    const double *, const double *, double *, double *, double *, double *,
    double *, double *, double *);
typedef int (*compute_cloud_fn)(
    int, int, int, int, int, const int *, const double *, double,
    const double *, const double *, double, const double *, int, int,
    double, double, const double *, double *, double *);
typedef const char *(*last_error_fn)(void);

static pthread_once_t load_once = PTHREAD_ONCE_INIT;
static void *cuda_library;
static available_fn cuda_available;
static compute_fn cuda_compute;
static compute_lw_fn cuda_compute_lw;
static compute_cloud_fn cuda_compute_cloud;
static last_error_fn cuda_last_error;
static char load_error[512];

static void load_cuda_library(void) {
  const char *path = getenv("OIFS_CUDA_RADIATION_LIBRARY");
  if (path == NULL || path[0] == '\0') {
    snprintf(load_error, sizeof(load_error),
             "OIFS_CUDA_RADIATION_LIBRARY is not set");
    return;
  }
  cuda_library = dlopen(path, RTLD_NOW | RTLD_LOCAL);
  if (cuda_library == NULL) {
    snprintf(load_error, sizeof(load_error), "%s", dlerror());
    return;
  }
  cuda_available = (available_fn)dlsym(cuda_library,
                                      "oifs_cuda_radiation_available");
  cuda_compute = (compute_fn)dlsym(cuda_library, "oifs_cuda_sw_compute_dp");
  cuda_compute_lw = (compute_lw_fn)dlsym(cuda_library, "oifs_cuda_lw_compute_dp");
  cuda_compute_cloud = (compute_cloud_fn)dlsym(
      cuda_library, "oifs_cuda_cloud_compute_dp");
  cuda_last_error = (last_error_fn)dlsym(
      cuda_library, "oifs_cuda_radiation_last_error");
  if (cuda_available == NULL || cuda_compute == NULL || cuda_compute_lw == NULL ||
      cuda_compute_cloud == NULL || cuda_last_error == NULL) {
    snprintf(load_error, sizeof(load_error),
             "CUDA radiation library has an incompatible API");
    dlclose(cuda_library);
    cuda_library = NULL;
    cuda_available = NULL;
    cuda_compute = NULL;
    cuda_compute_lw = NULL;
    cuda_compute_cloud = NULL;
    cuda_last_error = NULL;
  }
}

int oifs_cuda_bridge_cloud_compute_dp(
    int ng, int nlev, int ncol, int overlap_scheme, int is_beta_overlap,
    const int *seeds, const double *active, double frac_threshold,
    const double *cloud_fraction, const double *overlap_parameter,
    double decorrelation_scaling, const double *fractional_std,
    int pdf_ncdf, int pdf_nfsd, double pdf_fsd1,
    double pdf_inv_fsd_interval, const double *pdf_values,
    double *od_scaling, double *total_cloud_cover) {
  pthread_once(&load_once, load_cuda_library);
  if (cuda_compute_cloud == NULL) return -1;
  return cuda_compute_cloud(
      ng, nlev, ncol, overlap_scheme, is_beta_overlap, seeds, active,
      frac_threshold, cloud_fraction, overlap_parameter,
      decorrelation_scaling, fractional_std, pdf_ncdf, pdf_nfsd, pdf_fsd1,
      pdf_inv_fsd_interval, pdf_values, od_scaling, total_cloud_cover);
}

int oifs_cuda_bridge_lw_compute_dp(
    int ng, int nbands, int nlev, int ncol,
    int do_aerosol_scattering, int do_cloud_scattering, int do_derivatives,
    const double *od, const double *ssa, const double *asymmetry,
    const double *planck_hl, const double *emission, const double *albedo,
    const int *band_from_g, double cloud_fraction_threshold,
    const double *cloud_fraction, const double *total_cloud_cover,
    const double *od_scaling, const double *od_cloud,
    const double *ssa_cloud, const double *asymmetry_cloud,
    double *lw_up_clear, double *lw_dn_clear, double *lw_up, double *lw_dn,
    double *lw_dn_surf_clear_g, double *lw_dn_surf_g,
    double *lw_derivatives) {
  pthread_once(&load_once, load_cuda_library);
  if (cuda_compute_lw == NULL) return -1;
  return cuda_compute_lw(
      ng, nbands, nlev, ncol, do_aerosol_scattering, do_cloud_scattering,
      do_derivatives, od, ssa, asymmetry, planck_hl, emission, albedo,
      band_from_g, cloud_fraction_threshold, cloud_fraction,
      total_cloud_cover, od_scaling, od_cloud, ssa_cloud, asymmetry_cloud,
      lw_up_clear, lw_dn_clear, lw_up, lw_dn, lw_dn_surf_clear_g,
      lw_dn_surf_g, lw_derivatives);
}

int oifs_cuda_bridge_available(void) {
  pthread_once(&load_once, load_cuda_library);
  return cuda_available != NULL && cuda_available() != 0;
}

int oifs_cuda_bridge_sw_compute_dp(
    int ng, int nbands, int nlev, int ncol, int do_delta_scaling,
    const double *mu0, const double *od, const double *ssa,
    const double *asymmetry, const double *albedo_direct,
    const double *albedo_diffuse, const double *incoming_sw,
    const int *band_from_g, double cloud_fraction_threshold,
    const double *cloud_fraction, const double *total_cloud_cover,
    const double *od_scaling,
    const double *od_cloud, const double *ssa_cloud,
    const double *asymmetry_cloud,
    double *sw_up_clear, double *sw_dn_clear, double *sw_dn_direct_clear,
    double *sw_up, double *sw_dn, double *sw_dn_direct,
    double *sw_dn_diffuse_surf_clear_g, double *sw_dn_direct_surf_clear_g,
    double *sw_dn_diffuse_surf_g, double *sw_dn_direct_surf_g) {
  pthread_once(&load_once, load_cuda_library);
  if (cuda_compute == NULL) return -1;
  return cuda_compute(
      ng, nbands, nlev, ncol, do_delta_scaling, mu0, od, ssa, asymmetry,
      albedo_direct, albedo_diffuse, incoming_sw, band_from_g,
      cloud_fraction_threshold, cloud_fraction, total_cloud_cover, od_scaling,
      od_cloud, ssa_cloud, asymmetry_cloud, sw_up_clear, sw_dn_clear,
      sw_dn_direct_clear, sw_up, sw_dn, sw_dn_direct,
      sw_dn_diffuse_surf_clear_g, sw_dn_direct_surf_clear_g,
      sw_dn_diffuse_surf_g, sw_dn_direct_surf_g);
}

const char *oifs_cuda_bridge_last_error(void) {
  pthread_once(&load_once, load_cuda_library);
  if (cuda_last_error != NULL) return cuda_last_error();
  return load_error;
}
