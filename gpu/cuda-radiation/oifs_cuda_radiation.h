#pragma once

#ifdef __cplusplus
extern "C" {
#endif

int oifs_cuda_radiation_available(void);

int oifs_cuda_sw_compute_dp(
    int ng, int nbands, int nlev, int ncol, int do_delta_scaling,
    const double* mu0, const double* od, const double* ssa,
    const double* asymmetry, const double* albedo_direct,
    const double* albedo_diffuse, const double* incoming_sw,
    const int* band_from_g, double cloud_fraction_threshold,
    const double* cloud_fraction, const double* total_cloud_cover,
    const double* od_scaling,
    const double* od_cloud, const double* ssa_cloud,
    const double* asymmetry_cloud,
    double* sw_up_clear, double* sw_dn_clear,
    double* sw_dn_direct_clear, double* sw_up, double* sw_dn,
    double* sw_dn_direct, double* sw_dn_diffuse_surf_clear_g,
    double* sw_dn_direct_surf_clear_g, double* sw_dn_diffuse_surf_g,
    double* sw_dn_direct_surf_g);

void oifs_cuda_radiation_finalize(void);
const char* oifs_cuda_radiation_last_error(void);

#ifdef __cplusplus
}
#endif
