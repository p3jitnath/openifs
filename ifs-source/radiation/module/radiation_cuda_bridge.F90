module radiation_cuda_bridge

  use parkind1, only : jprb
  use, intrinsic :: iso_c_binding, only : c_double, c_int

  implicit none

  private
  public :: cuda_radiation_available, cuda_sw_compute, cuda_lw_compute

#if defined(OIFS_CUDA_RADIATION) && !defined(PARKIND1_SINGLE)
  interface
    integer(c_int) function oifs_cuda_radiation_available() &
         & bind(C,name='oifs_cuda_bridge_available')
      import :: c_int
    end function oifs_cuda_radiation_available

    integer(c_int) function oifs_cuda_sw_compute_dp(ng,nbands,nlev,ncol, &
         & do_delta_scaling,mu0,od,ssa,asymmetry,albedo_direct,albedo_diffuse, &
         & incoming_sw,band_from_g,cloud_fraction_threshold,cloud_fraction, &
         & total_cloud_cover,od_scaling,od_cloud,ssa_cloud,asymmetry_cloud, &
         & sw_up_clear,sw_dn_clear,sw_dn_direct_clear,sw_up,sw_dn, &
         & sw_dn_direct,sw_dn_diffuse_surf_clear_g, &
         & sw_dn_direct_surf_clear_g,sw_dn_diffuse_surf_g, &
         & sw_dn_direct_surf_g) &
         & bind(C,name='oifs_cuda_bridge_sw_compute_dp')
      import :: c_double, c_int
      integer(c_int), value :: ng, nbands, nlev, ncol, do_delta_scaling
      real(c_double), intent(in) :: mu0(*), od(*), ssa(*), asymmetry(*)
      real(c_double), intent(in) :: albedo_direct(*), albedo_diffuse(*), incoming_sw(*)
      integer(c_int), intent(in) :: band_from_g(*)
      real(c_double), value :: cloud_fraction_threshold
      real(c_double), intent(in) :: cloud_fraction(*), total_cloud_cover(*)
      real(c_double), intent(in) :: od_scaling(*)
      real(c_double), intent(in) :: od_cloud(*), ssa_cloud(*), asymmetry_cloud(*)
      real(c_double), intent(out) :: sw_up_clear(*), sw_dn_clear(*)
      real(c_double), intent(out) :: sw_dn_direct_clear(*), sw_up(*), sw_dn(*)
      real(c_double), intent(out) :: sw_dn_direct(*)
      real(c_double), intent(out) :: sw_dn_diffuse_surf_clear_g(*)
      real(c_double), intent(out) :: sw_dn_direct_surf_clear_g(*)
      real(c_double), intent(out) :: sw_dn_diffuse_surf_g(*), sw_dn_direct_surf_g(*)
    end function oifs_cuda_sw_compute_dp

    integer(c_int) function oifs_cuda_lw_compute_dp(ng,nbands,nlev,ncol, &
         & do_aerosol_scattering,do_cloud_scattering,do_derivatives, &
         & od,ssa,asymmetry,planck_hl,emission,albedo,band_from_g, &
         & cloud_fraction_threshold,cloud_fraction,total_cloud_cover, &
         & od_scaling,od_cloud,ssa_cloud,asymmetry_cloud,lw_up_clear, &
         & lw_dn_clear,lw_up,lw_dn,lw_dn_surf_clear_g,lw_dn_surf_g, &
         & lw_derivatives) bind(C,name='oifs_cuda_bridge_lw_compute_dp')
      import :: c_double, c_int
      integer(c_int), value :: ng, nbands, nlev, ncol
      integer(c_int), value :: do_aerosol_scattering, do_cloud_scattering
      integer(c_int), value :: do_derivatives
      real(c_double), intent(in) :: od(*), ssa(*), asymmetry(*), planck_hl(*)
      real(c_double), intent(in) :: emission(*), albedo(*)
      integer(c_int), intent(in) :: band_from_g(*)
      real(c_double), value :: cloud_fraction_threshold
      real(c_double), intent(in) :: cloud_fraction(*), total_cloud_cover(*)
      real(c_double), intent(in) :: od_scaling(*), od_cloud(*)
      real(c_double), intent(in) :: ssa_cloud(*), asymmetry_cloud(*)
      real(c_double), intent(out) :: lw_up_clear(*), lw_dn_clear(*)
      real(c_double), intent(out) :: lw_up(*), lw_dn(*)
      real(c_double), intent(out) :: lw_dn_surf_clear_g(*), lw_dn_surf_g(*)
      real(c_double), intent(out) :: lw_derivatives(*)
    end function oifs_cuda_lw_compute_dp
  end interface
#endif

contains

  logical function cuda_radiation_available()
#if defined(OIFS_CUDA_RADIATION) && !defined(PARKIND1_SINGLE)
    cuda_radiation_available = oifs_cuda_radiation_available() /= 0_c_int
#else
    cuda_radiation_available = .false.
#endif
  end function cuda_radiation_available

  integer function cuda_sw_compute(ng,nbands,nlev,ncol,do_delta_scaling, &
       & mu0,od,ssa,asymmetry,albedo_direct,albedo_diffuse,incoming_sw, &
       & band_from_g,cloud_fraction_threshold,cloud_fraction,total_cloud_cover, &
       & od_scaling,od_cloud,ssa_cloud,asymmetry_cloud,sw_up_clear, &
       & sw_dn_clear,sw_dn_direct_clear,sw_up,sw_dn,sw_dn_direct, &
       & sw_dn_diffuse_surf_clear_g,sw_dn_direct_surf_clear_g, &
       & sw_dn_diffuse_surf_g,sw_dn_direct_surf_g)
    integer, intent(in) :: ng, nbands, nlev, ncol
    logical, intent(in) :: do_delta_scaling
    real(jprb), intent(in) :: mu0(ncol)
    real(jprb), intent(in) :: od(ng,nlev,ncol), ssa(ng,nlev,ncol)
    real(jprb), intent(in) :: asymmetry(ng,nlev,ncol)
    real(jprb), intent(in) :: albedo_direct(ng,ncol), albedo_diffuse(ng,ncol)
    real(jprb), intent(in) :: incoming_sw(ng,ncol)
    integer, intent(in) :: band_from_g(ng)
    real(jprb), intent(in) :: cloud_fraction_threshold
    real(jprb), intent(in) :: cloud_fraction(ncol,nlev)
    real(jprb), intent(in) :: total_cloud_cover(ncol)
    real(jprb), intent(in) :: od_scaling(ng,nlev,ncol)
    real(jprb), intent(in) :: od_cloud(nbands,nlev,ncol)
    real(jprb), intent(in) :: ssa_cloud(nbands,nlev,ncol)
    real(jprb), intent(in) :: asymmetry_cloud(nbands,nlev,ncol)
    real(jprb), intent(out) :: sw_up_clear(ncol,nlev+1), sw_dn_clear(ncol,nlev+1)
    real(jprb), intent(out) :: sw_dn_direct_clear(ncol,nlev+1)
    real(jprb), intent(out) :: sw_up(ncol,nlev+1), sw_dn(ncol,nlev+1)
    real(jprb), intent(out) :: sw_dn_direct(ncol,nlev+1)
    real(jprb), intent(out) :: sw_dn_diffuse_surf_clear_g(ng,ncol)
    real(jprb), intent(out) :: sw_dn_direct_surf_clear_g(ng,ncol)
    real(jprb), intent(out) :: sw_dn_diffuse_surf_g(ng,ncol)
    real(jprb), intent(out) :: sw_dn_direct_surf_g(ng,ncol)

#if defined(OIFS_CUDA_RADIATION) && !defined(PARKIND1_SINGLE)
    cuda_sw_compute = oifs_cuda_sw_compute_dp(int(ng,c_int),int(nbands,c_int), &
         & int(nlev,c_int),int(ncol,c_int),merge(1_c_int,0_c_int,do_delta_scaling), &
         & mu0,od,ssa,asymmetry,albedo_direct,albedo_diffuse,incoming_sw, &
         & band_from_g,real(cloud_fraction_threshold,c_double),cloud_fraction, &
         & total_cloud_cover,od_scaling,od_cloud,ssa_cloud,asymmetry_cloud, &
         & sw_up_clear,sw_dn_clear,sw_dn_direct_clear,sw_up,sw_dn,sw_dn_direct, &
         & sw_dn_diffuse_surf_clear_g,sw_dn_direct_surf_clear_g, &
         & sw_dn_diffuse_surf_g,sw_dn_direct_surf_g)
#else
    cuda_sw_compute = -1
#endif
  end function cuda_sw_compute

  integer function cuda_lw_compute(ng,nbands,nlev,ncol, &
       & do_aerosol_scattering,do_cloud_scattering,do_derivatives, &
       & od,ssa,asymmetry,planck_hl,emission,albedo,band_from_g, &
       & cloud_fraction_threshold,cloud_fraction,total_cloud_cover, &
       & od_scaling,od_cloud,ssa_cloud,asymmetry_cloud,lw_up_clear, &
       & lw_dn_clear,lw_up,lw_dn,lw_dn_surf_clear_g,lw_dn_surf_g, &
       & lw_derivatives)
    integer, intent(in) :: ng, nbands, nlev, ncol
    logical, intent(in) :: do_aerosol_scattering, do_cloud_scattering
    logical, intent(in) :: do_derivatives
    real(jprb), intent(in) :: od(ng,nlev,ncol)
    real(jprb), intent(in) :: ssa(:,:,:), asymmetry(:,:,:)
    real(jprb), intent(in) :: planck_hl(ng,nlev+1,ncol)
    real(jprb), intent(in) :: emission(ng,ncol), albedo(ng,ncol)
    integer, intent(in) :: band_from_g(ng)
    real(jprb), intent(in) :: cloud_fraction_threshold
    real(jprb), intent(in) :: cloud_fraction(ncol,nlev)
    real(jprb), intent(in) :: total_cloud_cover(ncol)
    real(jprb), intent(in) :: od_scaling(ng,nlev,ncol)
    real(jprb), intent(in) :: od_cloud(nbands,nlev,ncol)
    real(jprb), intent(in) :: ssa_cloud(:,:,:), asymmetry_cloud(:,:,:)
    real(jprb), intent(out) :: lw_up_clear(ncol,nlev+1)
    real(jprb), intent(out) :: lw_dn_clear(ncol,nlev+1)
    real(jprb), intent(out) :: lw_up(ncol,nlev+1), lw_dn(ncol,nlev+1)
    real(jprb), intent(out) :: lw_dn_surf_clear_g(ng,ncol)
    real(jprb), intent(out) :: lw_dn_surf_g(ng,ncol)
    real(jprb), intent(out) :: lw_derivatives(ncol,nlev+1)

#if defined(OIFS_CUDA_RADIATION) && !defined(PARKIND1_SINGLE)
    cuda_lw_compute = oifs_cuda_lw_compute_dp(int(ng,c_int),int(nbands,c_int), &
         & int(nlev,c_int),int(ncol,c_int), &
         & merge(1_c_int,0_c_int,do_aerosol_scattering), &
         & merge(1_c_int,0_c_int,do_cloud_scattering), &
         & merge(1_c_int,0_c_int,do_derivatives),od,ssa,asymmetry,planck_hl, &
         & emission,albedo,band_from_g,real(cloud_fraction_threshold,c_double), &
         & cloud_fraction,total_cloud_cover,od_scaling,od_cloud,ssa_cloud, &
         & asymmetry_cloud,lw_up_clear,lw_dn_clear,lw_up,lw_dn, &
         & lw_dn_surf_clear_g,lw_dn_surf_g,lw_derivatives)
#else
    cuda_lw_compute = -1
#endif
  end function cuda_lw_compute

end module radiation_cuda_bridge
