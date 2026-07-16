! radiation_mcica_sw.F90 - Monte-Carlo Independent Column Approximation shortwave solver
!
! Copyright (C) 2015-2017 ECMWF
!
! Author:  Robin Hogan
! Email:   r.j.hogan@ecmwf.int
! License: see the COPYING file for details
!
! Modifications
!   2017-04-11  R. Hogan  Receive albedos at g-points
!   2017-04-22  R. Hogan  Store surface fluxes at all g-points
!   2017-10-23  R. Hogan  Renamed single-character variables

module radiation_mcica_sw

  implicit none

#ifdef OIFS_CUDA_RADIATION
  logical, save :: gpu_options_initialized = .false.
  logical, save :: gpu_enabled_option = .false.
  logical, save :: gpu_available_option = .false.
  logical, save :: gpu_validate_option = .false.
  logical, save :: gpu_cloud_enabled_option = .false.
  integer, save :: gpu_min_columns_option = 256
#endif

contains

  ! Provides elemental function "delta_eddington"
#include "radiation_delta_eddington.h"

  !---------------------------------------------------------------------
  ! Shortwave Monte Carlo Independent Column Approximation
  ! (McICA). This implementation performs a clear-sky and a cloudy-sky
  ! calculation, and then weights the two to get the all-sky fluxes
  ! according to the total cloud cover. This method reduces noise for
  ! low cloud cover situations, and exploits the clear-sky
  ! calculations that are usually performed for diagnostic purposes
  ! simultaneously. The cloud generator has been carefully written
  ! such that the stochastic cloud field satisfies the prescribed
  ! overlap parameter accounting for this weighting.
  subroutine solver_mcica_sw(nlev,istartcol,iendcol, &
       &  config, single_level, cloud, & 
       &  od, ssa, g, od_cloud, ssa_cloud, g_cloud, &
       &  albedo_direct, albedo_diffuse, incoming_sw, &
       &  flux)

    use parkind1, only           : jprb
    use yomhook,  only           : lhook, dr_hook, jphook
    use radiation_io,   only           : nulerr, radiation_abort
    use radiation_config, only         : config_type
    use radiation_single_level, only   : single_level_type
    use radiation_cloud, only          : cloud_type
    use radiation_flux, only           : flux_type
    use radiation_two_stream, only     : calc_two_stream_gammas_sw, &
         &                               calc_reflectance_transmittance_sw
    use radiation_adding_ica_sw, only  : adding_ica_sw
    use radiation_cloud_generator, only: cloud_generator
#ifdef OIFS_CUDA_RADIATION
    use radiation_cuda_bridge, only   : cuda_radiation_available, cuda_sw_compute, &
         &                              cuda_cloud_compute
#endif

    implicit none

    ! Inputs
    integer, intent(in) :: nlev               ! number of model levels
    integer, intent(in) :: istartcol, iendcol ! range of columns to process
    type(config_type),        intent(in) :: config
    type(single_level_type),  intent(in) :: single_level
    type(cloud_type),         intent(in) :: cloud

    ! Gas and aerosol optical depth, single-scattering albedo and
    ! asymmetry factor at each shortwave g-point
    real(jprb), intent(in), dimension(config%n_g_sw, nlev, istartcol:iendcol) :: &
         &  od, ssa, g

    ! Cloud and precipitation optical depth, single-scattering albedo and
    ! asymmetry factor in each shortwave band
    real(jprb), intent(in), dimension(config%n_bands_sw,nlev,istartcol:iendcol)   :: &
         &  od_cloud, ssa_cloud, g_cloud

    ! Direct and diffuse surface albedos, and the incoming shortwave
    ! flux into a plane perpendicular to the incoming radiation at
    ! top-of-atmosphere in each of the shortwave g points
    real(jprb), intent(in), dimension(config%n_g_sw,istartcol:iendcol) :: &
         &  albedo_direct, albedo_diffuse, incoming_sw

    ! Output
    type(flux_type), intent(inout):: flux

    ! Local variables

    ! Cosine of solar zenith angle
    real(jprb)                                 :: cos_sza

    ! Diffuse reflectance and transmittance for each layer in clear
    ! and all skies
    real(jprb), dimension(config%n_g_sw, nlev) :: ref_clear, trans_clear, reflectance, transmittance

    ! Fraction of direct beam scattered by a layer into the upwelling
    ! or downwelling diffuse streams, in clear and all skies
    real(jprb), dimension(config%n_g_sw, nlev) :: ref_dir_clear, trans_dir_diff_clear, ref_dir, trans_dir_diff

    ! Transmittance for the direct beam in clear and all skies
    real(jprb), dimension(config%n_g_sw, nlev) :: trans_dir_dir_clear, trans_dir_dir

#ifdef OIFS_CUDA_RADIATION
    real(jprb), allocatable, dimension(:,:,:) :: od_scaling_batch
    real(jprb), allocatable, dimension(:) :: total_cloud_cover_batch
    real(jprb), allocatable, dimension(:,:) :: gpu_sw_up_clear, gpu_sw_dn_clear
    real(jprb), allocatable, dimension(:,:) :: gpu_sw_dn_direct_clear
    real(jprb), allocatable, dimension(:,:) :: gpu_sw_up, gpu_sw_dn, gpu_sw_dn_direct
    real(jprb), allocatable, dimension(:,:) :: gpu_sw_dn_diffuse_surf_clear_g
    real(jprb), allocatable, dimension(:,:) :: gpu_sw_dn_direct_surf_clear_g
    real(jprb), allocatable, dimension(:,:) :: gpu_sw_dn_diffuse_surf_g
    real(jprb), allocatable, dimension(:,:) :: gpu_sw_dn_direct_surf_g
#endif

    ! Fluxes per g point
    real(jprb), dimension(config%n_g_sw, nlev+1) :: flux_up, flux_dn_diffuse, flux_dn_direct

    ! Combined gas+aerosol+cloud optical depth, single scattering
    ! albedo and asymmetry factor
    real(jprb), dimension(config%n_g_sw) :: od_total, ssa_total, g_total

    ! Two-stream coefficients
    real(jprb), dimension(config%n_g_sw) :: gamma1, gamma2, gamma3

    ! Optical depth scaling from the cloud generator, zero indicating
    ! clear skies
    real(jprb), dimension(config%n_g_sw,nlev) :: od_scaling

    ! Modified optical depth after McICA scaling to represent cloud
    ! inhomogeneity
    real(jprb), dimension(config%n_g_sw) :: od_cloud_new

    ! Total cloud cover output from the cloud generator
    real(jprb) :: total_cloud_cover

    ! Number of g points
    integer :: ng

    ! Loop indices for level and column
    integer :: jlev, jcol

#ifdef OIFS_CUDA_RADIATION
    logical :: use_gpu_sw, validate_gpu_sw
    character(len=32) :: gpu_radiation_env, gpu_min_columns_env, gpu_validate_env
    character(len=32) :: gpu_cloud_env
    integer :: gpu_min_columns, gpu_status, env_read_status
    real(jprb) :: gpu_error_ratio, gpu_cloud_error_ratio
#endif

    real(jphook) :: hook_handle

    if (lhook) call dr_hook('radiation_mcica_sw:solver_mcica_sw',0,hook_handle)

    if (.not. config%do_clear) then
      write(nulerr,'(a)') '*** Error: shortwave McICA requires clear-sky calculation to be performed'
      call radiation_abort()      
    end if

    ng = config%n_g_sw

#ifdef OIFS_CUDA_RADIATION
!$omp critical(oifs_cuda_radiation_env)
    if (.not. gpu_options_initialized) then
      gpu_radiation_env = ''
      call get_environment_variable('OIFS_GPU_RADIATION',gpu_radiation_env)
      gpu_enabled_option = trim(gpu_radiation_env) == '1'
      gpu_validate_env = ''
      call get_environment_variable('OIFS_GPU_VALIDATE',gpu_validate_env)
      gpu_validate_option = trim(gpu_validate_env) == '1'
      gpu_cloud_env = ''
      call get_environment_variable('OIFS_GPU_CLOUD',gpu_cloud_env)
      gpu_cloud_enabled_option = trim(gpu_cloud_env) == '1'
      gpu_min_columns_option = 256
      gpu_min_columns_env = ''
      call get_environment_variable('OIFS_GPU_MIN_COLUMNS',gpu_min_columns_env)
      if (len_trim(gpu_min_columns_env) > 0) then
        read(gpu_min_columns_env,*,iostat=env_read_status) gpu_min_columns
        if (env_read_status == 0) gpu_min_columns_option = max(1,gpu_min_columns)
      end if
      if (gpu_enabled_option) then
        gpu_available_option = cuda_radiation_available()
        if (.not. gpu_available_option) write(nulerr,'(a)') &
             & '*** CUDA shortwave radiation unavailable; using the CPU solver'
      end if
      gpu_options_initialized = .true.
    end if
    use_gpu_sw = gpu_enabled_option .and. gpu_available_option
    validate_gpu_sw = gpu_validate_option
    gpu_min_columns = gpu_min_columns_option
!$omp end critical(oifs_cuda_radiation_env)
    use_gpu_sw = use_gpu_sw .and. iendcol-istartcol+1 >= max(1,gpu_min_columns)

    if (use_gpu_sw) then
      allocate(od_scaling_batch(ng,nlev,istartcol:iendcol), &
           & total_cloud_cover_batch(istartcol:iendcol), &
           & gpu_sw_up_clear(istartcol:iendcol,nlev+1), &
           & gpu_sw_dn_clear(istartcol:iendcol,nlev+1), &
           & gpu_sw_dn_direct_clear(istartcol:iendcol,nlev+1), &
           & gpu_sw_up(istartcol:iendcol,nlev+1), &
           & gpu_sw_dn(istartcol:iendcol,nlev+1), &
           & gpu_sw_dn_direct(istartcol:iendcol,nlev+1), &
           & gpu_sw_dn_diffuse_surf_clear_g(ng,istartcol:iendcol), &
           & gpu_sw_dn_direct_surf_clear_g(ng,istartcol:iendcol), &
           & gpu_sw_dn_diffuse_surf_g(ng,istartcol:iendcol), &
           & gpu_sw_dn_direct_surf_g(ng,istartcol:iendcol))
      gpu_status = 7
      if (gpu_cloud_enabled_option) then
        gpu_status = cuda_cloud_compute(ng,nlev,iendcol-istartcol+1, &
             & config%i_overlap_scheme,config%use_beta_overlap, &
             & single_level%iseed(istartcol:iendcol), &
             & single_level%cos_sza(istartcol:iendcol), &
             & config%cloud_fraction_threshold,cloud%fraction(istartcol:iendcol,:), &
             & cloud%overlap_param(istartcol:iendcol,:), &
             & config%cloud_inhom_decorr_scaling, &
             & cloud%fractional_std(istartcol:iendcol,:),config%pdf_sampler%ncdf, &
             & config%pdf_sampler%nfsd,config%pdf_sampler%fsd1, &
             & config%pdf_sampler%inv_fsd_interval,config%pdf_sampler%val, &
             & od_scaling_batch,total_cloud_cover_batch)
      end if
      if (gpu_status == 7) then
        total_cloud_cover_batch = 0.0_jprb
        do jcol = istartcol,iendcol
          if (single_level%cos_sza(jcol) > 0.0_jprb) then
            call cloud_generator(ng,nlev,config%i_overlap_scheme, &
                 & single_level%iseed(jcol),config%cloud_fraction_threshold, &
                 & cloud%fraction(jcol,:),cloud%overlap_param(jcol,:), &
                 & config%cloud_inhom_decorr_scaling,cloud%fractional_std(jcol,:), &
                 & config%pdf_sampler,od_scaling_batch(:,:,jcol), &
                 & total_cloud_cover_batch(jcol),is_beta_overlap=config%use_beta_overlap)
          end if
        end do
      else if (gpu_status /= 0) then
        write(nulerr,'(a,i0)') '*** CUDA cloud generation failed with status ',gpu_status
        call radiation_abort()
      end if
      gpu_status = cuda_sw_compute(ng,config%n_bands_sw,nlev,iendcol-istartcol+1, &
           & config%do_sw_delta_scaling_with_gases, &
           & single_level%cos_sza(istartcol:iendcol),od,ssa,g,albedo_direct, &
           & albedo_diffuse,incoming_sw,config%i_band_from_reordered_g_sw, &
           & config%cloud_fraction_threshold,cloud%fraction(istartcol:iendcol,:), &
           & total_cloud_cover_batch,od_scaling_batch,od_cloud,ssa_cloud,g_cloud, &
           & gpu_sw_up_clear,gpu_sw_dn_clear,gpu_sw_dn_direct_clear, &
           & gpu_sw_up,gpu_sw_dn,gpu_sw_dn_direct, &
           & gpu_sw_dn_diffuse_surf_clear_g,gpu_sw_dn_direct_surf_clear_g, &
           & gpu_sw_dn_diffuse_surf_g,gpu_sw_dn_direct_surf_g)
      if (gpu_status /= 0) then
        write(nulerr,'(a,i0)') '*** CUDA shortwave radiation failed with status ',gpu_status
        call radiation_abort()
      end if
      if (.not. validate_gpu_sw) then
        flux%sw_up_clear(istartcol:iendcol,:) = gpu_sw_up_clear
        flux%sw_dn_clear(istartcol:iendcol,:) = gpu_sw_dn_clear
        flux%sw_up(istartcol:iendcol,:) = gpu_sw_up
        flux%sw_dn(istartcol:iendcol,:) = gpu_sw_dn
        if (allocated(flux%sw_dn_direct_clear)) then
          flux%sw_dn_direct_clear(istartcol:iendcol,:) = gpu_sw_dn_direct_clear
        end if
        if (allocated(flux%sw_dn_direct)) then
          flux%sw_dn_direct(istartcol:iendcol,:) = gpu_sw_dn_direct
        end if
        flux%sw_dn_diffuse_surf_clear_g(:,istartcol:iendcol) = gpu_sw_dn_diffuse_surf_clear_g
        flux%sw_dn_direct_surf_clear_g(:,istartcol:iendcol) = gpu_sw_dn_direct_surf_clear_g
        flux%sw_dn_diffuse_surf_g(:,istartcol:iendcol) = gpu_sw_dn_diffuse_surf_g
        flux%sw_dn_direct_surf_g(:,istartcol:iendcol) = gpu_sw_dn_direct_surf_g
        flux%cloud_cover_sw(istartcol:iendcol) = total_cloud_cover_batch
        if (lhook) call dr_hook('radiation_mcica_sw:solver_mcica_sw',1,hook_handle)
        return
      end if
    end if
#endif

    ! Loop through columns
#ifdef OIFS_CUDA_RADIATION
    gpu_cloud_error_ratio = 0.0_jprb
#endif
    do jcol = istartcol,iendcol
      ! Only perform calculation if sun above the horizon
      if (single_level%cos_sza(jcol) > 0.0_jprb) then
        cos_sza = single_level%cos_sza(jcol)

        ! Clear-sky calculation - first compute clear-sky reflectance,
        ! transmittance etc at each model level
        if (.not. config%do_sw_delta_scaling_with_gases) then
          ! Delta-Eddington scaling has already been performed to the
          ! aerosol part of od, ssa and g
          do jlev = 1,nlev
            call calc_two_stream_gammas_sw(ng, &
                 &  cos_sza, ssa(:,jlev,jcol), g(:,jlev,jcol), &
                 &  gamma1, gamma2, gamma3)
            call calc_reflectance_transmittance_sw(ng, &
                 &  cos_sza, od(:,jlev,jcol), ssa(:,jlev,jcol), &
                 &  gamma1, gamma2, gamma3, &
                 &  ref_clear(:,jlev), trans_clear(:,jlev), &
                 &  ref_dir_clear(:,jlev), trans_dir_diff_clear(:,jlev), &
                 &  trans_dir_dir_clear(:,jlev) )
          end do
        else
          ! Apply delta-Eddington scaling to the aerosol-gas mixture
          do jlev = 1,nlev
            od_total  =  od(:,jlev,jcol)
            ssa_total = ssa(:,jlev,jcol)
            g_total   =   g(:,jlev,jcol)
            call delta_eddington(od_total, ssa_total, g_total)
            call calc_two_stream_gammas_sw(ng, &
                 &  cos_sza, ssa_total, g_total, &
                 &  gamma1, gamma2, gamma3)
            call calc_reflectance_transmittance_sw(ng, &
                 &  cos_sza, od_total, ssa_total, &
                 &  gamma1, gamma2, gamma3, &
                 &  ref_clear(:,jlev), trans_clear(:,jlev), &
                 &  ref_dir_clear(:,jlev), trans_dir_diff_clear(:,jlev), &
                 &  trans_dir_dir_clear(:,jlev) )
          end do
        end if

        ! Use adding method to compute fluxes
        call adding_ica_sw(ng, nlev, incoming_sw(:,jcol), &
             &  albedo_diffuse(:,jcol), albedo_direct(:,jcol), spread(cos_sza,1,ng), &
             &  ref_clear, trans_clear, ref_dir_clear, trans_dir_diff_clear, &
             &  trans_dir_dir_clear, flux_up, flux_dn_diffuse, flux_dn_direct)
        
        ! Sum over g-points to compute and save clear-sky broadband
        ! fluxes
        flux%sw_up_clear(jcol,:) = sum(flux_up,1)
        if (allocated(flux%sw_dn_direct_clear)) then
          flux%sw_dn_direct_clear(jcol,:) &
               &  = sum(flux_dn_direct,1)
          flux%sw_dn_clear(jcol,:) = sum(flux_dn_diffuse,1) &
               &  + flux%sw_dn_direct_clear(jcol,:)
        else
          flux%sw_dn_clear(jcol,:) = sum(flux_dn_diffuse,1) &
               &  + sum(flux_dn_direct,1)
        end if
        ! Store spectral downwelling fluxes at surface
        flux%sw_dn_diffuse_surf_clear_g(:,jcol) = flux_dn_diffuse(:,nlev+1)
        flux%sw_dn_direct_surf_clear_g(:,jcol)  = flux_dn_direct(:,nlev+1)

        ! Do cloudy-sky calculation
        call cloud_generator(ng, nlev, config%i_overlap_scheme, &
             &  single_level%iseed(jcol), &
             &  config%cloud_fraction_threshold, &
             &  cloud%fraction(jcol,:), cloud%overlap_param(jcol,:), &
             &  config%cloud_inhom_decorr_scaling, cloud%fractional_std(jcol,:), &
             &  config%pdf_sampler, od_scaling, total_cloud_cover, &
             &  is_beta_overlap=config%use_beta_overlap)

#ifdef OIFS_CUDA_RADIATION
        if (use_gpu_sw .and. validate_gpu_sw) then
          if (total_cloud_cover >= config%cloud_fraction_threshold) then
            gpu_cloud_error_ratio = max(gpu_cloud_error_ratio,maxval(abs(od_scaling &
                 & -od_scaling_batch(:,:,jcol)) / (1.0e-13_jprb+2.0e-12_jprb &
                 & *max(abs(od_scaling),abs(od_scaling_batch(:,:,jcol))))))
          end if
          gpu_cloud_error_ratio = max(gpu_cloud_error_ratio,abs(total_cloud_cover &
               & -total_cloud_cover_batch(jcol)) / (1.0e-13_jprb+2.0e-12_jprb &
               & *max(abs(total_cloud_cover),abs(total_cloud_cover_batch(jcol)))))
        end if
#endif

        ! Store total cloud cover
        flux%cloud_cover_sw(jcol) = total_cloud_cover
        
        if (total_cloud_cover >= config%cloud_fraction_threshold) then
          ! Total-sky calculation
          do jlev = 1,nlev
            ! Compute combined gas+aerosol+cloud optical properties
            if (cloud%fraction(jcol,jlev) >= config%cloud_fraction_threshold) then
              od_cloud_new = od_scaling(:,jlev) &
                   &  * od_cloud(config%i_band_from_reordered_g_sw,jlev,jcol)
              od_total  = od(:,jlev,jcol) + od_cloud_new
              ssa_total = 0.0_jprb
              g_total   = 0.0_jprb
              where (od_total > 0.0_jprb)
                ssa_total = (ssa(:,jlev,jcol)*od(:,jlev,jcol) &
                     &     + ssa_cloud(config%i_band_from_reordered_g_sw,jlev,jcol) &
                     &     *  od_cloud_new) & 
                     &     / od_total
              end where
              where (ssa_total*od_total > 0.0_jprb)
                g_total = (g(:,jlev,jcol)*ssa(:,jlev,jcol)*od(:,jlev,jcol) &
                     &     +   g_cloud(config%i_band_from_reordered_g_sw,jlev,jcol) &
                     &     * ssa_cloud(config%i_band_from_reordered_g_sw,jlev,jcol) &
                     &     *  od_cloud_new) &
                     &     / (ssa_total*od_total)
              end where

              ! Apply delta-Eddington scaling to the cloud-aerosol-gas
              ! mixture
              if (config%do_sw_delta_scaling_with_gases) then
                call delta_eddington(od_total, ssa_total, g_total)
              end if

             ! Compute cloudy-sky reflectance, transmittance etc at
              ! each model level
              call calc_two_stream_gammas_sw(ng, &
                   &  cos_sza, ssa_total, g_total, &
                   &  gamma1, gamma2, gamma3)

              call calc_reflectance_transmittance_sw(ng, &
                   &  cos_sza, od_total, ssa_total, &
                   &  gamma1, gamma2, gamma3, &
                   &  reflectance(:,jlev), transmittance(:,jlev), &
                   &  ref_dir(:,jlev), trans_dir_diff(:,jlev), &
                   &  trans_dir_dir(:,jlev) )

            else
              ! Clear-sky layer: copy over clear-sky values
              reflectance(:,jlev) = ref_clear(:,jlev)
              transmittance(:,jlev) = trans_clear(:,jlev)
              ref_dir(:,jlev) = ref_dir_clear(:,jlev)
              trans_dir_diff(:,jlev) = trans_dir_diff_clear(:,jlev)
              trans_dir_dir(:,jlev) = trans_dir_dir_clear(:,jlev)
            end if
          end do
            
          ! Use adding method to compute fluxes for an overcast sky
          call adding_ica_sw(ng, nlev, incoming_sw(:,jcol), &
               &  albedo_diffuse(:,jcol), albedo_direct(:,jcol), spread(cos_sza,1,ng), &
               &  reflectance, transmittance, ref_dir, trans_dir_diff, &
               &  trans_dir_dir, flux_up, flux_dn_diffuse, flux_dn_direct)
          
          ! Store overcast broadband fluxes
          flux%sw_up(jcol,:) = sum(flux_up,1)
          if (allocated(flux%sw_dn_direct)) then
            flux%sw_dn_direct(jcol,:) = sum(flux_dn_direct,1)
            flux%sw_dn(jcol,:) = sum(flux_dn_diffuse,1) &
                 &  + flux%sw_dn_direct(jcol,:)
          else
            flux%sw_dn(jcol,:) = sum(flux_dn_diffuse,1) &
                 &  + sum(flux_dn_direct,1)
          end if

          ! Cloudy flux profiles currently assume completely overcast
          ! skies; perform weighted average with clear-sky profile
          flux%sw_up(jcol,:) =  total_cloud_cover *flux%sw_up(jcol,:) &
               &  + (1.0_jprb - total_cloud_cover)*flux%sw_up_clear(jcol,:)
          flux%sw_dn(jcol,:) =  total_cloud_cover *flux%sw_dn(jcol,:) &
               &  + (1.0_jprb - total_cloud_cover)*flux%sw_dn_clear(jcol,:)
          if (allocated(flux%sw_dn_direct)) then
            flux%sw_dn_direct(jcol,:) = total_cloud_cover *flux%sw_dn_direct(jcol,:) &
                 &  + (1.0_jprb - total_cloud_cover)*flux%sw_dn_direct_clear(jcol,:)
          end if
          ! Likewise for surface spectral fluxes
          flux%sw_dn_diffuse_surf_g(:,jcol) = flux_dn_diffuse(:,nlev+1)
          flux%sw_dn_direct_surf_g(:,jcol)  = flux_dn_direct(:,nlev+1)
          flux%sw_dn_diffuse_surf_g(:,jcol) = total_cloud_cover *flux%sw_dn_diffuse_surf_g(:,jcol) &
               &     + (1.0_jprb - total_cloud_cover)*flux%sw_dn_diffuse_surf_clear_g(:,jcol)
          flux%sw_dn_direct_surf_g(:,jcol) = total_cloud_cover *flux%sw_dn_direct_surf_g(:,jcol) &
               &     + (1.0_jprb - total_cloud_cover)*flux%sw_dn_direct_surf_clear_g(:,jcol)
          
        else
          ! No cloud in profile and clear-sky fluxes already
          ! calculated: copy them over
          flux%sw_up(jcol,:) = flux%sw_up_clear(jcol,:)
          flux%sw_dn(jcol,:) = flux%sw_dn_clear(jcol,:)
          if (allocated(flux%sw_dn_direct)) then
            flux%sw_dn_direct(jcol,:) = flux%sw_dn_direct_clear(jcol,:)
          end if
          flux%sw_dn_diffuse_surf_g(:,jcol) = flux%sw_dn_diffuse_surf_clear_g(:,jcol)
          flux%sw_dn_direct_surf_g(:,jcol)  = flux%sw_dn_direct_surf_clear_g(:,jcol)

        end if ! Cloud is present in profile

      else
        ! Set fluxes to zero if sun is below the horizon
        flux%sw_up(jcol,:) = 0.0_jprb
        flux%sw_dn(jcol,:) = 0.0_jprb
        if (allocated(flux%sw_dn_direct)) then
          flux%sw_dn_direct(jcol,:) = 0.0_jprb
        end if
        flux%sw_up_clear(jcol,:) = 0.0_jprb
        flux%sw_dn_clear(jcol,:) = 0.0_jprb
        if (allocated(flux%sw_dn_direct_clear)) then
          flux%sw_dn_direct_clear(jcol,:) = 0.0_jprb
        end if
        flux%sw_dn_diffuse_surf_g(:,jcol) = 0.0_jprb
        flux%sw_dn_direct_surf_g(:,jcol)  = 0.0_jprb
        flux%sw_dn_diffuse_surf_clear_g(:,jcol) = 0.0_jprb
        flux%sw_dn_direct_surf_clear_g(:,jcol)  = 0.0_jprb
      end if ! Sun above horizon

    end do ! Loop over columns

#ifdef OIFS_CUDA_RADIATION
    if (use_gpu_sw .and. validate_gpu_sw) then
      if (gpu_cloud_error_ratio > 1.0_jprb) then
        write(nulerr,'(a,es12.4)') &
             & '*** CUDA shortwave cloud validation failed; normalized max error = ', &
             & gpu_cloud_error_ratio
        call radiation_abort()
      end if
      gpu_error_ratio = 0.0_jprb
      gpu_error_ratio = max(gpu_error_ratio,maxval(abs(flux%sw_up_clear(istartcol:iendcol,:)-gpu_sw_up_clear) &
           & / (5.0e-8_jprb + 2.0e-12_jprb*max(abs(flux%sw_up_clear(istartcol:iendcol,:)),abs(gpu_sw_up_clear)))))
      gpu_error_ratio = max(gpu_error_ratio,maxval(abs(flux%sw_dn_clear(istartcol:iendcol,:)-gpu_sw_dn_clear) &
           & / (5.0e-8_jprb + 2.0e-12_jprb*max(abs(flux%sw_dn_clear(istartcol:iendcol,:)),abs(gpu_sw_dn_clear)))))
      gpu_error_ratio = max(gpu_error_ratio,maxval(abs(flux%sw_up(istartcol:iendcol,:)-gpu_sw_up) &
           & / (5.0e-8_jprb + 2.0e-12_jprb*max(abs(flux%sw_up(istartcol:iendcol,:)),abs(gpu_sw_up)))))
      gpu_error_ratio = max(gpu_error_ratio,maxval(abs(flux%sw_dn(istartcol:iendcol,:)-gpu_sw_dn) &
           & / (5.0e-8_jprb + 2.0e-12_jprb*max(abs(flux%sw_dn(istartcol:iendcol,:)),abs(gpu_sw_dn)))))
      gpu_error_ratio = max(gpu_error_ratio,maxval(abs(flux%sw_dn_diffuse_surf_clear_g(:,istartcol:iendcol) &
           & - gpu_sw_dn_diffuse_surf_clear_g) / (5.0e-8_jprb + 2.0e-12_jprb &
           & * max(abs(flux%sw_dn_diffuse_surf_clear_g(:,istartcol:iendcol)),abs(gpu_sw_dn_diffuse_surf_clear_g)))))
      gpu_error_ratio = max(gpu_error_ratio,maxval(abs(flux%sw_dn_direct_surf_clear_g(:,istartcol:iendcol) &
           & - gpu_sw_dn_direct_surf_clear_g) / (5.0e-8_jprb + 2.0e-12_jprb &
           & * max(abs(flux%sw_dn_direct_surf_clear_g(:,istartcol:iendcol)),abs(gpu_sw_dn_direct_surf_clear_g)))))
      gpu_error_ratio = max(gpu_error_ratio,maxval(abs(flux%sw_dn_diffuse_surf_g(:,istartcol:iendcol) &
           & - gpu_sw_dn_diffuse_surf_g) / (5.0e-8_jprb + 2.0e-12_jprb &
           & * max(abs(flux%sw_dn_diffuse_surf_g(:,istartcol:iendcol)),abs(gpu_sw_dn_diffuse_surf_g)))))
      gpu_error_ratio = max(gpu_error_ratio,maxval(abs(flux%sw_dn_direct_surf_g(:,istartcol:iendcol) &
           & - gpu_sw_dn_direct_surf_g) / (5.0e-8_jprb + 2.0e-12_jprb &
           & * max(abs(flux%sw_dn_direct_surf_g(:,istartcol:iendcol)),abs(gpu_sw_dn_direct_surf_g)))))
      gpu_error_ratio = max(gpu_error_ratio,maxval(abs(flux%cloud_cover_sw(istartcol:iendcol)-total_cloud_cover_batch) &
           & / (5.0e-13_jprb + 2.0e-12_jprb*max(abs(flux%cloud_cover_sw(istartcol:iendcol)), &
           & abs(total_cloud_cover_batch)))))
      if (allocated(flux%sw_dn_direct_clear)) then
        gpu_error_ratio = max(gpu_error_ratio,maxval(abs(flux%sw_dn_direct_clear(istartcol:iendcol,:) &
             & - gpu_sw_dn_direct_clear) / (5.0e-8_jprb + 2.0e-12_jprb &
             & * max(abs(flux%sw_dn_direct_clear(istartcol:iendcol,:)),abs(gpu_sw_dn_direct_clear)))))
      end if
      if (allocated(flux%sw_dn_direct)) then
        gpu_error_ratio = max(gpu_error_ratio,maxval(abs(flux%sw_dn_direct(istartcol:iendcol,:)-gpu_sw_dn_direct) &
             & / (5.0e-8_jprb + 2.0e-12_jprb*max(abs(flux%sw_dn_direct(istartcol:iendcol,:)), &
             & abs(gpu_sw_dn_direct)))))
      end if
      if (gpu_error_ratio > 1.0_jprb) then
        write(nulerr,'(a,es12.4)') '*** CUDA shortwave validation failed; normalized max error = ',gpu_error_ratio
        write(nulerr,'(a,es12.4)') '    sw_up_clear max abs = ', &
             & maxval(abs(flux%sw_up_clear(istartcol:iendcol,:)-gpu_sw_up_clear))
        write(nulerr,'(a,es12.4)') '    sw_dn_clear max abs = ', &
             & maxval(abs(flux%sw_dn_clear(istartcol:iendcol,:)-gpu_sw_dn_clear))
        write(nulerr,'(a,es12.4)') '    sw_up max abs = ', &
             & maxval(abs(flux%sw_up(istartcol:iendcol,:)-gpu_sw_up))
        write(nulerr,'(a,es12.4)') '    sw_dn max abs = ', &
             & maxval(abs(flux%sw_dn(istartcol:iendcol,:)-gpu_sw_dn))
        write(nulerr,'(a,es12.4)') '    sw_dn_diffuse_surf_clear_g max abs = ', &
             & maxval(abs(flux%sw_dn_diffuse_surf_clear_g(:,istartcol:iendcol)-gpu_sw_dn_diffuse_surf_clear_g))
        write(nulerr,'(a,es12.4)') '    sw_dn_direct_surf_clear_g max abs = ', &
             & maxval(abs(flux%sw_dn_direct_surf_clear_g(:,istartcol:iendcol)-gpu_sw_dn_direct_surf_clear_g))
        write(nulerr,'(a,es12.4)') '    sw_dn_diffuse_surf_g max abs = ', &
             & maxval(abs(flux%sw_dn_diffuse_surf_g(:,istartcol:iendcol)-gpu_sw_dn_diffuse_surf_g))
        write(nulerr,'(a,es12.4)') '    sw_dn_direct_surf_g max abs = ', &
             & maxval(abs(flux%sw_dn_direct_surf_g(:,istartcol:iendcol)-gpu_sw_dn_direct_surf_g))
        if (allocated(flux%sw_dn_direct_clear)) write(nulerr,'(a,es12.4)') &
             & '    sw_dn_direct_clear max abs = ', &
             & maxval(abs(flux%sw_dn_direct_clear(istartcol:iendcol,:)-gpu_sw_dn_direct_clear))
        if (allocated(flux%sw_dn_direct)) write(nulerr,'(a,es12.4)') &
             & '    sw_dn_direct max abs = ', &
             & maxval(abs(flux%sw_dn_direct(istartcol:iendcol,:)-gpu_sw_dn_direct))
        call radiation_abort()
      end if

      flux%sw_up_clear(istartcol:iendcol,:) = gpu_sw_up_clear
      flux%sw_dn_clear(istartcol:iendcol,:) = gpu_sw_dn_clear
      flux%sw_up(istartcol:iendcol,:) = gpu_sw_up
      flux%sw_dn(istartcol:iendcol,:) = gpu_sw_dn
      if (allocated(flux%sw_dn_direct_clear)) &
           & flux%sw_dn_direct_clear(istartcol:iendcol,:) = gpu_sw_dn_direct_clear
      if (allocated(flux%sw_dn_direct)) flux%sw_dn_direct(istartcol:iendcol,:) = gpu_sw_dn_direct
      flux%sw_dn_diffuse_surf_clear_g(:,istartcol:iendcol) = gpu_sw_dn_diffuse_surf_clear_g
      flux%sw_dn_direct_surf_clear_g(:,istartcol:iendcol) = gpu_sw_dn_direct_surf_clear_g
      flux%sw_dn_diffuse_surf_g(:,istartcol:iendcol) = gpu_sw_dn_diffuse_surf_g
      flux%sw_dn_direct_surf_g(:,istartcol:iendcol) = gpu_sw_dn_direct_surf_g
      flux%cloud_cover_sw(istartcol:iendcol) = total_cloud_cover_batch
    end if
#endif

    if (lhook) call dr_hook('radiation_mcica_sw:solver_mcica_sw',1,hook_handle)
  end subroutine solver_mcica_sw

end module radiation_mcica_sw
