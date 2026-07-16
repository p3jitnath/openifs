! radiation_mcica_lw.F90 - Monte-Carlo Independent Column Approximation longtwave solver
!
! Copyright (C) 2015-2017 ECMWF
!
! Author:  Robin Hogan
! Email:   r.j.hogan@ecmwf.int
! License: see the COPYING file for details
!
! Modifications
!   2017-04-11  R. Hogan  Receive emission/albedo rather than planck/emissivity
!   2017-04-22  R. Hogan  Store surface fluxes at all g-points
!   2017-07-12  R. Hogan  Call fast adding method if only clouds scatter
!   2017-10-23  R. Hogan  Renamed single-character variables

module radiation_mcica_lw

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

  !---------------------------------------------------------------------
  ! Longwave Monte Carlo Independent Column Approximation
  ! (McICA). This implementation performs a clear-sky and a cloudy-sky
  ! calculation, and then weights the two to get the all-sky fluxes
  ! according to the total cloud cover. This method reduces noise for
  ! low cloud cover situations, and exploits the clear-sky
  ! calculations that are usually performed for diagnostic purposes
  ! simultaneously. The cloud generator has been carefully written
  ! such that the stochastic cloud field satisfies the prescribed
  ! overlap parameter accounting for this weighting.
  subroutine solver_mcica_lw(nlev,istartcol,iendcol, &
       &  config, single_level, cloud, & 
       &  od, ssa, g, od_cloud, ssa_cloud, g_cloud, planck_hl, &
       &  emission, albedo, &
       &  flux)

    use parkind1, only           : jprb
    use yomhook,  only           : lhook, dr_hook, jphook
    use radiation_io,   only           : nulerr, radiation_abort
    use radiation_config, only         : config_type
    use radiation_single_level, only   : single_level_type
    use radiation_cloud, only          : cloud_type
    use radiation_flux, only           : flux_type
    use radiation_two_stream, only     : calc_two_stream_gammas_lw, &
         &                               calc_reflectance_transmittance_lw, &
         &                               calc_no_scattering_transmittance_lw
    use radiation_adding_ica_lw, only  : adding_ica_lw, fast_adding_ica_lw, &
         &                               calc_fluxes_no_scattering_lw
    use radiation_lw_derivatives, only : calc_lw_derivatives_ica, modify_lw_derivatives_ica
    use radiation_cloud_generator, only: cloud_generator
#ifdef OIFS_CUDA_RADIATION
    use radiation_cuda_bridge, only   : cuda_radiation_available, cuda_lw_compute, &
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
    ! asymmetry factor at each longwave g-point
    real(jprb), intent(in), dimension(config%n_g_lw, nlev, istartcol:iendcol) :: &
         &  od
    real(jprb), intent(in), dimension(config%n_g_lw_if_scattering, nlev, istartcol:iendcol) :: &
         &  ssa, g

    ! Cloud and precipitation optical depth, single-scattering albedo and
    ! asymmetry factor in each longwave band
    real(jprb), intent(in), dimension(config%n_bands_lw,nlev,istartcol:iendcol)   :: &
         &  od_cloud
    real(jprb), intent(in), dimension(config%n_bands_lw_if_scattering, &
         &  nlev,istartcol:iendcol) :: ssa_cloud, g_cloud

    ! Planck function at each half-level and the surface
    real(jprb), intent(in), dimension(config%n_g_lw,nlev+1,istartcol:iendcol) :: &
         &  planck_hl

    ! Emission (Planck*emissivity) and albedo (1-emissivity) at the
    ! surface at each longwave g-point
    real(jprb), intent(in), dimension(config%n_g_lw, istartcol:iendcol) :: emission, albedo

    ! Output
    type(flux_type), intent(inout):: flux

    ! Local variables

    ! Diffuse reflectance and transmittance for each layer in clear
    ! and all skies
    real(jprb), dimension(config%n_g_lw, nlev) :: ref_clear, trans_clear, reflectance, transmittance

    ! Emission by a layer into the upwelling or downwelling diffuse
    ! streams, in clear and all skies
    real(jprb), dimension(config%n_g_lw, nlev) :: source_up_clear, source_dn_clear, source_up, source_dn

    ! Fluxes per g point
    real(jprb), dimension(config%n_g_lw, nlev+1) :: flux_up, flux_dn
    real(jprb), dimension(config%n_g_lw, nlev+1) :: flux_up_clear, flux_dn_clear

#ifdef OIFS_CUDA_RADIATION
    real(jprb), allocatable, dimension(:,:,:) :: od_scaling_batch
    real(jprb), allocatable, dimension(:) :: total_cloud_cover_batch
    real(jprb), allocatable, dimension(:) :: cloud_active_batch
    integer, allocatable, dimension(:) :: cloud_seed_batch
    real(jprb), allocatable, dimension(:,:) :: gpu_lw_up_clear, gpu_lw_dn_clear
    real(jprb), allocatable, dimension(:,:) :: gpu_lw_up, gpu_lw_dn
    real(jprb), allocatable, dimension(:,:) :: gpu_lw_dn_surf_clear_g
    real(jprb), allocatable, dimension(:,:) :: gpu_lw_dn_surf_g
    real(jprb), allocatable, dimension(:,:) :: gpu_lw_derivatives
#endif

    ! Combined gas+aerosol+cloud optical depth, single scattering
    ! albedo and asymmetry factor
    real(jprb), dimension(config%n_g_lw) :: od_total, ssa_total, g_total

    ! Two-stream coefficients
    real(jprb), dimension(config%n_g_lw) :: gamma1, gamma2

    ! Optical depth scaling from the cloud generator, zero indicating
    ! clear skies
    real(jprb), dimension(config%n_g_lw,nlev) :: od_scaling

    ! Modified optical depth after McICA scaling to represent cloud
    ! inhomogeneity
    real(jprb), dimension(config%n_g_lw) :: od_cloud_new

    ! Total cloud cover output from the cloud generator
    real(jprb) :: total_cloud_cover

    ! Identify clear-sky layers
    logical :: is_clear_sky_layer(nlev)

    ! Index of the highest cloudy layer
    integer :: i_cloud_top

    ! Number of g points
    integer :: ng

    ! Loop indices for level and column
    integer :: jlev, jcol

#ifdef OIFS_CUDA_RADIATION
    logical :: use_gpu_lw, validate_gpu_lw
    character(len=32) :: gpu_radiation_env, gpu_min_columns_env, gpu_validate_env
    character(len=32) :: gpu_cloud_env
    integer :: gpu_min_columns, gpu_status, env_read_status
    real(jprb) :: gpu_error_ratio, gpu_cloud_error_ratio
    real(jprb) :: gpu_cloud_scaling_max_abs, gpu_cloud_cover_max_abs
#endif

    real(jphook) :: hook_handle

    if (lhook) call dr_hook('radiation_mcica_lw:solver_mcica_lw',0,hook_handle)

    if (.not. config%do_clear) then
      write(nulerr,'(a)') '*** Error: longwave McICA requires clear-sky calculation to be performed'
      call radiation_abort()      
    end if

    ng = config%n_g_lw

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
             & '*** CUDA longwave radiation unavailable; using the CPU solver'
      end if
      gpu_options_initialized = .true.
    end if
    use_gpu_lw = gpu_enabled_option .and. gpu_available_option
    validate_gpu_lw = gpu_validate_option
    gpu_min_columns = gpu_min_columns_option
!$omp end critical(oifs_cuda_radiation_env)
    use_gpu_lw = use_gpu_lw .and. iendcol-istartcol+1 >= max(1,gpu_min_columns)

    if (use_gpu_lw) then
      allocate(od_scaling_batch(ng,nlev,istartcol:iendcol), &
           & total_cloud_cover_batch(istartcol:iendcol), &
           & cloud_active_batch(istartcol:iendcol), &
           & cloud_seed_batch(istartcol:iendcol), &
           & gpu_lw_up_clear(istartcol:iendcol,nlev+1), &
           & gpu_lw_dn_clear(istartcol:iendcol,nlev+1), &
           & gpu_lw_up(istartcol:iendcol,nlev+1), &
           & gpu_lw_dn(istartcol:iendcol,nlev+1), &
           & gpu_lw_dn_surf_clear_g(ng,istartcol:iendcol), &
           & gpu_lw_dn_surf_g(ng,istartcol:iendcol), &
           & gpu_lw_derivatives(istartcol:iendcol,nlev+1))
      cloud_active_batch = 1.0_jprb
      cloud_seed_batch = single_level%iseed(istartcol:iendcol)+997
      gpu_status = 7
      if (gpu_cloud_enabled_option) then
        gpu_status = cuda_cloud_compute(ng,nlev,iendcol-istartcol+1, &
             & config%i_overlap_scheme,config%use_beta_overlap,.not. validate_gpu_lw, &
             & cloud_seed_batch, &
             & cloud_active_batch,config%cloud_fraction_threshold, &
             & cloud%fraction(istartcol:iendcol,:), &
             & cloud%overlap_param(istartcol:iendcol,:), &
             & config%cloud_inhom_decorr_scaling, &
             & cloud%fractional_std(istartcol:iendcol,:),config%pdf_sampler%ncdf, &
             & config%pdf_sampler%nfsd,config%pdf_sampler%fsd1, &
             & config%pdf_sampler%inv_fsd_interval,config%pdf_sampler%val, &
             & od_scaling_batch,total_cloud_cover_batch)
      end if
      if (gpu_status == 7) then
        do jcol = istartcol,iendcol
          call cloud_generator(ng,nlev,config%i_overlap_scheme, &
               & single_level%iseed(jcol)+997,config%cloud_fraction_threshold, &
               & cloud%fraction(jcol,:),cloud%overlap_param(jcol,:), &
               & config%cloud_inhom_decorr_scaling,cloud%fractional_std(jcol,:), &
               & config%pdf_sampler,od_scaling_batch(:,:,jcol), &
               & total_cloud_cover_batch(jcol),is_beta_overlap=config%use_beta_overlap)
        end do
      else if (gpu_status /= 0) then
        write(nulerr,'(a,i0)') '*** CUDA cloud generation failed with status ',gpu_status
        call radiation_abort()
      end if
      gpu_status = cuda_lw_compute(ng,config%n_bands_lw,nlev,iendcol-istartcol+1, &
           & config%do_lw_aerosol_scattering,config%do_lw_cloud_scattering, &
           & config%do_lw_derivatives,od,ssa,g,planck_hl,emission,albedo, &
           & config%i_band_from_reordered_g_lw,config%cloud_fraction_threshold, &
           & cloud%fraction(istartcol:iendcol,:),total_cloud_cover_batch, &
           & od_scaling_batch,od_cloud,ssa_cloud,g_cloud,gpu_lw_up_clear, &
           & gpu_lw_dn_clear,gpu_lw_up,gpu_lw_dn,gpu_lw_dn_surf_clear_g, &
           & gpu_lw_dn_surf_g,gpu_lw_derivatives)
      if (gpu_status /= 0) then
        write(nulerr,'(a,i0)') '*** CUDA longwave radiation failed with status ',gpu_status
        call radiation_abort()
      end if
      if (.not. validate_gpu_lw) then
        flux%lw_up_clear(istartcol:iendcol,:) = gpu_lw_up_clear
        flux%lw_dn_clear(istartcol:iendcol,:) = gpu_lw_dn_clear
        flux%lw_up(istartcol:iendcol,:) = gpu_lw_up
        flux%lw_dn(istartcol:iendcol,:) = gpu_lw_dn
        flux%lw_dn_surf_clear_g(:,istartcol:iendcol) = gpu_lw_dn_surf_clear_g
        flux%lw_dn_surf_g(:,istartcol:iendcol) = gpu_lw_dn_surf_g
        flux%cloud_cover_lw(istartcol:iendcol) = total_cloud_cover_batch
        if (config%do_lw_derivatives) &
             & flux%lw_derivatives(istartcol:iendcol,:) = gpu_lw_derivatives
        if (lhook) call dr_hook('radiation_mcica_lw:solver_mcica_lw',1,hook_handle)
        return
      end if
    end if
#endif

    ! Loop through columns
#ifdef OIFS_CUDA_RADIATION
    gpu_cloud_error_ratio = 0.0_jprb
    gpu_cloud_scaling_max_abs = 0.0_jprb
    gpu_cloud_cover_max_abs = 0.0_jprb
#endif
    do jcol = istartcol,iendcol

      ! Clear-sky calculation
      if (config%do_lw_aerosol_scattering) then
        ! Scattering case: first compute clear-sky reflectance,
        ! transmittance etc at each model level
        do jlev = 1,nlev
          ssa_total = ssa(:,jlev,jcol)
          g_total   = g(:,jlev,jcol)
          call calc_two_stream_gammas_lw(ng, ssa_total, g_total, &
               &  gamma1, gamma2)
          call calc_reflectance_transmittance_lw(ng, &
               &  od(:,jlev,jcol), gamma1, gamma2, &
               &  planck_hl(:,jlev,jcol), planck_hl(:,jlev+1,jcol), &
               &  ref_clear(:,jlev), trans_clear(:,jlev), &
               &  source_up_clear(:,jlev), source_dn_clear(:,jlev))
        end do
        ! Then use adding method to compute fluxes
        call adding_ica_lw(ng, nlev, &
             &  ref_clear, trans_clear, source_up_clear, source_dn_clear, &
             &  emission(:,jcol), albedo(:,jcol), &
             &  flux_up_clear, flux_dn_clear)
        
      else
        ! Non-scattering case: use simpler functions for
        ! transmission and emission
        do jlev = 1,nlev
          call calc_no_scattering_transmittance_lw(ng, od(:,jlev,jcol), &
               &  planck_hl(:,jlev,jcol), planck_hl(:,jlev+1, jcol), &
               &  trans_clear(:,jlev), source_up_clear(:,jlev), source_dn_clear(:,jlev))
        end do
        ! Simpler down-then-up method to compute fluxes
        call calc_fluxes_no_scattering_lw(ng, nlev, &
             &  trans_clear, source_up_clear, source_dn_clear, &
             &  emission(:,jcol), albedo(:,jcol), &
             &  flux_up_clear, flux_dn_clear)
        
        ! Ensure that clear-sky reflectance is zero since it may be
        ! used in cloudy-sky case
        ref_clear = 0.0_jprb
      end if

      ! Sum over g-points to compute broadband fluxes
      flux%lw_up_clear(jcol,:) = sum(flux_up_clear,1)
      flux%lw_dn_clear(jcol,:) = sum(flux_dn_clear,1)
      ! Store surface spectral downwelling fluxes
      flux%lw_dn_surf_clear_g(:,jcol) = flux_dn_clear(:,nlev+1)

      ! Do cloudy-sky calculation; add a prime number to the seed in
      ! the longwave
      call cloud_generator(ng, nlev, config%i_overlap_scheme, &
           &  single_level%iseed(jcol) + 997, &
           &  config%cloud_fraction_threshold, &
           &  cloud%fraction(jcol,:), cloud%overlap_param(jcol,:), &
           &  config%cloud_inhom_decorr_scaling, cloud%fractional_std(jcol,:), &
           &  config%pdf_sampler, od_scaling, total_cloud_cover, &
           &  is_beta_overlap=config%use_beta_overlap)

#ifdef OIFS_CUDA_RADIATION
      if (use_gpu_lw .and. validate_gpu_lw) then
        if (total_cloud_cover >= config%cloud_fraction_threshold) then
          gpu_cloud_error_ratio = max(gpu_cloud_error_ratio,maxval(abs(od_scaling &
               & -od_scaling_batch(:,:,jcol)) / (1.0e-13_jprb+2.0e-12_jprb &
               & *max(abs(od_scaling),abs(od_scaling_batch(:,:,jcol))))))
          gpu_cloud_scaling_max_abs = max(gpu_cloud_scaling_max_abs, &
               & maxval(abs(od_scaling-od_scaling_batch(:,:,jcol))))
        end if
        gpu_cloud_error_ratio = max(gpu_cloud_error_ratio,abs(total_cloud_cover &
             & -total_cloud_cover_batch(jcol)) / (1.0e-13_jprb+2.0e-12_jprb &
             & *max(abs(total_cloud_cover),abs(total_cloud_cover_batch(jcol)))))
        gpu_cloud_cover_max_abs = max(gpu_cloud_cover_max_abs, &
             & abs(total_cloud_cover-total_cloud_cover_batch(jcol)))
      end if
#endif
      
      ! Store total cloud cover
      flux%cloud_cover_lw(jcol) = total_cloud_cover
      
      if (total_cloud_cover >= config%cloud_fraction_threshold) then
        ! Total-sky calculation

        is_clear_sky_layer = .true.
        i_cloud_top = nlev+1
        do jlev = 1,nlev
          ! Compute combined gas+aerosol+cloud optical properties
          if (cloud%fraction(jcol,jlev) >= config%cloud_fraction_threshold) then
            is_clear_sky_layer(jlev) = .false.
            ! Get index to the first cloudy layer from the top
            if (i_cloud_top > jlev) then
              i_cloud_top = jlev
            end if

            od_cloud_new = od_scaling(:,jlev) &
                 &  * od_cloud(config%i_band_from_reordered_g_lw,jlev,jcol)
            od_total = od(:,jlev,jcol) + od_cloud_new
            ssa_total = 0.0_jprb
            g_total   = 0.0_jprb

            if (config%do_lw_cloud_scattering) then
              ! Scattering case: calculate reflectance and
              ! transmittance at each model level
              if (config%do_lw_aerosol_scattering) then
                where (od_total > 0.0_jprb)
                  ssa_total = (ssa(:,jlev,jcol)*od(:,jlev,jcol) &
                       &     + ssa_cloud(config%i_band_from_reordered_g_lw,jlev,jcol) &
                       &     *  od_cloud_new) & 
                       &     / od_total
                end where
                where (ssa_total*od_total > 0.0_jprb)
                  g_total = (g(:,jlev,jcol)*ssa(:,jlev,jcol)*od(:,jlev,jcol) &
                       &     +   g_cloud(config%i_band_from_reordered_g_lw,jlev,jcol) &
                       &     * ssa_cloud(config%i_band_from_reordered_g_lw,jlev,jcol) &
                       &     *  od_cloud_new) &
                       &     / (ssa_total*od_total)
                end where
              else
                where (od_total > 0.0_jprb)
                  ssa_total = ssa_cloud(config%i_band_from_reordered_g_lw,jlev,jcol) &
                       &     * od_cloud_new / od_total
                end where
                where (ssa_total*od_total > 0.0_jprb)
                  g_total = g_cloud(config%i_band_from_reordered_g_lw,jlev,jcol) &
                       &     * ssa_cloud(config%i_band_from_reordered_g_lw,jlev,jcol) &
                       &     *  od_cloud_new / (ssa_total*od_total)
                end where
              end if
            
              ! Compute cloudy-sky reflectance, transmittance etc at
              ! each model level
              call calc_two_stream_gammas_lw(ng, ssa_total, g_total, &
                   &  gamma1, gamma2)
              call calc_reflectance_transmittance_lw(ng, &
                   &  od_total, gamma1, gamma2, &
                   &  planck_hl(:,jlev,jcol), planck_hl(:,jlev+1,jcol), &
                   &  reflectance(:,jlev), transmittance(:,jlev), source_up(:,jlev), source_dn(:,jlev))
            else
              ! No-scattering case: use simpler functions for
              ! transmission and emission
              call calc_no_scattering_transmittance_lw(ng, od_total, &
                   &  planck_hl(:,jlev,jcol), planck_hl(:,jlev+1, jcol), &
                   &  transmittance(:,jlev), source_up(:,jlev), source_dn(:,jlev))
            end if

          else
            ! Clear-sky layer: copy over clear-sky values
            reflectance(:,jlev) = ref_clear(:,jlev)
            transmittance(:,jlev) = trans_clear(:,jlev)
            source_up(:,jlev) = source_up_clear(:,jlev)
            source_dn(:,jlev) = source_dn_clear(:,jlev)
          end if
        end do
        
        if (config%do_lw_aerosol_scattering) then
          ! Use adding method to compute fluxes for an overcast sky,
          ! allowing for scattering in all layers
          call adding_ica_lw(ng, nlev, reflectance, transmittance, source_up, source_dn, &
               &  emission(:,jcol), albedo(:,jcol), &
               &  flux_up, flux_dn)
        else if (config%do_lw_cloud_scattering) then
          ! Use adding method to compute fluxes but optimize for the
          ! presence of clear-sky layers
!          call adding_ica_lw(ng, nlev, reflectance, transmittance, source_up, source_dn, &
!               &  emission(:,jcol), albedo(:,jcol), &
!               &  flux_up, flux_dn)
          call fast_adding_ica_lw(ng, nlev, reflectance, transmittance, source_up, source_dn, &
               &  emission(:,jcol), albedo(:,jcol), &
               &  is_clear_sky_layer, i_cloud_top, flux_dn_clear, &
               &  flux_up, flux_dn)
        else
          ! Simpler down-then-up method to compute fluxes
          call calc_fluxes_no_scattering_lw(ng, nlev, &
               &  transmittance, source_up, source_dn, emission(:,jcol), albedo(:,jcol), &
               &  flux_up, flux_dn)
        end if
        
        ! Store overcast broadband fluxes
        flux%lw_up(jcol,:) = sum(flux_up,1)
        flux%lw_dn(jcol,:) = sum(flux_dn,1)

        ! Cloudy flux profiles currently assume completely overcast
        ! skies; perform weighted average with clear-sky profile
        flux%lw_up(jcol,:) =  total_cloud_cover *flux%lw_up(jcol,:) &
             &  + (1.0_jprb - total_cloud_cover)*flux%lw_up_clear(jcol,:)
        flux%lw_dn(jcol,:) =  total_cloud_cover *flux%lw_dn(jcol,:) &
             &  + (1.0_jprb - total_cloud_cover)*flux%lw_dn_clear(jcol,:)
        ! Store surface spectral downwelling fluxes
        flux%lw_dn_surf_g(:,jcol) = total_cloud_cover*flux_dn(:,nlev+1) &
             &  + (1.0_jprb - total_cloud_cover)*flux%lw_dn_surf_clear_g(:,jcol)

        ! Compute the longwave derivatives needed by Hogan and Bozzo
        ! (2015) approximate radiation update scheme
        if (config%do_lw_derivatives) then
          call calc_lw_derivatives_ica(ng, nlev, jcol, transmittance, flux_up(:,nlev+1), &
               &                       flux%lw_derivatives)
          if (total_cloud_cover < 1.0_jprb - config%cloud_fraction_threshold) then
            ! Modify the existing derivative with the contribution from the clear sky
            call modify_lw_derivatives_ica(ng, nlev, jcol, trans_clear, flux_up_clear(:,nlev+1), &
                 &                         1.0_jprb-total_cloud_cover, flux%lw_derivatives)
          end if
        end if

      else
        ! No cloud in profile and clear-sky fluxes already
        ! calculated: copy them over
        flux%lw_up(jcol,:) = flux%lw_up_clear(jcol,:)
        flux%lw_dn(jcol,:) = flux%lw_dn_clear(jcol,:)
        flux%lw_dn_surf_g(:,jcol) = flux%lw_dn_surf_clear_g(:,jcol)
        if (config%do_lw_derivatives) then
          call calc_lw_derivatives_ica(ng, nlev, jcol, trans_clear, flux_up_clear(:,nlev+1), &
               &                       flux%lw_derivatives)
 
        end if
      end if ! Cloud is present in profile
    end do

#ifdef OIFS_CUDA_RADIATION
    if (use_gpu_lw .and. validate_gpu_lw) then
      if (gpu_cloud_error_ratio > 1.0_jprb) then
        write(nulerr,'(a,es12.4)') &
             & '*** CUDA longwave cloud validation failed; normalized max error = ', &
             & gpu_cloud_error_ratio
        write(nulerr,'(a,es12.4)') '    od_scaling max abs = ',gpu_cloud_scaling_max_abs
        write(nulerr,'(a,es12.4)') '    cloud cover max abs = ',gpu_cloud_cover_max_abs
        call radiation_abort()
      end if
      gpu_error_ratio = 0.0_jprb
      gpu_error_ratio = max(gpu_error_ratio,maxval(abs(flux%lw_up_clear(istartcol:iendcol,:)-gpu_lw_up_clear) &
           & / (1.0e-5_jprb + 2.0e-12_jprb*max(abs(flux%lw_up_clear(istartcol:iendcol,:)),abs(gpu_lw_up_clear)))))
      gpu_error_ratio = max(gpu_error_ratio,maxval(abs(flux%lw_dn_clear(istartcol:iendcol,:)-gpu_lw_dn_clear) &
           & / (1.0e-5_jprb + 2.0e-12_jprb*max(abs(flux%lw_dn_clear(istartcol:iendcol,:)),abs(gpu_lw_dn_clear)))))
      gpu_error_ratio = max(gpu_error_ratio,maxval(abs(flux%lw_up(istartcol:iendcol,:)-gpu_lw_up) &
           & / (1.0e-5_jprb + 2.0e-12_jprb*max(abs(flux%lw_up(istartcol:iendcol,:)),abs(gpu_lw_up)))))
      gpu_error_ratio = max(gpu_error_ratio,maxval(abs(flux%lw_dn(istartcol:iendcol,:)-gpu_lw_dn) &
           & / (1.0e-5_jprb + 2.0e-12_jprb*max(abs(flux%lw_dn(istartcol:iendcol,:)),abs(gpu_lw_dn)))))
      gpu_error_ratio = max(gpu_error_ratio,maxval(abs(flux%lw_dn_surf_clear_g(:,istartcol:iendcol) &
           & -gpu_lw_dn_surf_clear_g) / (1.0e-5_jprb + 2.0e-12_jprb &
           & *max(abs(flux%lw_dn_surf_clear_g(:,istartcol:iendcol)),abs(gpu_lw_dn_surf_clear_g)))))
      gpu_error_ratio = max(gpu_error_ratio,maxval(abs(flux%lw_dn_surf_g(:,istartcol:iendcol) &
           & -gpu_lw_dn_surf_g) / (1.0e-5_jprb + 2.0e-12_jprb &
           & *max(abs(flux%lw_dn_surf_g(:,istartcol:iendcol)),abs(gpu_lw_dn_surf_g)))))
      gpu_error_ratio = max(gpu_error_ratio,maxval(abs(flux%cloud_cover_lw(istartcol:iendcol)-total_cloud_cover_batch) &
           & / (5.0e-13_jprb + 2.0e-12_jprb*max(abs(flux%cloud_cover_lw(istartcol:iendcol)), &
           & abs(total_cloud_cover_batch)))))
      if (config%do_lw_derivatives) then
        gpu_error_ratio = max(gpu_error_ratio,maxval(abs(flux%lw_derivatives(istartcol:iendcol,:)-gpu_lw_derivatives) &
             & / (1.0e-9_jprb + 2.0e-12_jprb*max(abs(flux%lw_derivatives(istartcol:iendcol,:)), &
             & abs(gpu_lw_derivatives)))))
      end if
      if (gpu_error_ratio > 1.0_jprb) then
        write(nulerr,'(a,es12.4)') '*** CUDA longwave validation failed; normalized max error = ',gpu_error_ratio
        write(nulerr,'(a,es12.4)') '    lw_up_clear max abs = ', &
             & maxval(abs(flux%lw_up_clear(istartcol:iendcol,:)-gpu_lw_up_clear))
        write(nulerr,'(a,es12.4)') '    lw_dn_clear max abs = ', &
             & maxval(abs(flux%lw_dn_clear(istartcol:iendcol,:)-gpu_lw_dn_clear))
        write(nulerr,'(a,es12.4)') '    lw_up max abs = ', &
             & maxval(abs(flux%lw_up(istartcol:iendcol,:)-gpu_lw_up))
        write(nulerr,'(a,es12.4)') '    lw_dn max abs = ', &
             & maxval(abs(flux%lw_dn(istartcol:iendcol,:)-gpu_lw_dn))
        write(nulerr,'(a,es12.4)') '    lw_dn_surf_clear_g max abs = ', &
             & maxval(abs(flux%lw_dn_surf_clear_g(:,istartcol:iendcol)-gpu_lw_dn_surf_clear_g))
        write(nulerr,'(a,es12.4)') '    lw_dn_surf_g max abs = ', &
             & maxval(abs(flux%lw_dn_surf_g(:,istartcol:iendcol)-gpu_lw_dn_surf_g))
        if (config%do_lw_derivatives) write(nulerr,'(a,es12.4)') &
             & '    lw_derivatives max abs = ', &
             & maxval(abs(flux%lw_derivatives(istartcol:iendcol,:)-gpu_lw_derivatives))
        call radiation_abort()
      end if
      flux%lw_up_clear(istartcol:iendcol,:) = gpu_lw_up_clear
      flux%lw_dn_clear(istartcol:iendcol,:) = gpu_lw_dn_clear
      flux%lw_up(istartcol:iendcol,:) = gpu_lw_up
      flux%lw_dn(istartcol:iendcol,:) = gpu_lw_dn
      flux%lw_dn_surf_clear_g(:,istartcol:iendcol) = gpu_lw_dn_surf_clear_g
      flux%lw_dn_surf_g(:,istartcol:iendcol) = gpu_lw_dn_surf_g
      flux%cloud_cover_lw(istartcol:iendcol) = total_cloud_cover_batch
      if (config%do_lw_derivatives) &
           & flux%lw_derivatives(istartcol:iendcol,:) = gpu_lw_derivatives
    end if
#endif

    if (lhook) call dr_hook('radiation_mcica_lw:solver_mcica_lw',1,hook_handle)
    
  end subroutine solver_mcica_lw

end module radiation_mcica_lw
