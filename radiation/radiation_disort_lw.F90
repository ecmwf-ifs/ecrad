! radiation_disort_lw.F90 - Longwave DISORT solver
!
! (C) Copyright 2022- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
!
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!
! Author:  Robin Hogan
! Email:   r.j.hogan@ecmwf.int
!
! Modifications

module radiation_disort_lw

  public

contains

  !---------------------------------------------------------------------
  ! Longwave disort solver, in which clouds are assumed to fill
  ! the gridbox horizontally
  subroutine solver_disort_lw(nlev,istartcol,iendcol, &
       &  config, cloud, & 
       &  od, ssa, g, od_cloud, ssa_cloud, pf_cloud, planck_hl, &
       &  emission, albedo, &
       &  flux, do_homogeneous)

    use parkind1, only           : jprb
    use yomhook,  only           : lhook, dr_hook

    use radiation_constants, only      : Pi
    use radiation_config, only         : config_type
    use radiation_cloud, only          : cloud_type
    use radiation_flux, only           : flux_type, indexed_sum_profile
    !use radiation_lw_derivatives, only : calc_lw_derivatives_ica
    use radiation_gaussian_quadrature, only : gaussian_quadrature, self_test, QuadratureName
    use radiation_cloud_cover, only    : cloud_cover
    use radiation_io, only             : nulout
    
    implicit none

    ! Number of angles for radiance calculation (not needed, so set to
    ! 1)
    integer, parameter :: nphi = 1
    
    ! Inputs
    integer, intent(in) :: nlev               ! number of model levels
    integer, intent(in) :: istartcol, iendcol ! range of columns to process
    type(config_type),        intent(in) :: config
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
         &  nlev,istartcol:iendcol) :: ssa_cloud
    real(jprb), intent(in), dimension(config%n_bands_lw_if_scattering, &
         &  nlev,istartcol:iendcol,config%n_pf_lw) :: pf_cloud

    ! Planck function at each half-level and the surface
    real(jprb), intent(in), dimension(config%n_g_lw,nlev+1,istartcol:iendcol) :: &
         &  planck_hl
  
    ! Emission (Planck*emissivity) and albedo (1-emissivity) at the
    ! surface at each longwave g-point
    real(jprb), intent(in), dimension(config%n_g_lw, istartcol:iendcol) &
         &  :: emission, albedo

    ! Do we homogenize the cloud to fill the gridbox (default false)?
    logical, optional, intent(in) :: do_homogeneous
    
    ! Output
    type(flux_type), intent(inout):: flux

    ! Local variables

    ! Optical depth and single scattering albedo
    real(jprb), dimension(nlev) :: od_in, ssa_in

    ! Coefficients of a Legendre decompositon of the phase function
    real(jprb), dimension(0:config%n_pf_lw,nlev) :: pf_mom

    ! Quadrature points for angular integration
    real(jprb) :: mu_quad(config%n_angles_per_hemisphere_lw*2)
    real(jprb) :: weight_quad(config%n_angles_per_hemisphere_lw*2)
    
    ! Dummy arguments to DISORT
    real(jprb) :: phi(nphi), mu(config%n_angles_per_hemisphere_lw*2)
    real(jprb) :: h_layer(0:nlev)
    real(jprb) :: rhoq(config%n_angles_per_hemisphere_lw, &
         &             0:config%n_angles_per_hemisphere_lw, &
         &             0:(config%n_angles_per_hemisphere_lw*2-1))
    real(jprb) :: rhou(config%n_angles_per_hemisphere_lw*2, &
         &             0:config%n_angles_per_hemisphere_lw, &
         &             0:(config%n_angles_per_hemisphere_lw*2-1))
    real(jprb) :: rho_accurate(config%n_angles_per_hemisphere_lw*2,nphi)
    real(jprb) :: bemst(config%n_angles_per_hemisphere_lw)
    real(jprb) :: emust(config%n_angles_per_hemisphere_lw*2)
    real(jprb) :: dfdt(nlev+1), uavg(nlev+1)
    real(jprb) :: uu(config%n_angles_per_hemisphere_lw*2,nlev+1,nphi)
    real(jprb) :: albmed(config%n_angles_per_hemisphere_lw*2)
    real(jprb) :: trnmed(config%n_angles_per_hemisphere_lw*2)
    real(jprb) :: flux_dn_dir(nlev+1), od_out(nlev+1)

    ! Planck functions per sterad (black-body flux divided by pi)
    real(jprb) :: planck_hl_per_sr(nlev+1), planck_surf_per_sr
    
    character*127 :: header
    
    ! Fluxes per g point, clear- and total-sky
    real(jprb), dimension(config%n_g_lw, nlev+1) :: flux_up_clear, flux_dn_clear
    real(jprb), dimension(config%n_g_lw, nlev+1) :: flux_up, flux_dn

    ! Temporary scattering optical depth
    real(jprb) :: scat_od_new(nlev)
    
    ! Is there any cloud in the profile?
    logical :: is_cloudy_profile

    ! Do we fill the cloud horizontally?
    logical :: do_homogeneous_local

    ! If not do_homogeneous_local, we scale the cloud optical depth by
    ! 1.0/cloud_cover to estimate the in-cloud value; note that the
    ! assumption is that all clouds have a cloud fraction equal to the
    ! cloud cover
    real(jprb) :: od_cloud_scaling

    ! Number of g points
    integer :: ng

    ! Loop indices for level and column
    integer :: jlev, jcol, jg, jcomp

    real(jprb) :: hook_handle

    if (lhook) call dr_hook('radiation_disort_lw:solver_disort_lw',0,hook_handle)

    if (present(do_homogeneous)) then
      do_homogeneous_local = do_homogeneous
    else
      do_homogeneous_local = .false.
    end if

!    call self_test(3)
    call gaussian_quadrature(config%i_quadrature_lw, config%n_angles_per_hemisphere_lw, &
         &  mu_quad, weight_quad, config%i_quadrature_power_lw)
    
    do jcomp = 1,config%n_angles_per_hemisphere_lw
      mu_quad    (jcomp + config%n_angles_per_hemisphere_lw) = -mu_quad    (jcomp)
      weight_quad(jcomp + config%n_angles_per_hemisphere_lw) =  weight_quad(jcomp)
    end do
    
    !if (config%iverbose >= 5) then
    !  write(nulout, *) 'Quadrature: ', QuadratureName(config%i_quadrature_lw), config%n_angles_per_hemisphere_lw, &
    !       &  ' mu=[', mu_quad, '] wt=[', weight_quad, ']'
    !end if

    ng = config%n_g_lw

    header = 'ecRad profile'
    !dummy = 0.5_jprb
    !dummy2 = 0.5_jprb
    !dummy3 = 0.5_jprb
    phi = 0.0_jprb
    mu = 0.0_jprb
    
    ! Loop through columns
    do jcol = istartcol,iendcol

      ! Is there any cloud in the profile?
      is_cloudy_profile = .false.
      do jlev = 1,nlev
        if (cloud%fraction(jcol,jlev) >= config%cloud_fraction_threshold) then
          is_cloudy_profile = .true.
          exit
        end if
      end do

      od_cloud_scaling = 1.0_jprb
      if (is_cloudy_profile) then
        if (.not. do_homogeneous_local) then
          ! Squash grid-box mean cloud properties in to the
          ! cloud-covered fraction of the gridbox
          flux%cloud_cover_lw(jcol) = cloud_cover(nlev, config%i_overlap_scheme, &
               &  cloud%fraction(jcol,:), cloud%overlap_param(jcol,:), config%use_beta_overlap)
          od_cloud_scaling = 1.0_jprb
        else
          ! Assumption is that clouds fill the box horizontally, so
          ! cloud cover is 1
          flux%cloud_cover_lw(jcol) = 1.0_jprb
        end if
      else
        flux%cloud_cover_lw(jcol) = 0.0_jprb
      end if

      ! Clear sky calculation
      do jg = 1,ng
        od_in = od(jg,:,jcol)
        pf_mom = 0.0_jprb;
        pf_mom(0,:) = 1.0_jprb
        if (config%do_lw_aerosol_scattering) then
          ssa_in = ssa(jg,:,jcol)
          ! Assume aerosols have a Henyey-Greenstein phase function
          ! expressed here in terms of Legendre polynomials
          do jcomp = 1,config%n_pf_lw
            ! Unnormalized coefficients
            !pf_mom(jcomp,:) = g(jg,:,jcol)**jcomp * (2.0_jprb*jcomp + 1.0_jprb)
            ! DISORT normalization
            pf_mom(jcomp,:) = g(jg,:,jcol)**jcomp
          end do
        else
          ssa_in = 0.0_jprb
        end if

        planck_hl_per_sr = planck_hl(jg,:,jcol) / Pi
        planck_surf_per_sr = emission(jg,jcol)/(Pi*(1.0_jprb-albedo(jg,jcol)))
        
        call disort(nlev, config%n_pf_lw, &
             &  config%n_angles_per_hemisphere_lw*2, &
             &  config%n_angles_per_hemisphere_lw*2, nphi, nlev+1, &
             &  .false., .false., 0, .true., spread(.false.,1,5), &
             &  .true., .true., .false., .false., od_in, &
             &  ssa_in, pf_mom, planck_hl_per_sr, -1.0_jprb, -1.0_jprb, &
             &  od_out, 1.0_jprb, 0.0_jprb, mu, phi, 0.0_jprb, &
             &  0.0_jprb, albedo(jg,jcol), planck_surf_per_sr, &
             &  0.0_jprb, 1.0_jprb, 6371.0_jprb, h_layer, &
             &  rhoq, rhou, rho_accurate, bemst, emust, &
             &  0.0_jprb, header, flux_dn_dir, flux_dn_clear(jg,:), flux_up_clear(jg,:), &
             &  dfdt, uavg, uu, albmed, trnmed, mu_quad, weight_quad)

        if (is_cloudy_profile) then
          ! Add cloud properties to existing profile
          if (config%do_lw_cloud_scattering) then
            scat_od_new = od_in*ssa_in + od_cloud_scaling &
                 &  * od_cloud(config%i_band_from_reordered_g_lw(jg),:,jcol) &
                 &  * ssa_cloud(config%i_band_from_reordered_g_lw(jg),:,jcol)

            if (config%do_lw_aerosol_scattering) then
              ! Aerosol scattering properties are already present -
              ! need to add cloud properties with a weighted average
              ! by scattering coefficient
              do jlev = 1,nlev
                pf_mom(1:,jlev) = (od_in(jlev)*ssa_in(jlev)*pf_mom(1:,jlev) &
                     &  + od_cloud_scaling &
                     &  * od_cloud(config%i_band_from_reordered_g_lw(jg),jlev,jcol) &
                     &  * ssa_cloud(config%i_band_from_reordered_g_lw(jg),jlev,jcol) &
                     &  * pf_cloud(config%i_band_from_reordered_g_lw(jg),jlev,jcol,:)) &
                     &  / max(scat_od_new(jlev), 1.0e-24_jprb)
              end do
            else
              ! Aerosol scattering properties are not present but
              ! cloud scattering properties are: copy over cloud
              ! properties
              pf_mom(1:,:) = transpose(pf_cloud(config%i_band_from_reordered_g_lw(jg),:,jcol,:))
            end if
                
            od_in = od_in + od_cloud_scaling &
                 &        * od_cloud(config%i_band_from_reordered_g_lw(jg),:,jcol)
            ssa_in = scat_od_new / max(od_in, 1.0e-24_jprb)
          else
            ! No cloud or aerosol scattering: simply add the cloud
            ! optical depth
            od_in = od_in + od_cloud_scaling &
                 &        * od_cloud(config%i_band_from_reordered_g_lw(jg),:,jcol)
          end if

          call disort(nlev, config%n_pf_lw, &
               &  config%n_angles_per_hemisphere_lw*2, &
               &  config%n_angles_per_hemisphere_lw*2, nphi, nlev+1, &
               &  .false., .false., 0, .true., spread(.false.,1,5), &
               &  .true., .true., .false., .false., od_in, &
               &  ssa_in, pf_mom, planck_hl_per_sr, -1.0_jprb, -1.0_jprb, &
               &  od_out, 1.0_jprb, 0.0_jprb, mu, phi, 0.0_jprb, &
               &  0.0_jprb, albedo(jg,jcol), planck_surf_per_sr, &
               &  0.0_jprb, 1.0_jprb, 6371.0_jprb, h_layer, &
               &  rhoq, rhou, rho_accurate, bemst, emust, &
               &  0.0_jprb, header, flux_dn_dir, flux_dn(jg,:), flux_up(jg,:), &
               &  dfdt, uavg, uu, albmed, trnmed, mu_quad, weight_quad)
        end if
        
      end do

      ! Sum over g-points to compute broadband fluxes
      flux%lw_up_clear(jcol,:) = sum(flux_up_clear,1)
      flux%lw_dn_clear(jcol,:) = sum(flux_dn_clear,1)
      ! Store surface spectral downwelling fluxes
      flux%lw_dn_surf_clear_g(:,jcol) = flux_dn_clear(:,nlev+1)
      
      ! Save the spectral fluxes if required
      if (config%do_save_spectral_flux) then
        call indexed_sum_profile(flux_up_clear, config%i_spec_from_reordered_g_lw, &
             &                   flux%lw_up_clear_band(:,jcol,:))
        call indexed_sum_profile(flux_dn_clear, config%i_spec_from_reordered_g_lw, &
             &                   flux%lw_dn_clear_band(:,jcol,:))
      end if

      if (is_cloudy_profile) then
        ! Sum over g-points to compute broadband fluxes
        flux%lw_up(jcol,:) = sum(flux_up,1)
        flux%lw_dn(jcol,:) = sum(flux_dn,1)
        ! Store surface spectral downwelling fluxes
        flux%lw_dn_surf_g(:,jcol) = flux_dn(:,nlev+1)
        
        ! Save the spectral fluxes if required
        if (config%do_save_spectral_flux) then
          call indexed_sum_profile(flux_up, config%i_spec_from_reordered_g_lw, &
               &                   flux%lw_up_band(:,jcol,:))
          call indexed_sum_profile(flux_dn, config%i_spec_from_reordered_g_lw, &
               &                   flux%lw_dn_band(:,jcol,:))
        end if

        if (.not. do_homogeneous_local) then
          ! Scale the all-sky fluxes to account for the cloud cover
          flux%lw_up(jcol,:) = flux%cloud_cover_lw(jcol) * flux%lw_up(jcol,:) &
               & + (1.0_jprb - flux%cloud_cover_lw(jcol))* flux%lw_up_clear(jcol,:)
          flux%lw_dn(jcol,:) = flux%cloud_cover_lw(jcol) * flux%lw_dn(jcol,:) &
               & + (1.0_jprb - flux%cloud_cover_lw(jcol))* flux%lw_dn_clear(jcol,:)
          flux%lw_dn_surf_g(:,jcol) = flux%cloud_cover_lw(jcol)*flux%lw_dn_surf_g(:,jcol) &
               & + (1.0_jprb - flux%cloud_cover_lw(jcol))* flux%lw_dn_surf_clear_g(:,jcol)
          if (config%do_save_spectral_flux) then
            flux%lw_up_band(:,jcol,:) = flux%cloud_cover_lw(jcol) * flux%lw_up_band(:,jcol,:) &
                 & + (1.0_jprb - flux%cloud_cover_lw(jcol)) * flux%lw_up_clear_band(:,jcol,:)
            flux%lw_dn_band(:,jcol,:) = flux%cloud_cover_lw(jcol) * flux%lw_dn_band(:,jcol,:) &
                 & + (1.0_jprb - flux%cloud_cover_lw(jcol)) * flux%lw_dn_clear_band(:,jcol,:)
          end if
        end if
        
      else
        ! Clear-sky profile: copy over clear-sky fluxes
        flux%lw_up(jcol,:) = flux%lw_up_clear(jcol,:)
        flux%lw_dn(jcol,:) = flux%lw_dn_clear(jcol,:)
        flux%lw_dn_surf_g(:,jcol) = flux_dn(:,nlev+1)
        if (config%do_save_spectral_flux) then
          flux%lw_up_band(:,jcol,:) = flux%lw_up_clear_band(:,jcol,:)
          flux%lw_dn_band(:,jcol,:) = flux%lw_dn_clear_band(:,jcol,:)
        end if
      end if
      
      if (config%do_lw_derivatives) then
        flux%lw_derivatives(jcol,:) = 1.0_jprb
      end if
      
    end do

    if (lhook) call dr_hook('radiation_disort_lw:solver_disort_lw',1,hook_handle)
    
  end subroutine solver_disort_lw


  !---------------------------------------------------------------------
  ! Longwave disort solver with no clouds
  subroutine solver_cloudless_disort_lw(nlev,istartcol,iendcol, &
       &  config, od, ssa, g, planck_hl, emission, albedo, &
       &  flux)

    use parkind1, only           : jprb
    use yomhook,  only           : lhook, dr_hook

    use radiation_constants, only      : Pi
    use radiation_config, only         : config_type
    use radiation_flux, only           : flux_type, indexed_sum_profile
    !use radiation_lw_derivatives, only : calc_lw_derivatives_ica
    use radiation_gaussian_quadrature, only : gaussian_quadrature, self_test
        
    implicit none

    ! Number of angles for radiance calculation (not needed, so set to
    ! 1)
    integer, parameter :: nphi = 1
    
    ! Inputs
    integer, intent(in) :: nlev               ! number of model levels
    integer, intent(in) :: istartcol, iendcol ! range of columns to process
    type(config_type),        intent(in) :: config

    ! Gas and aerosol optical depth, single-scattering albedo and
    ! asymmetry factor at each longwave g-point
    real(jprb), intent(in), dimension(config%n_g_lw, nlev, istartcol:iendcol) :: &
         &  od
    real(jprb), intent(in), dimension(config%n_g_lw_if_scattering, nlev, istartcol:iendcol) :: &
         &  ssa, g

    ! Planck function at each half-level and the surface
    real(jprb), intent(in), dimension(config%n_g_lw,nlev+1,istartcol:iendcol) :: &
         &  planck_hl
  
    ! Emission (Planck*emissivity) and albedo (1-emissivity) at the
    ! surface at each longwave g-point
    real(jprb), intent(in), dimension(config%n_g_lw, istartcol:iendcol) &
         &  :: emission, albedo

    ! Output
    type(flux_type), intent(inout):: flux

    ! Local variables

    ! Optical depth and single scattering albedo
    real(jprb), dimension(nlev) :: od_in, ssa_in

    ! Coefficients of a Legendre decompositon of the phase function
    real(jprb), dimension(0:config%n_pf_lw,nlev) :: pf_mom

    ! Quadrature points for angular integration
    real(jprb) :: mu_quad(config%n_angles_per_hemisphere_lw*2)
    real(jprb) :: weight_quad(config%n_angles_per_hemisphere_lw*2)
    
    ! Dummy arguments to DISORT
    real(jprb) :: phi(nphi), mu(config%n_angles_per_hemisphere_lw*2)
    real(jprb) :: h_layer(0:nlev)
    real(jprb) :: rhoq(config%n_angles_per_hemisphere_lw, &
         &             0:config%n_angles_per_hemisphere_lw, &
         &             0:(config%n_angles_per_hemisphere_lw*2-1))
    real(jprb) :: rhou(config%n_angles_per_hemisphere_lw*2, &
         &             0:config%n_angles_per_hemisphere_lw, &
         &             0:(config%n_angles_per_hemisphere_lw*2-1))
    real(jprb) :: rho_accurate(config%n_angles_per_hemisphere_lw*2,nphi)
    real(jprb) :: bemst(config%n_angles_per_hemisphere_lw)
    real(jprb) :: emust(config%n_angles_per_hemisphere_lw*2)
    real(jprb) :: dfdt(nlev+1), uavg(nlev+1)
    real(jprb) :: uu(config%n_angles_per_hemisphere_lw*2,nlev+1,nphi)
    real(jprb) :: albmed(config%n_angles_per_hemisphere_lw*2)
    real(jprb) :: trnmed(config%n_angles_per_hemisphere_lw*2)
    real(jprb) :: flux_dn_dir(nlev+1), od_out(nlev+1)

    ! Planck functions per sterad (black-body flux divided by pi)
    real(jprb) :: planck_hl_per_sr(nlev+1), planck_surf_per_sr
    
    character*127 :: header
    
    ! Fluxes per g point, clear- and total-sky
    real(jprb), dimension(config%n_g_lw, nlev+1) :: flux_up_clear, flux_dn_clear

    ! Temporary scattering optical depth
    real(jprb) :: scat_od_new(nlev)
    
    ! Number of g points
    integer :: ng

    ! Loop indices for level and column
    integer :: jlev, jcol, jg, jcomp

    real(jprb) :: hook_handle

    if (lhook) call dr_hook('radiation_disort_lw:solver_cloudless_disort_lw',0,hook_handle)

!    call self_test(3)
    call gaussian_quadrature(config%i_quadrature_lw, config%n_angles_per_hemisphere_lw, &
         &  mu_quad, weight_quad, config%i_quadrature_power_lw)

    do jcomp = 1,config%n_angles_per_hemisphere_lw
      mu_quad    (jcomp + config%n_angles_per_hemisphere_lw) = -mu_quad    (jcomp)
      weight_quad(jcomp + config%n_angles_per_hemisphere_lw) =  weight_quad(jcomp)
    end do
    
    ng = config%n_g_lw

    header = 'ecRad profile'
    !dummy = 0.5_jprb
    !dummy2 = 0.5_jprb
    !dummy3 = 0.5_jprb
    phi = 0.0_jprb
    mu = 0.0_jprb
    
    ! Loop through columns
    do jcol = istartcol,iendcol

      ! Clear sky calculation
      do jg = 1,ng
        od_in = od(jg,:,jcol)
        pf_mom = 0.0_jprb;
        pf_mom(0,:) = 1.0_jprb
        if (config%do_lw_aerosol_scattering) then
          ssa_in = ssa(jg,:,jcol)
          ! Assume aerosols have a Henyey-Greenstein phase function
          ! expressed here in terms of Legendre polynomials
          do jcomp = 1,config%n_pf_lw
            ! Unnormalized coefficients
            !pf_mom(jcomp,:) = g(jg,:,jcol)**jcomp * (2.0_jprb*jcomp + 1.0_jprb)
            ! DISORT normalization
            pf_mom(jcomp,:) = g(jg,:,jcol)**jcomp
          end do
        else
          ssa_in = 0.0_jprb
        end if

        planck_hl_per_sr = planck_hl(jg,:,jcol) / Pi
        planck_surf_per_sr = emission(jg,jcol)/(Pi*(1.0_jprb-albedo(jg,jcol)))
        
        call disort(nlev, config%n_pf_lw, &
             &  config%n_angles_per_hemisphere_lw*2, &
             &  config%n_angles_per_hemisphere_lw*2, nphi, nlev+1, &
             &  .false., .false., 0, .true., spread(.false.,1,5), &
             &  .true., .true., .false., .false., od_in, &
             &  ssa_in, pf_mom, planck_hl_per_sr, -1.0_jprb, -1.0_jprb, &
             &  od_out, 1.0_jprb, 0.0_jprb, mu, phi, 0.0_jprb, &
             &  0.0_jprb, albedo(jg,jcol), planck_surf_per_sr, &
             &  0.0_jprb, 1.0_jprb, 6371.0_jprb, h_layer, &
             &  rhoq, rhou, rho_accurate, bemst, emust, &
             &  0.0_jprb, header, flux_dn_dir, flux_dn_clear(jg,:), flux_up_clear(jg,:), &
             &  dfdt, uavg, uu, albmed, trnmed, mu_quad, weight_quad)

      end do

      ! Sum over g-points to compute broadband fluxes
      flux%lw_up(jcol,:) = sum(flux_up_clear,1)
      flux%lw_dn(jcol,:) = sum(flux_dn_clear,1)
      ! Store surface spectral downwelling fluxes
      flux%lw_dn_surf_g(:,jcol) = flux_dn_clear(:,nlev+1)
      
      ! Save the spectral fluxes if required
      if (config%do_save_spectral_flux) then
        call indexed_sum_profile(flux_up_clear, config%i_spec_from_reordered_g_lw, &
             &                   flux%lw_up_band(:,jcol,:))
        call indexed_sum_profile(flux_dn_clear, config%i_spec_from_reordered_g_lw, &
             &                   flux%lw_dn_band(:,jcol,:))
      end if

      if (config%do_clear) then
        ! Copy over clear-sky fluxes
        flux%lw_up_clear(jcol,:) = flux%lw_up(jcol,:)
        flux%lw_dn_clear(jcol,:) = flux%lw_dn(jcol,:)
        flux%lw_dn_surf_clear_g(:,jcol) = flux_dn_clear(:,nlev+1)
        if (config%do_save_spectral_flux) then
          flux%lw_up_clear_band(:,jcol,:) = flux%lw_up_band(:,jcol,:)
          flux%lw_dn_clear_band(:,jcol,:) = flux%lw_dn_band(:,jcol,:)
        end if
      end if
      
      if (config%do_lw_derivatives) then
        flux%lw_derivatives(jcol,:) = 1.0_jprb
      end if
      
    end do
    
    if (lhook) call dr_hook('radiation_disort_lw:solver_cloudless_disort_lw',1,hook_handle)
    
  end subroutine solver_cloudless_disort_lw

  
end module radiation_disort_lw
