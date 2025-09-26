! radiation_disort_sw.F90 - Shortwave DISORT solver
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

module radiation_disort_sw

  public

contains

  !---------------------------------------------------------------------
  ! Longwave disort solver, in which clouds are assumed to fill
  ! the gridbox horizontally
  subroutine solver_disort_sw(nlev,istartcol,iendcol, &
       &  config, single_level, cloud, & 
       &  od, ssa, g, od_cloud, ssa_cloud, pf_cloud, &
       &  albedo_direct, albedo_diffuse, incoming_sw, &
       &  flux)

    use parkind1, only           : jprb
    use yomhook,  only           : lhook, dr_hook, jphook

    use radiation_constants, only      : Pi
    use radiation_config, only         : config_type
    use radiation_single_level, only   : single_level_type
    use radiation_cloud, only          : cloud_type
    use radiation_flux, only           : flux_type, indexed_sum_profile, &
         &                               add_indexed_sum_profile
    use radiation_gaussian_quadrature, only : gaussian_quadrature, self_test

    implicit none

    ! Number of angles for radiance calculation (not needed, so set to
    ! 1)
    integer, parameter :: nphi = 1

    ! Inputs
    integer, intent(in) :: nlev               ! number of model levels
    integer, intent(in) :: istartcol, iendcol ! range of columns to process
    type(config_type),        intent(in) :: config
    type(single_level_type),  intent(in) :: single_level
    type(cloud_type),         intent(in) :: cloud

    ! Gas and aerosol optical depth, single-scattering albedo and
    ! asymmetry factor at each longwave g-point
    real(jprb), intent(in), dimension(config%n_g_sw, nlev, istartcol:iendcol) :: &
         &  od, ssa, g

    ! Cloud and precipitation optical depth, single-scattering albedo and
    ! asymmetry factor in each longwave band
    real(jprb), intent(in), dimension(config%n_bands_sw,nlev,istartcol:iendcol)   :: &
         &  od_cloud, ssa_cloud
    real(jprb), intent(in), dimension(config%n_bands_sw, &
         &  nlev,istartcol:iendcol,config%n_pf_sw) :: pf_cloud

    ! Direct and diffuse surface albedos, and the incoming shortwave
    ! flux into a plane perpendicular to the incoming radiation at
    ! top-of-atmosphere in each of the shortwave g points
    real(jprb), intent(in), dimension(config%n_g_sw, istartcol:iendcol) &
         &  :: albedo_direct, albedo_diffuse, incoming_sw

    ! Output
    type(flux_type), intent(inout):: flux

    ! Local variables

    ! Optical depth and single scattering albedo
    real(jprb), dimension(nlev) :: od_in, ssa_in

    ! Coefficients of a Legendre decompositon of the phase function
    real(jprb), dimension(0:config%n_pf_sw,nlev) :: pf_mom

    ! Quadrature points for angular integration
    real(jprb) :: mu_quad(config%n_angles_per_hemisphere_sw*2)
    real(jprb) :: weight_quad(config%n_angles_per_hemisphere_sw*2)

    ! Dummy arguments to DISORT
    real(jprb) :: phi(nphi), mu(config%n_angles_per_hemisphere_sw*2)
    real(jprb) :: h_layer(0:nlev)
    real(jprb) :: rhoq(config%n_angles_per_hemisphere_sw, &
         &             0:config%n_angles_per_hemisphere_sw, &
         &             0:(config%n_angles_per_hemisphere_sw*2-1))
    real(jprb) :: rhou(config%n_angles_per_hemisphere_sw*2, &
         &             0:config%n_angles_per_hemisphere_sw, &
         &             0:(config%n_angles_per_hemisphere_sw*2-1))
    real(jprb) :: rho_accurate(config%n_angles_per_hemisphere_sw*2,nphi)
    real(jprb) :: bemst(config%n_angles_per_hemisphere_sw)
    real(jprb) :: emust(config%n_angles_per_hemisphere_sw*2)
    real(jprb) :: dfdt(nlev+1), uavg(nlev+1)
    real(jprb) :: uu(config%n_angles_per_hemisphere_sw*2,nlev+1,nphi)
    real(jprb) :: albmed(config%n_angles_per_hemisphere_sw*2)
    real(jprb) :: trnmed(config%n_angles_per_hemisphere_sw*2)
    real(jprb) :: od_out(nlev+1)

    ! Planck functions per sterad (black-body flux divided by pi)
    real(jprb) :: planck_hl_per_sr(nlev+1), planck_surf_per_sr

    character*127 :: header

    ! Fluxes per g point, clear- and total-sky
    real(jprb), dimension(config%n_g_sw, nlev+1) :: flux_up_clear, flux_dn_dif_clear, flux_dn_dir_clear
    real(jprb), dimension(config%n_g_sw, nlev+1) :: flux_up, flux_dn_dif, flux_dn_dir

    ! Temporary scattering optical depth
    real(jprb) :: scat_od_new(nlev)

    ! Is there any cloud in the profile?
    logical :: is_cloudy_profile

    ! Number of g points
    integer :: ng

    ! Loop indices for level and column
    integer :: jlev, jcol, jg, jcomp

    real(jphook) :: hook_handle

    if (lhook) call dr_hook('radiation_disort_sw:solver_disort_sw',0,hook_handle)

    h_layer = 0.0_jprb
    
    !    call self_test(3)
    call gaussian_quadrature(config%i_quadrature_sw, config%n_angles_per_hemisphere_sw, &
         &  mu_quad, weight_quad, config%i_quadrature_power_sw)

    ! Duplicate for opposite hemisphere
    do jcomp = 1,config%n_angles_per_hemisphere_sw
      mu_quad    (jcomp + config%n_angles_per_hemisphere_sw) = -mu_quad    (jcomp)
      weight_quad(jcomp + config%n_angles_per_hemisphere_sw) =  weight_quad(jcomp)
    end do

    ng = config%n_g_sw

    header = 'ecRad profile'
    phi = 0.0_jprb
    mu = 0.0_jprb

    ! Loop through columns
    do jcol = istartcol,iendcol

      ! Only perform calculation if sun above the horizon
      if (single_level%cos_sza(jcol) > 0.0_jprb) then

        ! Is there any cloud in the profile?
        is_cloudy_profile = .false.
        do jlev = 1,nlev
          if (cloud%fraction(jcol,jlev) >= config%cloud_fraction_threshold) then
            is_cloudy_profile = .true.
            exit
          end if
        end do

        ! Clear sky calculation
        do jg = 1,ng
          od_in = od(jg,:,jcol)
          pf_mom = 0.0_jprb;
          pf_mom(0,:) = 1.0_jprb

          ssa_in = ssa(jg,:,jcol)
          ! Assume aerosols have a Henyey-Greenstein phase function
          ! expressed here in terms of Legendre polynomials
          do jcomp = 1,config%n_pf_sw
            ! Unnormalized coefficients
            !pf_mom(jcomp,:) = g(jg,:,jcol)**jcomp * (2.0_jprb*jcomp + 1.0_jprb)
            ! DISORT normalization
            pf_mom(jcomp,:) = g(jg,:,jcol)**jcomp
          end do

          planck_hl_per_sr = 0.0_jprb
          planck_surf_per_sr = 0.0_jprb

          call disort(nlev, config%n_pf_sw, &
               &  config%n_angles_per_hemisphere_sw*2, &
               &  config%n_angles_per_hemisphere_sw*2, nphi, nlev+1, &
               &  .false., .false., 0, .true., spread(.false.,1,5), &
               &  .false., .true., .false., .false., od_in, &
               &  ssa_in, pf_mom, planck_hl_per_sr, -1.0_jprb, -1.0_jprb, &
               &  od_out, single_level%cos_sza(jcol), 0.0_jprb, mu, phi, incoming_sw(jg,jcol), &
               &  0.0_jprb, albedo_diffuse(jg,jcol), planck_surf_per_sr, &
               &  0.0_jprb, 1.0_jprb, 6371.0_jprb, h_layer, &
               &  rhoq, rhou, rho_accurate, bemst, emust, &
               &  0.0_jprb, header, flux_dn_dir_clear(jg,:), flux_dn_dif_clear(jg,:), flux_up_clear(jg,:), &
               &  dfdt, uavg, uu, albmed, trnmed, mu_quad, weight_quad)

          if (is_cloudy_profile) then
            ! Add cloud properties to existing profile
            scat_od_new = od_in*ssa_in &
                 &  + od_cloud(config%i_band_from_reordered_g_sw(jg),:,jcol) &
                 &  * ssa_cloud(config%i_band_from_reordered_g_sw(jg),:,jcol)

            ! Aerosol scattering properties are already present - need
            ! to add cloud properties with a weighted average by
            ! scattering coefficient
            do jlev = 1,nlev
              pf_mom(1:,jlev) = (od_in(jlev)*ssa_in(jlev)*pf_mom(1:,jlev) &
                   &  + od_cloud(config%i_band_from_reordered_g_sw(jg),jlev,jcol) &
                   &  * ssa_cloud(config%i_band_from_reordered_g_sw(jg),jlev,jcol) &
                   &  * pf_cloud(config%i_band_from_reordered_g_sw(jg),jlev,jcol,:)) &
                   &  / max(scat_od_new(jlev), 1.0e-24_jprb)
            end do

            od_in = od_in + od_cloud(config%i_band_from_reordered_g_sw(jg),:,jcol)
            ssa_in = scat_od_new / max(od_in, 1.0e-24_jprb)

            call disort(nlev, config%n_pf_sw, &
                 &  config%n_angles_per_hemisphere_sw*2, &
                 &  config%n_angles_per_hemisphere_sw*2, nphi, nlev+1, &
                 &  .false., .false., 0, .true., spread(.false.,1,5), &
                 &  .false., .true., .false., .false., od_in, &
                 &  ssa_in, pf_mom, planck_hl_per_sr, -1.0_jprb, -1.0_jprb, &
                 &  od_out, single_level%cos_sza(jcol), 0.0_jprb, mu, phi, incoming_sw(jg,jcol), &
                 &  0.0_jprb, albedo_diffuse(jg,jcol), planck_surf_per_sr, &
                 &  0.0_jprb, 1.0_jprb, 6371.0_jprb, h_layer, &
                 &  rhoq, rhou, rho_accurate, bemst, emust, &
                 &  0.0_jprb, header, flux_dn_dir(jg,:), flux_dn_dif(jg,:), flux_up(jg,:), &
                 &  dfdt, uavg, uu, albmed, trnmed, mu_quad, weight_quad)
          end if

        end do

        ! Sum over g-points to compute broadband fluxes
        flux%sw_up_clear(jcol,:)     = sum(flux_up_clear,1)
        flux%sw_dn_direct_clear(jcol,:) = sum(flux_dn_dir_clear,1)
        flux%sw_dn_clear(jcol,:)     = sum(flux_dn_dif_clear,1) + flux%sw_dn_direct_clear(jcol,:)
        ! Store surface spectral downwelling fluxes
        flux%sw_dn_diffuse_surf_clear_g(:,jcol) = flux_dn_dif_clear(:,nlev+1)
        flux%sw_dn_direct_surf_clear_g(:,jcol)  = flux_dn_dir_clear(:,nlev+1)

        ! Save the spectral fluxes if required
        if (config%do_save_spectral_flux) then
          call indexed_sum_profile(flux_up_clear, config%i_spec_from_reordered_g_sw, &
               &                   flux%sw_up_clear_band(:,jcol,:))
          call indexed_sum_profile(flux_dn_dir_clear, config%i_spec_from_reordered_g_sw, &
               &                   flux%sw_dn_clear_band(:,jcol,:))
          if (allocated(flux%sw_dn_direct_clear_band)) then
            flux%sw_dn_direct_clear_band(:,jcol,:) &
                 &  = flux%sw_dn_clear_band(:,jcol,:)
          end if
          call add_indexed_sum_profile(flux_dn_dif_clear, &
               &                       config%i_spec_from_reordered_g_sw, &
               &                       flux%sw_dn_clear_band(:,jcol,:))
        end if

        if (is_cloudy_profile) then
          ! Sum over g-points to compute broadband fluxes
          flux%sw_up(jcol,:) = sum(flux_up,1)
          flux%sw_dn_direct(jcol,:) = sum(flux_dn_dir,1)
          flux%sw_dn(jcol,:) = sum(flux_dn_dif,1)+flux%sw_dn_direct(jcol,:)
          ! Store surface spectral downwelling fluxes
          flux%sw_dn_diffuse_surf_g(:,jcol) = flux_dn_dif(:,nlev+1)
          flux%sw_dn_direct_surf_g(:,jcol)  = flux_dn_dir(:,nlev+1)

          ! Save the spectral fluxes if required
          if (config%do_save_spectral_flux) then
            call indexed_sum_profile(flux_up, config%i_spec_from_reordered_g_sw, &
                 &                   flux%sw_up_band(:,jcol,:))
            call indexed_sum_profile(flux_dn_dir, config%i_spec_from_reordered_g_sw, &
                 &                   flux%sw_dn_band(:,jcol,:))
            if (allocated(flux%sw_dn_direct_band)) then
              flux%sw_dn_direct_band(:,jcol,:) &
                   &  = flux%sw_dn_band(:,jcol,:)
            end if
            call add_indexed_sum_profile(flux_dn_dif, &
                 &                       config%i_spec_from_reordered_g_sw, &
                 &                       flux%sw_dn_band(:,jcol,:))          
          end if
        else
          ! Clear-sky profile: copy over clear-sky fluxes
          flux%sw_up(jcol,:) = flux%sw_up_clear(jcol,:)
          flux%sw_dn(jcol,:) = flux%sw_dn_clear(jcol,:)
          flux%sw_dn_direct(jcol,:) = flux%sw_dn_direct_clear(jcol,:)
          flux%sw_dn_diffuse_surf_clear_g(:,jcol) = flux_dn_dif_clear(:,nlev+1)
          flux%sw_dn_direct_surf_clear_g(:,jcol) = flux_dn_dir_clear(:,nlev+1)
          if (config%do_save_spectral_flux) then
            flux%sw_up_band(:,jcol,:) = flux%sw_up_clear_band(:,jcol,:)
            flux%sw_dn_band(:,jcol,:) = flux%sw_dn_clear_band(:,jcol,:)
            flux%sw_dn_direct_band(:,jcol,:) = flux%sw_dn_direct_clear_band(:,jcol,:)
          end if
        end if

      else ! Sun below the horizon
        ! Set fluxes to zero if sun is below the horizon
        flux%sw_up(jcol,:) = 0.0_jprb
        flux%sw_dn(jcol,:) = 0.0_jprb
        if (allocated(flux%sw_dn_direct)) then
          flux%sw_dn_direct(jcol,:) = 0.0_jprb
        end if
        flux%sw_dn_diffuse_surf_g(:,jcol) = 0.0_jprb
        flux%sw_dn_direct_surf_g(:,jcol)  = 0.0_jprb

        if (config%do_clear) then
          flux%sw_up_clear(jcol,:) = 0.0_jprb
          flux%sw_dn_clear(jcol,:) = 0.0_jprb
          if (allocated(flux%sw_dn_direct_clear)) then
            flux%sw_dn_direct_clear(jcol,:) = 0.0_jprb
          end if
          flux%sw_dn_diffuse_surf_clear_g(:,jcol) = 0.0_jprb
          flux%sw_dn_direct_surf_clear_g(:,jcol)  = 0.0_jprb
        end if

        if (config%do_save_spectral_flux) then
          flux%sw_dn_band(:,jcol,:) = 0.0_jprb
          flux%sw_up_band(:,jcol,:) = 0.0_jprb
          if (allocated(flux%sw_dn_direct_band)) then
            flux%sw_dn_direct_band(:,jcol,:) = 0.0_jprb
          end if
          if (config%do_clear) then
            flux%sw_dn_clear_band(:,jcol,:) = 0.0_jprb
            flux%sw_up_clear_band(:,jcol,:) = 0.0_jprb
            if (allocated(flux%sw_dn_direct_clear_band)) then
              flux%sw_dn_direct_clear_band(:,jcol,:) = 0.0_jprb
            end if
          end if
        end if
      end if

    end do

    if (lhook) call dr_hook('radiation_disort_sw:solver_disort_sw',1,hook_handle)

  end subroutine solver_disort_sw

end module radiation_disort_sw
