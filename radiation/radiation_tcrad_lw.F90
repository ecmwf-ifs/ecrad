! radiation_tcrad_lw.F90 - Interface to TCRAD solver
!
! (C) Copyright 2021- ECMWF.
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

module radiation_tcrad_lw

  public

contains

  subroutine solver_tcrad_lw(nlev,istartcol,iendcol, &
       &  config, thermodynamics, cloud, & 
       &  od, ssa, g, od_cloud, ssa_cloud, g_cloud, planck_hl, &
       &  emission, albedo, &
       &  flux)

    use parkind1, only           : jprb
    use yomhook,  only           : lhook, dr_hook

!    use radiation_io, only             : nulout
    use radiation_config, only         : config_type, IPdfShapeGamma
    use radiation_thermodynamics, only : thermodynamics_type
    use radiation_cloud, only          : cloud_type
    use radiation_flux, only           : flux_type, indexed_sum_profile
    use tcrad_clear_sky_solver, only   : calc_clear_sky_flux
    use tcrad_2region_solver, only     : calc_flux_2region => calc_flux, &
         &                               calc_no_scattering_flux_2region => calc_no_scattering_flux
    use tcrad_3region_solver, only     : calc_flux_3region => calc_flux, &
         &                               calc_no_scattering_flux_3region => calc_no_scattering_flux
    use radiation_constants,  only : GasConstantDryAir, AccelDueToGravity

    implicit none

    real(jprb), parameter :: R_over_g = GasConstantDryAir / AccelDueToGravity

    ! Inputs
    integer, intent(in) :: nlev               ! number of model levels
    integer, intent(in) :: istartcol, iendcol ! range of columns to process
    type(config_type),        intent(in) :: config
    type(thermodynamics_type),intent(in) :: thermodynamics
    type(cloud_type),         intent(in) :: cloud

    ! Gas and aerosol optical depth of each layer at each longwave
    ! g-point
    real(jprb), intent(in), dimension(config%n_g_lw,nlev,istartcol:iendcol) :: od

    ! Gas and aerosol single-scattering albedo and asymmetry factor,
    ! only if longwave scattering by aerosols is to be represented
    real(jprb), intent(in), &
         &  dimension(config%n_g_lw_if_scattering,nlev,istartcol:iendcol) :: ssa, g

    ! Cloud and precipitation optical depth of each layer in each
    ! longwave band
    real(jprb), intent(in) :: od_cloud(config%n_bands_lw,nlev,istartcol:iendcol)

    ! Cloud and precipitation single-scattering albedo and asymmetry
    ! factor, only if longwave scattering by clouds is to be
    ! represented
    real(jprb), intent(in), dimension(config%n_bands_lw_if_scattering, &
         &                            nlev,istartcol:iendcol) :: ssa_cloud, g_cloud

    ! Planck function (emitted flux from a black body) at half levels
    ! and at the surface at each longwave g-point
    real(jprb), intent(in), dimension(config%n_g_lw,nlev+1,istartcol:iendcol) :: planck_hl

    ! Emission (Planck*emissivity) and albedo (1-emissivity) at the
    ! surface at each longwave g-point
    real(jprb), intent(in), dimension(config%n_g_lw, istartcol:iendcol) :: emission, albedo

    ! Output
    type(flux_type), intent(inout):: flux

    ! Local variables

    ! Regridded cloud properties
    real(jprb), dimension(config%n_g_lw,nlev) :: od_cloud_regrid, ssa_cloud_regrid, g_cloud_regrid

    ! Additional geometric properties needed for 3D effects: inverse
    ! of the cloud separation scale (Fielding et al. 20202) (m-1), and
    ! the layer thickness (m)
    real(jprb), dimension(nlev) :: inv_cloud_separation_scale, layer_thickness

    ! Spectral fluxes
    real(jprb), dimension(config%n_g_lw,nlev+1) :: flux_up, flux_dn

    real(jprb) :: hook_handle

    ! Column loop index
    integer :: jcol

    if (lhook) call dr_hook('radiation_tcrad_lw:solver_tcrad_lw',0,hook_handle)

    flux%lw_up(istartcol:iendcol,:) = 0.0_jprb
    flux%lw_dn(istartcol:iendcol,:) = 0.0_jprb

    if (.not. config%do_lw_cloud_scattering) then
      ! FIX we have allocated memory that is not used
      ssa_cloud_regrid = 0.0_jprb
      g_cloud_regrid   = 0.0_jprb
    end if

    do jcol = istartcol,iendcol

      if (config%do_clear) then
        call calc_clear_sky_flux(config%n_g_lw, nlev, emission(:,jcol), albedo(:,jcol), &
             &  planck_hl(:,:,jcol), od(:,:,jcol), flux_up, flux_dn, &
             &  n_angles_per_hem=config%n_angles_per_hemisphere_lw)

        flux%lw_up_clear(jcol,:) = sum(flux_up,1)
        flux%lw_dn_clear(jcol,:) = sum(flux_dn,1)

        if (config%do_save_spectral_flux) then
          call indexed_sum_profile(flux_up, &
               &           config%i_spec_from_reordered_g_lw, &
               &           flux%lw_up_clear_band(:,jcol,:))
          call indexed_sum_profile(flux_dn, &
               &           config%i_spec_from_reordered_g_lw, &
               &           flux%lw_dn_clear_band(:,jcol,:))
        end if

      end if

      ! If we do 3D effects then we need to provide the layer
      ! thickness and cloud separation scale
      if (config%do_3d_effects) then
        layer_thickness = R_over_g * 0.5_jprb*(thermodynamics%temperature_hl(jcol,1:nlev) &
             &                                +thermodynamics%temperature_hl(jcol,2:nlev+1)) &
             &               * log(thermodynamics%pressure_hl(jcol,2:nlev+1) &
             &                     /max(thermodynamics%pressure_hl(jcol,1:nlev),1.0e-2_jprb))
        inv_cloud_separation_scale = cloud%inv_cloud_effective_size(jcol,:) &
             &  * sqrt(cloud%fraction(jcol,:)*(1.0_jprb-cloud%fraction(jcol,:)))
      end if

      od_cloud_regrid    = od_cloud(config%i_band_from_reordered_g_lw,:,jcol)
      if (config%do_lw_cloud_scattering) then
        ssa_cloud_regrid = ssa_cloud(config%i_band_from_reordered_g_lw,:,jcol)
        g_cloud_regrid   = g_cloud(config%i_band_from_reordered_g_lw,:,jcol)
        if (config%nregions == 2) then
          if (config%do_3d_effects) then
            ! Two regions with scattering and 3D effects
            call calc_flux_2region(config%n_g_lw, nlev, emission(:,jcol), albedo(:,jcol), &
                 &         planck_hl(:,:,jcol), cloud%fraction(jcol,:), &
                 &         od(:,:,jcol), od_cloud_regrid, ssa_cloud_regrid, g_cloud_regrid, &
                 &         cloud%overlap_param(jcol,:), &
                 &         flux_up, flux_dn, &
                 &         n_angles_per_hem=config%n_angles_per_hemisphere_lw, &
                 &         layer_thickness=layer_thickness, inv_cloud_scale=inv_cloud_separation_scale, &
                 &         cloud_cover=flux%cloud_cover_lw(jcol))
          else
            ! Two regions with scattering but without 3D effects
            call calc_flux_2region(config%n_g_lw, nlev, emission(:,jcol), albedo(:,jcol), &
                 &         planck_hl(:,:,jcol), cloud%fraction(jcol,:), &
                 &         od(:,:,jcol), od_cloud_regrid, ssa_cloud_regrid, g_cloud_regrid, &
                 &         cloud%overlap_param(jcol,:), &
                 &         flux_up, flux_dn, &
                 &         n_angles_per_hem=config%n_angles_per_hemisphere_lw, &
                 &         cloud_cover=flux%cloud_cover_lw(jcol))
          end if
        else
          if (config%do_3d_effects) then
            ! Three regions with scattering and 3D effects
            call calc_flux_3region(config%n_g_lw, nlev, emission(:,jcol), albedo(:,jcol), &
                 &         planck_hl(:,:,jcol), cloud%fraction(jcol,:), cloud%fractional_std(jcol,:), &
                 &         od(:,:,jcol), od_cloud_regrid, ssa_cloud_regrid, g_cloud_regrid, &
                 &         cloud%overlap_param(jcol,:), &
                 &         flux_up, flux_dn, &
                 &         n_angles_per_hem=config%n_angles_per_hemisphere_lw, &
                 &         layer_thickness=layer_thickness, inv_cloud_scale=inv_cloud_separation_scale, &
                 &         cloud_cover=flux%cloud_cover_lw(jcol))
          else
            ! Three regions with scattering but without 3D effects
            call calc_flux_3region(config%n_g_lw, nlev, emission(:,jcol), albedo(:,jcol), &
                 &         planck_hl(:,:,jcol), cloud%fraction(jcol,:), cloud%fractional_std(jcol,:), &
                 &         od(:,:,jcol), od_cloud_regrid, ssa_cloud_regrid, g_cloud_regrid, &
                 &         cloud%overlap_param(jcol,:), &
                 &         flux_up, flux_dn, &
                 &         n_angles_per_hem=config%n_angles_per_hemisphere_lw, &
                 &         cloud_cover=flux%cloud_cover_lw(jcol))
          end if
        end if
      else
        ! Assume that the optical depth of clouds is the absorption
        ! optical depth
        if (config%nregions == 2) then
          ! Two regions and no scattering or 3D effects
          call calc_no_scattering_flux_2region(config%n_g_lw, nlev, emission(:,jcol), albedo(:,jcol), &
               &         planck_hl(:,:,jcol), cloud%fraction(jcol,:), &
               &         od(:,:,jcol), od_cloud_regrid, &
               &         cloud%overlap_param(jcol,:), &
               &         flux_up, flux_dn, &
               &         n_angles_per_hem=config%n_angles_per_hemisphere_lw, &
               &         do_3d_effects=config%do_3d_effects, &
               &         cloud_cover=flux%cloud_cover_lw(jcol))
        else
          ! Thre regions and no scattering or 3D effects
          call calc_no_scattering_flux_3region(config%n_g_lw, nlev, emission(:,jcol), albedo(:,jcol), &
               &         planck_hl(:,:,jcol), &
               &         cloud%fraction(jcol,:), cloud%fractional_std(jcol,:), &
               &         od(:,:,jcol), od_cloud_regrid, &
               &         cloud%overlap_param(jcol,:), &
               &         flux_up, flux_dn, &
               &         n_angles_per_hem=config%n_angles_per_hemisphere_lw, &
               &         do_3d_effects=config%do_3d_effects, &
               &         cloud_cover=flux%cloud_cover_lw(jcol))
        end if
      end if

      flux%lw_up(jcol,:) = sum(flux_up,1)
      flux%lw_dn(jcol,:) = sum(flux_dn,1)

      if (config%do_save_spectral_flux) then
        call indexed_sum_profile(flux_up, &
             &           config%i_spec_from_reordered_g_lw, &
             &           flux%lw_up_band(:,jcol,:))
        call indexed_sum_profile(flux_dn, &
             &           config%i_spec_from_reordered_g_lw, &
             &           flux%lw_dn_band(:,jcol,:))
      end if
     
      ! Store surface spectral downwelling fluxes
      flux%lw_dn_surf_g(:,jcol) = flux_dn(:,nlev+1)

    end do

    if (lhook) call dr_hook('radiation_tcrad_lw:solver_tcrad_lw',1,hook_handle)

  end subroutine solver_tcrad_lw


  subroutine radiance_solver_tcrad_lw(nlev,istartcol,iendcol, &
       &  config, thermodynamics, cloud, cos_sensor_zenith_angle, &
       &  od, ssa, g, od_cloud, ssa_cloud, g_cloud, planck_hl, &
       &  emission, albedo, flux)

    use parkind1, only           : jprb
    use yomhook,  only           : lhook, dr_hook

!    use radiation_io, only             : nulout
    use radiation_config, only         : config_type, IPdfShapeGamma
    use radiation_thermodynamics, only : thermodynamics_type
    use radiation_cloud, only          : cloud_type
    use radiation_flux, only           : flux_type, indexed_sum
    use radiation_cloudless_lw, only   : solver_cloudless_lw
    use tcrad_2region_solver, only     : calc_radiance_2region => calc_radiance, &
         &                               calc_no_scattering_radiance_2region => calc_no_scattering_radiance
    use tcrad_3region_solver, only     : calc_radiance_3region => calc_radiance, &
         &                               calc_no_scattering_radiance_3region => calc_no_scattering_radiance
    use radiation_constants,  only : GasConstantDryAir, AccelDueToGravity

    implicit none

    real(jprb), parameter :: R_over_g = GasConstantDryAir / AccelDueToGravity

    ! Inputs
    integer, intent(in) :: nlev               ! number of model levels
    integer, intent(in) :: istartcol, iendcol ! range of columns to process
    type(config_type),        intent(in) :: config
    type(thermodynamics_type),intent(in) :: thermodynamics
    type(cloud_type),         intent(in) :: cloud

    real(jprb), intent(in) :: cos_sensor_zenith_angle(istartcol:iendcol)

    ! Gas and aerosol optical depth of each layer at each longwave
    ! g-point
    real(jprb), intent(in), dimension(config%n_g_lw,nlev,istartcol:iendcol) :: od

    ! Gas and aerosol single-scattering albedo and asymmetry factor,
    ! only if longwave scattering by aerosols is to be represented
    real(jprb), intent(in), &
         &  dimension(config%n_g_lw_if_scattering,nlev,istartcol:iendcol) :: ssa, g

    ! Cloud and precipitation optical depth of each layer in each
    ! longwave band
    real(jprb), intent(in) :: od_cloud(config%n_bands_lw,nlev,istartcol:iendcol)

    ! Cloud and precipitation single-scattering albedo and asymmetry
    ! factor, only if longwave scattering by clouds is to be
    ! represented
    real(jprb), intent(in), dimension(config%n_bands_lw_if_scattering, &
         &                            nlev,istartcol:iendcol) :: ssa_cloud, g_cloud

    ! Planck function (emitted flux from a black body) at half levels
    ! and at the surface at each longwave g-point
    real(jprb), intent(in), dimension(config%n_g_lw,nlev+1,istartcol:iendcol) :: planck_hl

    ! Emission (Planck*emissivity) and albedo (1-emissivity) at the
    ! surface at each longwave g-point
    real(jprb), intent(in), dimension(config%n_g_lw, istartcol:iendcol) :: emission, albedo

    ! Output
    type(flux_type), intent(inout):: flux

    ! Local variables

    ! Regridded cloud properties
    real(jprb), dimension(config%n_g_lw,nlev) :: od_cloud_regrid, ssa_cloud_regrid, g_cloud_regrid

    ! Additional geometric properties needed for 3D effects: inverse
    ! of the cloud separation scale (Fielding et al. 20202) (m-1), and
    ! the layer thickness (m)
    real(jprb), dimension(nlev) :: inv_cloud_separation_scale, layer_thickness

    real(jprb), dimension(config%n_g_lw) :: spectral_radiance

    real(jprb) :: hook_handle

    integer :: jcol

    if (lhook) call dr_hook('radiation_tcrad_lw:radiance_solver_tcrad_lw',0,hook_handle)

    if (.not. config%do_lw_cloud_scattering) then
      ! FIX we have allocated memory that is not used
      ssa_cloud_regrid = 0.0_jprb
      g_cloud_regrid   = 0.0_jprb
    end if

    do jcol = istartcol,iendcol

      ! If we do 3D effects then we need to provide the layer
      ! thickness and cloud separation scale
      if (config%do_3d_effects) then
        layer_thickness = R_over_g * 0.5_jprb*(thermodynamics%temperature_hl(jcol,1:nlev) &
             &                                +thermodynamics%temperature_hl(jcol,2:nlev+1)) &
             &               * log(thermodynamics%pressure_hl(jcol,2:nlev+1) &
             &                     /max(thermodynamics%pressure_hl(jcol,1:nlev),1.0e-2_jprb))
        inv_cloud_separation_scale = cloud%inv_cloud_effective_size(jcol,:) &
             &  * sqrt(cloud%fraction(jcol,:)*(1.0_jprb-cloud%fraction(jcol,:)))
      end if

      od_cloud_regrid    = od_cloud(config%i_band_from_reordered_g_lw,:,jcol)
      if (config%do_lw_cloud_scattering) then
        ssa_cloud_regrid = ssa_cloud(config%i_band_from_reordered_g_lw,:,jcol)
        g_cloud_regrid   = g_cloud(config%i_band_from_reordered_g_lw,:,jcol)
        if (config%nregions == 2) then
          if (config%do_3d_effects) then
            call calc_radiance_2region(config%n_g_lw, nlev, emission(:,jcol), albedo(:,jcol), &
                 &         planck_hl(:,:,jcol), cloud%fraction(jcol,:), &
                 &         od(:,:,jcol), od_cloud_regrid, ssa_cloud_regrid, g_cloud_regrid, &
                 &         cloud%overlap_param(jcol,:), &
                 &         cos_sensor_zenith_angle(jcol), spectral_radiance, &
                 &         cloud_cover=flux%cloud_cover_lw(jcol), &
                 &         layer_thickness=layer_thickness, inv_cloud_scale=inv_cloud_separation_scale)
          else
            call calc_radiance_2region(config%n_g_lw, nlev, emission(:,jcol), albedo(:,jcol), &
                 &         planck_hl(:,:,jcol), cloud%fraction(jcol,:), &
                 &         od(:,:,jcol), od_cloud_regrid, ssa_cloud_regrid, g_cloud_regrid, &
                 &         cloud%overlap_param(jcol,:), &
                 &         cos_sensor_zenith_angle(jcol), spectral_radiance, &
                 &         cloud_cover=flux%cloud_cover_lw(jcol))
          end if
        else
          if (config%do_3d_effects) then
            call calc_radiance_3region(config%n_g_lw, nlev, emission(:,jcol), albedo(:,jcol), &
                 &         planck_hl(:,:,jcol), cloud%fraction(jcol,:), cloud%fractional_std(jcol,:), &
                 &         od(:,:,jcol), od_cloud_regrid, ssa_cloud_regrid, g_cloud_regrid, &
                 &         cloud%overlap_param(jcol,:), &
                 &         cos_sensor_zenith_angle(jcol), spectral_radiance, &
                 &         cloud_cover=flux%cloud_cover_lw(jcol), &
                 &         layer_thickness=layer_thickness, inv_cloud_scale=inv_cloud_separation_scale)
          else
            call calc_radiance_3region(config%n_g_lw, nlev, emission(:,jcol), albedo(:,jcol), &
                 &         planck_hl(:,:,jcol), cloud%fraction(jcol,:), cloud%fractional_std(jcol,:), &
                 &         od(:,:,jcol), od_cloud_regrid, ssa_cloud_regrid, g_cloud_regrid, &
                 &         cloud%overlap_param(jcol,:), &
                 &         cos_sensor_zenith_angle(jcol), spectral_radiance, &
                 &         cloud_cover=flux%cloud_cover_lw(jcol))
          end if
        end if
      else
        ! Assume that the optical depth of clouds is the absorption
        ! optical depth
        if (config%nregions == 2) then
          call calc_no_scattering_radiance_2region(config%n_g_lw, nlev, emission(:,jcol), albedo(:,jcol), &
               &         planck_hl(:,:,jcol), cloud%fraction(jcol,:), &
               &         od(:,:,jcol), od_cloud_regrid, &
               &         cloud%overlap_param(jcol,:), &
               &         cos_sensor_zenith_angle(jcol), spectral_radiance, &
               &         do_3d_effects=config%do_3d_effects, &
               &         cloud_cover=flux%cloud_cover_lw(jcol))
        else
          call calc_no_scattering_radiance_3region(config%n_g_lw, nlev, emission(:,jcol), albedo(:,jcol), &
               &         planck_hl(:,:,jcol), &
               &         cloud%fraction(jcol,:), cloud%fractional_std(jcol,:), &
               &         od(:,:,jcol), od_cloud_regrid, &
               &         cloud%overlap_param(jcol,:), &
               &         cos_sensor_zenith_angle(jcol), spectral_radiance, &
               &         do_3d_effects=config%do_3d_effects, &
               &         cloud_cover=flux%cloud_cover_lw(jcol))
        end if
      end if

      if (config%do_save_spectral_flux) then
        flux%lw_radiance_band(:,jcol) = spectral_radiance
      else
        call indexed_sum(spectral_radiance, &
             &           config%i_band_from_reordered_g_lw, &
             &           flux%lw_radiance_band(:,jcol))
      end if
     
    end do

    if (lhook) call dr_hook('radiation_tcrad_lw:radiance_solver_tcrad_lw',1,hook_handle)

  end subroutine radiance_solver_tcrad_lw

end module radiation_tcrad_lw
