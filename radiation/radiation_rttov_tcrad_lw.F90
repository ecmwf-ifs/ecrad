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

!  public

module radiation_rttov_tcrad_lw

  use iso_c_binding

contains

  subroutine radiance_solver_rttov_tcrad_lw(nlev,istartcol,iendcol, &
       &  config, thermodynamics, cloud, cos_sensor_zenith_angle, &
       &  od, ssa, g, od_cloud, ssa_cloud, g_cloud, planck_hl, &
       &  emission, albedo, flux)

    use parkind1, only           : jprb
    use yomhook,  only           : lhook, dr_hook, jphook

!    use radiation_io, only             : nulout
    use radiation_config, only         : config_type, IPdfShapeGamma
    use radiation_thermodynamics, only : thermodynamics_type
    use radiation_cloud, only          : cloud_type
    use radiation_flux, only           : flux_type, indexed_sum
    use radiation_cloudless_lw, only   : solver_cloudless_lw
    use tcrad, only                    : calc_radiance
    use tcrad_clear_sky, only          : calc_radiance_clear_sky
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

    real(jprb), dimension(config%n_g_lw) :: spectral_radiance, spectral_radiance_clear
    
    integer :: jcol

    real(jphook) :: hook_handle

    if (lhook) call dr_hook('radiation_tcrad_lw:radiance_solver_tcrad_lw',0,hook_handle)

    if (.not. config%do_lw_cloud_scattering) then
      ! FIX we have allocated memory that is not used
      ssa_cloud_regrid = 0.0_jprb
      g_cloud_regrid   = 0.0_jprb
    end if

    do jcol = istartcol,iendcol
      
      od_cloud_regrid    = od_cloud(config%i_band_from_reordered_g_lw,:,jcol)
      if (config%do_lw_cloud_scattering) then
        ssa_cloud_regrid = ssa_cloud(config%i_band_from_reordered_g_lw,:,jcol)
        g_cloud_regrid   = g_cloud(config%i_band_from_reordered_g_lw,:,jcol)
        
        call calc_radiance(config%n_g_lw, nlev,  &
             &         emission(:,jcol), albedo(:,jcol), &
             &         planck_hl(:,:,jcol), cloud%fraction(jcol,:), cloud%fractional_std(jcol,:), &
             &         od(:,:,jcol), od_cloud_regrid, ssa_cloud_regrid, g_cloud_regrid, &
             &         cloud%overlap_param(jcol,:), cos_sensor_zenith_angle(jcol), &
             &         spectral_radiance, &
             &         cloud_cover=flux%cloud_cover_lw(jcol), &
             &         do_specular_surface=config%do_specular_surface)

      end if

      if (config%do_save_spectral_flux) then
        flux%lw_radiance_band(:,jcol) = spectral_radiance
      else
        call indexed_sum(spectral_radiance, &
             &           config%i_band_from_reordered_g_lw, &
             &           flux%lw_radiance_band(:,jcol))
      end if

      if (config%do_clear) then
        call calc_radiance_clear_sky(config%n_g_lw, nlev,  &
             &         emission(:,jcol), albedo(:,jcol), &
             &         planck_hl(:,:,jcol), od(:,:,jcol), &
             &         cos_sensor_zenith_angle(jcol), &
             &         spectral_radiance_clear, &
             &         do_specular_surface=config%do_specular_surface)
        if (config%do_save_spectral_flux) then
          flux%lw_radiance_clear_band(:,jcol) = spectral_radiance_clear
        else
          call indexed_sum(spectral_radiance_clear, &
               &           config%i_band_from_reordered_g_lw, &
               &           flux%lw_radiance_clear_band(:,jcol))
        end if
      end if
      
    end do
    
    if (lhook) call dr_hook('radiation_tcrad_lw:radiance_solver_tcrad_lw',1,hook_handle)
    
  end subroutine radiance_solver_rttov_tcrad_lw
    
  subroutine radiance_solver_rttov_tcrad_lw_ad(nlev,istartcol,iendcol, &
       &  config, thermodynamics, cloud, cos_sensor_zenith_angle, &
       &  od, ssa, g, od_cloud, ssa_cloud, g_cloud, planck_hl, &
       &  emission, albedo, flux)

    use parkind1, only           : jprb
    use yomhook,  only           : lhook, dr_hook, jphook

!    use radiation_io, only             : nulout
    use radiation_config, only         : config_type, IPdfShapeGamma
    use radiation_thermodynamics, only : thermodynamics_type
    use radiation_cloud, only          : cloud_type
    use radiation_flux, only           : flux_type, indexed_sum
    use radiation_cloudless_lw, only   : solver_cloudless_lw
    use tcrad_ad, only                    : calc_radiance_ad
    use radiation_constants,  only : GasConstantDryAir, AccelDueToGravity

    implicit none

#include "tcrad.inc"

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


    ! Adjoints

    real(jprb), dimension(config%n_g_lw) :: spectral_radiance_ad, emission_ad, albedo_ad;
    real(jprb), dimension(config%n_g_lw,nlev+1) :: planck_hl_ad
    real(jprb), dimension(nlev) :: cloud_fraction_ad, fractional_std_ad
    real(jprb), dimension(config%n_g_lw,nlev) :: od_ad, od_cloud_ad, ssa_cloud_ad, g_cloud_ad
    real(jprb), dimension(nlev-1) :: overlap_param_ad
    
    integer :: jcol, jlev

    real(jphook) :: hook_handle

    if (lhook) call dr_hook('radiation_tcrad_lw:radiance_solver_tcrad_lw_ad',0,hook_handle)

    if (.not. config%do_lw_cloud_scattering) then
      ! FIX we have allocated memory that is not used
      ssa_cloud_regrid = 0.0_jprb
      g_cloud_regrid   = 0.0_jprb
    end if
    
    do jcol = istartcol,iendcol
      
      od_cloud_regrid    = od_cloud(config%i_band_from_reordered_g_lw,:,jcol)
      if (config%do_lw_cloud_scattering) then
        ssa_cloud_regrid = ssa_cloud(config%i_band_from_reordered_g_lw,:,jcol)
        g_cloud_regrid   = g_cloud(config%i_band_from_reordered_g_lw,:,jcol)

        
        spectral_radiance_ad = 0.0
        spectral_radiance_ad(1) = 1.0
        emission_ad=0.0
        albedo_ad=0.0
        planck_hl_ad=0.0
        cloud_fraction_ad=0.0
        od_ad=0.0
        od_cloud_ad=0.0
        ssa_cloud_ad=0.0
        g_cloud_ad=0.0
        call calc_radiance_ad(config%n_g_lw, nlev,  &
             &         emission(:,jcol), albedo(:,jcol), &
             &         planck_hl(:,:,jcol), cloud%fraction(jcol,:), cloud%fractional_std(jcol,:), &
             &         od(:,:,jcol), od_cloud_regrid, ssa_cloud_regrid, g_cloud_regrid, &
             &         cloud%overlap_param(jcol,:), cos_sensor_zenith_angle(jcol), &
             &         spectral_radiance,  &
             &         emission_ad, albedo_ad, planck_hl_ad, cloud_fraction_ad, fractional_std_ad, od_ad, od_cloud_ad, &
             &         ssa_cloud_ad, g_cloud_ad, &
             &         overlap_param_ad, spectral_radiance_ad, &
             &         cloud_cover=flux%cloud_cover_lw(jcol), do_specular_surface=config%do_specular_surface)

        ! Write to unit 101 the gradients of the radiance with respect to the inputs
        do jlev = 1,nlev
          write(101,*) jcol, jlev, spectral_radiance(1), planck_hl_ad(1,jlev), cloud_fraction_ad(jlev), od_ad(1,jlev), &
               &  od_cloud_ad(1,jlev), ssa_cloud_ad(1,jlev), g_cloud_ad(1,jlev)
        end do
        
      end if

      if (config%do_save_spectral_flux) then
        flux%lw_radiance_band(:,jcol) = spectral_radiance
      else
        call indexed_sum(spectral_radiance, &
             &           config%i_band_from_reordered_g_lw, &
             &           flux%lw_radiance_band(:,jcol))
      end if
      
    end do

    if (lhook) call dr_hook('radiation_tcrad_lw:radiance_solver_tcrad_lw_ad',1,hook_handle)

  end subroutine radiance_solver_rttov_tcrad_lw_ad
  
end module radiation_rttov_tcrad_lw
