! tcrad_clear_sky.F90 - Clear-sky no-scattering radiance model
!
! (C) Copyright 2014- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
!
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!
! This file provides the capability to compute thermal-infrared and
! microwave radiances in all-sky conditions using the Tripleclouds
! approximation for treating cloud heterogeneity.
!
! Author:  Robin Hogan
! Email:   r.j.hogan@ecmwf.int
!

module tcrad_clear_sky

contains
  
!---------------------------------------------------------------------
! Compute the TOA clear-sky radiance neglecting atmospheric
! scattering. Note that this is not very memory efficient, but that's
! so that it is easier to generate the adjoint from this code.
subroutine calc_radiance_clear_sky(nspec, nlev, surf_emission, surf_albedo, &
     &  planck_hl, od, mu, radiance, do_specular_surface)

  use parkind1, only : jprb, jpim
  use tcrad,    only : lw_diffusivity, OD_THRESH_2STREAM
  use yomhook,  only : lhook, dr_hook, jphook
  
  implicit none

  real(jprb), parameter :: PI          = acos(-1.0_jprb)
  real(jprb), parameter :: ONE_OVER_PI = 1.0_jprb / PI

  ! Inputs

  ! Number of spectral intervals and levels. Note that all
  ! level-dependent variables count down from top-of-atmosphere.
  integer(jpim), intent(in) :: nspec, nlev

  ! Surface upwards emission, in W m-2 (i.e. emissivity multiplied
  ! by Planck function at the surface skin temperature) integrated
  ! across each spectral interval, and albedo in the same intervals
  real(jprb), intent(in), dimension(nspec) :: surf_emission, surf_albedo

  ! Planck function integrated over each spectral interval at each
  ! half-level, in W m-2 (i.e. the flux emitted by a horizontal
  ! black-body surface)
  real(jprb), intent(in), dimension(nspec,nlev+1) :: planck_hl

  ! Layer optical depth of gas and aerosol
  real(jprb), intent(in), dimension(nspec,nlev) :: od

  ! Cosine of the sensor zenith angle
  real(jprb), intent(in) :: mu

  ! Outputs

  ! Radiances in W m sr-1
  real(jprb), intent(out), dimension(nspec) :: radiance

  ! Optional inputs

  ! Is scattering from the surface treated specularly?  Appropriate
  ! for microwave scattering by the sea surface. Default is NO.
  logical, intent(in), optional :: do_specular_surface

  ! Local variables

  ! Transmittance of each layer for the downward and upward paths. In
  ! the case of specular reflection they are the same and so not
  ! recomputed. "transmittance" is up in all cases and also down in
  ! the case of specular reflection when the path through the
  ! atmosphere is the same.
  real(jprb), dimension(nspec,nlev) :: transmittance, transmittance_dn

  ! Rate of emission up from the top or down through the base of
  ! each layer and region (W m-2)
  real(jprb), dimension(nspec,nlev) :: source_up, source_dn

  ! Profile of radiances in direction of sensor
  real(jprb), dimension(nspec, nlev+1) :: radiance_dn, radiance_up

  ! Cosine and secants of the downward and upward propagating radiances
  real(jprb) :: secant_dn, secant_up

  ! Working variables
  real(jprb) :: coeff, coeff_dn_base, coeff_dn_top, coeff_up_base, coeff_up_top

  integer(jpim) :: jlev, jspec
  
  ! Is scattering from the surface treated specularly?  Appropriate
  ! for microwave scattering by the sea surface.
  logical :: do_specular_surface_local

  real(jphook) :: hook_handle

  if (lhook) call dr_hook('tcrad_clear_sky:calc_radiance_clear_sky',0,hook_handle)

  ! By default we do not treat the surface as a specular reflector
  if (present(do_specular_surface)) then
    do_specular_surface_local = do_specular_surface
  else
    do_specular_surface_local = .false.
  endif

  secant_up = 1.0_jprb / mu
  if (do_specular_surface_local) then
    ! For specular reflectance at the surface we calculate the
    ! downward radiance at the same angle as the upward
    secant_dn = secant_up
  else
    ! For Lambertian reflectance at the surface we want the downward
    ! flux, so use an angle consistent with the two-stream
    ! approximation
    secant_dn = lw_diffusivity
  end if

  ! Zero top-of-atmosphere radiance
  radiance_dn(:,1) = 0.0_jprb

  if (do_specular_surface) then
    ! Loop down through the atmosphere computing the profile of
    ! transmittances and sources
    do jlev = 1,nlev
      ! No-scattering solution in clear-sky region: compute emission
      ! assuming the Planck function to vary linearly with optical depth
      ! within the layer (e.g. Wiscombe, JQSRT 1976).
      do jspec = 1,nspec
        coeff = secant_dn*od(jspec,jlev)
        if (od(jspec,jlev) > OD_THRESH_2STREAM) then
          transmittance(jspec,jlev) = exp(-coeff)
          coeff = (planck_hl(jspec,jlev+1)-planck_hl(jspec,jlev)) / coeff
          coeff_up_top  =  coeff + planck_hl(jspec,jlev)
          coeff_up_base =  coeff + planck_hl(jspec,jlev+1)
          coeff_dn_top  = -coeff + planck_hl(jspec,jlev)
          coeff_dn_base = -coeff + planck_hl(jspec,jlev+1)
          source_up(jspec,jlev) = coeff_up_top &
               &  - transmittance(jspec,jlev) * coeff_up_base
          source_dn(jspec,jlev) = coeff_dn_base &
               &  - transmittance(jspec,jlev) * coeff_dn_top
        else
          ! Linear limit at low optical depth
          transmittance(jspec,jlev) = 1.0_jprb - coeff
          source_up(jspec,jlev) = coeff * 0.5_jprb &
               &  * (planck_hl(jspec,jlev)+planck_hl(jspec,jlev+1))
          source_dn(jspec,jlev) = source_up(jspec,jlev)
        end if
      end do
    end do

    ! Calculate radiance profile downwards    
    do jlev = 1,nlev
      ! No-scattering solution in clear-sky region: compute emission
      ! assuming the Planck function to vary linearly with optical depth
      ! within the layer (e.g. Wiscombe, JQSRT 1976).
      radiance_dn(:,jlev+1) = transmittance(:,jlev)*radiance_dn(:,jlev) &
           &  + ONE_OVER_PI * source_dn(:,jlev)
    end do
    
  else
    ! Loop down through the atmosphere computing the profile of
    ! transmittances and sources
    do jlev = 1,nlev
      ! No-scattering solution in clear-sky region: compute emission
      ! assuming the Planck function to vary linearly with optical depth
      ! within the layer (e.g. Wiscombe, JQSRT 1976).
      do jspec = 1,nspec
        if (od(jspec,jlev) > OD_THRESH_2STREAM) then
          ! Downward part 
          coeff = secant_dn*od(jspec,jlev)
          transmittance_dn(jspec,jlev) = exp(-coeff)
          coeff = (planck_hl(jspec,jlev+1)-planck_hl(jspec,jlev)) / coeff
          !coeff_up_top  =  coeff + planck_hl(jspec,jlev)
          !coeff_up_base =  coeff + planck_hl(jspec,jlev+1)
          coeff_dn_top  = -coeff + planck_hl(jspec,jlev)
          coeff_dn_base = -coeff + planck_hl(jspec,jlev+1)
          source_dn(jspec,jlev) = coeff_dn_base &
               &  - transmittance_dn(jspec,jlev) * coeff_dn_top
          ! Upward part
          coeff = secant_up*od(jspec,jlev)
          transmittance(jspec,jlev) = exp(-coeff)
          coeff = (planck_hl(jspec,jlev+1)-planck_hl(jspec,jlev)) / coeff
          coeff_up_top  =  coeff + planck_hl(jspec,jlev)
          coeff_up_base =  coeff + planck_hl(jspec,jlev+1)
          !coeff_dn_top  = -coeff + planck_hl(jspec,jlev)
          !coeff_dn_base = -coeff + planck_hl(jspec,jlev+1)
          source_up(jspec,jlev) = coeff_up_top &
               &  - transmittance(jspec,jlev) * coeff_up_base
        else
          ! Linear limit at low optical depth: downward part
          coeff = secant_dn*od(jspec,jlev)
          transmittance_dn(jspec,jlev) = 1.0_jprb - coeff
          source_dn(jspec,jlev) = coeff * 0.5_jprb &
               &  * (planck_hl(jspec,jlev)+planck_hl(jspec,jlev+1))
          ! Upward part
          coeff = secant_up*od(jspec,jlev)
          transmittance(jspec,jlev) = 1.0_jprb - coeff
          source_up(jspec,jlev) = coeff * 0.5_jprb &
               &  * (planck_hl(jspec,jlev)+planck_hl(jspec,jlev+1))
        end if
      end do
    end do

    ! Calculate radiance profile downwards
    do jlev = 1,nlev
      ! No-scattering solution in clear-sky region: compute emission
      ! assuming the Planck function to vary linearly with optical depth
      ! within the layer (e.g. Wiscombe, JQSRT 1976).
      radiance_dn(:,jlev+1) = transmittance_dn(:,jlev)*radiance_dn(:,jlev) &
           &  + ONE_OVER_PI * source_dn(:,jlev)
    end do
  end if

  radiance_up(:,nlev+1) = surf_emission*ONE_OVER_PI + surf_albedo*radiance_dn(:,nlev+1)
  do jlev = nlev,1,-1
    radiance_up(:,jlev) = transmittance(:,jlev)*radiance_up(:,jlev+1) &
         &  + ONE_OVER_PI * source_up(:,jlev)
  end do
  
  radiance = radiance_up(:,1)
  
  if (lhook) call dr_hook('tcrad_clear_sky:calc_radiance_clear_sky',1,hook_handle)

end subroutine calc_radiance_clear_sky

end module tcrad_clear_sky
