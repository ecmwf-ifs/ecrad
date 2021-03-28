! tcrad_clear_sky_radiance.h - Calculate clear-sky radiances in TCRAD
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

module tcrad_clear_sky_radiance

contains

  !---------------------------------------------------------------------
  ! Compute clear-sky upward radiance profile by solving the
  ! Schwarzschild radiative transfer equation assuming the source term
  ! to vary linearly with optical depth in each layer. This routine adds
  ! to any existing radiance profile, useful if used as part of a
  ! multi-stream flux calculation, e.g. the delta-2-plus-4 method of Fu
  ! et al. (1997).
  subroutine calc_clear_sky_radiance_up(nspec, nlev, &
       &  weight, surf_up, &
       &  transmittance, source_up, radiance_up)

    use parkind1, only           : jpim, jprb
    use yomhook,  only           : lhook, dr_hook

    implicit none

    ! Inputs

    ! Number of spectral intervals and levels
    integer(jpim), intent(in) :: nspec, nlev

    ! Weight sources by this amount
    real(jprb), intent(in) :: weight

    ! Surface upwelling flux in W m-2
    real(jprb), intent(in),  dimension(nspec) :: surf_up

    ! Transmittance of each layer in the direction of the radiance; this
    ! does not include diffuse transmittance, i.e. rays that may be
    ! scattered as they pass through the layer.  We use an assumed-shape
    ! array because this argument may have been sliced from a
    ! multi-region array in which case the level dimension will not be
    ! contiguous.
    real(jprb), intent(in),  dimension(:,:) :: transmittance

    ! Upward source from the top of the layer in the direction of the
    ! radiance, which may include Planck emission, and scattering. We
    ! use an assumed-shape array because this argument may have been
    ! sliced from a multi-region array in which case the level dimension
    ! will not be contiguous.
    real(jprb), intent(in),  dimension(:,:) :: source_up ! (nspec,nlev)

    ! Output

    ! Upward radiance profile: note that we add to any existing radiance
    ! profile, useful when summing over multiple angles to get a flux
    real(jprb), intent(inout), dimension(nspec,nlev+1) :: radiance_up

    ! Local variables

    ! Radiance per region at base and top of each layer
    real(jprb), dimension(nspec) :: radiance_base, radiance_top

    ! Loop index for level
    integer(jpim) :: jlev

    real(jprb) :: hook_handle

    if (lhook) call dr_hook('tcrad:calc_clear_sky_radiance_up',0,hook_handle)

    ! Set surface upward radiance
    radiance_base = weight * surf_up

    ! Save radiance profile averaged over regions, adding to existing
    ! values
    radiance_up(:,nlev+1) = radiance_up(:,nlev+1) + radiance_base

    do jlev = nlev,1,-1
      ! Solution to Schwarzschild equation
      radiance_top = transmittance(:,jlev)*radiance_base + weight * source_up(:,jlev)
      ! Save radiances
      radiance_up(:,jlev) = radiance_up(:,jlev) + radiance_base
    end do

    if (lhook) call dr_hook('tcrad:calc_radiance_clear_sky_up',1,hook_handle)

  end subroutine calc_clear_sky_radiance_up


  !---------------------------------------------------------------------
  ! Compute clear-sky downward radiance profile by solving the
  ! Schwarzschild radiative transfer equation assuming the source term
  ! to vary linearly with optical depth in each layer. This routine adds
  ! to any existing radiance profile, useful if used as part of a
  ! multi-stream flux calculation, e.g. the delta-2-plus-4 method of Fu
  ! et al. (1997).
  subroutine calc_clear_sky_radiance_dn(nspec, nlev, &
       &  weight, transmittance, source_dn, radiance_dn)

    use parkind1, only           : jpim, jprb
    use yomhook,  only           : lhook, dr_hook

    implicit none

    ! Inputs

    ! Number of spectral intervals and levels
    integer(jpim), intent(in) :: nspec, nlev

    ! Weight sources by this amount
    real(jprb), intent(in) :: weight

    ! Transmittance of each layer in the direction of the radiance; this
    ! does not include diffuse transmittance, i.e. rays that may be
    ! scattered as they pass through the layer.  We use an assumed-shape
    ! array because this argument may have been sliced from a
    ! multi-region array in which case the level dimension will not be
    ! contiguous.
    real(jprb), intent(in),  dimension(:,:) :: transmittance ! (nspec,nlev)

    ! Down source from the base of the layer in the direction of the
    ! radiance, which may include Planck emission, and scattering.  We
    ! use an assumed-shape array because this argument may have been
    ! sliced from a multi-region array in which case the level dimension
    ! will not be contiguous.
    real(jprb), intent(in),  dimension(:,:) :: source_dn ! (nspec,nlev)

    ! Output

    ! Downward radiance profile: note that we add to any existing radiance
    ! profile, useful when summing over multiple angles to get a flux
    real(jprb), intent(inout), dimension(nspec,nlev+1) :: radiance_dn

    ! Local variables

    ! Radiance per region at base and top of each layer
    real(jprb), dimension(nspec) :: radiance_base, radiance_top

    ! Loop index for level
    integer(jpim) :: jlev

    real(jprb) :: hook_handle

    if (lhook) call dr_hook('tcrad:calc_ckear_sky_radiance_dn',0,hook_handle)

    ! Start with zero at TOA
    radiance_top = 0.0_jprb

    do jlev = 1,nlev
      ! Solution to Schwarzschild equation
      radiance_base = transmittance(:,jlev)*radiance_top + weight * source_dn(:,jlev)
      ! Save radiances
      radiance_dn(:,jlev+1) = radiance_dn(:,jlev+1) + radiance_top
    end do

    if (lhook) call dr_hook('tcrad:calc_clear_sky_radiance_dn',1,hook_handle)

  end subroutine calc_clear_sky_radiance_dn

end module tcrad_clear_sky_radiance

