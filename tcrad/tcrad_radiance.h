! tcrad_radiance.h - Calculate radiances in TCRAD -*- f90 -*-
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
! This file is included in the modules specifying the NREGION
! parameter (typically 2 or 3) which makes this routine use a
! "doubleclouds" or "tripleclouds" assumption.
!

!---------------------------------------------------------------------
! Compute upward radiance profile by solving the Schwarzschild
! radiative transfer equation in each region of each layer assuming
! the source term to vary linearly with optical depth in each
! layer. The overlap matrix is used to translate the radiances exiting
! the top of one layer into radiances in the base of the layer above,
! consistent with the Tripleclouds approximation. This routine adds to
! any existing radiance profile, useful if used as part of a
! multi-stream flux calculation, e.g. the delta-2-plus-4 method of Fu
! et al. (1997).
subroutine calc_radiance_up(nspec, nlev, &
     &  weight, surf_up, &
     &  transmittance, source_up, u_overlap, radiance_up)

  use parkind1, only           : jpim, jprb
  use yomhook,  only           : lhook, dr_hook
  
  implicit none
  
  ! Inputs

  ! Number of spectral intervals and levels
  integer(jpim), intent(in) :: nspec, nlev

  ! Weight sources by this amount
  real(jprb), intent(in) :: weight

  ! Surface upwelling flux in W m-2
  real(jprb), intent(in),  dimension(nspec,NREGION) :: surf_up

  ! Transmittance of each layer and region in the direction of the
  ! radiance; this does not include diffuse transmittance, i.e. rays
  ! that may be scattered as they pass through the layer
  real(jprb), intent(in),  dimension(nspec,NREGION,nlev) :: transmittance

  ! Upward source from the top of the layer in the direction of the
  ! radiance, which may include Planck emission, and scattering
  real(jprb), intent(in),  dimension(nspec,NREGION,nlev) :: source_up

  ! Upward overlap matrix - see Hogan et al. (JGR 2016) for definition
  real(jprb), intent(in),  dimension(NREGION,NREGION,nlev+1) :: u_overlap

  ! Output

  ! Upward radiance profile: note that we add to any existing radiance
  ! profile, useful when summing over multiple angles to get a flux
  real(jprb), intent(inout), dimension(nspec,nlev+1) :: radiance_up

  ! Local variables

  ! Radiance per region at base and top of each layer
  real(jprb), dimension(nspec,NREGION) :: radiance_base, radiance_top

  ! Loop indices for spectral interval, level and region
  integer(jpim) :: jspec, jlev, jreg

  real(jprb) :: hook_handle

  if (lhook) call dr_hook('tcrad:calc_radiance_up',0,hook_handle)

  ! Set surface upward radiance
  radiance_base = weight * surf_up

  ! Save radiance profile averaged over regions, adding to existing
  ! values
  radiance_up(:,nlev+1) = radiance_up(:,nlev+1) + sum(radiance_base,2)

  do jlev = nlev,1,-1
    ! Solution to Schwarzschild equation
    radiance_top = transmittance(:,:,jlev)*radiance_base + weight * source_up(:,:,jlev)
    ! Overlap rules to obtain radiances at base of the layer above
    radiance_base = singlemat_x_vec(nspec,u_overlap(:,:,jlev),radiance_top)
    ! Save radiances
    radiance_up(:,jlev) = radiance_up(:,jlev) + sum(radiance_base,2)
  end do

  if (lhook) call dr_hook('tcrad:calc_radiance_up',1,hook_handle)

end subroutine calc_radiance_up


!---------------------------------------------------------------------
! Compute downward radiance profile by solving the Schwarzschild
! radiative transfer equation in each region of each layer assuming
! the source term to vary linearly with optical depth in each
! layer. The overlap matrix is used to translate the radiances exiting
! the base of one layer into radiances in the top of the layer below,
! consistent with the Tripleclouds approximation. This routine adds to
! any existing radiance profile, useful if used as part of a
! multi-stream flux calculation, e.g. the delta-2-plus-4 method of Fu
! et al. (1997).
subroutine calc_radiance_dn(nspec, nlev, &
     &  weight, transmittance, source_dn, v_overlap, radiance_dn)

  use parkind1, only           : jpim, jprb
  use yomhook,  only           : lhook, dr_hook
  
  implicit none
  
  ! Inputs

  ! Number of spectral intervals and levels
  integer(jpim), intent(in) :: nspec, nlev

  ! Weight sources by this amount
  real(jprb), intent(in) :: weight

  ! Transmittance of each layer and region in the direction of the
  ! radiance; this does not include diffuse transmittance, i.e. rays
  ! that may be scattered as they pass through the layer
  real(jprb), intent(in),  dimension(nspec,NREGION,nlev) :: transmittance

  ! Down source from the base of the layer in the direction of the
  ! radiance, which may include Planck emission, and scattering
  real(jprb), intent(in),  dimension(nspec,NREGION,nlev) :: source_dn

  ! Downward overlap matrix - see Shonk and Hogan (2008) for definition
  real(jprb), intent(in),  dimension(NREGION,NREGION,nlev+1) :: v_overlap

  ! Output

  ! Downward radiance profile: note that we add to any existing radiance
  ! profile, useful when summing over multiple angles to get a flux
  real(jprb), intent(inout), dimension(nspec,nlev+1) :: radiance_dn

  ! Local variables

  ! Radiance per region at base and top of each layer
  real(jprb), dimension(nspec,NREGION) :: radiance_base, radiance_top

  ! Loop indices for spectral interval, level and region
  integer(jpim) :: jspec, jlev, jreg

  real(jprb) :: hook_handle

  if (lhook) call dr_hook('tcrad:calc_radiance_dn',0,hook_handle)

  ! Start with zero at TOA
  radiance_top = 0.0_jprb

  do jlev = 1,nlev
    ! Solution to Schwarzschild equation
    radiance_base = transmittance(:,:,jlev)*radiance_top + weight * source_dn(:,:,jlev)
    ! Overlap rules to obtain radiances at base of the layer above
    radiance_top = singlemat_x_vec(nspec,v_overlap(:,:,jlev+1),radiance_base)
    ! Save radiances
    radiance_dn(:,jlev+1) = radiance_dn(:,jlev+1) + sum(radiance_top,2)
  end do

  if (lhook) call dr_hook('tcrad:calc_radiance_dn',1,hook_handle)

end subroutine calc_radiance_dn
