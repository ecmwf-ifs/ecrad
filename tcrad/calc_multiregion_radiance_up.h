! calc_multiregion_radiance_up.h - Upward radiances in multi-region profile -*- f90 -*-
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

subroutine calc_multiregion_radiance_up(nspec, nlev, surf_emission, &
     &  transmittance, source_up, u_overlap, radiance_up)

  use parkind1, only           : jpim, jprb
  use yomhook,  only           : lhook, dr_hook
  
  implicit none
  
  ! Inputs

  ! Number of spectral intervals and levels
  integer(jpim), intent(in) :: nspec, nlev

  real(jprb), intent(in),  dimension(nspec) :: surf_emission

  real(jprb), intent(in),  dimension(nspec,NREGION,nlev) :: transmittance
  real(jprb), intent(in),  dimension(nspec,NREGION,nlev) :: source_up

  real(jprb), intent(in),  dimension(NREGION,NREGION,nlev) :: u_overlap

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

  if (lhook) call dr_hook('calc_multiregion_radiance_up',0,hook_handle)

  do jreg = 1,NREGION
    radiance_base(:,jreg) = u_overlap(jreg,1,nlev+1)*surf_emission(:)
  end do

  radiance_up(:,nlev+1) = radiance_up(:,nlev+1) + sum(radiance_base,2)

  do jlev = nlev,1,-1
    radiance_top = transmittance(:,:,jlev)*radiance_base + source_up(:,:,jlev)
    radiance_base = singlemat_x_vec(nspec,u_overlap(:,:,jlev),radiance_top)
    radiance_up(:,jlev) = radiance_up(:,jlev) + sum(radiance_base,2)
  end do

  if (lhook) call dr_hook('calc_multiregion_radiance_up',1,hook_handle)

end subroutine calc_multiregion_radiance_up
