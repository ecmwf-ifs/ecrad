! radiation_liquid_optics_nielsen.F90
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
! Author:  Robin Hogan
! Email:   r.j.hogan@ecmwf.int
! License: see the COPYING file for details
!

module radiation_liquid_optics_nielsen

  implicit none

  integer, parameter :: NLiqOpticsCoeffsNielsenSW = 8

contains

  !---------------------------------------------------------------------
  subroutine calc_liq_optics_nielsen(nb, coeff, lwp, re, od, scat_od, g)

    use parkind1, only : jprb
    !use yomhook,  only : lhook, dr_hook

    ! Number of bands
    integer, intent(in)  :: nb
    ! Coefficients read from a data file
    real(jprb), intent(in) :: coeff(:,:)
    ! Liquid water path (kg m-2) and effective radius (m)
    real(jprb), intent(in) :: lwp, re
    ! Total optical depth, scattering optical depth and asymmetry factor
    real(jprb), intent(out) :: od(nb), scat_od(nb), g(nb)

    ! Liquid water path in g m-2
    real(jprb) :: lwp_gm_2
    ! Effective radius in microns, and its inverse
    real(jprb) :: re_um, inv_re_um
    !real(jprb) :: hook_handle

    !if (lhook) call dr_hook('radiation_liquid_optics_socrates:calc_liq_optics_socrates',0,hook_handle)

    lwp_gm_2 = lwp * 1000.0_jprb
    re_um= re*1e6_jprb
    inv_re_um = 1.0_jprb / re_um

    od=lwp_gm_2 * coeff(1:nb,1)*(re_um**coeff(1:nb,2))
    scat_od=od*(coeff(1:nb,3)+coeff(1:nb,4)*re_um)
    g=coeff(1:nb,5)+coeff(1:nb,6)*re_um+(coeff(1:nb,7)*EXP(coeff(1:nb,8)*re_um))

    !if (lhook) call dr_hook('radiation_liquid_optics_socrates:calc_liq_optics_socrates',1,hook_handle)

  end subroutine calc_liq_optics_nielsen

end module radiation_liquid_optics_nielsen
