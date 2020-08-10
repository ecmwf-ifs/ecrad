! radiation_liquid_optics_socrates.F90 - SOCRATES method for parameterizing liquid droplet optics
!
! Copyright (C) 2014-2020 ECMWF
!
! Author:  Robin Hogan
! Email:   r.j.hogan@ecmwf.int
! License: see the COPYING file for details
!
! Modifications
!   2020-08-10  R. Hogan  Bounded re to be >=1.2um and <=50um

module radiation_liquid_optics_socrates

  use parkind1, only : jprb

  implicit none
  public

  ! SOCRATES (Edwards-Slingo) parameterizes info on the dependence of
  ! the scattering properties in each band on effective radius in
  ! terms of 16 coefficients
  integer, parameter :: NLiqOpticsCoeffsSOCRATES = 16

  ! Range of valid input effective radius, in microns
  real(jprb), parameter :: MinEffectiveRadius = 1.2e-6
  real(jprb), parameter :: MaxEffectiveRadius = 50.0e-6

contains

  !---------------------------------------------------------------------
  ! Compute liquid-droplet scattering properties using a
  ! parameterization consisting of Pade approximants from the
  ! SOCRATES (Edwards-Slingo) code
  subroutine calc_liq_optics_socrates(nb, coeff, lwp, re_in, od, scat_od, g)

    use parkind1, only : jprb
    !use yomhook,  only : lhook, dr_hook

    ! Number of bands
    integer, intent(in)  :: nb
    ! Coefficients read from a data file
    real(jprb), intent(in) :: coeff(:,:)
    ! Liquid water path (kg m-2) and effective radius (m)
    real(jprb), intent(in) :: lwp, re_in
    ! Total optical depth, scattering optical depth and asymmetry factor
    real(jprb), intent(out) :: od(nb), scat_od(nb), g(nb)

    ! Local effective radius (m), after applying bounds
    real(jprb) :: re

    !real(jprb) :: hook_handle

    !if (lhook) call dr_hook('radiation_liquid_optics_socrates:calc_liq_optics_socrates',0,hook_handle)

    ! Apply the bounds of validity to effective radius
    re = max(MinEffectiveRadius, min(re_in, MaxEffectiveRadius))

    od = lwp * (coeff(1:nb,1) + re*(coeff(1:nb,2) + re*coeff(1:nb,3))) &
         &  / (1.0_jprb + re*(coeff(1:nb,4) + re*(coeff(1:nb,5) &
         &  + re*coeff(1:nb,6))))
    scat_od = od * (1.0_jprb &
         &  - (coeff(1:nb,7) + re*(coeff(1:nb,8) + re*coeff(1:nb,9))) &
         &  / (1.0_jprb + re*(coeff(1:nb,10) + re*coeff(1:nb,11))))
    g = (coeff(1:nb,12) + re*(coeff(1:nb,13) + re*coeff(1:nb,14))) &
         &  / (1.0_jprb + re*(coeff(1:nb,15) + re*coeff(1:nb,16)))

    !if (lhook) call dr_hook('radiation_liquid_optics_socrates:calc_liq_optics_socrates',1,hook_handle)

  end subroutine calc_liq_optics_socrates

end module radiation_liquid_optics_socrates
