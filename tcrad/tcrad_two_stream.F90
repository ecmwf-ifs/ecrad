! tcrad_two_stream.F90 - Two-stream and related layer solutions for TCRAD package
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

module tcrad_two_stream

  use parkind1, only : jprb

  ! Elsasser's factor: the effective factor by which the zenith
  ! optical depth needs to be multiplied to account for longwave
  ! transmission at all angles through the atmosphere.  Alternatively
  ! think of acos(1/LW_DIFFUSIVITY) to be the effective zenith angle
  ! of longwave radiation.
  real(jprb), parameter :: LW_DIFFUSIVITY = 1.66_jprb

contains

  subroutine calc_reflectance_transmittance(nspec, nlev, nreg, &
       &  is_cloud_free_layer, planck_hl, od, ssa, asymmetry, &
       &  reflectance, transmittance, source_up, source_dn)
        
    use parkind1, only           : jpim, jprb
    use yomhook,  only           : lhook, dr_hook

    implicit none

    ! Inputs

    ! Number of spectral intervals, levels and regions
    integer(jpim), intent(in) :: nspec, nlev, nreg

    logical,    intent(in), dimension(nlev)         :: is_cloud_free_layer

    real(jprb), intent(in), dimension(nspec,nlev+1) :: planck_hl

    real(jprb), intent(in), dimension(nspec,nreg,nlev) :: od
    real(jprb), intent(in), dimension(nspec,2:nreg,nlev) :: ssa
    real(jprb), intent(in), dimension(nspec,nlev) :: asymmetry

    ! Outputs

    real(jprb), intent(out), dimension(nspec,nreg,nlev) :: reflectance, transmittance
    real(jprb), intent(out), dimension(nspec,nreg,nlev) :: source_up, source_dn

    real(jprb) :: gamma1, gamma2, coeff, coeff_dn_base, coeff_dn_top, coeff_up_base, coeff_up_top
    real(jprb) :: factor, exponential, exponential2, k_exponent, reftrans_factor

    integer(jpim) :: jspec, jreg, jlev

    real(jprb) :: hook_handle

    if (lhook) call dr_hook('multiregion_two_stream:calc_reflectance_transmittance',0,hook_handle)

    do jlev = 1,nlev
      do jreg = 1,nreg
        if (jreg == 1 .or. is_cloud_free_layer(jlev)) then
          ! No-scattering solution: compute upward and downward
          ! emission assuming the Planck function to vary linearly
          ! with optical depth within the layer (e.g. Wiscombe , JQSRT
          ! 1976).
          do jspec = 1,nspec
            reflectance(jspec,jreg,jlev) = 0.0_jprb
            if (od(jspec,jreg,jlev) > 1.0e-3_jprb) then
              coeff = LW_DIFFUSIVITY*od(jspec,jreg,jlev)
              transmittance(jspec,jreg,jlev) = exp(-coeff)
              coeff = (planck_hl(jspec,jlev+1)-planck_hl(jspec,jlev)) / coeff
              coeff_up_top  =  coeff + planck_hl(jspec,jlev)
              coeff_up_base =  coeff + planck_hl(jspec,jlev+1)
              coeff_dn_top  = -coeff + planck_hl(jspec,jlev)
              coeff_dn_base = -coeff + planck_hl(jspec,jlev+1)
              source_up(jspec,jreg,jlev) = coeff_up_top &
                   &  - transmittance(jspec,jreg,jlev) * coeff_up_base
              source_dn(jspec,jreg,jlev) = coeff_dn_base &
                   &  - transmittance(jspec,jreg,jlev) * coeff_dn_top
            else
              ! Linear limit at low optical depth
              coeff = LW_DIFFUSIVITY*od(jspec,jreg,jlev)
              transmittance(jspec,jreg,jlev) = 1.0_jprb - coeff
              source_up(jspec,jreg,jlev) = coeff * 0.5_jprb &
                   &  * (planck_hl(jspec,jlev)+planck_hl(jspec,jlev+1))
              source_dn(jspec,jreg,jlev) = source_up(jspec,jreg,jlev)
            end if
          end do
        else
          ! Scattering solution
          factor = (LW_DIFFUSIVITY * 0.5_jprb) * ssa(jspec,jreg,jlev)
          gamma1 = LW_DIFFUSIVITY - factor*(1.0_jprb + asymmetry(jspec,jlev))
          gamma2 = factor * (1.0_jprb - asymmetry(jspec,jlev))

          do jspec = 1,nspec
            k_exponent = sqrt(max((gamma1 - gamma2) * (gamma1 + gamma2), &
                 1.E-12_jprb)) ! Eq 18 of Meador & Weaver (1980)
            if (od(jspec,jreg,jlev) > 1.0e-3_jprb) then
              exponential = exp(-k_exponent*od(jspec,jreg,jlev))
              exponential2 = exponential*exponential
              reftrans_factor = 1.0 / (k_exponent + gamma1&
                   &  + (k_exponent - gamma1)*exponential2)
              ! Meador & Weaver (1980) Eq. 25
              reflectance(jspec,jreg,jlev) = gamma2 &
                   &  * (1.0_jprb - exponential2) * reftrans_factor
              ! Meador & Weaver (1980) Eq. 26
              transmittance(jspec,jreg,jlev) = 2.0_jprb * k_exponent &
                   &  * exponential * reftrans_factor
      
              ! Compute upward and downward emission assuming the
              ! Planck function to vary linearly with optical depth
              ! within the layer (e.g. Wiscombe , JQSRT 1976).

              ! Stackhouse and Stephens (JAS 1991) Eqs 5 & 12
              coeff = (planck_hl(jspec,jlev+1)-planck_hl(jspec,jlev)) &
                   / (od(jspec,jreg,jlev)*(gamma1+gamma2))
              coeff_up_top  =  coeff + planck_hl(jspec,jlev)
              coeff_up_base =  coeff + planck_hl(jspec,jlev+1)
              coeff_dn_top  = -coeff + planck_hl(jspec,jlev)
              coeff_dn_base = -coeff + planck_hl(jspec,jlev+1)
              source_up(jspec,jreg,jlev) = coeff_up_top &
                   &  - reflectance(jspec,jreg,jlev) * coeff_dn_top &
                   &  - transmittance(jspec,jreg,jlev) * coeff_up_base
              source_dn(jspec,jreg,jlev) = coeff_dn_base &
                   &  - reflectance(jspec,jreg,jlev) * coeff_up_base &
                   &  - transmittance(jspec,jreg,jlev) * coeff_dn_top
            else
              reflectance(jspec,jreg,jlev) = gamma2 * od(jspec,jreg,jlev)
              transmittance(jspec,jreg,jlev) = (1.0_jprb - k_exponent*od(jspec,jreg,jlev)) &
                   &  / (1.0_jprb + od(jspec,jreg,jlev)*(gamma1-k_exponent))
              source_up(jspec,jreg,jlev) = (1.0_jprb - reflectance(jspec,jreg,jlev) &
                   &  - transmittance(jspec,jreg,jlev)) &
                   &       * 0.5 * (planck_hl(jspec,jlev) + planck_hl(jspec,jlev+1))
              source_dn(jspec,jreg,jlev) = source_up(jspec,jreg,jlev)
            end if
          end do
        end if
      end do
    end do

    if (lhook) call dr_hook('multiregion_two_stream:calc_reflectance_transmittance',1,hook_handle)

  end subroutine calc_reflectance_transmittance

end module tcrad_two_stream
