! radiation_two_stream.F90 - Compute two-stream coefficients
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
!
! Modifications
!   2017-05-04  P Dueben/R Hogan  Use JPRD where double precision essential
!   2017-07-12  R Hogan  Optimized LW coeffs in low optical depth case
!   2017-07-26  R Hogan  Added calc_frac_scattered_diffuse_sw routine
!   2017-10-23  R Hogan  Renamed single-character variables
!   2020-12-xx  P Ukkonen Performance optimizations:
!   1) Faster single precision computations by adjusting minimum k value instead of using JPRD
!   2) Faster reformulated version of calc_reflectance_transmittance_sw from RTE (sw_two_stream)

module radiation_two_stream

  use parkind1, only : jprb, jprd

  implicit none
  public

  ! Elsasser's factor: the effective factor by which the zenith
  ! optical depth needs to be multiplied to account for longwave
  ! transmission at all angles through the atmosphere.  Alternatively
  ! think of acos(1/lw_diffusivity) to be the effective zenith angle
  ! of longwave radiation.
  real(jprd), parameter :: LwDiffusivity = 1.66_jprd

  ! Shortwave diffusivity factor assumes hemispheric isotropy, assumed
  ! by Zdunkowski's scheme and most others; note that for efficiency
  ! this parameter is not used in the calculation of the gamma values,
  ! but is used in the SPARTACUS solver.
  real(jprb), parameter :: SwDiffusivity = 2.00_jprb

  ! Make minimum k value depend on precision, allowing to avoid JPRD altogether
#ifdef SINGLE_PRECISION
  real(jprb), parameter :: k_min = 1.e-4_jprb
#else
  real(jprb), parameter :: k_min = 1.e-12_jprb
#endif

  ! The routines in this module can be called millions of times, so
  !calling Dr Hook for each one may be a significant overhead.
  !Uncomment the following to turn Dr Hook on.
!#define DO_DR_HOOK_TWO_STREAM

contains

#ifdef FAST_EXPONENTIAL
  !---------------------------------------------------------------------
  ! Fast exponential for negative arguments: a Pade approximant that
  ! doesn't go negative for negative arguments, applied to arg/8, and
  ! the result is then squared three times
  elemental function exp_fast(arg) result(ex)
    real(jprb), intent(in)  :: arg
    real(jprb)              :: ex
    ex = 1.0_jprb / (1.0_jprb + arg*(-0.125_jprb &
        + arg*(0.0078125_jprb - 0.000325520833333333_jprb * arg)))
    ex = ex*ex
    ex = ex*ex
    ex = ex*ex
  end function exp_fast
#else
#define exp_fast exp
#endif

  !---------------------------------------------------------------------
  ! Calculate the two-stream coefficients gamma1 and gamma2 for the
  ! longwave
  subroutine calc_two_stream_gammas_lw(ng, ssa, g, &
       &                               gamma1, gamma2)

#ifdef DO_DR_HOOK_TWO_STREAM
    use yomhook, only : lhook, dr_hook
#endif

    integer, intent(in) :: ng
    ! Sngle scattering albedo and asymmetry factor:
    real(jprb), intent(in),  dimension(:) :: ssa, g
    real(jprb), intent(out), dimension(:) :: gamma1, gamma2

    real(jprb) :: factor

    integer    :: jg

#ifdef DO_DR_HOOK_TWO_STREAM
    real(jprb) :: hook_handle

    if (lhook) call dr_hook('radiation_two_stream:calc_two_stream_gammas_lw',0,hook_handle)
#endif

    do jg = 1, ng
      ! Fu et al. (1997), Eq 2.9 and 2.10:
      !      gamma1(jg) = LwDiffusivity * (1.0_jprb - 0.5_jprb*ssa(jg) &
      !           &                    * (1.0_jprb + g(jg)))
      !      gamma2(jg) = LwDiffusivity * 0.5_jprb * ssa(jg) &
      !           &                    * (1.0_jprb - g(jg))
      ! Reduce number of multiplications
      factor = (LwDiffusivity * 0.5_jprb) * ssa(jg)
      gamma1(jg) = LwDiffusivity - factor*(1.0_jprb + g(jg))
      gamma2(jg) = factor * (1.0_jprb - g(jg))
    end do

#ifdef DO_DR_HOOK_TWO_STREAM
    if (lhook) call dr_hook('radiation_two_stream:calc_two_stream_gammas_lw',1,hook_handle)
#endif

  end subroutine calc_two_stream_gammas_lw


  !---------------------------------------------------------------------
  ! Calculate the two-stream coefficients gamma1-gamma4 in the
  ! shortwave
  subroutine calc_two_stream_gammas_sw(ng, mu0, ssa, g, &
       &                               gamma1, gamma2, gamma3)

#ifdef DO_DR_HOOK_TWO_STREAM
    use yomhook, only : lhook, dr_hook
#endif

    integer, intent(in) :: ng
    ! Cosine of solar zenith angle, single scattering albedo and
    ! asymmetry factor:
    real(jprb), intent(in)                :: mu0
    real(jprb), intent(in),  dimension(:) :: ssa, g
    real(jprb), intent(out), dimension(:) :: gamma1, gamma2, gamma3

    real(jprb) :: factor

    integer    :: jg

#ifdef DO_DR_HOOK_TWO_STREAM
    real(jprb) :: hook_handle

    if (lhook) call dr_hook('radiation_two_stream:calc_two_stream_gammas_sw',0,hook_handle)
#endif

    ! Zdunkowski "PIFM" (Zdunkowski et al., 1980; Contributions to
    ! Atmospheric Physics 53, 147-66)
    do jg = 1, ng
      !      gamma1(jg) = 2.0_jprb  - ssa(jg) * (1.25_jprb + 0.75_jprb*g(jg))
      !      gamma2(jg) = 0.75_jprb *(ssa(jg) * (1.0_jprb - g(jg)))
      !      gamma3(jg) = 0.5_jprb  - (0.75_jprb*mu0)*g(jg)
      ! Optimized version:
      factor = 0.75_jprb*g(jg)
      gamma1(jg) = 2.0_jprb  - ssa(jg) * (1.25_jprb + factor)
      gamma2(jg) = ssa(jg) * (0.75_jprb - factor)
      gamma3(jg) = 0.5_jprb  - mu0*factor
    end do

#ifdef DO_DR_HOOK_TWO_STREAM
    if (lhook) call dr_hook('radiation_two_stream:calc_two_stream_gammas_sw',1,hook_handle)
#endif

  end subroutine calc_two_stream_gammas_sw


  !---------------------------------------------------------------------
  ! Compute the longwave reflectance and transmittance to diffuse
  ! radiation using the Meador & Weaver formulas, as well as the
  ! upward flux at the top and the downward flux at the base of the
  ! layer due to emission from within the layer assuming a linear
  ! variation of Planck function within the layer.
  subroutine calc_reflectance_transmittance_lw(ng, &
       &    od, gamma1, gamma2, planck_top, planck_bot, &
       &    reflectance, transmittance, source_up, source_dn)

#ifdef DO_DR_HOOK_TWO_STREAM
    use yomhook, only : lhook, dr_hook
#endif

    integer, intent(in) :: ng

    ! Optical depth and single scattering albedo
    real(jprb), intent(in), dimension(ng) :: od

    ! The two transfer coefficients from the two-stream
    ! differentiatial equations (computed by
    ! calc_two_stream_gammas_lw)
    real(jprb), intent(in), dimension(ng) :: gamma1, gamma2

    ! The Planck terms (functions of temperature) at the top and
    ! bottom of the layer
    real(jprb), intent(in), dimension(ng) :: planck_top, planck_bot

    ! The diffuse reflectance and transmittance, i.e. the fraction of
    ! diffuse radiation incident on a layer from either top or bottom
    ! that is reflected back or transmitted through
    real(jprb), intent(out), dimension(ng) :: reflectance, transmittance

    ! The upward emission at the top of the layer and the downward
    ! emission at its base, due to emission from within the layer
    real(jprb), intent(out), dimension(ng) :: source_up, source_dn

    real(jprb) :: k_exponent, reftrans_factor
    real(jprb) :: exponential  ! = exp(-k_exponent*od)
    real(jprb) :: exponential2 ! = exp(-2*k_exponent*od)

    real(jprb) :: coeff, coeff_up_top, coeff_up_bot, coeff_dn_top, coeff_dn_bot

    integer :: jg

#ifdef DO_DR_HOOK_TWO_STREAM
    real(jprb) :: hook_handle

    if (lhook) call dr_hook('radiation_two_stream:calc_reflectance_transmittance_lw',0,hook_handle)
#endif

    do jg = 1, ng
      if (od(jg) > 1.0e-3_jprb) then
        k_exponent = sqrt(max((gamma1(jg) - gamma2(jg)) * (gamma1(jg) + gamma2(jg)), &
          k_min)) ! Eq 18 of Meador & Weaver (1980)
        exponential = exp_fast(-k_exponent*od(jg))
        exponential2 = exponential*exponential
        reftrans_factor = 1.0 / (k_exponent + gamma1(jg) + (k_exponent - gamma1(jg))*exponential2)
        ! Meador & Weaver (1980) Eq. 25
        reflectance(jg) = gamma2(jg) * (1.0_jprb - exponential2) * reftrans_factor
        ! Meador & Weaver (1980) Eq. 26
        transmittance(jg) = 2.0_jprb * k_exponent * exponential * reftrans_factor
      
        ! Compute upward and downward emission assuming the Planck
        ! function to vary linearly with optical depth within the layer
        ! (e.g. Wiscombe , JQSRT 1976).

        ! Stackhouse and Stephens (JAS 1991) Eqs 5 & 12
        coeff = (planck_bot(jg)-planck_top(jg)) / (od(jg)*(gamma1(jg)+gamma2(jg)))
        coeff_up_top  =  coeff + planck_top(jg)
        coeff_up_bot  =  coeff + planck_bot(jg)
        coeff_dn_top  = -coeff + planck_top(jg)
        coeff_dn_bot  = -coeff + planck_bot(jg)
        source_up(jg) =  coeff_up_top - reflectance(jg) * coeff_dn_top - transmittance(jg) * coeff_up_bot
        source_dn(jg) =  coeff_dn_bot - reflectance(jg) * coeff_up_bot - transmittance(jg) * coeff_dn_top
      else
        k_exponent = sqrt(max((gamma1(jg) - gamma2(jg)) * (gamma1(jg) + gamma2(jg)), &
             k_min)) ! Eq 18 of Meador & Weaver (1980)
        reflectance(jg) = gamma2(jg) * od(jg)
        transmittance(jg) = (1.0_jprb - k_exponent*od(jg)) / (1.0_jprb + od(jg)*(gamma1(jg)-k_exponent))
        source_up(jg) = (1.0_jprb - reflectance(jg) - transmittance(jg)) &
             &       * 0.5 * (planck_top(jg) + planck_bot(jg))
        source_dn(jg) = source_up(jg)
      end if
    end do
    
#ifdef DO_DR_HOOK_TWO_STREAM
    if (lhook) call dr_hook('radiation_two_stream:calc_reflectance_transmittance_lw',1,hook_handle)
#endif
  
  end subroutine calc_reflectance_transmittance_lw
  


  !---------------------------------------------------------------------
  ! As calc_reflectance_transmittance_lw but for an isothermal layer
  subroutine calc_reflectance_transmittance_isothermal_lw(ng, &
       &    od, gamma1, gamma2, planck, &
       &    reflectance, transmittance, source)

#ifdef DO_DR_HOOK_TWO_STREAM
    use yomhook, only : lhook, dr_hook
#endif

    integer, intent(in) :: ng

    ! Optical depth and single scattering albedo
    real(jprb), intent(in), dimension(ng) :: od

    ! The two transfer coefficients from the two-stream
    ! differentiatial equations (computed by
    ! calc_two_stream_gammas_lw)
    real(jprb), intent(in), dimension(ng) :: gamma1, gamma2

    ! The Planck terms (functions of temperature) constant through the
    ! layer
    real(jprb), intent(in), dimension(ng) :: planck

    ! The diffuse reflectance and transmittance, i.e. the fraction of
    ! diffuse radiation incident on a layer from either top or bottom
    ! that is reflected back or transmitted through
    real(jprb), intent(out), dimension(ng) :: reflectance, transmittance

    ! The upward emission at the top of the layer and the downward
    ! emission at its base, due to emission from within the layer
    real(jprb), intent(out), dimension(ng) :: source

    real(jprb) :: k_exponent, reftrans_factor
    real(jprb) :: exponential  ! = exp(-k_exponent*od)
    real(jprb) :: exponential2 ! = exp(-2*k_exponent*od)

    integer :: jg

#ifdef DO_DR_HOOK_TWO_STREAM
    real(jprb) :: hook_handle

    if (lhook) call dr_hook('radiation_two_stream:calc_reflectance_transmittance_isothermal_lw',0,hook_handle)
#endif

    do jg = 1, ng
      k_exponent = sqrt(max((gamma1(jg) - gamma2(jg)) * (gamma1(jg) + gamma2(jg)), &
           k_min)) ! Eq 18 of Meador & Weaver (1980)
      exponential = exp_fast(-k_exponent*od(jg))
      exponential2 = exponential*exponential
      reftrans_factor = 1.0 / (k_exponent + gamma1(jg) + (k_exponent - gamma1(jg))*exponential2)
      ! Meador & Weaver (1980) Eq. 25
      reflectance(jg) = gamma2(jg) * (1.0_jprb - exponential2) * reftrans_factor
      ! Meador & Weaver (1980) Eq. 26
      transmittance(jg) = 2.0_jprb * k_exponent * exponential * reftrans_factor
      
      ! Emissivity of layer is one minus reflectance minus
      ! transmittance, multiply by Planck function to get emitted
      ! ousrce
      source(jg) = planck(jg) * (1.0_jprb - reflectance(jg) - transmittance(jg))
    end do
    
#ifdef DO_DR_HOOK_TWO_STREAM
    if (lhook) call dr_hook('radiation_two_stream:calc_reflectance_transmittance_isothermal_lw',1,hook_handle)
#endif
  
  end subroutine calc_reflectance_transmittance_isothermal_lw
  


  !---------------------------------------------------------------------
  ! Compute the longwave transmittance to diffuse radiation in the
  ! no-scattering case, as well as the upward flux at the top and the
  ! downward flux at the base of the layer due to emission from within
  ! the layer assuming a linear variation of Planck function within
  ! the layer.
  subroutine calc_no_scattering_transmittance_lw(ng, &
       &    od, planck_top, planck_bot, transmittance, source_up, source_dn)

#ifdef DO_DR_HOOK_TWO_STREAM
    use yomhook, only : lhook, dr_hook
#endif

    integer, intent(in) :: ng

    ! Optical depth and single scattering albedo
    real(jprb), intent(in), dimension(ng) :: od

    ! The Planck terms (functions of temperature) at the top and
    ! bottom of the layer
    real(jprb), intent(in), dimension(ng) :: planck_top, planck_bot

    ! The diffuse transmittance, i.e. the fraction of diffuse
    ! radiation incident on a layer from either top or bottom that is
    ! reflected back or transmitted through
    real(jprb), intent(out), dimension(ng) :: transmittance

    ! The upward emission at the top of the layer and the downward
    ! emission at its base, due to emission from within the layer
    real(jprb), intent(out), dimension(ng) :: source_up, source_dn

    real(jprb) :: coeff, coeff_up_top, coeff_up_bot, coeff_dn_top, coeff_dn_bot !, planck_mean

    integer :: jg

#ifdef DO_DR_HOOK_TWO_STREAM
    real(jprb) :: hook_handle

    if (lhook) call dr_hook('radiation_two_stream:calc_no_scattering_transmittance_lw',0,hook_handle)
#endif

    do jg = 1, ng
      ! Compute upward and downward emission assuming the Planck
      ! function to vary linearly with optical depth within the layer
      ! (e.g. Wiscombe , JQSRT 1976).
      if (od(jg) > 1.0e-3) then
        ! Simplified from calc_reflectance_transmittance_lw above
        coeff = LwDiffusivity*od(jg)
        transmittance(jg) = exp_fast(-coeff)
        coeff = (planck_bot(jg)-planck_top(jg)) / coeff
        coeff_up_top  =  coeff + planck_top(jg)
        coeff_up_bot  =  coeff + planck_bot(jg)
        coeff_dn_top  = -coeff + planck_top(jg)
        coeff_dn_bot  = -coeff + planck_bot(jg)
        source_up(jg) =  coeff_up_top - transmittance(jg) * coeff_up_bot
        source_dn(jg) =  coeff_dn_bot - transmittance(jg) * coeff_dn_top
      else
        ! Linear limit at low optical depth
        coeff = LwDiffusivity*od(jg)
        transmittance(jg) = 1.0_jprb - coeff
        source_up(jg) = coeff * 0.5_jprb * (planck_top(jg)+planck_bot(jg))
        source_dn(jg) = source_up(jg)
      end if
    end do

    ! Method in the older IFS radiation scheme
    !    do j = 1, n
    !      coeff = od(jg) / (3.59712_jprb + od(jg))
    !      planck_mean = 0.5_jprb * (planck_top(jg) + planck_bot(jg))
    !      
    !      source_up(jg) = (1.0_jprb-transmittance(jg)) * (planck_mean + (planck_top(jg)    - planck_mean) * coeff)
    !      source_dn(jg) = (1.0_jprb-transmittance(jg)) * (planck_mean + (planck_bot(jg) - planck_mean) * coeff)
    !    end do

#ifdef DO_DR_HOOK_TWO_STREAM
    if (lhook) call dr_hook('radiation_two_stream:calc_no_scattering_transmittance_lw',1,hook_handle)
#endif

  end subroutine calc_no_scattering_transmittance_lw
   
   
  !---------------------------------------------------------------------
  ! Compute the shortwave reflectance and transmittance to diffuse
  ! radiation using the Meador & Weaver formulas, as well as the
  ! "direct" reflection and transmission, which really means the rate
  ! of transfer of direct solar radiation (into a plane perpendicular
  ! to the direct beam) into diffuse upward and downward streams at
  ! the top and bottom of the layer, respectively.  Finally,
  ! trans_dir_dir is the transmittance of the atmosphere to direct
  ! radiation with no scattering.
  subroutine calc_reflectance_transmittance_sw(ng, mu0, od, ssa, &
       &      gamma1, gamma2, gamma3, ref_diff, trans_diff, &
       &      ref_dir, trans_dir_diff, trans_dir_dir)
    
#ifdef DO_DR_HOOK_TWO_STREAM
    use yomhook, only : lhook, dr_hook
#endif

    integer, intent(in) :: ng

    ! Cosine of solar zenith angle
    real(jprb), intent(in) :: mu0

    ! Optical depth and single scattering albedo
    real(jprb), intent(in), dimension(ng) :: od, ssa

    ! The three transfer coefficients from the two-stream
    ! differentiatial equations (computed by calc_two_stream_gammas)
    real(jprb), intent(in), dimension(ng) :: gamma1, gamma2, gamma3

    ! The direct reflectance and transmittance, i.e. the fraction of
    ! incoming direct solar radiation incident at the top of a layer
    ! that is either reflected back (ref_dir) or scattered but
    ! transmitted through the layer to the base (trans_dir_diff)
    real(jprb), intent(out), dimension(ng) :: ref_dir, trans_dir_diff

    ! The diffuse reflectance and transmittance, i.e. the fraction of
    ! diffuse radiation incident on a layer from either top or bottom
    ! that is reflected back or transmitted through
    real(jprb), intent(out), dimension(ng) :: ref_diff, trans_diff

    ! Transmittance of the direct been with no scattering
    real(jprb), intent(out), dimension(ng) :: trans_dir_dir

    real(jprb) :: gamma4, alpha1, alpha2, k_exponent, reftrans_factor
    real(jprb) :: exponential0 ! = exp(-od/mu0)
    real(jprb) :: exponential  ! = exp(-k_exponent*od)
    real(jprb) :: exponential2 ! = exp(-2*k_exponent*od)
    real(jprb) :: k_mu0, k_gamma3, k_gamma4
    real(jprb) :: k_2_exponential, od_over_mu0

    integer :: jg

#ifdef DO_DR_HOOK_TWO_STREAM
    real(jprb) :: hook_handle

    if (lhook) call dr_hook('radiation_two_stream:calc_reflectance_transmittance_sw',0,hook_handle)
#endif

    do jg = 1, ng
      od_over_mu0 = max(od(jg) / mu0, 0.0_jprb)
      ! In the IFS this appears to be faster without testing the value
      ! of od_over_mu0:
      if (.true.) then
!      if (od_over_mu0 > 1.0e-6_jprb) then
        gamma4 = 1.0_jprb - gamma3(jg)
        alpha1 = gamma1(jg)*gamma4     + gamma2(jg)*gamma3(jg) ! Eq. 16
        alpha2 = gamma1(jg)*gamma3(jg) + gamma2(jg)*gamma4    ! Eq. 17
        
        ! Note that if the minimum value is reduced (e.g. to 1.0e-24)
        ! then noise starts to appear as a function of solar zenith
        ! angle
        k_exponent = sqrt(max((gamma1(jg) - gamma2(jg)) * (gamma1(jg) + gamma2(jg)), &
             &       k_min)) ! Eq 18
        k_mu0 = k_exponent*mu0
        k_gamma3 = k_exponent*gamma3(jg)
        k_gamma4 = k_exponent*gamma4
        ! Check for mu0 <= 0!
        exponential0 = exp_fast(-od_over_mu0)
        trans_dir_dir(jg) = exponential0
        exponential = exp_fast(-k_exponent*od(jg))
        
        exponential2 = exponential*exponential
        k_2_exponential = 2.0_jprb * k_exponent * exponential
        
        if (k_mu0 == 1.0_jprb) then
          k_mu0 = 1.0_jprb - 10.0_jprb*epsilon(1.0_jprb)
        end if
        
        reftrans_factor = 1.0_jprb / (k_exponent + gamma1(jg) + (k_exponent - gamma1(jg))*exponential2)
        
        ! Meador & Weaver (1980) Eq. 25
        ref_diff(jg) = gamma2(jg) * (1.0_jprb - exponential2) * reftrans_factor
        
        ! Meador & Weaver (1980) Eq. 26
        trans_diff(jg) = k_2_exponential * reftrans_factor
        
        ! Here we need mu0 even though it wasn't in Meador and Weaver
        ! because we are assuming the incoming direct flux is defined
        ! to be the flux into a plane perpendicular to the direction of
        ! the sun, not into a horizontal plane
        reftrans_factor = mu0 * ssa(jg) * reftrans_factor / (1.0_jprb - k_mu0*k_mu0)
        
        ! Meador & Weaver (1980) Eq. 14, multiplying top & bottom by
        ! exp(-k_exponent*od) in case of very high optical depths
        ref_dir(jg) = reftrans_factor &
             &  * ( (1.0_jprb - k_mu0) * (alpha2 + k_gamma3) &
             &     -(1.0_jprb + k_mu0) * (alpha2 - k_gamma3)*exponential2 &
             &     -k_2_exponential*(gamma3(jg) - alpha2*mu0)*exponential0)
        
        ! Meador & Weaver (1980) Eq. 15, multiplying top & bottom by
        ! exp(-k_exponent*od), minus the 1*exp(-od/mu0) term representing direct
        ! unscattered transmittance.  
        trans_dir_diff(jg) = reftrans_factor * ( k_2_exponential*(gamma4 + alpha1*mu0) &
            & - exponential0 &
            & * ( (1.0_jprb + k_mu0) * (alpha1 + k_gamma4) &
            &    -(1.0_jprb - k_mu0) * (alpha1 - k_gamma4) * exponential2) )

      else
        ! Low optical-depth limit; see equations 19, 20 and 27 from
        ! Meador & Weaver (1980)
        trans_diff(jg)     = 1.0_jprb - gamma1(jg) * od(jg)
        ref_diff(jg)       = gamma2(jg) * od(jg)
        trans_dir_diff(jg) = (1.0_jprb - gamma3(jg)) * ssa(jg) * od(jg)
        ref_dir(jg)        = gamma3(jg) * ssa(jg) * od(jg)
        trans_dir_dir(jg)  = 1.0_jprb - od_over_mu0
      end if
    end do
    
#ifdef DO_DR_HOOK_TWO_STREAM
    if (lhook) call dr_hook('radiation_two_stream:calc_reflectance_transmittance_sw',1,hook_handle)
#endif
 
  end subroutine calc_reflectance_transmittance_sw
  
  ! -------------------------------------------------------------------------------------------------
  !
  !   Two-stream procedure "sw_two_stream" adapted from RTE (Radiative Transfer for Energetics) 
  !   code, developed by Robert Pincus, and including some optimization by Peter Ukkonen 
  !   for RTE+RRTMGP-NN: https://github.com/peterukk/rte-rrtmgp-nn
  !   Functionally identical to calc_reflectance_transmittance_sw
  !
  ! -------------------------------------------------------------------------------------------------
  !
  ! Two-stream solutions to direct and diffuse reflectance and transmittance for a layer
  !    with optical depth od, single scattering albedo ssa, and asymmetery parameter g.
  !
  ! Equations are developed in Meador and Weaver, 1980,
  !    doi:10.1175/1520-0469(1980)037<0630:TSATRT>2.0.CO;2
  !
  pure subroutine calc_reflectance_transmittance_sw_opt(ngpt, mu0, od, ssa, gamma1 ,gamma2, gamma3, &
                                ref_dif, trans_dif, ref_dir, trans_dir_dif, trans_dir_dir)
    integer,                          intent(in)  :: ngpt
    real(jprb),                       intent(in)  :: mu0
    real(jprb), dimension(ngpt),      intent(in)  :: od, ssa
    real(jprb), dimension(ngpt),      intent(in)  :: gamma1, gamma2, gamma3
    real(jprb), dimension(ngpt),      intent(out) :: ref_dif, trans_dif, ref_dir, trans_dir_dif, trans_dir_dir
    ! -----------------------
    integer  :: i, j

    ! Variables used in Meador and Weaver
    real(jprb), dimension(ngpt) :: gamma4, alpha1, alpha2, k
    ! Ancillary variables
    real(jprb), dimension(ngpt) :: exp_minuskod, exp_minus2kod, RT_term
    real(jprb) :: k_gamma3, k_gamma4  ! Need to be in double precision?
    real(jprb) :: k_mu, k_mu2, mu0_inv

    mu0_inv = 1._jprb/mu0

    do i = 1, ngpt
      ! Zdunkowski Practical Improved Flux Method "PIFM"
      !  (Zdunkowski et al., 1980;  Contributions to Atmospheric Physics 53, 147-66)
      !
      ! gamma1(i)= (8._jprb - ssa(i) * (5._jprb + 3._jprb * g(i))) * .25_jprb
      ! gamma2(i)=  3._jprb *(ssa(i) * (1._jprb -         g(i))) * .25_jprb
      ! gamma3(i)= (2._jprb - 3._jprb * mu0 *              g(i) ) * .25_jprb
      gamma4(i)=  1._jprb - gamma3(i)

      alpha1(i) = gamma1(i) * gamma4(i) + gamma2(i) * gamma3(i)           ! Eq. 16
      alpha2(i) = gamma1(i) * gamma3(i) + gamma2(i) * gamma4(i)           ! Eq. 17

      k(i) = sqrt(max((gamma1(i) - gamma2(i)) * (gamma1(i) + gamma2(i)),  k_min))
    end do

    exp_minuskod(:) = exp_fast(-od(:)*k(:))
    !
    ! Diffuse reflection and transmission
    !
    do i = 1, ngpt
      exp_minus2kod(i)  = exp_minuskod(i) * exp_minuskod(i)

      ! Refactored to avoid rounding errors when k, gamma1 are of very different magnitudes
      RT_term(i) = 1._jprb / (k(i) * (1._jprb + exp_minus2kod(i)) + gamma1(i) * (1._jprb - exp_minus2kod(i)) )

      ! Equation 25
      ref_dif(i) = RT_term(i) * gamma2(i) * (1._jprb - exp_minus2kod(i))

      ! Equation 26
      trans_dif(i) = RT_term(i) * 2._jprb * k(i) * exp_minuskod(i)
      
    end do

    !
    ! Transmittance of direct, unscattered beam. Also used below
    !
    trans_dir_dir(:) = exp_fast(-od(:)*mu0_inv)
    !
    ! Direct reflect and transmission
    !
    do i = 1, ngpt
      k_mu     = k(i) * mu0
      k_mu2    = k_mu*k_mu
      k_gamma3 = k(i) * gamma3(i)
      k_gamma4 = k(i) * gamma4(i)
      !
      ! Equation 14, multiplying top and bottom by exp_fast(-k*od)
      !   and rearranging to avoid div by 0.      
      ! Here we need mu0 even though it wasn't in Meador and Weaver
      ! because we are assuming the incoming direct flux is defined
      ! to be the flux into a plane perpendicular to the direction of
      ! the sun, not into a horizontal plane   
      RT_term(i) =  mu0 * ssa(i) *  &
      RT_term(i) / merge(1._jprb - k_mu2, epsilon(1._jprb), abs(1._jprb - k_mu2) >= epsilon(1._jprb))
      !  --> divide by (1 - kmu2) when (1-kmu2)> eps, otherwise divide by eps

      ref_dir(i) = RT_term(i)  *                              &
              (   (1._jprb - k_mu) * (alpha2(i) + k_gamma3) -  &
                  (1._jprb + k_mu) * (alpha2(i) - k_gamma3) * exp_minus2kod(i) - &
            2.0_jprb * (k_gamma3 - alpha2(i) * k_mu)  * exp_minuskod (i) * trans_dir_dir(i)  )
      !
      ! Equation 15, multiplying top and bottom by exp(-k*od),
      !   multiplying through by exp(-od/mu0) to
      !   prefer underflow to overflow
      ! Omitting direct transmittance
      !
      trans_dir_dif(i) = -RT_term(i) *                                                                 &
                  ((1._jprb + k_mu) * (alpha1(i) + k_gamma4)                     * trans_dir_dir(i) - &
                    (1._jprb - k_mu) * (alpha1(i) - k_gamma4) * exp_minus2kod(i) * trans_dir_dir(i) - &
                    2.0_jprb * (k_gamma4 + alpha1(i) * k_mu)  * exp_minuskod (i))
                    
    end do
  
   end subroutine calc_reflectance_transmittance_sw_opt
  !---------------------------------------------------------------------
  ! As above but with height as a vertical coordinate rather than
  ! optical depth
  subroutine calc_reflectance_transmittance_z_sw(ng, mu0, depth, &
       &      gamma0, gamma1, gamma2, gamma3, gamma4, &
       &      ref_diff, trans_diff, ref_dir, trans_dir_diff, trans_dir_dir)
    
#ifdef DO_DR_HOOK_TWO_STREAM
    use yomhook, only : lhook, dr_hook
#endif

    integer, intent(in) :: ng

    ! Cosine of solar zenith angle
    real(jprb), intent(in) :: mu0

    ! Layer depth
    real(jprb), intent(in) :: depth

    ! The four transfer coefficients from the two-stream
    ! differentiatial equations
    real(jprb), intent(in), dimension(ng) :: gamma1, gamma2, gamma3, gamma4

    ! An additional coefficient for direct unscattered flux "Fdir"
    ! such that dFdir/dz = -gamma0*Fdir
    real(jprb), intent(in), dimension(ng) :: gamma0

    ! The direct reflectance and transmittance, i.e. the fraction of
    ! incoming direct solar radiation incident at the top of a layer
    ! that is either reflected back (ref_dir) or scattered but
    ! transmitted through the layer to the base (trans_dir_diff)
    real(jprb), intent(out), dimension(ng) :: ref_dir, trans_dir_diff

    ! The diffuse reflectance and transmittance, i.e. the fraction of
    ! diffuse radiation incident on a layer from either top or bottom
    ! that is reflected back or transmitted through
    real(jprb), intent(out), dimension(ng) :: ref_diff, trans_diff

    ! Transmittance of the direct been with no scattering
    real(jprb), intent(out), dimension(ng) :: trans_dir_dir

    real(jprb) :: alpha1, alpha2, k_exponent, reftrans_factor
    real(jprb) :: exponential0 ! = exp(-od/mu0)
    real(jprb) :: exponential  ! = exp(-k_exponent*od)
    real(jprb) :: exponential2 ! = exp(-2*k_exponent*od)
    real(jprb) :: k_mu0, k_gamma3, k_gamma4
    real(jprb) :: k_2_exponential, od_over_mu0

    integer :: jg

#ifdef DO_DR_HOOK_TWO_STREAM
    real(jprb) :: hook_handle

    if (lhook) call dr_hook('radiation_two_stream:calc_reflectance_transmittance_z_sw',0,hook_handle)
#endif

    do jg = 1, ng
      od_over_mu0 = max(gamma0(jg) * depth, 0.0_jprb)
      ! In the IFS this appears to be faster without testing the value
      ! of od_over_mu0:
      if (.true.) then
!      if (od_over_mu0 > 1.0e-6_jprb) then
        alpha1 = gamma1(jg)*gamma4(jg) + gamma2(jg)*gamma3(jg) ! Eq. 16
        alpha2 = gamma1(jg)*gamma3(jg) + gamma2(jg)*gamma4(jg) ! Eq. 17
        
        ! Note that if the minimum value is reduced (e.g. to 1.0e-24)
        ! then noise starts to appear as a function of solar zenith
        ! angle
        k_exponent = sqrt(max((gamma1(jg) - gamma2(jg)) * (gamma1(jg) + gamma2(jg)), &
             &      k_min)) ! Eq 18
        k_mu0 = k_exponent*mu0
        k_gamma3 = k_exponent*gamma3(jg)
        k_gamma4 = k_exponent*gamma4(jg)
        ! Check for mu0 <= 0!
        exponential0 = exp_fast(-od_over_mu0)
        trans_dir_dir(jg) = exponential0
        exponential = exp_fast(-k_exponent*depth)
        
        exponential2 = exponential*exponential
        k_2_exponential = 2.0_jprb * k_exponent * exponential
        
        if (k_mu0 == 1.0_jprb) then
          k_mu0 = 1.0_jprb - 10.0_jprb*epsilon(1.0_jprb)
        end if
        
        reftrans_factor = 1.0_jprb / (k_exponent + gamma1(jg) + (k_exponent - gamma1(jg))*exponential2)
        
        ! Meador & Weaver (1980) Eq. 25
        ref_diff(jg) = gamma2(jg) * (1.0_jprb - exponential2) * reftrans_factor
        
        ! Meador & Weaver (1980) Eq. 26
        trans_diff(jg) = k_2_exponential * reftrans_factor
        
        ! Here we need mu0 even though it wasn't in Meador and Weaver
        ! because we are assuming the incoming direct flux is defined
        ! to be the flux into a plane perpendicular to the direction of
        ! the sun, not into a horizontal plane
        reftrans_factor = mu0 * reftrans_factor / (1.0_jprb - k_mu0*k_mu0)
        
        ! Meador & Weaver (1980) Eq. 14, multiplying top & bottom by
        ! exp(-k_exponent*od) in case of very high optical depths
        ref_dir(jg) = reftrans_factor &
             &  * ( (1.0_jprb - k_mu0) * (alpha2 + k_gamma3) &
             &     -(1.0_jprb + k_mu0) * (alpha2 - k_gamma3)*exponential2 &
             &     -k_2_exponential*(gamma3(jg) - alpha2*mu0)*exponential0)
        
        ! Meador & Weaver (1980) Eq. 15, multiplying top & bottom by
        ! exp(-k_exponent*od), minus the 1*exp(-od/mu0) term representing direct
        ! unscattered transmittance.  
        trans_dir_diff(jg) = reftrans_factor * ( k_2_exponential*(gamma4(jg) + alpha1*mu0) &
            & - exponential0 &
            & * ( (1.0_jprb + k_mu0) * (alpha1 + k_gamma4) &
            &    -(1.0_jprb - k_mu0) * (alpha1 - k_gamma4) * exponential2) )

      else
        ! Low optical-depth limit; see equations 19, 20 and 27 from
        ! Meador & Weaver (1980)
        trans_diff(jg)     = 1.0_jprb - gamma1(jg) * depth
        ref_diff(jg)       = gamma2(jg) * depth
        trans_dir_diff(jg) = (1.0_jprb - gamma3(jg)) * depth
        ref_dir(jg)        = gamma3(jg) * depth
        trans_dir_dir(jg)  = 1.0_jprb - od_over_mu0
      end if
    end do
    
#ifdef DO_DR_HOOK_TWO_STREAM
    if (lhook) call dr_hook('radiation_two_stream:calc_reflectance_transmittance_z_sw',1,hook_handle)
#endif
 
  end subroutine calc_reflectance_transmittance_z_sw
  

  !---------------------------------------------------------------------
  ! Compute the fraction of shortwave transmitted diffuse radiation
  ! that is scattered during its transmission, used to compute
  ! entrapment in SPARTACUS
  subroutine calc_frac_scattered_diffuse_sw(ng, od, &
       &      gamma1, gamma2, frac_scat_diffuse)
    
#ifdef DO_DR_HOOK_TWO_STREAM
    use yomhook, only : lhook, dr_hook
#endif

    integer, intent(in) :: ng

    ! Optical depth
    real(jprb), intent(in), dimension(ng) :: od

    ! The first two transfer coefficients from the two-stream
    ! differentiatial equations (computed by calc_two_stream_gammas)
    real(jprb), intent(in), dimension(ng) :: gamma1, gamma2

    ! The fraction of shortwave transmitted diffuse radiation that is
    ! scattered during its transmission
    real(jprb), intent(out), dimension(ng) :: frac_scat_diffuse

    real(jprb) :: k_exponent, reftrans_factor
    real(jprb) :: exponential  ! = exp(-k_exponent*od)
    real(jprb) :: exponential2 ! = exp(-2*k_exponent*od)
    real(jprb) :: k_2_exponential
    integer :: jg

#ifdef DO_DR_HOOK_TWO_STREAM
    real(jprb) :: hook_handle

    if (lhook) call dr_hook('radiation_two_stream:calc_frac_scattered_diffuse_sw',0,hook_handle)
#endif

    do jg = 1, ng
      ! Note that if the minimum value is reduced (e.g. to 1.0e-24)
      ! then noise starts to appear as a function of solar zenith
      ! angle
      k_exponent = sqrt(max((gamma1(jg) - gamma2(jg)) * (gamma1(jg) + gamma2(jg)), &
           &       k_min)) ! Eq 18
      exponential = exp_fast(-k_exponent*od(jg))
      exponential2 = exponential*exponential
      k_2_exponential = 2.0_jprb * k_exponent * exponential
        
      reftrans_factor = 1.0_jprb / (k_exponent + gamma1(jg) + (k_exponent - gamma1(jg))*exponential2)
        
      ! Meador & Weaver (1980) Eq. 26.
      ! Until 1.1.8, used LwDiffusivity instead of 2.0, although the
      ! effect is very small
      !      frac_scat_diffuse(jg) = 1.0_jprb - min(1.0_jprb,exp_fast(-LwDiffusivity*od(jg)) &
      !           &  / max(1.0e-8_jprb, k_2_exponential * reftrans_factor))
      frac_scat_diffuse(jg) = 1.0_jprb &
           &  - min(1.0_jprb,exp_fast(-2.0_jprb*od(jg)) &
           &  / max(1.0e-8_jprb, k_2_exponential * reftrans_factor))
    end do
    
#ifdef DO_DR_HOOK_TWO_STREAM
    if (lhook) call dr_hook('radiation_two_stream:calc_frac_scattered_diffuse_sw',1,hook_handle)
#endif
 
  end subroutine calc_frac_scattered_diffuse_sw

end module radiation_two_stream
