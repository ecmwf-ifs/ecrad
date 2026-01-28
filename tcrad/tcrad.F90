! tcrad.F90 - Tripleclouds radiance model
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

module tcrad

  use parkind1, only : jpim, jprb

  implicit none
  public

  ! Three regions at each height, one clear and two cloudy
  integer(jpim), parameter :: NREGION = 3
  
  ! Two-stream scheme can be:
  !   Elsasser: radiance field assumed to be two pencil beams
  !     travelling at cosine-zenith angles (mu) of +/-1/1.66;
  !   Eddington: radiance field assumed to be L(mu)=L0+mu*L1.
  enum, bind(c)
    enumerator ITwoStreamElsasser, ITwoStreamEddington, ITwoStreamLegendre, &
         &  ITwoStreamHybrid, ITwoStreamScaledWiscombeGrams
  end enum

  ! Two stream scheme currently in use 
  integer(jpim) :: i_two_stream_scheme = ITwoStreamElsasser
  
  ! The effective factor by which the zenith optical depth needs to be
  ! multiplied to account for longwave transmission at all angles
  ! through the atmosphere.  Alternatively think of
  ! acos(1/lw_diffusivity) to be the effective zenith angle of
  ! longwave radiation.
  real(jprb) :: lw_diffusivity = 1.66_jprb ! Elsasser default

  ! Optionally have a different diffusivity in cloud, where radiation
  ! tends to be more isotropic
  real(jprb), private :: lw_diffusivity_cloud = 1.66_jprb
  
  ! To avoid division by near-zero values use simpler formulae in the
  ! low optical depth regime
#ifdef PARKIND1_SINGLE
  real(jprb), parameter :: OD_THRESH_2STREAM = 1.0e-4_jprb
  real(jprb), parameter :: OD_THRESH         = 1.0e-4_jprb
  real(jprb), parameter :: MIN_K_SQUARED     = 1.0e-6_jprb
#else
  real(jprb), parameter :: OD_THRESH_2STREAM = 1.0e-7_jprb
  real(jprb), parameter :: OD_THRESH         = 1.0e-7_jprb
  real(jprb), parameter :: MIN_K_SQUARED     = 1.0e-12_jprb
#endif

contains

!---------------------------------------------------------------------
! Set the two-stream scheme; this overwrites global module variables
! so should normally be called from outside a parallel block
subroutine set_two_stream_scheme(i_scheme)
    
  integer(jpim), intent(in) :: i_scheme

  ! Only overwrite global module variables if they need changing
  if (i_scheme == ITwoStreamEddington) then
    ! Toon et al. (1989), Table 1
    i_two_stream_scheme  = ITwoStreamEddington
    lw_diffusivity       = 2.0_jprb
    lw_diffusivity_cloud = 2.0_jprb
  else if (i_scheme == ITwoStreamElsasser) then
    ! Elsasser (1942): "Quadrature" from Toon et al.'s Table 1 but
      ! with alternative diffusivity
    i_two_stream_scheme  = ITwoStreamElsasser
    lw_diffusivity       = 1.66_jprb
    lw_diffusivity_cloud = 1.66_jprb
  else if (i_scheme == ITwoStreamLegendre) then
    ! Hemispheric mean from Toon et al.'s Table 1
    i_two_stream_scheme  = ITwoStreamLegendre
    lw_diffusivity       = 2.0_jprb
    lw_diffusivity_cloud = 2.0_jprb
  else if (i_scheme == ITwoStreamHybrid) then
    ! Hybrid between Elsasser (which is best for propagation in
    ! clear skies) and Legendre (which is arguably better for
    ! scattering in cloudy skies)
    i_two_stream_scheme  = ITwoStreamHybrid
    lw_diffusivity       = 1.66_jprb
    lw_diffusivity_cloud = 2.0_jprb
  else if (i_scheme == ITwoStreamScaledWiscombeGrams) then
    i_two_stream_scheme  = ITwoStreamScaledWiscombeGrams
    lw_diffusivity       = 1.66_jprb
    lw_diffusivity_cloud = 2.0_jprb
  end if
    
end subroutine set_two_stream_scheme


!---------------------------------------------------------------------
! Compute the longwave reflectance and transmittance to diffuse
! radiation using the Meador & Weaver (1980) two-stream formulas, as
! well as the upward flux at the top and the downward flux at the
! base of the layer due to emission from within the layer assuming a
! linear variation of Planck function within the layer.
subroutine calc_reflectance_transmittance(nspec, nlev, nreg, &
     &  region_fracs, planck_hl, od, ssa, asymmetry, &
     &  reflectance, transmittance, source_up, source_dn)
        
  use yomhook,  only           : lhook, dr_hook, jphook

  ! Inputs

  ! Number of spectral intervals, levels and regions
  integer(jpim), intent(in) :: nspec, nlev, nreg

  ! Fraction of the gridbox occupied by each region (summing to 1)
  ! at each level
  real(jprb), intent(in), dimension(nreg,nlev) :: region_fracs

  ! Planck function integrated over each spectral interval at each
  ! half-level, in W m-2 (i.e. the flux emitted by a horizontal
  ! black-body surface)
  real(jprb), intent(in), dimension(nspec,nlev+1) :: planck_hl

  ! Optical depth in each region and layer
  real(jprb), intent(in), dimension(nspec,nreg,nlev) :: od
  
  ! Single scattering albedo in each cloudy region
  real(jprb), intent(in), dimension(nspec,2:nreg,nlev) :: ssa

  ! Asymmetry factor of clouds
  real(jprb), intent(in), dimension(nspec,nlev) :: asymmetry

  ! Outputs

  ! Layer reflectance and transmittance
  real(jprb), intent(out), dimension(nspec,nreg,nlev) :: reflectance, transmittance

  ! The upward emission at the top of the layer and the downward
  ! emission at its base, due to emission from within the layer, in
  ! units of Watts of power per square metre of the entire gridbox,
  ! so emission is proportional to the size of each region
  real(jprb), intent(out), dimension(nspec,nreg,nlev) :: source_up, source_dn

  ! Two-stream exchange coefficients
  real(jprb) :: gamma1, gamma2

  ! Working variables
  real(jprb) :: coeff, coeff_dn_base, coeff_dn_top, coeff_up_base, coeff_up_top
  real(jprb) :: factor, exponential, exponential2, k_exponent, reftrans_factor

  ! Loop indices
  integer(jpim) :: jspec, jreg, jlev

  real(jphook) :: hook_handle

  if (lhook) call dr_hook('tcrad:calc_reflectance_transmittance',0,hook_handle)

  ! Set cloudy regions to default values
  reflectance(:,2:,:)   = 0.0_jprb
  transmittance(:,2:,:) = 1.0_jprb
  source_up(:,2:,:)     = 0.0_jprb
  source_dn(:,2:,:)     = 0.0_jprb

  do jlev = 1,nlev
    ! No-scattering solution in clear-sky region: compute upward and
    ! downward emission assuming the Planck function to vary
    ! linearly with optical depth within the layer (e.g. Wiscombe,
    ! JQSRT 1976).
    do jspec = 1,nspec
      reflectance(jspec,1,jlev) = 0.0_jprb
      if (od(jspec,1,jlev) > OD_THRESH_2STREAM) then
        coeff = lw_diffusivity*od(jspec,1,jlev)
        transmittance(jspec,1,jlev) = exp(-coeff)
        coeff = (planck_hl(jspec,jlev+1)-planck_hl(jspec,jlev)) / coeff
        coeff_up_top  =  coeff + planck_hl(jspec,jlev)
        coeff_up_base =  coeff + planck_hl(jspec,jlev+1)
        coeff_dn_top  = -coeff + planck_hl(jspec,jlev)
        coeff_dn_base = -coeff + planck_hl(jspec,jlev+1)
        source_up(jspec,1,jlev) = coeff_up_top &
             &  - transmittance(jspec,1,jlev) * coeff_up_base
        source_dn(jspec,1,jlev) = coeff_dn_base &
             &  - transmittance(jspec,1,jlev) * coeff_dn_top
      else
        ! Linear limit at low optical depth
        coeff = lw_diffusivity*od(jspec,1,jlev)
        transmittance(jspec,1,jlev) = 1.0_jprb - coeff
        source_up(jspec,1,jlev) = coeff * 0.5_jprb &
             &  * (planck_hl(jspec,jlev)+planck_hl(jspec,jlev+1))
        source_dn(jspec,1,jlev) = source_up(jspec,1,jlev)
      end if
      ! Scale the sources by the area fraction of the region
      source_up(jspec,1,jlev) = region_fracs(1,jlev)*source_up(jspec,1,jlev)
      source_dn(jspec,1,jlev) = region_fracs(1,jlev)*source_dn(jspec,1,jlev)
    end do

    if (region_fracs(1,jlev) < 1.0_jprb) then
      ! We have some cloud
      do jreg = 2,nreg
        ! Scattering solution
        do jspec = 1,nspec
          if (i_two_stream_scheme == ITwoStreamEddington) then
            ! See Meador & Weaver (1980), Table 1; Toon et al. (1989), Table 1
            gamma1 = 1.75_jprb - ssa(jspec,jreg,jlev)*(1.0_jprb + 0.75_jprb*asymmetry(jspec,jlev))
            gamma2 = ssa(jspec,jreg,jlev)*(1.0_jprb - 0.75_jprb*asymmetry(jspec,jlev)) - 0.25_jprb
          else if (i_two_stream_scheme == ITwoStreamScaledWiscombeGrams) then
            ! Wiscombe-Grams backscatter fraction applied to
            ! de-scaled asymmety factor
            factor = 0.5_jprb * (1.0_jprb-0.75_jprb*asymmetry(jspec,jlev) &
                 &                        /(1.0_jprb-asymmetry(jspec,jlev)))
            gamma1 = lw_diffusivity_cloud * (1.0_jprb - ssa(jspec,jreg,jlev)*(1.0_jprb-factor))
            gamma2 = lw_diffusivity_cloud * ssa(jspec,jreg,jlev) * factor
          else
            ! Both Elsasser and Legendre schemes use these formulae,
            ! but with different values for lw_diffusivity_cloud
            ! See Fu et al. (1997), Eqs. 2.9 and 2.10
            factor = (lw_diffusivity_cloud * 0.5_jprb) * ssa(jspec,jreg,jlev)
            gamma1 = lw_diffusivity_cloud - factor*(1.0_jprb + asymmetry(jspec,jlev))
            gamma2 = factor * (1.0_jprb - asymmetry(jspec,jlev))
          end if
          
          k_exponent = sqrt(max((gamma1 - gamma2) * (gamma1 + gamma2), &
               1.E-12_jprb)) ! Eq 18 of Meador & Weaver (1980)
          if (od(jspec,jreg,jlev) > OD_THRESH_2STREAM) then
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
            ! Low optical depth approximation
            reflectance(jspec,jreg,jlev) = gamma2 * od(jspec,jreg,jlev)
            transmittance(jspec,jreg,jlev) &
                 &  = (1.0_jprb - k_exponent*od(jspec,jreg,jlev)) &
                 &  / (1.0_jprb + od(jspec,jreg,jlev)*(gamma1-k_exponent))
            source_up(jspec,jreg,jlev) &
                 &  = (1.0_jprb - reflectance(jspec,jreg,jlev) &
                 &  - transmittance(jspec,jreg,jlev)) &
                 &    * 0.5 * (planck_hl(jspec,jlev) + planck_hl(jspec,jlev+1))
            source_dn(jspec,jreg,jlev) = source_up(jspec,jreg,jlev)
          end if
          ! Scale the sources by the area fraction of the region
          source_up(jspec,jreg,jlev) &
               &  = region_fracs(jreg,jlev)*source_up(jspec,jreg,jlev)
          source_dn(jspec,jreg,jlev) &
               &  = region_fracs(jreg,jlev)*source_dn(jspec,jreg,jlev)
        end do
      end do
    end if
  end do

  if (lhook) call dr_hook('tcrad:calc_reflectance_transmittance',1,hook_handle)

end subroutine calc_reflectance_transmittance


!---------------------------------------------------------------------
! Calculate the transmittance of each layer and region along a path
! with consine of zenith angle "mu", as well as (optionally) the
! emission up from the top of the layer and down through its base.

! WARNING: this routine does not work perfectly in single precision:
! if you run in an independent column configuration in which the
! high cloud-free stratosphere is being used as a cloudy region, the
! OD_THESH value becomes important, and the higher value for
! stability required in single precision is too high for accuracy.

subroutine calc_radiance_trans_source(nspec, nlev, nreg, &
     &  mu, region_fracs, planck_hl, od, ssa, asymmetry, &
     &  flux_up_base, flux_dn_top, &
     &  transmittance, source_up, source_dn)
    
  use yomhook,  only           : lhook, dr_hook, jphook

  ! Parameter
  real(jprb), parameter :: PI = acos(-1.0_jprb)

  ! Inputs

  ! Number of spectral intervals, levels and regions
  integer(jpim), intent(in) :: nspec, nlev, nreg

  ! Cosine of the zenith angle (positive)
  real(jprb) :: mu

  ! Fraction of the gridbox occupied by each region (summing to 1)
  ! at each level
  real(jprb), intent(in), dimension(nreg,nlev) :: region_fracs

  ! Planck function integrated over each spectral interval at each
  ! half-level, in W m-2 (i.e. the flux emitted by a horizontal
  ! black-body surface)
  real(jprb), intent(in), dimension(nspec,nlev+1) :: planck_hl

  ! Optical depth in each region and layer
  real(jprb), intent(in), dimension(nspec,nreg,nlev) :: od

  ! Single scattering albedo in each cloudy region
  real(jprb), intent(in), dimension(nspec,2:nreg,nlev) :: ssa

  ! Asymmetry factor of clouds
  real(jprb), intent(in), dimension(nspec,nlev) :: asymmetry

  ! Upward and downward fluxes at the top and base of each layer and
  ! region, in Watts of power per square metre of the entire
  ! gridbox, so the energy is scaled by the size of each region
  real(jprb), intent(in), dimension(nspec,nreg,nlev) :: flux_up_base, flux_dn_top
  
  ! Outputs

  ! Layer transmittance at the requested zenith angle
  real(jprb), intent(out), dimension(nspec,nreg,nlev) :: transmittance

  ! Optional outputs

  ! Source term up from the top of the layer or down from its base,
  ! in Watts of power per square metre of the entire gridbox, so the
  ! energy is scaled by the size of each region. Since the user may
  ! only require a radiance up or down, these output arguments are
  ! optional.
  real(jprb), intent(out), dimension(nspec,nreg,nlev), optional &
       &  :: source_up, source_dn
  
  ! Local variables
  
  ! Working variables in W m-2
  real(jprb), dimension(nspec,nreg) :: planck_top, planck_base
  real(jprb), dimension(nspec,nreg) :: source_top, source_base

  ! Other working variables
  real(jprb) :: secant, factor, coeff, gamma1, gamma2, k_exponent, rt_factor
  real(jprb) :: exponential, ssa_local, one_minus_kmu
  real(jprb) :: p_same, p_opposite, planck_prime, c1, c2, scaling1, scaling2
    
  ! Maximum number of active regions in a layer (1 in a cloud-free layer)
  integer(jpim) :: max_reg

  ! Loop indices for level, spectral interval and region
  integer(jpim) :: jlev, jspec, jreg

  real(jphook) :: hook_handle

  if (lhook) call dr_hook('tcrad:calc_radiance_trans_source',0,hook_handle)

  secant = 1.0_jprb / mu

  do jlev = 1,nlev

    if (region_fracs(1,jlev) < 1.0_jprb) then
      max_reg = nreg;
    else
      max_reg = 1;
    end if

    if (max_reg > 1) then
      ! Cloudy layer: scale the Planck terms by the region fraction
      ! and also by the single-scattering co-albedo
      planck_top(:,1) = planck_hl(:,jlev) * region_fracs(1,jlev)
      planck_top(:,2:nreg) = spread(planck_hl(:,jlev),2,nreg-1) &
           &  * (1.0_jprb - 0.0*ssa(:,2:nreg,jlev)) &
           &  * spread(region_fracs(2:nreg,jlev),1,nspec)
      planck_base(:,1) = planck_hl(:,jlev+1) * region_fracs(1,jlev)
      planck_base(:,2:nreg) = spread(planck_hl(:,jlev+1),2,nreg-1) &
           &  * (1.0_jprb - 0.0*ssa(:,2:nreg,jlev)) &
           &  * spread(region_fracs(2:nreg,jlev),1,nspec)
    else
      ! Clear layer
      max_reg = 1 ! Redundant
      planck_top(:,1)  = planck_hl(:,jlev)
      planck_base(:,1) = planck_hl(:,jlev+1)
    end if

    do jreg = 1,max_reg
      ! Transmittance of the layer to a beam of radiation
      transmittance(:,jreg,jlev) = exp(-od(:,jreg,jlev)*secant)
    end do
    
    ! Clear region
    if (present(source_up)) then
      do jspec = 1,nspec
        if (od(jspec,1,jlev) > OD_THRESH) then
          planck_prime = (planck_base(jspec,1)-planck_top(jspec,1)) / od(jspec,1,jlev)
          source_up(jspec,1,jlev) &
               &  = planck_top(jspec,1)-planck_base(jspec,1)*transmittance(jspec,1,jlev) &
               &     + planck_prime*mu*(1.0_jprb - transmittance(jspec,1,jlev))
        else
          ! At low optical depths the effective Planck function is
          ! half the top and bottom values, and we avoid the
          ! division by optical depth
          source_up(jspec,1,jlev) = od(jspec,1,jlev) * 0.5_jprb &
               &  * (planck_base(jspec,1)+planck_top(jspec,1)) / mu
        end if
      end do
    end if
    
    if (present(source_dn)) then
      do jspec = 1,nspec
        if (od(jspec,1,jlev) > OD_THRESH) then
          planck_prime = (planck_base(jspec,1)-planck_top(jspec,1)) / od(jspec,1,jlev)
          source_dn(jspec,1,jlev) &
               &  = planck_base(jspec,1)-planck_top(jspec,1)*transmittance(jspec,1,jlev) &
               &     - planck_prime*mu*(1.0_jprb - transmittance(jspec,1,jlev))
        else
          ! At low optical depths the effective Planck function is
          ! half the top and bottom values, and we avoid the
          ! division by optical depth
          source_dn(jspec,1,jlev) = od(jspec,1,jlev) * 0.5_jprb &
               &  * (planck_base(jspec,1)+planck_top(jspec,1)) / mu
        end if
      end do
    end if
    
    ! Cloudy regions
    do jreg = 2,max_reg
      do jspec = 1,nspec
        if (od(jspec,jreg,jlev) > OD_THRESH) then
          if (i_two_stream_scheme == ITwoStreamEddington) then
            ! See Meador & Weaver (1980), Table 1; Toon et al. (1989), Table 1
            gamma1 = 1.75_jprb - ssa(jspec,jreg,jlev)*(1.0_jprb + 0.75_jprb*asymmetry(jspec,jlev))
            gamma2 = ssa(jspec,jreg,jlev)*(1.0_jprb - 0.75_jprb*asymmetry(jspec,jlev)) - 0.25_jprb
          else if (i_two_stream_scheme == ITwoStreamScaledWiscombeGrams) then
            ! Wiscombe-Grams backscatter fraction applied to
            ! de-scaled asymmety factor
            factor = 0.5_jprb * (1.0_jprb-0.75_jprb*asymmetry(jspec,jlev)/(1.0_jprb-asymmetry(jspec,jlev)))
            gamma1 = lw_diffusivity_cloud * (1.0_jprb - ssa(jspec,jreg,jlev)*(1.0_jprb-factor))
            gamma2 = lw_diffusivity_cloud * ssa(jspec,jreg,jlev) * factor
          else
            ! See Fu et al. (1997), Eqs. 2.9 and 2.10; also
            ! "Quadrature" from Toon et al. (1989) Table 1 but with
            ! generalized diffusivity
            factor = (lw_diffusivity_cloud * 0.5_jprb) * ssa(jspec,jreg,jlev)
            gamma1 = lw_diffusivity_cloud - factor*(1.0_jprb + asymmetry(jspec,jlev))
            gamma2 = factor * (1.0_jprb - asymmetry(jspec,jlev))
          end if
          k_exponent = sqrt(max((gamma1 - gamma2) * (gamma1 + gamma2), &
               &  MIN_K_SQUARED)) ! Eq 18 of Meador & Weaver (1980)
          
          ! Phase functions from upwelling flux to upwelling radiance (or down to down)
          p_same     = 1.0_jprb + 3.0_jprb * asymmetry(jspec,jlev) * mu / lw_diffusivity_cloud
          ! Phase function from downwelling flux to upwelling radiance (or up to down)
          p_opposite = 1.0_jprb - 3.0_jprb * asymmetry(jspec,jlev) * mu / lw_diffusivity_cloud
          
          ! Compute the coefficients, c1 and c2, of the exponentials
          ! describing the analytic profile of upwelling and
          ! downwelling fluxes within the layer
          planck_prime = (planck_base(jspec,jreg)-planck_top(jspec,jreg)) / od(jspec,jreg,jlev)
          exponential = exp(-k_exponent*od(jspec,jreg,jlev))
          coeff = planck_prime / (gamma1+gamma2)
          rt_factor = 1.0_jprb / (k_exponent + gamma1 + (k_exponent-gamma1) &
               &                  *exponential*exponential)
          factor = exponential * gamma2 / (gamma1 + k_exponent)
          c1 = rt_factor * (flux_up_base(jspec,jreg,jlev) - factor*flux_dn_top(jspec,jreg,jlev) &
               &  - (planck_base(jspec,jreg)+coeff) + factor*(planck_top(jspec,jreg)-coeff))
          c2 = rt_factor * (flux_dn_top(jspec,jreg,jlev) - factor*flux_up_base(jspec,jreg,jlev) &
               &  -(planck_top(jspec,jreg)-coeff) + factor*(planck_base(jspec,jreg)+coeff))
          
          ! Scaling factors for the coefficients
          one_minus_kmu = 1.0_jprb - k_exponent*mu
          scaling1 = (exponential - transmittance(jspec,jreg,jlev)) &
               &  / merge(one_minus_kmu, epsilon(1.0_jprb), abs(one_minus_kmu)>epsilon(1.0_jprb))
          scaling2 = (1.0_jprb-exponential*transmittance(jspec,jreg,jlev))/(1.0_jprb+k_exponent*mu)

          if (present(source_up)) then
            ! Direct emission plus scattering from the part of the
            ! fluxes due to internal emission and having a linear
            ! structure
            source_up(jspec,jreg,jlev) &
                 &  = 0.5_jprb*ssa(jspec,jreg,jlev)*(1.0_jprb - transmittance(jspec,jreg,jlev)) &
                 &    * planck_prime*(p_same-p_opposite)/(gamma1+gamma2) &
                 &  + planck_top(jspec,jreg)-planck_base(jspec,jreg)*transmittance(jspec,jreg,jlev) &
                 &     + planck_prime*mu*(1.0_jprb - transmittance(jspec,jreg,jlev))
            
            ! Scattering from the exponential part of the flux,
            ! whether caused by external or internal sources
            source_up(jspec,jreg,jlev) = source_up(jspec,jreg,jlev) &
                 &  + 0.5_jprb*ssa(jspec,jreg,jlev) &
                 &  * (p_same     * ((gamma1+k_exponent)*scaling1*c1 +  gamma2*scaling2*c2) &
                 &    +p_opposite * ( gamma2*scaling1*c1             + (gamma1+k_exponent)*scaling2*c2))
          end if
          if (present(source_dn)) then
            ! Direct emission plus scattering from the part of the
            ! fluxes due to internal emission and having a linear
            ! structure
            source_dn(jspec,jreg,jlev) &
                 &  = - 0.5_jprb*ssa(jspec,jreg,jlev)*(1.0_jprb - transmittance(jspec,jreg,jlev)) &
                 &    * planck_prime*(p_same-p_opposite)/(gamma1+gamma2) &
                 &  + planck_base(jspec,jreg)-planck_top(jspec,jreg)*transmittance(jspec,jreg,jlev) &
                 &     - planck_prime*mu*(1.0_jprb - transmittance(jspec,jreg,jlev))
              
            ! Scattering from the exponential part of the flux,
            ! whether caused by external or internal sources
            source_dn(jspec,jreg,jlev) = source_dn(jspec,jreg,jlev) &
                 &  + 0.5_jprb*ssa(jspec,jreg,jlev) &
                 &  * (p_opposite * ((gamma1+k_exponent)*scaling2*c1 +  gamma2*scaling1*c2) &
                 &    +p_same     * ( gamma2*scaling2*c1             + (gamma1+k_exponent)*scaling1*c2))
          end if
            
        else
          ! Low optical depth approximation: emission only
          if (present(source_up)) then
            source_up(jspec,jreg,jlev) = od(jspec,jreg,jlev) &
                 &  * 0.5_jprb*(planck_base(jspec,jreg)+planck_top(jspec,jreg)) / mu
          end if
          if (present(source_dn)) then
            source_dn(jspec,jreg,jlev) = od(jspec,jreg,jlev) &
                 &  * 0.5_jprb*(planck_base(jspec,jreg)+planck_top(jspec,jreg)) / mu
          end if
        end if
          
      end do ! jspec
    end do ! jreg
  end do ! jlev
  
  if (lhook) call dr_hook('tcrad:calc_radiance_trans_source',1,hook_handle)
  
end subroutine calc_radiance_trans_source


!---------------------------------------------------------------------
! Compute the region fractions and the optical depth scalings for the
! optically "thick" and "thin" regions of a Tripleclouds
! representation of a sub-grid PDF of cloud optical depth. Following
! Shonk and Hogan (2008), the 16th percentile is used for the thin
! region, and the formulas estimate this for both lognormal and gamma
! distributions. However, an adjustment is needed for the gamma
! distribution at large fractional standard deviations.
subroutine calc_region_properties(nlev, &
     &  cloud_fraction, do_gamma, frac_std, &
     &  reg_fracs, od_scaling, cloud_fraction_threshold)

  use parkind1,     only : jprb
  use yomhook, only : lhook, dr_hook, jphook

  ! Minimum od_scaling in the case of a Gamma distribution
  real(jprb), parameter :: MinGammaODScaling = 0.025_jprb

  ! At large fractional standard deviations (FSDs), we cannot capture
  ! the behaviour of a gamma distribution with two equally weighted
  ! points; we need to weight the first ("lower") point more.  The
  ! weight of the first point is normally 0.5, but for FSDs between
  ! 1.5 and 3.725 it rises linearly to 0.9, and for higher FSD it is
  ! capped at this value.  The weight of the second point is one minus
  ! the first point.
  real(jprb), parameter :: MinLowerFrac      = 0.5_jprb
  real(jprb), parameter :: MaxLowerFrac      = 0.9_jprb
  real(jprb), parameter :: FSDAtMinLowerFrac = 1.5_jprb
  real(jprb), parameter :: FSDAtMaxLowerFrac = 3.725_jprb
  ! Between FSDAtMinLowerFrac and FSDAtMaxLowerFrac,
  ! LowerFrac=LowerFracFSDIntercept+FSD*LowerFracFSDGradient
  real(jprb), parameter :: LowerFracFSDGradient &
       &  = (MaxLowerFrac-MinLowerFrac) / (FSDAtMaxLowerFrac-FSDAtMinLowerFrac)
  real(jprb), parameter :: LowerFracFSDIntercept &
       &  = MinLowerFrac - FSDAtMinLowerFrac*LowerFracFSDGradient
  
  ! Number of levels and regions
  integer, intent(in) :: nlev
  
  ! Do we do a lognormal or gamma distribution?
  logical, intent(in) :: do_gamma

  ! Cloud fraction, i.e. the fraction of the gridbox assigned to all
  ! regions numbered 2 and above (region 1 is clear sky)
  real(jprb), intent(in), dimension(:)  :: cloud_fraction ! (nlev)

  ! Fractional standard deviation of in-cloud water content
  real(jprb), intent(in), dimension(:)  :: frac_std       ! (nlev)

  ! Fractional area coverage of each region
  real(jprb), intent(out) :: reg_fracs(NREGION,nlev)

  ! Optical depth scaling for the cloudy regions
  real(jprb), intent(out) :: od_scaling(2:NREGION,nlev)

  ! Regions smaller than this are ignored
  real(jprb), intent(in), optional :: cloud_fraction_threshold
  
  ! In case the user doesn't supply cloud_fraction_threshold we use
  ! a default value
  real(jprb) :: frac_threshold
  
  ! Loop indices
  integer :: jlev
  
  real(jphook) :: hook_handle

  if (lhook) call dr_hook('tcrad:calc_region_properties',0,hook_handle)

  if (present(cloud_fraction_threshold)) then
    frac_threshold = cloud_fraction_threshold
  else
    frac_threshold = 10.0_jprb * epsilon(1.0_jprb)
  end if
  
  ! Two cloudy regions with optical depth scaled by 1-x and 1+x.
  
  ! Simple version which fails when fractional_std >= 1:
  !od_scaling(2) = 1.0_jprb-cloud%fractional_std(jlev)
  ! According to Shonk and Hogan (2008), 1-FSD should correspond to
  ! the 16th percentile. 
  if (.not. do_gamma) then
    ! If we treat the distribution as a lognormal such that the
    ! equivalent Normal has a mean mu and standard deviation sigma,
    ! then the 16th percentile of the lognormal is very close to
    ! exp(mu-sigma).
    do jlev = 1,nlev
      if (cloud_fraction(jlev) < frac_threshold) then
        reg_fracs(1,jlev)    = 1.0_jprb
        reg_fracs(2:3,jlev)  = 0.0_jprb
        od_scaling(2:3,jlev) = 1.0_jprb
      else
        reg_fracs(1,jlev) = 1.0_jprb - cloud_fraction(jlev)
        reg_fracs(2:3,jlev) = cloud_fraction(jlev)*0.5_jprb
        od_scaling(2,jlev) &
             &  = exp(-sqrt(log(frac_std(jlev)**2+1))) &
             &  / sqrt(frac_std(jlev)**2+1)
        od_scaling(3,jlev) = 2.0_jprb - od_scaling(2,jlev)
      end if
    end do
  else
    ! If we treat the distribution as a gamma then the 16th percentile
    ! is close to the following.  Note that because it becomes
    ! vanishingly small for FSD >~ 2, we have a lower limit of 1/40,
    ! and for higher FSDs reduce the fractional cover of the denser
    ! region - see region_fractions routine below
    do jlev = 1,nlev
      if (cloud_fraction(jlev) < frac_threshold) then
        reg_fracs(1,jlev)    = 1.0_jprb
        reg_fracs(2:3,jlev)  = 0.0_jprb
        od_scaling(2:3,jlev) = 1.0_jprb
      else
        ! Fraction of the clear-sky region
        reg_fracs(1,jlev) = 1.0_jprb - cloud_fraction(jlev)
        ! Fraction and optical-depth scaling of the lower of the two
        ! cloudy regions: the fraction of the thicker and thinner
        ! cloudy regions are not necessarily of the same area,
        ! following the appendix of Hogan et al. (2019).
        reg_fracs(2,jlev) = cloud_fraction(jlev) &
             &  * max(MinLowerFrac, min(MaxLowerFrac, &
             &  LowerFracFSDIntercept + frac_std(jlev)*LowerFracFSDGradient))
        od_scaling(2,jlev) = MinGammaODScaling &
             &  + (1.0_jprb - MinGammaODScaling) &
             &    * exp(-frac_std(jlev)*(1.0_jprb + 0.5_jprb*frac_std(jlev) &
             &                     *(1.0_jprb+0.5_jprb*frac_std(jlev))))
        ! Fraction of the upper of the two cloudy regions
        reg_fracs(3,jlev) = 1.0_jprb-reg_fracs(1,jlev)-reg_fracs(2,jlev)
        ! Ensure conservation of the mean optical depth
        od_scaling(3,jlev) = (cloud_fraction(jlev) &
             &  -reg_fracs(2,jlev)*od_scaling(2,jlev)) / reg_fracs(3,jlev)
      end if
    end do
  end if
  
  if (lhook) call dr_hook('tcrad:calc_region_properties',1,hook_handle)

end subroutine calc_region_properties


!---------------------------------------------------------------------
! Treat A as an m-by-m square matrix and b as n NREGION-element vectors
! (with the n dimension varying fastest), and perform n matrix-vector
! multiplications
pure function singlemat_x_vec(n,A,b)

  use parkind1, only : jprb

  implicit none

  integer,    intent(in)                             :: n
  real(jprb), intent(in), dimension(NREGION,NREGION) :: A
  real(jprb), intent(in), dimension(n,NREGION)       :: b
  real(jprb),             dimension(n,NREGION)       :: singlemat_x_vec
  
  integer    :: j1, j2
  
  ! Array-wise assignment
  singlemat_x_vec = 0.0_jprb
  
  do j1 = 1,NREGION
    do j2 = 1,NREGION
      singlemat_x_vec(:,j1) = singlemat_x_vec(:,j1) + A(j1,j2)*b(:,j2)
    end do
  end do
  
end function singlemat_x_vec


!---------------------------------------------------------------------
! Calculate a matrix expressing the overlap of regions in adjacent
! layers, using the Hogan and Illingworth (2000) "alpha" overlap
! parameter, but allowing for the two cloudy regions in the
! Tripleclouds assumption to have different areas
pure function calc_alpha_overlap_matrix(op, op_inhom, &
     &  frac_upper, frac_lower) result(overlap_matrix)

  use parkind1, only : jprb
    
  ! Overlap parameter for cloud boundaries and for internal
  ! inhomogeneities
  real(jprb), intent(in) :: op, op_inhom

  ! Fraction of the gridbox occupied by each region in the upper and
  ! lower layers
  real(jprb), intent(in), dimension(NREGION) :: frac_upper, frac_lower

  ! Output overlap matrix
  real(jprb) :: overlap_matrix(NREGION,NREGION)

  ! Combined cloud cover of pair of layers
  real(jprb) :: pair_cloud_cover

  ! Cloud fraction of upper and lower layers
  real(jprb) :: cf_upper, cf_lower

  ! One divided by cloud fraction
  real(jprb) :: one_over_cf

  ! Fraction of domain with cloud in both layers
  real(jprb) :: frac_both

  cf_upper = sum(frac_upper(2:NREGION))
  cf_lower = sum(frac_lower(2:NREGION))

  pair_cloud_cover = op*max(cf_upper,cf_lower) &
       &  + (1.0_jprb - op) &
       &  * (cf_upper+cf_lower-cf_upper*cf_lower)
  
  ! Clear in both layers
  overlap_matrix(1,1) = 1.0_jprb - pair_cloud_cover

  ! Clear in upper layer, cloudy in lower layer
  one_over_cf = 1.0_jprb / max(cf_lower, 1.0e-6_jprb)
  overlap_matrix(1,2) = (pair_cloud_cover - cf_upper) &
       &              * frac_lower(2) * one_over_cf
  overlap_matrix(1,3) = (pair_cloud_cover - cf_upper) &
       &              * frac_lower(3) * one_over_cf
  ! Clear in lower layer, cloudy in upper layer
  one_over_cf = 1.0_jprb / max(cf_upper, 1.0e-6_jprb)
  overlap_matrix(2,1) = (pair_cloud_cover - cf_lower) &
       &              * frac_upper(2) * one_over_cf
  overlap_matrix(3,1) = (pair_cloud_cover - cf_lower) &
       &              * frac_upper(3) * one_over_cf
  ! Cloudy in both layers: frac_both is the fraction of the
  ! gridbox with cloud in both layers
  frac_both = cf_upper + cf_lower - pair_cloud_cover
  ! Treat low and high optical-depth regions within frac_both as one
  ! treats clear and cloudy skies in the whole domain; redefine the
  ! following variables treating the high optical-depth region as
  ! the cloud
  cf_upper = frac_upper(3) / max(cf_upper, 1.0e-6_jprb)
  cf_lower = frac_lower(3) / max(cf_lower, 1.0e-6_jprb)
  pair_cloud_cover = op_inhom*max(cf_upper,cf_lower) &
       &  + (1.0_jprb - op_inhom) &
       &  * (cf_upper+cf_lower-cf_upper*cf_lower)
  ! Assign overlaps for this 2x2 section of the 3x3 matrix as for
  ! the 2-region case above, but multiplied by frac_both
  overlap_matrix(2,2) = frac_both * (1.0_jprb - pair_cloud_cover)
  overlap_matrix(2,3) = frac_both * (pair_cloud_cover - cf_upper)
  overlap_matrix(3,2) = frac_both * (pair_cloud_cover - cf_lower)
  overlap_matrix(3,3) = frac_both * (cf_upper+cf_lower-pair_cloud_cover)

end function calc_alpha_overlap_matrix


!---------------------------------------------------------------------
! Compute the upward and downward overlap matrices u_overlap and
! v_overlap, respectively, where u_overlap is defined such that
! y=u_overlap*x, where x is a vector of upwelling fluxes in each
! region just below an interface, and y is a vector of upwelling
! fluxes in each region just above that interface. For nlev model
! levels there are nlev+1 interfaces including the ground and
! top-of-atmosphere, and so that is one of the dimensions of u_overlap
! and v_overlap.
subroutine calc_overlap_matrices(nlev, &
     &     region_fracs, overlap_param, u_overlap, v_overlap, &
     &     decorrelation_scaling, &
     &     cloud_fraction_threshold, cloud_cover)

  use parkind1,     only : jprb
  use yomhook, only : lhook, dr_hook, jphook

  ! Number of levels and regions
  integer,  intent(in) :: nlev

  ! Area fraction of each region: region 1 is clear sky, and 2+ are
  ! the cloudy regions (only one or two cloudy regions are supported)
  real(jprb), intent(in), dimension(1:NREGION,nlev)  :: region_fracs

  ! The overlap "alpha" overlap parameter of Hogan & Illingworth
  real(jprb), intent(in), dimension(:)  :: overlap_param  ! (nlev-1)

  ! Output overlap matrices
  real(jprb), intent(out), dimension(NREGION,NREGION,nlev+1) &
       &  :: u_overlap, v_overlap
  
  ! For regions 2 and above, the overlap decorrelation length for
  ! cloud boundaries is scaled by this amount to obtain the overlap
  ! decorrelation length for cloud inhomogeneities. Typically this
  ! number is 0.5, but if omitted it will be assumed to be one (same
  ! decorrelation for cloud boundaries and in-cloud inhomogeneities)
  real(jprb), intent(in), optional :: decorrelation_scaling

  ! Regions smaller than this are ignored
  real(jprb), intent(in), optional :: cloud_fraction_threshold

  ! The diagnosed cloud cover is an optional output
  real(jprb), intent(out), optional :: cloud_cover

  ! Loop indices for column, level, region and the regions in the
  ! upper and lower layers for an interface
  integer  :: jlev, jupper, jlower

  ! Overlap matrix (non-directional)
  real(jprb) :: overlap_matrix(NREGION,NREGION)

  ! Fraction of the gridbox occupied by each region in the upper and
  ! lower layers for an interface
  real(jprb) :: frac_upper(NREGION), frac_lower(NREGION)

  ! Beta overlap parameter for each region
  real(jprb) :: op(NREGION)

  ! In case the user doesn't supply cloud_fraction_threshold we use a
  ! default value
  real(jprb) :: frac_threshold

  ! The decorrelation scaling to use, in case decorrelation_scaling
  ! was not provided
  real(jprb) :: used_decorrelation_scaling

  real(jphook) :: hook_handle

  if (lhook) call dr_hook('tcrad:calc_overlap_matrices',0,hook_handle)

  if (present(decorrelation_scaling)) then
    used_decorrelation_scaling = decorrelation_scaling
  else
    used_decorrelation_scaling = 1.0_jprb
  end if
  
  if (present(cloud_fraction_threshold)) then
    frac_threshold = cloud_fraction_threshold
  else
    frac_threshold = 10.0_jprb * epsilon(1.0_jprb)
  end if
    
  ! Outer space is treated as one clear-sky region, so the fractions
  ! are assigned as such
  frac_upper(1) = 1.0_jprb
  frac_upper(2:NREGION) = 0.0_jprb
  
  ! Overlap parameter is irrelevant when there is only one region in
  ! the upper layer
  op = 1.0_jprb

  ! Loop down through the atmosphere, where jlev indexes each
  ! half-level starting at 1 for the top-of-atmosphere, as well as
  ! indexing each level starting at 1 for the top-most level.
  do jlev = 1,nlev+1
    ! Fraction of each region just below the interface
    if (jlev > nlev) then
      ! We are at the surface: treat as a single clear-sky region
      frac_lower(1) = 1.0_jprb
      frac_lower(2:NREGION) = 0.0_jprb
    else
      frac_lower = region_fracs(1:NREGION,jlev)
    end if
    
    ! Compute the overlap parameter of the interface just below the
    ! current full level
    if (jlev == 1 .or. jlev > nlev) then
      ! We are at the surface or top-of-atmosphere: overlap
      ! parameter is irrelevant
      op = 1.0_jprb
    else
      ! We are not at the surface
      op(1) = overlap_param(jlev-1)
      ! For cloudy regions, scale the cloud-boundary overlap
      ! parameter to obtain the cloud-inhomogeneity overlap
      ! parameter as follows
      if (op(1) >= 0.0_jprb) then
        op(2:NREGION) = op(1)**(1.0_jprb/used_decorrelation_scaling)
      else
        op(2:NREGION) = op(1)
      end if
    end if
     
    overlap_matrix = calc_alpha_overlap_matrix( &
         &  op(1), op(2), frac_upper, frac_lower)

    ! Convert to directional overlap matrices
    do jupper = 1,NREGION
      do jlower = 1,NREGION
        if (frac_lower(jlower) >= frac_threshold) then
          u_overlap(jupper,jlower,jlev) = overlap_matrix(jupper,jlower) &
               &  / frac_lower(jlower)
        else
          u_overlap(jupper,jlower,jlev) = 0.0_jprb
        end if
        if (frac_upper(jupper) >= frac_threshold) then
          v_overlap(jlower,jupper,jlev) = overlap_matrix(jupper,jlower) &
               &  / frac_upper(jupper)
        else
          v_overlap(jlower,jupper,jlev) = 0.0_jprb
        end if
      end do
    end do
    frac_upper = frac_lower
    
  end do ! levels
  
  ! Compute cloud cover from one of the directional overlap matrices
  if (present(cloud_cover)) then
    cloud_cover = 1.0_jprb - product(v_overlap(1,1,:))
  end if
  
  if (lhook) call dr_hook('tcrad:calc_overlap_matrices',1,hook_handle)
  
end subroutine calc_overlap_matrices


!-------------------------------------------------------------
! Compute fluxes at top and base of each layer and region from
! precomputed layer transmittance, reflectance and upward/downward
! layer sources, and precomputed upward and downward cloud overlap
! matrices, using the Tripleclouds two-stream method of Shonk and
! Hogan (2008).
subroutine calc_two_stream_flux(nspec, nlev, surf_emission, surf_albedo, &
     &  reflectance, transmittance, source_up, source_dn, &
     &  is_cloud_free_layer, u_overlap, v_overlap, &
     &  flux_up_base, flux_dn_base, flux_up_top, flux_dn_top)

  use parkind1, only           : jpim, jprb
  use yomhook, only : lhook, dr_hook, jphook
  
  implicit none
  
  ! Inputs

  ! Number of spectral intervals and levels
  integer(jpim), intent(in) :: nspec, nlev

  ! Surface upwards emission, in W m-2 (i.e. emissivity multiplied by
  ! Planck function at the surface skin temperature) integrated across
  ! each spectral interval, and albedo in the same intervals
  real(jprb), intent(in),  dimension(nspec) :: surf_emission, surf_albedo
  
  ! Reflectance and transmittance of each layer and region
  real(jprb), intent(in),  dimension(nspec,NREGION,nlev) :: reflectance
  real(jprb), intent(in),  dimension(nspec,NREGION,nlev) :: transmittance

  ! Rate of emission upward from the top of a layer and downwards from
  ! its base, due to emission within the layer, in Watts of power per
  ! square metre of the entire gridbox, so the energy is scaled by the
  ! size of each region
  real(jprb), intent(in),  dimension(nspec,NREGION,nlev) :: source_up
  real(jprb), intent(in),  dimension(nspec,NREGION,nlev) :: source_dn

  ! Where are the cloud free layers?  Includes dummy layers at 0 and
  ! nlev+1 which are also cloud free
  logical,    intent(in) :: is_cloud_free_layer(0:nlev+1)

  ! U and V overlap matrices, defined by Shonk and Hogan (2008)
  real(jprb), intent(in),  dimension(NREGION,NREGION,nlev+1) :: u_overlap, v_overlap
  
  ! Outputs

  ! Upwelling and downwelling fluxes at the top and base of each layer
  ! in the regions of that layer, in Watts of power per square metre
  ! of the entire gridbox, so the energy is scaled by the size of each
  ! region
  real(jprb), intent(out), dimension(nspec,NREGION,nlev) :: flux_up_base, flux_dn_base
  real(jprb), intent(out), dimension(nspec,NREGION,nlev) :: flux_up_top, flux_dn_top
  
  ! Local variables

  ! Total albedo of the atmosphere/surface just above a layer
  ! interface with respect to downwelling diffuse radiation at that
  ! interface, where level index = 1 corresponds to the
  ! top-of-atmosphere
  real(jprb), dimension(nspec, NREGION, nlev+1) :: total_albedo

  ! Upwelling radiation just above a layer interface due to emission
  ! below that interface, where level index = 1 corresponds to the
  ! top-of-atmosphere
  real(jprb), dimension(nspec, NREGION, nlev+1) :: total_source

  ! Total albedo and source of the atmosphere just below a layer interface
  real(jprb), dimension(nspec, NREGION) :: total_albedo_below, total_source_below

  ! Term in the adding method to replace divisions by multiplications
  real(jprb), dimension(nspec, NREGION) :: inv_denom

  ! Index of highest layer containing cloud
  integer(jpim) :: icloudtop

  ! Loop indices for spectral interval, level and region
  integer(jpim) :: jspec, jlev, jreg, jreg2

  real(jphook) :: hook_handle

  if (lhook) call dr_hook('tcrad:calc_two_stream_flux',0,hook_handle)

  ! --------------------------------------------------------
  ! Section 1: Prepare variables and arrays
  ! --------------------------------------------------------

  ! Find cloud top
  icloudtop = nlev
  do jlev = 1,nlev
    if (.not. is_cloud_free_layer(jlev)) then
      icloudtop = jlev
      exit
    end if
  end do

  ! --------------------------------------------------------
  ! Section 2: Clear-sky downwelling fluxes from TOA to cloud top
  ! --------------------------------------------------------
  ! Initialize outputs
  flux_up_base = 0.0_jprb
  flux_dn_base = 0.0_jprb
  flux_up_top  = 0.0_jprb
  flux_dn_top  = 0.0_jprb

  ! Downwelling fluxes from TOA to cloud top, neglecting clear-air
  ! scattering. Initial flux_dn_top in the clear region is already
  ! zero.
  do jlev = 1,icloudtop
    if (jlev > 1) then
      flux_dn_top(:,1,jlev) = flux_dn_base(:,1,jlev-1)
    end if
    flux_dn_base(:,1,jlev) = source_dn(:,1,jlev) &
         &                 + transmittance(:,1,jlev)*flux_dn_top(:,1,jlev)
  end do

  ! --------------------------------------------------------
  ! Section 3: Compute total sources and albedos up to cloud top
  ! --------------------------------------------------------
  total_albedo = 0.0_jprb
  total_source = 0.0_jprb
  ! Calculate the upwelling radiation emitted by the surface, and copy
  ! the surface albedo into total_albedo
  do jreg = 1,NREGION
    do jspec = 1,nspec
      total_source(jspec,jreg,nlev+1) = u_overlap(jreg,1,nlev+1)*surf_emission(jspec)
      total_albedo(jspec,jreg,nlev+1) = surf_albedo(jspec)
    end do
  end do

  ! Work up from the surface computing the total albedo of the
  ! atmosphere and the total upwelling due to emission below each
  ! level below using the Adding Method
  do jlev = nlev,icloudtop,-1
    total_albedo_below = 0.0_jprb

    if (is_cloud_free_layer(jlev)) then
      total_albedo_below = 0.0_jprb
      total_source_below = 0.0_jprb
      do jspec = 1,nspec
        inv_denom(jspec,1) = 1.0_jprb &
             &  / (1.0_jprb - total_albedo(jspec,1,jlev+1)*reflectance(jspec,1,jlev))
        total_albedo_below(jspec,1) = reflectance(jspec,1,jlev) &
             &  + transmittance(jspec,1,jlev)*transmittance(jspec,1,jlev)*total_albedo(jspec,1,jlev+1) &
             &  * inv_denom(jspec,1)
        total_source_below(jspec,1) = source_up(jspec,1,jlev) &
             &  + transmittance(jspec,1,jlev)*(total_source(jspec,1,jlev+1) &
             &  + total_albedo(jspec,1,jlev+1)*source_dn(jspec,1,jlev)) &
             &  * inv_denom(jspec,1)
      end do
    else
      inv_denom = 1.0_jprb / (1.0_jprb - total_albedo(:,:,jlev+1)*reflectance(:,:,jlev))
      total_albedo_below = reflectance(:,:,jlev) &
           &  + transmittance(:,:,jlev)*transmittance(:,:,jlev)*total_albedo(:,:,jlev+1) &
           &  * inv_denom
      total_source_below = source_up(:,:,jlev) &
           &  + transmittance(:,:,jlev)*(total_source(:,:,jlev+1) &
           &  + total_albedo(:,:,jlev+1)*source_dn(:,:,jlev)) &
           &  * inv_denom
    end if
    
    ! Account for cloud overlap when converting albedo below a layer
    ! interface to the equivalent values just above
    if (is_cloud_free_layer(jlev) .and. is_cloud_free_layer(jlev-1)) then
      total_albedo(:,:,jlev) = total_albedo_below(:,:)
      total_source(:,:,jlev) = total_source_below(:,:)
    else
      total_source(:,:,jlev) = singlemat_x_vec(nspec,&
           &  u_overlap(:,:,jlev), total_source_below)
      ! Use overlap matrix and exclude "anomalous" horizontal
      ! transport described by Shonk & Hogan (2008).  Therefore, the
      ! operation we perform is essentially diag(total_albedo) =
      ! matmul(transpose(v_overlap), diag(total_albedo_below)).
      do jreg = 1,NREGION
        do jreg2 = 1,NREGION
          total_albedo(:,jreg,jlev) &
               &  = total_albedo(:,jreg,jlev) &
               &  + total_albedo_below(:,jreg2) * v_overlap(jreg2,jreg,jlev)
        end do
      end do
      
    end if
        
  end do ! Reverse loop over levels

  ! --------------------------------------------------------
  ! Section 4: Compute fluxes up to top-of-atmosphere
  ! --------------------------------------------------------

  if (icloudtop > 1) then
    ! Compute the fluxes in the layer just above the highest cloud
    flux_up_base(:,1,icloudtop-1) = total_source(:,1,icloudtop) & 
         &  + total_albedo(:,1,icloudtop)*flux_dn_base(:,1,icloudtop-1)
    flux_up_top(:,1,icloudtop-1) = source_up(:,1,icloudtop-1) &
         &  + transmittance(:,1,icloudtop-1)*flux_up_base(:,1,icloudtop-1)
    ! Compute fluxes in remaining layers up to TOA
    do jlev = icloudtop-2,1,-1
      flux_up_base(:,1,jlev) = flux_up_top(:,1,jlev+1)
      flux_up_top(:,1,jlev) = source_up(:,1,jlev) &
           &                + transmittance(:,1,jlev)*flux_up_base(:,1,jlev)
    end do
  end if

  ! --------------------------------------------------------
  ! Section 5: Compute fluxes from cloud top down to surface
  ! --------------------------------------------------------

  ! Copy over downwelling spectral fluxes at top of first scattering
  ! layer, using overlap matrix to translate to the regions of the
  ! first layer of cloud
  if (icloudtop > 1) then
    do jreg = 1,NREGION
      flux_dn_top(:,jreg,icloudtop) &
           &  = v_overlap(jreg,1,icloudtop)*flux_dn_base(:,1,icloudtop-1)
    end do
    ! else the highest layer is cloudy, in which case
    ! flux_dn_top(:,:,jlev=1) is already set to zero
  end if


  ! Final loop back down through the atmosphere to compute fluxes
  do jlev = icloudtop,nlev

    if (is_cloud_free_layer(jlev)) then
      do jspec = 1,nspec
        flux_dn_base(jspec,1,jlev) = (transmittance(jspec,1,jlev)*flux_dn_top(jspec,1,jlev) &
             &  + reflectance(jspec,1,jlev)*total_source(jspec,1,jlev+1) &
             &  + source_dn(jspec,1,jlev) ) &
             &  / (1.0_jprb - reflectance(jspec,1,jlev)*total_albedo(jspec,1,jlev+1))
        flux_up_base(jspec,1,jlev) = total_source(jspec,1,jlev+1) &
             &  + flux_dn_base(jspec,1,jlev)*total_albedo(jspec,1,jlev+1)
        flux_up_top(jspec,1,jlev) = flux_up_base(jspec,1,jlev)*transmittance(jspec,1,jlev) &
             &  + source_up(jspec,1,jlev) + flux_dn_top(jspec,1,jlev)*reflectance(jspec,1,jlev)
      end do
    else
      flux_dn_base(:,:,jlev) = (transmittance(:,:,jlev)*flux_dn_top(:,:,jlev) &
           &     + reflectance(:,:,jlev)*total_source(:,:,jlev+1) + source_dn(:,:,jlev) ) &
           &  / (1.0_jprb - reflectance(:,:,jlev)*total_albedo(:,:,jlev+1))
      flux_up_base(:,:,jlev) = total_source(:,:,jlev+1) &
           &  + flux_dn_base(:,:,jlev)*total_albedo(:,:,jlev+1)
      flux_up_top(:,:,jlev) = flux_up_base(:,:,jlev)*transmittance(:,:,jlev) &
           &  + source_up(:,:,jlev) + flux_dn_top(:,:,jlev)*reflectance(:,:,jlev)
    end if
    
    if (jlev < nlev) then
      if (.not. (is_cloud_free_layer(jlev) .and. is_cloud_free_layer(jlev+1))) then
        ! Account for overlap rules in translating fluxes just above a
        ! layer interface to the values just below
        flux_dn_top(:,:,jlev+1) = singlemat_x_vec(nspec, &
             &  v_overlap(:,:,jlev+1), flux_dn_base(:,:,jlev))
      else 
        ! Two clear-sky layers: copy the fluxes
        flux_dn_top(:,1,jlev+1) = flux_dn_base(:,1,jlev)
      end if
    end if

  end do

  if (lhook) call dr_hook('tcrad:calc_two_stream_flux',1,hook_handle)

end subroutine calc_two_stream_flux


!------------------------------------------------------------------
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
  use yomhook, only : lhook, dr_hook, jphook
  
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

  ! Loop index for level
  integer(jpim) :: jlev

  real(jphook) :: hook_handle

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
  use yomhook, only : lhook, dr_hook, jphook
  
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

  ! Loop index for level
  integer(jpim) :: jlev

  real(jphook) :: hook_handle

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


!---------------------------------------------------------------------
! Compute the TOA or surface radiance including the effects of scattering
subroutine calc_radiance(nspec, nlev, surf_emission, surf_albedo, planck_hl, &
     &  cloud_fraction, fractional_std, &
     &  od_clear, od_cloud, ssa_cloud, asymmetry_cloud, &
     &  overlap_param, mu, radiance, cloud_cover, &
     &  do_specular_surface)

  use parkind1, only           : jpim, jprb
  use yomhook, only : lhook, dr_hook, jphook

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

  ! Profile of cloud fraction 
  real(jprb), intent(in), dimension(nlev) :: cloud_fraction

  ! Profile of the fractional standard deviation (i.e. standard
  ! deviation divided by mean) of the horizontal in-cloud water
  ! content (or cloud extinction coefficient) distribution.
  real(jprb), intent(in), dimension(nlev) :: fractional_std

  ! Layer optical depth of gas and aerosol
  real(jprb), intent(in), dimension(nspec,nlev) :: od_clear

  ! Layer optical depth of cloud averaged only over the cloudy part
  ! of the gridbox, single scattering albedo and asymmetry
  ! factor. If delta-Eddington scaling is required then this should
  ! already have been done.
  real(jprb), intent(in), dimension(nspec,nlev) :: od_cloud, &
       &  ssa_cloud, asymmetry_cloud

  ! Overlap parameter governing how clouds in adjacent layers are
  ! overlapped - this is the "alpha" of Hogan and Illingworth
  ! (2000). It is defined only between layers, hence the nlev-1
  ! elements.
  real(jprb), intent(in), dimension(nlev-1) :: overlap_param

  ! Cosine of the sensor zenith angle
  real(jprb), intent(in) :: mu

  ! Outputs

  ! Radiances in W m sr-1
  real(jprb), intent(out), dimension(nspec) :: radiance

  ! Return cloud cover computed from cloud fraction profile and
  ! overlap rules
  real(jprb), intent(out), optional :: cloud_cover

  ! Optional inputs

  ! Is scattering from the surface treated specularly?  Appropriate
  ! for microwave scattering by the sea surface. Default is NO.
  logical, intent(in), optional :: do_specular_surface

  ! Local variables

  ! Combined gas/aerosol/cloud optical depth in each region
  real(jprb), dimension(nspec,NREGION,nlev)   :: od

  ! Single scattering albedo of the cloudy regions (ssa=0 in the
  ! clear region)
  real(jprb), dimension(nspec,2:NREGION,nlev) :: ssa

  ! Reflectance and transmittance of each layer and region
  real(jprb), dimension(nspec,NREGION,nlev) :: reflectance, transmittance

  ! Rate of emission up from the top or down through the base of
  ! each layer and region (W m-2)
  real(jprb), dimension(nspec,NREGION,nlev) :: source_up, source_dn

  ! Which layers are cloud-free?  Dummy cloud-free layers are added
  ! above TOA (level 0) and below the ground (level nlev+1).
  logical :: is_cloud_free_layer(0:nlev+1)

  ! Upward and downward overlap matrices - see Hogan et al. (JGR
  ! 2016) for definitions
  real(jprb), dimension(NREGION,NREGION,nlev+1) :: u_overlap, v_overlap

  ! Upwelling and downwelling fluxes at the top and base of each
  ! layer in each region, in W m-2
  real(jprb), dimension(nspec,NREGION,nlev) :: flux_up_base, flux_dn_base
  real(jprb), dimension(nspec,NREGION,nlev) :: flux_up_top, flux_dn_top

  ! Profile of radiances in direction of sensor
  real(jprb), dimension(nspec, nlev+1) :: radiance_profile

  ! Flux up at the surface from radiance, used for specular reflection
  real(jprb), dimension(nspec,NREGION) :: flux_up_surface

  ! Cloud optical depth scaling in each cloudy region
  real(jprb) :: od_scaling(2:NREGION,nlev)

  ! Fractional area coverage of each region
  real(jprb) :: region_fracs(1:NREGION,nlev)

  ! Area of the vertical interface between each pair of regions,
  ! divided by the horizontal area of the domain. For 3 regions there
  ! are two areas: between regions 1 and 2 and between regions 2 and 3
  ! (regions 1 and 3 are assumed not to touch).
  real(jprb) :: region_edge_area(NREGION-1,nlev)

  ! Cloud fractions below this are ignored
  real(jprb), parameter :: cloud_fraction_threshold = 1.0e-6

  ! Is scattering from the surface treated specularly?  Appropriate
  ! for microwave scattering by the sea surface.
  logical :: do_specular_surface_local

  ! Loop indices for region
  integer(jpim) :: jreg

  real(jphook) :: hook_handle

  if (lhook) call dr_hook('tcrad:calc_radiance',0,hook_handle)

  ! By default we do not treat the surface as a specular reflector
  if (present(do_specular_surface)) then
    do_specular_surface_local = do_specular_surface
  else
    do_specular_surface_local = .false.
  endif
  
  ! Compute the wavelength-independent region fractions and
  ! optical-depth scalings
  call calc_region_properties(nlev, cloud_fraction, &
       &  .true., fractional_std, &
       &  region_fracs, &
       &  od_scaling, cloud_fraction_threshold)

  ! Compute wavelength-independent overlap matrices u_overlap and
  ! v_overlap
  call calc_overlap_matrices(nlev, &
       &  region_fracs, overlap_param, &
       &  u_overlap, v_overlap, &
       &  0.5_jprb, &
       &  cloud_fraction_threshold, &
       &  cloud_cover)

  ! Average gas and cloud properties noting that: (1) region 1 is
  ! cloud-free so we copy over the gas optical depth; (2) gases only
  ! absorb so the single scattering albedo (ssa) of region is
  ! implicitly zero and we don't even have an array dimension for
  ! it; (3) since the gases don't scatter, the asymmetry factor of
  ! the gas-cloud mixture is equal to the value for clouds,
  ! regardless of the optical depth scaling (od_scaling) so we
  ! simply use the asymmetry_cloud variable when calculating
  ! reflectance and transmittance.
  od(:,1,:) = od_clear
  do jreg = 2,NREGION
    od(:,jreg,:) = od_clear + od_cloud*spread(od_scaling(jreg,:),1,nspec)
    ssa(:,jreg,:) = ssa_cloud(:,:)*od_cloud(:,:) &
         &  * spread(od_scaling(jreg,:),1,nspec) / od(:,jreg,:)
  end do

  ! Identify cloud-free layers
  is_cloud_free_layer(0) = .true.
  is_cloud_free_layer(1:nlev) = (region_fracs(1,:) == 1.0_jprb)
  is_cloud_free_layer(nlev+1) = .true.

  ! Compute layer-wise properties
  call calc_reflectance_transmittance(nspec, nlev, NREGION, &
       &  region_fracs, planck_hl, od, ssa, asymmetry_cloud, &
       &  reflectance, transmittance, source_up, source_dn)

  ! Classic Tripleclouds method to compute flux profile
  call calc_two_stream_flux(nspec, nlev, surf_emission, surf_albedo, &
       &  reflectance, transmittance, source_up, source_dn, &
       &  is_cloud_free_layer, u_overlap, v_overlap, &
       &  flux_up_base, flux_dn_base, flux_up_top, flux_dn_top)

  ! calc_radiance_up/dn is additive so the radiance profile needs to
  ! be initialized to zero
  radiance_profile = 0.0_jprb

  if (do_specular_surface_local) then
    ! Upward directed radiance measured at top-of-atmosphere, but
    ! first a downward directed radiance which is specularly reflected
    ! from the surface

    ! Compute transmittance and source towards sensor and towards
    ! surface
    call calc_radiance_trans_source(nspec, nlev, NREGION, mu, &
         &  region_fracs, planck_hl, od, ssa, asymmetry_cloud, &
         &  flux_up_base, flux_dn_top, transmittance, &
         &  source_up=source_up, source_dn=source_dn)

    call calc_radiance_dn(nspec, nlev, &
         &  ONE_OVER_PI, transmittance, source_dn, v_overlap, radiance_profile)
    ! Reflect surface downward radiance upwards, spread into the
    ! regions of the lowest layer, and convert to a flux (with PI)
    flux_up_surface = spread(surf_emission + PI*surf_albedo*radiance_profile(:,nlev+1),2,NREGION) &
         &          * spread(region_fracs(:,nlev),1,nspec)
    call calc_radiance_up(nspec, nlev, &
         &  ONE_OVER_PI, flux_up_surface, &
         &  transmittance, source_up, u_overlap, radiance_profile)
    
    ! Extract top-of-atmosphere value
    radiance = radiance_profile(:,1)

  else
    ! Upward directed radiance measured at top-of-atmosphere,
    ! Lambertian surface

    ! Compute transmittance and source towards sensor
    call calc_radiance_trans_source(nspec, nlev, NREGION, mu, &
         &  region_fracs, planck_hl, od, ssa, asymmetry_cloud, &
         &  flux_up_base, flux_dn_top, transmittance, source_up=source_up)
    call calc_radiance_up(nspec, nlev, &
         &  ONE_OVER_PI, flux_up_base(:,:,nlev), &
         &  transmittance, source_up, u_overlap, radiance_profile)
    ! Extract top-of-atmosphere value
    radiance = radiance_profile(:,1)

  end if

  if (lhook) call dr_hook('tcrad:calc_radiance',1,hook_handle)

end subroutine calc_radiance

end module tcrad
