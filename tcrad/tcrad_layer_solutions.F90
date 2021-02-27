! tcrad_layer_solutions.F90 - Two-stream and related layer solutions for TCRAD package
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

module tcrad_layer_solutions

  use parkind1, only : jpim, jprb

  implicit none

  ! Elsasser's factor: the effective factor by which the zenith
  ! optical depth needs to be multiplied to account for longwave
  ! transmission at all angles through the atmosphere.  Alternatively
  ! think of acos(1/LW_DIFFUSIVITY) to be the effective zenith angle
  ! of longwave radiation.
  real(jprb), parameter :: LW_DIFFUSIVITY = 1.66_jprb

  integer(jpim), parameter :: MAX_GAUSS_LEGENDRE_POINTS = 4

contains

  !---------------------------------------------------------------------
  ! Return Gauss-Legendre quadrature points
  subroutine gauss_legendre(npoint, xpoint, weight)

    integer(jpim), intent(in)  :: npoint
    real(jprb),    intent(out) :: xpoint(:), weight(:)

    ! Set default values
    xpoint = 0.5_jprb
    weight = 0.0_jprb

    if (npoint <= 1) then
      xpoint(1) = 0.5_jprb
      weight(1) = 1.0_jprb
    else if (npoint == 2) then
      xpoint(1:2) = [0.211324865405187_jprb, 0.788675134594813_jprb]
      weight(1:2) = [0.5_jprb, 0.5_jprb]
    else if (npoint == 3) then
      xpoint(1:3) = [0.112701665379258_jprb, 0.5_jprb, &
           &         0.887298334620742_jprb]
      weight(1:3) = [0.277777777777777_jprb, 0.444444444444444_jprb, &
           &         0.277777777777777_jprb]
    else
      xpoint(1:4) = [0.0694318442029737_jprb, 0.330009478207572_jprb, &
           &         0.669990521792428_jprb,  0.930568155797026_jprb]
      weight(1:4) = [0.173927422568727_jprb, 0.326072577431273_jprb, &
           &         0.326072577431273_jprb, 0.173927422568727_jprb]
    end if

  end subroutine gauss_legendre


  !---------------------------------------------------------------------
  ! Compute the longwave reflectance and transmittance to diffuse
  ! radiation using the Meador & Weaver (1980) two-stream formulas, as
  ! well as the upward flux at the top and the downward flux at the
  ! base of the layer due to emission from within the layer assuming a
  ! linear variation of Planck function within the layer.
  subroutine calc_reflectance_transmittance(nspec, nlev, nreg, &
       &  region_fracs, planck_hl, od, ssa, asymmetry, &
       &  reflectance, transmittance, source_up, source_dn)
        
    use yomhook,  only           : lhook, dr_hook

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

    real(jprb) :: hook_handle

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
        if (od(jspec,1,jlev) > 1.0e-3_jprb) then
          coeff = LW_DIFFUSIVITY*od(jspec,1,jlev)
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
          coeff = LW_DIFFUSIVITY*od(jspec,1,jlev)
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
            factor = (LW_DIFFUSIVITY * 0.5_jprb) * ssa(jspec,jreg,jlev)
            gamma1 = LW_DIFFUSIVITY - factor*(1.0_jprb + asymmetry(jspec,jlev))
            gamma2 = factor * (1.0_jprb - asymmetry(jspec,jlev))

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
              ! Low optical depth approximation
              reflectance(jspec,jreg,jlev) = gamma2 * od(jspec,jreg,jlev)
              transmittance(jspec,jreg,jlev) = (1.0_jprb - k_exponent*od(jspec,jreg,jlev)) &
                   &  / (1.0_jprb + od(jspec,jreg,jlev)*(gamma1-k_exponent))
              source_up(jspec,jreg,jlev) = (1.0_jprb - reflectance(jspec,jreg,jlev) &
                   &  - transmittance(jspec,jreg,jlev)) &
                   &       * 0.5 * (planck_hl(jspec,jlev) + planck_hl(jspec,jlev+1))
              source_dn(jspec,jreg,jlev) = source_up(jspec,jreg,jlev)
            end if
            ! Scale the sources by the area fraction of the region
            source_up(jspec,jreg,jlev) = region_fracs(jreg,jlev)*source_up(jspec,jreg,jlev)
            source_dn(jspec,jreg,jlev) = region_fracs(jreg,jlev)*source_dn(jspec,jreg,jlev)
          end do
        end do
      end if
    end do

    if (lhook) call dr_hook('tcrad:calc_reflectance_transmittance',1,hook_handle)

  end subroutine calc_reflectance_transmittance


  !---------------------------------------------------------------------
  ! Compute the transmittance to a beam of radiation at a particular
  ! zenith angle cosine (mu), as well as optionally the source from
  ! the layer in that direction up and/or down. The latter includes
  ! both emission, and the source term from scattering using the
  ! incoming flux terms from a previous two-stream calculation.
  subroutine calc_radiance_source(nspec, nlev, nreg, &
       &  mu, region_fracs, planck_hl, od, ssa, asymmetry, &
       &  flux_up_base, flux_dn_base, flux_up_top, flux_dn_top, &
       &  transmittance, source_up, source_dn)
    
    use yomhook,  only           : lhook, dr_hook

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
    real(jprb), intent(in), dimension(nspec,nreg,nlev) :: flux_up_base, flux_dn_base
    real(jprb), intent(in), dimension(nspec,nreg,nlev) :: flux_up_top, flux_dn_top
  
    ! Outputs

    ! Layer transmittance at the requested zenith angle
    real(jprb), intent(out), dimension(nspec,nreg,nlev) :: transmittance

    ! Source term up from the top of the layer or down from its base,
    ! in Watts of power per square metre of the entire gridbox, so the
    ! energy is scaled by the size of each region. Since the user may
    ! only require a radiance up or down, these output arguments are
    ! optional.
    real(jprb), intent(out), dimension(nspec,nreg,nlev), optional &
         &  :: source_up, source_dn

    ! Working variables in W m-2
    real(jprb), dimension(nspec,nreg) :: planck_top, planck_base
    real(jprb), dimension(nspec,nreg) :: source_top, source_base

    ! Other working variables
    real(jprb) :: secant, factor, coeff

    ! Maximum number of active regions in a layer (1 in a cloud-free layer)
    integer(jpim) :: max_reg

    ! Loop indices for level, spectral interval and region
    integer(jpim) :: jlev, jspec, jreg

    real(jprb) :: hook_handle

    if (lhook) call dr_hook('tcrad:calc_radiance_source',0,hook_handle)

    secant = 1.0_jprb / mu

    ! 0.5: half the scattering goes up and half down
    factor = 0.5_jprb * 3.0_jprb * mu / LW_DIFFUSIVITY

    do jlev = 1,nlev

      if (region_fracs(1,jlev) < 1.0_jprb) then
        ! Cloudy layer: scale the Planck terms by the region fraction
        ! and also by the single-scattering co-albedo
        max_reg = nreg
        planck_top(:,1) = planck_hl(:,jlev) * region_fracs(1,jlev)
        planck_top(:,2:nreg) = spread(planck_hl(:,jlev),2,nreg-1) &
             &  * (1.0_jprb - ssa(:,2:nreg,jlev)) &
             &  * spread(region_fracs(2:nreg,jlev),1,nspec)
        planck_base(:,1) = planck_hl(:,jlev+1) * region_fracs(1,jlev)
        planck_base(:,2:nreg) = spread(planck_hl(:,jlev+1),2,nreg-1) &
             &  * (1.0_jprb - ssa(:,2:nreg,jlev)) &
             &  * spread(region_fracs(2:nreg,jlev),1,nspec)
        ! Compute transmittance in all regions
        transmittance(:,1:nreg,jlev) = exp(-od(:,1:nreg,jlev)*secant)
      else
        ! Clear layer
        max_reg = 1
        planck_top(:,1)  = planck_hl(:,jlev)
        planck_base(:,1) = planck_hl(:,jlev+1)
        ! Compute transmittance in active region
        transmittance(:,1,jlev) = exp(-od(:,1,jlev)*secant)
        transmittance(:,2:nreg,jlev) = 1.0_jprb
      end if

      if (present(source_up)) then
        ! Compute the rate of energy emitted or scattered in the
        ! upward direction mu at the top and base of the layer: first
        ! the Planck emission for all regions...
        source_top  = planck_top
        source_base = planck_base
        ! ...then scattering from the scattering source function, but
        ! only in cloudy regions
        if (max_reg > 1) then
          source_top(:,2:nreg) = source_top(:,2:nreg) &
               &  + ssa(:,2:nreg,jlev) &
               &  * (flux_up_top(:,2:nreg,jlev) &
               &     * (0.5_jprb + factor*spread(asymmetry(:,jlev),2,nreg-1)) &
               &    +flux_dn_top(:,2:nreg,jlev) &
               &     * (0.5_jprb - factor*spread(asymmetry(:,jlev),2,nreg-1)))
          source_base(:,2:nreg) = source_base(:,2:nreg) &
               &  + ssa(:,2:nreg,jlev) &
               &  * (flux_up_base(:,2:nreg,jlev) &
               &     * (0.5_jprb + factor*spread(asymmetry(:,jlev),2,nreg-1)) &
               &    +flux_dn_base(:,2:nreg,jlev) &
               &     * (0.5_jprb - factor*spread(asymmetry(:,jlev),2,nreg-1)))
        else
          source_up(:,2:nreg,jlev) = 0.0_jprb
        end if

        ! Compute the energy making it to the top of the layer
        ! accounting for attenuation within it
        do jreg = 1,max_reg
          do jspec = 1,nspec
            coeff = (source_base(jspec,jreg)-source_top(jspec,jreg)) &
                 &   * mu / od(jspec,jreg,jlev)
            source_up(jspec,jreg,jlev) = coeff + source_top(jspec,jreg) &
                 - transmittance(jspec,jreg,jlev) * (coeff + source_base(jspec,jreg))
          end do
        end do
      end if

      if (present(source_dn)) then
        ! Compute the rate of energy emitted or scattered in the
        ! downward direction mu at the top and base of the layer:
        ! first the Planck emission for all regions...
        source_top  = planck_top
        source_base = planck_base
        ! ...then scattering from the scattering source function, but
        ! only in cloudy regions
        if (max_reg > 1) then
          source_top(:,2:nreg) = source_top(:,2:nreg) &
               &  + ssa(:,2:nreg,jlev) &
               &  * (flux_up_top(:,2:nreg,jlev) &
               &     * (0.5_jprb - factor*spread(asymmetry(:,jlev),2,nreg-1)) &
               &    +flux_dn_top(:,2:nreg,jlev) &
               &     * (0.5_jprb + factor*spread(asymmetry(:,jlev),2,nreg-1)))
          source_base(:,2:nreg) = source_base(:,2:nreg) &
               &  + ssa(:,2:nreg,jlev) &
               &  * (flux_up_base(:,2:nreg,jlev) &
               &     * (0.5_jprb - factor*spread(asymmetry(:,jlev),2,nreg-1)) &
               &    +flux_dn_base(:,2:nreg,jlev) &
               &     * (0.5_jprb + factor*spread(asymmetry(:,jlev),2,nreg-1)))
        else
          source_dn(:,2:nreg,jlev) = 0.0_jprb
        end if
        ! Compute the energy making it to the top of the layer
        ! accounting for attenuation within it
        do jreg = 1,max_reg
          do jspec = 1,nspec
            coeff = (source_top(jspec,jreg)-source_base(jspec,jreg)) &
                 &   * mu / od(jspec,jreg,jlev)
            source_dn(jspec,jreg,jlev) = coeff + source_base(jspec,jreg) &
                 - transmittance(jspec,jreg,jlev) * (coeff + source_top(jspec,jreg))
          end do
        end do
      end if

    end do

    if (lhook) call dr_hook('tcrad:calc_radiance_source',1,hook_handle)

  end subroutine calc_radiance_source


  !---------------------------------------------------------------------
  ! Compute the transmittance to a beam of radiation at a particular
  ! zenith angle cosine (mu), as well as optionally the source from
  ! the layer in that direction up and/or down. The latter includes
  ! only emission, so is suitable to be used in a no-scattering
  ! radiance calculation.
  subroutine calc_no_scattering_radiance_source(nspec, nlev, nreg, &
       &  mu, region_fracs, planck_hl, od, &
       &  transmittance, source_up, source_dn)
      
    use yomhook,  only           : lhook, dr_hook

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
  
    ! Outputs

    ! Layer transmittance at the requested zenith angle
    real(jprb), intent(out), dimension(nspec,nreg,nlev) :: transmittance

    ! Source term up from the top of the layer or down from its base,
    ! in Watts of power per square metre of the entire gridbox, so the
    ! energy is scaled by the size of each region. Since the user may
    ! only require a radiance up or down, these output arguments are
    ! optional.
    real(jprb), intent(out), dimension(nspec,nreg,nlev), optional &
         &  :: source_up, source_dn

    ! Working variables in W m-2
    real(jprb) :: source_top, source_base

    ! Other working variables
    real(jprb) :: secant, coeff
   
    ! Maximum number of active regions in a layer (1 in a cloud-free layer)
    integer(jpim) :: max_reg

    ! Loop indices for level, spectral interval and region
    integer(jpim) :: jlev, jspec, jreg

    real(jprb) :: hook_handle

    if (lhook) call dr_hook('tcrad:calc_no_scattering_radiance_source',0,hook_handle)

    secant = 1.0_jprb / mu

    do jlev = 1,nlev

      if (region_fracs(1,jlev) < 1.0_jprb) then
        ! Cloudy layer: scale the Planck terms by the region fraction
        ! and also by the single-scattering co-albedo
        max_reg = nreg
        transmittance(:,1:max_reg,jlev) = exp(-od(:,1:max_reg,jlev)*secant)
      else
        ! Clear layer
        max_reg = 1
        transmittance(:,1,jlev) = exp(-od(:,1,jlev)*secant)
        transmittance(:,2:nreg,jlev) = 1.0_jprb
      end if

      if (present(source_up)) then
        ! Compute upward source from layer top due to Planck emission
        ! within the layer
        do jreg = 1,max_reg
          do jspec = 1,nspec
            source_top  = planck_hl(jspec,jlev)   * region_fracs(jreg,jlev)
            source_base = planck_hl(jspec,jlev+1) * region_fracs(jreg,jlev)
            coeff = (source_base - source_top) &
                 &   * mu / od(jspec,jreg,jlev)
            source_up(jspec,jreg,jlev) = coeff + source_top &
                 - transmittance(jspec,jreg,jlev) * (coeff + source_base)
          end do
        end do
        if (max_reg == 1) then
          ! In a clear layer, set cloudy values to zero
          source_up(:,2:nreg,jlev) = 0.0_jprb
        end if
      end if

      if (present(source_dn)) then
        ! Compute downward source from layer base due to Planck emission
        ! within the layer
        do jreg = 1,max_reg
          do jspec = 1,nspec
            source_top  = planck_hl(jspec,jlev)   * region_fracs(jreg,jlev)
            source_base = planck_hl(jspec,jlev+1) * region_fracs(jreg,jlev)
            coeff = (source_top - source_base) &
                 &   * mu / od(jspec,jreg,jlev)
            source_dn(jspec,jreg,jlev) = coeff + source_base &
                 - transmittance(jspec,jreg,jlev) * (coeff + source_top)
          end do
        end do
        if (max_reg == 1) then
          ! In a clear layer, set cloudy values to zero
          source_dn(:,2:nreg,jlev) = 0.0_jprb
        end if
      end if

    end do

    if (lhook) call dr_hook('tcrad:calc_no_scattering_radiance_source',1,hook_handle)

  end subroutine calc_no_scattering_radiance_source

end module tcrad_layer_solutions
