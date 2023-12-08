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
  public

  ! Two-stream scheme can be Elsasser: radiance field assumed to be
  ! two pencil beams travelling at cosine-zenith angles (mu) of
  ! +/-1/1.66; or Eddington: radiance field assumed to be
  ! L(mu)=L0+mu*L1.
  enum, bind(c)
    enumerator ITwoStreamElsasser, ITwoStreamEddington
  end enum

  ! Two stream scheme currently in use 
  integer(jpim) :: i_two_stream_scheme = ITwoStreamElsasser
  
  ! The effective factor by which the zenith optical depth needs to be
  ! multiplied to account for longwave transmission at all angles
  ! through the atmosphere.  Alternatively think of
  ! acos(1/lw_diffusivity) to be the effective zenith angle of
  ! longwave radiation.
  real(jprb) :: lw_diffusivity = 1.66_jprb ! Elsasser default

  ! To avoid division by near-zero values use simpler formulae in the
  ! low optical depth regime
  real(jprb), parameter :: OD_THRESH_2STREAM = 1.0e-3_jprb
  real(jprb), parameter :: OD_THRESH = 1.0e-3_jprb
!#warning(fixme)
  !real(jprb), parameter :: OD_THRESH = 1.0e-6_jprb

#ifdef PARKIND1_SINGLE
  real(jprb), parameter :: MIN_K_SQUARED = 1.0e-6_jprb
#else
  real(jprb), parameter :: MIN_K_SQUARED = 1.0e-12_jprb
#endif
  
  integer(jpim), parameter :: MAX_GAUSS_LEGENDRE_POINTS = 8

contains

  !---------------------------------------------------------------------
  ! Set the two-stream scheme; this overwrites global module variables
  ! so should normally be called from outside a parallel block
  subroutine set_two_stream_scheme(i_scheme)
    
    integer(jpim), intent(in) :: i_scheme

    ! Only overwrite global module variables if they need changing
    if (i_scheme == ITwoStreamEddington &
         .and. i_two_stream_scheme /= ITwoStreamEddington) then
      ! Toon et al. (1989), Table 1
      i_two_stream_scheme = ITwoStreamEddington
      lw_diffusivity = 2.0_jprb
    else if (i_two_stream_scheme /= ITwoStreamElsasser) then
      ! Elsasser (1942)
      i_two_stream_scheme = ITwoStreamElsasser
      lw_diffusivity = 1.66_jprb
    end if
    
  end subroutine set_two_stream_scheme
  
  
  !---------------------------------------------------------------------
  ! Return Gauss-Legendre quadrature points, or Gauss-Jacobi-5 for
  ! negative values
  subroutine gauss_legendre(npoint, xpoint, weight)

    integer(jpim), intent(in)  :: npoint
    real(jprb),    intent(out) :: xpoint(:), weight(:)

    ! Set default values
    xpoint = 0.5_jprb
    weight = 0.0_jprb

    if (npoint < 0) then
      call gauss_jacobi_5(-npoint, xpoint, weight)
    else if (npoint == 1 .or. npoint == 0) then
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
    else if (npoint == 4) then
      xpoint(1:4) = [0.0694318442029737_jprb, 0.330009478207572_jprb, &
           &         0.669990521792428_jprb,  0.930568155797026_jprb]
      weight(1:4) = [0.173927422568727_jprb, 0.326072577431273_jprb, &
           &         0.326072577431273_jprb, 0.173927422568727_jprb]
    else if (npoint == 5) then
      xpoint(1:5) = [0.9530899230_jprb, 0.7692346551_jprb, &
           &  0.5000000000_jprb, 0.2307653449_jprb, 0.0469100770_jprb]
      weight(1:5) = [0.1184634425_jprb, 0.2393143352_jprb, &
           &  0.2844444444_jprb, 0.2393143352_jprb, 0.1184634425_jprb]
    else if (npoint == 6) then
      xpoint(1:6) = [0.9662347571_jprb, 0.8306046932_jprb, 0.6193095930_jprb, &
           &         0.3806904070_jprb, 0.1693953068_jprb, 0.0337652429_jprb]
      weight(1:6) = [0.0856622462_jprb, 0.1803807865_jprb, 0.2339569673_jprb, &
           &         0.2339569673_jprb, 0.1803807865_jprb, 0.0856622462_jprb]
      
    else if (npoint == 7) then
      xpoint(1:7) = [0.9745539562_jprb, 0.8707655928_jprb, 0.7029225757_jprb, &
           &         0.5000000000_jprb, 0.2970774243_jprb, 0.1292344072_jprb, &
           &         0.0254460438_jprb]
      weight(1:7) = [0.0647424831_jprb, 0.1398526957_jprb, 0.1909150253_jprb, &
           &         0.2089795918_jprb, 0.1909150253_jprb, 0.1398526957_jprb, &
           &         0.0647424831_jprb]
    else
      xpoint(1:8) = [0.9801449282_jprb, 0.8983332387_jprb, 0.7627662050_jprb, &
           &         0.5917173212_jprb, 0.4082826788_jprb, 0.2372337950_jprb, &
           &         0.1016667613_jprb, 0.0198550718_jprb]
      weight(1:8) = [0.0506142681_jprb, 0.1111905172_jprb, 0.1568533229_jprb, &
           &         0.1813418917_jprb, 0.1813418917_jprb, 0.1568533229_jprb,&
           &         0.1111905172_jprb, 0.0506142681_jprb]
    end if

  end subroutine gauss_legendre

  
  !---------------------------------------------------------------------
  ! Return Gauss-Jacobi-5 quadrature points
  subroutine gauss_jacobi_5(npoint, xpoint, weight)

    integer(jpim), intent(in)  :: npoint
    real(jprb),    intent(out) :: xpoint(:), weight(:)

    if (npoint == 1) then
      xpoint(1) = 0.62973761_jprb
      weight(1) = 0.79398148_jprb
    else if (npoint == 2) then
      xpoint(1:2) = [0.25099074_jprb, 0.79084740_jprb]
      weight(1:2) = [0.45823479_jprb, 0.48680354_jprb]
    else if (npoint == 3) then
      xpoint(1:3) = [0.10249222_jprb, 0.44179603_jprb, 0.86337516_jprb]
      weight(1:3) = [0.21358706_jprb, 0.43864096_jprb, 0.32931125_jprb]
    else
      xpoint(1:4) = [0.04545867_jprb, 0.23223344_jprb, &
           &         0.57401988_jprb, 0.90307760_jprb]
      weight(1:4) = [0.10126647_jprb, 0.27681291_jprb, &
           &         0.37658817_jprb, 0.23801084_jprb]
    end if

  end subroutine gauss_jacobi_5


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
            if (i_two_stream_scheme == ITwoStreamElsasser) then
              ! See Fu et al. (1997), Eqs. 2.9 and 2.10
              factor = (lw_diffusivity * 0.5_jprb) * ssa(jspec,jreg,jlev)
              gamma1 = lw_diffusivity - factor*(1.0_jprb + asymmetry(jspec,jlev))
              gamma2 = factor * (1.0_jprb - asymmetry(jspec,jlev))
            else
              ! See Meador & Weaver (1980), Table 1; Toon et al. (1989), Table 1
              gamma1 = 1.75_jprb - ssa(jspec,jreg,jlev)*(1.0_jprb + 0.75_jprb*asymmetry(jspec,jlev))
              gamma2 = ssa(jspec,jreg,jlev)*(1.0_jprb - 0.75_jprb*asymmetry(jspec,jlev)) - 0.25_jprb
            end if

            k_exponent = sqrt(max((gamma1 - gamma2) * (gamma1 + gamma2), &
                 &  MIN_K_SQUARED)) ! Eq 18 of Meador & Weaver (1980)
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
  ! Compute the rates of emission and scattering into a particular
  ! zenith angle cosine (mu) at the top and base of each layer and
  ! region
  subroutine calc_radiance_rates(nspec, nlev, nreg, &
       &  mu, region_fracs, planck_hl, ssa, asymmetry, &
       &  flux_up_base, flux_dn_base, flux_up_top, flux_dn_top, &
       &  rate_up_top, rate_up_base, rate_dn_top, rate_dn_base)
    
    use yomhook,  only           : lhook, dr_hook, jphook

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

    ! Rate of emission/scattering upwards or downwards at the top or
    ! base of each layer and region, in Watts of power per square
    ! metre of the entire gridbox. These terms are available instead
    ! of source_up and source_dn for the 3D radiance calculation. Note
    ! that this is the rate of emission/scattering along the radiance
    ! direction per unit optical depth in that layer region, but since
    ! optical depth is dimensionless, these rates still have units of
    ! W m-2.
    real(jprb), intent(out), dimension(nspec,nreg,nlev), optional &
         &  :: rate_up_top, rate_up_base, rate_dn_top, rate_dn_base

    ! Local variables

    ! Working variables in W m-2
    real(jprb), dimension(nspec,nreg) :: planck_top, planck_base
    real(jprb), dimension(nspec,nreg) :: source_top, source_base

    ! Other working variables
    real(jprb) :: secant, factor, coeff

    ! Maximum number of active regions in a layer (1 in a cloud-free layer)
    integer(jpim) :: max_reg

    ! Loop indices for level, spectral interval and region
    integer(jpim) :: jlev, jspec, jreg

    real(jphook) :: hook_handle

    if (lhook) call dr_hook('tcrad:calc_radiance_rates',0,hook_handle)

    secant = 1.0_jprb / mu

    ! 0.5: half the scattering goes up and half down
    factor = 0.5_jprb * 3.0_jprb * mu / lw_diffusivity

    do jlev = 1,nlev

      ! CHECK / FIX: is the source that per unit optical depth in the
      ! vertical or along the direction to the sensor?

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
      else
        ! Clear layer
        max_reg = 1
        planck_top(:,1)  = planck_hl(:,jlev)
        planck_base(:,1) = planck_hl(:,jlev+1)
      end if

      if (present(rate_up_top) .and. present(rate_up_base)) then
        ! Compute the rate of energy emitted or scattered in the
        ! upward direction mu at the top and base of the layer: first
        ! the Planck emission for all regions...
        rate_up_top(:,:,jlev)  = planck_top
        rate_up_base(:,:,jlev) = planck_base
        ! ...then scattering from the scattering source function, but
        ! only in cloudy regions
        if (max_reg > 1) then
          rate_up_top(:,2:nreg,jlev) = rate_up_top(:,2:nreg,jlev) &
               &  + ssa(:,2:nreg,jlev) &
               &  * (flux_up_top(:,2:nreg,jlev) &
               &     * (0.5_jprb + factor*spread(asymmetry(:,jlev),2,nreg-1)) &
               &    +flux_dn_top(:,2:nreg,jlev) &
               &     * (0.5_jprb - factor*spread(asymmetry(:,jlev),2,nreg-1)))
          rate_up_base(:,2:nreg,jlev) = rate_up_base(:,2:nreg,jlev) &
               &  + ssa(:,2:nreg,jlev) &
               &  * (flux_up_base(:,2:nreg,jlev) &
               &     * (0.5_jprb + factor*spread(asymmetry(:,jlev),2,nreg-1)) &
               &    +flux_dn_base(:,2:nreg,jlev) &
               &     * (0.5_jprb - factor*spread(asymmetry(:,jlev),2,nreg-1)))
        end if
      end if

      if (present(rate_dn_top) .and. present(rate_dn_base)) then
        ! Compute the rate of energy emitted or scattered in the
        ! downward direction mu at the top and base of the layer:
        ! first the Planck emission for all regions...
        rate_dn_top(:,:,jlev)  = planck_top
        rate_dn_base(:,:,jlev) = planck_base
        ! ...then scattering from the scattering source function, but
        ! only in cloudy regions
        if (max_reg > 1) then
          rate_dn_top(:,2:nreg,jlev) = rate_dn_top(:,2:nreg,jlev) &
               &  + ssa(:,2:nreg,jlev) &
               &  * (flux_up_top(:,2:nreg,jlev) &
               &     * (0.5_jprb - factor*spread(asymmetry(:,jlev),2,nreg-1)) &
               &    +flux_dn_top(:,2:nreg,jlev) &
               &     * (0.5_jprb + factor*spread(asymmetry(:,jlev),2,nreg-1)))
          rate_dn_base(:,2:nreg,jlev) = rate_dn_base(:,2:nreg,jlev) &
               &  + ssa(:,2:nreg,jlev) &
               &  * (flux_up_base(:,2:nreg,jlev) &
               &     * (0.5_jprb - factor*spread(asymmetry(:,jlev),2,nreg-1)) &
               &    +flux_dn_base(:,2:nreg,jlev) &
               &     * (0.5_jprb + factor*spread(asymmetry(:,jlev),2,nreg-1)))
        end if
      end if

    end do

    if (lhook) call dr_hook('tcrad:calc_radiance_rates',1,hook_handle)

  end subroutine calc_radiance_rates

  
  !---------------------------------------------------------------------
  ! Compute the transmittance to a beam of radiation at a particular
  ! zenith angle cosine (mu), as well as optionally the source from
  ! the layer in that direction up and/or down. The latter includes
  ! only emission, so is suitable to be used in a no-scattering
  ! radiance calculation.
  subroutine calc_no_scattering_radiance_source(nspec, nlev, nreg, &
       &  mu, region_fracs, planck_hl, od, &
       &  transmittance, source_up, source_dn)
      
    use yomhook,  only           : lhook, dr_hook, jphook

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

    real(jphook) :: hook_handle

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
            if (od(jspec,jreg,jlev) > OD_THRESH) then
              coeff = (source_base - source_top) &
                   &   * mu / od(jspec,jreg,jlev)
              source_up(jspec,jreg,jlev) = coeff + source_top &
                   - transmittance(jspec,jreg,jlev) * (coeff + source_base)
            else
              source_up(jspec,jreg,jlev) = od(jspec,jreg,jlev) &
                   &  * 0.5_jprb * (source_top+source_base) / mu
            end if
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
            if (od(jspec,jreg,jlev) > OD_THRESH) then
              coeff = (source_top - source_base) &
                   &   * mu / od(jspec,jreg,jlev)
              source_dn(jspec,jreg,jlev) = coeff + source_base &
                   - transmittance(jspec,jreg,jlev) * (coeff + source_top)
            else
              source_dn(jspec,jreg,jlev) = od(jspec,jreg,jlev) &
                   &  * 0.5_jprb * (source_top+source_base) / mu
            end if
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


  !---------------------------------------------------------------------
  ! Compute the clear-sky transmittance to a beam of radiation at a
  ! particular zenith angle cosine (mu), as well as optionally the
  ! source from the layer in that direction up and/or down. The latter
  ! includes only emission, so is suitable to be used in a clear-sky
  ! no-scattering radiance calculation.
  subroutine calc_clear_sky_trans_source(nspec, nlev, &
       &  mu, planck_hl, od, &
       &  transmittance, source_up, source_dn)
      
    use yomhook,  only           : lhook, dr_hook, jphook

    ! Inputs

    ! Number of spectral intervals and levels
    integer(jpim), intent(in) :: nspec, nlev

    ! Cosine of the zenith angle (positive)
    real(jprb) :: mu

    ! Planck function integrated over each spectral interval at each
    ! half-level, in W m-2 (i.e. the flux emitted by a horizontal
    ! black-body surface)
    real(jprb), intent(in), dimension(nspec,nlev+1) :: planck_hl

    ! Optical depth in each layer
    real(jprb), intent(in), dimension(nspec,nlev) :: od
  
    ! Outputs

    ! Layer transmittance at the requested zenith angle
    real(jprb), intent(out), dimension(nspec,nlev) :: transmittance

    ! Source term up from the top of the layer or down from its base,
    ! in Watts of power per square metre of the entire gridbox, so the
    ! energy is scaled by the size of each region. Since the user may
    ! only require a radiance up or down, these output arguments are
    ! optional.
    real(jprb), intent(out), dimension(nspec,nlev), optional &
         &  :: source_up, source_dn

    ! Working variables
    real(jprb) :: secant, coeff
   
    ! Loop indices for level and spectral interval
    integer(jpim) :: jlev, jspec

    real(jphook) :: hook_handle

    if (lhook) call dr_hook('tcrad:calc_clear_sky_trans_source',0,hook_handle)

    secant = 1.0_jprb / mu

    do jlev = 1,nlev

      transmittance(:,jlev) = exp(-od(:,jlev)*secant)

      if (present(source_up)) then
        ! Compute upward source from layer top due to Planck emission
        ! within the layer
        do jspec = 1,nspec
          if (od(jspec,jlev) > OD_THRESH) then
            coeff = (planck_hl(jspec,jlev+1) - planck_hl(jspec,jlev)) &
                 &   * mu / od(jspec,jlev)
            source_up(jspec,jlev) = coeff + planck_hl(jspec,jlev) &
                 - transmittance(jspec,jlev) * (coeff + planck_hl(jspec,jlev+1))
          else
            source_up(jspec,jlev) = od(jspec,jlev) &
                 &  * 0.5_jprb * (planck_hl(jspec,jlev)+planck_hl(jspec,jlev+1)) / mu
          end if
        end do
      end if

      if (present(source_dn)) then
        ! Compute downward source from layer base due to Planck emission
        ! within the layer
        do jspec = 1,nspec
          if (od(jspec,jlev) > OD_THRESH) then
            coeff = (planck_hl(jspec,jlev) - planck_hl(jspec,jlev+1)) &
                 &   * mu / od(jspec,jlev)
            source_dn(jspec,jlev) = coeff + planck_hl(jspec,jlev+1) &
                 - transmittance(jspec,jlev) * (coeff + planck_hl(jspec,jlev))
          else
            source_dn(jspec,jlev) = od(jspec,jlev) &
                 &  * 0.5_jprb * (planck_hl(jspec,jlev)+planck_hl(jspec,jlev+1)) / mu
          end if
        end do
      end if

    end do

    if (lhook) call dr_hook('tcrad:calc_clear_sky_trans_source',1,hook_handle)

  end subroutine calc_clear_sky_trans_source


  !---------------------------------------------------------------------
  ! Calculate the transmittance of each layer and region along a path
  ! with consine of zenith angle "mu", as well as (optionally) the
  ! emission up from the top of the layer and down through its base,
  ! computed neglecting 3D effects
  subroutine calc_radiance_trans_source(nspec, nlev, nreg, &
             &  mu, region_fracs, od, transmittance, &
             &  rate_up_top, rate_up_base, rate_dn_top, rate_dn_base, &
             &  source_up, source_dn)

    use yomhook,  only           : lhook, dr_hook, jphook

    ! Inputs

    ! Number of spectral intervals, levels and regions
    integer(jpim), intent(in) :: nspec, nlev, nreg

    ! Cosine of the zenith angle (positive)
    real(jprb) :: mu

    ! Fraction of the gridbox occupied by each region (summing to 1)
    ! at each level
    real(jprb), intent(in), dimension(nreg,nlev) :: region_fracs

    ! Optical depth in each region and layer
    real(jprb), intent(in), dimension(nspec,nreg,nlev) :: od
  
    ! Outputs

    ! Layer transmittance at the requested zenith angle
    real(jprb), intent(out), dimension(nspec,nreg,nlev) :: transmittance

    ! Optional inputs

    ! Rate of emission/scattering upwards or downwards at the top or
    ! base of each layer and region, in Watts of power per square
    ! metre of the entire gridbox. These terms are available instead
    ! of source_up and source_dn for the 3D radiance calculation. Note
    ! that this is the rate of emission/scattering along the radiance
    ! direction per unit optical depth in that layer region, but since
    ! optical depth is dimensionless, these rates still have units of
    ! W m-2.
    real(jprb), intent(in), dimension(nspec,nreg,nlev), optional &
         &  :: rate_up_top, rate_up_base, rate_dn_top, rate_dn_base

    ! Optional outputs

    ! Source term up from the top of the layer or down from its base,
    ! in Watts of power per square metre of the entire gridbox, so the
    ! energy is scaled by the size of each region. Since the user may
    ! only require a radiance up or down, these output arguments are
    ! optional.
    real(jprb), intent(out), dimension(nspec,nreg,nlev), optional &
         &  :: source_up, source_dn

    ! Local variables

    ! Working variables
    real(jprb) :: secant, coeff

    ! Maximum number of active regions in layer
    integer(jpim) :: max_reg

    ! Loop indices
    integer(jpim) :: jlev, jreg, jspec

    real(jphook) :: hook_handle

    if (lhook) call dr_hook('tcrad:calc_radiance_trans_source',0,hook_handle)

    secant = 1.0_jprb / mu

    do jlev = 1,nlev

      if (region_fracs(1,jlev) < 1.0_jprb) then
        max_reg = nreg;
      else
        max_reg = 1;
      end if

      do jreg = 1,max_reg
        transmittance(:,jreg,jlev) = exp(-od(:,jreg,jlev)*secant)

        if (present(source_dn)) then
          do jspec = 1,nspec
            if (od(jspec,jreg,jlev) > OD_THRESH) then
              coeff = (rate_dn_top(jspec,jreg,jlev) - rate_dn_base(jspec,jreg,jlev)) &
                   &   * mu / od(jspec,jreg,jlev)
              source_dn(jspec,jreg,jlev) = coeff + rate_dn_base(jspec,jreg,jlev) &
                   - transmittance(jspec,jreg,jlev) * (coeff + rate_dn_top(jspec,jreg,jlev))
            else
              source_dn(jspec,jreg,jlev) = od(jspec,jreg,jlev) * 0.5_jprb &
                   &  * (rate_dn_top(jspec,jreg,jlev) + rate_dn_base(jspec,jreg,jlev)) / mu
            end if
          end do
        end if
        if (present(source_up)) then
          do jspec = 1,nspec
            if (od(jspec,jreg,jlev) > OD_THRESH) then
              coeff = (rate_up_base(jspec,jreg,jlev) - rate_up_top(jspec,jreg,jlev)) &
                   &   * mu / od(jspec,jreg,jlev)
              source_up(jspec,jreg,jlev) = coeff + rate_up_top(jspec,jreg,jlev) &
                   - transmittance(jspec,jreg,jlev) * (coeff + rate_up_base(jspec,jreg,jlev))
            else
              source_up(jspec,jreg,jlev) = od(jspec,jreg,jlev) * 0.5_jprb &
                   &  * (rate_up_top(jspec,jreg,jlev) + rate_up_base(jspec,jreg,jlev)) / mu
            end if
          end do
        end if
      end do

      if (max_reg == 1) then
        transmittance(:,2:,jlev) = 1.0_jprb
        if (present(source_dn)) then
          source_dn(:,2:,jlev) = 0.0_jprb
        end if
        if (present(source_up)) then
          source_up(:,2:,jlev) = 0.0_jprb
        end if
      end if
      
    end do

    if (lhook) call dr_hook('tcrad:calc_radiance_trans_source',1,hook_handle)

  end subroutine calc_radiance_trans_source


!#warning "calc_radiance_source is deprecated"

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
    
    use yomhook,  only           : lhook, dr_hook, jphook

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

    real(jphook) :: hook_handle

    if (lhook) call dr_hook('tcrad:calc_radiance_source',0,hook_handle)

    secant = 1.0_jprb / mu

    ! 0.5: half the scattering goes up and half down
    factor = 0.5_jprb * 3.0_jprb * mu / lw_diffusivity

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
            if (od(jspec,jreg,jlev) > OD_THRESH) then
              coeff = (source_base(jspec,jreg)-source_top(jspec,jreg)) &
                   &   * mu / od(jspec,jreg,jlev)
              source_up(jspec,jreg,jlev) = coeff + source_top(jspec,jreg) &
                   - transmittance(jspec,jreg,jlev) * (coeff + source_base(jspec,jreg))
            else
              source_up(jspec,jreg,jlev) = od(jspec,jreg,jlev) * 0.5_jprb &
                   &  * (source_base(jspec,jreg)+source_top(jspec,jreg)) / mu
            end if
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
            if (od(jspec,jreg,jlev) > OD_THRESH) then
              coeff = (source_top(jspec,jreg)-source_base(jspec,jreg)) &
                   &   * mu / od(jspec,jreg,jlev)
              source_dn(jspec,jreg,jlev) = coeff + source_base(jspec,jreg) &
                   - transmittance(jspec,jreg,jlev) * (coeff + source_top(jspec,jreg))
            else
              source_dn(jspec,jreg,jlev) = od(jspec,jreg,jlev) * 0.5_jprb &
                   &  * (source_base(jspec,jreg)+source_top(jspec,jreg)) / mu
            end if
          end do
        end do
      end if

    end do

    if (lhook) call dr_hook('tcrad:calc_radiance_source',1,hook_handle)

  end subroutine calc_radiance_source


  !---------------------------------------------------------------------
  ! Calculate the transmittance of each layer and region along a path
  ! with consine of zenith angle "mu", as well as (optionally) the
  ! emission up from the top of the layer and down through its base,
  ! computed neglecting 3D effects. This version uses an exact
  ! solution to the two-stream equations, rather than the
  ! calc_radiance_rates/calc_radiance_trans_source combination, which
  ! assume a linear variation of fluxes with optical depth.
  subroutine calc_radiance_trans_source_exact(nspec, nlev, nreg, &
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

    if (lhook) call dr_hook('tcrad:calc_radiance_trans_source_exact',0,hook_handle)

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
             &  * (1.0_jprb - 0.0_jprb*ssa(:,2:nreg,jlev)) &
             &  * spread(region_fracs(2:nreg,jlev),1,nspec)
        planck_base(:,1) = planck_hl(:,jlev+1) * region_fracs(1,jlev)
        planck_base(:,2:nreg) = spread(planck_hl(:,jlev+1),2,nreg-1) &
             &  * (1.0_jprb - 0.0_jprb*ssa(:,2:nreg,jlev)) &
             &  * spread(region_fracs(2:nreg,jlev),1,nspec)
      else
        ! Clear layer
        max_reg = 1
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
            if (i_two_stream_scheme == ITwoStreamElsasser) then
              ! See Fu et al. (1997), Eqs. 2.9 and 2.10
              factor = (lw_diffusivity * 0.5_jprb) * ssa(jspec,jreg,jlev)
              gamma1 = lw_diffusivity - factor*(1.0_jprb + asymmetry(jspec,jlev))
              gamma2 = factor * (1.0_jprb - asymmetry(jspec,jlev))
            else
              ! See Meador & Weaver (1980), Table 1; Toon et al. (1989), Table 1
              gamma1 = 1.75_jprb - ssa(jspec,jreg,jlev)*(1.0_jprb + 0.75_jprb*asymmetry(jspec,jlev))
              gamma2 = ssa(jspec,jreg,jlev)*(1.0_jprb - 0.75_jprb*asymmetry(jspec,jlev)) - 0.25_jprb
            end if
            k_exponent = sqrt(max((gamma1 - gamma2) * (gamma1 + gamma2), &
                 &  MIN_K_SQUARED)) ! Eq 18 of Meador & Weaver (1980)
              
            ! Phase functions from upwelling flux to upwelling radiance (or down to down)
            p_same     = 1.0_jprb + 3.0_jprb * asymmetry(jspec,jlev) * mu / lw_diffusivity
            ! Phase function from downwelling flux to upwelling radiance (or up to down)
            p_opposite = 1.0_jprb - 3.0_jprb * asymmetry(jspec,jlev) * mu / lw_diffusivity

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
              source_up(jspec,1,jlev) = ssa(jspec,jreg,jlev)*od(jspec,jreg,jlev) &
                   &  * 0.5_jprb*(planck_base(jspec,jreg)+planck_top(jspec,jreg)) / mu
            end if
            if (present(source_dn)) then
              source_dn(jspec,1,jlev) = ssa(jspec,jreg,jlev)*od(jspec,jreg,jlev) &
                   &  * 0.5_jprb*(planck_base(jspec,jreg)+planck_top(jspec,jreg)) / mu
            end if
          end if
          
        end do ! jspec
      end do ! jreg
    end do ! jlev

    if (lhook) call dr_hook('tcrad:calc_radiance_trans_source_exact',1,hook_handle)

  end subroutine calc_radiance_trans_source_exact

end module tcrad_layer_solutions
