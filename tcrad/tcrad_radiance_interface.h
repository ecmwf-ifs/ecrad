! tcrad_radiance_interface.F90 - Interface routines for TCRAD radiances -*- f90 -*-
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
!
! This file is included in modules specifying the NREGION parameter
! (typically 2 or 3) which makes this routine use a "doubleclouds" or
! "tripleclouds" assumption.
!

!---------------------------------------------------------------------
! Compute the TOA or surface radiance including the effects of scattering
subroutine calc_radiance(nspec, nlev, surf_emission, surf_albedo, planck_hl, &
     &  cloud_fraction, &
#if NUM_REGIONS == 3
     &  fractional_std, &
#endif
     &  od_clear, od_cloud, ssa_cloud, asymmetry_cloud, &
     &  overlap_param, mu, radiance, cloud_cover, &
     &  layer_thickness, inv_cloud_scale, do_specular_surface)

  use parkind1, only           : jpim, jprb
  use yomhook,  only           : lhook, dr_hook
  use tcrad_layer_solutions, only   : calc_reflectance_transmittance, &
       &  calc_radiance_rates, calc_radiance_trans_source, LW_DIFFUSIVITY

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

#if NUM_REGIONS == 3
  ! Profile of the fractional standard deviation (i.e. standard
  ! deviation divided by mean) of the horizontal in-cloud water
  ! content (or cloud extinction coefficient) distribution.
  real(jprb), intent(in), dimension(nlev) :: fractional_std
#endif

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

  ! If 3D effects are to be simulated we need the layer thickness in
  ! metres...
  real(jprb), intent(in), optional :: layer_thickness(nlev)

  ! ...and the cloud horizontal scale in metres, where we use the
  ! cloud separation scale defined by Fielding et al. (QJRMS 2020)
  real(jprb), intent(in), optional :: inv_cloud_scale(nlev)

  ! Is scattering from the surface treated specularly?  Appropriate
  ! for microwave scattering by the sea surface.
  logical, intent(in), optional :: do_specular_surface

  ! Local variables

  ! Combined gas/aerosol/cloud optical depth in each region
  real(jprb), dimension(nspec,NREGION,nlev)   :: od

  ! Single scattering albedo of the cloudy regions (ssa=0 in the
  ! clear region)
  real(jprb), dimension(nspec,2:NREGION,nlev) :: ssa

  ! Reflectance and transmittance of each layer and region
  real(jprb), dimension(nspec,NREGION,nlev) :: reflectance, transmittance

  ! Transmission matrix for 3D effects
  real(jprb), dimension(nspec,NREGION,NREGION,nlev) :: transmittance_mat

  ! Rate of emission up from the top or down through the base of
  ! each layer and region (W m-2)
  real(jprb), dimension(nspec,NREGION,nlev) :: source_up, source_dn

  ! Rate of emission/scattering in the direction of a particular
  ! radiance at the top and base of each layer and region, per unit
  ! optical depth (W m-2), used for 3D radiances
  real(jprb), dimension(nspec,NREGION,nlev) :: rate_up_top, rate_up_base
  real(jprb), dimension(nspec,NREGION,nlev) :: rate_dn_top, rate_dn_base

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

  ! Do we represent 3D effects in the radiance calculations?
  logical :: do_3d_effects

  ! Is scattering from the surface treated specularly?  Appropriate
  ! for microwave scattering by the sea surface.
  logical :: do_specular_surface_local

  ! Loop indices for region
  integer(jpim) :: jreg

  real(jprb) :: hook_handle

  if (lhook) call dr_hook('tcrad:calc_radiance',0,hook_handle)

  if (present(layer_thickness) .and. present(inv_cloud_scale)) then
    do_3d_effects = .true.
  else
    do_3d_effects = .false.
  end if

  if (present(do_specular_surface)) then
    do_specular_surface_local = do_specular_surface
  else
    do_specular_surface_local = .false.
  end if

  ! Compute the wavelength-independent region fractions and
  ! optical-depth scalings
  call calc_region_properties(nlev, cloud_fraction, &
#if NUM_REGIONS == 3
       &  .true., fractional_std, &
#endif
       &  region_fracs, &
       &  od_scaling, cloud_fraction_threshold)

  if (do_3d_effects) then
    call calc_region_edge_areas(nlev, region_fracs, layer_thickness, &
         &                      inv_cloud_scale, region_edge_area)
  end if

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

  if (mu >= 0.0_jprb .and. do_specular_surface_local) then
    ! Upward directed radiance measured at top-of-atmosphere, but
    ! first a downward directed radiance which is specularly reflected
    ! from the surface

    ! Compute transmittance and source towards sensor and towards
    ! surface
    call calc_radiance_rates(nspec, nlev, NREGION, mu, &
         &  region_fracs, planck_hl, ssa, asymmetry_cloud, &
         &  flux_up_base, flux_dn_base, flux_up_top, flux_dn_top, &
         &  rate_up_top=rate_up_top, rate_up_base=rate_up_base, &
         &  rate_dn_top=rate_dn_top, rate_dn_base=rate_dn_base)

    ! Compute surface radiance excluding 3D effects (not worth
    ! considering this detail in the downward pass)
    if (do_3d_effects) then
      call calc_radiance_trans_source_3d(nspec, nlev, &
           &  mu, region_fracs, region_edge_area, od, &
           &  transmittance_mat, &
           &  rate_up_top=rate_up_top, rate_up_base=rate_up_base, &
           &  rate_dn_top=rate_dn_top, rate_dn_base=rate_dn_base, &
           &  source_up=source_up, source_dn=source_dn)
      call calc_radiance_dn_3d(nspec, nlev, ONE_OVER_PI, &
           &  transmittance_mat, source_dn, &
           &  v_overlap, radiance_profile)
      ! Reflect surface downward radiance upwards, spread into the
      ! regions of the lowest layer, and convert to a flux (with PI)
      flux_up_surface = spread(surf_emission + PI*surf_albedo*radiance_profile(:,nlev+1),2,NREGION) &
           &          * spread(region_fracs(:,nlev),1,nspec)
      call calc_radiance_up_3d(nspec, nlev, ONE_OVER_PI, &
           &  flux_up_surface, &
           &  transmittance_mat, source_up, &
           &  u_overlap, radiance_profile)
    else
      call calc_radiance_trans_source(nspec, nlev, NREGION, mu, &
           &  region_fracs, od, transmittance, &
           &  rate_up_top=rate_up_top, rate_up_base=rate_up_base, &
           &  rate_dn_top=rate_dn_top, rate_dn_base=rate_dn_base, &
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
    end if
    ! Extract top-of-atmosphere value
    radiance = radiance_profile(:,1)

  else if (mu >= 0.0_jprb) then
    ! Upward directed radiance measured at top-of-atmosphere,
    ! Lambertian surface

    ! Compute transmittance and source towards sensor
    call calc_radiance_rates(nspec, nlev, NREGION, mu, &
         &  region_fracs, planck_hl, ssa, asymmetry_cloud, &
         &  flux_up_base, flux_dn_base, flux_up_top, flux_dn_top, &
         &  rate_up_top=rate_up_top, rate_up_base=rate_up_base)

    if (do_3d_effects) then
      call calc_radiance_trans_source_3d(nspec, nlev, &
           &  mu, region_fracs, region_edge_area, od, &
           &  transmittance_mat, &
           &  rate_up_top=rate_up_top, rate_up_base=rate_up_base, &
           &  source_up=source_up)
      call calc_radiance_up_3d(nspec, nlev, ONE_OVER_PI, &
           &  flux_up_base(:,:,nlev), &
           &  transmittance_mat, source_up, &
           &  u_overlap, radiance_profile)
    else
      call calc_radiance_trans_source(nspec, nlev, NREGION, mu, &
           &  region_fracs, od, transmittance, &
           &  rate_up_top=rate_up_top, rate_up_base=rate_up_base, &
           &  source_up=source_up)
      call calc_radiance_up(nspec, nlev, &
           &  ONE_OVER_PI, flux_up_base(:,:,nlev), &
           &  transmittance, source_up, u_overlap, radiance_profile)
    end if
    ! Extract top-of-atmosphere value
    radiance = radiance_profile(:,1)

  else
    ! Downward directed radiance measured at the surface

    ! Compute transmittance and source towards sensor
    call calc_radiance_rates(nspec, nlev, NREGION, -mu, &
         &  region_fracs, planck_hl, ssa, asymmetry_cloud, &
         &  flux_up_base, flux_dn_base, flux_up_top, flux_dn_top, &
         &  rate_dn_top=rate_dn_top, rate_dn_base=rate_dn_base)

    if (do_3d_effects) then
      call calc_radiance_trans_source_3d(nspec, nlev, &
           &  -mu, region_fracs, region_edge_area, od, &
           &  transmittance_mat, &
           &  rate_dn_top=rate_dn_top, rate_dn_base=rate_dn_base, &
           &  source_dn=source_dn)
      call calc_radiance_dn_3d(nspec, nlev, ONE_OVER_PI, &
           &  transmittance_mat, source_dn, &
           &  v_overlap, radiance_profile)
    else
      call calc_radiance_trans_source(nspec, nlev, NREGION, mu, &
           &  region_fracs, od, transmittance, &
           &  rate_dn_top=rate_dn_top, rate_dn_base=rate_dn_base, &
           &  source_dn=source_dn)
      call calc_radiance_dn(nspec, nlev, &
           &  ONE_OVER_PI, transmittance, source_dn, v_overlap, radiance_profile)
    end if
    ! Extract surface value
    radiance = radiance_profile(:,nlev+1)

  end if

  if (lhook) call dr_hook('tcrad:calc_radiance',1,hook_handle)

end subroutine calc_radiance


!---------------------------------------------------------------------
! Compute the TOA or surface radiance neglecting the effects of
! scattering
subroutine calc_no_scattering_radiance(nspec, nlev, surf_emission, surf_albedo, &
     &  planck_hl, cloud_fraction, &
#if NUM_REGIONS == 3
     &  fractional_std, &
#endif
     &  od_clear, od_cloud, overlap_param, mu, radiance, do_3d_effects, cloud_cover)

  use parkind1, only           : jpim, jprb
  use yomhook,  only           : lhook, dr_hook
  use tcrad_layer_solutions, only : calc_no_scattering_radiance_source, LW_DIFFUSIVITY

  implicit none

  real(jprb), parameter :: ONE_OVER_PI = 1.0_jprb / acos(-1.0_jprb)

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

#if NUM_REGIONS == 3
  ! Profile of the fractional standard deviation (i.e. standard
  ! deviation divided by mean) of the horizontal in-cloud water
  ! content (or cloud extinction coefficient) distribution.
  real(jprb), intent(in), dimension(nlev) :: fractional_std
#endif

  ! Layer optical depth of gas and aerosol
  real(jprb), intent(in), dimension(nspec,nlev) :: od_clear

  ! Layer optical depth of cloud averaged only over the cloudy part
  ! of the gridbox. If Chou scaling is required then this should
  ! have been done already.
  real(jprb), intent(in), dimension(nspec,nlev) :: od_cloud

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

  ! Optional inputs

  ! Do we represent 3D effects in the radiance calculations?
  logical,    intent(in), optional :: do_3d_effects

  ! Return cloud cover computed from cloud fraction profile and
  ! overlap rules
  real(jprb), intent(out), optional :: cloud_cover

  ! Local variables

  ! Combined gas/aerosol/cloud optical depth in each region
  real(jprb), dimension(nspec,NREGION,nlev)   :: od

  ! Transmittance of each layer and region
  real(jprb), dimension(nspec,NREGION,nlev) :: transmittance

  ! Rate of emission up from the top or down through the base of
  ! each layer and region (W m-2)
  real(jprb), dimension(nspec,NREGION,nlev) :: source_up, source_dn

  ! Which layers are cloud-free?  Dummy cloud-free layers are added
  ! above TOA (level 0) and below the ground (level nlev+1).
  logical :: is_cloud_free_layer(0:nlev+1)

  ! Upward and downward overlap matrices - see Hogan et al. (JGR
  ! 2016) for definitions
  real(jprb), dimension(NREGION,NREGION,nlev+1) :: u_overlap, v_overlap

  ! Surface upwelling flux (W m-2)
  real(jprb), dimension(nspec,NREGION) :: flux_up_surf

  ! Profile of radiances in direction of sensor
  real(jprb), dimension(nspec, nlev+1) :: radiance_profile

  ! Cloud optical depth scaling in each cloudy region
  real(jprb) :: od_scaling(2:NREGION,nlev)

  ! Fractional area coverage of each region
  real(jprb) :: region_fracs(1:NREGION,nlev)

  ! Cloud fractions below this are ignored
  real(jprb), parameter :: cloud_fraction_threshold = 1.0e-6

  ! Local versions of optional arguments
  logical :: do_3d_effects_local

  ! Loop indices for region
  integer(jpim) :: jreg

  real(jprb) :: hook_handle

  if (lhook) call dr_hook('tcrad:calc_no_scattering_radiance',0,hook_handle)

  ! Store local values for optional variables
  if (present(do_3d_effects)) then
    do_3d_effects_local = do_3d_effects
  else
    do_3d_effects_local = .false.
  end if

  ! Compute the wavelength-independent region fractions and
  ! optical-depth scalings
  call calc_region_properties(nlev, cloud_fraction, &
#if NUM_REGIONS == 3
       &  .true., fractional_std, &
#endif
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

  ! Average gas and cloud properties
  od(:,1,:) = od_clear
  do jreg = 2,NREGION
    od(:,jreg,:) = od_clear + od_cloud*spread(od_scaling(jreg,:),1,nspec)
  end do

  ! Identify cloud-free layers
  is_cloud_free_layer(0) = .true.
  is_cloud_free_layer(1:nlev) = (region_fracs(1,:) == 1.0_jprb)
  is_cloud_free_layer(nlev+1) = .true.

  ! calc_radiance_up/dn is additive so the radiance profile needs to
  ! be initialized to zero
  radiance_profile = 0.0_jprb

  if (mu >= 0.0_jprb) then
    ! Upward directed radiance measured at top-of-atmosphere

    ! In order to obtain the surface reflectance we first calculate
    ! the downward flux profile

    ! Compute transmittance and source towards sensor
    call calc_no_scattering_radiance_source(nspec, nlev, NREGION, &
         &  LW_DIFFUSIVITY, region_fracs, planck_hl, od, &
         &  transmittance, source_dn=source_dn)

    ! Downward radiances - in fact these are scaled to be fluxes
    call calc_radiance_dn(nspec, nlev, 1.0_jprb, &
         &  transmittance, source_dn, v_overlap, radiance_profile)

    ! Surface upwelling flux: note that the reflection term has no
    ! memory of the location of clouds above, which is "anomalous
    ! horizontal transport" in the terminology of Shonk & Hogan
    ! (2008), although in the longwave this effect ought to be small.
    flux_up_surf = spread(surf_emission + radiance_profile(:,nlev+1)*surf_albedo, &
         &                2, NREGION) &
         &       * spread(region_fracs(:,nlev), 1, nspec)

    ! Compute transmittance and source towards sensor
    call calc_no_scattering_radiance_source(nspec, nlev, NREGION, mu, &
         &  region_fracs, planck_hl, od, &
         &  transmittance, source_up=source_up)
    
    ! Compute radiance profile
    call calc_radiance_up(nspec, nlev, &
         &  ONE_OVER_PI, flux_up_surf, &
         &  transmittance, source_up, u_overlap, radiance_profile)

    ! Extract top-of-atmosphere value
    radiance = radiance_profile(:,1)

  else
    ! Downward directed radiance measured at the surface

    ! Compute transmittance and source towards sensor
    call calc_no_scattering_radiance_source(nspec, nlev, NREGION, -mu, &
         &  region_fracs, planck_hl, od, &
         &  transmittance, source_dn=source_dn)
    
    ! Compute radiance profile
    call calc_radiance_dn(nspec, nlev, &
         &  ONE_OVER_PI, transmittance, source_dn, v_overlap, radiance_profile)

    ! Extract surface value
    radiance = radiance_profile(:,nlev+1)

  end if

  if (lhook) call dr_hook('tcrad:calc_no_scattering_radiance',1,hook_handle)

end subroutine calc_no_scattering_radiance

