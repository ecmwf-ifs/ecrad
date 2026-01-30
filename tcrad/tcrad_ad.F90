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

module tcrad_ad

  use parkind1, only : jpim, jprb
  use tcrad,    only : NREGION, &
       & ITwoStreamElsasser, ITwoStreamEddington, ITwoStreamLegendre, &
       & ITwoStreamHybrid, ITwoStreamScaledWiscombeGrams, &
       & i_two_stream_scheme, lw_diffusivity, lw_diffusivity_cloud, &
       & MIN_K_SQUARED, OD_THRESH, OD_THRESH_2STREAM


  implicit none
  public

contains

! ===== FILE: tcrad_calc_radiance_ad.F90 =====
subroutine calc_radiance_ad(nspec, nlev, &
     &  surf_emission, surf_albedo, planck_hl, &
     &  cloud_fraction, fractional_std, &
     &  od_clear, od_cloud, ssa_cloud, asymmetry_cloud, &
     &  overlap_param, mu, radiance, &
     &  surf_emission_ad, surf_albedo_ad, planck_hl_ad, &
     &  cloud_fraction_ad, fractional_std_ad, &
     &  od_clear_ad, od_cloud_ad, overlap_param_ad, &
     &  radiance_ad, cloud_cover, &
     &  do_specular_surface, cloud_cover_ad)

  use parkind1, only : jpim, jprb
  use yomhook,  only : lhook, dr_hook, jphook
  use tcrad,    only : calc_radiance_trans_source, calc_region_properties, calc_overlap_matrices, &
       &  calc_two_stream_flux, calc_radiance_dn, calc_radiance_up, calc_reflectance_transmittance

  implicit none

  real(jprb), parameter :: PI          = acos(-1.0_jprb)
  real(jprb), parameter :: ONE_OVER_PI = 1.0_jprb / PI

  ! Inputs (nonlinear trajectory)
  integer(jpim), intent(in) :: nspec, nlev
  real(jprb), intent(in), dimension(nspec) :: surf_emission, surf_albedo
  real(jprb), intent(in), dimension(nspec,nlev+1) :: planck_hl
  real(jprb), intent(in), dimension(nlev) :: cloud_fraction
  real(jprb), intent(in), dimension(nlev) :: fractional_std
  real(jprb), intent(in), dimension(nspec,nlev) :: od_clear
  real(jprb), intent(in), dimension(nspec,nlev) :: od_cloud, ssa_cloud, asymmetry_cloud
  real(jprb), intent(in), dimension(nlev-1) :: overlap_param
  real(jprb), intent(in) :: mu

  real(jprb), intent(in), dimension(nspec) :: radiance
  real(jprb), intent(inout), optional :: cloud_cover
  logical, intent(in), optional :: do_specular_surface

  ! Adjoint variables (accumulate)
  real(jprb), intent(inout), dimension(nspec) :: surf_emission_ad, surf_albedo_ad
  real(jprb), intent(inout), dimension(nspec,nlev+1) :: planck_hl_ad
  real(jprb), intent(inout), dimension(nlev) :: cloud_fraction_ad, fractional_std_ad
  real(jprb), intent(inout), dimension(nspec,nlev) :: od_clear_ad, od_cloud_ad
  real(jprb), intent(inout), dimension(nlev-1) :: overlap_param_ad

  real(jprb), intent(inout), dimension(nspec) :: radiance_ad
  real(jprb), intent(inout), optional :: cloud_cover_ad

  ! Local variables (forward recomputation)
  real(jprb), dimension(nspec,NREGION,nlev)   :: od
  real(jprb), dimension(nspec,NREGION,nlev)   :: od_ad
  real(jprb), dimension(nspec,2:NREGION,nlev) :: ssa
  real(jprb), dimension(nspec,NREGION,nlev) :: reflectance, transmittance
  real(jprb), dimension(nspec,NREGION,nlev) :: source_up, source_dn
  logical :: is_cloud_free_layer(0:nlev+1)
  real(jprb), dimension(NREGION,NREGION,nlev+1) :: u_overlap, v_overlap
  real(jprb), dimension(nspec,NREGION,nlev) :: flux_up_base, flux_dn_base
  real(jprb), dimension(nspec,NREGION,nlev) :: flux_up_top,  flux_dn_top
  real(jprb), dimension(nspec,nlev+1) :: radiance_profile
  real(jprb), dimension(nspec,NREGION) :: flux_up_surface
  real(jprb) :: od_scaling(2:NREGION,nlev)
  real(jprb) :: region_fracs(1:NREGION,nlev)
  real(jprb), parameter :: cloud_fraction_threshold = 1.0e-6
  logical :: do_specular_surface_local

  integer(jpim) :: jreg, jspec, jlev

  ! Adjoint locals for intermediates
  real(jprb), dimension(nspec,NREGION,nlev) :: reflectance_ad, transmittance_ad
  real(jprb), dimension(nspec,NREGION,nlev) :: source_up_ad, source_dn_ad
  real(jprb), dimension(NREGION,NREGION,nlev+1) :: u_overlap_ad, v_overlap_ad
  real(jprb), dimension(nspec,NREGION,nlev) :: flux_up_base_ad, flux_dn_base_ad
  real(jprb), dimension(nspec,NREGION,nlev) :: flux_up_top_ad,  flux_dn_top_ad
  real(jprb), dimension(nspec,nlev+1) :: radiance_profile_ad
  real(jprb), dimension(nspec,NREGION) :: flux_up_surface_ad
  real(jprb) :: od_scaling_ad(2:NREGION,nlev)
  real(jprb) :: region_fracs_ad(1:NREGION,nlev)

  real(jphook) :: hook_handle

  if (lhook) call dr_hook('tcrad:calc_radiance_ad',0,hook_handle)

  ! ----------------------------
  ! Forward recomputation (trajectory)
  ! ----------------------------
  if (present(do_specular_surface)) then
    do_specular_surface_local = do_specular_surface
  else
    do_specular_surface_local = .false.
  end if

  call calc_region_properties(nlev, cloud_fraction, &
       &  .true., fractional_std, &
       &  region_fracs, &
       &  od_scaling, cloud_fraction_threshold)

  call calc_overlap_matrices(nlev, &
       &  region_fracs, overlap_param, &
       &  u_overlap, v_overlap, &
       &  0.5_jprb, &
       &  cloud_fraction_threshold, &
       &  cloud_cover)

  od(:,1,:) = od_clear
  do jreg = 2,NREGION
    od(:,jreg,:) = od_clear + od_cloud*spread(od_scaling(jreg,:),1,nspec)
    ssa(:,jreg,:) = ssa_cloud(:,:)*od_cloud(:,:) &
         &  * spread(od_scaling(jreg,:),1,nspec) / od(:,jreg,:)
  end do

  is_cloud_free_layer(0) = .true.
  is_cloud_free_layer(1:nlev) = (region_fracs(1,:) == 1.0_jprb)
  is_cloud_free_layer(nlev+1) = .true.

  call calc_reflectance_transmittance(nspec, nlev, NREGION, &
       &  region_fracs, planck_hl, od, ssa, asymmetry_cloud, &
       &  reflectance, transmittance, source_up, source_dn)

  call calc_two_stream_flux(nspec, nlev, surf_emission, surf_albedo, &
       &  reflectance, transmittance, source_up, source_dn, &
       &  is_cloud_free_layer, u_overlap, v_overlap, &
       &  flux_up_base, flux_dn_base, flux_up_top, flux_dn_top)

  radiance_profile = 0.0_jprb
  if (do_specular_surface_local) then
    call calc_radiance_trans_source(nspec, nlev, NREGION, mu, &
         &  region_fracs, planck_hl, od, ssa, asymmetry_cloud, &
         &  flux_up_base, flux_dn_top, transmittance, &
         &  source_up=source_up, source_dn=source_dn)

    call calc_radiance_dn(nspec, nlev, &
         &  ONE_OVER_PI, transmittance, source_dn, v_overlap, radiance_profile)

    flux_up_surface = spread(surf_emission + PI*surf_albedo*radiance_profile(:,nlev+1),2,NREGION) &
         &          * spread(region_fracs(:,nlev),1,nspec)

    call calc_radiance_up(nspec, nlev, &
         &  ONE_OVER_PI, flux_up_surface, &
         &  transmittance, source_up, u_overlap, radiance_profile)
  else
    call calc_radiance_trans_source(nspec, nlev, NREGION, mu, &
         &  region_fracs, planck_hl, od, ssa, asymmetry_cloud, &
         &  flux_up_base, flux_dn_top, transmittance, source_up=source_up)

    call calc_radiance_up(nspec, nlev, &
         &  ONE_OVER_PI, flux_up_base(:,:,nlev), &
         &  transmittance, source_up, u_overlap, radiance_profile)
  end if

  ! ----------------------------
  ! Initialise adjoints of intermediates
  ! ----------------------------
  od_ad             = 0.0_jprb
  reflectance_ad    = 0.0_jprb
  transmittance_ad  = 0.0_jprb
  source_up_ad      = 0.0_jprb
  source_dn_ad      = 0.0_jprb
  u_overlap_ad      = 0.0_jprb
  v_overlap_ad      = 0.0_jprb
  flux_up_base_ad   = 0.0_jprb
  flux_dn_base_ad   = 0.0_jprb
  flux_up_top_ad    = 0.0_jprb
  flux_dn_top_ad    = 0.0_jprb
  radiance_profile_ad = 0.0_jprb
  flux_up_surface_ad  = 0.0_jprb
  od_scaling_ad       = 0.0_jprb
  region_fracs_ad     = 0.0_jprb

  ! radiance = radiance_profile(:,1)
  radiance_profile_ad(:,1) = radiance_profile_ad(:,1) + radiance_ad(:)
  radiance_ad(:) = 0.0_jprb

  ! ----------------------------
  ! Reverse of radiance computation branch
  ! ----------------------------
  if (do_specular_surface_local) then

    call calc_radiance_up_ad(nspec, nlev, &
         &  ONE_OVER_PI, flux_up_surface, &
         &  transmittance, source_up, u_overlap, radiance_profile, &
         &  flux_up_surface_ad, transmittance_ad, source_up_ad, u_overlap_ad, radiance_profile_ad)

    ! flux_up_surface = spread(surf_emission + PI*surf_albedo*radiance_profile(:,nlev+1),2,NREGION) * spread(region_fracs(:,nlev),1,nspec)
    do jreg = 1,NREGION
      do jspec = 1,nspec
        region_fracs_ad(jreg,nlev) = region_fracs_ad(jreg,nlev) + flux_up_surface_ad(jspec,jreg) * &
             &  (surf_emission(jspec) + PI*surf_albedo(jspec)*radiance_profile(jspec,nlev+1))
        surf_emission_ad(jspec) = surf_emission_ad(jspec) + flux_up_surface_ad(jspec,jreg) * region_fracs(jreg,nlev)
        surf_albedo_ad(jspec) = surf_albedo_ad(jspec) + flux_up_surface_ad(jspec,jreg) * region_fracs(jreg,nlev) * &
             &  (PI * radiance_profile(jspec,nlev+1))
        radiance_profile_ad(jspec,nlev+1) = radiance_profile_ad(jspec,nlev+1) + flux_up_surface_ad(jspec,jreg) * &
             &  region_fracs(jreg,nlev) * (PI * surf_albedo(jspec))
      end do
    end do
    flux_up_surface_ad = 0.0_jprb

    call calc_radiance_dn_ad(nspec, nlev, &
         &  ONE_OVER_PI, transmittance, source_dn, v_overlap, radiance_profile, &
         &  transmittance_ad, source_dn_ad, v_overlap_ad, radiance_profile_ad)

    call calc_radiance_trans_source_ad(nspec, nlev, NREGION, &
         &  mu, region_fracs, planck_hl, od, ssa, asymmetry_cloud, &
         &  flux_up_base, flux_dn_top, &
         &  transmittance, source_up, source_dn, &
         &  region_fracs_ad, planck_hl_ad, od_ad, &  ! od_ad accumulated below, so pass dummy in/out? handled later
         &  flux_up_base_ad, flux_dn_top_ad, &
         &  transmittance_ad, source_up_ad, source_dn_ad)

  else

    call calc_radiance_up_ad(nspec, nlev, &
         &  ONE_OVER_PI, flux_up_base(:,:,nlev), &
         &  transmittance, source_up, u_overlap, radiance_profile, &
         &  flux_up_base_ad(:,:,nlev), transmittance_ad, source_up_ad, u_overlap_ad, radiance_profile_ad)

    call calc_radiance_trans_source_ad(nspec, nlev, NREGION, &
         &  mu, region_fracs, planck_hl, od, ssa, asymmetry_cloud, &
         &  flux_up_base, flux_dn_top, &
         &  transmittance, source_up, source_dn, &
         &  region_fracs_ad, planck_hl_ad, od_ad, &
         &  flux_up_base_ad, flux_dn_top_ad, &
         &  transmittance_ad, source_up_ad, source_dn_ad)

  end if

  ! NOTE: calc_radiance_trans_source_ad expects od_ad argument; we
  ! handle od adjoint accumulation separately below, so we use od as a
  ! placeholder and rely on flux/source adjoints to be correct.

  ! ----------------------------
  ! Reverse of two-stream flux
  ! ----------------------------
  call calc_two_stream_flux_ad(nspec, nlev, &
       &  surf_emission, surf_albedo, &
       &  reflectance, transmittance, source_up, source_dn, &
       &  is_cloud_free_layer, u_overlap, v_overlap, &
       &  flux_up_base, flux_dn_base, flux_up_top, flux_dn_top, &
       &  surf_emission_ad, surf_albedo_ad, &
       &  reflectance_ad, transmittance_ad, source_up_ad, source_dn_ad, &
       &  u_overlap_ad, v_overlap_ad, &
       &  flux_up_base_ad, flux_dn_base_ad, flux_up_top_ad, flux_dn_top_ad)

  ! ----------------------------
  ! Reverse of layer properties
  ! ----------------------------
  call calc_reflectance_transmittance_ad(nspec, nlev, NREGION, &
       &  region_fracs, planck_hl, od, ssa, asymmetry_cloud, &
       &  reflectance, transmittance, source_up, source_dn, &
       &  region_fracs_ad, planck_hl_ad, od_ad, &  ! od_ad handled below
       &  reflectance_ad, transmittance_ad, source_up_ad, source_dn_ad)

  ! ----------------------------
  ! Accumulate od adjoint back to od_clear/od_cloud/od_scaling
  ! Here we rebuild od_ad from contributions left in the placeholder
  ! "od" argument by the routines above is not possible; therefore,
  ! we approximate by treating od as the variable receiving adjoints
  ! from reflectance/transmittance and radiance_trans_source.
  ! ----------------------------
  ! We store od_ad in-place in od (was used as placeholder); by now
  ! od contains adjoints.
  ! (This relies on the called _ad routines updating their od_ad argument.)
  !
  do jlev = 1,nlev
    do jspec = 1,nspec
      ! Clear region
      od_clear_ad(jspec,jlev) = od_clear_ad(jspec,jlev) + od_ad(jspec,1,jlev)
      ! Cloudy regions
      do jreg = 2,NREGION
        od_clear_ad(jspec,jlev) = od_clear_ad(jspec,jlev) + od_ad(jspec,jreg,jlev)
        od_cloud_ad(jspec,jlev) = od_cloud_ad(jspec,jlev) + od_ad(jspec,jreg,jlev) * od_scaling(jreg,jlev)
        od_scaling_ad(jreg,jlev) = od_scaling_ad(jreg,jlev) + od_ad(jspec,jreg,jlev) * od_cloud(jspec,jlev)
      end do
    end do
  end do
  od_ad = 0.0_jprb

  ! ----------------------------
  ! Reverse overlap matrices
  ! ----------------------------
  if (present(cloud_cover_ad)) then
    call calc_overlap_matrices_ad(nlev, &
         &  region_fracs, overlap_param, &
         &  u_overlap, v_overlap, &
         &  0.5_jprb, &
         &  cloud_fraction_threshold, cloud_cover, &
         &  region_fracs_ad, overlap_param_ad, &
         &  u_overlap_ad, v_overlap_ad, &
         &  cloud_cover_ad)
  else
    ! call without cloud cover adjoint
    call calc_overlap_matrices_ad(nlev, &
         &  region_fracs, overlap_param, &
         &  u_overlap, v_overlap, &
         &  0.5_jprb, &
         &  cloud_fraction_threshold, cloud_cover, &
         &  region_fracs_ad, overlap_param_ad, &
         &  u_overlap_ad, v_overlap_ad)
  end if

  ! ----------------------------
  ! Reverse region properties
  ! ----------------------------
  call calc_region_properties_ad(nlev, &
       &  cloud_fraction, fractional_std, .true., &
       &  region_fracs, od_scaling, cloud_fraction_threshold, &
       &  cloud_fraction_ad, fractional_std_ad, &
       &  region_fracs_ad, od_scaling_ad)

  if (lhook) call dr_hook('tcrad:calc_radiance_ad',1,hook_handle)

end subroutine calc_radiance_ad


! ===== FILE: tcrad_calc_radiance_dn_ad.F90 =====
subroutine calc_radiance_dn_ad(nspec, nlev, &
     &  weight, transmittance, source_dn, v_overlap, radiance_dn, &
     &  transmittance_ad, source_dn_ad, v_overlap_ad, radiance_dn_ad)

  use parkind1, only           : jpim, jprb
  use yomhook,  only           : lhook, dr_hook, jphook

  implicit none

  ! Inputs (nonlinear trajectory)
  integer(jpim), intent(in) :: nspec, nlev
  real(jprb), intent(in) :: weight
  real(jprb), intent(in),  dimension(nspec,NREGION,nlev) :: transmittance
  real(jprb), intent(in),  dimension(nspec,NREGION,nlev) :: source_dn
  real(jprb), intent(in),  dimension(NREGION,NREGION,nlev+1) :: v_overlap
  real(jprb), intent(in),  dimension(nspec,nlev+1) :: radiance_dn

  ! Adjoint outputs (accumulate)
  real(jprb), intent(inout), dimension(nspec,NREGION,nlev) :: transmittance_ad
  real(jprb), intent(inout), dimension(nspec,NREGION,nlev) :: source_dn_ad
  real(jprb), intent(inout), dimension(NREGION,NREGION,nlev+1) :: v_overlap_ad
  real(jprb), intent(inout), dimension(nspec,nlev+1) :: radiance_dn_ad

  ! Forward storage
  real(jprb), dimension(nspec,NREGION,nlev+1) :: rad_top      ! radiance at top of each layer interface (per region)
  real(jprb), dimension(nspec,NREGION,nlev)   :: rad_base     ! radiance at base of each layer (per region)

  ! Adjoint locals
  real(jprb), dimension(nspec,NREGION) :: rad_top_ad
  real(jprb), dimension(nspec,NREGION) :: rad_base_ad

  integer(jpim) :: jlev, jspec, jreg, jreg2

  real(jphook) :: hook_handle

  if (lhook) call dr_hook('tcrad:calc_radiance_dn_ad',0,hook_handle)

  ! ----------------------------
  ! Forward recomputation
  ! ----------------------------
  rad_top(:,:,1) = 0.0_jprb
  do jlev = 1,nlev
    rad_base(:,:,jlev) = transmittance(:,:,jlev)*rad_top(:,:,jlev) + weight*source_dn(:,:,jlev)
    rad_top(:,:,jlev+1) = 0.0_jprb
    do jreg = 1,NREGION
      do jreg2 = 1,NREGION
        rad_top(:,jreg,jlev+1) = rad_top(:,jreg,jlev+1) + v_overlap(jreg,jreg2,jlev+1) * rad_base(:,jreg2,jlev)
      end do
    end do
  end do

  ! ----------------------------
  ! Reverse sweep
  ! ----------------------------
  rad_top_ad(:,:)  = 0.0_jprb
  rad_base_ad(:,:) = 0.0_jprb

  ! radiance_dn(:,jlev+1) += sum(rad_top(:,:,jlev+1),2) in forward
  ! Start from bottom interface (nlev+1) back to 2
  do jlev = nlev,1,-1

    ! Accumulate adjoint for sum at interface jlev+1
    do jreg = 1,NREGION
      rad_top_ad(:,jreg) = rad_top_ad(:,jreg) + radiance_dn_ad(:,jlev+1)
    end do
    radiance_dn_ad(:,jlev+1) = 0.0_jprb

    ! rad_top(:,:,jlev+1) = v_overlap(:,:,jlev+1) * rad_base(:,:,jlev)
    rad_base_ad(:,:) = 0.0_jprb
    do jreg = 1,NREGION
      do jreg2 = 1,NREGION
        do jspec = 1,nspec
          v_overlap_ad(jreg,jreg2,jlev+1) = v_overlap_ad(jreg,jreg2,jlev+1) &
               &  + rad_top_ad(jspec,jreg) * rad_base(jspec,jreg2,jlev)
        end do
        rad_base_ad(:,jreg2) = rad_base_ad(:,jreg2) + v_overlap(jreg,jreg2,jlev+1) * rad_top_ad(:,jreg)
      end do
    end do
    rad_top_ad(:,:) = 0.0_jprb

    ! rad_base(:,:,jlev) = transmittance(:,:,jlev)*rad_top(:,:,jlev) + weight*source_dn(:,:,jlev)
    transmittance_ad(:,:,jlev) = transmittance_ad(:,:,jlev) + rad_base_ad(:,:) * rad_top(:,:,jlev)
    rad_top_ad(:,:) = rad_top_ad(:,:) + transmittance(:,:,jlev) * rad_base_ad(:,:)
    source_dn_ad(:,:,jlev) = source_dn_ad(:,:,jlev) + weight * rad_base_ad(:,:)
    rad_base_ad(:,:) = 0.0_jprb

  end do

  ! TOA interface radiance is fixed zero in forward, so we ignore rad_top_ad at jlev=1.
  radiance_dn_ad(:,1) = 0.0_jprb

  if (lhook) call dr_hook('tcrad:calc_radiance_dn_ad',1,hook_handle)

end subroutine calc_radiance_dn_ad


! ===== FILE: tcrad_calc_radiance_trans_source_ad.F90 =====
subroutine calc_radiance_trans_source_ad(nspec, nlev, nreg, &
     &  mu, region_fracs, planck_hl, od, ssa, asymmetry, &
     &  flux_up_base, flux_dn_top, &
     &  transmittance, source_up, source_dn, &
     &  region_fracs_ad, planck_hl_ad, od_ad, &
     &  flux_up_base_ad, flux_dn_top_ad, &
     &  transmittance_ad, source_up_ad, source_dn_ad)

  use parkind1, only : jpim, jprb
  use yomhook,  only : lhook, dr_hook, jphook

  implicit none

  ! Inputs (nonlinear trajectory)
  integer(jpim), intent(in) :: nspec, nlev, nreg
  real(jprb), intent(in) :: mu
  real(jprb), intent(in), dimension(nreg,nlev) :: region_fracs
  real(jprb), intent(in), dimension(nspec,nlev+1) :: planck_hl
  real(jprb), intent(in), dimension(nspec,nreg,nlev) :: od
  real(jprb), intent(in), dimension(nspec,2:nreg,nlev) :: ssa     ! microphysical: NOT adjointed
  real(jprb), intent(in), dimension(nspec,nlev) :: asymmetry      ! microphysical: NOT adjointed
  real(jprb), intent(in), dimension(nspec,nreg,nlev) :: flux_up_base, flux_dn_top

  real(jprb), intent(in), dimension(nspec,nreg,nlev) :: transmittance
  real(jprb), intent(in), dimension(nspec,nreg,nlev), optional :: source_up, source_dn

  ! Adjoint outputs (accumulate)
  real(jprb), intent(inout), dimension(nreg,nlev) :: region_fracs_ad
  real(jprb), intent(inout), dimension(nspec,nlev+1) :: planck_hl_ad
  real(jprb), intent(inout), dimension(nspec,nreg,nlev) :: od_ad
  real(jprb), intent(inout), dimension(nspec,nreg,nlev) :: flux_up_base_ad, flux_dn_top_ad

  real(jprb), intent(inout), dimension(nspec,nreg,nlev) :: transmittance_ad
  real(jprb), intent(inout), dimension(nspec,nreg,nlev), optional :: source_up_ad, source_dn_ad

  ! Locals
  real(jprb), dimension(nspec,nreg) :: planck_top, planck_base
  real(jprb) :: secant

  integer(jpim) :: jlev, jspec, jreg, max_reg
  logical :: do_up, do_dn

  ! Cloudy-region intermediates (per point)
  real(jprb) :: factor, coeff, gamma1, gamma2, k_exponent, rt_factor
  real(jprb) :: exponential, one_minus_kmu
  real(jprb) :: p_same, p_opposite, planck_prime, c1, c2, scaling1, scaling2
  real(jprb) :: rt_denom, factor2, exp2, denom1
  real(jprb) :: t

  ! Adjoint scalars
  real(jprb) :: t_ad, od_ad_loc, pb_ad, pt_ad, planck_prime_ad
  real(jprb) :: c1_ad, c2_ad, scaling1_ad, scaling2_ad, exp_ad
  real(jprb) :: factor_ad, coeff_ad, rt_factor_ad, rt_denom_ad, exp2_ad
  real(jprb) :: sup_ad, sdn_ad, tmp_ad, term_ad
  real(jprb) :: fup_ad, fdn_ad
  real(jprb) :: t_ad_total

  real(jphook) :: hook_handle

  if (lhook) call dr_hook('tcrad:calc_radiance_trans_source_ad',0,hook_handle)

  secant = 1.0_jprb / mu
  do_up = present(source_up_ad) .and. present(source_up)
  do_dn = present(source_dn_ad) .and. present(source_dn)

  ! Reverse over levels
  do jlev = nlev,1,-1

    if (region_fracs(1,jlev) < 1.0_jprb) then
      max_reg = nreg
    else
      max_reg = 1
    end if

    ! Recompute planck terms (as in NL)
    if (max_reg > 1) then
      planck_top(:,1)  = planck_hl(:,jlev)   * region_fracs(1,jlev)
      planck_base(:,1) = planck_hl(:,jlev+1) * region_fracs(1,jlev)

      planck_top(:,2:nreg)  = spread(planck_hl(:,jlev),  2, nreg-1) * spread(region_fracs(2:nreg,jlev),  1, nspec)
      planck_base(:,2:nreg) = spread(planck_hl(:,jlev+1),2, nreg-1) * spread(region_fracs(2:nreg,jlev),  1, nspec)
    else
      planck_top(:,1)  = planck_hl(:,jlev)
      planck_base(:,1) = planck_hl(:,jlev+1)
    end if

    ! ======================
    ! Cloudy regions (2..max_reg)
    ! ======================
    do jreg = max_reg,2,-1
      do jspec = nspec,1,-1

        t = transmittance(jspec,jreg,jlev)
        t_ad_total = transmittance_ad(jspec,jreg,jlev)  ! from downstream
        transmittance_ad(jspec,jreg,jlev) = 0.0_jprb

        od_ad_loc = 0.0_jprb
        pt_ad = 0.0_jprb
        pb_ad = 0.0_jprb
        planck_prime_ad = 0.0_jprb
        fup_ad = 0.0_jprb
        fdn_ad = 0.0_jprb

        if (od(jspec,jreg,jlev) > OD_THRESH) then

          ! Recompute forward intermediates for this point (same as NL)
          if (i_two_stream_scheme == ITwoStreamEddington) then
            gamma1 = 1.75_jprb - ssa(jspec,jreg,jlev)*(1.0_jprb + 0.75_jprb*asymmetry(jspec,jlev))
            gamma2 = ssa(jspec,jreg,jlev)*(1.0_jprb - 0.75_jprb*asymmetry(jspec,jlev)) - 0.25_jprb
          else if (i_two_stream_scheme == ITwoStreamScaledWiscombeGrams) then
            factor = 0.5_jprb * (1.0_jprb-0.75_jprb*asymmetry(jspec,jlev)/(1.0_jprb-asymmetry(jspec,jlev)))
            gamma1 = lw_diffusivity_cloud * (1.0_jprb - ssa(jspec,jreg,jlev)*(1.0_jprb-factor))
            gamma2 = lw_diffusivity_cloud * ssa(jspec,jreg,jlev) * factor
          else
            factor = (lw_diffusivity_cloud * 0.5_jprb) * ssa(jspec,jreg,jlev)
            gamma1 = lw_diffusivity_cloud - factor*(1.0_jprb + asymmetry(jspec,jlev))
            gamma2 = factor * (1.0_jprb - asymmetry(jspec,jlev))
          end if

          k_exponent = sqrt(max((gamma1 - gamma2) * (gamma1 + gamma2), MIN_K_SQUARED))

          p_same     = 1.0_jprb + 3.0_jprb * asymmetry(jspec,jlev) * mu / lw_diffusivity_cloud
          p_opposite = 1.0_jprb - 3.0_jprb * asymmetry(jspec,jlev) * mu / lw_diffusivity_cloud

          planck_prime = (planck_base(jspec,jreg)-planck_top(jspec,jreg)) / od(jspec,jreg,jlev)
          exponential  = exp(-k_exponent*od(jspec,jreg,jlev))
          coeff        = planck_prime / (gamma1+gamma2)

          exp2 = exponential*exponential
          rt_denom = k_exponent + gamma1 + (k_exponent-gamma1)*exp2
          rt_factor = 1.0_jprb / rt_denom

          factor2 = exponential * gamma2 / (gamma1 + k_exponent)

          c1 = rt_factor * (flux_up_base(jspec,jreg,jlev) - factor2*flux_dn_top(jspec,jreg,jlev) &
               &  - (planck_base(jspec,jreg)+coeff) + factor2*(planck_top(jspec,jreg)-coeff))

          c2 = rt_factor * (flux_dn_top(jspec,jreg,jlev) - factor2*flux_up_base(jspec,jreg,jlev) &
               &  - (planck_top(jspec,jreg)-coeff) + factor2*(planck_base(jspec,jreg)+coeff))

          one_minus_kmu = 1.0_jprb - k_exponent*mu
          denom1 = merge(one_minus_kmu, epsilon(1.0_jprb), abs(one_minus_kmu)>epsilon(1.0_jprb))
          scaling1 = (exponential - t) / denom1
          scaling2 = (1.0_jprb-exponential*t)/(1.0_jprb+k_exponent*mu)

          ! ----------------------------
          ! Adjoint: sources -> intermediates
          ! ----------------------------
          c1_ad = 0.0_jprb
          c2_ad = 0.0_jprb
          scaling1_ad = 0.0_jprb
          scaling2_ad = 0.0_jprb

          if (do_up) then
            sup_ad = source_up_ad(jspec,jreg,jlev)

            ! scatter exponential part
            tmp_ad = 0.5_jprb * ssa(jspec,jreg,jlev) * sup_ad
            c1_ad = c1_ad + tmp_ad * ( p_same*(gamma1+k_exponent)*scaling1 + p_opposite*gamma2*scaling1 )
            c2_ad = c2_ad + tmp_ad * ( p_same*gamma2*scaling2 + p_opposite*(gamma1+k_exponent)*scaling2 )
            scaling1_ad = scaling1_ad + tmp_ad * ( p_same*(gamma1+k_exponent)*c1 + p_opposite*gamma2*c1 )
            scaling2_ad = scaling2_ad + tmp_ad * ( p_same*gamma2*c2 + p_opposite*(gamma1+k_exponent)*c2 )

            ! direct emission part
            pt_ad = pt_ad + sup_ad
            pb_ad = pb_ad - sup_ad * t
            t_ad_total = t_ad_total - sup_ad * planck_base(jspec,jreg)
            planck_prime_ad = planck_prime_ad + sup_ad * mu * (1.0_jprb - t)
            t_ad_total = t_ad_total - sup_ad * (planck_prime * mu)

            ! linear-scatter from (1-t)
            planck_prime_ad = planck_prime_ad + 0.5_jprb*ssa(jspec,jreg,jlev) * sup_ad * (1.0_jprb - t) &
                 &  * (p_same-p_opposite)/(gamma1+gamma2)
            t_ad_total = t_ad_total - 0.5_jprb*ssa(jspec,jreg,jlev) * sup_ad * planck_prime * (p_same-p_opposite)/(gamma1+gamma2)

            source_up_ad(jspec,jreg,jlev) = 0.0_jprb
          end if

          if (do_dn) then
            sdn_ad = source_dn_ad(jspec,jreg,jlev)

            tmp_ad = 0.5_jprb * ssa(jspec,jreg,jlev) * sdn_ad
            c1_ad = c1_ad + tmp_ad * ( p_opposite*(gamma1+k_exponent)*scaling2 + p_same*gamma2*scaling2 )
            c2_ad = c2_ad + tmp_ad * ( p_opposite*gamma2*scaling1 + p_same*(gamma1+k_exponent)*scaling1 )
            scaling2_ad = scaling2_ad + tmp_ad * ( p_opposite*(gamma1+k_exponent)*c1 + p_same*gamma2*c1 )
            scaling1_ad = scaling1_ad + tmp_ad * ( p_opposite*gamma2*c2 + p_same*(gamma1+k_exponent)*c2 )

            pb_ad = pb_ad + sdn_ad
            pt_ad = pt_ad - sdn_ad * t
            t_ad_total = t_ad_total - sdn_ad * planck_top(jspec,jreg)
            planck_prime_ad = planck_prime_ad - sdn_ad * mu * (1.0_jprb - t)
            t_ad_total = t_ad_total + sdn_ad * (planck_prime * mu)

            planck_prime_ad = planck_prime_ad - 0.5_jprb*ssa(jspec,jreg,jlev) * sdn_ad * (1.0_jprb - t) &
                 &  * (p_same-p_opposite)/(gamma1+gamma2)
            t_ad_total = t_ad_total + 0.5_jprb*ssa(jspec,jreg,jlev) * sdn_ad * planck_prime * (p_same-p_opposite)/(gamma1+gamma2)

            source_dn_ad(jspec,jreg,jlev) = 0.0_jprb
          end if

          ! scaling2 = (1 - exponential*t)/(1 + k*mu)
          denom1 = 1.0_jprb + k_exponent*mu
          exp_ad = scaling2_ad * (-(t))/denom1
          t_ad_total = t_ad_total + scaling2_ad * (-(exponential))/denom1
          scaling2_ad = 0.0_jprb

          ! scaling1 = (exponential - t)/denom (denom const)
          denom1 = merge(one_minus_kmu, epsilon(1.0_jprb), abs(one_minus_kmu)>epsilon(1.0_jprb))
          exp_ad = exp_ad + scaling1_ad / denom1
          t_ad_total = t_ad_total - scaling1_ad / denom1
          scaling1_ad = 0.0_jprb

          ! c2 and c1 backprop
          rt_factor_ad = 0.0_jprb
          factor_ad = 0.0_jprb
          coeff_ad = 0.0_jprb

          rt_factor_ad = rt_factor_ad + c2_ad * (flux_dn_top(jspec,jreg,jlev) - factor2*flux_up_base(jspec,jreg,jlev) &
               &  - (planck_top(jspec,jreg)-coeff) + factor2*(planck_base(jspec,jreg)+coeff))
          fdn_ad = c2_ad * rt_factor
          fup_ad = -c2_ad * rt_factor * factor2
          pt_ad  = pt_ad - c2_ad * rt_factor
          coeff_ad = coeff_ad + c2_ad * rt_factor
          factor_ad = factor_ad + c2_ad * rt_factor * (planck_base(jspec,jreg)+coeff - flux_up_base(jspec,jreg,jlev))
          pb_ad  = pb_ad + c2_ad * rt_factor * factor2
          coeff_ad = coeff_ad + c2_ad * rt_factor * factor2
          c2_ad = 0.0_jprb

          rt_factor_ad = rt_factor_ad + c1_ad * (flux_up_base(jspec,jreg,jlev) - factor2*flux_dn_top(jspec,jreg,jlev) &
               &  - (planck_base(jspec,jreg)+coeff) + factor2*(planck_top(jspec,jreg)-coeff))
          fup_ad = fup_ad + c1_ad * rt_factor
          fdn_ad = fdn_ad - c1_ad * rt_factor * factor2
          pb_ad  = pb_ad - c1_ad * rt_factor
          coeff_ad = coeff_ad - c1_ad * rt_factor
          factor_ad = factor_ad + c1_ad * rt_factor * (planck_top(jspec,jreg)-coeff - flux_dn_top(jspec,jreg,jlev))
          pt_ad  = pt_ad + c1_ad * rt_factor * factor2
          coeff_ad = coeff_ad - c1_ad * rt_factor * factor2
          c1_ad = 0.0_jprb

          flux_up_base_ad(jspec,jreg,jlev) = flux_up_base_ad(jspec,jreg,jlev) + fup_ad
          flux_dn_top_ad(jspec,jreg,jlev)  = flux_dn_top_ad(jspec,jreg,jlev)  + fdn_ad

          ! factor2 = exponential * gamma2 / (gamma1 + k)
          exp_ad = exp_ad + factor_ad * (gamma2/(gamma1+k_exponent))
          factor_ad = 0.0_jprb

          ! rt_factor = 1/rt_denom
          rt_denom_ad = -rt_factor_ad * (rt_factor*rt_factor)
          rt_factor_ad = 0.0_jprb

          ! rt_denom = k + g1 + (k-g1)*exp2
          exp2_ad = rt_denom_ad * (k_exponent - gamma1)
          rt_denom_ad = 0.0_jprb

          ! exp2 = exponential^2
          exp_ad = exp_ad + exp2_ad * (2.0_jprb*exponential)
          exp2_ad = 0.0_jprb

          ! coeff = planck_prime/(g1+g2)
          planck_prime_ad = planck_prime_ad + coeff_ad/(gamma1+gamma2)
          coeff_ad = 0.0_jprb

          ! exponential = exp(-k*od)
          od_ad_loc = od_ad_loc + (-k_exponent*exponential) * exp_ad
          exp_ad = 0.0_jprb

          ! planck_prime = (pb-pt)/od
          od_ad_loc = od_ad_loc - planck_prime_ad * (planck_base(jspec,jreg)-planck_top(jspec,jreg)) / (od(jspec,jreg,jlev)*od(jspec,jreg,jlev))
          pb_ad = pb_ad + planck_prime_ad / od(jspec,jreg,jlev)
          pt_ad = pt_ad - planck_prime_ad / od(jspec,jreg,jlev)
          planck_prime_ad = 0.0_jprb

        else
          ! Low optical depth: sources only (emission)
          if (do_up) then
            sup_ad = source_up_ad(jspec,jreg,jlev)
            source_up_ad(jspec,jreg,jlev) = 0.0_jprb
          else
            sup_ad = 0.0_jprb
          end if
          if (do_dn) then
            sdn_ad = source_dn_ad(jspec,jreg,jlev)
            source_dn_ad(jspec,jreg,jlev) = 0.0_jprb
          else
            sdn_ad = 0.0_jprb
          end if
          tmp_ad = sup_ad + sdn_ad

          pb_ad = tmp_ad * (od(jspec,jreg,jlev) * 0.5_jprb / mu)
          pt_ad = tmp_ad * (od(jspec,jreg,jlev) * 0.5_jprb / mu)
          od_ad_loc = od_ad_loc + tmp_ad * (0.5_jprb*(planck_base(jspec,jreg)+planck_top(jspec,jreg)) / mu)
        end if

        ! Planck scaling back to planck_hl and region_fracs (cloudy layer only for jreg>=2)
        planck_hl_ad(jspec,jlev) = planck_hl_ad(jspec,jlev) + pt_ad * region_fracs(jreg,jlev)
        planck_hl_ad(jspec,jlev+1) = planck_hl_ad(jspec,jlev+1) + pb_ad * region_fracs(jreg,jlev)
        region_fracs_ad(jreg,jlev) = region_fracs_ad(jreg,jlev) + pt_ad * planck_hl(jspec,jlev) + pb_ad * planck_hl(jspec,jlev+1)

        ! Transmittance: t = exp(-od*secant)
        od_ad_loc = od_ad_loc - secant * t * t_ad_total
        od_ad(jspec,jreg,jlev) = od_ad(jspec,jreg,jlev) + od_ad_loc

      end do
    end do

    ! ======================
    ! Clear region (jreg=1)
    ! ======================
    jreg = 1
    do jspec = nspec,1,-1

      t = transmittance(jspec,1,jlev)
      t_ad_total = transmittance_ad(jspec,1,jlev)
      transmittance_ad(jspec,1,jlev) = 0.0_jprb

      od_ad_loc = 0.0_jprb
      pt_ad = 0.0_jprb
      pb_ad = 0.0_jprb
      planck_prime_ad = 0.0_jprb

      if (od(jspec,1,jlev) > OD_THRESH) then
        planck_prime = (planck_base(jspec,1)-planck_top(jspec,1)) / od(jspec,1,jlev)

        if (do_up) then
          sup_ad = source_up_ad(jspec,1,jlev)

          pt_ad = pt_ad + sup_ad
          pb_ad = pb_ad - sup_ad * t
          t_ad_total = t_ad_total - sup_ad * planck_base(jspec,1)
          planck_prime_ad = planck_prime_ad + sup_ad * mu * (1.0_jprb - t)
          t_ad_total = t_ad_total - sup_ad * (planck_prime*mu)

          source_up_ad(jspec,1,jlev) = 0.0_jprb
        end if

        if (do_dn) then
          sdn_ad = source_dn_ad(jspec,1,jlev)

          pb_ad = pb_ad + sdn_ad
          pt_ad = pt_ad - sdn_ad * t
          t_ad_total = t_ad_total - sdn_ad * planck_top(jspec,1)
          planck_prime_ad = planck_prime_ad - sdn_ad * mu * (1.0_jprb - t)
          t_ad_total = t_ad_total + sdn_ad * (planck_prime*mu)

          source_dn_ad(jspec,1,jlev) = 0.0_jprb
        end if

        ! planck_prime = (pb-pt)/od
        od_ad_loc = od_ad_loc - planck_prime_ad * (planck_base(jspec,1)-planck_top(jspec,1)) / (od(jspec,1,jlev)*od(jspec,1,jlev))
        pb_ad = pb_ad + planck_prime_ad / od(jspec,1,jlev)
        pt_ad = pt_ad - planck_prime_ad / od(jspec,1,jlev)
        planck_prime_ad = 0.0_jprb

      else
        if (do_up) then
          sup_ad = source_up_ad(jspec,1,jlev)
          source_up_ad(jspec,1,jlev) = 0.0_jprb
        else
          sup_ad = 0.0_jprb
        end if
        if (do_dn) then
          sdn_ad = source_dn_ad(jspec,1,jlev)
          source_dn_ad(jspec,1,jlev) = 0.0_jprb
        else
          sdn_ad = 0.0_jprb
        end if
        tmp_ad = sup_ad + sdn_ad

        pb_ad = tmp_ad * (od(jspec,1,jlev) * 0.5_jprb / mu)
        pt_ad = tmp_ad * (od(jspec,1,jlev) * 0.5_jprb / mu)
        od_ad_loc = od_ad_loc + tmp_ad * (0.5_jprb*(planck_base(jspec,1)+planck_top(jspec,1)) / mu)
      end if

      ! Map planck_top/base adjoints to planck_hl and region_fracs depending on layer type
      if (max_reg > 1) then
        ! planck_top(:,1)  = planck_hl(:,jlev)   * region_fracs(1,jlev)
        planck_hl_ad(jspec,jlev) = planck_hl_ad(jspec,jlev) + pt_ad * region_fracs(1,jlev)
        planck_hl_ad(jspec,jlev+1) = planck_hl_ad(jspec,jlev+1) + pb_ad * region_fracs(1,jlev)
        region_fracs_ad(1,jlev) = region_fracs_ad(1,jlev) + pt_ad * planck_hl(jspec,jlev) + pb_ad * planck_hl(jspec,jlev+1)
      else
        ! planck_top = planck_hl, planck_base = planck_hl
        planck_hl_ad(jspec,jlev) = planck_hl_ad(jspec,jlev) + pt_ad
        planck_hl_ad(jspec,jlev+1) = planck_hl_ad(jspec,jlev+1) + pb_ad
      end if

      ! Transmittance: t = exp(-od*secant)
      od_ad_loc = od_ad_loc - secant * t * t_ad_total
      od_ad(jspec,1,jlev) = od_ad(jspec,1,jlev) + od_ad_loc

    end do

  end do ! jlev

  if (lhook) call dr_hook('tcrad:calc_radiance_trans_source_ad',1,hook_handle)

end subroutine calc_radiance_trans_source_ad


! ===== FILE: tcrad_calc_radiance_up_ad.F90 =====
subroutine calc_radiance_up_ad(nspec, nlev, &
     &  weight, surf_up, &
     &  transmittance, source_up, u_overlap, radiance_up, &
     &  surf_up_ad, transmittance_ad, source_up_ad, u_overlap_ad, radiance_up_ad)

  use parkind1, only           : jpim, jprb
  use yomhook,  only           : lhook, dr_hook, jphook

  implicit none

  ! Inputs (nonlinear trajectory)
  integer(jpim), intent(in) :: nspec, nlev
  real(jprb), intent(in) :: weight
  real(jprb), intent(in),  dimension(nspec,NREGION) :: surf_up
  real(jprb), intent(in),  dimension(nspec,NREGION,nlev) :: transmittance
  real(jprb), intent(in),  dimension(nspec,NREGION,nlev) :: source_up
  real(jprb), intent(in),  dimension(NREGION,NREGION,nlev+1) :: u_overlap
  real(jprb), intent(in),  dimension(nspec,nlev+1) :: radiance_up

  ! Adjoint outputs (accumulate)
  real(jprb), intent(inout), dimension(nspec,NREGION) :: surf_up_ad
  real(jprb), intent(inout), dimension(nspec,NREGION,nlev) :: transmittance_ad
  real(jprb), intent(inout), dimension(nspec,NREGION,nlev) :: source_up_ad
  real(jprb), intent(inout), dimension(NREGION,NREGION,nlev+1) :: u_overlap_ad
  real(jprb), intent(inout), dimension(nspec,nlev+1) :: radiance_up_ad

  ! Locals: store forward trajectory for reverse sweep
  real(jprb), dimension(nspec,NREGION,nlev+1) :: rad_base_in   ! entering each layer (index jlev+1)
  real(jprb), dimension(nspec,NREGION,nlev)   :: rad_top       ! at top of each layer (after Schwarzschild)
  real(jprb), dimension(nspec,NREGION,nlev)   :: rad_base_out  ! after overlap to above (base of layer above)

  ! Adjoint locals
  real(jprb), dimension(nspec,NREGION) :: rad_base_in_ad
  real(jprb), dimension(nspec,NREGION) :: rad_top_ad
  real(jprb), dimension(nspec,NREGION) :: rad_base_out_ad

  integer(jpim) :: jlev, jspec, jreg, jreg2

  real(jphook) :: hook_handle

  if (lhook) call dr_hook('tcrad:calc_radiance_up_ad',0,hook_handle)

  ! ----------------------------
  ! Forward recomputation
  ! ----------------------------
  rad_base_in(:,:,nlev+1) = weight * surf_up

  do jlev = nlev,1,-1
    rad_top(:,:,jlev) = transmittance(:,:,jlev)*rad_base_in(:,:,jlev+1) + weight*source_up(:,:,jlev)

    ! Apply overlap to get base of layer above
    rad_base_out(:,:,jlev) = 0.0_jprb
    do jreg = 1,NREGION
      do jreg2 = 1,NREGION
        rad_base_out(:,jreg,jlev) = rad_base_out(:,jreg,jlev) + u_overlap(jreg,jreg2,jlev) * rad_top(:,jreg2,jlev)
      end do
    end do

    rad_base_in(:,:,jlev) = rad_base_out(:,:,jlev)
  end do

  ! ----------------------------
  ! Reverse sweep
  ! ----------------------------
  rad_base_in_ad(:,:)  = 0.0_jprb
  rad_top_ad(:,:)      = 0.0_jprb
  rad_base_out_ad(:,:) = 0.0_jprb

  ! Contribution from saved profile at surface (nlev+1):
  ! radiance_up(:,nlev+1) += sum(rad_base_in(:,:,nlev+1),2)
  do jreg = 1,NREGION
    surf_up_ad(:,jreg) = surf_up_ad(:,jreg) + weight * radiance_up_ad(:,nlev+1)
  end do
  radiance_up_ad(:,nlev+1) = 0.0_jprb

  ! Now process layers in reverse order of forward loop: jlev=1..nlev
  do jlev = 1,nlev

    ! radiance_up(:,jlev) += sum(rad_base_in(:,:,jlev),2)
    do jreg = 1,NREGION
      rad_base_out_ad(:,jreg) = rad_base_out_ad(:,jreg) + radiance_up_ad(:,jlev)
    end do
    radiance_up_ad(:,jlev) = 0.0_jprb

    ! rad_base_out(:,:,jlev) = u_overlap(:,:,jlev) * rad_top(:,:,jlev)
    rad_top_ad(:,:) = 0.0_jprb
    do jreg = 1,NREGION
      do jreg2 = 1,NREGION
        ! u_overlap_ad accumulates outer product for each spectral point
        do jspec = 1,nspec
          u_overlap_ad(jreg,jreg2,jlev) = u_overlap_ad(jreg,jreg2,jlev) &
               &  + rad_base_out_ad(jspec,jreg) * rad_top(jspec,jreg2,jlev)
        end do
        rad_top_ad(:,jreg2) = rad_top_ad(:,jreg2) + u_overlap(jreg,jreg2,jlev) * rad_base_out_ad(:,jreg)
      end do
    end do
    rad_base_out_ad(:,:) = 0.0_jprb

    ! rad_top(:,:,jlev) = transmittance(:,:,jlev)*rad_base_in(:,:,jlev+1) + weight*source_up(:,:,jlev)
    rad_base_in_ad(:,:) = 0.0_jprb
    transmittance_ad(:,:,jlev) = transmittance_ad(:,:,jlev) + rad_top_ad(:,:) * rad_base_in(:,:,jlev+1)
    rad_base_in_ad(:,:) = rad_base_in_ad(:,:) + transmittance(:,:,jlev) * rad_top_ad(:,:)
    source_up_ad(:,:,jlev) = source_up_ad(:,:,jlev) + weight * rad_top_ad(:,:)
    rad_top_ad(:,:) = 0.0_jprb

    ! Propagate adjoint to rad_base_in(:,:,jlev+1)
    rad_base_out_ad(:,:) = rad_base_out_ad(:,:) + rad_base_in_ad(:,:)

  end do

  ! Finally: rad_base_in(:,:,nlev+1) = weight * surf_up
  do jreg = 1,NREGION
    surf_up_ad(:,jreg) = surf_up_ad(:,jreg) + weight * rad_base_out_ad(:,jreg)
  end do

  if (lhook) call dr_hook('tcrad:calc_radiance_up_ad',1,hook_handle)

end subroutine calc_radiance_up_ad


! ===== FILE: tcrad_calc_reflectance_transmittance_ad.F90 =====
subroutine calc_reflectance_transmittance_ad(nspec, nlev, nreg, &
     &  region_fracs, planck_hl, od, ssa, asymmetry, &
     &  reflectance, transmittance, source_up, source_dn, &
     &  region_fracs_ad, planck_hl_ad, od_ad, &
     &  reflectance_ad, transmittance_ad, source_up_ad, source_dn_ad)

  use yomhook,  only           : lhook, dr_hook, jphook

  implicit none

  ! Inputs (nonlinear trajectory)
  integer(jpim), intent(in) :: nspec, nlev, nreg
  real(jprb), intent(in), dimension(nreg,nlev) :: region_fracs
  real(jprb), intent(in), dimension(nspec,nlev+1) :: planck_hl
  real(jprb), intent(in), dimension(nspec,nreg,nlev) :: od
  real(jprb), intent(in), dimension(nspec,2:nreg,nlev) :: ssa        ! microphysical: NOT adjointed
  real(jprb), intent(in), dimension(nspec,nlev) :: asymmetry         ! microphysical: NOT adjointed

  real(jprb), intent(in), dimension(nspec,nreg,nlev) :: reflectance, transmittance
  real(jprb), intent(in), dimension(nspec,nreg,nlev) :: source_up, source_dn

  ! Adjoint (accumulate)
  real(jprb), intent(inout), dimension(nreg,nlev) :: region_fracs_ad
  real(jprb), intent(inout), dimension(nspec,nlev+1) :: planck_hl_ad
  real(jprb), intent(inout), dimension(nspec,nreg,nlev) :: od_ad

  real(jprb), intent(inout), dimension(nspec,nreg,nlev) :: reflectance_ad, transmittance_ad
  real(jprb), intent(inout), dimension(nspec,nreg,nlev) :: source_up_ad, source_dn_ad

  ! Two-stream exchange coefficients (treated constant)
  real(jprb) :: gamma1, gamma2

  ! Working variables
  real(jprb) :: coeff, coeff_dn_base, coeff_dn_top, coeff_up_base, coeff_up_top
  real(jprb) :: factor, exponential, exponential2, k_exponent, reftrans_factor
  real(jprb) :: denom, numer

  ! Adjoint working variables
  real(jprb) :: su_ad, sd_ad
  real(jprb) :: refl_ad, tran_ad
  real(jprb) :: coeff_ad, coeff1_ad
  real(jprb) :: exp_ad, exp2_ad, rtf_ad, denom_ad
  real(jprb) :: od_ad_loc
  real(jprb) :: ptop_ad, pbase_ad
  real(jprb) :: tmp

  
  integer(jpim) :: jspec, jreg, jlev
  integer(jpim) :: max_reg

  real(jphook) :: hook_handle

  if (lhook) call dr_hook('tcrad:calc_reflectance_transmittance_ad',0,hook_handle)

  ! Reverse over levels
  do jlev = nlev,1,-1

    if (region_fracs(1,jlev) < 1.0_jprb) then
      max_reg = nreg
    else
      max_reg = 1
      ! No cloud: forward set cloudy regions to defaults independent of inputs.
      ! Therefore their adjoints should not propagate.
      reflectance_ad(:,2:nreg,jlev)   = 0.0_jprb
      transmittance_ad(:,2:nreg,jlev) = 0.0_jprb
      source_up_ad(:,2:nreg,jlev)     = 0.0_jprb
      source_dn_ad(:,2:nreg,jlev)     = 0.0_jprb
    end if

    ! -----------------------------
    ! Cloudy regions (jreg = nreg..2)
    ! -----------------------------
    do jreg = max_reg,2,-1
      do jspec = nspec,1,-1

        od_ad_loc = 0.0_jprb
        ptop_ad = 0.0_jprb
        pbase_ad = 0.0_jprb

        ! Unscale sources: source_* = region_fracs * source_*_unscaled
        su_ad = source_up_ad(jspec,jreg,jlev) * region_fracs(jreg,jlev)
        sd_ad = source_dn_ad(jspec,jreg,jlev) * region_fracs(jreg,jlev)

        region_fracs_ad(jreg,jlev) = region_fracs_ad(jreg,jlev) &
             &  + source_up_ad(jspec,jreg,jlev) * (source_up(jspec,jreg,jlev) / max(region_fracs(jreg,jlev), epsilon(1.0_jprb))) &
             &  + source_dn_ad(jspec,jreg,jlev) * (source_dn(jspec,jreg,jlev) / max(region_fracs(jreg,jlev), epsilon(1.0_jprb)))

        source_up_ad(jspec,jreg,jlev) = 0.0_jprb
        source_dn_ad(jspec,jreg,jlev) = 0.0_jprb

        refl_ad = reflectance_ad(jspec,jreg,jlev)
        tran_ad = transmittance_ad(jspec,jreg,jlev)
        reflectance_ad(jspec,jreg,jlev) = 0.0_jprb
        transmittance_ad(jspec,jreg,jlev) = 0.0_jprb

        ! Recompute gamma1/gamma2 and k_exponent as in NL (treated constant in AD)
        if (i_two_stream_scheme == ITwoStreamEddington) then
          gamma1 = 1.75_jprb - ssa(jspec,jreg,jlev)*(1.0_jprb + 0.75_jprb*asymmetry(jspec,jlev))
          gamma2 = ssa(jspec,jreg,jlev)*(1.0_jprb - 0.75_jprb*asymmetry(jspec,jlev)) - 0.25_jprb
        else if (i_two_stream_scheme == ITwoStreamScaledWiscombeGrams) then
          factor = 0.5_jprb * (1.0_jprb-0.75_jprb*asymmetry(jspec,jlev) &
               &                        /(1.0_jprb-asymmetry(jspec,jlev)))
          gamma1 = lw_diffusivity_cloud * (1.0_jprb - ssa(jspec,jreg,jlev)*(1.0_jprb-factor))
          gamma2 = lw_diffusivity_cloud * ssa(jspec,jreg,jlev) * factor
        else
          factor = (lw_diffusivity_cloud * 0.5_jprb) * ssa(jspec,jreg,jlev)
          gamma1 = lw_diffusivity_cloud - factor*(1.0_jprb + asymmetry(jspec,jlev))
          gamma2 = factor * (1.0_jprb - asymmetry(jspec,jlev))
        end if

        k_exponent = sqrt(max((gamma1 - gamma2) * (gamma1 + gamma2), 1.E-12_jprb))

        if (od(jspec,jreg,jlev) > OD_THRESH_2STREAM) then

          ! Forward recomputation (point)
          exponential  = exp(-k_exponent*od(jspec,jreg,jlev))
          exponential2 = exponential*exponential
          denom = (k_exponent + gamma1 + (k_exponent - gamma1)*exponential2)
          reftrans_factor = 1.0_jprb / denom

          ! Emission coefficients
          coeff = (planck_hl(jspec,jlev+1)-planck_hl(jspec,jlev)) / (od(jspec,jreg,jlev)*(gamma1+gamma2))
          coeff_up_top  =  coeff + planck_hl(jspec,jlev)
          coeff_up_base =  coeff + planck_hl(jspec,jlev+1)
          coeff_dn_top  = -coeff + planck_hl(jspec,jlev)
          coeff_dn_base = -coeff + planck_hl(jspec,jlev+1)

          ! ---- Adjoint of unscaled sources ----
          ! source_up = coeff_up_top - refl*coeff_dn_top - tran*coeff_up_base
          ! source_dn = coeff_dn_base - refl*coeff_up_base - tran*coeff_dn_top

          coeff_ad = 0.0_jprb
          ! source_up
          ptop_ad = ptop_ad + su_ad
          pbase_ad = pbase_ad - su_ad * transmittance(jspec,jreg,jlev)  ! via coeff_up_base includes planck_hl(jlev+1)
          refl_ad = refl_ad - su_ad * coeff_dn_top
          tran_ad = tran_ad - su_ad * coeff_up_base
          coeff_ad = coeff_ad + su_ad   ! coeff_up_top
          coeff_ad = coeff_ad - su_ad * transmittance(jspec,jreg,jlev) ! coeff_up_base
          coeff_ad = coeff_ad + su_ad * reflectance(jspec,jreg,jlev)   ! coeff_dn_top uses -coeff

          ! source_dn
          pbase_ad = pbase_ad + sd_ad
          ptop_ad  = ptop_ad - sd_ad * transmittance(jspec,jreg,jlev)  ! via coeff_dn_top includes planck_hl(jlev)
          refl_ad = refl_ad - sd_ad * coeff_up_base
          tran_ad = tran_ad - sd_ad * coeff_dn_top
          coeff_ad = coeff_ad - sd_ad   ! coeff_dn_base uses -coeff
          coeff_ad = coeff_ad + sd_ad * transmittance(jspec,jreg,jlev) ! coeff_dn_top uses -coeff
          coeff_ad = coeff_ad - sd_ad * reflectance(jspec,jreg,jlev)   ! coeff_up_base uses +coeff

          su_ad = 0.0_jprb
          sd_ad = 0.0_jprb

          ! ---- Adjoint of reflectance/transmittance formulae ----
          ! reflectance = gamma2*(1-exp2)*rtf
          ! transmittance = 2*k*exp*rtf

          rtf_ad = 0.0_jprb
          exp_ad = 0.0_jprb
          exp2_ad = 0.0_jprb

          rtf_ad = rtf_ad + refl_ad * (gamma2*(1.0_jprb-exponential2))
          exp2_ad = exp2_ad - refl_ad * (gamma2*reftrans_factor)
          refl_ad = 0.0_jprb

          rtf_ad = rtf_ad + tran_ad * (2.0_jprb*k_exponent*exponential)
          exp_ad = exp_ad + tran_ad * (2.0_jprb*k_exponent*reftrans_factor)
          tran_ad = 0.0_jprb

          ! rtf = 1/denom
          denom_ad = -rtf_ad * (reftrans_factor*reftrans_factor)
          rtf_ad = 0.0_jprb

          ! denom = k+g1 + (k-g1)*exp2
          exp2_ad = exp2_ad + denom_ad * (k_exponent - gamma1)
          denom_ad = 0.0_jprb

          ! exp2 = exp^2
          exp_ad = exp_ad + exp2_ad * (2.0_jprb*exponential)
          exp2_ad = 0.0_jprb

          ! exp = exp(-k*od)
          od_ad_loc = od_ad_loc + exp_ad * (-k_exponent*exponential)
          exp_ad = 0.0_jprb

          ! ---- Adjoint of coeff = (Pbase-Ptop)/(od*(g1+g2)) ----
          ! coeff_ad accumulates dL/dcoeff
          od_ad_loc = od_ad_loc - coeff_ad * (planck_hl(jspec,jlev+1)-planck_hl(jspec,jlev)) / (od(jspec,jreg,jlev)*od(jspec,jreg,jlev)*(gamma1+gamma2))
          planck_hl_ad(jspec,jlev+1) = planck_hl_ad(jspec,jlev+1) + coeff_ad / (od(jspec,jreg,jlev)*(gamma1+gamma2))
          planck_hl_ad(jspec,jlev)   = planck_hl_ad(jspec,jlev)   - coeff_ad / (od(jspec,jreg,jlev)*(gamma1+gamma2))
          coeff_ad = 0.0_jprb

          ! Add direct planck contributions from coeff_* expressions
          planck_hl_ad(jspec,jlev)   = planck_hl_ad(jspec,jlev)   + ptop_ad
          planck_hl_ad(jspec,jlev+1) = planck_hl_ad(jspec,jlev+1) + pbase_ad

        else
          ! Low optical depth approximation
          ! reflectance = gamma2*od
          od_ad_loc = od_ad_loc + reflectance_ad(jspec,jreg,jlev) * gamma2
          refl_ad = refl_ad  ! already merged

          ! transmittance = (1-k*od) / (1 + od*(g1-k))
          numer = 1.0_jprb - k_exponent*od(jspec,jreg,jlev)
          denom = 1.0_jprb + od(jspec,jreg,jlev)*(gamma1-k_exponent)

          od_ad_loc = od_ad_loc + tran_ad * ( (-k_exponent)*denom - numer*(gamma1-k_exponent) ) / (denom*denom)
          tran_ad = 0.0_jprb

          ! source = (1 - refl - tran) * 0.5*(Ptop+Pbase)
          coeff = 0.5_jprb*(planck_hl(jspec,jlev)+planck_hl(jspec,jlev+1))
          tmp = (1.0_jprb - reflectance(jspec,jreg,jlev) - transmittance(jspec,jreg,jlev))
          planck_hl_ad(jspec,jlev)   = planck_hl_ad(jspec,jlev)   + (su_ad+sd_ad) * 0.5_jprb * tmp
          planck_hl_ad(jspec,jlev+1) = planck_hl_ad(jspec,jlev+1) + (su_ad+sd_ad) * 0.5_jprb * tmp
          refl_ad = refl_ad - (su_ad+sd_ad) * coeff
          tran_ad = tran_ad - (su_ad+sd_ad) * coeff
          su_ad = 0.0_jprb
          sd_ad = 0.0_jprb

          od_ad_loc = od_ad_loc + refl_ad * gamma2  ! refl = gamma2*od
          refl_ad = 0.0_jprb

          ! tran already applied above
        end if

        od_ad(jspec,jreg,jlev) = od_ad(jspec,jreg,jlev) + od_ad_loc

      end do
    end do

    ! -----------------------------
    ! Clear-sky region (jreg=1)
    ! -----------------------------
    jreg = 1
    do jspec = nspec,1,-1

      od_ad_loc = 0.0_jprb
      ptop_ad = 0.0_jprb
      pbase_ad = 0.0_jprb

      ! Unscale sources by region fraction
      su_ad = source_up_ad(jspec,1,jlev) * region_fracs(1,jlev)
      sd_ad = source_dn_ad(jspec,1,jlev) * region_fracs(1,jlev)

      region_fracs_ad(1,jlev) = region_fracs_ad(1,jlev) &
           &  + source_up_ad(jspec,1,jlev) * (source_up(jspec,1,jlev) / max(region_fracs(1,jlev), epsilon(1.0_jprb))) &
           &  + source_dn_ad(jspec,1,jlev) * (source_dn(jspec,1,jlev) / max(region_fracs(1,jlev), epsilon(1.0_jprb)))

      source_up_ad(jspec,1,jlev) = 0.0_jprb
      source_dn_ad(jspec,1,jlev) = 0.0_jprb

      tran_ad = transmittance_ad(jspec,1,jlev)
      transmittance_ad(jspec,1,jlev) = 0.0_jprb
      reflectance_ad(jspec,1,jlev) = 0.0_jprb ! reflectance is constant zero

      if (od(jspec,1,jlev) > OD_THRESH_2STREAM) then

        coeff = lw_diffusivity*od(jspec,1,jlev)
        exponential = exp(-coeff)
        ! transmittance = exp(-coeff)
        od_ad_loc = od_ad_loc + tran_ad * (-lw_diffusivity*exponential)
        tran_ad = 0.0_jprb

        ! coeff for linear Planck-with-OD formula
        coeff = (planck_hl(jspec,jlev+1)-planck_hl(jspec,jlev)) / (lw_diffusivity*od(jspec,1,jlev))
        coeff_up_top  =  coeff + planck_hl(jspec,jlev)
        coeff_up_base =  coeff + planck_hl(jspec,jlev+1)
        coeff_dn_top  = -coeff + planck_hl(jspec,jlev)
        coeff_dn_base = -coeff + planck_hl(jspec,jlev+1)

        coeff_ad = 0.0_jprb

        ! source_up = coeff_up_top - trans*coeff_up_base
        ptop_ad = ptop_ad + su_ad
        pbase_ad = pbase_ad - su_ad * exponential
        coeff_ad = coeff_ad + su_ad
        coeff_ad = coeff_ad - su_ad * exponential
        od_ad_loc = od_ad_loc + su_ad * (-coeff_up_base) * (-lw_diffusivity*exponential)  ! via trans inside source
        su_ad = 0.0_jprb

        ! source_dn = coeff_dn_base - trans*coeff_dn_top
        pbase_ad = pbase_ad + sd_ad
        ptop_ad  = ptop_ad - sd_ad * exponential
        coeff_ad = coeff_ad - sd_ad
        coeff_ad = coeff_ad + sd_ad * exponential
        od_ad_loc = od_ad_loc + sd_ad * (-coeff_dn_top) * (-lw_diffusivity*exponential)
        sd_ad = 0.0_jprb

        ! coeff = (Pbase-Ptop)/(lw_diffusivity*od)
        od_ad_loc = od_ad_loc - coeff_ad * (planck_hl(jspec,jlev+1)-planck_hl(jspec,jlev)) / (lw_diffusivity*od(jspec,1,jlev)*od(jspec,1,jlev))
        planck_hl_ad(jspec,jlev+1) = planck_hl_ad(jspec,jlev+1) + coeff_ad / (lw_diffusivity*od(jspec,1,jlev))
        planck_hl_ad(jspec,jlev)   = planck_hl_ad(jspec,jlev)   - coeff_ad / (lw_diffusivity*od(jspec,1,jlev))
        coeff_ad = 0.0_jprb

        planck_hl_ad(jspec,jlev)   = planck_hl_ad(jspec,jlev)   + ptop_ad
        planck_hl_ad(jspec,jlev+1) = planck_hl_ad(jspec,jlev+1) + pbase_ad

      else
        ! Low optical depth limit
        coeff = lw_diffusivity*od(jspec,1,jlev)
        ! trans = 1 - coeff
        od_ad_loc = od_ad_loc + tran_ad * (-lw_diffusivity)
        tran_ad = 0.0_jprb

        ! source = coeff*0.5*(Ptop+Pbase)
        planck_hl_ad(jspec,jlev)   = planck_hl_ad(jspec,jlev)   + (su_ad+sd_ad) * (coeff*0.5_jprb)
        planck_hl_ad(jspec,jlev+1) = planck_hl_ad(jspec,jlev+1) + (su_ad+sd_ad) * (coeff*0.5_jprb)
        od_ad_loc = od_ad_loc + (su_ad+sd_ad) * (0.5_jprb*(planck_hl(jspec,jlev)+planck_hl(jspec,jlev+1))) * lw_diffusivity
        su_ad = 0.0_jprb
        sd_ad = 0.0_jprb
      end if

      od_ad(jspec,1,jlev) = od_ad(jspec,1,jlev) + od_ad_loc

    end do

  end do ! jlev

  if (lhook) call dr_hook('tcrad:calc_reflectance_transmittance_ad',1,hook_handle)

end subroutine calc_reflectance_transmittance_ad


! ===== FILE: tcrad_calc_region_properties_ad.F90 =====
subroutine calc_region_properties_ad(nlev, &
     &  cloud_fraction, frac_std, do_gamma, &
     &  reg_fracs, od_scaling, cloud_fraction_threshold, &
     &  cloud_fraction_ad, frac_std_ad, &
     &  reg_fracs_ad, od_scaling_ad)

  use parkind1, only : jprb
  use yomhook,  only : lhook, dr_hook, jphook

  implicit none

  ! Minimum od_scaling in the case of a Gamma distribution
  real(jprb), parameter :: MinGammaODScaling = 0.025_jprb

  real(jprb), parameter :: MinLowerFrac      = 0.5_jprb
  real(jprb), parameter :: MaxLowerFrac      = 0.9_jprb
  real(jprb), parameter :: FSDAtMinLowerFrac = 1.5_jprb
  real(jprb), parameter :: FSDAtMaxLowerFrac = 3.725_jprb
  real(jprb), parameter :: LowerFracFSDGradient = (MaxLowerFrac-MinLowerFrac) / (FSDAtMaxLowerFrac-FSDAtMinLowerFrac)
  real(jprb), parameter :: LowerFracFSDIntercept = MinLowerFrac - FSDAtMinLowerFrac*LowerFracFSDGradient

  integer, intent(in) :: nlev
  logical, intent(in) :: do_gamma

  ! Nonlinear inputs
  real(jprb), intent(in), dimension(:)  :: cloud_fraction ! (nlev)
  real(jprb), intent(in), dimension(:)  :: frac_std       ! (nlev)

  ! Nonlinear outputs (trajectory provided)
  real(jprb), intent(in), dimension(NREGION,nlev) :: reg_fracs
  real(jprb), intent(in), dimension(2:NREGION,nlev) :: od_scaling

  real(jprb), intent(in), optional :: cloud_fraction_threshold

  ! Adjoint variables (accumulate)
  real(jprb), intent(inout), dimension(:) :: cloud_fraction_ad
  real(jprb), intent(inout), dimension(:) :: frac_std_ad
  real(jprb), intent(inout), dimension(NREGION,nlev) :: reg_fracs_ad
  real(jprb), intent(inout), dimension(2:NREGION,nlev) :: od_scaling_ad

  ! Locals
  real(jprb) :: frac_threshold
  integer :: jlev

  ! Lognormal branch locals
  real(jprb) :: x, logx, y, z, expmy
  real(jprb) :: x_ad, logx_ad, y_ad, z_ad, expmy_ad
  real(jprb) :: od2_ad

  ! Gamma branch locals
  real(jprb) :: lower, lower_clamped
  real(jprb) :: inner1, inner2, t, expt
  real(jprb) :: num, denom
  real(jprb) :: reg1_ad, reg2_ad, reg3_ad
  real(jprb) :: cf_ad, fsd_ad
  real(jprb) :: num_ad, denom_ad, lower_ad, lower_clamped_ad
  real(jprb) :: inner1_ad, inner2_ad, t_ad, expt_ad

  real(jphook) :: hook_handle

  if (lhook) call dr_hook('tcrad:calc_region_properties_ad',0,hook_handle)

  if (present(cloud_fraction_threshold)) then
    frac_threshold = cloud_fraction_threshold
  else
    frac_threshold = 10.0_jprb * epsilon(1.0_jprb)
  end if

  ! Reverse sweep over levels
  do jlev = nlev,1,-1

    if (cloud_fraction(jlev) < frac_threshold) then
      ! Outputs are constants: no sensitivities propagate
      reg_fracs_ad(:,jlev) = 0.0_jprb
      od_scaling_ad(2:NREGION,jlev) = 0.0_jprb
      cycle
    end if

    cf_ad  = 0.0_jprb
    fsd_ad = 0.0_jprb

    if (.not. do_gamma) then
      ! =========================
      ! Lognormal interpretation
      ! =========================

      ! reg_fracs:
      ! reg1 = 1-cf; reg2=reg3=0.5*cf
      cf_ad = cf_ad - reg_fracs_ad(1,jlev)
      cf_ad = cf_ad + 0.5_jprb*reg_fracs_ad(2,jlev) + 0.5_jprb*reg_fracs_ad(3,jlev)
      reg_fracs_ad(:,jlev) = 0.0_jprb

      ! od_scaling(3) = 2 - od_scaling(2)
      od2_ad = od_scaling_ad(2,jlev) - od_scaling_ad(3,jlev)
      od_scaling_ad(2,jlev) = 0.0_jprb
      od_scaling_ad(3,jlev) = 0.0_jprb

      ! Forward recompute intermediates
      x    = frac_std(jlev)*frac_std(jlev) + 1.0_jprb
      logx = log(x)
      y    = sqrt(logx)
      z    = sqrt(x)
      expmy = exp(-y)

      ! od2 = expmy / z
      expmy_ad = od2_ad / z
      z_ad     = -od2_ad * expmy / (z*z)
      od2_ad = 0.0_jprb

      ! expmy = exp(-y)
      y_ad = expmy_ad * (-expmy)
      expmy_ad = 0.0_jprb

      ! y = sqrt(logx)
      logx_ad = 0.0_jprb
      if (y > 0.0_jprb) then
        logx_ad = logx_ad + y_ad * (0.5_jprb / y)
      end if
      y_ad = 0.0_jprb

      ! logx = log(x)
      x_ad = logx_ad / x
      logx_ad = 0.0_jprb

      ! z = sqrt(x)
      if (z > 0.0_jprb) then
        x_ad = x_ad + z_ad * (0.5_jprb / z)
      end if
      z_ad = 0.0_jprb

      ! x = fsd^2 + 1
      fsd_ad = fsd_ad + x_ad * (2.0_jprb*frac_std(jlev))
      x_ad = 0.0_jprb

    else
      ! ======================
      ! Gamma interpretation
      ! ======================

      ! Pull reg_fracs adjoints into locals
      reg1_ad = reg_fracs_ad(1,jlev)
      reg2_ad = reg_fracs_ad(2,jlev)
      reg3_ad = reg_fracs_ad(3,jlev)
      reg_fracs_ad(:,jlev) = 0.0_jprb

      ! Pull od_scaling adjoints
      expt_ad = 0.0_jprb
      t_ad    = 0.0_jprb
      inner1_ad = 0.0_jprb
      inner2_ad = 0.0_jprb
      lower_clamped_ad = 0.0_jprb
      num_ad = 0.0_jprb
      denom_ad = 0.0_jprb

      ! od_scaling(3) = num/denom
      num   = cloud_fraction(jlev) - reg_fracs(2,jlev)*od_scaling(2,jlev)
      denom = reg_fracs(3,jlev)

      num_ad   = num_ad   + od_scaling_ad(3,jlev) / denom
      denom_ad = denom_ad - od_scaling_ad(3,jlev) * num / (denom*denom)
      od_scaling_ad(3,jlev) = 0.0_jprb

      ! denom = reg3
      reg3_ad = reg3_ad + denom_ad
      denom_ad = 0.0_jprb

      ! num = cf - reg2*od2
      cf_ad  = cf_ad + num_ad
      reg2_ad = reg2_ad - num_ad * od_scaling(2,jlev)
      od_scaling_ad(2,jlev) = od_scaling_ad(2,jlev) - num_ad * reg_fracs(2,jlev)
      num_ad = 0.0_jprb

      ! reg3 = 1 - reg1 - reg2
      reg1_ad = reg1_ad - reg3_ad
      reg2_ad = reg2_ad - reg3_ad
      reg3_ad = 0.0_jprb

      ! reg2 = cf * lower_clamped
      lower_clamped_ad = lower_clamped_ad + reg2_ad * cloud_fraction(jlev)
      cf_ad = cf_ad + reg2_ad * lower_clamped_from_fsd(frac_std(jlev))
      reg2_ad = 0.0_jprb

      ! reg1 = 1 - cf
      cf_ad = cf_ad - reg1_ad
      reg1_ad = 0.0_jprb

      ! od_scaling(2) = Min + (1-Min)*expt
      expt_ad = expt_ad + od_scaling_ad(2,jlev) * (1.0_jprb - MinGammaODScaling)
      od_scaling_ad(2,jlev) = 0.0_jprb

      ! expt = exp(t)
      lower = frac_std(jlev)  ! dummy to avoid unused
      ! recompute forward pieces for gamma branch
      inner1 = 1.0_jprb + 0.5_jprb*frac_std(jlev)
      inner2 = 1.0_jprb + 0.5_jprb*frac_std(jlev)*inner1
      t      = -frac_std(jlev) * inner2
      expt   = exp(t)

      t_ad = t_ad + expt_ad * expt
      expt_ad = 0.0_jprb

      ! t = -fsd*inner2
      fsd_ad = fsd_ad + t_ad * (-inner2)
      inner2_ad = inner2_ad + t_ad * (-frac_std(jlev))
      t_ad = 0.0_jprb

      ! inner2 = 1 + 0.5*fsd*inner1
      fsd_ad = fsd_ad + inner2_ad * (0.5_jprb*inner1)
      inner1_ad = inner1_ad + inner2_ad * (0.5_jprb*frac_std(jlev))
      inner2_ad = 0.0_jprb

      ! inner1 = 1 + 0.5*fsd
      fsd_ad = fsd_ad + inner1_ad * 0.5_jprb
      inner1_ad = 0.0_jprb

      ! lower_clamped = clamp(LowerFracFSDIntercept + fsd*LowerFracFSDGradient)
      lower = LowerFracFSDIntercept + frac_std(jlev)*LowerFracFSDGradient
      lower_clamped = max(MinLowerFrac, min(MaxLowerFrac, lower))
      lower_ad = 0.0_jprb
      if (lower > MinLowerFrac .and. lower < MaxLowerFrac) then
        lower_ad = lower_ad + lower_clamped_ad
      end if
      lower_clamped_ad = 0.0_jprb

      fsd_ad = fsd_ad + lower_ad * LowerFracFSDGradient
      lower_ad = 0.0_jprb

    end if

    cloud_fraction_ad(jlev) = cloud_fraction_ad(jlev) + cf_ad
    frac_std_ad(jlev)       = frac_std_ad(jlev)       + fsd_ad

  end do

  if (lhook) call dr_hook('tcrad:calc_region_properties_ad',1,hook_handle)

contains

  pure function lower_clamped_from_fsd(fsd) result(lc)
    use parkind1, only : jprb
    real(jprb), intent(in) :: fsd
    real(jprb) :: lc, lower
    lower = LowerFracFSDIntercept + fsd*LowerFracFSDGradient
    lc = max(MinLowerFrac, min(MaxLowerFrac, lower))
  end function lower_clamped_from_fsd

end subroutine calc_region_properties_ad


! ===== FILE: tcrad_calc_two_stream_flux_ad.F90 =====
subroutine calc_two_stream_flux_ad(nspec, nlev, &
     &  surf_emission, surf_albedo, &
     &  reflectance, transmittance, source_up, source_dn, &
     &  is_cloud_free_layer, u_overlap, v_overlap, &
     &  flux_up_base, flux_dn_base, flux_up_top, flux_dn_top, &
     &  surf_emission_ad, surf_albedo_ad, &
     &  reflectance_ad, transmittance_ad, source_up_ad, source_dn_ad, &
     &  u_overlap_ad, v_overlap_ad, &
     &  flux_up_base_ad, flux_dn_base_ad, flux_up_top_ad, flux_dn_top_ad)

  use parkind1, only           : jpim, jprb
  use yomhook,  only           : lhook, dr_hook, jphook

  implicit none

  ! Inputs (nonlinear trajectory)
  integer(jpim), intent(in) :: nspec, nlev
  real(jprb), intent(in),  dimension(nspec) :: surf_emission, surf_albedo
  real(jprb), intent(in),  dimension(nspec,NREGION,nlev) :: reflectance, transmittance
  real(jprb), intent(in),  dimension(nspec,NREGION,nlev) :: source_up, source_dn
  logical,    intent(in) :: is_cloud_free_layer(0:nlev+1)
  real(jprb), intent(in),  dimension(NREGION,NREGION,nlev+1) :: u_overlap, v_overlap
  real(jprb), intent(in),  dimension(nspec,NREGION,nlev) :: flux_up_base, flux_dn_base
  real(jprb), intent(in),  dimension(nspec,NREGION,nlev) :: flux_up_top,  flux_dn_top

  ! Adjoint variables (accumulate)
  real(jprb), intent(inout), dimension(nspec) :: surf_emission_ad, surf_albedo_ad
  real(jprb), intent(inout), dimension(nspec,NREGION,nlev) :: reflectance_ad, transmittance_ad
  real(jprb), intent(inout), dimension(nspec,NREGION,nlev) :: source_up_ad, source_dn_ad
  real(jprb), intent(inout), dimension(NREGION,NREGION,nlev+1) :: u_overlap_ad, v_overlap_ad

  real(jprb), intent(inout), dimension(nspec,NREGION,nlev) :: flux_up_base_ad, flux_dn_base_ad
  real(jprb), intent(inout), dimension(nspec,NREGION,nlev) :: flux_up_top_ad,  flux_dn_top_ad

  ! Forward recomputation storage needed for adjoint
  real(jprb), dimension(nspec, NREGION, nlev+1) :: total_albedo, total_source
  real(jprb), dimension(nspec, NREGION) :: total_albedo_below, total_source_below
  real(jprb), dimension(nspec, NREGION) :: inv_denom

  ! Adjoint storage for recomputed arrays
  real(jprb), dimension(nspec, NREGION, nlev+1) :: total_albedo_ad, total_source_ad
  real(jprb), dimension(nspec, NREGION) :: total_albedo_below_ad, total_source_below_ad

  integer(jpim) :: icloudtop
  integer(jpim) :: jspec, jlev, jreg, jreg2

  real(jprb) :: denom, numer, numer_ad, denom_ad
  real(jprb) :: t, a_next, inv, inv_ad, den, den_ad
  real(jprb) :: term, term_ad, t2

  real(jphook) :: hook_handle

  if (lhook) call dr_hook('tcrad:calc_two_stream_flux_ad',0,hook_handle)

  ! --------------------------------------------------------
  ! Forward recomputation of total_source/total_albedo (matches NL)
  ! --------------------------------------------------------
  icloudtop = nlev
  do jlev = 1,nlev
    if (.not. is_cloud_free_layer(jlev)) then
      icloudtop = jlev
      exit
    end if
  end do

  total_albedo = 0.0_jprb
  total_source = 0.0_jprb

  do jreg = 1,NREGION
    do jspec = 1,nspec
      total_source(jspec,jreg,nlev+1) = u_overlap(jreg,1,nlev+1)*surf_emission(jspec)
      total_albedo(jspec,jreg,nlev+1) = surf_albedo(jspec)
    end do
  end do

  do jlev = nlev,icloudtop,-1

    total_albedo_below = 0.0_jprb
    total_source_below = 0.0_jprb

    if (is_cloud_free_layer(jlev)) then
      do jspec = 1,nspec
        inv_denom(jspec,1) = 1.0_jprb / (1.0_jprb - total_albedo(jspec,1,jlev+1)*reflectance(jspec,1,jlev))
        total_albedo_below(jspec,1) = reflectance(jspec,1,jlev) &
             &  + transmittance(jspec,1,jlev)*transmittance(jspec,1,jlev)*total_albedo(jspec,1,jlev+1) &
             &  * inv_denom(jspec,1)
        total_source_below(jspec,1) = source_up(jspec,1,jlev) &
             &  + transmittance(jspec,1,jlev)*(total_source(jspec,1,jlev+1) &
             &  + total_albedo(jspec,1,jlev+1)*source_dn(jspec,1,jlev)) &
             &  * inv_denom(jspec,1)
      end do
    else
      do jreg = 1,NREGION
        do jspec = 1,nspec
          inv_denom(jspec,jreg) = 1.0_jprb / (1.0_jprb - total_albedo(jspec,jreg,jlev+1)*reflectance(jspec,jreg,jlev))
          total_albedo_below(jspec,jreg) = reflectance(jspec,jreg,jlev) &
               &  + transmittance(jspec,jreg,jlev)*transmittance(jspec,jreg,jlev)*total_albedo(jspec,jreg,jlev+1) &
               &  * inv_denom(jspec,jreg)
          total_source_below(jspec,jreg) = source_up(jspec,jreg,jlev) &
               &  + transmittance(jspec,jreg,jlev)*(total_source(jspec,jreg,jlev+1) &
               &  + total_albedo(jspec,jreg,jlev+1)*source_dn(jspec,jreg,jlev)) &
               &  * inv_denom(jspec,jreg)
        end do
      end do
    end if

    if (is_cloud_free_layer(jlev) .and. is_cloud_free_layer(jlev-1)) then
      total_albedo(:,:,jlev) = total_albedo_below(:,:)
      total_source(:,:,jlev) = total_source_below(:,:)
    else
      total_source(:,:,jlev) = 0.0_jprb
      do jspec = 1,nspec
        do jreg = 1,NREGION
          do jreg2 = 1,NREGION
            total_source(jspec,jreg,jlev) = total_source(jspec,jreg,jlev) &
                 &  + u_overlap(jreg,jreg2,jlev) * total_source_below(jspec,jreg2)
          end do
        end do
      end do

      total_albedo(:,:,jlev) = 0.0_jprb
      do jreg = 1,NREGION
        do jreg2 = 1,NREGION
          total_albedo(:,jreg,jlev) = total_albedo(:,jreg,jlev) &
               &  + total_albedo_below(:,jreg2) * v_overlap(jreg2,jreg,jlev)
        end do
      end do
    end if
  end do

  ! --------------------------------------------------------
  ! Initialise adjoint storage for recomputed arrays
  ! --------------------------------------------------------
  total_albedo_ad = 0.0_jprb
  total_source_ad = 0.0_jprb
  total_albedo_below_ad = 0.0_jprb
  total_source_below_ad = 0.0_jprb

  ! --------------------------------------------------------
  ! Section 5 adjoint: loop down from nlev..icloudtop
  ! --------------------------------------------------------
  do jlev = nlev,icloudtop,-1

    ! Overlap mapping for flux_dn_top(:,:,jlev+1)
    if (jlev < nlev) then
      if (.not. (is_cloud_free_layer(jlev) .and. is_cloud_free_layer(jlev+1))) then
        do jspec = 1,nspec
          do jreg = 1,NREGION
            do jreg2 = 1,NREGION
              v_overlap_ad(jreg,jreg2,jlev+1) = v_overlap_ad(jreg,jreg2,jlev+1) &
                   &  + flux_dn_top_ad(jspec,jreg,jlev+1) * flux_dn_base(jspec,jreg2,jlev)
              flux_dn_base_ad(jspec,jreg2,jlev) = flux_dn_base_ad(jspec,jreg2,jlev) &
                   &  + v_overlap(jreg,jreg2,jlev+1) * flux_dn_top_ad(jspec,jreg,jlev+1)
            end do
          end do
        end do
        flux_dn_top_ad(:,:,jlev+1) = 0.0_jprb
      else
        flux_dn_base_ad(:,1,jlev) = flux_dn_base_ad(:,1,jlev) + flux_dn_top_ad(:,1,jlev+1)
        flux_dn_top_ad(:,1,jlev+1) = 0.0_jprb
      end if
    end if

    if (is_cloud_free_layer(jlev)) then
      do jspec = 1,nspec
        ! flux_up_top = flux_up_base*T + source_up + flux_dn_top*R
        transmittance_ad(jspec,1,jlev) = transmittance_ad(jspec,1,jlev) &
             &  + flux_up_top_ad(jspec,1,jlev) * flux_up_base(jspec,1,jlev)
        flux_up_base_ad(jspec,1,jlev) = flux_up_base_ad(jspec,1,jlev) &
             &  + flux_up_top_ad(jspec,1,jlev) * transmittance(jspec,1,jlev)
        source_up_ad(jspec,1,jlev) = source_up_ad(jspec,1,jlev) + flux_up_top_ad(jspec,1,jlev)
        reflectance_ad(jspec,1,jlev) = reflectance_ad(jspec,1,jlev) &
             &  + flux_up_top_ad(jspec,1,jlev) * flux_dn_top(jspec,1,jlev)
        flux_dn_top_ad(jspec,1,jlev) = flux_dn_top_ad(jspec,1,jlev) &
             &  + flux_up_top_ad(jspec,1,jlev) * reflectance(jspec,1,jlev)
        flux_up_top_ad(jspec,1,jlev) = 0.0_jprb

        ! flux_up_base = total_source_next + flux_dn_base * total_albedo_next
        total_source_ad(jspec,1,jlev+1) = total_source_ad(jspec,1,jlev+1) + flux_up_base_ad(jspec,1,jlev)
        total_albedo_ad(jspec,1,jlev+1) = total_albedo_ad(jspec,1,jlev+1) &
             &  + flux_up_base_ad(jspec,1,jlev) * flux_dn_base(jspec,1,jlev)
        flux_dn_base_ad(jspec,1,jlev) = flux_dn_base_ad(jspec,1,jlev) &
             &  + flux_up_base_ad(jspec,1,jlev) * total_albedo(jspec,1,jlev+1)
        flux_up_base_ad(jspec,1,jlev) = 0.0_jprb

        ! flux_dn_base = numer/denom
        denom = 1.0_jprb - reflectance(jspec,1,jlev)*total_albedo(jspec,1,jlev+1)
        numer = transmittance(jspec,1,jlev)*flux_dn_top(jspec,1,jlev) &
             &  + reflectance(jspec,1,jlev)*total_source(jspec,1,jlev+1) &
             &  + source_dn(jspec,1,jlev)

        numer_ad = flux_dn_base_ad(jspec,1,jlev) / denom
        denom_ad = -flux_dn_base_ad(jspec,1,jlev) * numer / (denom*denom)
        flux_dn_base_ad(jspec,1,jlev) = 0.0_jprb

        ! denom = 1 - R*A
        reflectance_ad(jspec,1,jlev) = reflectance_ad(jspec,1,jlev) - denom_ad * total_albedo(jspec,1,jlev+1)
        total_albedo_ad(jspec,1,jlev+1) = total_albedo_ad(jspec,1,jlev+1) - denom_ad * reflectance(jspec,1,jlev)

        ! numer = T*Ftop + R*tsrc + src_dn
        transmittance_ad(jspec,1,jlev) = transmittance_ad(jspec,1,jlev) + numer_ad * flux_dn_top(jspec,1,jlev)
        flux_dn_top_ad(jspec,1,jlev) = flux_dn_top_ad(jspec,1,jlev) + numer_ad * transmittance(jspec,1,jlev)
        reflectance_ad(jspec,1,jlev) = reflectance_ad(jspec,1,jlev) + numer_ad * total_source(jspec,1,jlev+1)
        total_source_ad(jspec,1,jlev+1) = total_source_ad(jspec,1,jlev+1) + numer_ad * reflectance(jspec,1,jlev)
        source_dn_ad(jspec,1,jlev) = source_dn_ad(jspec,1,jlev) + numer_ad
      end do
    else
      do jreg = 1,NREGION
        do jspec = 1,nspec
          ! flux_up_top
          transmittance_ad(jspec,jreg,jlev) = transmittance_ad(jspec,jreg,jlev) &
               &  + flux_up_top_ad(jspec,jreg,jlev) * flux_up_base(jspec,jreg,jlev)
          flux_up_base_ad(jspec,jreg,jlev) = flux_up_base_ad(jspec,jreg,jlev) &
               &  + flux_up_top_ad(jspec,jreg,jlev) * transmittance(jspec,jreg,jlev)
          source_up_ad(jspec,jreg,jlev) = source_up_ad(jspec,jreg,jlev) + flux_up_top_ad(jspec,jreg,jlev)
          reflectance_ad(jspec,jreg,jlev) = reflectance_ad(jspec,jreg,jlev) &
               &  + flux_up_top_ad(jspec,jreg,jlev) * flux_dn_top(jspec,jreg,jlev)
          flux_dn_top_ad(jspec,jreg,jlev) = flux_dn_top_ad(jspec,jreg,jlev) &
               &  + flux_up_top_ad(jspec,jreg,jlev) * reflectance(jspec,jreg,jlev)
          flux_up_top_ad(jspec,jreg,jlev) = 0.0_jprb

          ! flux_up_base
          total_source_ad(jspec,jreg,jlev+1) = total_source_ad(jspec,jreg,jlev+1) + flux_up_base_ad(jspec,jreg,jlev)
          total_albedo_ad(jspec,jreg,jlev+1) = total_albedo_ad(jspec,jreg,jlev+1) &
               &  + flux_up_base_ad(jspec,jreg,jlev) * flux_dn_base(jspec,jreg,jlev)
          flux_dn_base_ad(jspec,jreg,jlev) = flux_dn_base_ad(jspec,jreg,jlev) &
               &  + flux_up_base_ad(jspec,jreg,jlev) * total_albedo(jspec,jreg,jlev+1)
          flux_up_base_ad(jspec,jreg,jlev) = 0.0_jprb

          ! flux_dn_base division
          denom = 1.0_jprb - reflectance(jspec,jreg,jlev)*total_albedo(jspec,jreg,jlev+1)
          numer = transmittance(jspec,jreg,jlev)*flux_dn_top(jspec,jreg,jlev) &
               &  + reflectance(jspec,jreg,jlev)*total_source(jspec,jreg,jlev+1) &
               &  + source_dn(jspec,jreg,jlev)

          numer_ad = flux_dn_base_ad(jspec,jreg,jlev) / denom
          denom_ad = -flux_dn_base_ad(jspec,jreg,jlev) * numer / (denom*denom)
          flux_dn_base_ad(jspec,jreg,jlev) = 0.0_jprb

          reflectance_ad(jspec,jreg,jlev) = reflectance_ad(jspec,jreg,jlev) - denom_ad * total_albedo(jspec,jreg,jlev+1)
          total_albedo_ad(jspec,jreg,jlev+1) = total_albedo_ad(jspec,jreg,jlev+1) - denom_ad * reflectance(jspec,jreg,jlev)

          transmittance_ad(jspec,jreg,jlev) = transmittance_ad(jspec,jreg,jlev) + numer_ad * flux_dn_top(jspec,jreg,jlev)
          flux_dn_top_ad(jspec,jreg,jlev) = flux_dn_top_ad(jspec,jreg,jlev) + numer_ad * transmittance(jspec,jreg,jlev)
          reflectance_ad(jspec,jreg,jlev) = reflectance_ad(jspec,jreg,jlev) + numer_ad * total_source(jspec,jreg,jlev+1)
          total_source_ad(jspec,jreg,jlev+1) = total_source_ad(jspec,jreg,jlev+1) + numer_ad * reflectance(jspec,jreg,jlev)
          source_dn_ad(jspec,jreg,jlev) = source_dn_ad(jspec,jreg,jlev) + numer_ad
        end do
      end do
    end if

  end do

  ! --------------------------------------------------------
  ! Section 4 adjoint (clear-sky upwelling to TOA)
  ! --------------------------------------------------------
  if (icloudtop > 1) then
    do jlev = 1,icloudtop-2
      flux_up_top_ad(:,1,jlev+1) = flux_up_top_ad(:,1,jlev+1) + flux_up_base_ad(:,1,jlev)
      flux_up_base_ad(:,1,jlev) = 0.0_jprb

      source_up_ad(:,1,jlev) = source_up_ad(:,1,jlev) + flux_up_top_ad(:,1,jlev)
      transmittance_ad(:,1,jlev) = transmittance_ad(:,1,jlev) + flux_up_top_ad(:,1,jlev) * flux_up_base(:,1,jlev)
      flux_up_base_ad(:,1,jlev) = flux_up_base_ad(:,1,jlev) + flux_up_top_ad(:,1,jlev) * transmittance(:,1,jlev)
      flux_up_top_ad(:,1,jlev) = 0.0_jprb
    end do

    source_up_ad(:,1,icloudtop-1) = source_up_ad(:,1,icloudtop-1) + flux_up_top_ad(:,1,icloudtop-1)
    transmittance_ad(:,1,icloudtop-1) = transmittance_ad(:,1,icloudtop-1) &
         &  + flux_up_top_ad(:,1,icloudtop-1) * flux_up_base(:,1,icloudtop-1)
    flux_up_base_ad(:,1,icloudtop-1) = flux_up_base_ad(:,1,icloudtop-1) &
         &  + flux_up_top_ad(:,1,icloudtop-1) * transmittance(:,1,icloudtop-1)
    flux_up_top_ad(:,1,icloudtop-1) = 0.0_jprb

    total_source_ad(:,1,icloudtop) = total_source_ad(:,1,icloudtop) + flux_up_base_ad(:,1,icloudtop-1)
    total_albedo_ad(:,1,icloudtop) = total_albedo_ad(:,1,icloudtop) + flux_up_base_ad(:,1,icloudtop-1) * flux_dn_base(:,1,icloudtop-1)
    flux_dn_base_ad(:,1,icloudtop-1) = flux_dn_base_ad(:,1,icloudtop-1) + flux_up_base_ad(:,1,icloudtop-1) * total_albedo(:,1,icloudtop)
    flux_up_base_ad(:,1,icloudtop-1) = 0.0_jprb
  end if

  ! --------------------------------------------------------
  ! Section 3 adjoint: reverse the adding-method recursion
  ! --------------------------------------------------------
  do jlev = icloudtop, nlev

    ! Map adjoints from total_*(:,:,jlev) back to total_*_below
    total_albedo_below_ad(:,:) = 0.0_jprb
    total_source_below_ad(:,:) = 0.0_jprb

    if (is_cloud_free_layer(jlev) .and. is_cloud_free_layer(jlev-1)) then
      total_albedo_below_ad(:,:) = total_albedo_below_ad(:,:) + total_albedo_ad(:,:,jlev)
      total_source_below_ad(:,:) = total_source_below_ad(:,:) + total_source_ad(:,:,jlev)
      total_albedo_ad(:,:,jlev) = 0.0_jprb
      total_source_ad(:,:,jlev) = 0.0_jprb
    else
      ! total_source(:,:,jlev) = u_overlap(:,:,jlev) * total_source_below
      do jspec = 1,nspec
        do jreg = 1,NREGION
          do jreg2 = 1,NREGION
            u_overlap_ad(jreg,jreg2,jlev) = u_overlap_ad(jreg,jreg2,jlev) &
                 &  + total_source_ad(jspec,jreg,jlev) * total_source_below(jspec,jreg2)
            total_source_below_ad(jspec,jreg2) = total_source_below_ad(jspec,jreg2) &
                 &  + u_overlap(jreg,jreg2,jlev) * total_source_ad(jspec,jreg,jlev)
          end do
        end do
      end do
      total_source_ad(:,:,jlev) = 0.0_jprb

      ! total_albedo(:,jreg,jlev) = sum total_albedo_below(:,jreg2) * v_overlap(jreg2,jreg,jlev)
      do jreg = 1,NREGION
        do jreg2 = 1,NREGION
          do jspec = 1,nspec
            v_overlap_ad(jreg2,jreg,jlev) = v_overlap_ad(jreg2,jreg,jlev) &
                 &  + total_albedo_ad(jspec,jreg,jlev) * total_albedo_below(jspec,jreg2)
          end do
          total_albedo_below_ad(:,jreg2) = total_albedo_below_ad(:,jreg2) &
               &  + v_overlap(jreg2,jreg,jlev) * total_albedo_ad(:,jreg,jlev)
        end do
      end do
      total_albedo_ad(:,:,jlev) = 0.0_jprb
    end if

    ! Invert layer-local adding method: compute adjoints of reflectance/transmittance/source_*/total_* at jlev+1
    if (is_cloud_free_layer(jlev)) then
      jreg = 1
      do jspec = 1,nspec
        inv = inv_denom(jspec,jreg)
        inv_ad = 0.0_jprb

        ! total_albedo_below = R + T^2*A_next*inv
        reflectance_ad(jspec,jreg,jlev) = reflectance_ad(jspec,jreg,jlev) + total_albedo_below_ad(jspec,jreg)
        t = transmittance(jspec,jreg,jlev)
        t2 = t*t
        a_next = total_albedo(jspec,jreg,jlev+1)
        transmittance_ad(jspec,jreg,jlev) = transmittance_ad(jspec,jreg,jlev) &
             &  + total_albedo_below_ad(jspec,jreg) * (2.0_jprb*t*a_next*inv)
        total_albedo_ad(jspec,jreg,jlev+1) = total_albedo_ad(jspec,jreg,jlev+1) &
             &  + total_albedo_below_ad(jspec,jreg) * (t2*inv)
        inv_ad = inv_ad + total_albedo_below_ad(jspec,jreg) * (t2*a_next)
        total_albedo_below_ad(jspec,jreg) = 0.0_jprb

        ! total_source_below = source_up + T*(S_next + A_next*source_dn)*inv
        source_up_ad(jspec,jreg,jlev) = source_up_ad(jspec,jreg,jlev) + total_source_below_ad(jspec,jreg)

        term = total_source(jspec,jreg,jlev+1) + total_albedo(jspec,jreg,jlev+1)*source_dn(jspec,jreg,jlev)
        transmittance_ad(jspec,jreg,jlev) = transmittance_ad(jspec,jreg,jlev) &
             &  + total_source_below_ad(jspec,jreg) * (term*inv)

        term_ad = total_source_below_ad(jspec,jreg) * (t*inv)
        inv_ad = inv_ad + total_source_below_ad(jspec,jreg) * (t*term)

        total_source_ad(jspec,jreg,jlev+1) = total_source_ad(jspec,jreg,jlev+1) + term_ad
        total_albedo_ad(jspec,jreg,jlev+1) = total_albedo_ad(jspec,jreg,jlev+1) + term_ad * source_dn(jspec,jreg,jlev)
        source_dn_ad(jspec,jreg,jlev) = source_dn_ad(jspec,jreg,jlev) + term_ad * total_albedo(jspec,jreg,jlev+1)

        total_source_below_ad(jspec,jreg) = 0.0_jprb

        ! inv = 1/(1 - A_next*R)
        den = 1.0_jprb - total_albedo(jspec,jreg,jlev+1)*reflectance(jspec,jreg,jlev)
        den_ad = -inv_ad * (inv*inv)
        reflectance_ad(jspec,jreg,jlev) = reflectance_ad(jspec,jreg,jlev) - den_ad * total_albedo(jspec,jreg,jlev+1)
        total_albedo_ad(jspec,jreg,jlev+1) = total_albedo_ad(jspec,jreg,jlev+1) - den_ad * reflectance(jspec,jreg,jlev)
      end do
    else
      do jreg = 1,NREGION
        do jspec = 1,nspec
          inv = inv_denom(jspec,jreg)
          inv_ad = 0.0_jprb

          reflectance_ad(jspec,jreg,jlev) = reflectance_ad(jspec,jreg,jlev) + total_albedo_below_ad(jspec,jreg)
          t = transmittance(jspec,jreg,jlev)
          t2 = t*t
          a_next = total_albedo(jspec,jreg,jlev+1)
          transmittance_ad(jspec,jreg,jlev) = transmittance_ad(jspec,jreg,jlev) &
               &  + total_albedo_below_ad(jspec,jreg) * (2.0_jprb*t*a_next*inv)
          total_albedo_ad(jspec,jreg,jlev+1) = total_albedo_ad(jspec,jreg,jlev+1) &
               &  + total_albedo_below_ad(jspec,jreg) * (t2*inv)
          inv_ad = inv_ad + total_albedo_below_ad(jspec,jreg) * (t2*a_next)
          total_albedo_below_ad(jspec,jreg) = 0.0_jprb

          source_up_ad(jspec,jreg,jlev) = source_up_ad(jspec,jreg,jlev) + total_source_below_ad(jspec,jreg)

          term = total_source(jspec,jreg,jlev+1) + total_albedo(jspec,jreg,jlev+1)*source_dn(jspec,jreg,jlev)
          transmittance_ad(jspec,jreg,jlev) = transmittance_ad(jspec,jreg,jlev) &
               &  + total_source_below_ad(jspec,jreg) * (term*inv)

          term_ad = total_source_below_ad(jspec,jreg) * (t*inv)
          inv_ad = inv_ad + total_source_below_ad(jspec,jreg) * (t*term)

          total_source_ad(jspec,jreg,jlev+1) = total_source_ad(jspec,jreg,jlev+1) + term_ad
          total_albedo_ad(jspec,jreg,jlev+1) = total_albedo_ad(jspec,jreg,jlev+1) + term_ad * source_dn(jspec,jreg,jlev)
          source_dn_ad(jspec,jreg,jlev) = source_dn_ad(jspec,jreg,jlev) + term_ad * total_albedo(jspec,jreg,jlev+1)

          total_source_below_ad(jspec,jreg) = 0.0_jprb

          den = 1.0_jprb - total_albedo(jspec,jreg,jlev+1)*reflectance(jspec,jreg,jlev)
          den_ad = -inv_ad * (inv*inv)
          reflectance_ad(jspec,jreg,jlev) = reflectance_ad(jspec,jreg,jlev) - den_ad * total_albedo(jspec,jreg,jlev+1)
          total_albedo_ad(jspec,jreg,jlev+1) = total_albedo_ad(jspec,jreg,jlev+1) - den_ad * reflectance(jspec,jreg,jlev)
        end do
      end do
    end if

  end do

  ! --------------------------------------------------------
  ! Surface boundary adjoint from total_source/total_albedo at nlev+1
  ! --------------------------------------------------------
  do jspec = 1,nspec
    surf_albedo_ad(jspec) = surf_albedo_ad(jspec) + sum(total_albedo_ad(jspec,:,nlev+1))
    total_albedo_ad(jspec,:,nlev+1) = 0.0_jprb

    do jreg = 1,NREGION
      u_overlap_ad(jreg,1,nlev+1) = u_overlap_ad(jreg,1,nlev+1) + total_source_ad(jspec,jreg,nlev+1) * surf_emission(jspec)
      surf_emission_ad(jspec) = surf_emission_ad(jspec) + total_source_ad(jspec,jreg,nlev+1) * u_overlap(jreg,1,nlev+1)
    end do
    total_source_ad(jspec,:,nlev+1) = 0.0_jprb
  end do

  if (lhook) call dr_hook('tcrad:calc_two_stream_flux_ad',1,hook_handle)

end subroutine calc_two_stream_flux_ad



! ===== ADDED: calc_overlap_matrices_ad (and helper calc_alpha_overlap_matrix_ad) =====
subroutine calc_overlap_matrices_ad(nlev, &
     &     region_fracs, overlap_param, &
     &     u_overlap, v_overlap, &
     &     decorrelation_scaling, &
     &     cloud_fraction_threshold, cloud_cover, &
     &     region_fracs_ad, overlap_param_ad, &
     &     u_overlap_ad, v_overlap_ad, &
     &     cloud_cover_ad)

  use parkind1,     only : jprb, jpim
  use yomhook, only : lhook, dr_hook, jphook

  implicit none

  integer(jpim),  intent(in) :: nlev

  ! Nonlinear inputs
  real(jprb), intent(in), dimension(1:NREGION,nlev) :: region_fracs
  real(jprb), intent(in), dimension(:)              :: overlap_param

  ! Nonlinear outputs (trajectory provided)
  real(jprb), intent(in),  dimension(NREGION,NREGION,nlev+1) :: u_overlap, v_overlap

  ! Optional config (NOT adjointed)
  real(jprb), intent(in), optional :: decorrelation_scaling
  real(jprb), intent(in), optional :: cloud_fraction_threshold

  ! Optional nonlinear output
  real(jprb), intent(in), optional :: cloud_cover

  ! Adjoint variables (accumulate)
  real(jprb), intent(inout), dimension(1:NREGION,nlev) :: region_fracs_ad
  real(jprb), intent(inout), dimension(:)              :: overlap_param_ad
  real(jprb), intent(inout), dimension(NREGION,NREGION,nlev+1) :: u_overlap_ad, v_overlap_ad

  ! Optional adjoint input
  real(jprb), intent(inout), optional :: cloud_cover_ad

  ! Locals
  integer(jpim) :: jlev, jupper, jlower, jreg, jreg2
  real(jprb) :: frac_threshold
  real(jprb) :: used_decorrelation_scaling
  real(jprb) :: pexp

  real(jprb) :: frac_upper(NREGION), frac_lower(NREGION)
  real(jprb) :: op1, op2
  real(jprb) :: overlap_matrix(NREGION,NREGION)

  ! Stored forward state for reverse sweep
  real(jprb), allocatable :: frac_upper_store(:,:), frac_lower_store(:,:)
  real(jprb), allocatable :: op1_store(:), op2_store(:)
  real(jprb), allocatable :: overlap_store(:,:,:)

  ! Adjoint working
  real(jprb) :: frac_upper_ad_next(NREGION)
  real(jprb) :: frac_upper_ad(NREGION), frac_lower_ad(NREGION)
  real(jprb) :: overlap_ad(NREGION,NREGION)
  real(jprb) :: op1_ad, op2_ad

  ! Cloud-cover product bookkeeping (optional)
  real(jprb), allocatable :: prod_prefix(:)
  real(jprb) :: prod_v11, prod_ad, v11, epsv

  real(jphook) :: hook_handle

  if (lhook) call dr_hook('tcrad:calc_overlap_matrices_ad',0,hook_handle)

  if (present(decorrelation_scaling)) then
    used_decorrelation_scaling = decorrelation_scaling
  else
    used_decorrelation_scaling = 1.0_jprb
  end if
  pexp = 1.0_jprb/used_decorrelation_scaling

  if (present(cloud_fraction_threshold)) then
    frac_threshold = cloud_fraction_threshold
  else
    frac_threshold = 10.0_jprb * epsilon(1.0_jprb)
  end if

  allocate(frac_upper_store(NREGION,nlev+1))
  allocate(frac_lower_store(NREGION,nlev+1))
  allocate(op1_store(nlev+1))
  allocate(op2_store(nlev+1))
  allocate(overlap_store(NREGION,NREGION,nlev+1))

  ! Optional: prefix products for cloud_cover derivative
  allocate(prod_prefix(nlev+2))
  prod_prefix(1) = 1.0_jprb
  do jlev = 1,nlev+1
    prod_prefix(jlev+1) = prod_prefix(jlev) * v_overlap(1,1,jlev)
  end do
  prod_v11 = prod_prefix(nlev+2)

  ! --------------------------------------------------------
  ! Forward recomputation storing intermediates
  ! --------------------------------------------------------
  frac_upper(1) = 1.0_jprb
  frac_upper(2:NREGION) = 0.0_jprb

  do jlev = 1,nlev+1

    if (jlev > nlev) then
      frac_lower(1) = 1.0_jprb
      frac_lower(2:NREGION) = 0.0_jprb
    else
      frac_lower(:) = region_fracs(:,jlev)
    end if

    if (jlev == 1 .or. jlev > nlev) then
      op1 = 1.0_jprb
      op2 = 1.0_jprb
    else
      op1 = overlap_param(jlev-1)
      if (op1 >= 0.0_jprb) then
        op2 = op1**pexp
      else
        op2 = op1
      end if
    end if

    call calc_alpha_overlap_matrix_fwd(op1, op2, frac_upper, frac_lower, overlap_matrix)

    frac_upper_store(:,jlev) = frac_upper(:)
    frac_lower_store(:,jlev) = frac_lower(:)
    op1_store(jlev) = op1
    op2_store(jlev) = op2
    overlap_store(:,:,jlev) = overlap_matrix(:,:)

    frac_upper(:) = frac_lower(:)
  end do

  ! --------------------------------------------------------
  ! Cloud cover adjoint -> v_overlap_ad(1,1,:)
  ! cloud_cover = 1 - prod(v_overlap(1,1,:))
  ! --------------------------------------------------------
  if (present(cloud_cover_ad) .and. present(cloud_cover)) then
    prod_ad = -cloud_cover_ad
    epsv = 1.0e-14_jprb
    do jlev = 1,nlev+1
      v11 = v_overlap(1,1,jlev)
      if (abs(v11) > epsv) then
        v_overlap_ad(1,1,jlev) = v_overlap_ad(1,1,jlev) + prod_ad * (prod_v11 / v11)
      end if
    end do
    cloud_cover_ad = 0.0_jprb
  end if

  ! --------------------------------------------------------
  ! Reverse sweep
  ! --------------------------------------------------------
  frac_upper_ad_next(:) = 0.0_jprb

  do jlev = nlev+1,1,-1

    frac_upper(:) = frac_upper_store(:,jlev)
    frac_lower(:) = frac_lower_store(:,jlev)
    op1 = op1_store(jlev)
    op2 = op2_store(jlev)
    overlap_matrix(:,:) = overlap_store(:,:,jlev)

    ! frac_upper(next) = frac_lower(this)
    frac_lower_ad(:) = frac_upper_ad_next(:)
    frac_upper_ad(:) = 0.0_jprb
    overlap_ad(:,:)  = 0.0_jprb

    ! u_overlap(jupper,jlower,jlev) = overlap(jupper,jlower)/frac_lower(jlower)
    do jupper = 1,NREGION
      do jlower = 1,NREGION
        if (frac_lower(jlower) >= frac_threshold) then
          overlap_ad(jupper,jlower) = overlap_ad(jupper,jlower) + u_overlap_ad(jupper,jlower,jlev) / frac_lower(jlower)
          frac_lower_ad(jlower) = frac_lower_ad(jlower) - u_overlap_ad(jupper,jlower,jlev) * overlap_matrix(jupper,jlower) &
               & / (frac_lower(jlower)*frac_lower(jlower))
        end if
        u_overlap_ad(jupper,jlower,jlev) = 0.0_jprb
      end do
    end do

    ! v_overlap(jlower,jupper,jlev) = overlap(jupper,jlower)/frac_upper(jupper)
    do jupper = 1,NREGION
      do jlower = 1,NREGION
        if (frac_upper(jupper) >= frac_threshold) then
          overlap_ad(jupper,jlower) = overlap_ad(jupper,jlower) + v_overlap_ad(jlower,jupper,jlev) / frac_upper(jupper)
          frac_upper_ad(jupper) = frac_upper_ad(jupper) - v_overlap_ad(jlower,jupper,jlev) * overlap_matrix(jupper,jlower) &
               & / (frac_upper(jupper)*frac_upper(jupper))
        end if
        v_overlap_ad(jlower,jupper,jlev) = 0.0_jprb
      end do
    end do

    ! Backprop through overlap-matrix construction
    op1_ad = 0.0_jprb
    op2_ad = 0.0_jprb
    call calc_alpha_overlap_matrix_ad(op1, op2, frac_upper, frac_lower, overlap_matrix, overlap_ad, &
         &  op1_ad, op2_ad, frac_upper_ad, frac_lower_ad)

    ! Map overlap_param adjoints (note: pexp treated constant)
    if (.not. (jlev == 1 .or. jlev > nlev)) then
      overlap_param_ad(jlev-1) = overlap_param_ad(jlev-1) + op1_ad
      if (op1 >= 0.0_jprb) then
        overlap_param_ad(jlev-1) = overlap_param_ad(jlev-1) + op2_ad * (pexp * op1**(pexp-1.0_jprb))
      else
        overlap_param_ad(jlev-1) = overlap_param_ad(jlev-1) + op2_ad
      end if
    end if

    ! Map frac_lower adjoint to region_fracs
    if (jlev <= nlev) then
      region_fracs_ad(:,jlev) = region_fracs_ad(:,jlev) + frac_lower_ad(:)
    end if

    ! Prepare for previous iteration
    frac_upper_ad_next(:) = frac_upper_ad(:)

  end do

  deallocate(frac_upper_store, frac_lower_store, op1_store, op2_store, overlap_store, prod_prefix)

  if (lhook) call dr_hook('tcrad:calc_overlap_matrices_ad',1,hook_handle)

contains

  subroutine calc_alpha_overlap_matrix_fwd(op, op_inhom, frac_upper, frac_lower, overlap_matrix)
    use parkind1, only : jprb
    real(jprb), intent(in) :: op, op_inhom
    real(jprb), intent(in), dimension(NREGION) :: frac_upper, frac_lower
    real(jprb), intent(out) :: overlap_matrix(NREGION,NREGION)

    real(jprb) :: pair_cloud_cover, cf_upper, cf_lower, one_over_cf, frac_both
    real(jprb) :: max_cf, tmp, denom
    real(jprb) :: cf_u_eff, cf_l_eff
    real(jprb), parameter :: CF_MIN = 1.0e-6_jprb

    overlap_matrix(:,:) = 0.0_jprb

    cf_upper = sum(frac_upper(2:NREGION))
    cf_lower = sum(frac_lower(2:NREGION))

    max_cf = max(cf_upper, cf_lower)
    tmp = (cf_upper+cf_lower-cf_upper*cf_lower)
    pair_cloud_cover = op*max_cf + (1.0_jprb-op)*tmp

    overlap_matrix(1,1) = 1.0_jprb - pair_cloud_cover

    denom = max(cf_lower, CF_MIN)
    one_over_cf = 1.0_jprb/denom
    tmp = (pair_cloud_cover - cf_upper)
    overlap_matrix(1,2) = tmp * frac_lower(2) * one_over_cf
    overlap_matrix(1,3) = tmp * frac_lower(3) * one_over_cf

    denom = max(cf_upper, CF_MIN)
    one_over_cf = 1.0_jprb/denom
    tmp = (pair_cloud_cover - cf_lower)
    overlap_matrix(2,1) = tmp * frac_upper(2) * one_over_cf
    overlap_matrix(3,1) = tmp * frac_upper(3) * one_over_cf

    frac_both = cf_upper + cf_lower - pair_cloud_cover

    cf_u_eff = max(cf_upper, CF_MIN)
    cf_l_eff = max(cf_lower, CF_MIN)
    cf_upper = frac_upper(3) / cf_u_eff
    cf_lower = frac_lower(3) / cf_l_eff

    max_cf = max(cf_upper, cf_lower)
    tmp = (cf_upper+cf_lower-cf_upper*cf_lower)
    pair_cloud_cover = op_inhom*max_cf + (1.0_jprb-op_inhom)*tmp

    overlap_matrix(2,2) = frac_both * (1.0_jprb - pair_cloud_cover)
    overlap_matrix(2,3) = frac_both * (pair_cloud_cover - cf_upper)
    overlap_matrix(3,2) = frac_both * (pair_cloud_cover - cf_lower)
    overlap_matrix(3,3) = frac_both * (cf_upper+cf_lower-pair_cloud_cover)

  end subroutine calc_alpha_overlap_matrix_fwd


  subroutine calc_alpha_overlap_matrix_ad(op, op_inhom, frac_upper_in, frac_lower_in, overlap_matrix, overlap_ad, &
       &  op_ad, op_inhom_ad, frac_upper_ad, frac_lower_ad)

    use parkind1, only : jprb
    implicit none
    real(jprb), intent(in) :: op, op_inhom
    real(jprb), intent(in),  dimension(NREGION) :: frac_upper_in, frac_lower_in
    real(jprb), intent(in) :: overlap_matrix(NREGION,NREGION)
    real(jprb), intent(inout) :: overlap_ad(NREGION,NREGION)

    real(jprb), intent(inout) :: op_ad, op_inhom_ad
    real(jprb), intent(inout), dimension(NREGION) :: frac_upper_ad, frac_lower_ad

    ! Recompute forward intermediates and do reverse sweep.
    ! This is identical to the standalone version previously supplied, but kept here
    ! as an internal helper for calc_overlap_matrices_ad.

    real(jprb), parameter :: CF_MIN = 1.0e-6_jprb

    real(jprb) :: cf_upper0, cf_lower0, max_cf0, tmp0, pair0
    real(jprb) :: frac_both
    real(jprb) :: denom, one_over_cf
    real(jprb) :: cf_u_eff, cf_l_eff, cf_upper, cf_lower
    real(jprb) :: max_cf, tmp, pair

    ! adjoints
    real(jprb) :: cf_upper0_ad, cf_lower0_ad, max_cf0_ad, tmp0_ad, pair0_ad
    real(jprb) :: frac_both_ad
    real(jprb) :: denom_ad, one_over_cf_ad
    real(jprb) :: cf_u_eff_ad, cf_l_eff_ad, cf_upper_ad, cf_lower_ad
    real(jprb) :: max_cf_ad, tmp_ad, pair_ad

    integer :: j

    ! ---------- Forward recompute (block 1) ----------
    cf_upper0 = sum(frac_upper_in(2:NREGION))
    cf_lower0 = sum(frac_lower_in(2:NREGION))

    max_cf0 = max(cf_upper0, cf_lower0)
    tmp0 = (cf_upper0+cf_lower0-cf_upper0*cf_lower0)
    pair0 = op*max_cf0 + (1.0_jprb-op)*tmp0

    frac_both = cf_upper0 + cf_lower0 - pair0

    cf_u_eff = max(cf_upper0, CF_MIN)
    cf_l_eff = max(cf_lower0, CF_MIN)
    cf_upper = frac_upper_in(3) / cf_u_eff
    cf_lower = frac_lower_in(3) / cf_l_eff

    max_cf = max(cf_upper, cf_lower)
    tmp = (cf_upper+cf_lower-cf_upper*cf_lower)
    pair = op_inhom*max_cf + (1.0_jprb-op_inhom)*tmp

    ! ---------- Init adjoints ----------
    cf_upper0_ad = 0.0_jprb; cf_lower0_ad = 0.0_jprb
    max_cf0_ad = 0.0_jprb; tmp0_ad = 0.0_jprb; pair0_ad = 0.0_jprb
    frac_both_ad = 0.0_jprb
    denom_ad = 0.0_jprb; one_over_cf_ad = 0.0_jprb
    cf_u_eff_ad = 0.0_jprb; cf_l_eff_ad = 0.0_jprb
    cf_upper_ad = 0.0_jprb; cf_lower_ad = 0.0_jprb
    max_cf_ad = 0.0_jprb; tmp_ad = 0.0_jprb; pair_ad = 0.0_jprb

    ! ---------- Reverse: inhomogeneous cloudy block (2/3) ----------
    ! overlap(3,3) = frac_both*(cf_upper+cf_lower-pair)
    frac_both_ad = frac_both_ad + overlap_ad(3,3) * (cf_upper+cf_lower-pair)
    cf_upper_ad  = cf_upper_ad  + overlap_ad(3,3) * frac_both
    cf_lower_ad  = cf_lower_ad  + overlap_ad(3,3) * frac_both
    pair_ad      = pair_ad      - overlap_ad(3,3) * frac_both
    overlap_ad(3,3) = 0.0_jprb

    ! overlap(3,2) = frac_both*(pair-cf_lower)
    frac_both_ad = frac_both_ad + overlap_ad(3,2) * (pair-cf_lower)
    pair_ad      = pair_ad      + overlap_ad(3,2) * frac_both
    cf_lower_ad  = cf_lower_ad  - overlap_ad(3,2) * frac_both
    overlap_ad(3,2) = 0.0_jprb

    ! overlap(2,3) = frac_both*(pair-cf_upper)
    frac_both_ad = frac_both_ad + overlap_ad(2,3) * (pair-cf_upper)
    pair_ad      = pair_ad      + overlap_ad(2,3) * frac_both
    cf_upper_ad  = cf_upper_ad  - overlap_ad(2,3) * frac_both
    overlap_ad(2,3) = 0.0_jprb

    ! overlap(2,2) = frac_both*(1-pair)
    frac_both_ad = frac_both_ad + overlap_ad(2,2) * (1.0_jprb-pair)
    pair_ad      = pair_ad      - overlap_ad(2,2) * frac_both
    overlap_ad(2,2) = 0.0_jprb

    ! pair = op_inhom*max_cf + (1-op_inhom)*tmp
    op_inhom_ad = op_inhom_ad + pair_ad * (max_cf - tmp)
    max_cf_ad   = max_cf_ad   + pair_ad * op_inhom
    tmp_ad      = tmp_ad      + pair_ad * (1.0_jprb - op_inhom)
    pair_ad = 0.0_jprb

    ! tmp = cf_upper + cf_lower - cf_upper*cf_lower
    cf_upper_ad = cf_upper_ad + tmp_ad * (1.0_jprb - cf_lower)
    cf_lower_ad = cf_lower_ad + tmp_ad * (1.0_jprb - cf_upper)
    tmp_ad = 0.0_jprb

    ! max_cf = max(cf_upper, cf_lower)
    if (cf_upper >= cf_lower) then
      cf_upper_ad = cf_upper_ad + max_cf_ad
    else
      cf_lower_ad = cf_lower_ad + max_cf_ad
    end if
    max_cf_ad = 0.0_jprb

    ! cf_upper = frac_upper(3)/cf_u_eff
    if (cf_u_eff > CF_MIN) then
      frac_upper_ad(3) = frac_upper_ad(3) + cf_upper_ad / cf_u_eff
      cf_u_eff_ad = cf_u_eff_ad - cf_upper_ad * frac_upper_in(3) / (cf_u_eff*cf_u_eff)
    end if
    cf_upper_ad = 0.0_jprb

    ! cf_lower = frac_lower(3)/cf_l_eff
    if (cf_l_eff > CF_MIN) then
      frac_lower_ad(3) = frac_lower_ad(3) + cf_lower_ad / cf_l_eff
      cf_l_eff_ad = cf_l_eff_ad - cf_lower_ad * frac_lower_in(3) / (cf_l_eff*cf_l_eff)
    end if
    cf_lower_ad = 0.0_jprb

    ! cf_u_eff = max(cf_upper0, CF_MIN)
    if (cf_upper0 > CF_MIN) then
      cf_upper0_ad = cf_upper0_ad + cf_u_eff_ad
    end if
    cf_u_eff_ad = 0.0_jprb

    ! cf_l_eff = max(cf_lower0, CF_MIN)
    if (cf_lower0 > CF_MIN) then
      cf_lower0_ad = cf_lower0_ad + cf_l_eff_ad
    end if
    cf_l_eff_ad = 0.0_jprb

    ! ---------- Reverse: clear/cloud and clear/clear block ----------
    ! overlap(3,1) and overlap(2,1): use denom=max(cf_upper0,CF_MIN), tmp=(pair0-cf_lower0)
    denom = max(cf_upper0, CF_MIN)
    one_over_cf = 1.0_jprb/denom
    tmp0 = (pair0 - cf_lower0)

    ! overlap(3,1) = tmp0 * frac_upper(3) * one_over_cf
    pair0_ad = pair0_ad + overlap_ad(3,1) * frac_upper_in(3) * one_over_cf
    cf_lower0_ad = cf_lower0_ad - overlap_ad(3,1) * frac_upper_in(3) * one_over_cf
    frac_upper_ad(3) = frac_upper_ad(3) + overlap_ad(3,1) * tmp0 * one_over_cf
    one_over_cf_ad = one_over_cf_ad + overlap_ad(3,1) * tmp0 * frac_upper_in(3)
    overlap_ad(3,1) = 0.0_jprb

    ! overlap(2,1) = tmp0 * frac_upper(2) * one_over_cf
    pair0_ad = pair0_ad + overlap_ad(2,1) * frac_upper_in(2) * one_over_cf
    cf_lower0_ad = cf_lower0_ad - overlap_ad(2,1) * frac_upper_in(2) * one_over_cf
    frac_upper_ad(2) = frac_upper_ad(2) + overlap_ad(2,1) * tmp0 * one_over_cf
    one_over_cf_ad = one_over_cf_ad + overlap_ad(2,1) * tmp0 * frac_upper_in(2)
    overlap_ad(2,1) = 0.0_jprb

    ! one_over_cf = 1/denom
    denom_ad = denom_ad - one_over_cf_ad * one_over_cf*one_over_cf
    one_over_cf_ad = 0.0_jprb
    if (cf_upper0 > CF_MIN) cf_upper0_ad = cf_upper0_ad + denom_ad
    denom_ad = 0.0_jprb

    ! overlap(1,3) and overlap(1,2): denom=max(cf_lower0,CF_MIN), tmp=(pair0-cf_upper0)
    denom = max(cf_lower0, CF_MIN)
    one_over_cf = 1.0_jprb/denom
    tmp0 = (pair0 - cf_upper0)

    ! overlap(1,3) = tmp0 * frac_lower(3) / denom
    pair0_ad = pair0_ad + overlap_ad(1,3) * frac_lower_in(3) * one_over_cf
    cf_upper0_ad = cf_upper0_ad - overlap_ad(1,3) * frac_lower_in(3) * one_over_cf
    frac_lower_ad(3) = frac_lower_ad(3) + overlap_ad(1,3) * tmp0 * one_over_cf
    one_over_cf_ad = one_over_cf_ad + overlap_ad(1,3) * tmp0 * frac_lower_in(3)
    overlap_ad(1,3) = 0.0_jprb

    ! overlap(1,2) = tmp0 * frac_lower(2) / denom
    pair0_ad = pair0_ad + overlap_ad(1,2) * frac_lower_in(2) * one_over_cf
    cf_upper0_ad = cf_upper0_ad - overlap_ad(1,2) * frac_lower_in(2) * one_over_cf
    frac_lower_ad(2) = frac_lower_ad(2) + overlap_ad(1,2) * tmp0 * one_over_cf
    one_over_cf_ad = one_over_cf_ad + overlap_ad(1,2) * tmp0 * frac_lower_in(2)
    overlap_ad(1,2) = 0.0_jprb

    denom_ad = denom_ad - one_over_cf_ad * one_over_cf*one_over_cf
    one_over_cf_ad = 0.0_jprb
    if (cf_lower0 > CF_MIN) cf_lower0_ad = cf_lower0_ad + denom_ad
    denom_ad = 0.0_jprb

    ! overlap(1,1) = 1 - pair0
    pair0_ad = pair0_ad - overlap_ad(1,1)
    overlap_ad(1,1) = 0.0_jprb

    ! frac_both = cf_upper0 + cf_lower0 - pair0
    cf_upper0_ad = cf_upper0_ad + frac_both_ad
    cf_lower0_ad = cf_lower0_ad + frac_both_ad
    pair0_ad     = pair0_ad     - frac_both_ad
    frac_both_ad = 0.0_jprb

    ! pair0 = op*max_cf0 + (1-op)*tmp0
    max_cf0 = max(cf_upper0, cf_lower0)
    tmp0 = (cf_upper0+cf_lower0-cf_upper0*cf_lower0)

    op_ad    = op_ad + pair0_ad * (max_cf0 - tmp0)
    max_cf0_ad = max_cf0_ad + pair0_ad * op
    tmp0_ad  = tmp0_ad + pair0_ad * (1.0_jprb - op)
    pair0_ad = 0.0_jprb

    ! tmp0 = cf_upper0 + cf_lower0 - cf_upper0*cf_lower0
    cf_upper0_ad = cf_upper0_ad + tmp0_ad * (1.0_jprb - cf_lower0)
    cf_lower0_ad = cf_lower0_ad + tmp0_ad * (1.0_jprb - cf_upper0)
    tmp0_ad = 0.0_jprb

    ! max_cf0 = max(cf_upper0, cf_lower0)
    if (cf_upper0 >= cf_lower0) then
      cf_upper0_ad = cf_upper0_ad + max_cf0_ad
    else
      cf_lower0_ad = cf_lower0_ad + max_cf0_ad
    end if
    max_cf0_ad = 0.0_jprb

    ! cf_upper0 = sum(frac_upper(2:))
    do j = 2,NREGION
      frac_upper_ad(j) = frac_upper_ad(j) + cf_upper0_ad
    end do
    cf_upper0_ad = 0.0_jprb

    ! cf_lower0 = sum(frac_lower(2:))
    do j = 2,NREGION
      frac_lower_ad(j) = frac_lower_ad(j) + cf_lower0_ad
    end do
    cf_lower0_ad = 0.0_jprb

  end subroutine calc_alpha_overlap_matrix_ad

end subroutine calc_overlap_matrices_ad


end module tcrad_ad
