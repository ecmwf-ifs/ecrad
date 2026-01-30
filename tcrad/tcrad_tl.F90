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

module tcrad_tl

  use parkind1, only : jpim, jprb
  use tcrad,    only : NREGION, &
       & ITwoStreamElsasser, ITwoStreamEddington, ITwoStreamLegendre, &
       & ITwoStreamHybrid, ITwoStreamScaledWiscombeGrams, &
       & i_two_stream_scheme, lw_diffusivity, lw_diffusivity_cloud, &
       & MIN_K_SQUARED, OD_THRESH, OD_THRESH_2STREAM

  implicit none
  public

contains

! ===== FILE: tcrad_calc_overlap_matrices_tl.F90 =====
function calc_alpha_overlap_matrix_tl(op, op_tl, op_inhom, op_inhom_tl, &
     &  frac_upper, frac_upper_tl, frac_lower, frac_lower_tl) result(overlap_matrix)

  use parkind1, only : jprb

  ! Overlap parameter for cloud boundaries and for internal inhomogeneities
  real(jprb), intent(in) :: op, op_tl, op_inhom, op_inhom_tl

  ! Fraction of the gridbox occupied by each region in the upper and lower layers
  real(jprb), intent(in), dimension(NREGION) :: frac_upper, frac_lower
  real(jprb), intent(in), dimension(NREGION) :: frac_upper_tl, frac_lower_tl

  ! Output overlap matrix: we return the NONLINEAR matrix; its TL is returned
  ! in the module variable overlap_matrix_tl via an internal saved array below.
  !
  ! NOTE: For ease of inclusion within existing code, this function returns the
  ! nonlinear overlap_matrix and we provide the TL in a companion argument-like
  ! saved variable overlap_matrix_tl. In the calling TL routine we immediately
  ! copy it out.
  real(jprb) :: overlap_matrix(NREGION,NREGION)

  ! Locals (nonlinear)
  real(jprb) :: pair_cloud_cover
  real(jprb) :: cf_upper, cf_lower
  real(jprb) :: one_over_cf
  real(jprb) :: frac_both
  real(jprb) :: cf_eff, cf_eff2, max_cf
  real(jprb) :: tmp, denom

  ! Locals (TL)
  real(jprb) :: pair_cloud_cover_tl
  real(jprb) :: cf_upper_tl, cf_lower_tl
  real(jprb) :: one_over_cf_tl
  real(jprb) :: frac_both_tl
  real(jprb) :: max_cf_tl
  real(jprb) :: tmp_tl, denom_tl
  real(jprb) :: overlap_matrix_tl(NREGION,NREGION)

  ! To avoid zero-division in TL, mirror the nonlinear safeguards
  real(jprb), parameter :: CF_MIN = 1.0e-6_jprb

  integer :: j

  overlap_matrix(:,:) = 0.0_jprb
  overlap_matrix_tl(:,:) = 0.0_jprb

  ! Cloud fraction of upper and lower layers
  cf_upper = sum(frac_upper(2:NREGION))
  cf_lower = sum(frac_lower(2:NREGION))
  cf_upper_tl = sum(frac_upper_tl(2:NREGION))
  cf_lower_tl = sum(frac_lower_tl(2:NREGION))

  max_cf = max(cf_upper, cf_lower)
  max_cf_tl = merge(cf_upper_tl, cf_lower_tl, cf_upper >= cf_lower)

  tmp = (cf_upper+cf_lower-cf_upper*cf_lower)
  tmp_tl = (cf_upper_tl+cf_lower_tl) - (cf_upper_tl*cf_lower + cf_upper*cf_lower_tl)

  pair_cloud_cover = op*max_cf + (1.0_jprb - op)*tmp
  pair_cloud_cover_tl = op_tl*max_cf + op*max_cf_tl + (-op_tl)*tmp + (1.0_jprb - op)*tmp_tl

  ! Clear in both layers
  overlap_matrix(1,1) = 1.0_jprb - pair_cloud_cover
  overlap_matrix_tl(1,1) = -pair_cloud_cover_tl

  ! Clear in upper layer, cloudy in lower layer
  cf_eff = max(cf_lower, CF_MIN)
  denom = cf_eff
  denom_tl = merge(cf_lower_tl, 0.0_jprb, cf_lower > CF_MIN)
  one_over_cf = 1.0_jprb / denom
  one_over_cf_tl = -one_over_cf*one_over_cf * denom_tl

  tmp = (pair_cloud_cover - cf_upper)
  tmp_tl = (pair_cloud_cover_tl - cf_upper_tl)

  overlap_matrix(1,2) = tmp * frac_lower(2) * one_over_cf
  overlap_matrix_tl(1,2) = tmp_tl * frac_lower(2) * one_over_cf &
       &                 + tmp * frac_lower_tl(2) * one_over_cf &
       &                 + tmp * frac_lower(2) * one_over_cf_tl
  overlap_matrix(1,3) = tmp * frac_lower(3) * one_over_cf
  overlap_matrix_tl(1,3) = tmp_tl * frac_lower(3) * one_over_cf &
       &                 + tmp * frac_lower_tl(3) * one_over_cf &
       &                 + tmp * frac_lower(3) * one_over_cf_tl

  ! Clear in lower layer, cloudy in upper layer
  cf_eff = max(cf_upper, CF_MIN)
  denom = cf_eff
  denom_tl = merge(cf_upper_tl, 0.0_jprb, cf_upper > CF_MIN)
  one_over_cf = 1.0_jprb / denom
  one_over_cf_tl = -one_over_cf*one_over_cf * denom_tl

  tmp = (pair_cloud_cover - cf_lower)
  tmp_tl = (pair_cloud_cover_tl - cf_lower_tl)

  overlap_matrix(2,1) = tmp * frac_upper(2) * one_over_cf
  overlap_matrix_tl(2,1) = tmp_tl * frac_upper(2) * one_over_cf &
       &                 + tmp * frac_upper_tl(2) * one_over_cf &
       &                 + tmp * frac_upper(2) * one_over_cf_tl
  overlap_matrix(3,1) = tmp * frac_upper(3) * one_over_cf
  overlap_matrix_tl(3,1) = tmp_tl * frac_upper(3) * one_over_cf &
       &                 + tmp * frac_upper_tl(3) * one_over_cf &
       &                 + tmp * frac_upper(3) * one_over_cf_tl

  ! Cloudy in both layers
  frac_both = cf_upper + cf_lower - pair_cloud_cover
  frac_both_tl = cf_upper_tl + cf_lower_tl - pair_cloud_cover_tl

  ! Redefine cloud fractions treating region 3 as the "cloud" within the cloudy part
  cf_eff  = max(cf_upper, CF_MIN)
  denom   = cf_eff
  denom_tl = merge(cf_upper_tl, 0.0_jprb, cf_upper > CF_MIN)
  cf_upper = frac_upper(3) / denom
  cf_upper_tl = (frac_upper_tl(3)*denom - frac_upper(3)*denom_tl) / (denom*denom)

  cf_eff2 = max(cf_lower, CF_MIN)
  denom   = cf_eff2
  denom_tl = merge(cf_lower_tl, 0.0_jprb, cf_lower > CF_MIN)
  cf_lower = frac_lower(3) / denom
  cf_lower_tl = (frac_lower_tl(3)*denom - frac_lower(3)*denom_tl) / (denom*denom)

  max_cf = max(cf_upper, cf_lower)
  max_cf_tl = merge(cf_upper_tl, cf_lower_tl, cf_upper >= cf_lower)

  tmp = (cf_upper+cf_lower-cf_upper*cf_lower)
  tmp_tl = (cf_upper_tl+cf_lower_tl) - (cf_upper_tl*cf_lower + cf_upper*cf_lower_tl)

  pair_cloud_cover = op_inhom*max_cf + (1.0_jprb - op_inhom)*tmp
  pair_cloud_cover_tl = op_inhom_tl*max_cf + op_inhom*max_cf_tl &
       &              + (-op_inhom_tl)*tmp + (1.0_jprb - op_inhom)*tmp_tl

  overlap_matrix(2,2) = frac_both * (1.0_jprb - pair_cloud_cover)
  overlap_matrix_tl(2,2) = frac_both_tl * (1.0_jprb - pair_cloud_cover) &
       &                 + frac_both * (-pair_cloud_cover_tl)

  overlap_matrix(2,3) = frac_both * (pair_cloud_cover - cf_upper)
  overlap_matrix_tl(2,3) = frac_both_tl * (pair_cloud_cover - cf_upper) &
       &                 + frac_both * (pair_cloud_cover_tl - cf_upper_tl)

  overlap_matrix(3,2) = frac_both * (pair_cloud_cover - cf_lower)
  overlap_matrix_tl(3,2) = frac_both_tl * (pair_cloud_cover - cf_lower) &
       &                 + frac_both * (pair_cloud_cover_tl - cf_lower_tl)

  overlap_matrix(3,3) = frac_both * (cf_upper+cf_lower-pair_cloud_cover)
  overlap_matrix_tl(3,3) = frac_both_tl * (cf_upper+cf_lower-pair_cloud_cover) &
       &                 + frac_both * ((cf_upper_tl+cf_lower_tl) - pair_cloud_cover_tl)

  ! Store TL in a module variable via host association workaround:
  ! In practice, the calling TL routine will include this function in the same
  ! file and will call the internal helper below to retrieve overlap_matrix_tl.
  call save_overlap_matrix_tl(overlap_matrix_tl)

contains
  subroutine save_overlap_matrix_tl(mat_tl)
    real(jprb), intent(in) :: mat_tl(NREGION,NREGION)
    ! Use a named common block to "return" TL without changing function signature
    common /tcrad_alpha_overlap_tl_common/ overlap_matrix_tl_common
    real(jprb) :: overlap_matrix_tl_common(NREGION,NREGION)
    overlap_matrix_tl_common(:,:) = mat_tl(:,:)
  end subroutine save_overlap_matrix_tl

end function calc_alpha_overlap_matrix_tl


subroutine calc_overlap_matrices_tl(nlev, &
     &     region_fracs, region_fracs_tl, overlap_param, overlap_param_tl, &
     &     u_overlap, u_overlap_tl, v_overlap, v_overlap_tl, &
     &     decorrelation_scaling, &
     &     cloud_fraction_threshold, cloud_cover, cloud_cover_tl)

  use parkind1,     only : jprb
  use yomhook, only : lhook, dr_hook, jphook

  integer,  intent(in) :: nlev

  ! Atmospheric/macrophysical inputs (TL w.r.t. these)
  real(jprb), intent(in), dimension(1:NREGION,nlev)  :: region_fracs
  real(jprb), intent(in), dimension(1:NREGION,nlev)  :: region_fracs_tl
  real(jprb), intent(in), dimension(:)  :: overlap_param
  real(jprb), intent(in), dimension(:)  :: overlap_param_tl

  ! Outputs (nonlinear + TL)
  real(jprb), intent(out), dimension(NREGION,NREGION,nlev+1) :: u_overlap, v_overlap
  real(jprb), intent(out), dimension(NREGION,NREGION,nlev+1) :: u_overlap_tl, v_overlap_tl

  ! Optional config (NOT TL'd)
  real(jprb), intent(in), optional :: decorrelation_scaling
  real(jprb), intent(in), optional :: cloud_fraction_threshold

  ! Optional outputs
  real(jprb), intent(out), optional :: cloud_cover
  real(jprb), intent(out), optional :: cloud_cover_tl

  integer  :: jlev, jupper, jlower

  real(jprb) :: overlap_matrix(NREGION,NREGION)
  real(jprb) :: overlap_matrix_tl(NREGION,NREGION)

  real(jprb) :: frac_upper(NREGION), frac_lower(NREGION)
  real(jprb) :: frac_upper_tl(NREGION), frac_lower_tl(NREGION)

  real(jprb) :: op(NREGION), op_tl(NREGION)

  real(jprb) :: frac_threshold
  real(jprb) :: used_decorrelation_scaling
  real(jprb) :: expo, expo_tl, pexp

  real(jprb) :: denom, denom_tl

  real(jprb) :: prod_v11, prod_v11_tl

  real(jphook) :: hook_handle

  ! Common block used to receive TL from calc_alpha_overlap_matrix_tl
  common /tcrad_alpha_overlap_tl_common/ overlap_matrix_tl_common
  real(jprb) :: overlap_matrix_tl_common(NREGION,NREGION)

  if (lhook) call dr_hook('tcrad:calc_overlap_matrices_tl',0,hook_handle)

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

  u_overlap(:,:,:) = 0.0_jprb
  v_overlap(:,:,:) = 0.0_jprb
  u_overlap_tl(:,:,:) = 0.0_jprb
  v_overlap_tl(:,:,:) = 0.0_jprb

  ! Top boundary: single clear-sky region
  frac_upper(1) = 1.0_jprb
  frac_upper(2:NREGION) = 0.0_jprb
  frac_upper_tl(:) = 0.0_jprb

  op(:) = 1.0_jprb
  op_tl(:) = 0.0_jprb

  do jlev = 1,nlev+1

    if (jlev > nlev) then
      frac_lower(1) = 1.0_jprb
      frac_lower(2:NREGION) = 0.0_jprb
      frac_lower_tl(:) = 0.0_jprb
    else
      frac_lower = region_fracs(1:NREGION,jlev)
      frac_lower_tl = region_fracs_tl(1:NREGION,jlev)
    end if

    if (jlev == 1 .or. jlev > nlev) then
      op(:) = 1.0_jprb
      op_tl(:) = 0.0_jprb
    else
      op(1) = overlap_param(jlev-1)
      op_tl(1) = overlap_param_tl(jlev-1)

      if (op(1) >= 0.0_jprb) then
        pexp = 1.0_jprb/used_decorrelation_scaling
        op(2:NREGION) = op(1)**pexp
        ! TL: d(op^p)/d(op)=p*op^(p-1)
        op_tl(2:NREGION) = pexp * op(1)**(pexp-1.0_jprb) * op_tl(1)
      else
        op(2:NREGION) = op(1)
        op_tl(2:NREGION) = op_tl(1)
      end if
    end if

    overlap_matrix = calc_alpha_overlap_matrix_tl(op(1), op_tl(1), op(2), op_tl(2), &
         &  frac_upper, frac_upper_tl, frac_lower, frac_lower_tl)

    overlap_matrix_tl(:,:) = overlap_matrix_tl_common(:,:)

    ! Convert to directional overlap matrices
    do jupper = 1,NREGION
      do jlower = 1,NREGION
        if (frac_lower(jlower) >= frac_threshold) then
          u_overlap(jupper,jlower,jlev) = overlap_matrix(jupper,jlower) / frac_lower(jlower)
          denom = frac_lower(jlower)
          denom_tl = frac_lower_tl(jlower)
          u_overlap_tl(jupper,jlower,jlev) = (overlap_matrix_tl(jupper,jlower)*denom &
               &                           - overlap_matrix(jupper,jlower)*denom_tl) / (denom*denom)
        else
          u_overlap(jupper,jlower,jlev) = 0.0_jprb
          u_overlap_tl(jupper,jlower,jlev) = 0.0_jprb
        end if
        if (frac_upper(jupper) >= frac_threshold) then
          v_overlap(jlower,jupper,jlev) = overlap_matrix(jupper,jlower) / frac_upper(jupper)
          denom = frac_upper(jupper)
          denom_tl = frac_upper_tl(jupper)
          v_overlap_tl(jlower,jupper,jlev) = (overlap_matrix_tl(jupper,jlower)*denom &
               &                           - overlap_matrix(jupper,jlower)*denom_tl) / (denom*denom)
        else
          v_overlap(jlower,jupper,jlev) = 0.0_jprb
          v_overlap_tl(jlower,jupper,jlev) = 0.0_jprb
        end if
      end do
    end do

    frac_upper = frac_lower
    frac_upper_tl = frac_lower_tl

  end do

  if (present(cloud_cover)) then
    prod_v11 = 1.0_jprb
    prod_v11_tl = 0.0_jprb
    do jlev = 1,nlev+1
      prod_v11_tl = prod_v11_tl * v_overlap(1,1,jlev) + prod_v11 * v_overlap_tl(1,1,jlev)
      prod_v11    = prod_v11    * v_overlap(1,1,jlev)
    end do
    cloud_cover = 1.0_jprb - prod_v11
    if (present(cloud_cover_tl)) then
      cloud_cover_tl = -prod_v11_tl
    end if
  else
    if (present(cloud_cover_tl)) then
      cloud_cover_tl = 0.0_jprb
    end if
  end if

  if (lhook) call dr_hook('tcrad:calc_overlap_matrices_tl',1,hook_handle)

end subroutine calc_overlap_matrices_tl


! ===== FILE: tcrad_calc_radiance_dn_tl.F90 =====
subroutine calc_radiance_dn_tl(nspec, nlev, &
     &  weight, &
     &  transmittance, transmittance_tl, &
     &  source_dn, source_dn_tl, &
     &  v_overlap, v_overlap_tl, &
     &  radiance_dn, radiance_dn_tl)

  use parkind1, only           : jpim, jprb
  use yomhook,  only           : lhook, dr_hook, jphook

  implicit none

  ! Inputs
  integer(jpim), intent(in) :: nspec, nlev
  real(jprb), intent(in) :: weight

  ! Atmospheric/macrophysical inputs (TL w.r.t. these)
  real(jprb), intent(in),  dimension(nspec,NREGION,nlev) :: transmittance
  real(jprb), intent(in),  dimension(nspec,NREGION,nlev) :: transmittance_tl

  real(jprb), intent(in),  dimension(nspec,NREGION,nlev) :: source_dn
  real(jprb), intent(in),  dimension(nspec,NREGION,nlev) :: source_dn_tl

  ! Overlap matrices (treated as atmospheric/macrophysical)
  real(jprb), intent(in),  dimension(NREGION,NREGION,nlev+1) :: v_overlap
  real(jprb), intent(in),  dimension(NREGION,NREGION,nlev+1) :: v_overlap_tl

  ! Outputs (accumulate into existing profiles)
  real(jprb), intent(inout), dimension(nspec,nlev+1) :: radiance_dn
  real(jprb), intent(inout), dimension(nspec,nlev+1) :: radiance_dn_tl

  ! Local variables
  real(jprb), dimension(nspec,NREGION) :: radiance_base, radiance_top
  real(jprb), dimension(nspec,NREGION) :: radiance_base_tl, radiance_top_tl

  integer(jpim) :: jlev, jspec, jreg, jreg2

  real(jphook) :: hook_handle

  if (lhook) call dr_hook('tcrad:calc_radiance_dn_tl',0,hook_handle)

  ! Start with zero at TOA
  radiance_top    = 0.0_jprb
  radiance_top_tl = 0.0_jprb

  do jlev = 1,nlev
    ! Schwarzschild equation in each region
    radiance_base = transmittance(:,:,jlev)*radiance_top + weight * source_dn(:,:,jlev)
    radiance_base_tl = transmittance_tl(:,:,jlev)*radiance_top &
         &           + transmittance(:,:,jlev)*radiance_top_tl &
         &           + weight * source_dn_tl(:,:,jlev)

    ! Apply overlap rules to obtain radiances at top of layer below:
    ! radiance_top = v_overlap(:,:,jlev+1) * radiance_base
    radiance_top(:,:)    = 0.0_jprb
    radiance_top_tl(:,:) = 0.0_jprb
    do jreg = 1,NREGION
      do jreg2 = 1,NREGION
        radiance_top(:,jreg) = radiance_top(:,jreg) &
             &  + v_overlap(jreg,jreg2,jlev+1) * radiance_base(:,jreg2)
        radiance_top_tl(:,jreg) = radiance_top_tl(:,jreg) &
             &  + v_overlap_tl(jreg,jreg2,jlev+1) * radiance_base(:,jreg2) &
             &  + v_overlap(jreg,jreg2,jlev+1) * radiance_base_tl(:,jreg2)
      end do
    end do

    ! Save radiances averaged over regions (add to existing)
    radiance_dn(:,jlev+1)    = radiance_dn(:,jlev+1)    + sum(radiance_top,2)
    radiance_dn_tl(:,jlev+1) = radiance_dn_tl(:,jlev+1) + sum(radiance_top_tl,2)
  end do

  if (lhook) call dr_hook('tcrad:calc_radiance_dn_tl',1,hook_handle)

end subroutine calc_radiance_dn_tl


! ===== FILE: tcrad_calc_radiance_tl.F90 =====
subroutine calc_radiance_tl(nspec, nlev, &
     &  surf_emission, surf_emission_tl, &
     &  surf_albedo,  surf_albedo_tl, &
     &  planck_hl,    planck_hl_tl, &
     &  cloud_fraction, cloud_fraction_tl, &
     &  fractional_std, fractional_std_tl, &
     &  od_clear, od_clear_tl, &
     &  od_cloud, od_cloud_tl, &
     &  ssa_cloud, asymmetry_cloud, &
     &  overlap_param, overlap_param_tl, &
     &  mu, &
     &  radiance, radiance_tl, &
     &  cloud_cover, cloud_cover_tl, &
     &  do_specular_surface)

  use parkind1, only           : jpim, jprb
  use yomhook,  only           : lhook, dr_hook, jphook

  implicit none

  real(jprb), parameter :: PI          = acos(-1.0_jprb)
  real(jprb), parameter :: ONE_OVER_PI = 1.0_jprb / PI

  ! Inputs
  integer(jpim), intent(in) :: nspec, nlev

  ! Atmospheric/macrophysical inputs (TL w.r.t. these)
  real(jprb), intent(in), dimension(nspec) :: surf_emission, surf_albedo
  real(jprb), intent(in), dimension(nspec) :: surf_emission_tl, surf_albedo_tl

  real(jprb), intent(in), dimension(nspec,nlev+1) :: planck_hl
  real(jprb), intent(in), dimension(nspec,nlev+1) :: planck_hl_tl

  real(jprb), intent(in), dimension(nlev) :: cloud_fraction, fractional_std
  real(jprb), intent(in), dimension(nlev) :: cloud_fraction_tl, fractional_std_tl

  real(jprb), intent(in), dimension(nspec,nlev) :: od_clear, od_cloud
  real(jprb), intent(in), dimension(nspec,nlev) :: od_clear_tl, od_cloud_tl

  ! Microphysical inputs (NO TL w.r.t. these)
  real(jprb), intent(in), dimension(nspec,nlev) :: ssa_cloud, asymmetry_cloud

  ! Overlap parameter (TL w.r.t. this)
  real(jprb), intent(in), dimension(nlev-1) :: overlap_param
  real(jprb), intent(in), dimension(nlev-1) :: overlap_param_tl

  ! Geometry (not TL'd here)
  real(jprb), intent(in) :: mu

  ! Outputs
  real(jprb), intent(out), dimension(nspec) :: radiance
  real(jprb), intent(out), dimension(nspec) :: radiance_tl

  real(jprb), intent(out), optional :: cloud_cover
  real(jprb), intent(out), optional :: cloud_cover_tl

  logical, intent(in), optional :: do_specular_surface

  ! Local variables (nonlinear)
  real(jprb), dimension(nspec,NREGION,nlev)   :: od
  real(jprb), dimension(nspec,NREGION,nlev)   :: od_tl

  real(jprb), dimension(nspec,2:NREGION,nlev) :: ssa

  real(jprb), dimension(nspec,NREGION,nlev) :: reflectance, transmittance
  real(jprb), dimension(nspec,NREGION,nlev) :: source_up, source_dn

  real(jprb), dimension(nspec,NREGION,nlev) :: reflectance_tl, transmittance_tl
  real(jprb), dimension(nspec,NREGION,nlev) :: source_up_tl, source_dn_tl

  logical :: is_cloud_free_layer(0:nlev+1)

  real(jprb), dimension(NREGION,NREGION,nlev+1) :: u_overlap, v_overlap
  real(jprb), dimension(NREGION,NREGION,nlev+1) :: u_overlap_tl, v_overlap_tl

  real(jprb), dimension(nspec,NREGION,nlev) :: flux_up_base, flux_dn_base
  real(jprb), dimension(nspec,NREGION,nlev) :: flux_up_top,  flux_dn_top

  real(jprb), dimension(nspec,NREGION,nlev) :: flux_up_base_tl, flux_dn_base_tl
  real(jprb), dimension(nspec,NREGION,nlev) :: flux_up_top_tl,  flux_dn_top_tl

  real(jprb), dimension(nspec, nlev+1) :: radiance_profile
  real(jprb), dimension(nspec, nlev+1) :: radiance_profile_tl

  real(jprb), dimension(nspec,NREGION) :: flux_up_surface
  real(jprb), dimension(nspec,NREGION) :: flux_up_surface_tl

  real(jprb) :: od_scaling(2:NREGION,nlev)
  real(jprb) :: od_scaling_tl(2:NREGION,nlev)

  real(jprb) :: region_fracs(1:NREGION,nlev)
  real(jprb) :: region_fracs_tl(1:NREGION,nlev)

  real(jprb), parameter :: cloud_fraction_threshold = 1.0e-6

  logical :: do_specular_surface_local

  integer(jpim) :: jreg, jlev, jspec

  real(jphook) :: hook_handle

  if (lhook) call dr_hook('tcrad:calc_radiance_tl',0,hook_handle)

  ! Default: Lambertian surface
  if (present(do_specular_surface)) then
    do_specular_surface_local = do_specular_surface
  else
    do_specular_surface_local = .false.
  endif

  ! --------------------------------------------------------
  ! Region properties (wavelength-independent)
  ! --------------------------------------------------------
  call calc_region_properties_tl(nlev, &
       &  cloud_fraction, cloud_fraction_tl, &
       &  .true., fractional_std, fractional_std_tl, &
       &  region_fracs, region_fracs_tl, &
       &  od_scaling, od_scaling_tl, cloud_fraction_threshold)

  ! --------------------------------------------------------
  ! Overlap matrices
  ! --------------------------------------------------------
  call calc_overlap_matrices_tl(nlev, &
       &  region_fracs, region_fracs_tl, &
       &  overlap_param, overlap_param_tl, &
       &  u_overlap, u_overlap_tl, v_overlap, v_overlap_tl, &
       &  0.5_jprb, &
       &  cloud_fraction_threshold, &
       &  cloud_cover, cloud_cover_tl)

  ! --------------------------------------------------------
  ! Combine gas/aerosol/cloud optical depth into regions
  ! --------------------------------------------------------
  od(:,1,:)    = od_clear
  od_tl(:,1,:) = od_clear_tl

  do jreg = 2,NREGION
    od(:,jreg,:) = od_clear + od_cloud*spread(od_scaling(jreg,:),1,nspec)
    od_tl(:,jreg,:) = od_clear_tl &
         &          + od_cloud_tl*spread(od_scaling(jreg,:),1,nspec) &
         &          + od_cloud*spread(od_scaling_tl(jreg,:),1,nspec)

    ! ssa is treated as microphysical (NO TL): compute nonlinear only
    ssa(:,jreg,:) = ssa_cloud(:,:)*od_cloud(:,:) &
         &  * spread(od_scaling(jreg,:),1,nspec) / od(:,jreg,:)
  end do

  ! Cloud-free layers
  is_cloud_free_layer(0) = .true.
  is_cloud_free_layer(1:nlev) = (region_fracs(1,:) == 1.0_jprb)
  is_cloud_free_layer(nlev+1) = .true.

  ! --------------------------------------------------------
  ! Layer-wise diffuse properties (two-stream)
  ! --------------------------------------------------------
  call calc_reflectance_transmittance_tl(nspec, nlev, NREGION, &
       &  region_fracs, region_fracs_tl, &
       &  planck_hl, planck_hl_tl, &
       &  od, od_tl, ssa, asymmetry_cloud, &
       &  reflectance, reflectance_tl, &
       &  transmittance, transmittance_tl, &
       &  source_up, source_up_tl, source_dn, source_dn_tl)

  ! --------------------------------------------------------
  ! Tripleclouds flux profile (adding method)
  ! --------------------------------------------------------
  call calc_two_stream_flux_tl(nspec, nlev, &
       &  surf_emission, surf_emission_tl, surf_albedo, surf_albedo_tl, &
       &  reflectance, reflectance_tl, transmittance, transmittance_tl, &
       &  source_up, source_up_tl, source_dn, source_dn_tl, &
       &  is_cloud_free_layer, u_overlap, u_overlap_tl, v_overlap, v_overlap_tl, &
       &  flux_up_base, flux_up_base_tl, flux_dn_base, flux_dn_base_tl, &
       &  flux_up_top, flux_up_top_tl, flux_dn_top, flux_dn_top_tl)

  ! --------------------------------------------------------
  ! Radiance profiles are additive, so initialize to zero
  ! --------------------------------------------------------
  radiance_profile    = 0.0_jprb
  radiance_profile_tl = 0.0_jprb

  if (do_specular_surface_local) then
    ! --------------------------------------------------------
    ! Specular surface branch: compute downward radiance then reflect
    ! --------------------------------------------------------

    call calc_radiance_trans_source_tl(nspec, nlev, NREGION, &
         &  mu, &
         &  region_fracs, region_fracs_tl, &
         &  planck_hl, planck_hl_tl, &
         &  od, od_tl, ssa, asymmetry_cloud, &
         &  flux_up_base, flux_up_base_tl, flux_dn_top, flux_dn_top_tl, &
         &  transmittance, transmittance_tl, &
         &  source_up, source_up_tl, source_dn, source_dn_tl)

    call calc_radiance_dn_tl(nspec, nlev, &
         &  ONE_OVER_PI, &
         &  transmittance, transmittance_tl, &
         &  source_dn, source_dn_tl, &
         &  v_overlap, v_overlap_tl, &
         &  radiance_profile, radiance_profile_tl)

    ! Reflect surface downward radiance upwards, spread into regions of lowest layer
    do jreg = 1,NREGION
      do jspec = 1,nspec
        flux_up_surface(jspec,jreg) = (surf_emission(jspec) + PI*surf_albedo(jspec)*radiance_profile(jspec,nlev+1)) &
             &                      * region_fracs(jreg,nlev)
        flux_up_surface_tl(jspec,jreg) = (surf_emission_tl(jspec) &
             &       + PI*(surf_albedo_tl(jspec)*radiance_profile(jspec,nlev+1) &
             &           + surf_albedo(jspec)*radiance_profile_tl(jspec,nlev+1))) * region_fracs(jreg,nlev) &
             &     + (surf_emission(jspec) + PI*surf_albedo(jspec)*radiance_profile(jspec,nlev+1)) &
             &       * region_fracs_tl(jreg,nlev)
      end do
    end do

    call calc_radiance_up_tl(nspec, nlev, &
         &  ONE_OVER_PI, &
         &  flux_up_surface, flux_up_surface_tl, &
         &  transmittance, transmittance_tl, &
         &  source_up, source_up_tl, &
         &  u_overlap, u_overlap_tl, &
         &  radiance_profile, radiance_profile_tl)

    radiance    = radiance_profile(:,1)
    radiance_tl = radiance_profile_tl(:,1)

  else
    ! --------------------------------------------------------
    ! Lambertian surface branch: upward radiance only
    ! --------------------------------------------------------

    call calc_radiance_trans_source_tl(nspec, nlev, NREGION, &
         &  mu, &
         &  region_fracs, region_fracs_tl, &
         &  planck_hl, planck_hl_tl, &
         &  od, od_tl, ssa, asymmetry_cloud, &
         &  flux_up_base, flux_up_base_tl, flux_dn_top, flux_dn_top_tl, &
         &  transmittance, transmittance_tl, &
         &  source_up, source_up_tl)

    call calc_radiance_up_tl(nspec, nlev, &
         &  ONE_OVER_PI, &
         &  flux_up_base(:,:,nlev), flux_up_base_tl(:,:,nlev), &
         &  transmittance, transmittance_tl, &
         &  source_up, source_up_tl, &
         &  u_overlap, u_overlap_tl, &
         &  radiance_profile, radiance_profile_tl)

    radiance    = radiance_profile(:,1)
    radiance_tl = radiance_profile_tl(:,1)

  end if

  if (lhook) call dr_hook('tcrad:calc_radiance_tl',1,hook_handle)

end subroutine calc_radiance_tl


! ===== FILE: tcrad_calc_radiance_trans_source_tl.F90 =====
subroutine calc_radiance_trans_source_tl(nspec, nlev, nreg, &
     &  mu, region_fracs, region_fracs_tl, planck_hl, planck_hl_tl, &
     &  od, od_tl, ssa, asymmetry, &
     &  flux_up_base, flux_up_base_tl, flux_dn_top, flux_dn_top_tl, &
     &  transmittance, transmittance_tl, &
     &  source_up, source_up_tl, source_dn, source_dn_tl)

  use yomhook,  only           : lhook, dr_hook, jphook
  
  ! Inputs
  integer(jpim), intent(in) :: nspec, nlev, nreg
  real(jprb), intent(in) :: mu

  ! Atmospheric / macrophysical inputs (TL w.r.t. these)
  real(jprb), intent(in),  dimension(nreg,nlev)   :: region_fracs
  real(jprb), intent(in),  dimension(nreg,nlev)   :: region_fracs_tl
  real(jprb), intent(in),  dimension(nspec,nlev+1):: planck_hl
  real(jprb), intent(in),  dimension(nspec,nlev+1):: planck_hl_tl
  real(jprb), intent(in),  dimension(nspec,nreg,nlev) :: od
  real(jprb), intent(in),  dimension(nspec,nreg,nlev) :: od_tl
  real(jprb), intent(in),  dimension(nspec,nreg,nlev) :: flux_up_base, flux_dn_top
  real(jprb), intent(in),  dimension(nspec,nreg,nlev) :: flux_up_base_tl, flux_dn_top_tl

  ! Microphysical inputs (NO TL w.r.t. these)
  real(jprb), intent(in),  dimension(nspec,2:nreg,nlev) :: ssa
  real(jprb), intent(in),  dimension(nspec,nlev) :: asymmetry

  ! Outputs (nonlinear + TL)
  real(jprb), intent(out), dimension(nspec,nreg,nlev) :: transmittance
  real(jprb), intent(out), dimension(nspec,nreg,nlev) :: transmittance_tl

  real(jprb), intent(out), dimension(nspec,nreg,nlev), optional :: source_up, source_dn
  real(jprb), intent(out), dimension(nspec,nreg,nlev), optional :: source_up_tl, source_dn_tl

  ! Local variables
  real(jprb), dimension(nspec,nreg) :: planck_top, planck_base
  real(jprb), dimension(nspec,nreg) :: planck_top_tl, planck_base_tl

  real(jprb) :: secant, factor, coeff, gamma1, gamma2, k_exponent, rt_factor
  real(jprb) :: exponential, one_minus_kmu
  real(jprb) :: p_same, p_opposite, planck_prime, c1, c2, scaling1, scaling2
  real(jprb) :: ssa_local, denom

  ! TL locals
  real(jprb) :: exponential_tl, planck_prime_tl, coeff_tl, rt_factor_tl
  real(jprb) :: factor_tl, c1_tl, c2_tl, scaling1_tl, scaling2_tl
  real(jprb) :: a, a_tl, d, d_tl

  integer(jpim) :: max_reg
  integer(jpim) :: jlev, jspec, jreg

  real(jphook) :: hook_handle

  if (lhook) call dr_hook('tcrad:calc_radiance_trans_source_tl',0,hook_handle)

  ! Geometry is not perturbed
  secant = 1.0_jprb / mu

  transmittance_tl(:,:,:) = 0.0_jprb
  if (present(source_up_tl)) source_up_tl(:,:,:) = 0.0_jprb
  if (present(source_dn_tl)) source_dn_tl(:,:,:) = 0.0_jprb

  do jlev = 1,nlev

    if (region_fracs(1,jlev) < 1.0_jprb) then
      max_reg = nreg
    else
      max_reg = 1
    end if

    if (max_reg > 1) then
      ! Cloudy layer: Planck terms are scaled by region fraction.
      planck_top(:,1) = planck_hl(:,jlev) * region_fracs(1,jlev)
      planck_top_tl(:,1) = planck_hl_tl(:,jlev) * region_fracs(1,jlev) &
           &            + planck_hl(:,jlev)    * region_fracs_tl(1,jlev)

      ! Note: the (1 - 0*ssa) factor is unity; no TL w.r.t. ssa.
      planck_top(:,2:nreg) = spread(planck_hl(:,jlev),2,nreg-1) &
           &  * spread(region_fracs(2:nreg,jlev),1,nspec)
      planck_top_tl(:,2:nreg) = spread(planck_hl_tl(:,jlev),2,nreg-1) &
           &  * spread(region_fracs(2:nreg,jlev),1,nspec) &
           &  + spread(planck_hl(:,jlev),2,nreg-1) &
           &  * spread(region_fracs_tl(2:nreg,jlev),1,nspec)

      planck_base(:,1) = planck_hl(:,jlev+1) * region_fracs(1,jlev)
      planck_base_tl(:,1) = planck_hl_tl(:,jlev+1) * region_fracs(1,jlev) &
           &             + planck_hl(:,jlev+1)    * region_fracs_tl(1,jlev)

      planck_base(:,2:nreg) = spread(planck_hl(:,jlev+1),2,nreg-1) &
           &  * spread(region_fracs(2:nreg,jlev),1,nspec)
      planck_base_tl(:,2:nreg) = spread(planck_hl_tl(:,jlev+1),2,nreg-1) &
           &  * spread(region_fracs(2:nreg,jlev),1,nspec) &
           &  + spread(planck_hl(:,jlev+1),2,nreg-1) &
           &  * spread(region_fracs_tl(2:nreg,jlev),1,nspec)
    else
      ! Clear layer
      planck_top(:,1)  = planck_hl(:,jlev)
      planck_top_tl(:,1)  = planck_hl_tl(:,jlev)
      planck_base(:,1) = planck_hl(:,jlev+1)
      planck_base_tl(:,1) = planck_hl_tl(:,jlev+1)
    end if

    do jreg = 1,max_reg
      transmittance(:,jreg,jlev) = exp(-od(:,jreg,jlev)*secant)
      transmittance_tl(:,jreg,jlev) = transmittance(:,jreg,jlev) * (-secant) * od_tl(:,jreg,jlev)
    end do

    !-----------------------
    ! Clear region (jreg=1)
    !-----------------------
    if (present(source_up)) then
      do jspec = 1,nspec
        if (od(jspec,1,jlev) > OD_THRESH) then
          planck_prime = (planck_base(jspec,1)-planck_top(jspec,1)) / od(jspec,1,jlev)
          planck_prime_tl = ((planck_base_tl(jspec,1)-planck_top_tl(jspec,1)) * od(jspec,1,jlev) &
               &           - (planck_base(jspec,1)-planck_top(jspec,1)) * od_tl(jspec,1,jlev)) &
               &          / (od(jspec,1,jlev)*od(jspec,1,jlev))

          source_up(jspec,1,jlev) = planck_top(jspec,1) &
               &  - planck_base(jspec,1)*transmittance(jspec,1,jlev) &
               &  + planck_prime*mu*(1.0_jprb - transmittance(jspec,1,jlev))

          if (present(source_up_tl)) then
            source_up_tl(jspec,1,jlev) = planck_top_tl(jspec,1) &
                 &  - (planck_base_tl(jspec,1)*transmittance(jspec,1,jlev) &
                 &     + planck_base(jspec,1)*transmittance_tl(jspec,1,jlev)) &
                 &  + planck_prime_tl*mu*(1.0_jprb - transmittance(jspec,1,jlev)) &
                 &  + planck_prime*mu*(-transmittance_tl(jspec,1,jlev))
          end if
        else
          source_up(jspec,1,jlev) = od(jspec,1,jlev) * 0.5_jprb &
               &  * (planck_base(jspec,1)+planck_top(jspec,1)) / mu
          if (present(source_up_tl)) then
            source_up_tl(jspec,1,jlev) = od_tl(jspec,1,jlev) * 0.5_jprb &
                 &  * (planck_base(jspec,1)+planck_top(jspec,1)) / mu &
                 &  + od(jspec,1,jlev) * 0.5_jprb &
                 &  * (planck_base_tl(jspec,1)+planck_top_tl(jspec,1)) / mu
          end if
        end if
      end do
    end if

    if (present(source_dn)) then
      do jspec = 1,nspec
        if (od(jspec,1,jlev) > OD_THRESH) then
          planck_prime = (planck_base(jspec,1)-planck_top(jspec,1)) / od(jspec,1,jlev)
          planck_prime_tl = ((planck_base_tl(jspec,1)-planck_top_tl(jspec,1)) * od(jspec,1,jlev) &
               &           - (planck_base(jspec,1)-planck_top(jspec,1)) * od_tl(jspec,1,jlev)) &
               &          / (od(jspec,1,jlev)*od(jspec,1,jlev))

          source_dn(jspec,1,jlev) = planck_base(jspec,1) &
               &  - planck_top(jspec,1)*transmittance(jspec,1,jlev) &
               &  - planck_prime*mu*(1.0_jprb - transmittance(jspec,1,jlev))

          if (present(source_dn_tl)) then
            source_dn_tl(jspec,1,jlev) = planck_base_tl(jspec,1) &
                 &  - (planck_top_tl(jspec,1)*transmittance(jspec,1,jlev) &
                 &     + planck_top(jspec,1)*transmittance_tl(jspec,1,jlev)) &
                 &  - planck_prime_tl*mu*(1.0_jprb - transmittance(jspec,1,jlev)) &
                 &  - planck_prime*mu*(-transmittance_tl(jspec,1,jlev))
          end if
        else
          source_dn(jspec,1,jlev) = od(jspec,1,jlev) * 0.5_jprb &
               &  * (planck_base(jspec,1)+planck_top(jspec,1)) / mu
          if (present(source_dn_tl)) then
            source_dn_tl(jspec,1,jlev) = od_tl(jspec,1,jlev) * 0.5_jprb &
                 &  * (planck_base(jspec,1)+planck_top(jspec,1)) / mu &
                 &  + od(jspec,1,jlev) * 0.5_jprb &
                 &  * (planck_base_tl(jspec,1)+planck_top_tl(jspec,1)) / mu
          end if
        end if
      end do
    end if

    !-----------------------
    ! Cloudy regions (jreg>=2)
    !-----------------------
    do jreg = 2,max_reg
      do jspec = 1,nspec

        if (od(jspec,jreg,jlev) > OD_THRESH) then

          ! gamma1, gamma2, k_exponent, p_same, p_opposite depend only on microphysics/geometry: TL = 0
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
          planck_prime_tl = ((planck_base_tl(jspec,jreg)-planck_top_tl(jspec,jreg)) * od(jspec,jreg,jlev) &
               &           - (planck_base(jspec,jreg)-planck_top(jspec,jreg)) * od_tl(jspec,jreg,jlev)) &
               &          / (od(jspec,jreg,jlev)*od(jspec,jreg,jlev))

          exponential = exp(-k_exponent*od(jspec,jreg,jlev))
          exponential_tl = exponential * (-k_exponent) * od_tl(jspec,jreg,jlev)

          coeff = planck_prime / (gamma1+gamma2)
          coeff_tl = planck_prime_tl / (gamma1+gamma2)

          d = (k_exponent + gamma1) + (k_exponent-gamma1) * exponential*exponential
          d_tl = (k_exponent-gamma1) * 2.0_jprb * exponential * exponential_tl
          rt_factor = 1.0_jprb / d
          rt_factor_tl = -rt_factor*rt_factor * d_tl

          factor = exponential * gamma2 / (gamma1 + k_exponent)
          factor_tl = exponential_tl * gamma2 / (gamma1 + k_exponent)

          a = flux_up_base(jspec,jreg,jlev) - factor*flux_dn_top(jspec,jreg,jlev) &
               &  - (planck_base(jspec,jreg)+coeff) + factor*(planck_top(jspec,jreg)-coeff)
          a_tl = flux_up_base_tl(jspec,jreg,jlev) &
               &  - (factor_tl*flux_dn_top(jspec,jreg,jlev) + factor*flux_dn_top_tl(jspec,jreg,jlev)) &
               &  - (planck_base_tl(jspec,jreg)+coeff_tl) &
               &  + (factor_tl*(planck_top(jspec,jreg)-coeff) + factor*(planck_top_tl(jspec,jreg)-coeff_tl))
          c1 = rt_factor * a
          c1_tl = rt_factor_tl * a + rt_factor * a_tl

          a = flux_dn_top(jspec,jreg,jlev) - factor*flux_up_base(jspec,jreg,jlev) &
               &  - (planck_top(jspec,jreg)-coeff) + factor*(planck_base(jspec,jreg)+coeff)
          a_tl = flux_dn_top_tl(jspec,jreg,jlev) &
               &  - (factor_tl*flux_up_base(jspec,jreg,jlev) + factor*flux_up_base_tl(jspec,jreg,jlev)) &
               &  - (planck_top_tl(jspec,jreg)-coeff_tl) &
               &  + (factor_tl*(planck_base(jspec,jreg)+coeff) + factor*(planck_base_tl(jspec,jreg)+coeff_tl))
          c2 = rt_factor * a
          c2_tl = rt_factor_tl * a + rt_factor * a_tl

          one_minus_kmu = 1.0_jprb - k_exponent*mu
          denom = merge(one_minus_kmu, epsilon(1.0_jprb), abs(one_minus_kmu)>epsilon(1.0_jprb))
          scaling1 = (exponential - transmittance(jspec,jreg,jlev)) / denom
          scaling1_tl = (exponential_tl - transmittance_tl(jspec,jreg,jlev)) / denom

          scaling2 = (1.0_jprb - exponential*transmittance(jspec,jreg,jlev)) / (1.0_jprb + k_exponent*mu)
          scaling2_tl = (-(exponential_tl*transmittance(jspec,jreg,jlev) + exponential*transmittance_tl(jspec,jreg,jlev))) &
               &      / (1.0_jprb + k_exponent*mu)

          if (present(source_up)) then
            source_up(jspec,jreg,jlev) = 0.5_jprb*ssa(jspec,jreg,jlev)*(1.0_jprb - transmittance(jspec,jreg,jlev)) &
                 &    * planck_prime*(p_same-p_opposite)/(gamma1+gamma2) &
                 &  + planck_top(jspec,jreg)-planck_base(jspec,jreg)*transmittance(jspec,jreg,jlev) &
                 &  + planck_prime*mu*(1.0_jprb - transmittance(jspec,jreg,jlev))

            source_up(jspec,jreg,jlev) = source_up(jspec,jreg,jlev) &
                 &  + 0.5_jprb*ssa(jspec,jreg,jlev) &
                 &  * (p_same     * ((gamma1+k_exponent)*scaling1*c1 +  gamma2*scaling2*c2) &
                 &    +p_opposite * ( gamma2*scaling1*c1             + (gamma1+k_exponent)*scaling2*c2))

            if (present(source_up_tl)) then
              ssa_local = ssa(jspec,jreg,jlev)

              ! TL of the "direct emission + linear structure scattering" part
              source_up_tl(jspec,jreg,jlev) = 0.5_jprb*ssa_local * ( &
                   &   (-transmittance_tl(jspec,jreg,jlev))*planck_prime &
                   &   + (1.0_jprb - transmittance(jspec,jreg,jlev))*planck_prime_tl ) &
                   &  * (p_same-p_opposite)/(gamma1+gamma2) &
                   &  + planck_top_tl(jspec,jreg) &
                   &  - (planck_base_tl(jspec,jreg)*transmittance(jspec,jreg,jlev) &
                   &     + planck_base(jspec,jreg)*transmittance_tl(jspec,jreg,jlev)) &
                   &  + planck_prime_tl*mu*(1.0_jprb - transmittance(jspec,jreg,jlev)) &
                   &  + planck_prime*mu*(-transmittance_tl(jspec,jreg,jlev))

              ! TL of the exponential scattering part
              source_up_tl(jspec,jreg,jlev) = source_up_tl(jspec,jreg,jlev) &
                   &  + 0.5_jprb*ssa_local * ( &
                   &  p_same * ( (gamma1+k_exponent) * (scaling1_tl*c1 + scaling1*c1_tl) &
                   &          +  gamma2           * (scaling2_tl*c2 + scaling2*c2_tl) ) &
                   & +p_opposite * ( gamma2 * (scaling1_tl*c1 + scaling1*c1_tl) &
                   &            + (gamma1+k_exponent) * (scaling2_tl*c2 + scaling2*c2_tl) ) )
            end if
          end if

          if (present(source_dn)) then
            source_dn(jspec,jreg,jlev) = - 0.5_jprb*ssa(jspec,jreg,jlev)*(1.0_jprb - transmittance(jspec,jreg,jlev)) &
                 &    * planck_prime*(p_same-p_opposite)/(gamma1+gamma2) &
                 &  + planck_base(jspec,jreg)-planck_top(jspec,jreg)*transmittance(jspec,jreg,jlev) &
                 &  - planck_prime*mu*(1.0_jprb - transmittance(jspec,jreg,jlev))

            source_dn(jspec,jreg,jlev) = source_dn(jspec,jreg,jlev) &
                 &  + 0.5_jprb*ssa(jspec,jreg,jlev) &
                 &  * (p_opposite * ((gamma1+k_exponent)*scaling2*c1 +  gamma2*scaling1*c2) &
                 &    +p_same     * ( gamma2*scaling2*c1             + (gamma1+k_exponent)*scaling1*c2))

            if (present(source_dn_tl)) then
              ssa_local = ssa(jspec,jreg,jlev)

              source_dn_tl(jspec,jreg,jlev) = -0.5_jprb*ssa_local * ( &
                   &   (-transmittance_tl(jspec,jreg,jlev))*planck_prime &
                   &   + (1.0_jprb - transmittance(jspec,jreg,jlev))*planck_prime_tl ) &
                   &  * (p_same-p_opposite)/(gamma1+gamma2) &
                   &  + planck_base_tl(jspec,jreg) &
                   &  - (planck_top_tl(jspec,jreg)*transmittance(jspec,jreg,jlev) &
                   &     + planck_top(jspec,jreg)*transmittance_tl(jspec,jreg,jlev)) &
                   &  - planck_prime_tl*mu*(1.0_jprb - transmittance(jspec,jreg,jlev)) &
                   &  - planck_prime*mu*(-transmittance_tl(jspec,jreg,jlev))

              source_dn_tl(jspec,jreg,jlev) = source_dn_tl(jspec,jreg,jlev) &
                   &  + 0.5_jprb*ssa_local * ( &
                   &  p_opposite * ( (gamma1+k_exponent) * (scaling2_tl*c1 + scaling2*c1_tl) &
                   &            +  gamma2           * (scaling1_tl*c2 + scaling1*c2_tl) ) &
                   & +p_same     * ( gamma2 * (scaling2_tl*c1 + scaling2*c1_tl) &
                   &            + (gamma1+k_exponent) * (scaling1_tl*c2 + scaling1*c2_tl) ) )
            end if
          end if

        else
          ! Low optical depth approximation: emission only
          if (present(source_up)) then
            source_up(jspec,jreg,jlev) = od(jspec,jreg,jlev) &
                 &  * 0.5_jprb*(planck_base(jspec,jreg)+planck_top(jspec,jreg)) / mu
            if (present(source_up_tl)) then
              source_up_tl(jspec,jreg,jlev) = od_tl(jspec,jreg,jlev) &
                   &  * 0.5_jprb*(planck_base(jspec,jreg)+planck_top(jspec,jreg)) / mu &
                   &  + od(jspec,jreg,jlev) &
                   &  * 0.5_jprb*(planck_base_tl(jspec,jreg)+planck_top_tl(jspec,jreg)) / mu
            end if
          end if
          if (present(source_dn)) then
            source_dn(jspec,jreg,jlev) = od(jspec,jreg,jlev) &
                 &  * 0.5_jprb*(planck_base(jspec,jreg)+planck_top(jspec,jreg)) / mu
            if (present(source_dn_tl)) then
              source_dn_tl(jspec,jreg,jlev) = od_tl(jspec,jreg,jlev) &
                   &  * 0.5_jprb*(planck_base(jspec,jreg)+planck_top(jspec,jreg)) / mu &
                   &  + od(jspec,jreg,jlev) &
                   &  * 0.5_jprb*(planck_base_tl(jspec,jreg)+planck_top_tl(jspec,jreg)) / mu
            end if
          end if
        end if

      end do
    end do

  end do

  if (lhook) call dr_hook('tcrad:calc_radiance_trans_source_tl',1,hook_handle)

end subroutine calc_radiance_trans_source_tl


! ===== FILE: tcrad_calc_radiance_up_tl.F90 =====
subroutine calc_radiance_up_tl(nspec, nlev, &
     &  weight, &
     &  surf_up, surf_up_tl, &
     &  transmittance, transmittance_tl, &
     &  source_up, source_up_tl, &
     &  u_overlap, u_overlap_tl, &
     &  radiance_up, radiance_up_tl)

  use parkind1, only           : jpim, jprb
  use yomhook,  only           : lhook, dr_hook, jphook

  implicit none

  ! Inputs
  integer(jpim), intent(in) :: nspec, nlev
  real(jprb), intent(in) :: weight

  ! Atmospheric/macrophysical inputs (TL w.r.t. these)
  real(jprb), intent(in),  dimension(nspec,NREGION) :: surf_up
  real(jprb), intent(in),  dimension(nspec,NREGION) :: surf_up_tl

  real(jprb), intent(in),  dimension(nspec,NREGION,nlev) :: transmittance
  real(jprb), intent(in),  dimension(nspec,NREGION,nlev) :: transmittance_tl

  real(jprb), intent(in),  dimension(nspec,NREGION,nlev) :: source_up
  real(jprb), intent(in),  dimension(nspec,NREGION,nlev) :: source_up_tl

  real(jprb), intent(in),  dimension(NREGION,NREGION,nlev+1) :: u_overlap
  real(jprb), intent(in),  dimension(NREGION,NREGION,nlev+1) :: u_overlap_tl

  ! Outputs (inout accumulation)
  real(jprb), intent(inout), dimension(nspec,nlev+1) :: radiance_up
  real(jprb), intent(inout), dimension(nspec,nlev+1) :: radiance_up_tl

  ! Local variables
  real(jprb), dimension(nspec,NREGION) :: radiance_base, radiance_top
  real(jprb), dimension(nspec,NREGION) :: radiance_base_tl, radiance_top_tl

  integer(jpim) :: jlev, jspec, jreg, jreg2

  real(jphook) :: hook_handle

  if (lhook) call dr_hook('tcrad:calc_radiance_up_tl',0,hook_handle)

  ! Surface upward radiance
  radiance_base    = weight * surf_up
  radiance_base_tl = weight * surf_up_tl

  ! Save radiance profile averaged over regions (add to existing)
  radiance_up(:,nlev+1)    = radiance_up(:,nlev+1)    + sum(radiance_base,2)
  radiance_up_tl(:,nlev+1) = radiance_up_tl(:,nlev+1) + sum(radiance_base_tl,2)

  do jlev = nlev,1,-1
    ! Schwarzschild solution in each region
    radiance_top = transmittance(:,:,jlev)*radiance_base + weight * source_up(:,:,jlev)
    radiance_top_tl = transmittance_tl(:,:,jlev)*radiance_base &
         &          + transmittance(:,:,jlev)*radiance_base_tl &
         &          + weight * source_up_tl(:,:,jlev)

    ! Apply overlap rules to obtain radiances at base of layer above:
    ! radiance_base = u_overlap(:,:,jlev) * radiance_top
    radiance_base(:,:)    = 0.0_jprb
    radiance_base_tl(:,:) = 0.0_jprb
    do jreg = 1,NREGION
      do jreg2 = 1,NREGION
        radiance_base(:,jreg) = radiance_base(:,jreg) &
             &  + u_overlap(jreg,jreg2,jlev) * radiance_top(:,jreg2)
        radiance_base_tl(:,jreg) = radiance_base_tl(:,jreg) &
             &  + u_overlap_tl(jreg,jreg2,jlev) * radiance_top(:,jreg2) &
             &  + u_overlap(jreg,jreg2,jlev) * radiance_top_tl(:,jreg2)
      end do
    end do

    ! Save radiances averaged over regions
    radiance_up(:,jlev)    = radiance_up(:,jlev)    + sum(radiance_base,2)
    radiance_up_tl(:,jlev) = radiance_up_tl(:,jlev) + sum(radiance_base_tl,2)
  end do

  if (lhook) call dr_hook('tcrad:calc_radiance_up_tl',1,hook_handle)

end subroutine calc_radiance_up_tl


! ===== FILE: tcrad_calc_reflectance_transmittance_tl.F90 =====
subroutine calc_reflectance_transmittance_tl(nspec, nlev, nreg, &
     &  region_fracs, region_fracs_tl, planck_hl, planck_hl_tl, &
     &  od, od_tl, ssa, asymmetry, &
     &  reflectance, reflectance_tl, transmittance, transmittance_tl, &
     &  source_up, source_up_tl, source_dn, source_dn_tl)

  use yomhook,  only           : lhook, dr_hook, jphook
  
  ! Inputs
  integer(jpim), intent(in) :: nspec, nlev, nreg

  ! Atmospheric/macrophysical inputs (TL w.r.t. these)
  real(jprb), intent(in), dimension(nreg,nlev) :: region_fracs
  real(jprb), intent(in), dimension(nreg,nlev) :: region_fracs_tl
  real(jprb), intent(in), dimension(nspec,nlev+1) :: planck_hl
  real(jprb), intent(in), dimension(nspec,nlev+1) :: planck_hl_tl
  real(jprb), intent(in), dimension(nspec,nreg,nlev) :: od
  real(jprb), intent(in), dimension(nspec,nreg,nlev) :: od_tl

  ! Microphysical inputs (NO TL w.r.t. these)
  real(jprb), intent(in), dimension(nspec,2:nreg,nlev) :: ssa
  real(jprb), intent(in), dimension(nspec,nlev) :: asymmetry

  ! Outputs (nonlinear + TL)
  real(jprb), intent(out), dimension(nspec,nreg,nlev) :: reflectance, transmittance
  real(jprb), intent(out), dimension(nspec,nreg,nlev) :: reflectance_tl, transmittance_tl
  real(jprb), intent(out), dimension(nspec,nreg,nlev) :: source_up, source_dn
  real(jprb), intent(out), dimension(nspec,nreg,nlev) :: source_up_tl, source_dn_tl

  ! Two-stream exchange coefficients
  real(jprb) :: gamma1, gamma2

  ! Working variables
  real(jprb) :: coeff, coeff_dn_base, coeff_dn_top, coeff_up_base, coeff_up_top
  real(jprb) :: factor, exponential, exponential2, k_exponent, reftrans_factor

  ! TL working variables
  real(jprb) :: coeff1, coeff1_tl
  real(jprb) :: exponential_tl, exponential2_tl, reftrans_factor_tl
  real(jprb) :: refl_tl, tran_tl
  real(jprb) :: denom, denom_tl

  ! Loop indices
  integer(jpim) :: jspec, jreg, jlev

  real(jphook) :: hook_handle

  if (lhook) call dr_hook('tcrad:calc_reflectance_transmittance_tl',0,hook_handle)

  ! Set cloudy regions to default values
  reflectance(:,2:,:)   = 0.0_jprb
  transmittance(:,2:,:) = 1.0_jprb
  source_up(:,2:,:)     = 0.0_jprb
  source_dn(:,2:,:)     = 0.0_jprb

  reflectance_tl(:,:,:)   = 0.0_jprb
  transmittance_tl(:,:,:) = 0.0_jprb
  source_up_tl(:,:,:)     = 0.0_jprb
  source_dn_tl(:,:,:)     = 0.0_jprb

  do jlev = 1,nlev
    ! Clear-sky region (no scattering)
    do jspec = 1,nspec
      reflectance(jspec,1,jlev) = 0.0_jprb
      reflectance_tl(jspec,1,jlev) = 0.0_jprb

      if (od(jspec,1,jlev) > OD_THRESH_2STREAM) then
        coeff1 = lw_diffusivity*od(jspec,1,jlev)
        coeff1_tl = lw_diffusivity*od_tl(jspec,1,jlev)

        transmittance(jspec,1,jlev) = exp(-coeff1)
        transmittance_tl(jspec,1,jlev) = transmittance(jspec,1,jlev) * (-coeff1_tl)

        coeff = (planck_hl(jspec,jlev+1)-planck_hl(jspec,jlev)) / coeff1
        coeff_up_top  =  coeff + planck_hl(jspec,jlev)
        coeff_up_base =  coeff + planck_hl(jspec,jlev+1)
        coeff_dn_top  = -coeff + planck_hl(jspec,jlev)
        coeff_dn_base = -coeff + planck_hl(jspec,jlev+1)

        ! TL of coeff
        denom = coeff1*coeff1
        denom_tl = 2.0_jprb*coeff1*coeff1_tl
        coeff1 = coeff  ! reuse variable name for clarity in TL below (coeff already computed)
        ! recompute TL explicitly with quotient rule:
        coeff1_tl = ((planck_hl_tl(jspec,jlev+1)-planck_hl_tl(jspec,jlev)) * (lw_diffusivity*od(jspec,1,jlev)) &
             &      - (planck_hl(jspec,jlev+1)-planck_hl(jspec,jlev)) * (lw_diffusivity*od_tl(jspec,1,jlev))) &
             &     / ((lw_diffusivity*od(jspec,1,jlev))*(lw_diffusivity*od(jspec,1,jlev)))

        ! Now coeff1 is TL of coeff (renamed for convenience)
        ! coeff_up_top_tl etc.
        source_up(jspec,1,jlev) = coeff_up_top &
             &  - transmittance(jspec,1,jlev) * coeff_up_base
        source_dn(jspec,1,jlev) = coeff_dn_base &
             &  - transmittance(jspec,1,jlev) * coeff_dn_top

        source_up_tl(jspec,1,jlev) = (coeff1_tl + planck_hl_tl(jspec,jlev)) &
             &  - (transmittance_tl(jspec,1,jlev) * coeff_up_base &
             &     + transmittance(jspec,1,jlev) * (coeff1_tl + planck_hl_tl(jspec,jlev+1)))
        source_dn_tl(jspec,1,jlev) = (-coeff1_tl + planck_hl_tl(jspec,jlev+1)) &
             &  - (transmittance_tl(jspec,1,jlev) * coeff_dn_top &
             &     + transmittance(jspec,1,jlev) * (-coeff1_tl + planck_hl_tl(jspec,jlev)))
      else
        ! Linear limit at low optical depth
        coeff = lw_diffusivity*od(jspec,1,jlev)
        transmittance(jspec,1,jlev) = 1.0_jprb - coeff
        transmittance_tl(jspec,1,jlev) = -lw_diffusivity*od_tl(jspec,1,jlev)

        source_up(jspec,1,jlev) = coeff * 0.5_jprb &
             &  * (planck_hl(jspec,jlev)+planck_hl(jspec,jlev+1))
        source_dn(jspec,1,jlev) = source_up(jspec,1,jlev)

        source_up_tl(jspec,1,jlev) = (lw_diffusivity*od_tl(jspec,1,jlev)) * 0.5_jprb &
             &  * (planck_hl(jspec,jlev)+planck_hl(jspec,jlev+1)) &
             &  + coeff * 0.5_jprb * (planck_hl_tl(jspec,jlev)+planck_hl_tl(jspec,jlev+1))
        source_dn_tl(jspec,1,jlev) = source_up_tl(jspec,1,jlev)
      end if

      ! Scale the sources by the area fraction of the region
      source_up_tl(jspec,1,jlev) = region_fracs_tl(1,jlev)*source_up(jspec,1,jlev) &
           &                     + region_fracs(1,jlev)  *source_up_tl(jspec,1,jlev)
      source_dn_tl(jspec,1,jlev) = region_fracs_tl(1,jlev)*source_dn(jspec,1,jlev) &
           &                     + region_fracs(1,jlev)  *source_dn_tl(jspec,1,jlev)

      source_up(jspec,1,jlev) = region_fracs(1,jlev)*source_up(jspec,1,jlev)
      source_dn(jspec,1,jlev) = region_fracs(1,jlev)*source_dn(jspec,1,jlev)
    end do

    if (region_fracs(1,jlev) < 1.0_jprb) then
      ! Cloudy regions
      do jreg = 2,nreg
        do jspec = 1,nspec

          ! gamma1/gamma2/k_exponent depend only on microphysics => TL=0
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

          k_exponent = sqrt(max((gamma1 - gamma2) * (gamma1 + gamma2), &
               1.E-12_jprb))

          if (od(jspec,jreg,jlev) > OD_THRESH_2STREAM) then
            exponential = exp(-k_exponent*od(jspec,jreg,jlev))
            exponential_tl = exponential * (-k_exponent) * od_tl(jspec,jreg,jlev)

            exponential2 = exponential*exponential
            exponential2_tl = 2.0_jprb*exponential*exponential_tl

            denom = (k_exponent + gamma1) + (k_exponent - gamma1)*exponential2
            denom_tl = (k_exponent - gamma1)*exponential2_tl
            reftrans_factor = 1.0_jprb / denom
            reftrans_factor_tl = -reftrans_factor*reftrans_factor * denom_tl

            reflectance(jspec,jreg,jlev) = gamma2 * (1.0_jprb - exponential2) * reftrans_factor
            refl_tl = gamma2 * ( (-exponential2_tl) * reftrans_factor &
                 &           + (1.0_jprb - exponential2) * reftrans_factor_tl )
            reflectance_tl(jspec,jreg,jlev) = refl_tl

            transmittance(jspec,jreg,jlev) = 2.0_jprb * k_exponent * exponential * reftrans_factor
            tran_tl = 2.0_jprb * k_exponent * (exponential_tl*reftrans_factor + exponential*reftrans_factor_tl)
            transmittance_tl(jspec,jreg,jlev) = tran_tl

            coeff = (planck_hl(jspec,jlev+1)-planck_hl(jspec,jlev)) &
                 / (od(jspec,jreg,jlev)*(gamma1+gamma2))
            coeff_up_top  =  coeff + planck_hl(jspec,jlev)
            coeff_up_base =  coeff + planck_hl(jspec,jlev+1)
            coeff_dn_top  = -coeff + planck_hl(jspec,jlev)
            coeff_dn_base = -coeff + planck_hl(jspec,jlev+1)

            ! TL of coeff (gamma1+gamma2 constant)
            denom = od(jspec,jreg,jlev)*(gamma1+gamma2)
            denom_tl = od_tl(jspec,jreg,jlev)*(gamma1+gamma2)
            coeff1_tl = ((planck_hl_tl(jspec,jlev+1)-planck_hl_tl(jspec,jlev))*denom &
                 &     - (planck_hl(jspec,jlev+1)-planck_hl(jspec,jlev))*denom_tl) &
                 &    / (denom*denom)

            source_up(jspec,jreg,jlev) = coeff_up_top &
                 &  - reflectance(jspec,jreg,jlev) * coeff_dn_top &
                 &  - transmittance(jspec,jreg,jlev) * coeff_up_base
            source_dn(jspec,jreg,jlev) = coeff_dn_base &
                 &  - reflectance(jspec,jreg,jlev) * coeff_up_base &
                 &  - transmittance(jspec,jreg,jlev) * coeff_dn_top

            source_up_tl(jspec,jreg,jlev) = (coeff1_tl + planck_hl_tl(jspec,jlev)) &
                 &  - (reflectance_tl(jspec,jreg,jlev)*coeff_dn_top &
                 &     + reflectance(jspec,jreg,jlev)*(-coeff1_tl + planck_hl_tl(jspec,jlev))) &
                 &  - (transmittance_tl(jspec,jreg,jlev)*coeff_up_base &
                 &     + transmittance(jspec,jreg,jlev)*(coeff1_tl + planck_hl_tl(jspec,jlev+1)))

            source_dn_tl(jspec,jreg,jlev) = (-coeff1_tl + planck_hl_tl(jspec,jlev+1)) &
                 &  - (reflectance_tl(jspec,jreg,jlev)*coeff_up_base &
                 &     + reflectance(jspec,jreg,jlev)*(coeff1_tl + planck_hl_tl(jspec,jlev+1))) &
                 &  - (transmittance_tl(jspec,jreg,jlev)*coeff_dn_top &
                 &     + transmittance(jspec,jreg,jlev)*(-coeff1_tl + planck_hl_tl(jspec,jlev)))
          else
            ! Low optical depth approximation
            reflectance(jspec,jreg,jlev) = gamma2 * od(jspec,jreg,jlev)
            reflectance_tl(jspec,jreg,jlev) = gamma2 * od_tl(jspec,jreg,jlev)

            transmittance(jspec,jreg,jlev) &
                 &  = (1.0_jprb - k_exponent*od(jspec,jreg,jlev)) &
                 &  / (1.0_jprb + od(jspec,jreg,jlev)*(gamma1-k_exponent))

            denom = 1.0_jprb + od(jspec,jreg,jlev)*(gamma1-k_exponent)
            denom_tl = od_tl(jspec,jreg,jlev)*(gamma1-k_exponent)
            transmittance_tl(jspec,jreg,jlev) = ( (-k_exponent*od_tl(jspec,jreg,jlev))*denom &
                 &  - (1.0_jprb - k_exponent*od(jspec,jreg,jlev))*denom_tl ) / (denom*denom)

            source_up(jspec,jreg,jlev) &
                 &  = (1.0_jprb - reflectance(jspec,jreg,jlev) &
                 &  - transmittance(jspec,jreg,jlev)) &
                 &    * 0.5_jprb * (planck_hl(jspec,jlev) + planck_hl(jspec,jlev+1))
            source_dn(jspec,jreg,jlev) = source_up(jspec,jreg,jlev)

            source_up_tl(jspec,jreg,jlev) = (-(reflectance_tl(jspec,jreg,jlev)+transmittance_tl(jspec,jreg,jlev))) &
                 &  * 0.5_jprb * (planck_hl(jspec,jlev) + planck_hl(jspec,jlev+1)) &
                 &  + (1.0_jprb - reflectance(jspec,jreg,jlev) - transmittance(jspec,jreg,jlev)) &
                 &  * 0.5_jprb * (planck_hl_tl(jspec,jlev) + planck_hl_tl(jspec,jlev+1))
            source_dn_tl(jspec,jreg,jlev) = source_up_tl(jspec,jreg,jlev)
          end if

          ! Scale the sources by the area fraction of the region
          source_up_tl(jspec,jreg,jlev) = region_fracs_tl(jreg,jlev)*source_up(jspec,jreg,jlev) &
               &                        + region_fracs(jreg,jlev)  *source_up_tl(jspec,jreg,jlev)
          source_dn_tl(jspec,jreg,jlev) = region_fracs_tl(jreg,jlev)*source_dn(jspec,jreg,jlev) &
               &                        + region_fracs(jreg,jlev)  *source_dn_tl(jspec,jreg,jlev)

          source_up(jspec,jreg,jlev) = region_fracs(jreg,jlev)*source_up(jspec,jreg,jlev)
          source_dn(jspec,jreg,jlev) = region_fracs(jreg,jlev)*source_dn(jspec,jreg,jlev)

        end do
      end do
    end if
  end do

  if (lhook) call dr_hook('tcrad:calc_reflectance_transmittance_tl',1,hook_handle)

end subroutine calc_reflectance_transmittance_tl


! ===== FILE: tcrad_calc_region_properties_tl.F90 =====
subroutine calc_region_properties_tl(nlev, &
     &  cloud_fraction, cloud_fraction_tl, do_gamma, frac_std, frac_std_tl, &
     &  reg_fracs, reg_fracs_tl, od_scaling, od_scaling_tl, cloud_fraction_threshold)

  use parkind1,     only : jprb
  use yomhook, only : lhook, dr_hook, jphook

  ! Minimum od_scaling in the case of a Gamma distribution
  real(jprb), parameter :: MinGammaODScaling = 0.025_jprb

  real(jprb), parameter :: MinLowerFrac      = 0.5_jprb
  real(jprb), parameter :: MaxLowerFrac      = 0.9_jprb
  real(jprb), parameter :: FSDAtMinLowerFrac = 1.5_jprb
  real(jprb), parameter :: FSDAtMaxLowerFrac = 3.725_jprb
  real(jprb), parameter :: LowerFracFSDGradient &
       &  = (MaxLowerFrac-MinLowerFrac) / (FSDAtMaxLowerFrac-FSDAtMinLowerFrac)
  real(jprb), parameter :: LowerFracFSDIntercept &
       &  = MinLowerFrac - FSDAtMinLowerFrac*LowerFracFSDGradient

  integer, intent(in) :: nlev
  logical, intent(in) :: do_gamma

  ! Atmospheric/macrophysical inputs (TL w.r.t. these)
  real(jprb), intent(in), dimension(:)  :: cloud_fraction ! (nlev)
  real(jprb), intent(in), dimension(:)  :: cloud_fraction_tl ! (nlev)
  real(jprb), intent(in), dimension(:)  :: frac_std       ! (nlev)
  real(jprb), intent(in), dimension(:)  :: frac_std_tl    ! (nlev)

  ! Outputs (nonlinear + TL)
  real(jprb), intent(out) :: reg_fracs(NREGION,nlev)
  real(jprb), intent(out) :: reg_fracs_tl(NREGION,nlev)
  real(jprb), intent(out) :: od_scaling(2:NREGION,nlev)
  real(jprb), intent(out) :: od_scaling_tl(2:NREGION,nlev)

  real(jprb), intent(in), optional :: cloud_fraction_threshold

  real(jprb) :: frac_threshold
  integer :: jlev
  real(jphook) :: hook_handle

  ! Locals for derivatives
  real(jprb) :: x, x_tl, logx, logx_tl, y, y_tl, z, z_tl
  real(jprb) :: expmy, expmy_tl
  real(jprb) :: lower, lower_clamped, lower_tl
  real(jprb) :: t, t_tl, inner1, inner1_tl, inner2, inner2_tl
  real(jprb) :: expt, expt_tl
  real(jprb) :: num, num_tl, denom, denom_tl

  if (lhook) call dr_hook('tcrad:calc_region_properties_tl',0,hook_handle)

  if (present(cloud_fraction_threshold)) then
    frac_threshold = cloud_fraction_threshold
  else
    frac_threshold = 10.0_jprb * epsilon(1.0_jprb)
  end if

  reg_fracs_tl(:,:)     = 0.0_jprb
  od_scaling_tl(:,:)    = 0.0_jprb

  if (.not. do_gamma) then
    ! Lognormal interpretation
    do jlev = 1,nlev
      if (cloud_fraction(jlev) < frac_threshold) then
        reg_fracs(1,jlev)    = 1.0_jprb
        reg_fracs(2:3,jlev)  = 0.0_jprb
        od_scaling(2:3,jlev) = 1.0_jprb

        reg_fracs_tl(:,jlev)   = 0.0_jprb
        od_scaling_tl(2:3,jlev)= 0.0_jprb
      else
        reg_fracs(1,jlev) = 1.0_jprb - cloud_fraction(jlev)
        reg_fracs(2:3,jlev) = cloud_fraction(jlev)*0.5_jprb

        reg_fracs_tl(1,jlev) = -cloud_fraction_tl(jlev)
        reg_fracs_tl(2,jlev) = 0.5_jprb*cloud_fraction_tl(jlev)
        reg_fracs_tl(3,jlev) = 0.5_jprb*cloud_fraction_tl(jlev)

        ! od_scaling(2) = exp(-sqrt(log(fsd^2+1))) / sqrt(fsd^2+1)
        x = frac_std(jlev)*frac_std(jlev) + 1.0_jprb
        x_tl = 2.0_jprb*frac_std(jlev)*frac_std_tl(jlev)

        logx = log(x)
        logx_tl = x_tl / x

        y = sqrt(logx)
        if (y > 0.0_jprb) then
          y_tl = 0.5_jprb * logx_tl / y
        else
          y_tl = 0.0_jprb
        end if

        z = sqrt(x)
        if (z > 0.0_jprb) then
          z_tl = 0.5_jprb * x_tl / z
        else
          z_tl = 0.0_jprb
        end if

        expmy = exp(-y)
        expmy_tl = expmy * (-y_tl)

        od_scaling(2,jlev) = expmy / z
        od_scaling_tl(2,jlev) = expmy_tl / z - expmy * z_tl / (z*z)

        od_scaling(3,jlev) = 2.0_jprb - od_scaling(2,jlev)
        od_scaling_tl(3,jlev) = -od_scaling_tl(2,jlev)
      end if
    end do
  else
    ! Gamma interpretation
    do jlev = 1,nlev
      if (cloud_fraction(jlev) < frac_threshold) then
        reg_fracs(1,jlev)    = 1.0_jprb
        reg_fracs(2:3,jlev)  = 0.0_jprb
        od_scaling(2:3,jlev) = 1.0_jprb

        reg_fracs_tl(:,jlev)   = 0.0_jprb
        od_scaling_tl(2:3,jlev)= 0.0_jprb
      else
        ! Clear-sky fraction
        reg_fracs(1,jlev) = 1.0_jprb - cloud_fraction(jlev)
        reg_fracs_tl(1,jlev) = -cloud_fraction_tl(jlev)

        ! Lower cloudy-region fraction (clamped function of FSD)
        lower = LowerFracFSDIntercept + frac_std(jlev)*LowerFracFSDGradient
        lower_clamped = max(MinLowerFrac, min(MaxLowerFrac, lower))

        if (lower > MinLowerFrac .and. lower < MaxLowerFrac) then
          lower_tl = frac_std_tl(jlev) * LowerFracFSDGradient
        else
          lower_tl = 0.0_jprb
        end if

        reg_fracs(2,jlev) = cloud_fraction(jlev) * lower_clamped
        reg_fracs_tl(2,jlev) = cloud_fraction_tl(jlev) * lower_clamped &
             &               + cloud_fraction(jlev) * lower_tl

        ! od_scaling(2) for gamma case
        inner1 = 1.0_jprb + 0.5_jprb*frac_std(jlev)
        inner1_tl = 0.5_jprb*frac_std_tl(jlev)

        inner2 = 1.0_jprb + 0.5_jprb*frac_std(jlev)*inner1
        inner2_tl = 0.5_jprb*(frac_std_tl(jlev)*inner1 + frac_std(jlev)*inner1_tl)

        t = -frac_std(jlev) * inner2
        t_tl = -(frac_std_tl(jlev)*inner2 + frac_std(jlev)*inner2_tl)

        expt = exp(t)
        expt_tl = expt * t_tl

        od_scaling(2,jlev) = MinGammaODScaling &
             &  + (1.0_jprb - MinGammaODScaling) * expt
        od_scaling_tl(2,jlev) = (1.0_jprb - MinGammaODScaling) * expt_tl

        ! Upper cloudy-region fraction
        reg_fracs(3,jlev) = 1.0_jprb - reg_fracs(1,jlev) - reg_fracs(2,jlev)
        reg_fracs_tl(3,jlev) = -reg_fracs_tl(1,jlev) - reg_fracs_tl(2,jlev)

        ! Ensure conservation of the mean optical depth
        num = cloud_fraction(jlev) - reg_fracs(2,jlev)*od_scaling(2,jlev)
        num_tl = cloud_fraction_tl(jlev) &
             &   - (reg_fracs_tl(2,jlev)*od_scaling(2,jlev) + reg_fracs(2,jlev)*od_scaling_tl(2,jlev))

        denom = reg_fracs(3,jlev)
        denom_tl = reg_fracs_tl(3,jlev)

        od_scaling(3,jlev) = num / denom
        od_scaling_tl(3,jlev) = (num_tl*denom - num*denom_tl) / (denom*denom)
      end if
    end do
  end if

  if (lhook) call dr_hook('tcrad:calc_region_properties_tl',1,hook_handle)

end subroutine calc_region_properties_tl


! ===== FILE: tcrad_calc_two_stream_flux_tl.F90 =====
subroutine calc_two_stream_flux_tl(nspec, nlev, &
     &  surf_emission, surf_emission_tl, surf_albedo, surf_albedo_tl, &
     &  reflectance, reflectance_tl, transmittance, transmittance_tl, &
     &  source_up, source_up_tl, source_dn, source_dn_tl, &
     &  is_cloud_free_layer, u_overlap, u_overlap_tl, v_overlap, v_overlap_tl, &
     &  flux_up_base, flux_up_base_tl, flux_dn_base, flux_dn_base_tl, &
     &  flux_up_top, flux_up_top_tl, flux_dn_top, flux_dn_top_tl)

  use parkind1, only           : jpim, jprb
  use yomhook,  only           : lhook, dr_hook, jphook

  implicit none

  ! Inputs
  integer(jpim), intent(in) :: nspec, nlev

  ! Atmospheric/macrophysical inputs (TL w.r.t. these)
  real(jprb), intent(in),  dimension(nspec) :: surf_emission, surf_albedo
  real(jprb), intent(in),  dimension(nspec) :: surf_emission_tl, surf_albedo_tl

  real(jprb), intent(in),  dimension(nspec,NREGION,nlev) :: reflectance, transmittance
  real(jprb), intent(in),  dimension(nspec,NREGION,nlev) :: source_up, source_dn

  real(jprb), intent(in),  dimension(nspec,NREGION,nlev) :: reflectance_tl, transmittance_tl
  real(jprb), intent(in),  dimension(nspec,NREGION,nlev) :: source_up_tl, source_dn_tl

  logical,    intent(in) :: is_cloud_free_layer(0:nlev+1)

  real(jprb), intent(in),  dimension(NREGION,NREGION,nlev+1) :: u_overlap, v_overlap
  real(jprb), intent(in),  dimension(NREGION,NREGION,nlev+1) :: u_overlap_tl, v_overlap_tl

  ! Outputs (nonlinear + TL)
  real(jprb), intent(out), dimension(nspec,NREGION,nlev) :: flux_up_base, flux_dn_base
  real(jprb), intent(out), dimension(nspec,NREGION,nlev) :: flux_up_top, flux_dn_top
  real(jprb), intent(out), dimension(nspec,NREGION,nlev) :: flux_up_base_tl, flux_dn_base_tl
  real(jprb), intent(out), dimension(nspec,NREGION,nlev) :: flux_up_top_tl, flux_dn_top_tl

  ! Local variables
  real(jprb), dimension(nspec, NREGION, nlev+1) :: total_albedo, total_source
  real(jprb), dimension(nspec, NREGION, nlev+1) :: total_albedo_tl, total_source_tl

  real(jprb), dimension(nspec, NREGION) :: total_albedo_below, total_source_below
  real(jprb), dimension(nspec, NREGION) :: total_albedo_below_tl, total_source_below_tl

  real(jprb), dimension(nspec, NREGION) :: inv_denom, inv_denom_tl

  integer(jpim) :: icloudtop
  integer(jpim) :: jspec, jlev, jreg, jreg2

  real(jprb) :: denom, denom_tl, numer, numer_tl
  real(jprb) :: inner, inner_tl

  real(jphook) :: hook_handle

  if (lhook) call dr_hook('tcrad:calc_two_stream_flux_tl',0,hook_handle)

  ! --------------------------------------------------------
  ! Section 1: Prepare variables and arrays
  ! --------------------------------------------------------

  icloudtop = nlev
  do jlev = 1,nlev
    if (.not. is_cloud_free_layer(jlev)) then
      icloudtop = jlev
      exit
    end if
  end do

  ! Initialize outputs
  flux_up_base    = 0.0_jprb
  flux_dn_base    = 0.0_jprb
  flux_up_top     = 0.0_jprb
  flux_dn_top     = 0.0_jprb

  flux_up_base_tl = 0.0_jprb
  flux_dn_base_tl = 0.0_jprb
  flux_up_top_tl  = 0.0_jprb
  flux_dn_top_tl  = 0.0_jprb

  total_albedo    = 0.0_jprb
  total_source    = 0.0_jprb
  total_albedo_tl = 0.0_jprb
  total_source_tl = 0.0_jprb

  ! --------------------------------------------------------
  ! Section 2: Clear-sky downwelling fluxes from TOA to cloud top
  ! --------------------------------------------------------
  do jlev = 1,icloudtop
    if (jlev > 1) then
      flux_dn_top(:,1,jlev)    = flux_dn_base(:,1,jlev-1)
      flux_dn_top_tl(:,1,jlev) = flux_dn_base_tl(:,1,jlev-1)
    end if

    flux_dn_base(:,1,jlev) = source_dn(:,1,jlev) &
         &                 + transmittance(:,1,jlev)*flux_dn_top(:,1,jlev)

    flux_dn_base_tl(:,1,jlev) = source_dn_tl(:,1,jlev) &
         &                    + transmittance_tl(:,1,jlev)*flux_dn_top(:,1,jlev) &
         &                    + transmittance(:,1,jlev)*flux_dn_top_tl(:,1,jlev)
  end do

  ! --------------------------------------------------------
  ! Section 3: Compute total sources and albedos up to cloud top
  ! --------------------------------------------------------

  ! Surface boundary
  do jreg = 1,NREGION
    do jspec = 1,nspec
      total_source(jspec,jreg,nlev+1) = u_overlap(jreg,1,nlev+1)*surf_emission(jspec)
      total_source_tl(jspec,jreg,nlev+1) = u_overlap_tl(jreg,1,nlev+1)*surf_emission(jspec) &
           &                           + u_overlap(jreg,1,nlev+1)*surf_emission_tl(jspec)

      total_albedo(jspec,jreg,nlev+1) = surf_albedo(jspec)
      total_albedo_tl(jspec,jreg,nlev+1) = surf_albedo_tl(jspec)
    end do
  end do

  ! Work up from the surface using the Adding Method
  do jlev = nlev,icloudtop,-1

    total_albedo_below    = 0.0_jprb
    total_source_below    = 0.0_jprb
    total_albedo_below_tl = 0.0_jprb
    total_source_below_tl = 0.0_jprb

    ! inv_denom = 1 / (1 - total_albedo_next * reflectance)
    denom = 0.0_jprb ! dummy init

    if (is_cloud_free_layer(jlev)) then
      do jspec = 1,nspec
        denom = 1.0_jprb - total_albedo(jspec,1,jlev+1)*reflectance(jspec,1,jlev)
        denom_tl = -( total_albedo_tl(jspec,1,jlev+1)*reflectance(jspec,1,jlev) &
             &       + total_albedo(jspec,1,jlev+1)*reflectance_tl(jspec,1,jlev) )
        inv_denom(jspec,1) = 1.0_jprb / denom
        inv_denom_tl(jspec,1) = -inv_denom(jspec,1)*inv_denom(jspec,1) * denom_tl

        total_albedo_below(jspec,1) = reflectance(jspec,1,jlev) &
             &  + transmittance(jspec,1,jlev)*transmittance(jspec,1,jlev)*total_albedo(jspec,1,jlev+1) &
             &  * inv_denom(jspec,1)

        total_albedo_below_tl(jspec,1) = reflectance_tl(jspec,1,jlev) &
             &  + (2.0_jprb*transmittance(jspec,1,jlev)*transmittance_tl(jspec,1,jlev) &
             &     * total_albedo(jspec,1,jlev+1) &
             &     + transmittance(jspec,1,jlev)*transmittance(jspec,1,jlev)*total_albedo_tl(jspec,1,jlev+1)) &
             &     * inv_denom(jspec,1) &
             &  + transmittance(jspec,1,jlev)*transmittance(jspec,1,jlev)*total_albedo(jspec,1,jlev+1) &
             &     * inv_denom_tl(jspec,1)

        inner = total_source(jspec,1,jlev+1) + total_albedo(jspec,1,jlev+1)*source_dn(jspec,1,jlev)
        inner_tl = total_source_tl(jspec,1,jlev+1) &
             &   + total_albedo_tl(jspec,1,jlev+1)*source_dn(jspec,1,jlev) &
             &   + total_albedo(jspec,1,jlev+1)*source_dn_tl(jspec,1,jlev)

        total_source_below(jspec,1) = source_up(jspec,1,jlev) &
             &  + transmittance(jspec,1,jlev) * inner * inv_denom(jspec,1)

        total_source_below_tl(jspec,1) = source_up_tl(jspec,1,jlev) &
             &  + transmittance_tl(jspec,1,jlev) * inner * inv_denom(jspec,1) &
             &  + transmittance(jspec,1,jlev) * inner_tl * inv_denom(jspec,1) &
             &  + transmittance(jspec,1,jlev) * inner * inv_denom_tl(jspec,1)
      end do
    else
      do jreg = 1,NREGION
        do jspec = 1,nspec
          denom = 1.0_jprb - total_albedo(jspec,jreg,jlev+1)*reflectance(jspec,jreg,jlev)
          denom_tl = -( total_albedo_tl(jspec,jreg,jlev+1)*reflectance(jspec,jreg,jlev) &
               &       + total_albedo(jspec,jreg,jlev+1)*reflectance_tl(jspec,jreg,jlev) )
          inv_denom(jspec,jreg) = 1.0_jprb / denom
          inv_denom_tl(jspec,jreg) = -inv_denom(jspec,jreg)*inv_denom(jspec,jreg) * denom_tl

          total_albedo_below(jspec,jreg) = reflectance(jspec,jreg,jlev) &
               &  + transmittance(jspec,jreg,jlev)*transmittance(jspec,jreg,jlev)*total_albedo(jspec,jreg,jlev+1) &
               &  * inv_denom(jspec,jreg)

          total_albedo_below_tl(jspec,jreg) = reflectance_tl(jspec,jreg,jlev) &
               &  + (2.0_jprb*transmittance(jspec,jreg,jlev)*transmittance_tl(jspec,jreg,jlev) &
               &     * total_albedo(jspec,jreg,jlev+1) &
               &     + transmittance(jspec,jreg,jlev)*transmittance(jspec,jreg,jlev)*total_albedo_tl(jspec,jreg,jlev+1)) &
               &     * inv_denom(jspec,jreg) &
               &  + transmittance(jspec,jreg,jlev)*transmittance(jspec,jreg,jlev)*total_albedo(jspec,jreg,jlev+1) &
               &     * inv_denom_tl(jspec,jreg)

          inner = total_source(jspec,jreg,jlev+1) + total_albedo(jspec,jreg,jlev+1)*source_dn(jspec,jreg,jlev)
          inner_tl = total_source_tl(jspec,jreg,jlev+1) &
               &   + total_albedo_tl(jspec,jreg,jlev+1)*source_dn(jspec,jreg,jlev) &
               &   + total_albedo(jspec,jreg,jlev+1)*source_dn_tl(jspec,jreg,jlev)

          total_source_below(jspec,jreg) = source_up(jspec,jreg,jlev) &
               &  + transmittance(jspec,jreg,jlev) * inner * inv_denom(jspec,jreg)

          total_source_below_tl(jspec,jreg) = source_up_tl(jspec,jreg,jlev) &
               &  + transmittance_tl(jspec,jreg,jlev) * inner * inv_denom(jspec,jreg) &
               &  + transmittance(jspec,jreg,jlev) * inner_tl * inv_denom(jspec,jreg) &
               &  + transmittance(jspec,jreg,jlev) * inner * inv_denom_tl(jspec,jreg)
        end do
      end do
    end if

    ! Convert "below interface" to "above interface" with overlap rules
    if (is_cloud_free_layer(jlev) .and. is_cloud_free_layer(jlev-1)) then
      total_albedo(:,:,jlev)    = total_albedo_below(:,:)
      total_source(:,:,jlev)    = total_source_below(:,:)
      total_albedo_tl(:,:,jlev) = total_albedo_below_tl(:,:)
      total_source_tl(:,:,jlev) = total_source_below_tl(:,:)
    else
      ! total_source(:, :, jlev) = u_overlap(:,:,jlev) * total_source_below
      total_source(:,:,jlev)    = 0.0_jprb
      total_source_tl(:,:,jlev) = 0.0_jprb
      do jspec = 1,nspec
        do jreg = 1,NREGION
          do jreg2 = 1,NREGION
            total_source(jspec,jreg,jlev) = total_source(jspec,jreg,jlev) &
                 &  + u_overlap(jreg,jreg2,jlev) * total_source_below(jspec,jreg2)
            total_source_tl(jspec,jreg,jlev) = total_source_tl(jspec,jreg,jlev) &
                 &  + u_overlap_tl(jreg,jreg2,jlev) * total_source_below(jspec,jreg2) &
                 &  + u_overlap(jreg,jreg2,jlev) * total_source_below_tl(jspec,jreg2)
          end do
        end do
      end do

      ! total_albedo(:,jreg,jlev) = sum_jreg2 total_albedo_below(:,jreg2) * v_overlap(jreg2,jreg,jlev)
      total_albedo(:,:,jlev)    = 0.0_jprb
      total_albedo_tl(:,:,jlev) = 0.0_jprb
      do jreg = 1,NREGION
        do jreg2 = 1,NREGION
          total_albedo(:,jreg,jlev) = total_albedo(:,jreg,jlev) &
               &  + total_albedo_below(:,jreg2) * v_overlap(jreg2,jreg,jlev)
          total_albedo_tl(:,jreg,jlev) = total_albedo_tl(:,jreg,jlev) &
               &  + total_albedo_below_tl(:,jreg2) * v_overlap(jreg2,jreg,jlev) &
               &  + total_albedo_below(:,jreg2) * v_overlap_tl(jreg2,jreg,jlev)
        end do
      end do
    end if

  end do

  ! --------------------------------------------------------
  ! Section 4: Compute upwelling fluxes up to top-of-atmosphere
  ! --------------------------------------------------------

  if (icloudtop > 1) then
    flux_up_base(:,1,icloudtop-1) = total_source(:,1,icloudtop) &
         &  + total_albedo(:,1,icloudtop)*flux_dn_base(:,1,icloudtop-1)

    flux_up_base_tl(:,1,icloudtop-1) = total_source_tl(:,1,icloudtop) &
         &  + total_albedo_tl(:,1,icloudtop)*flux_dn_base(:,1,icloudtop-1) &
         &  + total_albedo(:,1,icloudtop)*flux_dn_base_tl(:,1,icloudtop-1)

    flux_up_top(:,1,icloudtop-1) = source_up(:,1,icloudtop-1) &
         &  + transmittance(:,1,icloudtop-1)*flux_up_base(:,1,icloudtop-1)

    flux_up_top_tl(:,1,icloudtop-1) = source_up_tl(:,1,icloudtop-1) &
         &  + transmittance_tl(:,1,icloudtop-1)*flux_up_base(:,1,icloudtop-1) &
         &  + transmittance(:,1,icloudtop-1)*flux_up_base_tl(:,1,icloudtop-1)

    do jlev = icloudtop-2,1,-1
      flux_up_base(:,1,jlev)    = flux_up_top(:,1,jlev+1)
      flux_up_base_tl(:,1,jlev) = flux_up_top_tl(:,1,jlev+1)

      flux_up_top(:,1,jlev) = source_up(:,1,jlev) &
           &                + transmittance(:,1,jlev)*flux_up_base(:,1,jlev)

      flux_up_top_tl(:,1,jlev) = source_up_tl(:,1,jlev) &
           &                   + transmittance_tl(:,1,jlev)*flux_up_base(:,1,jlev) &
           &                   + transmittance(:,1,jlev)*flux_up_base_tl(:,1,jlev)
    end do
  end if

  ! --------------------------------------------------------
  ! Section 5: Compute fluxes from cloud top down to surface
  ! --------------------------------------------------------

  if (icloudtop > 1) then
    do jreg = 1,NREGION
      flux_dn_top(:,jreg,icloudtop) = v_overlap(jreg,1,icloudtop)*flux_dn_base(:,1,icloudtop-1)
      flux_dn_top_tl(:,jreg,icloudtop) = v_overlap_tl(jreg,1,icloudtop)*flux_dn_base(:,1,icloudtop-1) &
           &                          + v_overlap(jreg,1,icloudtop)*flux_dn_base_tl(:,1,icloudtop-1)
    end do
  end if

  do jlev = icloudtop,nlev

    if (is_cloud_free_layer(jlev)) then
      do jspec = 1,nspec

        denom = 1.0_jprb - reflectance(jspec,1,jlev)*total_albedo(jspec,1,jlev+1)
        denom_tl = -( reflectance_tl(jspec,1,jlev)*total_albedo(jspec,1,jlev+1) &
             &       + reflectance(jspec,1,jlev)*total_albedo_tl(jspec,1,jlev+1) )

        numer = transmittance(jspec,1,jlev)*flux_dn_top(jspec,1,jlev) &
             &  + reflectance(jspec,1,jlev)*total_source(jspec,1,jlev+1) &
             &  + source_dn(jspec,1,jlev)

        numer_tl = transmittance_tl(jspec,1,jlev)*flux_dn_top(jspec,1,jlev) &
             &    + transmittance(jspec,1,jlev)*flux_dn_top_tl(jspec,1,jlev) &
             &    + reflectance_tl(jspec,1,jlev)*total_source(jspec,1,jlev+1) &
             &    + reflectance(jspec,1,jlev)*total_source_tl(jspec,1,jlev+1) &
             &    + source_dn_tl(jspec,1,jlev)

        flux_dn_base(jspec,1,jlev) = numer / denom
        flux_dn_base_tl(jspec,1,jlev) = (numer_tl*denom - numer*denom_tl) / (denom*denom)

        flux_up_base(jspec,1,jlev) = total_source(jspec,1,jlev+1) &
             &  + flux_dn_base(jspec,1,jlev)*total_albedo(jspec,1,jlev+1)

        flux_up_base_tl(jspec,1,jlev) = total_source_tl(jspec,1,jlev+1) &
             &  + flux_dn_base_tl(jspec,1,jlev)*total_albedo(jspec,1,jlev+1) &
             &  + flux_dn_base(jspec,1,jlev)*total_albedo_tl(jspec,1,jlev+1)

        flux_up_top(jspec,1,jlev) = flux_up_base(jspec,1,jlev)*transmittance(jspec,1,jlev) &
             &  + source_up(jspec,1,jlev) + flux_dn_top(jspec,1,jlev)*reflectance(jspec,1,jlev)

        flux_up_top_tl(jspec,1,jlev) = flux_up_base_tl(jspec,1,jlev)*transmittance(jspec,1,jlev) &
             &  + flux_up_base(jspec,1,jlev)*transmittance_tl(jspec,1,jlev) &
             &  + source_up_tl(jspec,1,jlev) &
             &  + flux_dn_top_tl(jspec,1,jlev)*reflectance(jspec,1,jlev) &
             &  + flux_dn_top(jspec,1,jlev)*reflectance_tl(jspec,1,jlev)
      end do

      ! Other regions unused in clear layer (keep zero TL)
      flux_dn_base(:,2:NREGION,jlev)    = 0.0_jprb
      flux_dn_base_tl(:,2:NREGION,jlev) = 0.0_jprb
      flux_up_base(:,2:NREGION,jlev)    = 0.0_jprb
      flux_up_base_tl(:,2:NREGION,jlev) = 0.0_jprb
      flux_up_top(:,2:NREGION,jlev)     = 0.0_jprb
      flux_up_top_tl(:,2:NREGION,jlev)  = 0.0_jprb

    else

      do jreg = 1,NREGION
        do jspec = 1,nspec
          denom = 1.0_jprb - reflectance(jspec,jreg,jlev)*total_albedo(jspec,jreg,jlev+1)
          denom_tl = -( reflectance_tl(jspec,jreg,jlev)*total_albedo(jspec,jreg,jlev+1) &
               &       + reflectance(jspec,jreg,jlev)*total_albedo_tl(jspec,jreg,jlev+1) )

          numer = transmittance(jspec,jreg,jlev)*flux_dn_top(jspec,jreg,jlev) &
               &  + reflectance(jspec,jreg,jlev)*total_source(jspec,jreg,jlev+1) &
               &  + source_dn(jspec,jreg,jlev)

          numer_tl = transmittance_tl(jspec,jreg,jlev)*flux_dn_top(jspec,jreg,jlev) &
               &    + transmittance(jspec,jreg,jlev)*flux_dn_top_tl(jspec,jreg,jlev) &
               &    + reflectance_tl(jspec,jreg,jlev)*total_source(jspec,jreg,jlev+1) &
               &    + reflectance(jspec,jreg,jlev)*total_source_tl(jspec,jreg,jlev+1) &
               &    + source_dn_tl(jspec,jreg,jlev)

          flux_dn_base(jspec,jreg,jlev) = numer / denom
          flux_dn_base_tl(jspec,jreg,jlev) = (numer_tl*denom - numer*denom_tl) / (denom*denom)

          flux_up_base(jspec,jreg,jlev) = total_source(jspec,jreg,jlev+1) &
               &  + flux_dn_base(jspec,jreg,jlev)*total_albedo(jspec,jreg,jlev+1)

          flux_up_base_tl(jspec,jreg,jlev) = total_source_tl(jspec,jreg,jlev+1) &
               &  + flux_dn_base_tl(jspec,jreg,jlev)*total_albedo(jspec,jreg,jlev+1) &
               &  + flux_dn_base(jspec,jreg,jlev)*total_albedo_tl(jspec,jreg,jlev+1)

          flux_up_top(jspec,jreg,jlev) = flux_up_base(jspec,jreg,jlev)*transmittance(jspec,jreg,jlev) &
               &  + source_up(jspec,jreg,jlev) + flux_dn_top(jspec,jreg,jlev)*reflectance(jspec,jreg,jlev)

          flux_up_top_tl(jspec,jreg,jlev) = flux_up_base_tl(jspec,jreg,jlev)*transmittance(jspec,jreg,jlev) &
               &  + flux_up_base(jspec,jreg,jlev)*transmittance_tl(jspec,jreg,jlev) &
               &  + source_up_tl(jspec,jreg,jlev) &
               &  + flux_dn_top_tl(jspec,jreg,jlev)*reflectance(jspec,jreg,jlev) &
               &  + flux_dn_top(jspec,jreg,jlev)*reflectance_tl(jspec,jreg,jlev)
        end do
      end do

    end if

    if (jlev < nlev) then
      if (.not. (is_cloud_free_layer(jlev) .and. is_cloud_free_layer(jlev+1))) then
        flux_dn_top(:,:,jlev+1)    = 0.0_jprb
        flux_dn_top_tl(:,:,jlev+1) = 0.0_jprb
        do jspec = 1,nspec
          do jreg = 1,NREGION
            do jreg2 = 1,NREGION
              flux_dn_top(jspec,jreg,jlev+1) = flux_dn_top(jspec,jreg,jlev+1) &
                   &  + v_overlap(jreg,jreg2,jlev+1) * flux_dn_base(jspec,jreg2,jlev)
              flux_dn_top_tl(jspec,jreg,jlev+1) = flux_dn_top_tl(jspec,jreg,jlev+1) &
                   &  + v_overlap_tl(jreg,jreg2,jlev+1) * flux_dn_base(jspec,jreg2,jlev) &
                   &  + v_overlap(jreg,jreg2,jlev+1) * flux_dn_base_tl(jspec,jreg2,jlev)
            end do
          end do
        end do
      else
        flux_dn_top(:,1,jlev+1)    = flux_dn_base(:,1,jlev)
        flux_dn_top_tl(:,1,jlev+1) = flux_dn_base_tl(:,1,jlev)
        flux_dn_top(:,2:NREGION,jlev+1)    = 0.0_jprb
        flux_dn_top_tl(:,2:NREGION,jlev+1) = 0.0_jprb
      end if
    end if

  end do

  if (lhook) call dr_hook('tcrad:calc_two_stream_flux_tl',1,hook_handle)

end subroutine calc_two_stream_flux_tl


! ===== FILE: tcrad_set_two_stream_scheme_tl.F90 =====
subroutine set_two_stream_scheme_tl(i_scheme)

  integer(jpim), intent(in) :: i_scheme

  ! No tangent-linear variables: this routine only sets module parameters
  ! based on an integer switch (no dependence on atmospheric state).

  call set_two_stream_scheme(i_scheme)

end subroutine set_two_stream_scheme_tl




end module tcrad_tl
