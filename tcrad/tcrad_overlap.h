! tcrad_overlap.h - Calculate cloud overlap quantities -*- f90 -*-
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
  if (NREGION == 2) then
    ! Clear in upper layer, cloudy in lower layer
    overlap_matrix(1,2) = pair_cloud_cover - cf_upper
    ! Clear in lower layer, cloudy in upper layer
    overlap_matrix(2,1) = pair_cloud_cover - cf_lower
    ! Cloudy in both layers
    overlap_matrix(2,2) = cf_upper + cf_lower - pair_cloud_cover
  else
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
  end if

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
     &     region_fracs, overlap_param, u_overlap, v_overlap, decorrelation_scaling, &
     &     cloud_fraction_threshold, cloud_cover)

  use parkind1,     only : jprb
  use yomhook,      only : lhook, dr_hook

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

  real(jprb) :: hook_handle

  if (lhook) call dr_hook('calc_overlap_',0,hook_handle)

  if (present(decorrelation_scaling)) then
    used_decorrelation_scaling = decorrelation_scaling
  else
    used_decorrelation_scaling = 1.0_jprb
  end if
  
  if (present(cloud_fraction_threshold)) then
    frac_threshold = cloud_fraction_threshold
  else
    frac_threshold = 1.0e-20_jprb
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
  
  if (lhook) call dr_hook('calc_overlap_matrices',1,hook_handle)
  
end subroutine calc_overlap_matrices
  
