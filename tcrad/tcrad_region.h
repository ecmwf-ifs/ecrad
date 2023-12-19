! tcrad_region.h - Properties of horizontal regions in the Tripleclouds assumption -*- f90 -*-
!
! (C) Copyright 2016- ECMWF.
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

#if NUM_REGIONS == 2

!---------------------------------------------------------------------
! Compute the region fractions and the (trivial) optical depth scaling
! for the case of only two regions (one clear and one cloudy)
subroutine calc_region_properties(nlev, &
     &  cloud_fraction, reg_fracs, od_scaling, cloud_fraction_threshold)

  use parkind1,     only : jprb
  use yomhook,      only : lhook, dr_hook, jphook

  ! Number of levels
  integer, intent(in) :: nlev
  
  ! Cloud fraction, i.e. the fraction of the gridbox assigned to all
  ! regions numbered 2 and above (region 1 is clear sky)
  real(jprb), intent(in), dimension(:)  :: cloud_fraction ! (nlev)

  ! Fractional area coverage of each region
  real(jprb), intent(out) :: reg_fracs(NREGION,nlev)

  ! Optical depth scaling for the cloudy regions
  real(jprb), intent(out) :: od_scaling(2:NREGION,nlev)

  ! Regions smaller than this are ignored
  real(jprb), intent(in), optional :: cloud_fraction_threshold
  
  ! In case the user doesn't supply cloud_fraction_threshold we use
  ! a default value
  real(jprb) :: frac_threshold
  
  real(jphook) :: hook_handle

  if (lhook) call dr_hook('tcrad:calc_region_properties',0,hook_handle)

  if (present(cloud_fraction_threshold)) then
    frac_threshold = cloud_fraction_threshold
  else
    frac_threshold = 1.0e-20_jprb
  end if
  ! Only one clear-sky and one cloudy region: cloudy region is
  ! homogeneous
  reg_fracs(2,1:nlev)  = cloud_fraction
  reg_fracs(1,1:nlev)  = 1.0_jprb - reg_fracs(2,1:nlev)
  od_scaling(2,1:nlev) = 1.0_jprb

end subroutine calc_region_properties


#elif NUM_REGIONS == 3

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
  use yomhook,      only : lhook, dr_hook, jphook

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
    frac_threshold = 1.0e-20_jprb
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
        reg_fracs(1,jlev)       = 1.0_jprb
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

#else

#error "calc_region_properties only defined for 2 or 3 regions"

#endif

subroutine calc_region_edge_areas(nlev, region_fracs, layer_thickness, &
     &                            inv_cloud_scale, region_edge_area)

  use parkind1,     only : jprb
  use yomhook,      only : lhook, dr_hook, jphook

  ! Number of levels
  integer, intent(in) :: nlev
  
  ! Fractional area coverage of each region
  real(jprb), intent(in) :: region_fracs(NREGION,nlev)

  ! Layer thickness in metres
  real(jprb), intent(in) :: layer_thickness(nlev)

  ! Inverse of the cloud horizontal separation scale, in m-1, using
  ! the definition of Fielding et al. (QJRMS 2020)
  real(jprb), intent(in) :: inv_cloud_scale(nlev)

  ! Area of the vertical interface between each pair of regions,
  ! divided by the horizontal area of the domain. For 3 regions there
  ! are two areas: between regions 1 and 2 and between regions 2 and 3
  ! (regions 1 and 3 are assumed not to touch).
  real(jprb), intent(out) :: region_edge_area(NREGION-1,nlev)

  real(jphook) :: hook_handle

  if (lhook) call dr_hook('tcrad:calc_region_edge_areas',0,hook_handle)

  ! Eq. 3 of Fielding et al. (2020), noting that 1-region_fracs(1,:)
  ! is the cloud fraction
  region_edge_area(1,:) = 4.0_jprb * layer_thickness * inv_cloud_scale &
       &  * sqrt(region_fracs(1,:)*(1.0_jprb-region_fracs(1,:)))

#if NUM_REGIONS == 3

  ! Apply the same formula but treating the cloud fraction as the
  ! fractional coverage of optically thick cloud only (with fraction
  ! region_fracs(3,:))
  region_edge_area(2,:) = 4.0_jprb * layer_thickness * inv_cloud_scale &
       &  * sqrt(region_fracs(3,:)*(1.0_jprb-region_fracs(3,:)))

#endif

  if (lhook) call dr_hook('tcrad:calc_region_edge_areas',1,hook_handle)

end subroutine calc_region_edge_areas
