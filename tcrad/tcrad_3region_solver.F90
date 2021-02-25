! tcrad_3region_solver.F90 - Solve for longwave fluxes or radiances with the Tripleclouds assumption
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

module tcrad_3region_solver

  use parkind1, only : jpim
  use tcrad_two_stream, only : calc_reflectance_transmittance

  integer(jpim), parameter :: NREGION = 3

contains

#include "matrix_functions.h"

#include "calc_overlap_matrices.h"

#include "calc_multiregion_flux.h"

#include "calc_multiregion_radiance_up.h"


  !---------------------------------------------------------------------
  ! Compute the optical depth scalings for the optically "thick" and
  ! "thin" regions of a Tripleclouds representation of a sub-grid PDF
  ! of cloud optical depth. Following Shonk and Hogan (2008), the 16th
  ! percentile is used for the thin region, and the formulas estimate
  ! this for both lognormal and gamma distributions. However, an
  ! adjustment is needed for the gamma distribution at large
  ! fractional standard deviations.
  subroutine calc_region_properties(nlev, do_gamma, &
       &  cloud_fraction, frac_std, reg_fracs, od_scaling, cloud_fraction_threshold)

    use parkind1,     only : jprb
    use yomhook,      only : lhook, dr_hook

    ! Minimum od_scaling in the case of a Gamma distribution
    real(jprb), parameter :: MinGammaODScaling = 0.025_jprb

    ! At large fractional standard deviations (FSDs), we cannot
    ! capture the behaviour of a gamma distribution with two equally
    ! weighted points; we need to weight the first ("lower") point
    ! more.  The weight of the first point is normally 0.5, but for
    ! FSDs between 1.5 and 3.725 it rises linearly to 0.9, and for
    ! higher FSD it is capped at this value.  The weight of the second
    ! point is one minus the first point.
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
    real(jprb), intent(out) :: reg_fracs(1:NREGION,nlev)

    ! Optical depth scaling for the cloudy regions
    real(jprb), intent(out) :: od_scaling(2:NREGION,nlev)

    ! Regions smaller than this are ignored
    real(jprb), intent(in), optional :: cloud_fraction_threshold

    ! In case the user doesn't supply cloud_fraction_threshold we use
    ! a default value
    real(jprb) :: frac_threshold

    ! Loop indices
    integer :: jlev

    real(jprb) :: hook_handle

    if (lhook) call dr_hook('tripleclouds_solver:calc_region_properties',0,hook_handle)

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
      ! If we treat the distribution as a gamma then the 16th
      ! percentile is close to the following.  Note that because it
      ! becomes vanishingly small for FSD >~ 2, we have a lower limit
      ! of 1/40, and for higher FSDs reduce the fractional cover of
      ! the denser region - see region_fractions routine below
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

    if (lhook) call dr_hook('tripleclouds_solver:calc_region_properties',1,hook_handle)

  end subroutine calc_region_properties

  subroutine calc_flux(nspec, nlev, surf_emission, surf_albedo, planck_hl, &
       &  cloud_fraction, fractional_std, od_gas, od_cloud, ssa_cloud, asymmetry_cloud, &
       &  overlap_param, flux_up, flux_dn, n_stream_per_hem, do_3d)
    
    use parkind1, only           : jpim, jprb
    use yomhook,  only           : lhook, dr_hook
    use tcrad_two_stream, only   : calc_reflectance_transmittance

    implicit none

    ! Inputs

    ! Number of spectral intervals and levels
    integer(jpim), intent(in) :: nspec, nlev


    real(jprb), intent(in), dimension(nspec) :: surf_emission, surf_albedo

    real(jprb), intent(in), dimension(nspec,nlev+1) :: planck_hl

    real(jprb), intent(in), dimension(nlev)  :: cloud_fraction, fractional_std

    real(jprb), intent(in), dimension(nspec,nlev) :: od_gas, od_cloud, ssa_cloud, asymmetry_cloud

    real(jprb), intent(in), dimension(nlev-1) :: overlap_param

    real(jprb), intent(out), dimension(nspec,nlev+1) :: flux_up, flux_dn

    integer,    intent(in), optional :: n_stream_per_hem
    logical,    intent(in), optional :: do_3d

    real(jprb), dimension(nspec,NREGION,nlev)   :: od
    real(jprb), dimension(nspec,2:NREGION,nlev) :: ssa

    real(jprb), dimension(nspec,NREGION,nlev) :: reflectance, transmittance
    real(jprb), dimension(nspec,NREGION,nlev) :: source_up, source_dn
    
    logical :: is_cloud_free_layer(0:nlev+1)

    real(jprb), dimension(NREGION,NREGION,nlev+1) :: u_overlap, v_overlap
  
    real(jprb), dimension(nspec,NREGION,nlev) :: flux_up_base, flux_dn_base
    real(jprb), dimension(nspec,NREGION,nlev) :: flux_up_top, flux_dn_top

    real(jprb) :: od_scaling(2:NREGION,nlev)
    
    ! Fractional area coverage of each region
    real(jprb) :: region_fracs(1:NREGION,nlev)

    real(jprb), parameter :: cloud_fraction_threshold = 1.0e-6

    real(jprb) :: cloud_cover

    integer(jpim) :: n_stream_per_hem_local
    logical :: do_3d_local

    integer(jpim) :: jreg

    if (present(n_stream_per_hem)) then
      n_stream_per_hem_local = n_stream_per_hem
    else
      n_stream_per_hem_local = 1
    end if

    if (present(do_3d)) then
      do_3d_local = do_3d
    else
      do_3d_local = .false.
    end if


    ! Compute the wavelength-independent region fractions and
    ! optical-depth scalings
    call calc_region_properties(nlev, &
         &  .true., cloud_fraction, fractional_std, region_fracs, &
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
    od(:,1,:) = od_gas
    do jreg = 2,NREGION
      od(:,jreg,:) = od_gas + od_cloud*spread(od_scaling(jreg,:),1,nspec)
      ssa(:,jreg,:) = ssa_cloud(:,:)*od_cloud(:,:) &
           &  * spread(od_scaling(jreg,:),1,nspec) / od(:,jreg,:)
    end do

    is_cloud_free_layer(0) = .true.
    is_cloud_free_layer(1:nlev) = (region_fracs(1,:) == 1.0_jprb)
    is_cloud_free_layer(nlev+1) = .true.

    call calc_reflectance_transmittance(nspec, nlev, NREGION, &
         &  region_fracs, planck_hl, od, ssa, asymmetry_cloud, &
         &  reflectance, transmittance, source_up, source_dn)

    call calc_multiregion_flux(nspec, nlev, surf_emission, surf_albedo, &
         &  reflectance, transmittance, source_up, source_dn, &
         &  is_cloud_free_layer, u_overlap, v_overlap, &
         &  flux_up_base, flux_dn_base, flux_up_top, flux_dn_top)

    if (n_stream_per_hem_local > 1) then
      ! Fu et al. (1997) method: pass four beams through the
      ! atmosphere using the two-stream solution as the scattering
      ! source function

    else
      ! Simply take the existing two-stream fluxes
      flux_up(:,1:nlev) = sum(flux_up_top,2)
      flux_up(:,nlev+1) = sum(flux_up_base(:,:,nlev),2)
      flux_dn(:,1:nlev) = sum(flux_dn_top,2)
      flux_dn(:,nlev+1) = sum(flux_dn_base(:,:,nlev),2)
    end if

  end subroutine calc_flux

end module tcrad_3region_solver
