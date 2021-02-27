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

  ! By making the number of regions a parameter NREGION, the compiler
  ! can optimize better.  We also provide a preprocessor parameter as
  ! in one or two cases the operations on 2x2 and 3x3 matrices
  ! (e.g. matrix exponentials) are completely different algorithms.
#define NUM_REGIONS 3
  integer(jpim), parameter :: NREGION = NUM_REGIONS

contains

! The following header files define routines and functions that make
! use of the NREGION parameter

#include "tcrad_region.h"

#include "tcrad_matrix.h"

#include "tcrad_overlap.h"

#include "tcrad_flux.h"

#include "tcrad_radiance.h"

#undef NUM_REGIONS

  subroutine calc_flux(nspec, nlev, surf_emission, surf_albedo, planck_hl, &
       &  cloud_fraction, fractional_std, od_gas, od_cloud, ssa_cloud, asymmetry_cloud, &
       &  overlap_param, flux_up, flux_dn, n_stream_per_hem, do_3d)
    
    use parkind1, only           : jpim, jprb
    use yomhook,  only           : lhook, dr_hook
    use tcrad_two_stream, only   : calc_reflectance_transmittance, &
         &  calc_radiance_source, LW_DIFFUSIVITY

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

    real(jprb) :: cloud_cover, weight

    real(jprb), dimension(3) :: mu_list, weight_list

    integer(jpim) :: n_stream_per_hem_local
    logical :: do_3d_local

    integer(jpim) :: jreg, jstream

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

    if (n_stream_per_hem_local > 0) then
      ! Fu et al. (1997) method: pass N beams through the
      ! atmosphere using the two-stream solution as the scattering
      ! source function

      if (n_stream_per_hem_local == 1) then
        weight_list(1) = 1;
        mu_list = 1.0_jprb / LW_DIFFUSIVITY
      else if (n_stream_per_hem_local == 2) then
        weight_list(1:2) = [0.5_jprb, 0.5_jprb];
        mu_list(1:2) = [0.5_jprb-0.5_jprb/sqrt(3.0_jprb), &
             &          0.5_jprb+0.5_jprb/sqrt(3.0_jprb)]
      else
        n_stream_per_hem_local = 3 ! Maximum
        weight_list(1:3) = [5.0_jprb/18.0_jprb, 4.0_jprb/9.0_jprb, 5.0_jprb/18.0_jprb]
        mu_list(1:3) = [0.5_jprb-0.5_jprb*sqrt(0.6_jprb), 0.5_jprb, &
             &          0.5_jprb+0.5_jprb*sqrt(0.6_jprb)]
      end if

      flux_up = 0.0_jprb
      flux_dn = 0.0_jprb
      do jstream = 1,n_stream_per_hem_local
        weight = weight_list(jstream)*mu_list(jstream) &
             &  / sum(weight_list(1:n_stream_per_hem_local) &
             &          * mu_list(1:n_stream_per_hem_local))
        call calc_radiance_source(nspec, nlev, NREGION, &
             &  mu_list(jstream), &
             &  region_fracs, planck_hl, od, ssa, asymmetry_cloud, &
             &  flux_up_base, flux_dn_base, flux_up_top, flux_dn_top, &
             &  transmittance, source_up, source_dn)
        call calc_multiregion_radiance_dn(nspec, nlev, &
             &  weight, &
             &  transmittance, source_dn, v_overlap, flux_dn)
        call calc_multiregion_radiance_up(nspec, nlev, &
             &  weight, flux_up_base(:,:,nlev), &
             &  transmittance, source_up, u_overlap, flux_up)
      end do

    else
      ! Simply take the existing two-stream fluxes
      flux_up(:,1:nlev) = sum(flux_up_top,2)
      flux_up(:,nlev+1) = sum(flux_up_base(:,:,nlev),2)
      flux_dn(:,1:nlev) = sum(flux_dn_top,2)
      flux_dn(:,nlev+1) = sum(flux_dn_base(:,:,nlev),2)
    end if

  end subroutine calc_flux

  subroutine calc_no_scattering_flux(nspec, nlev, surf_emission, surf_albedo, planck_hl, &
       &  cloud_fraction, fractional_std, od_gas, od_cloud, &
       &  overlap_param, flux_up, flux_dn, n_stream_per_hem, do_3d)
    
    use parkind1, only           : jpim, jprb
    use yomhook,  only           : lhook, dr_hook
    use tcrad_two_stream, only   : calc_reflectance_transmittance, &
         &  calc_radiance_source, calc_no_scattering_radiance_source, &
         &  LW_DIFFUSIVITY

    implicit none

    ! Inputs

    ! Number of spectral intervals and levels
    integer(jpim), intent(in) :: nspec, nlev


    real(jprb), intent(in), dimension(nspec) :: surf_emission, surf_albedo

    real(jprb), intent(in), dimension(nspec,nlev+1) :: planck_hl

    real(jprb), intent(in), dimension(nlev)  :: cloud_fraction, fractional_std

    real(jprb), intent(in), dimension(nspec,nlev) :: od_gas, od_cloud

    real(jprb), intent(in), dimension(nlev-1) :: overlap_param

    real(jprb), intent(out), dimension(nspec,nlev+1) :: flux_up, flux_dn

    integer,    intent(in), optional :: n_stream_per_hem
    logical,    intent(in), optional :: do_3d

    real(jprb), dimension(nspec,NREGION,nlev)   :: od

    real(jprb), dimension(nspec,NREGION,nlev) :: transmittance
    real(jprb), dimension(nspec,NREGION,nlev) :: source_up, source_dn
    
    logical :: is_cloud_free_layer(0:nlev+1)

    real(jprb), dimension(NREGION,NREGION,nlev+1) :: u_overlap, v_overlap
  
    real(jprb), dimension(nspec,NREGION) :: flux_up_surf

    real(jprb) :: od_scaling(2:NREGION,nlev)
    
    ! Fractional area coverage of each region
    real(jprb) :: region_fracs(1:NREGION,nlev)

    real(jprb), parameter :: cloud_fraction_threshold = 1.0e-6

    real(jprb) :: cloud_cover

    real(jprb), dimension(3) :: mu_list, weight_list

    integer(jpim) :: n_stream_per_hem_local
    logical :: do_3d_local

    integer(jpim) :: jreg, jstream

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

    ! Average gas and cloud properties
    od(:,1,:) = od_gas
    do jreg = 2,NREGION
      od(:,jreg,:) = od_gas + od_cloud*spread(od_scaling(jreg,:),1,nspec)
    end do

    is_cloud_free_layer(0) = .true.
    is_cloud_free_layer(1:nlev) = (region_fracs(1,:) == 1.0_jprb)
    is_cloud_free_layer(nlev+1) = .true.

    if (n_stream_per_hem_local == 0) then
      n_stream_per_hem_local = 1
    end if
    ! Fu et al. (1997) method: pass N beams through the
    ! atmosphere using the two-stream solution as the scattering
    ! source function

    if (n_stream_per_hem_local == 1) then
      weight_list(1) = 1;
      mu_list = 1.0_jprb / LW_DIFFUSIVITY
    else if (n_stream_per_hem_local == 2) then
      weight_list(1:2) = [0.5_jprb, 0.5_jprb];
      mu_list(1:2) = [0.5_jprb-0.5_jprb/sqrt(3.0_jprb), &
           &          0.5_jprb+0.5_jprb/sqrt(3.0_jprb)]
    else
      n_stream_per_hem_local = 3 ! Maximum
      weight_list(1:3) = [5.0_jprb/18.0_jprb, 4.0_jprb/9.0_jprb, 5.0_jprb/18.0_jprb]
      mu_list(1:3) = [0.5_jprb-0.5_jprb*sqrt(0.6_jprb), 0.5_jprb, &
           &          0.5_jprb+0.5_jprb*sqrt(0.6_jprb)]
    end if
    
    flux_up = 0.0_jprb
    flux_dn = 0.0_jprb
    flux_up_surf = spread(surf_emission,2,NREGION)*spread(region_fracs(:,nlev),1,nspec)
    do jstream = 1,n_stream_per_hem_local
      call calc_no_scattering_radiance_source(nspec, nlev, NREGION, &
           &  mu_list(jstream), &
           &  region_fracs, planck_hl, od,  &
           &  transmittance, source_up, source_dn)
      call calc_multiregion_radiance_dn(nspec, nlev, &
           &  weight_list(jstream), &
           &  transmittance, source_dn, v_overlap, flux_dn)
      call calc_multiregion_radiance_up(nspec, nlev, &
           &  weight_list(jstream), flux_up_surf, &
           &  transmittance, source_up, u_overlap, flux_up)
    end do

  end subroutine calc_no_scattering_flux

end module tcrad_3region_solver
