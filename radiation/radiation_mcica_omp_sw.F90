! radiation_mcica_omp_sw.F90 - Monte-Carlo Independent Column Approximation shortwave solver
!
! (C) Copyright 2015- ECMWF.
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
! Modifications
!   2017-04-11  R. Hogan  Receive albedos at g-points
!   2017-04-22  R. Hogan  Store surface fluxes at all g-points
!   2017-10-23  R. Hogan  Renamed single-character variables

#include "ecrad_config.h"

module radiation_mcica_omp_sw

  use radiation_io, only : nulout
  public
  private :: solver_mcica_omp_sw_impl

contains

  ! Provides elemental function "delta_eddington"
#include "radiation_delta_eddington.h"

  !---------------------------------------------------------------------
  ! Shortwave Monte Carlo Independent Column Approximation
  ! (McICA). This implementation performs a clear-sky and a cloudy-sky
  ! calculation, and then weights the two to get the all-sky fluxes
  ! according to the total cloud cover. This method reduces noise for
  ! low cloud cover situations, and exploits the clear-sky
  ! calculations that are usually performed for diagnostic purposes
  ! simultaneously. The cloud generator has been carefully written
  ! such that the stochastic cloud field satisfies the prescribed
  ! overlap parameter accounting for this weighting.
  subroutine solver_mcica_omp_sw(nlev,istartcol,iendcol, &
       &  config, single_level, cloud, &
       &  od, ssa, g, od_cloud, ssa_cloud, g_cloud, &
       &  albedo_direct, albedo_diffuse, incoming_sw, &
       &  flux)

    use parkind1, only           : jprb
    use yomhook,  only           : lhook, dr_hook, jphook
    use radiation_io,   only           : nulerr, radiation_abort
    use radiation_config, only         : config_type
    use radiation_single_level, only   : single_level_type
    use radiation_cloud, only          : cloud_type
    use radiation_flux, only           : flux_type
#ifdef HAVE_ROCTX
    use roctx_profiling, only: roctxstartrange, roctxendrange
    use iso_c_binding, only: c_null_char
#endif

    implicit none

    integer, intent(in) :: nlev
    integer, intent(in) :: istartcol, iendcol
    type(config_type),        intent(in) :: config
    type(single_level_type),  intent(in) :: single_level
    type(cloud_type),         intent(in) :: cloud
    real(jprb), intent(in), dimension(config%n_g_sw, nlev, istartcol:iendcol) :: &
         &  od, ssa, g
    real(jprb), intent(in), dimension(config%n_bands_sw,nlev,istartcol:iendcol)   :: &
         &  od_cloud, ssa_cloud, g_cloud
    real(jprb), intent(in), dimension(config%n_g_sw,istartcol:iendcol) :: &
         &  albedo_direct, albedo_diffuse, incoming_sw
    type(flux_type), intent(inout):: flux

    real(jphook) :: hook_handle

#ifdef HAVE_ROCTX
    call roctxStartRange("radiation::mcica_omp_sw"//c_null_char)
#endif
    if (lhook) call dr_hook('radiation_mcica_omp_sw:solver_mcica_omp_sw',0,hook_handle)

    if (.not. config%do_clear) then
      write(nulerr,'(a)') '*** Error: shortwave McICA OMP requires clear-sky calculation to be performed'
      call radiation_abort()
    end if

    if (allocated(flux%sw_dn_direct) .and. allocated(flux%sw_dn_direct_clear)) then
      call solver_mcica_omp_sw_impl(nlev, size(cloud%fraction,1), istartcol, iendcol, config%n_g_sw, &
           & config%pdf_sampler%ncdf, config%pdf_sampler%nfsd, config%n_bands_sw, &
           & config%use_beta_overlap, &
           & config%do_sw_delta_scaling_with_gases, &
           & config%cloud_fraction_threshold, config%cloud_inhom_decorr_scaling, &
           & config%pdf_sampler%fsd1, config%pdf_sampler%inv_fsd_interval, &
           & config%i_band_from_reordered_g_sw, config%pdf_sampler%val, &
           & single_level%iseed, single_level%cos_sza, &
           & cloud%fraction, cloud%fractional_std, cloud%overlap_param, &
           & od, ssa, g, od_cloud, ssa_cloud, g_cloud, &
           & albedo_direct, albedo_diffuse, incoming_sw, &
           & flux%cloud_cover_sw, flux%sw_dn_diffuse_surf_clear_g, flux%sw_dn_direct_surf_clear_g, &
           & flux%sw_dn_diffuse_surf_g, flux%sw_dn_direct_surf_g, &
           & flux%sw_up_clear, flux%sw_dn_clear, flux%sw_up, flux%sw_dn, &
           & flux%sw_dn_direct, flux%sw_dn_direct_clear)
    else
      call solver_mcica_omp_sw_impl(nlev, size(cloud%fraction,1), istartcol, iendcol, config%n_g_sw, &
           & config%pdf_sampler%ncdf, config%pdf_sampler%nfsd, config%n_bands_sw, &
           & config%use_beta_overlap, &
           & config%do_sw_delta_scaling_with_gases, &
           & config%cloud_fraction_threshold, config%cloud_inhom_decorr_scaling, &
           & config%pdf_sampler%fsd1, config%pdf_sampler%inv_fsd_interval, &
           & config%i_band_from_reordered_g_sw, config%pdf_sampler%val, &
           & single_level%iseed, single_level%cos_sza, &
           & cloud%fraction, cloud%fractional_std, cloud%overlap_param, &
           & od, ssa, g, od_cloud, ssa_cloud, g_cloud, &
           & albedo_direct, albedo_diffuse, incoming_sw, &
           & flux%cloud_cover_sw, flux%sw_dn_diffuse_surf_clear_g, flux%sw_dn_direct_surf_clear_g, &
           & flux%sw_dn_diffuse_surf_g, flux%sw_dn_direct_surf_g, &
           & flux%sw_up_clear, flux%sw_dn_clear, flux%sw_up, flux%sw_dn)
    end if

#ifdef HAVE_ROCTX
    call roctxEndRange
#endif
    if (lhook) call dr_hook('radiation_mcica_omp_sw:solver_mcica_omp_sw',1,hook_handle)

  end subroutine solver_mcica_omp_sw

  subroutine solver_mcica_omp_sw_impl(nlev, ncol, istartcol, iendcol, ng, ncdf, nfsd, n_bands, &
       & use_beta_overlap, do_sw_delta_scaling_with_gases, &
       & cloud_fraction_threshold, cloud_inhom_decorr_scaling, fsd1, inv_fsd_interval, &
       & i_band_from_reordered_g_sw, pdf_val, iseed, cos_sza_col, &
       & cloud_fraction, cloud_fractional_std, cloud_overlap_param, &
       & od, ssa, g, od_cloud, ssa_cloud, g_cloud, &
       & albedo_direct, albedo_diffuse, incoming_sw, &
       & cloud_cover_sw, sw_dn_diffuse_surf_clear_g, sw_dn_direct_surf_clear_g, &
       & sw_dn_diffuse_surf_g, sw_dn_direct_surf_g, &
       & sw_up_clear, sw_dn_clear, sw_up, sw_dn, sw_dn_direct, sw_dn_direct_clear)

    use parkind1, only           : jprb
    use radiation_two_stream, only     : calc_two_stream_gammas_sw_single_band_omp, &
         &                               calc_reflectance_transmittance_sw_single_band_omp, &
         &                               calc_ref_trans_sw_omp, calc_ref_trans_sw_single_level_omp
    use radiation_adding_ica_sw, only  : adding_ica_sw_omp
    use radiation_cloud_generator_acc, only: cloud_generator_omp
    use radiation_cloud_cover, only    : beta2alpha, MaxCloudFrac

    implicit none

    integer, intent(in) :: nlev, ncol, istartcol, iendcol, ng, ncdf, nfsd, n_bands
    logical, intent(in) :: use_beta_overlap, do_sw_delta_scaling_with_gases
    real(jprb), intent(in) :: cloud_fraction_threshold, cloud_inhom_decorr_scaling
    real(jprb), intent(in) :: fsd1, inv_fsd_interval
    integer, intent(in) :: i_band_from_reordered_g_sw(ng)
    real(jprb), intent(in) :: pdf_val(ncdf,nfsd)
    integer, intent(in) :: iseed(ncol)
    real(jprb), intent(in) :: cos_sza_col(ncol)
    real(jprb), intent(in) :: cloud_fraction(ncol,nlev)
    real(jprb), intent(in) :: cloud_fractional_std(ncol,nlev)
    real(jprb), intent(in) :: cloud_overlap_param(ncol,nlev-1)
    real(jprb), intent(in) :: od(ng,nlev,istartcol:iendcol)
    real(jprb), intent(in) :: ssa(ng,nlev,istartcol:iendcol)
    real(jprb), intent(in) :: g(ng,nlev,istartcol:iendcol)
    real(jprb), intent(in) :: od_cloud(n_bands,nlev,istartcol:iendcol)
    real(jprb), intent(in) :: ssa_cloud(n_bands,nlev,istartcol:iendcol)
    real(jprb), intent(in) :: g_cloud(n_bands,nlev,istartcol:iendcol)
    real(jprb), intent(in) :: albedo_direct(ng,istartcol:iendcol)
    real(jprb), intent(in) :: albedo_diffuse(ng,istartcol:iendcol)
    real(jprb), intent(in) :: incoming_sw(ng,istartcol:iendcol)
    real(jprb), intent(inout) :: cloud_cover_sw(ncol)
    real(jprb), intent(inout) :: sw_dn_diffuse_surf_clear_g(ng,ncol)
    real(jprb), intent(inout) :: sw_dn_direct_surf_clear_g(ng,ncol)
    real(jprb), intent(inout) :: sw_dn_diffuse_surf_g(ng,ncol)
    real(jprb), intent(inout) :: sw_dn_direct_surf_g(ng,ncol)
    real(jprb), intent(inout) :: sw_up_clear(ncol,nlev+1)
    real(jprb), intent(inout) :: sw_dn_clear(ncol,nlev+1)
    real(jprb), intent(inout) :: sw_up(ncol,nlev+1)
    real(jprb), intent(inout) :: sw_dn(ncol,nlev+1)
    real(jprb), intent(inout), optional :: sw_dn_direct(ncol,nlev+1)
    real(jprb), intent(inout), optional :: sw_dn_direct_clear(ncol,nlev+1)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Local variables : Mapped into Global Memory
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! Fluxes per g point
    real(jprb), dimension(ng, nlev+1, istartcol:iendcol) :: flux_up, flux_dn_diffuse, flux_dn_direct
    real(jprb), dimension(ng, nlev+1, istartcol:iendcol) :: flux_up_clear, flux_dn_diffuse_clear, flux_dn_direct_clear

    ! workaround that allows inling of cloud generator
    real(jprb), dimension(ncdf, nfsd) :: sample_val

    ! copies to increase performance
    real(jprb), dimension(nlev, istartcol:iendcol) :: frac, frac_std
    real(jprb), dimension(nlev-1, istartcol:iendcol) :: overlap_param

    real(jprb), dimension(nlev, istartcol:iendcol) :: cum_cloud_cover
    real(jprb), dimension(nlev-1, istartcol:iendcol) :: pair_cloud_cover

    ! Cumulative product needed in computation of total_cloud_cover
    real(jprb) :: cum_product(istartcol:iendcol)

    ! First and last cloudy layers
    integer :: ibegin(istartcol:iendcol), iend(istartcol:iendcol)
 
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Local variables : Was stack and is now in Global Memory
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! Transmittance for the direct beam in clear and all skies
    real(jprb), dimension(ng, nlev, istartcol:iendcol) :: ref_clear,  trans_clear, ref_dir_clear
    !real(jprb), dimension(ng, nlev, istartcol:iendcol) :: trans_dir_diff_clear, trans_dir_dir_clear

    ! Fraction of direct beam scattered by a layer into the upwelling
    ! or downwelling diffuse streams, in clear and all skies
    real(jprb), dimension(ng, nlev, istartcol:iendcol) :: ref_dir, trans_dir_diff

    ! Transmittance for the direct beam in clear and all skies
    real(jprb), dimension(ng, nlev, istartcol:iendcol) :: trans_dir_dir

    ! Diffuse reflectance and transmittance for each layer in clear
    ! and all skies
    real(jprb), dimension(ng, nlev, istartcol:iendcol) :: reflectance, transmittance

    ! Optical depth scaling from the cloud generator, zero indicating
    ! clear skies
    real(jprb), dimension(ng,nlev,istartcol:iendcol) :: od_scaling

    ! Temporary working array
    real(jprb), dimension(ng,nlev+1,istartcol:iendcol) :: tmp_work_source

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Local variables : Stack
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    ! Cosine of solar zenith angle
    real(jprb)                                 :: cos_sza

    ! Combined gas+aerosol+cloud optical depth, single scattering
    ! albedo and asymmetry factor
    real(jprb) :: od_total, ssa_total, g_total

    ! Combined scattering optical depth
    real(jprb) :: scat_od

    ! Two-stream coefficients
    real(jprb) :: gamma1, gamma2, gamma3

    ! Modified optical depth after McICA scaling to represent cloud
    ! inhomogeneity
    real(jprb) :: od_cloud_new

    ! "Alpha" overlap parameter
    real(jprb) :: overlap_alpha

    ! Auxiliary for more efficient summation
    real(jprb) :: sum_up, sum_dn_direct, sum_dn_diffuse

    ! Loop indices for level, column and g point
    integer :: jlev, jcol, jg
    !real(jprb)  totalMem

!$OMP TARGET ENTER DATA MAP(ALLOC: ref_clear,  trans_clear, ref_dir_clear, &
    !$OMP&   od_scaling, tmp_work_source,&
    !$OMP&   ref_dir, trans_dir_diff, trans_dir_dir, reflectance, transmittance)
    !$OMP TARGET ENTER DATA MAP(ALLOC: flux_up, flux_dn_diffuse, flux_dn_direct, &
    !$OMP             flux_up_clear, flux_dn_diffuse_clear, flux_dn_direct_clear, &
    !$OMP             sample_val, frac, frac_std, overlap_param, &
    !$OMP             cum_cloud_cover, pair_cloud_cover, cum_product, ibegin, iend)
#if defined(__amdflang__)
    !$OMP TARGET DATA MAP(PRESENT, ALLOC: od, ssa, g, od_cloud, ssa_cloud, g_cloud, &
    !$OMP             albedo_direct, albedo_diffuse, incoming_sw, pdf_val, iseed, cos_sza_col, &
    !$OMP             cloud_fraction, cloud_fractional_std, cloud_overlap_param, &
    !$OMP             i_band_from_reordered_g_sw, cloud_cover_sw, &
    !$OMP             sw_dn_diffuse_surf_clear_g, sw_dn_direct_surf_clear_g, &
    !$OMP             sw_dn_diffuse_surf_g, sw_dn_direct_surf_g, &
    !$OMP             sw_up_clear, sw_dn_clear, sw_up, sw_dn)
#endif


    !totalMem = 9*ng * nlev
    !totalMem = totalMem+7*(ng)*(nlev+1) !flux_up,dn,up_clear,dn_clear,source
    !totalmem = totalMem*(iendcol-istartcol)*SIZEOF((real(jprb)))/1.e9
    !write(nulout,'(a,a,i0,a,g0.5)') __FILE__, " : LINE = ", __LINE__, " total_memory=",totalMem

    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO COLLAPSE(2)
    do jlev = 1,nfsd
      do jcol = 1,ncdf
        sample_val(jcol,jlev) = pdf_val(jcol,jlev)
      end do
    end do
    !$OMP END TARGET TEAMS DISTRIBUTE PARALLEL DO

    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO COLLAPSE(2)
    do jcol = istartcol,iendcol
      do jlev = 1, nlev
        frac(jlev, jcol) = cloud_fraction(jcol,jlev)
        frac_std(jlev, jcol) = cloud_fractional_std(jcol,jlev)
      end do
    end do
    !$OMP END TARGET TEAMS DISTRIBUTE PARALLEL DO

    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO COLLAPSE(2)
    do jcol = istartcol,iendcol
      do jlev = 1, nlev-1
        overlap_param(jlev, jcol) = cloud_overlap_param(jcol,jlev)
      end do
    end do
    !$OMP END TARGET TEAMS DISTRIBUTE PARALLEL DO

    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO COLLAPSE(2) PRIVATE(overlap_alpha)
    do jcol = istartcol,iendcol
      !Only perform calculation if sun above the horizon
      !---------------------------------------------------------------------
      ! manual inline from cum_cloud_cover_exp_ran >>>>>>>>>>>>>>>>>>>>>>>>
      ! Loop to compute total cloud cover and the cumulative cloud cover
      ! down to the base of each layer
      do jlev = 1,nlev-1
        if (cos_sza_col(jcol) > 0.0_jprb ) then
          ! Convert to "alpha" overlap parameter if necessary
          if (use_beta_overlap) then
            overlap_alpha = beta2alpha(overlap_param(jlev,jcol), &
                  &                     frac(jlev,jcol), frac(jlev+1,jcol))
          else
            overlap_alpha = overlap_param(jlev,jcol)
          end if
          ! Compute the combined cloud cover of layers jlev and jlev+1
          pair_cloud_cover(jlev, jcol) = overlap_alpha*max(frac(jlev,jcol),frac(jlev+1,jcol)) &
                &  + (1.0_jprb - overlap_alpha) &
                &  * (frac(jlev,jcol)+frac(jlev+1,jcol)-frac(jlev,jcol)*frac(jlev+1,jcol))
        end if
      end do
    end do
    !$OMP END TARGET TEAMS DISTRIBUTE PARALLEL DO

    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO
    do jcol = istartcol,iendcol
      !Only perform calculation if sun above the horizon
      if (cos_sza_col(jcol) > 0.0_jprb ) then
        cum_cloud_cover(1, jcol) = frac(1,jcol)
        cum_product(jcol) = 1.0_jprb - frac(1,jcol)
        do jlev = 1,nlev-1
          if (frac(jlev,jcol) >= MaxCloudFrac) then
            ! Cloud cover has reached one
            cum_product(jcol) = 0.0_jprb
          else
            cum_product(jcol) = cum_product(jcol) * (1.0_jprb - pair_cloud_cover(jlev, jcol)) &
                  &  / (1.0_jprb - frac(jlev,jcol))
          end if
          cum_cloud_cover(jlev+1, jcol) = 1.0_jprb - cum_product(jcol)
        end do
        cloud_cover_sw(jcol) = cum_cloud_cover(nlev,jcol);
        if (cloud_cover_sw(jcol) < cloud_fraction_threshold) then
          ! Treat column as clear sky: calling function therefore will not
          ! use od_scaling so we don't need to calculate it
          cloud_cover_sw(jcol) = 0.0_jprb
        end if
      end if
    end do
    !$OMP END TARGET TEAMS DISTRIBUTE PARALLEL DO

    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO
    do jcol = istartcol,iendcol
      !Only perform calculation if sun above the horizon
      if (cos_sza_col(jcol) > 0.0_jprb .and. cloud_cover_sw(jcol) >= cloud_fraction_threshold) then
        ! Cloud is present: need to calculate od_scaling
        ! Find range of cloudy layers
        ibegin(jcol) = nlev
        do jlev = 1, nlev
          if( frac(jlev,jcol) > 0.0_jprb ) then
            ibegin(jcol) = min(jlev, ibegin(jcol))
          end if
        end do

        iend(jcol) = ibegin(jcol)
        do jlev = ibegin(jcol)+1,nlev
          if (frac(jlev,jcol) > 0.0_jprb) then
            iend(jcol) = max(jlev, iend(jcol))
          end if
        end do
      end if
    end do
    !$OMP END TARGET TEAMS DISTRIBUTE PARALLEL DO

    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO COLLAPSE(2) PRIVATE(jcol, jg) FIRSTPRIVATE(istartcol, iendcol, ng, nlev) &
    !$OMP& THREAD_LIMIT(512)
    do jcol = istartcol,iendcol
       do jg = 1, ng
          ! Do cloudy-sky calculation
          call cloud_generator_omp(jg, ng, nlev, &
               &  iseed(jcol) + 0, & ! Workaround for nvhpc-24.1
               &  cloud_fraction_threshold, &
               &  frac(:,jcol), overlap_param(:,jcol), &
               &  cloud_inhom_decorr_scaling, frac_std(:,jcol), &
               &  ncdf, nfsd, &
               &  fsd1, inv_fsd_interval, &
               &  sample_val, &
               &  od_scaling(:,:,jcol), cloud_cover_sw(jcol)+0.0_jprb, & ! Workaround for nvhpc-24.1
               &  ibegin(jcol), iend(jcol), &
               &  cum_cloud_cover=cum_cloud_cover(:,jcol), &
               &  pair_cloud_cover=pair_cloud_cover(:,jcol))
       enddo
    enddo

    ! Split the former combined kernel so the always-on clear-sky path is not
    ! compiled together with cloudy two-stream/adding.
    ! Revert: restore radiation_mcica_omp_sw.F90.pre_split_l404
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO COLLAPSE(2) PRIVATE(cos_sza, gamma1, gamma2, gamma3, od_total, &
    !$OMP& ssa_total, g_total, jcol, jg, jlev) FIRSTPRIVATE(istartcol, iendcol, ng, nlev) THREAD_LIMIT(1024)
    do jcol = istartcol,iendcol
       do jg = 1, ng
          ! Only perform calculation if sun above the horizon
          if (cos_sza_col(jcol) > 0.0_jprb) then
             cos_sza = cos_sza_col(jcol)

             ! Clear-sky calculation - first compute clear-sky reflectance,
             ! transmittance etc at each model level
             if (.not. do_sw_delta_scaling_with_gases) then
                ! Delta-Eddington scaling has already been performed to the
                ! aerosol part of od, ssa and g
                call calc_ref_trans_sw_omp(jg, ng, nlev, &
                     &  cos_sza, od(:,:,jcol), ssa(:,:,jcol), g(:,:,jcol), &
                     &  ref_clear(:,:,jcol), trans_clear(:,:,jcol), &
                     &  ref_dir_clear(:,:,jcol), trans_dir_diff(:,:,jcol), &
                     &  trans_dir_dir(:,:,jcol))
             else
                ! Apply delta-Eddington scaling to the aerosol-gas mixture
                do jlev = 1,nlev
                   od_total  =  od(jg,jlev,jcol)
                   ssa_total = ssa(jg,jlev,jcol)
                   g_total   =   g(jg,jlev,jcol)
                   call delta_eddington(od_total, ssa_total, g_total)
                   call calc_two_stream_gammas_sw_single_band_omp(jg, &
                        &  cos_sza, ssa_total, g_total, &
                        &  gamma1, gamma2, gamma3)

                   call calc_reflectance_transmittance_sw_single_band_omp(jg, ng, &
                        &  cos_sza, od_total, ssa_total, &
                        &  gamma1, gamma2, gamma3, &
                        &  ref_clear(:,jlev,jcol), trans_clear(:,jlev,jcol), &
                        &  ref_dir_clear(:,jlev,jcol), trans_dir_diff(:,jlev,jcol), &
                        &  trans_dir_dir(:,jlev,jcol) )
                end do
             end if

             ! Use adding method to compute fluxes
             call adding_ica_sw_omp(jg, ng, nlev, incoming_sw(:,jcol), &
                  &  albedo_diffuse(:,jcol), albedo_direct(:,jcol), cos_sza, &
                  &  ref_clear(:,:,jcol), trans_clear(:,:,jcol), ref_dir_clear(:,:,jcol), trans_dir_diff(:,:,jcol), &
                  &  trans_dir_dir(:,:,jcol), flux_up(:,:,jcol), flux_dn_diffuse(:,:,jcol), flux_dn_direct(:,:,jcol), &
                  &  flux_up(:,:,jcol), flux_dn_diffuse(:,:,jcol), source=tmp_work_source(:,:,jcol))

             ! save temporarily clear-sky broadband fluxes
             do jlev = 1,nlev+1
                flux_up_clear(jg,jlev,jcol) = flux_up(jg,jlev,jcol)
                flux_dn_direct_clear(jg,jlev,jcol) = flux_dn_direct(jg,jlev,jcol)
                flux_dn_diffuse_clear(jg,jlev,jcol) = flux_dn_diffuse(jg,jlev,jcol)
             end do
             
             ! Store spectral downwelling fluxes at surface
             sw_dn_diffuse_surf_clear_g(jg,jcol) = flux_dn_diffuse(jg,nlev+1,jcol)
             sw_dn_direct_surf_clear_g(jg,jcol)  = flux_dn_direct(jg,nlev+1,jcol)
          else
             sw_dn_diffuse_surf_g(jg,jcol) = 0.0_jprb
             sw_dn_direct_surf_g(jg,jcol)  = 0.0_jprb
             sw_dn_diffuse_surf_clear_g(jg,jcol) = 0.0_jprb
             sw_dn_direct_surf_clear_g(jg,jcol)  = 0.0_jprb
          end if ! Sun above horizon
       end do !loop over spectral bands
    end do ! Loop over columns
    !$OMP END TARGET TEAMS DISTRIBUTE PARALLEL DO

    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO COLLAPSE(2) PRIVATE(cos_sza, od_cloud_new, od_total, &
    !$OMP& ssa_total, g_total, scat_od, jcol, jg, jlev) FIRSTPRIVATE(istartcol, iendcol, ng, nlev) THREAD_LIMIT(1024)
    do jcol = istartcol,iendcol
       do jg = 1, ng
          if (cos_sza_col(jcol) > 0.0_jprb) then
             cos_sza = cos_sza_col(jcol)

             if (cloud_cover_sw(jcol) >= cloud_fraction_threshold) then
                ! Total-sky calculation
                do jlev = 1,nlev
                   ! Compute combined gas+aerosol+cloud optical properties
                   if (frac(jlev,jcol) >= cloud_fraction_threshold) then
                      od_cloud_new = od_scaling(jg,jlev,jcol) &
                           &  * od_cloud(i_band_from_reordered_g_sw(jg),jlev,jcol)
                      od_total  = od(jg,jlev,jcol) + od_cloud_new
                      ssa_total = 0.0_jprb
                      g_total   = 0.0_jprb
                         
                      ! In single precision we need to protect against the
                      ! case that od_total > 0.0 and ssa_total > 0.0 but
                      ! od_total*ssa_total == 0 due to underflow
                      if (od_total > 0.0_jprb) then
                         scat_od = ssa(jg,jlev,jcol)*od(jg,jlev,jcol) &
                              &     + ssa_cloud(i_band_from_reordered_g_sw(jg),jlev,jcol) &
                              &     *  od_cloud_new
                         ssa_total = scat_od / od_total
                         if (scat_od > 0.0_jprb) then
                            g_total = (g(jg,jlev,jcol)*ssa(jg,jlev,jcol)*od(jg,jlev,jcol) &
                                 &     +   g_cloud(i_band_from_reordered_g_sw(jg),jlev,jcol) &
                                 &     * ssa_cloud(i_band_from_reordered_g_sw(jg),jlev,jcol) &
                                 &     *  od_cloud_new) &
                                 &     / scat_od
                         end if
                      end if
                      ! Apply delta-Eddington scaling to the cloud-aerosol-gas
                      ! mixture
                      if (do_sw_delta_scaling_with_gases) then
                         call delta_eddington(od_total, ssa_total, g_total)
                      end if

                      ! Compute cloudy-sky reflectance, transmittance etc at
                      ! each model level
                      call calc_ref_trans_sw_single_level_omp(jg, ng, &
                           &  cos_sza, od_total, ssa_total, g_total, &
                           &  reflectance(:,jlev,jcol), transmittance(:,jlev,jcol), &
                           &  ref_dir(:,jlev,jcol), trans_dir_diff(:,jlev,jcol), &
                           &  trans_dir_dir(:,jlev,jcol))
                   else
                      ! Clear-sky layer: copy over clear-sky values
                      reflectance(jg,jlev,jcol) = ref_clear(jg,jlev,jcol)
                      transmittance(jg,jlev,jcol) = trans_clear(jg,jlev,jcol)
                      ref_dir(jg,jlev,jcol) = ref_dir_clear(jg,jlev,jcol)
                      !trans_dir_diff(jg,jlev,jcol) = trans_dir_diff_clear(jg,jlev,jcol)
                      !trans_dir_dir(jg,jlev,jcol) = trans_dir_dir_clear(jg,jlev,jcol)
                   end if
                end do

                ! Use adding method to compute fluxes for an overcast sky
                call adding_ica_sw_omp(jg, ng, nlev, incoming_sw(:,jcol), &
                     &  albedo_diffuse(:,jcol), albedo_direct(:,jcol), cos_sza, &
                     &  reflectance(:,:,jcol), transmittance(:,:,jcol), ref_dir(:,:,jcol), trans_dir_diff(:,:,jcol), &
                     &  trans_dir_dir(:,:,jcol), flux_up(:,:,jcol), flux_dn_diffuse(:,:,jcol), flux_dn_direct(:,:,jcol), &
                     &  flux_up(:,:,jcol), flux_dn_diffuse(:,:,jcol), source=tmp_work_source(:,:,jcol))
                
                ! Likewise for surface spectral fluxes
                sw_dn_diffuse_surf_g(jg,jcol) = flux_dn_diffuse(jg,nlev+1,jcol)
                sw_dn_direct_surf_g(jg,jcol)  = flux_dn_direct(jg,nlev+1,jcol)
                sw_dn_diffuse_surf_g(jg,jcol) = cloud_cover_sw(jcol) *sw_dn_diffuse_surf_g(jg,jcol) &
                     &     + (1.0_jprb - cloud_cover_sw(jcol))*sw_dn_diffuse_surf_clear_g(jg,jcol)
                sw_dn_direct_surf_g(jg,jcol) = cloud_cover_sw(jcol) *sw_dn_direct_surf_g(jg,jcol) &
                     &     + (1.0_jprb - cloud_cover_sw(jcol))*sw_dn_direct_surf_clear_g(jg,jcol)
                
             else
                ! No cloud in profile and clear-sky fluxes already
                ! calculated: copy them over
                sw_dn_diffuse_surf_g(jg,jcol) = sw_dn_diffuse_surf_clear_g(jg,jcol)
                sw_dn_direct_surf_g(jg,jcol)  = sw_dn_direct_surf_clear_g(jg,jcol)

             end if ! Cloud is present in profile
          else
             sw_dn_diffuse_surf_g(jg,jcol) = 0.0_jprb
             sw_dn_direct_surf_g(jg,jcol)  = 0.0_jprb
             sw_dn_diffuse_surf_clear_g(jg,jcol) = 0.0_jprb
             sw_dn_direct_surf_clear_g(jg,jcol)  = 0.0_jprb
          end if ! Sun above horizon
       end do !loop over spectral bands
    end do ! Loop over columns
    !$OMP END TARGET TEAMS DISTRIBUTE PARALLEL DO

    ! Loop through columns
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO COLLAPSE(2) PRIVATE(cos_sza, sum_dn_diffuse, sum_dn_direct, sum_up) THREAD_LIMIT(16)
    do jcol = istartcol,iendcol
       do jlev = 1, nlev+1

        ! Only perform calculation if sun above the horizon
        if (cos_sza_col(jcol) > 0.0_jprb) then

          ! Sum over g-points to compute and save clear-sky broadband
          ! fluxes
          sum_up = 0.0_jprb
          sum_dn_direct = 0.0_jprb
          sum_dn_diffuse = 0.0_jprb

          do jg = 1,ng
             sum_up = sum_up + flux_up_clear(jg,jlev,jcol)
             sum_dn_direct = sum_dn_direct + flux_dn_direct_clear(jg,jlev,jcol)
             sum_dn_diffuse = sum_dn_diffuse + flux_dn_diffuse_clear(jg,jlev,jcol)
          end do
          sw_up_clear(jcol,jlev) = sum_up
          if (present(sw_dn_direct_clear)) then
             sw_dn_direct_clear(jcol,jlev) = sum_dn_direct
          end if
          sw_dn_clear(jcol,jlev) = sum_dn_diffuse + sum_dn_direct

          if (cloud_cover_sw(jcol) >= cloud_fraction_threshold) then
             ! Store overcast broadband fluxes
             sum_up = 0.0_jprb
             sum_dn_direct = 0.0_jprb
             sum_dn_diffuse = 0.0_jprb
             do jg = 1,ng
                sum_up = sum_up + flux_up(jg,jlev,jcol)
                sum_dn_direct = sum_dn_direct + flux_dn_direct(jg,jlev,jcol)
                sum_dn_diffuse = sum_dn_diffuse + flux_dn_diffuse(jg,jlev,jcol)
             end do
             sw_up(jcol,jlev) = sum_up
             if (present(sw_dn_direct)) then
                sw_dn_direct(jcol,jlev) = sum_dn_direct
             end if
             sw_dn(jcol,jlev) = sum_dn_diffuse + sum_dn_direct
          end if
       end if ! Sun above horizon
      end do ! Loop over columns
    end do
    !$OMP END TARGET TEAMS DISTRIBUTE PARALLEL DO

    ! Loop through columns
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO COLLAPSE(2) PRIVATE(cos_sza)
    do jlev = 1, nlev+1
      do jcol = istartcol,iendcol
        ! Only perform calculation if sun above the horizon
        if (cos_sza_col(jcol) > 0.0_jprb) then
          cos_sza = cos_sza_col(jcol)

          if (cloud_cover_sw(jcol) >= cloud_fraction_threshold) then
            ! Cloudy flux profiles currently assume completely overcast
            ! skies; perform weighted average with clear-sky profile
            sw_up(jcol,jlev) =  cloud_cover_sw(jcol) *sw_up(jcol,jlev) &
                &  + (1.0_jprb - cloud_cover_sw(jcol))*sw_up_clear(jcol,jlev)
            sw_dn(jcol,jlev) =  cloud_cover_sw(jcol) *sw_dn(jcol,jlev) &
                &  + (1.0_jprb - cloud_cover_sw(jcol))*sw_dn_clear(jcol,jlev)
            if (present(sw_dn_direct)) then
              sw_dn_direct(jcol,jlev) = cloud_cover_sw(jcol) *sw_dn_direct(jcol,jlev) &
                  &  + (1.0_jprb - cloud_cover_sw(jcol))*sw_dn_direct_clear(jcol,jlev)
            end if

          else
            ! No cloud in profile and clear-sky fluxes already
            ! calculated: copy them over
            sw_up(jcol,jlev) = sw_up_clear(jcol,jlev)
            sw_dn(jcol,jlev) = sw_dn_clear(jcol,jlev)
            if (present(sw_dn_direct)) then
              sw_dn_direct(jcol,jlev) = sw_dn_direct_clear(jcol,jlev)
            end if

          end if ! Cloud is present in profile
        else
          ! Set fluxes to zero if sun is below the horizon
          sw_up(jcol,jlev) = 0.0_jprb
          sw_dn(jcol,jlev) = 0.0_jprb
          if (present(sw_dn_direct)) then
            sw_dn_direct(jcol,jlev) = 0.0_jprb
          end if
          sw_up_clear(jcol,jlev) = 0.0_jprb
          sw_dn_clear(jcol,jlev) = 0.0_jprb
          if (present(sw_dn_direct_clear)) then
            sw_dn_direct_clear(jcol,jlev) = 0.0_jprb
          end if
        end if ! Sun above horizon
      end do ! Loop over columns
    end do ! Loop over levels
    !$OMP END TARGET TEAMS DISTRIBUTE PARALLEL DO

    !$OMP TARGET EXIT DATA MAP(DELETE: flux_up, flux_dn_diffuse, flux_dn_direct, &
    !$OMP             flux_up_clear, flux_dn_diffuse_clear, flux_dn_direct_clear, &
    !$OMP             sample_val, frac, frac_std, overlap_param, &
    !$OMP             cum_cloud_cover, pair_cloud_cover, cum_product, ibegin, iend)

    !$OMP TARGET EXIT DATA MAP(DELETE: ref_clear,  trans_clear, ref_dir_clear, &
    !$OMP&   od_scaling, tmp_work_source,&
    !$OMP&   ref_dir, trans_dir_diff, trans_dir_dir, reflectance, transmittance)

#if defined(__amdflang__)
    !$OMP END TARGET DATA
#endif

  end subroutine solver_mcica_omp_sw_impl


end module radiation_mcica_omp_sw
