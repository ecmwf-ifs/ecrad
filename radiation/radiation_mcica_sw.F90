! radiation_mcica_sw.F90 - Monte-Carlo Independent Column Approximation shortwave solver
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

module radiation_mcica_sw

  public

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
  subroutine solver_mcica_sw(nlev,istartcol,iendcol, &
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
    use radiation_two_stream, only     : calc_ref_trans_sw
    use radiation_adding_ica_sw, only  : adding_ica_sw
    use radiation_cloud_generator, only: cloud_generator

    implicit none

    ! Inputs
    integer, intent(in) :: nlev               ! number of model levels
    integer, intent(in) :: istartcol, iendcol ! range of columns to process
    type(config_type),        intent(in) :: config
    type(single_level_type),  intent(in) :: single_level
    type(cloud_type),         intent(in) :: cloud

    ! Gas and aerosol optical depth, single-scattering albedo and
    ! asymmetry factor at each shortwave g-point
    real(jprb), intent(in), dimension(config%n_g_sw, nlev, istartcol:iendcol) :: &
         &  od, ssa, g

    ! Cloud and precipitation optical depth, single-scattering albedo and
    ! asymmetry factor in each shortwave band
    real(jprb), intent(in), dimension(config%n_bands_sw,nlev,istartcol:iendcol)   :: &
         &  od_cloud, ssa_cloud, g_cloud

    ! Direct and diffuse surface albedos, and the incoming shortwave
    ! flux into a plane perpendicular to the incoming radiation at
    ! top-of-atmosphere in each of the shortwave g points
    real(jprb), intent(in), dimension(config%n_g_sw,istartcol:iendcol) :: &
         &  albedo_direct, albedo_diffuse, incoming_sw

    ! Output
    type(flux_type), intent(inout):: flux

    ! Local variables

    ! Cosine of solar zenith angle
    real(jprb)                                 :: cos_sza

    ! Diffuse reflectance and transmittance for each layer in clear
    ! and all skies
    real(jprb), dimension(config%n_g_sw, nlev) :: ref_clear, trans_clear, reflectance, transmittance

    ! Fraction of direct beam scattered by a layer into the upwelling
    ! or downwelling diffuse streams, in clear and all skies
    real(jprb), dimension(config%n_g_sw, nlev) :: ref_dir_clear, trans_dir_diff_clear, ref_dir, trans_dir_diff

    ! Transmittance for the direct beam in clear and all skies
    real(jprb), dimension(config%n_g_sw, nlev) :: trans_dir_dir_clear, trans_dir_dir

    ! Fluxes per g point
    real(jprb), dimension(config%n_g_sw, nlev+1) :: flux_up, flux_dn_diffuse, flux_dn_direct

    ! Combined gas+aerosol+cloud optical depth, single scattering
    ! albedo and asymmetry factor
    real(jprb), dimension(config%n_g_sw) :: od_total, ssa_total, g_total

    ! Combined scattering optical depth
    real(jprb) :: scat_od

    ! Optical depth scaling from the cloud generator, zero indicating
    ! clear skies
    real(jprb), dimension(config%n_g_sw,nlev) :: od_scaling

    ! Modified optical depth after McICA scaling to represent cloud
    ! inhomogeneity
    real(jprb), dimension(config%n_g_sw) :: od_cloud_new

    ! Total cloud cover output from the cloud generator
    real(jprb) :: total_cloud_cover

    ! Temporary storage for more efficient summation
#ifdef DWD_REDUCTION_OPTIMIZATIONS
    real(jprb), dimension(nlev+1,3) :: sum_aux
#else
    real(jprb) :: sum_up, sum_dn_diff, sum_dn_dir
#endif

    ! Number of g points
    integer :: ng

    ! Loop indices for level, column and g point
    integer :: jlev, jcol, jg

    real(jphook) :: hook_handle

    if (lhook) call dr_hook('radiation_mcica_sw:solver_mcica_sw',0,hook_handle)

    if (.not. config%do_clear) then
      write(nulerr,'(a)') '*** Error: shortwave McICA requires clear-sky calculation to be performed'
      call radiation_abort()      
    end if

    ng = config%n_g_sw

    ! Loop through columns
    do jcol = istartcol,iendcol
      ! Only perform calculation if sun above the horizon
      if (single_level%cos_sza(jcol) > 0.0_jprb) then
        cos_sza = single_level%cos_sza(jcol)

        ! Clear-sky calculation - first compute clear-sky reflectance,
        ! transmittance etc at each model level
        if (.not. config%do_sw_delta_scaling_with_gases) then
          ! Delta-Eddington scaling has already been performed to the
          ! aerosol part of od, ssa and g
          call calc_ref_trans_sw(ng*nlev, &
               &  cos_sza, od(:,:,jcol), ssa(:,:,jcol), g(:,:,jcol), &
               &  ref_clear, trans_clear, &
               &  ref_dir_clear, trans_dir_diff_clear, &
               &  trans_dir_dir_clear)
        else
          ! Apply delta-Eddington scaling to the aerosol-gas mixture
          do jlev = 1,nlev
            od_total  =  od(:,jlev,jcol)
            ssa_total = ssa(:,jlev,jcol)
            g_total   =   g(:,jlev,jcol)
            call delta_eddington(od_total, ssa_total, g_total)
            call calc_ref_trans_sw(ng, &
                 &  cos_sza, od_total, ssa_total, g_total, &
                 &  ref_clear(:,jlev), trans_clear(:,jlev), &
                 &  ref_dir_clear(:,jlev), trans_dir_diff_clear(:,jlev), &
                 &  trans_dir_dir_clear(:,jlev) )
          end do
        end if

        ! Use adding method to compute fluxes
        call adding_ica_sw(ng, nlev, incoming_sw(:,jcol), &
             &  albedo_diffuse(:,jcol), albedo_direct(:,jcol), spread(cos_sza,1,ng), &
             &  ref_clear, trans_clear, ref_dir_clear, trans_dir_diff_clear, &
             &  trans_dir_dir_clear, flux_up, flux_dn_diffuse, flux_dn_direct)
        
        ! Sum over g-points to compute and save clear-sky broadband
        ! fluxes. Note that the built-in "sum" function is very slow,
        ! and before being replaced by the alternatives below
        ! accounted for around 40% of the total cost of this routine.
#ifdef DWD_REDUCTION_OPTIMIZATIONS
        ! Optimized summation for the NEC architecture
        sum_aux(:,:) = 0.0_jprb
        do jg = 1,ng
          do jlev = 1,nlev+1
            sum_aux(jlev,1) = sum_aux(jlev,1) + flux_up(jg,jlev)
            sum_aux(jlev,2) = sum_aux(jlev,2) + flux_dn_direct(jg,jlev)
            sum_aux(jlev,3) = sum_aux(jlev,3) + flux_dn_diffuse(jg,jlev)
          end do
        end do
        flux%sw_up_clear(jcol,:) = sum_aux(:,1)
        flux%sw_dn_clear(jcol,:) = sum_aux(:,2) + sum_aux(:,3)
        if (allocated(flux%sw_dn_direct_clear)) then
          flux%sw_dn_direct_clear(jcol,:) = sum_aux(:,2)
        end if
#else
        ! Optimized summation for the x86-64 architecture
        do jlev = 1,nlev+1
          sum_up      = 0.0_jprb
          sum_dn_diff = 0.0_jprb
          sum_dn_dir  = 0.0_jprb
          !$omp simd reduction(+:sum_up, sum_dn_diff, sum_dn_dir)
          do jg = 1,ng
            sum_up      = sum_up      + flux_up(jg,jlev)
            sum_dn_diff = sum_dn_diff + flux_dn_diffuse(jg,jlev)
            sum_dn_dir  = sum_dn_dir  + flux_dn_direct(jg,jlev)
          end do
          flux%sw_up_clear(jcol,jlev) = sum_up
          flux%sw_dn_clear(jcol,jlev) = sum_dn_diff + sum_dn_dir
          if (allocated(flux%sw_dn_direct_clear)) then
            flux%sw_dn_direct_clear(jcol,jlev) = sum_dn_dir
          end if
        end do
#endif
        
        ! Store spectral downwelling fluxes at surface
        do jg = 1,ng
          flux%sw_dn_diffuse_surf_clear_g(jg,jcol) = flux_dn_diffuse(jg,nlev+1)
          flux%sw_dn_direct_surf_clear_g(jg,jcol)  = flux_dn_direct(jg,nlev+1)
        end do

        ! Do cloudy-sky calculation
        call cloud_generator(ng, nlev, config%i_overlap_scheme, &
             &  single_level%iseed(jcol), &
             &  config%cloud_fraction_threshold, &
             &  cloud%fraction(jcol,:), cloud%overlap_param(jcol,:), &
             &  config%cloud_inhom_decorr_scaling, cloud%fractional_std(jcol,:), &
             &  config%pdf_sampler, od_scaling, total_cloud_cover, &
             &  use_beta_overlap=config%use_beta_overlap, &
             &  use_vectorizable_generator=config%use_vectorizable_generator)

        ! Store total cloud cover
        flux%cloud_cover_sw(jcol) = total_cloud_cover
        
        if (total_cloud_cover >= config%cloud_fraction_threshold) then
          ! Total-sky calculation
          do jlev = 1,nlev
            ! Compute combined gas+aerosol+cloud optical properties
            if (cloud%fraction(jcol,jlev) >= config%cloud_fraction_threshold) then
              do jg = 1,ng
                od_cloud_new(jg) = od_scaling(jg,jlev) &
                   &  * od_cloud(config%i_band_from_reordered_g_sw(jg),jlev,jcol)
                od_total(jg)  = od(jg,jlev,jcol) + od_cloud_new(jg)
                ssa_total(jg) = 0.0_jprb
                g_total(jg)   = 0.0_jprb

                ! In single precision we need to protect against the
                ! case that od_total > 0.0 and ssa_total > 0.0 but
                ! od_total*ssa_total == 0 due to underflow
                if (od_total(jg) > 0.0_jprb) then
                  scat_od = ssa(jg,jlev,jcol)*od(jg,jlev,jcol) &
                       &     + ssa_cloud(config%i_band_from_reordered_g_sw(jg),jlev,jcol) &
                       &     *  od_cloud_new(jg)
                  ssa_total(jg) = scat_od / od_total(jg)
                  if (scat_od > 0.0_jprb) then
                    g_total(jg) = (g(jg,jlev,jcol)*ssa(jg,jlev,jcol)*od(jg,jlev,jcol) &
                         &     +   g_cloud(config%i_band_from_reordered_g_sw(jg),jlev,jcol) &
                         &     * ssa_cloud(config%i_band_from_reordered_g_sw(jg),jlev,jcol) &
                         &     *  od_cloud_new(jg)) &
                         &     / scat_od
                  end if
                end if
              end do

              ! Apply delta-Eddington scaling to the cloud-aerosol-gas
              ! mixture
              if (config%do_sw_delta_scaling_with_gases) then
                call delta_eddington(od_total, ssa_total, g_total)
              end if

              ! Compute cloudy-sky reflectance, transmittance etc at
              ! each model level
              call calc_ref_trans_sw(ng, &
                   &  cos_sza, od_total, ssa_total, g_total, &
                   &  reflectance(:,jlev), transmittance(:,jlev), &
                   &  ref_dir(:,jlev), trans_dir_diff(:,jlev), &
                   &  trans_dir_dir(:,jlev))
              
            else
              ! Clear-sky layer: copy over clear-sky values
              do jg = 1,ng
                reflectance(jg,jlev) = ref_clear(jg,jlev)
                transmittance(jg,jlev) = trans_clear(jg,jlev)
                ref_dir(jg,jlev) = ref_dir_clear(jg,jlev)
                trans_dir_diff(jg,jlev) = trans_dir_diff_clear(jg,jlev)
                trans_dir_dir(jg,jlev) = trans_dir_dir_clear(jg,jlev)
              end do
            end if
          end do
            
          ! Use adding method to compute fluxes for an overcast sky
          call adding_ica_sw(ng, nlev, incoming_sw(:,jcol), &
               &  albedo_diffuse(:,jcol), albedo_direct(:,jcol), spread(cos_sza,1,ng), &
               &  reflectance, transmittance, ref_dir, trans_dir_diff, &
               &  trans_dir_dir, flux_up, flux_dn_diffuse, flux_dn_direct)
          
          ! Store overcast broadband fluxes
#ifdef DWD_REDUCTION_OPTIMIZATIONS
          sum_aux(:,:) = 0.0_jprb
          do jg = 1,ng
            do jlev = 1,nlev+1
              sum_aux(jlev,1) = sum_aux(jlev,1) + flux_up(jg,jlev)
              sum_aux(jlev,2) = sum_aux(jlev,2) + flux_dn_direct(jg,jlev)
              sum_aux(jlev,3) = sum_aux(jlev,3) + flux_dn_diffuse(jg,jlev)
            end do
          end do
          flux%sw_up(jcol,:) = sum_aux(:,1)
          flux%sw_dn(jcol,:) = sum_aux(:,2) + sum_aux(:,3)
          if (allocated(flux%sw_dn_direct)) then
            flux%sw_dn_direct(jcol,:) = sum_aux(:,2)
          end if
#else
          do jlev = 1,nlev+1
            sum_up      = 0.0_jprb
            sum_dn_diff = 0.0_jprb
            sum_dn_dir  = 0.0_jprb
            !$omp simd reduction(+:sum_up, sum_dn_diff, sum_dn_dir)
            do jg = 1,ng
              sum_up      = sum_up      + flux_up(jg,jlev)
              sum_dn_diff = sum_dn_diff + flux_dn_diffuse(jg,jlev)
              sum_dn_dir  = sum_dn_dir  + flux_dn_direct(jg,jlev)
            end do
            flux%sw_up(jcol,jlev) = sum_up
            flux%sw_dn(jcol,jlev) = sum_dn_diff + sum_dn_dir
            if (allocated(flux%sw_dn_direct)) then
              flux%sw_dn_direct(jcol,jlev) = sum_dn_dir
            end if
          end do
#endif
          
          ! Cloudy flux profiles currently assume completely overcast
          ! skies; perform weighted average with clear-sky profile
          do jlev = 1, nlev+1
            flux%sw_up(jcol,jlev) =  total_cloud_cover *flux%sw_up(jcol,jlev) &
                 &     + (1.0_jprb - total_cloud_cover)*flux%sw_up_clear(jcol,jlev)
            flux%sw_dn(jcol,jlev) =  total_cloud_cover *flux%sw_dn(jcol,jlev) &
                 &     + (1.0_jprb - total_cloud_cover)*flux%sw_dn_clear(jcol,jlev)
            if (allocated(flux%sw_dn_direct)) then
              flux%sw_dn_direct(jcol,jlev) = total_cloud_cover *flux%sw_dn_direct(jcol,jlev) &
                   &  + (1.0_jprb - total_cloud_cover)*flux%sw_dn_direct_clear(jcol,jlev)
            end if
          end do
          ! Likewise for surface spectral fluxes
          do jg = 1,ng
            flux%sw_dn_diffuse_surf_g(jg,jcol) = flux_dn_diffuse(jg,nlev+1)
            flux%sw_dn_direct_surf_g(jg,jcol)  = flux_dn_direct(jg,nlev+1)
            flux%sw_dn_diffuse_surf_g(jg,jcol) = total_cloud_cover *flux%sw_dn_diffuse_surf_g(jg,jcol) &
                 &                 + (1.0_jprb - total_cloud_cover)*flux%sw_dn_diffuse_surf_clear_g(jg,jcol)
            flux%sw_dn_direct_surf_g(jg,jcol)  = total_cloud_cover *flux%sw_dn_direct_surf_g(jg,jcol) &
                 &                 + (1.0_jprb - total_cloud_cover)*flux%sw_dn_direct_surf_clear_g(jg,jcol)
          end do

        else
          ! No cloud in profile and clear-sky fluxes already
          ! calculated: copy them over
          do jlev = 1, nlev+1
            flux%sw_up(jcol,jlev) = flux%sw_up_clear(jcol,jlev)
            flux%sw_dn(jcol,jlev) = flux%sw_dn_clear(jcol,jlev)
            if (allocated(flux%sw_dn_direct)) then
              flux%sw_dn_direct(jcol,jlev) = flux%sw_dn_direct_clear(jcol,jlev)
            end if
          end do
          do jg = 1,ng
            flux%sw_dn_diffuse_surf_g(jg,jcol) = flux%sw_dn_diffuse_surf_clear_g(jg,jcol)
            flux%sw_dn_direct_surf_g(jg,jcol)  = flux%sw_dn_direct_surf_clear_g(jg,jcol)
          end do

        end if ! Cloud is present in profile

      else
        ! Set fluxes to zero if sun is below the horizon
        do jlev = 1, nlev+1
          flux%sw_up(jcol,jlev) = 0.0_jprb
          flux%sw_dn(jcol,jlev) = 0.0_jprb
          if (allocated(flux%sw_dn_direct)) then
            flux%sw_dn_direct(jcol,jlev) = 0.0_jprb
          end if
          flux%sw_up_clear(jcol,jlev) = 0.0_jprb
          flux%sw_dn_clear(jcol,jlev) = 0.0_jprb
          if (allocated(flux%sw_dn_direct_clear)) then
            flux%sw_dn_direct_clear(jcol,jlev) = 0.0_jprb
          end if
        end do
        do jg = 1,ng
          flux%sw_dn_diffuse_surf_g(jg,jcol) = 0.0_jprb
          flux%sw_dn_direct_surf_g(jg,jcol)  = 0.0_jprb
          flux%sw_dn_diffuse_surf_clear_g(jg,jcol) = 0.0_jprb
          flux%sw_dn_direct_surf_clear_g(jg,jcol)  = 0.0_jprb
        end do
      end if ! Sun above horizon

    end do ! Loop over columns

    if (lhook) call dr_hook('radiation_mcica_sw:solver_mcica_sw',1,hook_handle)
    
  end subroutine solver_mcica_sw

end module radiation_mcica_sw
