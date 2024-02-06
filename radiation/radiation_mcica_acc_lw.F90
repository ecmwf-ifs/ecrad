! radiation_mcica_acc_lw.F90 - Monte-Carlo Independent Column Approximation longtwave solver
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
!   2017-04-11  R. Hogan  Receive emission/albedo rather than planck/emissivity
!   2017-04-22  R. Hogan  Store surface fluxes at all g-points
!   2017-07-12  R. Hogan  Call fast adding method if only clouds scatter
!   2017-10-23  R. Hogan  Renamed single-character variables

#include "ecrad_config.h"

module radiation_mcica_acc_lw

  public

contains

  !---------------------------------------------------------------------
  ! Longwave Monte Carlo Independent Column Approximation
  ! (McICA). This implementation performs a clear-sky and a cloudy-sky
  ! calculation, and then weights the two to get the all-sky fluxes
  ! according to the total cloud cover. This method reduces noise for
  ! low cloud cover situations, and exploits the clear-sky
  ! calculations that are usually performed for diagnostic purposes
  ! simultaneously. The cloud generator has been carefully written
  ! such that the stochastic cloud field satisfies the prescribed
  ! overlap parameter accounting for this weighting.
  subroutine solver_mcica_acc_lw(nlev,istartcol,iendcol, &
       &  config, single_level, cloud, &
       &  od, ssa, g, od_cloud, ssa_cloud, g_cloud, planck_hl, &
       &  emission, albedo, &
       &  flux)

    use parkind1, only           : jprb
    use yomhook,  only           : lhook, dr_hook, jphook

    use radiation_io,   only           : nulerr, radiation_abort
    use radiation_config, only         : config_type
    use radiation_single_level, only   : single_level_type
    use radiation_cloud, only          : cloud_type
    use radiation_flux, only           : flux_type
    use radiation_two_stream, only     : calc_ref_trans_lw, &
         &                               calc_no_scattering_transmittance_lw
    use radiation_adding_ica_lw, only  : adding_ica_lw, fast_adding_ica_lw, &
         &                               calc_fluxes_no_scattering_lw
    use radiation_lw_derivatives, only : calc_lw_derivatives_ica, modify_lw_derivatives_ica
    use radiation_cloud_generator_acc, only: cloud_generator_acc
    use radiation_cloud_cover, only    : beta2alpha, MaxCloudFrac

    implicit none

    ! Inputs
    integer, intent(in) :: nlev               ! number of model levels
    integer, intent(in) :: istartcol, iendcol ! range of columns to process
    type(config_type),        intent(in) :: config
    type(single_level_type),  intent(in) :: single_level
    type(cloud_type),         intent(in) :: cloud

    ! Gas and aerosol optical depth, single-scattering albedo and
    ! asymmetry factor at each longwave g-point
    real(jprb), intent(in), dimension(config%n_g_lw, nlev, istartcol:iendcol) :: &
         &  od
    real(jprb), intent(in), dimension(config%n_g_lw_if_scattering, nlev, istartcol:iendcol) :: &
         &  ssa, g

    ! Cloud and precipitation optical depth, single-scattering albedo and
    ! asymmetry factor in each longwave band
    real(jprb), intent(in), dimension(config%n_bands_lw,nlev,istartcol:iendcol)   :: &
         &  od_cloud
    real(jprb), intent(in), dimension(config%n_bands_lw_if_scattering, &
         &  nlev,istartcol:iendcol) :: ssa_cloud, g_cloud

    ! Planck function at each half-level and the surface
    real(jprb), intent(in), dimension(config%n_g_lw,nlev+1,istartcol:iendcol) :: &
         &  planck_hl

    ! Emission (Planck*emissivity) and albedo (1-emissivity) at the
    ! surface at each longwave g-point
    real(jprb), intent(in), dimension(config%n_g_lw, istartcol:iendcol) :: emission, albedo

    ! Output
    type(flux_type), intent(inout):: flux

    ! Local variables

    ! Diffuse reflectance and transmittance for each layer in clear
    ! and all skies
    real(jprb), dimension(config%n_g_lw, nlev) :: ref_clear, trans_clear, reflectance, transmittance

    ! Emission by a layer into the upwelling or downwelling diffuse
    ! streams, in clear and all skies
    real(jprb), dimension(config%n_g_lw, nlev) :: source_up_clear, source_dn_clear, source_up, source_dn

    ! Fluxes per g point
    real(jprb), dimension(config%n_g_lw, nlev+1, istartcol:iendcol) :: flux_up, flux_dn
    real(jprb), dimension(config%n_g_lw, nlev+1, istartcol:iendcol) :: flux_up_clear, flux_dn_clear

    ! Combined gas+aerosol+cloud optical depth, single scattering
    ! albedo and asymmetry factor
    real(jprb), dimension(config%n_g_lw) :: od_total, ssa_total, g_total

    ! Combined scattering optical depth
    real(jprb) :: scat_od, scat_od_total(config%n_g_lw)

    ! Optical depth scaling from the cloud generator, zero indicating
    ! clear skies
    real(jprb), dimension(config%n_g_lw,nlev) :: od_scaling

    ! Modified optical depth after McICA scaling to represent cloud
    ! inhomogeneity
    real(jprb), dimension(config%n_g_lw) :: od_cloud_new

    ! Total cloud cover output from the cloud generator
    real(jprb) :: total_cloud_cover

    ! Identify clear-sky layers
    logical :: is_clear_sky_layer(nlev, istartcol:iendcol)

    ! workaround that allows inling of cloud generator
    real(jprb), dimension(config%pdf_sampler%ncdf, config%pdf_sampler%nfsd)  :: sample_val

    ! temporary arrays to increase performance
    real(jprb), dimension(nlev, istartcol:iendcol) :: frac, frac_std
    real(jprb), dimension(nlev-1, istartcol:iendcol) :: overlap_param

    ! temporary arrays
    real(jprb), dimension(nlev, istartcol:iendcol) :: cum_cloud_cover
    real(jprb), dimension(nlev-1, istartcol:iendcol) :: pair_cloud_cover

    ! "Alpha" overlap parameter
    real(jprb) :: overlap_alpha

    ! Cumulative product needed in computation of total_cloud_cover
    real(jprb) :: cum_product(istartcol:iendcol)

    ! First and last cloudy layers
    integer :: ibegin(istartcol:iendcol), iend(istartcol:iendcol)

    ! Temporary working array
    real(jprb), dimension(config%n_g_lw,nlev+1) :: tmp_work_albedo, &
      &                                            tmp_work_source
    real(jprb), dimension(config%n_g_lw,nlev) :: tmp_work_inv_denominator
    ! Temporary storage for more efficient summation
    real(jprb) :: sum_up, sum_dn, sum_up_clr, sum_dn_clr

    ! Index of the highest cloudy layer
    integer :: i_cloud_top

    ! Number of g points
    integer :: ng

    ! Loop indices for level, column and g point
    integer :: jlev, jcol, jg

    real(jphook) :: hook_handle

    if (lhook) call dr_hook('radiation_mcica_acc_lw:solver_mcica_acc_lw',0,hook_handle)

    if (.not. config%do_clear) then
      write(nulerr,'(a)') '*** Error: longwave McICA ACC requires clear-sky calculation to be performed'
      call radiation_abort()
    end if

    ng = config%n_g_lw

    !$ACC DATA CREATE(flux_up, flux_dn, flux_up_clear, flux_dn_clear, &
    !$ACC             is_clear_sky_layer, &
    !$ACC             sample_val, frac, frac_std, overlap_param, cum_cloud_cover, &
    !$ACC             pair_cloud_cover, cum_product, ibegin, iend) &
    !$ACC     PRESENT(config, single_level, cloud, od, ssa, g, od_cloud, ssa_cloud, &
    !$ACC             g_cloud, planck_hl, emission, albedo, flux)

    !$ACC PARALLEL DEFAULT(NONE) ASYNC(1)
    !$ACC LOOP GANG VECTOR COLLAPSE(2)
    do jlev = 1,config%pdf_sampler%nfsd
      do jcol = 1,config%pdf_sampler%ncdf
        sample_val(jcol,jlev) = config%pdf_sampler%val(jcol,jlev)
      end do
    end do
    !$ACC END PARALLEL

    !$ACC PARALLEL DEFAULT(NONE) ASYNC(1)
    !$ACC LOOP GANG VECTOR COLLAPSE(2)
    do jcol = istartcol,iendcol
      do jlev = 1, nlev
        frac(jlev, jcol) = cloud%fraction(jcol,jlev)
        frac_std(jlev, jcol) = cloud%fractional_std(jcol,jlev)
        is_clear_sky_layer(jlev,jcol) = .true.
      end do
    end do
    !$ACC END PARALLEL

    !$ACC PARALLEL DEFAULT(NONE) ASYNC(1)
    !$ACC LOOP GANG VECTOR COLLAPSE(2)
    do jcol = istartcol,iendcol
      do jlev = 1, nlev-1
        overlap_param(jlev, jcol) = cloud%overlap_param(jcol,jlev)
      end do
    end do
    !$ACC END PARALLEL


    !$ACC PARALLEL DEFAULT(NONE) ASYNC(1)
    !$ACC LOOP GANG VECTOR COLLAPSE(2) PRIVATE(overlap_alpha)
    do jcol = istartcol,iendcol
      !---------------------------------------------------------------------
      ! manual inline from cum_cloud_cover_exp_ran >>>>>>>>>>>>>>>>>>>>>>>>
      ! Loop to compute total cloud cover and the cumulative cloud cover
      ! down to the base of each layer
      do jlev = 1,nlev-1
        ! Convert to "alpha" overlap parameter if necessary
        if (config%use_beta_overlap) then
          overlap_alpha = beta2alpha(overlap_param(jlev,jcol), &
                &                     frac(jlev,jcol), frac(jlev+1,jcol))
        else
          overlap_alpha = overlap_param(jlev,jcol)
        end if

        ! Compute the combined cloud cover of layers jlev and jlev+1
        pair_cloud_cover(jlev, jcol) = overlap_alpha*max(frac(jlev,jcol),frac(jlev+1,jcol)) &
              &  + (1.0_jprb - overlap_alpha) &
              &  * (frac(jlev,jcol)+frac(jlev+1,jcol)-frac(jlev,jcol)*frac(jlev+1,jcol))
      end do
    end do
    !$ACC END PARALLEL

    !$ACC PARALLEL DEFAULT(NONE) ASYNC(1)
    !$ACC LOOP GANG VECTOR
    do jcol = istartcol,iendcol
      cum_cloud_cover(1, jcol) = frac(1,jcol)
      cum_product(jcol) = 1.0_jprb - frac(1,jcol)
      !$ACC LOOP SEQ
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
      flux%cloud_cover_lw(jcol) = cum_cloud_cover(nlev,jcol);
      if (flux%cloud_cover_lw(jcol) < config%cloud_fraction_threshold) then
        ! Treat column as clear sky: calling function therefore will not
        ! use od_scaling so we don't need to calculate it
        flux%cloud_cover_lw(jcol) = 0.0_jprb
      end if
    end do
    !$ACC END PARALLEL

    !$ACC PARALLEL DEFAULT(NONE) ASYNC(1)
    !$ACC LOOP GANG VECTOR
    do jcol = istartcol,iendcol
      if (flux%cloud_cover_lw(jcol) >= config%cloud_fraction_threshold) then
        ! Cloud is present: need to calculate od_scaling
        ! Find range of cloudy layers
        ibegin(jcol) = nlev
        !$ACC LOOP SEQ
        do jlev = 1, nlev
          if( frac(jlev,jcol) > 0.0_jprb ) then
            ibegin(jcol) = min(jlev, ibegin(jcol))
          end if
        end do

        iend(jcol) = ibegin(jcol)
        !$ACC LOOP SEQ
        do jlev = ibegin(jcol)+1,nlev
          if (frac(jlev,jcol) > 0.0_jprb) then
            iend(jcol) = max(jlev, iend(jcol))
          end if
        end do
      end if
    end do
    !$ACC END PARALLEL


    ! Loop through columns
    !$ACC PARALLEL DEFAULT(NONE) ASYNC(1)  &
    !$ACC   NUM_GANGS(iendcol-istartcol+1) NUM_WORKERS((config%n_g_lw-1)/32+1) VECTOR_LENGTH(32)
    !$ACC LOOP GANG PRIVATE(g_total, i_cloud_top, od_cloud_new, od_scaling, od_total, ref_clear, &
    !$ACC   reflectance, tmp_work_inv_denominator, tmp_work_albedo, tmp_work_source, &
    !$ACC   scat_od_total, source_dn, source_dn_clear, source_up, source_up_clear, ssa_total, &
    !$ACC   trans_clear, transmittance)
    do jcol = istartcol,iendcol

      ! Clear-sky calculation
#ifndef _OPENACC
      if (config%do_lw_aerosol_scattering) then
        ! Scattering case: first compute clear-sky reflectance,
        ! transmittance etc at each model level
        call calc_ref_trans_lw(ng*nlev, &
             &  od(:,:,jcol), ssa(:,:,jcol), g(:,:,jcol), &
             &  planck_hl(:,1:jlev,jcol), planck_hl(:,2:jlev+1,jcol), &
             &  ref_clear, trans_clear, &
             &  source_up_clear, source_dn_clear)
        ! Then use adding method to compute fluxes
        call adding_ica_lw(ng, nlev, &
             &  ref_clear, trans_clear, source_up_clear, source_dn_clear, &
             &  emission(:,jcol), albedo(:,jcol), &
             &  flux_up_clear(:,:,jcol), flux_dn_clear(:,:,jcol))

      else
#endif
        ! Non-scattering case: use simpler functions for
        ! transmission and emission
        call calc_no_scattering_transmittance_lw(ng*nlev, od(:,:,jcol), &
               &  planck_hl(:,1:nlev,jcol), planck_hl(:,2:nlev+1, jcol), &
               &  trans_clear, source_up_clear, &
               &  source_dn_clear)
        ! Simpler down-then-up method to compute fluxes
        call calc_fluxes_no_scattering_lw(ng, nlev, &
             &  trans_clear, source_up_clear, source_dn_clear, &
             &  emission(:,jcol), albedo(:,jcol), &
             &  flux_up_clear(:,:,jcol), flux_dn_clear(:,:,jcol))

        ! Ensure that clear-sky reflectance is zero since it may be
        ! used in cloudy-sky case
        !$ACC LOOP SEQ
        do jlev = 1,nlev
          !$ACC LOOP WORKER VECTOR
          do jg = 1,ng
            ref_clear(jg,jlev) = 0.0_jprb
          end do
        end do
#ifndef _OPENACC
      end if
#endif

      ! Store surface spectral downwelling fluxes
      !$ACC LOOP WORKER VECTOR
      do jg = 1,ng
        flux%lw_dn_surf_clear_g(jg,jcol) = flux_dn_clear(jg,nlev+1,jcol)
      end do

      ! Do cloudy-sky calculation; add a prime number to the seed in
      ! the longwave
      call cloud_generator_acc(ng, nlev, &
           &  single_level%iseed(jcol) + 997, &
           &  config%cloud_fraction_threshold, &
           &  frac(:,jcol), overlap_param(:,jcol), &
           &  config%cloud_inhom_decorr_scaling, frac_std(:,jcol), &
           &  config%pdf_sampler%ncdf, config%pdf_sampler%nfsd, &
           &  config%pdf_sampler%fsd1, config%pdf_sampler%inv_fsd_interval, &
           &  sample_val, &
           &  od_scaling, flux%cloud_cover_lw(jcol)+0.0_jprb, & ! Workaround for nvhpc-24.1
           &  ibegin(jcol), iend(jcol), &
           &  cum_cloud_cover=cum_cloud_cover(:,jcol), &
           &  pair_cloud_cover=pair_cloud_cover(:,jcol))

      if (flux%cloud_cover_lw(jcol) >= config%cloud_fraction_threshold) then
        ! Total-sky calculation

        i_cloud_top = nlev+1
        !$ACC LOOP SEQ
        do jlev = 1,nlev
          ! Compute combined gas+aerosol+cloud optical properties
          if (frac(jlev,jcol) >= config%cloud_fraction_threshold) then
            is_clear_sky_layer(jlev,jcol) = .false.
            ! Get index to the first cloudy layer from the top
            if (i_cloud_top > jlev) then
              i_cloud_top = jlev
            end if

            !$ACC LOOP WORKER VECTOR
            do jg = 1,ng
              od_cloud_new(jg) = od_scaling(jg,jlev) &
                 &  * od_cloud(config%i_band_from_reordered_g_lw(jg),jlev,jcol)
              od_total(jg)  = od(jg,jlev,jcol) + od_cloud_new(jg)
              ssa_total(jg) = 0.0_jprb
              g_total(jg)   = 0.0_jprb
            end do

            if (config%do_lw_cloud_scattering) then
              ! Scattering case: calculate reflectance and
              ! transmittance at each model level

#ifndef _OPENACC
              if (config%do_lw_aerosol_scattering) then
                ! In single precision we need to protect against the
                ! case that od_total > 0.0 and ssa_total > 0.0 but
                ! od_total*ssa_total == 0 due to underflow
                do jg = 1,ng
                  if (od_total(jg) > 0.0_jprb) then
                    scat_od_total(jg) = ssa(jg,jlev,jcol)*od(jg,jlev,jcol) &
                     &     + ssa_cloud(config%i_band_from_reordered_g_lw(jg),jlev,jcol) &
                     &     *  od_cloud_new(jg)
                    ssa_total(jg) = scat_od_total(jg) / od_total(jg)

                    if (scat_od_total(jg) > 0.0_jprb) then
                      g_total(jg) = (g(jg,jlev,jcol)*ssa(jg,jlev,jcol)*od(jg,jlev,jcol) &
                         &     +   g_cloud(config%i_band_from_reordered_g_lw(jg),jlev,jcol) &
                         &     * ssa_cloud(config%i_band_from_reordered_g_lw(jg),jlev,jcol) &
                         &     *  od_cloud_new(jg)) &
                         &     / scat_od_total(jg)
                    end if
                  end if
                end do

              else
#endif

                !$ACC LOOP WORKER VECTOR PRIVATE(scat_od)
                do jg = 1,ng
                  if (od_total(jg) > 0.0_jprb) then
                    scat_od = ssa_cloud(config%i_band_from_reordered_g_lw(jg),jlev,jcol) &
                         &     * od_cloud_new(jg)
                    ssa_total(jg) = scat_od / od_total(jg)
                    if (scat_od > 0.0_jprb) then
                      g_total(jg) = g_cloud(config%i_band_from_reordered_g_lw(jg),jlev,jcol) &
                           &     * ssa_cloud(config%i_band_from_reordered_g_lw(jg),jlev,jcol) &
                           &     *  od_cloud_new(jg) / scat_od
                    end if
                  end if
                end do

#ifndef _OPENACC
              end if
#endif

              ! Compute cloudy-sky reflectance, transmittance etc at
              ! each model level
              call calc_ref_trans_lw(ng, &
                   &  od_total, ssa_total, g_total, &
                   &  planck_hl(:,jlev,jcol), planck_hl(:,jlev+1,jcol), &
                   &  reflectance(:,jlev), transmittance(:,jlev), &
                   &  source_up(:,jlev), source_dn(:,jlev))
            else
              ! No-scattering case: use simpler functions for
              ! transmission and emission
              call calc_no_scattering_transmittance_lw(ng, od_total, &
                   &  planck_hl(:,jlev,jcol), planck_hl(:,jlev+1, jcol), &
                   &  transmittance(:,jlev), source_up(:,jlev), &
                   &  source_dn(:,jlev))
            end if

          else
            ! Clear-sky layer: copy over clear-sky values
            !$ACC LOOP WORKER VECTOR
            do jg = 1,ng
              reflectance(jg,jlev) = ref_clear(jg,jlev)
              transmittance(jg,jlev) = trans_clear(jg,jlev)
              source_up(jg,jlev) = source_up_clear(jg,jlev)
              source_dn(jg,jlev) = source_dn_clear(jg,jlev)
            end do
          end if
        end do

#ifndef _OPENACC
        if (config%do_lw_aerosol_scattering) then
          ! Use adding method to compute fluxes for an overcast sky,
          ! allowing for scattering in all layers
          call adding_ica_lw(ng, nlev, reflectance, transmittance, source_up, source_dn, &
               &  emission(:,jcol), albedo(:,jcol), &
               &  flux_up(:,:,jcol), flux_dn(:,:,jcol))
        else if (config%do_lw_cloud_scattering) then
#else
        if(config%do_lw_cloud_scattering) then
#endif
          ! Use adding method to compute fluxes but optimize for the
          ! presence of clear-sky layers
          call fast_adding_ica_lw(ng, nlev, reflectance, transmittance, source_up, source_dn, &
               &  emission(:,jcol), albedo(:,jcol), &
               &  is_clear_sky_layer(:,jcol), i_cloud_top, flux_dn_clear(:,:,jcol), &
               &  flux_up(:,:,jcol), flux_dn(:,:,jcol), &
               &  albedo=tmp_work_albedo, &
               &  source=tmp_work_source, &
               &  inv_denominator=tmp_work_inv_denominator)
        else
          ! Simpler down-then-up method to compute fluxes
          call calc_fluxes_no_scattering_lw(ng, nlev, &
               &  transmittance, source_up, source_dn, emission(:,jcol), albedo(:,jcol), &
               &  flux_up(:,:,jcol), flux_dn(:,:,jcol))
        end if

        ! Cloudy flux profiles currently assume completely overcast
        ! skies; perform weighted average with clear-sky profile
        ! Store surface spectral downwelling fluxes
        !$ACC LOOP WORKER VECTOR
        do jg = 1,ng
          flux%lw_dn_surf_g(jg,jcol) = flux%cloud_cover_lw(jcol)*flux_dn(jg,nlev+1,jcol) &
              &  + (1.0_jprb - flux%cloud_cover_lw(jcol))*flux%lw_dn_surf_clear_g(jg,jcol)
        end do

        ! Compute the longwave derivatives needed by Hogan and Bozzo
        ! (2015) approximate radiation update scheme
        if (config%do_lw_derivatives) then
          call calc_lw_derivatives_ica(ng, nlev, jcol, transmittance, flux_up(:,nlev+1,jcol), &
               &                       flux%lw_derivatives)
          if (flux%cloud_cover_lw(jcol) < 1.0_jprb - config%cloud_fraction_threshold) then
            ! Modify the existing derivative with the contribution from the clear sky
            call modify_lw_derivatives_ica(ng, nlev, jcol, trans_clear, flux_up_clear(:,nlev+1,jcol), &
                 &                         1.0_jprb-flux%cloud_cover_lw(jcol), flux%lw_derivatives)
          end if
        end if

      else
        ! No cloud in profile and clear-sky fluxes already
        ! calculated: copy them over
        !$ACC LOOP WORKER VECTOR
        do jg = 1,ng
          flux%lw_dn_surf_g(jg,jcol) = flux%lw_dn_surf_clear_g(jg,jcol)
        end do
        if (config%do_lw_derivatives) then
          call calc_lw_derivatives_ica(ng, nlev, jcol, trans_clear, flux_up_clear(:,nlev+1,jcol), &
               &                       flux%lw_derivatives)

        end if
      end if ! Cloud is present in profile
    end do
    !$ACC END PARALLEL

    ! Loop through columns
    !$ACC PARALLEL DEFAULT(NONE) ASYNC(1) &
    !$ACC   NUM_GANGS((iendcol-istartcol+1)*(nlev+1)) NUM_WORKERS(1) VECTOR_LENGTH(32*((config%n_g_lw-1)/32+1))
    !$ACC LOOP GANG COLLAPSE(2) PRIVATE(total_cloud_cover, sum_up, sum_dn, sum_up_clr, sum_dn_clr)
    do jcol = istartcol,iendcol
      do jlev = 1,nlev+1

        sum_up_clr = 0._jprb
        sum_dn_clr = 0._jprb
        !$ACC LOOP VECTOR REDUCTION(+:sum_up_clr, sum_dn_clr)
        do jg = 1,ng
          sum_up_clr = sum_up_clr + flux_up_clear(jg,jlev,jcol)
          sum_dn_clr = sum_dn_clr + flux_dn_clear(jg,jlev,jcol)
        end do
        flux%lw_up_clear(jcol,jlev) = sum_up_clr
        flux%lw_dn_clear(jcol,jlev) = sum_dn_clr

        total_cloud_cover = flux%cloud_cover_lw(jcol)
        if (total_cloud_cover >= config%cloud_fraction_threshold) then

          ! Store overcast broadband fluxes
          sum_up = 0._jprb
          sum_dn = 0._jprb
          !$ACC LOOP VECTOR REDUCTION(+:sum_up, sum_dn)
          do jg = 1,ng
            sum_up = sum_up + flux_up(jg,jlev,jcol)
            sum_dn = sum_dn + flux_dn(jg,jlev,jcol)
          end do
          flux%lw_up(jcol,jlev) = total_cloud_cover*sum_up + (1.0_jprb - total_cloud_cover)*sum_up_clr
          flux%lw_dn(jcol,jlev) = total_cloud_cover*sum_dn + (1.0_jprb - total_cloud_cover)*sum_dn_clr

        else

          flux%lw_up(jcol,jlev) = sum_up_clr
          flux%lw_dn(jcol,jlev) = sum_dn_clr

        end if
      end do
    end do
    !$ACC END PARALLEL

    !$ACC WAIT
    !$ACC END DATA

    if (lhook) call dr_hook('radiation_mcica_acc_lw:solver_mcica_acc_lw',1,hook_handle)

  end subroutine solver_mcica_acc_lw

end module radiation_mcica_acc_lw
