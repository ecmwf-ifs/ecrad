! radiation_mcica_lw.F90 - Monte-Carlo Independent Column Approximation longtwave solver
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

module radiation_mcica_lw

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
  subroutine solver_mcica_lw(nlev,istartcol,iendcol, &
       &  config, single_level, cloud, & 
       &  od, ssa, g, od_cloud, ssa_cloud, g_cloud, planck_hl, &
       &  emission, albedo, &
       &  flux)

    use parkind1, only           : jprb
    use yomhook,  only           : lhook, dr_hook

    use radiation_io,   only           : nulerr, radiation_abort
    use radiation_config, only         : config_type
    use radiation_single_level, only   : single_level_type
    use radiation_cloud, only          : cloud_type
    use radiation_flux, only           : flux_type
    use radiation_two_stream, only     : calc_two_stream_gammas_lw, &
         &                               calc_reflectance_transmittance_lw, &
         &                               calc_no_scattering_transmittance_lw
    use radiation_adding_ica_lw, only  : adding_ica_lw, fast_adding_ica_lw, &
         &                               calc_fluxes_no_scattering_lw
    use radiation_lw_derivatives, only : calc_lw_derivatives_ica, modify_lw_derivatives_ica
    use radiation_cloud_generator, only: cloud_generator
    use random_numbers_mix, only       : randomnumberstream

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
    real(jprb), dimension(config%n_g_lw, nlev, istartcol:iendcol) :: ref_clear, trans_clear, reflectance, transmittance

    ! Emission by a layer into the upwelling or downwelling diffuse
    ! streams, in clear and all skies
    real(jprb), dimension(config%n_g_lw, nlev, istartcol:iendcol) :: source_up_clear, source_dn_clear, source_up, source_dn

    ! Fluxes per g point
    real(jprb), dimension(config%n_g_lw, nlev+1, istartcol:iendcol) :: flux_up, flux_dn
    real(jprb), dimension(config%n_g_lw, nlev+1, istartcol:iendcol) :: flux_up_clear, flux_dn_clear

    ! Combined gas+aerosol+cloud optical depth, single scattering
    ! albedo and asymmetry factor
    real(jprb), dimension(config%n_g_lw,istartcol:iendcol) :: od_total, ssa_total, g_total

    ! Combined scattering optical depth
    real(jprb) :: scat_od, scat_od_total(config%n_g_lw,istartcol:iendcol)

    ! Two-stream coefficients
    real(jprb), dimension(config%n_g_lw,istartcol:iendcol) :: gamma1, gamma2

    ! Optical depth scaling from the cloud generator, zero indicating
    ! clear skies
    real(jprb), dimension(config%n_g_lw,nlev,istartcol:iendcol) :: od_scaling

    ! Modified optical depth after McICA scaling to represent cloud
    ! inhomogeneity
    real(jprb), dimension(config%n_g_lw, istartcol:iendcol) :: od_cloud_new

    ! Total cloud cover output from the cloud generator
    real(jprb) :: total_cloud_cover

    ! Identify clear-sky layers
    logical :: is_clear_sky_layer(nlev, istartcol:iendcol)

    ! Temporary working array
    real(jprb), dimension(nlev, istartcol:iendcol) :: tmp_work_nlev
    real(jprb), dimension(config%n_g_lw, istartcol:iendcol) :: tmp_work_ng
    real(jprb), dimension(nlev-1, istartcol:iendcol) :: tmp_work_nlevm1, &
      &                                                 tmp_work_nlevm2, &
      &                                                 tmp_work_nlevm3
    real(jprb), dimension(config%n_g_lw,nlev+1,istartcol:iendcol) :: tmp_work_ngnlevp1, &
      &                                                              tmp_work_ngnlevp2
    real(jprb), dimension(config%n_g_lw,nlev,istartcol:iendcol) :: tmp_work_ngnlev1, &
      &                                                            tmp_work_ngnlev2, &
      &                                                            tmp_work_ngnlev3
    real(jprb), dimension(999, istartcol:iendcol) :: tmp_work_jpwarmup_lf
    type(randomnumberstream), dimension(istartcol:iendcol) :: random_stream

    ! Index of the highest cloudy layer
    integer :: i_cloud_top

    ! Number of g points
    integer :: ng

    ! Loop indices for level, column and g point
    integer :: jlev, jcol, jg

    real(jprb) :: hook_handle

    if (lhook) call dr_hook('radiation_mcica_lw:solver_mcica_lw',0,hook_handle)

    if (.not. config%do_clear) then
      write(nulerr,'(a)') '*** Error: longwave McICA requires clear-sky calculation to be performed'
      call radiation_abort()      
    end if

    ng = config%n_g_lw

    ! Loop through columns
    do jcol = istartcol,iendcol

      ! Clear-sky calculation
      if (config%do_lw_aerosol_scattering) then
        ! Scattering case: first compute clear-sky reflectance,
        ! transmittance etc at each model level
        do jlev = 1,nlev
          ssa_total(:,jcol) = ssa(:,jlev,jcol)
          g_total(:,jcol)   = g(:,jlev,jcol)
          call calc_two_stream_gammas_lw(ng, ssa_total(:,jcol), g_total(:,jcol), &
               &  gamma1(:,jcol), gamma2(:,jcol))
          call calc_reflectance_transmittance_lw(ng, &
               &  od(:,jlev,jcol), gamma1(:,jcol), gamma2(:,jcol), &
               &  planck_hl(:,jlev,jcol), planck_hl(:,jlev+1,jcol), &
               &  ref_clear(:,jlev,jcol), trans_clear(:,jlev,jcol), &
               &  source_up_clear(:,jlev,jcol), source_dn_clear(:,jlev,jcol))
        end do
        ! Then use adding method to compute fluxes
        call adding_ica_lw(ng, nlev, &
             &  ref_clear(:,:,jcol), trans_clear(:,:,jcol), source_up_clear(:,:,jcol), source_dn_clear(:,:,jcol), &
             &  emission(:,jcol), albedo(:,jcol), &
             &  flux_up_clear(:,:,jcol), flux_dn_clear(:,:,jcol))
        
      else
        ! Non-scattering case: use simpler functions for
        ! transmission and emission
        do jlev = 1,nlev
          call calc_no_scattering_transmittance_lw(ng, od(:,jlev,jcol), &
               &  planck_hl(:,jlev,jcol), planck_hl(:,jlev+1, jcol), &
               &  trans_clear(:,jlev,jcol), source_up_clear(:,jlev,jcol), source_dn_clear(:,jlev,jcol))
        end do
        ! Simpler down-then-up method to compute fluxes
        call calc_fluxes_no_scattering_lw(ng, nlev, &
             &  trans_clear(:,:,jcol), source_up_clear(:,:,jcol), source_dn_clear(:,:,jcol), &
             &  emission(:,jcol), albedo(:,jcol), &
             &  flux_up_clear(:,:,jcol), flux_dn_clear(:,:,jcol))
        
        ! Ensure that clear-sky reflectance is zero since it may be
        ! used in cloudy-sky case
        ref_clear(:,:,jcol) = 0.0_jprb
      end if

      ! Sum over g-points to compute broadband fluxes
      flux%lw_up_clear(jcol,:) = sum(flux_up_clear(:,:,jcol),1)
      flux%lw_dn_clear(jcol,:) = sum(flux_dn_clear(:,:,jcol),1)
      ! Store surface spectral downwelling fluxes
      flux%lw_dn_surf_clear_g(:,jcol) = flux_dn_clear(:,nlev+1,jcol)

      ! Do cloudy-sky calculation; add a prime number to the seed in
      ! the longwave
      call cloud_generator(ng, nlev, config%i_overlap_scheme, &
           &  single_level%iseed(jcol) + 997, &
           &  config%cloud_fraction_threshold, &
           &  cloud%fraction(jcol,:), cloud%overlap_param(jcol,:), &
           &  config%cloud_inhom_decorr_scaling, cloud%fractional_std(jcol,:), &
           &  config%pdf_sampler, od_scaling(:,:,jcol), total_cloud_cover, &
           &  cum_cloud_cover=tmp_work_nlev(:,jcol), &
           &  rand_top=tmp_work_ng(:,jcol), &
           &  overlap_param_inhom=tmp_work_nlevm1(:,jcol), &
           &  pair_cloud_cover=tmp_work_nlevm2(:,jcol), &
           &  overhang=tmp_work_nlevm3(:,jcol), &
           &  tmp_work_ngnlev1=tmp_work_ngnlev1(:,:,jcol), &
           &  tmp_work_ngnlev2=tmp_work_ngnlev2(:,:,jcol), &
           &  tmp_work_ngnlev3=tmp_work_ngnlev3(:,:,jcol), &
           &  tmp_work_jpwarmup_lfg=tmp_work_jpwarmup_lf(:,jcol), &
           &  random_stream=random_stream(jcol), &
           &  use_beta_overlap=config%use_beta_overlap, &
           &  use_vectorizable_generator=config%use_vectorizable_generator)
      
      ! Store total cloud cover
      flux%cloud_cover_lw(jcol) = total_cloud_cover
      
      if (total_cloud_cover >= config%cloud_fraction_threshold) then
        ! Total-sky calculation

        is_clear_sky_layer(:,jcol) = .true.
        i_cloud_top = nlev+1
        do jlev = 1,nlev
          ! Compute combined gas+aerosol+cloud optical properties
          if (cloud%fraction(jcol,jlev) >= config%cloud_fraction_threshold) then
            is_clear_sky_layer(jlev,jcol) = .false.
            ! Get index to the first cloudy layer from the top
            if (i_cloud_top > jlev) then
              i_cloud_top = jlev
            end if

            do jg = 1,ng
              od_cloud_new(jg,jcol) = od_scaling(jg,jlev,jcol) &
                 &  * od_cloud(config%i_band_from_reordered_g_lw(jg),jlev,jcol)
              od_total(jg,jcol)  = od(jg,jlev,jcol) + od_cloud_new(jg,jcol)
              ssa_total(jg,jcol) = 0.0_jprb
              g_total(jg,jcol)   = 0.0_jprb
            end do

            if (config%do_lw_cloud_scattering) then
              ! Scattering case: calculate reflectance and
              ! transmittance at each model level

              if (config%do_lw_aerosol_scattering) then
                ! In single precision we need to protect against the
                ! case that od_total > 0.0 and ssa_total > 0.0 but
                ! od_total*ssa_total == 0 due to underflow
                do jg = 1,ng
                  if (od_total(jg,jcol) > 0.0_jprb) then
                    scat_od_total(jg,jcol) = ssa(jg,jlev,jcol)*od(jg,jlev,jcol) &
                     &     + ssa_cloud(config%i_band_from_reordered_g_lw(jg),jlev,jcol) &
                     &     *  od_cloud_new(jg,jcol)
                    ssa_total(jg,jcol) = scat_od_total(jg,jcol) / od_total(jg,jcol)

                    if (scat_od_total(jg,jcol) > 0.0_jprb) then
                      g_total(jg,jcol) = (g(jg,jlev,jcol)*ssa(jg,jlev,jcol)*od(jg,jlev,jcol) &
                         &     +   g_cloud(config%i_band_from_reordered_g_lw(jg),jlev,jcol) &
                         &     * ssa_cloud(config%i_band_from_reordered_g_lw(jg),jlev,jcol) &
                         &     *  od_cloud_new(jg,jcol)) &
                         &     / scat_od_total(jg,jcol)
                    end if
                  end if
                end do

              else

                do jg = 1,ng
                  if (od_total(jg,jcol) > 0.0_jprb) then
                    scat_od = ssa_cloud(config%i_band_from_reordered_g_lw(jg),jlev,jcol) &
                         &     * od_cloud_new(jg,jcol)
                    ssa_total(jg,jcol) = scat_od / od_total(jg,jcol)
                    if (scat_od > 0.0_jprb) then
                      g_total(jg,jcol) = g_cloud(config%i_band_from_reordered_g_lw(jg),jlev,jcol) &
                           &     * ssa_cloud(config%i_band_from_reordered_g_lw(jg),jlev,jcol) &
                           &     *  od_cloud_new(jg,jcol) / scat_od
                    end if
                  end if
                end do

              end if
            
              ! Compute cloudy-sky reflectance, transmittance etc at
              ! each model level
              call calc_two_stream_gammas_lw(ng, ssa_total(:,jcol), g_total(:,jcol), &
                   &  gamma1(:,jcol), gamma2(:,jcol))
              call calc_reflectance_transmittance_lw(ng, &
                   &  od_total(:,jcol), gamma1(:,jcol), gamma2(:,jcol), &
                   &  planck_hl(:,jlev,jcol), planck_hl(:,jlev+1,jcol), &
                   &  reflectance(:,jlev,jcol), transmittance(:,jlev,jcol), source_up(:,jlev,jcol), source_dn(:,jlev,jcol))
            else
              ! No-scattering case: use simpler functions for
              ! transmission and emission
              call calc_no_scattering_transmittance_lw(ng, od_total(:,jcol), &
                   &  planck_hl(:,jlev,jcol), planck_hl(:,jlev+1, jcol), &
                   &  transmittance(:,jlev,jcol), source_up(:,jlev,jcol), source_dn(:,jlev,jcol))
            end if

          else
            ! Clear-sky layer: copy over clear-sky values
            reflectance(:,jlev,jcol) = ref_clear(:,jlev,jcol)
            transmittance(:,jlev,jcol) = trans_clear(:,jlev,jcol)
            source_up(:,jlev,jcol) = source_up_clear(:,jlev,jcol)
            source_dn(:,jlev,jcol) = source_dn_clear(:,jlev,jcol)
          end if
        end do
        
        if (config%do_lw_aerosol_scattering) then
          ! Use adding method to compute fluxes for an overcast sky,
          ! allowing for scattering in all layers
          call adding_ica_lw(ng, nlev, reflectance(:,:,jcol), transmittance(:,:,jcol), source_up(:,:,jcol), source_dn(:,:,jcol), &
               &  emission(:,jcol), albedo(:,jcol), &
               &  flux_up(:,:,jcol), flux_dn(:,:,jcol))
        else if (config%do_lw_cloud_scattering) then
          ! Use adding method to compute fluxes but optimize for the
          ! presence of clear-sky layers
!          call adding_ica_lw(ng, nlev, reflectance(:,:,jcol), transmittance(:,:,jcol), source_up(:,:,jcol), source_dn(:,:,jcol), &
!               &  emission(:,jcol), albedo(:,jcol), &
!               &  flux_up(:,:,jcol), flux_dn(:,:,jcol))
          call fast_adding_ica_lw(ng, nlev, reflectance(:,:,jcol), transmittance(:,:,jcol), source_up(:,:,jcol), source_dn(:,:,jcol), &
               &  emission(:,jcol), albedo(:,jcol), &
               &  is_clear_sky_layer(:,jcol), i_cloud_top, flux_dn_clear(:,:,jcol), &
               &  flux_up(:,:,jcol), flux_dn(:,:,jcol), &
               &  albedo=tmp_work_ngnlevp1(:,:,jcol), &
               &  source=tmp_work_ngnlevp2(:,:,jcol), &
               &  inv_denominator=tmp_work_ngnlev1(:,:,jcol))
        else
          ! Simpler down-then-up method to compute fluxes
          call calc_fluxes_no_scattering_lw(ng, nlev, &
               &  transmittance(:,:,jcol), source_up(:,:,jcol), source_dn(:,:,jcol), emission(:,jcol), albedo(:,jcol), &
               &  flux_up(:,:,jcol), flux_dn(:,:,jcol))
        end if
        
        ! Store overcast broadband fluxes
        flux%lw_up(jcol,:) = sum(flux_up(:,:,jcol),1)
        flux%lw_dn(jcol,:) = sum(flux_dn(:,:,jcol),1)

        ! Cloudy flux profiles currently assume completely overcast
        ! skies; perform weighted average with clear-sky profile
        flux%lw_up(jcol,:) =  total_cloud_cover *flux%lw_up(jcol,:) &
             &  + (1.0_jprb - total_cloud_cover)*flux%lw_up_clear(jcol,:)
        flux%lw_dn(jcol,:) =  total_cloud_cover *flux%lw_dn(jcol,:) &
             &  + (1.0_jprb - total_cloud_cover)*flux%lw_dn_clear(jcol,:)
        ! Store surface spectral downwelling fluxes
        flux%lw_dn_surf_g(:,jcol) = total_cloud_cover*flux_dn(:,nlev+1,jcol) &
             &  + (1.0_jprb - total_cloud_cover)*flux%lw_dn_surf_clear_g(:,jcol)

        ! Compute the longwave derivatives needed by Hogan and Bozzo
        ! (2015) approximate radiation update scheme
        if (config%do_lw_derivatives) then
          call calc_lw_derivatives_ica(ng, nlev, jcol, transmittance(:,:,jcol), flux_up(:,nlev+1,jcol), &
               &                       flux%lw_derivatives)
          if (total_cloud_cover < 1.0_jprb - config%cloud_fraction_threshold) then
            ! Modify the existing derivative with the contribution from the clear sky
            call modify_lw_derivatives_ica(ng, nlev, jcol, trans_clear(:,:,jcol), flux_up_clear(:,nlev+1,jcol), &
                 &                         1.0_jprb-total_cloud_cover, flux%lw_derivatives)
          end if
        end if

      else
        ! No cloud in profile and clear-sky fluxes already
        ! calculated: copy them over
        flux%lw_up(jcol,:) = flux%lw_up_clear(jcol,:)
        flux%lw_dn(jcol,:) = flux%lw_dn_clear(jcol,:)
        flux%lw_dn_surf_g(:,jcol) = flux%lw_dn_surf_clear_g(:,jcol)
        if (config%do_lw_derivatives) then
          call calc_lw_derivatives_ica(ng, nlev, jcol, trans_clear(:,:,jcol), flux_up_clear(:,nlev+1,jcol), &
               &                       flux%lw_derivatives)
 
        end if
      end if ! Cloud is present in profile
    end do

    if (lhook) call dr_hook('radiation_mcica_lw:solver_mcica_lw',1,hook_handle)
    
  end subroutine solver_mcica_lw

end module radiation_mcica_lw
