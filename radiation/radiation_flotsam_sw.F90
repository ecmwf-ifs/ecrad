! radiation_flotsam_sw.F90 - Interface to FLOTSAM solar radiance model
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

module radiation_flotsam_sw

  public

contains

  ! This module contains the shortwave "FLOTSAM" solver

  subroutine radiance_solver_flotsam_sw(nlev,istartcol,iendcol, &
       &  config, single_level, cloud, & 
       &  od, ssa, g, od_cloud, ssa_cloud, g_cloud, &
       &  albedo_direct, albedo_diffuse, incoming_sw, &
       &  flux, use_stochastic_columns)

    use parkind1, only           : jprb
    use yomhook,  only           : lhook, dr_hook

    use radiation_config, only         : config_type, IPdfShapeGamma
    use radiation_single_level, only   : single_level_type
    use radiation_cloud, only          : cloud_type
    use radiation_cloud_cover, only    : cloud_cover
    use radiation_flux, only           : flux_type
    use radiation_optimal_columns, only: optimal_columns
    use radiation_cloud_generator, only: cloud_generator

    implicit none

#include <flotsam.inc>

    ! Inputs
    integer, intent(in) :: nlev               ! number of model levels
    integer, intent(in) :: istartcol, iendcol ! range of columns to process
    type(config_type),        intent(in) :: config
    type(single_level_type),  intent(in) :: single_level
    type(cloud_type),         intent(in) :: cloud

    ! Gas and aerosol optical depth, single-scattering albedo and
    ! asymmetry factor at each shortwave g-point
    real(jprb), intent(in), dimension(config%n_g_sw,nlev,istartcol:iendcol) :: &
         &  od, ssa, g

    ! Cloud and precipitation optical depth, single-scattering albedo and
    ! asymmetry factor in each shortwave band
    real(jprb), intent(in), dimension(config%n_bands_sw,nlev,istartcol:iendcol) :: &
         &  od_cloud, ssa_cloud, g_cloud

    ! Direct and diffuse surface albedos, and the incoming shortwave
    ! flux into a plane perpendicular to the incoming radiation at
    ! top-of-atmosphere in each of the shortwave g points
    real(jprb), intent(in), dimension(config%n_g_sw,istartcol:iendcol) :: &
         &  albedo_direct, albedo_diffuse, incoming_sw

    ! Do we use stochastic cloud generator or optimal columns?
    logical, optional, intent(in) :: use_stochastic_columns

    ! Output
    type(flux_type), intent(inout):: flux

    ! Local variables

    ! Absorption and scattering gas optical depths
    real(jprb), dimension(nlev,config%n_g_sw) :: &
         &  od_abs, od_rayleigh

    ! We compute total cloud cover, then assume clouds in each layer
    ! have a cloud fraction equal to the total cloud cover. This means
    ! rescaling the cloud optical depth
    real(jprb), dimension(nlev) :: od_cloud_scaled

    ! Fractional standard deviation of the total-column optical depth
    ! in the cloudy part of the gridbox
    real(jprb) :: cloudy_fsd_od(config%n_g_sw)

    integer :: ind(nlev)

    ! Solar incoming in each band
    real(jprb) :: weight(config%n_g_sw)

    real(jprb) :: albedo(FLOTSAM_NUM_ALBEDO_COMPONENTS)

    ! Weights and optical depths used by the optimal columns scheme
    !real(jprb) :: weight_sub(config%n_bands_sw,config%n_cloudy_subcolumns_sw)
    !real(jprb) :: od_cloud_sub(config%n_bands_sw,nlev,config%n_cloudy_subcolumns_sw)
    real(jprb), allocatable :: weight_sub_oc(:,:)
    real(jprb), allocatable :: od_cloud_sub_oc(:,:,:)

    ! Weight and optical depth scalings used by the stochastic columns scheme
    real(jprb), allocatable :: od_scaling_sub_sc(:,:) ! (config%n_cloudy_subcolumns_sw,nlev)
    real(jprb) :: od_cloud_sub_sc(nlev)
    real(jprb) :: weight_sub_sc
    
    ! Radiance from one subcolumn
    real(jprb) :: radiance_band_sub

    real(jprb) :: total_cloud_cover

    ! Phase function and components
    real(jprb) :: pf(nlev)
    real(jprb), allocatable :: pf_components(:,:)

    logical :: use_stochastic_columns_local

    ! Index to FLOTSAM band ID
    integer :: iband, istatus, ig

    ! Number of albedos
    integer :: nalbedo

    integer :: i

    ! Number of phase-function components
    integer :: n_pf_components

    ! Number of g points in band
    integer :: ng

    ! Loop counters
    integer :: jcol, jg, jband, jsub, jlev

    real(jprb) :: hook_handle

    if (lhook) call dr_hook('radiation_flotsam_sw:solver_flotsam_sw',0,hook_handle)

    ind = [(i,i=0,nlev-1)]

    iband = flotsam_new_band_profile()

    n_pf_components = flotsam_n_phase_function_components()
    allocate(pf_components(n_pf_components, nlev))

    pf = 1.0_jprb
    pf_components = 0.0_jprb
    pf_components(1,:) = 1.0_jprb

    ! Allocate relevant arrays depending on whether stochastic or
    ! optimal columns are being used
    if (present(use_stochastic_columns)) then
      use_stochastic_columns_local = use_stochastic_columns
    else
      use_stochastic_columns_local = .false.
    end if
    if (use_stochastic_columns_local) then
      allocate(od_scaling_sub_sc(config%n_cloudy_subcolumns_sw,nlev))
    else
      allocate(weight_sub_oc(config%n_bands_sw,config%n_cloudy_subcolumns_sw))
      allocate(od_cloud_sub_oc(config%n_bands_sw,nlev,config%n_cloudy_subcolumns_sw))
    end if

    do jcol = istartcol,iendcol

      if (single_level%cos_sensor_zenith_angle(jcol) > 0.0_jprb &
           &  .and. single_level%cos_sza(jcol) > 0.0_jprb) then

        istatus = flotsam_set_geometry(iband, single_level%cos_sza(jcol), &
             &                      single_level%cos_sensor_zenith_angle(jcol), &
             &  single_level%solar_azimuth_angle(jcol) - single_level%sensor_azimuth_angle(jcol))

        if (use_stochastic_columns_local) then
          call cloud_generator(config%n_cloudy_subcolumns_sw, nlev, config%i_overlap_scheme, &
               &  single_level%iseed(jcol), &
               &  config%cloud_fraction_threshold, &
               &  cloud%fraction(jcol,:), cloud%overlap_param(jcol,:), &
               &  config%cloud_inhom_decorr_scaling, cloud%fractional_std(jcol,:), &
               &  config%pdf_sampler, od_scaling_sub_sc, total_cloud_cover, &
               &  use_beta_overlap=config%use_beta_overlap, &
               &  use_vectorizable_generator=config%use_vectorizable_generator)

          ! Each cloudy subcolumn is weighted equally
          weight_sub_sc = total_cloud_cover / config%n_cloudy_subcolumns_sw
        else
          call optimal_columns(config%n_bands_sw, config%n_cloudy_subcolumns_sw, nlev, config%cloud_fraction_threshold, &
               &  cloud%fraction(jcol,:), cloud%overlap_param(jcol,:), cloud%fractional_std(jcol,:), &
               &  od_cloud(:,:,jcol), weight_sub_oc, od_cloud_sub_oc, total_cloud_cover, &
               &  cloudy_fsd_od=cloudy_fsd_od)
        end if

        flux%cloud_cover_sw(jcol) = total_cloud_cover

        do jband = 1,config%n_bands_sw
          
          ! Copy gas scattering properties
          ng = 0
          albedo = 0.0_jprb
          do jg = 1,config%n_g_sw
            if (config%i_band_from_reordered_g_sw(jg) == jband) then
              ng = ng + 1
              weight(ng)        = incoming_sw(jg,jcol)
              od_abs(:,ng)      = od(jg,:,jcol) * (1.0_jprb - ssa(jg,:,jcol))
              od_rayleigh(:,ng) = od(jg,:,jcol) * ssa(jg,:,jcol)
              ! Albedo is simply a weighted average
              albedo = albedo + albedo_diffuse(jg,jcol)*weight(ng)
            end if
          end do
          albedo = albedo / sum(weight(1:ng))
          weight(1:ng) = weight(1:ng) / sum(weight(1:ng))
          
          ! If we are over sea we optionally use the Cox-Munk model
          ! for ocean reflectance accounting for sun glint
          nalbedo = 1
          if (allocated(single_level%sea_fraction) &
               &  .and. allocated(config%i_ocean_reflectance_id_sw)) then
            if (config%i_ocean_reflectance_id_sw(jband) >= 0 &
                 &  .and. single_level%sea_fraction(jcol) > 0.5_jprb) then
              nalbedo = 4
              istatus = flotsam_get_ocean_albedo_components( &
                   &  config%i_ocean_reflectance_id_sw(jband), single_level%cos_sza(jcol), &
                   &  single_level%cos_sensor_zenith_angle(jcol), &
                   &  single_level%solar_azimuth_angle(jcol) - single_level%sensor_azimuth_angle(jcol), &
                   &  3.0_jprb, albedo)
            end if
          end if

          istatus = flotsam_init_band_profile_direct(iband, ng, nlev, weight, &
               &  od_rayleigh, od_abs)
          ! Clear sky
          istatus = flotsam_reflectance(iband, nalbedo, albedo, 0, ind, od_cloud(jband,:,jcol), &
               &  ssa_cloud(jband,:,jcol), pf, pf_components, flux%sw_radiance_clear_band(jband,jcol))

          if (flux%cloud_cover_sw(jcol) > 0.0_jprb) then
            ! Scale clear radiance down according to the cloud cover
            flux%sw_radiance_band(jband,jcol) &
                 &  = (1.0_jprb - flux%cloud_cover_sw(jcol)) * flux%sw_radiance_clear_band(jband,jcol)

            ! Sum over cloudy subcolumns
            do jsub = 1,config%n_cloudy_subcolumns_sw
              if (use_stochastic_columns_local) then
                od_cloud_sub_sc = od_cloud(jband,:,jcol) * od_scaling_sub_sc(jsub,:)
                istatus = flotsam_reflectance(iband, nalbedo, albedo, nlev, ind, od_cloud_sub_sc, &
                     &  ssa_cloud(jband,:,jcol), pf, pf_components, radiance_band_sub)
                flux%sw_radiance_band(jband,jcol) = flux%sw_radiance_band(jband,jcol) &
                     &  + weight_sub_sc*radiance_band_sub
              else
                istatus = flotsam_reflectance(iband, nalbedo, albedo, nlev, ind, od_cloud_sub_oc(jband,:,jsub), &
                     &  ssa_cloud(jband,:,jcol), pf, pf_components, radiance_band_sub)
                flux%sw_radiance_band(jband,jcol) = flux%sw_radiance_band(jband,jcol) &
                     &  + weight_sub_oc(jband,jsub)*radiance_band_sub
              end if
            end do

          else
            flux%sw_radiance_band(jband,jcol) = flux%sw_radiance_clear_band(jband,jcol)            
          end if
          
        end do

      else
        ! Either the sun or the sensor is below the horizon
        flux%sw_radiance_band(:,jcol) = 0.0_jprb
        flux%sw_radiance_clear_band(:,jcol) = 0.0_jprb
      end if

    end do

    istatus = flotsam_free_band_profile(iband)

    deallocate(pf_components)
    
    if (lhook) call dr_hook('radiation_flotsam_sw:solver_flotsam_sw',1,hook_handle)

  end subroutine radiance_solver_flotsam_sw

  subroutine allocate_ocean_reflectance_model(config)

    use parkind1, only           : jprb
    use yomhook,  only           : lhook, dr_hook
    use radiation_config, only   : config_type
    use radiation_io, only       : nulout
    
    implicit none

#include <flotsam.inc>

    type(config_type), intent(inout) :: config

    integer, parameter :: nwind = 10

    integer :: jband, jwind

    real(jprb) :: wind(nwind)

    real(jprb) :: wavelength

    real(jprb) :: hook_handle

    if (lhook) call dr_hook('radiation_flotsam_sw:allocate_ocean_reflectance_model',0,hook_handle)

    allocate(config%i_ocean_reflectance_id_sw(config%n_bands_sw))

    do jwind = 1,nwind
      wind(jwind) = jwind*2.0_jprb - 1.0_jprb
    end do

    do jband = 1,config%n_bands_sw

      wavelength = 0.01_jprb / (0.5_jprb*(config%gas_optics_sw%spectral_def%wavenumber1_band(jband) &
           &  + config%gas_optics_sw%spectral_def%wavenumber1_band(jband)))
      
      if (config%iverbosesetup >= 2) then
        write(nulout, '(a,f10.1,a)') 'Initializing ocean reflectance model at wavelength of ', &
             &  wavelength*1.0e9, ' nm'
      end if

      config%i_ocean_reflectance_id_sw(jband) &
           &  = flotsam_new_ocean_brdf(wavelength, nwind, wind);
      ! FIX: check ID is valid

    end do

    if (lhook) call dr_hook('radiation_flotsam_sw:allocate_ocean_reflectance_model',1,hook_handle)

  end subroutine allocate_ocean_reflectance_model

  pure function correlation(ns, v1, v2) 

    use parkind1, only : jprb

    integer, intent(in) :: ns
    real(jprb), intent(in) :: v1(ns), v2(ns)
    real(jprb) :: correlation

    real(jprb) :: mean1, mean2, cov, var1, var2

    mean1 = sum(v1) / ns
    mean2 = sum(v2) / ns

    cov = sum((v1-mean1)*(v2-mean2)) / ns
    var1 = sum((v1-mean1)**2) / ns
    var2 = sum((v2-mean2)**2) / ns
    correlation = cov / max(1.0e-12_jprb, sqrt(var1*var2))

  end function correlation

end module radiation_flotsam_sw
