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

  ! This module contains just one subroutine, the shortwave
  ! "FLOTSAM" solver

  subroutine radiance_solver_flotsam_sw(nlev,istartcol,iendcol, &
       &  config, single_level, cloud, & 
       &  od, ssa, g, od_cloud, ssa_cloud, g_cloud, &
       &  albedo_direct, albedo_diffuse, incoming_sw, &
       &  flux)

    use parkind1, only           : jprb
    use yomhook,  only           : lhook, dr_hook

    use radiation_config, only         : config_type, IPdfShapeGamma
    use radiation_single_level, only   : single_level_type
    use radiation_cloud, only          : cloud_type
    use radiation_cloud_cover, only    : cloud_cover
    use radiation_flux, only           : flux_type

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
!    real(jprb), intent(in), dimension(istartcol:iendcol,nlev,config%n_g_sw) :: &
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

    integer :: ind(nlev)

    real(jprb) :: weight(config%n_g_sw)

    real(jprb) :: albedo(1)

    ! Phase function and components
    real(jprb) :: pf(nlev)
    real(jprb), allocatable :: pf_components(:,:)

    ! Index to FLOTSAM band ID
    integer :: iband, istatus, ig

    integer :: i

    ! Number of phase-function components
    integer :: n_pf_components

    ! Number of g points in band
    integer :: ng

    ! Loop counter
    integer :: jcol, jg, jband

    real(jprb) :: hook_handle

    if (lhook) call dr_hook('radiation_flotsam_sw:solver_flotsam_sw',0,hook_handle)

    ind = [(i,i=0,nlev-1)]

    iband = flotsam_new_band_profile()

    n_pf_components = flotsam_n_phase_function_components()
    allocate(pf_components(n_pf_components, nlev))

    pf = 1.0_jprb
    pf_components = 0.0_jprb
    pf_components(1,:) = 1.0_jprb

    do jcol = istartcol,iendcol

      if (single_level%cos_sensor_zenith_angle(jcol) > 0.0_jprb &
           &  .and. single_level%cos_sza(jcol) > 0.0_jprb) then

        flux%cloud_cover_sw(jcol) = cloud_cover(nlev, config%i_overlap_scheme, &
             &  cloud%fraction(jcol,:), cloud%overlap_param(jcol,:), &
             &  config%use_beta_overlap)

        istatus = flotsam_set_geometry(iband, single_level%cos_sza(jcol), &
             &                      single_level%cos_sensor_zenith_angle(jcol), &
             &  single_level%solar_azimuth_angle(jcol) - single_level%sensor_azimuth_angle(jcol))
        
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
          
          istatus = flotsam_init_band_profile_direct(iband, ng, nlev, weight, &
               &  od_rayleigh, od_abs)
          ! Clear sky
          istatus = flotsam_reflectance(iband, 1, albedo, 0, ind, od_cloud(jband,:,jcol), &
               &  ssa_cloud(jband,:,jcol), pf, pf_components, flux%sw_radiance_clear_band(jband,jcol))

          if (flux%cloud_cover_sw(jcol) > 0.0_jprb) then
            ! Rescale cloud optical depth; the od_cloud variable
            od_cloud_scaled = cloud%fraction(jcol,:) * od_cloud(jband,:,jcol) &
                 &  * (1.0_jprb / flux%cloud_cover_sw(jcol))
            ! Cloudy sky
            istatus = flotsam_reflectance(iband, 1, albedo, nlev, ind, od_cloud_scaled, &
                 &  ssa_cloud(jband,:,jcol), pf, pf_components, flux%sw_radiance_band(jband,jcol))
            ! Weighted average between clear and cloudy
            flux%sw_radiance_band(jband,jcol) &
                 &  = flux%cloud_cover_sw(jcol) * flux%sw_radiance_band(jband,jcol) &
                 &  + (1.0_jprb-flux%cloud_cover_sw(jcol)) * flux%sw_radiance_clear_band(jband,jcol)
          else
            flux%sw_radiance_band(jband,jcol) = flux%sw_radiance_clear_band(jband,jcol)            
          end if
          
        end do
      
      else
        ! Either the sun or the sensor is below he horizon
        flux%sw_radiance_band(:,jcol) = 0.0_jprb
        flux%sw_radiance_clear_band(:,jcol) = 0.0_jprb
      end if

    end do

    istatus = flotsam_free_band_profile(iband)

    deallocate(pf_components)
    
    if (lhook) call dr_hook('radiation_flotsam_sw:solver_flotsam_sw',1,hook_handle)

  end subroutine radiance_solver_flotsam_sw

end module radiation_flotsam_sw
