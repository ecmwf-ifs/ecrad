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
    use radiation_regions, only        : calc_region_properties
    use radiation_overlap, only        : calc_overlap_matrices
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

    ! Index to FLOTSAM band ID
    integer :: iband

    ! Loop counter
    integer :: jcol, jg

    real(jprb) :: hook_handle

    if (lhook) call dr_hook('radiation_flotsam_sw:solver_flotsam_sw',0,hook_handle)

    iband = flotsam_new_band_profile()

    do jcol = istartcol,iendcol

    flotsam_set_geometry(iband, single_level%cos_solar_zenith_angle(jcol), &
         &                      single_level%cos_sensor_zenith_angle(jcol), &
         &  single_level%solar_azimuth_angle(jcol) - single_level%sensor_zenith_angle(jcol))
    
    do jg = 1,config%n_band_sw


    end do

    flotsam_free_band_profile(iband)
    
    if (lhook) call dr_hook('radiation_flotsam_sw:solver_flotsam_sw',1,hook_handle)

  end subroutine radiance_solver_flotsam_sw

end module radiation_flotsam_sw
