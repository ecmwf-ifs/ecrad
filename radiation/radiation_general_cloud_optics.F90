! radiation_general_cloud_optics.F90 - Computing generalized cloud optical properties
!
! Copyright (C) 2020 ECMWF
!
! Author:  Robin Hogan
! Email:   r.j.hogan@ecmwf.int
! License: see the COPYING file for details
!

module radiation_general_cloud_optics

  implicit none

  public
  
contains

  ! Provides elemental function "delta_eddington_scat_od"
#include "radiation_delta_eddington.h"


  !---------------------------------------------------------------------
  ! Load cloud scattering data; this subroutine delegates to one
  ! in radiation_general_cloud_optics_data.F90

  subroutine setup_general_cloud_optics(config)

    use parkind1,         only : jprb
    use yomhook,          only : lhook, dr_hook

    use radiation_io,     only : nulout, nulerr, radiation_abort
    use radiation_config, only : config_type, NMaxCloudTypes

    real(jprb), parameter :: SolarReferenceTemperature = 5777.0_jprb ! K
    real(jprb), parameter :: TerrestrialReferenceTemperature = 273.15_jprb ! K

    type(config_type), intent(inout) :: config

    character(len=511) :: file_name

    integer :: jtype ! loop index
    integer :: strlen

    real(jprb) :: hook_handle

    if (lhook) call dr_hook('radiation_general_cloud_optics:setup_general_cloud_optics',0,hook_handle)

    ! Count number of cloud types
    config%n_cloud_types = 0
    do jtype = 1,NMaxCloudTypes
      if (len_trim(config%cloud_type_name(jtype)) > 0) then
        config%n_cloud_types = jtype
      else
        exit
      end if
    end do

    if (config%n_cloud_types == 0) then
      config%cloud_type_name(1) = "mie_liquid"
      config%cloud_type_name(2) = "baum-general-habit-mixture_ice"
      config%n_cloud_types = 2
    end if

    ! Allocate structures
    if (config%do_sw) then
      allocate(config%cloud_optics_sw(config%n_cloud_types))
    end if

    if (config%do_lw) then
      allocate(config%cloud_optics_lw(config%n_cloud_types))
    end if

    ! Load cloud optics data
    do jtype = 1,config%n_cloud_types
      if (config%cloud_type_name(jtype)(1:1) == '/') then
        file_name = trim(config%cloud_type_name(jtype))
      else
        strlen = len_trim(config%cloud_type_name(jtype))
        if (config%cloud_type_name(jtype)(strlen-2:strlen) == ".nc") then
          file_name = trim(config%directory_name) &
               &  // '/' // trim(config%cloud_type_name(jtype))
        else
          file_name = trim(config%directory_name) &
               &  // '/' // trim(config%cloud_type_name(jtype)) &
               &  // '_scattering.nc'
        end if
      end if

      if (config%do_sw) then
        if (config%iverbosesetup >= 2) then
          write(nulout,'(a,i0,a)') 'Shortwave cloud type ', jtype, ':'
        end if
        call config%cloud_optics_sw(jtype)%setup(file_name, &
             &  config%gas_optics_sw%spectral_def, &
             &  use_thick_averaging=.true., &
             &  weighting_temperature=SolarReferenceTemperature, &
             &  iverbose=config%iverbosesetup)
      end if

      if (config%do_lw) then
        if (config%iverbosesetup >= 2) then
          write(nulout,'(a,i0,a)') 'Longwave cloud type ', jtype, ':'
        end if
        call config%cloud_optics_lw(jtype)%setup(file_name, &
             &  config%gas_optics_lw%spectral_def, &
             &  use_thick_averaging=.true., &
             &  weighting_temperature=TerrestrialReferenceTemperature, &
             &  iverbose=config%iverbosesetup)
      end if

    end do

    if (lhook) call dr_hook('radiation_general_cloud_optics:setup_general_cloud_optics',1,hook_handle)

  end subroutine setup_general_cloud_optics

  !---------------------------------------------------------------------
  ! Compute cloud optical properties
  subroutine general_cloud_optics(nlev,istartcol,iendcol, &
       &  config, thermodynamics, cloud, & 
       &  od_lw_cloud, ssa_lw_cloud, g_lw_cloud, &
       &  od_sw_cloud, ssa_sw_cloud, g_sw_cloud)

    use parkind1, only           : jprb
    use yomhook,  only           : lhook, dr_hook

    use radiation_io,     only : nulout, nulerr, radiation_abort
    use radiation_config, only : config_type
    use radiation_thermodynamics, only    : thermodynamics_type
    use radiation_cloud, only             : cloud_type
    use radiation_constants, only         : AccelDueToGravity
    !use radiation_general_cloud_optics_data, only : general_cloud_optics_type

    integer, intent(in) :: nlev               ! number of model levels
    integer, intent(in) :: istartcol, iendcol ! range of columns to process
    type(config_type), intent(in), target :: config
    type(thermodynamics_type),intent(in)  :: thermodynamics
    type(cloud_type),   intent(in)        :: cloud

    ! Layer optical depth, single scattering albedo and g factor of
    ! clouds in each longwave band, where the latter two
    ! variables are only defined if cloud longwave scattering is
    ! enabled (otherwise both are treated as zero).
    real(jprb), dimension(config%n_bands_lw,nlev,istartcol:iendcol), intent(out) :: &
         &   od_lw_cloud
    real(jprb), dimension(config%n_bands_lw_if_scattering,nlev,istartcol:iendcol), &
         &   intent(out) :: ssa_lw_cloud, g_lw_cloud

    ! Layer optical depth, single scattering albedo and g factor of
    ! clouds in each shortwave band
    real(jprb), dimension(config%n_bands_sw,nlev,istartcol:iendcol), intent(out) :: &
         &   od_sw_cloud, ssa_sw_cloud, g_sw_cloud

    integer :: jtype

    real(jprb) :: hook_handle

    if (lhook) call dr_hook('radiation_general_cloud_optics:general_cloud_optics',0,hook_handle)

    if (config%iverbose >= 2) then
      write(nulout,'(a)') 'Computing cloud absorption/scattering properties'
    end if

    ! Array-wise assignment
    od_lw_cloud  = 0.0_jprb
    od_sw_cloud  = 0.0_jprb
    ssa_sw_cloud = 0.0_jprb
    g_sw_cloud   = 0.0_jprb
    if (config%do_lw_cloud_scattering) then
      ssa_lw_cloud = 0.0_jprb
      g_lw_cloud   = 0.0_jprb
    end if

    if (config%do_lw) then
      ! For the moment, we use od_lw_cloud, ssa_lw_cloud and
      ! g_lw_cloud as containers for mass extinction coefficient, mass
      ! scattering coefficient and mass scattering coefficient x
      ! asymmetry factor, then scale after
      do jtype = 1,config%n_cloud_types
        if (config%do_lw_cloud_scattering) then
          call config%cloud_optics_lw(jtype)%add_optical_properties(config%n_bands_lw, nlev, &
               &  iendcol+1-istartcol, cloud%mixing_ratio(istartcol:iendcol,:,jtype), &
               &  cloud%effective_radius(istartcol:iendcol,:,jtype), od_lw_cloud, ssa_lw_cloud, g_lw_cloud)
        else
          call config%cloud_optics_lw(jtype)%add_optical_properties(config%n_bands_lw, nlev, &
               &  iendcol+1-istartcol, cloud%mixing_ratio(istartcol:iendcol,:,jtype), &
               &  cloud%effective_radius(istartcol:iendcol,:,jtype), od_lw_cloud)
        end if
      end do

      if (config%do_lw_cloud_scattering) then
        ! Note that original cloud optics does not do delta-Eddington
        ! scaling for liquid clouds in longwave
        call delta_eddington_extensive(od_lw_cloud, ssa_lw_cloud, g_lw_cloud)

        ! Scale to get asymmetry factor and single scattering albedo
        g_lw_cloud   = g_lw_cloud   / max(ssa_lw_cloud, 1.0e-15_jprb)
        ssa_lw_cloud = ssa_lw_cloud / max(od_lw_cloud,  1.0e-15_jprb)
      end if

      ! Scale mass extinction coefficient to obtain optical depth
      if (config%is_homogeneous) then
        od_lw_cloud = od_lw_cloud &
             &  * spread(transpose(thermodynamics%pressure_hl(istartcol:iendcol, 2:nlev+1) &
             &             -thermodynamics%pressure_hl(istartcol:iendcol, 1:nlev)) &
             &  * (1.0_jprb / AccelDueToGravity),1,config%n_bands_lw)
      else
        od_lw_cloud = od_lw_cloud &
             &  * spread(transpose((thermodynamics%pressure_hl(istartcol:iendcol, 2:nlev+1) &
             &              -thermodynamics%pressure_hl(istartcol:iendcol, 1:nlev)) &
             &    / max(config%cloud_fraction_threshold, cloud%fraction(istartcol:iendcol,:))) &
             &  * (1.0_jprb / AccelDueToGravity),1,config%n_bands_lw)
      end if
      
    end if

    if (config%do_sw) then
      ! For the moment, we use od_sw_cloud, ssa_sw_cloud and
      ! g_sw_cloud as containers for mass extinction coefficient, mass
      ! scattering coefficient and mass scattering coefficient x
      ! asymmetry factor, then scale after
      do jtype =1,config%n_cloud_types
        call config%cloud_optics_sw(jtype)%add_optical_properties(config%n_bands_sw, nlev, &
             &  iendcol+1-istartcol, cloud%mixing_ratio(istartcol:iendcol,:,jtype), &
             &  cloud%effective_radius(istartcol:iendcol,:,jtype), od_sw_cloud, ssa_sw_cloud, g_sw_cloud)
      end do

      if (.not. config%do_sw_delta_scaling_with_gases) then
        call delta_eddington_extensive(od_sw_cloud, ssa_sw_cloud, g_sw_cloud)
      end if

      ! Scale to get asymmetry factor and single scattering albedo
      g_sw_cloud   = g_sw_cloud   / max(ssa_sw_cloud, 1.0e-15_jprb)
      ssa_sw_cloud = ssa_sw_cloud / max(od_sw_cloud,  1.0e-15_jprb)

      ! Scale mass extinction coefficient to obtain optical depth
      if (config%is_homogeneous) then
        od_sw_cloud = od_sw_cloud &
             &  * spread(transpose(thermodynamics%pressure_hl(istartcol:iendcol, 2:nlev+1) &
             &             -thermodynamics%pressure_hl(istartcol:iendcol, 1:nlev)) &
             &  * (1.0_jprb / AccelDueToGravity),1,config%n_bands_sw)
      else
        od_sw_cloud = od_sw_cloud &
             &  * spread(transpose((thermodynamics%pressure_hl(istartcol:iendcol, 2:nlev+1) &
             &              -thermodynamics%pressure_hl(istartcol:iendcol, 1:nlev)) &
             &    / max(config%cloud_fraction_threshold, cloud%fraction(istartcol:iendcol,:))) &
             &  * (1.0_jprb / AccelDueToGravity),1,config%n_bands_sw)
      end if
      
    end if

    if (lhook) call dr_hook('radiation_general_cloud_optics:general_cloud_optics',1,hook_handle)

  end subroutine general_cloud_optics

end module radiation_general_cloud_optics
