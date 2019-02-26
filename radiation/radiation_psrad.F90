! radiation_psrad.F90 - Interface to Pincus-Stevens RRTMG code
!
! Copyright (C) 2014-2015 ECMWF
!
! Author:  Robin Hogan
! Email:   r.j.hogan@ecmwf.int
! License: see the COPYING file for details
!

module radiation_psrad

  implicit none

contains

  ! Call Pincus-Stevens radiation code, which currently has spectral
  ! sampling disabled so is a McICA implementation using RRTMG for gas
  ! absorption and homogeneous clouds
  subroutine psrad(ncol, nlev, istartcol, iendcol, config, &
       &  single_level, thermodynamics, gas, cloud, aerosol, flux)

    use parkind1, only : jprb

    use radiation_config,         only : config_type
    use radiation_single_level,   only : single_level_type
    use radiation_thermodynamics, only : thermodynamics_type
    use radiation_gas,            only : gas_type, &
         &  IMassMixingRatio, IVolumeMixingRatio, NMaxGases, &
         &  IH2O, ICO2, IO3, IN2O, ICO, ICH4, IO2, ICFC11, ICFC12
    use radiation_cloud,          only : cloud_type
    use radiation_aerosol,        only : aerosol_type
    use radiation_flux,           only : flux_type

    use mo_radiation, only             : setup_radiation, rrtm_interface

    ! Inputs
    integer, intent(in) :: ncol       ! number of columns
    integer, intent(in) :: nlev       ! number of model levels
    integer, intent(in) :: istartcol, iendcol    ! range of columns to process
    type(config_type),        intent(in) :: config
    type(single_level_type),  intent(in) :: single_level
    type(thermodynamics_type),intent(in) :: thermodynamics
    type(gas_type),           intent(inout) :: gas
    type(cloud_type),         intent(in) :: cloud
    type(aerosol_type),       intent(in) :: aerosol
    ! Output
    type(flux_type),          intent(inout):: flux

    ! Pressure at half levels in hPa
    real(jprb), dimension(ncol, nlev+1) :: pressure_hl_hPa

    ! Surface pressure
    real(jprb), dimension(ncol) :: pressure_surf_hPa

    ! Temperature and pressure on full levels
    real(jprb), dimension(ncol, nlev) :: pressure_fl_hPa, temperature_fl

    ! Effective radius in microns
    real(jprb), dimension(ncol, nlev) :: re_liq_um, re_ice_um

    if (thermodynamics%pressure_hl(1,2) < thermodynamics%pressure_hl(1,1)) then
      ! Arrays are arranged in order of decreasing pressure /
      ! increasing height: need to reverse them
      call psrad_reverse(ncol, nlev,istartcol,iendcol, config, &
           &  single_level, thermodynamics, gas, cloud, aerosol, flux)
    else
      ! Arrays arranged in order of increasing pressure / decreasing
      ! height: proceed normally
      call setup_radiation(trim(config%directory_name))

      ! Convert to hPa
      pressure_hl_hPa(istartcol:iendcol,:) &
           &  = thermodynamics%pressure_hl(istartcol:iendcol,:)*0.01
      pressure_surf_hPa(istartcol:iendcol) &
           &  = pressure_hl_hPa(istartcol:iendcol,nlev+1)

      ! Pressure and temperature on full levels are just taken to be
      ! the average of the values on half levels
      pressure_fl_hPa(istartcol:iendcol,1:nlev) = 0.5_jprb &
           &  * (pressure_hl_hPa(istartcol:iendcol,1:nlev) &
           &  + pressure_hl_hPa(istartcol:iendcol,2:nlev+1))
      temperature_fl(istartcol:iendcol,1:nlev) = 0.5_jprb &
           &  * (thermodynamics%temperature_hl(istartcol:iendcol,1:nlev) &
           &  + thermodynamics%temperature_hl(istartcol:iendcol,2:nlev+1))

      ! Convert effective radius from metres to microns
      re_liq_um = cloud%re_liq * 1.0e6_jprb
      re_ice_um = cloud%re_ice * 1.0e6_jprb

      call gas%set_units(IMassMixingRatio,   igas=IH2O)
      call gas%set_units(IMassMixingRatio,   igas=IO3)
      call gas%set_units(IVolumeMixingRatio, igas=ICO2)
      call gas%set_units(IVolumeMixingRatio, igas=ICH4)
      call gas%set_units(IVolumeMixingRatio, igas=IN2O)
      call gas%set_units(IVolumeMixingRatio, igas=ICFC11)
      call gas%set_units(IVolumeMixingRatio, igas=ICFC12)
      call gas%set_units(IVolumeMixingRatio, igas=IO2)

      ! Call Pincus-Stevens code
      call rrtm_interface(ncol, iendcol, nlev, &
           &  single_level%lw_emissivity(:,1), single_level%sw_albedo(:,1), &
           &  single_level%cos_sza, pressure_fl_hPa, &
           &  pressure_hl_hPa, pressure_surf_hPa, &
           &  temperature_fl, thermodynamics%temperature_hl, &
           &  single_level%skin_temperature, &
           &  gas%mixing_ratio(:,:,IH2O), gas%mixing_ratio(:,:,IH2O), &
           &  cloud%q_liq, cloud%q_ice, re_liq_um, re_ice_um, &
           &  cloud%fraction, gas%mixing_ratio(:,:,IO3), &
           &  gas%mixing_ratio(:,:,ICO2), gas%mixing_ratio(:,:,ICH4), &
           &  gas%mixing_ratio(:,:,IN2O), gas%mixing_ratio(:,:,ICFC11), &
           &  gas%mixing_ratio(:,:,ICFC12), gas%mixing_ratio(:,:,IO2), &
           &  flux%lw_dn, flux%lw_up, flux%sw_dn, flux%sw_up, &
           &  flux%lw_dn_clear, flux%lw_up_clear, &
           &  flux%sw_dn_clear, flux%sw_up_clear)
    end if

    ! Direct fluxes not set by PSrad so need to make sure they are
    ! zero
    if (allocated(flux%sw_dn_direct)) then
      flux%sw_dn_direct = 0.0_jprb
    end if
    if (allocated(flux%sw_dn_direct_clear)) then
      flux%sw_dn_direct_clear = 0.0_jprb
    end if

  end subroutine psrad


  ! This is almost entirely a copy-and-paste from radiation_reverse
  subroutine psrad_reverse(ncol, nlev,istartcol,iendcol, config, &
       &  single_level, thermodynamics, gas, cloud, aerosol, flux)

    use parkind1, only : jprb

    use radiation_config,         only : config_type
    use radiation_single_level,   only : single_level_type
    use radiation_thermodynamics, only : thermodynamics_type
    use radiation_gas,            only : gas_type, IO2
    use radiation_cloud,          only : cloud_type
    use radiation_aerosol,        only : aerosol_type
    use radiation_flux,           only : flux_type

    implicit none

    ! Inputs
    integer, intent(in) :: ncol               ! number of columns
    integer, intent(in) :: nlev               ! number of model levels
    integer, intent(in) :: istartcol, iendcol ! range of columns to process
    type(config_type),        intent(in) :: config
    type(single_level_type),  intent(in) :: single_level
    type(thermodynamics_type),intent(in) :: thermodynamics
    type(gas_type),           intent(in) :: gas
    type(cloud_type),         intent(in) :: cloud
    type(aerosol_type),       intent(in) :: aerosol
    ! Output
    type(flux_type),          intent(inout):: flux

    ! Reversed data structures
    type(thermodynamics_type) :: thermodynamics_rev
    type(gas_type)            :: gas_rev
    type(cloud_type)          :: cloud_rev
    type(aerosol_type)        :: aerosol_rev
    type(flux_type)           :: flux_rev

    ! Allocate reversed arrays
    call thermodynamics_rev%allocate(ncol, nlev)
    call cloud_rev%allocate(ncol, nlev)
    call flux_rev%allocate(config, istartcol, iendcol, nlev)

    ! Fill reversed thermodynamic arrays
    thermodynamics_rev%pressure_hl(istartcol:iendcol,:) &
         &  = thermodynamics%pressure_hl(istartcol:iendcol, nlev+1:1:-1)
    thermodynamics_rev%temperature_hl(istartcol:iendcol,:) &
         &  = thermodynamics%temperature_hl(istartcol:iendcol, nlev+1:1:-1)

    ! Fill reversed gas arrays
    call gas%reverse(istartcol, iendcol, gas_rev)

    ! Fill reversed cloud arrays
    cloud_rev%q_liq(istartcol:iendcol,:) &
         &  = cloud%q_liq(istartcol:iendcol,nlev:1:-1)
    cloud_rev%re_liq(istartcol:iendcol,:) &
         &  = cloud%re_liq(istartcol:iendcol,nlev:1:-1)
    cloud_rev%q_ice(istartcol:iendcol,:) &
         &  = cloud%q_ice(istartcol:iendcol,nlev:1:-1)
    cloud_rev%re_ice(istartcol:iendcol,:) &
         &  = cloud%re_ice(istartcol:iendcol,nlev:1:-1)
    cloud_rev%fraction(istartcol:iendcol,:) &
         &  = cloud%fraction(istartcol:iendcol,nlev:1:-1)
    cloud_rev%overlap_param(istartcol:iendcol,:) &
         &  = cloud%overlap_param(istartcol:iendcol,nlev-1:1:-1)
    if (allocated(cloud%fractional_std)) then
      cloud_rev%fractional_std(istartcol:iendcol,:) &
           &  = cloud%fractional_std(istartcol:iendcol,nlev:1:-1)
    else
      cloud_rev%fractional_std(istartcol:iendcol,:) = 0.0_jprb       
    end if
    if (allocated(cloud%inv_cloud_effective_size)) then
      cloud_rev%inv_cloud_effective_size(istartcol:iendcol,:) &
           &  = cloud%inv_cloud_effective_size(istartcol:iendcol,nlev:1:-1)
    else
      cloud_rev%inv_cloud_effective_size(istartcol:iendcol,:) = 0.0_jprb       
    end if

    ! Run radiation scheme on reversed profiles
    call psrad(ncol, nlev, istartcol, iendcol, &
         &  config, single_level, thermodynamics_rev, gas_rev, &
         &  cloud_rev, aerosol_rev, flux_rev)

    ! Reorder fluxes
    if (allocated(flux%lw_up)) then
      flux%lw_up(istartcol:iendcol,:) &
           &  = flux_rev%lw_up(istartcol:iendcol,nlev+1:1:-1)
      flux%lw_dn(istartcol:iendcol,:) &
           &  = flux_rev%lw_dn(istartcol:iendcol,nlev+1:1:-1)
      if (allocated(flux%lw_up_clear)) then
        flux%lw_up_clear(istartcol:iendcol,:) &
             &  = flux_rev%lw_up_clear(istartcol:iendcol,nlev+1:1:-1)
        flux%lw_dn_clear(istartcol:iendcol,:) &
             &  = flux_rev%lw_dn_clear(istartcol:iendcol,nlev+1:1:-1)
      end if
    end if
    if (allocated(flux%sw_up)) then
      flux%sw_up(istartcol:iendcol,:) &
           &  = flux_rev%sw_up(istartcol:iendcol,nlev+1:1:-1)
      flux%sw_dn(istartcol:iendcol,:) &
           &  = flux_rev%sw_dn(istartcol:iendcol,nlev+1:1:-1)
      if (allocated(flux%sw_dn_direct)) then
        flux%sw_dn_direct(istartcol:iendcol,:) &
             &  = flux_rev%sw_dn_direct(istartcol:iendcol,nlev+1:1:-1)
      end if
      if (allocated(flux%sw_up_clear)) then
        flux%sw_up_clear(istartcol:iendcol,:) &
             &  = flux_rev%sw_up_clear(istartcol:iendcol,nlev+1:1:-1)
        flux%sw_dn_clear(istartcol:iendcol,:) &
             &  = flux_rev%sw_dn_clear(istartcol:iendcol,nlev+1:1:-1)
        if (allocated(flux%sw_dn_direct_clear)) then
          flux%sw_dn_direct_clear(istartcol:iendcol,:) &
               &  = flux_rev%sw_dn_direct_clear(istartcol:iendcol,nlev+1:1:-1)
        end if
      end if
    end if

    ! Deallocate reversed arrays
    call thermodynamics_rev%deallocate
    call gas_rev%deallocate
    call cloud_rev%deallocate
    call flux_rev%deallocate

  end subroutine psrad_reverse

end module radiation_psrad
