! ecrad_ifs_driver.F90 - Configure and run the IFS-style configuration of ECRAD radiation scheme
!
! (C) Copyright 2014- ECMWF.
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
! ECRAD is the radiation scheme used in the ECMWF Integrated
! Forecasting System in cycle 43R3 and later. Several solvers are
! available, including McICA, Tripleclouds and SPARTACUS (the Speedy
! Algorithm for Radiative Transfer through Cloud Sides, a modification
! of the two-stream formulation of shortwave and longwave radiative
! transfer to account for 3D radiative effects). Gas optical
! properties are provided by the RRTM-G gas optics scheme.

! This program takes three arguments:
! 1) Namelist file to configure the radiation calculation
! 2) Name of a NetCDF file containing one or more atmospheric profiles
! 3) Name of output NetCDF file

program ecrad_ifs_driver

  ! --------------------------------------------------------
  ! Section 1: Declarations
  ! --------------------------------------------------------
  use parkind1,                 only : jprb ! Working precision

  use radiation_config,         only : config_type
  use radiation_single_level,   only : single_level_type
  use radiation_thermodynamics, only : thermodynamics_type
  use radiation_gas,            only : gas_type
  use radiation_cloud,          only : cloud_type
  use radiation_aerosol,        only : aerosol_type
  use radiation_flux,           only : flux_type
  use ecrad_driver_config,      only : driver_config_type
  use type_model,               only : model

  use ecrad_standalone,         only : ecrad_standalone_setup, &
       &   ecrad_standalone_save_output, ecrad_standalone_run
  use ecrad_ifs,                only : ifs_config_type, ecrad_ifs_setup, &
       &   ecrad_ifs_interpolate_in, ecrad_ifs_run, ecrad_ifs_interpolate_out, &
       &   ecrad_ifs_validate

  implicit none

  ! Derived types for the inputs to the radiation scheme
  type(config_type)         :: config
  type(single_level_type)   :: single_level
  type(thermodynamics_type) :: thermodynamics
  type(gas_type)            :: gas
  type(cloud_type)          :: cloud
  type(aerosol_type)        :: aerosol

  ! Configuration specific to this driver
  type(driver_config_type)  :: driver_config

  ! Derived type containing outputs from the radiation scheme
  type(flux_type)           :: flux_out
  type(flux_type)           :: flux_ref

  ! Model type and IFS config values
  type(model)               :: ydmodel
  type(ifs_config_type)     :: ifs_config

  integer :: ncol, nlev         ! Number of columns and levels
  !integer :: istartcol, iendcol ! Range of columns to process

  ! Name of file names specified on command line
  character(len=512) :: input_file_name, output_file_name, nml_file_name
  integer            :: istatus ! Result of command_argument_count

  ! monolithic IFS data structure to pass to radiation scheme
  real(kind=jprb), allocatable :: zrgp(:,:,:)

  ! Seed for random number generator
  integer, allocatable :: iseed(:,:)

  ! Empty thermodynamics type to store pressure_hl for output at the end
  type(thermodynamics_type)    :: thermodynamics_out

  ! Per-block flux data structure to validate outputs
  type(flux_type), allocatable :: flux_blocks(:)

  logical :: is_validating


  ! --------------------------------------------------------
  ! Section 2: Configure
  ! --------------------------------------------------------

  ! Check program called with correct number of arguments
  if (command_argument_count() < 3) then
    stop 'Usage: ecrad config.nam input_file.nc output_file.nc [output_surface_file.nc]'
  end if

  ! Use namelist to configure the radiation calculation
  call get_command_argument(1, nml_file_name, status=istatus)
  if (istatus /= 0) then
    stop 'Failed to read name of namelist file as string of length < 512'
  end if

  ! Get NetCDF input file name
  call get_command_argument(2, input_file_name, status=istatus)
  if (istatus /= 0) then
    stop 'Failed to read name of input NetCDF file as string of length < 512'
  end if

  ! Get NetCDF output file name
  call get_command_argument(3, output_file_name, status=istatus)
  if (istatus /= 0) then
    stop 'Failed to read name of output NetCDF file as string of length < 512'
  end if

  ! --------------------------------------------------------
  ! Section 3: Read input data file
  ! --------------------------------------------------------

  call ecrad_standalone_setup( &
    & nml_file_name, input_file_name, driver_config, config, &
    & single_level, thermodynamics, gas, cloud, aerosol, &
    & ncol, nlev )

  call ecrad_ifs_setup(nml_file_name, driver_config, config, ydmodel, ncol)

  call ecrad_ifs_interpolate_in( &
   & driver_config, ifs_config, ncol, nlev, &
   & single_level, thermodynamics, gas, cloud, aerosol, &
   & ydmodel, zrgp, iseed, thermodynamics_out )

  ! --------------------------------------------------------
  ! Section 4: Call radiation scheme
  ! --------------------------------------------------------

  ! Call standalone scheme to generate reference results
  call ecrad_standalone_run( &
   & ncol, nlev, driver_config, config, &
   & single_level, thermodynamics, gas, cloud, aerosol, flux_ref )

  ! Deallocate input data structures
  call single_level%deallocate
  call thermodynamics%deallocate
  call gas%deallocate
  call cloud%deallocate
  call aerosol%deallocate

  call ecrad_ifs_run( &
   & driver_config, ifs_config, ncol, nlev, &
   & ydmodel, iseed, zrgp, flux_blocks )

  call ecrad_ifs_interpolate_out( &
   & driver_config, ncol, nlev, ydmodel, flux_blocks, flux_out )

  call ecrad_standalone_save_output( &
    & output_file_name, driver_config, ydmodel%yrml_phy_rad%yradiation%rad_config, &
    & thermodynamics_out, flux_out )

  is_validating = ecrad_ifs_validate(config, flux_ref, flux_out, ncol, nlev)

  call flux_out%deallocate
  deallocate(thermodynamics_out%pressure_hl)

  if (.not. is_validating) then
    error stop
  endif

end program ecrad_ifs_driver
