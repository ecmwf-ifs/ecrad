! ecrad_driver.F90 - Driver for offline ECRAD radiation scheme
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

program ecrad_driver

  ! --------------------------------------------------------
  ! Section 1: Declarations
  ! --------------------------------------------------------

  use radiation_config,         only : config_type
  use radiation_single_level,   only : single_level_type
  use radiation_thermodynamics, only : thermodynamics_type
  use radiation_gas,            only : gas_type
  use radiation_cloud,          only : cloud_type
  use radiation_aerosol,        only : aerosol_type
  use ecrad_driver_config,      only : driver_config_type
  use radiation_flux,           only : flux_type

  use ecrad_standalone, only: ecrad_standalone_setup, ecrad_standalone_run, &
                            & ecrad_standalone_save_output

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
  type(flux_type)           :: flux

  integer :: ncol, nlev         ! Number of columns and levels

  ! Name of file names specified on command line
  character(len=512) :: nml_file_name, input_file_name, output_file_name
  integer            :: istatus ! Result of get_command_argument


  ! --------------------------------------------------------
  ! Section 2: Configure
  ! --------------------------------------------------------

  ! Check program called with correct number of arguments
  if (command_argument_count() < 3) then
    stop 'Usage: ecrad config.nam input_file.nc output_file.nc'
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

  ! Run the standalone radiation scheme
  call ecrad_standalone_setup( &
    & nml_file_name, input_file_name, driver_config, config, &
    & single_level, thermodynamics, gas, cloud, aerosol, &
    & ncol, nlev )

  ! --------------------------------------------------------
  ! Section 4: Call radiation scheme
  ! --------------------------------------------------------

  call ecrad_standalone_run( &
    & ncol, nlev, driver_config, config, &
    & single_level, thermodynamics, gas, cloud, aerosol, &
    & flux )

  ! --------------------------------------------------------
  ! Section 5: Check and save output
  ! --------------------------------------------------------

  call ecrad_standalone_save_output( &
    & output_file_name, driver_config, config, &
    & thermodynamics, flux )

end program ecrad_driver
