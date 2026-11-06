! (C) Copyright 2026- ECMWF.
! Licensed under the Apache Licence Version 2.0; see COPYING.

program missing_gases
  use parkind1, only : jprb, jpim
  use easy_netcdf, only : netcdf_file
  use radiation_io, only : radiation_abort
  use ecrad_driver_config, only : driver_config_type
  use ecrad_driver_read_input, only : read_input
  use radiation_config, only : config_type
  use radiation_single_level, only : single_level_type
  use radiation_thermodynamics, only : thermodynamics_type
  use radiation_gas, only : gas_type, IH2O, ICO, INO2
  use radiation_cloud, only : cloud_type
  use radiation_aerosol, only : aerosol_type
  implicit none

  type(netcdf_file) :: file
  type(config_type) :: config
  type(driver_config_type) :: driver_config
  type(single_level_type) :: single_level
  type(thermodynamics_type) :: thermodynamics
  type(gas_type) :: gas
  type(cloud_type) :: cloud
  type(aerosol_type) :: aerosol
  character(len=1024) :: config_name, input_name
  integer(kind=jpim) :: ncol, nlev

  call get_command_argument(1, config_name)
  call get_command_argument(2, input_name)
  call config%read(trim(config_name))
  call driver_config%read(trim(config_name))
  call file%open(trim(input_name), iverbose=0)
  call file%transpose_matrices(.true.)
  call read_input(file, config, driver_config, ncol, nlev, &
       &          single_level, thermodynamics, gas, cloud, aerosol)
  call file%close()

  if (gas%is_present(ICO) .or. gas%is_present(INO2)) call radiation_abort('Unexpected gas in input')
  if (any(gas%mixing_ratio(:,:,ICO) /= 0.0_jprb)) call radiation_abort('Missing CO is not zero')
  if (any(gas%mixing_ratio(:,:,INO2) /= 0.0_jprb)) call radiation_abort('Missing NO2 is not zero')
  if (.not. gas%is_present(IH2O)) call radiation_abort('Missing water vapour in input')
  if (.not. any(gas%mixing_ratio(:,:,IH2O) > 0.0_jprb)) call radiation_abort('Water vapour is not positive')
end program missing_gases
