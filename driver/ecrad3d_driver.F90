! ecrad3d_driver.F90 - Driver for 3D version of offline ECRAD radiation scheme
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
! Forecasting System in cycle 43R3 and later. This version uses an
! efficient method to represent 3D radiation transport between
! columns.
!
! This program takes three arguments:
! 1) Namelist file to configure the radiation calculation
! 2) Name of a NetCDF file containing one or more atmospheric profiles
! 3) Name of output NetCDF file

program ecrad3d_driver

  ! --------------------------------------------------------
  ! Section 1: Declarations
  ! --------------------------------------------------------
  use parkind1,                 only : jprb, jprd ! Working/double precision

  use radiation_io,             only : nulout
  use radiation_interface,      only : setup_radiation, radiation, set_gas_units
  use radiation_config,         only : config_type, ISolverPOMART3D, ISolverPOMART3DTICA
  use radiation_single_level,   only : single_level_type
  use radiation_thermodynamics, only : thermodynamics_type
  use radiation_gas,            only : gas_type, &
       &   IVolumeMixingRatio, IMassMixingRatio, &
       &   IH2O, ICO2, IO3, IN2O, ICO, ICH4, IO2, ICFC11, ICFC12, &
       &   IHCFC22, ICCl4, GasName, GasLowerCaseName, NMaxGases
  use radiation_cloud,          only : cloud_type
  use radiation_aerosol,        only : aerosol_type
  use ecrad3d_geometry,         only : geometry_type
  use radiation_flux,           only : flux_type
  use radiation_save,           only : save_fluxes, save_inputs, save_radiances
  use radiation_general_cloud_optics, only : save_general_cloud_optics
  use ecrad3d_save,             only : save_fluxes_3d    => save_fluxes, &
       &                               save_radiances_3d => save_radiances
  use ecrad_driver_config,      only : driver_config_type
  use ecrad3d_driver_read_input,only : read_input
  use ecrad3d_interface,        only : ecrad3d
  use easy_netcdf
  use print_matrix_mod,         only : print_matrix
  
  implicit none

  ! The NetCDF file containing the input profiles
  type(netcdf_file)         :: file

  ! Derived types for the inputs to the radiation scheme
  type(config_type)         :: config
  type(single_level_type)   :: single_level
  type(thermodynamics_type) :: thermodynamics
  type(gas_type)            :: gas
  type(cloud_type)          :: cloud
  type(aerosol_type)        :: aerosol
  type(geometry_type)       :: geometry

  ! Configuration specific to this driver
  type(driver_config_type)  :: driver_config

  ! Derived type containing outputs from the radiation scheme
  type(flux_type)           :: flux

  integer :: ncol, nlev         ! Number of columns and levels
  integer :: istartcol, iendcol ! Range of columns to process

  ! Name of file names specified on command line
  character(len=512) :: file_name
  integer            :: istatus ! Result of command_argument_count

  ! For parallel processing of multiple blocks
  integer :: jblock, nblock ! Block loop index and number

#ifndef NO_OPENMP
  ! OpenMP functions
  integer, external :: omp_get_thread_num
  real(kind=jprd), external :: omp_get_wtime
  ! Start/stop time in seconds
  real(kind=jprd) :: tstart, tstop
#endif

  ! Loop index for repeats (for benchmarking)
  integer :: jrepeat

  ! Are any variables out of bounds?
  logical :: is_out_of_bounds

  ! Do we use the explicit 3D solver, which requires a different
  ! interface?
  logical :: use_ecrad3d_solver
  
  ! --------------------------------------------------------
  ! Section 2: Configure
  ! --------------------------------------------------------

  ! Check program called with correct number of arguments
  if (command_argument_count() < 3) then
    stop 'Usage: ecrad3d config.nam input_file.nc output_file.nc'
  end if

  ! Use namelist to configure the radiation calculation
  call get_command_argument(1, file_name, status=istatus)
  if (istatus /= 0) then
    error stop 'Failed to read name of namelist file as string of length < 512'
  end if

  ! Read "radiation" namelist into radiation configuration type
  call config%read(file_name=file_name)

  ! Read "radiation_driver" namelist into radiation driver config type
  call driver_config%read(file_name)

  if (driver_config%iverbose >= 2) then
    write(nulout,'(a)') '------------------------- OFFLINE ECRAD3D RADIATION SCHEME -------------------------'
    write(nulout,'(a)') 'Copyright (C) 2014- ECMWF'
    write(nulout,'(a)') 'Contact: Robin Hogan (r.j.hogan@ecmwf.int)'
#ifdef PARKIND1_SINGLE
    write(nulout,'(a)') 'Floating-point precision: single'
#else
    write(nulout,'(a)') 'Floating-point precision: double'
#endif
    call config%print(driver_config%iverbose)
  end if

  ! If we use the explicit 3D solver then a different ecRad interface
  ! must be called
  if ((config%do_sw .and. (config%i_solver_sw == ISolverPOMART3D .or. config%i_solver_sw == ISolverPOMART3DTICA)) &
       & .or. (config%do_lw .and. config%i_solver_lw == ISolverPOMART3D)) then
    if (config%do_sw .and. config%do_lw .and. config%i_solver_sw /= config%i_solver_lw) then
      error stop 'The POMART3D solver must be selected for both SW and LW, or neither'
    end if
    use_ecrad3d_solver = .true.
  else
    use_ecrad3d_solver = .false.    
  end if
  
  ! Setup the radiation scheme: load the coefficients for gas and
  ! cloud optics
  call setup_radiation(config)

  if (driver_config%do_save_aerosol_optics) then
    call config%aerosol_optics%save('aerosol_optics.nc', iverbose=driver_config%iverbose)
  end if

  if (driver_config%do_save_cloud_optics .and. config%use_general_cloud_optics) then
    call save_general_cloud_optics(config, 'hydrometeor_optics', iverbose=driver_config%iverbose)
  end if

  ! --------------------------------------------------------
  ! Section 3: Read input data file
  ! --------------------------------------------------------

  ! Get NetCDF input file name
  call get_command_argument(2, file_name, status=istatus)
  if (istatus /= 0) then
    error stop 'Failed to read name of input NetCDF file as string of length < 512'
  end if

  ! Open the file and configure the way it is read
  call file%open(trim(file_name), iverbose=driver_config%iverbose)

  ! Get NetCDF output file name
  call get_command_argument(3, file_name, status=istatus)
  if (istatus /= 0) then
    error stop 'Failed to read name of output NetCDF file as string of length < 512'
  end if

  ! 2D arrays are assumed to be stored in the file with height varying
  ! more rapidly than column index. Specifying "true" here transposes
  ! all 2D arrays so that the column index varies fastest within the
  ! program.
  call file%transpose_matrices(.true.)

  ! Read input variables from NetCDF file
  call read_input(file, config, driver_config, ncol, nlev, &
       &          geometry, single_level, thermodynamics, &
       &          gas, cloud, aerosol)

  ! Close input file
  call file%close()

  ! Set first and last columns to process
  if (driver_config%iendcol < 1 .or. driver_config%iendcol > ncol) then
    driver_config%iendcol = ncol
  end if

  if (driver_config%istartcol > driver_config%iendcol) then
    write(nulout,'(a,i0,a,i0,a,i0,a)') '*** Error: requested column range (', &
         &  driver_config%istartcol, &
         &  ' to ', driver_config%iendcol, ') is out of the range in the data (1 to ', &
         &  ncol, ')'
    stop 1
  end if
  
  ! Store inputs
  if (driver_config%do_save_inputs) then
    call save_inputs('inputs.nc', config, single_level, thermodynamics, &
         &                gas, cloud, aerosol, &
         &                lat=spread(0.0_jprb,1,ncol), &
         &                lon=spread(0.0_jprb,1,ncol), &
         &                iverbose=driver_config%iverbose)
  end if

  ! --------------------------------------------------------
  ! Section 4: Call radiation scheme
  ! --------------------------------------------------------

  ! Ensure the units of the gas mixing ratios are what is required
  ! by the gas absorption model
  call set_gas_units(config, gas)

  ! Compute saturation with respect to liquid (needed for aerosol
  ! hydration) call
  call thermodynamics%calc_saturation_wrt_liquid(driver_config%istartcol,driver_config%iendcol)

  ! Check inputs are within physical bounds, printing message if not
  is_out_of_bounds =     gas%out_of_physical_bounds(driver_config%istartcol, driver_config%iendcol, &
       &                                            driver_config%do_correct_unphysical_inputs) &
       & .or.   single_level%out_of_physical_bounds(driver_config%istartcol, driver_config%iendcol, &
       &                                            driver_config%do_correct_unphysical_inputs) &
       & .or. thermodynamics%out_of_physical_bounds(driver_config%istartcol, driver_config%iendcol, &
       &                                            driver_config%do_correct_unphysical_inputs) &
       & .or.          cloud%out_of_physical_bounds(driver_config%istartcol, driver_config%iendcol, &
       &                                            driver_config%do_correct_unphysical_inputs) &
       & .or.        aerosol%out_of_physical_bounds(driver_config%istartcol, driver_config%iendcol, &
       &                                            driver_config%do_correct_unphysical_inputs) 
  
  ! Allocate memory for the flux profiles, which may include arrays
  ! of dimension n_bands_sw/n_bands_lw, so must be called after
  ! setup_radiation
  if (config%do_radiances) then
    call flux%allocate_radiances_only(config, 1, ncol)
  else
    call flux%allocate(config, 1, ncol, nlev, do_divergence=use_ecrad3d_solver)
  end if
  
  if (driver_config%iverbose >= 2) then
    write(nulout,'(a)')  'Performing radiative transfer calculations'
  end if
  
  ! Option of repeating calculation multiple time for more accurate
  ! profiling
#ifndef NO_OPENMP
  tstart = omp_get_wtime() 
#endif
  do jrepeat = 1,driver_config%nrepeat

    if (use_ecrad3d_solver) then
      ! Call the ECRAD3D radiation scheme
      call ecrad3d(config, geometry, single_level, thermodynamics, gas, cloud, aerosol, &
           &       flux)
    else
      
      if (driver_config%do_parallel) then
        ! Run radiation scheme over blocks of columns in parallel
      
        ! Compute number of blocks to process
        nblock = (driver_config%iendcol - driver_config%istartcol &
             &  + driver_config%nblocksize) / driver_config%nblocksize
        
        !$OMP PARALLEL DO PRIVATE(istartcol, iendcol) SCHEDULE(RUNTIME)
        do jblock = 1, nblock
          ! Specify the range of columns to process.
          istartcol = (jblock-1) * driver_config%nblocksize &
               &    + driver_config%istartcol
          iendcol = min(istartcol + driver_config%nblocksize - 1, &
               &        driver_config%iendcol)
          
          if (driver_config%iverbose >= 3) then
#ifndef NO_OPENMP
            write(nulout,'(a,i0,a,i0,a,i0)')  'Thread ', omp_get_thread_num(), &
                 &  ' processing columns ', istartcol, '-', iendcol
#else
            write(nulout,'(a,i0,a,i0)')  'Processing columns ', istartcol, '-', iendcol
#endif
          end if
        
          ! Call the ECRAD radiation scheme
          call radiation(ncol, nlev, istartcol, iendcol, config, &
               &  single_level, thermodynamics, gas, cloud, aerosol, flux)
          
        end do
        !$OMP END PARALLEL DO
        
      end if
    end if
  end do

#ifndef NO_OPENMP
  tstop = omp_get_wtime()
  write(nulout, '(a,g12.5,a)') 'Time elapsed in radiative transfer: ', tstop-tstart, ' seconds'
#endif

  ! --------------------------------------------------------
  ! Section 5: Check and save output
  ! --------------------------------------------------------

  is_out_of_bounds = flux%out_of_physical_bounds(driver_config%istartcol, driver_config%iendcol)

  ! Save the fluxes or radiances in the output file
  if (config%do_radiances) then
    ! Store radiances (or brightness temperatures) in the output file
    if (driver_config%do_save_brightness_temperature) then
      call flux%convert_radiance_to_brightness_temperature(ncol, config%gas_optics_lw)
    end if

    if (.not. geometry%is_3d) then
      call save_radiances(file_name, config, single_level, flux, &
           &   iverbose=driver_config%iverbose, is_hdf5_file=driver_config%do_write_hdf5, &
           &   experiment_name=driver_config%experiment_name)

    else
      call save_radiances_3d(file_name, config, geometry, single_level, flux, &
           &   iverbose=driver_config%iverbose, is_hdf5_file=driver_config%do_write_hdf5, &
           &   experiment_name=driver_config%experiment_name)
    end if
  else
    ! Store fluxes
    if (.not. geometry%is_3d) then
      call save_fluxes(file_name, config, thermodynamics, flux, &
           &   iverbose=driver_config%iverbose, is_hdf5_file=driver_config%do_write_hdf5, &
           &   experiment_name=driver_config%experiment_name, &
           &   is_double_precision=driver_config%do_write_double_precision)
    else
      call save_fluxes_3d(file_name, config, geometry, thermodynamics, flux, &
           &   iverbose=driver_config%iverbose, is_hdf5_file=driver_config%do_write_hdf5, &
           &   experiment_name=driver_config%experiment_name, &
           &   is_double_precision=driver_config%do_write_double_precision)
    end if
  end if
  
  if (driver_config%iverbose >= 2) then
    write(nulout,'(a)') '------------------------------------------------------------------------------------'
  end if

end program ecrad3d_driver
