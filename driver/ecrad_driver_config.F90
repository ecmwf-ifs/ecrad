! ecrad_driver_config.F90 - Configure driver for offline ecRad radiation scheme
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

module ecrad_driver_config

  use parkind1,                      only : jprb

  implicit none

  public

  ! Max length of "experiment" global attribute
  integer, parameter :: NMaxStringLength = 2000

  ! Maximum number of spectral diagnostics
  integer, parameter :: NMaxSpectralDiag = 256
  
  type driver_config_type

     ! Parallel settings
     logical :: do_parallel
     integer :: nblocksize ! Number of columns processed at once

     ! Override values from the radiation_override namelist (mostly
     ! related to clouds): these will override any values in the
     ! NetCDF data file (or scale them)
     real(jprb) :: fractional_std_override
     real(jprb) :: overlap_decorr_length_override
     real(jprb) :: high_inv_effective_size_override   = -1.0_jprb ! m-1
     real(jprb) :: middle_inv_effective_size_override = -1.0_jprb ! m-1
     real(jprb) :: low_inv_effective_size_override    = -1.0_jprb ! m-1
     real(jprb) :: effective_size_scaling
     real(jprb) :: sw_albedo_override
     real(jprb) :: lw_emissivity_override
     real(jprb) :: q_liq_scaling, q_ice_scaling
     real(jprb) :: cloud_fraction_scaling
     real(jprb) :: overlap_decorr_length_scaling
     real(jprb) :: skin_temperature_override ! K
     real(jprb) :: solar_irradiance_override ! W m-2
     real(jprb) :: solar_cycle_multiplier_override
     real(jprb) :: cos_sza_override
     real(jprb) :: cloud_inhom_separation_factor  = 1.0_jprb
     real(jprb) :: cloud_separation_scale_surface = -1.0_jprb
     real(jprb) :: cloud_separation_scale_toa     = -1.0_jprb
     real(jprb) :: cloud_separation_scale_power   = 1.0_jprb
     real(jprb) :: h2o_scaling    = 1.0_jprb
     real(jprb) :: co2_scaling    = 1.0_jprb
     real(jprb) :: o3_scaling     = 1.0_jprb
     real(jprb) :: co_scaling     = 1.0_jprb
     real(jprb) :: ch4_scaling    = 1.0_jprb
     real(jprb) :: n2o_scaling    = 1.0_jprb
     real(jprb) :: o2_scaling     = 1.0_jprb
     real(jprb) :: cfc11_scaling  = 1.0_jprb
     real(jprb) :: cfc12_scaling  = 1.0_jprb
     real(jprb) :: hcfc22_scaling = 1.0_jprb
     real(jprb) :: ccl4_scaling   = 1.0_jprb
     real(jprb) :: no2_scaling    = 1.0_jprb

     ! Optional monotonically increasing wavelength bounds (m) for
     ! shortwave spectral flux diagnostics, to be written to
     ! sw_diagnostic_file_name
     real(jprb) :: sw_diag_wavelength_bound(NMaxSpectralDiag+1) = -1.0_jprb

     ! Name of file to write shortwave spectral diagnostics to, but
     ! only if sw_diag_wavelength_bound is populated via the namelist
     character(len=NMaxStringLength) :: sw_diag_file_name = 'sw_diagnostics.nc'

     ! Number of shortwave diagnostics, worked out from first
     ! unassigned value in sw_diag_wavelength_bound after reading
     ! namelist
     integer :: n_sw_diag
     
     ! Volume mixing ratios (m3 m-3) in model layers (or equivalently
     ! mole fractions (mol mol-1)) are typically stored in the input
     ! file with a name like co2_vmr, but the suffix can be overridden
     ! by the user
     character(len=32) :: vmr_suffix_str = '_vmr'

     ! Process a limited number of columns (iendcol=0 indicates to
     ! process from istartcol up to the end)
     integer :: istartcol, iendcol

     ! Save inputs in "inputs.nc"
     logical :: do_save_inputs
     
     ! Save aerosol optical properties to "aerosol_optics.nc"
     logical :: do_save_aerosol_optics

     ! Save aerosol optical properties to "hydrometeor_optics*.nc"
     logical :: do_save_cloud_optics

     ! Save only net and surface/TOA fluxes, rather than up and down
     logical :: do_save_net_fluxes
     
     ! Do we ignore the inv_inhom_effective_size variable and instead
     ! assume the scale of cloud inhomogeneities is the same as the
     ! scale of the clouds themselves?
     logical :: do_ignore_inhom_effective_size = .false.

     ! Number of repeats (for benchmarking)
     integer :: nrepeat

     ! Do we correct unphysical inputs (e.g. negative gas concentrations)?
     logical :: do_correct_unphysical_inputs = .false.

     ! Do we write NetCDF4/HDF5 file format, needed for very large
     ! files?
     logical :: do_write_hdf5 = .false.

     ! Do we write fluxes in double precision?
     logical :: do_write_double_precision = .false.

     ! Name of the experiment, to save in output file
     character(len=NMaxStringLength) :: experiment_name = ''

     ! Control verbosity in driver routine: 0=none (no output to
     ! standard output; write to standard error only if an error
     ! occurs), 1=warning, 2=info, 3=progress, 4=detailed, 5=debug
     integer :: iverbose

   contains
     procedure :: read => read_config_from_namelist

  end type driver_config_type

contains

  !---------------------------------------------------------------------
  ! This subroutine reads configuration data from a namelist file, and
  ! anything that is not in the namelists will be set to default
  ! values. If optional output argument "is_success" is present, then on
  ! error (e.g. missing file) it will be set to .false.; if this
  ! argument is missing then on error the program will be aborted.
  subroutine read_config_from_namelist(this, file_name, is_success)

    use radiation_io, only : nulerr, radiation_abort

    class(driver_config_type), intent(inout) :: this
    character(*), intent(in)          :: file_name
    logical, intent(out), optional    :: is_success

    integer :: iosopen ! Status after calling open

    ! Override and scaling values
    real(jprb) :: fractional_std
    real(jprb) :: overlap_decorr_length
    real(jprb) :: inv_effective_size
    real(jprb) :: high_inv_effective_size
    real(jprb) :: middle_inv_effective_size
    real(jprb) :: low_inv_effective_size
    real(jprb) :: effective_size_scaling
    real(jprb) :: sw_albedo
    real(jprb) :: lw_emissivity
    real(jprb) :: q_liquid_scaling, q_ice_scaling
    real(jprb) :: cloud_fraction_scaling
    real(jprb) :: overlap_decorr_length_scaling
    real(jprb) :: skin_temperature
    real(jprb) :: cos_solar_zenith_angle
    real(jprb) :: solar_irradiance_override
    real(jprb) :: solar_cycle_multiplier_override
    real(jprb) :: cloud_inhom_separation_factor
    real(jprb) :: cloud_separation_scale_surface
    real(jprb) :: cloud_separation_scale_toa
    real(jprb) :: cloud_separation_scale_power
    real(jprb) :: h2o_scaling   
    real(jprb) :: co2_scaling   
    real(jprb) :: o3_scaling    
    real(jprb) :: co_scaling    
    real(jprb) :: ch4_scaling   
    real(jprb) :: n2o_scaling
    real(jprb) :: o2_scaling    
    real(jprb) :: cfc11_scaling 
    real(jprb) :: cfc12_scaling 
    real(jprb) :: hcfc22_scaling
    real(jprb) :: ccl4_scaling  
    real(jprb) :: no2_scaling   
    real(jprb) :: sw_diag_wavelength_bound(NMaxSpectralDiag+1)
    character(len=NMaxStringLength) :: sw_diag_file_name
    character(len=32) :: vmr_suffix_str
    character(len=NMaxStringLength) :: experiment_name

    ! Parallel settings
    logical :: do_parallel
    integer :: nblocksize

    logical :: do_save_inputs, do_save_aerosol_optics, do_save_net_fluxes, &
         &  do_save_cloud_optics, do_ignore_inhom_effective_size, &
         &  do_correct_unphysical_inputs, do_write_hdf5, &
         &  do_write_double_precision
    integer :: nrepeat

    ! Process a limited number of columns (iendcol=0 indicates to
    ! process from istartcol up to the end)
    integer :: istartcol, iendcol

    ! Verbosity
    integer :: iverbose

    ! Are we going to override the effective size?
    logical :: do_override_eff_size

    ! Loop index
    integer :: jdiag
    
    namelist /radiation_driver/ fractional_std, &
         &  overlap_decorr_length, inv_effective_size, sw_albedo, &
         &  high_inv_effective_size, middle_inv_effective_size, &
         &  low_inv_effective_size, cloud_inhom_separation_factor, &
         &  effective_size_scaling, cos_solar_zenith_angle, &
         &  lw_emissivity, q_liquid_scaling, q_ice_scaling, &
         &  istartcol, iendcol, solar_irradiance_override, &
         &  solar_cycle_multiplier_override, &
         &  cloud_fraction_scaling, overlap_decorr_length_scaling, &
         &  skin_temperature, do_parallel, nblocksize, iverbose, &
         &  nrepeat, do_save_inputs, do_ignore_inhom_effective_size, &
         &  do_save_aerosol_optics, do_save_net_fluxes, do_save_cloud_optics, &
         &  cloud_separation_scale_toa, cloud_separation_scale_surface, &
         &  cloud_separation_scale_power, do_correct_unphysical_inputs, &
         &  do_write_hdf5, h2o_scaling, co2_scaling, o3_scaling, co_scaling, &
         &  ch4_scaling, o2_scaling, cfc11_scaling, cfc12_scaling, &
         &  hcfc22_scaling, no2_scaling, n2o_scaling, ccl4_scaling, &
         &  vmr_suffix_str, experiment_name, do_write_double_precision, &
         &  sw_diag_wavelength_bound, sw_diag_file_name

    ! Default values
    do_parallel = .true.
    do_save_inputs = .false.
    do_save_aerosol_optics = .false.
    do_save_cloud_optics = .false.
    do_save_net_fluxes = .false.
    do_ignore_inhom_effective_size = .false.
    nblocksize = 8

    ! Negative values indicate no override will take place
    fractional_std = -1.0_jprb
    overlap_decorr_length = -1.0_jprb
    inv_effective_size = -1.0_jprb
    high_inv_effective_size = -1.0_jprb
    middle_inv_effective_size = -1.0_jprb
    low_inv_effective_size = -1.0_jprb
    effective_size_scaling = -1.0_jprb
    sw_albedo = -1.0_jprb
    lw_emissivity = -1.0_jprb
    q_liquid_scaling = -1.0_jprb
    q_ice_scaling = -1.0_jprb
    cloud_fraction_scaling = -1.0_jprb
    overlap_decorr_length_scaling = -1.0_jprb
    skin_temperature = -1.0_jprb
    cos_solar_zenith_angle = -1.0_jprb
    solar_irradiance_override = -1.0_jprb
    solar_cycle_multiplier_override = -2.0e6_jprb
    cloud_inhom_separation_factor = 1.0_jprb
    cloud_separation_scale_toa = -1.0_jprb
    cloud_separation_scale_surface = -1.0_jprb
    cloud_separation_scale_power = 1.0_jprb
    h2o_scaling    = 1.0_jprb
    co2_scaling    = 1.0_jprb
    o3_scaling     = 1.0_jprb
    co_scaling     = 1.0_jprb
    ch4_scaling    = 1.0_jprb
    n2o_scaling    = 1.0_jprb
    o2_scaling     = 1.0_jprb
    cfc11_scaling  = 1.0_jprb
    cfc12_scaling  = 1.0_jprb
    hcfc22_scaling = 1.0_jprb
    ccl4_scaling   = 1.0_jprb
    no2_scaling    = 1.0_jprb
    vmr_suffix_str = '_vmr';
    iverbose = 2 ! Default verbosity is "warning"
    istartcol = 0
    iendcol = 0
    nrepeat = 1
    do_correct_unphysical_inputs = .false.
    do_write_hdf5 = .false.
    do_write_double_precision = .false.
    experiment_name = ''
    sw_diag_wavelength_bound = this%sw_diag_wavelength_bound
    sw_diag_file_name = this%sw_diag_file_name
    
    ! Open the namelist file and read the radiation_driver namelist
    open(unit=10, iostat=iosopen, file=trim(file_name))
    if (iosopen /= 0) then
      ! An error occurred
      if (present(is_success)) then
        is_success = .false.
        ! We now continue the subroutine so that the default values
        ! are placed in the config structure
      else
        write(nulerr,'(a,a,a)') '*** Error: namelist file "', &
             &                trim(file_name), '" not found'
        call radiation_abort('Driver configuration error')
      end if
    else
      ! Read the radiation_driver namelist, noting that it is not an
      ! error if this namelist is not present, provided all the required
      ! variables are present in the NetCDF data file instead
      read(unit=10, nml=radiation_driver)
      close(unit=10)
    end if

    ! Copy namelist data into configuration object
    this%do_parallel = do_parallel
    this%do_save_inputs = do_save_inputs
    this%do_save_aerosol_optics = do_save_aerosol_optics
    this%do_save_cloud_optics = do_save_cloud_optics
    this%do_save_net_fluxes = do_save_net_fluxes
    this%do_ignore_inhom_effective_size = do_ignore_inhom_effective_size
    this%nblocksize = nblocksize
    this%iverbose = iverbose
    this%nrepeat = nrepeat
    if (istartcol < 1) then
      this%istartcol = 1
    else
      this%istartcol = istartcol
    end if
    if (iendcol < 1) then
      this%iendcol = 0
    else
      this%iendcol = iendcol
    end if

    ! Set override values
    this%fractional_std_override = fractional_std
    this%overlap_decorr_length_override = overlap_decorr_length

    do_override_eff_size = .false.
    if (inv_effective_size >= 0.0_jprb) then
      this%high_inv_effective_size_override = inv_effective_size
      this%middle_inv_effective_size_override = inv_effective_size
      this%low_inv_effective_size_override = inv_effective_size
    end if
    if (high_inv_effective_size >= 0.0_jprb) then
      this%high_inv_effective_size_override = high_inv_effective_size
      do_override_eff_size = .true.
    end if
    if (middle_inv_effective_size >= 0.0_jprb) then
      this%middle_inv_effective_size_override = middle_inv_effective_size
      do_override_eff_size = .true.
    end if
    if (low_inv_effective_size >= 0.0_jprb) then
      this%low_inv_effective_size_override = low_inv_effective_size
      do_override_eff_size = .true.
    end if

    if (do_override_eff_size &
         &  .and. (this%high_inv_effective_size_override < 0.0_jprb &
              .or. this%middle_inv_effective_size_override < 0.0_jprb &
              .or. this%low_inv_effective_size_override < 0.0_jprb)) then
      write(nulerr,'(a)') '*** Error: inverse effective cloud size not specified for high, middle and low clouds"'
      call radiation_abort('Driver configuration error')
    end if

    this%effective_size_scaling = effective_size_scaling
    this%sw_albedo_override = sw_albedo
    this%lw_emissivity_override = lw_emissivity
    this%q_liq_scaling = q_liquid_scaling
    this%q_ice_scaling = q_ice_scaling
    this%cloud_fraction_scaling = cloud_fraction_scaling
    this%overlap_decorr_length_scaling = overlap_decorr_length_scaling
    this%skin_temperature_override = skin_temperature
    this%cos_sza_override = cos_solar_zenith_angle
    this%solar_irradiance_override = solar_irradiance_override
    this%solar_cycle_multiplier_override = solar_cycle_multiplier_override
    this%cloud_inhom_separation_factor = cloud_inhom_separation_factor
    this%cloud_separation_scale_toa = cloud_separation_scale_toa
    this%cloud_separation_scale_surface = cloud_separation_scale_surface
    this%cloud_separation_scale_power = cloud_separation_scale_power
    this%do_correct_unphysical_inputs = do_correct_unphysical_inputs
    this%do_write_hdf5  = do_write_hdf5
    this%do_write_double_precision = do_write_double_precision
    this%h2o_scaling    = h2o_scaling
    this%co2_scaling    = co2_scaling
    this%o3_scaling     = o3_scaling
    this%co_scaling     = co_scaling
    this%ch4_scaling    = ch4_scaling
    this%n2o_scaling    = n2o_scaling
    this%o2_scaling     = o2_scaling
    this%cfc11_scaling  = cfc11_scaling
    this%cfc12_scaling  = cfc12_scaling
    this%hcfc22_scaling = hcfc22_scaling
    this%ccl4_scaling   = ccl4_scaling
    this%no2_scaling    = no2_scaling
    this%vmr_suffix_str = trim(vmr_suffix_str)
    this%experiment_name= trim(experiment_name)
    
    this%sw_diag_file_name = trim(sw_diag_file_name)
    this%sw_diag_wavelength_bound = sw_diag_wavelength_bound
    ! Work out number of shortwave diagnostics from first negative
    ! wavelength bound, noting that the number of diagnostics is one
    ! fewer than the number of valid bounds
    do jdiag = 0,NMaxSpectralDiag
      if (this%sw_diag_wavelength_bound(jdiag+1) < 0.0_jprb) then
        this%n_sw_diag = max(0,jdiag-1)
        exit
      end if
    end do
    
  end subroutine read_config_from_namelist

end module ecrad_driver_config
