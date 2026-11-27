! ecrad4py.F90 - C-interoperable interface to ecRad for the pyecrad Python package
!
! (C) Copyright 2025- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
!
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!
! Author:  Sébastien Riette
! Email:   sebastien.riette@meteo.fr
!
! This module provides procedures with C binding that are compiled
! into a shared library (libecrad4py) and called from Python by the
! pyecrad package (see pyecrad/ecrad4py.py):
!
!   ecrad4py_setup  - read the "radiation" namelist and set up the
!                     scheme (load gas, cloud and aerosol optics);
!                     returns an opaque handle to the configuration
!   ecrad4py_free   - release a configuration created by ecrad4py_setup
!   ecrad4py_run    - run the radiation scheme on a set of columns
!   ecrad4py_get_*  - return the values of some ecRad constants
!
! The sequence of operations in ecrad4py_run follows the offline
! driver (driver/ecrad_driver.F90 and ecrad_driver_read_input.F90),
! so that the same inputs give the same fluxes. Floating-point
! arguments are always passed in double precision; they are converted
! to the working precision (jprb) when copied into the ecRad derived
! types, so the library may be built in either precision.
!
! Arrays are dimensioned as in the ecRad derived types, i.e. with the
! column index varying fastest: (ncol,nlev) for layer quantities,
! (ncol,nlev+1) for half-level quantities and (ncol,nlev-1) for the
! cloud overlap parameter. All quantities are in SI units.
!
! Each configuration is allocated on the heap and its address is
! returned to the caller, so there is no mutable module data. A given
! configuration may be used by several concurrent calls to
! ecrad4py_run, but ecrad4py_setup and ecrad4py_free must not be
! called concurrently with ecrad4py_run on the same handle.
!
! Errors are handled by radiation_abort, which terminates the calling
! process; the Python layer therefore checks its inputs beforehand.

module ecrad4py

  use, intrinsic :: iso_c_binding, only : c_ptr, c_null_ptr, c_loc, c_f_pointer, &
       &                                  c_associated, c_char, c_int64_t, c_double
  use radiation_config,            only : config_type

  implicit none

  public

  ! Length of the character arrays used to pass file and directory
  ! names
  integer, parameter :: NStringLength = 512

contains

  !---------------------------------------------------------------------
  ! Convert a blank-padded C character array to a Fortran string
  function c_chars_to_string(c_chars) result(str)

    character(kind=c_char), intent(in) :: c_chars(:)
    character(len=size(c_chars))       :: str

    integer :: jc

    do jc = 1,size(c_chars)
      str(jc:jc) = c_chars(jc)
    end do

  end function c_chars_to_string


  !---------------------------------------------------------------------
  ! Read the "radiation" namelist from namelist_file_name and set up
  ! the radiation scheme. If directory_name is not blank, it overrides
  ! the directory_name namelist entry, i.e. the location of the data
  ! files. A handle to the configuration is returned in config_handle.
  subroutine ecrad4py_setup(namelist_file_name, directory_name, config_handle) &
       &  bind(C, name='ecrad4py_setup')

    use radiation_interface, only : setup_radiation
    use radiation_io,        only : nulerr, radiation_abort
    use yomhook,             only : dr_hook_init
#ifdef HAVE_FIAT
    use mpl_module,          only : mpl_init
#endif

    character(kind=c_char), intent(in)  :: namelist_file_name(NStringLength)
    character(kind=c_char), intent(in)  :: directory_name(NStringLength)
    type(c_ptr),            intent(out) :: config_handle

    type(config_type), pointer          :: config
    character(len=NStringLength)        :: file_name_str, directory_name_str
    logical                             :: is_success

    ! Initialise MPL and DR_HOOK as the offline driver does (the RRTMG
    ! setup queries MPL, which aborts if MPL is not initialised). Both
    ! calls do nothing if the initialisation has already been done, for
    ! instance by a previous configuration.
#ifdef HAVE_FIAT
    call mpl_init(ldinfo=.false.)
#endif
    call dr_hook_init()

    file_name_str      = c_chars_to_string(namelist_file_name)
    directory_name_str = c_chars_to_string(directory_name)

    allocate(config)

    ! Read "radiation" namelist into radiation configuration type
    call config%read(file_name=trim(file_name_str), is_success=is_success)
    if (.not. is_success) then
      write(nulerr,'(a,a)') '*** Error reading the "radiation" namelist from ', &
           &  trim(file_name_str)
      call radiation_abort()
    end if

    if (len_trim(directory_name_str) > 0) then
      if (len_trim(directory_name_str) > len(config%directory_name)) then
        call radiation_abort('ecrad4py_setup: directory name too long')
      end if
      config%directory_name = trim(directory_name_str)
    end if

    ! Setup the radiation scheme: load the coefficients for gas,
    ! cloud and aerosol optics
    call setup_radiation(config)

    config_handle = c_loc(config)

  end subroutine ecrad4py_setup


  !---------------------------------------------------------------------
  ! Release a configuration created by ecrad4py_setup
  subroutine ecrad4py_free(config_handle) bind(C, name='ecrad4py_free')

    type(c_ptr), intent(inout) :: config_handle

    type(config_type), pointer :: config

    if (c_associated(config_handle)) then
      call c_f_pointer(config_handle, config)
      deallocate(config)
      config_handle = c_null_ptr
    end if

  end subroutine ecrad4py_free


  !---------------------------------------------------------------------
  ! Value used by ecRad to denote gas concentrations expressed as
  ! volume mixing ratio
  function ecrad4py_get_ivolumemixingratio() result(iunits) &
       &  bind(C, name='ecrad4py_get_ivolumemixingratio')

    use radiation_gas, only : IVolumeMixingRatio

    integer(c_int64_t) :: iunits

    iunits = int(IVolumeMixingRatio, c_int64_t)

  end function ecrad4py_get_ivolumemixingratio


  !---------------------------------------------------------------------
  ! Value used by ecRad to denote gas concentrations expressed as mass
  ! mixing ratio
  function ecrad4py_get_imassmixingratio() result(iunits) &
       &  bind(C, name='ecrad4py_get_imassmixingratio')

    use radiation_gas, only : IMassMixingRatio

    integer(c_int64_t) :: iunits

    iunits = int(IMassMixingRatio, c_int64_t)

  end function ecrad4py_get_imassmixingratio


  !---------------------------------------------------------------------
  ! Return 1 if the library was compiled in single precision (jprb
  ! is 4 bytes), 0 otherwise
  function ecrad4py_is_single_precision() result(is_single) &
       &  bind(C, name='ecrad4py_is_single_precision')

    use parkind1, only : jprb, jprd

    integer(c_int64_t) :: is_single

    if (kind(1.0_jprb) == kind(1.0_jprd)) then
      is_single = 0
    else
      is_single = 1
    end if

  end function ecrad4py_is_single_precision


  !---------------------------------------------------------------------
  ! Return some configuration settings, as 0 (false) or 1 (true) for
  ! logicals, so that the caller can check its inputs before calling
  ! ecrad4py_run
  subroutine ecrad4py_get_config_info(config_handle, ido_lw, ido_sw, ido_clouds, &
       &  iuse_aerosols, n_aerosol_types, iuse_spartacus) &
       &  bind(C, name='ecrad4py_get_config_info')

    use radiation_config, only : ISolverSPARTACUS
    use radiation_io,     only : radiation_abort

    type(c_ptr),        intent(in)  :: config_handle
    integer(c_int64_t), intent(out) :: ido_lw, ido_sw, ido_clouds, iuse_aerosols
    integer(c_int64_t), intent(out) :: n_aerosol_types, iuse_spartacus

    type(config_type), pointer :: config

    if (.not. c_associated(config_handle)) then
      call radiation_abort('ecrad4py_get_config_info: configuration has not been set up')
    end if
    call c_f_pointer(config_handle, config)

    ido_lw          = merge(1_c_int64_t, 0_c_int64_t, config%do_lw)
    ido_sw          = merge(1_c_int64_t, 0_c_int64_t, config%do_sw)
    ido_clouds      = merge(1_c_int64_t, 0_c_int64_t, config%do_clouds)
    iuse_aerosols   = merge(1_c_int64_t, 0_c_int64_t, config%use_aerosols)
    n_aerosol_types = int(config%n_aerosol_types, c_int64_t)
    iuse_spartacus  = merge(1_c_int64_t, 0_c_int64_t, &
         &  config%i_solver_sw == ISolverSPARTACUS .or. config%i_solver_lw == ISolverSPARTACUS)

  end subroutine ecrad4py_get_config_info


  !---------------------------------------------------------------------
  ! Run the radiation scheme on ncol columns of nlev layers using the
  ! configuration config_handle. Columns are processed in blocks of
  ! nblocksize columns, in parallel if OpenMP is enabled; if
  ! nblocksize <= 0, all columns are processed in a single call.
  !
  ! Optional arguments may be absent (null pointers from C):
  !  - sw_albedo_direct: if absent, sw_albedo is used for both direct
  !    and diffuse radiation
  !  - iseed: if absent, iseed(jcol)=jcol as in the offline driver
  !  - overlap_param: if absent, computed from the overlap
  !    decorrelation length overlap_decorr_length (m)
  !  - fractional_std: if absent, set to zero
  !  - inv_cloud_effective_size: required by the SPARTACUS solver
  !  - inv_inhom_effective_size: optional for the SPARTACUS solver
  !  - aerosol_mmr: required if config%use_aerosols
  !  - gases: each gas is considered only if both its units and its
  !    mixing ratio are present
  !
  ! Fluxes that are not computed with this configuration (e.g. clear
  ! sky fluxes if config%do_clear is false) are set to NaN.
  subroutine ecrad4py_run(config_handle, ncol, nlev, nblocksize, &
       &  pressure_hl, temperature_hl, solar_irradiance, &
       &  spectral_solar_cycle_multiplier, cos_solar_zenith_angle, &
       &  cloud_fraction, q_liquid, re_liquid, q_ice, re_ice, &
       &  skin_temperature, nalbedobands, sw_albedo, nemissivitybands, lw_emissivity, &
       &  overlap_decorr_length, sw_albedo_direct, iseed, overlap_param, fractional_std, &
       &  inv_cloud_effective_size, inv_inhom_effective_size, &
       &  naerosoltypes, aerosol_mmr, &
       &  h2o_units, h2o, co2_units, co2, o3_units, o3, n2o_units, n2o, &
       &  co_units, co, ch4_units, ch4, o2_units, o2, cfc11_units, cfc11, &
       &  cfc12_units, cfc12, hcfc22_units, hcfc22, ccl4_units, ccl4, no2_units, no2, &
       &  lw_up, lw_dn, lw_up_clear, lw_dn_clear, cloud_cover_lw, &
       &  sw_up, sw_dn, sw_up_clear, sw_dn_clear, &
       &  sw_dn_direct, sw_dn_direct_clear, cloud_cover_sw) &
       &  bind(C, name='ecrad4py_run')

    use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan
    use parkind1,                 only : jprb
    use radiation_interface,      only : radiation, set_gas_units
    use radiation_single_level,   only : single_level_type
    use radiation_thermodynamics, only : thermodynamics_type
    use radiation_cloud,          only : cloud_type
    use radiation_aerosol,        only : aerosol_type
    use radiation_flux,           only : flux_type
    use radiation_config,         only : ISolverSPARTACUS
    use radiation_gas,            only : gas_type, IMassMixingRatio, &
         &  IH2O, ICO2, IO3, IN2O, ICO, ICH4, IO2, ICFC11, ICFC12, IHCFC22, ICCl4, INO2
    use radiation_io,             only : radiation_abort

    ! Configuration created by ecrad4py_setup
    type(c_ptr),        intent(in) :: config_handle
    ! Number of columns and layers, and number of columns per block
    integer(c_int64_t), intent(in) :: ncol, nlev, nblocksize
    ! Pressure (Pa) and temperature (K) at half levels
    real(c_double),     intent(in) :: pressure_hl(ncol,nlev+1)
    real(c_double),     intent(in) :: temperature_hl(ncol,nlev+1)
    ! Total solar irradiance at Earth (W m-2) and scaling of the
    ! spectral solar cycle (+1 solar maximum, -1 minimum, 0 mean)
    real(c_double),     intent(in) :: solar_irradiance
    real(c_double),     intent(in) :: spectral_solar_cycle_multiplier
    ! Cosine of the solar zenith angle
    real(c_double),     intent(in) :: cos_solar_zenith_angle(ncol)
    ! Cloud fraction (0-1), liquid and ice mass mixing ratios (kg/kg)
    ! and effective radii (m)
    real(c_double),     intent(in) :: cloud_fraction(ncol,nlev)
    real(c_double),     intent(in) :: q_liquid(ncol,nlev), re_liquid(ncol,nlev)
    real(c_double),     intent(in) :: q_ice(ncol,nlev), re_ice(ncol,nlev)
    ! Skin temperature (K)
    real(c_double),     intent(in) :: skin_temperature(ncol)
    ! Shortwave albedo in each albedo band
    integer(c_int64_t), intent(in) :: nalbedobands
    real(c_double),     intent(in) :: sw_albedo(ncol,nalbedobands)
    ! Longwave emissivity in each emissivity band
    integer(c_int64_t), intent(in) :: nemissivitybands
    real(c_double),     intent(in) :: lw_emissivity(ncol,nemissivitybands)
    ! Overlap decorrelation length (m), used if overlap_param is absent
    real(c_double),     intent(in) :: overlap_decorr_length
    ! Optional direct shortwave albedo
    real(c_double), optional, intent(in) :: sw_albedo_direct(ncol,nalbedobands)
    ! Seed for the random number generator of McICA
    integer(c_int64_t), optional, intent(in) :: iseed(ncol)
    ! Cloud overlap parameter between adjacent layers
    real(c_double), optional, intent(in) :: overlap_param(ncol,nlev-1)
    ! Fractional standard deviation of in-cloud water content
    real(c_double), optional, intent(in) :: fractional_std(ncol,nlev)
    ! Inverse of the cloud and inhomogeneity effective sizes (m-1),
    ! needed by SPARTACUS
    real(c_double), optional, intent(in) :: inv_cloud_effective_size(ncol,nlev)
    real(c_double), optional, intent(in) :: inv_inhom_effective_size(ncol,nlev)
    ! Aerosol mass mixing ratios (kg/kg)
    integer(c_int64_t), intent(in) :: naerosoltypes
    real(c_double), optional, intent(in) :: aerosol_mmr(ncol,nlev,naerosoltypes)
    ! Gas units (IMassMixingRatio or IVolumeMixingRatio) and mixing
    ! ratios
    integer(c_int64_t), optional, intent(in) :: h2o_units, co2_units, o3_units, &
         &  n2o_units, co_units, ch4_units, o2_units, cfc11_units, cfc12_units, &
         &  hcfc22_units, ccl4_units, no2_units
    real(c_double), optional, intent(in), dimension(ncol,nlev) :: h2o, co2, o3, &
         &  n2o, co, ch4, o2, cfc11, cfc12, hcfc22, ccl4, no2
    ! Output longwave and shortwave fluxes (W m-2) at half levels and
    ! total cloud cover diagnosed by each solver
    real(c_double), intent(out) :: lw_up(ncol,nlev+1), lw_dn(ncol,nlev+1)
    real(c_double), intent(out) :: lw_up_clear(ncol,nlev+1), lw_dn_clear(ncol,nlev+1)
    real(c_double), intent(out) :: cloud_cover_lw(ncol)
    real(c_double), intent(out) :: sw_up(ncol,nlev+1), sw_dn(ncol,nlev+1)
    real(c_double), intent(out) :: sw_up_clear(ncol,nlev+1), sw_dn_clear(ncol,nlev+1)
    real(c_double), intent(out) :: sw_dn_direct(ncol,nlev+1), sw_dn_direct_clear(ncol,nlev+1)
    real(c_double), intent(out) :: cloud_cover_sw(ncol)

    type(config_type), pointer :: config
    type(single_level_type)    :: single_level
    type(thermodynamics_type)  :: thermodynamics
    type(cloud_type)           :: cloud
    type(gas_type)             :: gas
    type(aerosol_type)         :: aerosol
    type(flux_type)            :: flux

    ! Number of columns, layers and columns per block as default integers
    integer :: ncolumn, nlayer, nblockcol
    integer :: nblock, jblock, istartcol, iendcol, jcol

    real(c_double) :: nan

    if (.not. c_associated(config_handle)) then
      call radiation_abort('ecrad4py_run: configuration has not been set up')
    end if
    call c_f_pointer(config_handle, config)

    ncolumn = int(ncol)
    nlayer  = int(nlev)

    ! --------------------------------------------------------
    ! Thermodynamics
    ! --------------------------------------------------------

    call thermodynamics%allocate(ncolumn, nlayer)
    thermodynamics%pressure_hl    = real(pressure_hl, jprb)
    thermodynamics%temperature_hl = real(temperature_hl, jprb)

    ! --------------------------------------------------------
    ! Sun and surface
    ! --------------------------------------------------------

    call single_level%allocate(ncolumn, int(nalbedobands), int(nemissivitybands), &
         &  use_sw_albedo_direct=present(sw_albedo_direct), is_simple_surface=.true.)
    single_level%solar_irradiance = real(solar_irradiance, jprb)
    single_level%spectral_solar_cycle_multiplier = real(spectral_solar_cycle_multiplier, jprb)
    single_level%cos_sza          = real(cos_solar_zenith_angle, jprb)
    single_level%skin_temperature = real(skin_temperature, jprb)
    single_level%sw_albedo        = real(sw_albedo, jprb)
    if (present(sw_albedo_direct)) then
      single_level%sw_albedo_direct = real(sw_albedo_direct, jprb)
    end if
    single_level%lw_emissivity    = real(lw_emissivity, jprb)

    ! Seed for the random number generator in McICA: as in the
    ! offline driver, iseed(jcol)=jcol unless provided
    if (present(iseed)) then
      single_level%iseed = int(iseed)
    else
      do jcol = 1,ncolumn
        single_level%iseed(jcol) = jcol
      end do
    end if

    ! --------------------------------------------------------
    ! Clouds
    ! --------------------------------------------------------

    if (config%do_clouds) then
      ! Two hydrometeor types: (1) liquid and (2) ice, accessible via
      ! cloud%q_liq, cloud%q_ice, cloud%re_liq and cloud%re_ice
      call cloud%allocate(ncolumn, nlayer, use_inhom_effective_size=.false.)
      cloud%fraction = real(cloud_fraction, jprb)
      cloud%q_liq    = real(q_liquid, jprb)
      cloud%re_liq   = real(re_liquid, jprb)
      cloud%q_ice    = real(q_ice, jprb)
      cloud%re_ice   = real(re_ice, jprb)

      if (present(overlap_param)) then
        cloud%overlap_param = real(overlap_param, jprb)
      else
        call cloud%set_overlap_param(thermodynamics, real(overlap_decorr_length, jprb))
      end if

      if (present(fractional_std)) then
        cloud%fractional_std = real(fractional_std, jprb)
      else
        cloud%fractional_std = 0.0_jprb
      end if

      ! Cloud properties needed by SPARTACUS
      if (config%i_solver_sw == ISolverSPARTACUS &
           &  .or. config%i_solver_lw == ISolverSPARTACUS) then
        if (.not. present(inv_cloud_effective_size)) then
          call radiation_abort('ecrad4py_run: inv_cloud_effective_size is required by SPARTACUS')
        end if
        allocate(cloud%inv_cloud_effective_size(ncolumn,nlayer))
        cloud%inv_cloud_effective_size = real(inv_cloud_effective_size, jprb)
        ! If absent, the effective size of in-cloud inhomogeneities is
        ! taken equal to the cloud effective size, as in the offline
        ! driver
        if (present(inv_inhom_effective_size)) then
          allocate(cloud%inv_inhom_effective_size(ncolumn,nlayer))
          cloud%inv_inhom_effective_size = real(inv_inhom_effective_size, jprb)
        end if
      end if
    end if

    ! --------------------------------------------------------
    ! Aerosols
    ! --------------------------------------------------------

    if (config%use_aerosols) then
      if (.not. present(aerosol_mmr)) then
        call radiation_abort('ecrad4py_run: aerosol_mmr is required when use_aerosols is true')
      end if
      call aerosol%allocate(ncolumn, 1, nlayer, int(naerosoltypes))
      aerosol%mixing_ratio = real(aerosol_mmr, jprb)
    end if

    ! --------------------------------------------------------
    ! Gases
    ! --------------------------------------------------------

    call gas%allocate(ncolumn, nlayer)
    call put_optional_gas(IH2O,    h2o_units,    h2o)
    call put_optional_gas(ICO2,    co2_units,    co2)
    call put_optional_gas(IO3,     o3_units,     o3)
    call put_optional_gas(IN2O,    n2o_units,    n2o)
    call put_optional_gas(ICO,     co_units,     co)
    call put_optional_gas(ICH4,    ch4_units,    ch4)
    call put_optional_gas(IO2,     o2_units,     o2)
    call put_optional_gas(ICFC11,  cfc11_units,  cfc11)
    call put_optional_gas(ICFC12,  cfc12_units,  cfc12)
    call put_optional_gas(IHCFC22, hcfc22_units, hcfc22)
    call put_optional_gas(ICCl4,   ccl4_units,   ccl4)
    call put_optional_gas(INO2,    no2_units,    no2)

    ! Ensure the units of the gas mixing ratios are what is required
    ! by the gas absorption model
#ifdef BITIDENTITY_TESTING
    ! Follow the same conversion path as the offline driver so that
    ! the results can be compared bit for bit
    call gas%set_units(IMassMixingRatio)
#endif
    call set_gas_units(config, gas)

    ! Compute saturation with respect to liquid (needed for aerosol
    ! hydration)
    call thermodynamics%calc_saturation_wrt_liquid(1, ncolumn)

    ! --------------------------------------------------------
    ! Call radiation scheme
    ! --------------------------------------------------------

    ! Allocate memory for the flux profiles, which must be done after
    ! setup_radiation
    call flux%allocate(config, 1, ncolumn, nlayer)

    if (nblocksize <= 0) then
      call radiation(ncolumn, nlayer, 1, ncolumn, config, &
           &  single_level, thermodynamics, gas, cloud, aerosol, flux)
    else
      nblockcol = int(nblocksize)
      nblock = (ncolumn - 1 + nblockcol) / nblockcol
      !$OMP PARALLEL DO PRIVATE(istartcol, iendcol) SCHEDULE(RUNTIME)
      do jblock = 1,nblock
        istartcol = (jblock-1) * nblockcol + 1
        iendcol   = min(istartcol + nblockcol - 1, ncolumn)
        call radiation(ncolumn, nlayer, istartcol, iendcol, config, &
             &  single_level, thermodynamics, gas, cloud, aerosol, flux)
      end do
      !$OMP END PARALLEL DO
    end if

    ! --------------------------------------------------------
    ! Copy outputs, using NaN for quantities not computed
    ! --------------------------------------------------------

    nan = ieee_value(nan, ieee_quiet_nan)

    call copy_flux(flux%lw_up,       lw_up)
    call copy_flux(flux%lw_dn,       lw_dn)
    call copy_flux(flux%lw_up_clear, lw_up_clear)
    call copy_flux(flux%lw_dn_clear, lw_dn_clear)
    call copy_flux(flux%sw_up,       sw_up)
    call copy_flux(flux%sw_dn,       sw_dn)
    call copy_flux(flux%sw_up_clear, sw_up_clear)
    call copy_flux(flux%sw_dn_clear, sw_dn_clear)
    call copy_flux(flux%sw_dn_direct,       sw_dn_direct)
    call copy_flux(flux%sw_dn_direct_clear, sw_dn_direct_clear)

    ! Cloud cover is only meaningful if clouds are represented, and
    ! is not written by the offline driver otherwise
    cloud_cover_lw = nan
    cloud_cover_sw = nan
    if (config%do_clouds) then
      if (config%do_lw .and. associated(flux%cloud_cover_lw)) then
        cloud_cover_lw = real(flux%cloud_cover_lw, c_double)
      end if
      if (config%do_sw .and. associated(flux%cloud_cover_sw)) then
        cloud_cover_sw = real(flux%cloud_cover_sw, c_double)
      end if
    end if

    call flux%deallocate()
    call gas%deallocate()
    if (config%use_aerosols) then
      call aerosol%deallocate()
    end if
    if (config%do_clouds) then
      call cloud%deallocate()
    end if
    call single_level%deallocate()
    call thermodynamics%deallocate()

  contains

    ! Store a gas if both its units and mixing ratio are present
    subroutine put_optional_gas(igas, iunits, mixing_ratio)

      integer,                      intent(in) :: igas
      integer(c_int64_t), optional, intent(in) :: iunits
      real(c_double),     optional, intent(in) :: mixing_ratio(:,:)

      if (present(iunits) .and. present(mixing_ratio)) then
        call gas%put(igas, int(iunits), mixing_ratio)
      end if

    end subroutine put_optional_gas

    ! Copy a flux profile if it has been computed, or set it to NaN
    subroutine copy_flux(flux_in, flux_out)

      real(jprb), pointer, intent(in)  :: flux_in(:,:)
      real(c_double),      intent(out) :: flux_out(:,:)

      if (associated(flux_in)) then
        flux_out = real(flux_in, c_double)
      else
        flux_out = nan
      end if

    end subroutine copy_flux

  end subroutine ecrad4py_run

end module ecrad4py
