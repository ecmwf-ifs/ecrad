! ecrad4py.F90 - Driver for python ECRAD radiation scheme
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
! Author:  Sébastien Riette
! Email:   sebastien.riette@meteo.fr
!
! This module contains two subroutine:
! 1) setup with namelist filereading (to be called once)
! 2) run to run the scheme on input profile
!
! The strucuture of the code is taken from the offline driver

module ecrad4py

  use radiation_config,         only : config_type

  implicit none

  type(config_type)         :: config

  contains

    subroutine setup(namelist_file_name, directory_name)
      use radiation_interface,      only : setup_radiation
      character(len=512), intent(in) :: namelist_file_name
      character(len=511), intent(in) :: directory_name

      ! Read "radiation" namelist into radiation configuration type
      call config%read(file_name=namelist_file_name)
      config%directory_name = directory_name

      ! Setup the radiation scheme: load the coefficients for gas and
      ! cloud optics, currently from RRTMG
      call setup_radiation(config)
    end subroutine setup

    subroutine run(ncol, nlev, pressure_hl, temperature_hl, solar_irradiance, &
                  &spectral_solar_cycle_multiplier, &
                  &cos_solar_zenith_angle, cloud_fraction, fractional_std, &
                  &q_liquid, re_liquid, q_ice, re_ice, iseed, overlap_param, &
                  &skin_temperature, nalbedobands, sw_albedo, sw_albedo_direct, &
                  &nemissivitygpoints, lw_emissivity, q, o3, &
                  &lw_up, lw_dn, sw_up, sw_dn)

      use parkind1,                 only : jprb, jprd ! Working/double precision

      use radiation_interface,      only : radiation, set_gas_units
      use radiation_thermodynamics, only : thermodynamics_type
      use radiation_single_level,   only : single_level_type
      use radiation_cloud,          only : cloud_type
      use radiation_config,         only : ISolverSPARTACUS, IGasModelMonochromatic
      use radiation_aerosol,        only : aerosol_type
      use radiation_gas,            only : gas_type, &
                                         & IVolumeMixingRatio, IMassMixingRatio, &
                                         & IH2O, ICO2, IO3, IN2O, ICO, ICH4, IO2, ICFC11, ICFC12, &
                                         & IHCFC22, ICCl4, INO2, GasName, GasLowerCaseName, NMaxGases
      use radiation_flux,           only : flux_type

      integer, intent(in) :: ncol, nlev
      real(kind=jprd), dimension(ncol, nlev+1), intent(in) :: pressure_hl ! pressure (Pa) on half-levels
      real(kind=jprd), dimension(ncol, nlev+1), intent(in) :: temperature_hl ! temperature (K) on half-levels
      real(kind=jprd), intent(in) :: solar_irradiance ! solar irradiance (W m-2)
      real(kind=jprd), intent(in) :: spectral_solar_cycle_multiplier ! +1.0 solar max, -1.0 min, 0.0 mean solar spectrum
      real(kind=jprd), dimension(ncol), intent(in) :: cos_solar_zenith_angle ! cosine of the solar zenith angle
      real(kind=jprd), dimension(ncol, nlev), intent(in) :: cloud_fraction ! cloud fraction
      real(kind=jprd), dimension(ncol, nlev), intent(in) :: fractional_std ! fractional standard deviation of in-cloud water content
      real(kind=jprd), dimension(ncol, nlev), intent(in) :: q_liquid ! liquid specific content (kg/kg)
      real(kind=jprd), dimension(ncol, nlev), intent(in) :: re_liquid ! liquid effective radius (m)
      real(kind=jprd), dimension(ncol, nlev), intent(in) :: q_ice ! ice specific content (kg/kg)
      real(kind=jprd), dimension(ncol, nlev), intent(in) :: re_ice ! ice effective radius (m)
      integer, dimension(ncol), intent(in) :: iseed ! Seed for random number generator in McICA
      real(kind=jprd), dimension(ncol, nlev-1) :: overlap_param ! overlap of cloud boundaries
      real(kind=jprd), dimension(ncol), intent(in) :: skin_temperature ! Ts (K)
      integer, intent(in) :: nalbedobands
      real(kind=jprd), dimension(ncol, nalbedobands), intent(in) :: sw_albedo
      real(kind=jprd), dimension(ncol, nalbedobands), intent(in) :: sw_albedo_direct
      integer, intent(in) :: nemissivitygpoints
      real(kind=jprd), dimension(ncol, nemissivitygpoints), intent(in) :: lw_emissivity
      real(kind=jprd), dimension(ncol, nlev), intent(in) :: q ! vapour mixing ratio
      real(kind=jprd), dimension(ncol, nlev), intent(in) :: o3 ! o3 mixing ratio
      real(kind=jprd), dimension(ncol, nlev+1), intent(out) :: lw_up
      real(kind=jprd), dimension(ncol, nlev+1), intent(out) :: lw_dn
      real(kind=jprd), dimension(ncol, nlev+1), intent(out) :: sw_up
      real(kind=jprd), dimension(ncol, nlev+1), intent(out) :: sw_dn

      type(single_level_type) :: single_level
      type(thermodynamics_type) :: thermodynamics
      type(cloud_type), target  :: cloud
      type(gas_type)            :: gas
      type(aerosol_type)        :: aerosol
      type(flux_type)           :: flux

      character(40)             :: gas_var_name
      integer                   :: jgas

      ! Pressure and temperature (SI units) are on half-levels, i.e. of length (ncol,nlev+1)
      thermodynamics%pressure_hl = pressure_hl
      thermodynamics%temperature_hl = temperature_hl

      ! --------------------------------------------------------
      ! Sun
      ! --------------------------------------------------------

      single_level%solar_irradiance = solar_irradiance
      single_level%spectral_solar_cycle_multiplier = spectral_solar_cycle_multiplier
      single_level%cos_sza = cos_solar_zenith_angle

      ! --------------------------------------------------------
      ! Cloud 
      ! --------------------------------------------------------

      if (config%do_clouds) then
        cloud%fraction = cloud_fraction
        cloud%fractional_std = fractional_std
        cloud%ntype=2
        allocate(cloud%mixing_ratio(ncol, nlev, cloud%ntype))
        allocate(cloud%effective_radius(ncol, nlev, cloud%ntype))

        cloud%mixing_ratio(:,:,1) = q_liquid
        cloud%mixing_ratio(:,:,2) = q_ice
        cloud%effective_radius(:,:,1) = re_liquid
        cloud%effective_radius(:,:,2) = re_ice
        cloud%q_liq  => cloud%mixing_ratio(:,:,1)
        cloud%q_ice  => cloud%mixing_ratio(:,:,2)
        cloud%re_liq => cloud%effective_radius(:,:,1)
        cloud%re_ice => cloud%effective_radius(:,:,2)

        call single_level%init_seed_simple(1,ncol) ! Simple initialization of the seeds for the Monte Carlo scheme
        single_level%iseed = iseed

        cloud%overlap_param = overlap_param

        ! --------------------------------------------------------
        ! Cloud properties needed by SPARTACUS
        ! --------------------------------------------------------
        if (config%i_solver_sw == ISolverSPARTACUS .or. config%i_solver_lw == ISolverSPARTACUS) then
          print*, '*** Error: SPARTACUS solver not implemented'
          stop
        endif ! Spartacus
      endif ! Cloud

      ! --------------------------------------------------------
      ! Surface properties
      ! --------------------------------------------------------

      single_level%is_simple_surface = .true.

      single_level%skin_temperature = skin_temperature
      single_level%sw_albedo = sw_albedo
      single_level%sw_albedo_direct = sw_albedo_direct
      single_level%lw_emissivity = lw_emissivity

      ! --------------------------------------------------------
      ! Aerosol and gas concentrations
      ! --------------------------------------------------------

      if (config%use_aerosols) then
        print*, '*** Error: aerosols not implemented'
        call flush()
        stop
      endif

      ! Water vapour and ozone are always in terms of mass mixing ratio
      ! (kg/kg) and always 2D arrays with dimensions (ncol,nlev), unlike
      ! other gases (see below)

      call gas%allocate(ncol, nlev)

      ! Loop through all radiatively important gases
      do jgas = 1, NMaxGases
        if (jgas == IH2O) then
          call gas%put(IH2O, IMassMixingRatio, q)
        else if (jgas == IO3) then
          call gas%put(IO3, IMassMixingRatio, o3)
        else
          gas_var_name = trim(GasLowerCaseName(jgas))
          !irank = file%get_rank(trim(gas_var_name))
          !if (irank == 0) then
          !  ! Store this as a well-mixed gas
          !  call gas%put_well_mixed(jgas, IVolumeMixingRatio, well_mixed_gas_vmr)
          !else if (irank == 2) then
          !  call gas%put(jgas, IVolumeMixingRatio, gas_mr)
          !else if (irank > 0) then
          !  write(nulout,'(a,a,a)')  '***  Error: ', trim(gas_var_name), ' does not have 0 or 2 dimensions'
          !  stop
          !end if
        end if
      end do

      ! --------------------------------------------------------
      ! Call radiation scheme
      ! --------------------------------------------------------
    
      ! Ensure the units of the gas mixing ratios are what is required
      ! by the gas absorption model
      call set_gas_units(config, gas)
    
      ! Compute saturation with respect to liquid (needed for aerosol hydration) call...
      call thermodynamics%calc_saturation_wrt_liquid(1, ncol)
    
      ! Allocate memory for the flux profiles, which may include arrays
      ! of dimension n_bands_sw/n_bands_lw, so must be called after
      ! setup_radiation
      call flux%allocate(config, 1, ncol, nlev)
        
      ! Call the ECRAD radiation scheme
      call radiation(ncol, nlev, 1, ncol, &
           &  config, single_level, thermodynamics, gas, cloud, aerosol, flux)
    
      ! --------------------------------------------------------
      ! Output
      ! --------------------------------------------------------
      if (config%i_gas_model_lw == IGasModelMonochromatic .and. &
         &config%mono_lw_wavelength > 0.0_jprb) then
          print*, '*** Error: lw flux unit is W m-3'
          call flush()
          stop
      endif
    
      if (config%do_lw) then
        lw_up = flux%lw_up
        lw_dn = flux%lw_dn
      else
        lw_up = 0.
        lw_dn = 0.
      endif
      if (config%do_sw) then
         sw_up = flux%sw_up
         sw_dn = flux%sw_dn
      else
        sw_up = 0.
        sw_dn = 0.
      endif
      if (config%do_clouds) then
        deallocate(cloud%mixing_ratio)
        deallocate(cloud%effective_radius)
      endif
      call gas%deallocate()
      call flux%deallocate()
    
    end subroutine run
end module ecrad4py
