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

    subroutine setup(namelist_file_name, directory_name) bind(C, name='setup_')
      use radiation_interface,      only : setup_radiation
      use, intrinsic :: iso_c_binding, only: c_char
      character(kind=c_char), dimension(512), intent(in) :: namelist_file_name
      character(kind=c_char), dimension(511), intent(in) :: directory_name
      character(len=size(namelist_file_name)) :: string_namelist_file_name
      character(len=size(directory_name)) :: string_directory_name
      integer :: jk

      do jk=1, size(namelist_file_name)
        string_namelist_file_name(jk:jk)=namelist_file_name(jk)
      enddo
      do jk=1, size(directory_name)
        string_directory_name(jk:jk)=directory_name(jk)
      enddo

      ! Read "radiation" namelist into radiation configuration type
      call config%read(file_name=string_namelist_file_name)
      config%directory_name = string_directory_name

      ! Setup the radiation scheme: load the coefficients for gas and
      ! cloud optics, currently from RRTMG
      call setup_radiation(config)
    end subroutine setup

    function get_IVolumeMixingRatio() bind(C, name='get_IVolumeMixingRatio_')
      use radiation_gas,               only : IVolumeMixingRatio
      use, intrinsic :: iso_c_binding, only: c_int64_t
      integer(kind=c_int64_t) :: get_IVolumeMixingRatio
      get_IVolumeMixingRatio=int(IVolumeMixingRatio, kind(get_IVolumeMixingRatio))
    end function get_IVolumeMixingRatio

    function get_IMassMixingRatio() bind(C, name='get_IMassMixingRatio_')
      use radiation_gas,               only : IMassMixingRatio
      use, intrinsic :: iso_c_binding, only: c_int64_t
      integer(kind=c_int64_t) :: get_IMassMixingRatio
      get_IMassMixingRatio=int(IMassMixingRatio, kind(get_IMassMixingRatio))
    end function get_IMassMixingRatio

    subroutine run(ncol, nlev, nblocksize, pressure_hl, temperature_hl, solar_irradiance, &
                  &spectral_solar_cycle_multiplier, &
                  &cos_solar_zenith_angle, cloud_fraction, fractional_std, &
                  &q_liquid, re_liquid, q_ice, re_ice, iseed, overlap_param, &
                  &skin_temperature, nalbedobands, sw_albedo, sw_albedo_direct, &
                  &nemissivitygpoints, lw_emissivity, &
                  &q_unit, q, co2_unit, co2, o3_unit, o3, n2o_unit, n2o, &
                  &co_unit, co, ch4_unit, ch4, o2_unit, o2, cfc11_unit, cfc11, &
                  &cfc12_unit, cfc12, hcfc22_unit, hcfc22, ccl4_unit, ccl4, no2_unit, no2, &
                  &naerosols, aerosols, inv_cloud_effective_size, inv_inhom_effective_size, &
                  &lw_up, lw_dn, lw_up_clear, lw_dn_clear, cloud_cover_lw, &
                  &sw_up, sw_dn, sw_up_clear, sw_dn_clear, cloud_cover_sw) bind(C, name='run_')

      use, intrinsic :: iso_c_binding, only: c_int64_t, c_double
      use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
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
      use radiation_io,             only : radiation_abort

      integer(kind=c_int64_t), intent(in) :: ncol, nlev, nblocksize
      real(kind=c_double), dimension(ncol, nlev+1), intent(in) :: pressure_hl ! pressure (Pa) on half-levels
      real(kind=c_double), dimension(ncol, nlev+1), intent(in) :: temperature_hl ! temperature (K) on half-levels
      real(kind=c_double), intent(in) :: solar_irradiance ! solar irradiance (W m-2)
      real(kind=c_double), intent(in) :: spectral_solar_cycle_multiplier ! +1.0 solar max, -1.0 min, 0.0 mean solar spectrum
      real(kind=c_double), dimension(ncol), intent(in) :: cos_solar_zenith_angle ! cosine of the solar zenith angle
      real(kind=c_double), dimension(ncol, nlev), intent(in) :: cloud_fraction ! cloud fraction
      real(kind=c_double), dimension(ncol, nlev), intent(in) :: fractional_std ! fractional standard deviation of in-cloud water content
      real(kind=c_double), dimension(ncol, nlev), intent(in) :: q_liquid ! liquid specific content (kg/kg)
      real(kind=c_double), dimension(ncol, nlev), intent(in) :: re_liquid ! liquid effective radius (m)
      real(kind=c_double), dimension(ncol, nlev), intent(in) :: q_ice ! ice specific content (kg/kg)
      real(kind=c_double), dimension(ncol, nlev), intent(in) :: re_ice ! ice effective radius (m)
      integer(kind=c_int64_t), dimension(ncol), intent(in) :: iseed ! Seed for random number generator in McICA
      real(kind=c_double), dimension(ncol, nlev-1) :: overlap_param ! overlap of cloud boundaries
      real(kind=c_double), dimension(ncol), intent(in) :: skin_temperature ! Ts (K)
      integer(kind=c_int64_t), intent(in) :: nalbedobands
      real(kind=c_double), dimension(ncol, nalbedobands), intent(in) :: sw_albedo
      real(kind=c_double), dimension(ncol, nalbedobands), optional, intent(in) :: sw_albedo_direct
      integer(kind=c_int64_t), intent(in) :: nemissivitygpoints
      real(kind=c_double), dimension(ncol, nemissivitygpoints), intent(in) :: lw_emissivity
      integer(kind=c_int64_t), optional, intent(in) :: q_unit ! vapour unit: mass or volume mixing ratio
      real(kind=c_double), dimension(ncol, nlev), optional, intent(in) :: q ! vapour mass/volume mixing ratio value
      integer(kind=c_int64_t), optional, intent(in) :: co2_unit
      real(kind=c_double), dimension(ncol, nlev), optional, intent(in) ::co2
      integer(kind=c_int64_t), optional, intent(in) :: o3_unit
      real(kind=c_double), dimension(ncol, nlev), optional, intent(in) :: o3
      integer(kind=c_int64_t), optional, intent(in) :: n2o_unit
      real(kind=c_double), dimension(ncol, nlev), optional, intent(in) :: n2o
      integer(kind=c_int64_t), optional, intent(in) :: co_unit
      real(kind=c_double), dimension(ncol, nlev), optional, intent(in) :: co
      integer(kind=c_int64_t), optional, intent(in) :: ch4_unit
      real(kind=c_double), dimension(ncol, nlev), optional, intent(in) :: ch4
      integer(kind=c_int64_t), optional, intent(in) :: o2_unit
      real(kind=c_double), dimension(ncol, nlev), optional, intent(in) :: o2
      integer(kind=c_int64_t), optional, intent(in) :: cfc11_unit
      real(kind=c_double), dimension(ncol, nlev), optional, intent(in) :: cfc11
      integer(kind=c_int64_t), optional, intent(in) :: cfc12_unit
      real(kind=c_double), dimension(ncol, nlev), optional, intent(in) :: cfc12
      integer(kind=c_int64_t), optional, intent(in) :: hcfc22_unit
      real(kind=c_double), dimension(ncol, nlev), optional, intent(in) :: hcfc22
      integer(kind=c_int64_t), optional, intent(in) :: ccl4_unit
      real(kind=c_double), dimension(ncol, nlev), optional, intent(in) :: ccl4
      integer(kind=c_int64_t), optional, intent(in) :: no2_unit
      real(kind=c_double), dimension(ncol, nlev), optional, intent(in) :: no2
      integer(kind=c_int64_t), intent(in) :: naerosols
      real(kind=c_double), dimension(ncol, nlev, naerosols), optional, intent(in) :: aerosols
      real(kind=c_double), dimension(ncol, nlev), optional, intent(in) :: inv_cloud_effective_size
      real(kind=c_double), dimension(ncol, nlev), optional, intent(in) :: inv_inhom_effective_size
      real(kind=c_double), dimension(ncol, nlev+1), intent(out) :: lw_up
      real(kind=c_double), dimension(ncol, nlev+1), intent(out) :: lw_dn
      real(kind=c_double), dimension(ncol, nlev+1), intent(out) :: lw_up_clear
      real(kind=c_double), dimension(ncol, nlev+1), intent(out) :: lw_dn_clear
      real(kind=c_double), dimension(ncol), intent(out) :: cloud_cover_lw
      real(kind=c_double), dimension(ncol, nlev+1), intent(out) :: sw_up
      real(kind=c_double), dimension(ncol, nlev+1), intent(out) :: sw_dn
      real(kind=c_double), dimension(ncol, nlev+1), intent(out) :: sw_up_clear
      real(kind=c_double), dimension(ncol, nlev+1), intent(out) :: sw_dn_clear
      real(kind=c_double), dimension(ncol), intent(out) :: cloud_cover_sw

      type(single_level_type) :: single_level
      type(thermodynamics_type) :: thermodynamics
      type(cloud_type), target  :: cloud
      type(gas_type)            :: gas
      type(aerosol_type)        :: aerosol
      type(flux_type)           :: flux

      character(40)             :: gas_var_name
      integer                   :: jgas
      integer                   :: jblock, nblock, istartcol, iendcol

      ! Pressure and temperature (SI units) are on half-levels, i.e. of length (ncol, nlev+1)
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
        allocate(cloud%mixing_ratio(int(ncol), nlev, cloud%ntype))
        allocate(cloud%effective_radius(int(ncol), nlev, cloud%ntype))

        cloud%mixing_ratio(:,:,1) = q_liquid
        cloud%mixing_ratio(:,:,2) = q_ice
        cloud%effective_radius(:,:,1) = re_liquid
        cloud%effective_radius(:,:,2) = re_ice
        cloud%q_liq  => cloud%mixing_ratio(:,:,1)
        cloud%q_ice  => cloud%mixing_ratio(:,:,2)
        cloud%re_liq => cloud%effective_radius(:,:,1)
        cloud%re_ice => cloud%effective_radius(:,:,2)

        call single_level%init_seed_simple(1, int(ncol)) ! Simple initialization of the seeds for the Monte Carlo scheme
        single_level%iseed = int(iseed)

        cloud%overlap_param = overlap_param

        ! --------------------------------------------------------
        ! Cloud properties needed by SPARTACUS
        ! --------------------------------------------------------
        if (config%i_solver_sw == ISolverSPARTACUS .or. config%i_solver_lw == ISolverSPARTACUS) then
          if (present(inv_cloud_effective_size)) then
            cloud%inv_cloud_effective_size = inv_cloud_effective_size
            if (present(inv_inhom_effective_size)) then
              cloud%inv_inhom_effective_size = inv_inhom_effective_size
            endif
          else
            call radiation_abort('inv_cloud_effective_size array absent with SPARTACUS solver')
          endif
        endif ! Spartacus
      endif ! Cloud

      ! --------------------------------------------------------
      ! Surface properties
      ! --------------------------------------------------------

      single_level%is_simple_surface = .true.

      single_level%skin_temperature = skin_temperature
      single_level%sw_albedo = sw_albedo
      if(present(sw_albedo_direct)) then
        single_level%sw_albedo_direct = sw_albedo_direct
      endif
      single_level%lw_emissivity = lw_emissivity

      ! --------------------------------------------------------
      ! Aerosol and gas concentrations
      ! --------------------------------------------------------

      if (config%use_aerosols) then
        if (present(aerosols)) then
          aerosol%mixing_ratio = aerosols
        else
          call radiation_abort('aerosols array absent with config%use_aerosols==.true.')
        endif
      endif

      ! Water vapour and ozone are always in terms of mass mixing ratio
      ! (kg/kg) and always 2D arrays with dimensions (ncol, nlev), unlike
      ! other gases (see below)

      call gas%allocate(int(ncol), int(nlev))

      ! Loop through all radiatively important gases
      do jgas = 1, NMaxGases
        if (jgas == IH2O .and. present(q_unit) .and. present(q)) then
          call gas%put(IH2O, int(q_unit), q)
        else if (jgas == ICO2 .and. present(co2_unit) .and. present(co2)) then
          call gas%put(ICO2, int(co2_unit), co2)
        else if (jgas == IO3 .and. present(o3_unit) .and. present(o3)) then
          call gas%put(IO3, int(o3_unit), o3)
        else if (jgas == IN2O .and. present(n2o_unit) .and. present(n2o)) then
          call gas%put(IN2O, int(n2o_unit), n2o)
        else if (jgas == ICO .and. present(co_unit) .and. present(co)) then
          call gas%put(ICO, int(co_unit), co)
        else if (jgas == ICH4 .and. present(ch4_unit) .and. present(ch4)) then
          call gas%put(ICH4, int(ch4_unit), ch4)
        else if (jgas == IO2 .and. present(o2_unit) .and. present(o2)) then
          call gas%put(IO2, int(o2_unit), o2)
        else if (jgas == ICFC11 .and. present(cfc11_unit) .and. present(cfc11)) then
          call gas%put(ICFC11, int(cfc11_unit), cfc11)
        else if (jgas == ICFC12 .and. present(cfc12_unit) .and. present(cfc12)) then
          call gas%put(ICFC12, int(cfc12_unit), cfc12)
        else if (jgas == IHCFC22 .and. present(hcfc22_unit) .and. present(hcfc22)) then
          call gas%put(IHCFC22, int(hcfc22_unit), hcfc22)
        else if (jgas == ICCl4 .and. present(ccl4_unit) .and. present(ccl4)) then
          call gas%put(ICCl4, int(ccl4_unit), ccl4)
        else if (jgas == INO2 .and. present(no2_unit) .and. present(no2)) then
          call gas%put(INO2, int(no2_unit), no2)
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
      call thermodynamics%calc_saturation_wrt_liquid(1, int(ncol))
    
      ! Allocate memory for the flux profiles, which may include arrays
      ! of dimension n_bands_sw/n_bands_lw, so must be called after
      ! setup_radiation
      call flux%allocate(config, 1, int(ncol), int(nlev))
        
      ! Call the ECRAD radiation scheme
      if (nblocksize <= 0) then
        call radiation(int(ncol), int(nlev), 1, int(ncol), &                              
             &  config, single_level, thermodynamics, gas, cloud, aerosol, flux)
      else
        nblock = (int(ncol) - 1 + int(nblocksize)) / int(nblocksize)
        !$OMP PARALLEL DO PRIVATE(istartcol, iendcol) SCHEDULE(RUNTIME)
        do jblock = 1, nblock
          istartcol = (jblock - 1) * int(nblocksize) + 1
          iendcol = min(istartcol + int(nblocksize) - 1, int(ncol))
          call radiation(int(ncol), int(nlev), istartcol, iendcol, &
               &  config, single_level, thermodynamics, gas, cloud, aerosol, flux)
        enddo
      endif
    
      ! --------------------------------------------------------
      ! Output
      ! --------------------------------------------------------
      if (config%i_gas_model_lw == IGasModelMonochromatic .and. &
         &config%mono_lw_wavelength > 0.0_jprb) then
          call radiation_abort('lw flux unit is W m-3')
      endif
    
      if (config%do_lw) then
        lw_up = flux%lw_up
        lw_dn = flux%lw_dn
        if (config%do_clear) then
          lw_up_clear = flux%lw_up_clear
          lw_dn_clear = flux%lw_dn_clear
        else
          lw_up_clear = ieee_value(lw_up_clear, ieee_quiet_nan)
          lw_dn_clear = ieee_value(lw_dn_clear, ieee_quiet_nan)
        endif
        if (config%do_clouds) then
          cloud_cover_lw = flux%cloud_cover_lw
        else
          cloud_cover_lw = ieee_value(cloud_cover_lw, ieee_quiet_nan)
        endif
      else
        lw_up = ieee_value(lw_up, ieee_quiet_nan)
        lw_dn = ieee_value(lw_dn, ieee_quiet_nan)
        lw_up_clear = ieee_value(lw_up_clear, ieee_quiet_nan)
        lw_dn_clear = ieee_value(lw_dn_clear, ieee_quiet_nan)
        cloud_cover_lw = ieee_value(cloud_cover_lw, ieee_quiet_nan)
      endif
      if (config%do_sw) then
         sw_up = flux%sw_up
         sw_dn = flux%sw_dn
         if (config%do_clear) then
           sw_up_clear = flux%sw_up_clear
           sw_dn_clear = flux%sw_dn_clear
         else
           sw_up_clear = ieee_value(sw_up_clear, ieee_quiet_nan)
           sw_dn_clear = ieee_value(sw_dn_clear, ieee_quiet_nan)
         endif
         if (config%do_clouds) then
          cloud_cover_sw = flux%cloud_cover_sw
        else
          cloud_cover_sw = ieee_value(cloud_cover_sw, ieee_quiet_nan)
        endif
      else
        sw_up = ieee_value(sw_up, ieee_quiet_nan)
        sw_dn = ieee_value(sw_dn, ieee_quiet_nan)
        sw_up_clear = ieee_value(sw_up_clear, ieee_quiet_nan)
        sw_dn_clear = ieee_value(sw_dn_clear, ieee_quiet_nan)
        cloud_cover_sw = ieee_value(cloud_cover_sw, ieee_quiet_nan)
      endif
      if (config%do_clouds) then
        deallocate(cloud%mixing_ratio)
        deallocate(cloud%effective_radius)
      endif
      call gas%deallocate()
      call flux%deallocate()
    
    end subroutine run
end module ecrad4py
