! ecrad_ifs_driver_blocked.F90 - Driver for offline ECRAD radiation scheme
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
! 1) Namelist file to configure the radiation calculation, but note
!    that only the radiation_config group is read
! 2) Name of a NetCDF file containing one or more atmospheric profiles
! 3) Name of output NetCDF file
!
! This version uses the infrastructure of the IFS, such as computing
! effective radius and cloud overlap from latitude and other
! variables. To configure ecRad in this version you need to edit
! ifs/yoerad.F90 in the ecRad package, but these options can be
! overridden with the "radiation" namelist. This file requires the
! input data to have compatible settings, e.g. the right number of
! aerosol variables, and surface albedo/emissivity bands; a test file
! satisfying this requirement is test/ifs/ecrad_meridian.nc in the
! ecRad package.
!
! Note that the purpose of this file is simply to demonstrate the use
! of the setup_radiation_scheme and radiation_scheme routines as well
! as the use of a blocked memory layout to improve cache efficiency;
! all the rest is using the offline ecRad driver containers to read
! a NetCDF file to memory and pass it into these routines.

program ecrad_ifs_driver

  ! --------------------------------------------------------
  ! Section 1: Declarations
  ! --------------------------------------------------------
  use parkind1,                 only : jprb, jprd ! Working/double precision
  use yomhook,                  only : dr_hook_init
#ifdef HAVE_FIAT
  use mpl_module,               only : mpl_init, mpl_end
#endif

  use radiation_io,             only : nulout
  use radiation_single_level,   only : single_level_type
  use radiation_thermodynamics, only : thermodynamics_type

  use radiation_gas,            only : gas_type, IMassMixingRatio, &
       &   IH2O, ICO2, IO3, IN2O, INO2, ICO, ICH4, IO2, ICFC11, ICFC12, &
       &   IHCFC22, ICCl4
  use radiation_cloud,          only : cloud_type
  use radiation_aerosol,        only : aerosol_type
  use radiation_flux,           only : flux_type
  use radiation_save,           only : save_net_fluxes
  use radiation_setup,          only : tradiation, setup_radiation_scheme
  use radiation_constants,      only : Pi
  use ecrad_driver_config,      only : driver_config_type
  use ecrad_driver_read_input,  only : read_input
  use easy_netcdf
  use ifs_blocking
  use iso_fortran_env,          only : int64
#ifdef HAVE_NVTX
  use nvtx
#endif
#ifdef HAVE_ROCTX
  use roctx_profiling, only: roctxstartrange, roctxendrange
  use iso_c_binding, only: c_null_char
#endif

  implicit none

#include "radiation_scheme.intfb.h"

  integer(kind=int64)    :: count_rate,t(4)

  ! The NetCDF file containing the input profiles
  type(netcdf_file)         :: file

  ! Configuration for the radiation scheme, IFS style
  type(tradiation)          :: yradiation

  ! Derived types for the inputs to the radiation scheme
  type(single_level_type)   :: single_level
  type(thermodynamics_type) :: thermodynamics
  type(gas_type)            :: gas
  type(cloud_type)          :: cloud
  type(aerosol_type)        :: aerosol

  ! Configuration specific to this driver
  type(driver_config_type)  :: driver_config

  ! Derived type containing outputs from the radiation scheme
  type(flux_type)           :: flux

  ! Additional arrays passed to radiation_scheme
  real(jprb), allocatable, dimension(:) :: ccn_land, ccn_sea, sin_latitude, longitude_rad, land_frac
  real(jprb), allocatable, dimension(:,:) :: pressure_fl, temperature_fl
  real(jprb), allocatable, dimension(:) :: flux_sw_direct_normal, flux_uv, flux_par, flux_par_clear, &
       &  emissivity_out
  real(jprb), allocatable, dimension(:,:) :: flux_diffuse_band, flux_direct_band

  ! Bespoke data types to set-up the blocked memory layout
  type(ifs_config_type)        :: ifs_config
  real(kind=jprb), allocatable :: zrgp(:,:,:) ! monolithic IFS data structure  
#ifdef BITIDENTITY_TESTING
  integer, allocatable         :: iseed(:,:) ! Seed for random number generator
  integer, allocatable         :: iseed_single_block(:) ! Seed for random number generator
#endif

  integer :: ncol, nlev         ! Number of columns and levels
  integer :: nproma             ! block size

  ! Name of file names specified on command line
  character(len=512) :: file_name
  integer            :: istatus ! Result of command_argument_count

  real(kind=jprd) :: total_dt, kernel_dt

  ! For demonstration of get_sw_weights later on
  ! Ultraviolet weightings
  !integer    :: nweight_uv
  !integer    :: iband_uv(100)
  !real(jprb) :: weight_uv(100)
  ! Photosynthetically active radiation weightings
  !integer    :: nweight_par
  !integer    :: iband_par(100)
  !real(jprb) :: weight_par(100)

  ! Loop index for repeats (for benchmarking)
  integer :: jrepeat

  ! Loop index
  integer :: jrl, ibeg, iend, il, ib
  integer :: ngpblks
  
  ! Are any variables out of bounds?
  logical :: is_out_of_bounds

  real(kind=jprb), allocatable :: PMU0(:)                        !(KLON, ) ! Cosine of solar zenith ang
  real(kind=jprb), allocatable :: PTEMPERATURE_SKIN(:)           !(KLON) ! (K)
  real(kind=jprb), allocatable :: PALBEDO_DIF(:,:)               !(KLON,YRADIATION%YRERAD%NSW)
  real(kind=jprb), allocatable :: PALBEDO_DIR(:,:)               !(KLON,YRADIATION%YRERAD%NSW)
  real(kind=jprb), allocatable :: PSPECTRALEMISS(:,:)            !(KLON,YRADIATION%YRERAD%NLWEMISS)
  real(kind=jprb), allocatable :: PCCN_LAND(:)!(KLON)
  real(kind=jprb), allocatable :: PCCN_SEA(:)!(KLON)
  real(kind=jprb), allocatable :: PGELAM(:)                      !(KLON)
  real(kind=jprb), allocatable :: PGEMU(:)                       !(KLON)
  real(kind=jprb), allocatable :: PLAND_SEA_MASK(:)              !((KLON)
  real(kind=jprb), allocatable :: PPRESSURE(:,:)                 !((KLON,KLEV)    ! (Pa)
  real(kind=jprb), allocatable :: PTEMPERATURE(:,:)              !(KLON,KLEV) ! (K)
  real(kind=jprb), allocatable :: PPRESSURE_H(:,:)               !(KLON,KLEV+1)    ! (Pa)
  real(kind=jprb), allocatable :: PTEMPERATURE_H(:,:)            !(KLON,KLEV+1) ! (K)
  real(kind=jprb), allocatable :: PQ(:,:)!(KLON,KLEV)
  real(kind=jprb), allocatable :: PCO2(:,:)!(KLON,KLEV)
  real(kind=jprb), allocatable :: PCH4(:,:)!(KLON,KLEV)
  real(kind=jprb), allocatable :: PN2O(:,:)!(KLON,KLEV)
  real(kind=jprb), allocatable :: PNO2(:,:)!(KLON,KLEV)
  real(kind=jprb), allocatable :: PCFC11(:,:)!(KLON,KLEV)
  real(kind=jprb), allocatable :: PCFC12(:,:)!(KLON,KLEV)
  real(kind=jprb), allocatable :: PHCFC22(:,:)!(KLON,KLEV)
  real(kind=jprb), allocatable :: PCCL4(:,:)!(KLON,KLEV)
  real(kind=jprb), allocatable :: PO3(:,:)!(KLON,KLEV)
  real(kind=jprb), allocatable :: PCLOUD_FRAC(:,:)!(KLON,KLEV)
  real(kind=jprb), allocatable :: PQ_LIQUID(:,:)!(KLON,KLEV)
  real(kind=jprb), allocatable :: PQ_ICE(:,:)!(KLON,KLEV)
  real(kind=jprb), allocatable :: PQ_RAIN(:,:)!(KLON,KLEV)
  real(kind=jprb), allocatable :: PQ_SNOW(:,:)!(KLON,KLEV)
  real(kind=jprb), allocatable :: PAEROSOL_OLD(:,:,:)!(KLON,6,KLEV)
  real(kind=jprb), allocatable :: PAEROSOL(:,:,:)!(KLON,KLEV,KAEROSOL)
  ! OUT
  real(kind=jprb), allocatable :: PFLUX_SW(:,:)!(KLON,KLEV+1)
  real(kind=jprb), allocatable :: PFLUX_LW(:,:)!(KLON,KLEV+1)
  real(kind=jprb), allocatable :: PFLUX_SW_CLEAR(:,:)!(KLON,KLEV+1)
  real(kind=jprb), allocatable :: PFLUX_LW_CLEAR(:,:)!(KLON,KLEV+1)
  real(kind=jprb), allocatable :: PFLUX_SW_DN(:)!(KLON)
  real(kind=jprb), allocatable :: PFLUX_LW_DN(:)!(KLON)
  real(kind=jprb), allocatable :: PFLUX_SW_DN_CLEAR(:)!(KLON)
  real(kind=jprb), allocatable :: PFLUX_LW_DN_CLEAR(:)!(KLON)
  real(kind=jprb), allocatable :: PFLUX_DIR(:)!(KLON)
  real(kind=jprb), allocatable :: PFLUX_DIR_CLEAR(:)!(KLON)
  real(kind=jprb), allocatable :: PFLUX_DIR_INTO_SUN(:)!(KLON)
  real(kind=jprb), allocatable :: PFLUX_UV(:)!(KLON)
  real(kind=jprb), allocatable :: PFLUX_PAR(:)!(KLON)
  real(kind=jprb), allocatable :: PFLUX_PAR_CLEAR(:)!(KLON)
  real(kind=jprb), allocatable :: PFLUX_SW_DN_TOA(:)!(KLON)
  real(kind=jprb), allocatable :: PEMIS_OUT(:)!(KLON)
  real(kind=jprb), allocatable :: PLWDERIVATIVE(:,:)!(KLON,KLEV+1)
  real(kind=jprb), allocatable :: PSWDIFFUSEBAND(:,:)!(KLON,YRADIATION%YRERAD%NSW)
  real(kind=jprb), allocatable :: PSWDIRECTBAND(:,:)!(KLON,YRADIATION%YRERAD%NSW)
  real(kind=jprb), allocatable :: PRE_LIQ(:,:)!(KLON, KLEV)
  real(kind=jprb), allocatable :: PRE_ICE(:,:)!(KLON, KLEV)
  real(kind=jprb), allocatable :: PCLOUD_OVERLAP(:,:)!(KLON, KLEV-1)

!  integer    :: iband(20), nweights
!  real(jprb) :: weight(20)


  call system_clock(count_rate=count_rate)

  ! Initialise MPI if not done yet
#ifdef HAVE_FIAT
  call mpl_init
#endif

  call dr_hook_init()

  ! --------------------------------------------------------
  ! Section 2: Configure
  ! --------------------------------------------------------

  ! Check program called with correct number of arguments
  if (command_argument_count() < 3) then
    stop 'Usage: ecrad config.nam input_file.nc output_file.nc'
  end if

  ! Use namelist to configure the radiation calculation
  call get_command_argument(1, file_name, status=istatus)
  if (istatus /= 0) then
    stop 'Failed to read name of namelist file as string of length < 512'
  end if

  ! Read "radiation_driver" namelist into radiation driver config type
  call driver_config%read(file_name)
  nproma = driver_config%nblocksize

  if (driver_config%iverbose >= 2) then
    write(nulout,'(a)') '-------------------------- OFFLINE ECRAD RADIATION SCHEME --------------------------'
    write(nulout,'(a)') 'Copyright (C) 2014- ECMWF'
    write(nulout,'(a)') 'Contact: Robin Hogan (r.j.hogan@ecmwf.int)'
#ifdef PARKIND1_SINGLE
    write(nulout,'(a)') 'Floating-point precision: single'
#else
    write(nulout,'(a)') 'Floating-point precision: double'
#endif
  end if

  ! Albedo/emissivity intervals may be specified like this
  !call config%define_sw_albedo_intervals(6, &
  !     &  [0.25e-6_jprb, 0.44e-6_jprb, 0.69e-6_jprb, &
  !     &     1.19_jprb, 2.38e-6_jprb], [1,2,3,4,5,6], &
  !     &   do_nearest=.false.)
  !call config%define_lw_emiss_intervals(3, &
  !     &  [8.0e-6_jprb, 13.0e-6_jprb], [1,2,1], &
  !     &   do_nearest=.false.)

  ! If monochromatic aerosol properties are required, then the
  ! wavelengths can be specified (in metres) as follows - these can be
  ! whatever you like for the general aerosol optics, but must match
  ! the monochromatic values in the aerosol input file for the older
  ! aerosol optics
  !call config%set_aerosol_wavelength_mono( &
  !     &  [3.4e-07_jprb, 3.55e-07_jprb, 3.8e-07_jprb, 4.0e-07_jprb, 4.4e-07_jprb, &
  !     &   4.69e-07_jprb, 5.0e-07_jprb, 5.32e-07_jprb, 5.5e-07_jprb, 6.45e-07_jprb, &
  !     &   6.7e-07_jprb, 8.0e-07_jprb, 8.58e-07_jprb, 8.65e-07_jprb, 1.02e-06_jprb, &
  !     &   1.064e-06_jprb, 1.24e-06_jprb, 1.64e-06_jprb, 2.13e-06_jprb, 1.0e-05_jprb])

  call yradiation%rad_config%read(file_name=file_name)

  ! Setup aerosols
  if (yradiation%rad_config%use_aerosols) then
    yradiation%yrerad%naermacc = 1 ! MACC-derived aerosol climatology on a NMCLAT x NMCLON grid
  else
    yradiation%yrerad%naermacc = 0
  endif

  ! Setup the radiation scheme: load the coefficients for gas and
  ! cloud optics, currently from RRTMG
  call setup_radiation_scheme(yradiation, .true., file_name=file_name)
  ! Or call without specifying the namelist filename, in which case
  ! the default settings are from yoerad.F90
  !call setup_radiation_scheme(yradiation, .true.)

  ! Demonstration of how to get weights for UV and PAR fluxes
  !if (config%do_sw) then
  !  call config%get_sw_weights(0.2e-6_jprb, 0.4415e-6_jprb,&
  !       &  nweight_uv, iband_uv, weight_uv,&
  !       &  'ultraviolet')
  !  call config%get_sw_weights(0.4e-6_jprb, 0.7e-6_jprb,&
  !       &  nweight_par, iband_par, weight_par,&
  !       &  'photosynthetically active radiation, PAR')
  !end if

  ! --------------------------------------------------------
  ! Section 3: Read input data file
  ! --------------------------------------------------------

  ! Get NetCDF input file name
  call get_command_argument(2, file_name, status=istatus)
  if (istatus /= 0) then
    stop 'Failed to read name of input NetCDF file as string of length < 512'
  end if

  ! Open the file and configure the way it is read
  call file%open(trim(file_name), iverbose=driver_config%iverbose)

  ! Get NetCDF output file name
  call get_command_argument(3, file_name, status=istatus)
  if (istatus /= 0) then
    stop 'Failed to read name of output NetCDF file as string of length < 512'
  end if

  ! 2D arrays are assumed to be stored in the file with height varying
  ! more rapidly than column index. Specifying "true" here transposes
  ! all 2D arrays so that the column index varies fastest within the
  ! program.
  call file%transpose_matrices(.true.)

  ! Read input variables from NetCDF file, noting that cloud overlap
  ! and effective radius are ignored
  call read_input(file, yradiation%rad_config, driver_config, ncol, nlev, &
       &          single_level, thermodynamics, &
       &          gas, cloud, aerosol)

  ! Latitude is used for cloud overlap and ice effective radius
  if (file%exists('lat')) then
    call file%get('lat', sin_latitude)
    sin_latitude = sin(sin_latitude * Pi/180.0_jprb)
  else
    allocate(sin_latitude(ncol))
    sin_latitude = 0.0_jprb
  end if

  if (file%exists('lon')) then
    call file%get('lon', longitude_rad)
    longitude_rad = longitude_rad * Pi/180.0_jprb
  else
    allocate(longitude_rad(ncol))
    longitude_rad = 0.0_jprb
  end if

  ! Close input file
  call file%close()

  ! Convert gas units to mass-mixing ratio
  call gas%set_units(gas, IMassMixingRatio, lacc=.false.)

  ! Compute seed from skin temperature residual
  !  single_level%iseed = int(1.0e9*(single_level%skin_temperature &
  !       &                            -int(single_level%skin_temperature)))

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

  ! --------------------------------------------------------
  ! Section 4: Call radiation scheme
  ! --------------------------------------------------------

  ! Compute saturation with respect to liquid (needed for aerosol
  ! hydration) call
  call thermodynamics%calc_saturation_wrt_liquid(thermodynamics, driver_config%istartcol,driver_config%iendcol)

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
  call flux%allocate(yradiation%rad_config, 1, ncol, nlev)

  ! set relevant fluxes to zero
  flux%lw_up(:,:) = 0._jprb
  flux%lw_dn(:,:) = 0._jprb
  flux%sw_up(:,:) = 0._jprb
  flux%sw_dn(:,:) = 0._jprb
  flux%sw_dn_direct(:,:) = 0._jprb
  flux%lw_up_clear(:,:) = 0._jprb
  flux%lw_dn_clear(:,:) = 0._jprb
  flux%sw_up_clear(:,:) = 0._jprb
  flux%sw_dn_clear(:,:) = 0._jprb
  flux%sw_dn_direct_clear(:,:) = 0._jprb

  flux%lw_dn_surf_canopy(:,:) = 0._jprb
  flux%sw_dn_diffuse_surf_canopy(:,:) = 0._jprb
  flux%sw_dn_direct_surf_canopy(:,:) = 0._jprb
  flux%lw_derivatives(:,:) = 0._jprb

  ! Allocate memory for additional arrays
  allocate(ccn_land(ncol))
  allocate(ccn_sea(ncol))
  allocate(land_frac(ncol))
  allocate(pressure_fl(ncol,nlev))
  allocate(temperature_fl(ncol,nlev))
  allocate(flux_sw_direct_normal(ncol))
  allocate(flux_uv(ncol))
  allocate(flux_par(ncol))
  allocate(flux_par_clear(ncol))
  allocate(emissivity_out(ncol))
  allocate(flux_diffuse_band(ncol,yradiation%yrerad%nsw))
  allocate(flux_direct_band(ncol,yradiation%yrerad%nsw))

  ccn_land = yradiation%yrerad%rccnlnd
  ccn_sea = yradiation%yrerad%rccnsea
  pressure_fl = 0.5_jprb * (thermodynamics%pressure_hl(:,1:nlev)+thermodynamics%pressure_hl(:,2:nlev+1))
  temperature_fl = 0.5_jprb * (thermodynamics%temperature_hl(:,1:nlev)+thermodynamics%temperature_hl(:,2:nlev+1))

  ! --------------------------------------------------------
  ! Section 4a: Reshuffle into blocked memory layout
  ! --------------------------------------------------------

  call ifs_setup_indices(driver_config, ifs_config, yradiation, nlev)

  call ifs_copy_inputs_to_blocked(driver_config, ifs_config, yradiation,&
        & ncol, nlev, single_level, thermodynamics, gas, cloud, aerosol,&
        & sin_latitude, longitude_rad, land_frac, pressure_fl, temperature_fl,&
        & zrgp &
#ifdef BITIDENTITY_TESTING
        &, iseed=iseed &
#endif
        & )

  nproma=driver_config%nblocksize        ! nproma size
  ngpblks=(ncol-1)/nproma+1              ! number of column blocks
  allocate(PMU0(nproma))                        !(KLON, ) ! Cosine of solar zenith ang
  allocate(PTEMPERATURE_SKIN(nproma))           !(KLON) ! (K)
  allocate(PALBEDO_DIF(nproma,yradiation%yrerad%nsw))               !(KLON,YRADIATION%YRERAD%NSW)
  allocate(PALBEDO_DIR(nproma,yradiation%yrerad%nsw))               !(KLON,YRADIATION%YRERAD%NSW)
  allocate(PSPECTRALEMISS(nproma,yradiation%yrerad%nlwemiss))            !(KLON,YRADIATION%YRERAD%NLWEMISS)
  allocate(PGELAM(nproma))                      !(KLON)
  allocate(PGEMU(nproma))                       !(KLON)
  allocate(PLAND_SEA_MASK(nproma))              !((KLON)
  allocate(PPRESSURE(nproma,nlev))          !((KLON,KLEV)    ! (Pa)
  allocate(PTEMPERATURE(nproma,nlev))                    !(KLON,KLEV) ! (K)
  allocate(PPRESSURE_H(nproma,nlev+1))                     !(KLON,KLEV+1)    ! (Pa)
  allocate(PTEMPERATURE_H(nproma,nlev+1))            !(KLON,KLEV+1) ! (K)
  allocate(PQ(nproma,nlev))!(KLON,KLEV)
  allocate(PCO2(nproma,nlev))!(KLON,KLEV)
  allocate(PCH4(nproma,nlev))!(KLON,KLEV)
  allocate(PN2O(nproma,nlev))!(KLON,KLEV)
  allocate(PNO2(nproma,nlev))!(KLON,KLEV)
  allocate(PCFC11(nproma,nlev))!(KLON,KLEV)
  allocate(PCFC12(nproma,nlev))!(KLON,KLEV)
  allocate(PHCFC22(nproma,nlev))!(KLON,KLEV)
  allocate(PCCL4(nproma,nlev))!(KLON,KLEV)
  allocate(PO3(nproma,nlev))!(KLON,KLEV)
  allocate(PCLOUD_FRAC(nproma,nlev))!(KLON,KLEV)
  allocate(PQ_LIQUID(nproma,nlev))!(KLON,KLEV)
  allocate(PQ_ICE(nproma,nlev))!(KLON,KLEV)
  allocate(PQ_RAIN(nproma,nlev))!(KLON,KLEV)
  allocate(PQ_SNOW(nproma,nlev))!(KLON,KLEV)
  allocate(PAEROSOL_OLD(nproma,6,nlev))!(KLON,6,KLEV)
  allocate(PAEROSOL(nproma,nlev,yradiation%rad_config%n_aerosol_types))!(KLON,KLEV,KAEROSOL)
  allocate(PCCN_LAND(nproma))!(KLON)
  allocate(PCCN_SEA(nproma))!(KLON)
#ifdef BITIDENTITY_TESTING
  ! OPTIONAL IN
  allocate(PRE_LIQ(nproma,nlev))!(KLON, KLEV)
  allocate(PRE_ICE(nproma,nlev))!(KLON, KLEV)
  allocate(PCLOUD_OVERLAP(nproma,nlev-1))!(KLON, KLEV-1)
  allocate(iseed_single_block(nproma))!(KLON, KLEV-1)
#endif
  ! OUT
  allocate(PFLUX_SW(nproma,nlev+1))!(KLON,KLEV+1)
  allocate(PFLUX_LW(nproma,nlev+1))!(KLON,KLEV+1)
  allocate(PFLUX_SW_CLEAR(nproma,nlev+1))!(KLON,KLEV+1)
  allocate(PFLUX_LW_CLEAR(nproma,nlev+1))!(KLON,KLEV+1)
  allocate(PFLUX_SW_DN(nproma))!(KLON)
  allocate(PFLUX_LW_DN(nproma))!(KLON)
  allocate(PFLUX_SW_DN_CLEAR(nproma))!(KLON)
  allocate(PFLUX_LW_DN_CLEAR(nproma))!(KLON)
  allocate(PFLUX_DIR(nproma))!(KLON)
  allocate(PFLUX_DIR_CLEAR(nproma))!(KLON)
  allocate(PFLUX_DIR_INTO_SUN(nproma))!(KLON)
  allocate(PFLUX_UV(nproma))!(KLON)
  allocate(PFLUX_PAR(nproma))!(KLON)
  allocate(PFLUX_PAR_CLEAR(nproma))!(KLON)
  allocate(PFLUX_SW_DN_TOA(nproma))!(KLON)
  allocate(PEMIS_OUT(nproma))!(KLON)
  allocate(PLWDERIVATIVE(nproma,nlev+1))!(KLON,KLEV+1)
  allocate(PSWDIFFUSEBAND(nproma,yradiation%yrerad%nsw))!(KLON,YRADIATION%YRERAD%NSW)
  allocate(PSWDIRECTBAND(nproma,yradiation%yrerad%nsw))!(KLON,YRADIATION%YRERAD%NSW)

  ! --------------------------------------------------------
  ! Section 4b: Call radiation_scheme with blocked memory data
  ! --------------------------------------------------------
  if (driver_config%iverbose >= 2) then
    write(nulout,'(a)')  'Performing radiative transfer calculations in new'
  end if

  ! Option of repeating calculation multiple time for more accurate
  ! profiling
  do jrepeat = 1,driver_config%nrepeat
     total_dt=0
     kernel_dt=0
  
#ifdef HAVE_NVTX
    call nvtxStartRange("ecrad_it")
#endif
#ifdef HAVE_ROCTX
    call roctxStartRange("ecrad_it"//c_null_char)
#endif  

#if defined(OMPGPU)
    !$OMP TARGET ENTER DATA MAP(ALLOC: PMU0, PTEMPERATURE_SKIN, &
    !$OMP& PALBEDO_DIF, PALBEDO_DIR, &
    !$OMP& PSPECTRALEMISS, &
    !$OMP& PCCN_LAND, PCCN_SEA, &
    !$OMP& PGELAM, PGEMU, PLAND_SEA_MASK, &
    !$OMP& PPRESSURE, PTEMPERATURE, &
    !$OMP& PPRESSURE_H, PTEMPERATURE_H, &
    !$OMP& PQ, PCO2, PCH4, PN2O, PNO2, &
    !$OMP& PCFC11, PCFC12, PHCFC22, PCCL4, PO3, &
    !$OMP& PCLOUD_FRAC, PQ_LIQUID, PQ_ICE, PQ_RAIN, PQ_SNOW, &
    !$OMP& PAEROSOL_OLD, PAEROSOL, &
    !$OMP& PFLUX_SW, PFLUX_LW, &
    !$OMP& PFLUX_SW_CLEAR, PFLUX_LW_CLEAR, &
    !$OMP& PFLUX_SW_DN, PFLUX_LW_DN, PFLUX_SW_DN_CLEAR, PFLUX_LW_DN_CLEAR, &
    !$OMP& PFLUX_DIR, PFLUX_DIR_CLEAR, PFLUX_DIR_INTO_SUN, &
    !$OMP& PFLUX_UV, PFLUX_PAR, PFLUX_PAR_CLEAR, &
    !$OMP& PFLUX_SW_DN_TOA, PEMIS_OUT, PLWDERIVATIVE, &
    !$OMP& PSWDIFFUSEBAND, PSWDIRECTBAND)
#ifdef BITIDENTITY_TESTING
    !$OMP TARGET ENTER DATA MAP(ALLOC: iseed_single_block,PRE_LIQ,PRE_ICE,PCLOUD_OVERLAP)
#endif ! BITIDENTITY_TESTING
#endif ! OMPGPU


#if defined(_OPENACC)
    !$ACC DATA CREATE(PMU0, PTEMPERATURE_SKIN, &
    !$ACC& PALBEDO_DIF, PALBEDO_DIR, &
    !$ACC& PSPECTRALEMISS, &
    !$ACC& PCCN_LAND, PCCN_SEA, &
    !$ACC& PGELAM, PGEMU, PLAND_SEA_MASK, &
    !$ACC& PPRESSURE, PTEMPERATURE, &
    !$ACC& PPRESSURE_H, PTEMPERATURE_H, &
    !$ACC& PQ, PCO2, PCH4, PN2O, PNO2, &
    !$ACC& PCFC11, PCFC12, PHCFC22, PCCL4, PO3, &
    !$ACC& PCLOUD_FRAC, PQ_LIQUID, PQ_ICE, PQ_RAIN, PQ_SNOW, &
    !$ACC& PAEROSOL_OLD, PAEROSOL, &
    !$ACC& PFLUX_SW, PFLUX_LW, &
    !$ACC& PFLUX_SW_CLEAR, PFLUX_LW_CLEAR, &
    !$ACC& PFLUX_SW_DN, PFLUX_LW_DN, PFLUX_SW_DN_CLEAR, PFLUX_LW_DN_CLEAR, &
    !$ACC& PFLUX_DIR, PFLUX_DIR_CLEAR, PFLUX_DIR_INTO_SUN, &
    !$ACC& PFLUX_UV, PFLUX_PAR, PFLUX_PAR_CLEAR, &
    !$ACC& PFLUX_SW_DN_TOA, PEMIS_OUT, PLWDERIVATIVE, &
    !$ACC& PSWDIFFUSEBAND, PSWDIRECTBAND)
#ifdef BITIDENTITY_TESTING
    !$ACC DATA CREATE(iseed_single_block,PRE_LIQ,PRE_ICE,PCLOUD_OVERLAP)
#endif ! BITIDENTITY_TESTING
#endif ! _OPENACC

!    if (driver_config%do_parallel) then
      ! Run radiation scheme over blocks of columns in parallel

      total_dt=0
      kernel_dt=0
      do jrl=1,ncol,nproma
        ibeg=jrl
        iend=min(ibeg+nproma-1,ncol)
        il=iend-ibeg+1
        ib=(jrl-1)/nproma+1

#ifdef HAVE_NVTX
        call nvtxStartRange("ifs_copy_inputs_to_blocked_new")
#endif
#ifdef HAVE_ROCTX
        call roctxStartRange("ifs_copy_inputs_to_blocked_new"//c_null_char)
#endif  

        ! Copy to the small array
        call ifs_copy_inputs_to_blocked_new(driver_config, ifs_config, yradiation,&
             &  single_level, ncol, nlev, jrl, zrgp, &
             &  PMU0, PTEMPERATURE_SKIN, PALBEDO_DIF, PALBEDO_DIR, &
             &  PSPECTRALEMISS, &
             &  PCCN_LAND, PCCN_SEA, &
             &  PGELAM, PGEMU, PLAND_SEA_MASK, &
             &  PPRESSURE, PTEMPERATURE, &
             &  PPRESSURE_H, PTEMPERATURE_H, &
             &  PQ, PCO2, PCH4, PN2O, PNO2, PCFC11, PCFC12, PHCFC22, PCCL4, PO3, &
             &  PCLOUD_FRAC, PQ_LIQUID, PQ_ICE, PQ_RAIN, PQ_SNOW, &
             &  PAEROSOL_OLD, PAEROSOL, &
             &  PFLUX_SW, PFLUX_LW, PFLUX_SW_CLEAR, PFLUX_LW_CLEAR, &
             &  PFLUX_SW_DN, PFLUX_LW_DN, PFLUX_SW_DN_CLEAR, PFLUX_LW_DN_CLEAR, &
             &  PFLUX_DIR, PFLUX_DIR_CLEAR, PFLUX_DIR_INTO_SUN, &
             &  PFLUX_UV, PFLUX_PAR, PFLUX_PAR_CLEAR, &
             &  PFLUX_SW_DN_TOA, PEMIS_OUT, PLWDERIVATIVE, &
             &  PSWDIFFUSEBAND, PSWDIRECTBAND, &
             ! OPTIONAL ARGUMENTS for bit-identical results in tests
#ifdef BITIDENTITY_TESTING
             &  iseed=iseed_single_block, PRE_LIQ=PRE_LIQ, PRE_ICE=PRE_ICE, PCLOUD_OVERLAP=PCLOUD_OVERLAP)
#else
             & )
#endif
             
#ifdef HAVE_NVTX
        call nvtxEndRange
        call nvtxStartRange("copy_to_device")
#endif
#ifdef HAVE_ROCTX
        call roctxEndRange
        call roctxStartRange("copy_to_device"//c_null_char)
#endif  

        call system_clock (count=t(1))

#if defined(OMPGPU)
        !$OMP TARGET UPDATE TO(PMU0, PTEMPERATURE_SKIN, &
        !$OMP& PALBEDO_DIF, PALBEDO_DIR, &
        !$OMP& PSPECTRALEMISS, &
        !$OMP& PCCN_LAND, PCCN_SEA, &
        !$OMP& PGELAM, PGEMU, PLAND_SEA_MASK, &
        !$OMP& PPRESSURE, PTEMPERATURE, &
        !$OMP& PPRESSURE_H, PTEMPERATURE_H, &
        !$OMP& PQ, PCO2, PCH4, PN2O, PNO2, &
        !$OMP& PCFC11, PCFC12, PHCFC22, PCCL4, PO3, &
        !$OMP& PCLOUD_FRAC, PQ_LIQUID, PQ_ICE, PQ_RAIN, PQ_SNOW, &
        !$OMP& PAEROSOL_OLD, PAEROSOL)
#ifdef BITIDENTITY_TESTING
        !$OMP TARGET UPDATE TO(iseed_single_block,PRE_LIQ,PRE_ICE,PCLOUD_OVERLAP)
#endif ! BITIDENTITY_TESTING
#endif ! OMPGPU

#if defined(_OPENACC)
        !$ACC UPDATE DEVICE(PMU0, PTEMPERATURE_SKIN, &
        !$ACC& PALBEDO_DIF, PALBEDO_DIR, &
        !$ACC& PSPECTRALEMISS, &
        !$ACC& PCCN_LAND, PCCN_SEA, &
        !$ACC& PGELAM, PGEMU, PLAND_SEA_MASK, &
        !$ACC& PPRESSURE, PTEMPERATURE, &
        !$ACC& PPRESSURE_H, PTEMPERATURE_H, &
        !$ACC& PQ, PCO2, PCH4, PN2O, PNO2, &
        !$ACC& PCFC11, PCFC12, PHCFC22, PCCL4, PO3, &
        !$ACC& PCLOUD_FRAC, PQ_LIQUID, PQ_ICE, PQ_RAIN, PQ_SNOW, &
        !$ACC& PAEROSOL_OLD, PAEROSOL)
#ifdef BITIDENTITY_TESTING
        !$ACC UPDATE DEVICE(iseed_single_block,PRE_LIQ,PRE_ICE,PCLOUD_OVERLAP)
#endif ! BITIDENTITY_TESTING
#endif ! _OPENACC

#ifdef HAVE_NVTX
        call nvtxEndRange
        call nvtxStartRange("radiation_scheme")
#endif
#ifdef HAVE_ROCTX
        call roctxEndRange
        call roctxStartRange("radiation_scheme"//c_null_char)
#endif  
        call system_clock (count=t(3))

        ! Call the ECRAD radiation scheme
        call radiation_scheme &
             & (yradiation, &
             &  1, il, nproma, &                       ! startcol, endcol, ncol
             &  nlev, size(aerosol%mixing_ratio,3), &    ! nlev, naerosols
             &  single_level%solar_irradiance, &                               ! solar_irrad
             ! array inputs
             & PMU0, PTEMPERATURE_SKIN, PALBEDO_DIF, PALBEDO_DIR, &
             &  PSPECTRALEMISS, &
             &  PCCN_LAND, PCCN_SEA, &
             &  PGELAM, PGEMU, PLAND_SEA_MASK, &
             &  PPRESSURE, PTEMPERATURE, &
             &  PPRESSURE_H, PTEMPERATURE_H, &
             &  PQ, PCO2, PCH4, PN2O, PNO2, &
             &  PCFC11, PCFC12, PHCFC22, PCCL4, PO3, &
             &  PCLOUD_FRAC, PQ_LIQUID, PQ_ICE, PQ_RAIN, PQ_SNOW, &
             &  PAEROSOL_OLD, PAEROSOL, &
             &  PFLUX_SW, PFLUX_LW, PFLUX_SW_CLEAR, PFLUX_LW_CLEAR, &
             &  PFLUX_SW_DN, PFLUX_LW_DN, PFLUX_SW_DN_CLEAR, PFLUX_LW_DN_CLEAR, &
             &  PFLUX_DIR, PFLUX_DIR_CLEAR, PFLUX_DIR_INTO_SUN, &
             &  PFLUX_UV, PFLUX_PAR, PFLUX_PAR_CLEAR, &
             &  PFLUX_SW_DN_TOA, PEMIS_OUT, PLWDERIVATIVE, &
             &  PSWDIFFUSEBAND, PSWDIRECTBAND, &
             ! OPTIONAL ARGUMENTS for bit-identical results in tests
#ifdef BITIDENTITY_TESTING
            ! To validate results against standalone ecrad, we overwrite effective
            ! radii, cloud overlap and seed with input values
             &  PRE_LIQ=PRE_LIQ, &
             &  PRE_ICE=PRE_ICE, &
             &  PCLOUD_OVERLAP=PCLOUD_OVERLAP, &
             &  iseed=iseed_single_block &
#endif
             & )

#ifdef HAVE_NVTX
        call nvtxEndRange
        call nvtxStartRange("copy_from_device")
#endif
#ifdef HAVE_ROCTX
        call roctxEndRange
        call roctxStartRange("copy_from_device"//c_null_char)
#endif  
        call system_clock (count=t(4))

#if defined(OMPGPU)
        !$OMP TARGET UPDATE FROM(PFLUX_SW, PFLUX_LW, &
        !$OMP&  PFLUX_SW_CLEAR, PFLUX_LW_CLEAR, &
        !$OMP&  PFLUX_SW_DN, PFLUX_LW_DN, PFLUX_SW_DN_CLEAR, PFLUX_LW_DN_CLEAR, &
        !$OMP&  PFLUX_DIR, PFLUX_DIR_CLEAR, PFLUX_DIR_INTO_SUN, &
        !$OMP&  PFLUX_UV, PFLUX_PAR, PFLUX_PAR_CLEAR, &
        !$OMP&  PFLUX_SW_DN_TOA, PEMIS_OUT, PLWDERIVATIVE, &
        !$OMP&  PSWDIFFUSEBAND, PSWDIRECTBAND)
#endif

#if defined(_OPENACC)
        !$ACC UPDATE HOST(PFLUX_SW, PFLUX_LW, &
        !$ACC&  PFLUX_SW_CLEAR, PFLUX_LW_CLEAR, &
        !$ACC&  PFLUX_SW_DN, PFLUX_LW_DN, PFLUX_SW_DN_CLEAR, PFLUX_LW_DN_CLEAR, &
        !$ACC&  PFLUX_DIR, PFLUX_DIR_CLEAR, PFLUX_DIR_INTO_SUN, &
        !$ACC&  PFLUX_UV, PFLUX_PAR, PFLUX_PAR_CLEAR, &
        !$ACC&  PFLUX_SW_DN_TOA, PEMIS_OUT, PLWDERIVATIVE, &
        !$ACC&  PSWDIFFUSEBAND, PSWDIRECTBAND)
#endif

        call system_clock (count=t(2))
        total_dt = total_dt + (t(2)-t(1))/dble(count_rate)
        kernel_dt = kernel_dt + (t(4)-t(3))/dble(count_rate)

#ifdef HAVE_NVTX
        call nvtxEndRange
        call nvtxStartRange("ifs_copy_inputs_from_blocked_new")
#endif
#ifdef HAVE_ROCTX
        call roctxEndRange
        call roctxStartRange("ifs_copy_inputs_from_blocked_new"//c_null_char)
#endif  

        call ifs_copy_fluxes_from_blocked_new(driver_config, ifs_config, yradiation, ncol, nlev, jrl, &
             &  PFLUX_SW, PFLUX_LW, PFLUX_SW_CLEAR, PFLUX_LW_CLEAR, &
             &  PFLUX_SW_DN, PFLUX_LW_DN, PFLUX_SW_DN_CLEAR, PFLUX_LW_DN_CLEAR, &
             &  PFLUX_DIR, PFLUX_DIR_CLEAR, PFLUX_DIR_INTO_SUN, &
             &  PFLUX_UV, PFLUX_PAR, PFLUX_PAR_CLEAR, &
             &  PFLUX_SW_DN_TOA, PEMIS_OUT, PLWDERIVATIVE, &
             &  PSWDIFFUSEBAND, PSWDIRECTBAND, zrgp)

#ifdef HAVE_NVTX
        call nvtxEndRange
#endif
#ifdef HAVE_ROCTX
        call roctxEndRange
#endif  

     end do

     write(nulout, '(a,g12.5,a)') 'time elapsed in radiative transfer: ', &
          &                         (total_dt), ' seconds'
     write(nulout, '(a,g12.5,a)') 'time elapsed in radiative transfer kernel: ', &
          &                         (kernel_dt), ' seconds'
     write(nulout, '(a,i0)') 'Columns/s : ', int((ncol)/(total_dt))
     write(nulout, '(a,i0)') 'Columns/s : ', int((ncol)/(kernel_dt))

#if defined(OMPGPU)
     !$OMP TARGET EXIT DATA MAP(DELETE: PMU0, PTEMPERATURE_SKIN, &
     !$OMP& PALBEDO_DIF, PALBEDO_DIR, &
     !$OMP& PSPECTRALEMISS, &
     !$OMP& PCCN_LAND, PCCN_SEA, &
     !$OMP& PGELAM, PGEMU, PLAND_SEA_MASK, &
     !$OMP& PPRESSURE, PTEMPERATURE, &
     !$OMP& PPRESSURE_H, PTEMPERATURE_H, &
     !$OMP& PQ, PCO2, PCH4, PN2O, PNO2, &
     !$OMP& PCFC11, PCFC12, PHCFC22, PCCL4, PO3, &
     !$OMP& PCLOUD_FRAC, PQ_LIQUID, PQ_ICE, PQ_RAIN, PQ_SNOW, &
     !$OMP& PAEROSOL_OLD, PAEROSOL, &
     !$OMP& PFLUX_SW, PFLUX_LW, &
     !$OMP& PFLUX_SW_CLEAR, PFLUX_LW_CLEAR, &
     !$OMP& PFLUX_SW_DN, PFLUX_LW_DN, PFLUX_SW_DN_CLEAR, PFLUX_LW_DN_CLEAR, &
     !$OMP& PFLUX_DIR, PFLUX_DIR_CLEAR, PFLUX_DIR_INTO_SUN, &
     !$OMP& PFLUX_UV, PFLUX_PAR, PFLUX_PAR_CLEAR, &
     !$OMP& PFLUX_SW_DN_TOA, PEMIS_OUT, PLWDERIVATIVE, &
     !$OMP& PSWDIFFUSEBAND, PSWDIRECTBAND) NOWAIT
#ifdef BITIDENTITY_TESTING
     !$OMP TARGET EXIT DATA MAP(DELETE: iseed_single_block, PRE_LIQ, PRE_ICE, PCLOUD_OVERLAP)
#endif ! BITIDENTITY_TESTING
#endif ! OMPGPU

#if defined(_OPENACC)
     !$ACC END DATA
#ifdef BITIDENTITY_TESTING
     !$ACC END DATA
#endif ! BITIDENTITY_TESTING
#endif ! _OPENACC

#ifdef HAVE_NVTX
     call nvtxEndRange
#endif
#ifdef HAVE_ROCTX
     call roctxEndRange
#endif  

  end do

  ! --------------------------------------------------------
  ! Section 4c: Copy fluxes from blocked memory data
  ! --------------------------------------------------------
  call ifs_copy_fluxes_from_blocked(driver_config, ifs_config, yradiation, ncol, nlev,&
          & zrgp, flux, flux_sw_direct_normal, flux_uv, flux_par, flux_par_clear, &
          & emissivity_out, flux_diffuse_band, flux_direct_band)

  ! "up" fluxes are actually net fluxes at this point - we modify the
  ! upwelling flux so that net=dn-up, while the TOA and surface
  ! downwelling fluxes are correct.
  flux%sw_up = -flux%sw_up
  flux%sw_up(:,1) = flux%sw_up(:,1)+flux%sw_dn(:,1)
  flux%sw_up(:,nlev+1) = flux%sw_up(:,nlev+1)+flux%sw_dn(:,nlev+1)

  flux%lw_up = -flux%lw_up
  flux%lw_up(:,1) = flux%lw_up(:,1)+flux%lw_dn(:,1)
  flux%lw_up(:,nlev+1) = flux%lw_up(:,nlev+1)+flux%lw_dn(:,nlev+1)

  flux%sw_up_clear = -flux%sw_up_clear
  flux%sw_up_clear(:,1) = flux%sw_up_clear(:,1)+flux%sw_dn_clear(:,1)
  flux%sw_up_clear(:,nlev+1) = flux%sw_up_clear(:,nlev+1)+flux%sw_dn_clear(:,nlev+1)

  flux%lw_up_clear = -flux%lw_up_clear
  flux%lw_up_clear(:,1) = flux%lw_up_clear(:,1)+flux%lw_dn_clear(:,1)
  flux%lw_up_clear(:,nlev+1) = flux%lw_up_clear(:,nlev+1)+flux%lw_dn_clear(:,nlev+1)

  ! Deallocate
  if (allocated(PMU0)) deallocate(PMU0)
  if (allocated(PTEMPERATURE_SKIN)) deallocate(PTEMPERATURE_SKIN)
  if (allocated(PALBEDO_DIF)) deallocate(PALBEDO_DIF)
  if (allocated(PALBEDO_DIR)) deallocate(PALBEDO_DIR)
  if (allocated(PSPECTRALEMISS)) deallocate(PSPECTRALEMISS)
  if (allocated(PGELAM)) deallocate(PGELAM)
  if (allocated(PGEMU)) deallocate(PGEMU)
  if (allocated(PLAND_SEA_MASK)) deallocate(PLAND_SEA_MASK)
  if (allocated(PPRESSURE)) deallocate(PPRESSURE)
  if (allocated(PTEMPERATURE)) deallocate(PTEMPERATURE)
  if (allocated(PPRESSURE_H)) deallocate(PPRESSURE_H)
  if (allocated(PTEMPERATURE_H)) deallocate(PTEMPERATURE_H)
  if (allocated(PQ)) deallocate(PQ)
  if (allocated(PCO2)) deallocate(PCO2)
  if (allocated(PCH4)) deallocate(PCH4)
  if (allocated(PN2O)) deallocate(PN2O)
  if (allocated(PNO2)) deallocate(PNO2)
  if (allocated(PCFC11)) deallocate(PCFC11)
  if (allocated(PCFC12)) deallocate(PCFC12)
  if (allocated(PHCFC22)) deallocate(PHCFC22)
  if (allocated(PCCL4)) deallocate(PCCL4)
  if (allocated(PO3)) deallocate(PO3)
  if (allocated(PCLOUD_FRAC)) deallocate(PCLOUD_FRAC)
  if (allocated(PQ_LIQUID)) deallocate(PQ_LIQUID)
  if (allocated(PQ_ICE)) deallocate(PQ_ICE)
  if (allocated(PQ_RAIN)) deallocate(PQ_RAIN)
  if (allocated(PQ_SNOW)) deallocate(PQ_SNOW)
  if (allocated(PAEROSOL_OLD)) deallocate(PAEROSOL_OLD)
  if (allocated(PAEROSOL)) deallocate(PAEROSOL)
  if (allocated(PCCN_LAND)) deallocate(PCCN_LAND)
  if (allocated(PCCN_SEA)) deallocate(PCCN_SEA)
  ! OPTIONAL IN
#ifdef BITIDENTITY_TESTING
  if (allocated(PRE_LIQ)) deallocate(PRE_LIQ)
  if (allocated(PRE_ICE)) deallocate(PRE_ICE)
  if (allocated(PCLOUD_OVERLAP)) deallocate(PCLOUD_OVERLAP)
  if (allocated(iseed_single_block)) deallocate(iseed_single_block)
#endif
  ! OUT
  if (allocated(PFLUX_SW)) deallocate(PFLUX_SW)
  if (allocated(PFLUX_LW)) deallocate(PFLUX_LW)
  if (allocated(PFLUX_SW_CLEAR)) deallocate(PFLUX_SW_CLEAR)
  if (allocated(PFLUX_LW_CLEAR)) deallocate(PFLUX_LW_CLEAR)
  if (allocated(PFLUX_SW_DN)) deallocate(PFLUX_SW_DN)
  if (allocated(PFLUX_LW_DN)) deallocate(PFLUX_LW_DN)
  if (allocated(PFLUX_SW_DN_CLEAR)) deallocate(PFLUX_SW_DN_CLEAR)
  if (allocated(PFLUX_LW_DN_CLEAR)) deallocate(PFLUX_LW_DN_CLEAR)
  if (allocated(PFLUX_DIR)) deallocate(PFLUX_DIR)
  if (allocated(PFLUX_DIR_CLEAR)) deallocate(PFLUX_DIR_CLEAR)
  if (allocated(PFLUX_DIR_INTO_SUN)) deallocate(PFLUX_DIR_INTO_SUN)
  if (allocated(PFLUX_UV)) deallocate(PFLUX_UV)
  if (allocated(PFLUX_PAR)) deallocate(PFLUX_PAR)
  if (allocated(PFLUX_PAR_CLEAR)) deallocate(PFLUX_PAR_CLEAR)
  if (allocated(PFLUX_SW_DN_TOA)) deallocate(PFLUX_SW_DN_TOA)
  if (allocated(PEMIS_OUT)) deallocate(PEMIS_OUT)
  if (allocated(PLWDERIVATIVE)) deallocate(PLWDERIVATIVE)
  if (allocated(PSWDIFFUSEBAND)) deallocate(PSWDIFFUSEBAND)
  if (allocated(PSWDIRECTBAND)) deallocate(PSWDIRECTBAND)
  if (allocated(zrgp)) deallocate(zrgp)
  if (allocated(iseed)) deallocate(iseed)
  
  ! --------------------------------------------------------
  ! Section 5: Check and save output
  ! --------------------------------------------------------

  ! This is unreliable because only the net fluxes are valid:
  !is_out_of_bounds = flux%out_of_physical_bounds(driver_config%istartcol, driver_config%iendcol)

  ! Store the fluxes in the output file
  yradiation%rad_config%do_surface_sw_spectral_flux = .false.
  yradiation%rad_config%do_canopy_fluxes_sw = .false.
  yradiation%rad_config%do_canopy_fluxes_lw = .false.

  call save_net_fluxes(file_name, yradiation%rad_config, thermodynamics, flux, &
       &   iverbose=driver_config%iverbose, is_hdf5_file=driver_config%do_write_hdf5, &
       &   experiment_name=driver_config%experiment_name, &
       &   is_double_precision=driver_config%do_write_double_precision)

  if (driver_config%iverbose >= 2) then
    write(nulout,'(a)') '------------------------------------------------------------------------------------'
  end if

  ! Finalise MPI if not done yet
#ifdef HAVE_FIAT
  call mpl_end(ldmeminfo=.false.)
#endif

end program ecrad_ifs_driver
