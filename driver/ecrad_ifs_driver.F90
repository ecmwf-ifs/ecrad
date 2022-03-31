! ecrad_ifs_driver.F90 - IFS driver for offline ECRAD radiation scheme
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
  use parkind1,                 only : jprb, jprd, jpim ! Working/double precision

  use radiation_io,             only : nulout
  use radiation_interface,      only : setup_radiation, radiation, set_gas_units
  ! use radiation_config,         only : config_type
  use radiation_single_level,   only : single_level_type
!  use radsurf_properties,       only : surface_type, print_surface_representation
!  use radsurf_intermediate,     only : surface_intermediate_type
!  use radsurf_flux,             only : surface_flux_type
!  use radsurf_save,             only : save_surface_fluxes
  use radiation_thermodynamics, only : thermodynamics_type
  use radiation_gas,            only : gas_type, &
       &   IVolumeMixingRatio, IMassMixingRatio, &
       &   IH2O, ICO2, IO3, IN2O, ICO, ICH4, IO2, ICFC11, ICFC12, &
       &   IHCFC22, ICCl4, GasName, GasLowerCaseName, NMaxGases
  use radiation_cloud,          only : cloud_type
  use radiation_aerosol,        only : aerosol_type
  use radiation_flux,           only : flux_type
  use radiation_save,           only : save_fluxes, save_inputs
  use radiation_setup,          only : tcompo, setup_radiation_scheme
  use ecrad_driver_config,      only : driver_config_type
  use ecrad_driver_read_input,  only : read_input
  use easy_netcdf
  use type_model,               only : model

  implicit none

  ! The NetCDF file containing the input profiles
  type(netcdf_file)         :: file

  ! Derived type for the surface inputs
!  type(surface_type)        :: surface

  ! Derived types for the inputs to the radiation scheme
  ! type(config_type)         :: config
  type(single_level_type)   :: single_level
  type(thermodynamics_type) :: thermodynamics
  type(gas_type)            :: gas
  type(cloud_type)          :: cloud
  type(aerosol_type)        :: aerosol

  ! Derived types for the surface fluxes
!  type(surface_intermediate_type) :: surface_intermediate
!  type(surface_flux_type)         :: surface_flux

  ! Configuration specific to this driver
  type(driver_config_type)  :: driver_config

  ! Derived type containing outputs from the radiation scheme
  type(flux_type)           :: flux

  ! Model type
  type(model)               :: ydmodel

  ! Dummy COMPO type
  type(tcompo)              :: ydcompo

  integer :: ncol, nlev         ! Number of columns and levels
  integer :: istartcol, iendcol ! Range of columns to process

  ! Name of file names specified on command line
  character(len=512) :: file_name
  integer            :: istatus ! Result of command_argument_count

  ! For parallel processing of multiple blocks
  integer :: jblock, nblock ! Block loop index and number
  integer, external :: omp_get_thread_num
  double precision, external :: omp_get_wtime

  ! Loop index for repeats (for benchmarking)
  integer :: jrepeat

  ! Are any variables out of bounds?
  logical :: is_out_of_bounds

  ! Are we using a complex surface representation stored in the
  ! "surface" structure?
!  logical :: is_complex_surface

!  integer    :: iband(20), nweights
!  real(jprb) :: weight(20)

  ! Start/stop time in seconds
  real(kind=jprd) :: tstart, tstop

  ! monolithic IFS data structure to pass to radiation scheme
  real(kind=jprb), allocatable :: zrgp(:,:,:)

  ! latitude, longitute
  real(kind=jprb), allocatable :: lat(:), lon(:)

  ! solar irradiance
  real(kind=jprb) :: zrii0
  ! number of column blocks
  integer :: ngpblks

  integer :: ifldsin, ifldsout, ifldstot

  integer :: nproma
  
  integer :: inext, iinbeg, iinend, ioutbeg, ioutend, igi, imu0, iamu0, iemiss, its, islm, iccnl,    &
       &     ibas, itop, igelam, igemu, iclon, islon, iald, ialp, iti, ipr, iqs, iwv, iclc, ilwa,    & 
       &     iiwa, iswa, irwa, irra, idp, ioz, iecpo3, ihpr, iaprs, ihti, iaero, ifrsod, icdir,      &
       &     ifrted, ifrsodc, ifrtedc, iemit, isudu, iuvdf, iparf, iparcf, itincf, ifdir, ifdif,     &
       &     ilwderivative, iswdirectband, iswdiffuseband, ifrso, iswfc, ifrth, ilwfc, iaer,         &
       &     iich4, iin2o, ino2, ic11, ic12, igix, iico2, iccno, ic22, icl4  

  integer :: ire_liq, ire_ice

  integer :: nlwout, nsw, nlwemiss, iradaer, naermacc

  logical :: lrayfm, lepo3ra, lchem_o3rad, ldiagforcing, lapproxlwupdate, lapproxswupdate, lldebug

  integer :: jrl, ibeg, iend, ib, il, ifld, jlev, jaer, joff, jalb, jemiss


  ! --------------------------------------------------------
  ! Section 2: Configure
  ! --------------------------------------------------------

  ! Check program called with correct number of arguments
  if (command_argument_count() < 3) then
    stop 'Usage: ecrad config.nam input_file.nc output_file.nc [output_surface_file.nc]'
  end if

  ! Use namelist to configure the radiation calculation
  call get_command_argument(1, file_name, status=istatus)
  if (istatus /= 0) then
    stop 'Failed to read name of namelist file as string of length < 512'
  end if

  ! Read "radiation" namelist into radiation configuration type
  ! call config%read(file_name=file_name)

  associate( &
    & yderdi=>ydmodel%yrml_phy_rad%yrerdi, &
    & ydeaeratm=>ydmodel%yrml_phy_rad%yreaeratm, &
    & ydephy=>ydmodel%yrml_phy_ec%yrephy, &
    & yderad=>ydmodel%yrml_phy_rad%yrerad, &
    & ydradiation=>ydmodel%yrml_phy_rad%yradiation &
  )
    ! Values from IFS:
    ! TODO: fill this from ecrad inputs

    ydephy%nalbedoscheme = 2
    ! ydeaeratm values???
    yderad%naermacc = 1  ! MACC aerosol climatology
    yderad%lapproxlwupdate = .true.  ! Hogan and Bozzo (2015) approx longwave updates
    yderad%lapproxswupdate = .true.  ! Hogan and Bozzo (2015) approx shortwave updates
    yderad%nliqopt = 4  ! 2 - SLINGO, 4 - SOCRATES
    yderad%niceopt = 3  ! 3 - ICEMODELFU, 4 - IDEMODELBARAN
    yderad%lfu_lw_ice_optics_bug = .false.
    yderad%nlwsolver = 0  ! 0 - McICA, 1 - SPARTACUS, 2 - SPARTACUS 3D, 3 - Tripleclouds
    yderad%nswsolver = 0  ! 0 - McICA, 1 - SPARTACUS, 2 - SPARTACUS 3D, 3 - Tripleclouds
    yderad%nlwscattering = 1  ! 1 - cloud scattering, 2 - cloud + aerosol scattering
    yderad%ncloudoverlap = 3  ! 1 - MAXIMUMRANDOM, 2 - EXPONENTIAL, 3 - EXPONENTIALRANDOM
    yderad%nsolarspectrum = 0
    yderad%ndecolat = 2
    yderad%nlwout = 1
    yderad%nlwemiss = 2
    yderad%nsw = 6

    ydradiation%rad_config%do_setup_ifsrrtm = .true.  ! Make sure ifsrrtm gets initialized

    call setup_radiation_scheme(yderdi, ydeaeratm, ydcompo, ydephy, yderad, &
                                & ydradiation, ldoutput=.true., file_name=file_name)
  end associate

  ! Read "radiation_driver" namelist into radiation driver config type
  call driver_config%read(file_name)

  if (driver_config%iverbose >= 2) then
    write(nulout,'(a)') '-------------------------- OFFLINE ECRAD RADIATION SCHEME --------------------------'
    write(nulout,'(a)') 'Copyright (C) 2014- ECMWF'
    write(nulout,'(a)') 'Contact: Robin Hogan (r.j.hogan@ecmwf.int)'
#ifdef SINGLE_PRECISION
    write(nulout,'(a)') 'Floating-point precision: single'
#else
    write(nulout,'(a)') 'Floating-point precision: double'
#endif
    call ydmodel%yrml_phy_rad%yradiation%rad_config%print(driver_config%iverbose)
  end if

  ! Albedo/emissivity intervals may be specified like this
  !call config%define_sw_albedo_intervals(6, &
  !     &  [0.25e-6_jprb, 0.44e-6_jprb, 0.69e-6_jprb, &
  !     &     1.19_jprb, 2.38e-6_jprb], [1,2,3,4,5,6], &
  !     &   do_nearest=.false.)
  !call config%define_lw_emiss_intervals(3, &
  !     &  [8.0e-6_jprb, 13.0e-6_jprb], [1,2,1], &
  !     &   do_nearest=.false.)

  ! Setup the radiation scheme: load the coefficients for gas and
  ! cloud optics, currently from RRTMG
  !call setup_radiation(config)

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

  ! Read input variables from NetCDF file
  call read_input(file, ydmodel%yrml_phy_rad%yradiation%rad_config, driver_config, ncol, nlev, &
!       &          is_complex_surface, surface, &
       &          single_level, thermodynamics, &
       &          gas, cloud, aerosol, &
       &          lat=lat, lon=lon)

!  if (is_complex_surface) then
!    config%do_canopy_fluxes_sw = .true.
!    config%do_canopy_fluxes_lw = .true.
!  end if

  ! Close input file
  call file%close()

!  if (is_complex_surface .and. driver_config%iverbose >= 2) then
!    call print_surface_representation(surface%i_representation)
!  end if

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
  
  ! Store inputs
  ! if (driver_config%do_save_inputs) then
  !   call save_inputs('inputs.nc', ydmodel%yrml_phy_rad%yradiation%rad_config, single_level, thermodynamics, &
  !        &                gas, cloud, aerosol, &
  !        &                lat=spread(0.0_jprb,1,ncol), &
  !        &                lon=spread(0.0_jprb,1,ncol), &
  !        &                iverbose=driver_config%iverbose)
  ! end if

  ! --------------------------------------------------------
  ! Section 4: Call radiation scheme
  ! --------------------------------------------------------

  ! Ensure the units of the gas mixing ratios are what is required
  ! by the gas absorption model
  !call set_gas_units(ydmodel%yrml_phy_rad%yradiation%rad_config, gas)

  ! Compute saturation with respect to liquid (needed for aerosol
  ! hydration) call
  !call thermodynamics%calc_saturation_wrt_liquid(driver_config%istartcol,driver_config%iendcol)

!  if (is_complex_surface) then
!    call surface_intermediate%allocate(driver_config%istartcol, driver_config%iendcol, &
!         &                             config, surface)
!    call surface_flux%allocate(config, driver_config%istartcol, driver_config%iendcol, &
!         &                     surface%i_representation)
!  end if

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
  !!!call flux%allocate(config, 1, ncol, nlev)

  !
  !
  !    IFS SETUP
  !
  !

  zrii0=single_level%solar_irradiance
  lldebug=(driver_config%iverbose>4)     ! debug 
  nproma=driver_config%nblocksize        ! nproma size
  ngpblks=(ncol-1)/nproma+1              ! number of column blocks

  ! number of emissivity intervals, = 2 in IFS
  if ( ydmodel%yrml_phy_rad%yradiation%rad_config%use_canopy_full_spectrum_lw ) then
    nlwemiss = ydmodel%yrml_phy_rad%yradiation%rad_config%n_g_lw
  else
    nlwemiss=2!surface%nemissbands ???
  endif

  ! number of shortwave spectral intervals, = 6 in IFS
  if (ydmodel%yrml_phy_rad%yradiation%rad_config%use_canopy_full_spectrum_sw) then
    nsw = ydmodel%yrml_phy_rad%yradiation%rad_config%n_g_sw  
  else
    nsw = 6!surface%nalbedobands    ???
  endif

  nlwout=1                           ! number of spectral intervals to pass LW downwelling flux
  lrayfm=.false.                     ! key for calling radiation scheme from MF 
  lepo3ra=.false.                    ! true if prognostic o-zone is passed to radiation
  lchem_o3rad=.false.                ! use chemistry o3 in radiation
  iradaer=ydmodel%yrml_phy_rad%yradiation%rad_config%n_aerosol_types     ! number of aerosol species
  naermacc=1                         ! MACC-derived aerosol climatology on a NMCLAT x NMCLON grid
  ldiagforcing=.false.               ! write input ozone, ghg and aerosol forcing to 3D fields
  lapproxlwupdate=.true.
  lapproxswupdate=.true.

  !
  ! RADINTG
  !
 
  !  INITIALISE INDICES FOR VARIABLE

  ! INDRAD is a CONTAIN'd function (at end of this routine)
 
  inext  =1
  iinbeg =1                        ! start of input variables
  igi    =indrad(inext,1,lldebug)
  imu0   =indrad(inext,1,.true.)
  iamu0  =indrad(inext,1,.true.)
  iemiss =indrad(inext,nlwemiss,.true.)
  its    =indrad(inext,1,.true.)
  islm   =indrad(inext,1,.true.)
  iccnl  =indrad(inext,1,.true.)
  iccno  =indrad(inext,1,.true.)
  ibas   =indrad(inext,1,.true.)
  itop   =indrad(inext,1,.true.)
  igelam =indrad(inext,1,.true.)
  igemu  =indrad(inext,1,.true.)
  iclon  =indrad(inext,1,.true.)
  islon  =indrad(inext,1,.true.)
  iald   =indrad(inext,nsw,.true.)
  ialp   =indrad(inext,nsw,.true.)
  iti    =indrad(inext,nlev,.true.)
  ipr    =indrad(inext,nlev,.true.)
  iqs    =indrad(inext,nlev,.true.)
  iwv    =indrad(inext,nlev,.true.)
  iclc   =indrad(inext,nlev,.true.)
  ilwa   =indrad(inext,nlev,.true.)
  iiwa   =indrad(inext,nlev,.true.)
  iswa   =indrad(inext,nlev,.true.)
  irwa   =indrad(inext,nlev,.true.)
  irra   =indrad(inext,nlev,.true.)
  idp    =indrad(inext,nlev,.true.)
  ioz    =indrad(inext,nlev,lrayfm)
  iecpo3 =indrad(inext,nlev ,.false.)
  ihpr   =indrad(inext,nlev+1,.true.) ! not used in ecrad
  iaprs  =indrad(inext,nlev+1,.true.)
  ihti   =indrad(inext,nlev+1,.true.)
  iaero  =indrad(inext,iradaer*nlev,.false.)

  iinend =inext-1                  ! end of input variables

  ioutbeg=inext                    ! start of output variables
  if (naermacc == 1) iaero = indrad(inext,iradaer*nlev,ldiagforcing)
  ifrsod =indrad(inext,1,.true.)
  ifrted =indrad(inext,nlwout,.true.)
  ifrsodc=indrad(inext,1,.true.)
  ifrtedc=indrad(inext,1,.true.)
  iemit  =indrad(inext,1,.true.)
  isudu  =indrad(inext,1,.true.)
  iuvdf  =indrad(inext,1,.true.)
  iparf  =indrad(inext,1,.true.)
  iparcf =indrad(inext,1,.true.)
  itincf =indrad(inext,1,.true.)
  ifdir  =indrad(inext,1,.true.)
  ifdif  =indrad(inext,1,.true.)
  icdir  =indrad(inext,1,.true.)
  ilwderivative =indrad(inext,nlev+1, lapproxlwupdate)
  iswdirectband =indrad(inext,nsw,lapproxswupdate)
  iswdiffuseband=indrad(inext,nsw,lapproxswupdate)
  ifrso  =indrad(inext,nlev+1,.true.)
  iswfc  =indrad(inext,nlev+1,.true.)
  ifrth  =indrad(inext,nlev+1,.true.)
  ilwfc  =indrad(inext,nlev+1,.true.)
  iaer   =indrad(inext,6*nlev,ldiagforcing)
  ioz    =indrad(inext,nlev,ldiagforcing)
  iico2  =indrad(inext,nlev,ldiagforcing)
  iich4  =indrad(inext,nlev,ldiagforcing)
  iin2o  =indrad(inext,nlev,ldiagforcing)
  ino2   =indrad(inext,nlev,ldiagforcing)
  ic11   =indrad(inext,nlev,ldiagforcing)
  ic12   =indrad(inext,nlev,ldiagforcing)
  ic22   =indrad(inext,nlev,ldiagforcing)
  icl4   =indrad(inext,nlev,ldiagforcing)
  igix   =indrad(inext,1,lldebug)

  ioutend=inext-1                  ! end of output variables

                                 ! start of local variables
  if(.not.ldiagforcing) then
    if (naermacc == 1)  iaero = indrad(inext,iradaer*nlev,.true.)
    iaer   =indrad(inext,nlev*6,.true.)
    ioz    =indrad(inext,nlev,.not.lrayfm)
    iico2  =indrad(inext,nlev,.true.)
    iich4  =indrad(inext,nlev,.true.)
    iin2o  =indrad(inext,nlev,.true.)
    ino2   =indrad(inext,nlev,.true.)
    ic11   =indrad(inext,nlev,.true.)
    ic12   =indrad(inext,nlev,.true.)
    ic22   =indrad(inext,nlev,.true.)
    icl4   =indrad(inext,nlev,.true.)
  endif
                                 ! end of local variables

                                  ! start of standalone inputs workaround variables
  ire_liq =indrad(inext,nlev,.true.)
  ire_ice =indrad(inext,nlev,.true.)
                                  ! end of standalone inputs workaround variables

  ifldsin = iinend - iinbeg +1
  ifldsout= ioutend-ioutbeg +1
  ifldstot= inext  - 1

  if( lldebug )then
    write(nulout,'("imu0   =",i0)')imu0
    write(nulout,'("iamu0  =",i0)')iamu0
    write(nulout,'("iemiss =",i0)')iemiss
    write(nulout,'("its    =",i0)')its
    write(nulout,'("islm   =",i0)')islm
    write(nulout,'("iccnl  =",i0)')iccnl
    write(nulout,'("iccno  =",i0)')iccno
    write(nulout,'("ibas   =",i0)')ibas
    write(nulout,'("itop   =",i0)')itop
    write(nulout,'("igelam =",i0)')igelam
    write(nulout,'("igemu  =",i0)')igemu
    write(nulout,'("iclon  =",i0)')iclon
    write(nulout,'("islon  =",i0)')islon
    write(nulout,'("iald   =",i0)')iald
    write(nulout,'("ialp   =",i0)')ialp
    write(nulout,'("iti    =",i0)')iti
    write(nulout,'("ipr    =",i0)')ipr
    write(nulout,'("iqs    =",i0)')iqs
    write(nulout,'("iwv    =",i0)')iwv
    write(nulout,'("iclc   =",i0)')iclc
    write(nulout,'("ilwa   =",i0)')ilwa
    write(nulout,'("iiwa   =",i0)')iiwa
    write(nulout,'("iswa   =",i0)')iswa
    write(nulout,'("irwa   =",i0)')irwa
    write(nulout,'("irra   =",i0)')irra
    write(nulout,'("idp    =",i0)')idp
    write(nulout,'("ioz    =",i0)')ioz
    write(nulout,'("iecpo3 =",i0)')iecpo3
    write(nulout,'("ihpr   =",i0)')ihpr
    write(nulout,'("iaprs  =",i0)')iaprs
    write(nulout,'("ihti   =",i0)')ihti
    write(nulout,'("ifrsod =",i0)')ifrsod
    write(nulout,'("ifrted=",i0)')ifrted 
    write(nulout,'("ifrsodc=",i0)')ifrsodc 
    write(nulout,'("ifrtedc=",i0)')ifrtedc 
    write(nulout,'("iemit  =",i0)')iemit
    write(nulout,'("isudu  =",i0)')isudu
    write(nulout,'("iuvdf  =",i0)')iuvdf
    write(nulout,'("iparf  =",i0)')iparf
    write(nulout,'("iparcf =",i0)')iparcf
    write(nulout,'("itincf =",i0)')itincf
    write(nulout,'("ifdir  =",i0)')ifdir
    write(nulout,'("ifdif  =",i0)')ifdif
    write(nulout,'("icdir  =",i0)')icdir
    write(nulout,'("ilwderivative  =",i0)')ilwderivative
    write(nulout,'("iswdirectband  =",i0)')iswdirectband
    write(nulout,'("iswdiffuseband =",i0)')iswdiffuseband
    write(nulout,'("ifrso  =",i0)')ifrso
    write(nulout,'("iswfc  =",i0)')iswfc
    write(nulout,'("ifrth  =",i0)')ifrth
    write(nulout,'("ilwfc  =",i0)')ilwfc
    write(nulout,'("igi    =",i0)')igi
    write(nulout,'("iaer   =",i0)')iaer
    write(nulout,'("iaero  =",i0)')iaero
    write(nulout,'("iico2  =",i0)')ico2
    write(nulout,'("iich4  =",i0)')ich4
    write(nulout,'("iin2o  =",i0)')in2o
    write(nulout,'("ino2   =",i0)')ino2
    write(nulout,'("ic11   =",i0)')ic11
    write(nulout,'("ic12   =",i0)')ic12
    write(nulout,'("ic22   =",i0)')ic22
    write(nulout,'("icl4   =",i0)')icl4
    write(nulout,'("ifldsin=",i0)')ifldsin
    write(nulout,'("ifldsout=",i0)')ifldsout
    write(nulout,'("ifldstot=",i0)')ifldstot
    write(nulout,'("ire_liq=",i0)')ire_liq
    write(nulout,'("ire_ice=",i0)')ire_ice
  endif
  
  ! Allocate blocked data structure
  allocate(zrgp(nproma,ifldstot,ngpblks))

  ! First touch
  !$OMP PARALLEL DO SCHEDULE(RUNTIME)&
  !$OMP&PRIVATE(IB,IFLD)
  do ib=1,ngpblks
    do ifld=1,ifldstot
      zrgp(:,ifld,ib) = 0._jprb
    enddo
  enddo
  !$OMP END PARALLEL DO
  !
  !
  !  END IFS SETUP
  !
  !

  ! REPLACED ich4 with iich4 due to clash
  ! REPLACED in2o with iin2o due to clash
  ! REPLACED ico2 with iico2 due to clash

  !  INPUT LOOP

  !$OMP PARALLEL DO SCHEDULE(RUNTIME)&
  !$OMP&PRIVATE(JRL,IBEG,IEND,IL,IB,JAER,JOFF,JLEV,JALB)
  do jrl=1,ncol,nproma

    ibeg=jrl
    iend=min(ibeg+nproma-1,ncol)
    il=iend-ibeg+1
    ib=(jrl-1)/nproma+1

    !* RADINTG:  3.      PREPARE INPUT ARRAYS

    ! zrgp(1:il,imu0,ib)  = ???
    zrgp(1:il,iamu0,ib)  =  single_level%cos_sza(ibeg:iend)   ! cosine of solar zenith ang (mu0)

    do jemiss=1,nlwemiss
      zrgp(1:il,iemiss+jemiss-1,ib)  =  single_level%lw_emissivity(ibeg:iend,jemiss)
    enddo

    zrgp(1:il,its,ib)      = single_level%skin_temperature(ibeg:iend)  ! skin temperature
    zrgp(1:il,islm,ib)     = 0.0_jprb ! Not in NetCDF inputs, see re_liq,re_ice workaround! ! land-sea mask
    zrgp(1:il,iccnl,ib)    = 0.0_jprb ! Not in NetCDF inputs, see re_liq,re_ice workaround! ! CCN over land
    zrgp(1:il,iccno,ib)    = 0.0_jprb ! Not in NetCDF inputs, see re_liq,re_ice workaround! ! CCN over sea
    ! zrgp(1:il,ibas,ib)     = ???
    ! zrgp(1:il,itop,ib)     = ???
    zrgp(1:il,igelam,ib)   = lon(ibeg:iend) ! longitude
    zrgp(1:il,igemu,ib)    = sin(lat(ibeg:iend)) ! sine of latitude
    ! zrgp(1:il,iclon,ib)    = ???
    ! zrgp(1:il,islon,ib)    = ???

    do jalb=1,nsw
      zrgp(1:il,iald+jalb-1,ib)  =  single_level%sw_albedo(ibeg:iend,jalb)
      zrgp(1:il,ialp+jalb-1,ib)  =  single_level%sw_albedo_direct(ibeg:iend,jalb)
    enddo
    
    do jlev=1,nlev
      zrgp(1:il,iti+jlev-1,ib)   = 0.0_jprb ! Not in NetCDF inputs, see re_liq,re_ice workaround and disabled SATUR ! full level temperature
      zrgp(1:il,ipr+jlev-1,ib)   = 0.0_jprb ! Not in NetCDF inputs, see re_liq,re_ice workaround and disabled SATUR ! full level pressure
      ! zrgp(1:il,iqs+jlev-1,ib)   = ???
    enddo

    do jlev=1,nlev
      zrgp(1:il,iwv+jlev-1,ib)   = gas%mixing_ratio(ibeg:iend,jlev,IH2O)
      zrgp(1:il,iclc+jlev-1,ib)  = cloud%fraction(ibeg:iend,jlev)
      zrgp(1:il,ilwa+jlev-1,ib)  = cloud%q_liq(ibeg:iend,jlev)
      zrgp(1:il,iiwa+jlev-1,ib)  = cloud%q_ice(ibeg:iend,jlev)
      zrgp(1:il,iswa+jlev-1,ib)  = 0._jprb  ! snow
      zrgp(1:il,irwa+jlev-1,ib)  = 0._jprb  ! rain

      ! zrgp(1:il,irra+jlev-1,ib)  = ???
      ! zrgp(1:il,idp+jlev-1,ib)   = ???
      ! zrgp(1:il,ifsd+jlev-1,ib)   = ???
      zrgp(1:il,ioz+jlev-1,ib)  = gas%mixing_ratio(ibeg:iend,jlev,IO3) ! extra treatment in rad scheme, to check
      ! zrgp(1:il,iecpo3+jlev-1,ib) = ???
    enddo

    zrgp(1:il,iaer,ib)  =  0._jprb ! old aerosol, not used
    joff=iaero 
    do jaer=1,iradaer
      do jlev=1,nlev
        zrgp(1:il,joff,ib) = aerosol%mixing_ratio(ibeg:iend,jlev,jaer)
        joff=joff+1
      enddo
    enddo

    do jlev=1,nlev+1
      ! zrgp(1:il,ihpr+jlev-1,ib)  = ???
      zrgp(1:il,iaprs+jlev-1,ib) = thermodynamics%pressure_hl(ibeg:iend,jlev)
      zrgp(1:il,ihti+jlev-1,ib)  = thermodynamics%temperature_hl(ibeg:iend,jlev)
    enddo

    ! -- by default, globally averaged concentrations (mmr)
    do jlev=1,nlev
      zrgp(1:il,iico2+jlev-1,ib) = gas%mixing_ratio(ibeg:iend,jlev,ICO2)
      zrgp(1:il,iich4+jlev-1,ib) = gas%mixing_ratio(ibeg:iend,jlev,ICH4)
      zrgp(1:il,iin2o+jlev-1,ib) = gas%mixing_ratio(ibeg:iend,jlev,IN2O)
      zrgp(1:il,ino2+jlev-1,ib)  = 0._jprb ! not really used?

      zrgp(1:il,ic11+jlev-1,ib)  = gas%mixing_ratio(ibeg:iend,jlev,ICFC11)
      zrgp(1:il,ic12+jlev-1,ib)  = gas%mixing_ratio(ibeg:iend,jlev,ICFC12)
      zrgp(1:il,ic22+jlev-1,ib)  = gas%mixing_ratio(ibeg:iend,jlev,IHCFC22)
      zrgp(1:il,icl4+jlev-1,ib)  = gas%mixing_ratio(ibeg:iend,jlev,ICCL4)
    enddo

    ! local workaround variables for standalone input files
    do jlev=1,nlev
      zrgp(1:il,ire_liq+jlev-1,ib) = cloud%re_liq(ibeg:iend,jlev)
      zrgp(1:il,ire_ice+jlev-1,ib) = cloud%re_ice(ibeg:iend,jlev)
    enddo
  enddo
  !$OMP END PARALLEL DO
 
  if (driver_config%iverbose >= 2) then
    write(nulout,'(a)')  'Performing radiative transfer calculations'
  end if
  
  ! Option of repeating calculation multiple time for more accurate
  ! profiling
  do jrepeat = 1,driver_config%nrepeat
      

    !$OMP PARALLEL DO SCHEDULE(RUNTIME)&
    !$OMP&PRIVATE(JRL,IBEG,IEND,IL,IB)
    do jrl=1,ncol,nproma

      ibeg=jrl
      iend=min(ibeg+nproma-1,ncol)
      il=iend-ibeg+1
      ib=(jrl-1)/nproma+1

      call radiation_scheme &
       & (ydmodel, &
       &  1, il, nproma, &                       ! startcol, endcol, ncol
       &  nlev, iradaer, &                       ! nlev, naerosols
       &  zrii0, &                               ! solar_irrad
       &  zrgp(1,iamu0,ib), zrgp(1,its,ib), &    ! mu0, skintemp
       &  zrgp(1,iald,ib), zrgp(1,ialp,ib), &    ! albedo_dif, albedo_dir
       &  zrgp(1,iemiss,ib), &                   ! spectral emissivity
       &  zrgp(1,iccnl,ib), zrgp(1,iccno,ib) ,&  ! CCN concentration, land and sea
       &  zrgp(1,igelam,ib),zrgp(1,igemu,ib), &  ! longitude, sine of latitude
       &  zrgp(1,islm,ib), &                     ! land sea mask
       &  zrgp(1,ipr,ib),   zrgp(1,iti,ib),  &   ! full level pressure and temperature
       &  zrgp(1,iaprs,ib), zrgp(1,ihti,ib), &   ! half-level pressure and temperature
       &  zrgp(1,iwv,ib),   zrgp(1,iico2,ib), zrgp(1,iich4,ib),zrgp(1,iin2o,ib), &
       &  zrgp(1,ino2,ib),  zrgp(1,ic11,ib),  zrgp(1,ic12,ib), zrgp(1,ic22,ib), &
       &  zrgp(1,icl4,ib),  zrgp(1,ioz,ib), &
       &  zrgp(1,iclc,ib),  zrgp(1,ilwa,ib),  zrgp(1,iiwa,ib), zrgp(1,irwa,ib), &
       &  zrgp(1,iswa,ib), &
       &  zrgp(1,iaer,ib),  zrgp(1,iaero,ib), &
       ! flux outputs
       &  zrgp(1,ifrso,ib), zrgp(1,ifrth,ib), zrgp(1,iswfc,ib),zrgp(1,ilwfc,ib),&
       &  zrgp(1,ifrsod,ib),zrgp(1,ifrted,ib), &
       &  zrgp(1,ifrsodc,ib),zrgp(1,ifrtedc,ib),& 
       &  zrgp(1,ifdir,ib), zrgp(1,icdir,ib), zrgp(1,isudu,ib), &
       &  zrgp(1,iuvdf,ib), zrgp(1,iparf,ib), &
       &  zrgp(1,iparcf,ib),zrgp(1,itincf,ib), &
       &  zrgp(1,iemit,ib) ,zrgp(1,ilwderivative,ib), &
       &  zrgp(1,iswdiffuseband,ib), zrgp(1,iswdirectband,ib),&
       ! workaround variables
       &  zrgp(1,ire_liq,ib), zrgp(1,ire_ice,ib), single_level%iseed(ibeg:iend))
    enddo
    !$OMP END PARALLEL DO

  ! ! Store inputs
  ! if (driver_config%do_save_inputs) then
  !   call save_inputs('inputs.nc', ydmodel%yrml_phy_rad%yradiation%rad_config, single_level, thermodynamics, &
  !        &                gas, cloud, aerosol, &
  !        &                lat=spread(0.0_jprb,1,ncol), &
  !        &                lon=spread(0.0_jprb,1,ncol), &
  !        &                iverbose=driver_config%iverbose)
  ! end if

 
      ! Run radiation scheme over blocks of columns in parallel
      
      ! Compute number of blocks to process
      !nblock = (driver_config%iendcol - driver_config%istartcol &
      !     &  + driver_config%nblocksize) / driver_config%nblocksize
     
      !tstart = omp_get_wtime() 
      !!$OMP PARALLEL DO PRIVATE(istartcol, iendcol) SCHEDULE(RUNTIME)
      !do jblock = 1, nblock
      !  ! Specify the range of columns to process.
      !  istartcol = (jblock-1) * driver_config%nblocksize &
      !       &    + driver_config%istartcol
      !  iendcol = min(istartcol + driver_config%nblocksize - 1, &
      !       &        driver_config%iendcol)
      !    
      !  if (driver_config%iverbose >= 3) then
      !    write(nulout,'(a,i0,a,i0,a,i0)')  'Thread ', omp_get_thread_num(), &
      !         &  ' processing columns ', istartcol, '-', iendcol
      !  end if
      !  
      !  if (is_complex_surface) then
      !    call surface_intermediate%calc_boundary_conditions(driver_config%istartcol, &
      !         &  driver_config%iendcol, config, surface, thermodynamics, gas, single_level)
      !  end if
      !  
      !  ! Call the ECRAD radiation scheme
      !  call radiation(ncol, nlev, istartcol, iendcol, config, &
      !       &  single_level, thermodynamics, gas, cloud, aerosol, flux)
        
      !  if (is_complex_surface) then
      !    call surface_intermediate%partition_fluxes(driver_config%istartcol, &
      !         &  driver_config%iendcol, config, surface, flux, surface_flux)
      !  end if
      !  
      !end do
      !!$OMP END PARALLEL DO
      !tstop = omp_get_wtime()
      !write(nulout, '(a,g11.5,a)') 'Time elapsed in radiative transfer: ', tstop-tstart, ' seconds'
      
    
  end do

  deallocate(zrgp)

  ! --------------------------------------------------------
  ! Section 5: Check and save output
  ! --------------------------------------------------------

  is_out_of_bounds = flux%out_of_physical_bounds(driver_config%istartcol, driver_config%iendcol)

  ! Store the fluxes in the output file
  call save_fluxes(file_name, ydmodel%yrml_phy_rad%yradiation%rad_config, thermodynamics, flux, &
       &   iverbose=driver_config%iverbose, is_hdf5_file=driver_config%do_write_hdf5, &
       &   experiment_name=driver_config%experiment_name)
    
!  if (is_complex_surface) then
!    ! Get NetCDF output file name for surface
!    call get_command_argument(4, file_name, status=istatus)
!    if (istatus /= 0) then
!      write(nulout,'(a)') 'Warning: file name for surface-flux outputs not provided'
!    else
!      call save_surface_fluxes(file_name, config, surface_flux, iverbose=driver_config%iverbose)
!    end if
!  end if

  if (driver_config%iverbose >= 2) then
    write(nulout,'(a)') '------------------------------------------------------------------------------------'
  end if

contains

integer(kind=jpim) function indrad(knext,kflds,lduse)
integer(kind=jpim), intent(inout) :: knext
integer(kind=jpim), intent(in) :: kflds
logical, intent(in) :: lduse

if( lduse ) then
  indrad=knext
  knext=knext+kflds
else
  indrad=-99999999
endif

end function indrad

end program ecrad_ifs_driver
