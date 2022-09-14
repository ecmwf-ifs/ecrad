! ecrad_ifs.F90 - Methods to set-up and run ecRAD in IFS configuration
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

module ecrad_ifs

  use parkind1,                 only : jprb, jprd, jpim ! Working/double precision, integer type

  implicit none

  public

  type :: ifs_config_type
    ! Some stray config values
    logical         :: lrayfm  ! key for calling radiation scheme from MF
    real(kind=jprb) :: zrii0   ! solar irradiance
    ! Offsets in ZRGP
    integer :: igi, imu0, iamu0, iemiss, its, islm, iccnl,    &
        &     ibas, itop, igelam, igemu, iclon, islon, iald, ialp, iti, ipr, iqs, iwv, iclc, ilwa,    &
        &     iiwa, iswa, irwa, irra, idp, ioz, iecpo3, ihpr, iaprs, ihti, iaero, ifrsod, icdir,      &
        &     ifrted, ifrsodc, ifrtedc, iemit, isudu, iuvdf, iparf, iparcf, itincf, ifdir, ifdif,     &
        &     ilwderivative, iswdirectband, iswdiffuseband, ifrso, iswfc, ifrth, ilwfc, iaer,         &
        &     iich4, iin2o, ino2, ic11, ic12, igix, iico2, iccno, ic22, icl4
    integer :: ire_liq, ire_ice, ioverlap
  end type ifs_config_type

contains

subroutine ecrad_ifs_setup(nml_file_name, driver_config, config, ydmodel, ncol)

  use radiation_io,             only : nulout
  use radiation_config,         only : config_type, &
       &                       ISolverMcICA, ISolverSpartacus, &
       &                       ISolverTripleclouds, &
       &                       IOverlapExponential, IOverlapMaximumRandom, &
       &                       IOverlapExponentialRandom
  use radiation_setup,          only : tcompo, setup_radiation_scheme
  use ecrad_driver_config,      only : driver_config_type
  use type_model,               only : model

  implicit none

  ! Name of namelist file name
  character(len=512), intent(in)          :: nml_file_name

  ! Configuration specific to this driver
  type(driver_config_type), intent(in)     :: driver_config

  ! Derived types for the inputs to the radiation scheme
  type(config_type), intent(in)            :: config

  ! Model type
  type(model), intent(inout)               :: ydmodel

  integer, intent(in) :: ncol         ! Number of columns

  ! Dummy COMPO type
  type(tcompo)              :: ydcompo

  ! number of column blocks, block size
  integer :: ngpblks, nproma
  logical :: lldebug

  associate( &
    & ydephy=>ydmodel%yrml_phy_ec%yrephy, &
    & yderad=>ydmodel%yrml_phy_rad%yrerad )

    !  -------------------------------------------------------
    !
    !    IFS SETUP  -  EXCERPT FROM SUECRAD AND RADINTG
    !
    !  -------------------------------------------------------

    lldebug=(driver_config%iverbose>4)     ! debug
    nproma=driver_config%nblocksize        ! nproma size
    ngpblks=(ncol-1)/nproma+1              ! number of column blocks

    write(nulout,'("ncol = ",i0,", nproma = ",i0,", ngpblks = ",i0)') ncol, nproma, ngpblks

    !
    ! Values from IFS  ! Unclear, ask Robin!
    !
    ydephy%nemissscheme = 0
    ydephy%nalbedoscheme = 2

    !
    ! Hard-coded in SUECRAD  ! Do they need to be configurable? Ask Robin!
    !
    yderad%lapproxlwupdate = .true.  ! Hogan and Bozzo (2015) approx longwave updates
    yderad%lapproxswupdate = .true.  ! Hogan and Bozzo (2015) approx shortwave updates
    yderad%nliqopt = 4  ! 2 - SLINGO, 4 - SOCRATES
    yderad%niceopt = 3  ! 3 - ICEMODELFU, 4 - IDEMODELBARAN
    yderad%nsolarspectrum = 0
    yderad%ndecolat = 2  ! DECORRELATION LENGTH FOR CF AND CW, 0: SPECIFIED INDEPENDENT OF LATITUDE, 1: SHONK-HOGAN, 2: IMPROVED
    yderad%rcloud_frac_std = 1.0_jprb
    yderad%ldiagforcing = .false.  ! write input ozone, ghg and aerosol forcing to 3D fields

    !
    ! Reconstruct IFS configuration values from ecRad's config
    !
    yderad%lfu_lw_ice_optics_bug = config%do_fu_lw_ice_optics_bug

    ! 0 - McICA, 1 - SPARTACUS, 2 - SPARTACUS 3D, 3 - Tripleclouds
    if(config%i_solver_lw == ISolverMcICA) then
      yderad%nlwsolver = 0
    elseif(config%i_solver_lw == ISolverSpartacus) then
      yderad%nlwsolver = 1
    elseif(config%i_solver_lw == ISolverTripleclouds) then
      yderad%nlwsolver = 3
    else
      call abor1('Unknown longwave solver type')
    endif

    ! 0 - McICA, 1 - SPARTACUS, 2 - SPARTACUS 3D, 3 - Tripleclouds
    if(config%i_solver_sw == ISolverMcICA) then
      yderad%nswsolver = 0
    elseif(config%i_solver_sw == ISolverSpartacus) then
      yderad%nswsolver = 1
    elseif(config%i_solver_sw == ISolverTripleclouds) then
      yderad%nswsolver = 3
    else
      call abor1('Unknown shortwave solver type')
    endif

    ! 1 - cloud scattering, 2 - cloud + aerosol scattering
    if(config%do_lw_cloud_scattering) then
      if(config%do_lw_aerosol_scattering) then
        yderad%nlwscattering = 2
      else
        yderad%nlwscattering = 1
      endif
    else
      yderad%nlwscattering = 0
    endif

    ! 1 - MAXIMUMRANDOM, 2 - EXPONENTIAL, 3 - EXPONENTIALRANDOM
    if(config%i_overlap_scheme == IOverlapMaximumRandom) then
      yderad%ncloudoverlap = 1
    elseif(config%i_overlap_scheme == IOverlapExponential) then
      yderad%ncloudoverlap = 2
    elseif(config%i_overlap_scheme == IOverlapExponentialRandom) then
      yderad%ncloudoverlap = 3
    else
      call abor1('Unknown cloud overlap scheme')
    endif

    ! number of shortwave spectral intervals, = 6 in IFS
    if (config%use_canopy_full_spectrum_sw) then
      yderad%nsw = config%n_g_sw
    else
      yderad%nsw = 6
    endif

    if (config%use_aerosols) then
      yderad%naermacc = 1  ! MACC-derived aerosol climatology on a NMCLAT x NMCLON grid
    else
      yderad%naermacc = 0
    endif

    !
    ! SUECRAD
    !

    ! Number of longwave surface emissivity intervals to use
    IF (YDEPHY%NEMISSSCHEME == 1) THEN
      ! We do a more accurate mapping for emissivity if NEMISSSCHEME==1.
      ! See susrad_mod.F90 for values for different surface types.
      YDERAD%NLWEMISS = 6
      IF (YDERAD%LAPPROXLWUPDATE) THEN
        ! Pass the same number of longwave downwelling surface spectral
        ! fluxes from ecRad to RADHEATN so that longwave approximate
        ! update scheme can be as accurate as possible
        YDERAD%NLWOUT = 6
      ELSE
        YDERAD%NLWOUT = 1
      ENDIF
      ! Create a spectral Planck look-up table, used by RADHEATN.  Note
      ! that this routine makes use of the length of its third argument.
      ! The following wavelength bounds (metres) match the RRTM band
      ! boundaries.
      CALL YDERAD%YSPECTPLANCK%INIT(6, &
          &  [ 8.4746E-6_JPRB, 10.2041E-6_JPRB, 12.1951E-6_JPRB, 15.8730E-6_JPRB, 28.5714E-6_JPRB ], &
          &  [ 1,2,3,4,5,6 ])
    ELSEIF (YDEPHY%NEMISSSCHEME == 0) THEN
      ! Traditional approach: one value of emissivty for parts of the
      ! spectrum on either side of the infrared atmospheric window
      ! (PEMIR), and one value for the window itself (PEMIW)
      YDERAD%NLWEMISS = 2
      ! ...and the longwave approximate update scheme uses a single
      ! broadband emissivity
      YDERAD%NLWOUT   = 1
      ! Create a spectral Planck look-up table, used by RADHEATN.  Note
      ! that this routine makes use of the length of its third argument.
      ! The wavelength bounds (metres) allow for the first emissivity to
      ! represent values outside the infrared atmospheric window, and the
      ! second emissivity to represent values within it.
      CALL YDERAD%YSPECTPLANCK%INIT(2, [ 8.0E-6_JPRB, 13.0E-6_JPRB ], [ 1,2,1 ])
    ELSE
      CALL ABOR1('RADIATION_SETUP: NEMISSSCHEME must be 0 or 1')
    ENDIF

    call setup_radiation_scheme( &
      & ydmodel%yrml_phy_rad%yrerdi, ydmodel%yrml_phy_rad%yreaeratm, ydcompo, &
      & ydephy, yderad, ydmodel%yrml_phy_rad%yradiation, &
      & ldoutput=.true., file_name=nml_file_name)

  end associate

end subroutine ecrad_ifs_setup

subroutine ecrad_ifs_interpolate_in ( &
  & driver_config, ifs_config, ncol, nlev, &
  & single_level, thermodynamics, gas, cloud, aerosol, &
  & ydmodel, zrgp, iseed, thermodynamics_out)

  use radiation_io,             only : nulout
  use radiation_config,         only : config_type
  use radiation_single_level,   only : single_level_type
  use radiation_thermodynamics, only : thermodynamics_type
  use radiation_gas,            only : gas_type, IMassMixingRatio, &
       &   IH2O, ICO2, IO3, IN2O, ICH4, ICFC11, ICFC12, IHCFC22, ICCL4
  use radiation_cloud,          only : cloud_type
  use radiation_aerosol,        only : aerosol_type
  use ecrad_driver_config,      only : driver_config_type
  use type_model,               only : model

  implicit none

  ! Configuration specific to this driver
  type(driver_config_type), intent(in)     :: driver_config

  type(ifs_config_type), intent(inout)     :: ifs_config

  ! Derived types for the inputs to the radiation scheme
  type(single_level_type)   :: single_level
  type(thermodynamics_type) :: thermodynamics
  type(gas_type)            :: gas
  type(cloud_type)          :: cloud
  type(aerosol_type)        :: aerosol

  ! Model type
  type(model), intent(inout)               :: ydmodel

  integer, intent(in) :: ncol, nlev         ! Number of columns and levels

  ! monolithic IFS data structure to pass to radiation scheme
  real(kind=jprb), intent(out), allocatable :: zrgp(:,:,:)

  ! Seed for random number generator
  integer, intent(out), allocatable         :: iseed(:,:)

  ! Empty thermodynamics type to store pressure_hl for output at the end
  type(thermodynamics_type), intent(inout)  :: thermodynamics_out

  ! number of column blocks, block size
  integer :: ngpblks, nproma

  integer :: ifldsin, ifldsout, ifldstot, inext, iinbeg, iinend, ioutbeg, ioutend
  integer :: jrl, ibeg, iend, il, ib, ifld, jemiss, jalb, jlev, joff, jaer

  logical :: lldebug


  ! Extract some config values
  lldebug=(driver_config%iverbose>4)     ! debug
  ifs_config%lrayfm=.false.
  ifs_config%zrii0=single_level%solar_irradiance
  nproma=driver_config%nblocksize        ! nproma size
  ngpblks=(ncol-1)/nproma+1              ! number of column blocks

  associate( &
    & yderad=>ydmodel%yrml_phy_rad%yrerad, &
    & rad_config=>ydmodel%yrml_phy_rad%yradiation%rad_config )
     !
    ! RADINTG
    !

    !  INITIALISE INDICES FOR VARIABLE

    ! INDRAD is a CONTAIN'd function (now a module function)

    inext  =1
    iinbeg =1                        ! start of input variables
    ifs_config%igi    =indrad(inext,1,lldebug)
    ifs_config%imu0   =indrad(inext,1,.true.)
    ifs_config%iamu0  =indrad(inext,1,.true.)
    ifs_config%iemiss =indrad(inext,yderad%nlwemiss,.true.)
    ifs_config%its    =indrad(inext,1,.true.)
    ifs_config%islm   =indrad(inext,1,.true.)
    ifs_config%iccnl  =indrad(inext,1,.true.)
    ifs_config%iccno  =indrad(inext,1,.true.)
    ifs_config%ibas   =indrad(inext,1,.true.)
    ifs_config%itop   =indrad(inext,1,.true.)
    ifs_config%igelam =indrad(inext,1,.true.)
    ifs_config%igemu  =indrad(inext,1,.true.)
    ifs_config%iclon  =indrad(inext,1,.true.)
    ifs_config%islon  =indrad(inext,1,.true.)
    ifs_config%iald   =indrad(inext,yderad%nsw,.true.)
    ifs_config%ialp   =indrad(inext,yderad%nsw,.true.)
    ifs_config%iti    =indrad(inext,nlev,.true.)
    ifs_config%ipr    =indrad(inext,nlev,.true.)
    ifs_config%iqs    =indrad(inext,nlev,.true.)
    ifs_config%iwv    =indrad(inext,nlev,.true.)
    ifs_config%iclc   =indrad(inext,nlev,.true.)
    ifs_config%ilwa   =indrad(inext,nlev,.true.)
    ifs_config%iiwa   =indrad(inext,nlev,.true.)
    ifs_config%iswa   =indrad(inext,nlev,.true.)
    ifs_config%irwa   =indrad(inext,nlev,.true.)
    ifs_config%irra   =indrad(inext,nlev,.true.)
    ifs_config%idp    =indrad(inext,nlev,.true.)
    ifs_config%ioz    =indrad(inext,nlev,ifs_config%lrayfm)
    ifs_config%iecpo3 =indrad(inext,nlev ,.false.)
    ifs_config%ihpr   =indrad(inext,nlev+1,.true.) ! not used in ecrad
    ifs_config%iaprs  =indrad(inext,nlev+1,.true.)
    ifs_config%ihti   =indrad(inext,nlev+1,.true.)
    ifs_config%iaero  =indrad(inext,rad_config%n_aerosol_types*nlev,.false.)

    iinend =inext-1                  ! end of input variables

    ioutbeg=inext                    ! start of output variables
    if (yderad%naermacc == 1) ifs_config%iaero = indrad(inext,rad_config%n_aerosol_types*nlev,yderad%ldiagforcing)
    ifs_config%ifrsod =indrad(inext,1,.true.)
    ifs_config%ifrted =indrad(inext,yderad%nlwout,.true.)
    ifs_config%ifrsodc=indrad(inext,1,.true.)
    ifs_config%ifrtedc=indrad(inext,1,.true.)
    ifs_config%iemit  =indrad(inext,1,.true.)
    ifs_config%isudu  =indrad(inext,1,.true.)
    ifs_config%iuvdf  =indrad(inext,1,.true.)
    ifs_config%iparf  =indrad(inext,1,.true.)
    ifs_config%iparcf =indrad(inext,1,.true.)
    ifs_config%itincf =indrad(inext,1,.true.)
    ifs_config%ifdir  =indrad(inext,1,.true.)
    ifs_config%ifdif  =indrad(inext,1,.true.)
    ifs_config%icdir  =indrad(inext,1,.true.)
    ifs_config%ilwderivative =indrad(inext,nlev+1, yderad%lapproxlwupdate)
    ifs_config%iswdirectband =indrad(inext,yderad%nsw,yderad%lapproxswupdate)
    ifs_config%iswdiffuseband=indrad(inext,yderad%nsw,yderad%lapproxswupdate)
    ifs_config%ifrso  =indrad(inext,nlev+1,.true.)
    ifs_config%iswfc  =indrad(inext,nlev+1,.true.)
    ifs_config%ifrth  =indrad(inext,nlev+1,.true.)
    ifs_config%ilwfc  =indrad(inext,nlev+1,.true.)
    ifs_config%iaer   =indrad(inext,6*nlev,yderad%ldiagforcing)
    ifs_config%ioz    =indrad(inext,nlev,yderad%ldiagforcing)
    ifs_config%iico2  =indrad(inext,nlev,yderad%ldiagforcing)
    ifs_config%iich4  =indrad(inext,nlev,yderad%ldiagforcing)
    ifs_config%iin2o  =indrad(inext,nlev,yderad%ldiagforcing)
    ifs_config%ino2   =indrad(inext,nlev,yderad%ldiagforcing)
    ifs_config%ic11   =indrad(inext,nlev,yderad%ldiagforcing)
    ifs_config%ic12   =indrad(inext,nlev,yderad%ldiagforcing)
    ifs_config%ic22   =indrad(inext,nlev,yderad%ldiagforcing)
    ifs_config%icl4   =indrad(inext,nlev,yderad%ldiagforcing)
    ifs_config%igix   =indrad(inext,1,lldebug)

    ioutend=inext-1                  ! end of output variables

                                  ! start of local variables
    if(.not.yderad%ldiagforcing) then
      if (yderad%naermacc == 1)  ifs_config%iaero = indrad(inext,rad_config%n_aerosol_types*nlev,.true.)
      ifs_config%iaer   =indrad(inext,nlev*6,.true.)
      ifs_config%ioz    =indrad(inext,nlev,.not.ifs_config%lrayfm)
      ifs_config%iico2  =indrad(inext,nlev,.true.)
      ifs_config%iich4  =indrad(inext,nlev,.true.)
      ifs_config%iin2o  =indrad(inext,nlev,.true.)
      ifs_config%ino2   =indrad(inext,nlev,.true.)
      ifs_config%ic11   =indrad(inext,nlev,.true.)
      ifs_config%ic12   =indrad(inext,nlev,.true.)
      ifs_config%ic22   =indrad(inext,nlev,.true.)
      ifs_config%icl4   =indrad(inext,nlev,.true.)
    endif
                                  ! end of local variables

                                    ! start of standalone inputs workaround variables
    ifs_config%ire_liq =indrad(inext,nlev,.true.)
    ifs_config%ire_ice =indrad(inext,nlev,.true.)
    ifs_config%ioverlap =indrad(inext,nlev-1,.true.)
                                    ! end of standalone inputs workaround variables

    ifldsin = iinend - iinbeg +1
    ifldsout= ioutend-ioutbeg +1
    ifldstot= inext  - 1

    if( lldebug )then
      call ydmodel%PRINT("ECRAD")
      write(nulout,'("imu0   =",i0)')ifs_config%imu0
      write(nulout,'("iamu0  =",i0)')ifs_config%iamu0
      write(nulout,'("iemiss =",i0)')ifs_config%iemiss
      write(nulout,'("its    =",i0)')ifs_config%its
      write(nulout,'("islm   =",i0)')ifs_config%islm
      write(nulout,'("iccnl  =",i0)')ifs_config%iccnl
      write(nulout,'("iccno  =",i0)')ifs_config%iccno
      write(nulout,'("ibas   =",i0)')ifs_config%ibas
      write(nulout,'("itop   =",i0)')ifs_config%itop
      write(nulout,'("igelam =",i0)')ifs_config%igelam
      write(nulout,'("igemu  =",i0)')ifs_config%igemu
      write(nulout,'("iclon  =",i0)')ifs_config%iclon
      write(nulout,'("islon  =",i0)')ifs_config%islon
      write(nulout,'("iald   =",i0)')ifs_config%iald
      write(nulout,'("ialp   =",i0)')ifs_config%ialp
      write(nulout,'("iti    =",i0)')ifs_config%iti
      write(nulout,'("ipr    =",i0)')ifs_config%ipr
      write(nulout,'("iqs    =",i0)')ifs_config%iqs
      write(nulout,'("iwv    =",i0)')ifs_config%iwv
      write(nulout,'("iclc   =",i0)')ifs_config%iclc
      write(nulout,'("ilwa   =",i0)')ifs_config%ilwa
      write(nulout,'("iiwa   =",i0)')ifs_config%iiwa
      write(nulout,'("iswa   =",i0)')ifs_config%iswa
      write(nulout,'("irwa   =",i0)')ifs_config%irwa
      write(nulout,'("irra   =",i0)')ifs_config%irra
      write(nulout,'("idp    =",i0)')ifs_config%idp
      write(nulout,'("ioz    =",i0)')ifs_config%ioz
      write(nulout,'("iecpo3 =",i0)')ifs_config%iecpo3
      write(nulout,'("ihpr   =",i0)')ifs_config%ihpr
      write(nulout,'("iaprs  =",i0)')ifs_config%iaprs
      write(nulout,'("ihti   =",i0)')ifs_config%ihti
      write(nulout,'("ifrsod =",i0)')ifs_config%ifrsod
      write(nulout,'("ifrted =",i0)')ifs_config%ifrted
      write(nulout,'("ifrsodc=",i0)')ifs_config%ifrsodc
      write(nulout,'("ifrtedc=",i0)')ifs_config%ifrtedc
      write(nulout,'("iemit  =",i0)')ifs_config%iemit
      write(nulout,'("isudu  =",i0)')ifs_config%isudu
      write(nulout,'("iuvdf  =",i0)')ifs_config%iuvdf
      write(nulout,'("iparf  =",i0)')ifs_config%iparf
      write(nulout,'("iparcf =",i0)')ifs_config%iparcf
      write(nulout,'("itincf =",i0)')ifs_config%itincf
      write(nulout,'("ifdir  =",i0)')ifs_config%ifdir
      write(nulout,'("ifdif  =",i0)')ifs_config%ifdif
      write(nulout,'("icdir  =",i0)')ifs_config%icdir
      write(nulout,'("ilwderivative  =",i0)')ifs_config%ilwderivative
      write(nulout,'("iswdirectband  =",i0)')ifs_config%iswdirectband
      write(nulout,'("iswdiffuseband =",i0)')ifs_config%iswdiffuseband
      write(nulout,'("ifrso  =",i0)')ifs_config%ifrso
      write(nulout,'("iswfc  =",i0)')ifs_config%iswfc
      write(nulout,'("ifrth  =",i0)')ifs_config%ifrth
      write(nulout,'("ilwfc  =",i0)')ifs_config%ilwfc
      write(nulout,'("igi    =",i0)')ifs_config%igi
      write(nulout,'("iaer   =",i0)')ifs_config%iaer
      write(nulout,'("iaero  =",i0)')ifs_config%iaero
      write(nulout,'("iico2  =",i0)')ifs_config%iico2
      write(nulout,'("iich4  =",i0)')ifs_config%iich4
      write(nulout,'("iin2o  =",i0)')ifs_config%iin2o
      write(nulout,'("ino2   =",i0)')ifs_config%ino2
      write(nulout,'("ic11   =",i0)')ifs_config%ic11
      write(nulout,'("ic12   =",i0)')ifs_config%ic12
      write(nulout,'("ic22   =",i0)')ifs_config%ic22
      write(nulout,'("icl4   =",i0)')ifs_config%icl4
      write(nulout,'("ire_liq=",i0)')ifs_config%ire_liq
      write(nulout,'("ire_ice=",i0)')ifs_config%ire_ice
      write(nulout,'("ioverlap=",i0)')ifs_config%ioverlap
      write(nulout,'("ifldsin =",i0)')ifldsin
      write(nulout,'("ifldsout=",i0)')ifldsout
      write(nulout,'("ifldstot=",i0)')ifldstot
    endif

    ! Allocate blocked data structure
    allocate(zrgp(nproma,ifldstot,ngpblks))
    allocate(iseed(nproma,ngpblks))
    allocate(thermodynamics_out%pressure_hl(ncol,nlev+1))

    ! First touch
    !$OMP PARALLEL DO SCHEDULE(RUNTIME)&
    !$OMP&PRIVATE(IB,IFLD)
    do ib=1,ngpblks
      do ifld=1,ifldstot
        zrgp(:,ifld,ib) = 0._jprb
      enddo
      iseed(:,ib) = 0
    enddo
    !$OMP END PARALLEL DO

    !  -------------------------------------------------------
    !
    !  END IFS SETUP
    !
    !  -------------------------------------------------------

    ! REPLACED ich4 with iich4 due to clash
    ! REPLACED in2o with iin2o due to clash
    ! REPLACED ico2 with iico2 due to clash

    !  -------------------------------------------------------
    !
    !  INPUT LOOP
    !
    !  -------------------------------------------------------

    !$OMP PARALLEL DO SCHEDULE(RUNTIME)&
    !$OMP&PRIVATE(JRL,IBEG,IEND,IL,IB,JAER,JOFF,JLEV,JALB)
    do jrl=1,ncol,nproma

      ibeg=jrl
      iend=min(ibeg+nproma-1,ncol)
      il=iend-ibeg+1
      ib=(jrl-1)/nproma+1

      !* RADINTG:  3.      PREPARE INPUT ARRAYS

      ! zrgp(1:il,imu0,ib)  = ???
      zrgp(1:il,ifs_config%iamu0,ib)  =  single_level%cos_sza(ibeg:iend)   ! cosine of solar zenith ang (mu0)

      do jemiss=1,yderad%nlwemiss
        zrgp(1:il,ifs_config%iemiss+jemiss-1,ib)  =  single_level%lw_emissivity(ibeg:iend,jemiss)
      enddo

      zrgp(1:il,ifs_config%its,ib)      = single_level%skin_temperature(ibeg:iend)  ! skin temperature
      zrgp(1:il,ifs_config%islm,ib)     = 0.0_jprb ! Not in NetCDF inputs, see re_liq,re_ice workaround! ! land-sea mask
      zrgp(1:il,ifs_config%iccnl,ib)    = 0.0_jprb ! Not in NetCDF inputs, see re_liq,re_ice workaround! ! CCN over land
      zrgp(1:il,ifs_config%iccno,ib)    = 0.0_jprb ! Not in NetCDF inputs, see re_liq,re_ice workaround! ! CCN over sea
      ! zrgp(1:il,ibas,ib)     = ???
      ! zrgp(1:il,itop,ib)     = ???
      zrgp(1:il,ifs_config%igelam,ib)   = 0.0_jprb ! lon(ibeg:iend) ! longitude
      zrgp(1:il,ifs_config%igemu,ib)    = 0.0_jprb ! sin(lat(ibeg:iend)) ! sine of latitude
      ! zrgp(1:il,iclon,ib)    = ???
      ! zrgp(1:il,islon,ib)    = ???

      do jalb=1,yderad%nsw
        zrgp(1:il,ifs_config%iald+jalb-1,ib)  =  single_level%sw_albedo(ibeg:iend,jalb)
        zrgp(1:il,ifs_config%ialp+jalb-1,ib)  =  single_level%sw_albedo_direct(ibeg:iend,jalb)
      enddo

      do jlev=1,nlev
        zrgp(1:il,ifs_config%iti+jlev-1,ib)   = 0.0_jprb ! Not in NetCDF inputs, see re_liq,re_ice workaround and disabled SATUR ! full level temperature
        zrgp(1:il,ifs_config%ipr+jlev-1,ib)   = 0.0_jprb ! Not in NetCDF inputs, see re_liq,re_ice workaround and disabled SATUR ! full level pressure
        ! zrgp(1:il,iqs+jlev-1,ib)   = ???
      enddo

      do jlev=1,nlev
        zrgp(1:il,ifs_config%iwv+jlev-1,ib)   = gas%mixing_ratio(ibeg:iend,jlev,IH2O) ! this is already in MassMixingRatio units
        zrgp(1:il,ifs_config%iclc+jlev-1,ib)  = cloud%fraction(ibeg:iend,jlev)
        zrgp(1:il,ifs_config%ilwa+jlev-1,ib)  = cloud%q_liq(ibeg:iend,jlev)
        zrgp(1:il,ifs_config%iiwa+jlev-1,ib)  = cloud%q_ice(ibeg:iend,jlev)
        zrgp(1:il,ifs_config%iswa+jlev-1,ib)  = 0._jprb  ! snow
        zrgp(1:il,ifs_config%irwa+jlev-1,ib)  = 0._jprb  ! rain

        ! zrgp(1:il,irra+jlev-1,ib)  = ???
        ! zrgp(1:il,idp+jlev-1,ib)   = ???
        ! zrgp(1:il,ifsd+jlev-1,ib)   = ???
        ! zrgp(1:il,iecpo3+jlev-1,ib) = ???
      enddo

      zrgp(1:il,ifs_config%iaer,ib)  =  0._jprb ! old aerosol, not used
      if (yderad%naermacc == 1) then
        joff=ifs_config%iaero
        do jaer=1,rad_config%n_aerosol_types
          do jlev=1,nlev
            zrgp(1:il,joff,ib) = aerosol%mixing_ratio(ibeg:iend,jlev,jaer)
            joff=joff+1
          enddo
        enddo
      endif

      do jlev=1,nlev+1
        ! zrgp(1:il,ihpr+jlev-1,ib)  = ???
        zrgp(1:il,ifs_config%iaprs+jlev-1,ib) = thermodynamics%pressure_hl(ibeg:iend,jlev)
        zrgp(1:il,ifs_config%ihti+jlev-1,ib)  = thermodynamics%temperature_hl(ibeg:iend,jlev)
      enddo

      ! -- by default, globally averaged concentrations (mmr)
      call gas%get(ICO2, IMassMixingRatio, zrgp(1:il,ifs_config%iico2:ifs_config%iico2+nlev-1,ib), istartcol=ibeg)
      call gas%get(ICH4, IMassMixingRatio, zrgp(1:il,ifs_config%iich4:ifs_config%iich4+nlev-1,ib), istartcol=ibeg)
      call gas%get(IN2O, IMassMixingRatio, zrgp(1:il,ifs_config%iin2o:ifs_config%iin2o+nlev-1,ib), istartcol=ibeg)

      call gas%get(ICFC11, IMassMixingRatio, zrgp(1:il,ifs_config%ic11:ifs_config%ic11+nlev-1,ib), istartcol=ibeg)
      call gas%get(ICFC12, IMassMixingRatio, zrgp(1:il,ifs_config%ic12:ifs_config%ic12+nlev-1,ib), istartcol=ibeg)
      call gas%get(IHCFC22,IMassMixingRatio, zrgp(1:il,ifs_config%ic22:ifs_config%ic22+nlev-1,ib), istartcol=ibeg)
      call gas%get(ICCL4,  IMassMixingRatio, zrgp(1:il,ifs_config%icl4:ifs_config%icl4+nlev-1,ib), istartcol=ibeg)

      call gas%get(IO3, IMassMixingRatio, zrgp(1:il,ifs_config%ioz:ifs_config%ioz+nlev-1,ib), istartcol=ibeg)
      ! convert ozone kg/kg to Pa*kg/kg
      do jlev=1,nlev
        zrgp(1:il,ifs_config%ioz+jlev-1,ib)  = zrgp(1:il,ifs_config%ioz+jlev-1,ib) &
              &                       * (thermodynamics%pressure_hl(ibeg:iend,jlev+1) &
              &                         - thermodynamics%pressure_hl(ibeg:iend,jlev))
      enddo

      ! local workaround variables for standalone input files
      do jlev=1,nlev
        ! missing full-level temperature and pressure as well as land-sea-mask
        zrgp(1:il,ifs_config%ire_liq+jlev-1,ib) = cloud%re_liq(ibeg:iend,jlev)
        zrgp(1:il,ifs_config%ire_ice+jlev-1,ib) = cloud%re_ice(ibeg:iend,jlev)
      enddo
      do jlev=1,nlev-1
        ! for the love of it, I can't figure this one out. Probably to do with
        ! my crude approach of setting PGEMU?
        zrgp(1:il,ifs_config%ioverlap+jlev-1,ib) = cloud%overlap_param(ibeg:iend,jlev)
      enddo
      iseed(1:il,ib) = single_level%iseed(ibeg:iend)
    enddo
    !$OMP END PARALLEL DO

    ! Store pressure for output
    thermodynamics_out%pressure_hl(:,:) = thermodynamics%pressure_hl(:,:)

  end associate

end subroutine ecrad_ifs_interpolate_in


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


subroutine ecrad_ifs_run (driver_config, ifs_config, ncol, nlev, ydmodel, iseed, zrgp, flux_blocks)

  use radiation_io,             only : nulout
  use radiation_flux,           only : flux_type
  use ecrad_driver_config,      only : driver_config_type
  use type_model,               only : model

  implicit none

#include "radiation_scheme.intfb.h"

  ! Configuration specific to this driver
  type(driver_config_type), intent(in)        :: driver_config

  type(ifs_config_type), intent(in)           :: ifs_config

  integer, intent(in) :: ncol, nlev         ! Number of columns and levels

  ! Model type
  type(model), intent(inout)                  :: ydmodel

  ! Seed for random number generator
  integer, allocatable, intent(inout)         :: iseed(:,:)

  ! monolithic IFS data structure to pass to radiation scheme
  real(kind=jprb), intent(inout), allocatable :: zrgp(:,:,:)

  ! Per-block flux data structure to validate outputs
  type(flux_type), intent(inout), allocatable :: flux_blocks(:)

  ! number of column blocks, block size
  integer :: ngpblks, nproma

  ! For parallel processing of multiple blocks
  !integer :: jblock, nblock ! Block loop index and number
  integer, external :: omp_get_thread_num
  double precision, external :: omp_get_wtime

  ! Start/stop time in seconds
  real(kind=jprd) :: tstart, tstop

  ! Loop index for repeats (for benchmarking)
  integer :: jrepeat

  integer :: jrl, ibeg, iend, ib, il, naerosols

  nproma=driver_config%nblocksize        ! nproma size
  ngpblks=(ncol-1)/nproma+1              ! number of column blocks
  naerosols=ydmodel%yrml_phy_rad%yradiation%rad_config%n_aerosol_types

  allocate(flux_blocks(ngpblks))

  ! Option of repeating calculation multiple time for more accurate
  ! profiling
  do jrepeat = 1,driver_config%nrepeat

    tstart = omp_get_wtime()

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
      &  nlev, naerosols, &    ! nlev, naerosols
      &  ifs_config%zrii0, &                               ! solar_irrad
      &  zrgp(1,ifs_config%iamu0,ib), zrgp(1,ifs_config%its,ib), &    ! mu0, skintemp
      &  zrgp(1,ifs_config%iald,ib) , zrgp(1,ifs_config%ialp,ib), &    ! albedo_dif, albedo_dir
      &  zrgp(1,ifs_config%iemiss,ib), &                   ! spectral emissivity
      &  zrgp(1,ifs_config%iccnl,ib), zrgp(1,ifs_config%iccno,ib) ,&  ! CCN concentration, land and sea
      &  zrgp(1,ifs_config%igelam,ib),zrgp(1,ifs_config%igemu,ib), &  ! longitude, sine of latitude
      &  zrgp(1,ifs_config%islm,ib), &                     ! land sea mask
      &  zrgp(1,ifs_config%ipr,ib),   zrgp(1,ifs_config%iti,ib),  &   ! full level pressure and temperature
      &  zrgp(1,ifs_config%iaprs,ib), zrgp(1,ifs_config%ihti,ib), &   ! half-level pressure and temperature
      &  zrgp(1,ifs_config%iwv,ib),   zrgp(1,ifs_config%iico2,ib), &
      &  zrgp(1,ifs_config%iich4,ib), zrgp(1,ifs_config%iin2o,ib), &
      &  zrgp(1,ifs_config%ino2,ib),  zrgp(1,ifs_config%ic11,ib), &
      &  zrgp(1,ifs_config%ic12,ib),  zrgp(1,ifs_config%ic22,ib), &
      &  zrgp(1,ifs_config%icl4,ib),  zrgp(1,ifs_config%ioz,ib), &
      &  zrgp(1,ifs_config%iclc,ib),  zrgp(1,ifs_config%ilwa,ib), &
      &  zrgp(1,ifs_config%iiwa,ib),  zrgp(1,ifs_config%irwa,ib), &
      &  zrgp(1,ifs_config%iswa,ib), &
      &  zrgp(1,ifs_config%iaer,ib),  zrgp(1,ifs_config%iaero,ib), &
      ! flux outputs
      &  zrgp(1,ifs_config%ifrso,ib), zrgp(1,ifs_config%ifrth,ib), &
      &  zrgp(1,ifs_config%iswfc,ib), zrgp(1,ifs_config%ilwfc,ib),&
      &  zrgp(1,ifs_config%ifrsod,ib),zrgp(1,ifs_config%ifrted,ib), &
      &  zrgp(1,ifs_config%ifrsodc,ib),zrgp(1,ifs_config%ifrtedc,ib),&
      &  zrgp(1,ifs_config%ifdir,ib), zrgp(1,ifs_config%icdir,ib), &
      &  zrgp(1,ifs_config%isudu,ib), &
      &  zrgp(1,ifs_config%iuvdf,ib), zrgp(1,ifs_config%iparf,ib), &
      &  zrgp(1,ifs_config%iparcf,ib),zrgp(1,ifs_config%itincf,ib), &
      &  zrgp(1,ifs_config%iemit,ib) ,zrgp(1,ifs_config%ilwderivative,ib), &
      &  zrgp(1,ifs_config%iswdiffuseband,ib), zrgp(1,ifs_config%iswdirectband,ib),&
      ! workaround variables
      &  zrgp(1,ifs_config%ire_liq,ib), zrgp(1,ifs_config%ire_ice,ib), iseed(1,ib),&
      &  zrgp(1,ifs_config%ioverlap,ib), flux_blocks(ib))
    enddo
    !$OMP END PARALLEL DO

    tstop = omp_get_wtime()
    write(nulout, '(a,g11.5,a)') 'Time elapsed in radiative transfer: ', tstop-tstart, ' seconds'
  end do

  deallocate(zrgp)
  deallocate(iseed)

end subroutine ecrad_ifs_run

subroutine ecrad_ifs_interpolate_out (driver_config, ncol, nlev, ydmodel, flux_blocks, flux_out)

  use radiation_flux,           only : flux_type
  use ecrad_driver_config,      only : driver_config_type
  use type_model,               only : model

  implicit none

  ! Configuration specific to this driver
  type(driver_config_type), intent(in)     :: driver_config

  ! Model type
  type(model), intent(in)                  :: ydmodel

  ! Per-block flux data structure to validate outputs
  type(flux_type), allocatable, intent(inout) :: flux_blocks(:)

  ! Derived type containing outputs from the radiation scheme
  type(flux_type), intent(inout)              :: flux_out

  integer, intent(in) :: ncol, nlev         ! Number of columns and levels

  integer :: nproma

  integer :: jrl, ibeg, iend, il, ib, jalb, jlev, joff, jaer

  nproma=driver_config%nblocksize        ! nproma size

  associate( rad_config=>ydmodel%yrml_phy_rad%yradiation%rad_config )

    !  -------------------------------------------------------
    !
    !  OUTPUT LOOP
    !
    !  -------------------------------------------------------

    ! Allocate memory for the flux profiles, which may include arrays
    ! of dimension n_bands_sw/n_bands_lw, so must be called after
    ! setup_radiation
    call flux_out%allocate(rad_config, 1, ncol, nlev)

    !$OMP PARALLEL DO SCHEDULE(RUNTIME)&
    !$OMP&PRIVATE(JRL,IBEG,IEND,IL,IB,JAER,JOFF,JLEV,JALB)
    do jrl=1,ncol,nproma

      ibeg=jrl
      iend=min(ibeg+nproma-1,ncol)
      il=iend-ibeg+1
      ib=(jrl-1)/nproma+1

      if (rad_config%do_lw) then
        flux_out%lw_up(ibeg:iend,:) = flux_blocks(ib)%lw_up(1:il,:)
        flux_out%lw_dn(ibeg:iend,:) = flux_blocks(ib)%lw_dn(1:il,:)

        if (rad_config%do_clear) then
          flux_out%lw_up_clear(ibeg:iend,:) = flux_blocks(ib)%lw_up_clear(1:il,:)
          flux_out%lw_dn_clear(ibeg:iend,:) = flux_blocks(ib)%lw_dn_clear(1:il,:)
        endif

        if (rad_config%do_lw_derivatives) then
          flux_out%lw_derivatives(ibeg:iend,:) = flux_blocks(ib)%lw_derivatives(1:il,:)
        endif

        if (rad_config%do_save_spectral_flux) then
          flux_out%lw_up_band(:,ibeg:iend,:) = flux_blocks(ib)%lw_up_band(:,1:il,:)
          flux_out%lw_dn_band(:,ibeg:iend,:) = flux_blocks(ib)%lw_dn_band(:,1:il,:)
          if (rad_config%do_clear) then
            flux_out%lw_up_clear_band(:,ibeg:iend,:) = flux_blocks(ib)%lw_up_clear_band(:,1:il,:)
            flux_out%lw_dn_clear_band(:,ibeg:iend,:) = flux_blocks(ib)%lw_dn_clear_band(:,1:il,:)
          endif
        endif

        if (rad_config%do_canopy_fluxes_lw) then
          flux_out%lw_dn_surf_canopy(:,ibeg:iend) = flux_blocks(ib)%lw_dn_surf_canopy(:,1:il)
        endif
      endif

      if (rad_config%do_sw) then
        flux_out%sw_up(ibeg:iend,:) = flux_blocks(ib)%sw_up(1:il,:)
        flux_out%sw_dn(ibeg:iend,:) = flux_blocks(ib)%sw_dn(1:il,:)
        if (rad_config%do_sw_direct) then
          flux_out%sw_dn_direct(ibeg:iend,:) = flux_blocks(ib)%sw_dn_direct(1:il,:)
        endif

        if (rad_config%do_clear) then
          flux_out%sw_up_clear(ibeg:iend,:) = flux_blocks(ib)%sw_up_clear(1:il,:)
          flux_out%sw_dn_clear(ibeg:iend,:) = flux_blocks(ib)%sw_dn_clear(1:il,:)
          if (rad_config%do_sw_direct) then
            flux_out%sw_dn_direct_clear(ibeg:iend,:) = flux_blocks(ib)%sw_dn_direct_clear(1:il,:)
          endif
        endif

        if (rad_config%do_save_spectral_flux) then
          flux_out%sw_up_band(:,ibeg:iend,:) = flux_blocks(ib)%sw_up_band(:,1:il,:)
          flux_out%sw_dn_band(:,ibeg:iend,:) = flux_blocks(ib)%sw_dn_band(:,1:il,:)

          if (rad_config%do_sw_direct) then
            flux_out%sw_dn_direct_band(:,ibeg:iend,:) = flux_blocks(ib)%sw_dn_direct_band(:,1:il,:)
          endif

          if (rad_config%do_clear) then
            flux_out%sw_up_clear_band(:,ibeg:iend,:) = flux_blocks(ib)%sw_up_clear_band(:,1:il,:)
            flux_out%sw_dn_clear_band(:,ibeg:iend,:) = flux_blocks(ib)%sw_dn_clear_band(:,1:il,:)
            if (rad_config%do_sw_direct) then
              flux_out%sw_dn_direct_clear_band(:,ibeg:iend,:) = flux_blocks(ib)%sw_dn_direct_band(:,1:il,:)
            endif
          endif

        else if (rad_config%do_surface_sw_spectral_flux) then
          flux_out%sw_dn_surf_band(:,ibeg:iend) = flux_blocks(ib)%sw_dn_surf_band(:,1:il)
          flux_out%sw_dn_direct_surf_band(:,ibeg:iend) = flux_blocks(ib)%sw_dn_direct_surf_band(:,1:il)

          if (rad_config%do_clear) then
            flux_out%sw_dn_surf_clear_band(:,ibeg:iend) = flux_blocks(ib)%sw_dn_surf_clear_band(:,1:il)
            flux_out%sw_dn_direct_surf_clear_band(:,ibeg:iend) = flux_blocks(ib)%sw_dn_direct_surf_clear_band(:,1:il)
          endif
        endif

        if (rad_config%do_canopy_fluxes_sw) then
          flux_out%sw_dn_direct_surf_canopy(:,ibeg:iend) = flux_blocks(ib)%sw_dn_direct_surf_canopy(:,1:il)
          flux_out%sw_dn_diffuse_surf_canopy(:,ibeg:iend) = flux_blocks(ib)%sw_dn_diffuse_surf_canopy(:,1:il)
        endif
      endif

      if (rad_config%do_lw .and. rad_config%do_clouds) then
        flux_out%cloud_cover_lw(ibeg:iend) = flux_blocks(ib)%cloud_cover_lw(1:il)
      endif
      if (rad_config%do_sw .and. rad_config%do_clouds) then
        flux_out%cloud_cover_sw(ibeg:iend) = flux_blocks(ib)%cloud_cover_sw(1:il)
      endif

      call flux_blocks(ib)%deallocate
    end do

  end associate

  deallocate(flux_blocks)
end subroutine ecrad_ifs_interpolate_out

end module ecrad_ifs
