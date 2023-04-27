MODULE RADIATION_SETUP

! RADIATION_SETUP - Setting up modular radiation scheme
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
! PURPOSE
! -------
!   The modular radiation scheme is contained in a separate
!   library. SETUP_RADIATION_SCHEME in this module sets up a small
!   derived type that contains the configuration object for the
!   radiation scheme, plus a small number of additional variables
!   needed for its implemenation in the IFS.
!
! INTERFACE
! ---------
!   SETUP_RADIATION_SCHEME is called from SUECRAD.  The radiation
!   scheme is actually run using the RADIATION_SCHEME routine (not in
!   this module).
!
! AUTHOR
! ------
!   Robin Hogan, ECMWF
!   Original: 2015-09-16
!
! MODIFICATIONS
! -------------
!   2017-03-03  R. Hogan   Put global variables in TRADIATION derived type
!   2017-11-17  S. Remy    Add Nitrates and SOA if NAERMACC=0
!   2017-11-28  R. Hogan   Delta scaling applied to particles only
!   2018-01-11  R. Hogan   Capability to scale solar spectrum in each band
!   2018-04-20  A. Bozzo   Added capability to read in aerosol optical properties
!                          at selected wavelengths
!   2019-01-21  R. Hogan   Explicit albedo and emissivity spectral definitions
!                          leading to smarter weighting in ecRad
!

!-----------------------------------------------------------------------

USE PARKIND1,         ONLY :   JPRB,JPIM
USE radiation_config, ONLY :   config_type, &
       &                       ISolverMcICA, ISolverSpartacus, &
       &                       ISolverTripleclouds, ISolverCloudless, &
       &                       ILiquidModelSlingo, ILiquidModelSOCRATES, &
       &                       IIceModelFu, IIceModelBaran, &
       &                       IOverlapExponential, IOverlapMaximumRandom, &
       &                       IOverlapExponentialRandom
USE YOERAD, ONLY : TERAD

IMPLICIT NONE

SAVE

! Background aerosol is specified in an ugly way: using the old Tegen
! fields that are in terms of optical depth, and converted to mass
! mixing ratio via the relevant mass-extinction coefficient. The
! following are the indices to the aerosol types used to describe
! tropospheric and stratospheric background aerosol.
INTEGER(KIND=JPIM), PARAMETER :: ITYPE_TROP_BG_AER = 8 ! hydrophobic organic
INTEGER(KIND=JPIM), PARAMETER :: ITYPE_STRAT_BG_AER=12 ! non-absorbing sulphate

! This derived type contains configuration information for the
! radiation scheme plus a few additional variables and parameters
! needed for the IFS interface to it
TYPE :: TRADIATION

  ! Configuration for wider aspects of the radiation scheme
  TYPE(TERAD) :: YRERAD

  ! Configuration information for the ecRad radiation scheme
  type(config_type)  :: rad_config

  ! Ultraviolet weightings
  INTEGER(KIND=JPIM) :: NWEIGHT_UV
  INTEGER(KIND=JPIM) :: IBAND_UV(100)
  REAL(KIND=JPRB)    :: WEIGHT_UV(100)
  ! Photosynthetically active radiation weightings
  INTEGER(KIND=JPIM) :: NWEIGHT_PAR
  INTEGER(KIND=JPIM) :: IBAND_PAR(100)
  REAL(KIND=JPRB)    :: WEIGHT_PAR(100)
  ! Mass-extinction coefficient (m2 kg-1) of tropospheric and
  ! stratospheric background aerosol at 550 nm
  REAL(KIND=JPRB)    :: TROP_BG_AER_MASS_EXT
  REAL(KIND=JPRB)    :: STRAT_BG_AER_MASS_EXT

END TYPE TRADIATION

! Dummy type
TYPE :: TCOMPO
  LOGICAL :: LAERNITRATE = .false.
  LOGICAL :: LAERSOA = .false.
END TYPE TCOMPO

CONTAINS

  ! This routine copies information between the IFS radiation
  ! configuration (stored mostly in YDERAD) and the radiation
  ! configuration of the modular radiation scheme (stored in
  ! PRADIATION%rad_config).  The optional input logical LDOUTPUT
  ! controls whether to print lots of information during the setup
  ! stage (default is no).
  SUBROUTINE SETUP_RADIATION_SCHEME(PRADIATION,LDOUTPUT,FILE_NAME)

    USE YOMHOOK,  ONLY : LHOOK, DR_HOOK, JPHOOK
    USE YOMLUN,   ONLY : NULOUT, NULERR
    !USE YOESRTWN, ONLY : NMPSRTM
    USE YOERAD,   ONLY : TERAD
    USE YOEPHY,   ONLY : TEPHY
    !USE YOMCOMPO, ONLY : TCOMPO

    USE RADIATION_INTERFACE,      ONLY : SETUP_RADIATION
    USE RADIATION_AEROSOL_OPTICS, ONLY : DRY_AEROSOL_MASS_EXTINCTION

    ! Radiation configuration information
    TYPE(TCOMPO) :: YDCOMPO
    TYPE(TRADIATION)  ,INTENT(INOUT), TARGET  :: PRADIATION
    CHARACTER(LEN=512),INTENT(IN), OPTIONAL   :: FILE_NAME

    ! Whether or not to print out information on the radiation scheme
    ! configuration
    LOGICAL, INTENT(IN), OPTIONAL :: LDOUTPUT

    ! Verbosity of configuration information 0=none, 1=warning,
    ! 2=info, 3=progress, 4=detailed, 5=debug
    INTEGER(KIND=JPIM) :: IVERBOSESETUP
    !INTEGER(KIND=JPIM) :: ISTAT

    ! Data directory name
    CHARACTER(LEN=256) :: CL_DATA_DIR

    ! Arrays to avoid temporaries
    REAL(KIND=JPRB)    :: ZWAVBOUND(15)
    INTEGER(KIND=JPIM) :: IBAND(16)

    ! Do we use the nearest ecRad band to the albedo/emissivity
    ! intervals, or a more intelligent weighting?
    LOGICAL :: LL_DO_NEAREST_SW_ALBEDO, LL_DO_NEAREST_LW_EMISS

    REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!#include "posname.intfb.h"
#include "abor1.intfb.h"

    IF (LHOOK) CALL DR_HOOK('RADIATION_SETUP:SETUP_RADIATION_SCHEME',0,ZHOOK_HANDLE)

    ! *** GENERAL SETUP ***
    ASSOCIATE(YDERAD=>PRADIATION%YRERAD)
    ASSOCIATE(RAD_CONFIG=>PRADIATION%RAD_CONFIG,&
              & LAERNITRATE=>YDCOMPO%LAERNITRATE, &
              & LAERSOA=>YDCOMPO%LAERSOA, &
              & YSPECTPLANCK=>YDERAD%YSPECTPLANCK)

    ! Configure verbosity of setup of radiation scheme
    IVERBOSESETUP = 4 ! Provide plenty of information
    IF (PRESENT(LDOUTPUT)) THEN
      IF (.NOT. LDOUTPUT) THEN
        IVERBOSESETUP = 1 ! Warnings and errors only
      ENDIF
    ENDIF
    RAD_CONFIG%IVERBOSESETUP = IVERBOSESETUP

    IF (IVERBOSESETUP > 1) THEN
      WRITE(NULOUT,'(a)') '-------------------------------------------------------------------------------'
      WRITE(NULOUT,'(a)') 'RADIATION_SETUP: ecRad 1.5'
    ENDIF

    ! Normal operation of the radiation scheme displays only errors
    ! and warnings
    RAD_CONFIG%IVERBOSE = 1

    ! Read data directory name from the DATA environment variable
    CALL GETENV("DATA", CL_DATA_DIR)
    IF (CL_DATA_DIR /= " ") THEN
      RAD_CONFIG%DIRECTORY_NAME = TRIM(CL_DATA_DIR) // "/ifsdata"
    ELSE
      ! If DATA not present, use the current directory
      RAD_CONFIG%DIRECTORY_NAME = "."
    ENDIF

    ! Do we do Hogan and Bozzo (2015) approximate longwave updates?
    RAD_CONFIG%DO_LW_DERIVATIVES = YDERAD%LAPPROXLWUPDATE

    ! If we are to perform Hogan and Bozzo (2015) approximate
    ! shortwave updates then we need the downwelling direct and
    ! diffuse shortwave fluxes at the surface in each albedo spectral
    ! interval
    RAD_CONFIG%DO_CANOPY_FLUXES_SW = YDERAD%LAPPROXSWUPDATE

    ! If we are to perform approximate longwave updates and we are
    ! using the new 6-interval longwave emissivity scheme then we need
    ! ecRad to compute the downwelling surface longwave fluxes in each
    ! emissivity spectral interval
    IF (YDERAD%NLWOUT > 1) THEN
      RAD_CONFIG%DO_CANOPY_FLUXES_LW = .TRUE.
    ENDIF

    ! Surface spectral fluxes are needed for UV and PAR calculations
    RAD_CONFIG%DO_SURFACE_SW_SPECTRAL_FLUX = .TRUE.


    ! *** SETUP GAS OPTICS ***

    ! Assume IFS has already set-up RRTM, so the setup_radiation
    ! routine below does not have to
    !RAD_CONFIG%DO_SETUP_IFSRRTM = .FALSE.


    ! *** SETUP CLOUD OPTICS ***

    ! Setup liquid optics
    IF (YDERAD%NLIQOPT == 2) THEN
      RAD_CONFIG%I_LIQ_MODEL = ILIQUIDMODELSLINGO
    ELSEIF (YDERAD%NLIQOPT == 4) THEN
      RAD_CONFIG%I_LIQ_MODEL = ILIQUIDMODELSOCRATES
    ELSE
      WRITE(NULERR,'(a,i0)') '*** Error: Unavailable liquid optics model in modular radiation scheme: NLIQOPT=', &
           &  YDERAD%NLIQOPT
      CALL ABOR1('RADIATION_SETUP: error interpreting NLIQOPT')
    ENDIF

    ! Setup ice optics
    IF (YDERAD%NICEOPT == 3) THEN
      RAD_CONFIG%I_ICE_MODEL = IICEMODELFU
      IF (YDERAD%LFU_LW_ICE_OPTICS_BUG) THEN
        RAD_CONFIG%DO_FU_LW_ICE_OPTICS_BUG = .TRUE.
      ENDIF
    ELSEIF (YDERAD%NICEOPT == 4) THEN
      RAD_CONFIG%I_ICE_MODEL = IICEMODELBARAN
    ELSE
      WRITE(NULERR,'(a,i0)') '*** Error: Unavailable ice optics model in modular radiation scheme: NICEOPT=', &
           &  YDERAD%NICEOPT
!!      CALL ABOR1('RADIATION_SETUP: error interpreting NICEOPT')   !db fix
    ENDIF

    ! For consistency with earlier versions of the IFS radiation
    ! scheme, until 45R1 we performed shortwave delta-Eddington
    ! scaling after the merge of the cloud, aerosol and gas optical
    ! properties.  Setting this to "false" does the scaling on the
    ! cloud and aerosol properties separately before merging with
    ! gases, which is more physically appropriate. The impact is very
    ! small (see item 6 of table 2 of Technical Memo 787).
    RAD_CONFIG%DO_SW_DELTA_SCALING_WITH_GASES = .FALSE.


    ! *** SETUP AEROSOLS ***

    RAD_CONFIG%USE_AEROSOLS = .TRUE.

    ! If monochromatic aerosol properties are available they will be
    ! read in automatically so the following is not needed
    !IF (YDEAERATM%LAERRAD) RAD_CONFIG%AEROSOL_OPTICS%READ_MONOCHROMATIC_OPTICS=.TRUE.

    IF (YDERAD%NAERMACC == 1) THEN
      ! Using MACC climatology or prognostic aerosol variables - in
      ! this case the aerosol optics file will be chosen automatically

      ! 12 IFS aerosol classes: 1-3 Sea salt, 4-6 Boucher desert dust,
      ! 7 hydrophilic organics, 8 hydrophobic organics, 9&10
      ! hydrophobic black carbon, 11 ammonium sulphate, 12 inactive
      ! SO2
      RAD_CONFIG%N_AEROSOL_TYPES = 12

      ! Indices to the aerosol optical properties in
      ! aerosol_ifs_rrtm_*.nc, for each class, where negative numbers
      ! index hydrophilic aerosol types and positive numbers index
      ! hydrophobic aerosol types
      RAD_CONFIG%I_AEROSOL_TYPE_MAP = 0 ! There can be up to 256 types
      RAD_CONFIG%I_AEROSOL_TYPE_MAP(1:12) = (/&
           &  -1,&! Sea salt, size bin 1 (OPAC)
           &  -2,&! Sea salt, size bin 2 (OPAC)
           &  -3,&! Sea salt, size bin 3 (OPAC)
           &  7,&! Desert dust, size bin 1 (Woodward 2001)
           &  8,&! Desert dust, size bin 2 (Woodward 2001)
           &  9,&! Desert dust, size bin 3 (Woodward 2001)
           &  -4,&! Hydrophilic organic matter (Hess, OPAC)
           &  10,&! Hydrophobic organic matter (Hess, OPAC)
           &  11,&! Black carbon (Hess, OPAC)
           &  11,&! Black carbon (Hess, OPAC)
           &  -5,&! Ammonium sulphate (GACP)
           &  14 /)  ! Stratospheric sulphate (GACP) [ climatology only ]

      ! Background aerosol mass-extinction coefficients are obtained
      ! after the configuration files have been read - see later in
      ! this routine.

      ! The default aerosol optics file is the following - please
      ! update here, not in radiation/module/radiation_config.F90
      RAD_CONFIG%AEROSOL_OPTICS_OVERRIDE_FILE_NAME = 'aerosol_ifs_rrtm_46R1_with_NI_AM.nc'

    ELSE
      ! Using Tegen climatology
      RAD_CONFIG%N_AEROSOL_TYPES = 6
      RAD_CONFIG%I_AEROSOL_TYPE_MAP = 0 ! There can be up to 256 types
      RAD_CONFIG%I_AEROSOL_TYPE_MAP(1:6) = (/&
           &  1,&! Continental background
           &  2,&! Maritime
           &  3,&! Desert
           &  4,&! Urban
           &  5,&! Volcanic active
           &  6 /)  ! Stratospheric background

      ! Manually set the aerosol optics file name (the directory will
      ! be added automatically)
      RAD_CONFIG%AEROSOL_OPTICS_OVERRIDE_FILE_NAME = 'aerosol_ifs_rrtm_tegen.nc'
    ENDIF

    ! *** SETUP SOLVER ***

    ! 3D effects are off by default
    RAD_CONFIG%DO_3D_EFFECTS = .FALSE.

    ! Select longwave solver
    SELECT CASE (YDERAD%NLWSOLVER)
    CASE(0)
      RAD_CONFIG%I_SOLVER_LW = ISOLVERMCICA
    CASE(1)
      RAD_CONFIG%I_SOLVER_LW = ISOLVERSPARTACUS
    CASE(2)
      RAD_CONFIG%I_SOLVER_LW = ISOLVERSPARTACUS
      RAD_CONFIG%DO_3D_EFFECTS = .TRUE.
    CASE(3)
      RAD_CONFIG%I_SOLVER_LW = ISOLVERTRIPLECLOUDS
    CASE(4)
      RAD_CONFIG%I_SOLVER_LW = ISOLVERCLOUDLESS
    CASE DEFAULT
      WRITE(NULERR,'(a,i0)') '*** Error: Unknown value for NLWSOLVER: ', YDERAD%NLWSOLVER
      CALL ABOR1('RADIATION_SETUP: error interpreting NLWSOLVER')
    END SELECT

    ! Select shortwave solver
    SELECT CASE (YDERAD%NSWSOLVER)
    CASE(0)
      RAD_CONFIG%I_SOLVER_SW = ISOLVERMCICA
    CASE(1)
      RAD_CONFIG%I_SOLVER_SW = ISOLVERSPARTACUS
      RAD_CONFIG%DO_3D_EFFECTS = .FALSE.
      IF (YDERAD%NLWSOLVER == 2) THEN
        CALL ABOR1('RADIATION_SETUP: cannot represent 3D effects in LW but not SW')
      ENDIF
    CASE(2)
      RAD_CONFIG%I_SOLVER_SW = ISOLVERSPARTACUS
      RAD_CONFIG%DO_3D_EFFECTS = .TRUE.
      IF (YDERAD%NLWSOLVER == 1) THEN
        CALL ABOR1('RADIATION_SETUP: cannot represent 3D effects in SW but not LW')
      ENDIF
    CASE(3)
      RAD_CONFIG%I_SOLVER_SW = ISOLVERTRIPLECLOUDS
    CASE(4)
      RAD_CONFIG%I_SOLVER_SW = ISOLVERCLOUDLESS
    CASE DEFAULT
      WRITE(NULERR,'(a,i0)') '*** Error: Unknown value for NSWSOLVER: ', YDERAD%NSWSOLVER
      CALL ABOR1('RADIATION_SETUP: error interpreting NSWSOLVER')
    END SELECT

    ! For stability the cloud effective size can't be too small in
    ! SPARTACUS
    RAD_CONFIG%MIN_CLOUD_EFFECTIVE_SIZE = 500.0_JPRB

    ! SPARTACUS solver requires delta scaling to be done separately
    ! for clouds & aerosols
    IF (RAD_CONFIG%I_SOLVER_SW == ISOLVERSPARTACUS) THEN
      RAD_CONFIG%DO_SW_DELTA_SCALING_WITH_GASES = .FALSE.
    ENDIF

    ! Do we represent longwave scattering?
    RAD_CONFIG%DO_LW_CLOUD_SCATTERING = .FALSE.
    RAD_CONFIG%DO_LW_AEROSOL_SCATTERING = .FALSE.
    SELECT CASE (YDERAD%NLWSCATTERING)
    CASE(1)
      RAD_CONFIG%DO_LW_CLOUD_SCATTERING = .TRUE.
    CASE(2)
      RAD_CONFIG%DO_LW_CLOUD_SCATTERING = .TRUE.
      IF (YDERAD%NAERMACC > 0) THEN
        ! Tegen climatology omits data required to do longwave
        ! scattering by aerosols, so only turn this on with a more
        ! recent scattering database
        RAD_CONFIG%DO_LW_AEROSOL_SCATTERING = .TRUE.
      ENDIF
    END SELECT

    SELECT CASE (YDERAD%NCLOUDOVERLAP)
    CASE (1)
      RAD_CONFIG%I_OVERLAP_SCHEME = IOVERLAPMAXIMUMRANDOM
    CASE (2)
      ! Use Exponential-Exponential cloud overlap to match original IFS
      ! implementation of Raisanen cloud generator
      RAD_CONFIG%I_OVERLAP_SCHEME = IOVERLAPEXPONENTIAL
    CASE (3)
      RAD_CONFIG%I_OVERLAP_SCHEME = IOVERLAPEXPONENTIALRANDOM
    CASE DEFAULT
      WRITE(NULERR,'(a,i0)') '*** Error: Unknown value for NCLOUDOVERLAP: ', YDERAD%NCLOUDOVERLAP
      CALL ABOR1('RADIATION_SETUP: error interpreting NCLOUDOVERLAP')
    END SELECT

    ! Change cloud overlap to exponential-random if Tripleclouds or
    ! SPARTACUS selected as both the shortwave and longwave solvers
    IF (RAD_CONFIG%I_OVERLAP_SCHEME /= IOVERLAPEXPONENTIALRANDOM &
         & .AND. (     RAD_CONFIG%I_SOLVER_SW == ISOLVERTRIPLECLOUDS &
         &        .OR. RAD_CONFIG%I_SOLVER_LW == ISOLVERTRIPLECLOUDS &
         &        .OR. RAD_CONFIG%I_SOLVER_SW == ISOLVERSPARTACUS &
         &        .OR. RAD_CONFIG%I_SOLVER_LW == ISOLVERSPARTACUS)) THEN
      IF (RAD_CONFIG%I_SOLVER_SW == RAD_CONFIG%I_SOLVER_LW) THEN
        WRITE(NULOUT,'(a)') 'Warning: Tripleclouds/SPARTACUS solver selected so changing cloud overlap to Exp-Ran'
        RAD_CONFIG%I_OVERLAP_SCHEME = IOVERLAPEXPONENTIALRANDOM
      ELSE
        ! If the solvers are not the same and exponential-random has
        ! not been selected then abort
        WRITE(NULERR,'(a)') '*** Error: Tripleclouds and SPARTACUS solvers can only simulate exponential-random overlap'
        CALL ABOR1('RADIATION_SETUP: Cloud overlap incompatible with solver')
      ENDIF

      ! For additional stability in SPARTACUS solver it helps if the
      ! cloud fraction threshold is higher than the default of 1.0e-6
      ! used for McICA; this is done for Tripleclouds too so that it
      ! is a good control for SPARTACUS.
      RAD_CONFIG%CLOUD_FRACTION_THRESHOLD = 2.5E-5_JPRB
    ENDIF


    ! Number of longwave surface emissivity intervals to use:
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
    CALL YDERAD%YSPECTPLANCK%INIT(2, [ 8.0E-6_JPRB, 13.0E-6_JPRB ], &
         &  [ 1,2,1 ])
      
    ! Populate the mapping between the 14 RRTM shortwave bands and the
    ! 6 albedo inputs.
    YDERAD%NSW = 6
    ZWAVBOUND(1:5) = [ 0.25e-6_jprb, 0.44e-6_jprb, 0.69e-6_jprb, &
         &             1.19e-6_jprb, 2.38e-6_jprb ]
    IBAND(1:6)  = [ 1,2,3,4,5,6 ]
    ! If NALBEDOSCHEME==2 then we are using the 6-component MODIS
    ! albedo climatology, and a weighted average is used to compute
    ! the albedos in each ecRad spectral band. If NALBEDOSCHEME==3
    ! then we use the diffuse part of the 4 components but still with
    ! a weighted average. Otherwise the older behaviour is followed:
    ! the nearest albedo interval to each band is selected, resulting
    ! in a discrete mapping that matches the one in YOESRTWN:NMPSRTM.
    ! Note that this tends to bias albedo high because there is a lot
    ! of energy around the interface between the UV-Vis and Near-IR
    ! channels, so this should be close to the 0.7 microns intended by
    ! the MODIS dataset, not shifted to the nearest RRTM band boundary
    ! at 0.625 microns.
    LL_DO_NEAREST_SW_ALBEDO = .FALSE.
    CALL RAD_CONFIG%DEFINE_SW_ALBEDO_INTERVALS(YDERAD%NSW, ZWAVBOUND, IBAND, &
         &  DO_NEAREST=LL_DO_NEAREST_SW_ALBEDO)

    ! Likewise between the 16 RRTM longwave bands and the NLWEMISS
    ! emissivity inputs - these are defined in suecrad.F90.
    LL_DO_NEAREST_LW_EMISS = .TRUE.
    CALL RAD_CONFIG%DEFINE_LW_EMISS_INTERVALS(UBOUND(YSPECTPLANCK%INTERVAL_MAP,1), &
         &  YSPECTPLANCK%WAVLEN_BOUND, YSPECTPLANCK%INTERVAL_MAP, &
         &  DO_NEAREST=LL_DO_NEAREST_LW_EMISS)

    ! Do we scale the incoming solar radiation in each band?
    IF (YDERAD%NSOLARSPECTRUM == 1) THEN
      IF (RAD_CONFIG%N_BANDS_SW /= 14) THEN
        WRITE(NULERR,'(a)') '*** Error: Shortwave must have 14 bands to apply spectral scaling'
        CALL ABOR1('RADIATION_SETUP: Shortwave must have 14 bands to apply spectral scaling')
      ELSE
        RAD_CONFIG%USE_SPECTRAL_SOLAR_SCALING = .TRUE.
      ENDIF
    ENDIF

    ! *** IMPLEMENT SETTINGS ***

    ! For advanced configuration, the configuration data for the
    ! "radiation" project can specified directly in the namelist.
    ! However, the variable naming convention is not consistent with
    ! the rest of the IFS.  For basic configuration there are specific
    ! variables in the NAERAD namelist available in the YDERAD
    ! structure.
    !CALL POSNAME(NULNAM, 'RADIATION', ISTAT)
    !SELECT CASE (ISTAT)
    !  CASE(0)
    !    CALL RAD_CONFIG%READ(UNIT=NULNAM)
    !  CASE(1)
    !    WRITE(NULOUT,'(a)') 'Namelist RADIATION not found, using settings from NAERAD only'
    !  CASE DEFAULT
    !    CALL ABOR1('RADIATION_SETUP: error reading RADIATION section of namelist file')
    !END SELECT
    IF (PRESENT(FILE_NAME)) THEN
      CALL RAD_CONFIG%READ(FILE_NAME=FILE_NAME)
    ENDIF

    ! Print configuration
    IF (IVERBOSESETUP > 1) THEN
      WRITE(NULOUT,'(a)') 'Radiation scheme settings:'
      CALL RAD_CONFIG%PRINT(IVERBOSE=IVERBOSESETUP)
    ENDIF

    ! Use configuration data to set-up radiation scheme, including
    ! reading scattering datafiles
    CALL SETUP_RADIATION(RAD_CONFIG)

    ! Get spectral weightings for UV and PAR
    CALL RAD_CONFIG%GET_SW_WEIGHTS(0.2E-6_JPRB, 0.4415E-6_JPRB,&
         &  PRADIATION%NWEIGHT_UV, PRADIATION%IBAND_UV, PRADIATION%WEIGHT_UV,&
         &  'ultraviolet')
    CALL RAD_CONFIG%GET_SW_WEIGHTS(0.4E-6_JPRB, 0.7E-6_JPRB,&
         &  PRADIATION%NWEIGHT_PAR, PRADIATION%IBAND_PAR, PRADIATION%WEIGHT_PAR,&
         &  'photosynthetically active radiation, PAR')

    IF (YDERAD%NAERMACC > 0) THEN
      ! With the MACC aerosol climatology we need to add in the
      ! background aerosol afterwards using the Tegen arrays.  In this
      ! case we first configure the background aerosol mass-extinction
      ! coefficient at 550 nm, which corresponds to the 10th RRTMG
      ! shortwave band.
      PRADIATION%TROP_BG_AER_MASS_EXT  = DRY_AEROSOL_MASS_EXTINCTION(RAD_CONFIG,&
           &                                   ITYPE_TROP_BG_AER, 550.0E-9_JPRB)
      PRADIATION%STRAT_BG_AER_MASS_EXT = DRY_AEROSOL_MASS_EXTINCTION(RAD_CONFIG,&
           &                                   ITYPE_STRAT_BG_AER, 550.0E-9_JPRB)

      WRITE(NULOUT,'(a,i0)') 'Tropospheric background uses aerosol type ',&
           &                 ITYPE_TROP_BG_AER
      WRITE(NULOUT,'(a,i0)') 'Stratospheric background uses aerosol type ',&
           &                 ITYPE_STRAT_BG_AER
    ELSE
      PRADIATION%TROP_BG_AER_MASS_EXT  = 0.0_JPRB
      PRADIATION%STRAT_BG_AER_MASS_EXT = 0.0_JPRB
    ENDIF

    IF (IVERBOSESETUP > 1) THEN
      WRITE(NULOUT,'(a)') '-------------------------------------------------------------------------------'
    ENDIF

    END ASSOCIATE
    END ASSOCIATE

    IF (LHOOK) CALL DR_HOOK('RADIATION_SETUP:SETUP_RADIATION_SCHEME',1,ZHOOK_HANDLE)

  END SUBROUTINE SETUP_RADIATION_SCHEME


END MODULE RADIATION_SETUP
