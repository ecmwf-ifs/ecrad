MODULE RADIATION_SETUP

! RADIATION_SETUP - Setting up modular radiation scheme
!
! PURPOSE
! -------
!   The modular radiation scheme is contained in a separate
!   library. SETUP_RADIATION_SCHEME in this module sets up a small
!   number of global variables needed to store the information for it.
!
!   Lower case is used for variables and types taken from the
!   radiation library
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
!
!-----------------------------------------------------------------------

  USE PARKIND1,         ONLY : JPRB
  USE radiation_config, ONLY : config_type, &
       &                       ISolverMcICA, ISolverSpartacus, &
       &                       ILiquidModelSlingo, ILiquidModelSOCRATES, &
       &                       IIceModelFu, IIceModelBaran, &
       &                       IOverlapExponential

  IMPLICIT NONE

  ! Store configuration information for the radiation scheme in a
  ! global variable
  type(config_type) :: rad_config

  ! Ultraviolet weightings
  INTEGER         :: NWEIGHT_UV
  INTEGER         :: IBAND_UV(100)
  REAL(KIND=JPRB) :: WEIGHT_UV(100)
  ! Photosynthetically active radiation weightings
  INTEGER         :: NWEIGHT_PAR
  INTEGER         :: IBAND_PAR(100)
  REAL(KIND=JPRB) :: WEIGHT_PAR(100)

  ! Background aerosol is specified in an ugly way: using the old
  ! Tegen fields that are in terms of optical depth, and converted to
  ! mass mixing ratio via the relevant mass-extinction coefficient
  INTEGER, PARAMETER :: ITYPE_TROP_BG_AER = 8 ! hydrophobic organic
  INTEGER, PARAMETER :: ITYPE_STRAT_BG_AER=12 ! non-absorbing sulphate
  REAL(KIND=JPRB)    :: TROP_BG_AER_MASS_EXT
  REAL(KIND=JPRB)    :: STRAT_BG_AER_MASS_EXT

CONTAINS

  ! This routine copies information between the IFS radiation
  ! configuration (stored in global variables) and the radiation
  ! configuration of the modular radiation scheme (stored in
  ! rad_config).  The optional input logical LOUTPUT controls whether
  ! to print lots of information during the setup stage (default is
  ! no).
  SUBROUTINE SETUP_RADIATION_SCHEME(LOUTPUT)

    USE YOMHOOK,  ONLY : LHOOK, DR_HOOK
    USE YOMLUN,   ONLY : NULNAM, NULOUT, NULERR
    USE YOESRTWN, ONLY : NMPSRTM
    USE YOERAD,   ONLY : YRERAD

    USE radiation_interface,      ONLY : setup_radiation
    USE radiation_aerosol_optics, ONLY : dry_aerosol_sw_mass_extinction

#include "posname.intfb.h"

    ! Whether or not to provide information on the radiation scheme
    ! configuration
    LOGICAL, INTENT(IN), OPTIONAL :: LOUTPUT

    ! Verbosity of configuration information 0=none, 1=warning,
    ! 2=info, 3=progress, 4=detailed, 5=debug
    INTEGER :: IVERBOSESETUP
    INTEGER :: ISTAT

    REAL(KIND=JPRB) :: ZHOOK_HANDLE

    IF (LHOOK) CALL DR_HOOK('RADIATION_SETUP:SETUP_RADIATION_SCHEME',0,ZHOOK_HANDLE)

    ! *** GENERAL SETUP ***

    ! Configure verbosity of setup of radiation scheme
    IVERBOSESETUP = 4 ! Provide plenty of information
    IF (PRESENT(LOUTPUT)) THEN
      IF (.NOT. LOUTPUT) THEN
        IVERBOSESETUP = 1 ! Warnings and errors only
      ENDIF
    ENDIF
    rad_config%iverbosesetup = IVERBOSESETUP

    IF (IVERBOSESETUP > 1) THEN
      WRITE(NULOUT,'(a)') '-------------------------------------------------------------------------------'
      WRITE(NULOUT,'(a)') 'RADIATION_SETUP'
    ENDIF

    ! Normal operation of the radiation scheme displays only errors
    ! and warnings
    rad_config%iverbose = 1

    ! For the time being, ensure a valid default directory name
    rad_config%directory_name = '/home/rd/parr/radiation_data'

    ! Do we do Hogan and Bozzo (2014) approximate longwave updates?
    rad_config%do_lw_derivatives = YRERAD%LAPPROXLWUPDATE

    ! Surface spectral fluxes are needed for spectral shortwave albedo
    ! calculation
    rad_config%do_surface_sw_spectral_flux = .TRUE.


    ! *** SETUP GAS OPTICS ***

    ! Assume IFS has already set-up RRTM, so the setup_radiation
    ! routine below does not have to
    rad_config%do_setup_ifsrrtm = .FALSE.


    ! *** SETUP CLOUD OPTICS ***

    ! Setup liquid optics
    IF (YRERAD%NLIQOPT == 2) THEN
      rad_config%i_liq_model = ILiquidModelSlingo
    ELSEIF (YRERAD%NLIQOPT == 3) THEN
      rad_config%i_liq_model = ILiquidModelSOCRATES
    ELSE
      WRITE(NULERR,'(a,i0)') 'Unavailable liquid optics model in modular radiation scheme: NLIQOPT=', &
           &  YRERAD%NLIQOPT
      CALL ABOR1('RADIATION_SETUP: error interpreting NLIQOPT')   
    ENDIF

    ! Setup ice optics
    IF (YRERAD%NICEOPT == 3) THEN
      rad_config%i_ice_model = IIceModelFu
    ELSEIF (YRERAD%NICEOPT == 4) THEN
      rad_config%i_ice_model = IIceModelBaran
    ELSE
      WRITE(NULERR,'(a,i0)') 'Unavailable ice optics model in modular radiation scheme: NICEOPT=', &
           &  YRERAD%NICEOPT
      CALL ABOR1('RADIATION_SETUP: error interpreting NICEOPT')   
    ENDIF

    ! For consistency with earlier versions of the IFS radiation
    ! scheme, we perform shortwave delta-Eddington scaling *after* the
    ! merge of the cloud, aerosol and gas optical properties.  Set
    ! this to "false" to do the scaling on the cloud and aerosol
    ! properties separately before merging with gases. Note that this
    ! is not compatible with the SPARTACUS solver.
    rad_config%do_sw_delta_scaling_with_gases = .TRUE.

    ! Use Exponential-Exponential cloud overlap to match original IFS
    ! implementation of Raisanen cloud generator
    rad_config%i_overlap_scheme = IOverlapExponential


    ! *** SETUP AEROSOLS ***

    rad_config%use_aerosols = .TRUE.

    IF (YRERAD%NAERMACC > 0) THEN
      ! Using MACC climatology - in this case the aerosol optics file
      ! will be chosen automatically

      ! 12 IFS aerosol classes: 1-3 Sea salt, 4-6 Boucher desert dust,
      ! 7 hydrophilic organics, 8 hydrophobic organics, 9&10
      ! hydrophobic black carbon, 11 ammonium sulphate, 12 inactive
      ! SO2
      rad_config%n_aerosol_types = 12

      ! Indices to the aerosol optical properties in
      ! aerosol_ifs_rrtm_*.nc, for each class, where negative numbers
      ! index hydrophilic aerosol types and positive numbers index
      ! hydrophobic aerosol types
      rad_config%i_aerosol_type_map = 0 ! There can be up to 256 types
      rad_config%i_aerosol_type_map(1:12) = (/ &
           &  -1, &  ! Sea salt, size bin 1 (OPAC)
           &  -2, &  ! Sea salt, size bin 2 (OPAC)
           &  -3, &  ! Sea salt, size bin 3 (OPAC)
           &   7, &  ! Desert dust, size bin 1 (Woodward 2001)
           &   8, &  ! Desert dust, size bin 2 (Woodward 2001)
           &   9, &  ! Desert dust, size bin 3 (Woodward 2001)
           &  -4, &  ! Hydrophilic organic matter (OPAC)
           &  10, &  ! Hydrophobic organic matter (OPAC)
           &  11, &  ! Black carbon (Boucher)
           &  11, &  ! Black carbon (Boucher)
           &  -5, &  ! Ammonium sulphate (OPAC)
           &  14 /)  ! Stratospheric sulphate (hand edited from OPAC)

      ! Background aerosol mass-extinction coefficients are obtained
      ! after the configuration files have been read - see later in
      ! this routine.

    ELSE
      ! Using Tegen climatology
      rad_config%n_aerosol_types = 6
      rad_config%i_aerosol_type_map = 0 ! There can be up to 256 types
      rad_config%i_aerosol_type_map(1:6) = (/ &
           &  1, &  ! Continental background
           &  2, &  ! Maritime
           &  3, &  ! Desert
           &  4, &  ! Urban
           &  5, &  ! Volcanic active
           &  6 /)  ! Stratospheric background

      ! Manually set the aerosol optics file name (the directory will
      ! be added automatically)
      rad_config%aerosol_optics_override_file_name = 'aerosol_ifs_rrtm_tegen.nc'
    ENDIF

    ! *** SETUP SOLVER ***

    ! 3D effects are off by default
    rad_config%do_3d_effects = .FALSE.

    ! Select longwave solver
    SELECT CASE (YRERAD%NLWSOLVER)
    CASE(0)
      rad_config%i_solver_lw = ISolverMcICA
    CASE(1)
      rad_config%i_solver_lw = ISolverSpartacus
    CASE(2)
      rad_config%i_solver_lw = ISolverSpartacus
      rad_config%do_3d_effects = .TRUE.
    CASE DEFAULT
      WRITE(NULERR,'(a,i0)') 'Unknown value for NLWSOLVER: ', YRERAD%NLWSOLVER
      CALL ABOR1('RADIATION_SETUP: error interpreting NLWSOLVER')
    END SELECT

    ! Select shortwave solver
    SELECT CASE (YRERAD%NSWSOLVER)
    CASE(0)
      rad_config%i_solver_sw = ISolverMcICA
    CASE(1)
      rad_config%i_solver_sw = ISolverSpartacus
      rad_config%do_3d_effects = .FALSE.
      IF (YRERAD%NLWSOLVER == 2) THEN
        CALL ABOR1('RADIATION_SETUP: cannot represent 3D effects in LW but not SW')
      ENDIF
    CASE(2)
      rad_config%i_solver_sw = ISolverSpartacus
      rad_config%do_3d_effects = .TRUE.
      IF (YRERAD%NLWSOLVER == 1) THEN
        CALL ABOR1('RADIATION_SETUP: cannot represent 3D effects in SW but not LW')
      ENDIF
    CASE DEFAULT
      WRITE(NULERR,'(a,i0)') 'Unknown value for NSWSOLVER: ', YRERAD%NSWSOLVER
      CALL ABOR1('RADIATION_SETUP: error interpreting NSWSOLVER')
    END SELECT

    ! SPARTACUS solver requires delta scaling to be done separately
    ! for clouds & aerosols
    IF (rad_config%i_solver_sw == ISolverSpartacus) THEN
      rad_config%do_sw_delta_scaling_with_gases = .FALSE.
    ENDIF

    ! Do we represent longwave scattering?
    rad_config%do_lw_cloud_scattering = .FALSE.
    rad_config%do_lw_aerosol_scattering = .FALSE.
    SELECT CASE (YRERAD%NLWSCATTERING)
    CASE(1)
      rad_config%do_lw_cloud_scattering = .TRUE.
    CASE(2)
      rad_config%do_lw_cloud_scattering = .TRUE.
      IF (YRERAD%NAERMACC > 0) THEN
        ! Tegen climatology omits data required to do longwave
        ! scattering by aerosols, so only turn this on with a more
        ! recent scattering database
        rad_config%do_lw_aerosol_scattering = .TRUE.
      ENDIF
    END SELECT


    ! *** IMPLEMENT SETTINGS ***

    ! For advanced configuration, the configuration data for the
    ! "radiation" project can specified directly in the namelist.
    ! However, the variable naming convention is not consistent with
    ! the rest of the IFS.  For basic configuration there are specific
    ! variables in the NAERAD namelist available in the YRERAD
    ! structure.
    CALL POSNAME(NULNAM, 'RADIATION', ISTAT)
    SELECT CASE (ISTAT)
      CASE(0)
        CALL rad_config%read(unit=NULNAM)
      CASE(1)
        WRITE(NULOUT,'(a)') 'Namelist RADIATION not found, using settings from NAERAD only'
      CASE DEFAULT
        CALL ABOR1('RADIATION_SETUP: error reading RADIATION section of namelist file')
    END SELECT

    ! Print configuration
    IF (IVERBOSESETUP > 1) THEN
      WRITE(NULOUT,'(a)') 'Radiation scheme settings:'
      CALL rad_config%print(IVERBOSE=IVERBOSESETUP)
    ENDIF

    ! Use configuration data to set-up radiation scheme, including
    ! reading scattering datafiles
    CALL setup_radiation(rad_config)

    ! Populate the mapping between the 14 RRTM shortwave bands and the
    ! 6 albedo inputs. The mapping according to the stated wavelength
    ! ranges of the 6-band model does not match the hard-wired mapping
    ! in NMPSRTM, but only the hard-wired values produce sensible
    ! results...
    ! Note that NMPSRTM(:)=(/  6, 6, 5, 5, 5, 5, 5, 4, 4, 3, 2, 2, 1, 6 /)
    rad_config%i_albedo_from_band_sw = NMPSRTM
    !    call rad_config%define_sw_albedo_intervals(6, &
    !         &  (/ 0.25e-6_jprb, 0.44e-6_jprb, 1.19e-6_jprb, &
    !         &     2.38e-6_jprb, 4.00e-6_jprb /),  (/ 1,2,3,4,5,6 /))
    
    ! Likewise between the 16 RRTM longwave bands and the 2 emissivity
    ! inputs (info taken from rrtm_ecrt_140gp_mcica.F90) representing
    ! outside and inside the window region of the spectrum
    ! rad_config%i_emiss_from_band_lw = (/ 1,1,1,1,1,2,2,2,1,1,1,1,1,1,1,1 /)
    call rad_config%define_lw_emiss_intervals(3, &
         &  (/ 8.0e-6_jprb,13.0e-6_jprb /),  (/ 1,2,1 /))

    ! Get spectral weightings for UV and PAR
    call rad_config%get_sw_weights(0.2e-6_jprb, 0.4415e-6_jprb, &
         &  NWEIGHT_UV, IBAND_UV, WEIGHT_UV, 'ultraviolet')
    call rad_config%get_sw_weights(0.4e-6_jprb, 0.7e-6_jprb, &
         &  NWEIGHT_PAR, IBAND_PAR, WEIGHT_PAR, &
         &  'photosynthetically active radiation, PAR')

    IF (YRERAD%NAERMACC > 0) THEN
      ! With the MACC aerosol climatology we need to add in the
      ! background aerosol afterwards using the Tegen arrays.  In this
      ! case we first configure the background aerosol mass-extinction
      ! coefficient at 550 nm, which corresponds to the 10th RRTMG
      ! shortwave band.
      TROP_BG_AER_MASS_EXT  = dry_aerosol_sw_mass_extinction(rad_config, &
           &                                   ITYPE_TROP_BG_AER, 10)
      STRAT_BG_AER_MASS_EXT = dry_aerosol_sw_mass_extinction(rad_config, &
           &                                   ITYPE_STRAT_BG_AER, 10)
      
      WRITE(NULOUT,'(a,i0)') 'Tropospheric bacground uses aerosol type ', &
           &                 ITYPE_TROP_BG_AER
      WRITE(NULOUT,'(a,i0)') 'Stratospheric bacground uses aerosol type ', &
           &                 ITYPE_STRAT_BG_AER
    ENDIF
      
    IF (IVERBOSESETUP > 1) THEN
      WRITE(NULOUT,'(a)') '-------------------------------------------------------------------------------'
    ENDIF

    IF (LHOOK) CALL DR_HOOK('RADIATION_SETUP:SETUP_RADIATION_SCHEME',1,ZHOOK_HANDLE)

  END SUBROUTINE SETUP_RADIATION_SCHEME

END MODULE RADIATION_SETUP
