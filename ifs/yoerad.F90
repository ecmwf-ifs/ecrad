MODULE YOERAD

USE PARKIND1,            ONLY : JPIM, JPRB
USE YOE_SPECTRAL_PLANCK, ONLY : TSPECTRALPLANCK

IMPLICIT NONE

SAVE

!     ------------------------------------------------------------------
!*    ** *YOERAD* - CONTROL OPTIONS FOR RADIATION CONFIGURATION
!     ------------------------------------------------------------------

! Here we have hard-coded the options that configure the radiation
! scheme; in the IFS they are set in suecrad.F90 with the option to
! override them using the NAERAD namelist
TYPE :: TERAD
  INTEGER(KIND=JPIM) :: NSW = 6
  INTEGER(KIND=JPIM) :: NLWEMISS = 2
  INTEGER(KIND=JPIM) :: NICEOPT = 3
  INTEGER(KIND=JPIM) :: NLIQOPT = 4
  INTEGER(KIND=JPIM) :: NRADIP = 3
  INTEGER(KIND=JPIM) :: NRADLP = 2
  INTEGER(KIND=JPIM) :: NLWOUT = 1
  INTEGER(KIND=JPIM) :: NDECOLAT = 2
  INTEGER(KIND=JPIM) :: NMINICE = 1
  INTEGER(KIND=JPIM) :: NAERMACC = 1
  INTEGER(KIND=JPIM) :: NMCVAR = 12
  INTEGER(KIND=JPIM) :: NLWSCATTERING = 1
  INTEGER(KIND=JPIM) :: NSWSOLVER = 3 ! Tripleclouds
  INTEGER(KIND=JPIM) :: NLWSOLVER = 3 ! Tripleclouds
  INTEGER(KIND=JPIM) :: NSOLARSPECTRUM = 0
  INTEGER(KIND=JPIM) :: NDUMPBADINPUTS = 0
  INTEGER(KIND=JPIM) :: NDUMPINPUTS = 0
  INTEGER(KIND=JPIM) :: NCLOUDOVERLAP = 3
  REAL(KIND=JPRB)    :: RCLOUD_FRAC_STD = 1.0_JPRB
  REAL(KIND=JPRB)    :: RCLOUD_SEPARATION_SCALE_TOA = 14000.0_JPRB
  REAL(KIND=JPRB)    :: RCLOUD_SEPARATION_SCALE_SURF = 2500.0_JPRB
  LOGICAL :: LFU_LW_ICE_OPTICS_BUG = .FALSE.
  LOGICAL :: LDIAGFORCING = .FALSE.
  LOGICAL :: LAPPROXLWUPDATE = .TRUE.
  LOGICAL :: LAPPROXSWUPDATE = .FALSE.
  LOGICAL :: LCCNL = .TRUE.
  LOGICAL :: LCCNO = .TRUE.
  REAL(KIND=JPRB) :: RCCNLND = 900.0_JPRB
  REAL(KIND=JPRB) :: RCCNSEA = 50.0_JPRB
  REAL(KIND=JPRB) :: RRE2DE = 0.64952_JPRB
  REAL(KIND=JPRB) :: RMINICE = 60.0_JPRB

  ! Look-up table for Planck function in emissivity intervals
  TYPE(TSPECTRALPLANCK) :: YSPECTPLANCK

END TYPE TERAD
!============================================================================

TYPE(TERAD), POINTER :: YRERAD => NULL()


!        * E.C.M.W.F. PHYSICS PACKAGE *

!     J.-J. MORCRETTE       E.C.M.W.F.      89/07/14
! Modifications
!    R J Hogan 20 May  2014: Added LApproxLwUpdate
!    R J Hogan 19 June 2014: Added LApproxSwUpdate
!    R J Hogan 19 Nov  2014: Added LCentredTimeSZA
!    R J Hogan 15 Apr  2015: Added LMannersSwUpdate
!    R J Hogan 24 Apr  2015: Added LAverageSZA
!    R J Hogan 18 Sept 2015: Added LUsePre2017Rad (was LUsePre2015Rad)
!    R J Hogan 1  Mar  2016: Added NLwSolver, NSwSolver, NLwScattering
!    A Bozzo      Feb  2017: Added logical to enable 3D aerosol climatology
!    R J Hogan 9  Mar  2018: Added NDUMPBADINPUTS
!    R J Hogan 22 Jan  2019: Added NLWEMISS, NCLOUDOVERLAP, NDUMPINPUTS
!    R J Hogan 4  Feb  2019: Added NLWOUT
!    R J Hogan 5  Feb  2019: Added YSPECTPLANCK
!    R J Hogan 11 Mar  2019: Added CGHG*FILE, CSOLARIRRADIANCEFILE

!  NAME     TYPE     PURPOSE
!  ----  :  ----   : ---------------------------------------------------
! LERAD1H: LOGICAL : .T. TO ALLOW MORE FREQUENT RADIATION CALCULATIONS
!                  : DURING FIRST N HOURS OF FORECAST
! NLNGR1H: INTEGER : NUMBER FORECAST HOURS DURING WHICH MORE FREQUENT
!                    RADIATION CALCULATIONS ARE REQUIRED
! LEPO3RA: LOGICAL : .T. IF PROGNOSTIC OZONE (EC) IS PASSED TO RADIATION
! NAER   : INTEGER : CONFIGURATION INDEX FOR AEROSOLS
! NMODE  : INTEGER : CONFIGURATION FOR RADIATION CODE: FLUX VS. RADIANCE
! NOZOCL : INTEGER : CHOICE OF OZONE CLIMATOLOGY (0 old, 1 new)
! NRADFR : INTEGER : FREQUENCY OF FULL RADIATION COMPUTATIONS
!                    IF(NRADFR.GT.0): RAD EVERY 'NRADFR' TIME-STEPS
!                    IF(NRADFR.LT.0): RAD EVERY '-NRADFR' HOURS
! NRADPFR: INTEGER : PRINT FREQUENCY FOR RAD.STATISTICS (in RAD.T.STEPS)
! NRADPLA: INTEGER : PRINT RAD.STATISTICS EVERY 'NRADPLA' ROWS
! NRADINT: INTEGER : RADIATION INTERPOLATION METHOD
!                  : 1 = SPECTRAL TRANSFORM INTERPOLATION
!                  : 2 =  4 POINT HORIZONTAL INTERPOLATION
!                  : 3 = 12 POINT HORIZONTAL INTERPOLATION
! NRADRES: INTEGER : RADIATION GRID SPECTRAL RESOLUTION
! NRADNFR: INTEGER : NORMAL   FREQUENCY OF RADIATION STEPS
! NRADSFR: INTEGER : START-UP FREQUENCY OF RADIATION STEPS
! NRADE1H: INTEGER : START-UP FREQUENCY OF RADIATION STEPS FOR EPS
! NRADE3H: INTEGER : SUBSEQUENT FREQUENCY OF RADIATION STEPS FOR EPS
! NRADELG: INTEGER : LENGTH IN HOURS DURING WHICH THE FREQUENCY OF RADIATION IS INCREASED FOR EPS
! NOVLP  : INTEGER : CLOUD OVERLAP CONFIGURATION IN PRE-MCRAD/ECRAD SCHEME
!                  : 1 = Max-rand (Geleyn & Hollingsworth)
!                  : 2 = Maximum
!                  : 3 = Random
! NRPROMA: INTEGER : VECTOR LENGTH FOR RADIATION CALCULATIONS
! NSW    : INTEGER : NUMBER OF SHORTWAVE SPECTRAL INTERVALS
! NSWNL  : INTEGER : NUMBER OF SHORTWAVE SPECTRAL INTERVALS IN NL MODEL
! NSWTL  : INTEGER : NUMBER OF SHORTWAVE SPECTRAL INTERVALS IN TL MODEL
! NTSW   : INTEGER : MAXIMUM POSSIBLE NUMBER OF SW SPECTRAL INTERVALS 
! NUV    : INTEGER : NUMBER OF UV SPECTRAL INTERVALS FOR THE UV PROCESSOR   
! LOPTRPROMA:LOGICAL: .T. NRPROMA will be optimised
!                   : .F. NRPROMA will not be optimised (forced
!                   :         by negative NRPROMA in namelist)

! NRADIP : INTEGER : INDEX FOR DIAGNOSIS OF ICE CLOUD EFFECTIVE RADIUS
!          0 = fixed at 40 microns
!          1 = Liou & Ou (1994) capped between 40-130 microns
!          2 = Liou & Ou but capped between 30 and 60 microns
!          3 = Sun & Rikus (1999) revised by Sun (2001)
! NRADLP : INTEGER : INDEX FOR DIAGNOSIS OF LIQ. CLOUD EFFECTIVE RADIUS
!          0 = ERA-15 function of pressure
!          1 = 10 microns over land, 13 microns over sea
!          2 = Martin_et_al (1994) in terms of land-sea number conc
!          3 = Linked to prognostic aerosols
! NICEOPT: INTEGER : INDEX FOR ICE CLOUD OPTICAL PROPERTIES
!          0 = SW Ebert-Curry, LW Smith & Shi (1992)
!          1 = SW Ebert-Curry, LW Ebert-Curry (1992)
!          2 = SW & LW Fu-Liou (1993)
!          3 = SW Fu (1996) LW Fu et al. (1998) + Chou et al. (1999) LW scatt approx
!   the following only available in newer modular radiation scheme:
!          4 = SW/LW Baran data fitted versus ice mixing ratio
! NLIQOPT: INTEGER : INDEX FOR LIQUID WATER CLOUD OPTICAL PROPERTIES
!          0 = SW Fouquart (1991) LW Smith-Shi (1992) YF/SmSh 
!          1 = SW Slingo (1989) LW Savijarvi (1997)
!          2 = SW Slingo (1989) LW Lindner-Li (2000)
!   the following only available in RADLSW, not RADLSWR:
!          3 = SW Nielsen       LW Smith-Shi
!   the following only available in newer modular radiation scheme:
!          4 = SW/LW SOCRATES scheme
!
! LONEWSW: LOGICAL : .T. IF NEW SW CODE IS ACTIVE
! LECSRAD: LOGICAL : .T. IF CLEAR-SKY RADIATION IS ARCHIVED AS PEXTR2
! NCSRADF: INTEGER : 1 IF ACCUMULATED, 2 IF INSTANTANEOUS
! LRRTM  : LOGICAL : .T. IF RRTM140MR IS USED FOR LW RADIATION TRANSFER

! LHVOLCA: LOGICAL : .T. IF USING HISTORICAL VOLCANIC AEROSOLS 
! LNEWAER: LOGICAL : .T. IF AEROSOL MONTHLY DISTRIBUTIONS ARE USED
! LNOTROAER:LOGICAL: .T. IF NO TROPOSPHERIC AEROSOLS
! CRTABLEDIR: CHAR : IF NRADINT > 0 SPECIFIES DIRECTORY PATH FOR RADIATION
!                  : GRID RTABLE NAMELIST
! CRTABLEFIL: CHAR : IF NRADINT > 0 SPECIFIES FILE NAME OF RADIATION 
!                  : GRID RTABLE NAMELIST
! LRAYL  : LOGICAL : .T. NEW RAYLEIGH FOR SW-6 VERSION

! RAOVLP : REAL    : COEFFICIENTS FOR ALPHA1 FACTOR IN HOGAN & 
! RBOVLP : REAL    : ILLINGWORTH's PARAMETRIZATION

! LCCNL  : LOGICAL : .T. IF CCN CONCENTRATION OVER LAND IS DIAGNOSED
! LCCNO  : LOGICAL : .T. IF CCN CONCENTRATION OVER OCEAN IS DIAGNOSED
! RCCNLND: REAL    : NUMBER CONCENTRATION (CM-3) OF CCNs OVER LAND
! RCCNSEA: REAL    : NUMBER CONCENTRATION (CM-3) OF CCNs OVER SEA

! LDIFFC : LOGICAL : .T. IF SAVIJARVI'S DIFFUSIVITY CORRECTION IS ON

! NINHOM : INTEGER : 0 IF NO INHOMOGENEITY SCALING EFFECT 
!                    1 IF SIMPLE 0.7 SCALING
!                    2 IF BARKER, 3 IF CAIRNS ET AL.
! RLWINHF: REAL    : INHOMOG. SCALING FACTOR FOR CLOUD LW OPTICAL THICKNESS
! RSWINHF: REAL    : INHOMOG. SCALING FACTOR FOR CLOUD SW OPTICAL THICKNESS

! NPERTAER : INTERGER : PERCENTAGE OF PERTURBATION FOR AEROSOL   
! NPERTOZONE : INTEGER : PERCENTAGE OF PERTURBATION FOR OZONE 
! NHINCSOL : INTEGER :  0: Total Solar Irradiance (TSI) fixed at 1366.0 W m-2
!                       1: Deprecated - use default
!                       2: Deprecated - use default
!                       3: Deprecated (was CMIP5) - use default
!                       4: TSI from CMIP6 NetCDF (default), or override with CSOLARIRRADIANCEFILE
! NSWWVCONTINUUM : INTEGER : 0 MT_CKD2.5 (SRTM default WV continuum)
!                            1 CAVIAR continuum (Shine et al. 2016)
! LECO2VAR: LOGICAL: .T. IF ERA-40/AMIP2 VARIABILITY OF GHG IS ON (ignored)
! LHGHG  : LOGICAL : .T. IF VARIABILITY OF GREENHOUSE GASES (INCLUDING CO2) IS ON
! N.B.: LHGHG supersedes LECO2VAR and allows using better specification of trace gases
! NSCEN  : INTEGER : 21st CENTURY SCENARIO FOR GHG (1=A1B, 2=A2, 3=B1)
! RRe2De : REAL    : CONVERSION FACTOR BETWWEN EFFECTIVE RADIUS AND PARTICLE SIZE
! RMINICE: REAL    : MINIMUM SIZE FOR ICE PARTICLES (um)
!                    FOR ICE
! NMINICE: INTEGER : 1-6 MINIMUM ICE PARTICLE SIZE DEPENDS ON LATITUDE, 0=INDEPENDENT OF LATITUDE
! NDECOLAT:INTEGER : DECORRELATION LENGTH FOR CF AND CW 
!                     0: SPECIFIED INDEPENDENT OF LATITUDE, 1: SHONK-HOGAN, 2: IMPROVED
! NMCICA : INTEGER :  0: NO McICA
!                     1: McICA w maximum-random in cloud generator
!                     2: McICA w generalized overlap in cloud generator
! LESO4HIS: LOGICAL:.T.: Use historical/projected SO4 data per decade and month
! NGHGRAD: INTEGER : configuration of 3D GHG climatologies accounted for in radiation
!                     0: global values
!                     1: CO2       2: CH4    3: N2O    4: NO2    5:CFC11   6:CFC12
!                    12: CO2+CH4  13: CO2+CH4+N2O     
!                    16: CO2+CH4+N2O+CFC11+CFC12
! LETRACGMS: LOGICAL : F=Cariolle climatol. T=GEMS-derived clim for CO2, CH4, O3
! LAERCLIM : LOGICAL : .T. for output of the climatological aerosol optical depth at 550 nm
! LAERVISI : LOGICAL : .T. for output of the visibility (from diagnsotic or prognostic aerosols)
! NVOLCVERT: INTEGER : Vertical distribution of volcanic aerosol
!                       0: original profile, diagnosed from T
!                       1: original profile, but upper boundary at 10hPa
!                       2: lower boundary diagnosed from ozone, upper boundary at 10hPa
! LVOLCSPEC: LOGICAL : T for specified volcanic aerosol
! LVOLCDAMP: LOGICAL : T for damping of specified volcanic aerosol from initial value
! RVOLCSPEC: REAL    : Specified volcanic aerosol (total optical depth) in NH/Tropics/SH
! RNs                : derived from Avogadro
! RSIGAIR: invariant terms in expression of Rayleigh scattering cross-section
! NREDGSW  : INTEGER : 0 full resolution for RRTM_SW (224)
!                      1 ECMWF High resolution model configuration (_SW: 112)
!                      2 ECMWF EPS configuration (_SW: 56)
! NREDGLW  : INTEGER : 0 full resolution for RRTM_LW (256)
!                      1 ECMWF High resolution model configuration (_LW: 140)
!                      2 ECMWF EPS configuration (_LW: 70)
! LDIAGFORCING : LOGICAL : T Write input ozone, ghg and aerosol forcing to 3D fields 
!                            To be used for diagnostics only; do not use in production runs
! NAERMACC : INTEGER : MACC-derived aerosol climatology on a NMCLAT x NMCLON grid
! RAESHxx  : REAL    : parameters related to scale height of MACC-derived aerosol climatology
! CVDAExx  : REAL    : scale heights of MACC-derived aerosol climatology
! LAERADJDU: LOGICAL : T adjust MACC-derived DU climatology
! LAERADCLI: LOGICAL : T if radiation uses the MACC-derived aerosol climatology
! LApproxLwUpdate : LOGICAL : Update the longwave upwelling flux every
!                             timestep/gridpoint using the stored rate
!                             of change of the fluxes with respect to
!                             the surface upwelling longwave flux
! LApproxSwUpdate : LOGICAL : Update the shortwave upwelling flux
!                             every gridpoint to account for the local
!                             value of surface albedo
! LMannersSwUpdate: LOGICAL : Update the shortwave flux every timestep
!                             using Manners et al. (2009) correction
!                             for solar zenith angle change
! LCentredTimeSZA : LOGICAL : Compute solar zenith angle in radiation
!                             scheme half way between calls to
!                             radiation scheme (rather than previous
!                             behaviour, which is half way between
!                             calls plus half a model timestep)
! LAverageSZA     : LOGICAL : Compute an averaged solar zenith angle
!                             across the time interval required
!                             (either a model timestep or a radiation
!                             timestep). Should be used with 
!                             LCentredTimeSZA=TRUE.
! LUsePre2017Rad  : LOGICAL : Use the pre-2017 radiation scheme, rather 
!                             than the modular scheme contained in the
!                             separate "radiation" library.  Note that
!                             the radiation library may make use of the
!                             pre-2017 RRTM-G gas optics.
! RCloud_Frac_Std : REAL    : Cloud water content horizontal fractional
!                             standard deviation in a gridbox
! LInterpInCloudMean : LOGICAL : When interpolating model fields to
!                             radiation grid, interpolate in-cloud
!                             mean water contents?  Better
!                             conservation achieved by interpolating
!                             gridbox-means.
! CGHGCLIMFILE : STRING :     Location of greenhouse gas climatology file,
!                             or empty if the default is to be used. If it
!                             starts with "." or "/" then a relative path
!                             is assumed, otherwise the default directory.
! CGHGTIMESERIESFILE:STRING : Location of greenhouse gas timeseries
!                             file, or empty if it is to be worked out
!                             from the NSCEN, YOECMIP%NGHGCMIP and
!                             YOECMIP%NRCP variables. If it starts
!                             with "." or "/" then a relative path is
!                             assumed, otherwise the default
!                             directory.
! CSOLARIRRADIANCEFILE:STRING:Location of Total Solar Irradiance file,
!                             or empty if the default is to be
!                             used. If it starts with "." or "/" then
!                             a relative path is assumed, otherwise
!                             the default directory.
! NLWEMISS      : INTEGER :   Number of emissivity spectral intervals, set 
!                             according to the value of NEMISSSCHEME; traditionally
!                             this has always been 2: outside the IR window and within
! NLWOUT        : INTEGER :   Number of spectral intervals to pass LW downwelling flux
!                             to RADHEATN; traditionally this was 1, but this led
!                             to errors with LAPPROXLWUPDATE=TRUE, which updated 
!                             fluxes using a single broadband emissivity. Now we can
!                             do approximate updates using full spectral emissivity.
! ------------------------------------------------------------------
! THE FOLLOWING ARE ONLY USED FOR THE ECRAD SCHEME (LUsePre2017Rad = .FALSE.)
! NLwScattering   : INTEGER : 0: No longwave scattering
!                             1: Longwave scattering by clouds only
!                             2: Longwave scattering by clouds and aerosols
! NSwSolver       : INTEGER :
! NLwSolver       : INTEGER : 0: McICA
!                             1: SPARTACUS 1D
!                             2: SPARTACUS 3D
!                             3: TripleClouds
! LFU_LW_ICE_OPTICS_BUG : LOGICAL : Continue to use bug in Fu LW ice
!                             optics whereby single scattering albedo is
!                             one minus what it should be
! NSOLARSPECTRUM : INTEGER :  0: Kurucz
!                             1: Coddington et al. (BAMS 2016)
! NDUMPBADINPUTS : INTEGER :  0: Warn only if fluxes out of physical bounds
!                             n: Write netcdf file of bad inputs up to n times per task
!                            -n: Abort if fluxes ever out of physical bounds
! NDUMPINPUTS    : INTEGER :  0: Do nothing
!                             n: Write netcdf file of all inputs up to n times per task
! NCLOUDOVERLAP  : INTEGER :  Cloud overlap scheme
!                             1: Maximum-random
!                             2: Exponential-exponential (the actual behaviour of McRad)
!                             3: Exponential-random (only option for Tripleclouds and SPARTACUS)
! RCLOUD_SEPARATION_SCALE_TOA, RCLOUD_SEPARATION_SCALE_SURF : REAL
!                 Cloud horizontal length scale, in metres, used to
!                 compute rate of horizontal exchange of radiation
!                 between clouds and clear skies in SPARTACUS solver
! ------------------------------------------------------------------
! KMODTS : INTEGER   : (A Bozzo) switch for different radiative transfer schemes for UV 
!                       = 0 Fouquart&Bonnel adapted by Morcrette and Arola
!                       = 1 eddington (joseph et al., 1976)
!                       = 2 pifm (zdunkowski et al., 1980)
!                       = 3 discrete ordinates (liou, 1973)
!     ------------------------------------------------------------------
! TRBKG : REAL tropospheric background OD@550nm for aerosol climatology.
!                  default for Tegen climatology was 0.03
! STBKG : REAL stratospheric background OD@550nm for aerosol climatology.
!     ------------------------------------------------------------------
! LDUSEASON : LOGICAL enables a monthly-varying scale height for the 
!                     dust aerosol climatology
! LAER3D : LOGICAL : to enable aerosol climatology in 3D


END MODULE YOERAD
