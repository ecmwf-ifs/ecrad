SUBROUTINE RADIATION_SCHEME &
     & (YRADIATION, KIDIA, KFDIA, KLON, KLEV, KAEROSOL, &
     &  PSOLAR_IRRADIANCE, &
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
     &  PRE_LIQ, PRE_ICE, ISEED, PCLOUD_OVERLAP)

! RADIATION_SCHEME - Interface to modular radiation scheme
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
!   library. This routine puts the the IFS arrays into appropriate
!   objects, computing the additional data that is required, and sends
!   it to the radiation scheme.  It returns net fluxes and surface
!   flux components needed by the rest of the model.
!
!   Lower case is used for variables and types taken from the
!   radiation library
!
! INTERFACE
! ---------
!    RADIATION_SCHEME is called from RADLSWR. The
!    SETUP_RADIATION_SCHEME routine (in the RADIATION_SETUP module)
!    populates the YRADIATION object, and should have been run first.
!
! AUTHOR
! ------
!   Robin Hogan, ECMWF
!   Original: 2015-09-16
!
! MODIFICATIONS
! -------------
!   2017-03-03  R. Hogan  Read configuration data from YRADIATION object
!   2017-05-11  R. Hogan  Pass KIDIA,KFDIA to get_layer_mass
!   2018-01-11  R. Hogan  Capability to scale solar spectrum in each band
!   2017-11-11  M. Ahlgrimm add variable FSD for cloud heterogeneity
!   2017-11-29  R. Hogan  Check fluxes in physical bounds
!   2019-01-22  R. Hogan  Use fluxes in albedo bands from ecRad
!   2019-01-23  R. Hogan  Spectral longwave emissivity in NLWEMISS bands
!   2019-02-04  R. Hogan  Pass out surface longwave downwelling in each emissivity interval
!   2019-02-07  R. Hogan  SPARTACUS cloud size from PARAM_CLOUD_EFFECTIVE_SEPARATION_ETA
!   2020-10-12  M. Leutbecher SPP abstraction
!   2021-08-26  R. Hogan  Added "true" Coddington solar spectrum
!   2022-07-08  R. El Khatib Contribution to the encapsulation of YOMCST and YOETHF
!   2023-01, M. Michou : modifications so that chemistry interacts with ecRad
!   2023-09, P. Nabat : add AEROCPL for TACTIC aerosol-radiation coupling
!
!-----------------------------------------------------------------------

! Modules from ifs or ifsaux libraries
USE PARKIND1       , ONLY : JPIM, JPRB, JPRD
USE YOMHOOK        , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOMCST         , ONLY : RPI, RSIGMA
USE YOMLUN_ECRAD   , ONLY : NULERR, NULOUT
USE RADIATION_SETUP, ONLY : ITYPE_TROP_BG_AER, ITYPE_STRAT_BG_AER, TRADIATION

! Modules from ecRad radiation library
USE RADIATION_CONFIG,         ONLY : ISOLVERSPARTACUS, CONFIG_TYPE
USE RADIATION_SINGLE_LEVEL,   ONLY : SINGLE_LEVEL_TYPE
USE RADIATION_THERMODYNAMICS, ONLY : THERMODYNAMICS_TYPE
USE RADIATION_GAS,            ONLY : GAS_TYPE,&
     &                               IMASSMIXINGRATIO, IVOLUMEMIXINGRATIO,&
     &                               IH2O, ICO2, ICH4, IN2O, ICFC11, ICFC12, IHCFC22, ICCL4, IO3, IO2
USE RADIATION_CLOUD,          ONLY : CLOUD_TYPE
USE RADIATION_AEROSOL,        ONLY : AEROSOL_TYPE
USE RADIATION_FLUX,           ONLY : FLUX_TYPE
USE RADIATION_INTERFACE,      ONLY : RADIATION, SET_GAS_UNITS
USE RADIATION_SAVE,           ONLY : SAVE_INPUTS, SAVE_FLUXES
#ifdef HAVE_NVTX
USE NVTX
#endif

IMPLICIT NONE

! INPUT ARGUMENTS

TYPE(TRADIATION), INTENT(IN),TARGET   :: YRADIATION

! *** Array dimensions and ranges
INTEGER(KIND=JPIM),INTENT(IN)   :: KIDIA    ! Start column to process
INTEGER(KIND=JPIM),INTENT(IN)   :: KFDIA    ! End column to process
INTEGER(KIND=JPIM),INTENT(IN)   :: KLON     ! Number of columns
INTEGER(KIND=JPIM),INTENT(IN)   :: KLEV     ! Number of levels
INTEGER(KIND=JPIM),INTENT(IN)   :: KAEROSOL ! Number of aerosol types

! *** Single-level fields
REAL(KIND=JPRB),   INTENT(IN) :: PSOLAR_IRRADIANCE ! (W m-2)
REAL(KIND=JPRB),   INTENT(IN) :: PMU0(KLON) ! Cosine of solar zenith ang
REAL(KIND=JPRB),   INTENT(IN) :: PTEMPERATURE_SKIN(KLON) ! (K)
! Diffuse and direct components of surface shortwave albedo
REAL(KIND=JPRB),   INTENT(IN) :: PALBEDO_DIF(KLON,YRADIATION%YRERAD%NSW)
REAL(KIND=JPRB),   INTENT(IN) :: PALBEDO_DIR(KLON,YRADIATION%YRERAD%NSW)
! Longwave spectral emissivity
REAL(KIND=JPRB),   INTENT(IN) :: PSPECTRALEMISS(KLON,YRADIATION%YRERAD%NLWEMISS)
! Longitude (radians), sine of latitude
REAL(KIND=JPRB),   INTENT(IN) :: PGELAM(KLON)
REAL(KIND=JPRB),   INTENT(IN) :: PGEMU(KLON)
! Land-sea mask
REAL(KIND=JPRB),   INTENT(IN) :: PLAND_SEA_MASK(KLON)

! *** Variables on full levels
REAL(KIND=JPRB),   INTENT(IN) :: PPRESSURE(KLON,KLEV)    ! (Pa)
REAL(KIND=JPRB),   INTENT(IN) :: PTEMPERATURE(KLON,KLEV) ! (K)
! *** Variables on half levels
REAL(KIND=JPRB),   INTENT(IN) :: PPRESSURE_H(KLON,KLEV+1)    ! (Pa)
REAL(KIND=JPRB),   INTENT(IN) :: PTEMPERATURE_H(KLON,KLEV+1) ! (K)

! *** Gas mass mixing ratios on full levels
REAL(KIND=JPRB),   INTENT(IN) :: PQ(KLON,KLEV)
REAL(KIND=JPRB),   INTENT(IN) :: PCO2(KLON,KLEV)
REAL(KIND=JPRB),   INTENT(IN) :: PCH4(KLON,KLEV)
REAL(KIND=JPRB),   INTENT(IN) :: PN2O(KLON,KLEV)
REAL(KIND=JPRB),   INTENT(IN) :: PNO2(KLON,KLEV)
REAL(KIND=JPRB),   INTENT(IN) :: PCFC11(KLON,KLEV)
REAL(KIND=JPRB),   INTENT(IN) :: PCFC12(KLON,KLEV)
REAL(KIND=JPRB),   INTENT(IN) :: PHCFC22(KLON,KLEV)
REAL(KIND=JPRB),   INTENT(IN) :: PCCL4(KLON,KLEV)
REAL(KIND=JPRB),   INTENT(IN) :: PO3(KLON,KLEV)

! *** Cloud fraction and hydrometeor mass mixing ratios
REAL(KIND=JPRB),   INTENT(IN) :: PCLOUD_FRAC(KLON,KLEV)
REAL(KIND=JPRB),   INTENT(IN) :: PQ_LIQUID(KLON,KLEV)
REAL(KIND=JPRB),   INTENT(IN) :: PQ_ICE(KLON,KLEV)
REAL(KIND=JPRB),   INTENT(IN) :: PQ_RAIN(KLON,KLEV)
REAL(KIND=JPRB),   INTENT(IN) :: PQ_SNOW(KLON,KLEV)

! *** Aerosol mass mixing ratios
REAL(KIND=JPRB),   INTENT(IN) :: PAEROSOL_OLD(KLON,6,KLEV)
REAL(KIND=JPRB),   INTENT(IN) :: PAEROSOL(KLON,KLEV,KAEROSOL)

REAL(KIND=JPRB),   INTENT(IN) :: PCCN_LAND(KLON)
REAL(KIND=JPRB),   INTENT(IN) :: PCCN_SEA(KLON)

! OUTPUT ARGUMENTS

! *** Net fluxes on half-levels (W m-2)
REAL(KIND=JPRB),  INTENT(INOUT) :: PFLUX_SW(KLON,KLEV+1)
REAL(KIND=JPRB),  INTENT(INOUT) :: PFLUX_LW(KLON,KLEV+1)
REAL(KIND=JPRB),  INTENT(INOUT) :: PFLUX_SW_CLEAR(KLON,KLEV+1)
REAL(KIND=JPRB),  INTENT(INOUT) :: PFLUX_LW_CLEAR(KLON,KLEV+1)

! *** Surface flux components (W m-2)
REAL(KIND=JPRB),  INTENT(INOUT) :: PFLUX_SW_DN(KLON)
REAL(KIND=JPRB),  INTENT(INOUT) :: PFLUX_LW_DN(KLON)
REAL(KIND=JPRB),  INTENT(INOUT) :: PFLUX_SW_DN_CLEAR(KLON)
REAL(KIND=JPRB),  INTENT(INOUT) :: PFLUX_LW_DN_CLEAR(KLON)
! Direct component of surface flux into horizontal plane
REAL(KIND=JPRB),  INTENT(INOUT) :: PFLUX_DIR(KLON)
REAL(KIND=JPRB),  INTENT(INOUT) :: PFLUX_DIR_CLEAR(KLON)
! As PFLUX_DIR but into a plane perpendicular to the sun
REAL(KIND=JPRB),  INTENT(INOUT) :: PFLUX_DIR_INTO_SUN(KLON)

! *** Ultraviolet and photosynthetically active radiation (W m-2)
REAL(KIND=JPRB),  INTENT(INOUT) :: PFLUX_UV(KLON)
REAL(KIND=JPRB),  INTENT(INOUT) :: PFLUX_PAR(KLON)
REAL(KIND=JPRB),  INTENT(INOUT) :: PFLUX_PAR_CLEAR(KLON)

! *** Other single-level diagnostics
! Top-of-atmosphere incident solar flux (W m-2)
REAL(KIND=JPRB),  INTENT(INOUT) :: PFLUX_SW_DN_TOA(KLON)
! Diagnosed longwave surface emissivity across the whole spectrum
REAL(KIND=JPRB),  INTENT(INOUT) :: PEMIS_OUT(KLON)

! Partial derivative of total-sky longwave upward flux at each level
! with respect to upward flux at surface, used to correct heating
! rates at gridpoints/timesteps between calls to the full radiation
! scheme.  Note that this version uses the convention of level index
! increasing downwards, unlike the local variable ZLwDerivative that
! is returned from the LW radiation scheme.
REAL(KIND=JPRB),  INTENT(INOUT) :: PLWDERIVATIVE(KLON,KLEV+1)

! Surface diffuse and direct downwelling shortwave flux in each
! shortwave albedo band, used in RADINTG to update the surface fluxes
! accounting for high-resolution albedo information
REAL(KIND=JPRB),  INTENT(INOUT) :: PSWDIFFUSEBAND(KLON,YRADIATION%YRERAD%NSW)
REAL(KIND=JPRB),  INTENT(INOUT) :: PSWDIRECTBAND (KLON,YRADIATION%YRERAD%NSW)

! Optional input arguments (Added for validating against ecrad standalone!)
REAL(KIND=JPRB), INTENT(IN), OPTIONAL :: PRE_LIQ(KLON, KLEV)
REAL(KIND=JPRB), INTENT(IN), OPTIONAL :: PRE_ICE(KLON, KLEV)
INTEGER,         INTENT(IN), OPTIONAL :: ISEED(KLON)
REAL(KIND=JPRB), INTENT(IN), OPTIONAL :: PCLOUD_OVERLAP(KLON, KLEV-1)

! Enable GPU code paths
LOGICAL                                :: LLACC

! LOCAL VARIABLES
TYPE(SINGLE_LEVEL_TYPE)   :: SINGLE_LEVEL
TYPE(THERMODYNAMICS_TYPE) :: THERMODYNAMICS
TYPE(GAS_TYPE)            :: GAS
TYPE(CLOUD_TYPE)          :: YLCLOUD
TYPE(AEROSOL_TYPE)        :: AEROSOL
TYPE(FLUX_TYPE)           :: FLUX

REAL(KIND=JPRB)           :: ZSIGMA

! Cloud effective radii in microns
REAL(KIND=JPRB)           :: ZRE_LIQUID_UM(KLON,KLEV)
REAL(KIND=JPRB)           :: ZRE_ICE_UM(KLON,KLEV)

! Cloud overlap decorrelation length for cloud boundaries in km
REAL(KIND=JPRB)           :: ZDECORR_LEN_KM(KLON)

! Ratio of cloud overlap decorrelation length for cloud water
! inhomogeneities to that for cloud boundaries (typically 0.5)
REAL(KIND=JPRB)           :: ZDECORR_LEN_RATIO

! The surface net longwave flux if the surface was a black body, used
! to compute the effective broadband surface emissivity
REAL(KIND=JPRB)           :: ZBLACK_BODY_NET_LW

! Layer mass in kg m-2
REAL(KIND=JPRB)           :: ZLAYER_MASS(KIDIA:KFDIA,KLEV)

! Time integers
! INTEGER(KIND=JPIM) :: ITIM, IDAY

! Loop indices
INTEGER(KIND=JPIM) :: JLON, JLEV, JBAND, JAER

! Have any fluxes been returned that are out of a physically
! reasonable range? This integer stores the number of blocks of fluxes
! that have contained a bad value so far, for this task.  NetCDF files
! will be written up to the value of NAERAD:NDUMPBADINPUTS.
INTEGER(KIND=JPIM), SAVE :: N_BAD_FLUXES = 0

! For debugging it can be useful to save input profiles and output
! fluxes without the condition that the fluxes are out of a reasonable
! range. NetCDF files will be written up to the value of
! NAERAD:NDUMPINPUTS.
INTEGER(KIND=JPIM), SAVE :: N_OUTPUT_FLUXES = 0

! NetCDF file name in case of bad fluxes
CHARACTER(LEN=512) :: CL_FILE_NAME

! Old Tegen aerosol scheme diagnosed by other logicals being false:
! the logic setting this needs to be the same as in
! su_radiation_scheme.F90
LOGICAL :: LL_USE_TEGEN_AEROSOLS

TYPE (CONFIG_TYPE), POINTER :: RAD_CONFIG

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE


! Dummy from MPL_MYRANK_MOD
INTEGER(KIND=JPIM) :: MPL_MYRANK
MPL_MYRANK() = 1

! Import time functions for iseed calculation
#include "fcttim.func.h"

#include "liquid_effective_radius.intfb.h"
#include "ice_effective_radius.intfb.h"
#include "cloud_overlap_decorr_len.intfb.h"
!#include "satur.intfb.h"
!#include "abor1.intfb.h"

IF (LHOOK) CALL DR_HOOK('RADIATION_SCHEME',0,ZHOOK_HANDLE)

! RAD_CONFIG is modified in this routine
! - either use a private object
! - or move the modification outside this routine
! For now use pointer laundering

RAD_CONFIG => YRADIATION%RAD_CONFIG

ASSOCIATE(YRERAD    =>YRADIATION%YRERAD, &
     &    NWEIGHT_UV=>YRADIATION%NWEIGHT_UV, &
     &    IBAND_UV  =>YRADIATION%IBAND_UV(:), &
     &    WEIGHT_UV =>YRADIATION%WEIGHT_UV(:), &
     &    NWEIGHT_PAR=>YRADIATION%NWEIGHT_PAR, &
     &    IBAND_PAR =>YRADIATION%IBAND_PAR(:), &
     &    WEIGHT_PAR=>YRADIATION%WEIGHT_PAR(:), &
     &    TROP_BG_AER_MASS_EXT=>YRADIATION%TROP_BG_AER_MASS_EXT, &
     &    STRAT_BG_AER_MASS_EXT=>YRADIATION%STRAT_BG_AER_MASS_EXT)
! ASSOCIATE(YDRADIATION=>YRADIATION%YRADIATION, &
!      &    YRERAD=>YRADIATION%YRERAD, &
!      &    YDERDI=>YRADIATION%YRERDI)
! ASSOCIATE(NWEIGHT_UV=>YDRADIATION%NWEIGHT_UV, &
!      &    IBAND_UV  =>YDRADIATION%IBAND_UV(:), &
!      &    WEIGHT_UV =>YDRADIATION%WEIGHT_UV(:), &
!      &    NWEIGHT_PAR=>YDRADIATION%NWEIGHT_PAR, &
!      &    IBAND_PAR =>YDRADIATION%IBAND_PAR(:), &
!      &    WEIGHT_PAR=>YDRADIATION%WEIGHT_PAR(:), &
!      &    TROP_BG_AER_MASS_EXT=>YDRADIATION%TROP_BG_AER_MASS_EXT, &
!      &    STRAT_BG_AER_MASS_EXT=>YDRADIATION%STRAT_BG_AER_MASS_EXT, &
!      &    LAERRADCPL=>YDEAERATM%LAERRADCPL, &
!      &    NAERMACC=>YRERAD%NAERMACC, &
!      &    NOZOCL=>YRERAD%NOZOCL,&
!      &    NTSW=>YRERAD%NTSW, NTLW=>YRERAD%NTLW)

! Create local copy of RSIGMA to ensure value is offloaded to GPU
ZSIGMA = RSIGMA

! Allocate memory in radiation objects
#ifdef HAVE_NVTX
call nvtxStartRange("allocate")
#endif
CALL SINGLE_LEVEL%ALLOCATE(KLON, YRERAD%NSW, YRERAD%NLWEMISS, &
     &                     USE_SW_ALBEDO_DIRECT=.TRUE.)
CALL THERMODYNAMICS%ALLOCATE(KLON, KLEV, USE_H2O_SAT=.TRUE.)
CALL GAS%ALLOCATE(KLON, KLEV)
CALL YLCLOUD%ALLOCATE(KLON, KLEV)
IF (YRERAD%NAERMACC == 1) THEN
  LL_USE_TEGEN_AEROSOLS = .FALSE.
  CALL AEROSOL%ALLOCATE(KLON, 1, KLEV, KAEROSOL) ! Prognostic or CAMS/MACC aerosols
ELSE
  LL_USE_TEGEN_AEROSOLS = .TRUE.
  WRITE(NULOUT,'(A)') 'RADIATION_SCHEME: Warning, assuming Tegen aerosols'
  CALL AEROSOL%ALLOCATE(KLON, 1, KLEV, 6) ! Tegen climatology
ENDIF
CALL FLUX%ALLOCATE(RAD_CONFIG, 1, KLON, KLEV)
#ifdef HAVE_NVTX
call nvtxEndRange
#endif

#if defined(_OPENACC) || defined(OMPGPU)
LLACC = .TRUE.
#else
LLACC = .FALSE.
#endif

!$ACC DATA &
!$ACC COPYIN(YRADIATION, YRERAD, RAD_CONFIG, SINGLE_LEVEL, THERMODYNAMICS, GAS, AEROSOL, YLCLOUD, FLUX) &
!$ACC ASYNC(1) IF(LLACC)
!$OMP TARGET ENTER DATA MAP(TO:YRADIATION, YRERAD, RAD_CONFIG, SINGLE_LEVEL, THERMODYNAMICS, GAS, AEROSOL, YLCLOUD, FLUX) IF(LLACC)
IF (LLACC) THEN
  CALL YRADIATION%RAD_CONFIG%CREATE_DEVICE(YRADIATION%RAD_CONFIG)
  ! CALL RAD_CONFIG%UPDATE_DEVICE()
  CALL SINGLE_LEVEL%CREATE_DEVICE(SINGLE_LEVEL)
  CALL THERMODYNAMICS%CREATE_DEVICE(THERMODYNAMICS)
  CALL GAS%CREATE_DEVICE(GAS)
  CALL AEROSOL%CREATE_DEVICE(AEROSOL)
  CALL YLCLOUD%CREATE_DEVICE(YLCLOUD)
  CALL FLUX%CREATE_DEVICE(FLUX)
ENDIF

!$ACC DATA &
!$ACC CREATE(ZRE_LIQUID_UM, ZRE_ICE_UM, ZDECORR_LEN_KM, ZLAYER_MASS) &
!$ACC PRESENT(PMU0, PTEMPERATURE_SKIN, PALBEDO_DIF, PALBEDO_DIR, PSPECTRALEMISS, &
!$ACC         PCCN_LAND, PCCN_SEA, PGEMU, PLAND_SEA_MASK, PPRESSURE, PTEMPERATURE, &
!$ACC         PPRESSURE_H, PTEMPERATURE_H, &
!$ACC         PCLOUD_FRAC, PQ_LIQUID, PQ_ICE, PQ_RAIN, PQ_SNOW, &
!$ACC         PAEROSOL_OLD, PAEROSOL, &
!$ACC         PFLUX_SW, PFLUX_LW, PFLUX_SW_CLEAR, PFLUX_LW_CLEAR, &
!$ACC         PFLUX_SW_DN, PFLUX_LW_DN, PFLUX_SW_DN_CLEAR, PFLUX_LW_DN_CLEAR, &
!$ACC         PFLUX_DIR, PFLUX_DIR_CLEAR, PFLUX_DIR_INTO_SUN, &
!$ACC         PFLUX_UV, PFLUX_PAR, PFLUX_PAR_CLEAR, PFLUX_SW_DN_TOA, &
!$ACC         PEMIS_OUT, PLWDERIVATIVE) &
!$ACC NO_CREATE(PSWDIRECTBAND, PSWDIFFUSEBAND, PRE_LIQ) ASYNC(1) IF(LLACC)

!$OMP TARGET ENTER DATA MAP(ALLOC:ZRE_LIQUID_UM, ZRE_ICE_UM, ZDECORR_LEN_KM,ZLAYER_MASS) IF(LLACC)

#ifdef DEBUG_WARNING
write(nulout,'(a,a,a,i0,a)') "    ", __FILE__, " : LINE = ", __LINE__, " WARNING : Checking the PRESENT status of subarrays does not work. This needs to be fixed! PJM 4-9-2025"
#endif
!!$OMP TARGET DATA MAP(PRESENT,ALLOC: &
!!$OMP&         PMU0, PTEMPERATURE_SKIN, PALBEDO_DIF, PALBEDO_DIR, PSPECTRALEMISS, &
!!$OMP&         PCCN_LAND, PCCN_SEA, PGEMU, PLAND_SEA_MASK, PPRESSURE, PTEMPERATURE, &
!!$OMP&         PPRESSURE_H, PTEMPERATURE_H, &
!!$OMP&         PCLOUD_FRAC, PQ_LIQUID, PQ_ICE, PQ_RAIN, PQ_SNOW, &
!!$OMP&         PAEROSOL_OLD, PAEROSOL, &
!!$OMP&         PFLUX_SW, PFLUX_LW, PFLUX_SW_CLEAR, PFLUX_LW_CLEAR, &
!!$OMP&         PFLUX_SW_DN, PFLUX_LW_DN, PFLUX_SW_DN_CLEAR, PFLUX_LW_DN_CLEAR, &
!!$OMP&         PFLUX_DIR, PFLUX_DIR_CLEAR, PFLUX_DIR_INTO_SUN, &
!!$OMP&         PFLUX_UV, PFLUX_PAR, PFLUX_PAR_CLEAR, PFLUX_SW_DN_TOA, &
!!$OMP&         PEMIS_OUT, PLWDERIVATIVE)
!!$OMP TARGET DATA MAP(PRESENT,ALLOC:PSWDIRECTBAND, PSWDIFFUSEBAND, PRE_LIQ)

#ifdef HAVE_NVTX
call nvtxStartRange("thermodynamics setup")
#endif
! Set thermodynamic profiles: simply copy over the half-level
! pressure and temperature
!$ACC PARALLEL DEFAULT(NONE) ASYNC(1) IF(LLACC) &
!$ACC   PRESENT(THERMODYNAMICS, THERMODYNAMICS%PRESSURE_HL, THERMODYNAMICS%TEMPERATURE_HL)
!$ACC LOOP GANG VECTOR COLLAPSE (2)
!$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO COLLAPSE(2) IF(LLACC)
DO JLEV = 1,KLEV+1
  DO JLON = KIDIA,KFDIA
    THERMODYNAMICS%PRESSURE_HL   (JLON,JLEV) = PPRESSURE_H   (JLON,JLEV)
    IF (JLEV /= KLEV + 1) THEN
      THERMODYNAMICS%TEMPERATURE_HL(JLON,JLEV) = PTEMPERATURE_H(JLON,JLEV)
    ELSE
      ! IFS currently sets the half-level temperature at the surface to be
      ! equal to the skin temperature. The radiation scheme takes as input
      ! only the half-level temperatures and assumes the Planck function to
      ! vary linearly in optical depth between half levels. In the lowest
      ! atmospheric layer, where the atmospheric temperature can be much
      ! cooler than the skin temperature, this can lead to significant
      ! differences between the effective temperature of this lowest layer
      ! and the true value in the model.
      !
      ! We may approximate the temperature profile in the lowest model level
      ! as piecewise linear between the top of the layer T[k-1/2], the
      ! centre of the layer T[k] and the base of the layer Tskin.  The mean
      ! temperature of the layer is then 0.25*T[k-1/2] + 0.5*T[k] +
      ! 0.25*Tskin, which can be achieved by setting the atmospheric
      ! temperature at the half-level corresponding to the surface as
      ! follows:
      THERMODYNAMICS%TEMPERATURE_HL(JLON,KLEV+1)&
          &  = PTEMPERATURE(JLON,KLEV)&
          &  + 0.5_JPRB * (PTEMPERATURE_H(JLON,KLEV+1)&
          &               -PTEMPERATURE_H(JLON,KLEV))

      ! Alternatively we respect the model's atmospheric temperature in the
      ! lowest model level by setting the temperature at the lowest
      ! half-level such that the mean temperature of the layer is correct:
      !THERMODYNAMICS%TEMPERATURE_HL(JLON,KLEV+1) &
      !     &  = 2.0_JPRB * PTEMPERATURE(JLON,KLEV) &
      !     &             - PTEMPERATURE_H(JLON,KLEV)
    ENDIF
  ENDDO
ENDDO
!$ACC END PARALLEL
!$OMP END TARGET TEAMS DISTRIBUTE PARALLEL DO

! Alternatively we respect the model's atmospheric temperature in the
! lowest model level by setting the temperature at the lowest
! half-level such that the mean temperature of the layer is correct:
!thermodynamics%temperature_hl(KIDIA:KFDIA,KLEV+1) &
!     &  = 2.0_JPRB * PTEMPERATURE(KIDIA:KFDIA,KLEV) &
!     &             - PTEMPERATURE_H(KIDIA:KFDIA,KLEV)

! Compute saturation specific humidity, used to hydrate aerosols. The
! "2" for the last argument indicates that the routine is not being
! called from within the convection scheme.
!CALL SATUR(KIDIA, KFDIA, KLON, 1, KLEV, .false., &
!     &  PPRESSURE, PTEMPERATURE, THERMODYNAMICS%H2O_SAT_LIQ, 2)
! Alternative approximate version using temperature and pressure from
! the thermodynamics structure
CALL thermodynamics%calc_saturation_wrt_liquid(thermodynamics, KIDIA, KFDIA, LACC=LLACC)

#ifdef HAVE_NVTX
call nvtxEndRange
call nvtxStartRange("single level setup")
#endif

! Set single-level fileds
SINGLE_LEVEL%SOLAR_IRRADIANCE              = PSOLAR_IRRADIANCE
!$ACC UPDATE DEVICE(SINGLE_LEVEL%SOLAR_IRRADIANCE) ASYNC(1) IF(LLACC)
!$OMP TARGET UPDATE TO(SINGLE_LEVEL%SOLAR_IRRADIANCE) IF(LLACC)

!$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(LLACC)
!$ACC LOOP GANG VECTOR
!$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO IF(LLACC)
DO JLON = KIDIA,KFDIA
  SINGLE_LEVEL%COS_SZA(JLON)          = PMU0(JLON)
  SINGLE_LEVEL%SKIN_TEMPERATURE(JLON) = PTEMPERATURE_SKIN(JLON)
ENDDO
!$OMP END TARGET TEAMS DISTRIBUTE PARALLEL DO

!$ACC LOOP GANG VECTOR COLLAPSE(2)
!$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO COLLAPSE(2) IF(LLACC)
DO JBAND = 1,YRERAD%NSW
  DO JLON = KIDIA,KFDIA
    SINGLE_LEVEL%SW_ALBEDO(JLON,JBAND)      = PALBEDO_DIF(JLON,JBAND)
    SINGLE_LEVEL%SW_ALBEDO_DIRECT(JLON,JBAND)=PALBEDO_DIR(JLON,JBAND)
  ENDDO
ENDDO
!$OMP END TARGET TEAMS DISTRIBUTE PARALLEL DO

!$ACC LOOP GANG VECTOR COLLAPSE(2)
!$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO COLLAPSE(2) IF(LLACC)
DO JBAND = 1,YRERAD%NLWEMISS
  DO JLON = KIDIA,KFDIA
    ! Spectral longwave emissivity
    SINGLE_LEVEL%LW_EMISSIVITY(JLON,JBAND)  = PSPECTRALEMISS(JLON,JBAND)
  ENDDO
ENDDO
!$ACC END PARALLEL
!$OMP END TARGET TEAMS DISTRIBUTE PARALLEL DO

! ! Do we adjust the solar spectrum for time in the solar cycle? Only
! ! available with the ecCKD gas optics scheme.
! IF (YRERAD%LSPECTRALSOLARCYCLE) THEN
!   ! Solar cycles are counted from solar minima, so whole numbers lead to a
!   ! multiplier of minus one (solar minimum), values ending in .5 lead
!   ! to a multiplier of plus one (solar maximum)
!   SINGLE_LEVEL%SPECTRAL_SOLAR_CYCLE_MULTIPLIER = YDERDI%RSOLARCYCLEMULT
! ELSE
  SINGLE_LEVEL%SPECTRAL_SOLAR_CYCLE_MULTIPLIER = 0.0_JPRB
! ENDIF
!$ACC UPDATE DEVICE(SINGLE_LEVEL%SPECTRAL_SOLAR_CYCLE_MULTIPLIER) IF(LLACC)
!$OMP TARGET UPDATE TO(SINGLE_LEVEL%SPECTRAL_SOLAR_CYCLE_MULTIPLIER) IF(LLACC)

! Create the relevant seed from date and time get the starting day
! and number of minutes since start
! IDAY = NDD(NINDAT)
! ITIM = NINT(NSTEP * YDMODEL%YRML_GCONF%YRRIP%TSTEP / 60.0_JPRB)
! DO JLON = KIDIA, KFDIA
!   ! This method gives a unique value for roughly every 1-km square
!   ! on the globe and every minute.  ASIN(PGEMU)*60 gives rough
!   ! latitude in degrees, which we multiply by 100 to give a unique
!   ! value for roughly every km. PGELAM*60*100 gives a unique number
!   ! for roughly every km of longitude around the equator, which we
!   ! multiply by 180*100 so there is no overlap with the latitude
!   ! values.  The result can be contained in a 32-byte integer (but
!   ! since random numbers are generated with the help of integer
!   ! overflow, it should not matter if the number did overflow).
!   SINGLE_LEVEL%ISEED(JLON) = ITIM + IDAY &
!        &  +  NINT(PGELAM(JLON)*108000000.0_JPRD &
!        &          + ASIN(PGEMU(JLON))*6000.0_JPRD)
! ENDDO

! Simple initialization of the seeds for the Monte Carlo scheme
call single_level%init_seed_simple(single_level, kidia, kfdia, lacc=.true.)

! Added for bit-identity validation against ecrad standalone:
! Overwrite seed with user-specified values
if (present(iseed)) then
  !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
  !$ACC LOOP GANG VECTOR
  !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO
  DO JLON = KIDIA,KFDIA
    single_level%iseed(jlon) = iseed(jlon)
  ENDDO
  !$ACC END PARALLEL
  !$OMP END TARGET TEAMS DISTRIBUTE PARALLEL DO
end if

! Set the solar spectrum scaling, if required
IF (YRERAD%NSOLARSPECTRUM /= 0) THEN
  ALLOCATE(SINGLE_LEVEL%SPECTRAL_SOLAR_SCALING(RAD_CONFIG%N_BANDS_SW))
  ! RRTMG uses the old Kurucz solar spectrum. The following scalings
  ! adjust it to match more recent measured spectra.
  IF (YRERAD%NSOLARSPECTRUM == 1) THEN
    ! The Whole Heliosphere Interval (WHI) 2008 reference spectrum for
    ! solar minimum conditions in 2008:
    ! https://lasp.colorado.edu/lisird/data/whi_ref_spectra This
    ! spectrum only extends to wavelengths of 2.4 microns, so a Kurucz
    ! is assumed to be correct for longer wavelengths. (Note that in
    ! previous cycles this was incorrectly labelled as the Coddington
    ! spectrum, which is below.)
    SINGLE_LEVEL%SPECTRAL_SOLAR_SCALING &
         &  = [ 1.0000_JPRB,  1.0000_JPRB,  1.0000_JPRB,  1.0478_JPRB, &
         &      1.0404_JPRB,  1.0317_JPRB,  1.0231_JPRB,  1.0054_JPRB, &
         &      0.98413_JPRB, 0.99863_JPRB, 0.99907_JPRB, 0.90589_JPRB, &
         &      0.92213_JPRB, 1.0000_JPRB ]
  ELSE
    ! The average of the last 33 years (3 solar cycles) of the
    ! Coddington et al. (BAMS, 2016) climate data record, which covers
    ! the entire spectrum
    SINGLE_LEVEL%SPECTRAL_SOLAR_SCALING &
         &  = [ 0.99892_JPRB, 0.99625_JPRB, 1.00822_JPRB, 1.01587_JPRB, &
         &      1.01898_JPRB, 1.01044_JPRB, 1.08441_JPRB, 0.99398_JPRB, &
         &      1.00553_JPRB, 0.99533_JPRB, 1.01509_JPRB, 0.92331_JPRB, &
         &      0.92681_JPRB, 0.99749_JPRB ]
  ENDIF
  !$ACC ENTER DATA COPYIN(SINGLE_LEVEL%SPECTRAL_SOLAR_SCALING) ASYNC(1) IF(LLACC)
  !$OMP TARGET UPDATE TO(SINGLE_LEVEL%SPECTRAL_SOLAR_SCALING) IF(LLACC)
ENDIF

#ifdef HAVE_NVTX
call nvtxEndRange
call nvtxStartRange("cloud setup")
#endif

! Set cloud fields
!$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(LLACC)
!$ACC LOOP GANG VECTOR COLLAPSE(2)
!$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO COLLAPSE(2) IF(LLACC)
DO JLEV = 1,KLEV
  DO JLON = KIDIA,KFDIA
    YLCLOUD%Q_LIQ(JLON,JLEV)    = PQ_LIQUID(JLON,JLEV)
    YLCLOUD%Q_ICE(JLON,JLEV)    = PQ_ICE(JLON,JLEV) + PQ_SNOW(JLON,JLEV)
    YLCLOUD%FRACTION(JLON,JLEV) = PCLOUD_FRAC(JLON,JLEV)
  ENDDO
ENDDO
!$ACC END PARALLEL
!$OMP END TARGET TEAMS DISTRIBUTE PARALLEL DO

! Compute effective radii and convert to metres
IF (PRESENT(PRE_LIQ)) THEN
  !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
  !$ACC LOOP GANG VECTOR COLLAPSE(2)
  !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO COLLAPSE(2)
  DO JLEV = 1,KLEV
    DO JLON = KIDIA,KFDIA
      YLCLOUD%RE_LIQ(JLON,JLEV) = PRE_LIQ(JLON,JLEV)
    ENDDO
  ENDDO
  !$ACC END PARALLEL
  !$OMP END TARGET TEAMS DISTRIBUTE PARALLEL DO
ELSE
  CALL LIQUID_EFFECTIVE_RADIUS(YRERAD, &
     &  KIDIA, KFDIA, KLON, KLEV, &
     &  PPRESSURE, PTEMPERATURE, PCLOUD_FRAC, PQ_LIQUID, PQ_RAIN, &
     &  PLAND_SEA_MASK, PCCN_LAND, PCCN_SEA, &
     &  ZRE_LIQUID_UM, LACC=LLACC)
  !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(LLACC)
  !$ACC LOOP GANG VECTOR COLLAPSE(2)
  !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO COLLAPSE(2) IF (LLACC)
  DO JLEV = 1,KLEV
    DO JLON = KIDIA,KFDIA
      YLCLOUD%RE_LIQ(JLON,JLEV) = ZRE_LIQUID_UM(JLON,JLEV) * 1.0E-6_JPRB
    ENDDO
  ENDDO
  !$ACC END PARALLEL
  !$OMP END TARGET TEAMS DISTRIBUTE PARALLEL DO
ENDIF

IF (PRESENT(PRE_ICE)) THEN
  !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
  !$ACC LOOP GANG VECTOR COLLAPSE(2)
  !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO COLLAPSE(2)
  DO JLEV = 1,KLEV
    DO JLON = KIDIA,KFDIA
      YLCLOUD%RE_ICE(JLON,JLEV) = PRE_ICE(JLON,JLEV)
    ENDDO
  ENDDO
  !$ACC END PARALLEL
  !$OMP END TARGET TEAMS DISTRIBUTE PARALLEL DO
ELSE
  CALL ICE_EFFECTIVE_RADIUS(YRERAD, KIDIA, KFDIA, KLON, KLEV, &
      &  PPRESSURE, PTEMPERATURE, PCLOUD_FRAC, PQ_ICE, PQ_SNOW, PGEMU, &
      &  ZRE_ICE_UM, LACC=LLACC)
  !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(LLACC)
  !$ACC LOOP GANG VECTOR COLLAPSE(2)
  !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO COLLAPSE(2) IF (LLACC)
  DO JLEV = 1,KLEV
    DO JLON = KIDIA,KFDIA
      YLCLOUD%RE_ICE(JLON,JLEV) = ZRE_ICE_UM(JLON,JLEV) * 1.0E-6_JPRB
    ENDDO
  ENDDO
  !$ACC END PARALLEL
  !$OMP END TARGET TEAMS DISTRIBUTE PARALLEL DO
ENDIF

! Get the cloud overlap decorrelation length (for cloud boundaries),
! in km, according to the parameterization specified by NDECOLAT,
! and insert into the "cloud" object. Also get the ratio of
! decorrelation lengths for cloud water content inhomogeneities and
! cloud boundaries, and set it in the "rad_config" object.
CALL CLOUD_OVERLAP_DECORR_LEN(KIDIA,KFDIA,KLON, &
     &  PGEMU,YRERAD%NDECOLAT, &
     &  PDECORR_LEN_EDGES_KM=ZDECORR_LEN_KM, PDECORR_LEN_RATIO=ZDECORR_LEN_RATIO, &
     &  LACC=LLACC)


! Compute cloud overlap parameter from decorrelation length
RAD_CONFIG%CLOUD_INHOM_DECORR_SCALING = ZDECORR_LEN_RATIO
!$ACC UPDATE DEVICE(RAD_CONFIG%CLOUD_INHOM_DECORR_SCALING) ASYNC(1) IF(LLACC)
!$OMP TARGET UPDATE TO(RAD_CONFIG%CLOUD_INHOM_DECORR_SCALING) IF(LLACC)

!$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(LLACC)
!$ACC LOOP GANG VECTOR
!$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO IF (LLACC)
DO JLON = KIDIA,KFDIA
  ZDECORR_LEN_KM(JLON) = 1000.0_JPRB*ZDECORR_LEN_KM(JLON)
ENDDO
!$ACC END PARALLEL
!$OMP END TARGET TEAMS DISTRIBUTE PARALLEL DO
CALL YLCLOUD%SET_OVERLAP_PARAM(YLCLOUD, THERMODYNAMICS,&
      &                       ZDECORR_LEN_KM,&
      &                       ISTARTCOL=KIDIA, IENDCOL=KFDIA, &
      &                       LACC=LLACC)
! Added for bit-identity validation against ecrad standalone:
! Overwrite overlap param with provided value
if(present(PCLOUD_OVERLAP)) then
  !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
  !$ACC LOOP GANG VECTOR COLLAPSE(2)
  !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO COLLAPSE(2)
  DO JLEV = 1,KLEV-1
    DO JLON = KIDIA,KFDIA
      YLCLOUD%OVERLAP_PARAM(JLON,JLEV) = PCLOUD_OVERLAP(JLON,JLEV)
    ENDDO
  ENDDO
  !$ACC END PARALLEL
  !$OMP END TARGET TEAMS DISTRIBUTE PARALLEL DO
endif

! Cloud water content fractional standard deviation is configurable
! from namelist NAERAD but must be globally constant. Before it was
! hard coded at 1.0.
CALL YLCLOUD%CREATE_FRACTIONAL_STD(YLCLOUD, KLON, KLEV, YRERAD%RCLOUD_FRAC_STD, LACC=LLACC)


IF (         RAD_CONFIG%I_SOLVER_LW == ISOLVERSPARTACUS &
     &  .OR. RAD_CONFIG%I_SOLVER_SW == ISOLVERSPARTACUS) THEN
  ! We are using the SPARTACUS solver so need to specify cloud scale,
  ! and use Mark Fielding's parameterization based on ARM data
  CALL YLCLOUD%PARAM_CLOUD_EFFECTIVE_SEPARATION_ETA(YLCLOUD, KLON, KLEV, &
       &  PPRESSURE_H, YRERAD%RCLOUD_SEPARATION_SCALE_SURF, &
       &  YRERAD%RCLOUD_SEPARATION_SCALE_TOA, 3.5_JPRB, 0.75_JPRB, &
       &  KIDIA, KFDIA)
ENDIF

! Compute the dry mass of each layer neglecting humidity effects, in
! kg m-2, needed to scale some of the aerosol inputs
CALL THERMODYNAMICS%GET_LAYER_MASS(THERMODYNAMICS, KIDIA,KFDIA,ZLAYER_MASS, LACC=LLACC)

#ifdef HAVE_NVTX
call nvtxEndRange
call nvtxStartRange("aerosol setup")
#endif

! Copy over aerosol mass mixing ratio
IF (.NOT. LL_USE_TEGEN_AEROSOLS) THEN
  ! MACC aerosol from climatology or prognostic aerosol variables -
  ! this is already in mass mixing ratio units with the required array
  ! orientation so we can copy it over directly
  ! AB need to cap the minimum mass mixing ratio/AOD to avoid instability
  ! in case of negative values in input
  !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(LLACC)
  !$ACC LOOP GANG VECTOR COLLAPSE(3)
  !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO COLLAPSE(3) IF(LLACC)
  DO JAER = 1,KAEROSOL
    DO JLEV = 1,KLEV
      DO JLON = KIDIA,KFDIA
        AEROSOL%MIXING_RATIO(JLON,JLEV,JAER) = MAX(PAEROSOL(JLON,JLEV,JAER),0.0_JPRB)
      ENDDO
    ENDDO
  ENDDO
  !$ACC END PARALLEL
  !$OMP END TARGET TEAMS DISTRIBUTE PARALLEL DO

  ! Add the tropospheric and stratospheric backgrounds contained in the
  ! old Tegen arrays - this is very ugly!
  IF (TROP_BG_AER_MASS_EXT > 0.0_JPRB) THEN
    !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(LLACC)
    !$ACC LOOP GANG VECTOR COLLAPSE(2)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO COLLAPSE(2) IF(LLACC)
    DO JLEV = 1,KLEV
      DO JLON = KIDIA,KFDIA
        AEROSOL%MIXING_RATIO(JLON,JLEV,ITYPE_TROP_BG_AER)&
            &  = AEROSOL%MIXING_RATIO(JLON,JLEV,ITYPE_TROP_BG_AER)&
            &  + PAEROSOL_OLD(JLON,1,JLEV)&
            &  / (ZLAYER_MASS(JLON,JLEV) * TROP_BG_AER_MASS_EXT)
      ENDDO
    ENDDO
    !$ACC END PARALLEL
    !$OMP END TARGET TEAMS DISTRIBUTE PARALLEL DO
  ENDIF
  IF (STRAT_BG_AER_MASS_EXT > 0.0_JPRB) THEN
    !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(LLACC)
    !$ACC LOOP GANG VECTOR COLLAPSE(2)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO COLLAPSE(2) IF(LLACC)
    DO JLEV = 1,KLEV
      DO JLON = KIDIA,KFDIA
        AEROSOL%MIXING_RATIO(JLON,JLEV,ITYPE_STRAT_BG_AER)&
            &  = AEROSOL%MIXING_RATIO(JLON,JLEV,ITYPE_STRAT_BG_AER)&
            &  + PAEROSOL_OLD(JLON,6,JLEV)&
            &  / (ZLAYER_MASS(JLON,JLEV) * STRAT_BG_AER_MASS_EXT)
      ENDDO
    ENDDO
    !$ACC END PARALLEL
    !$OMP END TARGET TEAMS DISTRIBUTE PARALLEL DO
  ENDIF
ELSE

  ! Tegen aerosol climatology - the array PAEROSOL_OLD contains the
  ! 550-nm optical depth in each layer. The optics data file
  ! aerosol_ifs_rrtm_tegen.nc does not contain mass extinction
  ! coefficient, but a scaling factor that the 550-nm optical depth
  ! should be multiplied by to obtain the optical depth in each
  ! spectral band.  Therefore, in order for the units to work out, we
  ! need to divide by the layer mass (in kg m-2) to obtain the 550-nm
  ! cross-section per unit mass of dry air (so in m2 kg-1).  We also
  ! need to permute the array.
  !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(LLACC)
  !$ACC LOOP GANG VECTOR COLLAPSE(3)
  !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO COLLAPSE(3) IF(LLACC)
  DO JLEV = 1,KLEV
    DO JAER = 1,6
      DO JLON=KIDIA,KFDIA
        AEROSOL%MIXING_RATIO(JLON,JLEV,JAER)&
          &  = PAEROSOL_OLD(JLON,JAER,JLEV)&
          &  / ZLAYER_MASS(JLON,JLEV)
      ENDDO
    ENDDO
  ENDDO
  !$ACC END PARALLEL
  !$OMP END TARGET TEAMS DISTRIBUTE PARALLEL DO

ENDIF

#ifdef HAVE_NVTX
call nvtxEndRange
call nvtxStartRange("gas setup")
#endif

! Insert gas mixing ratios (mol mol-1 or kg kg-1)
CALL GAS%PUT(GAS, IH2O,    IMASSMIXINGRATIO, PQ, LACC=LLACC)
CALL GAS%PUT(GAS, ICO2,    IMASSMIXINGRATIO, PCO2, LACC=LLACC)
CALL GAS%PUT(GAS, ICH4,    IMASSMIXINGRATIO, PCH4, LACC=LLACC)
CALL GAS%PUT(GAS, IN2O,    IMASSMIXINGRATIO, PN2O, LACC=LLACC)
CALL GAS%PUT(GAS, ICFC11,  IMASSMIXINGRATIO, PCFC11, LACC=LLACC)
CALL GAS%PUT(GAS, ICFC12,  IMASSMIXINGRATIO, PCFC12, LACC=LLACC)
CALL GAS%PUT(GAS, IHCFC22, IMASSMIXINGRATIO, PHCFC22, LACC=LLACC)
CALL GAS%PUT(GAS, ICCL4,   IMASSMIXINGRATIO, PCCL4, LACC=LLACC)
CALL GAS%PUT(GAS, IO3,     IMASSMIXINGRATIO, PO3, LACC=LLACC)
CALL GAS%PUT_WELL_MIXED(GAS, IO2, IVOLUMEMIXINGRATIO, 0.20944_JPRB, LACC=LLACC)

! Ensure the units of the gas mixing ratios are what is required by
! the gas absorption model
CALL SET_GAS_UNITS(RAD_CONFIG, GAS, LACC=LLACC)

#ifdef HAVE_NVTX
call nvtxEndRange
#endif

! Call radiation scheme
#ifdef HAVE_NVTX
call nvtxStartRange("radiation")
#endif

CALL RADIATION(KLON, KLEV, KIDIA, KFDIA, RAD_CONFIG,&
     &  SINGLE_LEVEL, THERMODYNAMICS, GAS, YLCLOUD, AEROSOL, FLUX)

#ifdef HAVE_NVTX
call nvtxEndRange
#endif

! Check fluxes are within physical bounds
IF (YRERAD%NDUMPBADINPUTS /= 0 &
     &  .AND. (N_BAD_FLUXES == 0 .OR. N_BAD_FLUXES < YRERAD%NDUMPBADINPUTS)) THEN
  IF (FLUX%OUT_OF_PHYSICAL_BOUNDS(KIDIA,KFDIA, LACC=LLACC)) THEN
#ifndef OMPGPU
!$OMP CRITICAL
#endif
#if defined(_OPENACC) || defined(OMPGPU)
    IF(LLACC) THEN
      !$ACC WAIT
      CALL SINGLE_LEVEL%UPDATE_HOST(SINGLE_LEVEL)
      CALL THERMODYNAMICS%UPDATE_HOST(THERMODYNAMICS)
      CALL GAS%UPDATE_HOST(GAS)
      CALL YLCLOUD%UPDATE_HOST(YLCLOUD)
      CALL AEROSOL%UPDATE_HOST(AEROSOL)
      CALL FLUX%UPDATE_HOST(FLUX)
      !$ACC UPDATE HOST(PGEMU, PGELAM)
      !$ACC WAIT
      !$OMP TARGET UPDATE FROM(PGEMU, PGELAM)
    ENDIF
#endif
    N_BAD_FLUXES = N_BAD_FLUXES+1
    WRITE(CL_FILE_NAME, '(A,I0,A,I0,A)') '/home/parr/ifs_dump/bad_inputs_', &
         &  MPL_MYRANK(), '_', N_BAD_FLUXES, '.nc'
    WRITE(NULERR,*) '  Writing ', TRIM(CL_FILE_NAME)
    ! Implicit assumption that KFDIA==KLON
    CALL SAVE_INPUTS(TRIM(CL_FILE_NAME), RAD_CONFIG, SINGLE_LEVEL, &
         &  THERMODYNAMICS, GAS, YLCLOUD, AEROSOL, &
         &  LAT=ASIN(PGEMU)*180.0_JPRB/RPI, LON=PGELAM*180.0_JPRB/RPI, IVERBOSE=3)
    WRITE(CL_FILE_NAME, '(A,I0,A,I0,A)') '/home/parr/ifs_dump/bad_outputs_', &
         &  MPL_MYRANK(), '_', N_BAD_FLUXES, '.nc'
    WRITE(NULERR,*) '  Writing ', TRIM(CL_FILE_NAME)
    CALL SAVE_FLUXES(TRIM(CL_FILE_NAME), RAD_CONFIG, THERMODYNAMICS, FLUX, IVERBOSE=3)
    IF (YRERAD%NDUMPBADINPUTS < 0) THEN
      ! Abort on the first set of bad fluxes
      CALL ABOR1("RADIATION_SCHEME: ABORT DUE TO FLUXES OUT OF PHYSICAL BOUNDS")
    ENDIF
#ifndef OMPGPU
!$OMP END CRITICAL
#endif
  ENDIF
ENDIF

! For debugging, do we store a certain number of inputs and outputs
! regardless of whether bad fluxes have been detected?
IF (N_OUTPUT_FLUXES < YRERAD%NDUMPINPUTS) THEN
#ifndef OMPGPU
!$OMP CRITICAL
#endif
#if defined(_OPENACC) || defined(OMPGPU)
    IF(LLACC) THEN
      !$ACC WAIT
      CALL SINGLE_LEVEL%UPDATE_HOST(SINGLE_LEVEL)
      CALL THERMODYNAMICS%UPDATE_HOST(THERMODYNAMICS)
      CALL GAS%UPDATE_HOST(GAS)
      CALL YLCLOUD%UPDATE_HOST(YLCLOUD)
      CALL AEROSOL%UPDATE_HOST(AEROSOL)
      !$ACC UPDATE HOST(PGEMU, PGELAM)
      !$ACC WAIT
      !$OMP TARGET UPDATE FROM(PGEMU, PGELAM)
    ENDIF
#endif
  N_OUTPUT_FLUXES = N_OUTPUT_FLUXES+1
  WRITE(CL_FILE_NAME, '(A,I0,A,I0,A)') '/home/parr/ifs_dump/inputs_', &
       &  MPL_MYRANK(), '_', N_OUTPUT_FLUXES, '.nc'
  WRITE(NULERR,*) '  Writing ', TRIM(CL_FILE_NAME)
  ! Implicit assumption that KFDIA==KLON
  CALL SAVE_INPUTS(TRIM(CL_FILE_NAME), RAD_CONFIG, SINGLE_LEVEL, &
       &  THERMODYNAMICS, GAS, YLCLOUD, AEROSOL, &
       &  LAT=ASIN(PGEMU)*180.0_JPRB/RPI, LON=PGELAM*180.0_JPRB/RPI, IVERBOSE=3)
  WRITE(CL_FILE_NAME, '(A,I0,A,I0,A)') '/home/parr/ifs_dump/outputs_', &
       &  MPL_MYRANK(), '_', N_OUTPUT_FLUXES, '.nc'
  WRITE(NULERR,*) '  Writing ', TRIM(CL_FILE_NAME)
  CALL SAVE_FLUXES(TRIM(CL_FILE_NAME), RAD_CONFIG, THERMODYNAMICS, FLUX, IVERBOSE=3)
#ifndef OMPGPU
!$OMP END CRITICAL
#endif
ENDIF

#ifdef HAVE_NVTX
call nvtxStartRange("compute fluxes")
#endif
! Compute required output fluxes
! First the net fluxes
!$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(LLACC)
!$ACC LOOP GANG VECTOR COLLAPSE(2)
!$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO COLLAPSE(2) IF(LLACC)
DO JLEV=1,KLEV+1
  DO JLON=KIDIA,KFDIA
    PFLUX_SW(JLON,JLEV) = FLUX%SW_DN(JLON,JLEV) - FLUX%SW_UP(JLON,JLEV)
    PFLUX_LW(JLON,JLEV) = FLUX%LW_DN(JLON,JLEV) - FLUX%LW_UP(JLON,JLEV)
    PFLUX_SW_CLEAR(JLON,JLEV)&
        &  = FLUX%SW_DN_CLEAR(JLON,JLEV) - FLUX%SW_UP_CLEAR(JLON,JLEV)
    PFLUX_LW_CLEAR(JLON,JLEV)&
        &  = FLUX%LW_DN_CLEAR(JLON,JLEV) - FLUX%LW_UP_CLEAR(JLON,JLEV)
  ENDDO
ENDDO
!$OMP END TARGET TEAMS DISTRIBUTE PARALLEL DO

!$ACC LOOP GANG VECTOR PRIVATE(ZBLACK_BODY_NET_LW)
!$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO PRIVATE(ZBLACK_BODY_NET_LW) FIRSTPRIVATE(ZSIGMA)
DO JLON=KIDIA,KFDIA
  ! Now the surface fluxes
  PFLUX_SW_DN      (JLON) = FLUX%SW_DN             (JLON,KLEV+1)
  ! Broadband longwave flux
  PFLUX_LW_DN      (JLON) = FLUX%LW_DN             (JLON,KLEV+1)
  PFLUX_SW_DN_CLEAR(JLON) = FLUX%SW_DN_CLEAR       (JLON,KLEV+1)
  PFLUX_LW_DN_CLEAR(JLON) = FLUX%LW_DN_CLEAR       (JLON,KLEV+1)
  PFLUX_DIR        (JLON) = FLUX%SW_DN_DIRECT      (JLON,KLEV+1)
  PFLUX_DIR_CLEAR  (JLON) = FLUX%SW_DN_DIRECT_CLEAR(JLON,KLEV+1)
  PFLUX_DIR_INTO_SUN(JLON) = 0.0_JPRB
  IF (PMU0(JLON) > EPSILON(1.0_JPRB)) THEN
    PFLUX_DIR_INTO_SUN(JLON) = PFLUX_DIR(JLON) / PMU0(JLON)
  ENDIF
  ! Top-of-atmosphere downwelling flux
  PFLUX_SW_DN_TOA(JLON) = FLUX%SW_DN(JLON,1)

  ! Compute UV fluxes as weighted sum of appropriate shortwave bands
  PFLUX_UV       (JLON) = 0.0_JPRB
  !$ACC LOOP SEQ
  DO JBAND = 1,NWEIGHT_UV
    PFLUX_UV(JLON) = PFLUX_UV(JLON) + WEIGHT_UV(JBAND)&
        &  * FLUX%SW_DN_SURF_BAND(IBAND_UV(JBAND),JLON)
  ENDDO

  ! Compute photosynthetically active radiation similarly
  PFLUX_PAR      (JLON) = 0.0_JPRB
  PFLUX_PAR_CLEAR(JLON) = 0.0_JPRB
  !$ACC LOOP SEQ
  DO JBAND = 1,NWEIGHT_PAR
    PFLUX_PAR(JLON) = PFLUX_PAR(JLON) + WEIGHT_PAR(JBAND)&
        &  * FLUX%SW_DN_SURF_BAND(IBAND_PAR(JBAND),JLON)
    PFLUX_PAR_CLEAR(JLON) = PFLUX_PAR_CLEAR(JLON)&
        &  + WEIGHT_PAR(JBAND)&
        &  * FLUX%SW_DN_SURF_CLEAR_BAND(IBAND_PAR(JBAND),JLON)
  ENDDO

  ! Compute effective broadband emissivity. This is only approximate -
  ! due to spectral variations in emissivity, it is not in general
  ! possible to provide a broadband emissivity that can reproduce the
  ! upwelling surface flux given the downwelling flux and the skin
  ! temperature.
  ZBLACK_BODY_NET_LW = PFLUX_LW_DN(JLON) &
      &  - ZSIGMA*PTEMPERATURE_SKIN(JLON)**4
  PEMIS_OUT(JLON) = PSPECTRALEMISS(JLON,1) ! Default value
  IF (ABS(ZBLACK_BODY_NET_LW) > 1.0E-5) THEN
    ! This calculation can go outside the range of any individual
    ! spectral emissivity value, so needs to be capped
    PEMIS_OUT(JLON) = MAX(0.8_JPRB, MIN(0.99_JPRB, PFLUX_LW(JLON,KLEV+1) / ZBLACK_BODY_NET_LW))
  ENDIF
ENDDO
!$OMP END TARGET TEAMS DISTRIBUTE PARALLEL DO

! Copy longwave derivatives
IF (YRERAD%LAPPROXLWUPDATE) THEN
  !$ACC LOOP GANG VECTOR COLLAPSE(2)
  !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO COLLAPSE(2)
  DO JLEV=1,KLEV+1
    DO JLON=KIDIA,KFDIA
        PLWDERIVATIVE(JLON,JLEV) = FLUX%LW_DERIVATIVES(JLON,JLEV)
    ENDDO
  ENDDO
  !$OMP END TARGET TEAMS DISTRIBUTE PARALLEL DO
ENDIF

! Store the shortwave downwelling fluxes in each albedo band
IF (YRERAD%LAPPROXSWUPDATE) THEN
  !$ACC LOOP GANG VECTOR COLLAPSE(2)
  !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO COLLAPSE(2)
  DO JBAND=1,YRERAD%NSW
    DO JLON=KIDIA,KFDIA
      PSWDIFFUSEBAND(JLON,JBAND) = FLUX%SW_DN_DIFFUSE_SURF_CANOPY(JBAND,JLON)
      PSWDIRECTBAND (JLON,JBAND) = FLUX%SW_DN_DIRECT_SURF_CANOPY (JBAND,JLON)
    ENDDO
  ENDDO
  !$OMP END TARGET TEAMS DISTRIBUTE PARALLEL DO
ENDIF
!$ACC END PARALLEL

#ifdef HAVE_NVTX
call nvtxEndRange
call nvtxStartRange("cleanup")
#endif

!$ACC END DATA
!$OMP TARGET EXIT DATA MAP(DELETE:ZRE_LIQUID_UM, ZRE_ICE_UM, ZDECORR_LEN_KM,ZLAYER_MASS)
IF (LLACC) THEN
  CALL FLUX%DELETE_DEVICE(FLUX)
  CALL YLCLOUD%DELETE_DEVICE(YLCLOUD)
  CALL AEROSOL%DELETE_DEVICE(AEROSOL)
  CALL GAS%DELETE_DEVICE(GAS)
  CALL THERMODYNAMICS%DELETE_DEVICE(THERMODYNAMICS)
  CALL SINGLE_LEVEL%DELETE_DEVICE(SINGLE_LEVEL)
  CALL YRADIATION%RAD_CONFIG%DELETE_DEVICE(YRADIATION%RAD_CONFIG)
ENDIF

!$ACC END DATA
!$ACC WAIT
!$OMP TARGET EXIT DATA MAP(DELETE:YRADIATION, YRERAD, RAD_CONFIG, SINGLE_LEVEL, THERMODYNAMICS, GAS, AEROSOL, YLCLOUD, FLUX)

CALL SINGLE_LEVEL%DEALLOCATE
CALL THERMODYNAMICS%DEALLOCATE
CALL GAS%DEALLOCATE
CALL YLCLOUD%DEALLOCATE
CALL AEROSOL%DEALLOCATE
CALL FLUX%DEALLOCATE

#ifdef HAVE_NVTX
call nvtxEndRange
#endif

END ASSOCIATE

IF (LHOOK) CALL DR_HOOK('RADIATION_SCHEME',1,ZHOOK_HANDLE)

END SUBROUTINE RADIATION_SCHEME
