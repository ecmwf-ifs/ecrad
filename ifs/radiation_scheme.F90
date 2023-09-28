SUBROUTINE RADIATION_SCHEME &
     & (YRADIATION,KIDIA, KFDIA, KLON, KLEV, KAEROSOL, &
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
     &  PSWDIFFUSEBAND, PSWDIRECTBAND)

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
!
!-----------------------------------------------------------------------

! Modules from ifs or ifsaux libraries
USE PARKIND1       , ONLY : JPIM, JPRB, JPRD
USE YOMHOOK        , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOMCST         , ONLY : RPI, RSIGMA ! Stefan-Boltzmann constant
USE YOMLUN         , ONLY : NULERR
USE RADIATION_SETUP, ONLY : ITYPE_TROP_BG_AER, ITYPE_STRAT_BG_AER, TRADIATION

! Modules from ecRad radiation library
USE RADIATION_CONFIG,         ONLY : ISOLVERSPARTACUS
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

IMPLICIT NONE

! INPUT ARGUMENTS

TYPE(TRADIATION), INTENT(IN)    :: YRADIATION

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
REAL(KIND=JPRB),  INTENT(OUT) :: PFLUX_SW(KLON,KLEV+1)
REAL(KIND=JPRB),  INTENT(OUT) :: PFLUX_LW(KLON,KLEV+1)
REAL(KIND=JPRB),  INTENT(OUT) :: PFLUX_SW_CLEAR(KLON,KLEV+1)
REAL(KIND=JPRB),  INTENT(OUT) :: PFLUX_LW_CLEAR(KLON,KLEV+1)

! *** Surface flux components (W m-2)
REAL(KIND=JPRB),  INTENT(OUT) :: PFLUX_SW_DN(KLON)
REAL(KIND=JPRB),  INTENT(OUT) :: PFLUX_LW_DN(KLON)
REAL(KIND=JPRB),  INTENT(OUT) :: PFLUX_SW_DN_CLEAR(KLON)
REAL(KIND=JPRB),  INTENT(OUT) :: PFLUX_LW_DN_CLEAR(KLON)
! Direct component of surface flux into horizontal plane
REAL(KIND=JPRB),  INTENT(OUT) :: PFLUX_DIR(KLON)
REAL(KIND=JPRB),  INTENT(OUT) :: PFLUX_DIR_CLEAR(KLON)
! As PFLUX_DIR but into a plane perpendicular to the sun
REAL(KIND=JPRB),  INTENT(OUT) :: PFLUX_DIR_INTO_SUN(KLON)

! *** Ultraviolet and photosynthetically active radiation (W m-2)
REAL(KIND=JPRB),  INTENT(OUT) :: PFLUX_UV(KLON)
REAL(KIND=JPRB),  INTENT(OUT) :: PFLUX_PAR(KLON)
REAL(KIND=JPRB),  INTENT(OUT) :: PFLUX_PAR_CLEAR(KLON)

! *** Other single-level diagnostics
! Top-of-atmosphere incident solar flux (W m-2)
REAL(KIND=JPRB),  INTENT(OUT) :: PFLUX_SW_DN_TOA(KLON)
! Diagnosed longwave surface emissivity across the whole spectrum
REAL(KIND=JPRB),  INTENT(OUT) :: PEMIS_OUT(KLON)

! Partial derivative of total-sky longwave upward flux at each level
! with respect to upward flux at surface, used to correct heating
! rates at gridpoints/timesteps between calls to the full radiation
! scheme.  Note that this version uses the convention of level index
! increasing downwards, unlike the local variable ZLwDerivative that
! is returned from the LW radiation scheme.
REAL(KIND=JPRB),  INTENT(OUT) :: PLWDERIVATIVE(KLON,KLEV+1)

! Surface diffuse and direct downwelling shortwave flux in each
! shortwave albedo band, used in RADINTG to update the surface fluxes
! accounting for high-resolution albedo information
REAL(KIND=JPRB),  INTENT(OUT) :: PSWDIFFUSEBAND(KLON,YRADIATION%YRERAD%NSW)
REAL(KIND=JPRB),  INTENT(OUT) :: PSWDIRECTBAND (KLON,YRADIATION%YRERAD%NSW)

! LOCAL VARIABLES
TYPE(SINGLE_LEVEL_TYPE)   :: SINGLE_LEVEL
TYPE(THERMODYNAMICS_TYPE) :: THERMODYNAMICS
TYPE(GAS_TYPE)            :: GAS
TYPE(CLOUD_TYPE)          :: YLCLOUD
TYPE(AEROSOL_TYPE)        :: AEROSOL
TYPE(FLUX_TYPE)           :: FLUX

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
REAL(KIND=JPRB)           :: ZBLACK_BODY_NET_LW(KIDIA:KFDIA)

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

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

! Dummy from YOMCT3
! INTEGER(KIND=JPIM) :: NSTEP = 0

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

ASSOCIATE(YRERAD    =>YRADIATION%YRERAD, &
     &    RAD_CONFIG=>YRADIATION%RAD_CONFIG, &
     &    NWEIGHT_UV=>YRADIATION%NWEIGHT_UV, &
     &    IBAND_UV  =>YRADIATION%IBAND_UV(:), &
     &    WEIGHT_UV =>YRADIATION%WEIGHT_UV(:), &
     &    NWEIGHT_PAR=>YRADIATION%NWEIGHT_PAR, &
     &    IBAND_PAR =>YRADIATION%IBAND_PAR(:), &
     &    WEIGHT_PAR=>YRADIATION%WEIGHT_PAR(:), &
     &    TROP_BG_AER_MASS_EXT=>YRADIATION%TROP_BG_AER_MASS_EXT, &
     &    STRAT_BG_AER_MASS_EXT=>YRADIATION%STRAT_BG_AER_MASS_EXT)
! Allocate memory in radiation objects
CALL SINGLE_LEVEL%ALLOCATE(KLON, YRERAD%NSW, YRERAD%NLWEMISS, &
     &                     USE_SW_ALBEDO_DIRECT=.TRUE.)
CALL THERMODYNAMICS%ALLOCATE(KLON, KLEV, USE_H2O_SAT=.TRUE.)
CALL GAS%ALLOCATE(KLON, KLEV)
CALL YLCLOUD%ALLOCATE(KLON, KLEV)
IF (YRERAD%NAERMACC == 1) THEN
  CALL AEROSOL%ALLOCATE(KLON, 1, KLEV, KAEROSOL) ! MACC aerosols
ELSE
  CALL AEROSOL%ALLOCATE(KLON, 1, KLEV, 6) ! Tegen climatology
ENDIF
CALL FLUX%ALLOCATE(RAD_CONFIG, 1, KLON, KLEV)

! Set thermodynamic profiles: simply copy over the half-level
! pressure and temperature
THERMODYNAMICS%PRESSURE_HL   (KIDIA:KFDIA,:) = PPRESSURE_H   (KIDIA:KFDIA,:)
THERMODYNAMICS%TEMPERATURE_HL(KIDIA:KFDIA,:) = PTEMPERATURE_H(KIDIA:KFDIA,:)

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
THERMODYNAMICS%TEMPERATURE_HL(KIDIA:KFDIA,KLEV+1)&
     &  = PTEMPERATURE(KIDIA:KFDIA,KLEV)&
     &  + 0.5_JPRB * (PTEMPERATURE_H(KIDIA:KFDIA,KLEV+1)&
     &               -PTEMPERATURE_H(KIDIA:KFDIA,KLEV))

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
CALL thermodynamics%calc_saturation_wrt_liquid(KIDIA, KFDIA)

! Set single-level fileds
SINGLE_LEVEL%SOLAR_IRRADIANCE              = PSOLAR_IRRADIANCE
SINGLE_LEVEL%COS_SZA(KIDIA:KFDIA)          = PMU0(KIDIA:KFDIA)
SINGLE_LEVEL%SKIN_TEMPERATURE(KIDIA:KFDIA) = PTEMPERATURE_SKIN(KIDIA:KFDIA)
SINGLE_LEVEL%SW_ALBEDO(KIDIA:KFDIA,:)      = PALBEDO_DIF(KIDIA:KFDIA,:)
SINGLE_LEVEL%SW_ALBEDO_DIRECT(KIDIA:KFDIA,:)=PALBEDO_DIR(KIDIA:KFDIA,:)
! Spectral longwave emissivity
SINGLE_LEVEL%LW_EMISSIVITY(KIDIA:KFDIA,:)  = PSPECTRALEMISS(KIDIA:KFDIA,:)

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
call single_level%init_seed_simple(kidia, kfdia)

! Set the solar spectrum scaling, if required
IF (YRERAD%NSOLARSPECTRUM == 1) THEN
  ALLOCATE(SINGLE_LEVEL%SPECTRAL_SOLAR_SCALING(RAD_CONFIG%N_BANDS_SW))
  ! Ratio of SORCE (Coddington et al. 2016) and Kurucz solar spectra
  SINGLE_LEVEL%SPECTRAL_SOLAR_SCALING &
       &  = (/  1.0, 1.0, 1.0, 1.0478, 1.0404, 1.0317, 1.0231, &
       &        1.0054, 0.98413, 0.99863, 0.99907, 0.90589, 0.92213, 1.0 /)
ENDIF

! Set cloud fields
YLCLOUD%Q_LIQ(KIDIA:KFDIA,:)    = PQ_LIQUID(KIDIA:KFDIA,:)
YLCLOUD%Q_ICE(KIDIA:KFDIA,:)    = PQ_ICE(KIDIA:KFDIA,:) + PQ_SNOW(KIDIA:KFDIA,:)
YLCLOUD%FRACTION(KIDIA:KFDIA,:) = PCLOUD_FRAC(KIDIA:KFDIA,:)

! Compute effective radii and convert to metres
CALL LIQUID_EFFECTIVE_RADIUS(YRERAD, &
     &  KIDIA, KFDIA, KLON, KLEV, &
     &  PPRESSURE, PTEMPERATURE, PCLOUD_FRAC, PQ_LIQUID, PQ_RAIN, &
     &  PLAND_SEA_MASK, PCCN_LAND, PCCN_SEA, &
     &  ZRE_LIQUID_UM) !, PPERT=PPERT)
YLCLOUD%RE_LIQ(KIDIA:KFDIA,:) = ZRE_LIQUID_UM(KIDIA:KFDIA,:) * 1.0E-6_JPRB

CALL ICE_EFFECTIVE_RADIUS(YRERAD, KIDIA, KFDIA, KLON, KLEV, &
     &  PPRESSURE, PTEMPERATURE, PCLOUD_FRAC, PQ_ICE, PQ_SNOW, PGEMU, &
     &  ZRE_ICE_UM) !, PPERT=PPERT)
YLCLOUD%RE_ICE(KIDIA:KFDIA,:) = ZRE_ICE_UM(KIDIA:KFDIA,:) * 1.0E-6_JPRB

! Get the cloud overlap decorrelation length (for cloud boundaries),
! in km, according to the parameterization specified by NDECOLAT,
! and insert into the "cloud" object. Also get the ratio of
! decorrelation lengths for cloud water content inhomogeneities and
! cloud boundaries, and set it in the "rad_config" object.
CALL CLOUD_OVERLAP_DECORR_LEN(KIDIA,KFDIA,KLON, &
     &  PGEMU,YRERAD%NDECOLAT, &
     &  PDECORR_LEN_EDGES_KM=ZDECORR_LEN_KM, PDECORR_LEN_RATIO=ZDECORR_LEN_RATIO)

! Compute cloud overlap parameter from decorrelation length
!RAD_CONFIG%CLOUD_INHOM_DECORR_SCALING = ZDECORR_LEN_RATIO
DO JLON = KIDIA,KFDIA
  CALL YLCLOUD%SET_OVERLAP_PARAM(THERMODYNAMICS,&
      &                       ZDECORR_LEN_KM(JLON)*1000.0_JPRB,&
      &                       ISTARTCOL=JLON, IENDCOL=JLON)
ENDDO
! Or we can call the routine on all columns at once
!CALL YLCLOUD%SET_OVERLAP_PARAM(THERMODYNAMICS,&
!     &                       ZDECORR_LEN_KM(KIDIA:KFDIA)*1000.0_JPRB,&
!     &                       ISTARTCOL=KIDIA, IENDCOL=KFDIA)

! Cloud water content fractional standard deviation is configurable
! from namelist NAERAD but must be globally constant. Before it was
! hard coded at 1.0.
CALL YLCLOUD%CREATE_FRACTIONAL_STD(KLON, KLEV, YRERAD%RCLOUD_FRAC_STD)


IF (         RAD_CONFIG%I_SOLVER_LW == ISOLVERSPARTACUS &
     &  .OR. RAD_CONFIG%I_SOLVER_SW == ISOLVERSPARTACUS) THEN
  ! We are using the SPARTACUS solver so need to specify cloud scale,
  ! and use Mark Fielding's parameterization based on ARM data
  CALL YLCLOUD%PARAM_CLOUD_EFFECTIVE_SEPARATION_ETA(KLON, KLEV, &
       &  PPRESSURE_H, YRERAD%RCLOUD_SEPARATION_SCALE_SURF, &
       &  YRERAD%RCLOUD_SEPARATION_SCALE_TOA, 3.5_JPRB, 0.75_JPRB, &
       &  KIDIA, KFDIA)
ENDIF

! Compute the dry mass of each layer neglecting humidity effects, in
! kg m-2, needed to scale some of the aerosol inputs
CALL THERMODYNAMICS%GET_LAYER_MASS(KIDIA,KFDIA,ZLAYER_MASS)

! Copy over aerosol mass mixing ratio
IF (YRERAD%NAERMACC == 1) THEN


  ! MACC aerosol from climatology or prognostic aerosol variables -
  ! this is already in mass mixing ratio units with the required array
  ! orientation so we can copy it over directly
  ! AB need to cap the minimum mass mixing ratio/AOD to avoid instability
  ! in case of negative values in input
  DO JAER = 1,KAEROSOL
    DO JLEV = 1,KLEV
      DO JLON = KIDIA,KFDIA
        AEROSOL%MIXING_RATIO(JLON,JLEV,JAER) = MAX(PAEROSOL(JLON,JLEV,JAER),0.0_JPRB)
      ENDDO
    ENDDO
  ENDDO

  IF (YRERAD%NAERMACC == 1) THEN
    ! Add the tropospheric and stratospheric backgrounds contained in the
    ! old Tegen arrays - this is very ugly!
    IF (TROP_BG_AER_MASS_EXT > 0.0_JPRB) THEN
      AEROSOL%MIXING_RATIO(KIDIA:KFDIA,:,ITYPE_TROP_BG_AER)&
           &  = AEROSOL%MIXING_RATIO(KIDIA:KFDIA,:,ITYPE_TROP_BG_AER)&
           &  + PAEROSOL_OLD(KIDIA:KFDIA,1,:)&
           &  / (ZLAYER_MASS * TROP_BG_AER_MASS_EXT)
    ENDIF
    IF (STRAT_BG_AER_MASS_EXT > 0.0_JPRB) THEN
      AEROSOL%MIXING_RATIO(KIDIA:KFDIA,:,ITYPE_STRAT_BG_AER)&
           &  = AEROSOL%MIXING_RATIO(KIDIA:KFDIA,:,ITYPE_STRAT_BG_AER)&
           &  + PAEROSOL_OLD(KIDIA:KFDIA,6,:)&
           &  / (ZLAYER_MASS * STRAT_BG_AER_MASS_EXT)
    ENDIF
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
  DO JLEV = 1,KLEV
    DO JAER = 1,6
      AEROSOL%MIXING_RATIO(KIDIA:KFDIA,JLEV,JAER)&
         &  = PAEROSOL_OLD(KIDIA:KFDIA,JAER,JLEV)&
         &  / ZLAYER_MASS(KIDIA:KFDIA,JLEV)
    ENDDO
  ENDDO

ENDIF

! Insert gas mixing ratios
CALL GAS%PUT(IH2O,    IMASSMIXINGRATIO, PQ)
CALL GAS%PUT(ICO2,    IMASSMIXINGRATIO, PCO2)
CALL GAS%PUT(ICH4,    IMASSMIXINGRATIO, PCH4)
CALL GAS%PUT(IN2O,    IMASSMIXINGRATIO, PN2O)
CALL GAS%PUT(ICFC11,  IMASSMIXINGRATIO, PCFC11)
CALL GAS%PUT(ICFC12,  IMASSMIXINGRATIO, PCFC12)
CALL GAS%PUT(IHCFC22, IMASSMIXINGRATIO, PHCFC22)
CALL GAS%PUT(ICCL4,   IMASSMIXINGRATIO, PCCL4)
CALL GAS%PUT(IO3,     IMASSMIXINGRATIO, PO3)
CALL GAS%PUT_WELL_MIXED(IO2, IVOLUMEMIXINGRATIO, 0.20944_JPRB)

! Ensure the units of the gas mixing ratios are what is required by
! the gas absorption model
CALL SET_GAS_UNITS(RAD_CONFIG, GAS)

!call save_inputs('inputs_ifs.nc', rad_config, single_level, thermodynamics, &
!     &           gas, ylcloud, aerosol, &
!     &           lat=spread(0.0_jprb,1,klon), &
!     &           lon=spread(0.0_jprb,1,klon), &
!     &           iverbose=2)

! Call radiation scheme
CALL RADIATION(KLON, KLEV, KIDIA, KFDIA, RAD_CONFIG,&
     &  SINGLE_LEVEL, THERMODYNAMICS, GAS, YLCLOUD, AEROSOL, FLUX)

! Check fluxes are within physical bounds
IF (YRERAD%NDUMPBADINPUTS /= 0 &
     &  .AND. (N_BAD_FLUXES == 0 .OR. N_BAD_FLUXES < YRERAD%NDUMPBADINPUTS)) THEN
  IF (FLUX%OUT_OF_PHYSICAL_BOUNDS(KIDIA,KFDIA)) THEN
!$OMP CRITICAL
    N_BAD_FLUXES = N_BAD_FLUXES+1
    WRITE(CL_FILE_NAME, '(A,I0,A,I0,A)') '/home/parr/ifs_dump/bad_inputs_', &
         &  MPL_MYRANK(), '_', N_BAD_FLUXES, '.nc'
    WRITE(NULERR,*) '  Writing ', TRIM(CL_FILE_NAME)
    ! Implicit assumption that KFDIA==KLON
    CALL SAVE_INPUTS(TRIM(CL_FILE_NAME), RAD_CONFIG, SINGLE_LEVEL, &
         &  THERMODYNAMICS, GAS, YLCLOUD, AEROSOL, &
         &  LAT=ASIN(PGEMU)*180.0/RPI, LON=PGELAM*180.0/RPI, IVERBOSE=3)
    WRITE(CL_FILE_NAME, '(A,I0,A,I0,A)') '/home/parr/ifs_dump/bad_outputs_', &
         &  MPL_MYRANK(), '_', N_BAD_FLUXES, '.nc'
    WRITE(NULERR,*) '  Writing ', TRIM(CL_FILE_NAME)
    CALL SAVE_FLUXES(TRIM(CL_FILE_NAME), RAD_CONFIG, THERMODYNAMICS, FLUX, IVERBOSE=3)
    IF (YRERAD%NDUMPBADINPUTS < 0) THEN
      ! Abort on the first set of bad fluxes
      CALL ABOR1("RADIATION_SCHEME: ABORT DUE TO FLUXES OUT OF PHYSICAL BOUNDS")
    ENDIF
!$OMP END CRITICAL
  ENDIF
ENDIF

! For debugging, do we store a certain number of inputs and outputs
! regardless of whether bad fluxes have been detected?
IF (N_OUTPUT_FLUXES < YRERAD%NDUMPINPUTS) THEN
!$OMP CRITICAL
  N_OUTPUT_FLUXES = N_OUTPUT_FLUXES+1
  WRITE(CL_FILE_NAME, '(A,I0,A,I0,A)') '/home/parr/ifs_dump/inputs_', &
       &  MPL_MYRANK(), '_', N_OUTPUT_FLUXES, '.nc'
  WRITE(NULERR,*) '  Writing ', TRIM(CL_FILE_NAME)
  ! Implicit assumption that KFDIA==KLON
  CALL SAVE_INPUTS(TRIM(CL_FILE_NAME), RAD_CONFIG, SINGLE_LEVEL, &
       &  THERMODYNAMICS, GAS, YLCLOUD, AEROSOL, &
       &  LAT=ASIN(PGEMU)*180.0/RPI, LON=PGELAM*180.0/RPI, IVERBOSE=3)
  WRITE(CL_FILE_NAME, '(A,I0,A,I0,A)') '/home/parr/ifs_dump/outputs_', &
       &  MPL_MYRANK(), '_', N_OUTPUT_FLUXES, '.nc'
  WRITE(NULERR,*) '  Writing ', TRIM(CL_FILE_NAME)
  CALL SAVE_FLUXES(TRIM(CL_FILE_NAME), RAD_CONFIG, THERMODYNAMICS, FLUX, IVERBOSE=3)
!$OMP END CRITICAL
ENDIF

! Compute required output fluxes
! First the net fluxes
PFLUX_SW(KIDIA:KFDIA,:) = FLUX%SW_DN(KIDIA:KFDIA,:) - FLUX%SW_UP(KIDIA:KFDIA,:)
PFLUX_LW(KIDIA:KFDIA,:) = FLUX%LW_DN(KIDIA:KFDIA,:) - FLUX%LW_UP(KIDIA:KFDIA,:)
PFLUX_SW_CLEAR(KIDIA:KFDIA,:)&
     &  = FLUX%SW_DN_CLEAR(KIDIA:KFDIA,:) - FLUX%SW_UP_CLEAR(KIDIA:KFDIA,:)
PFLUX_LW_CLEAR(KIDIA:KFDIA,:)&
     &  = FLUX%LW_DN_CLEAR(KIDIA:KFDIA,:) - FLUX%LW_UP_CLEAR(KIDIA:KFDIA,:)
! Now the surface fluxes
PFLUX_SW_DN      (KIDIA:KFDIA) = FLUX%SW_DN             (KIDIA:KFDIA,KLEV+1)
PFLUX_LW_DN      (KIDIA:KFDIA) = FLUX%LW_DN             (KIDIA:KFDIA,KLEV+1)
PFLUX_SW_DN_CLEAR(KIDIA:KFDIA) = FLUX%SW_DN_CLEAR       (KIDIA:KFDIA,KLEV+1)
PFLUX_LW_DN_CLEAR(KIDIA:KFDIA) = FLUX%LW_DN_CLEAR       (KIDIA:KFDIA,KLEV+1)
PFLUX_DIR        (KIDIA:KFDIA) = FLUX%SW_DN_DIRECT      (KIDIA:KFDIA,KLEV+1)
PFLUX_DIR_CLEAR  (KIDIA:KFDIA) = FLUX%SW_DN_DIRECT_CLEAR(KIDIA:KFDIA,KLEV+1)
PFLUX_DIR_INTO_SUN(KIDIA:KFDIA) = 0.0_JPRB
WHERE (PMU0(KIDIA:KFDIA) > EPSILON(1.0_JPRB))
  PFLUX_DIR_INTO_SUN(KIDIA:KFDIA) = PFLUX_DIR(KIDIA:KFDIA) / PMU0(KIDIA:KFDIA)
ENDWHERE
! Top-of-atmosphere downwelling flux
PFLUX_SW_DN_TOA(KIDIA:KFDIA) = FLUX%SW_DN(KIDIA:KFDIA,1)

! Compute UV fluxes as weighted sum of appropriate shortwave bands
PFLUX_UV       (KIDIA:KFDIA) = 0.0_JPRB
DO JBAND = 1,NWEIGHT_UV
!DEC$ IVDEP
  PFLUX_UV(KIDIA:KFDIA) = PFLUX_UV(KIDIA:KFDIA) + WEIGHT_UV(JBAND)&
       &  * FLUX%SW_DN_SURF_BAND(IBAND_UV(JBAND),KIDIA:KFDIA)
ENDDO

! Compute photosynthetically active radiation similarly
PFLUX_PAR      (KIDIA:KFDIA) = 0.0_JPRB
PFLUX_PAR_CLEAR(KIDIA:KFDIA) = 0.0_JPRB
DO JBAND = 1,NWEIGHT_PAR
!DEC$ IVDEP
  PFLUX_PAR(KIDIA:KFDIA) = PFLUX_PAR(KIDIA:KFDIA) + WEIGHT_PAR(JBAND)&
       &  * FLUX%SW_DN_SURF_BAND(IBAND_PAR(JBAND),KIDIA:KFDIA)
!DEC$ IVDEP
  PFLUX_PAR_CLEAR(KIDIA:KFDIA) = PFLUX_PAR_CLEAR(KIDIA:KFDIA)&
       &  + WEIGHT_PAR(JBAND)&
       &  * FLUX%SW_DN_SURF_CLEAR_BAND(IBAND_PAR(JBAND),KIDIA:KFDIA)
ENDDO

! Compute effective broadband emissivity. This is only approximate -
! due to spectral variations in emissivity, it is not in general
! possible to provide a broadband emissivity that can reproduce the
! upwelling surface flux given the downwelling flux and the skin
! temperature.
ZBLACK_BODY_NET_LW = PFLUX_LW_DN(KIDIA:KFDIA) &
     &  - RSIGMA*PTEMPERATURE_SKIN(KIDIA:KFDIA)**4
PEMIS_OUT(KIDIA:KFDIA) = PSPECTRALEMISS(KIDIA:KFDIA,1) ! Default value
WHERE (ABS(ZBLACK_BODY_NET_LW) > 1.0E-5)
  ! This calculation can go outside the range of any individual
  ! spectral emissivity value, so needs to be capped
  PEMIS_OUT(KIDIA:KFDIA) = MAX(0.8_JPRB, MIN(0.99_JPRB, PFLUX_LW(KIDIA:KFDIA,KLEV+1) / ZBLACK_BODY_NET_LW))
ENDWHERE

! Copy longwave derivatives
IF (YRERAD%LAPPROXLWUPDATE) THEN
  PLWDERIVATIVE(KIDIA:KFDIA,:) = FLUX%LW_DERIVATIVES(KIDIA:KFDIA,:)
ENDIF

! Store the shortwave downwelling fluxes in each albedo band
IF (YRERAD%LAPPROXSWUPDATE) THEN
  PSWDIFFUSEBAND(KIDIA:KFDIA,:) = TRANSPOSE(FLUX%SW_DN_DIFFUSE_SURF_CANOPY(:,KIDIA:KFDIA))
  PSWDIRECTBAND (KIDIA:KFDIA,:) = TRANSPOSE(FLUX%SW_DN_DIRECT_SURF_CANOPY (:,KIDIA:KFDIA))
ENDIF

CALL SINGLE_LEVEL%DEALLOCATE
CALL THERMODYNAMICS%DEALLOCATE
CALL GAS%DEALLOCATE
CALL YLCLOUD%DEALLOCATE
CALL AEROSOL%DEALLOCATE
CALL FLUX%DEALLOCATE

END ASSOCIATE

IF (LHOOK) CALL DR_HOOK('RADIATION_SCHEME',1,ZHOOK_HANDLE)

END SUBROUTINE RADIATION_SCHEME
