SUBROUTINE RADIATION_SCHEME &
     & (KIDIA, KFDIA, KLON, KLEV, KAEROSOL, &
     &  PSOLAR_IRRADIANCE, &
     &  PMU0, PTEMPERATURE_SKIN, PALBEDO_DIF, PALBEDO_DIR, &
     &  PEMIS, PEMIS_WINDOW, &
     &  PCCN_LAND, PCCN_SEA, &
     &  PGELAM, PGEMU, PLAND_SEA_MASK, &
     &  PPRESSURE, PTEMPERATURE, &
     &  PPRESSURE_H, PTEMPERATURE_H, &
     &  PQ, PCO2, PCH4, PN2O, PNO2, PCFC11, PCFC12, PHCFC22, PCCL4, PO3_DP, &
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
!    should have been run first.
!
! AUTHOR
! ------
!   Robin Hogan, ECMWF
!   Original: 2015-09-16
!
! MODIFICATIONS
! -------------
!
! TO DO
! -----
!
!-----------------------------------------------------------------------

! Modules from ifs or ifsaux libraries
USE PARKIND1 , ONLY : JPIM, JPRB
USE YOMHOOK  , ONLY : LHOOK, DR_HOOK
USE YOERAD   , ONLY : YRERAD
USE RADIATION_SETUP, ONLY : rad_config, &
     &  NWEIGHT_UV,  IBAND_UV,  WEIGHT_UV, &
     &  NWEIGHT_PAR, IBAND_PAR, WEIGHT_PAR, &
     &  ITYPE_TROP_BG_AER,  TROP_BG_AER_MASS_EXT, &
     &  ITYPE_STRAT_BG_AER, STRAT_BG_AER_MASS_EXT
USE YOMRIP0  , ONLY : NINDAT
USE YOMCT3   , ONLY : NSTEP
USE YOMRIP   , ONLY : YRRIP
USE YOMCST   , ONLY : RSIGMA ! Stefan-Boltzmann constant

! Modules from radiation library
USE radiation_single_level,   ONLY : single_level_type
USE radiation_thermodynamics, ONLY : thermodynamics_type
USE radiation_gas
USE radiation_cloud,          ONLY : cloud_type
USE radiation_aerosol,        ONLY : aerosol_type
USE radiation_flux,           ONLY : flux_type
USE radiation_interface,      ONLY : radiation, set_gas_units
USE radiation_save,           ONLY : save_inputs

IMPLICIT NONE

! INPUT ARGUMENTS

! *** Array dimensions and ranges
INTEGER(KIND=JPIM),INTENT(IN) :: KIDIA    ! Start column to process
INTEGER(KIND=JPIM),INTENT(IN) :: KFDIA    ! End column to process
INTEGER(KIND=JPIM),INTENT(IN) :: KLON     ! Number of columns
INTEGER(KIND=JPIM),INTENT(IN) :: KLEV     ! Number of levels
INTEGER(KIND=JPIM),INTENT(IN) :: KAEROSOL ! Number of aerosol types

! *** Single-level fields
REAL(KIND=JPRB),   INTENT(IN) :: PSOLAR_IRRADIANCE ! (W m-2)
REAL(KIND=JPRB),   INTENT(IN) :: PMU0(KLON) ! Cosine of solar zenith ang
REAL(KIND=JPRB),   INTENT(IN) :: PTEMPERATURE_SKIN(KLON) ! (K)
! Diffuse and direct components of surface shortwave albedo
REAL(KIND=JPRB),   INTENT(IN) :: PALBEDO_DIF(KLON,YRERAD%NSW)
REAL(KIND=JPRB),   INTENT(IN) :: PALBEDO_DIR(KLON,YRERAD%NSW)
! Longwave emissivity outside and inside the window region
REAL(KIND=JPRB),   INTENT(IN) :: PEMIS(KLON)
REAL(KIND=JPRB),   INTENT(IN) :: PEMIS_WINDOW(KLON)
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
REAL(KIND=JPRB),   INTENT(IN) :: PO3_DP(KLON,KLEV) ! (Pa*kg/kg) !

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
REAL(KIND=JPRB),  INTENT(OUT) :: PSWDIFFUSEBAND(KLON,YRERAD%NSW)
REAL(KIND=JPRB),  INTENT(OUT) :: PSWDIRECTBAND (KLON,YRERAD%NSW)

! LOCAL VARIABLES
TYPE(single_level_type)   :: single_level
TYPE(thermodynamics_type) :: thermodynamics
TYPE(gas_type)            :: gas
TYPE(cloud_type)          :: cloud
TYPE(aerosol_type)        :: aerosol
TYPE(flux_type)           :: flux

! Mass mixing ratio of ozone (kg/kg)
REAL(KIND=JPRB)           :: ZO3(KLON,KLEV)

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
INTEGER :: ITIM, IDAY

! Loop indices
INTEGER :: JLON, JLEV, JBAND, JB_ALBEDO, JAER

REAL(KIND=JPRB) :: ZHOOK_HANDLE

! Import time functions for iseed calculation
#include "fcttim.func.h"

#include "liquid_effective_radius.intfb.h"
#include "ice_effective_radius.intfb.h"
#include "cloud_overlap_decorr_len.intfb.h"
#include "satur.intfb.h"

IF (LHOOK) CALL DR_HOOK('RADIATION_SCHEME',0,ZHOOK_HANDLE)

! Allocate memory in radiation objects
CALL single_level%allocate(KLON, YRERAD%NSW, 2, &
     &                     use_sw_albedo_direct=.TRUE.)
CALL thermodynamics%allocate(KLON, KLEV, use_h2o_sat=.true.)
CALL gas%allocate(KLON, KLEV)
CALL cloud%allocate(KLON, KLEV)
IF (YRERAD%NAERMACC > 0) THEN
  CALL aerosol%allocate(KLON, 1, KLEV, KAEROSOL) ! MACC climatology
ELSE
  CALL aerosol%allocate(KLON, 1, KLEV, 6) ! Tegen climatology
ENDIF
CALL flux%allocate(rad_config, 1, KLON, KLEV)

! Set thermodynamic profiles: simply copy over the half-level
! pressure and temperature
thermodynamics%pressure_hl   (KIDIA:KFDIA,:) = PPRESSURE_H   (KIDIA:KFDIA,:)
thermodynamics%temperature_hl(KIDIA:KFDIA,:) = PTEMPERATURE_H(KIDIA:KFDIA,:)

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
thermodynamics%temperature_hl(KIDIA:KFDIA,KLEV+1) &
     &  = PTEMPERATURE(KIDIA:KFDIA,KLEV) &
     &  + 0.5_JPRB * (PTEMPERATURE_H(KIDIA:KFDIA,KLEV+1) &
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
CALL SATUR(KIDIA, KFDIA, KLON, 1, KLEV, &
     &  PPRESSURE, PTEMPERATURE, thermodynamics%h2o_sat_liq, 2)  
! Alternative approximate version using temperature and pressure from
! the thermodynamics structure
!CALL thermodynamics%calc_saturation_wrt_liquid(KIDIA, KFDIA)

! Set single-level fileds
single_level%solar_irradiance              = PSOLAR_IRRADIANCE
single_level%cos_sza(KIDIA:KFDIA)          = PMU0(KIDIA:KFDIA)
single_level%skin_temperature(KIDIA:KFDIA) = PTEMPERATURE_SKIN(KIDIA:KFDIA)
single_level%sw_albedo(KIDIA:KFDIA,:)      = PALBEDO_DIF(KIDIA:KFDIA,:)
single_level%sw_albedo_direct(KIDIA:KFDIA,:)=PALBEDO_DIR(KIDIA:KFDIA,:)
! Longwave emissivity is in two bands
single_level%lw_emissivity(KIDIA:KFDIA,1)  = PEMIS(KIDIA:KFDIA)
single_level%lw_emissivity(KIDIA:KFDIA,2)  = PEMIS_WINDOW(KIDIA:KFDIA)

! Create the relevant seed from date and time get the starting day
! and number of minutes since start
IDAY = NDD(NINDAT)
ITIM = NINT(NSTEP * YRRIP%TSTEP / 60.0_JPRB)
DO JLON = KIDIA, KFDIA
  ! This method gives a unique value for roughly every 1-km square
  ! on the globe and every minute.  ASIN(PGEMU)*60 gives rough
  ! latitude in degrees, which we multiply by 100 to give a unique
  ! value for roughly every km. PGELAM*60*100 gives a unique number
  ! for roughly every km of longitude around the equator, which we
  ! multiply by 180*100 so there is no overlap with the latitude
  ! values.  The result can be contained in a 32-byte integer (but
  ! since random numbers are generated with the help of integer
  ! overflow, it should not matter if the number did overflow).
  single_level%iseed(JLON) = ITIM + IDAY & 
       &  +  NINT(PGELAM(JLON)*108000000.0_JPRB &
       &          + ASIN(PGEMU(JLON))*6000.0_JPRB)
ENDDO

! Set cloud fields
cloud%q_liq(KIDIA:KFDIA,:)    = PQ_LIQUID(KIDIA:KFDIA,:)
cloud%q_ice(KIDIA:KFDIA,:)    = PQ_ICE(KIDIA:KFDIA,:) + PQ_SNOW(KIDIA:KFDIA,:)
cloud%fraction(KIDIA:KFDIA,:) = PCLOUD_FRAC(KIDIA:KFDIA,:)

! Compute effective radii and convert to metres
CALL LIQUID_EFFECTIVE_RADIUS(KIDIA, KFDIA, KLON, KLEV, &
     &  PPRESSURE, PTEMPERATURE, PCLOUD_FRAC, PQ_LIQUID, PQ_RAIN, &
     &  PLAND_SEA_MASK, PCCN_LAND, PCCN_SEA, &
     &  ZRE_LIQUID_UM)
cloud%re_liq(KIDIA:KFDIA,:) = ZRE_LIQUID_UM(KIDIA:KFDIA,:) * 1.0e-6_JPRB

CALL ICE_EFFECTIVE_RADIUS(KIDIA, KFDIA, KLON, KLEV, &
     &  PPRESSURE, PTEMPERATURE, PCLOUD_FRAC, PQ_ICE, PQ_SNOW, PGEMU, &
     &  ZRE_ICE_UM)
cloud%re_ice(KIDIA:KFDIA,:) = ZRE_ICE_UM(KIDIA:KFDIA,:) * 1.0e-6_JPRB

! Get the cloud overlap decorrelation length (for cloud boundaries),
! in km, according to the parameterization specified by NDECOLAT,
! and insert into the "cloud" object. Also get the ratio of
! decorrelation lengths for cloud water content inhomogeneities and
! cloud boundaries, and set it in the "rad_config" object.
CALL CLOUD_OVERLAP_DECORR_LEN(KIDIA, KFDIA, KLON, PGEMU, YRERAD%NDECOLAT, &
     &    ZDECORR_LEN_KM, PDECORR_LEN_RATIO=ZDECORR_LEN_RATIO)
rad_config%cloud_inhom_decorr_scaling = ZDECORR_LEN_RATIO
DO JLON = KIDIA,KFDIA
  CALL cloud%set_overlap_param(thermodynamics, &
       &                       ZDECORR_LEN_KM(JLON)*1000.0_JPRB, &
       &                       istartcol=JLON, iendcol=JLON)
ENDDO

! Cloud water content fractional standard deviation is configurable
! from namelist NAERAD but must be globally constant. Before it was
! hard coded at 1.0.
CALL cloud%create_fractional_std(KLON, KLEV, YRERAD%RCLOUD_FRAC_STD)

! By default mid and high cloud effective size is 10 km
CALL cloud%create_inv_cloud_effective_size(KLON,KLEV,1.0_JPRB/10000.0_JPRB)
! But for boundary clouds (eta > 0.8) we set it to 1 km
DO JLEV = 1,KLEV
  DO JLON = KIDIA,KFDIA
    IF (PPRESSURE(JLON,JLEV) > 0.8_JPRB * PPRESSURE_H(JLON,KLEV+1)) THEN
      cloud%inv_cloud_effective_size(JLON,JLEV) = 1.0e-3_JPRB
    ENDIF
  ENDDO
ENDDO


! Compute the dry mass of each layer neglecting humidity effects, in
! kg m-2, needed to scale some of the aerosol inputs
CALL thermodynamics%get_layer_mass(ZLAYER_MASS)

! Copy over aerosol mass mixing ratio
IF (YRERAD%NAERMACC > 0) THEN

  ! MACC aerosol climatology - this is already in mass mixing ratio
  ! units with the required array orientation so we can copy it over
  ! directly
  aerosol%mixing_ratio(KIDIA:KFDIA,:,:) = PAEROSOL(KIDIA:KFDIA,:,:)

  ! Add the tropospheric and stratospheric backgrounds contained in the
  ! old Tegen arrays - this is very ugly!
  IF (TROP_BG_AER_MASS_EXT > 0.0_JPRB) THEN
    aerosol%mixing_ratio(KIDIA:KFDIA,:,ITYPE_TROP_BG_AER) &
         &  = aerosol%mixing_ratio(KIDIA:KFDIA,:,ITYPE_TROP_BG_AER) &
         &  + PAEROSOL_OLD(KIDIA:KFDIA,1,:) &
         &  / (ZLAYER_MASS * TROP_BG_AER_MASS_EXT)
  ENDIF
  IF (STRAT_BG_AER_MASS_EXT > 0.0_JPRB) THEN
    aerosol%mixing_ratio(KIDIA:KFDIA,:,ITYPE_STRAT_BG_AER) &
         &  = aerosol%mixing_ratio(KIDIA:KFDIA,:,ITYPE_STRAT_BG_AER) &
         &  + PAEROSOL_OLD(KIDIA:KFDIA,6,:) &
         &  / (ZLAYER_MASS * STRAT_BG_AER_MASS_EXT)
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
      aerosol%mixing_ratio(KIDIA:KFDIA,JLEV,JAER) &
         &  = PAEROSOL_OLD(KIDIA:KFDIA,JAER,JLEV) &
         &  / ZLAYER_MASS(KIDIA:KFDIA,JLEV)
    ENDDO
  ENDDO

ENDIF


! Convert ozone Pa*kg/kg to kg/kg
DO JLEV = 1,KLEV
  DO JLON = KIDIA,KFDIA
    ZO3(JLON,JLEV) = PO3_DP(JLON,JLEV) &
         &         / (PPRESSURE_H(JLON,JLEV+1)-PPRESSURE_H(JLON,JLEV))
  ENDDO
ENDDO

! Insert gas mixing ratios
CALL gas%put(IH2O,    IMassMixingRatio, PQ)
CALL gas%put(ICO2,    IMassMixingRatio, PCO2)
CALL gas%put(ICH4,    IMassMixingRatio, PCH4)
CALL gas%put(IN2O,    IMassMixingRatio, PN2O)
CALL gas%put(ICFC11,  IMassMixingRatio, PCFC11)
CALL gas%put(ICFC12,  IMassMixingRatio, PCFC12)
CALL gas%put(IHCFC22, IMassMixingRatio, PHCFC22)
CALL gas%put(ICCL4,   IMassMixingRatio, PCCL4)
CALL gas%put(IO3,     IMassMixingRatio, ZO3)
CALL gas%put_well_mixed(IO2, IVolumeMixingRatio, 0.20944_JPRB)

! Ensure the units of the gas mixing ratios are what is required by
! the gas absorption model
call set_gas_units(rad_config, gas)

! Call radiation scheme
CALL radiation(KLON, KLEV, KIDIA, KFDIA, rad_config, &
     &  single_level, thermodynamics, gas, cloud, aerosol, flux)

! Compute required output fluxes
! First the net fluxes
PFLUX_SW(KIDIA:KFDIA,:) = flux%sw_dn(KIDIA:KFDIA,:) - flux%sw_up(KIDIA:KFDIA,:)
PFLUX_LW(KIDIA:KFDIA,:) = flux%lw_dn(KIDIA:KFDIA,:) - flux%lw_up(KIDIA:KFDIA,:)
PFLUX_SW_CLEAR(KIDIA:KFDIA,:) &
     &  = flux%sw_dn_clear(KIDIA:KFDIA,:) - flux%sw_up_clear(KIDIA:KFDIA,:)
PFLUX_LW_CLEAR(KIDIA:KFDIA,:) &
     &  = flux%lw_dn_clear(KIDIA:KFDIA,:) - flux%lw_up_clear(KIDIA:KFDIA,:)
! Now the surface fluxes
PFLUX_SW_DN      (KIDIA:KFDIA) = flux%sw_dn             (KIDIA:KFDIA,KLEV+1)
PFLUX_LW_DN      (KIDIA:KFDIA) = flux%lw_dn             (KIDIA:KFDIA,KLEV+1)
PFLUX_SW_DN_CLEAR(KIDIA:KFDIA) = flux%sw_dn_clear       (KIDIA:KFDIA,KLEV+1)
PFLUX_LW_DN_CLEAR(KIDIA:KFDIA) = flux%lw_dn_clear       (KIDIA:KFDIA,KLEV+1)
PFLUX_DIR        (KIDIA:KFDIA) = flux%sw_dn_direct      (KIDIA:KFDIA,KLEV+1)
PFLUX_DIR_CLEAR  (KIDIA:KFDIA) = flux%sw_dn_direct_clear(KIDIA:KFDIA,KLEV+1)
PFLUX_DIR_INTO_SUN(KIDIA:KFDIA) = 0.0_JPRB
WHERE (PMU0(KIDIA:KFDIA) > EPSILON(1.0_JPRB))
  PFLUX_DIR_INTO_SUN(KIDIA:KFDIA) = PFLUX_DIR(KIDIA:KFDIA) / PMU0(KIDIA:KFDIA)
END WHERE
! Top-of-atmosphere downwelling flux
PFLUX_SW_DN_TOA(KIDIA:KFDIA) = flux%sw_dn(KIDIA:KFDIA,1)

! Compute UV fluxes as weighted sum of appropriate shortwave bands
PFLUX_UV       (KIDIA:KFDIA) = 0.0_JPRB
DO JBAND = 1,NWEIGHT_UV
  PFLUX_UV(KIDIA:KFDIA) = PFLUX_UV(KIDIA:KFDIA) + WEIGHT_UV(JBAND) &
       &  * flux%sw_dn_surf_band(IBAND_UV(JBAND),KIDIA:KFDIA)
ENDDO

! Compute photosynthetically active radiation similarly
PFLUX_PAR      (KIDIA:KFDIA) = 0.0_JPRB
PFLUX_PAR_CLEAR(KIDIA:KFDIA) = 0.0_JPRB
DO JBAND = 1,NWEIGHT_PAR
  PFLUX_PAR(KIDIA:KFDIA) = PFLUX_PAR(KIDIA:KFDIA) + WEIGHT_PAR(JBAND) &
       &  * flux%sw_dn_surf_band(IBAND_PAR(JBAND),KIDIA:KFDIA)
  PFLUX_PAR_CLEAR(KIDIA:KFDIA) = PFLUX_PAR_CLEAR(KIDIA:KFDIA) &
       &  + WEIGHT_PAR(JBAND) &
       &  * flux%sw_dn_surf_clear_band(IBAND_PAR(JBAND),KIDIA:KFDIA)
ENDDO

! Compute effective broadband emissivity
ZBLACK_BODY_NET_LW = PFLUX_LW_DN(KIDIA:KFDIA) &
     &  - RSIGMA*PTEMPERATURE_SKIN(KIDIA:KFDIA)**4
PEMIS_OUT(KIDIA:KFDIA) = PEMIS(KIDIA:KFDIA)
WHERE (ABS(ZBLACK_BODY_NET_LW) > 1.0E-5) 
  PEMIS_OUT(KIDIA:KFDIA) = PFLUX_LW(KIDIA:KFDIA,KLEV+1) / ZBLACK_BODY_NET_LW
END WHERE

! Copy longwave derivatives
IF (YRERAD%LAPPROXLWUPDATE) THEN
  PLWDERIVATIVE(KIDIA:KFDIA,:) = flux%lw_derivatives(KIDIA:KFDIA,:)
END IF

! Store the shortwave downwelling fluxes in each albedo band
IF (YRERAD%LAPPROXSWUPDATE) THEN
  PSWDIFFUSEBAND(KIDIA:KFDIA,:) = 0.0_JPRB
  PSWDIRECTBAND (KIDIA:KFDIA,:) = 0.0_JPRB
  DO JBAND = 1,rad_config%n_bands_sw
    JB_ALBEDO = rad_config%i_albedo_from_band_sw(JBAND)
    DO JLON = KIDIA,KFDIA
      PSWDIFFUSEBAND(JLON,JB_ALBEDO) = PSWDIFFUSEBAND(JLON,JB_ALBEDO) &
           &  + flux%sw_dn_surf_band(JBAND,JLON) &
           &  - flux%sw_dn_direct_surf_band(JBAND,JLON)
      PSWDIRECTBAND(JLON,JB_ALBEDO)  = PSWDIRECTBAND(JLON,JB_ALBEDO) &
           &  + flux%sw_dn_direct_surf_band(JBAND,JLON)
    ENDDO
  ENDDO
ENDIF

CALL single_level%deallocate
CALL thermodynamics%deallocate
CALL gas%deallocate
CALL cloud%deallocate
CALL aerosol%deallocate
CALL flux%deallocate

IF (LHOOK) CALL DR_HOOK('RADIATION_SCHEME',1,ZHOOK_HANDLE)

END SUBROUTINE RADIATION_SCHEME
