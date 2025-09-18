! (C) Copyright 2005- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
!
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!
SUBROUTINE LIQUID_EFFECTIVE_RADIUS &
     & (YDERAD, KIDIA, KFDIA, KLON, KLEV, &
     &  PPRESSURE, PTEMPERATURE, PCLOUD_FRAC, PQ_LIQ, PQ_RAIN, &
     &  PLAND_FRAC, PCCN_LAND, PCCN_SEA, &
     &  PRE_UM, LACC)

! LIQUID_EFFECTIVE_RADIUS
!
! PURPOSE
! -------
!   Calculate effective radius of liquid clouds
!
! AUTHOR
! ------
!   Robin Hogan, ECMWF (using code extracted from radlswr.F90)
!   Original: 2015-09-24
!
! MODIFICATIONS
! -------------
!
!   M Leutbecher (Oct 2020) SPP abstraction
!   L Descamps   (Jul 2023) 3 hard-coded values transformed in YDERAD parameters
!   for MF EPS Random parameters
! -------------------------------------------------------------------

USE PARKIND1 , ONLY : JPIM, JPRB
USE YOMHOOK  , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOERAD   , ONLY : TERAD
USE YOERDU   , ONLY : REPLOG, REPSCW
USE YOMLUN_ECRAD  , ONLY : NULERR
USE YOMCST   , ONLY : RD, RPI

! -------------------------------------------------------------------

IMPLICIT NONE

! INPUT ARGUMENTS

! *** Array dimensions and ranges
TYPE(TERAD)       ,INTENT(IN) :: YDERAD
INTEGER(KIND=JPIM),INTENT(IN) :: KIDIA    ! Start column to process
INTEGER(KIND=JPIM),INTENT(IN) :: KFDIA    ! End column to process
INTEGER(KIND=JPIM),INTENT(IN) :: KLON     ! Number of columns
INTEGER(KIND=JPIM),INTENT(IN) :: KLEV     ! Number of levels

! *** Variables on model levels
REAL(KIND=JPRB),   INTENT(IN) :: PPRESSURE(KLON,KLEV)    ! (Pa)
REAL(KIND=JPRB),   INTENT(IN) :: PTEMPERATURE(KLON,KLEV) ! (K)
REAL(KIND=JPRB),   INTENT(IN) :: PCLOUD_FRAC(KLON,KLEV)
REAL(KIND=JPRB),   INTENT(IN) :: PQ_LIQ(KLON,KLEV)       ! (kg/kg)
REAL(KIND=JPRB),   INTENT(IN) :: PQ_RAIN(KLON,KLEV)      ! (kg/kg)

! *** Single-level variables
REAL(KIND=JPRB),   INTENT(IN) :: PLAND_FRAC(KLON)        ! 1=land, 0=sea
REAL(KIND=JPRB),   INTENT(IN) :: PCCN_LAND(KLON)
REAL(KIND=JPRB),   INTENT(IN) :: PCCN_SEA(KLON)

! OUTPUT ARGUMENT
! Effective radius
REAL(KIND=JPRB),  INTENT(OUT) :: PRE_UM(KLON,KLEV) ! (microns)

! OPTIONAL INPUT ARGUMENT

LOGICAL,          INTENT(IN), OPTIONAL :: LACC

! PARAMETERS

! Minimum and maximum effective radius, in microns
REAL(KIND=JPRB), PARAMETER :: PP_MIN_RE_UM =  4.0_JPRB
REAL(KIND=JPRB), PARAMETER :: PP_MAX_RE_UM = 30.0_JPRB

! LOCAL VARIABLES
INTEGER(KIND=JPIM) :: IRADLP ! ID of effective radius scheme to use
REAL(KIND=JPRB) :: ZCCN    ! CCN concentration (units?)

REAL(KIND=JPRB) :: ZSPECTRAL_DISPERSION
REAL(KIND=JPRB) :: ZNTOT_CM3 ! Number conc in cm-3
REAL(KIND=JPRB) :: ZRE_CUBED
REAL(KIND=JPRB) :: ZLWC_GM3, ZRWC_GM3 ! In-cloud liquid, rain content in g m-3
REAL(KIND=JPRB) :: ZAIR_DENSITY_GM3   ! Air density in g m-3
REAL(KIND=JPRB) :: ZRAIN_RATIO        ! Ratio of rain to liquid water content
REAL(KIND=JPRB) :: ZWOOD_FACTOR, ZRATIO

INTEGER(KIND=JPIM) :: JL, JK

LOGICAL :: LLACC

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

! -------------------------------------------------------------------

#include "abor1.intfb.h"

! -------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('LIQUID_EFFECTIVE_RADIUS',0,ZHOOK_HANDLE)

! -------------------------------------------------------------------

! Liquid optical properties from namelist / suecrad.F90
IRADLP = YDERAD%NRADLP

LLACC = .FALSE.
IF (PRESENT(LACC)) LLACC = LACC

SELECT CASE(IRADLP)
CASE(0)
  ! Very old parameterization as a function of pressure, used in ERA-15
#if defined(_OPENACC)
  !$ACC KERNELS DEFAULT(NONE) PRESENT(PRE_UM, PPRESSURE) ASYNC(1) IF(LLACC)
  PRE_UM(KIDIA:KFDIA,:) = 10.0_JPRB&
       &  + (100000.0_JPRB-PPRESSURE(KIDIA:KFDIA,:))*3.5_JPRB
  !$ACC END KERNELS
#endif
#if defined(OMPGPU)
  !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO IF(LLACC)
  DO JL = KIDIA,KFDIA
     PRE_UM(JL,:) = 10.0_JPRB&
          &  + (100000.0_JPRB-PPRESSURE(JL,:))*3.5_JPRB
  ENDDO
  !$OMP END TARGET TEAMS DISTRIBUTE PARALLEL DO
#endif

CASE(1)
  ! Simple distinction between land (10um) and ocean (13um) by Zhang
  ! and Rossow
  !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO IF(LLACC)
  !$ACC PARALLEL DEFAULT(NONE) PRESENT(PRE_UM, PLAND_FRAC) ASYNC(1) IF(LLACC)
  !$ACC LOOP GANG VECTOR
  DO JL = KIDIA,KFDIA
    IF (PLAND_FRAC(JL) < 0.5_JPRB) THEN
      PRE_UM(JL,:) = 13.0_JPRB
    ELSE
      PRE_UM(JL,:) = 10.0_JPRB
    ENDIF
  ENDDO
  !$ACC END PARALLEL
  !$OMP END TARGET TEAMS DISTRIBUTE PARALLEL DO

CASE(2)
  ! Martin et al. (JAS 1994)
  !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO PRIVATE(ZCCN, ZSPECTRAL_DISPERSION, ZNTOT_CM3, ZRATIO, &
  !$OMP&    ZAIR_DENSITY_GM3, ZLWC_GM3, ZRWC_GM3, ZRAIN_RATIO, ZWOOD_FACTOR, ZRE_CUBED) IF(LLACC)
  !$ACC PARALLEL DEFAULT(NONE) PRESENT(PLAND_FRAC, YDERAD, PCCN_SEA, PCCN_LAND, &
  !$ACC     PQ_LIQ, PQ_RAIN, PPRESSURE, PTEMPERATURE, PCLOUD_FRAC, PRE_UM) ASYNC(1) IF(LLACC)
  !$ACC LOOP GANG VECTOR PRIVATE(ZCCN, ZSPECTRAL_DISPERSION, ZNTOT_CM3, ZRATIO, &
  !$ACC     ZAIR_DENSITY_GM3, ZLWC_GM3, ZRWC_GM3, ZRAIN_RATIO, ZWOOD_FACTOR, ZRE_CUBED)
  DO JL = KIDIA,KFDIA
    ! First compute the cloud droplet concentration
    IF (PLAND_FRAC(JL) < 0.5_JPRB) THEN
      ! Sea case
      IF (YDERAD%LCCNO) THEN
        ZCCN = PCCN_SEA(JL)
      ELSE
        ZCCN = YDERAD%RCCNSEA
      ENDIF
     ZSPECTRAL_DISPERSION = 0.77_JPRB
      ! ZSPECTRAL_DISPERSION = YDERAD%RRATSEA
      ! Cloud droplet concentration in cm-3 (activated CCN) over
      ! ocean
      ZNTOT_CM3 = -1.15E-03_JPRB*ZCCN*ZCCN + 0.963_JPRB*ZCCN + 5.30_JPRB
    ELSE
      ! Land case
      IF (YDERAD%LCCNL) THEN
        ZCCN=PCCN_LAND(JL)
      ELSE
        ZCCN=YDERAD%RCCNLND
      ENDIF
     ZSPECTRAL_DISPERSION = 0.69_JPRB
      ! ZSPECTRAL_DISPERSION = YDERAD%RRATLAND
      ! Cloud droplet concentration in cm-3 (activated CCN) over
      ! land
      ZNTOT_CM3 = -2.10E-04_JPRB*ZCCN*ZCCN + 0.568_JPRB*ZCCN - 27.9_JPRB
    ENDIF

   ZRATIO = (0.222_JPRB/ZSPECTRAL_DISPERSION)**0.333_JPRB
    ! ZRATIO = (YDERAD%RRATDRI/ZSPECTRAL_DISPERSION)**0.333_JPRB

    !$ACC LOOP SEQ
    DO JK = 1,KLEV

      ! Only consider cloudy regions
      IF (PCLOUD_FRAC(JL,JK) >= 0.001_JPRB&
           &  .AND. (PQ_LIQ(JL,JK)+PQ_RAIN(JL,JK)) > 0.0_JPRB) THEN

        ! Compute liquid and rain water contents
        ZAIR_DENSITY_GM3 = 1000.0_JPRB * PPRESSURE(JL,JK)&
             &           / (RD*PTEMPERATURE(JL,JK))
        ! In-cloud mean water contents found by dividing by cloud
        ! fraction
        ZLWC_GM3 = ZAIR_DENSITY_GM3 * PQ_LIQ(JL,JK)  / PCLOUD_FRAC(JL,JK)
        ZRWC_GM3 = ZAIR_DENSITY_GM3 * PQ_RAIN(JL,JK) / PCLOUD_FRAC(JL,JK)

        ! Wood's (2000, eq. 19) adjustment to Martin et al's
        ! parameterization
        IF (ZLWC_GM3 > REPSCW) THEN
          ZRAIN_RATIO = ZRWC_GM3 / ZLWC_GM3
          ZWOOD_FACTOR = ((1.0_JPRB + ZRAIN_RATIO)**0.666_JPRB)&
               &     / (1.0_JPRB + 0.2_JPRB * ZRATIO*ZRAIN_RATIO)
        ELSE
          ZWOOD_FACTOR = 1.0_JPRB
        ENDIF

        ! g m-3 and cm-3 units cancel out with density of water
        ! 10^6/(1000*1000); need a factor of 10^6 to convert to
        ! microns and cubed root is factor of 100 which appears in
        ! equation below
        ZRE_CUBED = (3.0_JPRB * (ZLWC_GM3 + ZRWC_GM3))&
             &    / (4.0_JPRB*RPI*ZNTOT_CM3*ZSPECTRAL_DISPERSION)
        IF (ZRE_CUBED > REPLOG) THEN
          PRE_UM(JL,JK) = ZWOOD_FACTOR*100.0_JPRB*EXP(0.333_JPRB*LOG(ZRE_CUBED))
          ! Make sure effective radius is bounded in range 4-30 microns
          PRE_UM(JL,JK) = MAX(PP_MIN_RE_UM, MIN(PRE_UM(JL,JK), PP_MAX_RE_UM))
        ELSE
          PRE_UM(JL,JK) = PP_MIN_RE_UM
        ENDIF

      ELSE
        ! Cloud fraction or liquid+rain water content too low to
        ! consider this a cloud
        PRE_UM(JL,JK) = PP_MIN_RE_UM

      ENDIF

    ENDDO

  ENDDO
  !$ACC END PARALLEL
  !$OMP END TARGET TEAMS DISTRIBUTE PARALLEL DO

CASE DEFAULT
  WRITE(NULERR,'(A,I0,A)') 'LIQUID EFFECTIVE RADIUS OPTION IRADLP=',IRADLP,' NOT AVAILABLE'
  CALL ABOR1('ERROR IN LIQUID_EFFECTIVE_RADIUS')
END SELECT

! -------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('LIQUID_EFFECTIVE_RADIUS',1,ZHOOK_HANDLE)

END SUBROUTINE LIQUID_EFFECTIVE_RADIUS
