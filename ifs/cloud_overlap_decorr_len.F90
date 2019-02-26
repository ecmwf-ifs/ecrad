SUBROUTINE CLOUD_OVERLAP_DECORR_LEN &
     & (KIDIA, KFDIA, KLON, PGEMU, NDECOLAT, &
     &  PDECORR_LEN_EDGES_KM, PDECORR_LEN_WATER_KM, PDECORR_LEN_RATIO)

! CLOUD_OVERLAP_DECORR_LEN
!
! PURPOSE
! -------
!   Calculate the cloud overlap decorrelation length as a function of
!   latitude for use in the radiation scheme
!
! INTERFACE
! ---------
!   CLOUD_OVERLAP_DECORR_LEN is called from RADLSWR and RADIATION_SCHEME
!
! AUTHOR
! ------
!   Robin Hogan, ECMWF (using code extracted from radlswr.F90)
!   Original: 2016-02-16
!
! MODIFICATIONS
! -------------
!
! -------------------------------------------------------------------

USE PARKIND1 , ONLY : JPIM, JPRB
USE YOMHOOK  , ONLY : LHOOK, DR_HOOK
USE YOMCST   , ONLY : RPI
USE YOECLD   , ONLY : YRECLD

! -------------------------------------------------------------------

IMPLICIT NONE

! INPUT ARGUMENTS

! *** Array dimensions and ranges
INTEGER(KIND=JPIM),INTENT(IN) :: KIDIA    ! Start column to process
INTEGER(KIND=JPIM),INTENT(IN) :: KFDIA    ! End column to process
INTEGER(KIND=JPIM),INTENT(IN) :: KLON     ! Number of columns

! *** Configuration variable controlling the overlap scheme
INTEGER(KIND=JPIM),INTENT(IN) :: NDECOLAT

! *** Single-level variables 
REAL(KIND=JPRB),   INTENT(IN) :: PGEMU(KLON) ! Sine of latitude

! OUTPUT ARGUMENTS

! *** Decorrelation lengths for cloud edges and cloud water content,
! *** in km
REAL(KIND=JPRB), INTENT(OUT)           :: PDECORR_LEN_EDGES_KM(KLON)
REAL(KIND=JPRB), INTENT(OUT), OPTIONAL :: PDECORR_LEN_WATER_KM(KLON)
  
! Ratio of water-content to cloud-edge decorrelation lengths
REAL(KIND=JPRB), INTENT(OUT), OPTIONAL :: PDECORR_LEN_RATIO

! LOCAL VARIABLES

REAL(KIND=JPRB) :: ZRADIANS_TO_DEGREES, ZABS_LAT_DEG, ZCOS_LAT

INTEGER(KIND=JPIM) :: JL

REAL(KIND=JPRB) :: ZHOOK_HANDLE

! -------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('CLOUD_OVERLAP_DECORR_LEN',0,ZHOOK_HANDLE)
  
! -------------------------------------------------------------------

IF (NDECOLAT == 0) THEN

  ! Decorrelation lengths are constant values
  PDECORR_LEN_EDGES_KM(KIDIA:KFDIA) = YRECLD%RDECORR_CF
  IF (PRESENT(PDECORR_LEN_WATER_KM)) THEN
    PDECORR_LEN_WATER_KM(KIDIA:KFDIA) = YRECLD%RDECORR_CW
  ENDIF
  IF (PRESENT(PDECORR_LEN_RATIO)) THEN
    PDECORR_LEN_RATIO = YRECLD%RDECORR_CW / YRECLD%RDECORR_CF
  ENDIF

ELSE

  ZRADIANS_TO_DEGREES = 180.0_JPRB / RPI

  IF (NDECOLAT == 1) THEN
    ! Shonk et al. (2010) Eq. 13 formula
    DO JL = KIDIA,KFDIA
      ZABS_LAT_DEG = ABS(ASIN(PGEMU(JL)) * ZRADIANS_TO_DEGREES)
      PDECORR_LEN_EDGES_KM(JL) = 2.899_JPRB - 0.02759_JPRB * ZABS_LAT_DEG
    ENDDO
  ELSE ! NDECOLAT == 2
    DO JL = KIDIA,KFDIA
      ! Shonk et al. (2010) but smoothed over the equator
      ZCOS_LAT = COS(ASIN(PGEMU(JL)))
      PDECORR_LEN_EDGES_KM(JL) = 0.75_JPRB + 2.149_JPRB * ZCOS_LAT*ZCOS_LAT
    ENDDO
  ENDIF

  ! Both NDECOLAT = 1 and 2 assume that the decorrelation length for
  ! cloud water content is half that for cloud edges
  IF (PRESENT(PDECORR_LEN_WATER_KM)) THEN
    PDECORR_LEN_WATER_KM(KIDIA:KFDIA) = PDECORR_LEN_EDGES_KM(KIDIA:KFDIA) * 0.5_JPRB
  ENDIF

  IF (PRESENT(PDECORR_LEN_RATIO)) THEN
    PDECORR_LEN_RATIO = 0.5_JPRB
  ENDIF

ENDIF

! -------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('CLOUD_OVERLAP_DECORR_LEN',1,ZHOOK_HANDLE)

END SUBROUTINE CLOUD_OVERLAP_DECORR_LEN
