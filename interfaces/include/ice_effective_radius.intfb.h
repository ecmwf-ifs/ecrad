interface
SUBROUTINE ICE_EFFECTIVE_RADIUS &
 & (YDERAD,KIDIA, KFDIA, KLON, KLEV, &
 & PPRESSURE, PTEMPERATURE, PCLOUD_FRAC, PQ_ICE, PQ_SNOW, PGEMU, &
 & PRE_UM)
USE PARKIND1 , ONLY : JPIM, JPRB
USE YOERAD , ONLY : TERAD
TYPE(TERAD) ,INTENT(INOUT):: YDERAD
INTEGER(KIND=JPIM),INTENT(IN) :: KIDIA
INTEGER(KIND=JPIM),INTENT(IN) :: KFDIA
INTEGER(KIND=JPIM),INTENT(IN) :: KLON
INTEGER(KIND=JPIM),INTENT(IN) :: KLEV
REAL(KIND=JPRB), INTENT(IN) :: PPRESSURE(KLON,KLEV)
REAL(KIND=JPRB), INTENT(IN) :: PTEMPERATURE(KLON,KLEV)
REAL(KIND=JPRB), INTENT(IN) :: PCLOUD_FRAC(KLON,KLEV)
REAL(KIND=JPRB), INTENT(IN) :: PQ_ICE(KLON,KLEV)
REAL(KIND=JPRB), INTENT(IN) :: PQ_SNOW(KLON,KLEV)
REAL(KIND=JPRB), INTENT(IN) :: PGEMU(KLON)
REAL(KIND=JPRB), INTENT(OUT) :: PRE_UM(KLON,KLEV)
END SUBROUTINE ICE_EFFECTIVE_RADIUS
end interface