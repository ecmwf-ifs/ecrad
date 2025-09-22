SUBROUTINE SRTM_TAUMOL26 &
 & ( KIDIA   , KFDIA    , KLEV,&
 & P_COLMOL  ,K_LAYTROP,&
 & P_SFLUXZEN, P_TAUG   , P_TAUR    , PRMU0,   &
 & laytrop_min, laytrop_max)

!     Written by Eli J. Mlawer, Atmospheric & Environmental Research.

!     BAND 26:  22650-29000 cm-1 (low - nothing; high - nothing)

!      PARAMETER (MG=16, MXLAY=203, NBANDS=14)

! Modifications
!        M.Hamrud      01-Oct-2003 CY28 Cleaning

!     JJMorcrette 2003-02-24 adapted to ECMWF environment
!        D.Salmond  31-Oct-2007 Vector version in the style of RRTM from Meteo France & NEC
!     JJMorcrette 20110610 Flexible configuration for number of g-points

USE PARKIND1 , ONLY : JPIM, JPRB
USE PARSRTM  , ONLY : JPG
USE YOESRTM  , ONLY : NG26
USE YOESRTA26, ONLY : SFLUXREFC, RAYLC

IMPLICIT NONE

!-- Output
INTEGER(KIND=JPIM),INTENT(IN)    :: KIDIA, KFDIA
INTEGER(KIND=JPIM),INTENT(IN)    :: KLEV
REAL(KIND=JPRB)   ,INTENT(IN)    :: P_COLMOL(KIDIA:KFDIA,KLEV)
INTEGER(KIND=JPIM),INTENT(IN)    :: K_LAYTROP(KIDIA:KFDIA)

REAL(KIND=JPRB)   ,INTENT(INOUT) :: P_SFLUXZEN(KIDIA:KFDIA,JPG)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: P_TAUG(KIDIA:KFDIA,KLEV,JPG)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: P_TAUR(KIDIA:KFDIA,KLEV,JPG)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PRMU0(KIDIA:KFDIA)
INTEGER(KIND=JPIM), OPTIONAL, INTENT(INOUT) :: laytrop_min, laytrop_max
!- from AER
!- from INTFAC
!- from INTIND
!- from PRECISE
!- from PROFDATA
!- from SELF
INTEGER(KIND=JPIM) :: IG, I_LAY, I_LAYSOLFR(KIDIA:KFDIA), I_NLAYERS, IPLON
INTEGER(KIND=JPIM) :: llaytrop_min, llaytrop_max

#include "rrtm_utils.intfb.h"

    if (present(laytrop_min) .AND. present(laytrop_max)) then
       llaytrop_min = laytrop_min
       llaytrop_max = laytrop_max
    else
       CALL COMPUTE_LAYTROP_MIN_MAX(KIDIA, KFDIA, K_LAYTROP, llaytrop_min, llaytrop_max)
    endif

    !$ACC DATA CREATE(i_laysolfr) &
    !$ACC     PRESENT(p_colmol, k_laytrop, p_sfluxzen, p_taug, p_taur, prmu0)

    i_nlayers = klev
    !$ACC PARALLEL DEFAULT(NONE) ASYNC(1)
    !$ACC LOOP GANG VECTOR
    DO iplon = KIDIA, KFDIA
      i_laysolfr(iplon) = k_laytrop(iplon)
    ENDDO
    !$ACC END PARALLEL

    !$ACC WAIT
    !$ACC PARALLEL DEFAULT(NONE) ASYNC(1)
    !$ACC LOOP GANG VECTOR COLLAPSE(3)
    DO i_lay = 1, llaytrop_min
       DO iplon = KIDIA, KFDIA
!$NEC unroll(NG26)
         DO ig = 1 , ng26
           IF(i_lay == i_laysolfr(iplon)) p_sfluxzen(iplon,ig)=sfluxrefc(ig)
           p_taug(iplon,i_lay,ig) = 0.0_JPRB
           p_taur(iplon,i_lay,ig) = p_colmol(iplon,i_lay) * raylc(ig)
         ENDDO
       ENDDO
    ENDDO
    !$ACC END PARALLEL

    !$ACC PARALLEL DEFAULT(NONE) ASYNC(1)
    !$ACC LOOP GANG VECTOR COLLAPSE(2)
    DO i_lay = llaytrop_min+1, llaytrop_max
       DO iplon = KIDIA, KFDIA
          IF (i_lay <= k_laytrop(iplon)) THEN
!$NEC unroll(NG26)
            !$ACC LOOP SEQ
            DO ig = 1 , ng26
              IF(i_lay == i_laysolfr(iplon)) p_sfluxzen(iplon,ig)=sfluxrefc(ig)
              p_taug(iplon,i_lay,ig) = 0.0_JPRB
              p_taur(iplon,i_lay,ig) = p_colmol(iplon,i_lay) * raylc(ig)
            ENDDO
          ELSE
!$NEC unroll(NG26)
            !$ACC LOOP SEQ
            DO ig = 1 , ng26
              p_taug(iplon,i_lay,ig) = 0.0_JPRB
              p_taur(iplon,i_lay,ig) = p_colmol(iplon,i_lay) * raylc(ig)
            ENDDO
          ENDIF
       ENDDO
    ENDDO
    !$ACC END PARALLEL

    !$ACC PARALLEL DEFAULT(NONE) ASYNC(1)
    !$ACC LOOP GANG VECTOR COLLAPSE(3)
    DO ig = 1 , ng26
       DO i_lay = llaytrop_max+1, i_nlayers
         DO iplon = KIDIA, KFDIA
           p_taug(iplon,i_lay,ig) = 0.0_JPRB
           p_taur(iplon,i_lay,ig) = p_colmol(iplon,i_lay) * raylc(ig)
         ENDDO
       ENDDO
    ENDDO
    !$ACC END PARALLEL

    !$ACC WAIT
    !$ACC END DATA

END SUBROUTINE SRTM_TAUMOL26
