SUBROUTINE SRTM_TAUMOL27 &
 & ( KIDIA   , KFDIA    , KLEV,&
 & P_FAC00   , P_FAC01  , P_FAC10   , P_FAC11,&
 & K_JP      , K_JT     , K_JT1,&
 & P_COLMOL  , P_COLO3,&
 & K_LAYTROP,&
 & P_SFLUXZEN, P_TAUG   , P_TAUR    , PRMU0   &
 & )  

!     Written by Eli J. Mlawer, Atmospheric & Environmental Research.

!     BAND 27:  29000-38000 cm-1 (low - O3; high - O3)

! Modifications
!        M.Hamrud      01-Oct-2003 CY28 Cleaning

!     JJMorcrette 2003-02-24 adapted to ECMWF environment
!        D.Salmond  31-Oct-2007 Vector version in the style of RRTM from Meteo France & NEC
!     JJMorcrette 20110610 Flexible configuration for number of g-points

USE PARKIND1 , ONLY : JPIM, JPRB
USE YOMHOOK  , ONLY : LHOOK, DR_HOOK, JPHOOK
USE PARSRTM  , ONLY : JPG
USE YOESRTM  , ONLY : NG27
USE YOESRTA27, ONLY : ABSA, ABSB, SFLUXREFC, RAYLC, LAYREFFR, SCALEKUR  
USE YOESRTWN , ONLY : NSPA, NSPB

IMPLICIT NONE

!-- Output
INTEGER(KIND=JPIM),INTENT(IN)    :: KIDIA, KFDIA 
INTEGER(KIND=JPIM),INTENT(IN)    :: KLEV 
REAL(KIND=JPRB)   ,INTENT(IN)    :: P_FAC00(KIDIA:KFDIA,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: P_FAC01(KIDIA:KFDIA,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: P_FAC10(KIDIA:KFDIA,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: P_FAC11(KIDIA:KFDIA,KLEV) 
INTEGER(KIND=JPIM),INTENT(IN)    :: K_JP(KIDIA:KFDIA,KLEV) 
INTEGER(KIND=JPIM),INTENT(IN)    :: K_JT(KIDIA:KFDIA,KLEV) 
INTEGER(KIND=JPIM),INTENT(IN)    :: K_JT1(KIDIA:KFDIA,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: P_COLMOL(KIDIA:KFDIA,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: P_COLO3(KIDIA:KFDIA,KLEV) 
INTEGER(KIND=JPIM),INTENT(IN)    :: K_LAYTROP(KIDIA:KFDIA) 

REAL(KIND=JPRB)   ,INTENT(OUT)   :: P_SFLUXZEN(KIDIA:KFDIA,JPG) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: P_TAUG(KIDIA:KFDIA,KLEV,JPG) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: P_TAUR(KIDIA:KFDIA,KLEV,JPG) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PRMU0(KIDIA:KFDIA)
!- from INTFAC      
!- from INTIND
!- from PRECISE             
!- from PROFDATA             
!- from SELF             
INTEGER(KIND=JPIM) :: IG, IND0, IND1, I_LAY, I_LAYSOLFR(KIDIA:KFDIA), I_NLAYERS, IPLON

REAL(KIND=JPRB) :: Z_TAURAY  
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

IF (LHOOK) CALL DR_HOOK('SRTM_TAUMOL27',0,ZHOOK_HANDLE)

I_NLAYERS = KLEV

!     Compute the optical depth by interpolating in ln(pressure), 
!     temperature, and appropriate species.  Below LAYTROP, the water
!     vapor self-continuum is interpolated (in temperature) separately.  
DO I_LAY = 1, I_NLAYERS
  DO IPLON = KIDIA, KFDIA
    IF (PRMU0(IPLON) > 0.0_JPRB) THEN
      IF (I_LAY <= K_LAYTROP(IPLON)) THEN
        IND0 = ((K_JP(IPLON,I_LAY)-1)*5+(K_JT(IPLON,I_LAY)-1))*NSPA(27) + 1
        IND1 = (K_JP(IPLON,I_LAY)*5+(K_JT1(IPLON,I_LAY)-1))*NSPA(27) + 1

        !  DO IG = 1, NG(27)
!CDIR UNROLL=NG27
        DO IG = 1 , NG27
          Z_TAURAY = P_COLMOL(IPLON,I_LAY) * RAYLC(IG)
          P_TAUG(IPLON,I_LAY,IG) = P_COLO3(IPLON,I_LAY) * &
           & (P_FAC00(IPLON,I_LAY) * ABSA(IND0,IG) + &
           & P_FAC10(IPLON,I_LAY) * ABSA(IND0+1,IG) + &
           & P_FAC01(IPLON,I_LAY) * ABSA(IND1,IG) + &
           & P_FAC11(IPLON,I_LAY) * ABSA(IND1+1,IG))  
          !     &          + TAURAY
          !    SSA(LAY,IG) = TAURAY/TAUG(LAY,IG)
          P_TAUR(IPLON,I_LAY,IG) = Z_TAURAY
        ENDDO
      ENDIF
    ENDIF
  ENDDO
ENDDO

I_LAYSOLFR(:) = I_NLAYERS

DO I_LAY = 1, I_NLAYERS
  DO IPLON = KIDIA, KFDIA
    IF (PRMU0(IPLON) > 0.0_JPRB) THEN
      IF (I_LAY >= K_LAYTROP(IPLON)+1) THEN
        IF (K_JP(IPLON,I_LAY-1) < LAYREFFR .AND. K_JP(IPLON,I_LAY) >= LAYREFFR) &
         & I_LAYSOLFR(IPLON) = I_LAY  
        IND0 = ((K_JP(IPLON,I_LAY)-13)*5+(K_JT(IPLON,I_LAY)-1))*NSPB(27) + 1
        IND1 = ((K_JP(IPLON,I_LAY)-12)*5+(K_JT1(IPLON,I_LAY)-1))*NSPB(27) + 1

        !  DO IG = 1, NG(27)
!CDIR UNROLL=NG27
        DO IG = 1 , NG27
          Z_TAURAY = P_COLMOL(IPLON,I_LAY) * RAYLC(IG)
          P_TAUG(IPLON,I_LAY,IG) = P_COLO3(IPLON,I_LAY) * &
           & (P_FAC00(IPLON,I_LAY) * ABSB(IND0,IG) + &
           & P_FAC10(IPLON,I_LAY) * ABSB(IND0+1,IG) + &
           & P_FAC01(IPLON,I_LAY) * ABSB(IND1,IG) +  &
           & P_FAC11(IPLON,I_LAY) * ABSB(IND1+1,IG))  
          !     &          + TAURAY
          !    SSA(LAY,IG) = TAURAY/TAUG(LAY,IG)
          IF (I_LAY == I_LAYSOLFR(IPLON)) P_SFLUXZEN(IPLON,IG) = SCALEKUR * SFLUXREFC(IG) 
          P_TAUR(IPLON,I_LAY,IG) = Z_TAURAY
        ENDDO
      ENDIF
    ENDIF
  ENDDO
ENDDO

!-----------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('SRTM_TAUMOL27',1,ZHOOK_HANDLE)

END SUBROUTINE SRTM_TAUMOL27
