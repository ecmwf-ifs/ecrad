SUBROUTINE SRTM_TAUMOL18 &
 & ( KIDIA   , KFDIA    , KLEV,&
 & P_FAC00   , P_FAC01  , P_FAC10   , P_FAC11,&
 & K_JP      , K_JT     , K_JT1     , P_ONEMINUS,&
 & P_COLH2O  , P_COLCH4 , P_COLMOL,&
 & K_LAYTROP , P_SELFFAC, P_SELFFRAC, K_INDSELF  , P_FORFAC, P_FORFRAC, K_INDFOR,&
 & P_SFLUXZEN, P_TAUG   , P_TAUR    , PRMU0   &
 & )  

!     Written by Eli J. Mlawer, Atmospheric & Environmental Research.

!     BAND 18:  4000-4650 cm-1 (low - H2O,CH4; high - CH4)

! Modifications
!        M.Hamrud      01-Oct-2003 CY28 Cleaning

!     JJMorcrette 2003-02-24 adapted to ECMWF environment
!        D.Salmond  31-Oct-2007 Vector version in the style of RRTM from Meteo France & NEC
!     JJMorcrette 20110610 Flexible configuration for number of g-points

USE PARKIND1 , ONLY : JPIM, JPRB
USE YOMHOOK  , ONLY : LHOOK, DR_HOOK, JPHOOK
USE PARSRTM  , ONLY : JPG
USE YOESRTM  , ONLY : NG18
USE YOESRTA18, ONLY : ABSA, ABSB, FORREFC, SELFREFC, SFLUXREFC, RAYL, LAYREFFR, STRRAT  
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
REAL(KIND=JPRB)   ,INTENT(IN)    :: P_ONEMINUS(KIDIA:KFDIA) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: P_COLH2O(KIDIA:KFDIA,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: P_COLCH4(KIDIA:KFDIA,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: P_COLMOL(KIDIA:KFDIA,KLEV) 
INTEGER(KIND=JPIM),INTENT(IN)    :: K_LAYTROP(KIDIA:KFDIA) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: P_SELFFAC(KIDIA:KFDIA,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: P_SELFFRAC(KIDIA:KFDIA,KLEV) 
INTEGER(KIND=JPIM),INTENT(IN)    :: K_INDSELF(KIDIA:KFDIA,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: P_FORFAC(KIDIA:KFDIA,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: P_FORFRAC(KIDIA:KFDIA,KLEV) 
INTEGER(KIND=JPIM),INTENT(IN)    :: K_INDFOR(KIDIA:KFDIA,KLEV)

REAL(KIND=JPRB)   ,INTENT(OUT)   :: P_SFLUXZEN(JPG,KIDIA:KFDIA) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: P_TAUG(JPG,KIDIA:KFDIA,KLEV) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: P_TAUR(JPG,KIDIA:KFDIA,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PRMU0(KIDIA:KFDIA)
!- from INTFAC      
!- from INTIND
!- from PRECISE             
!- from PROFDATA             
!- from SELF             
INTEGER(KIND=JPIM) :: IG, IND0, IND1, INDS, INDF, JS, I_LAY, I_LAYSOLFR(KIDIA:KFDIA), I_NLAYERS, IPLON

INTEGER(KIND=JPIM) :: I_LAY_NEXT

REAL(KIND=JPRB) :: Z_FAC000, Z_FAC001, Z_FAC010, Z_FAC011, Z_FAC100, Z_FAC101,&
 & Z_FAC110, Z_FAC111, Z_FS, Z_SPECCOMB, Z_SPECMULT, Z_SPECPARM, &
 & Z_TAURAY  
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

IF (LHOOK) CALL DR_HOOK('SRTM_TAUMOL18',0,ZHOOK_HANDLE)
I_NLAYERS = KLEV

!     Compute the optical depth by interpolating in ln(pressure), 
!     temperature, and appropriate species.  Below LAYTROP, the water
!     vapor self-continuum is interpolated (in temperature) separately.  

I_LAYSOLFR(KIDIA:KFDIA) = K_LAYTROP(KIDIA:KFDIA)

DO I_LAY = 1, I_NLAYERS
  I_LAY_NEXT = MIN(I_NLAYERS, I_LAY+1)
  DO IPLON = KIDIA, KFDIA
    IF (PRMU0(IPLON) > 0.0_JPRB) THEN
      IF (I_LAY <= K_LAYTROP(IPLON)) THEN
        IF (K_JP(IPLON,I_LAY) < LAYREFFR .AND. K_JP(IPLON,I_LAY_NEXT) >= LAYREFFR) &
         & I_LAYSOLFR(IPLON) = MIN(I_LAY+1,K_LAYTROP(IPLON))  
        Z_SPECCOMB = P_COLH2O(IPLON,I_LAY) + STRRAT*P_COLCH4(IPLON,I_LAY)
        Z_SPECPARM = P_COLH2O(IPLON,I_LAY)/Z_SPECCOMB 
        IF (Z_SPECPARM >= P_ONEMINUS(IPLON)) Z_SPECPARM = P_ONEMINUS(IPLON)
        Z_SPECMULT = 8.*(Z_SPECPARM)
        JS = 1 + INT(Z_SPECMULT)
        Z_FS = MOD(Z_SPECMULT, 1.0_JPRB )
        ! Z_FAC000 = (1. - Z_FS) * P_FAC00(I_LAY)
        ! Z_FAC010 = (1. - Z_FS) * P_FAC10(I_LAY)
        ! Z_FAC100 = Z_FS * P_FAC00(I_LAY)
        ! Z_FAC110 = Z_FS * P_FAC10(I_LAY)
        ! Z_FAC001 = (1. - Z_FS) * P_FAC01(I_LAY)
        ! Z_FAC011 = (1. - Z_FS) * P_FAC11(I_LAY)
        ! Z_FAC101 = Z_FS * P_FAC01(I_LAY)
        ! Z_FAC111 = Z_FS * P_FAC11(I_LAY)
        IND0 = ((K_JP(IPLON,I_LAY)-1)*5+(K_JT(IPLON,I_LAY)-1))*NSPA(18) + JS
        IND1 = (K_JP(IPLON,I_LAY)*5+(K_JT1(IPLON,I_LAY)-1))*NSPA(18) + JS
        INDS = K_INDSELF(IPLON,I_LAY)
        INDF = K_INDFOR(IPLON,I_LAY)
        Z_TAURAY = P_COLMOL(IPLON,I_LAY) * RAYL

        !  DO IG = 1, NG(18)
!CDIR UNROLL=NG18
        DO IG = 1, NG18
          P_TAUG(IG,IPLON,I_LAY) = Z_SPECCOMB * &
           !    & (Z_FAC000 * ABSA(IG,IND0) + &
           !    & Z_FAC100 * ABSA(IG,IND0+1) + &
           !    & Z_FAC010 * ABSA(IG,IND0+9) + &
           !    & Z_FAC110 * ABSA(IG,IND0+10) + &
           !    & Z_FAC001 * ABSA(IG,IND1) + &
           !    & Z_FAC101 * ABSA(IG,IND1+1) + &
           !    & Z_FAC011 * ABSA(IG,IND1+9) + &
           !    & Z_FAC111 * ABSA(IG,IND1+10)) + &
           & (&
           & (1. - Z_FS) * ( ABSA(IG,IND0) * P_FAC00(IPLON,I_LAY) + &
           &                 ABSA(IG,IND0+9) * P_FAC10(IPLON,I_LAY) + &
           &                 ABSA(IG,IND1) * P_FAC01(IPLON,I_LAY) + &
           &                 ABSA(IG,IND1+9) * P_FAC11(IPLON,I_LAY) ) + &
           & Z_FS        * ( ABSA(IG,IND0+1) * P_FAC00(IPLON,I_LAY) + &
           &                 ABSA(IG,IND0+10) * P_FAC10(IPLON,I_LAY) + &
           &                 ABSA(IG,IND1+1) * P_FAC01(IPLON,I_LAY) + &
           &                 ABSA(IG,IND1+10) * P_FAC11(IPLON,I_LAY) ) &
           & ) + &
           & P_COLH2O(IPLON,I_LAY) * &
           & (P_SELFFAC(IPLON,I_LAY) * (SELFREFC(IG,INDS) + &
           & P_SELFFRAC(IPLON,I_LAY) * &
           & (SELFREFC(IG,INDS+1) - SELFREFC(IG,INDS))) + &
           & P_FORFAC(IPLON,I_LAY) * (FORREFC(IG,INDF) + &
           & P_FORFRAC(IPLON,I_LAY) * &
           & (FORREFC(IG,INDF+1) - FORREFC(IG,INDF))))   
          !     &           + TAURAY
          !    SSA(LAY,IG) = TAURAY/TAUG(LAY,IG)
          IF (I_LAY == I_LAYSOLFR(IPLON)) P_SFLUXZEN(IG,IPLON) = SFLUXREFC(IG,JS)  &
           & + Z_FS * (SFLUXREFC(IG,JS+1) - SFLUXREFC(IG,JS))  
          P_TAUR(IG,IPLON,I_LAY) = Z_TAURAY
        ENDDO
      ENDIF
    ENDIF
  ENDDO
ENDDO

DO I_LAY = 1, I_NLAYERS
  DO IPLON = KIDIA, KFDIA
    IF (PRMU0(IPLON) > 0.0_JPRB) THEN
      IF (I_LAY >= K_LAYTROP(IPLON)+1) THEN
        IND0 = ((K_JP(IPLON,I_LAY)-13)*5+(K_JT(IPLON,I_LAY)-1))*NSPB(18) + 1
        IND1 = ((K_JP(IPLON,I_LAY)-12)*5+(K_JT1(IPLON,I_LAY)-1))*NSPB(18) + 1
        Z_TAURAY = P_COLMOL(IPLON,I_LAY) * RAYL

        !  DO IG = 1, NG(18)
!CDIR UNROLL=NG18
        DO IG = 1, NG18
          P_TAUG(IG,IPLON,I_LAY) = P_COLCH4(IPLON,I_LAY) * &
           & (P_FAC00(IPLON,I_LAY) * ABSB(IG,IND0) + &
           & P_FAC10(IPLON,I_LAY) * ABSB(IG,IND0+1) + &
           & P_FAC01(IPLON,I_LAY) * ABSB(IG,IND1) +       &
           & P_FAC11(IPLON,I_LAY) * ABSB(IG,IND1+1))   
          !     &           + TAURAY
          !    SSA(LAY,IG) = TAURAY/TAUG(LAY,IG)
          P_TAUR(IG,IPLON,I_LAY) = Z_TAURAY
        ENDDO
      ENDIF
    ENDIF
  ENDDO
ENDDO

!-----------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('SRTM_TAUMOL18',1,ZHOOK_HANDLE)

END SUBROUTINE SRTM_TAUMOL18
