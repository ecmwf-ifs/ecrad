SUBROUTINE SRTM_TAUMOL27 &
 & ( KIDIA   , KFDIA    , KLEV,&
 & P_FAC00   , P_FAC01  , P_FAC10   , P_FAC11,&
 & K_JP      , K_JT     , K_JT1,&
 & P_COLMOL  , P_COLO3,&
 & K_LAYTROP,&
 & P_SFLUXZEN, P_TAUG   , P_TAUR    , PRMU0,   &
 & laytrop_min, laytrop_max)

!     Written by Eli J. Mlawer, Atmospheric & Environmental Research.

!     BAND 27:  29000-38000 cm-1 (low - O3; high - O3)

! Modifications
!        M.Hamrud      01-Oct-2003 CY28 Cleaning

!     JJMorcrette 2003-02-24 adapted to ECMWF environment
!        D.Salmond  31-Oct-2007 Vector version in the style of RRTM from Meteo France & NEC
!     JJMorcrette 20110610 Flexible configuration for number of g-points

USE PARKIND1 , ONLY : JPIM, JPRB
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

REAL(KIND=JPRB)   ,INTENT(INOUT) :: P_SFLUXZEN(KIDIA:KFDIA,JPG)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: P_TAUG(KIDIA:KFDIA,KLEV,JPG)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: P_TAUR(KIDIA:KFDIA,KLEV,JPG)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PRMU0(KIDIA:KFDIA)
INTEGER(KIND=JPIM),INTENT(IN)    :: laytrop_min, laytrop_max
!- from INTFAC
!- from INTIND
!- from PRECISE
!- from PROFDATA
!- from SELF
INTEGER(KIND=JPIM) :: IG, IND0, IND1, I_LAY, I_LAYSOLFR(KIDIA:KFDIA), I_NLAYERS, IPLON
REAL(KIND=JPRB) ::  &
 & Z_TAURAY

    !$ACC DATA CREATE(i_laysolfr) &
    !$ACC   PRESENT(P_FAC00, P_FAC01, P_FAC10, P_FAC11, K_JP, K_JT, K_JT1, &
    !$ACC           P_COLMOL, P_COLO3, K_LAYTROP, P_SFLUXZEN, P_TAUG, P_TAUR, &
    !$ACC           PRMU0)

    i_nlayers = klev
    !$ACC PARALLEL DEFAULT(NONE) ASYNC(1)
    !$ACC LOOP GANG(STATIC:1) VECTOR
    DO iplon = KIDIA,KFDIA
      i_laysolfr(iplon) = i_nlayers
    ENDDO
    !$ACC END PARALLEL

    !$ACC WAIT
    !$ACC PARALLEL DEFAULT(NONE) ASYNC(1)
    !$ACC LOOP GANG VECTOR COLLAPSE(2) PRIVATE(ind0, ind1)
    DO i_lay = 1, laytrop_min
       DO iplon = KIDIA, KFDIA
         ind0 = ((k_jp(iplon,i_lay)-1)*5+(k_jt(iplon,i_lay)-1))*nspa(27) + 1
         ind1 = (k_jp(iplon,i_lay)*5+(k_jt1(iplon,i_lay)-1))*nspa(27) + 1
         !$ACC LOOP SEQ PRIVATE(z_tauray)
!$NEC unroll(NG27)
         DO ig = 1 , ng27
           z_tauray = p_colmol(iplon,i_lay) * raylc(ig)
           p_taug(iplon,i_lay,ig) = p_colo3(iplon,i_lay) * &
                & (p_fac00(iplon,i_lay) * absa(ind0,ig)  + &
                & p_fac10(iplon,i_lay) * absa(ind0+1,ig) + &
                & p_fac01(iplon,i_lay) * absa(ind1,ig) +   &
                & p_fac11(iplon,i_lay) * absa(ind1+1,ig))
           p_taur(iplon,i_lay,ig) = z_tauray
         ENDDO
       ENDDO
    ENDDO
    !$ACC END PARALLEL

    !$ACC PARALLEL DEFAULT(NONE) ASYNC(1)
    !$ACC LOOP GANG VECTOR COLLAPSE(2) PRIVATE(ind0, ind1)
    DO i_lay = laytrop_min+1, laytrop_max
       DO iplon = KIDIA, KFDIA
          IF (i_lay <= k_laytrop(iplon)) THEN
            ind0 = ((k_jp(iplon,i_lay)-1)*5+(k_jt(iplon,i_lay)-1))*nspa(27) + 1
            ind1 = (k_jp(iplon,i_lay)*5+(k_jt1(iplon,i_lay)-1))*nspa(27) + 1
!$NEC unroll(NG27)
            !$ACC LOOP SEQ PRIVATE(z_tauray)
            DO ig = 1 , ng27
              z_tauray = p_colmol(iplon,i_lay) * raylc(ig)
              p_taug(iplon,i_lay,ig) = p_colo3(iplon,i_lay) * &
                   & (p_fac00(iplon,i_lay) * absa(ind0,ig)  + &
                   & p_fac10(iplon,i_lay) * absa(ind0+1,ig) + &
                   & p_fac01(iplon,i_lay) * absa(ind1,ig) +   &
                   & p_fac11(iplon,i_lay) * absa(ind1+1,ig))
              p_taur(iplon,i_lay,ig) = z_tauray
            ENDDO
          ELSE
            IF (k_jp(iplon,i_lay-1) < layreffr &
                 &    .AND. k_jp(iplon,i_lay) >= layreffr) i_laysolfr(iplon) = i_lay
            ind0 = ((k_jp(iplon,i_lay)-13)*5+(k_jt(iplon,i_lay)-1))*nspb(27) + 1
            ind1 = ((k_jp(iplon,i_lay)-12)*5+(k_jt1(iplon,i_lay)-1))*nspb(27)+ 1
!$NEC unroll(NG27)
            !$ACC LOOP SEQ PRIVATE(z_tauray)
            DO ig = 1 , ng27
              z_tauray = p_colmol(iplon,i_lay) * raylc(ig)
              p_taug(iplon,i_lay,ig) = p_colo3(iplon,i_lay) * &
                   & (p_fac00(iplon,i_lay) * absb(ind0,ig)  + &
                   & p_fac10(iplon,i_lay) * absb(ind0+1,ig) + &
                   & p_fac01(iplon,i_lay) * absb(ind1,ig) +   &
                   & p_fac11(iplon,i_lay) * absb(ind1+1,ig))
              IF (i_lay == i_laysolfr(iplon)) &
                   & p_sfluxzen(iplon,ig) = scalekur * sfluxrefc(ig)
              p_taur(iplon,i_lay,ig) = z_tauray
            ENDDO
          ENDIF
       ENDDO
    ENDDO
    !$ACC END PARALLEL

    !$ACC PARALLEL DEFAULT(NONE) ASYNC(1)
    !$ACC LOOP GANG VECTOR COLLAPSE(2) PRIVATE(ind0, ind1)
    DO i_lay = laytrop_max+1, i_nlayers
       DO iplon = KIDIA, KFDIA
         IF (k_jp(iplon,i_lay-1) < layreffr &
              &    .AND. k_jp(iplon,i_lay) >= layreffr) i_laysolfr(iplon) = i_lay
         ind0 = ((k_jp(iplon,i_lay)-13)*5+(k_jt(iplon,i_lay)-1))*nspb(27) + 1
         ind1 = ((k_jp(iplon,i_lay)-12)*5+(k_jt1(iplon,i_lay)-1))*nspb(27)+ 1
         !$ACC LOOP SEQ PRIVATE(z_tauray)
!$NEC unroll(NG27)
         DO ig = 1 , ng27
           z_tauray = p_colmol(iplon,i_lay) * raylc(ig)
           p_taug(iplon,i_lay,ig) = p_colo3(iplon,i_lay) * &
                & (p_fac00(iplon,i_lay) * absb(ind0,ig)  + &
                & p_fac10(iplon,i_lay) * absb(ind0+1,ig) + &
                & p_fac01(iplon,i_lay) * absb(ind1,ig) +   &
                & p_fac11(iplon,i_lay) * absb(ind1+1,ig))
           IF (i_lay == i_laysolfr(iplon)) &
                & p_sfluxzen(iplon,ig) = scalekur * sfluxrefc(ig)
           p_taur(iplon,i_lay,ig) = z_tauray
         ENDDO
       ENDDO
    ENDDO
    !$ACC END PARALLEL

    !$ACC WAIT
    !$ACC END DATA

END SUBROUTINE SRTM_TAUMOL27
