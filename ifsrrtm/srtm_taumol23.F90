SUBROUTINE SRTM_TAUMOL23 &
 & ( KIDIA   , KFDIA    , KLEV,&
 & P_FAC00   , P_FAC01  , P_FAC10   , P_FAC11,&
 & K_JP      , K_JT     , K_JT1,&
 & P_COLH2O  , P_COLMOL,&
 & K_LAYTROP , P_SELFFAC, P_SELFFRAC, K_INDSELF  , P_FORFAC, P_FORFRAC, K_INDFOR,&
 & P_SFLUXZEN, P_TAUG   , P_TAUR    , PRMU0,   &
 & laytrop_min, laytrop_max)

!     Written by Eli J. Mlawer, Atmospheric & Environmental Research.

!     BAND 23:  8050-12850 cm-1 (low - H2O; high - nothing)

! Modifications
!        M.Hamrud      01-Oct-2003 CY28 Cleaning

!     JJMorcrette 2003-02-24 adapted to ECMWF environment
!        D.Salmond  31-Oct-2007 Vector version in the style of RRTM from Meteo France & NEC
!     JJMorcrette 20110610 Flexible configuration for number of g-points

USE PARKIND1 , ONLY : JPIM, JPRB
USE PARSRTM  , ONLY : JPG
USE YOESRTM  , ONLY : NG23
USE YOESRTA23, ONLY : ABSA, FORREFC, SELFREFC, SFLUXREFC, RAYLC, LAYREFFR, GIVFAC
USE YOESRTWN , ONLY : NSPA

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
REAL(KIND=JPRB)   ,INTENT(IN)    :: P_COLH2O(KIDIA:KFDIA,KLEV)
REAL(KIND=JPRB)   ,INTENT(IN)    :: P_COLMOL(KIDIA:KFDIA,KLEV)
INTEGER(KIND=JPIM),INTENT(IN)    :: K_LAYTROP(KIDIA:KFDIA)
REAL(KIND=JPRB)   ,INTENT(IN)    :: P_SELFFAC(KIDIA:KFDIA,KLEV)
REAL(KIND=JPRB)   ,INTENT(IN)    :: P_SELFFRAC(KIDIA:KFDIA,KLEV)
INTEGER(KIND=JPIM),INTENT(IN)    :: K_INDSELF(KIDIA:KFDIA,KLEV)
REAL(KIND=JPRB)   ,INTENT(IN)    :: P_FORFAC(KIDIA:KFDIA,KLEV)
REAL(KIND=JPRB)   ,INTENT(IN)    :: P_FORFRAC(KIDIA:KFDIA,KLEV)
INTEGER(KIND=JPIM),INTENT(IN)    :: K_INDFOR(KIDIA:KFDIA,KLEV)

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
INTEGER(KIND=JPIM) :: IG, IND0, IND1, INDS, INDF, I_LAY, I_LAYSOLFR(KIDIA:KFDIA), I_NLAYERS, IPLON
INTEGER(KIND=JPIM) :: I_LAY_NEXT

REAL(KIND=JPRB) ::  &
 & Z_TAURAY

    !$ACC DATA CREATE(i_laysolfr) &
    !$ACC     PRESENT(p_fac00, p_fac01, p_fac10, p_fac11, k_jp, k_jt, k_jt1, &
    !$ACC             p_colh2o, p_colmol, k_laytrop, p_selffac, p_selffrac, &
    !$ACC             k_indself, p_forfac, p_forfrac, k_indfor, p_sfluxzen, &
    !$ACC             p_taug, p_taur , prmu0)

    i_nlayers = klev
    !$ACC WAIT
    !$ACC PARALLEL DEFAULT(NONE) ASYNC(1)
    !$ACC LOOP GANG(STATIC:1) VECTOR
    DO iplon = KIDIA, KFDIA
      i_laysolfr(iplon) = k_laytrop(iplon)
    ENDDO

    !$ACC LOOP SEQ
    DO i_lay = 1, laytrop_min
      !$ACC LOOP GANG(STATIC:1) VECTOR PRIVATE(ind0, ind1, inds, indf)
       DO iplon = KIDIA, KFDIA
         IF (k_jp(iplon,i_lay) < layreffr                            &
              &    .AND. k_jp(iplon,i_lay+1) >= layreffr)            &
              &    i_laysolfr(iplon) = MIN(i_lay+1,k_laytrop(iplon))
         ind0 = ((k_jp(iplon,i_lay)-1)*5+(k_jt(iplon,i_lay)-1))*nspa(23) + 1
         ind1 = (k_jp(iplon,i_lay)*5+(k_jt1(iplon,i_lay)-1))*nspa(23) + 1
         inds = k_indself(iplon,i_lay)
         indf = k_indfor(iplon,i_lay)

         !$ACC LOOP SEQ PRIVATE(z_tauray)
!$NEC unroll(NG23)
         DO ig = 1 , ng23
           z_tauray = p_colmol(iplon,i_lay) * raylc(ig)
           p_taug(iplon,i_lay,ig) = p_colh2o(iplon,i_lay) *         &
                & (givfac * (p_fac00(iplon,i_lay) * absa(ind0,ig) + &
                & p_fac10(iplon,i_lay) * absa(ind0+1,ig) +          &
                & p_fac01(iplon,i_lay) * absa(ind1,ig) +            &
                & p_fac11(iplon,i_lay) * absa(ind1+1,ig)) +         &
                & p_selffac(iplon,i_lay) * (selfrefc(inds,ig) +     &
                & p_selffrac(iplon,i_lay) *                         &
                & (selfrefc(inds+1,ig) - selfrefc(inds,ig))) +      &
                & p_forfac(iplon,i_lay) * (forrefc(indf,ig) +       &
                & p_forfrac(iplon,i_lay) *                          &
                & (forrefc(indf+1,ig) - forrefc(indf,ig))))
           IF (i_lay == i_laysolfr(iplon)) &
                p_sfluxzen(iplon,ig) = sfluxrefc(ig)
           p_taur(iplon,i_lay,ig) = z_tauray
         ENDDO
       ENDDO
    ENDDO

    !$ACC LOOP SEQ
    DO i_lay = laytrop_min+1, laytrop_max
       !$ACC LOOP GANG(STATIC:1) VECTOR PRIVATE(ind0, ind1, inds, indf)
       DO iplon = KIDIA, KFDIA
          IF (i_lay <= k_laytrop(iplon)) THEN
            IF (k_jp(iplon,i_lay) < layreffr                            &
                 &    .AND. k_jp(iplon,i_lay+1) >= layreffr)            &
                 &    i_laysolfr(iplon) = MIN(i_lay+1,k_laytrop(iplon))
            ind0 = ((k_jp(iplon,i_lay)-1)*5+(k_jt(iplon,i_lay)-1))*nspa(23) + 1
            ind1 = (k_jp(iplon,i_lay)*5+(k_jt1(iplon,i_lay)-1))*nspa(23) + 1
            inds = k_indself(iplon,i_lay)
            indf = k_indfor(iplon,i_lay)

!$NEC unroll(NG23)
            !$ACC LOOP SEQ PRIVATE(z_tauray)
            DO ig = 1 , ng23
              z_tauray = p_colmol(iplon,i_lay) * raylc(ig)
              p_taug(iplon,i_lay,ig) = p_colh2o(iplon,i_lay) *         &
                   & (givfac * (p_fac00(iplon,i_lay) * absa(ind0,ig) + &
                   & p_fac10(iplon,i_lay) * absa(ind0+1,ig) +          &
                   & p_fac01(iplon,i_lay) * absa(ind1,ig) +            &
                   & p_fac11(iplon,i_lay) * absa(ind1+1,ig)) +         &
                   & p_selffac(iplon,i_lay) * (selfrefc(inds,ig) +     &
                   & p_selffrac(iplon,i_lay) *                         &
                   & (selfrefc(inds+1,ig) - selfrefc(inds,ig))) +      &
                   & p_forfac(iplon,i_lay) * (forrefc(indf,ig) +       &
                   & p_forfrac(iplon,i_lay) *                          &
                   & (forrefc(indf+1,ig) - forrefc(indf,ig))))
              IF (i_lay == i_laysolfr(iplon)) &
                   p_sfluxzen(iplon,ig) = sfluxrefc(ig)
              p_taur(iplon,i_lay,ig) = z_tauray
            ENDDO
          ELSE
!$NEC unroll(NG23)
            !$ACC LOOP SEQ
            DO ig = 1 , ng23
              p_taug(iplon,i_lay,ig) = 0.0_JPRB
              p_taur(iplon,i_lay,ig) = p_colmol(iplon,i_lay) * raylc(ig)
            ENDDO
          ENDIF
       ENDDO
    ENDDO

    !$ACC LOOP SEQ
    DO ig = 1 , ng23
      !$ACC LOOP SEQ
      DO i_lay = laytrop_max+1, i_nlayers
        !$ACC LOOP GANG(STATIC:1) VECTOR
        DO iplon = KIDIA, KFDIA
          p_taug(iplon,i_lay,ig) = 0.0_JPRB
          p_taur(iplon,i_lay,ig) = p_colmol(iplon,i_lay) * raylc(ig)
        ENDDO
      ENDDO
    ENDDO
    !$ACC END PARALLEL

    !$ACC WAIT
    !$ACC END DATA

END SUBROUTINE SRTM_TAUMOL23
