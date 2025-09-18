SUBROUTINE SRTM_TAUMOL18 &
 & ( KIDIA   , KFDIA    , KLEV,&
 & P_FAC00   , P_FAC01  , P_FAC10   , P_FAC11,&
 & K_JP      , K_JT     , K_JT1     , P_ONEMINUS,&
 & P_COLH2O  , P_COLCH4 , P_COLMOL,&
 & K_LAYTROP , P_SELFFAC, P_SELFFRAC, K_INDSELF  , P_FORFAC, P_FORFRAC, K_INDFOR,&
 & P_SFLUXZEN, P_TAUG   , P_TAUR    , PRMU0,   &
 & laytrop_min, laytrop_max)

!     Written by Eli J. Mlawer, Atmospheric & Environmental Research.

!     BAND 18:  4000-4650 cm-1 (low - H2O,CH4; high - CH4)

! Modifications
!        M.Hamrud      01-Oct-2003 CY28 Cleaning

!     JJMorcrette 2003-02-24 adapted to ECMWF environment
!        D.Salmond  31-Oct-2007 Vector version in the style of RRTM from Meteo France & NEC
!     JJMorcrette 20110610 Flexible configuration for number of g-points

USE PARKIND1 , ONLY : JPIM, JPRB
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

REAL(KIND=JPRB)   ,INTENT(INOUT) :: P_SFLUXZEN(KIDIA:KFDIA,JPG)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: P_TAUG(KIDIA:KFDIA,KLEV,JPG)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: P_TAUR(KIDIA:KFDIA,KLEV,JPG)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PRMU0(KIDIA:KFDIA)
INTEGER(KIND=JPIM), OPTIONAL, INTENT(INOUT) :: laytrop_min, laytrop_max
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
INTEGER(KIND=JPIM) :: llaytrop_min, llaytrop_max

#include "rrtm_utils.intfb.h"

    if (present(laytrop_min) .AND. present(laytrop_max)) then
       llaytrop_min = laytrop_min
       llaytrop_max = laytrop_max
    else
       CALL COMPUTE_LAYTROP_MIN_MAX(KIDIA, KFDIA, K_LAYTROP, llaytrop_min, llaytrop_max)
    endif

    !$ACC DATA CREATE(I_LAYSOLFR) &
    !$ACC     PRESENT(P_FAC00, P_FAC01, P_FAC10, P_FAC11, K_JP, K_JT, K_JT1, &
    !$ACC             P_ONEMINUS, P_COLH2O, P_COLCH4, P_COLMOL, K_LAYTROP, &
    !$ACC             P_SELFFAC, P_SELFFRAC, K_INDSELF, P_FORFAC, P_FORFRAC, &
    !$ACC             K_INDFOR, P_SFLUXZEN, P_TAUG, P_TAUR, PRMU0)
    !$OMP TARGET ENTER DATA MAP(ALLOC: I_LAYSOLFR)
    !$OMP TARGET DATA MAP(PRESENT, ALLOC: P_FAC00, P_FAC01, P_FAC10, P_FAC11, K_JP, K_JT, K_JT1, &
    !$OMP             P_ONEMINUS, P_COLH2O, P_COLCH4, P_COLMOL, K_LAYTROP, &
    !$OMP             P_SELFFAC, P_SELFFRAC, K_INDSELF, P_FORFAC, P_FORFRAC, &
    !$OMP             K_INDFOR, P_SFLUXZEN, P_TAUG, P_TAUR, PRMU0)

    i_nlayers = klev
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO
    !$ACC WAIT
    !$ACC PARALLEL DEFAULT(NONE) ASYNC(1)
    !$ACC LOOP GANG(STATIC:1) VECTOR
    DO iplon = KIDIA, KFDIA
      i_laysolfr(iplon) = k_laytrop(iplon)
    ENDDO
    !$OMP END TARGET TEAMS DISTRIBUTE PARALLEL DO


    !$ACC LOOP SEQ
    DO i_lay = 1, llaytrop_min
      !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO PRIVATE(ind0, ind1, inds, indf, js, z_fs, &
      !$OMP   z_speccomb, z_specmult, z_specparm, z_tauray)
      !$ACC LOOP GANG(STATIC:1) VECTOR PRIVATE(ind0, ind1, inds, indf, js, z_fs, &
      !$ACC   z_speccomb, z_specmult, z_specparm, z_tauray)
      DO iplon = KIDIA, KFDIA
         IF (k_jp(iplon,i_lay) < layreffr                            &
              &    .AND. k_jp(iplon,i_lay+1) >= layreffr)            &
              &    i_laysolfr(iplon) = MIN(i_lay+1,k_laytrop(iplon))
         z_speccomb = p_colh2o(iplon,i_lay) + strrat*p_colch4(iplon,i_lay)
         z_specparm = p_colh2o(iplon,i_lay)/z_speccomb
         z_specparm = MIN(p_oneminus(iplon),z_specparm)
         z_specmult = 8._JPRB*(z_specparm)
         js = 1 + INT(z_specmult)
         z_fs = z_specmult - AINT(z_specmult)
         ind0 = ((k_jp(iplon,i_lay)-1)*5+(k_jt(iplon,i_lay)-1))*nspa(18) + js
         ind1 = (k_jp(iplon,i_lay)*5+(k_jt1(iplon,i_lay)-1))*nspa(18) + js
         inds = k_indself(iplon,i_lay)
         indf = k_indfor(iplon,i_lay)
         z_tauray = p_colmol(iplon,i_lay) * rayl

         !$ACC LOOP SEQ
!$NEC unroll(NG18)
         DO ig = 1, ng18
           p_taug(iplon,i_lay,ig) = z_speccomb *                            &
                & (                                                         &
                & (1._JPRB- z_fs) * ( absa(ind0,ig) * p_fac00(iplon,i_lay) +    &
                &                 absa(ind0+9,ig) * p_fac10(iplon,i_lay) +  &
                &                 absa(ind1,ig) * p_fac01(iplon,i_lay) +    &
                &                 absa(ind1+9,ig) * p_fac11(iplon,i_lay) )+ &
                & z_fs        * ( absa(ind0+1,ig) * p_fac00(iplon,i_lay) +  &
                &                 absa(ind0+10,ig) * p_fac10(iplon,i_lay) + &
                &                 absa(ind1+1,ig) * p_fac01(iplon,i_lay) +  &
                &                 absa(ind1+10,ig) * p_fac11(iplon,i_lay) ) &
                & ) +                                                       &
                & p_colh2o(iplon,i_lay) *                                   &
                & (p_selffac(iplon,i_lay) * (selfrefc(inds,ig) +            &
                & p_selffrac(iplon,i_lay) *                                 &
                & (selfrefc(inds+1,ig) - selfrefc(inds,ig))) +              &
                & p_forfac(iplon,i_lay) * (forrefc(indf,ig) +               &
                & p_forfrac(iplon,i_lay) *                                  &
                & (forrefc(indf+1,ig) - forrefc(indf,ig))))
           IF (i_lay == i_laysolfr(iplon))                           &
                &    p_sfluxzen(iplon,ig) = sfluxrefc(ig,js)         &
                &    + z_fs * (sfluxrefc(ig,js+1) - sfluxrefc(ig,js))
           p_taur(iplon,i_lay,ig) = z_tauray
         ENDDO
      ENDDO
      !$OMP END TARGET TEAMS DISTRIBUTE PARALLEL DO
    ENDDO

    !$ACC LOOP SEQ
    DO i_lay = llaytrop_min+1, llaytrop_max
       !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO PRIVATE(ind0, ind1, inds, indf, js, z_fs, z_speccomb, &
       !$OMP   z_specmult, z_specparm, z_tauray)
       !$ACC LOOP GANG(STATIC:1) VECTOR PRIVATE(ind0, ind1, inds, indf, js, z_fs, z_speccomb, z_specmult, z_specparm, &
       !$ACC   z_tauray)
       DO iplon = KIDIA, KFDIA
          IF (i_lay <= k_laytrop(iplon)) THEN
            IF (k_jp(iplon,i_lay) < layreffr                            &
                 &    .AND. k_jp(iplon,i_lay+1) >= layreffr)            &
                 &    i_laysolfr(iplon) = MIN(i_lay+1,k_laytrop(iplon))
            z_speccomb = p_colh2o(iplon,i_lay) + strrat*p_colch4(iplon,i_lay)
            z_specparm = p_colh2o(iplon,i_lay)/z_speccomb
            z_specparm = MIN(p_oneminus(iplon),z_specparm)
            z_specmult = 8._JPRB*(z_specparm)
            js = 1 + INT(z_specmult)
            z_fs = z_specmult - AINT(z_specmult)
            ind0 = ((k_jp(iplon,i_lay)-1)*5+(k_jt(iplon,i_lay)-1))*nspa(18) + js
            ind1 = (k_jp(iplon,i_lay)*5+(k_jt1(iplon,i_lay)-1))*nspa(18) + js
            inds = k_indself(iplon,i_lay)
            indf = k_indfor(iplon,i_lay)
            z_tauray = p_colmol(iplon,i_lay) * rayl

!$NEC unroll(NG18)
            !$ACC LOOP SEQ
            DO ig = 1, ng18
              p_taug(iplon,i_lay,ig) = z_speccomb *                            &
                   & (                                                         &
                   & (1._JPRB- z_fs) * ( absa(ind0,ig) * p_fac00(iplon,i_lay) +    &
                   &                 absa(ind0+9,ig) * p_fac10(iplon,i_lay) +  &
                   &                 absa(ind1,ig) * p_fac01(iplon,i_lay) +    &
                   &                 absa(ind1+9,ig) * p_fac11(iplon,i_lay) )+ &
                   & z_fs        * ( absa(ind0+1,ig) * p_fac00(iplon,i_lay) +  &
                   &                 absa(ind0+10,ig) * p_fac10(iplon,i_lay) + &
                   &                 absa(ind1+1,ig) * p_fac01(iplon,i_lay) +  &
                   &                 absa(ind1+10,ig) * p_fac11(iplon,i_lay) ) &
                   & ) +                                                       &
                   & p_colh2o(iplon,i_lay) *                                   &
                   & (p_selffac(iplon,i_lay) * (selfrefc(inds,ig) +            &
                   & p_selffrac(iplon,i_lay) *                                 &
                   & (selfrefc(inds+1,ig) - selfrefc(inds,ig))) +              &
                   & p_forfac(iplon,i_lay) * (forrefc(indf,ig) +               &
                   & p_forfrac(iplon,i_lay) *                                  &
                   & (forrefc(indf+1,ig) - forrefc(indf,ig))))
              IF (i_lay == i_laysolfr(iplon))                           &
                   &    p_sfluxzen(iplon,ig) = sfluxrefc(ig,js)         &
                   &    + z_fs * (sfluxrefc(ig,js+1) - sfluxrefc(ig,js))
              p_taur(iplon,i_lay,ig) = z_tauray
            ENDDO
          ELSE
            ind0 = ((k_jp(iplon,i_lay)-13)*5+(k_jt(iplon,i_lay)-1))*nspb(18) + 1
            ind1 = ((k_jp(iplon,i_lay)-12)*5+(k_jt1(iplon,i_lay)-1))*nspb(18)+ 1
            z_tauray = p_colmol(iplon,i_lay) * rayl

!$NEC unroll(NG18)
            !$ACC LOOP SEQ
            DO ig = 1, ng18
              p_taug(iplon,i_lay,ig) = p_colch4(iplon,i_lay) * &
                   & (p_fac00(iplon,i_lay) * absb(ind0,ig)  +  &
                   & p_fac10(iplon,i_lay) * absb(ind0+1,ig) +  &
                   & p_fac01(iplon,i_lay) * absb(ind1,ig)   +  &
                   & p_fac11(iplon,i_lay) * absb(ind1+1,ig))
              p_taur(iplon,i_lay,ig) = z_tauray
            ENDDO
          ENDIF
       ENDDO
       !$OMP END TARGET TEAMS DISTRIBUTE PARALLEL DO
    ENDDO

    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO COLLAPSE(2) PRIVATE(ind0, ind1, z_tauray)
    !$ACC LOOP SEQ
    DO i_lay = llaytrop_max+1, i_nlayers
      !$ACC LOOP GANG(STATIC:1) VECTOR PRIVATE(ind0, ind1, z_tauray)
       DO iplon = KIDIA, KFDIA
         ind0 = ((k_jp(iplon,i_lay)-13)*5+(k_jt(iplon,i_lay)-1))*nspb(18) + 1
         ind1 = ((k_jp(iplon,i_lay)-12)*5+(k_jt1(iplon,i_lay)-1))*nspb(18)+ 1
         z_tauray = p_colmol(iplon,i_lay) * rayl

         !$ACC LOOP SEQ
!$NEC unroll(NG18)
         DO ig = 1, ng18
           p_taug(iplon,i_lay,ig) = p_colch4(iplon,i_lay) * &
                & (p_fac00(iplon,i_lay) * absb(ind0,ig)  +  &
                & p_fac10(iplon,i_lay) * absb(ind0+1,ig) +  &
                & p_fac01(iplon,i_lay) * absb(ind1,ig)   +  &
                & p_fac11(iplon,i_lay) * absb(ind1+1,ig))
           p_taur(iplon,i_lay,ig) = z_tauray
         ENDDO
       ENDDO
    ENDDO
    !$ACC END PARALLEL
    !$OMP END TARGET TEAMS DISTRIBUTE PARALLEL DO

  !$ACC WAIT
  !$ACC END DATA

  !$OMP TARGET EXIT DATA MAP(DELETE: I_LAYSOLFR)
  !$OMP END TARGET DATA

END SUBROUTINE SRTM_TAUMOL18
