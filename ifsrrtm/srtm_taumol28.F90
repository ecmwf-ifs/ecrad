SUBROUTINE SRTM_TAUMOL28 &
 & ( KIDIA   , KFDIA    , KLEV,&
 & P_FAC00   , P_FAC01  , P_FAC10   , P_FAC11,&
 & K_JP      , K_JT     , K_JT1     , P_ONEMINUS,&
 & P_COLMOL  , P_COLO2  , P_COLO3,&
 & K_LAYTROP,&
 & P_SFLUXZEN, P_TAUG   , P_TAUR    , PRMU0,   &
 & laytrop_min, laytrop_max)

!     Written by Eli J. Mlawer, Atmospheric & Environmental Research.

!     BAND 28:  38000-50000 cm-1 (low - O3,O2; high - O3,O2)

! Modifications
!        M.Hamrud      01-Oct-2003 CY28 Cleaning

!     JJMorcrette 2003-02-24 adapted to ECMWF environment
!        D.Salmond  31-Oct-2007 Vector version in the style of RRTM from Meteo France & NEC
!     JJMorcrette 20010610 Flexible configuration for number of g-points

USE PARKIND1 , ONLY : JPIM, JPRB
USE PARSRTM  , ONLY : JPG
USE YOESRTM  , ONLY : NG28
USE YOESRTA28, ONLY : ABSA, ABSB, SFLUXREFC, RAYL, LAYREFFR, STRRAT
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
REAL(KIND=JPRB)   ,INTENT(IN)    :: P_COLMOL(KIDIA:KFDIA,KLEV)
REAL(KIND=JPRB)   ,INTENT(IN)    :: P_COLO2(KIDIA:KFDIA,KLEV)
REAL(KIND=JPRB)   ,INTENT(IN)    :: P_COLO3(KIDIA:KFDIA,KLEV)
INTEGER(KIND=JPIM),INTENT(IN)    :: K_LAYTROP(KIDIA:KFDIA)

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
INTEGER(KIND=JPIM) :: IG, IND0, IND1, JS, I_LAY, I_LAYSOLFR(KIDIA:KFDIA), I_NLAYERS, IPLON
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

    !$ACC DATA CREATE(i_laysolfr) &
    !$ACC     PRESENT(P_FAC00, P_FAC01, P_FAC10, P_FAC11, K_JP, K_JT, K_JT1, &
    !$ACC             P_ONEMINUS, P_COLMOL, P_COLO2, P_COLO3, K_LAYTROP, &
    !$ACC             P_SFLUXZEN, P_TAUG, P_TAUR, PRMU0)

    i_nlayers = klev
    !$ACC PARALLEL DEFAULT(NONE) ASYNC(1)
    !$ACC LOOP GANG(STATIC:1) VECTOR
    DO iplon = KIDIA,KFDIA
      i_laysolfr(iplon) = i_nlayers
    ENDDO
    !$ACC END PARALLEL

    !$ACC WAIT
    !$ACC PARALLEL DEFAULT(NONE) ASYNC(1)
    !$ACC LOOP GANG VECTOR COLLAPSE(2) PRIVATE(ind0, ind1, js, z_fs, z_speccomb, z_specmult, z_specparm, z_tauray)
    DO i_lay = 1, llaytrop_min
       DO iplon = KIDIA, KFDIA
         z_speccomb = p_colo3(iplon,i_lay) + strrat*p_colo2(iplon,i_lay)
         z_specparm = p_colo3(iplon,i_lay)/z_speccomb
         z_specparm = MIN(p_oneminus(iplon),z_specparm)
         z_specmult = 8._JPRB*(z_specparm)
         js = 1 + INT(z_specmult)
         z_fs = z_specmult - AINT(z_specmult)
         ind0 = ((k_jp(iplon,i_lay)-1)*5+(k_jt(iplon,i_lay)-1))*nspa(28) + js
         ind1 = (k_jp(iplon,i_lay)*5+(k_jt1(iplon,i_lay)-1))*nspa(28) + js
         z_tauray = p_colmol(iplon,i_lay) * rayl

         !$ACC LOOP SEQ
!$NEC unroll(NG28)
         DO ig = 1 , ng28
           p_taug(iplon,i_lay,ig) = z_speccomb * &
                & (&
                & (1._JPRB- z_fs) * ( absa(ind0,ig) * p_fac00(iplon,i_lay) +    &
                &                 absa(ind0+9,ig) * p_fac10(iplon,i_lay) +  &
                &                 absa(ind1,ig) * p_fac01(iplon,i_lay) +    &
                &                 absa(ind1+9,ig) * p_fac11(iplon,i_lay) )+ &
                & z_fs        * ( absa(ind0+1,ig) * p_fac00(iplon,i_lay) +  &
                &                 absa(ind0+10,ig) * p_fac10(iplon,i_lay) + &
                &                 absa(ind1+1,ig) * p_fac01(iplon,i_lay) +  &
                &                 absa(ind1+10,ig) * p_fac11(iplon,i_lay) ) &
                & )
           p_taur(iplon,i_lay,ig) = z_tauray
         ENDDO
       ENDDO
    ENDDO
    !$ACC END PARALLEL

    !$ACC PARALLEL DEFAULT(NONE) ASYNC(1)
    !$ACC LOOP GANG VECTOR COLLAPSE(2) PRIVATE(ind0, ind1, js, z_fs, z_speccomb, z_specmult, z_specparm, z_tauray)
    DO i_lay = llaytrop_min+1, llaytrop_max
       DO iplon = KIDIA, KFDIA
          IF (i_lay <= k_laytrop(iplon)) THEN
            z_speccomb = p_colo3(iplon,i_lay) + strrat*p_colo2(iplon,i_lay)
            z_specparm = p_colo3(iplon,i_lay)/z_speccomb
            z_specparm = MIN(p_oneminus(iplon),z_specparm)
            z_specmult = 8._JPRB*(z_specparm)
            js = 1 + INT(z_specmult)
            z_fs = z_specmult - AINT(z_specmult)
            ind0 = ((k_jp(iplon,i_lay)-1)*5+(k_jt(iplon,i_lay)-1))*nspa(28) + js
            ind1 = (k_jp(iplon,i_lay)*5+(k_jt1(iplon,i_lay)-1))*nspa(28) + js
            z_tauray = p_colmol(iplon,i_lay) * rayl

!$NEC unroll(NG28)
            !$ACC LOOP SEQ
            DO ig = 1 , ng28
              p_taug(iplon,i_lay,ig) = z_speccomb * &
                   & (&
                   & (1._JPRB- z_fs) * ( absa(ind0,ig) * p_fac00(iplon,i_lay) +    &
                   &                 absa(ind0+9,ig) * p_fac10(iplon,i_lay) +  &
                   &                 absa(ind1,ig) * p_fac01(iplon,i_lay) +    &
                   &                 absa(ind1+9,ig) * p_fac11(iplon,i_lay) )+ &
                   & z_fs        * ( absa(ind0+1,ig) * p_fac00(iplon,i_lay) +  &
                   &                 absa(ind0+10,ig) * p_fac10(iplon,i_lay) + &
                   &                 absa(ind1+1,ig) * p_fac01(iplon,i_lay) +  &
                   &                 absa(ind1+10,ig) * p_fac11(iplon,i_lay) ) &
                   & )
              p_taur(iplon,i_lay,ig) = z_tauray
            ENDDO
          ELSE
            IF (k_jp(iplon,i_lay-1) < layreffr &
                 &    .AND. k_jp(iplon,i_lay) >= layreffr) i_laysolfr(iplon) = i_lay
            z_speccomb = p_colo3(iplon,i_lay) + strrat*p_colo2(iplon,i_lay)
            z_specparm = p_colo3(iplon,i_lay)/z_speccomb
            z_specparm = MIN(p_oneminus(iplon),z_specparm)
            z_specmult = 4._JPRB*(z_specparm)
            js = 1 + INT(z_specmult)
            z_fs = z_specmult - AINT(z_specmult)
            ind0 = ((k_jp(iplon,i_lay)-13)*5+(k_jt(iplon,i_lay)-1))*nspb(28)+ js
            ind1 = ((k_jp(iplon,i_lay)-12)*5+(k_jt1(iplon,i_lay)-1))*nspb(28)+js
            z_tauray = p_colmol(iplon,i_lay) * rayl

!$NEC unroll(NG28)
            !$ACC LOOP SEQ
            DO ig = 1 , ng28
              p_taug(iplon,i_lay,ig) = z_speccomb * &
                   & (&
                   & (1._JPRB- z_fs) * ( absb(ind0,ig) * p_fac00(iplon,i_lay) +    &
                   &                 absb(ind0+5,ig) * p_fac10(iplon,i_lay) +  &
                   &                 absb(ind1,ig) * p_fac01(iplon,i_lay) +    &
                   &                 absb(ind1+5,ig) * p_fac11(iplon,i_lay) )+ &
                   & z_fs        * ( absb(ind0+1,ig) * p_fac00(iplon,i_lay) +  &
                   &                 absb(ind0+6,ig) * p_fac10(iplon,i_lay) +  &
                   &                 absb(ind1+1,ig) * p_fac01(iplon,i_lay) +  &
                   &                 absb(ind1+6,ig) * p_fac11(iplon,i_lay) )  &
                   & )
              IF (i_lay == i_laysolfr(iplon)) p_sfluxzen(iplon,ig) = sfluxrefc(ig,js) &
                   & + z_fs * (sfluxrefc(ig,js+1) - sfluxrefc(ig,js))
              p_taur(iplon,i_lay,ig) = z_tauray
            ENDDO
          ENDIF
       ENDDO
    ENDDO
    !$ACC END PARALLEL


    !$ACC PARALLEL DEFAULT(NONE) ASYNC(1)
    !$ACC LOOP GANG VECTOR COLLAPSE(2) PRIVATE(ind0, ind1, js, z_fs, z_speccomb, z_specmult, z_specparm, z_tauray)
    DO i_lay = llaytrop_max+1, i_nlayers
       DO iplon = KIDIA, KFDIA
         IF (k_jp(iplon,i_lay-1) < layreffr &
              &    .AND. k_jp(iplon,i_lay) >= layreffr) i_laysolfr(iplon) = i_lay
         z_speccomb = p_colo3(iplon,i_lay) + strrat*p_colo2(iplon,i_lay)
         z_specparm = p_colo3(iplon,i_lay)/z_speccomb
         z_specparm = MIN(p_oneminus(iplon),z_specparm)
         z_specmult = 4._JPRB*(z_specparm)
         js = 1 + INT(z_specmult)
         z_fs = z_specmult - AINT(z_specmult)
         ind0 = ((k_jp(iplon,i_lay)-13)*5+(k_jt(iplon,i_lay)-1))*nspb(28)+ js
         ind1 = ((k_jp(iplon,i_lay)-12)*5+(k_jt1(iplon,i_lay)-1))*nspb(28)+js
         z_tauray = p_colmol(iplon,i_lay) * rayl

         !$ACC LOOP SEQ
!$NEC unroll(NG28)
         DO ig = 1 , ng28
           p_taug(iplon,i_lay,ig) = z_speccomb * &
                & (&
                & (1._JPRB- z_fs) * ( absb(ind0,ig) * p_fac00(iplon,i_lay) +    &
                &                 absb(ind0+5,ig) * p_fac10(iplon,i_lay) +  &
                &                 absb(ind1,ig) * p_fac01(iplon,i_lay) +    &
                &                 absb(ind1+5,ig) * p_fac11(iplon,i_lay) )+ &
                & z_fs        * ( absb(ind0+1,ig) * p_fac00(iplon,i_lay) +  &
                &                 absb(ind0+6,ig) * p_fac10(iplon,i_lay) +  &
                &                 absb(ind1+1,ig) * p_fac01(iplon,i_lay) +  &
                &                 absb(ind1+6,ig) * p_fac11(iplon,i_lay) )  &
                & )
           IF (i_lay == i_laysolfr(iplon)) p_sfluxzen(iplon,ig) = sfluxrefc(ig,js) &
                & + z_fs * (sfluxrefc(ig,js+1) - sfluxrefc(ig,js))
           p_taur(iplon,i_lay,ig) = z_tauray
         ENDDO
       ENDDO
    ENDDO
    !$ACC END PARALLEL

    !$ACC WAIT
    !$ACC END DATA

END SUBROUTINE SRTM_TAUMOL28
