SUBROUTINE SRTM_TAUMOL25 &
 & ( KIDIA   , KFDIA    , KLEV,&
 & P_FAC00   , P_FAC01  , P_FAC10   , P_FAC11,&
 & K_JP      , K_JT     , K_JT1,&
 & P_COLH2O  , P_COLMOL , P_COLO3,&
 & K_LAYTROP,&
 & P_SFLUXZEN, P_TAUG   , P_TAUR    , PRMU0   &
 & )

!     Written by Eli J. Mlawer, Atmospheric & Environmental Research.

!     BAND 25:  16000-22650 cm-1 (low - H2O; high - nothing)

!      PARAMETER (MG=16, MXLAY=203, NBANDS=14)

! Modifications
!        M.Hamrud      01-Oct-2003 CY28 Cleaning

!     JJMorcrette 2003-02-24 adapted to ECMWF environment
!        D.Salmond  31-Oct-2007 Vector version in the style of RRTM from Meteo France & NEC
!     JJMorcrette 20110610 Flexible configuration for number of g-points

USE PARKIND1 , ONLY : JPIM, JPRB
USE PARSRTM  , ONLY : JPG
USE YOESRTM  , ONLY : NG25
USE YOESRTA25, ONLY : ABSA, SFLUXREFC, ABSO3AC, ABSO3BC, RAYLC, LAYREFFR
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
REAL(KIND=JPRB)   ,INTENT(IN)    :: P_COLO3(KIDIA:KFDIA,KLEV)
INTEGER(KIND=JPIM),INTENT(IN)    :: K_LAYTROP(KIDIA:KFDIA)

REAL(KIND=JPRB)   ,INTENT(INOUT) :: P_SFLUXZEN(KIDIA:KFDIA,JPG)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: P_TAUG(KIDIA:KFDIA,KLEV,JPG)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: P_TAUR(KIDIA:KFDIA,KLEV,JPG)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PRMU0(KIDIA:KFDIA)
!- from INTFAC
!- from INTIND
!- from PRECISE
!- from PROFDATA
!- from SELF
INTEGER(KIND=JPIM) :: IG, IND0, IND1, I_LAY, I_LAYSOLFR(KIDIA:KFDIA), I_NLAYERS, IPLON
INTEGER(KIND=JPIM) :: laytrop_min, laytrop_max
INTEGER(KIND=JPIM) :: I_LAY_NEXT

REAL(KIND=JPRB) ::  &
 & Z_TAURAY

    !$ACC DATA CREATE(I_LAYSOLFR) &
    !$ACC     PRESENT(P_FAC00, P_FAC01, P_FAC10, P_FAC11, K_JP, K_JT, K_JT1, &
    !$ACC             P_COLH2O, P_COLMOL, P_COLO3, K_LAYTROP, P_SFLUXZEN, &
    !$ACC             P_TAUG, P_TAUR, PRMU0)

#ifndef _OPENACC
    laytrop_min = MINVAL(k_laytrop(KIDIA:KFDIA))
    laytrop_max = MAXVAL(k_laytrop(KIDIA:KFDIA))
#else
    laytrop_min = HUGE(laytrop_min)
    laytrop_max = -HUGE(laytrop_max)
    !$ACC PARALLEL DEFAULT(NONE) ASYNC(1)
    !$ACC LOOP GANG VECTOR REDUCTION(min:laytrop_min) REDUCTION(max:laytrop_max)
    do iplon = KIDIA,KFDIA
      laytrop_min = MIN(laytrop_min, k_laytrop(iplon))
      laytrop_max = MAX(laytrop_max, k_laytrop(iplon))
    end do
    !$ACC END PARALLEL
#endif

    i_nlayers = klev
    !$ACC WAIT
    !$ACC PARALLEL DEFAULT(NONE) ASYNC(1)
    !$ACC LOOP GANG(STATIC:1) VECTOR
    DO iplon = KIDIA, KFDIA
      i_laysolfr(iplon) = k_laytrop(iplon)
    ENDDO

    !$ACC LOOP SEQ
    DO i_lay = 1, laytrop_min
      !$ACC LOOP GANG(STATIC:1) VECTOR PRIVATE(ind0, ind1)
       DO iplon = KIDIA, KFDIA
         IF (k_jp(iplon,i_lay) < layreffr .AND.   &
              &    k_jp(iplon,i_lay+1) >= layreffr) &
              &    i_laysolfr(iplon) = MIN(i_lay+1,k_laytrop(iplon))
         ind0 = ((k_jp(iplon,i_lay)-1)*5+(k_jt(iplon,i_lay)-1))*nspa(25) + 1
         ind1 = (k_jp(iplon,i_lay)*5+(k_jt1(iplon,i_lay)-1))*nspa(25) + 1
         !$ACC LOOP SEQ PRIVATE(z_tauray)
!$NEC unroll(NG25)
         DO ig = 1 , ng25
           z_tauray = p_colmol(iplon,i_lay) * raylc(ig)
           p_taug(iplon,i_lay,ig) = p_colh2o(iplon,i_lay) * &
                & (p_fac00(iplon,i_lay) * absa(ind0,ig)   + &
                & p_fac10(iplon,i_lay) * absa(ind0+1,ig)  + &
                & p_fac01(iplon,i_lay) * absa(ind1,ig)    + &
                & p_fac11(iplon,i_lay) * absa(ind1+1,ig)) + &
                & p_colo3(iplon,i_lay) * abso3ac(ig)
           IF(i_lay == i_laysolfr(iplon)) p_sfluxzen(iplon,ig)=sfluxrefc(ig)
           p_taur(iplon,i_lay,ig) = z_tauray
         ENDDO
       ENDDO
    ENDDO

    !$ACC LOOP SEQ
    DO i_lay = laytrop_min+1, laytrop_max
      !$ACC LOOP GANG(STATIC:1) VECTOR PRIVATE(ind0, ind1)
       DO iplon = KIDIA, KFDIA
          IF (i_lay <= k_laytrop(iplon)) THEN
            IF (k_jp(iplon,i_lay) < layreffr .AND.   &
                 &    k_jp(iplon,i_lay+1) >= layreffr) &
                 &    i_laysolfr(iplon) = MIN(i_lay+1,k_laytrop(iplon))
            ind0 = ((k_jp(iplon,i_lay)-1)*5+(k_jt(iplon,i_lay)-1))*nspa(25) + 1
            ind1 = (k_jp(iplon,i_lay)*5+(k_jt1(iplon,i_lay)-1))*nspa(25) + 1
!$NEC unroll(NG25)
            !$ACC LOOP SEQ PRIVATE(z_tauray)
            DO ig = 1 , ng25
              z_tauray = p_colmol(iplon,i_lay) * raylc(ig)
              p_taug(iplon,i_lay,ig) = p_colh2o(iplon,i_lay) * &
                   & (p_fac00(iplon,i_lay) * absa(ind0,ig)   + &
                   & p_fac10(iplon,i_lay) * absa(ind0+1,ig)  + &
                   & p_fac01(iplon,i_lay) * absa(ind1,ig)    + &
                   & p_fac11(iplon,i_lay) * absa(ind1+1,ig)) + &
                   & p_colo3(iplon,i_lay) * abso3ac(ig)
              IF(i_lay == i_laysolfr(iplon)) p_sfluxzen(iplon,ig)=sfluxrefc(ig)
              p_taur(iplon,i_lay,ig) = z_tauray
            ENDDO
          ELSE
!$NEC unroll(NG25)
            !$ACC LOOP SEQ PRIVATE(z_tauray)
            DO ig = 1 , ng25
              z_tauray = p_colmol(iplon,i_lay) * raylc(ig)
              p_taug(iplon,i_lay,ig) = p_colo3(iplon,i_lay) * abso3bc(ig)
              p_taur(iplon,i_lay,ig) = z_tauray
            ENDDO
          ENDIF
       ENDDO
    ENDDO

    !$ACC LOOP SEQ
    DO ig = 1 , ng25
      !$ACC LOOP SEQ
      DO i_lay = laytrop_max+1, i_nlayers
        !$ACC LOOP GANG(STATIC:1) VECTOR PRIVATE(z_tauray)
        DO iplon = KIDIA, KFDIA
          z_tauray = p_colmol(iplon,i_lay) * raylc(ig)
          p_taug(iplon,i_lay,ig) = p_colo3(iplon,i_lay) * abso3bc(ig)
          p_taur(iplon,i_lay,ig) = z_tauray
        ENDDO
      ENDDO
    ENDDO
    !$ACC END PARALLEL

    !$ACC WAIT
    !$ACC END DATA

END SUBROUTINE SRTM_TAUMOL25
