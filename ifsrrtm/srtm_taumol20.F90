SUBROUTINE SRTM_TAUMOL20 &
 & ( KIDIA   , KFDIA    , KLEV,&
 & P_FAC00   , P_FAC01  , P_FAC10   , P_FAC11,&
 & K_JP      , K_JT     , K_JT1,&
 & P_COLH2O  , P_COLCH4 , P_COLMOL,&
 & K_LAYTROP , P_SELFFAC, P_SELFFRAC, K_INDSELF  , P_FORFAC, P_FORFRAC, K_INDFOR,&
 & P_SFLUXZEN, P_TAUG   , P_TAUR    , PRMU0   &
 & )  

!     Written by Eli J. Mlawer, Atmospheric & Environmental Research.

!     BAND 20:  5150-6150 cm-1 (low - H2O; high - H2O)

! Modifications
!        M.Hamrud      01-Oct-2003 CY28 Cleaning

!     JJMorcrette 2003-02-24 adapted to ECMWF environment
!        D.Salmond  31-Oct-2007 Vector version in the style of RRTM from Meteo France & NEC
!     JJMorcrette 20110610 Flexible configuration for number of g-points

USE PARKIND1 , ONLY : JPIM, JPRB
USE YOMHOOK  , ONLY : LHOOK, DR_HOOK
USE PARSRTM  , ONLY : JPG
USE YOESRTM  , ONLY : NG20
USE YOESRTA20, ONLY : ABSA, ABSB, FORREFC, SELFREFC, SFLUXREFC, ABSCH4C, RAYL, LAYREFFR  
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
!- from INTFAC      
!- from INTIND
!- from PRECISE             
!- from PROFDATA             
!- from SELF             
INTEGER(KIND=JPIM) :: IG, IND0, IND1, INDS, INDF, I_LAY, I_LAYSOLFR(KIDIA:KFDIA), I_NLAYERS, IPLON
INTEGER(KIND=JPIM) :: laytrop_min, laytrop_max
INTEGER(KIND=JPIM) :: I_LAY_NEXT

REAL(KIND=JPRB) ::  &
 & Z_TAURAY  
REAL(KIND=JPRB) :: ZHOOK_HANDLE

    laytrop_min = MINVAL(k_laytrop(KIDIA:KFDIA))
    laytrop_max = MAXVAL(k_laytrop(KIDIA:KFDIA))

    i_nlayers = klev
    i_laysolfr(KIDIA:KFDIA) = k_laytrop(KIDIA:KFDIA)

    DO i_lay = 1, laytrop_min
       DO iplon = KIDIA, KFDIA
         IF (k_jp(iplon,i_lay) < layreffr                           &
              &    .AND. k_jp(iplon,i_lay+1) >= layreffr)           &
              &    i_laysolfr(iplon) = MIN(i_lay+1,k_laytrop(iplon))
         ind0 = ((k_jp(iplon,i_lay)-1)*5+(k_jt(iplon,i_lay)-1))*nspa(20) + 1
         ind1 = (k_jp(iplon,i_lay)*5+(k_jt1(iplon,i_lay)-1))*nspa(20) + 1
         inds = k_indself(iplon,i_lay)
         indf = k_indfor(iplon,i_lay)
         z_tauray = p_colmol(iplon,i_lay) * rayl
!$NEC unroll(NG20)
         DO ig = 1 , ng20
           p_taug(iplon,i_lay,ig) = p_colh2o(iplon,i_lay) *      &
                & ((p_fac00(iplon,i_lay) * absa(ind0,ig) +       &
                & p_fac10(iplon,i_lay) * absa(ind0+1,ig) +       &
                & p_fac01(iplon,i_lay) * absa(ind1,ig) +         &
                & p_fac11(iplon,i_lay) * absa(ind1+1,ig)) +      &
                & p_selffac(iplon,i_lay) * (selfrefc(inds,ig) +  &
                & p_selffrac(iplon,i_lay) *                      &
                & (selfrefc(inds+1,ig) - selfrefc(inds,ig))) +   &
                & p_forfac(iplon,i_lay) * (forrefc(indf,ig) +    &
                & p_forfrac(iplon,i_lay) *                       &
                & (forrefc(indf+1,ig) - forrefc(indf,ig))))      &
                & + p_colch4(iplon,i_lay) * absch4c(ig)
           p_taur(iplon,i_lay,ig) = z_tauray
           IF(i_lay == i_laysolfr(iplon)) p_sfluxzen(iplon,ig)=sfluxrefc(ig)
         ENDDO
       ENDDO
    ENDDO

    DO i_lay = laytrop_min+1, laytrop_max
       DO iplon = KIDIA, KFDIA
          IF (i_lay <= k_laytrop(iplon)) THEN
            IF (k_jp(iplon,i_lay) < layreffr                           &
                 &    .AND. k_jp(iplon,i_lay+1) >= layreffr)           &
                 &    i_laysolfr(iplon) = MIN(i_lay+1,k_laytrop(iplon))
            ind0 = ((k_jp(iplon,i_lay)-1)*5+(k_jt(iplon,i_lay)-1))*nspa(20) + 1
            ind1 = (k_jp(iplon,i_lay)*5+(k_jt1(iplon,i_lay)-1))*nspa(20) + 1
            inds = k_indself(iplon,i_lay)
            indf = k_indfor(iplon,i_lay)
            z_tauray = p_colmol(iplon,i_lay) * rayl
!$NEC unroll(NG20)
            DO ig = 1 , ng20
              p_taug(iplon,i_lay,ig) = p_colh2o(iplon,i_lay) *      &
                   & ((p_fac00(iplon,i_lay) * absa(ind0,ig) +       &
                   & p_fac10(iplon,i_lay) * absa(ind0+1,ig) +       &
                   & p_fac01(iplon,i_lay) * absa(ind1,ig) +         &
                   & p_fac11(iplon,i_lay) * absa(ind1+1,ig)) +      &
                   & p_selffac(iplon,i_lay) * (selfrefc(inds,ig) +  &
                   & p_selffrac(iplon,i_lay) *                      &
                   & (selfrefc(inds+1,ig) - selfrefc(inds,ig))) +   &
                   & p_forfac(iplon,i_lay) * (forrefc(indf,ig) +    &
                   & p_forfrac(iplon,i_lay) *                       &
                   & (forrefc(indf+1,ig) - forrefc(indf,ig))))      &
                   & + p_colch4(iplon,i_lay) * absch4c(ig)
              p_taur(iplon,i_lay,ig) = z_tauray
              IF(i_lay == i_laysolfr(iplon)) p_sfluxzen(iplon,ig)=sfluxrefc(ig)
            ENDDO
          ELSE
            ind0 = ((k_jp(iplon,i_lay)-13)*5+(k_jt(iplon,i_lay)-1))*nspb(20) + 1
            ind1 = ((k_jp(iplon,i_lay)-12)*5+(k_jt1(iplon,i_lay)-1))*nspb(20) +1
            indf = k_indfor(iplon,i_lay)
            z_tauray = p_colmol(iplon,i_lay) * rayl
!$NEC unroll(NG20)
            DO ig = 1 , ng20
              p_taug(iplon,i_lay,ig) = p_colh2o(iplon,i_lay) *   &
                   & (p_fac00(iplon,i_lay) * absb(ind0,ig) +     &
                   & p_fac10(iplon,i_lay) * absb(ind0+1,ig) +    &
                   & p_fac01(iplon,i_lay) * absb(ind1,ig) +      &
                   & p_fac11(iplon,i_lay) * absb(ind1+1,ig) +    &
                   & p_forfac(iplon,i_lay) * (forrefc(indf,ig) + &
                   & p_forfrac(iplon,i_lay) *                    &
                   & (forrefc(indf+1,ig) - forrefc(indf,ig)))) + &
                   & p_colch4(iplon,i_lay) * absch4c(ig)
              p_taur(iplon,i_lay,ig) = z_tauray
            ENDDO
          ENDIF
       ENDDO
    ENDDO

    DO i_lay = laytrop_max+1, i_nlayers
       DO iplon = KIDIA, KFDIA
         ind0 = ((k_jp(iplon,i_lay)-13)*5+(k_jt(iplon,i_lay)-1))*nspb(20) + 1
         ind1 = ((k_jp(iplon,i_lay)-12)*5+(k_jt1(iplon,i_lay)-1))*nspb(20) +1
         indf = k_indfor(iplon,i_lay)
         z_tauray = p_colmol(iplon,i_lay) * rayl
!$NEC unroll(NG20)
         DO ig = 1 , ng20
           p_taug(iplon,i_lay,ig) = p_colh2o(iplon,i_lay) *   &
                & (p_fac00(iplon,i_lay) * absb(ind0,ig) +     &
                & p_fac10(iplon,i_lay) * absb(ind0+1,ig) +    &
                & p_fac01(iplon,i_lay) * absb(ind1,ig) +      &
                & p_fac11(iplon,i_lay) * absb(ind1+1,ig) +    &
                & p_forfac(iplon,i_lay) * (forrefc(indf,ig) + &
                & p_forfrac(iplon,i_lay) *                    &
                & (forrefc(indf+1,ig) - forrefc(indf,ig)))) + &
                & p_colch4(iplon,i_lay) * absch4c(ig)
           p_taur(iplon,i_lay,ig) = z_tauray
         ENDDO
       ENDDO
    ENDDO


END SUBROUTINE SRTM_TAUMOL20
