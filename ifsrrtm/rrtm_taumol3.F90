!----------------------------------------------------------------------------
SUBROUTINE RRTM_TAUMOL3 (KIDIA,KFDIA,KLEV,taug,&
 & P_TAUAERL,FAC00,FAC01,FAC10,FAC11,FORFAC,FORFRAC,INDFOR,JP,JT,jt1,ONEMINUS,&
 & COLH2O,COLCO2,COLN2O,COLDRY,LAYTROP,SELFFAC,SELFFRAC,INDSELF,FRACS, &
 & RAT_H2OCO2, RAT_H2OCO2_1,MINORFRAC,INDMINOR)

!     BAND 3:  500-630 cm-1 (low - H2O,CO2; high - H2O,CO2)

!     AUTHOR.
!     -------
!      JJMorcrette, ECMWF

!     MODIFICATIONS.
!     --------------
!      M.Hamrud      01-Oct-2003 CY28 Cleaning
!      NEC           25-Oct-2007 Optimisations
!      JJMorcrette 20110613 flexible number of g-points
!      ABozzo 20130517 updated to rrtmg_lw_v4.85:
!     band 3:  500-630 cm-1 (low key - h2o,co2; low minor - n2o)
!                           (high key - h2o,co2; high minor - n2o)
! ---------------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM     ,JPRB

USE PARRRTM  , ONLY : JPBAND
USE YOERRTM  , ONLY : JPGPT  ,NG3   ,NGS2
USE YOERRTWN , ONLY : NSPA   ,NSPB
USE YOERRTA3 , ONLY : ABSA   ,ABSB   ,FRACREFA, FRACREFB,&
 & FORREF   ,SELFREF , KA_MN2O ,  KB_MN2O
USE YOERRTRF, ONLY : CHI_MLS

IMPLICIT NONE

INTEGER(KIND=JPIM),INTENT(IN)    :: KIDIA
INTEGER(KIND=JPIM),INTENT(IN)    :: KFDIA
INTEGER(KIND=JPIM),INTENT(IN)    :: KLEV
REAL(KIND=JPRB)   ,INTENT(INOUT) :: taug(KIDIA:KFDIA,JPGPT,KLEV)
REAL(KIND=JPRB)   ,INTENT(IN)    :: P_TAUAERL(KIDIA:KFDIA,KLEV,JPBAND)
REAL(KIND=JPRB)   ,INTENT(IN)    :: FAC00(KIDIA:KFDIA,KLEV)
REAL(KIND=JPRB)   ,INTENT(IN)    :: FAC01(KIDIA:KFDIA,KLEV)
REAL(KIND=JPRB)   ,INTENT(IN)    :: FAC10(KIDIA:KFDIA,KLEV)
REAL(KIND=JPRB)   ,INTENT(IN)    :: FAC11(KIDIA:KFDIA,KLEV)
REAL(KIND=JPRB)   ,INTENT(IN)    :: FORFAC(KIDIA:KFDIA,KLEV)
INTEGER(KIND=JPIM),INTENT(IN)    :: JP(KIDIA:KFDIA,KLEV)
INTEGER(KIND=JPIM),INTENT(IN)    :: JT(KIDIA:KFDIA,KLEV)
INTEGER(KIND=JPIM),INTENT(IN)    :: jt1(KIDIA:KFDIA,KLEV)
REAL(KIND=JPRB)   ,INTENT(IN)    :: ONEMINUS
REAL(KIND=JPRB)   ,INTENT(IN)    :: COLH2O(KIDIA:KFDIA,KLEV)
REAL(KIND=JPRB)   ,INTENT(IN)    :: COLCO2(KIDIA:KFDIA,KLEV)
REAL(KIND=JPRB)   ,INTENT(IN)    :: COLN2O(KIDIA:KFDIA,KLEV)
REAL(KIND=JPRB)   ,INTENT(IN)    :: COLDRY(KIDIA:KFDIA,KLEV)
INTEGER(KIND=JPIM),INTENT(IN)    :: LAYTROP(KIDIA:KFDIA)
REAL(KIND=JPRB)   ,INTENT(IN)    :: SELFFAC(KIDIA:KFDIA,KLEV)
REAL(KIND=JPRB)   ,INTENT(IN)    :: SELFFRAC(KIDIA:KFDIA,KLEV)
INTEGER(KIND=JPIM),INTENT(IN)    :: INDSELF(KIDIA:KFDIA,KLEV)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: FRACS(KIDIA:KFDIA,JPGPT,KLEV)

REAL(KIND=JPRB)   ,INTENT(IN)   :: RAT_H2OCO2(KIDIA:KFDIA,KLEV)
REAL(KIND=JPRB)   ,INTENT(IN)   :: RAT_H2OCO2_1(KIDIA:KFDIA,KLEV)
INTEGER(KIND=JPIM),INTENT(IN)   :: INDFOR(KIDIA:KFDIA,KLEV)
REAL(KIND=JPRB)   ,INTENT(IN)   :: FORFRAC(KIDIA:KFDIA,KLEV)
REAL(KIND=JPRB)   ,INTENT(IN)   :: MINORFRAC(KIDIA:KFDIA,KLEV)
INTEGER(KIND=JPIM),INTENT(IN)   :: INDMINOR(KIDIA:KFDIA,KLEV)
! ---------------------------------------------------------------------------

REAL(KIND=JPRB) :: SPECCOMB,SPECCOMB1,SPECCOMB_MN2O, &
& SPECCOMB_PLANCK
REAL(KIND=JPRB) :: REFRAT_PLANCK_A, REFRAT_PLANCK_B, REFRAT_M_A, REFRAT_M_B

INTEGER(KIND=JPIM) :: IND0,IND1,INDS,INDF,INDM
INTEGER(KIND=JPIM) :: IG, JS, LAY, JS1,JMN2O,JPL

REAL(KIND=JPRB) :: FS, SPECMULT, SPECPARM,  &
 & FS1, SPECMULT1, SPECPARM1,   &
 & FMN2O, SPECMULT_MN2O, SPECPARM_MN2O,   &
 & fpl, SPECMULT_PLANCK, SPECPARM_PLANCK

REAL(KIND=JPRB) :: ADJFAC,ADJCOLN2O,RATN2O,CHI_N2O

 REAL(KIND=JPRB) ::  FAC000, FAC100, FAC200,&
 & FAC010, FAC110, FAC210, &
 & FAC001, FAC101, FAC201, &
 & FAC011, FAC111, FAC211
REAL(KIND=JPRB) :: P, P4, FK0, FK1, FK2
REAL(KIND=JPRB) :: TAUFOR,TAUSELF,N2OM1,N2OM2,ABSN2O,TAU_MAJOR(ng3),TAU_MAJOR1(ng3)

    !     local integer arrays
    INTEGER(KIND=JPIM) :: laytrop_min, laytrop_max
    integer(KIND=JPIM) :: ixc(KLEV), ixlow(KFDIA,KLEV), ixhigh(KFDIA,KLEV)
    INTEGER(KIND=JPIM) :: ich, icl, ixc0, ixp, jc, jl

#define MOD1(x) ((x) - AINT((x)))

    !$ACC DATA PRESENT(taug, P_TAUAERL, FAC00, FAC01, FAC10, FAC11, FORFAC, JP, &
    !$ACC             JT, jt1, COLH2O, COLCO2, COLN2O, COLDRY, LAYTROP, &
    !$ACC             SELFFAC, SELFFRAC, INDSELF, FRACS, RAT_H2OCO2, &
    !$ACC             RAT_H2OCO2_1, INDFOR, FORFRAC, MINORFRAC, INDMINOR)

#ifndef _OPENACC
    laytrop_min = MINVAL(laytrop)
    laytrop_max = MAXVAL(laytrop)

    ixlow  = 0
    ixhigh = 0
    ixc    = 0

    ! create index lists for mixed layers
    do lay = laytrop_min+1, laytrop_max
      icl = 0
      ich = 0
      do jc = KIDIA, KFDIA
        if ( lay <= laytrop(jc) ) then
          icl = icl + 1
          ixlow(icl,lay) = jc
        else
          ich = ich + 1
          ixhigh(ich,lay) = jc
        endif
      enddo
      ixc(lay) = icl
    enddo
#else
    laytrop_min = HUGE(laytrop_min)
    laytrop_max = -HUGE(laytrop_max)
    !$ACC PARALLEL DEFAULT(NONE) ASYNC(1)
    !$ACC LOOP GANG VECTOR REDUCTION(min:laytrop_min) REDUCTION(max:laytrop_max)
    do jc = KIDIA,KFDIA
      laytrop_min = MIN(laytrop_min, laytrop(jc))
      laytrop_max = MAX(laytrop_max, laytrop(jc))
    end do
    !$ACC END PARALLEL
#endif

!     Compute the optical depth by interpolating in ln(pressure),
!     temperature, and appropriate species.  Below LAYTROP, the water
!     vapor self-continuum is interpolated (in temperature) separately.


! Minor gas mapping levels:
!     lower - n2o, p = 706.272 mbar, t = 278.94 k
!     upper - n2o, p = 95.58 mbar, t = 215.7 k

!  P = 212.725 mb
      REFRAT_PLANCK_A = CHI_MLS(1,9)/CHI_MLS(2,9)

!  P = 95.58 mb
      REFRAT_PLANCK_B = CHI_MLS(1,13)/CHI_MLS(2,13)

!  P = 706.270mb
      REFRAT_M_A = CHI_MLS(1,3)/CHI_MLS(2,3)

!  P = 95.58 mb
      REFRAT_M_B = CHI_MLS(1,13)/CHI_MLS(2,13)

      ! Lower atmosphere loop
      !$ACC WAIT
      !$ACC PARALLEL DEFAULT(NONE) ASYNC(1)
      !$ACC LOOP GANG VECTOR COLLAPSE(2) PRIVATE(speccomb, specparm, specmult, js, fs, speccomb1, specparm1, specmult1, &
      !$ACC   js1, fs1, speccomb_mn2o, specparm_mn2o, specmult_mn2o, jmn2o, fmn2o, chi_n2o, ratn2o, adjfac, adjcoln2o, &
      !$ACC   speccomb_planck, specparm_planck, specmult_planck, jpl, fpl, ind0, ind1, inds, indf, indm, p, p4, fk0, &
      !$ACC   fk1, fk2, fac000, fac100, fac200, fac010, fac110, fac210, fac001, fac101, fac201, fac011, fac111, &
      !$ACC   fac211, tau_major, tau_major1)
      do lay = 1, laytrop_min
        do jl = KIDIA, KFDIA

          speccomb = colh2o(jl,lay) + rat_h2oco2(jl,lay)*colco2(jl,lay)
          specparm = MIN(colh2o(jl,lay)/speccomb,oneminus)
          specmult = 8._JPRB*(specparm)
          js = 1 + int(specmult)
          fs = MOD1(specmult)

          speccomb1 = colh2o(jl,lay) + rat_h2oco2_1(jl,lay)*colco2(jl,lay)
          specparm1 = MIN(colh2o(jl,lay)/speccomb1,oneminus)
          specmult1 = 8._JPRB*(specparm1)
          js1 = 1 + int(specmult1)
          fs1 = MOD1(specmult1)

          speccomb_mn2o = colh2o(jl,lay) + refrat_m_a*colco2(jl,lay)
          specparm_mn2o = MIN(colh2o(jl,lay)/speccomb_mn2o,oneminus)
          specmult_mn2o = 8._JPRB*specparm_mn2o
          jmn2o = 1 + int(specmult_mn2o)
          fmn2o = MOD1(specmult_mn2o)
          !  In atmospheres where the amount of N2O is too great to be considered
          !  a minor species, adjust the column amount of N2O by an empirical factor
          !  to obtain the proper contribution.
          chi_n2o = coln2o(jl,lay)/coldry(jl,lay)
          ratn2o = 1.e20_JPRB*chi_n2o/chi_mls(4,jp(jl,lay)+1)
          if (ratn2o .gt. 1.5_JPRB) then
            adjfac = 0.5_JPRB+(ratn2o-0.5_JPRB)**0.65_JPRB
            adjcoln2o = adjfac*chi_mls(4,jp(jl,lay)+1)*coldry(jl,lay)*1.e-20_JPRB
          else
            adjcoln2o = coln2o(jl,lay)
          endif

          speccomb_planck = colh2o(jl,lay)+refrat_planck_a*colco2(jl,lay)
          specparm_planck = MIN(colh2o(jl,lay)/speccomb_planck,oneminus)
          specmult_planck = 8._JPRB*specparm_planck
          jpl = 1 + int(specmult_planck)
          fpl = MOD1(specmult_planck)

          ind0 = ((jp(jl,lay)-1)*5+(jt(jl,lay)-1))*nspa(3) + js
          ind1 = (jp(jl,lay)*5+(jt1(jl,lay)-1))*nspa(3) + js1
          inds = indself(jl,lay)
          indf = indfor(jl,lay)
          indm = indminor(jl,lay)

          if (specparm .lt. 0.125_JPRB) then
            p = fs - 1._JPRB
            p4 = p**4
            fk0 = p4
            fk1 = 1._JPRB - p - 2.0_JPRB*p4
            fk2 = p + p4
            fac000 = fk0*fac00(jl,lay)
            fac100 = fk1*fac00(jl,lay)
            fac200 = fk2*fac00(jl,lay)
            fac010 = fk0*fac10(jl,lay)
            fac110 = fk1*fac10(jl,lay)
            fac210 = fk2*fac10(jl,lay)
          else if (specparm .gt. 0.875_JPRB) then
            p = -fs
            p4 = p**4
            fk0 = p4
            fk1 = 1._JPRB - p - 2.0_JPRB*p4
            fk2 = p + p4
            fac000 = fk0*fac00(jl,lay)
            fac100 = fk1*fac00(jl,lay)
            fac200 = fk2*fac00(jl,lay)
            fac010 = fk0*fac10(jl,lay)
            fac110 = fk1*fac10(jl,lay)
            fac210 = fk2*fac10(jl,lay)
          else
            fac000 = (1._JPRB - fs) * fac00(jl,lay)
            fac010 = (1._JPRB - fs) * fac10(jl,lay)
            fac100 = fs * fac00(jl,lay)
            fac110 = fs * fac10(jl,lay)
            fac200 = 0._JPRB
            fac210 = 0._JPRB
          endif

          if (specparm1 .lt. 0.125_JPRB) then
            p = fs1 - 1._JPRB
            p4 = p**4
            fk0 = p4
            fk1 = 1._JPRB - p - 2.0_JPRB*p4
            fk2 = p + p4
            fac001 = fk0*fac01(jl,lay)
            fac101 = fk1*fac01(jl,lay)
            fac201 = fk2*fac01(jl,lay)
            fac011 = fk0*fac11(jl,lay)
            fac111 = fk1*fac11(jl,lay)
            fac211 = fk2*fac11(jl,lay)
          else if (specparm1 .gt. 0.875_JPRB) then
            p = -fs1
            p4 = p**4
            fk0 = p4
            fk1 = 1._JPRB - p - 2.0_JPRB*p4
            fk2 = p + p4
            fac001 = fk0*fac01(jl,lay)
            fac101 = fk1*fac01(jl,lay)
            fac201 = fk2*fac01(jl,lay)
            fac011 = fk0*fac11(jl,lay)
            fac111 = fk1*fac11(jl,lay)
            fac211 = fk2*fac11(jl,lay)
          else
            fac001 = (1._JPRB - fs1) * fac01(jl,lay)
            fac011 = (1._JPRB - fs1) * fac11(jl,lay)
            fac101 = fs1 * fac01(jl,lay)
            fac111 = fs1 * fac11(jl,lay)
            fac201 = 0._JPRB
            fac211 = 0._JPRB
          endif

          if (specparm .lt. 0.125_JPRB) then
!$NEC unroll(NG3)
            tau_major(1:ng3) = speccomb *    &
             (fac000 * absa(ind0,1:ng3)    + &
              fac100 * absa(ind0+1,1:ng3)  + &
              fac200 * absa(ind0+2,1:ng3)  + &
              fac010 * absa(ind0+9,1:ng3)  + &
              fac110 * absa(ind0+10,1:ng3) + &
              fac210 * absa(ind0+11,1:ng3))
          else if (specparm .gt. 0.875_JPRB) then
!$NEC unroll(NG3)
            tau_major(1:ng3) = speccomb *   &
             (fac200 * absa(ind0-1,1:ng3) + &
              fac100 * absa(ind0,1:ng3)   + &
              fac000 * absa(ind0+1,1:ng3) + &
              fac210 * absa(ind0+8,1:ng3) + &
              fac110 * absa(ind0+9,1:ng3) + &
              fac010 * absa(ind0+10,1:ng3))
          else
!$NEC unroll(NG3)
            tau_major(1:ng3) = speccomb *   &
             (fac000 * absa(ind0,1:ng3)   + &
              fac100 * absa(ind0+1,1:ng3) + &
              fac010 * absa(ind0+9,1:ng3) + &
              fac110 * absa(ind0+10,1:ng3))
          endif

          if (specparm1 .lt. 0.125_JPRB) then
!$NEC unroll(NG3)
            tau_major1(1:ng3) = speccomb1 *  &
             (fac001 * absa(ind1,1:ng3)    + &
              fac101 * absa(ind1+1,1:ng3)  + &
              fac201 * absa(ind1+2,1:ng3)  + &
              fac011 * absa(ind1+9,1:ng3)  + &
              fac111 * absa(ind1+10,1:ng3) + &
              fac211 * absa(ind1+11,1:ng3))
          else if (specparm1 .gt. 0.875_JPRB) then
!$NEC unroll(NG3)
            tau_major1(1:ng3) = speccomb1 * &
             (fac201 * absa(ind1-1,1:ng3) + &
              fac101 * absa(ind1,1:ng3)   + &
              fac001 * absa(ind1+1,1:ng3) + &
              fac211 * absa(ind1+8,1:ng3) + &
              fac111 * absa(ind1+9,1:ng3) + &
              fac011 * absa(ind1+10,1:ng3))
          else
!$NEC unroll(NG3)
            tau_major1(1:ng3) = speccomb1 * &
             (fac001 * absa(ind1,1:ng3)   + &
              fac101 * absa(ind1+1,1:ng3) + &
              fac011 * absa(ind1+9,1:ng3) + &
              fac111 * absa(ind1+10,1:ng3))
          endif

          !$ACC LOOP SEQ PRIVATE(tauself, taufor, n2om1, n2om2, absn2o)
!$NEC unroll(NG3)
          do ig = 1, ng3
            tauself = selffac(jl,lay)* (selfref(inds,ig) + selffrac(jl,lay) * &
                 (selfref(inds+1,ig) - selfref(inds,ig)))
            taufor = forfac(jl,lay) * (forref(indf,ig) + forfrac(jl,lay) * &
                 (forref(indf+1,ig) - forref(indf,ig)))
            n2om1 = ka_mn2o(jmn2o,indm,ig) + fmn2o * &
                 (ka_mn2o(jmn2o+1,indm,ig) - ka_mn2o(jmn2o,indm,ig))
            n2om2 = ka_mn2o(jmn2o,indm+1,ig) + fmn2o * &
                 (ka_mn2o(jmn2o+1,indm+1,ig) - ka_mn2o(jmn2o,indm+1,ig))
            absn2o = n2om1 + minorfrac(jl,lay) * (n2om2 - n2om1)

            taug(jl,ngs2+ig,lay) = tau_major(ig) + tau_major1(ig) &
                 + tauself + taufor &
                 + adjcoln2o*absn2o
            fracs(jl,ngs2+ig,lay) = fracrefa(ig,jpl) + fpl * &
                 (fracrefa(ig,jpl+1)-fracrefa(ig,jpl))
          enddo
        enddo

      enddo
      !$ACC END PARALLEL

      ! Upper atmosphere loop
      !$ACC PARALLEL DEFAULT(NONE) ASYNC(1)
      !$ACC LOOP GANG VECTOR COLLAPSE(2) PRIVATE(speccomb, specparm, specmult, js, fs, speccomb1, specparm1, specmult1, &
      !$ACC   js1, fs1, speccomb_mn2o, specparm_mn2o, specmult_mn2o, jmn2o, fmn2o, chi_n2o, ratn2o, adjfac, adjcoln2o, &
      !$ACC   speccomb_planck, specparm_planck, specmult_planck, jpl, fpl, ind0, ind1, indf, indm, fac000, fac100, &
      !$ACC   fac010, fac110, fac001, fac101, fac011, fac111)
      do lay = laytrop_max+1, KLEV
        do jl = KIDIA, KFDIA

          speccomb = colh2o(jl,lay) + rat_h2oco2(jl,lay)*colco2(jl,lay)
          specparm = MIN(colh2o(jl,lay)/speccomb,oneminus)
          specmult = 4._JPRB*(specparm)
          js = 1 + int(specmult)
          fs = MOD1(specmult)

          speccomb1 = colh2o(jl,lay) + rat_h2oco2_1(jl,lay)*colco2(jl,lay)
          specparm1 = MIN(colh2o(jl,lay)/speccomb1,oneminus)
          specmult1 = 4._JPRB*(specparm1)
          js1 = 1 + int(specmult1)
          fs1 = MOD1(specmult1)

          fac000 = (1._JPRB - fs) * fac00(jl,lay)
          fac010 = (1._JPRB - fs) * fac10(jl,lay)
          fac100 = fs * fac00(jl,lay)
          fac110 = fs * fac10(jl,lay)
          fac001 = (1._JPRB - fs1) * fac01(jl,lay)
          fac011 = (1._JPRB - fs1) * fac11(jl,lay)
          fac101 = fs1 * fac01(jl,lay)
          fac111 = fs1 * fac11(jl,lay)

          speccomb_mn2o = colh2o(jl,lay) + refrat_m_b*colco2(jl,lay)
          specparm_mn2o = MIN(colh2o(jl,lay)/speccomb_mn2o,oneminus)
          specmult_mn2o = 4._JPRB*specparm_mn2o
          jmn2o = 1 + int(specmult_mn2o)
          fmn2o = MOD1(specmult_mn2o)
          !  In atmospheres where the amount of N2O is too great to be considered
          !  a minor species, adjust the column amount of N2O by an empirical factor
          !  to obtain the proper contribution.
          chi_n2o = coln2o(jl,lay)/coldry(jl,lay)
          ratn2o = 1.e20_JPRB*chi_n2o/chi_mls(4,jp(jl,lay)+1)
          if (ratn2o .gt. 1.5_JPRB) then
            adjfac = 0.5_JPRB+(ratn2o-0.5_JPRB)**0.65_JPRB
            adjcoln2o = adjfac*chi_mls(4,jp(jl,lay)+1)*coldry(jl,lay)*1.e-20_JPRB
          else
            adjcoln2o = coln2o(jl,lay)
          endif

          speccomb_planck = colh2o(jl,lay)+refrat_planck_b*colco2(jl,lay)
          specparm_planck = MIN(colh2o(jl,lay)/speccomb_planck,oneminus)
          specmult_planck = 4._JPRB*specparm_planck
          jpl= 1 + int(specmult_planck)
          fpl = MOD1(specmult_planck)

          ind0 = ((jp(jl,lay)-13)*5+(jt(jl,lay)-1))*nspb(3) + js
          ind1 = ((jp(jl,lay)-12)*5+(jt1(jl,lay)-1))*nspb(3) + js1
          indf = indfor(jl,lay)
          indm = indminor(jl,lay)
          !$ACC LOOP SEQ PRIVATE(taufor, n2om1, n2om2, absn2o)
!$NEC unroll(NG3)
          do ig = 1, ng3
            taufor = forfac(jl,lay) * (forref(indf,ig) + &
                 forfrac(jl,lay) * (forref(indf+1,ig) - forref(indf,ig)))
            n2om1 = kb_mn2o(jmn2o,indm,ig) + fmn2o * &
                 (kb_mn2o(jmn2o+1,indm,ig)-kb_mn2o(jmn2o,indm,ig))
            n2om2 = kb_mn2o(jmn2o,indm+1,ig) + fmn2o * &
                 (kb_mn2o(jmn2o+1,indm+1,ig)-kb_mn2o(jmn2o,indm+1,ig))
            absn2o = n2om1 + minorfrac(jl,lay) * (n2om2 - n2om1)
            taug(jl,ngs2+ig,lay) = speccomb * &
                 (fac000 * absb(ind0,ig) + &
                 fac100 * absb(ind0+1,ig) + &
                 fac010 * absb(ind0+5,ig) + &
                 fac110 * absb(ind0+6,ig)) &
                 + speccomb1 * &
                 (fac001 * absb(ind1,ig) +  &
                 fac101 * absb(ind1+1,ig) + &
                 fac011 * absb(ind1+5,ig) + &
                 fac111 * absb(ind1+6,ig))  &
                 + taufor &
                 + adjcoln2o*absn2o
            fracs(jl,ngs2+ig,lay) = fracrefb(ig,jpl) + fpl * &
                 (fracrefb(ig,jpl+1)-fracrefb(ig,jpl))
          enddo
        enddo

      enddo
      !$ACC END PARALLEL

      IF (laytrop_max /= laytrop_min) THEN
        ! Mixed loop
        ! Lower atmosphere part
        !$ACC PARALLEL DEFAULT(NONE) ASYNC(1)
        !$ACC LOOP GANG VECTOR COLLAPSE(2) PRIVATE(speccomb, specparm, specmult, js, fs, speccomb1, specparm1, &
        !$ACC   specmult1, js1, fs1, speccomb_mn2o, specparm_mn2o, specmult_mn2o, jmn2o, fmn2o, chi_n2o, ratn2o, adjfac,&
        !$ACC   adjcoln2o, speccomb_planck, specparm_planck, specmult_planck, jpl, fpl, ind0, ind1, inds, indf, indm, p,&
        !$ACC   p4, fk0, &
        !$ACC   fk1, fk2, fac000, fac100, fac200, fac010, fac110, fac210, fac001, fac101, fac201, &
        !$ACC   fac011, fac111, fac211, tau_major, tau_major1)
        do lay = laytrop_min+1, laytrop_max
#ifdef _OPENACC
          do jl = KIDIA, KFDIA
            if ( lay <= laytrop(jl) ) then
#else

          ixc0 = ixc(lay)

!$NEC ivdep
          do ixp = 1, ixc0
            jl = ixlow(ixp,lay)
#endif

            speccomb = colh2o(jl,lay) + rat_h2oco2(jl,lay)*colco2(jl,lay)
            specparm = MIN(colh2o(jl,lay)/speccomb,oneminus)
            specmult = 8._JPRB*(specparm)
            js = 1 + int(specmult)
            fs = MOD1(specmult)

            speccomb1 = colh2o(jl,lay) + rat_h2oco2_1(jl,lay)*colco2(jl,lay)
            specparm1 = MIN(colh2o(jl,lay)/speccomb1,oneminus)
            specmult1 = 8._JPRB*(specparm1)
            js1 = 1 + int(specmult1)
            fs1 = MOD1(specmult1)

            speccomb_mn2o = colh2o(jl,lay) + refrat_m_a*colco2(jl,lay)
            specparm_mn2o = MIN(colh2o(jl,lay)/speccomb_mn2o,oneminus)
            specmult_mn2o = 8._JPRB*specparm_mn2o
            jmn2o = 1 + int(specmult_mn2o)
            fmn2o = MOD1(specmult_mn2o)
            !  In atmospheres where the amount of N2O is too great to be considered
            !  a minor species, adjust the column amount of N2O by an empirical factor
            !  to obtain the proper contribution.
            chi_n2o = coln2o(jl,lay)/coldry(jl,lay)
            ratn2o = 1.e20_JPRB*chi_n2o/chi_mls(4,jp(jl,lay)+1)
            if (ratn2o .gt. 1.5_JPRB) then
              adjfac = 0.5_JPRB+(ratn2o-0.5_JPRB)**0.65_JPRB
              adjcoln2o = adjfac*chi_mls(4,jp(jl,lay)+1)*coldry(jl,lay)*1.e-20_JPRB
            else
              adjcoln2o = coln2o(jl,lay)
            endif

            speccomb_planck = colh2o(jl,lay)+refrat_planck_a*colco2(jl,lay)
            specparm_planck = MIN(colh2o(jl,lay)/speccomb_planck,oneminus)
            specmult_planck = 8._JPRB*specparm_planck
            jpl = 1 + int(specmult_planck)
            fpl = MOD1(specmult_planck)

            ind0 = ((jp(jl,lay)-1)*5+(jt(jl,lay)-1))*nspa(3) + js
            ind1 = (jp(jl,lay)*5+(jt1(jl,lay)-1))*nspa(3) + js1
            inds = indself(jl,lay)
            indf = indfor(jl,lay)
            indm = indminor(jl,lay)

            if (specparm .lt. 0.125_JPRB) then
              p = fs - 1._JPRB
              p4 = p**4
              fk0 = p4
              fk1 = 1._JPRB - p - 2.0_JPRB*p4
              fk2 = p + p4
              fac000 = fk0*fac00(jl,lay)
              fac100 = fk1*fac00(jl,lay)
              fac200 = fk2*fac00(jl,lay)
              fac010 = fk0*fac10(jl,lay)
              fac110 = fk1*fac10(jl,lay)
              fac210 = fk2*fac10(jl,lay)
            else if (specparm .gt. 0.875_JPRB) then
              p = -fs
              p4 = p**4
              fk0 = p4
              fk1 = 1._JPRB - p - 2.0_JPRB*p4
              fk2 = p + p4
              fac000 = fk0*fac00(jl,lay)
              fac100 = fk1*fac00(jl,lay)
              fac200 = fk2*fac00(jl,lay)
              fac010 = fk0*fac10(jl,lay)
              fac110 = fk1*fac10(jl,lay)
              fac210 = fk2*fac10(jl,lay)
            else
              fac000 = (1._JPRB - fs) * fac00(jl,lay)
              fac010 = (1._JPRB - fs) * fac10(jl,lay)
              fac100 = fs * fac00(jl,lay)
              fac110 = fs * fac10(jl,lay)
              fac200 = 0._JPRB
              fac210 = 0._JPRB
            endif

            if (specparm1 .lt. 0.125_JPRB) then
              p = fs1 - 1._JPRB
              p4 = p**4
              fk0 = p4
              fk1 = 1._JPRB - p - 2.0_JPRB*p4
              fk2 = p + p4
              fac001 = fk0*fac01(jl,lay)
              fac101 = fk1*fac01(jl,lay)
              fac201 = fk2*fac01(jl,lay)
              fac011 = fk0*fac11(jl,lay)
              fac111 = fk1*fac11(jl,lay)
              fac211 = fk2*fac11(jl,lay)
            else if (specparm1 .gt. 0.875_JPRB) then
              p = -fs1
              p4 = p**4
              fk0 = p4
              fk1 = 1._JPRB - p - 2.0_JPRB*p4
              fk2 = p + p4
              fac001 = fk0*fac01(jl,lay)
              fac101 = fk1*fac01(jl,lay)
              fac201 = fk2*fac01(jl,lay)
              fac011 = fk0*fac11(jl,lay)
              fac111 = fk1*fac11(jl,lay)
              fac211 = fk2*fac11(jl,lay)
            else
              fac001 = (1._JPRB - fs1) * fac01(jl,lay)
              fac011 = (1._JPRB - fs1) * fac11(jl,lay)
              fac101 = fs1 * fac01(jl,lay)
              fac111 = fs1 * fac11(jl,lay)
              fac201 = 0._JPRB
              fac211 = 0._JPRB
            endif

            if (specparm .lt. 0.125_JPRB) then
!$NEC unroll(NG3)
              tau_major(1:ng3) = speccomb *    &
              (fac000 * absa(ind0,1:ng3)    + &
                fac100 * absa(ind0+1,1:ng3)  + &
                fac200 * absa(ind0+2,1:ng3)  + &
                fac010 * absa(ind0+9,1:ng3)  + &
                fac110 * absa(ind0+10,1:ng3) + &
                fac210 * absa(ind0+11,1:ng3))
            else if (specparm .gt. 0.875_JPRB) then
!$NEC unroll(NG3)
              tau_major(1:ng3) = speccomb *   &
              (fac200 * absa(ind0-1,1:ng3) + &
                fac100 * absa(ind0,1:ng3)   + &
                fac000 * absa(ind0+1,1:ng3) + &
                fac210 * absa(ind0+8,1:ng3) + &
                fac110 * absa(ind0+9,1:ng3) + &
                fac010 * absa(ind0+10,1:ng3))
            else
!$NEC unroll(NG3)
              tau_major(1:ng3) = speccomb *   &
              (fac000 * absa(ind0,1:ng3)   + &
                fac100 * absa(ind0+1,1:ng3) + &
                fac010 * absa(ind0+9,1:ng3) + &
                fac110 * absa(ind0+10,1:ng3))
            endif

            if (specparm1 .lt. 0.125_JPRB) then
!$NEC unroll(NG3)
              tau_major1(1:ng3) = speccomb1 *  &
              (fac001 * absa(ind1,1:ng3)    + &
                fac101 * absa(ind1+1,1:ng3)  + &
                fac201 * absa(ind1+2,1:ng3)  + &
                fac011 * absa(ind1+9,1:ng3)  + &
                fac111 * absa(ind1+10,1:ng3) + &
                fac211 * absa(ind1+11,1:ng3))
            else if (specparm1 .gt. 0.875_JPRB) then
!$NEC unroll(NG3)
              tau_major1(1:ng3) = speccomb1 * &
              (fac201 * absa(ind1-1,1:ng3) + &
                fac101 * absa(ind1,1:ng3)   + &
                fac001 * absa(ind1+1,1:ng3) + &
                fac211 * absa(ind1+8,1:ng3) + &
                fac111 * absa(ind1+9,1:ng3) + &
                fac011 * absa(ind1+10,1:ng3))
            else
!$NEC unroll(NG3)
              tau_major1(1:ng3) = speccomb1 * &
              (fac001 * absa(ind1,1:ng3)   + &
                fac101 * absa(ind1+1,1:ng3) + &
                fac011 * absa(ind1+9,1:ng3) + &
                fac111 * absa(ind1+10,1:ng3))
            endif

!$NEC unroll(NG3)
            !$ACC LOOP SEQ PRIVATE(tauself, taufor, n2om1, n2om2, absn2o)
            do ig = 1, ng3
              tauself = selffac(jl,lay)* (selfref(inds,ig) + selffrac(jl,lay) * &
                  (selfref(inds+1,ig) - selfref(inds,ig)))
              taufor = forfac(jl,lay) * (forref(indf,ig) + forfrac(jl,lay) * &
                  (forref(indf+1,ig) - forref(indf,ig)))
              n2om1 = ka_mn2o(jmn2o,indm,ig) + fmn2o * &
                  (ka_mn2o(jmn2o+1,indm,ig) - ka_mn2o(jmn2o,indm,ig))
              n2om2 = ka_mn2o(jmn2o,indm+1,ig) + fmn2o * &
                  (ka_mn2o(jmn2o+1,indm+1,ig) - ka_mn2o(jmn2o,indm+1,ig))
              absn2o = n2om1 + minorfrac(jl,lay) * (n2om2 - n2om1)

              taug(jl,ngs2+ig,lay) = tau_major(ig) + tau_major1(ig) &
                  + tauself + taufor &
                  + adjcoln2o*absn2o
              fracs(jl,ngs2+ig,lay) = fracrefa(ig,jpl) + fpl * &
                  (fracrefa(ig,jpl+1)-fracrefa(ig,jpl))
            enddo
#ifdef _OPENACC
         else
#else
          enddo

          ! Upper atmosphere part
          ixc0 = KFDIA - KIDIA + 1 - ixc0
!$NEC ivdep
          do ixp = 1, ixc0
            jl = ixhigh(ixp,lay)
#endif

            speccomb = colh2o(jl,lay) + rat_h2oco2(jl,lay)*colco2(jl,lay)
            specparm = MIN(colh2o(jl,lay)/speccomb,oneminus)
            specmult = 4._JPRB*(specparm)
            js = 1 + int(specmult)
            fs = MOD1(specmult)

            speccomb1 = colh2o(jl,lay) + rat_h2oco2_1(jl,lay)*colco2(jl,lay)
            specparm1 = MIN(colh2o(jl,lay)/speccomb1,oneminus)
            specmult1 = 4._JPRB*(specparm1)
            js1 = 1 + int(specmult1)
            fs1 = MOD1(specmult1)

            fac000 = (1._JPRB - fs) * fac00(jl,lay)
            fac010 = (1._JPRB - fs) * fac10(jl,lay)
            fac100 = fs * fac00(jl,lay)
            fac110 = fs * fac10(jl,lay)
            fac001 = (1._JPRB - fs1) * fac01(jl,lay)
            fac011 = (1._JPRB - fs1) * fac11(jl,lay)
            fac101 = fs1 * fac01(jl,lay)
            fac111 = fs1 * fac11(jl,lay)

            speccomb_mn2o = colh2o(jl,lay) + refrat_m_b*colco2(jl,lay)
            specparm_mn2o = MIN(colh2o(jl,lay)/speccomb_mn2o,oneminus)
            specmult_mn2o = 4._JPRB*specparm_mn2o
            jmn2o = 1 + int(specmult_mn2o)
            fmn2o = MOD1(specmult_mn2o)
            !  In atmospheres where the amount of N2O is too great to be considered
            !  a minor species, adjust the column amount of N2O by an empirical factor
            !  to obtain the proper contribution.
            chi_n2o = coln2o(jl,lay)/coldry(jl,lay)
            ratn2o = 1.e20_JPRB*chi_n2o/chi_mls(4,jp(jl,lay)+1)
            if (ratn2o .gt. 1.5_JPRB) then
              adjfac = 0.5_JPRB+(ratn2o-0.5_JPRB)**0.65_JPRB
              adjcoln2o = adjfac*chi_mls(4,jp(jl,lay)+1)*coldry(jl,lay)*1.e-20_JPRB
            else
              adjcoln2o = coln2o(jl,lay)
            endif

            speccomb_planck = colh2o(jl,lay)+refrat_planck_b*colco2(jl,lay)
            specparm_planck = MIN(colh2o(jl,lay)/speccomb_planck,oneminus)
            specmult_planck = 4._JPRB*specparm_planck
            jpl= 1 + int(specmult_planck)
            fpl = MOD1(specmult_planck)

            ind0 = ((jp(jl,lay)-13)*5+(jt(jl,lay)-1))*nspb(3) + js
            ind1 = ((jp(jl,lay)-12)*5+(jt1(jl,lay)-1))*nspb(3) + js1
            indf = indfor(jl,lay)
            indm = indminor(jl,lay)
!$NEC unroll(NG3)
            !$ACC LOOP SEQ PRIVATE(taufor, n2om1, n2om2, absn2o)
            do ig = 1, ng3
              taufor = forfac(jl,lay) * (forref(indf,ig) + &
                  forfrac(jl,lay) * (forref(indf+1,ig) - forref(indf,ig)))
              n2om1 = kb_mn2o(jmn2o,indm,ig) + fmn2o * &
                  (kb_mn2o(jmn2o+1,indm,ig)-kb_mn2o(jmn2o,indm,ig))
              n2om2 = kb_mn2o(jmn2o,indm+1,ig) + fmn2o * &
                  (kb_mn2o(jmn2o+1,indm+1,ig)-kb_mn2o(jmn2o,indm+1,ig))
              absn2o = n2om1 + minorfrac(jl,lay) * (n2om2 - n2om1)
              taug(jl,ngs2+ig,lay) = speccomb * &
                  (fac000 * absb(ind0,ig) + &
                  fac100 * absb(ind0+1,ig) + &
                  fac010 * absb(ind0+5,ig) + &
                  fac110 * absb(ind0+6,ig)) &
                  + speccomb1 * &
                  (fac001 * absb(ind1,ig) +  &
                  fac101 * absb(ind1+1,ig) + &
                  fac011 * absb(ind1+5,ig) + &
                  fac111 * absb(ind1+6,ig))  &
                  + taufor &
                  + adjcoln2o*absn2o
              fracs(jl,ngs2+ig,lay) = fracrefb(ig,jpl) + fpl * &
                  (fracrefb(ig,jpl+1)-fracrefb(ig,jpl))
            enddo
#ifdef _OPENACC
           endif
#endif
          enddo

        enddo
        !$ACC END PARALLEL

      ENDIF

      !$ACC END DATA

END SUBROUTINE RRTM_TAUMOL3
