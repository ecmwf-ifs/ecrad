!----------------------------------------------------------------------------
SUBROUTINE RRTM_TAUMOL15 (KIDIA,KFDIA,KLEV,taug,&
 & P_TAUAERL,fac00,fac01,fac10,fac11,forfac,forfrac,indfor,jp,jt,jt1,oneminus,&
 & colh2o,colco2,coln2o,laytrop,selffac,selffrac,indself,fracs, &
 & rat_n2oco2, rat_n2oco2_1,minorfrac,indminor,scaleminor,colbrd,laytrop_min,laytrop_max)

!     BAND 15:  2380-2600 cm-1 (low - N2O,CO2; high - nothing)

!     AUTHOR.
!     -------
!      JJMorcrette, ECMWF

!     MODIFICATIONS.
!     --------------
!      M.Hamrud      01-Oct-2003 CY28 Cleaning
!      NEC           25-Oct-2007 Optimisations
!      JJMorcrette 20110613 flexible number of g-points
!     ABozzo 2001306 updated to rrtmg v4.85
!     band 15:  2380-2600 cm-1 (low - n2o,co2; low minor - n2)
!                              (high - nothing)
! ---------------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM     ,JPRB

USE PARRRTM  , ONLY : JPBAND
USE YOERRTM  , ONLY : JPGPT  ,NGS14  ,NG15
USE YOERRTWN , ONLY : NSPA
USE YOERRTA15, ONLY : ABSA   ,KA_MN2,FRACREFA,SELFREF,FORREF
USE YOERRTRF, ONLY : CHI_MLS

IMPLICIT NONE

INTEGER(KIND=JPIM),INTENT(IN)    :: KIDIA
INTEGER(KIND=JPIM),INTENT(IN)    :: KFDIA
INTEGER(KIND=JPIM),INTENT(IN)    :: KLEV
REAL(KIND=JPRB)   ,INTENT(INOUT) :: taug(KIDIA:KFDIA,JPGPT,KLEV)
REAL(KIND=JPRB)   ,INTENT(IN)    :: P_TAUAERL(KIDIA:KFDIA,KLEV,JPBAND)
REAL(KIND=JPRB)   ,INTENT(IN)    :: fac00(KIDIA:KFDIA,KLEV)
REAL(KIND=JPRB)   ,INTENT(IN)    :: fac01(KIDIA:KFDIA,KLEV)
REAL(KIND=JPRB)   ,INTENT(IN)    :: fac10(KIDIA:KFDIA,KLEV)
REAL(KIND=JPRB)   ,INTENT(IN)    :: fac11(KIDIA:KFDIA,KLEV)
INTEGER(KIND=JPIM),INTENT(IN)    :: jp(KIDIA:KFDIA,KLEV)
INTEGER(KIND=JPIM),INTENT(IN)    :: jt(KIDIA:KFDIA,KLEV)
INTEGER(KIND=JPIM),INTENT(IN)    :: jt1(KIDIA:KFDIA,KLEV)
REAL(KIND=JPRB)   ,INTENT(IN)    :: oneminus
REAL(KIND=JPRB)   ,INTENT(IN)    :: colh2o(KIDIA:KFDIA,KLEV)
REAL(KIND=JPRB)   ,INTENT(IN)    :: colco2(KIDIA:KFDIA,KLEV)
REAL(KIND=JPRB)   ,INTENT(IN)    :: coln2o(KIDIA:KFDIA,KLEV)
INTEGER(KIND=JPIM),INTENT(IN)    :: laytrop(KIDIA:KFDIA)
REAL(KIND=JPRB)   ,INTENT(IN)    :: selffac(KIDIA:KFDIA,KLEV)
REAL(KIND=JPRB)   ,INTENT(IN)    :: selffrac(KIDIA:KFDIA,KLEV)
INTEGER(KIND=JPIM),INTENT(IN)    :: indself(KIDIA:KFDIA,KLEV)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: fracs(KIDIA:KFDIA,JPGPT,KLEV)

REAL(KIND=JPRB)   ,INTENT(IN)   :: rat_n2oco2(KIDIA:KFDIA,KLEV)
REAL(KIND=JPRB)   ,INTENT(IN)   :: rat_n2oco2_1(KIDIA:KFDIA,KLEV)
INTEGER(KIND=JPIM),INTENT(IN)   :: indfor(KIDIA:KFDIA,KLEV)
REAL(KIND=JPRB)   ,INTENT(IN)   :: forfac(KIDIA:KFDIA,KLEV)
REAL(KIND=JPRB)   ,INTENT(IN)   :: forfrac(KIDIA:KFDIA,KLEV)
REAL(KIND=JPRB)   ,INTENT(IN)   :: minorfrac(KIDIA:KFDIA,KLEV)
INTEGER(KIND=JPIM),INTENT(IN)   :: indminor(KIDIA:KFDIA,KLEV)
REAL(KIND=JPRB)   ,INTENT(IN)    :: scaleminor(KIDIA:KFDIA,KLEV)
REAL(KIND=JPRB)   ,INTENT(IN)    :: colbrd(KIDIA:KFDIA,KLEV)
INTEGER(KIND=JPIM), OPTIONAL, INTENT(INOUT) :: laytrop_min, laytrop_max
! ---------------------------------------------------------------------------

INTEGER(KIND=JPIM) :: IG, IND0, IND1, INDS,INDF,INDM, JS,JS1,JPL,JMN2, lay
REAL(KIND=JPRB) :: refrat_planck_a, refrat_m_a
REAL(KIND=JPRB) :: taufor,tauself,tau_major(ng15),tau_major1(ng15), n2m1, n2m2, taun2,scalen2
REAL(KIND=JPRB) ::  fac000, fac100, fac200,&
 & fac010, fac110, fac210, &
 & fac001, fac101, fac201, &
 & fac011, fac111, fac211
REAL(KIND=JPRB) :: p, p4, fk0, fk1, fk2

REAL(KIND=JPRB) :: fs, specmult, specparm,speccomb,  &
& fs1, specmult1, specparm1,speccomb1, &
& fmn2, specmult_mn2, specparm_mn2,speccomb_mn2, &
& fpl, specmult_planck, specparm_planck,speccomb_planck
    !     local integer arrays
    integer(KIND=JPIM) :: ixc(KLEV), ixlow(KFDIA,KLEV), ixhigh(KFDIA,KLEV)
    INTEGER(KIND=JPIM) :: ich, icl, ixc0, ixp, jc, jl
INTEGER(KIND=JPIM) :: llaytrop_min, llaytrop_max

#define MOD1(x) ((x) - AINT((x)))

#include "rrtm_utils.intfb.h"

    if (present(laytrop_min) .AND. present(laytrop_max)) then
       llaytrop_min = laytrop_min
       llaytrop_max = laytrop_max
    else
       CALL COMPUTE_LAYTROP_MIN_MAX(KIDIA, KFDIA, LAYTROP, llaytrop_min, llaytrop_max)
    endif

    !$ACC DATA PRESENT(taug, P_TAUAERL, fac00, fac01, fac10, fac11, jp, jt, jt1, &
    !$ACC             colh2o,  colco2,  coln2o,  laytrop, selffac, selffrac, &
    !$ACC             indself, fracs, rat_n2oco2, rat_n2oco2_1, indfor, forfac, &
    !$ACC             forfrac, minorfrac, indminor, scaleminor, colbrd)

#if !defined(_OPENACC) && !defined(OMPGPU)
    ixlow  = 0
    ixhigh = 0
    ixc    = 0

    ! create index lists for mixed layers
    do lay = llaytrop_min+1, llaytrop_max
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
#endif

      ! Minor gas mapping level :
      !     Lower - Nitrogen Continuum, P = 1053., T = 294.

      ! Calculate reference ratio to be used in calculation of Planck
      ! fraction in lower atmosphere.
      ! P = 1053. mb (Level 1)
      refrat_planck_a = chi_mls(4,1)/chi_mls(2,1)

      ! P = 1053.
      refrat_m_a = chi_mls(4,1)/chi_mls(2,1)

      ! Compute the optical depth by interpolating in ln(pressure),
      ! temperature, and appropriate species.  Below laytrop, the water
      ! vapor self-continuum and foreign continuum is interpolated
      ! (in temperature) separately.

      ! Lower atmosphere loop
      !$ACC WAIT
      !$ACC PARALLEL DEFAULT(NONE) ASYNC(1)
      !$ACC LOOP GANG VECTOR COLLAPSE(2) &
      !$ACC   PRIVATE(ind0, ind1, inds, indf, indm, &
      !$ACC           js, js1, jpl, jmn2, &
      !$ACC           fac000, fac100, fac200, &
      !$ACC           fac010, fac110, fac210, &
      !$ACC           fac001, fac101, fac201, &
      !$ACC           fac011, fac111, fac211, &
      !$ACC           p, p4, fk0, fk1, fk2, scalen2, &
      !$ACC           fs, specmult, specparm,speccomb,  &
      !$ACC           fs1, specmult1, specparm1,speccomb1, &
      !$ACC           fmn2, specmult_mn2, specparm_mn2,speccomb_mn2, &
      !$ACC           fpl, specmult_planck, &
      !$ACC           specparm_planck,speccomb_planck, &
      !$ACC           tau_major, tau_major1)
      do lay = 1, llaytrop_min
        do jl = KIDIA, KFDIA

          speccomb = coln2o(jl,lay) + rat_n2oco2(jl,lay)*colco2(jl,lay)
          specparm = MIN(coln2o(jl,lay)/speccomb,oneminus)
          specmult = 8._JPRB*(specparm)
          js = 1 + int(specmult)
          fs = MOD1(specmult)

          speccomb1 = coln2o(jl,lay) + rat_n2oco2_1(jl,lay)*colco2(jl,lay)
          specparm1 = MIN(coln2o(jl,lay)/speccomb1,oneminus)
          specmult1 = 8._JPRB*(specparm1)
          js1 = 1 + int(specmult1)
          fs1 = MOD1(specmult1)

          speccomb_mn2 = coln2o(jl,lay) + refrat_m_a*colco2(jl,lay)
          specparm_mn2 = MIN(coln2o(jl,lay)/speccomb_mn2,oneminus)
          specmult_mn2 = 8._JPRB*specparm_mn2
          jmn2 = 1 + int(specmult_mn2)
          fmn2 = MOD1(specmult_mn2)

          speccomb_planck = coln2o(jl,lay)+refrat_planck_a*colco2(jl,lay)
          specparm_planck = MIN(coln2o(jl,lay)/speccomb_planck,oneminus)
          specmult_planck = 8._JPRB*specparm_planck
          jpl = 1 + int(specmult_planck)
          fpl = MOD1(specmult_planck)

          ind0 = ((jp(jl,lay)-1)*5+(jt(jl,lay)-1))*nspa(15) + js
          ind1 = (jp(jl,lay)*5+(jt1(jl,lay)-1))*nspa(15) + js1
          inds = indself(jl,lay)
          indf = indfor(jl,lay)
          indm = indminor(jl,lay)

          scalen2 = colbrd(jl,lay)*scaleminor(jl,lay)

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
!$NEC unroll(NG15)
            tau_major(1:ng15) = speccomb *    &
             (fac000 * absa(ind0,1:ng15)    + &
              fac100 * absa(ind0+1,1:ng15)  + &
              fac200 * absa(ind0+2,1:ng15)  + &
              fac010 * absa(ind0+9,1:ng15)  + &
              fac110 * absa(ind0+10,1:ng15) + &
              fac210 * absa(ind0+11,1:ng15))
          else if (specparm .gt. 0.875_JPRB) then
!$NEC unroll(NG15)
            tau_major(1:ng15) = speccomb *   &
             (fac200 * absa(ind0-1,1:ng15) + &
              fac100 * absa(ind0,1:ng15)   + &
              fac000 * absa(ind0+1,1:ng15) + &
              fac210 * absa(ind0+8,1:ng15) + &
              fac110 * absa(ind0+9,1:ng15) + &
              fac010 * absa(ind0+10,1:ng15))
          else
!$NEC unroll(NG15)
            tau_major(1:ng15) = speccomb *   &
             (fac000 * absa(ind0,1:ng15)   + &
              fac100 * absa(ind0+1,1:ng15) + &
              fac010 * absa(ind0+9,1:ng15) + &
              fac110 * absa(ind0+10,1:ng15))
          endif

          if (specparm1 .lt. 0.125_JPRB) then
!$NEC unroll(NG15)
            tau_major1(1:ng15) = speccomb1 *  &
             (fac001 * absa(ind1,1:ng15)    + &
              fac101 * absa(ind1+1,1:ng15)  + &
              fac201 * absa(ind1+2,1:ng15)  + &
              fac011 * absa(ind1+9,1:ng15)  + &
              fac111 * absa(ind1+10,1:ng15) + &
              fac211 * absa(ind1+11,1:ng15))
          else if (specparm1 .gt. 0.875_JPRB) then
!$NEC unroll(NG15)
            tau_major1(1:ng15) = speccomb1 * &
             (fac201 * absa(ind1-1,1:ng15) + &
              fac101 * absa(ind1,1:ng15)   + &
              fac001 * absa(ind1+1,1:ng15) + &
              fac211 * absa(ind1+8,1:ng15) + &
              fac111 * absa(ind1+9,1:ng15) + &
              fac011 * absa(ind1+10,1:ng15))
          else
!$NEC unroll(NG15)
            tau_major1(1:ng15) = speccomb1 * &
             (fac001 * absa(ind1,1:ng15)   + &
              fac101 * absa(ind1+1,1:ng15) + &
              fac011 * absa(ind1+9,1:ng15) + &
              fac111 * absa(ind1+10,1:ng15))
          endif

          !$ACC LOOP SEQ PRIVATE(taufor, tauself, n2m1, n2m2, taun2)
!$NEC unroll(NG15)
          do ig = 1, ng15
            tauself = selffac(jl,lay)* (selfref(inds,ig) + selffrac(jl,lay) * &
                 (selfref(inds+1,ig) - selfref(inds,ig)))
            taufor =  forfac(jl,lay) * (forref(indf,ig) + forfrac(jl,lay) * &
                 (forref(indf+1,ig) - forref(indf,ig)))
            n2m1 = ka_mn2(jmn2,indm,ig) + fmn2 * &
                 (ka_mn2(jmn2+1,indm,ig) - ka_mn2(jmn2,indm,ig))
            n2m2 = ka_mn2(jmn2,indm+1,ig) + fmn2 * &
                 (ka_mn2(jmn2+1,indm+1,ig) - ka_mn2(jmn2,indm+1,ig))
            taun2 = scalen2 * (n2m1 + minorfrac(jl,lay) * (n2m2 - n2m1))

            taug(jl,ngs14+ig,lay) = tau_major(ig) + tau_major1(ig) &
                 + tauself + taufor &
                 + taun2
            fracs(jl,ngs14+ig,lay) = fracrefa(ig,jpl) + fpl * &
                 (fracrefa(ig,jpl+1)-fracrefa(ig,jpl))
          enddo
        enddo
      enddo
      !$ACC END PARALLEL

      ! Upper atmosphere loop
      !$ACC PARALLEL DEFAULT(NONE) ASYNC(1)
      !$ACC LOOP GANG VECTOR COLLAPSE(3)
      do ig = 1, ng15
        do lay = llaytrop_max+1, KLEV
          do jl = KIDIA, KFDIA
            taug(jl,ngs14+ig,lay) = 0.0_JPRB
            fracs(jl,ngs14+ig,lay) = 0.0_JPRB
          enddo
        enddo
      enddo
      !$ACC END PARALLEL

      IF (llaytrop_max /= llaytrop_min) THEN
        ! Mixed loop
        ! Lower atmosphere part
        !$ACC PARALLEL DEFAULT(NONE) ASYNC(1)
        !$ACC LOOP GANG VECTOR COLLAPSE(2) PRIVATE(speccomb, specparm, specmult, js, fs, speccomb1, specparm1, &
        !$ACC   specmult1, js1, fs1, speccomb_mn2, specparm_mn2, specmult_mn2, jmn2, fmn2, scalen2, speccomb_planck, &
        !$ACC   specparm_planck, specmult_planck, jpl, fpl, ind0, ind1, inds, indf, indm, p, p4, fk0, fk1, fk2, &
        !$ACC   fac000, fac100, fac200, fac010, fac110, fac210, fac001, fac101, fac201, fac011, fac111, fac211, &
        !$ACC   tau_major, tau_major1)
        do lay = llaytrop_min+1, llaytrop_max
#ifdef _OPENACC
          do jl = KIDIA, KFDIA
            if ( lay <= laytrop(jl) ) then
#else
          ixc0 = ixc(lay)

!$NEC ivdep
          do ixp = 1, ixc0
            jl = ixlow(ixp,lay)
#endif

            speccomb = coln2o(jl,lay) + rat_n2oco2(jl,lay)*colco2(jl,lay)
            specparm = MIN(coln2o(jl,lay)/speccomb,oneminus)
            specmult = 8._JPRB*(specparm)
            js = 1 + int(specmult)
            fs = MOD1(specmult)

            speccomb1 = coln2o(jl,lay) + rat_n2oco2_1(jl,lay)*colco2(jl,lay)
            specparm1 = MIN(coln2o(jl,lay)/speccomb1,oneminus)
            specmult1 = 8._JPRB*(specparm1)
            js1 = 1 + int(specmult1)
            fs1 = MOD1(specmult1)

            speccomb_mn2 = coln2o(jl,lay) + refrat_m_a*colco2(jl,lay)
            specparm_mn2 = MIN(coln2o(jl,lay)/speccomb_mn2,oneminus)
            specmult_mn2 = 8._JPRB*specparm_mn2
            jmn2 = 1 + int(specmult_mn2)
            fmn2 = MOD1(specmult_mn2)

            speccomb_planck = coln2o(jl,lay)+refrat_planck_a*colco2(jl,lay)
            specparm_planck = MIN(coln2o(jl,lay)/speccomb_planck,oneminus)
            specmult_planck = 8._JPRB*specparm_planck
            jpl = 1 + int(specmult_planck)
            fpl = MOD1(specmult_planck)

            ind0 = ((jp(jl,lay)-1)*5+(jt(jl,lay)-1))*nspa(15) + js
            ind1 = (jp(jl,lay)*5+(jt1(jl,lay)-1))*nspa(15) + js1
            inds = indself(jl,lay)
            indf = indfor(jl,lay)
            indm = indminor(jl,lay)

            scalen2 = colbrd(jl,lay)*scaleminor(jl,lay)

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
!$NEC unroll(NG15)
              tau_major(1:ng15) = speccomb *    &
              (fac000 * absa(ind0,1:ng15)    + &
                fac100 * absa(ind0+1,1:ng15)  + &
                fac200 * absa(ind0+2,1:ng15)  + &
                fac010 * absa(ind0+9,1:ng15)  + &
                fac110 * absa(ind0+10,1:ng15) + &
                fac210 * absa(ind0+11,1:ng15))
            else if (specparm .gt. 0.875_JPRB) then
!$NEC unroll(NG15)
              tau_major(1:ng15) = speccomb *   &
              (fac200 * absa(ind0-1,1:ng15) + &
                fac100 * absa(ind0,1:ng15)   + &
                fac000 * absa(ind0+1,1:ng15) + &
                fac210 * absa(ind0+8,1:ng15) + &
                fac110 * absa(ind0+9,1:ng15) + &
                fac010 * absa(ind0+10,1:ng15))
            else
!$NEC unroll(NG15)
              tau_major(1:ng15) = speccomb *   &
              (fac000 * absa(ind0,1:ng15)   + &
                fac100 * absa(ind0+1,1:ng15) + &
                fac010 * absa(ind0+9,1:ng15) + &
                fac110 * absa(ind0+10,1:ng15))
            endif

            if (specparm1 .lt. 0.125_JPRB) then
!$NEC unroll(NG15)
              tau_major1(1:ng15) = speccomb1 *  &
              (fac001 * absa(ind1,1:ng15)    + &
                fac101 * absa(ind1+1,1:ng15)  + &
                fac201 * absa(ind1+2,1:ng15)  + &
                fac011 * absa(ind1+9,1:ng15)  + &
                fac111 * absa(ind1+10,1:ng15) + &
                fac211 * absa(ind1+11,1:ng15))
            else if (specparm1 .gt. 0.875_JPRB) then
!$NEC unroll(NG15)
              tau_major1(1:ng15) = speccomb1 * &
              (fac201 * absa(ind1-1,1:ng15) + &
                fac101 * absa(ind1,1:ng15)   + &
                fac001 * absa(ind1+1,1:ng15) + &
                fac211 * absa(ind1+8,1:ng15) + &
                fac111 * absa(ind1+9,1:ng15) + &
                fac011 * absa(ind1+10,1:ng15))
            else
!$NEC unroll(NG15)
              tau_major1(1:ng15) = speccomb1 * &
              (fac001 * absa(ind1,1:ng15)   + &
                fac101 * absa(ind1+1,1:ng15) + &
                fac011 * absa(ind1+9,1:ng15) + &
                fac111 * absa(ind1+10,1:ng15))
            endif

!$NEC unroll(NG15)
            !$ACC LOOP SEQ PRIVATE(tauself, taufor, n2m1, n2m2, taun2)
            do ig = 1, ng15
              tauself = selffac(jl,lay)* (selfref(inds,ig) + selffrac(jl,lay) * &
                  (selfref(inds+1,ig) - selfref(inds,ig)))
              taufor =  forfac(jl,lay) * (forref(indf,ig) + forfrac(jl,lay) * &
                  (forref(indf+1,ig) - forref(indf,ig)))
              n2m1 = ka_mn2(jmn2,indm,ig) + fmn2 * &
                  (ka_mn2(jmn2+1,indm,ig) - ka_mn2(jmn2,indm,ig))
              n2m2 = ka_mn2(jmn2,indm+1,ig) + fmn2 * &
                  (ka_mn2(jmn2+1,indm+1,ig) - ka_mn2(jmn2,indm+1,ig))
              taun2 = scalen2 * (n2m1 + minorfrac(jl,lay) * (n2m2 - n2m1))

              taug(jl,ngs14+ig,lay) = tau_major(ig) + tau_major1(ig) &
                  + tauself + taufor &
                  + taun2
              fracs(jl,ngs14+ig,lay) = fracrefa(ig,jpl) + fpl * &
                  (fracrefa(ig,jpl+1)-fracrefa(ig,jpl))
            enddo
#ifdef _OPENACC
         else
#else
          enddo

          ! Upper atmosphere part
          ixc0 = KFDIA - KIDIA + 1 - ixc0
#endif

          !$ACC LOOP SEQ
          do ig = 1, ng15
#ifndef _OPENACC
!$NEC ivdep
            do ixp = 1, ixc0
              jl = ixhigh(ixp,lay)
#endif

              taug(jl,ngs14+ig,lay) = 0.0_JPRB
              fracs(jl,ngs14+ig,lay) = 0.0_JPRB
            enddo
#ifdef _OPENACC
           endif
#endif
          enddo

        enddo
        !$ACC END PARALLEL

      ENDIF

      !$ACC END DATA

END SUBROUTINE RRTM_TAUMOL15
