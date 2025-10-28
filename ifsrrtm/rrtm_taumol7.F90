!----------------------------------------------------------------------------
SUBROUTINE RRTM_TAUMOL7 (KIDIA,KFDIA,KLEV,taug,&
 & P_TAUAERL,fac00,fac01,fac10,fac11,forfac,forfrac,indfor,jp,jt,jt1,oneminus,&
 & colh2o,colo3,colco2,coldry,laytrop,selffac,selffrac,indself,fracs, &
 & rat_h2oo3, rat_h2oo3_1,minorfrac,indminor,laytrop_min,laytrop_max)

!     BAND 7:  980-1080 cm-1 (low - H2O,O3; high - O3)

!     AUTHOR.
!     -------
!      JJMorcrette, ECMWF

!     MODIFICATIONS.
!     --------------
!      M.Hamrud      01-Oct-2003 CY28 Cleaning
!      NEC           25-Oct-2007 Optimisations
!      JJMorcrette 20110613 flexible number of g-points
!      ABozzo 201306 updated to rrtmg v4.85
!     band 7:  980-1080 cm-1 (low key - h2o,o3; low minor - co2)
!                            (high key - o3; high minor - co2)
! ---------------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM     ,JPRB

USE PARRRTM  , ONLY : JPBAND
USE YOERRTM  , ONLY : JPGPT  ,NG7   ,NGS6
USE YOERRTWN , ONLY : NSPA   ,NSPB
USE YOERRTA7 , ONLY : ABSA   ,ABSB   ,KA_MCO2,KB_MCO2 ,FRACREFA ,FRACREFB,SELFREF,FORREF
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
REAL(KIND=JPRB)   ,INTENT(IN)    :: colo3(KIDIA:KFDIA,KLEV)
REAL(KIND=JPRB)   ,INTENT(IN)    :: colco2(KIDIA:KFDIA,KLEV)
REAL(KIND=JPRB)   ,INTENT(IN)    :: coldry(KIDIA:KFDIA,KLEV)
INTEGER(KIND=JPIM),INTENT(IN)    :: laytrop(KIDIA:KFDIA)
REAL(KIND=JPRB)   ,INTENT(IN)    :: selffac(KIDIA:KFDIA,KLEV)
REAL(KIND=JPRB)   ,INTENT(IN)    :: selffrac(KIDIA:KFDIA,KLEV)
INTEGER(KIND=JPIM),INTENT(IN)    :: indself(KIDIA:KFDIA,KLEV)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: fracs(KIDIA:KFDIA,JPGPT,KLEV)

REAL(KIND=JPRB)   ,INTENT(IN)   :: rat_h2oo3(KIDIA:KFDIA,KLEV)
REAL(KIND=JPRB)   ,INTENT(IN)   :: rat_h2oo3_1(KIDIA:KFDIA,KLEV)
INTEGER(KIND=JPIM),INTENT(IN)   :: indfor(KIDIA:KFDIA,KLEV)
REAL(KIND=JPRB)   ,INTENT(IN)   :: forfrac(KIDIA:KFDIA,KLEV)
REAL(KIND=JPRB)   ,INTENT(IN)   :: forfac(KIDIA:KFDIA,KLEV)
REAL(KIND=JPRB)   ,INTENT(IN)   :: minorfrac(KIDIA:KFDIA,KLEV)
INTEGER(KIND=JPIM),INTENT(IN)   :: indminor(KIDIA:KFDIA,KLEV)
INTEGER(KIND=JPIM), OPTIONAL, INTENT(INOUT) :: laytrop_min, laytrop_max

! ---------------------------------------------------------------------------

REAL(KIND=JPRB) :: speccomb,speccomb1, &
& speccomb_mco2, speccomb_planck
INTEGER(KIND=JPIM) :: ind0,ind1,inds,indf,indm

INTEGER(KIND=JPIM) :: IG, JS, lay, JS1, JPL, JMCO2

REAL(KIND=JPRB) :: refrat_planck_a, refrat_m_a
REAL(KIND=JPRB) :: chi_co2, ratco2, adjfac, adjcolco2
REAL(KIND=JPRB) ::  fac000, fac100, fac200,&
 & fac010, fac110, fac210, &
 & fac001, fac101, fac201, &
 & fac011, fac111, fac211
REAL(KIND=JPRB) :: p, p4, fk0, fk1, fk2
REAL(KIND=JPRB) :: taufor,tauself,tau_major(ng7),tau_major1(ng7), co2m1, co2m2, absco2


REAL(KIND=JPRB) :: fs, specmult, specparm,  &
& fs1, specmult1, specparm1, &
& fpl, specmult_PLANCK, specparm_PLANCK, &
& fmco2, specmult_mco2, specparm_mco2
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
    !$ACC             colh2o, colo3, colco2, coldry, laytrop, selffac, selffrac, indself, fracs, &
    !$ACC             rat_h2oo3, rat_h2oo3_1, indfor, forfrac, forfac, &
    !$ACC             minorfrac, indminor)
#ifndef __NVCOMPILER
    !$OMP TARGET DATA MAP(PRESENT, ALLOC:taug, P_TAUAERL, fac00, fac01, fac10, fac11, jp, jt, jt1, &
    !$OMP             colh2o, colo3, colco2, coldry, laytrop, selffac, selffrac, indself, fracs, &
    !$OMP             rat_h2oo3, rat_h2oo3_1, indfor, forfrac, forfac, &
    !$OMP             minorfrac, indminor)
#endif

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


      ! P = 706.2620 mb
      refrat_planck_a = chi_mls(1,3)/chi_mls(3,3)

      ! P = 706.2720 mb
      refrat_m_a = chi_mls(1,3)/chi_mls(3,3)

      ! Compute the optical depth by interpolating in ln(pressure),
      ! temperature, and appropriate species.  Below laytrop, the water
      ! vapor self-continuum and foreign continuum is interpolated
      ! (in temperature) separately.

      ! Lower atmosphere loop
      !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO COLLAPSE(2) PRIVATE(speccomb, speccomb1, speccomb_mco2, speccomb_planck, ind0, ind1, &
      !$OMP   inds, indf, indm, js, js1, jpl, jmco2, chi_co2, ratco2, adjfac, adjcolco2, fac000, fac100, fac200, &
      !$OMP   fac010, fac110, fac210, fac001, fac101, fac201, fac011, fac111, fac211, p, p4, fk0, fk1, fk2, fs, &
      !$OMP   specmult, specparm, fs1, specmult1, specparm1, fpl, specmult_PLANCK, specparm_PLANCK, fmco2, &
      !$OMP   specmult_mco2, specparm_mco2, tau_major, tau_major1, taufor, tauself, co2m1, co2m2, absco2) THREAD_LIMIT(128)
      !$ACC WAIT
      !$ACC PARALLEL DEFAULT(NONE) ASYNC(1)
      !$ACC LOOP GANG VECTOR COLLAPSE(2) PRIVATE(speccomb, speccomb1, speccomb_mco2, speccomb_planck, ind0, ind1, &
      !$ACC   inds, indf, indm, js, js1, jpl, jmco2, chi_co2, ratco2, adjfac, adjcolco2, fac000, fac100, fac200, &
      !$ACC   fac010, fac110, fac210, fac001, fac101, fac201, fac011, fac111, fac211, p, p4, fk0, fk1, fk2, fs, &
      !$ACC   specmult, specparm, fs1, specmult1, specparm1, fpl, specmult_PLANCK, specparm_PLANCK, fmco2, &
      !$ACC   specmult_mco2, specparm_mco2, tau_major, tau_major1)
      do lay = 1, llaytrop_min
        do jl = KIDIA, KFDIA

          speccomb = colh2o(jl,lay) + rat_h2oo3(jl,lay)*colo3(jl,lay)
          specparm = MIN(colh2o(jl,lay)/speccomb,oneminus)
          specmult = 8._JPRB*(specparm)
          js = 1 + int(specmult)
          fs = MOD1(specmult)

          speccomb1 = colh2o(jl,lay) + rat_h2oo3_1(jl,lay)*colo3(jl,lay)
          specparm1 = MIN(colh2o(jl,lay)/speccomb1,oneminus)
          specmult1 = 8._JPRB*(specparm1)
          js1 = 1 + int(specmult1)
          fs1 = MOD1(specmult1)

          speccomb_mco2 = colh2o(jl,lay) + refrat_m_a*colo3(jl,lay)
          specparm_mco2 = MIN(colh2o(jl,lay)/speccomb_mco2,oneminus)
          specmult_mco2 = 8._JPRB*specparm_mco2

          jmco2 = 1 + int(specmult_mco2)
          fmco2 = MOD1(specmult_mco2)

          !  In atmospheres where the amount of CO2 is too great to be considered
          !  a minor species, adjust the column amount of CO2 by an empirical factor
          !  to obtain the proper contribution.
          chi_co2 = colco2(jl,lay)/(coldry(jl,lay))
          ratco2 = 1.e20_JPRB*chi_co2/chi_mls(2,jp(jl,lay)+1)
          if (ratco2 .gt. 3.0_JPRB) then
            adjfac = 3.0_JPRB+(ratco2-3.0_JPRB)**0.79_JPRB
            adjcolco2 = adjfac*chi_mls(2,jp(jl,lay)+1)*coldry(jl,lay)*1.e-20_JPRB
          else
            adjcolco2 = colco2(jl,lay)
          endif

          speccomb_planck = colh2o(jl,lay)+refrat_planck_a*colo3(jl,lay)
          specparm_planck = MIN(colh2o(jl,lay)/speccomb_planck,oneminus)
          specmult_planck = 8._JPRB*specparm_planck
          jpl = 1 + int(specmult_planck)
          fpl = MOD1(specmult_planck)

          ind0 = ((jp(jl,lay)-1)*5+(jt(jl,lay)-1))*nspa(7) + js
          ind1 = (jp(jl,lay)*5+(jt1(jl,lay)-1))*nspa(7) + js1
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
#if defined(__amdflang__) && defined(OMPGPU)
            do ig = 1, ng7
            tau_major(ig) = speccomb *    &
             (fac000 * absa(ind0,ig)    + &
              fac100 * absa(ind0+1,ig)  + &
              fac200 * absa(ind0+2,ig)  + &
              fac010 * absa(ind0+9,ig)  + &
              fac110 * absa(ind0+10,ig) + &
              fac210 * absa(ind0+11,ig))
            end do
#else
!$NEC unroll(NG7)
            tau_major(1:ng7) = speccomb *    &
             (fac000 * absa(ind0,1:ng7)    + &
              fac100 * absa(ind0+1,1:ng7)  + &
              fac200 * absa(ind0+2,1:ng7)  + &
              fac010 * absa(ind0+9,1:ng7)  + &
              fac110 * absa(ind0+10,1:ng7) + &
              fac210 * absa(ind0+11,1:ng7))
#endif
          else if (specparm .gt. 0.875_JPRB) then
#if defined(__amdflang__) && defined(OMPGPU)
            do ig = 1, ng7
            tau_major(ig) = speccomb *   &
             (fac200 * absa(ind0-1,ig) + &
              fac100 * absa(ind0,ig)   + &
              fac000 * absa(ind0+1,ig) + &
              fac210 * absa(ind0+8,ig) + &
              fac110 * absa(ind0+9,ig) + &
              fac010 * absa(ind0+10,ig))
            end do
#else
!$NEC unroll(NG7)
            tau_major(1:ng7) = speccomb *   &
             (fac200 * absa(ind0-1,1:ng7) + &
              fac100 * absa(ind0,1:ng7)   + &
              fac000 * absa(ind0+1,1:ng7) + &
              fac210 * absa(ind0+8,1:ng7) + &
              fac110 * absa(ind0+9,1:ng7) + &
              fac010 * absa(ind0+10,1:ng7))
#endif
          else
#if defined(__amdflang__) && defined(OMPGPU)
            do ig = 1, ng7
            tau_major(ig) = speccomb *   &
             (fac000 * absa(ind0,ig)   + &
              fac100 * absa(ind0+1,ig) + &
              fac010 * absa(ind0+9,ig) + &
              fac110 * absa(ind0+10,ig))
            end do
#else
!$NEC unroll(NG7)
            tau_major(1:ng7) = speccomb *   &
             (fac000 * absa(ind0,1:ng7)   + &
              fac100 * absa(ind0+1,1:ng7) + &
              fac010 * absa(ind0+9,1:ng7) + &
              fac110 * absa(ind0+10,1:ng7))
#endif
          endif

          if (specparm1 .lt. 0.125_JPRB) then
#if defined(__amdflang__) && defined(OMPGPU)
            do ig = 1, ng7
            tau_major1(ig) = speccomb1 *  &
             (fac001 * absa(ind1,ig)    + &
              fac101 * absa(ind1+1,ig)  + &
              fac201 * absa(ind1+2,ig)  + &
              fac011 * absa(ind1+9,ig)  + &
              fac111 * absa(ind1+10,ig) + &
              fac211 * absa(ind1+11,ig))
            end do
#else
!$NEC unroll(NG7)
            tau_major1(1:ng7) = speccomb1 *  &
             (fac001 * absa(ind1,1:ng7)    + &
              fac101 * absa(ind1+1,1:ng7)  + &
              fac201 * absa(ind1+2,1:ng7)  + &
              fac011 * absa(ind1+9,1:ng7)  + &
              fac111 * absa(ind1+10,1:ng7) + &
              fac211 * absa(ind1+11,1:ng7))
#endif
          else if (specparm1 .gt. 0.875_JPRB) then
#if defined(__amdflang__) && defined(OMPGPU)
            do ig = 1, ng7
            tau_major1(ig) = speccomb1 * &
             (fac201 * absa(ind1-1,ig) + &
              fac101 * absa(ind1,ig)   + &
              fac001 * absa(ind1+1,ig) + &
              fac211 * absa(ind1+8,ig) + &
              fac111 * absa(ind1+9,ig) + &
              fac011 * absa(ind1+10,ig))
            end do
#else
!$NEC unroll(NG7)
            tau_major1(1:ng7) = speccomb1 * &
             (fac201 * absa(ind1-1,1:ng7) + &
              fac101 * absa(ind1,1:ng7)   + &
              fac001 * absa(ind1+1,1:ng7) + &
              fac211 * absa(ind1+8,1:ng7) + &
              fac111 * absa(ind1+9,1:ng7) + &
              fac011 * absa(ind1+10,1:ng7))
#endif
          else
#if defined(__amdflang__) && defined(OMPGPU)
            do ig = 1, ng7
            tau_major1(ig) = speccomb1 * &
             (fac001 * absa(ind1,ig)   + &
              fac101 * absa(ind1+1,ig) + &
              fac011 * absa(ind1+9,ig) + &
              fac111 * absa(ind1+10,ig))
            end do
#else
!$NEC unroll(NG7)
            tau_major1(1:ng7) = speccomb1 * &
             (fac001 * absa(ind1,1:ng7)   + &
              fac101 * absa(ind1+1,1:ng7) + &
              fac011 * absa(ind1+9,1:ng7) + &
              fac111 * absa(ind1+10,1:ng7))
#endif
          endif

          !$ACC LOOP SEQ PRIVATE(taufor, tauself, co2m1, co2m2, absco2)
!$NEC unroll(NG7)
          do ig = 1, ng7
            tauself = selffac(jl,lay)* (selfref(inds,ig) + selffrac(jl,lay) * &
                 (selfref(inds+1,ig) - selfref(inds,ig)))
            taufor = forfac(jl,lay) * (forref(indf,ig) + forfrac(jl,lay) * &
                 (forref(indf+1,ig) - forref(indf,ig)))
            co2m1 = ka_mco2(jmco2,indm,ig) + fmco2 * &
                 (ka_mco2(jmco2+1,indm,ig) - ka_mco2(jmco2,indm,ig))
            co2m2 = ka_mco2(jmco2,indm+1,ig) + fmco2 * &
                 (ka_mco2(jmco2+1,indm+1,ig) - ka_mco2(jmco2,indm+1,ig))
            absco2 = co2m1 + minorfrac(jl,lay) * (co2m2 - co2m1)

            taug(jl,ngs6+ig,lay) = tau_major(ig) + tau_major1(ig) &
                 + tauself + taufor &
                 + adjcolco2*absco2
            fracs(jl,ngs6+ig,lay) = fracrefa(ig,jpl) + fpl * &
                 (fracrefa(ig,jpl+1)-fracrefa(ig,jpl))
          enddo
        enddo

      enddo
      !$ACC END PARALLEL
      !$OMP END TARGET TEAMS DISTRIBUTE PARALLEL DO

      ! Upper atmosphere loop
      !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO COLLAPSE(2) PRIVATE(ind0, ind1, indm, chi_co2, ratco2, adjfac, adjcolco2, absco2)
      !$ACC PARALLEL DEFAULT(NONE) ASYNC(1)
      !$ACC LOOP GANG VECTOR COLLAPSE(2) PRIVATE(ind0, ind1, indm, chi_co2, ratco2, adjfac, adjcolco2)
      do lay = llaytrop_max+1, KLEV
        do jl = KIDIA, KFDIA

          !  In atmospheres where the amount of CO2 is too great to be considered
          !  a minor species, adjust the column amount of CO2 by an empirical factor
          !  to obtain the proper contribution.
          chi_co2 = colco2(jl,lay)/(coldry(jl,lay))
          ratco2 = 1.e20_JPRB*chi_co2/chi_mls(2,jp(jl,lay)+1)
          if (ratco2 .gt. 3.0_JPRB) then
            adjfac = 2.0_JPRB+(ratco2-2.0_JPRB)**0.79_JPRB
            adjcolco2 = adjfac*chi_mls(2,jp(jl,lay)+1)*coldry(jl,lay)*1.e-20_JPRB
          else
            adjcolco2 = colco2(jl,lay)
          endif

          ind0 = ((jp(jl,lay)-13)*5+(jt(jl,lay)-1))*nspb(7) + 1
          ind1 = ((jp(jl,lay)-12)*5+(jt1(jl,lay)-1))*nspb(7) + 1
          indm = indminor(jl,lay)
          !$ACC LOOP SEQ PRIVATE(absco2)
!$NEC unroll(NG7)
          do ig = 1, ng7
            absco2 = kb_mco2(indm,ig) + minorfrac(jl,lay) * &
                 (kb_mco2(indm+1,ig) - kb_mco2(indm,ig))
            taug(jl,ngs6+ig,lay) = colo3(jl,lay) * &
                 (fac00(jl,lay) * absb(ind0,ig) + &
                 fac10(jl,lay) * absb(ind0+1,ig) + &
                 fac01(jl,lay) * absb(ind1,ig) + &
                 fac11(jl,lay) * absb(ind1+1,ig)) &
                 + adjcolco2 * absco2
            fracs(jl,ngs6+ig,lay) = fracrefb(ig)
          enddo
        enddo
      enddo
      !$ACC END PARALLEL
      !$ACC WAIT
      !$OMP END TARGET TEAMS DISTRIBUTE PARALLEL DO

      ! Empirical modification to code to improve stratospheric cooling rates
      ! for o3.  Revised to apply weighting for g-point reduction in this band.
      !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO COLLAPSE(2)
      !$ACC PARALLEL DEFAULT(NONE) ASYNC(1)
      !$ACC LOOP GANG VECTOR COLLAPSE(2)
      do lay = llaytrop_max+1, KLEV
        do jl = KIDIA, KFDIA
          taug(jl,ngs6+6,lay)=taug(jl,ngs6+6,lay)*0.92_JPRB
          taug(jl,ngs6+7,lay)=taug(jl,ngs6+7,lay)*0.88_JPRB
          taug(jl,ngs6+8,lay)=taug(jl,ngs6+8,lay)*1.07_JPRB
          taug(jl,ngs6+9,lay)=taug(jl,ngs6+9,lay)*1.1_JPRB
          taug(jl,ngs6+10,lay)=taug(jl,ngs6+10,lay)*0.99_JPRB
          taug(jl,ngs6+11,lay)=taug(jl,ngs6+11,lay)*0.855_JPRB
        enddo
      enddo
      !$ACC END PARALLEL
      !$OMP END TARGET TEAMS DISTRIBUTE PARALLEL DO

      IF (llaytrop_max /= llaytrop_min) THEN
        ! Mixed loop
        ! Lower atmosphere part
        !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO COLLAPSE(2) PRIVATE(speccomb, specparm, specmult, js, fs, speccomb1, specparm1, &
        !$OMP   specmult1, js1, fs1, speccomb_mco2, specparm_mco2, specmult_mco2, jmco2, fmco2, adjfac, chi_co2, &
        !$OMP   ratco2, adjcolco2, speccomb_planck, specparm_planck, specmult_planck, jpl, fpl, ind0, ind1, inds, &
        !$OMP   indf, indm, p, p4, fk0, fk1, fk2, fac000, fac100, fac200, fac010, fac110, fac210, fac001, fac101, &
        !$OMP   fac201, fac011, fac111, fac211, tau_major, tau_major1, tauself, taufor, co2m1, co2m2, absco2)
        !$ACC PARALLEL DEFAULT(NONE) ASYNC(1)
        !$ACC LOOP GANG VECTOR COLLAPSE(2) PRIVATE(speccomb, specparm, specmult, js, fs, speccomb1, specparm1, &
        !$ACC   specmult1, js1, fs1, speccomb_mco2, specparm_mco2, specmult_mco2, jmco2, fmco2, adjfac, chi_co2, &
        !$ACC   ratco2, adjcolco2, speccomb_planck, specparm_planck, specmult_planck, jpl, fpl, ind0, ind1, inds, &
        !$ACC   indf, indm, p, p4, fk0, fk1, fk2, fac000, fac100, fac200, fac010, fac110, fac210, fac001, fac101, &
        !$ACC   fac201, fac011, fac111, fac211, tau_major, tau_major1)
        do lay = llaytrop_min+1, llaytrop_max
#if defined(_OPENACC) || defined(OMPGPU)
          do jl = KIDIA, KFDIA
            if ( lay <= laytrop(jl) ) then
#else
          ixc0 = ixc(lay)

!$NEC ivdep
          do ixp = 1, ixc0
            jl = ixlow(ixp,lay)
#endif

            speccomb = colh2o(jl,lay) + rat_h2oo3(jl,lay)*colo3(jl,lay)
            specparm = MIN(colh2o(jl,lay)/speccomb,oneminus)
            specmult = 8._JPRB*(specparm)
            js = 1 + int(specmult)
            fs = MOD1(specmult)

            speccomb1 = colh2o(jl,lay) + rat_h2oo3_1(jl,lay)*colo3(jl,lay)
            specparm1 = MIN(colh2o(jl,lay)/speccomb1,oneminus)
            specmult1 = 8._JPRB*(specparm1)
            js1 = 1 + int(specmult1)
            fs1 = MOD1(specmult1)

            speccomb_mco2 = colh2o(jl,lay) + refrat_m_a*colo3(jl,lay)
            specparm_mco2 = MIN(colh2o(jl,lay)/speccomb_mco2,oneminus)
            specmult_mco2 = 8._JPRB*specparm_mco2

            jmco2 = 1 + int(specmult_mco2)
            fmco2 = MOD1(specmult_mco2)

            !  In atmospheres where the amount of CO2 is too great to be considered
            !  a minor species, adjust the column amount of CO2 by an empirical factor
            !  to obtain the proper contribution.
            chi_co2 = colco2(jl,lay)/(coldry(jl,lay))
            ratco2 = 1.e20_JPRB*chi_co2/chi_mls(2,jp(jl,lay)+1)
            if (ratco2 .gt. 3.0_JPRB) then
              adjfac = 3.0_JPRB+(ratco2-3.0_JPRB)**0.79_JPRB
              adjcolco2 = adjfac*chi_mls(2,jp(jl,lay)+1)*coldry(jl,lay)*1.e-20_JPRB
            else
              adjcolco2 = colco2(jl,lay)
            endif

            speccomb_planck = colh2o(jl,lay)+refrat_planck_a*colo3(jl,lay)
            specparm_planck = MIN(colh2o(jl,lay)/speccomb_planck,oneminus)
            specmult_planck = 8._JPRB*specparm_planck
            jpl = 1 + int(specmult_planck)
            fpl = MOD1(specmult_planck)

            ind0 = ((jp(jl,lay)-1)*5+(jt(jl,lay)-1))*nspa(7) + js
            ind1 = (jp(jl,lay)*5+(jt1(jl,lay)-1))*nspa(7) + js1
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
#if defined(__amdflang__) && defined(OMPGPU)
            do ig = 1, ng7
              tau_major(ig) = speccomb *    &
              (fac000 * absa(ind0,ig)    + &
                fac100 * absa(ind0+1,ig)  + &
                fac200 * absa(ind0+2,ig)  + &
                fac010 * absa(ind0+9,ig)  + &
                fac110 * absa(ind0+10,ig) + &
                fac210 * absa(ind0+11,ig))
            end do
#else
!$NEC unroll(NG7)
              tau_major(1:ng7) = speccomb *    &
              (fac000 * absa(ind0,1:ng7)    + &
                fac100 * absa(ind0+1,1:ng7)  + &
                fac200 * absa(ind0+2,1:ng7)  + &
                fac010 * absa(ind0+9,1:ng7)  + &
                fac110 * absa(ind0+10,1:ng7) + &
                fac210 * absa(ind0+11,1:ng7))
#endif
            else if (specparm .gt. 0.875_JPRB) then
#if defined(__amdflang__) && defined(OMPGPU)
            do ig = 1, ng7
              tau_major(ig) = speccomb *   &
              (fac200 * absa(ind0-1,ig) + &
                fac100 * absa(ind0,ig)   + &
                fac000 * absa(ind0+1,ig) + &
                fac210 * absa(ind0+8,ig) + &
                fac110 * absa(ind0+9,ig) + &
                fac010 * absa(ind0+10,ig))
            end do
#else
!$NEC unroll(NG7)
              tau_major(1:ng7) = speccomb *   &
              (fac200 * absa(ind0-1,1:ng7) + &
                fac100 * absa(ind0,1:ng7)   + &
                fac000 * absa(ind0+1,1:ng7) + &
                fac210 * absa(ind0+8,1:ng7) + &
                fac110 * absa(ind0+9,1:ng7) + &
                fac010 * absa(ind0+10,1:ng7))
#endif
            else
#if defined(__amdflang__) && defined(OMPGPU)
            do ig = 1, ng7
              tau_major(ig) = speccomb *   &
              (fac000 * absa(ind0,ig)   + &
                fac100 * absa(ind0+1,ig) + &
                fac010 * absa(ind0+9,ig) + &
                fac110 * absa(ind0+10,ig))
            end do
#else
!$NEC unroll(NG7)
              tau_major(1:ng7) = speccomb *   &
              (fac000 * absa(ind0,1:ng7)   + &
                fac100 * absa(ind0+1,1:ng7) + &
                fac010 * absa(ind0+9,1:ng7) + &
                fac110 * absa(ind0+10,1:ng7))
#endif
            endif

            if (specparm1 .lt. 0.125_JPRB) then
#if defined(__amdflang__) && defined(OMPGPU)
            do ig = 1, ng7
              tau_major1(ig) = speccomb1 *  &
              (fac001 * absa(ind1,ig)    + &
                fac101 * absa(ind1+1,ig)  + &
                fac201 * absa(ind1+2,ig)  + &
                fac011 * absa(ind1+9,ig)  + &
                fac111 * absa(ind1+10,ig) + &
                fac211 * absa(ind1+11,ig))
            end do
#else
!$NEC unroll(NG7)
              tau_major1(1:ng7) = speccomb1 *  &
              (fac001 * absa(ind1,1:ng7)    + &
                fac101 * absa(ind1+1,1:ng7)  + &
                fac201 * absa(ind1+2,1:ng7)  + &
                fac011 * absa(ind1+9,1:ng7)  + &
                fac111 * absa(ind1+10,1:ng7) + &
                fac211 * absa(ind1+11,1:ng7))
#endif
            else if (specparm1 .gt. 0.875_JPRB) then
#if defined(__amdflang__) && defined(OMPGPU)
            do ig = 1, ng7
              tau_major1(ig) = speccomb1 * &
              (fac201 * absa(ind1-1,ig) + &
                fac101 * absa(ind1,ig)   + &
                fac001 * absa(ind1+1,ig) + &
                fac211 * absa(ind1+8,ig) + &
                fac111 * absa(ind1+9,ig) + &
                fac011 * absa(ind1+10,ig))
            end do
#else
!$NEC unroll(NG7)
              tau_major1(1:ng7) = speccomb1 * &
              (fac201 * absa(ind1-1,1:ng7) + &
                fac101 * absa(ind1,1:ng7)   + &
                fac001 * absa(ind1+1,1:ng7) + &
                fac211 * absa(ind1+8,1:ng7) + &
                fac111 * absa(ind1+9,1:ng7) + &
                fac011 * absa(ind1+10,1:ng7))
#endif
            else
#if defined(__amdflang__) && defined(OMPGPU)
            do ig = 1, ng7
              tau_major1(ig) = speccomb1 * &
              (fac001 * absa(ind1,ig)   + &
                fac101 * absa(ind1+1,ig) + &
                fac011 * absa(ind1+9,ig) + &
                fac111 * absa(ind1+10,ig))
            end do
#else
!$NEC unroll(NG7)
              tau_major1(1:ng7) = speccomb1 * &
              (fac001 * absa(ind1,1:ng7)   + &
                fac101 * absa(ind1+1,1:ng7) + &
                fac011 * absa(ind1+9,1:ng7) + &
                fac111 * absa(ind1+10,1:ng7))
#endif
            endif

!$NEC unroll(NG7)
            !$ACC LOOP SEQ PRIVATE(tauself, taufor, co2m1, co2m2, absco2)
            do ig = 1, ng7
              tauself = selffac(jl,lay)* (selfref(inds,ig) + selffrac(jl,lay) * &
                  (selfref(inds+1,ig) - selfref(inds,ig)))
              taufor = forfac(jl,lay) * (forref(indf,ig) + forfrac(jl,lay) * &
                  (forref(indf+1,ig) - forref(indf,ig)))
              co2m1 = ka_mco2(jmco2,indm,ig) + fmco2 * &
                  (ka_mco2(jmco2+1,indm,ig) - ka_mco2(jmco2,indm,ig))
              co2m2 = ka_mco2(jmco2,indm+1,ig) + fmco2 * &
                  (ka_mco2(jmco2+1,indm+1,ig) - ka_mco2(jmco2,indm+1,ig))
              absco2 = co2m1 + minorfrac(jl,lay) * (co2m2 - co2m1)

              taug(jl,ngs6+ig,lay) = tau_major(ig) + tau_major1(ig) &
                  + tauself + taufor &
                  + adjcolco2*absco2
              fracs(jl,ngs6+ig,lay) = fracrefa(ig,jpl) + fpl * &
                  (fracrefa(ig,jpl+1)-fracrefa(ig,jpl))
            enddo
#if defined(_OPENACC) || defined(OMPGPU)
         else
#else
          enddo

          ! Upper atmosphere part
          ixc0 = KFDIA - KIDIA + 1 - ixc0
!$NEC ivdep
          do ixp = 1, ixc0
            jl = ixhigh(ixp,lay)
#endif

            !  In atmospheres where the amount of CO2 is too great to be considered
            !  a minor species, adjust the column amount of CO2 by an empirical factor
            !  to obtain the proper contribution.
            chi_co2 = colco2(jl,lay)/(coldry(jl,lay))
            ratco2 = 1.e20_JPRB*chi_co2/chi_mls(2,jp(jl,lay)+1)
            if (ratco2 .gt. 3.0_JPRB) then
              adjfac = 2.0_JPRB+(ratco2-2.0_JPRB)**0.79_JPRB
              adjcolco2 = adjfac*chi_mls(2,jp(jl,lay)+1)*coldry(jl,lay)*1.e-20_JPRB
            else
              adjcolco2 = colco2(jl,lay)
            endif

            ind0 = ((jp(jl,lay)-13)*5+(jt(jl,lay)-1))*nspb(7) + 1
            ind1 = ((jp(jl,lay)-12)*5+(jt1(jl,lay)-1))*nspb(7) + 1
            indm = indminor(jl,lay)
!$NEC unroll(NG7)
            !$ACC LOOP SEQ PRIVATE (absco2)
            do ig = 1, ng7
              absco2 = kb_mco2(indm,ig) + minorfrac(jl,lay) * &
                  (kb_mco2(indm+1,ig) - kb_mco2(indm,ig))
              taug(jl,ngs6+ig,lay) = colo3(jl,lay) * &
                  (fac00(jl,lay) * absb(ind0,ig) + &
                  fac10(jl,lay) * absb(ind0+1,ig) + &
                  fac01(jl,lay) * absb(ind1,ig) + &
                  fac11(jl,lay) * absb(ind1+1,ig)) &
                  + adjcolco2 * absco2
              fracs(jl,ngs6+ig,lay) = fracrefb(ig)
            enddo
#if defined(_OPENACC) || defined(OMPGPU)
#else
          enddo

          ! Empirical modification to code to improve stratospheric cooling rates
          ! for o3.  Revised to apply weighting for g-point reduction in this band.

!$NEC ivdep
          do ixp = 1, ixc0
            jl = ixhigh(ixp,lay)
#endif
            taug(jl,ngs6+6,lay)=taug(jl,ngs6+6,lay)*0.92_JPRB
            taug(jl,ngs6+7,lay)=taug(jl,ngs6+7,lay)*0.88_JPRB
            taug(jl,ngs6+8,lay)=taug(jl,ngs6+8,lay)*1.07_JPRB
            taug(jl,ngs6+9,lay)=taug(jl,ngs6+9,lay)*1.1_JPRB
            taug(jl,ngs6+10,lay)=taug(jl,ngs6+10,lay)*0.99_JPRB
            taug(jl,ngs6+11,lay)=taug(jl,ngs6+11,lay)*0.855_JPRB
#if defined(_OPENACC) || defined(OMPGPU)
           endif
#endif
          enddo

        enddo
        !$ACC END PARALLEL
        !$OMP END TARGET TEAMS DISTRIBUTE PARALLEL DO

      ENDIF

      !$ACC END DATA
#ifndef __NVCOMPILER
      !$OMP END TARGET DATA
#endif

END SUBROUTINE RRTM_TAUMOL7
