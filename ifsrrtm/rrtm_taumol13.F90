!----------------------------------------------------------------------------
SUBROUTINE RRTM_TAUMOL13 (KIDIA,KFDIA,KLEV,taug,&
 & P_TAUAERL,fac00,fac01,fac10,fac11,forfac,forfrac,indfor,jp,jt,jt1,oneminus,&
 & colh2o,coln2o,colco2,colo3,coldry,laytrop,selffac,selffrac,indself,fracs, &
 & rat_h2on2o, rat_h2on2o_1,minorfrac,indminor,laytrop_min,laytrop_max)

!     BAND 13:  2080-2250 cm-1 (low - H2O,N2O; high - nothing)

!     AUTHOR.
!     -------
!      JJMorcrette, ECMWF

  !     MODIFICATIONS.
!     --------------
!      M.Hamrud      01-Oct-2003 CY28 Cleaning
!      NEC           25-Oct-2007 Optimisations
!      JJMorcrette 20110613 flexible number of g-points
!      ABozzo 201306 updated to rrtmg v4.85
!     band 13:  2080-2250 cm-1 (low key - h2o,n2o; high minor - o3 minor)
! ---------------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM     ,JPRB

USE PARRRTM  , ONLY : JPBAND
USE YOERRTM  , ONLY : JPGPT  ,NG13  ,NGS12
USE YOERRTWN , ONLY : NSPA
USE YOERRTA13, ONLY : ABSA   ,FRACREFA,FRACREFB,SELFREF,FORREF,KA_MCO2, KA_MCO, KB_MO3
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
REAL(KIND=JPRB)   ,INTENT(IN)    :: coln2o(KIDIA:KFDIA,KLEV)
REAL(KIND=JPRB)   ,INTENT(IN)    :: colco2(KIDIA:KFDIA,KLEV)
REAL(KIND=JPRB)   ,INTENT(IN)    :: colo3(KIDIA:KFDIA,KLEV)
REAL(KIND=JPRB)   ,INTENT(IN)    :: coldry(KIDIA:KFDIA,KLEV)
INTEGER(KIND=JPIM),INTENT(IN)    :: laytrop(KIDIA:KFDIA)
REAL(KIND=JPRB)   ,INTENT(IN)    :: selffac(KIDIA:KFDIA,KLEV)
REAL(KIND=JPRB)   ,INTENT(IN)    :: selffrac(KIDIA:KFDIA,KLEV)
INTEGER(KIND=JPIM),INTENT(IN)    :: indself(KIDIA:KFDIA,KLEV)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: fracs(KIDIA:KFDIA,JPGPT,KLEV)

REAL(KIND=JPRB)   ,INTENT(IN)   :: rat_h2on2o(KIDIA:KFDIA,KLEV)
REAL(KIND=JPRB)   ,INTENT(IN)   :: rat_h2on2o_1(KIDIA:KFDIA,KLEV)
INTEGER(KIND=JPIM),INTENT(IN)   :: indfor(KIDIA:KFDIA,KLEV)
REAL(KIND=JPRB)   ,INTENT(IN)   :: forfac(KIDIA:KFDIA,KLEV)
REAL(KIND=JPRB)   ,INTENT(IN)   :: forfrac(KIDIA:KFDIA,KLEV)
REAL(KIND=JPRB)   ,INTENT(IN)   :: minorfrac(KIDIA:KFDIA,KLEV)
INTEGER(KIND=JPIM),INTENT(IN)   :: indminor(KIDIA:KFDIA,KLEV)
INTEGER(KIND=JPIM), OPTIONAL, INTENT(INOUT) :: laytrop_min, laytrop_max

! ---------------------------------------------------------------------------


REAL(KIND=JPRB) :: speccomb,speccomb1, speccomb_planck, &
                   & speccomb_mco2, speccomb_mco
INTEGER(KIND=JPIM) :: ind0,ind1,inds,indf,indm
INTEGER(KIND=JPIM) :: IG, JS, lay, JS1, JPL, JMCO2, JMCO

REAL(KIND=JPRB) :: refrat_planck_a, refrat_m_a, refrat_m_a3
REAL(KIND=JPRB) ::  fac000, fac100, fac200,&
 & fac010, fac110, fac210, &
 & fac001, fac101, fac201, &
 & fac011, fac111, fac211
REAL(KIND=JPRB) :: p, p4, fk0, fk1, fk2

REAL(KIND=JPRB) :: taufor,tauself,tau_major(ng13),tau_major1(ng13), co2m1, co2m2, absco2
REAL(KIND=JPRB) :: com1, com2, absco, abso3
REAL(KIND=JPRB) :: chi_co2, ratco2, adjfac, adjcolco2

REAL(KIND=JPRB) :: fs, specmult, specparm,  &
& fs1, specmult1, specparm1, &
& fmco2, specmult_mco2, specparm_mco2, &
& fmco , specmult_mco , specparm_mco , &
& fpl, specmult_planck, specparm_planck

REAL(KIND=JPRB)   :: colco(KIDIA:KFDIA,KLEV) !left =0 for now,not passed from rrtm_gasbas1a

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
    !$ACC             colh2o, coln2o, colco2, colo3, coldry, laytrop, selffac, &
    !$ACC             selffrac, indself, fracs, rat_h2on2o, rat_h2on2o_1, &
    !$ACC             indfor, forfac, forfrac, minorfrac, indminor) &
    !$ACC       CREATE(colco)
#ifndef __NVCOMPILER
    !$OMP TARGET DATA MAP(PRESENT, ALLOC: taug, P_TAUAERL, fac00, fac01, fac10, fac11, jp, jt, jt1, &
    !$OMP             colh2o, coln2o, colco2, colo3, coldry, laytrop, selffac, &
    !$OMP             selffrac, indself, fracs, rat_h2on2o, rat_h2on2o_1, &
    !$OMP             indfor, forfac, forfrac, minorfrac, indminor)
#endif
    !$OMP TARGET ENTER DATA MAP(ALLOC: colco)

    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO COLLAPSE(2)
    !$ACC PARALLEL DEFAULT(NONE) ASYNC(1)
    !$ACC LOOP GANG VECTOR COLLAPSE(2)
    do lay = 1,KLEV
      do jc = KIDIA, KFDIA
        colco(jc,lay) = 0._JPRB
      end do
    end do
    !$ACC END PARALLEL
    !$OMP END TARGET TEAMS DISTRIBUTE PARALLEL DO

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

      ! P = 473.420 mb (Level 5)
      refrat_planck_a = chi_mls(1,5)/chi_mls(4,5)

      ! P = 1053. (Level 1)
      refrat_m_a = chi_mls(1,1)/chi_mls(4,1)

      ! P = 706. (Level 3)
      refrat_m_a3 = chi_mls(1,3)/chi_mls(4,3)

      ! Compute the optical depth by interpolating in ln(pressure),
      ! temperature, and appropriate species.  Below laytrop, the water
      ! vapor self-continuum and foreign continuum is interpolated
      ! (in temperature) separately.

      ! Lower atmosphere loop
      !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO COLLAPSE(2) PRIVATE(speccomb, speccomb1, speccomb_planck, speccomb_mco2, speccomb_mco, &
      !$OMP   ind0, ind1, &
      !$OMP   inds, indf, indm, js, js1, jpl, jmco2, jmco, fac000, fac100, fac200, fac010, &
      !$OMP   fac110, fac210, fac001, fac101, fac201, fac011, fac111, fac211, p, p4, fk0, fk1, fk2, chi_co2, ratco2, &
      !$OMP   adjfac, adjcolco2, fs, specmult, specparm, fs1, specmult1, specparm1, fmco2, specmult_mco2, &
      !$OMP   specparm_mco2, fmco, specmult_mco, specparm_mco, fpl, specmult_planck, specparm_planck, tau_major, &
      !$OMP   tau_major1, taufor, tauself, co2m1, co2m2, absco2, &
      !$OMP   com1, com2, absco) !THREAD_LIMIT(128)
      !$ACC WAIT
      !$ACC PARALLEL DEFAULT(NONE) ASYNC(1)
      !$ACC LOOP GANG VECTOR COLLAPSE(2) PRIVATE(speccomb, speccomb1, speccomb_planck, speccomb_mco2, speccomb_mco, &
      !$ACC   ind0, ind1, &
      !$ACC   inds, indf, indm, js, js1, jpl, jmco2, jmco, fac000, fac100, fac200, fac010, &
      !$ACC   fac110, fac210, fac001, fac101, fac201, fac011, fac111, fac211, p, p4, fk0, fk1, fk2, chi_co2, ratco2, &
      !$ACC   adjfac, adjcolco2, fs, specmult, specparm, fs1, specmult1, specparm1, fmco2, specmult_mco2, &
      !$ACC   specparm_mco2, fmco, specmult_mco, specparm_mco, fpl, specmult_planck, specparm_planck, tau_major, &
      !$ACC   tau_major1)
      do lay = 1, llaytrop_min
        do jl = KIDIA, KFDIA

          speccomb = colh2o(jl,lay) + rat_h2on2o(jl,lay)*coln2o(jl,lay)
          specparm = MIN(colh2o(jl,lay)/speccomb,oneminus)
          specmult = 8._JPRB*(specparm)
          js = 1 + int(specmult)
          fs = MOD1(specmult)

          speccomb1 = colh2o(jl,lay) + rat_h2on2o_1(jl,lay)*coln2o(jl,lay)
          specparm1 = MIN(colh2o(jl,lay)/speccomb1,oneminus)
          specmult1 = 8._JPRB*(specparm1)
          js1 = 1 + int(specmult1)
          fs1 = MOD1(specmult1)

          speccomb_mco2 = colh2o(jl,lay) + refrat_m_a*coln2o(jl,lay)
          specparm_mco2 = MIN(colh2o(jl,lay)/speccomb_mco2,oneminus)
          specmult_mco2 = 8._JPRB*specparm_mco2
          jmco2 = 1 + int(specmult_mco2)
          fmco2 = MOD1(specmult_mco2)

          !  In atmospheres where the amount of CO2 is too great to be considered
          !  a minor species, adjust the column amount of CO2 by an empirical factor
          !  to obtain the proper contribution.
          chi_co2 = colco2(jl,lay)/(coldry(jl,lay))
          ratco2 = 1.e20_JPRB*chi_co2/3.55e-4_JPRB
          if (ratco2 .gt. 3.0_JPRB) then
            adjfac = 2.0_JPRB+(ratco2-2.0_JPRB)**0.68_JPRB
            adjcolco2 = adjfac*3.55e-4_JPRB*coldry(jl,lay)*1.e-20_JPRB
          else
            adjcolco2 = colco2(jl,lay)
          endif

          speccomb_mco = colh2o(jl,lay) + refrat_m_a3*coln2o(jl,lay)
          specparm_mco = MIN(colh2o(jl,lay)/speccomb_mco,oneminus)
          specmult_mco = 8._JPRB*specparm_mco
          jmco = 1 + int(specmult_mco)
          fmco = MOD1(specmult_mco)

          speccomb_planck = colh2o(jl,lay)+refrat_planck_a*coln2o(jl,lay)
          specparm_planck = MIN(colh2o(jl,lay)/speccomb_planck,oneminus)
          specmult_planck = 8._JPRB*specparm_planck
          jpl = 1 + int(specmult_planck)
          fpl = MOD1(specmult_planck)

          ind0 = ((jp(jl,lay)-1)*5+(jt(jl,lay)-1))*nspa(13) + js
          ind1 = (jp(jl,lay)*5+(jt1(jl,lay)-1))*nspa(13) + js1
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
            do ig = 1, ng13
            tau_major(ig) = speccomb *    &
             (fac000 * absa(ind0,ig)    + &
              fac100 * absa(ind0+1,ig)  + &
              fac200 * absa(ind0+2,ig)  + &
              fac010 * absa(ind0+9,ig)  + &
              fac110 * absa(ind0+10,ig) + &
              fac210 * absa(ind0+11,ig))
            end do
#else
!$NEC unroll(NG13)
            tau_major(1:ng13) = speccomb *    &
             (fac000 * absa(ind0,1:ng13)    + &
              fac100 * absa(ind0+1,1:ng13)  + &
              fac200 * absa(ind0+2,1:ng13)  + &
              fac010 * absa(ind0+9,1:ng13)  + &
              fac110 * absa(ind0+10,1:ng13) + &
              fac210 * absa(ind0+11,1:ng13))
#endif
          else if (specparm .gt. 0.875_JPRB) then
#if defined(__amdflang__) && defined(OMPGPU)
            do ig = 1, ng13
            tau_major(ig) = speccomb *   &
             (fac200 * absa(ind0-1,ig) + &
              fac100 * absa(ind0,ig)   + &
              fac000 * absa(ind0+1,ig) + &
              fac210 * absa(ind0+8,ig) + &
              fac110 * absa(ind0+9,ig) + &
              fac010 * absa(ind0+10,ig))
            end do
#else
!$NEC unroll(NG13)
            tau_major(1:ng13) = speccomb *   &
             (fac200 * absa(ind0-1,1:ng13) + &
              fac100 * absa(ind0,1:ng13)   + &
              fac000 * absa(ind0+1,1:ng13) + &
              fac210 * absa(ind0+8,1:ng13) + &
              fac110 * absa(ind0+9,1:ng13) + &
              fac010 * absa(ind0+10,1:ng13))
#endif
          else
#if defined(__amdflang__) && defined(OMPGPU)
            do ig = 1, ng13
            tau_major(ig) = speccomb *   &
             (fac000 * absa(ind0,ig)   + &
              fac100 * absa(ind0+1,ig) + &
              fac010 * absa(ind0+9,ig) + &
              fac110 * absa(ind0+10,ig))
            end do
#else
!$NEC unroll(NG13)
            tau_major(1:ng13) = speccomb *   &
             (fac000 * absa(ind0,1:ng13)   + &
              fac100 * absa(ind0+1,1:ng13) + &
              fac010 * absa(ind0+9,1:ng13) + &
              fac110 * absa(ind0+10,1:ng13))
#endif
          endif

          if (specparm1 .lt. 0.125_JPRB) then
#if defined(__amdflang__) && defined(OMPGPU)
            do ig = 1, ng13
            tau_major1(ig) = speccomb1 *  &
             (fac001 * absa(ind1,ig)    + &
              fac101 * absa(ind1+1,ig)  + &
              fac201 * absa(ind1+2,ig)  + &
              fac011 * absa(ind1+9,ig)  + &
              fac111 * absa(ind1+10,ig) + &
              fac211 * absa(ind1+11,ig))
            end do
#else
!$NEC unroll(NG13)
            tau_major1(1:ng13) = speccomb1 *  &
             (fac001 * absa(ind1,1:ng13)    + &
              fac101 * absa(ind1+1,1:ng13)  + &
              fac201 * absa(ind1+2,1:ng13)  + &
              fac011 * absa(ind1+9,1:ng13)  + &
              fac111 * absa(ind1+10,1:ng13) + &
              fac211 * absa(ind1+11,1:ng13))
#endif
          else if (specparm1 .gt. 0.875_JPRB) then
#if defined(__amdflang__) && defined(OMPGPU)
            do ig = 1, ng13
            tau_major1(ig) = speccomb1 * &
             (fac201 * absa(ind1-1,ig) + &
              fac101 * absa(ind1,ig)   + &
              fac001 * absa(ind1+1,ig) + &
              fac211 * absa(ind1+8,ig) + &
              fac111 * absa(ind1+9,ig) + &
              fac011 * absa(ind1+10,ig))
            end do
#else
!$NEC unroll(NG13)
            tau_major1(1:ng13) = speccomb1 * &
             (fac201 * absa(ind1-1,1:ng13) + &
              fac101 * absa(ind1,1:ng13)   + &
              fac001 * absa(ind1+1,1:ng13) + &
              fac211 * absa(ind1+8,1:ng13) + &
              fac111 * absa(ind1+9,1:ng13) + &
              fac011 * absa(ind1+10,1:ng13))
#endif
          else
#if defined(__amdflang__) && defined(OMPGPU)
            do ig = 1, ng13
            tau_major1(ig) = speccomb1 * &
             (fac001 * absa(ind1,ig)   + &
              fac101 * absa(ind1+1,ig) + &
              fac011 * absa(ind1+9,ig) + &
              fac111 * absa(ind1+10,ig))
            end do
#else
!$NEC unroll(NG13)
            tau_major1(1:ng13) = speccomb1 * &
             (fac001 * absa(ind1,1:ng13)   + &
              fac101 * absa(ind1+1,1:ng13) + &
              fac011 * absa(ind1+9,1:ng13) + &
              fac111 * absa(ind1+10,1:ng13))
#endif
          endif

          !$ACC LOOP SEQ PRIVATE(taufor, tauself, co2m1, co2m2, absco2, &
          !$ACC   com1, com2, absco)
!$NEC unroll(NG13)
          do ig = 1, ng13
            tauself = selffac(jl,lay)* (selfref(inds,ig) + selffrac(jl,lay) * &
                 (selfref(inds+1,ig) - selfref(inds,ig)))
            taufor = forfac(jl,lay) * (forref(indf,ig) + forfrac(jl,lay) * &
                 (forref(indf+1,ig) - forref(indf,ig)))
            co2m1 = ka_mco2(jmco2,indm,ig) + fmco2 * &
                 (ka_mco2(jmco2+1,indm,ig) - ka_mco2(jmco2,indm,ig))
            co2m2 = ka_mco2(jmco2,indm+1,ig) + fmco2 * &
                 (ka_mco2(jmco2+1,indm+1,ig) - ka_mco2(jmco2,indm+1,ig))
            absco2 = co2m1 + minorfrac(jl,lay) * (co2m2 - co2m1)
            com1 = ka_mco(jmco,indm,ig) + fmco * &
                 (ka_mco(jmco+1,indm,ig) - ka_mco(jmco,indm,ig))
            com2 = ka_mco(jmco,indm+1,ig) + fmco * &
                 (ka_mco(jmco+1,indm+1,ig) - ka_mco(jmco,indm+1,ig))
            absco = com1 + minorfrac(jl,lay) * (com2 - com1)

            taug(jl,ngs12+ig,lay) = tau_major(ig) + tau_major1(ig) &
                 + tauself + taufor &
                 + adjcolco2*absco2 &
                 + colco(jl,lay)*absco
            fracs(jl,ngs12+ig,lay) = fracrefa(ig,jpl) + fpl * &
                 (fracrefa(ig,jpl+1)-fracrefa(ig,jpl))
          enddo
        enddo

      enddo
      !$ACC END PARALLEL
      !$OMP END TARGET TEAMS DISTRIBUTE PARALLEL DO

      ! Upper atmosphere loop
      !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO COLLAPSE(2) PRIVATE(indm, abso3)
      !$ACC PARALLEL DEFAULT(NONE) ASYNC(1)
      !$ACC LOOP GANG VECTOR COLLAPSE(2) PRIVATE(indm)
      do lay = llaytrop_max+1, KLEV
        do jl = KIDIA, KFDIA

          indm = indminor(jl,lay)
          !$ACC LOOP SEQ PRIVATE(abso3)
!$NEC unroll(NG13)
          do ig = 1, ng13
            abso3 = kb_mo3(indm,ig) + minorfrac(jl,lay) * &
                 (kb_mo3(indm+1,ig) - kb_mo3(indm,ig))
            taug(jl,ngs12+ig,lay) = colo3(jl,lay)*abso3
            fracs(jl,ngs12+ig,lay) =  fracrefb(ig)
          enddo
        enddo

      enddo
      !$ACC END PARALLEL
      !$OMP END TARGET TEAMS DISTRIBUTE PARALLEL DO

      IF (llaytrop_max /= llaytrop_min) THEN
        ! Mixed loop
        ! Lower atmosphere part
        !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO COLLAPSE(2) PRIVATE(speccomb, specparm, specmult, js, fs, speccomb1, specparm1, &
        !$OMP   specmult1, js1, fs1, speccomb_mco2, specparm_mco2, specmult_mco2, jmco2, fmco2, chi_co2, ratco2, adjfac, &
        !$OMP   adjcolco2, speccomb_mco, specparm_mco, specmult_mco, jmco, fmco, speccomb_planck, specparm_planck, &
        !$OMP   specmult_planck, jpl, fpl, ind0, ind1, inds, indf, indm, p, p4, fk0, fk1, fk2, fac000, fac100, fac200, &
        !$OMP   fac010, fac110, fac210, fac001, fac101, fac201, fac011, fac111, fac211, tau_major, tau_major1, &
        !$OMP tauself, taufor, co2m1, absco2, com1, com2, absco, abso3)
        !$ACC PARALLEL DEFAULT(NONE) ASYNC(1)
        !$ACC LOOP GANG VECTOR COLLAPSE(2) PRIVATE(speccomb, specparm, specmult, js, fs, speccomb1, specparm1, &
        !$ACC   specmult1, js1, fs1, speccomb_mco2, specparm_mco2, specmult_mco2, jmco2, fmco2, chi_co2, ratco2, adjfac, &
        !$ACC   adjcolco2, speccomb_mco, specparm_mco, specmult_mco, jmco, fmco, speccomb_planck, specparm_planck, &
        !$ACC   specmult_planck, jpl, fpl, ind0, ind1, inds, indf, indm, p, p4, fk0, fk1, fk2, fac000, fac100, fac200, &
        !$ACC   fac010, fac110, fac210, fac001, fac101, fac201, fac011, fac111, fac211, tau_major, tau_major1)
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

            speccomb = colh2o(jl,lay) + rat_h2on2o(jl,lay)*coln2o(jl,lay)
            specparm = MIN(colh2o(jl,lay)/speccomb,oneminus)
            specmult = 8._JPRB*(specparm)
            js = 1 + int(specmult)
            fs = MOD1(specmult)

            speccomb1 = colh2o(jl,lay) + rat_h2on2o_1(jl,lay)*coln2o(jl,lay)
            specparm1 = MIN(colh2o(jl,lay)/speccomb1,oneminus)
            specmult1 = 8._JPRB*(specparm1)
            js1 = 1 + int(specmult1)
            fs1 = MOD1(specmult1)

            speccomb_mco2 = colh2o(jl,lay) + refrat_m_a*coln2o(jl,lay)
            specparm_mco2 = MIN(colh2o(jl,lay)/speccomb_mco2,oneminus)
            specmult_mco2 = 8._JPRB*specparm_mco2
            jmco2 = 1 + int(specmult_mco2)
            fmco2 = MOD1(specmult_mco2)

            !  In atmospheres where the amount of CO2 is too great to be considered
            !  a minor species, adjust the column amount of CO2 by an empirical factor
            !  to obtain the proper contribution.
            chi_co2 = colco2(jl,lay)/(coldry(jl,lay))
            ratco2 = 1.e20_JPRB*chi_co2/3.55e-4_JPRB
            if (ratco2 .gt. 3.0_JPRB) then
              adjfac = 2.0_JPRB+(ratco2-2.0_JPRB)**0.68_JPRB
              adjcolco2 = adjfac*3.55e-4_JPRB*coldry(jl,lay)*1.e-20_JPRB
            else
              adjcolco2 = colco2(jl,lay)
            endif

            speccomb_mco = colh2o(jl,lay) + refrat_m_a3*coln2o(jl,lay)
            specparm_mco = MIN(colh2o(jl,lay)/speccomb_mco,oneminus)
            specmult_mco = 8._JPRB*specparm_mco
            jmco = 1 + int(specmult_mco)
            fmco = MOD1(specmult_mco)

            speccomb_planck = colh2o(jl,lay)+refrat_planck_a*coln2o(jl,lay)
            specparm_planck = MIN(colh2o(jl,lay)/speccomb_planck,oneminus)
            specmult_planck = 8._JPRB*specparm_planck
            jpl = 1 + int(specmult_planck)
            fpl = MOD1(specmult_planck)

            ind0 = ((jp(jl,lay)-1)*5+(jt(jl,lay)-1))*nspa(13) + js
            ind1 = (jp(jl,lay)*5+(jt1(jl,lay)-1))*nspa(13) + js1
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
            do ig = 1, ng13
              tau_major(ig) = speccomb *    &
              (fac000 * absa(ind0,ig)    + &
                fac100 * absa(ind0+1,ig)  + &
                fac200 * absa(ind0+2,ig)  + &
                fac010 * absa(ind0+9,ig)  + &
                fac110 * absa(ind0+10,ig) + &
                fac210 * absa(ind0+11,ig))
            end do
#else
!$NEC unroll(NG13)
              tau_major(1:ng13) = speccomb *    &
              (fac000 * absa(ind0,1:ng13)    + &
                fac100 * absa(ind0+1,1:ng13)  + &
                fac200 * absa(ind0+2,1:ng13)  + &
                fac010 * absa(ind0+9,1:ng13)  + &
                fac110 * absa(ind0+10,1:ng13) + &
                fac210 * absa(ind0+11,1:ng13))
#endif
            else if (specparm .gt. 0.875_JPRB) then
#if defined(__amdflang__) && defined(OMPGPU)
            do ig = 1, ng13
              tau_major(ig) = speccomb *   &
              (fac200 * absa(ind0-1,ig) + &
                fac100 * absa(ind0,ig)   + &
                fac000 * absa(ind0+1,ig) + &
                fac210 * absa(ind0+8,ig) + &
                fac110 * absa(ind0+9,ig) + &
                fac010 * absa(ind0+10,ig))
            end do
#else
!$NEC unroll(NG13)
              tau_major(1:ng13) = speccomb *   &
              (fac200 * absa(ind0-1,1:ng13) + &
                fac100 * absa(ind0,1:ng13)   + &
                fac000 * absa(ind0+1,1:ng13) + &
                fac210 * absa(ind0+8,1:ng13) + &
                fac110 * absa(ind0+9,1:ng13) + &
                fac010 * absa(ind0+10,1:ng13))
#endif
            else
#if defined(__amdflang__) && defined(OMPGPU)
            do ig = 1, ng13
              tau_major(ig) = speccomb *   &
              (fac000 * absa(ind0,ig)   + &
                fac100 * absa(ind0+1,ig) + &
                fac010 * absa(ind0+9,ig) + &
                fac110 * absa(ind0+10,ig))
            end do
#else
!$NEC unroll(NG13)
              tau_major(1:ng13) = speccomb *   &
              (fac000 * absa(ind0,1:ng13)   + &
                fac100 * absa(ind0+1,1:ng13) + &
                fac010 * absa(ind0+9,1:ng13) + &
                fac110 * absa(ind0+10,1:ng13))
#endif
            endif

            if (specparm1 .lt. 0.125_JPRB) then
#if defined(__amdflang__) && defined(OMPGPU)
            do ig = 1, ng13
              tau_major1(ig) = speccomb1 *  &
              (fac001 * absa(ind1,ig)    + &
                fac101 * absa(ind1+1,ig)  + &
                fac201 * absa(ind1+2,ig)  + &
                fac011 * absa(ind1+9,ig)  + &
                fac111 * absa(ind1+10,ig) + &
                fac211 * absa(ind1+11,ig))
            end do
#else
!$NEC unroll(NG13)
              tau_major1(1:ng13) = speccomb1 *  &
              (fac001 * absa(ind1,1:ng13)    + &
                fac101 * absa(ind1+1,1:ng13)  + &
                fac201 * absa(ind1+2,1:ng13)  + &
                fac011 * absa(ind1+9,1:ng13)  + &
                fac111 * absa(ind1+10,1:ng13) + &
                fac211 * absa(ind1+11,1:ng13))
#endif
            else if (specparm1 .gt. 0.875_JPRB) then
#if defined(__amdflang__) && defined(OMPGPU)
            do ig = 1, ng13
              tau_major1(ig) = speccomb1 * &
              (fac201 * absa(ind1-1,ig) + &
                fac101 * absa(ind1,ig)   + &
                fac001 * absa(ind1+1,ig) + &
                fac211 * absa(ind1+8,ig) + &
                fac111 * absa(ind1+9,ig) + &
                fac011 * absa(ind1+10,ig))
            end do
#else
!$NEC unroll(NG13)
              tau_major1(1:ng13) = speccomb1 * &
              (fac201 * absa(ind1-1,1:ng13) + &
                fac101 * absa(ind1,1:ng13)   + &
                fac001 * absa(ind1+1,1:ng13) + &
                fac211 * absa(ind1+8,1:ng13) + &
                fac111 * absa(ind1+9,1:ng13) + &
                fac011 * absa(ind1+10,1:ng13))
#endif
            else
#if defined(__amdflang__) && defined(OMPGPU)
            do ig = 1, ng13
              tau_major1(ig) = speccomb1 * &
              (fac001 * absa(ind1,ig)   + &
                fac101 * absa(ind1+1,ig) + &
                fac011 * absa(ind1+9,ig) + &
                fac111 * absa(ind1+10,ig))
            end do
#else
!$NEC unroll(NG13)
              tau_major1(1:ng13) = speccomb1 * &
              (fac001 * absa(ind1,1:ng13)   + &
                fac101 * absa(ind1+1,1:ng13) + &
                fac011 * absa(ind1+9,1:ng13) + &
                fac111 * absa(ind1+10,1:ng13))
#endif
            endif

!$NEC unroll(NG13)
            !$ACC LOOP SEQ PRIVATE(tauself, taufor, co2m1, absco2, com1, com2, absco)
            do ig = 1, ng13
              tauself = selffac(jl,lay)* (selfref(inds,ig) + selffrac(jl,lay) * &
                  (selfref(inds+1,ig) - selfref(inds,ig)))
              taufor = forfac(jl,lay) * (forref(indf,ig) + forfrac(jl,lay) * &
                  (forref(indf+1,ig) - forref(indf,ig)))
              co2m1 = ka_mco2(jmco2,indm,ig) + fmco2 * &
                  (ka_mco2(jmco2+1,indm,ig) - ka_mco2(jmco2,indm,ig))
              co2m2 = ka_mco2(jmco2,indm+1,ig) + fmco2 * &
                  (ka_mco2(jmco2+1,indm+1,ig) - ka_mco2(jmco2,indm+1,ig))
              absco2 = co2m1 + minorfrac(jl,lay) * (co2m2 - co2m1)
              com1 = ka_mco(jmco,indm,ig) + fmco * &
                  (ka_mco(jmco+1,indm,ig) - ka_mco(jmco,indm,ig))
              com2 = ka_mco(jmco,indm+1,ig) + fmco * &
                  (ka_mco(jmco+1,indm+1,ig) - ka_mco(jmco,indm+1,ig))
              absco = com1 + minorfrac(jl,lay) * (com2 - com1)

              taug(jl,ngs12+ig,lay) = tau_major(ig) + tau_major1(ig) &
                  + tauself + taufor &
                  + adjcolco2*absco2 &
                  + colco(jl,lay)*absco
              fracs(jl,ngs12+ig,lay) = fracrefa(ig,jpl) + fpl * &
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

            indm = indminor(jl,lay)
!$NEC unroll(NG13)
            !$ACC LOOP SEQ PRIVATE(abso3)
            do ig = 1, ng13
              abso3 = kb_mo3(indm,ig) + minorfrac(jl,lay) * &
                  (kb_mo3(indm+1,ig) - kb_mo3(indm,ig))
              taug(jl,ngs12+ig,lay) = colo3(jl,lay)*abso3
              fracs(jl,ngs12+ig,lay) =  fracrefb(ig)
            enddo
#if defined(_OPENACC) || defined(OMPGPU)
           endif
#endif
          enddo

        enddo
        !$ACC END PARALLEL
        !$OMP END TARGET TEAMS DISTRIBUTE PARALLEL DO

      ENDIF

      !$ACC WAIT
      !$ACC END DATA
      !$OMP TARGET EXIT DATA MAP(DELETE:colco)
#ifndef __NVCOMPILER
      !$OMP END TARGET DATA
#endif

END SUBROUTINE RRTM_TAUMOL13
