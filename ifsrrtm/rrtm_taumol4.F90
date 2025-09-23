!----------------------------------------------------------------------------
SUBROUTINE RRTM_TAUMOL4 (KIDIA,KFDIA,KLEV,taug,&
 & P_TAUAERL,fac00,fac01,fac10,fac11,forfac,forfrac,indfor,jp,jt,jt1,oneminus,&
 & colh2o,colco2,colo3,laytrop,selffac,selffrac,indself,fracs, &
 & rat_h2oco2, rat_h2oco2_1, rat_o3co2, rat_o3co2_1,laytrop_min,laytrop_max)

!     BAND 4:  630-700 cm-1 (low - H2O,CO2; high - O3,CO2)

!     AUTHOR.
!     -------
!      JJMorcrette, ECMWF

!     MODIFICATIONS.
!     --------------
!      M.Hamrud      01-Oct-2003 CY28 Cleaning
!      NEC           25-Oct-2007 Optimisations
!      JJMorcrette 20110613 flexible number of g-points
!      ABozzo 201306 updated to rrtmg v4.85
! ---------------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM     ,JPRB

USE PARRRTM  , ONLY : JPBAND
USE YOERRTM  , ONLY : JPGPT  ,NG4   ,NGS3
USE YOERRTWN , ONLY : NSPA   ,NSPB
USE YOERRTA4 , ONLY : ABSA   ,ABSB   ,FRACREFA, FRACREFB,SELFREF,FORREF
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
REAL(KIND=JPRB)   ,INTENT(IN)    :: colo3(KIDIA:KFDIA,KLEV)
INTEGER(KIND=JPIM),INTENT(IN)    :: laytrop(KIDIA:KFDIA)
REAL(KIND=JPRB)   ,INTENT(IN)    :: selffac(KIDIA:KFDIA,KLEV)
REAL(KIND=JPRB)   ,INTENT(IN)    :: selffrac(KIDIA:KFDIA,KLEV)
INTEGER(KIND=JPIM),INTENT(IN)    :: indself(KIDIA:KFDIA,KLEV)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: fracs(KIDIA:KFDIA,JPGPT,KLEV)

REAL(KIND=JPRB)   ,INTENT(IN)   :: rat_h2oco2(KIDIA:KFDIA,KLEV)
REAL(KIND=JPRB)   ,INTENT(IN)   :: rat_h2oco2_1(KIDIA:KFDIA,KLEV)
REAL(KIND=JPRB)   ,INTENT(IN)   :: rat_o3co2(KIDIA:KFDIA,KLEV)
REAL(KIND=JPRB)   ,INTENT(IN)   :: rat_o3co2_1(KIDIA:KFDIA,KLEV)
INTEGER(KIND=JPIM),INTENT(IN)   :: indfor(KIDIA:KFDIA,KLEV)
REAL(KIND=JPRB)   ,INTENT(IN)   :: forfac(KIDIA:KFDIA,KLEV)
REAL(KIND=JPRB)   ,INTENT(IN)   :: forfrac(KIDIA:KFDIA,KLEV)
INTEGER(KIND=JPIM), OPTIONAL, INTENT(INOUT) :: laytrop_min, laytrop_max
! ---------------------------------------------------------------------------

REAL(KIND=JPRB) :: speccomb,speccomb1, speccomb_planck
INTEGER(KIND=JPIM) :: ind0,ind1,inds,indf

INTEGER(KIND=JPIM) :: IG, JS, lay, JS1, JPL


REAL(KIND=JPRB) :: refrat_planck_a, refrat_planck_b
 REAL(KIND=JPRB) ::  fac000, fac100, fac200,&
 & fac010, fac110, fac210, &
 & fac001, fac101, fac201, &
 & fac011, fac111, fac211
REAL(KIND=JPRB) :: p, p4, fk0, fk1, fk2
REAL(KIND=JPRB) :: taufor,tauself,tau_major(ng4),tau_major1(ng4)
REAL(KIND=JPRB) :: fs, specmult, specparm,  &
 & fs1, specmult1, specparm1, &
 & fpl, specmult_PLANCK, specparm_PLANCK

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
!$ACC             colh2o, colco2, colo3, laytrop, selffac, selffrac, indself, &
!$ACC             fracs, rat_h2oco2, rat_h2oco2_1, rat_o3co2, rat_o3co2_1, &
!$ACC             indfor, forfac, forfrac)
!$OMP TARGET DATA MAP(PRESENT, ALLOC: taug, P_TAUAERL, fac00, fac01, fac10, fac11, jp, jt, jt1, &
!$OMP             colh2o, colco2, colo3, laytrop, selffac, selffrac, indself, &
!$OMP             fracs, rat_h2oco2, rat_h2oco2_1, rat_o3co2, rat_o3co2_1, &
!$OMP             indfor, forfac, forfrac)

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


!     Compute the optical depth by interpolating in ln(pressure),
!     temperature, and appropriate species.  Below LAYTROP, the water
!     vapor self-continuum is interpolated (in temperature) separately.

      ! P =   142.5940 mb
      refrat_planck_a = chi_mls(1,11)/chi_mls(2,11)

      ! P = 95.58350 mb
      refrat_planck_b = chi_mls(3,13)/chi_mls(2,13)

      ! Compute the optical depth by interpolating in ln(pressure) and
      ! temperature, and appropriate species.  Below laytrop, the water
      ! vapor self-continuum and foreign continuum is interpolated (in temperature)
      ! separately.

      ! Lower atmosphere loop
      !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO COLLAPSE(2) PRIVATE(speccomb,speccomb1, speccomb_planck, ind0, ind1, inds, indf, js, js1, &
      !$OMP   jpl, fac000, fac100, fac200, fac010, fac110, fac210, fac001, fac101, fac201, fac011, fac111, fac211, p, &
      !$OMP   p4, fk0, fk1, fk2, fs, specmult, specparm, fs1, specmult1, specparm1, fpl, specmult_PLANCK, &
      !$OMP   specparm_PLANCK, tau_major, tau_major1, taufor, tauself) THREAD_LIMIT(128)
      !$ACC WAIT
      !$ACC PARALLEL DEFAULT(NONE) ASYNC(1)
      !$ACC LOOP GANG VECTOR COLLAPSE(2) PRIVATE(speccomb,speccomb1, speccomb_planck, ind0, ind1, inds, indf, js, js1, &
      !$ACC   jpl, fac000, fac100, fac200, fac010, fac110, fac210, fac001, fac101, fac201, fac011, fac111, fac211, p, &
      !$ACC   p4, fk0, fk1, fk2, fs, specmult, specparm, fs1, specmult1, specparm1, fpl, specmult_PLANCK, &
      !$ACC   specparm_PLANCK, tau_major, tau_major1)
      do lay = 1, llaytrop_min
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

          speccomb_planck = colh2o(jl,lay)+refrat_planck_a*colco2(jl,lay)
          specparm_planck = MIN(colh2o(jl,lay)/speccomb_planck,oneminus)
          specmult_planck = 8._JPRB*specparm_planck
          jpl= 1 + int(specmult_planck)
          fpl = MOD1(specmult_planck)

          ind0 = ((jp(jl,lay)-1)*5+(jt(jl,lay)-1))*nspa(4) + js
          ind1 = (jp(jl,lay)*5+(jt1(jl,lay)-1))*nspa(4) + js1
          inds = indself(jl,lay)
          indf = indfor(jl,lay)

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
            do ig = 1, ng4
            tau_major(ig) = speccomb *    &
             (fac000 * absa(ind0,ig)    + &
              fac100 * absa(ind0+1,ig)  + &
              fac200 * absa(ind0+2,ig)  + &
              fac010 * absa(ind0+9,ig)  + &
              fac110 * absa(ind0+10,ig) + &
              fac210 * absa(ind0+11,ig))
            enddo
#else
!$NEC unroll(NG4)
            tau_major(1:ng4) = speccomb *    &
             (fac000 * absa(ind0,1:ng4)    + &
              fac100 * absa(ind0+1,1:ng4)  + &
              fac200 * absa(ind0+2,1:ng4)  + &
              fac010 * absa(ind0+9,1:ng4)  + &
              fac110 * absa(ind0+10,1:ng4) + &
              fac210 * absa(ind0+11,1:ng4))
#endif
          else if (specparm .gt. 0.875_JPRB) then
#if defined(__amdflang__) && defined(OMPGPU)
            do ig = 1, ng4
            tau_major(ig) = speccomb *   &
             (fac200 * absa(ind0-1,ig) + &
              fac100 * absa(ind0,ig)   + &
              fac000 * absa(ind0+1,ig) + &
              fac210 * absa(ind0+8,ig) + &
              fac110 * absa(ind0+9,ig) + &
              fac010 * absa(ind0+10,ig))
            enddo
#else
!$NEC unroll(NG4)
            tau_major(1:ng4) = speccomb *   &
             (fac200 * absa(ind0-1,1:ng4) + &
              fac100 * absa(ind0,1:ng4)   + &
              fac000 * absa(ind0+1,1:ng4) + &
              fac210 * absa(ind0+8,1:ng4) + &
              fac110 * absa(ind0+9,1:ng4) + &
              fac010 * absa(ind0+10,1:ng4))
#endif
          else
#if defined(__amdflang__) && defined(OMPGPU)
            do ig = 1, ng4
             tau_major(ig) = speccomb *   &
              (fac000 * absa(ind0,ig)   + &
               fac100 * absa(ind0+1,ig) + &
               fac010 * absa(ind0+9,ig) + &
               fac110 * absa(ind0+10,ig))
            enddo
#else
!$NEC unroll(NG4)
             tau_major(1:ng4) = speccomb *   &
              (fac000 * absa(ind0,1:ng4)   + &
               fac100 * absa(ind0+1,1:ng4) + &
               fac010 * absa(ind0+9,1:ng4) + &
               fac110 * absa(ind0+10,1:ng4))
#endif
          endif

          if (specparm1 .lt. 0.125_JPRB) then
#if defined(__amdflang__) && defined(OMPGPU)
            do ig = 1, ng4
            tau_major1(ig) = speccomb1 *  &
             (fac001 * absa(ind1,ig)    + &
              fac101 * absa(ind1+1,ig)  + &
              fac201 * absa(ind1+2,ig)  + &
              fac011 * absa(ind1+9,ig)  + &
              fac111 * absa(ind1+10,ig) + &
              fac211 * absa(ind1+11,ig))
            enddo
#else
!$NEC unroll(NG4)
            tau_major1(1:ng4) = speccomb1 *  &
             (fac001 * absa(ind1,1:ng4)    + &
              fac101 * absa(ind1+1,1:ng4)  + &
              fac201 * absa(ind1+2,1:ng4)  + &
              fac011 * absa(ind1+9,1:ng4)  + &
              fac111 * absa(ind1+10,1:ng4) + &
              fac211 * absa(ind1+11,1:ng4))
#endif
          else if (specparm1 .gt. 0.875_JPRB) then
#if defined(__amdflang__) && defined(OMPGPU)
            do ig = 1, ng4
            tau_major1(ig) = speccomb1 * &
             (fac201 * absa(ind1-1,ig) + &
              fac101 * absa(ind1,ig)   + &
              fac001 * absa(ind1+1,ig) + &
              fac211 * absa(ind1+8,ig) + &
              fac111 * absa(ind1+9,ig) + &
              fac011 * absa(ind1+10,ig))
            enddo
#else
!$NEC unroll(NG4)
            tau_major1(1:ng4) = speccomb1 * &
             (fac201 * absa(ind1-1,1:ng4) + &
              fac101 * absa(ind1,1:ng4)   + &
              fac001 * absa(ind1+1,1:ng4) + &
              fac211 * absa(ind1+8,1:ng4) + &
              fac111 * absa(ind1+9,1:ng4) + &
              fac011 * absa(ind1+10,1:ng4))
#endif
          else
#if defined(__amdflang__) && defined(OMPGPU)
            do ig = 1, ng4
            tau_major1(ig) = speccomb1 * &
             (fac001 * absa(ind1,ig)   + &
              fac101 * absa(ind1+1,ig) + &
              fac011 * absa(ind1+9,ig) + &
              fac111 * absa(ind1+10,ig))
            enddo
#else
!$NEC unroll(NG4)
            tau_major1(1:ng4) = speccomb1 * &
             (fac001 * absa(ind1,1:ng4)   + &
              fac101 * absa(ind1+1,1:ng4) + &
              fac011 * absa(ind1+9,1:ng4) + &
              fac111 * absa(ind1+10,1:ng4))
#endif
          endif

          !$ACC LOOP SEQ PRIVATE(taufor, tauself)
!$NEC unroll(NG4)
          do ig = 1, ng4
            tauself = selffac(jl,lay)* (selfref(inds,ig) + selffrac(jl,lay) * &
                 (selfref(inds+1,ig) - selfref(inds,ig)))
            taufor =  forfac(jl,lay) * (forref(indf,ig) + forfrac(jl,lay) * &
                 (forref(indf+1,ig) - forref(indf,ig)))

            taug(jl,ngs3+ig,lay) = tau_major(ig) + tau_major1(ig) &
                 + tauself + taufor
            fracs(jl,ngs3+ig,lay) = fracrefa(ig,jpl) + fpl * &
                 (fracrefa(ig,jpl+1)-fracrefa(ig,jpl))
          enddo
        enddo

      enddo
      !$ACC END PARALLEL
      !$OMP END TARGET TEAMS DISTRIBUTE PARALLEL DO

      ! Upper atmosphere loop
      !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO COLLAPSE(2) PRIVATE(speccomb, speccomb1, speccomb_planck, ind0, ind1, js, js1, jpl, &
      !$OMP   fac000, fac100, fac010, fac110, fac001, fac101, fac011, fac111, p4, fs, specmult, specparm, fs1, &
      !$OMP   specmult1, specparm1, fpl, specmult_PLANCK, specparm_PLANCK)
      !$ACC PARALLEL DEFAULT(NONE) ASYNC(1)
      !$ACC LOOP GANG VECTOR COLLAPSE(2) PRIVATE(speccomb, speccomb1, speccomb_planck, ind0, ind1, js, js1, jpl, &
      !$ACC   fac000, fac100, fac010, fac110, fac001, fac101, fac011, fac111, p4, fs, specmult, specparm, fs1, &
      !$ACC   specmult1, specparm1, fpl, specmult_PLANCK, specparm_PLANCK)
      do lay = llaytrop_max+1, KLEV
        do jl = KIDIA, KFDIA

          speccomb = colo3(jl,lay) + rat_o3co2(jl,lay)*colco2(jl,lay)
          specparm = MIN(colo3(jl,lay)/speccomb,oneminus)
          specmult = 4._JPRB*(specparm)
          js = 1 + int(specmult)
          fs = MOD1(specmult)

          speccomb1 = colo3(jl,lay) + rat_o3co2_1(jl,lay)*colco2(jl,lay)
          specparm1 = MIN(colo3(jl,lay)/speccomb1,oneminus)
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

          speccomb_planck = colo3(jl,lay)+refrat_planck_b*colco2(jl,lay)
          specparm_planck = MIN(colo3(jl,lay)/speccomb_planck,oneminus)
          specmult_planck = 4._JPRB*specparm_planck
          jpl= 1 + int(specmult_planck)
          fpl = MOD1(specmult_planck)

          ind0 = ((jp(jl,lay)-13)*5+(jt(jl,lay)-1))*nspb(4) + js
          ind1 = ((jp(jl,lay)-12)*5+(jt1(jl,lay)-1))*nspb(4) + js1
          !$ACC LOOP SEQ
!$NEC unroll(NG4)
          do ig = 1, ng4
            taug(jl,ngs3+ig,lay) =  speccomb * &
                 (fac000 * absb(ind0,ig) + &
                 fac100 * absb(ind0+1,ig) + &
                 fac010 * absb(ind0+5,ig) + &
                 fac110 * absb(ind0+6,ig)) &
                 + speccomb1 * &
                 (fac001 * absb(ind1,ig) +  &
                 fac101 * absb(ind1+1,ig) + &
                 fac011 * absb(ind1+5,ig) + &
                 fac111 * absb(ind1+6,ig))
            fracs(jl,ngs3+ig,lay) = fracrefb(ig,jpl) + fpl * &
                 (fracrefb(ig,jpl+1)-fracrefb(ig,jpl))
          enddo
        enddo
      enddo
      !$ACC END PARALLEL
      !$ACC WAIT
      !$OMP END TARGET TEAMS DISTRIBUTE PARALLEL DO

      ! Empirical modification to code to improve stratospheric cooling rates
      ! for co2.  Revised to apply weighting for g-point reduction in this band.
      !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO COLLAPSE(2)
      !$ACC PARALLEL DEFAULT(NONE) ASYNC(1)
      !$ACC LOOP GANG VECTOR COLLAPSE(2)
      do lay = llaytrop_max+1, KLEV
        do jl = KIDIA, KFDIA
          taug(jl,ngs3+8,lay)=taug(jl,ngs3+8,lay)*0.92_JPRB
          taug(jl,ngs3+9,lay)=taug(jl,ngs3+9,lay)*0.88_JPRB
          taug(jl,ngs3+10,lay)=taug(jl,ngs3+10,lay)*1.07_JPRB
          taug(jl,ngs3+11,lay)=taug(jl,ngs3+11,lay)*1.1_JPRB
          taug(jl,ngs3+12,lay)=taug(jl,ngs3+12,lay)*0.99_JPRB
          taug(jl,ngs3+13,lay)=taug(jl,ngs3+13,lay)*0.88_JPRB
          taug(jl,ngs3+14,lay)=taug(jl,ngs3+14,lay)*0.943_JPRB
        enddo
      enddo
      !$ACC END PARALLEL
      !$OMP END TARGET TEAMS DISTRIBUTE PARALLEL DO

      IF (llaytrop_max /= llaytrop_min) THEN
        ! Mixed loop
        ! Lower atmosphere part
        !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO COLLAPSE(2) PRIVATE(speccomb, specparm, specmult, js, fs, speccomb1, specparm1, &
        !$OMP   specmult1, js1, fs1, speccomb_planck, specparm_planck, specmult_planck, jpl, fpl, ind0, ind1, inds, &
        !$OMP   indf, p, p4, fk0, fk1, fk2, fac000, fac100, fac200, fac010, fac110, fac210, fac001, fac101, fac201, &
        !$OMP   fac011, fac111, fac211, tau_major, tau_major1, tauself, taufor)
        !$ACC PARALLEL DEFAULT(NONE) ASYNC(1)
        !$ACC LOOP GANG VECTOR COLLAPSE(2) PRIVATE(speccomb, specparm, specmult, js, fs, speccomb1, specparm1, &
        !$ACC   specmult1, js1, fs1, speccomb_planck, specparm_planck, specmult_planck, jpl, fpl, ind0, ind1, inds, &
        !$ACC   indf, p, p4, fk0, fk1, fk2, fac000, fac100, fac200, fac010, fac110, fac210, fac001, fac101, fac201, &
        !$ACC   fac011, fac111, fac211, tau_major, tau_major1)
        DO lay = llaytrop_min+1, llaytrop_max
#if defined(_OPENACC) || defined(OMPGPU)
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

            speccomb_planck = colh2o(jl,lay)+refrat_planck_a*colco2(jl,lay)
            specparm_planck = MIN(colh2o(jl,lay)/speccomb_planck,oneminus)
            specmult_planck = 8._JPRB*specparm_planck
            jpl= 1 + int(specmult_planck)
            fpl = MOD1(specmult_planck)

            ind0 = ((jp(jl,lay)-1)*5+(jt(jl,lay)-1))*nspa(4) + js
            ind1 = (jp(jl,lay)*5+(jt1(jl,lay)-1))*nspa(4) + js1
            inds = indself(jl,lay)
            indf = indfor(jl,lay)

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
            do ig = 1, ng4
              tau_major(ig) = speccomb *    &
              (fac000 * absa(ind0,ig)    + &
                fac100 * absa(ind0+1,ig)  + &
                fac200 * absa(ind0+2,ig)  + &
                fac010 * absa(ind0+9,ig)  + &
                fac110 * absa(ind0+10,ig) + &
                fac210 * absa(ind0+11,ig))
            enddo
#else
!$NEC unroll(NG4)
              tau_major(1:ng4) = speccomb *    &
              (fac000 * absa(ind0,1:ng4)    + &
                fac100 * absa(ind0+1,1:ng4)  + &
                fac200 * absa(ind0+2,1:ng4)  + &
                fac010 * absa(ind0+9,1:ng4)  + &
                fac110 * absa(ind0+10,1:ng4) + &
                fac210 * absa(ind0+11,1:ng4))
#endif
            else if (specparm .gt. 0.875_JPRB) then
#if defined(__amdflang__) && defined(OMPGPU)
            do ig = 1, ng4
              tau_major(ig) = speccomb *   &
              (fac200 * absa(ind0-1,ig) + &
                fac100 * absa(ind0,ig)   + &
                fac000 * absa(ind0+1,ig) + &
                fac210 * absa(ind0+8,ig) + &
                fac110 * absa(ind0+9,ig) + &
                fac010 * absa(ind0+10,ig))
            enddo
#else
!$NEC unroll(NG4)
              tau_major(1:ng4) = speccomb *   &
              (fac200 * absa(ind0-1,1:ng4) + &
                fac100 * absa(ind0,1:ng4)   + &
                fac000 * absa(ind0+1,1:ng4) + &
                fac210 * absa(ind0+8,1:ng4) + &
                fac110 * absa(ind0+9,1:ng4) + &
                fac010 * absa(ind0+10,1:ng4))
#endif
            else
#if defined(__amdflang__) && defined(OMPGPU)
            do ig = 1, ng4
              tau_major(ig) = speccomb *   &
                (fac000 * absa(ind0,ig)   + &
                fac100 * absa(ind0+1,ig) + &
                fac010 * absa(ind0+9,ig) + &
                fac110 * absa(ind0+10,ig))
            enddo
#else
!$NEC unroll(NG4)
              tau_major(1:ng4) = speccomb *   &
                (fac000 * absa(ind0,1:ng4)   + &
                fac100 * absa(ind0+1,1:ng4) + &
                fac010 * absa(ind0+9,1:ng4) + &
                fac110 * absa(ind0+10,1:ng4))
#endif
            endif

            if (specparm1 .lt. 0.125_JPRB) then
#if defined(__amdflang__) && defined(OMPGPU)
            do ig = 1, ng4
              tau_major1(ig) = speccomb1 *  &
              (fac001 * absa(ind1,ig)    + &
                fac101 * absa(ind1+1,ig)  + &
                fac201 * absa(ind1+2,ig)  + &
                fac011 * absa(ind1+9,ig)  + &
                fac111 * absa(ind1+10,ig) + &
                fac211 * absa(ind1+11,ig))
            enddo
#else
!$NEC unroll(NG4)
              tau_major1(1:ng4) = speccomb1 *  &
              (fac001 * absa(ind1,1:ng4)    + &
                fac101 * absa(ind1+1,1:ng4)  + &
                fac201 * absa(ind1+2,1:ng4)  + &
                fac011 * absa(ind1+9,1:ng4)  + &
                fac111 * absa(ind1+10,1:ng4) + &
                fac211 * absa(ind1+11,1:ng4))
#endif
            else if (specparm1 .gt. 0.875_JPRB) then
#if defined(__amdflang__) && defined(OMPGPU)
            do ig = 1, ng4
              tau_major1(ig) = speccomb1 * &
              (fac201 * absa(ind1-1,ig) + &
                fac101 * absa(ind1,ig)   + &
                fac001 * absa(ind1+1,ig) + &
                fac211 * absa(ind1+8,ig) + &
                fac111 * absa(ind1+9,ig) + &
                fac011 * absa(ind1+10,ig))
            enddo
#else
!$NEC unroll(NG4)
              tau_major1(1:ng4) = speccomb1 * &
              (fac201 * absa(ind1-1,1:ng4) + &
                fac101 * absa(ind1,1:ng4)   + &
                fac001 * absa(ind1+1,1:ng4) + &
                fac211 * absa(ind1+8,1:ng4) + &
                fac111 * absa(ind1+9,1:ng4) + &
                fac011 * absa(ind1+10,1:ng4))
#endif
            else
#if defined(__amdflang__) && defined(OMPGPU)
            do ig = 1, ng4
              tau_major1(ig) = speccomb1 * &
              (fac001 * absa(ind1,ig)   + &
                fac101 * absa(ind1+1,ig) + &
                fac011 * absa(ind1+9,ig) + &
                fac111 * absa(ind1+10,ig))
            enddo
#else
!$NEC unroll(NG4)
              tau_major1(1:ng4) = speccomb1 * &
              (fac001 * absa(ind1,1:ng4)   + &
                fac101 * absa(ind1+1,1:ng4) + &
                fac011 * absa(ind1+9,1:ng4) + &
                fac111 * absa(ind1+10,1:ng4))
#endif
            endif

!$NEC unroll(NG4)
            !$ACC LOOP SEQ PRIVATE(tauself, taufor)
            do ig = 1, ng4
              tauself = selffac(jl,lay)* (selfref(inds,ig) + selffrac(jl,lay) * &
                  (selfref(inds+1,ig) - selfref(inds,ig)))
              taufor =  forfac(jl,lay) * (forref(indf,ig) + forfrac(jl,lay) * &
                  (forref(indf+1,ig) - forref(indf,ig)))

              taug(jl,ngs3+ig,lay) = tau_major(ig) + tau_major1(ig) &
                  + tauself + taufor
              fracs(jl,ngs3+ig,lay) = fracrefa(ig,jpl) + fpl * &
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

            speccomb = colo3(jl,lay) + rat_o3co2(jl,lay)*colco2(jl,lay)
            specparm = MIN(colo3(jl,lay)/speccomb,oneminus)
            specmult = 4._JPRB*(specparm)
            js = 1 + int(specmult)
            fs = MOD1(specmult)

            speccomb1 = colo3(jl,lay) + rat_o3co2_1(jl,lay)*colco2(jl,lay)
            specparm1 = MIN(colo3(jl,lay)/speccomb1,oneminus)
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

            speccomb_planck = colo3(jl,lay)+refrat_planck_b*colco2(jl,lay)
            specparm_planck = MIN(colo3(jl,lay)/speccomb_planck,oneminus)
            specmult_planck = 4._JPRB*specparm_planck
            jpl= 1 + int(specmult_planck)
            fpl = MOD1(specmult_planck)

            ind0 = ((jp(jl,lay)-13)*5+(jt(jl,lay)-1))*nspb(4) + js
            ind1 = ((jp(jl,lay)-12)*5+(jt1(jl,lay)-1))*nspb(4) + js1
!$NEC unroll(NG4)
            !$ACC LOOP SEQ
            do ig = 1, ng4
              taug(jl,ngs3+ig,lay) =  speccomb * &
                  (fac000 * absb(ind0,ig) + &
                  fac100 * absb(ind0+1,ig) + &
                  fac010 * absb(ind0+5,ig) + &
                  fac110 * absb(ind0+6,ig)) &
                  + speccomb1 * &
                  (fac001 * absb(ind1,ig) +  &
                  fac101 * absb(ind1+1,ig) + &
                  fac011 * absb(ind1+5,ig) + &
                  fac111 * absb(ind1+6,ig))
              fracs(jl,ngs3+ig,lay) = fracrefb(ig,jpl) + fpl * &
                  (fracrefb(ig,jpl+1)-fracrefb(ig,jpl))
            enddo
#if defined(_OPENACC) || defined(OMPGPU)
#else
          enddo

          ! Empirical modification to code to improve stratospheric cooling rates
          ! for co2.  Revised to apply weighting for g-point reduction in this band.
!$NEC ivdep
          do ixp = 1, ixc0
            jl = ixhigh(ixp,lay)
#endif
            taug(jl,ngs3+8,lay)=taug(jl,ngs3+8,lay)*0.92_JPRB
            taug(jl,ngs3+9,lay)=taug(jl,ngs3+9,lay)*0.88_JPRB
            taug(jl,ngs3+10,lay)=taug(jl,ngs3+10,lay)*1.07_JPRB
            taug(jl,ngs3+11,lay)=taug(jl,ngs3+11,lay)*1.1_JPRB
            taug(jl,ngs3+12,lay)=taug(jl,ngs3+12,lay)*0.99_JPRB
            taug(jl,ngs3+13,lay)=taug(jl,ngs3+13,lay)*0.88_JPRB
            taug(jl,ngs3+14,lay)=taug(jl,ngs3+14,lay)*0.943_JPRB
#if defined(_OPENACC) || defined(OMPGPU)
           endif
#endif
          enddo

        ENDDO
        !$ACC END PARALLEL
        !$OMP END TARGET TEAMS DISTRIBUTE PARALLEL DO

      END IF

      !$ACC END DATA
      !$OMP END TARGET DATA

END SUBROUTINE RRTM_TAUMOL4
