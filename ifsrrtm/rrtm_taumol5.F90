!----------------------------------------------------------------------------
SUBROUTINE RRTM_TAUMOL5 (KIDIA,KFDIA,KLEV,taug,wx,&
 & P_TAUAERL,fac00,fac01,fac10,fac11,forfac,forfrac,indfor,jp,jt,jt1,oneminus,&
 & colh2o,colco2, colo3,laytrop,selffac,selffrac,indself,fracs, &
 & rat_h2oco2, rat_h2oco2_1, rat_o3co2, rat_o3co2_1,minorfrac,indminor)

!     BAND 5:  700-820 cm-1 (low - H2O,CO2; high - O3,CO2)

!     AUTHOR.
!     -------
!      JJMorcrette, ECMWF

!     MODIFICATIONS.
!     --------------
!      M.Hamrud      01-Oct-2003 CY28 Cleaning
!      NEC           25-Oct-2007 Optimisations
!      JJMorcrette 20110613 flexible number of g-points
!      ABozzo 201306 updated to rrtmg v4.85
!      band 5:  700-820 cm-1 (low key - h2o,co2; low minor - o3, ccl4)
!                            (high key - o3,co2)
! ---------------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM     ,JPRB

USE PARRRTM  , ONLY : JPBAND ,JPXSEC
USE YOERRTM  , ONLY : JPGPT  ,NG5    ,NGS4
USE YOERRTWN , ONLY : NSPA   ,NSPB
USE YOERRTA5 , ONLY : ABSA   ,ABSB   ,CCL4   , FRACREFA, FRACREFB,SELFREF,FORREF, &
 & KA_MO3
USE YOERRTRF, ONLY : CHI_MLS

IMPLICIT NONE

INTEGER(KIND=JPIM),INTENT(IN)    :: KIDIA
INTEGER(KIND=JPIM),INTENT(IN)    :: KFDIA
INTEGER(KIND=JPIM),INTENT(IN)    :: KLEV
REAL(KIND=JPRB)   ,INTENT(INOUT) :: taug(KIDIA:KFDIA,JPGPT,KLEV)
REAL(KIND=JPRB)   ,INTENT(IN)    :: wx(KIDIA:KFDIA,JPXSEC,KLEV) ! Amount of trace gases
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
REAL(KIND=JPRB)   ,INTENT(IN)   :: forfrac(KIDIA:KFDIA,KLEV)
REAL(KIND=JPRB)   ,INTENT(IN)   :: forfac(KIDIA:KFDIA,KLEV)
REAL(KIND=JPRB)   ,INTENT(IN)   :: minorfrac(KIDIA:KFDIA,KLEV)
INTEGER(KIND=JPIM),INTENT(IN)   :: indminor(KIDIA:KFDIA,KLEV)
! ---------------------------------------------------------------------------

REAL(KIND=JPRB) :: speccomb,speccomb1, &
& speccomb_mo3, speccomb_planck
INTEGER(KIND=JPIM) :: ind0,ind1,inds,indf,indm

INTEGER(KIND=JPIM) :: IG, JS, lay, JS1, JPL, JMO3

REAL(KIND=JPRB) :: refrat_planck_a, refrat_planck_b,refrat_m_a

REAL(KIND=JPRB) ::  fac000, fac100, fac200,&
 & fac010, fac110, fac210, &
 & fac001, fac101, fac201, &
 & fac011, fac111, fac211
REAL(KIND=JPRB) :: p, p4, fk0, fk1, fk2
REAL(KIND=JPRB) :: taufor,tauself,tau_major(ng5),tau_major1(ng5), o3m1, o3m2, abso3

REAL(KIND=JPRB) :: fs, specmult, specparm, &
& fs1, specmult1, specparm1, &
& fpl, specmult_PLANCK, specparm_PLANCK, &
& fmo3, specmult_MO3, specparm_MO3
    !     local integer arrays
    INTEGER(KIND=JPIM) :: laytrop_min, laytrop_max
    integer(KIND=JPIM) :: ixc(KLEV), ixlow(KFDIA,KLEV), ixhigh(KFDIA,KLEV)
    INTEGER(KIND=JPIM) :: ich, icl, ixc0, ixp, jc, jl

#define MOD1(x) ((x) - AINT((x)))

    !$ACC DATA PRESENT(taug, wx, P_TAUAERL, fac00, fac01, fac10, fac11, jp, jt, &
    !$ACC             jt1, colh2o, colco2, colo3, laytrop, selffac, selffrac, &
    !$ACC             indself, fracs, rat_h2oco2, rat_h2oco2_1, rat_o3co2, &
    !$ACC             rat_o3co2_1, indfor, forfrac, forfac, minorfrac, indminor)

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

! Minor gas mapping level :
!     lower - o3, p = 317.34 mbar, t = 240.77 k
!     lower - ccl4

!     Compute the optical depth by interpolating in ln(pressure),
!     temperature, and appropriate species.  Below LAYTROP, the water
!     vapor self-continuum and foreign continuum is
!     interpolated (in temperature) separately.

      ! P = 473.420 mb
      refrat_planck_a = chi_mls(1,5)/chi_mls(2,5)

      ! P = 0.2369 mb
      refrat_planck_b = chi_mls(3,43)/chi_mls(2,43)

      ! P = 317.3480
      refrat_m_a = chi_mls(1,7)/chi_mls(2,7)

      ! Compute the optical depth by interpolating in ln(pressure) and
      ! temperature, and appropriate species.  Below laytrop, the
      ! water vapor self-continuum and foreign continuum is
      ! interpolated (in temperature) separately.

      ! Lower atmosphere loop
      !$ACC WAIT
      !$ACC PARALLEL DEFAULT(NONE) ASYNC(1)
      !$ACC LOOP GANG VECTOR COLLAPSE(2) PRIVATE(speccomb,speccomb1, speccomb_mo3, speccomb_planck, ind0, ind1, inds, &
      !$ACC   indf, indm, js, js1, jpl, jmo3, fac000, fac100, fac200, fac010, fac110, fac210, fac001, fac101, fac201, &
      !$ACC   fac011, fac111, fac211, p, p4, fk0, fk1, fk2, fs, specmult, specparm, fs1, specmult1, specparm1, &
      !$ACC   fpl, specmult_PLANCK, specparm_PLANCK, fmo3, specmult_MO3, specparm_MO3, tau_major, tau_major1)
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

          speccomb_mo3 = colh2o(jl,lay) + refrat_m_a*colco2(jl,lay)
          specparm_mo3 = MIN(colh2o(jl,lay)/speccomb_mo3,oneminus)
          specmult_mo3 = 8._JPRB*specparm_mo3
          jmo3 = 1 + int(specmult_mo3)
          fmo3 = MOD1(specmult_mo3)

          speccomb_planck = colh2o(jl,lay)+refrat_planck_a*colco2(jl,lay)
          specparm_planck = MIN(colh2o(jl,lay)/speccomb_planck,oneminus)
          specmult_planck = 8._JPRB*specparm_planck
          jpl= 1 + int(specmult_planck)
          fpl = MOD1(specmult_planck)

          ind0 = ((jp(jl,lay)-1)*5+(jt(jl,lay)-1))*nspa(5) + js
          ind1 = (jp(jl,lay)*5+(jt1(jl,lay)-1))*nspa(5) + js1
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
!$NEC unroll(NG5)
            tau_major(1:ng5) = speccomb *    &
             (fac000 * absa(ind0,1:ng5)    + &
              fac100 * absa(ind0+1,1:ng5)  + &
              fac200 * absa(ind0+2,1:ng5)  + &
              fac010 * absa(ind0+9,1:ng5)  + &
              fac110 * absa(ind0+10,1:ng5) + &
              fac210 * absa(ind0+11,1:ng5))
          else if (specparm .gt. 0.875_JPRB) then
!$NEC unroll(NG5)
            tau_major(1:ng5) = speccomb *   &
             (fac200 * absa(ind0-1,1:ng5) + &
              fac100 * absa(ind0,1:ng5)   + &
              fac000 * absa(ind0+1,1:ng5) + &
              fac210 * absa(ind0+8,1:ng5) + &
              fac110 * absa(ind0+9,1:ng5) + &
              fac010 * absa(ind0+10,1:ng5))
          else
!$NEC unroll(NG5)
            tau_major(1:ng5) = speccomb *   &
             (fac000 * absa(ind0,1:ng5)   + &
              fac100 * absa(ind0+1,1:ng5) + &
              fac010 * absa(ind0+9,1:ng5) + &
              fac110 * absa(ind0+10,1:ng5))
          endif

          if (specparm1 .lt. 0.125_JPRB) then
!$NEC unroll(NG5)
            tau_major1(1:ng5) = speccomb1 *  &
             (fac001 * absa(ind1,1:ng5)    + &
              fac101 * absa(ind1+1,1:ng5)  + &
              fac201 * absa(ind1+2,1:ng5)  + &
              fac011 * absa(ind1+9,1:ng5)  + &
              fac111 * absa(ind1+10,1:ng5) + &
              fac211 * absa(ind1+11,1:ng5))
          else if (specparm1 .gt. 0.875_JPRB) then
!$NEC unroll(NG5)
            tau_major1(1:ng5) = speccomb1 * &
             (fac201 * absa(ind1-1,1:ng5) + &
              fac101 * absa(ind1,1:ng5)   + &
              fac001 * absa(ind1+1,1:ng5) + &
              fac211 * absa(ind1+8,1:ng5) + &
              fac111 * absa(ind1+9,1:ng5) + &
              fac011 * absa(ind1+10,1:ng5))
          else
!$NEC unroll(NG5)
            tau_major1(1:ng5) = speccomb1 * &
             (fac001 * absa(ind1,1:ng5)   + &
              fac101 * absa(ind1+1,1:ng5) + &
              fac011 * absa(ind1+9,1:ng5) + &
              fac111 * absa(ind1+10,1:ng5))
          endif

          !$ACC LOOP SEQ PRIVATE(taufor,tauself, o3m1, o3m2, abso3)
!$NEC unroll(NG5)
          do ig = 1, ng5
            tauself = selffac(jl,lay) * (selfref(inds,ig) + selffrac(jl,lay) * &
                 (selfref(inds+1,ig) - selfref(inds,ig)))
            taufor =  forfac(jl,lay) * (forref(indf,ig) + forfrac(jl,lay) * &
                 (forref(indf+1,ig) - forref(indf,ig)))
            o3m1 = ka_mo3(jmo3,indm,ig) + fmo3 * &
                 (ka_mo3(jmo3+1,indm,ig)-ka_mo3(jmo3,indm,ig))
            o3m2 = ka_mo3(jmo3,indm+1,ig) + fmo3 * &
                 (ka_mo3(jmo3+1,indm+1,ig)-ka_mo3(jmo3,indm+1,ig))
            abso3 = o3m1 + minorfrac(jl,lay)*(o3m2-o3m1)

            taug(jl,ngs4+ig,lay) = tau_major(ig) + tau_major1(ig) &
                 + tauself + taufor &
                 + abso3*colo3(jl,lay) &
                 + wx(jl,1,lay) * ccl4(ig)
            fracs(jl,ngs4+ig,lay) = fracrefa(ig,jpl) + fpl * &
                 (fracrefa(ig,jpl+1)-fracrefa(ig,jpl))
          enddo
        enddo

      enddo
      !$ACC END PARALLEL

      ! Upper atmosphere loop
      !$ACC PARALLEL DEFAULT(NONE) ASYNC(1)
      !$ACC LOOP GANG VECTOR COLLAPSE(2) PRIVATE(speccomb, speccomb1, speccomb_planck, ind0, ind1, indf, indm, js, &
      !$ACC   js1, jpl, fac000, fac100, fac010, fac110, fac001, fac101, fac011, fac111, fs, specmult, specparm, fs1, &
      !$ACC   specmult1, specparm1, fpl, specmult_PLANCK, specparm_PLANCK)
      do lay = laytrop_max+1, KLEV
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

          ind0 = ((jp(jl,lay)-13)*5+(jt(jl,lay)-1))*nspb(5) + js
          ind1 = ((jp(jl,lay)-12)*5+(jt1(jl,lay)-1))*nspb(5) + js1
          !$ACC LOOP SEQ
!$NEC unroll(NG5)
          do ig = 1, ng5
            taug(jl,ngs4+ig,lay) = speccomb * &
                 (fac000 * absb(ind0,ig) + &
                 fac100 * absb(ind0+1,ig) + &
                 fac010 * absb(ind0+5,ig) + &
                 fac110 * absb(ind0+6,ig)) &
                 + speccomb1 * &
                 (fac001 * absb(ind1,ig) + &
                 fac101 * absb(ind1+1,ig) + &
                 fac011 * absb(ind1+5,ig) + &
                 fac111 * absb(ind1+6,ig))  &
                 + wx(jl,1,lay) * ccl4(ig)
            fracs(jl,ngs4+ig,lay) = fracrefb(ig,jpl) + fpl * &
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
        !$ACC   specmult1, js1, fs1, speccomb_mo3, specparm_mo3, specmult_mo3, jmo3, fmo3, speccomb_planck, &
        !$ACC   specparm_planck, specmult_planck, jpl, fpl, ind0, ind1, inds, indf, indm, p, p4, fk0, fac000, fac100, &
        !$ACC   fac010, fac110, fac001, fac101, fac011, fac111, fac200, fac201, fac210, fac211, tau_major, tau_major1)
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

            speccomb_mo3 = colh2o(jl,lay) + refrat_m_a*colco2(jl,lay)
            specparm_mo3 = MIN(colh2o(jl,lay)/speccomb_mo3,oneminus)
            specmult_mo3 = 8._JPRB*specparm_mo3
            jmo3 = 1 + int(specmult_mo3)
            fmo3 = MOD1(specmult_mo3)

            speccomb_planck = colh2o(jl,lay)+refrat_planck_a*colco2(jl,lay)
            specparm_planck = MIN(colh2o(jl,lay)/speccomb_planck,oneminus)
            specmult_planck = 8._JPRB*specparm_planck
            jpl= 1 + int(specmult_planck)
            fpl = MOD1(specmult_planck)

            ind0 = ((jp(jl,lay)-1)*5+(jt(jl,lay)-1))*nspa(5) + js
            ind1 = (jp(jl,lay)*5+(jt1(jl,lay)-1))*nspa(5) + js1
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
!$NEC unroll(NG5)
              tau_major(1:ng5) = speccomb *    &
              (fac000 * absa(ind0,1:ng5)    + &
                fac100 * absa(ind0+1,1:ng5)  + &
                fac200 * absa(ind0+2,1:ng5)  + &
                fac010 * absa(ind0+9,1:ng5)  + &
                fac110 * absa(ind0+10,1:ng5) + &
                fac210 * absa(ind0+11,1:ng5))
            else if (specparm .gt. 0.875_JPRB) then
!$NEC unroll(NG5)
              tau_major(1:ng5) = speccomb *   &
              (fac200 * absa(ind0-1,1:ng5) + &
                fac100 * absa(ind0,1:ng5)   + &
                fac000 * absa(ind0+1,1:ng5) + &
                fac210 * absa(ind0+8,1:ng5) + &
                fac110 * absa(ind0+9,1:ng5) + &
                fac010 * absa(ind0+10,1:ng5))
            else
!$NEC unroll(NG5)
              tau_major(1:ng5) = speccomb *   &
              (fac000 * absa(ind0,1:ng5)   + &
                fac100 * absa(ind0+1,1:ng5) + &
                fac010 * absa(ind0+9,1:ng5) + &
                fac110 * absa(ind0+10,1:ng5))
            endif

            if (specparm1 .lt. 0.125_JPRB) then
!$NEC unroll(NG5)
              tau_major1(1:ng5) = speccomb1 *  &
              (fac001 * absa(ind1,1:ng5)    + &
                fac101 * absa(ind1+1,1:ng5)  + &
                fac201 * absa(ind1+2,1:ng5)  + &
                fac011 * absa(ind1+9,1:ng5)  + &
                fac111 * absa(ind1+10,1:ng5) + &
                fac211 * absa(ind1+11,1:ng5))
            else if (specparm1 .gt. 0.875_JPRB) then
!$NEC unroll(NG5)
              tau_major1(1:ng5) = speccomb1 * &
              (fac201 * absa(ind1-1,1:ng5) + &
                fac101 * absa(ind1,1:ng5)   + &
                fac001 * absa(ind1+1,1:ng5) + &
                fac211 * absa(ind1+8,1:ng5) + &
                fac111 * absa(ind1+9,1:ng5) + &
                fac011 * absa(ind1+10,1:ng5))
            else
!$NEC unroll(NG5)
              tau_major1(1:ng5) = speccomb1 * &
              (fac001 * absa(ind1,1:ng5)   + &
                fac101 * absa(ind1+1,1:ng5) + &
                fac011 * absa(ind1+9,1:ng5) + &
                fac111 * absa(ind1+10,1:ng5))
            endif

!$NEC unroll(NG5)
            !$ACC LOOP SEQ PRIVATE(tauself, taufor, o3m1, o3m2, abso3)
            do ig = 1, ng5
              tauself = selffac(jl,lay) * (selfref(inds,ig) + selffrac(jl,lay) * &
                  (selfref(inds+1,ig) - selfref(inds,ig)))
              taufor =  forfac(jl,lay) * (forref(indf,ig) + forfrac(jl,lay) * &
                  (forref(indf+1,ig) - forref(indf,ig)))
              o3m1 = ka_mo3(jmo3,indm,ig) + fmo3 * &
                  (ka_mo3(jmo3+1,indm,ig)-ka_mo3(jmo3,indm,ig))
              o3m2 = ka_mo3(jmo3,indm+1,ig) + fmo3 * &
                  (ka_mo3(jmo3+1,indm+1,ig)-ka_mo3(jmo3,indm+1,ig))
              abso3 = o3m1 + minorfrac(jl,lay)*(o3m2-o3m1)

              taug(jl,ngs4+ig,lay) = tau_major(ig) + tau_major1(ig) &
                  + tauself + taufor &
                  + abso3*colo3(jl,lay) &
                  + wx(jl,1,lay) * ccl4(ig)
              fracs(jl,ngs4+ig,lay) = fracrefa(ig,jpl) + fpl * &
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

            ind0 = ((jp(jl,lay)-13)*5+(jt(jl,lay)-1))*nspb(5) + js
            ind1 = ((jp(jl,lay)-12)*5+(jt1(jl,lay)-1))*nspb(5) + js1
!$NEC unroll(NG5)
            !$ACC LOOP SEQ
            do ig = 1, ng5
              taug(jl,ngs4+ig,lay) = speccomb * &
                  (fac000 * absb(ind0,ig) + &
                  fac100 * absb(ind0+1,ig) + &
                  fac010 * absb(ind0+5,ig) + &
                  fac110 * absb(ind0+6,ig)) &
                  + speccomb1 * &
                  (fac001 * absb(ind1,ig) + &
                  fac101 * absb(ind1+1,ig) + &
                  fac011 * absb(ind1+5,ig) + &
                  fac111 * absb(ind1+6,ig))  &
                  + wx(jl,1,lay) * ccl4(ig)
              fracs(jl,ngs4+ig,lay) = fracrefb(ig,jpl) + fpl * &
                  (fracrefb(ig,jpl+1)-fracrefb(ig,jpl))
            enddo
#ifdef _OPENACC
           endif
#endif
          enddo

        enddo
        !$ACC END PARALLEL

    END IF

    !$ACC END DATA

END SUBROUTINE RRTM_TAUMOL5
