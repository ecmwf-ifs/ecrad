!----------------------------------------------------------------------------
SUBROUTINE RRTM_TAUMOL9 (KIDIA,KFDIA,KLEV,taug,&
 & P_TAUAERL,fac00,fac01,fac10,fac11,forfac,forfrac,indfor,jp,jt,jt1,oneminus,&
 & colh2o,coln2o,colch4,coldry,laytrop,K_LAYSWTCH,K_LAYLOW,selffac,selffrac,indself,fracs, &
 & rat_h2och4,rat_h2och4_1,minorfrac,indminor,laytrop_min,laytrop_max)

!     BAND 9:  1180-1390 cm-1 (low - H2O,CH4; high - CH4)

!     AUTHOR.
!     -------
!      JJMorcrette, ECMWF

!     MODIFICATIONS.
!     --------------
!      M.Hamrud      01-Oct-2003 CY28 Cleaning
!      NEC           25-Oct-2007 Optimisations
!      JJMorcrette 20110613 flexible number of g-points
!      201306 ABozzo updated to rrtmg v4.85
!     band 9:  1180-1390 cm-1 (low key - h2o,ch4; low minor - n2o)
!                             (high key - ch4; high minor - n2o)
! ---------------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM     ,JPRB

USE PARRRTM  , ONLY : JPBAND
USE YOERRTM  , ONLY : JPGPT  ,NG9   ,NGS8
USE YOERRTWN , ONLY :      NSPA   ,NSPB
USE YOERRTA9 , ONLY : ABSA   ,ABSB   ,FRACREFA, FRACREFB,SELFREF,FORREF,KA_MN2O, KB_MN2O
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
REAL(KIND=JPRB)   ,INTENT(IN)    :: colch4(KIDIA:KFDIA,KLEV)
REAL(KIND=JPRB)   ,INTENT(IN)    :: coldry(KIDIA:KFDIA,KLEV)
INTEGER(KIND=JPIM),INTENT(IN)    :: laytrop(KIDIA:KFDIA)
INTEGER(KIND=JPIM),INTENT(IN)    :: K_LAYSWTCH(KIDIA:KFDIA)
INTEGER(KIND=JPIM),INTENT(IN)    :: K_LAYLOW(KIDIA:KFDIA)
REAL(KIND=JPRB)   ,INTENT(IN)    :: selffac(KIDIA:KFDIA,KLEV)
REAL(KIND=JPRB)   ,INTENT(IN)    :: selffrac(KIDIA:KFDIA,KLEV)
INTEGER(KIND=JPIM),INTENT(IN)    :: indself(KIDIA:KFDIA,KLEV)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: fracs(KIDIA:KFDIA,JPGPT,KLEV)

REAL(KIND=JPRB)   ,INTENT(IN)   :: rat_h2och4(KIDIA:KFDIA,KLEV)
REAL(KIND=JPRB)   ,INTENT(IN)   :: rat_h2och4_1(KIDIA:KFDIA,KLEV)
INTEGER(KIND=JPIM),INTENT(IN)   :: indfor(KIDIA:KFDIA,KLEV)
REAL(KIND=JPRB)   ,INTENT(IN)   :: forfac(KIDIA:KFDIA,KLEV)
REAL(KIND=JPRB)   ,INTENT(IN)   :: forfrac(KIDIA:KFDIA,KLEV)
REAL(KIND=JPRB)   ,INTENT(IN)   :: minorfrac(KIDIA:KFDIA,KLEV)
INTEGER(KIND=JPIM),INTENT(IN)   :: indminor(KIDIA:KFDIA,KLEV)
INTEGER(KIND=JPIM),INTENT(IN)   :: laytrop_min, laytrop_max



! ---------------------------------------------------------------------------

INTEGER(KIND=JPIM) :: ind0,ind1,inds,indf,indm
INTEGER(KIND=JPIM) :: IG, JS, lay, JS1, JMN2O, JPL

REAL(KIND=JPRB) :: speccomb,speccomb1,speccomb_mn2o, &
& speccomb_planck
REAL(KIND=JPRB) :: refrat_planck_a, refrat_m_a

REAL(KIND=JPRB) :: fs, specmult, specparm,  &
 & fs1, specmult1, specparm1,   &
 & fmn2o, specmult_mn2o, specparm_mn2o,   &
 & fpl, specmult_planck, specparm_planck

REAL(KIND=JPRB) :: adjfac,adjcoln2o,ratn2o,chi_n2o
REAL(KIND=JPRB) ::  fac000, fac100, fac200,&
 & fac010, fac110, fac210, &
 & fac001, fac101, fac201, &
 & fac011, fac111, fac211
REAL(KIND=JPRB) :: p, p4, fk0, fk1, fk2
REAL(KIND=JPRB) :: taufor,tauself,n2om1,n2om2,absn2o,tau_major(ng9),tau_major1(ng9)

    !     local integer arrays
    integer(KIND=JPIM) :: ixc(KLEV), ixlow(KFDIA,KLEV), ixhigh(KFDIA,KLEV)
    INTEGER(KIND=JPIM) :: ich, icl, ixc0, ixp, jc, jl

#define MOD1(x) ((x) - AINT((x)))

    !$ACC DATA PRESENT(taug, P_TAUAERL, fac00, fac01, fac10, fac11, jp, jt, jt1, &
    !$ACC             colh2o, coln2o, colch4, coldry, laytrop, K_LAYSWTCH, &
    !$ACC             K_LAYLOW, selffac, selffrac, indself, fracs, rat_h2och4, &
    !$ACC             rat_h2och4_1, indfor, forfac, forfrac, minorfrac, indminor)
    !$OMP TARGET DATA MAP(PRESENT, ALLOC: taug, P_TAUAERL, fac00, fac01, fac10, fac11, jp, jt, jt1, &
    !$OMP             colh2o, coln2o, colch4, coldry, laytrop, K_LAYSWTCH, &
    !$OMP             K_LAYLOW, selffac, selffrac, indself, fracs, rat_h2och4, &
    !$OMP             rat_h2och4_1, indfor, forfac, forfrac, minorfrac, indminor)

#if !defined(_OPENACC) && !defined(OMPGPU)
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
#endif

      ! P = 212 mb
      refrat_planck_a = chi_mls(1,9)/chi_mls(6,9)

      ! P = 706.272 mb
      refrat_m_a = chi_mls(1,3)/chi_mls(6,3)

      ! Compute the optical depth by interpolating in ln(pressure),
      ! temperature, and appropriate species.  Below laytrop, the water
      ! vapor self-continuum and foreign continuum is interpolated
      ! (in temperature) separately.

      ! Lower atmosphere loop
      !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO COLLAPSE(2) PRIVATE(ind0, ind1, inds, indf, indm, js, js1, jmn2o, jpl, speccomb, &
      !$OMP   speccomb1, speccomb_mn2o, speccomb_planck, fs, specmult, specparm, fs1, specmult1, specparm1, fmn2o, &
      !$OMP   specmult_mn2o, specparm_mn2o, fpl, specmult_planck, specparm_planck, adjfac, adjcoln2o, ratn2o, &
      !$OMP   chi_n2o, fac000, fac100, fac200, fac010, fac110, fac210, fac001, fac101, fac201, fac011, fac111, &
      !$OMP   fac211, p, p4, fk0, fk1, fk2, tau_major, tau_major1, taufor, tauself, n2om1, n2om2, absn2o) THREAD_LIMIT(128)
      !$ACC WAIT
      !$ACC PARALLEL DEFAULT(NONE) ASYNC(1)
      !$ACC LOOP GANG VECTOR COLLAPSE(2) PRIVATE(ind0, ind1, inds, indf, indm, js, js1, jmn2o, jpl, speccomb, &
      !$ACC   speccomb1, speccomb_mn2o, speccomb_planck, fs, specmult, specparm, fs1, specmult1, specparm1, fmn2o, &
      !$ACC   specmult_mn2o, specparm_mn2o, fpl, specmult_planck, specparm_planck, adjfac, adjcoln2o, ratn2o, &
      !$ACC   chi_n2o, fac000, fac100, fac200, fac010, fac110, fac210, fac001, fac101, fac201, fac011, fac111, &
      !$ACC   fac211, p, p4, fk0, fk1, fk2, tau_major, tau_major1)
      do lay = 1, laytrop_min
        do jl = KIDIA, KFDIA

          speccomb = colh2o(jl,lay) + rat_h2och4(jl,lay)*colch4(jl,lay)
          specparm = MIN(colh2o(jl,lay)/speccomb,oneminus)
          specmult = 8._JPRB*(specparm)
          js = 1 + int(specmult)
          fs = MOD1(specmult)

          speccomb1 = colh2o(jl,lay) + rat_h2och4_1(jl,lay)*colch4(jl,lay)
          specparm1 = MIN(colh2o(jl,lay)/speccomb1,oneminus)
          specmult1 = 8._JPRB*(specparm1)
          js1 = 1 + int(specmult1)
          fs1 = MOD1(specmult1)

          speccomb_mn2o = colh2o(jl,lay) + refrat_m_a*colch4(jl,lay)
          specparm_mn2o = MIN(colh2o(jl,lay)/speccomb_mn2o,oneminus)
          specmult_mn2o = 8._JPRB*specparm_mn2o
          jmn2o = 1 + int(specmult_mn2o)
          fmn2o = MOD1(specmult_mn2o)

          !  In atmospheres where the amount of N2O is too great to be considered
          !  a minor species, adjust the column amount of N2O by an empirical factor
          !  to obtain the proper contribution.
          chi_n2o = coln2o(jl,lay)/(coldry(jl,lay))
          ratn2o = 1.e20_JPRB*chi_n2o/chi_mls(4,jp(jl,lay)+1)
          if (ratn2o .gt. 1.5_JPRB) then
            adjfac = 0.5_JPRB+(ratn2o-0.5_JPRB)**0.65_JPRB
            adjcoln2o = adjfac*chi_mls(4,jp(jl,lay)+1)*coldry(jl,lay)*1.e-20_JPRB
          else
            adjcoln2o = coln2o(jl,lay)
          endif

          speccomb_planck = colh2o(jl,lay)+refrat_planck_a*colch4(jl,lay)
          specparm_planck = MIN(colh2o(jl,lay)/speccomb_planck,oneminus)
          specmult_planck = 8._JPRB*specparm_planck
          jpl= 1 + int(specmult_planck)
          fpl = MOD1(specmult_planck)

          ind0 = ((jp(jl,lay)-1)*5+(jt(jl,lay)-1))*nspa(9) + js
          ind1 = (jp(jl,lay)*5+(jt1(jl,lay)-1))*nspa(9) + js1
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
!$NEC unroll(NG9)
            tau_major(1:ng9) = speccomb *    &
             (fac000 * absa(ind0,1:ng9)    + &
              fac100 * absa(ind0+1,1:ng9)  + &
              fac200 * absa(ind0+2,1:ng9)  + &
              fac010 * absa(ind0+9,1:ng9)  + &
              fac110 * absa(ind0+10,1:ng9) + &
              fac210 * absa(ind0+11,1:ng9))
          else if (specparm .gt. 0.875_JPRB) then
!$NEC unroll(NG9)
            tau_major(1:ng9) = speccomb *   &
             (fac200 * absa(ind0-1,1:ng9) + &
              fac100 * absa(ind0,1:ng9)   + &
              fac000 * absa(ind0+1,1:ng9) + &
              fac210 * absa(ind0+8,1:ng9) + &
              fac110 * absa(ind0+9,1:ng9) + &
              fac010 * absa(ind0+10,1:ng9))
          else
!$NEC unroll(NG9)
            tau_major(1:ng9) = speccomb *   &
             (fac000 * absa(ind0,1:ng9)   + &
              fac100 * absa(ind0+1,1:ng9) + &
              fac010 * absa(ind0+9,1:ng9) + &
              fac110 * absa(ind0+10,1:ng9))
          endif

          if (specparm1 .lt. 0.125_JPRB) then
!$NEC unroll(NG9)
            tau_major1(1:ng9) = speccomb1 *  &
             (fac001 * absa(ind1,1:ng9)    + &
              fac101 * absa(ind1+1,1:ng9)  + &
              fac201 * absa(ind1+2,1:ng9)  + &
              fac011 * absa(ind1+9,1:ng9)  + &
              fac111 * absa(ind1+10,1:ng9) + &
              fac211 * absa(ind1+11,1:ng9))
          else if (specparm1 .gt. 0.875_JPRB) then
!$NEC unroll(NG9)
            tau_major1(1:ng9) = speccomb1 * &
             (fac201 * absa(ind1-1,1:ng9) + &
              fac101 * absa(ind1,1:ng9)   + &
              fac001 * absa(ind1+1,1:ng9) + &
              fac211 * absa(ind1+8,1:ng9) + &
              fac111 * absa(ind1+9,1:ng9) + &
              fac011 * absa(ind1+10,1:ng9))
          else
!$NEC unroll(NG9)
            tau_major1(1:ng9) = speccomb1 * &
             (fac001 * absa(ind1,1:ng9)   + &
              fac101 * absa(ind1+1,1:ng9) + &
              fac011 * absa(ind1+9,1:ng9) + &
              fac111 * absa(ind1+10,1:ng9))
          endif

          !$ACC LOOP SEQ PRIVATE(taufor, tauself, n2om1, n2om2, absn2o)
!$NEC unroll(NG9)
          do ig = 1, ng9
            tauself = selffac(jl,lay)* (selfref(inds,ig) + selffrac(jl,lay) * &
                 (selfref(inds+1,ig) - selfref(inds,ig)))
            taufor = forfac(jl,lay) * (forref(indf,ig) + forfrac(jl,lay) * &
                 (forref(indf+1,ig) - forref(indf,ig)))
            n2om1 = ka_mn2o(jmn2o,indm,ig) + fmn2o * &
                 (ka_mn2o(jmn2o+1,indm,ig) - ka_mn2o(jmn2o,indm,ig))
            n2om2 = ka_mn2o(jmn2o,indm+1,ig) + fmn2o * &
                 (ka_mn2o(jmn2o+1,indm+1,ig) - ka_mn2o(jmn2o,indm+1,ig))
            absn2o = n2om1 + minorfrac(jl,lay) * (n2om2 - n2om1)

            taug(jl,ngs8+ig,lay) = tau_major(ig) + tau_major1(ig) &
                 + tauself + taufor &
                 + adjcoln2o*absn2o
            fracs(jl,ngs8+ig,lay) = fracrefa(ig,jpl) + fpl * &
                 (fracrefa(ig,jpl+1)-fracrefa(ig,jpl))
          enddo
        enddo

      enddo
      !$ACC END PARALLEL
      !$OMP END TARGET TEAMS DISTRIBUTE PARALLEL DO

      ! Upper atmosphere loop
      !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO COLLAPSE(2) PRIVATE(ind0, ind1, indm, adjfac, adjcoln2o, ratn2o, chi_n2o, absn2o)
      !$ACC PARALLEL DEFAULT(NONE) ASYNC(1)
      !$ACC LOOP GANG VECTOR COLLAPSE(2) PRIVATE(ind0, ind1, indm, adjfac, adjcoln2o, ratn2o, chi_n2o)
      do lay = laytrop_max+1, KLEV
        do jl = KIDIA, KFDIA

          !  In atmospheres where the amount of N2O is too great to be considered
          !  a minor species, adjust the column amount of N2O by an empirical factor
          !  to obtain the proper contribution.
          chi_n2o = coln2o(jl,lay)/(coldry(jl,lay))
          ratn2o = 1.e20_JPRB*chi_n2o/chi_mls(4,jp(jl,lay)+1)
          if (ratn2o .gt. 1.5_JPRB) then
            adjfac = 0.5_JPRB+(ratn2o-0.5_JPRB)**0.65_JPRB
            adjcoln2o = adjfac*chi_mls(4,jp(jl,lay)+1)*coldry(jl,lay)*1.e-20_JPRB
          else
            adjcoln2o = coln2o(jl,lay)
          endif

          ind0 = ((jp(jl,lay)-13)*5+(jt(jl,lay)-1))*nspb(9) + 1
          ind1 = ((jp(jl,lay)-12)*5+(jt1(jl,lay)-1))*nspb(9) + 1
          indm = indminor(jl,lay)
          !$ACC LOOP SEQ PRIVATE(absn2o)
!$NEC unroll(NG9)
          do ig = 1, ng9
            absn2o = kb_mn2o(indm,ig) + minorfrac(jl,lay) * &
                 (kb_mn2o(indm+1,ig) - kb_mn2o(indm,ig))
            taug(jl,ngs8+ig,lay) = colch4(jl,lay) * &
                 (fac00(jl,lay) * absb(ind0,ig) + &
                 fac10(jl,lay) * absb(ind0+1,ig) + &
                 fac01(jl,lay) * absb(ind1,ig) +  &
                 fac11(jl,lay) * absb(ind1+1,ig)) &
                 + adjcoln2o*absn2o
            fracs(jl,ngs8+ig,lay) = fracrefb(ig)
          enddo
        enddo

      enddo
      !$ACC END PARALLEL
      !$OMP END TARGET TEAMS DISTRIBUTE PARALLEL DO

      IF (laytrop_max /= laytrop_min) THEN
        ! Mixed loop
        ! Lower atmosphere part
        !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO COLLAPSE(2) PRIVATE(speccomb, specparm, specmult, js, fs, speccomb1, specparm1, &
        !$OMP   specmult1, js1, fs1, speccomb_mn2o, specparm_mn2o, specmult_mn2o, jmn2o, fmn2o, chi_n2o, ratn2o, &
        !$OMP   adjfac, adjcoln2o, speccomb_planck, specparm_planck, specmult_planck, jpl, fpl, ind0, ind1, inds, &
        !$OMP   indf, indm, p, p4, fk0, fk1, fk2, fac000, fac100, fac200, fac010, fac110, fac210, fac001, fac101, &
        !$OMP   fac201, fac011, fac111, fac211, tau_major, tau_major1, tauself, taufor, n2om1, n2om2, absn2o)
        !$ACC PARALLEL DEFAULT(NONE) ASYNC(1)
        !$ACC LOOP GANG VECTOR COLLAPSE(2) PRIVATE(speccomb, specparm, specmult, js, fs, speccomb1, specparm1, &
        !$ACC   specmult1, js1, fs1, speccomb_mn2o, specparm_mn2o, specmult_mn2o, jmn2o, fmn2o, chi_n2o, ratn2o, &
        !$ACC   adjfac, adjcoln2o, speccomb_planck, specparm_planck, specmult_planck, jpl, fpl, ind0, ind1, inds, &
        !$ACC   indf, indm, p, p4, fk0, fk1, fk2, fac000, fac100, fac200, fac010, fac110, fac210, fac001, fac101, &
        !$ACC   fac201, fac011, fac111, fac211, tau_major, tau_major1)
        do lay = laytrop_min+1, laytrop_max
#if defined(_OPENACC) || defined(OMPGPU)
          do jl = KIDIA, KFDIA
            if ( lay <= laytrop(jl) ) then
#else
          ixc0 = ixc(lay)

!$NEC ivdep
          do ixp = 1, ixc0
            jl = ixlow(ixp,lay)
#endif

            speccomb = colh2o(jl,lay) + rat_h2och4(jl,lay)*colch4(jl,lay)
            specparm = MIN(colh2o(jl,lay)/speccomb,oneminus)
            specmult = 8._JPRB*(specparm)
            js = 1 + int(specmult)
            fs = MOD1(specmult)

            speccomb1 = colh2o(jl,lay) + rat_h2och4_1(jl,lay)*colch4(jl,lay)
            specparm1 = MIN(colh2o(jl,lay)/speccomb1,oneminus)
            specmult1 = 8._JPRB*(specparm1)
            js1 = 1 + int(specmult1)
            fs1 = MOD1(specmult1)

            speccomb_mn2o = colh2o(jl,lay) + refrat_m_a*colch4(jl,lay)
            specparm_mn2o = MIN(colh2o(jl,lay)/speccomb_mn2o,oneminus)
            specmult_mn2o = 8._JPRB*specparm_mn2o
            jmn2o = 1 + int(specmult_mn2o)
            fmn2o = MOD1(specmult_mn2o)

            !  In atmospheres where the amount of N2O is too great to be considered
            !  a minor species, adjust the column amount of N2O by an empirical factor
            !  to obtain the proper contribution.
            chi_n2o = coln2o(jl,lay)/(coldry(jl,lay))
            ratn2o = 1.e20_JPRB*chi_n2o/chi_mls(4,jp(jl,lay)+1)
            if (ratn2o .gt. 1.5_JPRB) then
              adjfac = 0.5_JPRB+(ratn2o-0.5_JPRB)**0.65_JPRB
              adjcoln2o = adjfac*chi_mls(4,jp(jl,lay)+1)*coldry(jl,lay)*1.e-20_JPRB
            else
              adjcoln2o = coln2o(jl,lay)
            endif

            speccomb_planck = colh2o(jl,lay)+refrat_planck_a*colch4(jl,lay)
            specparm_planck = MIN(colh2o(jl,lay)/speccomb_planck,oneminus)
            specmult_planck = 8._JPRB*specparm_planck
            jpl= 1 + int(specmult_planck)
            fpl = MOD1(specmult_planck)

            ind0 = ((jp(jl,lay)-1)*5+(jt(jl,lay)-1))*nspa(9) + js
            ind1 = (jp(jl,lay)*5+(jt1(jl,lay)-1))*nspa(9) + js1
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
!$NEC unroll(NG9)
              tau_major(1:ng9) = speccomb *    &
              (fac000 * absa(ind0,1:ng9)    + &
                fac100 * absa(ind0+1,1:ng9)  + &
                fac200 * absa(ind0+2,1:ng9)  + &
                fac010 * absa(ind0+9,1:ng9)  + &
                fac110 * absa(ind0+10,1:ng9) + &
                fac210 * absa(ind0+11,1:ng9))
            else if (specparm .gt. 0.875_JPRB) then
!$NEC unroll(NG9)
              tau_major(1:ng9) = speccomb *   &
              (fac200 * absa(ind0-1,1:ng9) + &
                fac100 * absa(ind0,1:ng9)   + &
                fac000 * absa(ind0+1,1:ng9) + &
                fac210 * absa(ind0+8,1:ng9) + &
                fac110 * absa(ind0+9,1:ng9) + &
                fac010 * absa(ind0+10,1:ng9))
            else
!$NEC unroll(NG9)
              tau_major(1:ng9) = speccomb *   &
              (fac000 * absa(ind0,1:ng9)   + &
                fac100 * absa(ind0+1,1:ng9) + &
                fac010 * absa(ind0+9,1:ng9) + &
                fac110 * absa(ind0+10,1:ng9))
            endif

            if (specparm1 .lt. 0.125_JPRB) then
!$NEC unroll(NG9)
              tau_major1(1:ng9) = speccomb1 *  &
              (fac001 * absa(ind1,1:ng9)    + &
                fac101 * absa(ind1+1,1:ng9)  + &
                fac201 * absa(ind1+2,1:ng9)  + &
                fac011 * absa(ind1+9,1:ng9)  + &
                fac111 * absa(ind1+10,1:ng9) + &
                fac211 * absa(ind1+11,1:ng9))
            else if (specparm1 .gt. 0.875_JPRB) then
!$NEC unroll(NG9)
              tau_major1(1:ng9) = speccomb1 * &
              (fac201 * absa(ind1-1,1:ng9) + &
                fac101 * absa(ind1,1:ng9)   + &
                fac001 * absa(ind1+1,1:ng9) + &
                fac211 * absa(ind1+8,1:ng9) + &
                fac111 * absa(ind1+9,1:ng9) + &
                fac011 * absa(ind1+10,1:ng9))
            else
!$NEC unroll(NG9)
              tau_major1(1:ng9) = speccomb1 * &
              (fac001 * absa(ind1,1:ng9)   + &
                fac101 * absa(ind1+1,1:ng9) + &
                fac011 * absa(ind1+9,1:ng9) + &
                fac111 * absa(ind1+10,1:ng9))
            endif

!$NEC unroll(NG9)
            !$ACC LOOP SEQ PRIVATE(tauself, taufor, n2om1, n2om2, absn2o)
            do ig = 1, ng9
              tauself = selffac(jl,lay)* (selfref(inds,ig) + selffrac(jl,lay) * &
                  (selfref(inds+1,ig) - selfref(inds,ig)))
              taufor = forfac(jl,lay) * (forref(indf,ig) + forfrac(jl,lay) * &
                  (forref(indf+1,ig) - forref(indf,ig)))
              n2om1 = ka_mn2o(jmn2o,indm,ig) + fmn2o * &
                  (ka_mn2o(jmn2o+1,indm,ig) - ka_mn2o(jmn2o,indm,ig))
              n2om2 = ka_mn2o(jmn2o,indm+1,ig) + fmn2o * &
                  (ka_mn2o(jmn2o+1,indm+1,ig) - ka_mn2o(jmn2o,indm+1,ig))
              absn2o = n2om1 + minorfrac(jl,lay) * (n2om2 - n2om1)

              taug(jl,ngs8+ig,lay) = tau_major(ig) + tau_major1(ig) &
                  + tauself + taufor &
                  + adjcoln2o*absn2o
              fracs(jl,ngs8+ig,lay) = fracrefa(ig,jpl) + fpl * &
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

            !  In atmospheres where the amount of N2O is too great to be considered
            !  a minor species, adjust the column amount of N2O by an empirical factor
            !  to obtain the proper contribution.
            chi_n2o = coln2o(jl,lay)/(coldry(jl,lay))
            ratn2o = 1.e20_JPRB*chi_n2o/chi_mls(4,jp(jl,lay)+1)
            if (ratn2o .gt. 1.5_JPRB) then
              adjfac = 0.5_JPRB+(ratn2o-0.5_JPRB)**0.65_JPRB
              adjcoln2o = adjfac*chi_mls(4,jp(jl,lay)+1)*coldry(jl,lay)*1.e-20_JPRB
            else
              adjcoln2o = coln2o(jl,lay)
            endif

            ind0 = ((jp(jl,lay)-13)*5+(jt(jl,lay)-1))*nspb(9) + 1
            ind1 = ((jp(jl,lay)-12)*5+(jt1(jl,lay)-1))*nspb(9) + 1
            indm = indminor(jl,lay)
!$NEC unroll(NG9)
            !$ACC LOOP SEQ PRIVATE(absn2o)
            do ig = 1, ng9
              absn2o = kb_mn2o(indm,ig) + minorfrac(jl,lay) * &
                  (kb_mn2o(indm+1,ig) - kb_mn2o(indm,ig))
              taug(jl,ngs8+ig,lay) = colch4(jl,lay) * &
                  (fac00(jl,lay) * absb(ind0,ig) + &
                  fac10(jl,lay) * absb(ind0+1,ig) + &
                  fac01(jl,lay) * absb(ind1,ig) +  &
                  fac11(jl,lay) * absb(ind1+1,ig)) &
                  + adjcoln2o*absn2o
              fracs(jl,ngs8+ig,lay) = fracrefb(ig)
            enddo
#if defined(_OPENACC) || defined(OMPGPU)
           endif
#endif
          enddo

        enddo
        !$ACC END PARALLEL
        !$OMP END TARGET TEAMS DISTRIBUTE PARALLEL DO

      ENDIF

      !$ACC END DATA
      !$OMP END TARGET DATA

END SUBROUTINE RRTM_TAUMOL9
