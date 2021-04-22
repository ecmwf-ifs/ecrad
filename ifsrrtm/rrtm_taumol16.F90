!----------------------------------------------------------------------------
SUBROUTINE RRTM_TAUMOL16 (KIDIA,KFDIA,KLEV,taug,&
 & P_TAUAERL,fac00,fac01,fac10,fac11,forfac,forfrac,indfor,jp,jt,jt1,oneminus,&
 & colh2o,colch4,laytrop,selffac,selffrac,indself,fracs, &
 & rat_h2och4,rat_h2och4_1)  

!     BAND 16:  2600-3000 cm-1 (low - H2O,CH4; high - nothing)

!     AUTHOR.
!     -------
!      JJMorcrette, ECMWF

!     MODIFICATIONS.
!     --------------
!      M.Hamrud      01-Oct-2003 CY28 Cleaning
!      NEC           25-Oct-2007 Optimisations
!      JJMorcrette 20110613 flexible number of g-points
!      ABozzo 201306 updated to rrtmg v4.85
!      band 16:  2600-3250 cm-1 (low key- h2o,ch4; high key - ch4)
! ---------------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK

USE PARRRTM  , ONLY : JPBAND
USE YOERRTM  , ONLY : JPGPT  ,NGS15  ,NG16
USE YOERRTWN , ONLY : NSPA,NSPB   
USE YOERRTA16, ONLY : ABSA,ABSB,FRACREFA,FRACREFB,SELFREF,FORREF 
USE YOERRTRF, ONLY : CHI_MLS

IMPLICIT NONE

INTEGER(KIND=JPIM),INTENT(IN)    :: KIDIA
INTEGER(KIND=JPIM),INTENT(IN)    :: KFDIA
INTEGER(KIND=JPIM),INTENT(IN)    :: KLEV 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: taug(KIDIA:KFDIA,JPGPT,KLEV) 
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
REAL(KIND=JPRB)   ,INTENT(IN)    :: colch4(KIDIA:KFDIA,KLEV) 
INTEGER(KIND=JPIM),INTENT(IN)    :: laytrop(KIDIA:KFDIA) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: selffac(KIDIA:KFDIA,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: selffrac(KIDIA:KFDIA,KLEV) 
INTEGER(KIND=JPIM),INTENT(IN)    :: indself(KIDIA:KFDIA,KLEV) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: fracs(KIDIA:KFDIA,JPGPT,KLEV) 

REAL(KIND=JPRB)   ,INTENT(IN)   :: rat_h2och4(KIDIA:KFDIA,KLEV)
REAL(KIND=JPRB)   ,INTENT(IN)   :: rat_h2och4_1(KIDIA:KFDIA,KLEV)
INTEGER(KIND=JPIM),INTENT(IN)   :: indfor(KIDIA:KFDIA,KLEV)
REAL(KIND=JPRB)   ,INTENT(IN)   :: forfac(KIDIA:KFDIA,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)   :: forfrac(KIDIA:KFDIA,KLEV) 

! ---------------------------------------------------------------------------

INTEGER(KIND=JPIM) :: IG, IND0, IND1, INDS,INDF, JS,JS1,JPL,lay

REAL(KIND=JPRB) ::  fac000, fac100, fac200,&
 & fac010, fac110, fac210, &
 & fac001, fac101, fac201, &
 & fac011, fac111, fac211
REAL(KIND=JPRB) :: p, p4, fk0, fk1, fk2

REAL(KIND=JPRB) :: refrat_planck_a
REAL(KIND=JPRB) :: taufor,tauself,tau_major(ng16),tau_major1(ng16)
REAL(KIND=JPRB) :: fs, specmult, specparm,speccomb,  &
& fs1, specmult1, specparm1,speccomb1, &
& fpl, specmult_planck, specparm_planck,speccomb_planck
    !     local integer arrays
    INTEGER(KIND=JPIM) :: laytrop_min, laytrop_max
    integer(KIND=JPIM) :: ixc(KLEV), ixlow(KFDIA,KLEV), ixhigh(KFDIA,KLEV)
    INTEGER(KIND=JPIM) :: ich, icl, ixc0, ixp, jc, jl

#define MOD1(x) ((x) - AINT((x)))

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


      ! P = 387. mb (Level 6)
      refrat_planck_a = chi_mls(1,6)/chi_mls(6,6)

      ! Compute the optical depth by interpolating in ln(pressure),
      ! temperature,and appropriate species.  Below laytrop, the water
      ! vapor self-continuum and foreign continuum is interpolated
      ! (in temperature) separately.

      ! Lower atmosphere loop
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

          speccomb_planck = colh2o(jl,lay)+refrat_planck_a*colch4(jl,lay)
          specparm_planck = MIN(colh2o(jl,lay)/speccomb_planck,oneminus)
          specmult_planck = 8._JPRB*specparm_planck
          jpl = 1 + int(specmult_planck)
          fpl = MOD1(specmult_planck)

          ind0 = ((jp(jl,lay)-1)*5+(jt(jl,lay)-1))*nspa(16) + js
          ind1 = (jp(jl,lay)*5+(jt1(jl,lay)-1))*nspa(16) + js1
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
!$NEC unroll(NG16)
            tau_major(1:ng16) = speccomb *    &
             (fac000 * absa(ind0,1:ng16)    + &
              fac100 * absa(ind0+1,1:ng16)  + &
              fac200 * absa(ind0+2,1:ng16)  + &
              fac010 * absa(ind0+9,1:ng16)  + &
              fac110 * absa(ind0+10,1:ng16) + &
              fac210 * absa(ind0+11,1:ng16))
          else if (specparm .gt. 0.875_JPRB) then
!$NEC unroll(NG16)
            tau_major(1:ng16) = speccomb *   &
             (fac200 * absa(ind0-1,1:ng16) + &
              fac100 * absa(ind0,1:ng16)   + &
              fac000 * absa(ind0+1,1:ng16) + &
              fac210 * absa(ind0+8,1:ng16) + &
              fac110 * absa(ind0+9,1:ng16) + &
              fac010 * absa(ind0+10,1:ng16))
          else
!$NEC unroll(NG16)
             tau_major(1:ng16) = speccomb *   &
              (fac000 * absa(ind0,1:ng16)   + &
               fac100 * absa(ind0+1,1:ng16) + &
               fac010 * absa(ind0+9,1:ng16) + &
               fac110 * absa(ind0+10,1:ng16))
          endif

          if (specparm1 .lt. 0.125_JPRB) then
!$NEC unroll(NG16)
             tau_major1(1:ng16) = speccomb1 *  &
              (fac001 * absa(ind1,1:ng16)    + &
               fac101 * absa(ind1+1,1:ng16)  + &
               fac201 * absa(ind1+2,1:ng16)  + &
               fac011 * absa(ind1+9,1:ng16)  + &
               fac111 * absa(ind1+10,1:ng16) + &
               fac211 * absa(ind1+11,1:ng16))
          else if (specparm1 .gt. 0.875_JPRB) then
!$NEC unroll(NG16)
            tau_major1(1:ng16) = speccomb1 * &
             (fac201 * absa(ind1-1,1:ng16) + &
              fac101 * absa(ind1,1:ng16)   + &
              fac001 * absa(ind1+1,1:ng16) + &
              fac211 * absa(ind1+8,1:ng16) + &
              fac111 * absa(ind1+9,1:ng16) + &
              fac011 * absa(ind1+10,1:ng16))
          else
!$NEC unroll(NG16)
            tau_major1(1:ng16) = speccomb1 * &
             (fac001 * absa(ind1,1:ng16)   + &
              fac101 * absa(ind1+1,1:ng16) + &
              fac011 * absa(ind1+9,1:ng16) + &
              fac111 * absa(ind1+10,1:ng16))
          endif

!$NEC unroll(NG16)
          do ig = 1, ng16
            tauself = selffac(jl,lay)* (selfref(inds,ig) + selffrac(jl,lay) * &
                 (selfref(inds+1,ig) - selfref(inds,ig)))
            taufor =  forfac(jl,lay) * (forref(indf,ig) + forfrac(jl,lay) * &
                 (forref(indf+1,ig) - forref(indf,ig)))

            taug(jl,ngs15+ig,lay) = tau_major(ig) + tau_major1(ig) &
                 + tauself + taufor
            fracs(jl,ngs15+ig,lay) = fracrefa(ig,jpl) + fpl * &
                 (fracrefa(ig,jpl+1)-fracrefa(ig,jpl))
          enddo
        enddo

      enddo

      ! Upper atmosphere loop
      do lay = laytrop_max+1, KLEV
        do jl = KIDIA, KFDIA

          ind0 = ((jp(jl,lay)-13)*5+(jt(jl,lay)-1))*nspb(16) + 1
          ind1 = ((jp(jl,lay)-12)*5+(jt1(jl,lay)-1))*nspb(16) + 1
!$NEC unroll(NG16)
          do ig = 1, ng16
            taug(jl,ngs15+ig,lay) = colch4(jl,lay) * &
                 (fac00(jl,lay) * absb(ind0,ig) + &
                 fac10(jl,lay) * absb(ind0+1,ig) + &
                 fac01(jl,lay) * absb(ind1,ig) + &
                 fac11(jl,lay) * absb(ind1+1,ig))
            fracs(jl,ngs15+ig,lay) = fracrefb(ig)
          enddo
        enddo

      enddo

      IF (laytrop_max == laytrop_min) RETURN
      ! Mixed loop
      ! Lower atmosphere part
      do lay = laytrop_min+1, laytrop_max
        ixc0 = ixc(lay)

!$NEC ivdep
        do ixp = 1, ixc0
          jl = ixlow(ixp,lay)

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

          speccomb_planck = colh2o(jl,lay)+refrat_planck_a*colch4(jl,lay)
          specparm_planck = MIN(colh2o(jl,lay)/speccomb_planck,oneminus)
          specmult_planck = 8._JPRB*specparm_planck
          jpl = 1 + int(specmult_planck)
          fpl = MOD1(specmult_planck)

          ind0 = ((jp(jl,lay)-1)*5+(jt(jl,lay)-1))*nspa(16) + js
          ind1 = (jp(jl,lay)*5+(jt1(jl,lay)-1))*nspa(16) + js1
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
!$NEC unroll(NG16)
            tau_major(1:ng16) = speccomb *    &
             (fac000 * absa(ind0,1:ng16)    + &
              fac100 * absa(ind0+1,1:ng16)  + &
              fac200 * absa(ind0+2,1:ng16)  + &
              fac010 * absa(ind0+9,1:ng16)  + &
              fac110 * absa(ind0+10,1:ng16) + &
              fac210 * absa(ind0+11,1:ng16))
          else if (specparm .gt. 0.875_JPRB) then
!$NEC unroll(NG16)
            tau_major(1:ng16) = speccomb *   &
             (fac200 * absa(ind0-1,1:ng16) + &
              fac100 * absa(ind0,1:ng16)   + &
              fac000 * absa(ind0+1,1:ng16) + &
              fac210 * absa(ind0+8,1:ng16) + &
              fac110 * absa(ind0+9,1:ng16) + &
              fac010 * absa(ind0+10,1:ng16))
          else
!$NEC unroll(NG16)
             tau_major(1:ng16) = speccomb *   &
              (fac000 * absa(ind0,1:ng16)   + &
               fac100 * absa(ind0+1,1:ng16) + &
               fac010 * absa(ind0+9,1:ng16) + &
               fac110 * absa(ind0+10,1:ng16))
          endif

          if (specparm1 .lt. 0.125_JPRB) then
!$NEC unroll(NG16)
             tau_major1(1:ng16) = speccomb1 *  &
              (fac001 * absa(ind1,1:ng16)    + &
               fac101 * absa(ind1+1,1:ng16)  + &
               fac201 * absa(ind1+2,1:ng16)  + &
               fac011 * absa(ind1+9,1:ng16)  + &
               fac111 * absa(ind1+10,1:ng16) + &
               fac211 * absa(ind1+11,1:ng16))
          else if (specparm1 .gt. 0.875_JPRB) then
!$NEC unroll(NG16)
            tau_major1(1:ng16) = speccomb1 * &
             (fac201 * absa(ind1-1,1:ng16) + &
              fac101 * absa(ind1,1:ng16)   + &
              fac001 * absa(ind1+1,1:ng16) + &
              fac211 * absa(ind1+8,1:ng16) + &
              fac111 * absa(ind1+9,1:ng16) + &
              fac011 * absa(ind1+10,1:ng16))
          else
!$NEC unroll(NG16)
            tau_major1(1:ng16) = speccomb1 * &
             (fac001 * absa(ind1,1:ng16)   + &
              fac101 * absa(ind1+1,1:ng16) + &
              fac011 * absa(ind1+9,1:ng16) + &
              fac111 * absa(ind1+10,1:ng16))
          endif

!$NEC unroll(NG16)
          do ig = 1, ng16
            tauself = selffac(jl,lay)* (selfref(inds,ig) + selffrac(jl,lay) * &
                 (selfref(inds+1,ig) - selfref(inds,ig)))
            taufor =  forfac(jl,lay) * (forref(indf,ig) + forfrac(jl,lay) * &
                 (forref(indf+1,ig) - forref(indf,ig)))

            taug(jl,ngs15+ig,lay) = tau_major(ig) + tau_major1(ig) &
                 + tauself + taufor
            fracs(jl,ngs15+ig,lay) = fracrefa(ig,jpl) + fpl * &
                 (fracrefa(ig,jpl+1)-fracrefa(ig,jpl))
          enddo
        enddo

        ! Upper atmosphere part
        ixc0 = KFDIA - KIDIA + 1 - ixc0
!$NEC ivdep
        do ixp = 1, ixc0
          jl = ixhigh(ixp,lay)

          ind0 = ((jp(jl,lay)-13)*5+(jt(jl,lay)-1))*nspb(16) + 1
          ind1 = ((jp(jl,lay)-12)*5+(jt1(jl,lay)-1))*nspb(16) + 1
!$NEC unroll(NG16)
          do ig = 1, ng16
            taug(jl,ngs15+ig,lay) = colch4(jl,lay) * &
                 (fac00(jl,lay) * absb(ind0,ig) + &
                 fac10(jl,lay) * absb(ind0+1,ig) + &
                 fac01(jl,lay) * absb(ind1,ig) + &
                 fac11(jl,lay) * absb(ind1+1,ig))
            fracs(jl,ngs15+ig,lay) = fracrefb(ig)
          enddo
        enddo

      enddo


END SUBROUTINE RRTM_TAUMOL16
