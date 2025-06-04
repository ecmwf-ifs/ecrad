!*******************************************************************************
SUBROUTINE RRTM_TAUMOL8 (KIDIA,KFDIA,KLEV,taug,wx,&
 & P_TAUAERL,fac00,fac01,fac10,fac11,forfac,forfrac,indfor,jp,jt,jt1,&
 & colh2o,colo3,coln2o,colco2,coldry,laytrop,selffac,selffrac,indself,fracs, &
 & minorfrac,indminor)

!     BAND 8:  1080-1180 cm-1 (low (i.e.>~300mb) - H2O; high - O3)

!     AUTHOR.
!     -------
!      JJMorcrette, ECMWF

!     MODIFICATIONS.
!     --------------
!      M.Hamrud      01-Oct-2003 CY28 Cleaning
!      NEC           25-Oct-2007 Optimisations
!      JJMorcrette 20110613 flexible number of g-points
!      ABozzo 201306 updated to rrtmg v4.85
!     band 8:  1080-1180 cm-1 (low key - h2o; low minor - co2,o3,n2o)
!                             (high key - o3; high minor - co2, n2o)
! ---------------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM     ,JPRB

USE PARRRTM  , ONLY : JPBAND ,JPXSEC
USE YOERRTM  , ONLY : JPGPT  ,NG8   ,NGS7
USE YOERRTWN , ONLY : NSPA   ,NSPB
USE YOERRTA8 , ONLY : ABSA   ,ABSB   ,FRACREFA, FRACREFB,SELFREF,KA_MCO2 ,KB_MCO2  ,&
 & KA_MN2O , KB_MN2O,KA_MO3,CFC12  ,CFC22ADJ,FORREF
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
REAL(KIND=JPRB)   ,INTENT(IN)    :: colh2o(KIDIA:KFDIA,KLEV)
REAL(KIND=JPRB)   ,INTENT(IN)    :: colo3(KIDIA:KFDIA,KLEV)
REAL(KIND=JPRB)   ,INTENT(IN)    :: coln2o(KIDIA:KFDIA,KLEV)
REAL(KIND=JPRB)   ,INTENT(IN)    :: colco2(KIDIA:KFDIA,KLEV)
REAL(KIND=JPRB)   ,INTENT(IN)    :: coldry(KIDIA:KFDIA,KLEV)
INTEGER(KIND=JPIM),INTENT(IN)    :: laytrop(KIDIA:KFDIA)
REAL(KIND=JPRB)   ,INTENT(IN)    :: selffac(KIDIA:KFDIA,KLEV)
REAL(KIND=JPRB)   ,INTENT(IN)    :: selffrac(KIDIA:KFDIA,KLEV)
INTEGER(KIND=JPIM),INTENT(IN)    :: indself(KIDIA:KFDIA,KLEV)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: fracs(KIDIA:KFDIA,JPGPT,KLEV)

INTEGER(KIND=JPIM),INTENT(IN)   :: indfor(KIDIA:KFDIA,KLEV)
REAL(KIND=JPRB)   ,INTENT(IN)   :: forfrac(KIDIA:KFDIA,KLEV)
REAL(KIND=JPRB)   ,INTENT(IN)   :: forfac(KIDIA:KFDIA,KLEV)
REAL(KIND=JPRB)   ,INTENT(IN)   :: minorfrac(KIDIA:KFDIA,KLEV)
INTEGER(KIND=JPIM),INTENT(IN)   :: indminor(KIDIA:KFDIA,KLEV)

! ---------------------------------------------------------------------------

INTEGER(KIND=JPIM) :: ind0,ind1,inds,indf,indm

INTEGER(KIND=JPIM) :: IG, lay

REAL(KIND=JPRB) :: chi_co2, ratco2, adjfac, adjcolco2
REAL(KIND=JPRB) :: taufor,tauself, abso3, absco2, absn2o
    !     local integer arrays
    INTEGER(KIND=JPIM) :: laytrop_min, laytrop_max
    integer(KIND=JPIM) :: ixc(KLEV), ixlow(KFDIA,KLEV), ixhigh(KFDIA,KLEV)
    INTEGER(KIND=JPIM) :: ich, icl, ixc0, ixp, jc, jl

    !$ACC DATA PRESENT(taug, wx, P_TAUAERL, fac00, fac01, fac10, fac11, jp, jt, &
    !$ACC             jt1, colh2o, colo3, coln2o, colco2, coldry, laytrop, &
    !$ACC             selffac, selffrac, indself, fracs, indfor, forfrac, forfac, &
    !$ACC             minorfrac, indminor)

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

! Minor gas mapping level:
!     lower - co2, p = 1053.63 mb, t = 294.2 k
!     lower - o3,  p = 317.348 mb, t = 240.77 k
!     lower - n2o, p = 706.2720 mb, t= 278.94 k
!     lower - cfc12,cfc11
!     upper - co2, p = 35.1632 mb, t = 223.28 k
!     upper - n2o, p = 8.716e-2 mb, t = 226.03 k

! Compute the optical depth by interpolating in ln(pressure) and
! temperature, and appropriate species.  Below laytrop, the water vapor
! self-continuum and foreign continuum is interpolated (in temperature)
! separately.

      ! Lower atmosphere loop
      !$ACC WAIT
      !$ACC PARALLEL DEFAULT(NONE) ASYNC(1)
      !$ACC LOOP GANG VECTOR COLLAPSE(2) PRIVATE(ind0, ind1, inds, indf, indm, chi_co2, ratco2, adjfac, adjcolco2)
      do lay = 1, laytrop_min
        do jl = KIDIA, KFDIA

          !  In atmospheres where the amount of CO2 is too great to be considered
          !  a minor species, adjust the column amount of CO2 by an empirical factor
          !  to obtain the proper contribution.
          chi_co2 = colco2(jl,lay)/(coldry(jl,lay))
          ratco2 = 1.e20_JPRB*chi_co2/chi_mls(2,jp(jl,lay)+1)
          if (ratco2 .gt. 3.0_JPRB) then
            adjfac = 2.0_JPRB+(ratco2-2.0_JPRB)**0.65_JPRB
            adjcolco2 = adjfac*chi_mls(2,jp(jl,lay)+1)*coldry(jl,lay)*1.e-20_JPRB
          else
            adjcolco2 = colco2(jl,lay)
          endif

          ind0 = ((jp(jl,lay)-1)*5+(jt(jl,lay)-1))*nspa(8) + 1
          ind1 = (jp(jl,lay)*5+(jt1(jl,lay)-1))*nspa(8) + 1
          inds = indself(jl,lay)
          indf = indfor(jl,lay)
          indm = indminor(jl,lay)
          !$ACC LOOP SEQ PRIVATE(taufor, tauself, abso3, absco2, absn2o)
!$NEC unroll(NG8)
          do ig = 1, ng8
            tauself = selffac(jl,lay) * (selfref(inds,ig) + selffrac(jl,lay) * &
                 (selfref(inds+1,ig) - selfref(inds,ig)))
            taufor = forfac(jl,lay) * (forref(indf,ig) + forfrac(jl,lay) * &
                 (forref(indf+1,ig) - forref(indf,ig)))
            absco2 =  (ka_mco2(indm,ig) + minorfrac(jl,lay) * &
                 (ka_mco2(indm+1,ig) - ka_mco2(indm,ig)))
            abso3 =  (ka_mo3(indm,ig) + minorfrac(jl,lay) * &
                 (ka_mo3(indm+1,ig) - ka_mo3(indm,ig)))
            absn2o =  (ka_mn2o(indm,ig) + minorfrac(jl,lay) * &
                 (ka_mn2o(indm+1,ig) - ka_mn2o(indm,ig)))
            taug(jl,ngs7+ig,lay) = colh2o(jl,lay) * &
                 (fac00(jl,lay) * absa(ind0,ig) + &
                 fac10(jl,lay) * absa(ind0+1,ig) + &
                 fac01(jl,lay) * absa(ind1,ig) +  &
                 fac11(jl,lay) * absa(ind1+1,ig)) &
                 + tauself + taufor &
                 + adjcolco2*absco2 &
                 + colo3(jl,lay) * abso3 &
                 + coln2o(jl,lay) * absn2o &
                 + wx(jl,3,lay) * cfc12(ig) &
                 + wx(jl,4,lay) * cfc22adj(ig)
            fracs(jl,ngs7+ig,lay) = fracrefa(ig)
          enddo
        enddo

      enddo
      !$ACC END PARALLEL

      ! Upper atmosphere loop
      !$ACC PARALLEL DEFAULT(NONE) ASYNC(1)
      !$ACC LOOP GANG VECTOR COLLAPSE(2) PRIVATE(ind0, ind1, indm, chi_co2, ratco2, adjfac, adjcolco2)
      do lay = laytrop_max+1, KLEV
        do jl = KIDIA, KFDIA

          !  In atmospheres where the amount of CO2 is too great to be considered
          !  a minor species, adjust the column amount of CO2 by an empirical factor
          !  to obtain the proper contribution.
          chi_co2 = colco2(jl,lay)/coldry(jl,lay)
          ratco2 = 1.e20_JPRB*chi_co2/chi_mls(2,jp(jl,lay)+1)
          if (ratco2 .gt. 3.0_JPRB) then
            adjfac = 2.0_JPRB+(ratco2-2.0_JPRB)**0.65_JPRB
            adjcolco2 = adjfac*chi_mls(2,jp(jl,lay)+1) * coldry(jl,lay)*1.e-20_JPRB
          else
            adjcolco2 = colco2(jl,lay)
          endif

          ind0 = ((jp(jl,lay)-13)*5+(jt(jl,lay)-1))*nspb(8) + 1
          ind1 = ((jp(jl,lay)-12)*5+(jt1(jl,lay)-1))*nspb(8) + 1
          indm = indminor(jl,lay)
          !$ACC LOOP SEQ PRIVATE(absco2, absn2o)
!$NEC unroll(NG8)
          do ig = 1, ng8
            absco2 =  (kb_mco2(indm,ig) + minorfrac(jl,lay) * &
                 (kb_mco2(indm+1,ig) - kb_mco2(indm,ig)))
            absn2o =  (kb_mn2o(indm,ig) + minorfrac(jl,lay) * &
                 (kb_mn2o(indm+1,ig) - kb_mn2o(indm,ig)))
            taug(jl,ngs7+ig,lay) = colo3(jl,lay) * &
                 (fac00(jl,lay) * absb(ind0,ig) + &
                 fac10(jl,lay) * absb(ind0+1,ig) + &
                 fac01(jl,lay) * absb(ind1,ig) + &
                 fac11(jl,lay) * absb(ind1+1,ig)) &
                 + adjcolco2*absco2 &
                 + coln2o(jl,lay)*absn2o &
                 + wx(jl,3,lay) * cfc12(ig) &
                 + wx(jl,4,lay) * cfc22adj(ig)
            fracs(jl,ngs7+ig,lay) = fracrefb(ig)
          enddo
        enddo

      enddo
      !$ACC END PARALLEL

      IF (laytrop_max /= laytrop_min) THEN
        ! Mixed loop
        ! Lower atmosphere part
        !$ACC PARALLEL DEFAULT(NONE) ASYNC(1)
        !$ACC LOOP GANG VECTOR COLLAPSE(2) PRIVATE(chi_co2, ratco2, adjfac, adjcolco2, ind0, ind1, inds, indf, indm)
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

            !  In atmospheres where the amount of CO2 is too great to be considered
            !  a minor species, adjust the column amount of CO2 by an empirical factor
            !  to obtain the proper contribution.
            chi_co2 = colco2(jl,lay)/(coldry(jl,lay))
            ratco2 = 1.e20_JPRB*chi_co2/chi_mls(2,jp(jl,lay)+1)
            if (ratco2 .gt. 3.0_JPRB) then
              adjfac = 2.0_JPRB+(ratco2-2.0_JPRB)**0.65_JPRB
              adjcolco2 = adjfac*chi_mls(2,jp(jl,lay)+1)*coldry(jl,lay)*1.e-20_JPRB
            else
              adjcolco2 = colco2(jl,lay)
            endif

            ind0 = ((jp(jl,lay)-1)*5+(jt(jl,lay)-1))*nspa(8) + 1
            ind1 = (jp(jl,lay)*5+(jt1(jl,lay)-1))*nspa(8) + 1
            inds = indself(jl,lay)
            indf = indfor(jl,lay)
            indm = indminor(jl,lay)
!$NEC unroll(NG8)
            !$ACC LOOP SEQ PRIVATE(tauself, taufor, absco2, abso3, absn2o)
            do ig = 1, ng8
              tauself = selffac(jl,lay) * (selfref(inds,ig) + selffrac(jl,lay) * &
                  (selfref(inds+1,ig) - selfref(inds,ig)))
              taufor = forfac(jl,lay) * (forref(indf,ig) + forfrac(jl,lay) * &
                  (forref(indf+1,ig) - forref(indf,ig)))
              absco2 =  (ka_mco2(indm,ig) + minorfrac(jl,lay) * &
                  (ka_mco2(indm+1,ig) - ka_mco2(indm,ig)))
              abso3 =  (ka_mo3(indm,ig) + minorfrac(jl,lay) * &
                  (ka_mo3(indm+1,ig) - ka_mo3(indm,ig)))
              absn2o =  (ka_mn2o(indm,ig) + minorfrac(jl,lay) * &
                  (ka_mn2o(indm+1,ig) - ka_mn2o(indm,ig)))
              taug(jl,ngs7+ig,lay) = colh2o(jl,lay) * &
                  (fac00(jl,lay) * absa(ind0,ig) + &
                  fac10(jl,lay) * absa(ind0+1,ig) + &
                  fac01(jl,lay) * absa(ind1,ig) +  &
                  fac11(jl,lay) * absa(ind1+1,ig)) &
                  + tauself + taufor &
                  + adjcolco2*absco2 &
                  + colo3(jl,lay) * abso3 &
                  + coln2o(jl,lay) * absn2o &
                  + wx(jl,3,lay) * cfc12(ig) &
                  + wx(jl,4,lay) * cfc22adj(ig)
              fracs(jl,ngs7+ig,lay) = fracrefa(ig)
            enddo
#ifdef _OPENACC
         else
#else
          enddo

          ! Upper atmosphere loop
          ixc0 = KFDIA - KIDIA + 1 - ixc0
!$NEC ivdep
          do ixp = 1, ixc0
            jl = ixhigh(ixp,lay)
#endif

            !  In atmospheres where the amount of CO2 is too great to be considered
            !  a minor species, adjust the column amount of CO2 by an empirical factor
            !  to obtain the proper contribution.
            chi_co2 = colco2(jl,lay)/coldry(jl,lay)
            ratco2 = 1.e20_JPRB*chi_co2/chi_mls(2,jp(jl,lay)+1)
            if (ratco2 .gt. 3.0_JPRB) then
              adjfac = 2.0_JPRB+(ratco2-2.0_JPRB)**0.65_JPRB
              adjcolco2 = adjfac*chi_mls(2,jp(jl,lay)+1) * coldry(jl,lay)*1.e-20_JPRB
            else
              adjcolco2 = colco2(jl,lay)
            endif

            ind0 = ((jp(jl,lay)-13)*5+(jt(jl,lay)-1))*nspb(8) + 1
            ind1 = ((jp(jl,lay)-12)*5+(jt1(jl,lay)-1))*nspb(8) + 1
            indm = indminor(jl,lay)
!$NEC unroll(NG8)
            !$ACC LOOP SEQ PRIVATE(absco2, absn2o)
            do ig = 1, ng8
              absco2 =  (kb_mco2(indm,ig) + minorfrac(jl,lay) * &
                  (kb_mco2(indm+1,ig) - kb_mco2(indm,ig)))
              absn2o =  (kb_mn2o(indm,ig) + minorfrac(jl,lay) * &
                  (kb_mn2o(indm+1,ig) - kb_mn2o(indm,ig)))
              taug(jl,ngs7+ig,lay) = colo3(jl,lay) * &
                  (fac00(jl,lay) * absb(ind0,ig) + &
                  fac10(jl,lay) * absb(ind0+1,ig) + &
                  fac01(jl,lay) * absb(ind1,ig) + &
                  fac11(jl,lay) * absb(ind1+1,ig)) &
                  + adjcolco2*absco2 &
                  + coln2o(jl,lay)*absn2o &
                  + wx(jl,3,lay) * cfc12(ig) &
                  + wx(jl,4,lay) * cfc22adj(ig)
              fracs(jl,ngs7+ig,lay) = fracrefb(ig)
            enddo
#ifdef _OPENACC
           endif
#endif
          enddo

        enddo
        !$ACC END PARALLEL

      ENDIF

      !$ACC END DATA

END SUBROUTINE RRTM_TAUMOL8
