!******************************************************************************
SUBROUTINE RRTM_TAUMOL11 (KIDIA,KFDIA,KLEV,taug,&
 & P_TAUAERL,fac00,fac01,fac10,fac11,forfac,forfrac,indfor,jp,jt,jt1,&
 & colh2o,colo2,laytrop,selffac,selffrac,indself,fracs,minorfrac,indminor,scaleminor)  

!     BAND 11:  1480-1800 cm-1 (low - H2O; high - H2O)

!     AUTHOR.
!     -------
!      JJMorcrette, ECMWF

!     MODIFICATIONS.
!     --------------
!      M.Hamrud      01-Oct-2003 CY28 Cleaning
!      NEC           25-Oct-2007 Optimisations
!      JJMorcrette 20110613 flexible number of g-points
!      ABozzo updated to rrtmg v4.85
!     band 11:  1480-1800 cm-1 (low - h2o; low minor - o2)
!                              (high key - h2o; high minor - o2)
! ---------------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK

USE PARRRTM  , ONLY : JPBAND
USE YOERRTM  , ONLY : JPGPT  ,NG11  ,NGS10
USE YOERRTWN , ONLY :      NSPA   ,NSPB
USE YOERRTA11, ONLY : ABSA   ,ABSB   ,FRACREFA, FRACREFB,SELFREF,FORREF, &
                     & KA_MO2, KB_MO2

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
REAL(KIND=JPRB)   ,INTENT(IN)    :: colh2o(KIDIA:KFDIA,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: colo2(KIDIA:KFDIA,KLEV) 
INTEGER(KIND=JPIM),INTENT(IN)    :: laytrop(KIDIA:KFDIA) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: selffac(KIDIA:KFDIA,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: selffrac(KIDIA:KFDIA,KLEV) 
INTEGER(KIND=JPIM),INTENT(IN)    :: indself(KIDIA:KFDIA,KLEV) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: fracs(KIDIA:KFDIA,JPGPT,KLEV) 

INTEGER(KIND=JPIM),INTENT(IN)   :: indfor(KIDIA:KFDIA,KLEV)
REAL(KIND=JPRB)   ,INTENT(IN)   :: forfac(KIDIA:KFDIA,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)   :: forfrac(KIDIA:KFDIA,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)   :: minorfrac(KIDIA:KFDIA,KLEV)
INTEGER(KIND=JPIM),INTENT(IN)   :: indminor(KIDIA:KFDIA,KLEV)
REAL(KIND=JPRB)   ,INTENT(IN)   :: scaleminor(KIDIA:KFDIA,KLEV)
! ---------------------------------------------------------------------------

INTEGER(KIND=JPIM) :: ind0,ind1
INTEGER(KIND=JPIM) :: inds,indf,indm

INTEGER(KIND=JPIM) :: IG, lay
REAL(KIND=JPRB) :: taufor,tauself,scaleO2, tauO2
    !     local integer arrays
    INTEGER(KIND=JPIM) :: laytrop_min, laytrop_max
    integer(KIND=JPIM) :: ixc(KLEV), ixlow(KFDIA,KLEV), ixhigh(KFDIA,KLEV)
    INTEGER(KIND=JPIM) :: ich, icl, ixc0, ixp, jc, jl

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

! Minor gas mapping level :
!     lower - o2, p = 706.2720 mbar, t = 278.94 k
!     upper - o2, p = 4.758820 mbarm t = 250.85 k

!     Compute the optical depth by interpolating in ln(pressure) and 
!     temperature.  Below LAYTROP, the water vapor self-continuum and foreign continuum 
!     is interpolated (in temperature) separately.
  
      ! Lower atmosphere loop
      do lay = 1, laytrop_min
        do jl = KIDIA, KFDIA

          ind0 = ((jp(jl,lay)-1)*5+(jt(jl,lay)-1))*nspa(11) + 1
          ind1 = (jp(jl,lay)*5+(jt1(jl,lay)-1))*nspa(11) + 1
          inds = indself(jl,lay)
          indf = indfor(jl,lay)
          indm = indminor(jl,lay)
          scaleo2 = colo2(jl,lay)*scaleminor(jl,lay)
!$NEC unroll(NG11)
          do ig = 1, ng11
            tauself = selffac(jl,lay) * (selfref(inds,ig) + selffrac(jl,lay) * &
                 (selfref(inds+1,ig) - selfref(inds,ig)))
            taufor = forfac(jl,lay) * (forref(indf,ig) + forfrac(jl,lay) * &
                 (forref(indf+1,ig) - forref(indf,ig)))
            tauo2 =  scaleo2 * (ka_mo2(indm,ig) + minorfrac(jl,lay) * &
                 (ka_mo2(indm+1,ig) - ka_mo2(indm,ig)))
            taug(jl,ngs10+ig,lay) = colh2o(jl,lay) * &
                 (fac00(jl,lay) * absa(ind0,ig) + &
                 fac10(jl,lay) * absa(ind0+1,ig) + &
                 fac01(jl,lay) * absa(ind1,ig) + &
                 fac11(jl,lay) * absa(ind1+1,ig)) &
                 + tauself + taufor &
                 + tauo2
            fracs(jl,ngs10+ig,lay) = fracrefa(ig)
          enddo
        enddo

      enddo

      ! Upper atmosphere loop
      do lay = laytrop_max+1, KLEV
        do jl = KIDIA, KFDIA

          ind0 = ((jp(jl,lay)-13)*5+(jt(jl,lay)-1))*nspb(11) + 1
          ind1 = ((jp(jl,lay)-12)*5+(jt1(jl,lay)-1))*nspb(11) + 1
          indf = indfor(jl,lay)
          indm = indminor(jl,lay)
          scaleo2 = colo2(jl,lay)*scaleminor(jl,lay)
!$NEC unroll(NG11)
          do ig = 1, ng11
            taufor = forfac(jl,lay) * (forref(indf,ig) + forfrac(jl,lay) * &
                 (forref(indf+1,ig) - forref(indf,ig)))
            tauo2 =  scaleo2 * (kb_mo2(indm,ig) + minorfrac(jl,lay) * &
                 (kb_mo2(indm+1,ig) - kb_mo2(indm,ig)))
            taug(jl,ngs10+ig,lay) = colh2o(jl,lay) * &
                 (fac00(jl,lay) * absb(ind0,ig) + &
                 fac10(jl,lay) * absb(ind0+1,ig) + &
                 fac01(jl,lay) * absb(ind1,ig) + &
                 fac11(jl,lay) * absb(ind1+1,ig))  &
                 + taufor &
                 + tauo2
            fracs(jl,ngs10+ig,lay) = fracrefb(ig)
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

          ind0 = ((jp(jl,lay)-1)*5+(jt(jl,lay)-1))*nspa(11) + 1
          ind1 = (jp(jl,lay)*5+(jt1(jl,lay)-1))*nspa(11) + 1
          inds = indself(jl,lay)
          indf = indfor(jl,lay)
          indm = indminor(jl,lay)
          scaleo2 = colo2(jl,lay)*scaleminor(jl,lay)
!$NEC unroll(NG11)
          do ig = 1, ng11
            tauself = selffac(jl,lay) * (selfref(inds,ig) + selffrac(jl,lay) * &
                 (selfref(inds+1,ig) - selfref(inds,ig)))
            taufor = forfac(jl,lay) * (forref(indf,ig) + forfrac(jl,lay) * &
                 (forref(indf+1,ig) - forref(indf,ig)))
            tauo2 =  scaleo2 * (ka_mo2(indm,ig) + minorfrac(jl,lay) * &
                 (ka_mo2(indm+1,ig) - ka_mo2(indm,ig)))
            taug(jl,ngs10+ig,lay) = colh2o(jl,lay) * &
                 (fac00(jl,lay) * absa(ind0,ig) + &
                 fac10(jl,lay) * absa(ind0+1,ig) + &
                 fac01(jl,lay) * absa(ind1,ig) + &
                 fac11(jl,lay) * absa(ind1+1,ig)) &
                 + tauself + taufor &
                 + tauo2
            fracs(jl,ngs10+ig,lay) = fracrefa(ig)
          enddo
        enddo

        ! Upper atmosphere part
        ixc0 = KFDIA - KIDIA + 1 - ixc0
!$NEC ivdep
        do ixp = 1, ixc0
          jl = ixhigh(ixp,lay)

          ind0 = ((jp(jl,lay)-13)*5+(jt(jl,lay)-1))*nspb(11) + 1
          ind1 = ((jp(jl,lay)-12)*5+(jt1(jl,lay)-1))*nspb(11) + 1
          indf = indfor(jl,lay)
          indm = indminor(jl,lay)
          scaleo2 = colo2(jl,lay)*scaleminor(jl,lay)
!$NEC unroll(NG11)
          do ig = 1, ng11
            taufor = forfac(jl,lay) * (forref(indf,ig) + forfrac(jl,lay) * &
                 (forref(indf+1,ig) - forref(indf,ig)))
            tauo2 =  scaleo2 * (kb_mo2(indm,ig) + minorfrac(jl,lay) * &
                 (kb_mo2(indm+1,ig) - kb_mo2(indm,ig)))
            taug(jl,ngs10+ig,lay) = colh2o(jl,lay) * &
                 (fac00(jl,lay) * absb(ind0,ig) + &
                 fac10(jl,lay) * absb(ind0+1,ig) + &
                 fac01(jl,lay) * absb(ind1,ig) + &
                 fac11(jl,lay) * absb(ind1+1,ig))  &
                 + taufor &
                 + tauo2
            fracs(jl,ngs10+ig,lay) = fracrefb(ig)
          enddo
        enddo

      enddo

END SUBROUTINE RRTM_TAUMOL11
