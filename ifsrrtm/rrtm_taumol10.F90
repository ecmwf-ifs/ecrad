!*******************************************************************************
SUBROUTINE RRTM_TAUMOL10 (KIDIA,KFDIA,KLEV,taug,&
 & P_TAUAERL,fac00,fac01,fac10,fac11,forfac,forfrac,indfor,jp,jt,jt1,&
 & colh2o,laytrop,selffac,selffrac,indself,fracs)  

!     BAND 10:  1390-1480 cm-1 (low - H2O; high - H2O)

!     AUTHOR.
!     -------
!     JJMorcrette, ECMWF

!     MODIFICATIONS.
!     --------------
!      M.Hamrud      01-Oct-2003 CY28 Cleaning
!      NEC           25-Oct-2007 Optimisations
!      JJMorcrette 20110613 flexible number of g-points
!      ABozzo 201306 updated to rrtmg v4.85
! ---------------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK

USE PARRRTM  , ONLY : JPBAND
USE YOERRTM  , ONLY : JPGPT  ,NG10   ,NGS9
USE YOERRTWN , ONLY : NSPA   ,NSPB
USE YOERRTA10, ONLY : ABSA   ,ABSB   ,FRACREFA, FRACREFB, FORREF   ,SELFREF

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
INTEGER(KIND=JPIM),INTENT(IN)    :: laytrop(KIDIA:KFDIA) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: fracs(KIDIA:KFDIA,JPGPT,KLEV) 

REAL(KIND=JPRB)   ,INTENT(IN)    :: selffac(KIDIA:KFDIA,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: selffrac(KIDIA:KFDIA,KLEV) 
INTEGER(KIND=JPIM),INTENT(IN)    :: indself(KIDIA:KFDIA,KLEV) 
INTEGER(KIND=JPIM),INTENT(IN)    :: indfor(KIDIA:KFDIA,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: forfrac(KIDIA:KFDIA,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: forfac(KIDIA:KFDIA,KLEV) 
! ---------------------------------------------------------------------------

INTEGER(KIND=JPIM) :: ind0,ind1
INTEGER(KIND=JPIM) :: inds,indf

INTEGER(KIND=JPIM) :: IG, lay
REAL(KIND=JPRB) :: taufor,tauself
    !     local integer arrays
    INTEGER(KIND=JPIM) :: laytrop_min, laytrop_max
    integer(KIND=JPIM) :: ixc(KLEV), ixlow(KFDIA,KLEV), ixhigh(KFDIA,KLEV)
    INTEGER(KIND=JPIM) :: ich, icl, ixc0, ixp, jc, jl

    !$ACC DATA PRESENT(taug, P_TAUAERL, fac00, fac01, fac10, fac11, jp, jt, jt1, &
    !$ACC             colh2o, laytrop, fracs, selffac, selffrac, indself, &
    !$ACC             indfor, forfrac, forfac)

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

!     Compute the optical depth by interpolating in ln(pressure) and 
!     temperature.  

      ! Lower atmosphere loop
      !$ACC WAIT
      !$ACC PARALLEL DEFAULT(NONE) ASYNC(1)
      !$ACC LOOP GANG VECTOR COLLAPSE(2) PRIVATE(ind0, ind1, inds, indf)
      do lay = 1, laytrop_min
        do jl = KIDIA, KFDIA

          ind0 = ((jp(jl,lay)-1)*5+(jt(jl,lay)-1))*nspa(10) + 1
          ind1 = (jp(jl,lay)*5+(jt1(jl,lay)-1))*nspa(10) + 1
          inds = indself(jl,lay)
          indf = indfor(jl,lay)
          !$ACC LOOP SEQ PRIVATE(tauself, taufor)
!$NEC unroll(NG10)
          do ig = 1, ng10
            tauself = selffac(jl,lay) * (selfref(inds,ig) + selffrac(jl,lay) * &
                 (selfref(inds+1,ig) - selfref(inds,ig)))
            taufor = forfac(jl,lay) * (forref(indf,ig) + forfrac(jl,lay) * &
                 (forref(indf+1,ig) - forref(indf,ig)))
            taug(jl,ngs9+ig,lay) = colh2o(jl,lay) * &
                 (fac00(jl,lay) * absa(ind0,ig) + &
                 fac10(jl,lay) * absa(ind0+1,ig) + &
                 fac01(jl,lay) * absa(ind1,ig) + &
                 fac11(jl,lay) * absa(ind1+1,ig))  &
                 + tauself + taufor
            fracs(jl,ngs9+ig,lay) = fracrefa(ig)
          enddo
        enddo

      enddo
      !$ACC END PARALLEL

      ! Upper atmosphere loop
      !$ACC PARALLEL DEFAULT(NONE) ASYNC(1)
      !$ACC LOOP GANG VECTOR COLLAPSE(2) PRIVATE(ind0, ind1, indf)
      do lay = laytrop_max+1, KLEV
        do jl = KIDIA, KFDIA

          ind0 = ((jp(jl,lay)-13)*5+(jt(jl,lay)-1))*nspb(10) + 1
          ind1 = ((jp(jl,lay)-12)*5+(jt1(jl,lay)-1))*nspb(10) + 1
          indf = indfor(jl,lay)
          !$ACC LOOP SEQ PRIVATE(taufor)
!$NEC unroll(NG10)
          do ig = 1, ng10
            taufor = forfac(jl,lay) * (forref(indf,ig) + forfrac(jl,lay) * &
                 (forref(indf+1,ig) - forref(indf,ig)))
            taug(jl,ngs9+ig,lay) = colh2o(jl,lay) * &
                 (fac00(jl,lay) * absb(ind0,ig) + &
                 fac10(jl,lay) * absb(ind0+1,ig) + &
                 fac01(jl,lay) * absb(ind1,ig) +  &
                 fac11(jl,lay) * absb(ind1+1,ig)) &
                 + taufor
            fracs(jl,ngs9+ig,lay) = fracrefb(ig)
          enddo
        enddo

      enddo
      !$ACC END PARALLEL

      IF (laytrop_max /= laytrop_min) THEN
        ! Mixed loop
        ! Lower atmosphere part
        !$ACC PARALLEL DEFAULT(NONE) ASYNC(1)
        !$ACC LOOP GANG VECTOR COLLAPSE(2) PRIVATE(ind0, ind1, inds, indf)
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

            ind0 = ((jp(jl,lay)-1)*5+(jt(jl,lay)-1))*nspa(10) + 1
            ind1 = (jp(jl,lay)*5+(jt1(jl,lay)-1))*nspa(10) + 1
            inds = indself(jl,lay)
            indf = indfor(jl,lay)
!$NEC unroll(NG10)
            !$ACC LOOP SEQ PRIVATE(tauself, taufor)
            do ig = 1, ng10
              tauself = selffac(jl,lay) * (selfref(inds,ig) + selffrac(jl,lay) * &
                  (selfref(inds+1,ig) - selfref(inds,ig)))
              taufor = forfac(jl,lay) * (forref(indf,ig) + forfrac(jl,lay) * &
                  (forref(indf+1,ig) - forref(indf,ig)))
              taug(jl,ngs9+ig,lay) = colh2o(jl,lay) * &
                  (fac00(jl,lay) * absa(ind0,ig) + &
                  fac10(jl,lay) * absa(ind0+1,ig) + &
                  fac01(jl,lay) * absa(ind1,ig) + &
                  fac11(jl,lay) * absa(ind1+1,ig))  &
                  + tauself + taufor
              fracs(jl,ngs9+ig,lay) = fracrefa(ig)
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

            ind0 = ((jp(jl,lay)-13)*5+(jt(jl,lay)-1))*nspb(10) + 1
            ind1 = ((jp(jl,lay)-12)*5+(jt1(jl,lay)-1))*nspb(10) + 1
            indf = indfor(jl,lay)
!$NEC unroll(NG10)
            !$ACC LOOP SEQ PRIVATE (taufor)
            do ig = 1, ng10
              taufor = forfac(jl,lay) * (forref(indf,ig) + forfrac(jl,lay) * &
                  (forref(indf+1,ig) - forref(indf,ig)))
              taug(jl,ngs9+ig,lay) = colh2o(jl,lay) * &
                  (fac00(jl,lay) * absb(ind0,ig) + &
                  fac10(jl,lay) * absb(ind0+1,ig) + &
                  fac01(jl,lay) * absb(ind1,ig) +  &
                  fac11(jl,lay) * absb(ind1+1,ig)) &
                  + taufor
              fracs(jl,ngs9+ig,lay) = fracrefb(ig)
            enddo
#ifdef _OPENACC
           endif
#endif
          enddo

        enddo
        !$ACC END PARALLEL

      ENDIF

      !$ACC END DATA

END SUBROUTINE RRTM_TAUMOL10
