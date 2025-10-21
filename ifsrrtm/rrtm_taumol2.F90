!----------------------------------------------------------------------------
SUBROUTINE RRTM_TAUMOL2 (KIDIA,KFDIA,KLEV,taug,PAVEL,coldry,&
 & P_TAUAERL,fac00,fac01,fac10,fac11,forfac,forfrac,indfor,jp,jt,jt1,&
 & colh2o,laytrop,selffac,selffrac,indself,fracs,laytrop_min,laytrop_max)

!     BAND 2:  250-500 cm-1 (low - H2O; high - H2O)

!     AUTHOR.
!     -------
!      JJMorcrette, ECMWF

!     MODIFICATIONS.
!     --------------
!      M.Hamrud      01-Oct-2003 CY28 Cleaning
!      NEC           25-Oct-2007 Optimisations
!      JJMorcrette 20110613 flexible number of g-points
!      ABozzo 201305 updated to rrtmg_lw_v4.85:
!*********
!     band 2:  350-500 cm-1 (low key - h2o; high key - h2o)
!
!     note: previous version of rrtm band 2:
!           250 - 500 cm-1 (low - h2o; high - h2o)
!
! ---------------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM     ,JPRB

USE PARRRTM  , ONLY : JPBAND
USE YOERRTM  , ONLY : JPGPT  ,NG2   ,NGS1
USE YOERRTWN , ONLY : NSPA   ,NSPB
USE YOERRTA2 , ONLY : ABSA   ,ABSB   ,FRACREFA, FRACREFB,&
 & FORREF   ,SELFREF

IMPLICIT NONE

INTEGER(KIND=JPIM),INTENT(IN)    :: KIDIA
INTEGER(KIND=JPIM),INTENT(IN)    :: KFDIA
INTEGER(KIND=JPIM),INTENT(IN)    :: KLEV
REAL(KIND=JPRB)   ,INTENT(INOUT) :: taug(KIDIA:KFDIA,JPGPT,KLEV)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PAVEL(KIDIA:KFDIA,KLEV) ! Layer pressures (hPa)
REAL(KIND=JPRB)   ,INTENT(IN)    :: coldry(KIDIA:KFDIA,KLEV)
REAL(KIND=JPRB)   ,INTENT(IN)    :: P_TAUAERL(KIDIA:KFDIA,KLEV,JPBAND)
REAL(KIND=JPRB)   ,INTENT(IN)    :: fac00(KIDIA:KFDIA,KLEV)
REAL(KIND=JPRB)   ,INTENT(IN)    :: fac01(KIDIA:KFDIA,KLEV)
REAL(KIND=JPRB)   ,INTENT(IN)    :: fac10(KIDIA:KFDIA,KLEV)
REAL(KIND=JPRB)   ,INTENT(IN)    :: fac11(KIDIA:KFDIA,KLEV)
REAL(KIND=JPRB)   ,INTENT(IN)    :: forfrac(KIDIA:KFDIA,KLEV)
REAL(KIND=JPRB)   ,INTENT(IN)    :: forfac(KIDIA:KFDIA,KLEV)
INTEGER(KIND=JPIM),INTENT(IN)    :: jp(KIDIA:KFDIA,KLEV)
INTEGER(KIND=JPIM),INTENT(IN)    :: jt(KIDIA:KFDIA,KLEV)
INTEGER(KIND=JPIM),INTENT(IN)    :: jt1(KIDIA:KFDIA,KLEV)
REAL(KIND=JPRB)   ,INTENT(IN)    :: colh2o(KIDIA:KFDIA,KLEV)
INTEGER(KIND=JPIM),INTENT(IN)    :: laytrop(KIDIA:KFDIA)
REAL(KIND=JPRB)   ,INTENT(IN)    :: selffac(KIDIA:KFDIA,KLEV)
REAL(KIND=JPRB)   ,INTENT(IN)    :: selffrac(KIDIA:KFDIA,KLEV)
INTEGER(KIND=JPIM),INTENT(IN)    :: indself(KIDIA:KFDIA,KLEV)
INTEGER(KIND=JPIM),INTENT(IN)    :: indfor(KIDIA:KFDIA,KLEV)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: fracs(KIDIA:KFDIA,JPGPT,KLEV)
INTEGER(KIND=JPIM), OPTIONAL, INTENT(INOUT) :: laytrop_min, laytrop_max

! ---------------------------------------------------------------------------

INTEGER(KIND=JPIM) :: ind0,ind1,inds, indf

INTEGER(KIND=JPIM) :: IG, lay

REAL(KIND=JPRB) :: taufor,tauself,corradj,pp
    !     local integer arrays
    integer(KIND=JPIM) :: ixc(KLEV), ixlow(KFDIA,KLEV), ixhigh(KFDIA,KLEV)
    INTEGER(KIND=JPIM) :: ich, icl, ixc0, ixp, jc, jl
INTEGER(KIND=JPIM) :: llaytrop_min, llaytrop_max

#include "rrtm_utils.intfb.h"

    if (present(laytrop_min) .AND. present(laytrop_max)) then
       llaytrop_min = laytrop_min
       llaytrop_max = laytrop_max
    else
       CALL COMPUTE_LAYTROP_MIN_MAX(KIDIA, KFDIA, LAYTROP, llaytrop_min, llaytrop_max)
    endif

    !$ACC DATA PRESENT(taug, PAVEL, coldry, P_TAUAERL, fac00, fac01, fac10, &
    !$ACC             fac11, forfrac, forfac, jp, jt, jt1, colh2o, laytrop, &
    !$ACC             selffac, selffrac, indself, indfor, fracs)
#ifndef __NVCOMPILER
    !$OMP TARGET DATA MAP(PRESENT, ALLOC: taug, PAVEL, coldry, P_TAUAERL, fac00, fac01, fac10, &
    !$OMP             fac11, forfrac, forfac, jp, jt, jt1, colh2o, laytrop, &
    !$OMP             selffac, selffrac, indself, indfor, fracs)
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


!     Compute the optical depth by interpolating in ln(pressure) and
!     temperature.  Below LAYTROP, the water vapor self-continuum is
!     interpolated (in temperature) separately.

      ! Lower atmosphere loop
      !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO COLLAPSE(2) PRIVATE(ind0, ind1, inds, indf, pp, corradj, tauself, taufor)
      !$ACC WAIT
      !$ACC PARALLEL DEFAULT(NONE) ASYNC(1)
      !$ACC LOOP GANG VECTOR COLLAPSE(2) PRIVATE(ind0, ind1, inds, indf, pp, corradj)
      do lay = 1, llaytrop_min
        do jl = KIDIA, KFDIA

          ind0 = ((jp(jl,lay)-1)*5+(jt(jl,lay)-1))*nspa(2) + 1
          ind1 = (jp(jl,lay)*5+(jt1(jl,lay)-1))*nspa(2) + 1
          inds = indself(jl,lay)
          indf = indfor(jl,lay)
          pp = pavel(jl,lay)
          corradj = 1._JPRB - .05_JPRB * (pp - 100._JPRB) / 900._JPRB
          !$ACC LOOP SEQ PRIVATE(tauself, taufor)
!$NEC unroll(NG2)
          do ig = 1, ng2
            tauself = selffac(jl,lay) * (selfref(inds,ig) + selffrac(jl,lay) * &
                 (selfref(inds+1,ig) - selfref(inds,ig)))
            taufor =  forfac(jl,lay) * (forref(indf,ig) + forfrac(jl,lay) * &
                 (forref(indf+1,ig) - forref(indf,ig)))
            taug(jl,ngs1+ig,lay) = corradj * (colh2o(jl,lay) * &
                 (fac00(jl,lay) * absa(ind0,ig) + &
                 fac10(jl,lay) * absa(ind0+1,ig) + &
                 fac01(jl,lay) * absa(ind1,ig) + &
                 fac11(jl,lay) * absa(ind1+1,ig)) &
                 + tauself + taufor)
            fracs(jl,ngs1+ig,lay) = fracrefa(ig)
          enddo
        enddo
      enddo
      !$ACC END PARALLEL
      !$OMP END TARGET TEAMS DISTRIBUTE PARALLEL DO

      ! Upper atmosphere loop
      !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO COLLAPSE(2) PRIVATE(ind0, ind1, indf, taufor)
      !$ACC PARALLEL DEFAULT(NONE) ASYNC(1)
      !$ACC LOOP GANG VECTOR COLLAPSE(2) PRIVATE(ind0, ind1, indf)
      do lay = llaytrop_max+1, KLEV
        do jl = KIDIA, KFDIA

          ind0 = ((jp(jl,lay)-13)*5+(jt(jl,lay)-1))*nspb(2) + 1
          ind1 = ((jp(jl,lay)-12)*5+(jt1(jl,lay)-1))*nspb(2) + 1
          indf = indfor(jl,lay)
          !$ACC LOOP SEQ PRIVATE(taufor)
!$NEC unroll(NG2)
          do ig = 1, ng2
            taufor =  forfac(jl,lay) * (forref(indf,ig) + &
                 forfrac(jl,lay) * (forref(indf+1,ig) - forref(indf,ig)))
            taug(jl,ngs1+ig,lay) = colh2o(jl,lay) * &
                 (fac00(jl,lay) * absb(ind0,ig) + &
                 fac10(jl,lay) * absb(ind0+1,ig) + &
                 fac01(jl,lay) * absb(ind1,ig) + &
                 fac11(jl,lay) * absb(ind1+1,ig)) &
                 + taufor
            fracs(jl,ngs1+ig,lay) = fracrefb(ig)
          enddo
        enddo

      enddo
      !$ACC END PARALLEL
      !$OMP END TARGET TEAMS DISTRIBUTE PARALLEL DO

      IF (llaytrop_max /= llaytrop_min) THEN
        ! Mixed loop
        ! Lower atmosphere part
        !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO COLLAPSE(2) PRIVATE(ind0, ind1, inds, indf, pp, corradj, tauself, taufor)
        !$ACC PARALLEL DEFAULT(NONE) ASYNC(1)
        !$ACC LOOP GANG VECTOR COLLAPSE(2) PRIVATE(ind0, ind1, inds, indf, pp, corradj)
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

            ind0 = ((jp(jl,lay)-1)*5+(jt(jl,lay)-1))*nspa(2) + 1
            ind1 = (jp(jl,lay)*5+(jt1(jl,lay)-1))*nspa(2) + 1
            inds = indself(jl,lay)
            indf = indfor(jl,lay)
            pp = pavel(jl,lay)
            corradj = 1._JPRB - .05_JPRB * (pp - 100._JPRB) / 900._JPRB
!$NEC unroll(NG2)
            !$ACC LOOP SEQ PRIVATE(tauself, taufor)
            do ig = 1, ng2
              tauself = selffac(jl,lay) * (selfref(inds,ig) + selffrac(jl,lay) * &
                  (selfref(inds+1,ig) - selfref(inds,ig)))
              taufor =  forfac(jl,lay) * (forref(indf,ig) + forfrac(jl,lay) * &
                  (forref(indf+1,ig) - forref(indf,ig)))
              taug(jl,ngs1+ig,lay) = corradj * (colh2o(jl,lay) * &
                 (fac00(jl,lay) * absa(ind0,ig) + &
                 fac10(jl,lay) * absa(ind0+1,ig) + &
                 fac01(jl,lay) * absa(ind1,ig) + &
                 fac11(jl,lay) * absa(ind1+1,ig)) &
                 + tauself + taufor)
              fracs(jl,ngs1+ig,lay) = fracrefa(ig)
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

            ind0 = ((jp(jl,lay)-13)*5+(jt(jl,lay)-1))*nspb(2) + 1
            ind1 = ((jp(jl,lay)-12)*5+(jt1(jl,lay)-1))*nspb(2) + 1
            indf = indfor(jl,lay)
!$NEC unroll(NG2)
            !$ACC LOOP SEQ PRIVATE(taufor)
            do ig = 1, ng2
              taufor =  forfac(jl,lay) * (forref(indf,ig) + &
                  forfrac(jl,lay) * (forref(indf+1,ig) - forref(indf,ig)))
              taug(jl,ngs1+ig,lay) = colh2o(jl,lay) * &
                  (fac00(jl,lay) * absb(ind0,ig) + &
                  fac10(jl,lay) * absb(ind0+1,ig) + &
                  fac01(jl,lay) * absb(ind1,ig) + &
                  fac11(jl,lay) * absb(ind1+1,ig)) &
                  + taufor
              fracs(jl,ngs1+ig,lay) = fracrefb(ig)
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
#ifndef __NVCOMPILER
    !$OMP END TARGET DATA
#endif

END SUBROUTINE RRTM_TAUMOL2
