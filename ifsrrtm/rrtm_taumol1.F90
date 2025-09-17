SUBROUTINE RRTM_TAUMOL1 (KIDIA,KFDIA,KLEV,taug,PAVEL,&
 & P_TAUAERL,FAC00,FAC01,FAC10,FAC11,FORFAC,FORFRAC,INDFOR,JP,JT,jt1,&
 & COLH2O,LAYTROP,SELFFAC,SELFFRAC,INDSELF,fracs,MINORFRAC,INDMINOR,&
 & SCALEMINORN2,COLBRD,laytrop_min,laytrop_max)

!******************************************************************************
!                                                                             *
!                  Optical depths developed for the                           *
!                                                                             *
!                RAPID RADIATIVE TRANSFER MODEL (RRTM)                        *
!                                                                             *
!            ATMOSPHERIC AND ENVIRONMENTAL RESEARCH, INC.                     *
!                        840 MEMORIAL DRIVE                                   *
!                        CAMBRIDGE, MA 02139                                  *
!                                                                             *
!                           ELI J. MLAWER                                     *
!                         STEVEN J. TAUBMAN                                   *
!                         SHEPARD A. CLOUGH                                   *
!                                                                             *
!                       email:  mlawer@aer.com                                *
!                                                                             *
!        The authors wish to acknowledge the contributions of the             *
!        following people:  Patrick D. Brown, Michael J. Iacono,              *
!        Ronald E. Farren, Luke Chen, Robert Bergstrom.                       *
!                                                                             *
!******************************************************************************
!     TAUMOL                                                                  *
!                                                                             *
!     This file contains the subroutines TAUGBn (where n goes from            *
!     1 to 16).  TAUGBn calculates the optical depths and Planck fractions    *
!     per g-value and layer for band n.                                       *
!                                                                             *
!  Output:  optical depths (unitless)                                         *
!           fractions needed to compute Planck functions at every layer       *
!               and g-value                                                   *
!                                                                             *
!     COMMON /TAUGCOM/  TAUG(MXLAY,MG)                                        *
!     COMMON /PLANKG/   FRACS(MXLAY,MG)                                       *
!                                                                             *
!  Input                                                                      *
!                                                                             *
!     COMMON /FEATURES/ NG(NBANDS),NSPA(NBANDS),NSPB(NBANDS)                  *
!     COMMON /PRECISE/  ONEMINUS                                              *
!     COMMON /PROFILE/  KLEV,PAVEL(MXLAY),TAVEL(MXLAY),                    *
!    &                  PZ(0:MXLAY),TZ(0:MXLAY),TBOUND                        *
!     COMMON /PROFDATA/ LAYTROP,LAYSWTCH,LAYLOW,                              *
!    &                  COLH2O(MXLAY),COLCO2(MXLAY),                          *
!    &                  COLO3(MXLAY),COLN2O(MXLAY),COLCH4(MXLAY),             *
!    &                  COLO2(MXLAY),CO2MULT(MXLAY)                           *
!     COMMON /INTFAC/   FAC00(MXLAY),FAC01(MXLAY),                            *
!    &                  FAC10(MXLAY),FAC11(MXLAY)                             *
!     COMMON /INTIND/   JP(MXLAY),JT(KIDIA:KFDIA,MXLAY),jt1(KIDIA:KFDIA,MXLAY)                        *
!     COMMON /SELF/     SELFFAC(MXLAY), SELFFRAC(MXLAY), INDSELF(KIDIA:KFDIA,MXLAY)       *
!                                                                             *
!     Description:                                                            *
!     NG(IBAND) - number of g-values in band IBAND                            *
!     NSPA(IBAND) - for the lower atmosphere, the number of reference         *
!                   atmospheres that are stored for band IBAND per            *
!                   pressure level and temperature.  Each of these            *
!                   atmospheres has different relative amounts of the         *
!                   key species for the band (i.e. different binary           *
!                   species parameters).                                      *
!     NSPB(IBAND) - same for upper atmosphere                                 *
!     ONEMINUS - since problems are caused in some cases by interpolation     *
!                parameters equal to or greater than 1, for these cases       *
!                these parameters are set to this value, slightly < 1.        *
!     PAVEL - layer pressures (mb)                                            *
!     TAVEL - layer temperatures (degrees K)                                  *
!     PZ - level pressures (mb)                                               *
!     TZ - level temperatures (degrees K)                                     *
!     LAYTROP - layer at which switch is made from one combination of         *
!               key species to another                                        *
!     COLH2O, COLCO2, COLO3, COLN2O, COLCH4 - column amounts of water         *
!               vapor,carbon dioxide, ozone, nitrous ozide, methane,          *
!               respectively (molecules/cm**2)                                *
!     CO2MULT - for bands in which carbon dioxide is implemented as a         *
!               trace species, this is the factor used to multiply the        *
!               band's average CO2 absorption coefficient to get the added    *
!               contribution to the optical depth relative to 355 ppm.        *
!     FACij(lay) - for layer lay, these are factors that are needed to        *
!                  compute the interpolation factors that multiply the        *
!                  appropriate reference k-values.  A value of 0 (1) for      *
!                  i,j indicates that the corresponding factor multiplies     *
!                  reference k-value for the lower (higher) of the two        *
!                  appropriate temperatures, and altitudes, respectively.     *
!     JP - the index of the lower (in altitude) of the two appropriate        *
!          reference pressure levels needed for interpolation                 *
!     JT, jt1 - the indices of the lower of the two appropriate reference     *
!               temperatures needed for interpolation (for pressure           *
!               levels JP and JP+1, respectively)                             *
!     SELFFAC - scale factor needed to water vapor self-continuum, equals     *
!               (water vapor density)/(atmospheric density at 296K and        *
!               1013 mb)                                                      *
!     SELFFRAC - factor needed for temperature interpolation of reference     *
!                water vapor self-continuum data                              *
!     INDSELF - index of the lower of the two appropriate reference           *
!               temperatures needed for the self-continuum interpolation      *
!                                                                             *
!  Data input                                                                 *
!     COMMON /Kn/ KA(NSPA(n),5,13,MG), KB(NSPB(n),5,13:59,MG), SELFREF(10,MG) *
!        (note:  n is the band number)                                        *
!                                                                             *
!     Description:                                                            *
!     KA - k-values for low reference atmospheres (no water vapor             *
!          self-continuum) (units: cm**2/molecule)                            *
!     KB - k-values for high reference atmospheres (all sources)              *
!          (units: cm**2/molecule)                                            *
!     SELFREF - k-values for water vapor self-continuum for reference         *
!               atmospheres (used below LAYTROP)                              *
!               (units: cm**2/molecule)                                       *
!                                                                             *
!     DIMENSION ABSA(65*NSPA(n),MG), ABSB(235*NSPB(n),MG)                     *
!     EQUIVALENCE (KA,ABSA),(KB,ABSB)                                         *
!                                                                             *
!******************************************************************************

!     BAND 1:  10-250 cm-1 (low - H2O; high - H2O)

!     AUTHOR.
!     -------
!      JJMorcrette, ECMWF, from
!      Eli J. Mlawer, Atmospheric & Environmental Research.
!      (Revised by Michael J. Iacono, Atmospheric & Environmental Research.)

!     MODIFICATIONS.
!     --------------
!      D Salmond   2000-05-15 speed-up
!      JJMorcrette 2000-05-17 speed-up
!      M.Hamrud      01-Oct-2003 CY28 Cleaning
!      NEC           25-Oct-2007 Optimisations
!      JJMorcrette 20110613 flexible number of g-points
!      ABozzo 200130517 updated to rrtmg_lw_v4.85:
!*********
!     band 1:  10-350 cm-1 (low key - h2o; low minor - n2)
!                          (high key - h2o; high minor - n2)
!
!     note: previous versions of rrtm band 1:
!           10-250 cm-1 (low - h2o; high - h2o)
! ---------------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM     ,JPRB

USE PARRRTM  , ONLY : JPBAND
USE YOERRTM  , ONLY : JPGPT  ,NG1
USE YOERRTWN , ONLY : NSPA   ,NSPB
USE YOERRTA1 , ONLY : ABSA   ,ABSB   ,FRACREFA, FRACREFB,&
 & FORREF   ,SELFREF,  KA_MN2, KB_MN2

IMPLICIT NONE

INTEGER(KIND=JPIM),INTENT(IN)    :: KIDIA
INTEGER(KIND=JPIM),INTENT(IN)    :: KFDIA
INTEGER(KIND=JPIM),INTENT(IN)    :: KLEV
REAL(KIND=JPRB)   ,INTENT(IN)    :: PAVEL(KIDIA:KFDIA,KLEV) ! Layer pressures (hPa)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: taug(KIDIA:KFDIA,JPGPT,KLEV)
REAL(KIND=JPRB)   ,INTENT(IN)    :: P_TAUAERL(KIDIA:KFDIA,KLEV,JPBAND)
REAL(KIND=JPRB)   ,INTENT(IN)    :: FAC00(KIDIA:KFDIA,KLEV)
REAL(KIND=JPRB)   ,INTENT(IN)    :: FAC01(KIDIA:KFDIA,KLEV)
REAL(KIND=JPRB)   ,INTENT(IN)    :: FAC10(KIDIA:KFDIA,KLEV)
REAL(KIND=JPRB)   ,INTENT(IN)    :: FAC11(KIDIA:KFDIA,KLEV)
REAL(KIND=JPRB)   ,INTENT(IN)    :: FORFAC(KIDIA:KFDIA,KLEV)
REAL(KIND=JPRB)   ,INTENT(IN)    :: FORFRAC(KIDIA:KFDIA,KLEV)
INTEGER(KIND=JPIM),INTENT(IN)    :: JP(KIDIA:KFDIA,KLEV)
INTEGER(KIND=JPIM),INTENT(IN)    :: JT(KIDIA:KFDIA,KLEV)
INTEGER(KIND=JPIM),INTENT(IN)    :: jt1(KIDIA:KFDIA,KLEV)
REAL(KIND=JPRB)   ,INTENT(IN)    :: COLH2O(KIDIA:KFDIA,KLEV)
INTEGER(KIND=JPIM),INTENT(IN)    :: LAYTROP(KIDIA:KFDIA)
REAL(KIND=JPRB)   ,INTENT(IN)    :: SELFFAC(KIDIA:KFDIA,KLEV)
REAL(KIND=JPRB)   ,INTENT(IN)    :: SELFFRAC(KIDIA:KFDIA,KLEV)
REAL(KIND=JPRB)   ,INTENT(IN)    :: MINORFRAC(KIDIA:KFDIA,KLEV)
INTEGER(KIND=JPIM),INTENT(IN)    :: INDSELF(KIDIA:KFDIA,KLEV)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: fracs(KIDIA:KFDIA,JPGPT,KLEV)

INTEGER(KIND=JPIM),INTENT(IN)    :: INDFOR(KIDIA:KFDIA,KLEV)
INTEGER(KIND=JPIM),INTENT(IN)    :: INDMINOR(KIDIA:KFDIA,KLEV)
REAL(KIND=JPRB)   ,INTENT(IN)    :: SCALEMINORN2(KIDIA:KFDIA,KLEV)
REAL(KIND=JPRB)   ,INTENT(IN)    :: COLBRD(KIDIA:KFDIA,KLEV)
INTEGER(KIND=JPIM),INTENT(IN)    :: laytrop_min, laytrop_max
! ---------------------------------------------------------------------------

INTEGER(KIND=JPIM) :: IND0,IND1,INDS
INTEGER(KIND=JPIM) :: INDF,INDM

INTEGER(KIND=JPIM) :: IG, LAY
REAL(KIND=JPRB) :: taufor,tauself,corradj,pp,scalen2, taun2
    !     local integer arrays
    integer(KIND=JPIM) :: ixc(KLEV), ixlow(KFDIA,KLEV), ixhigh(KFDIA,KLEV)
    INTEGER(KIND=JPIM) :: ich, icl, ixc0, ixp, jc, jl

    !$ACC DATA PRESENT(PAVEL, taug, P_TAUAERL, FAC00, FAC01, FAC10, FAC11, &
    !$ACC             FORFAC, FORFRAC, JP, JT, jt1, COLH2O, LAYTROP, SELFFAC, &
    !$ACC             SELFFRAC, MINORFRAC, INDSELF, fracs, &
    !$ACC             INDFOR, INDMINOR, SCALEMINORN2, COLBRD)

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


! Minor gas mapping levels:
!     lower - n2, p = 142.5490 mbar, t = 215.70 k
!     upper - n2, p = 142.5490 mbar, t = 215.70 k

!     Compute the optical depth by interpolating in ln(pressure) and
!     temperature.  Below LAYTROP, the water vapor self-continuum and
!     foreign continuum is interpolated (in temperature) separately.

      ! Lower atmosphere loop
      !$ACC WAIT
      !$ACC PARALLEL DEFAULT(NONE) ASYNC(1)
      !$ACC LOOP GANG VECTOR COLLAPSE(2) PRIVATE(ind0, ind1, inds, indf, indm, pp, corradj)
      do lay = 1, laytrop_min
        do jl = KIDIA, KFDIA

          ind0 = ((jp(jl,lay)-1)*5+(jt(jl,lay)-1))*nspa(1) + 1
          ind1 = (jp(jl,lay)*5+(jt1(jl,lay)-1))*nspa(1) + 1
          inds = indself(jl,lay)
          indf = indfor(jl,lay)
          indm = indminor(jl,lay)
          pp = pavel(jl,lay)
          corradj =  1.0_JPRB
          if (pp .lt. 250._JPRB) then
            corradj = 1._JPRB - 0.15_JPRB * (250._JPRB-pp) / 154.4_JPRB
          endif

          scalen2 = colbrd(jl,lay) * scaleminorn2(jl,lay)
          !$ACC LOOP SEQ PRIVATE(tauself, taufor, taun2)
!$NEC unroll(NG1)
          do ig = 1, ng1
            tauself = selffac(jl,lay) * (selfref(inds,ig) + selffrac(jl,lay) * &
                 (selfref(inds+1,ig) - selfref(inds,ig)))
            taufor =  forfac(jl,lay) * (forref(indf,ig) + forfrac(jl,lay) * &
                 (forref(indf+1,ig) -  forref(indf,ig)))
            taun2 = scalen2*(ka_mn2(indm,ig) + &
                 minorfrac(jl,lay) * (ka_mn2(indm+1,ig) - ka_mn2(indm,ig)))
            taug(jl,ig,lay) = corradj * (colh2o(jl,lay) * &
                 (fac00(jl,lay) * absa(ind0,ig) + &
                 fac10(jl,lay) * absa(ind0+1,ig) + &
                 fac01(jl,lay) * absa(ind1,ig) + &
                 fac11(jl,lay) * absa(ind1+1,ig)) &
                 + tauself + taufor + taun2)
            fracs(jl,ig,lay) = fracrefa(ig)
          enddo
        enddo
      enddo
      !$ACC END PARALLEL

      ! Upper atmosphere loop
      !$ACC PARALLEL DEFAULT(NONE) ASYNC(1)
      !$ACC LOOP GANG VECTOR COLLAPSE(2) PRIVATE(ind0, ind1, indf, indm, pp, corradj, scalen2)
      do lay = laytrop_max+1, KLEV
        do jl = KIDIA, KFDIA

          ind0 = ((jp(jl,lay)-13)*5+(jt(jl,lay)-1))*nspb(1) + 1
          ind1 = ((jp(jl,lay)-12)*5+(jt1(jl,lay)-1))*nspb(1) + 1
          indf = indfor(jl,lay)
          indm = indminor(jl,lay)
          pp = pavel(jl,lay)
          corradj =  1._JPRB - 0.15_JPRB * (pp / 95.6_JPRB)

          scalen2 = colbrd(jl,lay) * scaleminorn2(jl,lay)
          !$ACC LOOP SEQ PRIVATE(taufor, taun2)
!$NEC unroll(NG1)
          do ig = 1, ng1
            taufor = forfac(jl,lay) * (forref(indf,ig) + &
                 forfrac(jl,lay) * (forref(indf+1,ig) - forref(indf,ig)))
            taun2 = scalen2*(kb_mn2(indm,ig) + &
                 minorfrac(jl,lay) * (kb_mn2(indm+1,ig) - kb_mn2(indm,ig)))
            taug(jl,ig,lay) = corradj * (colh2o(jl,lay) * &
                 (fac00(jl,lay) * absb(ind0,ig) + &
                 fac10(jl,lay) * absb(ind0+1,ig) + &
                 fac01(jl,lay) * absb(ind1,ig) + &
                 fac11(jl,lay) * absb(ind1+1,ig)) &
                 + taufor + taun2)
            fracs(jl,ig,lay) = fracrefb(ig)
          enddo
        enddo
      enddo
      !$ACC END PARALLEL

      IF (laytrop_max /= laytrop_min) THEN
        ! Mixed loop
        ! Lower atmosphere part
        !$ACC PARALLEL DEFAULT(NONE) ASYNC(1)
        !$ACC LOOP GANG VECTOR COLLAPSE(2) PRIVATE(ind0, ind1, inds, indf, indm, pp, corradj, scalen2)
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

            ind0 = ((jp(jl,lay)-1)*5+(jt(jl,lay)-1))*nspa(1) + 1
            ind1 = (jp(jl,lay)*5+(jt1(jl,lay)-1))*nspa(1) + 1
            inds = indself(jl,lay)
            indf = indfor(jl,lay)
            indm = indminor(jl,lay)
            pp = pavel(jl,lay)
            corradj =  1.0_JPRB
            if (pp .lt. 250._JPRB) then
              corradj = 1._JPRB - 0.15_JPRB * (250._JPRB-pp) / 154.4_JPRB
            endif

            scalen2 = colbrd(jl,lay) * scaleminorn2(jl,lay)
!$NEC unroll(NG1)
            !$ACC LOOP SEQ PRIVATE(tauself, taufor, taun2)
            do ig = 1, ng1
              tauself = selffac(jl,lay) * (selfref(inds,ig) + selffrac(jl,lay) * &
                  (selfref(inds+1,ig) - selfref(inds,ig)))
              taufor =  forfac(jl,lay) * (forref(indf,ig) + forfrac(jl,lay) * &
                  (forref(indf+1,ig) -  forref(indf,ig)))
              taun2 = scalen2*(ka_mn2(indm,ig) + &
                  minorfrac(jl,lay) * (ka_mn2(indm+1,ig) - ka_mn2(indm,ig)))
              taug(jl,ig,lay) = corradj * (colh2o(jl,lay) * &
                  (fac00(jl,lay) * absa(ind0,ig) + &
                  fac10(jl,lay) * absa(ind0+1,ig) + &
                  fac01(jl,lay) * absa(ind1,ig) + &
                  fac11(jl,lay) * absa(ind1+1,ig)) &
                  + tauself + taufor + taun2)
              fracs(jl,ig,lay) = fracrefa(ig)
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

            ind0 = ((jp(jl,lay)-13)*5+(jt(jl,lay)-1))*nspb(1) + 1
            ind1 = ((jp(jl,lay)-12)*5+(jt1(jl,lay)-1))*nspb(1) + 1
            indf = indfor(jl,lay)
            indm = indminor(jl,lay)
            pp = pavel(jl,lay)
            corradj =  1._JPRB - 0.15_JPRB * (pp / 95.6_JPRB)

            scalen2 = colbrd(jl,lay) * scaleminorn2(jl,lay)
!$NEC unroll(NG1)
            !$ACC LOOP SEQ PRIVATE(taufor, taun2)
            do ig = 1, ng1
              taufor = forfac(jl,lay) * (forref(indf,ig) + &
                  forfrac(jl,lay) * (forref(indf+1,ig) - forref(indf,ig)))
              taun2 = scalen2*(kb_mn2(indm,ig) + &
                  minorfrac(jl,lay) * (kb_mn2(indm+1,ig) - kb_mn2(indm,ig)))
              taug(jl,ig,lay) = corradj * (colh2o(jl,lay) * &
                  (fac00(jl,lay) * absb(ind0,ig) + &
                  fac10(jl,lay) * absb(ind0+1,ig) + &
                  fac01(jl,lay) * absb(ind1,ig) + &
                  fac11(jl,lay) * absb(ind1+1,ig)) &
                  + taufor + taun2)
              fracs(jl,ig,lay) = fracrefb(ig)
            enddo
#ifdef _OPENACC
           endif
#endif
          enddo

        enddo
        !$ACC END PARALLEL

      ENDIF

      !$ACC END DATA

END SUBROUTINE RRTM_TAUMOL1
