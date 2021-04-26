!----------------------------------------------------------------------------
SUBROUTINE RRTM_TAUMOL6 (KIDIA,KFDIA,KLEV,taug,wx,&
 & P_TAUAERL,fac00,fac01,fac10,fac11,forfac,forfrac,indfor,jp,jt,jt1,&
 & colh2o,colco2,coldry,laytrop,selffac,selffrac,indself,fracs,minorfrac,indminor)  

!     BAND 6:  820-980 cm-1 (low - H2O; high - nothing)

!     AUTHOR.
!     -------
!      JJMorcrette, ECMWF

!     MODIFICATIONS.
!     --------------
!      M.Hamrud      01-Oct-2003 CY28 Cleaning
!      NEC           25-Oct-2007 Optimisations
!      JJMorcrette 20110613 flexible number of g-points
!      ABozzo 201306 updated to rrtmg v4.85
!     band 6:  820-980 cm-1 (low key - h2o; low minor - co2)
!                           (high key - nothing; high minor - cfc11, cfc12)
! ---------------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK

USE PARRRTM  , ONLY : JPBAND ,JPXSEC
USE YOERRTM  , ONLY : JPGPT  ,NG6   ,NGS5
USE YOERRTWN , ONLY : NSPA   
USE YOERRTA6 , ONLY : ABSA   ,KA_MCO2 ,CFC11ADJ , CFC12  ,&
 & FRACREFA,SELFREF,FORREF  
USE YOERRTRF, ONLY : CHI_MLS

IMPLICIT NONE

INTEGER(KIND=JPIM),INTENT(IN)    :: KIDIA
INTEGER(KIND=JPIM),INTENT(IN)    :: KFDIA
INTEGER(KIND=JPIM),INTENT(IN)    :: KLEV 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: taug(KIDIA:KFDIA,JPGPT,KLEV) 
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
REAL(KIND=JPRB)   ,INTENT(IN)    :: colco2(KIDIA:KFDIA,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: coldry(KIDIA:KFDIA,KLEV) 
INTEGER(KIND=JPIM),INTENT(IN)    :: laytrop(KIDIA:KFDIA) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: selffac(KIDIA:KFDIA,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: selffrac(KIDIA:KFDIA,KLEV) 
INTEGER(KIND=JPIM),INTENT(IN)    :: indself(KIDIA:KFDIA,KLEV) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: fracs(KIDIA:KFDIA,JPGPT,KLEV) 

INTEGER(KIND=JPIM),INTENT(IN)   :: indfor(KIDIA:KFDIA,KLEV)
REAL(KIND=JPRB)   ,INTENT(IN)   :: forfac(KIDIA:KFDIA,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)   :: forfrac(KIDIA:KFDIA,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)   :: minorfrac(KIDIA:KFDIA,KLEV)
INTEGER(KIND=JPIM),INTENT(IN)   :: indminor(KIDIA:KFDIA,KLEV)

! ---------------------------------------------------------------------------

INTEGER(KIND=JPIM) :: ind0,ind1,inds,indf,indm

INTEGER(KIND=JPIM) :: IG, lay

REAL(KIND=JPRB) :: adjfac,adjcolco2,ratco2,chi_co2
REAL(KIND=JPRB) :: taufor,tauself,absco2
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

! Minor gas mapping level:
!     lower - co2, p = 706.2720 mb, t = 294.2 k
!     upper - cfc11, cfc12


!     Compute the optical depth by interpolating in ln(pressure) and
!     temperature. The water vapor self- and foreign- continuum is interpolated
!     (in temperature) separately.  

      ! Lower atmosphere loop
      do lay = 1, laytrop_min
        do jl = KIDIA, KFDIA

          ! In atmospheres where the amount of CO2 is too great to be considered
          ! a minor species, adjust the column amount of CO2 by an empirical factor
          ! to obtain the proper contribution.
          chi_co2 = colco2(jl,lay)/(coldry(jl,lay))
          ratco2 = 1.e20_JPRB*chi_co2/chi_mls(2,jp(jl,lay)+1)
          if (ratco2 .gt. 3.0_JPRB) then
            adjfac = 2.0_JPRB+(ratco2-2.0_JPRB)**0.77_JPRB
            adjcolco2 = adjfac*chi_mls(2,jp(jl,lay)+1)*coldry(jl,lay)*1.e-20_JPRB
          else
            adjcolco2 = colco2(jl,lay)
          endif

          ind0 = ((jp(jl,lay)-1)*5+(jt(jl,lay)-1))*nspa(6) + 1
          ind1 = (jp(jl,lay)*5+(jt1(jl,lay)-1))*nspa(6) + 1
          inds = indself(jl,lay)
          indf = indfor(jl,lay)
          indm = indminor(jl,lay)
!$NEC unroll(NG6)
          do ig = 1, ng6
            tauself = selffac(jl,lay) * (selfref(inds,ig) + selffrac(jl,lay) * &
                 (selfref(inds+1,ig) - selfref(inds,ig)))
            taufor =  forfac(jl,lay) * (forref(indf,ig) + forfrac(jl,lay) * &
               (forref(indf+1,ig) - forref(indf,ig)))
            absco2 =  (ka_mco2(indm,ig) + minorfrac(jl,lay) * &
                 (ka_mco2(indm+1,ig) - ka_mco2(indm,ig)))
            taug(jl,ngs5+ig,lay) = colh2o(jl,lay) * &
                 (fac00(jl,lay) * absa(ind0,ig) + &
                 fac10(jl,lay) * absa(ind0+1,ig) + &
                 fac01(jl,lay) * absa(ind1,ig) +  &
                 fac11(jl,lay) * absa(ind1+1,ig))  &
                 + tauself + taufor &
                 + adjcolco2 * absco2 &
                 + wx(jl,2,lay) * cfc11adj(ig) &
                 + wx(jl,3,lay) * cfc12(ig)
            fracs(jl,ngs5+ig,lay) = fracrefa(ig)
          enddo
        enddo
      enddo

      ! Upper atmosphere loop
      ! Nothing important goes on above laytrop in this band.
      do ig = 1, ng6
        do lay = laytrop_max+1, KLEV
          do jl = KIDIA, KFDIA
            taug(jl,ngs5+ig,lay) = 0.0_JPRB &
                 + wx(jl,2,lay) * cfc11adj(ig) &
                 + wx(jl,3,lay) * cfc12(ig)
            fracs(jl,ngs5+ig,lay) = fracrefa(ig)
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

          ! In atmospheres where the amount of CO2 is too great to be considered
          ! a minor species, adjust the column amount of CO2 by an empirical factor
          ! to obtain the proper contribution.
          chi_co2 = colco2(jl,lay)/(coldry(jl,lay))
          ratco2 = 1.e20_JPRB*chi_co2/chi_mls(2,jp(jl,lay)+1)
          if (ratco2 .gt. 3.0_JPRB) then
            adjfac = 2.0_JPRB+(ratco2-2.0_JPRB)**0.77_JPRB
            adjcolco2 = adjfac*chi_mls(2,jp(jl,lay)+1)*coldry(jl,lay)*1.e-20_JPRB
          else
            adjcolco2 = colco2(jl,lay)
          endif

          ind0 = ((jp(jl,lay)-1)*5+(jt(jl,lay)-1))*nspa(6) + 1
          ind1 = (jp(jl,lay)*5+(jt1(jl,lay)-1))*nspa(6) + 1
          inds = indself(jl,lay)
          indf = indfor(jl,lay)
          indm = indminor(jl,lay)
!$NEC unroll(NG6)
          do ig = 1, ng6
            tauself = selffac(jl,lay) * (selfref(inds,ig) + selffrac(jl,lay) * &
                 (selfref(inds+1,ig) - selfref(inds,ig)))
            taufor =  forfac(jl,lay) * (forref(indf,ig) + forfrac(jl,lay) * &
               (forref(indf+1,ig) - forref(indf,ig)))
            absco2 =  (ka_mco2(indm,ig) + minorfrac(jl,lay) * &
                 (ka_mco2(indm+1,ig) - ka_mco2(indm,ig)))
            taug(jl,ngs5+ig,lay) = colh2o(jl,lay) * &
                 (fac00(jl,lay) * absa(ind0,ig) + &
                 fac10(jl,lay) * absa(ind0+1,ig) + &
                 fac01(jl,lay) * absa(ind1,ig) +  &
                 fac11(jl,lay) * absa(ind1+1,ig))  &
                 + tauself + taufor &
                 + adjcolco2 * absco2 &
                 + wx(jl,2,lay) * cfc11adj(ig) &
                 + wx(jl,3,lay) * cfc12(ig)
            fracs(jl,ngs5+ig,lay) = fracrefa(ig)
          enddo
        enddo

        ! Upper atmosphere part
        ! Nothing important goes on above laytrop in this band.
        ixc0 = KFDIA - KIDIA + 1 - ixc0

        do ig = 1, ng6
!$NEC ivdep
          do ixp = 1, ixc0
            jl = ixhigh(ixp,lay)
            taug(jl,ngs5+ig,lay) = 0.0_JPRB &
                 + wx(jl,2,lay) * cfc11adj(ig) &
                 + wx(jl,3,lay) * cfc12(ig)
            fracs(jl,ngs5+ig,lay) = fracrefa(ig)
          enddo
        enddo

      enddo

END SUBROUTINE RRTM_TAUMOL6
