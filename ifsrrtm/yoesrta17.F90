MODULE YOESRTA17

USE PARKIND1  ,ONLY : JPIM     ,JPRB,JPRD

IMPLICIT NONE

PUBLIC

SAVE

!     -----------------------------------------------------------------
!*    ** *YOESRTA17* - SRTM COEFFICIENTS FOR INTERVAL 17
!     BAND 17:  3250-4000 cm-1 (low - H2O,CO2; high - H2O,CO2)
!     -----------------------------------------------------------------

INTEGER(KIND=JPIM), PARAMETER ::  JPG = 16, NG17 = 16, NGS16=16

REAL(KIND=JPRB) :: KA(JPG,9,5,13)  
REAL(KIND=JPRB) :: KB(JPG,5,5,13:59)
REAL(KIND=JPRD) :: KA_D(JPG,9,5,13)  
REAL(KIND=JPRD) :: KB_D(JPG,5,5,13:59)
REAL(KIND=JPRB) :: SELFREF(JPG,10)  ,FORREF(JPG,4)
REAL(KIND=JPRB) :: SFLUXREF(JPG,5)
REAL(KIND=JPRB) :: RAYL              ,STRRAT
INTEGER(KIND=JPIM) :: LAYREFFR

REAL(KIND=JPRB) :: KAC(NG17,9,5,13),ABSA(NG17,585)
REAL(KIND=JPRB) :: KBC(NG17,5,5,13:59),ABSB(NG17,1175)
REAL(KIND=JPRB) :: SELFREFC(NG17,10),FORREFC(NG17,4)
REAL(KIND=JPRB) :: SFLUXREFC(NG17,5)

!EQUIVALENCE (KA(1,1,1,1),ABSA(1,1)), (KB(1,1,13,1),ABSB(1,1))
EQUIVALENCE (KAC(1,1,1,1),ABSA(1,1)), (KBC(1,1,1,13),ABSB(1,1))
!     -----------------------------------------------------------------
!        * E.C.M.W.F. PHYSICS PACKAGE ** RRTM SW RADIATION **

!     J.-J. MORCRETTE       E.C.M.W.F.      02/10/29
!     M. J. IACONO          AER             12/09/03

!  NAME     TYPE     PURPOSE
!  ----   : ----   : ---------------------------------------------------
! KA      : REAL     absorption coefficient of major absorber
! KB      : REAL     absorption coefficient of secondary absorber
! SELFREF : REAL     self brodening coefficient for water vapour
! FORREF  : REAL     foreign broadening coefficient for water vapour
! SFLUXREF: REAL     Incident solar radiation in the spectral interval
! RAYL    : REAL     Rayleigh scattering parameter
! STRRAT  : REAL     weighting factor for the transition between tropospheric 
!                    and stratospheric computations
! LAYREFFR: INTEGER  reference level for the transition
! KAC     : REAL     Reduced g-point array for KA
! KBC     : REAL     Reduced g-point array for KB
! SELFREFC: REAL     Reduced g-point array for SELFREF
! FORREFC : REAL     Reduced g-point array for FORREF
!SFLUXREFC: REAL     Reduced g-point array for SFLUXREF
!     -----------------------------------------------------------------
END MODULE YOESRTA17

