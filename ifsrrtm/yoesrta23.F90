MODULE YOESRTA23

USE PARKIND1  ,ONLY : JPIM     ,JPRB,JPRD

IMPLICIT NONE

SAVE

!     -----------------------------------------------------------------
!*    ** *YOESRTA23* - SRTM COEFFICIENTS FOR INTERVAL 23
!     BAND 23:  8050-12850 cm-1 (low - H2O; high - nothing)
!     -----------------------------------------------------------------

INTEGER(KIND=JPIM), PARAMETER :: JPG = 16, NG23 = 16

REAL(KIND=JPRB) :: KA(5,13,JPG)   
REAL(KIND=JPRD) :: KA_D(5,13,JPG)   
REAL(KIND=JPRB) :: SELFREF(10,JPG),FORREF(3,JPG)
REAL(KIND=JPRB) :: SFLUXREF(JPG)  ,RAYL(JPG)
REAL(KIND=JPRB) :: GIVFAC
INTEGER(KIND=JPIM) :: LAYREFFR

REAL(KIND=JPRB) :: KAC(5,13,NG23)   ,ABSA(65,NG23)
REAL(KIND=JPRB) :: SELFREFC(10,NG23),FORREFC(3,NG23)
REAL(KIND=JPRB) :: SFLUXREFC(NG23)  ,RAYLC(NG23)

!EQUIVALENCE (KA(1,1,1),ABSA(1,1))
EQUIVALENCE (KAC(1,1,1),ABSA(1,1))

!     -----------------------------------------------------------------
!        * E.C.M.W.F. PHYSICS PACKAGE ** RRTM SW RADIATION **

!     J.-J. MORCRETTE       E.C.M.W.F.      02/10/29
!     M. J. IACONO          AER             12/09/03

!  NAME     TYPE     PURPOSE
!  ----   : ----   : ---------------------------------------------------
! LAYREFFR: INTEGER
! KA      : REAL     absorption coefficient of major absorber
! SELFREF : REAL     self brodening coefficient for water vapour
! FORREF  : REAL     foreign broadening coefficient for water vapour
! SFLUXREF: REAL     Incident solar radiation in the spectral interval
! RAYL    : REAL     Rayleigh scattering parameter
! GIVFAC  : REAL     adjustment factor
! LAYREFFR: INTEGER  reference level for the transition
! KAC     : REAL     Reduced g-point array for KA
! SELFREFC: REAL     Reduced g-point array for SELFREF
! FORREFC : REAL     Reduced g-point array for FORREF
!SFLUXREFC: REAL     Reduced g-point array for SFLUXREF
! RAYLC   : REAL     Reduced g-point array for RAYL
!     -----------------------------------------------------------------
END MODULE YOESRTA23

