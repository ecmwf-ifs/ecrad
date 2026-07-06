MODULE YOESRTA24

USE PARKIND1  ,ONLY : JPIM     ,JPRB,JPRD

IMPLICIT NONE

PUBLIC

SAVE

!     -----------------------------------------------------------------
!*    ** *YOESRTA24* - SRTM COEFFICIENTS FOR INTERVAL 24
!     BAND 24: 12850-16000 cm-1 (low - H2O,O2; high - O2)
!     -----------------------------------------------------------------

INTEGER(KIND=JPIM), PARAMETER :: JPG = 16, NG24 = 16

REAL(KIND=JPRB) :: KA(JPG,9,5,13) 
REAL(KIND=JPRB) :: KB(JPG,5,13:59)
REAL(KIND=JPRD) :: KA_D(JPG,9,5,13) 
REAL(KIND=JPRD) :: KB_D(JPG,5,13:59)
REAL(KIND=JPRB) :: SELFREF(JPG,10),FORREF(JPG,3)
REAL(KIND=JPRB) :: SFLUXREF(JPG,9)
REAL(KIND=JPRB) :: ABSO3A(JPG), ABSO3B(JPG), RAYLA(JPG,9), RAYLB(JPG)
REAL(KIND=JPRB) :: STRRAT
INTEGER(KIND=JPIM) :: LAYREFFR

REAL(KIND=JPRB) :: KAC(NG24,9,5,13) ,ABSA(NG24,585)
REAL(KIND=JPRB) :: KBC(NG24,5,13:59),ABSB(NG24,235)
REAL(KIND=JPRB) :: SELFREFC(NG24,10),FORREFC(NG24,3)
REAL(KIND=JPRB) :: SFLUXREFC(NG24,9)
REAL(KIND=JPRB) :: ABSO3AC(NG24), ABSO3BC(NG24), RAYLAC(NG24,9), RAYLBC(NG24)

!EQUIVALENCE (KA(1,1,1,1),ABSA(1,1)), (KB(1,13,1),ABSB(1,1))
EQUIVALENCE (KAC(1,1,1,1),ABSA(1,1)), (KBC(1,1,13),ABSB(1,1))

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
! ABSO3A  : REAL     O3 absorption coefficient in first part of band
! ABSO3B  : REAL     O3 absorption coefficient in second part of band
! RAYLA   : REAL     Rayleigh scattering parameter in first part of band
! RAYLB   : REAL     Rayleigh scattering parameter in second part of band   
! STRRAT  : REAL     weighting factor for the transition between tropospheric 
!                    and stratospheric computations
! LAYREFFR: INTEGER  reference level for the transition
! KAC     : REAL     Reduced g-point array for KA
! KBC     : REAL     Reduced g-point array for KB
! SELFREFC: REAL     Reduced g-point array for SELFREF
! FORREFC : REAL     Reduced g-point array for FORREF
!SFLUXREFC: REAL     Reduced g-point array for SFLUXREF
! ABSO3AC : REAL     Reduced g-point array for ABSO3A
! ABSO3BC : REAL     Reduced g-point array for ABSO3B
! RAYLAC  : REAL     Reduced g-point array for RAYLA
! RAYLBC  : REAL     Reduced g-point array for RAYLB
!     -----------------------------------------------------------------
END MODULE YOESRTA24

