MODULE YOESRTA27

USE PARKIND1  ,ONLY : JPIM     ,JPRB,JPRD

IMPLICIT NONE

PUBLIC

SAVE

!     -----------------------------------------------------------------
!*    ** *YOESRTA27* - SRTM COEFFICIENTS FOR INTERVAL 27
!     BAND 27: 29000-38000 cm-1 (low - O3; high - O3)
!     -----------------------------------------------------------------

INTEGER(KIND=JPIM), PARAMETER :: JPG = 16, NG27 = 16

REAL(KIND=JPRB) :: KA(JPG,5,13)   
REAL(KIND=JPRB) :: KB(JPG,5,13:59)
REAL(KIND=JPRD) :: KA_D(JPG,5,13)   
REAL(KIND=JPRD) :: KB_D(JPG,5,13:59)
REAL(KIND=JPRB) :: SFLUXREF(JPG)  ,RAYL(JPG)
REAL(KIND=JPRB) :: SCALEKUR
INTEGER(KIND=JPIM) :: LAYREFFR

REAL(KIND=JPRB) :: KAC(NG27,5,13)   ,ABSA(NG27,65)
REAL(KIND=JPRB) :: KBC(NG27,5,13:59),ABSB(NG27,235)
REAL(KIND=JPRB) :: SFLUXREFC(NG27)  ,RAYLC(NG27)

!EQUIVALENCE (KA(1,1,1),ABSA(1,1)), (KB(1,13,1),ABSB(1,1))
EQUIVALENCE (KAC(1,1,1),ABSA(1,1)), (KBC(1,1,13),ABSB(1,1))

!     -----------------------------------------------------------------
!        * E.C.M.W.F. PHYSICS PACKAGE ** RRTM SW RADIATION **

!     J.-J. MORCRETTE       E.C.M.W.F.      02/10/29
!     M. J. IACONO          AER             12/09/03

!  NAME     TYPE     PURPOSE
!  ----   : ----   : ---------------------------------------------------
! FRACREFA: REAL
! KA      : REAL     absorption coefficient of major absorber
! KB      : REAL     absorption coefficient of secondary absorber
! SFLUXREF: REAL     Incident solar radiation in the spectral interval
! RAYL    : REAL     Rayleigh scattering parameter
! SCALEKUR: REAL     weighting factor to account for Kurucz's correction to solar spectrum 
! LAYREFFR: INTEGER  reference level for the transition
! KAC     : REAL     Reduced g-point array for KA
! KBC     : REAL     Reduced g-point array for KB
!SFLUXREFC: REAL     Reduced g-point array for SFLUXREF
! RAYLC   : REAL     Reduced g-point array for RAYL
!     -----------------------------------------------------------------
END MODULE YOESRTA27

