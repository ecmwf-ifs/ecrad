MODULE YOESRTA28

USE PARKIND1  ,ONLY : JPIM     ,JPRB,JPRD

IMPLICIT NONE

SAVE

!     -----------------------------------------------------------------
!*    ** *YOESRTA28* - SRTM COEFFICIENTS FOR INTERVAL 28
!     BAND 28: 38000-50000 cm-1 (low - O3, O2; high - O3, O2)
!     -----------------------------------------------------------------

INTEGER(KIND=JPIM), PARAMETER :: JPG = 16, NG28 = 16

REAL(KIND=JPRB) :: KA(9,5,13,JPG)   
REAL(KIND=JPRB) :: KB(5,5,13:59,JPG)
REAL(KIND=JPRD) :: KA_D(9,5,13,JPG)   
REAL(KIND=JPRD) :: KB_D(5,5,13:59,JPG)
REAL(KIND=JPRB) :: SFLUXREF(JPG,5)
REAL(KIND=JPRB) :: RAYL              ,STRRAT
INTEGER(KIND=JPIM) :: LAYREFFR

REAL(KIND=JPRB) :: KAC(9,5,13,NG28)   ,ABSA(585,NG28)
REAL(KIND=JPRB) :: KBC(5,5,13:59,NG28),ABSB(1175,NG28)
REAL(KIND=JPRB) :: SFLUXREFC(NG28,5)

!EQUIVALENCE (KA(1,1,1,1),ABSA(1,1)), (KB(1,1,13,1),ABSB(1,1))
EQUIVALENCE (KAC(1,1,1,1),ABSA(1,1)), (KBC(1,1,13,1),ABSB(1,1))

!     -----------------------------------------------------------------
!        * E.C.M.W.F. PHYSICS PACKAGE ** RRTM SW RADIATION **

!     J.-J. MORCRETTE       E.C.M.W.F.      02/10/29
!     M. J. IACONO          AER             12/09/03

!  NAME     TYPE     PURPOSE
!  ----   : ----   : ---------------------------------------------------
! KA      : REAL     absorption coefficient of major absorber
! KB      : REAL     absorption coefficient of secondary absorber
! SFLUXREF: REAL     Incident solar radiation in the spectral interval
! RAYL    : REAL     Rayleigh scattering parameter
! STRRAT  : REAL     weighting factor for the transition between tropospheric 
!                    and stratospheric computations
! LAYREFFR: INTEGER  reference level for the transition
! KAC     : REAL     Reduced g-point array for KA
! KBC     : REAL     Reduced g-point array for KB
!SFLUXREFC: REAL     Reduced g-point array for SFLUXREF
!     -----------------------------------------------------------------
END MODULE YOESRTA28

