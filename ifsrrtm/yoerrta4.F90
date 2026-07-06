MODULE YOERRTA4

USE PARKIND1  ,ONLY : JPIM     ,JPRB

IMPLICIT NONE

PUBLIC

SAVE

!     -----------------------------------------------------------------
!*    ** *YOERRTA4* - RRTM COEFFICIENTS FOR INTERVAL 5
!     BAND 4:  630-700 cm-1 (low - H2O,CO2; high - O3,CO2)
!     201306 ABozzo updated to rrtmg v4.85
!     -----------------------------------------------------------------

INTEGER(KIND=JPIM), PARAMETER :: NG4  = 14

REAL(KIND=JPRB) :: FRACREFA(NG4,9)  ,FRACREFB(NG4,5)
REAL(KIND=JPRB) :: KA(NG4,9,5,13)   ,ABSA(NG4,585)
REAL(KIND=JPRB) :: KB(NG4,5,5,13:59),ABSB(NG4,1175)
REAL(KIND=JPRB) :: SELFREF(NG4,10)  ,FORREF(NG4,4)

EQUIVALENCE (KA(1,1,1,1),ABSA(1,1)),(KB(1,1,1,13),ABSB(1,1))

!     ------------------------------------------------------------------
!        * E.C.M.W.F. PHYSICS PACKAGE ** RRTM LW RADIATION **

!     J.-J. MORCRETTE       E.C.M.W.F.      98/07/14

!  NAME     TYPE     PURPOSE
!  ----   : ----   : ---------------------------------------------------
! ABSA    : REAL     absorption coefficient of major absorber for M reference tropospheric 
!                    pressures and N reference tropospheric temperatures 
! ABSB    : REAL     absorption coefficient of secondary absorber for M reference stratospheric
!                    pressures and N reference stratospheric temperatures 
! FRACREFA: REAL     distance from r and T reference tabulated points (troposphere)
! FRACREFB: REAL     distance from r and T reference tabulated points (stratosphere)
! KA      : REAL     absorption coefficient of major absorber (equiv. to ABSA)   
! KB      : REAL     absorption coefficient of secondary absorber (equiv. to ABSB)   
! SELFREF : REAL     self broadening coefficient for water vapour
! STRRATn : REAL     weighting factors for the transition between tropospheric 
!                    and stratospheric computations
!     -----------------------------------------------------------------
END MODULE YOERRTA4
