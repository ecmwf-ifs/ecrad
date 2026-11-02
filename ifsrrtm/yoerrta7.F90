MODULE YOERRTA7

USE PARKIND1  ,ONLY : JPIM     ,JPRB

IMPLICIT NONE

PUBLIC

SAVE

!     -----------------------------------------------------------------
!*    ** *YOERRTA7* - RRTM COEFFICIENTS FOR INTERVAL 7
!     BAND 7:  980-1080 cm-1 (low - H2O,O3; high - O3)
!     ABozzo updated to rrtmg v4.85
!     -----------------------------------------------------------------

INTEGER(KIND=JPIM), PARAMETER :: NG7  = 12

REAL(KIND=JPRB) :: FRACREFA(NG7,9)

REAL(KIND=JPRB) , DIMENSION(NG7) :: FRACREFB
REAL(KIND=JPRB) :: KA(NG7,9,5,13) ,ABSA(NG7,585)
REAL(KIND=JPRB) :: KB(NG7,5,13:59),ABSB(NG7,235)
REAL(KIND=JPRB) :: SELFREF(NG7,10)
REAL(KIND=JPRB) :: KA_MCO2(NG7,9,19)
REAL(KIND=JPRB) :: KB_MCO2(NG7,19)
REAL(KIND=JPRB) :: FORREF(NG7,4)

EQUIVALENCE (KA(1,1,1,1),ABSA(1,1)),(KB(1,1,13),ABSB(1,1))

!     -----------------------------------------------------------------
!        * E.C.M.W.F. PHYSICS PACKAGE *

!     J.-J. MORCRETTE       E.C.M.W.F.      98/07/14

!  NAME     TYPE     PURPOSE
!  ----   : ----   : --------------------------------------------------- 
! ABSA    : REAL     absorption coefficient of major absorber for M reference tropospheric 
!                    pressures and N reference tropospheric temperatures 
! ABSB    : REAL     absorption coefficient of secondary absorber for M reference stratospheric
!                    pressures and N reference stratospheric temperatures 
! ABSCO2  : REAL     absorption coefficient for CO2
! FRACREFA: REAL     distance from r and T reference tabulated points (troposphere)
! FRACREFB: REAL     distance from r and T reference tabulated points (stratosphere)
! KA      : REAL     absorption coefficient of major absorber (equiv. to ABSA)   
! KB      : REAL     absorption coefficient of secondary absorber (equiv. to ABSB)   
! SELFREF : REAL     self broadening coefficient for water vapour
! STRRAT  : REAL     weighting factors for the transition between tropospheric 
!                    and stratospheric computations
!     -----------------------------------------------------------------
END MODULE YOERRTA7
