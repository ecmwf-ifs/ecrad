MODULE YOERRTO5

USE PARKIND1  ,ONLY : JPIM     ,JPRB,JPRD

IMPLICIT NONE

PUBLIC

SAVE

!     -----------------------------------------------------------------
!*    ** *YOERRTO5* - RRTM ORIGINAL COEFFICIENTS FOR INTERVAL 5
!     BAND 5:  700-820 cm-1 (low - H2O,CO2; high - O3,CO2)
!     ABozzo 201306 updated to rrtmg v4.85
!     -----------------------------------------------------------------

INTEGER(KIND=JPIM), PARAMETER :: NO5  = 16

REAL(KIND=JPRB) :: FRACREFAO(NO5,9) ,FRACREFBO(NO5,5)

REAL(KIND=JPRB) , DIMENSION(NO5) :: CCL4O

REAL(KIND=JPRB) :: KAO_MO3(NO5,9,19)
REAL(KIND=JPRB) :: KAO(NO5,9,5,13)
REAL(KIND=JPRB) :: KBO(NO5,5,5,13:59)
REAL(KIND=JPRD) :: KAO_D(NO5,9,5,13)
REAL(KIND=JPRD) :: KBO_D(NO5,5,5,13:59)
REAL(KIND=JPRB) :: SELFREFO(NO5,10)
REAL(KIND=JPRB) :: FORREFO(NO5,4)


!     -----------------------------------------------------------------
!        * E.C.M.W.F. PHYSICS PACKAGE ** RRTM LW RADIATION **

!     J.-J. MORCRETTE       E.C.M.W.F.      98/07/14

!  NAME     TYPE     PURPOSE
!  ----   : ----   : ---------------------------------------------------
! FRACREFA: REAL    
! FRACREFB: REAL
! KA      : REAL     
! KB      : REAL     
! SELFREF : REAL     
!     -----------------------------------------------------------------
END MODULE YOERRTO5
