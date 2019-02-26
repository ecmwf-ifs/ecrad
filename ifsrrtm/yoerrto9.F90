MODULE YOERRTO9

USE PARKIND1  ,ONLY : JPIM     ,JPRB,JPRD

IMPLICIT NONE

SAVE

!     -----------------------------------------------------------------
!*    ** *YOERRTO9* - RRTM ORIGINAL COEFFICIENTS FOR INTERVAL 9
!     BAND 9:  1180-1390 cm-1 (low - H2O,CH4; high - CH4)
!     ABozzo 201306 updated to rrtmg v4.85
!     -----------------------------------------------------------------

INTEGER(KIND=JPIM), PARAMETER :: NO9  = 16

REAL(KIND=JPRB) :: FRACREFAO(NO9,9)

REAL(KIND=JPRB) , DIMENSION(NO9) :: FRACREFBO
! 48 = 3*NO9      


REAL(KIND=JPRB) :: KAO(9,5,13,NO9)
REAL(KIND=JPRB) :: KBO(5,13:59,NO9)
REAL(KIND=JPRD) :: KAO_D(9,5,13,NO9)
REAL(KIND=JPRD) :: KBO_D(5,13:59,NO9)
REAL(KIND=JPRB) :: KAO_MN2O(9,19,NO9)
REAL(KIND=JPRB) :: KBO_MN2O(19,NO9)
REAL(KIND=JPRB) :: SELFREFO(10,NO9)
REAL(KIND=JPRB) :: FORREFO(4,NO9)
!     -----------------------------------------------------------------
!        * E.C.M.W.F. PHYSICS PACKAGE ** RRTM LW RADIATION **

!     J.-J. MORCRETTE       E.C.M.W.F.      98/07/14

!  NAME     TYPE     PURPOSE
!  ----   : ----   : ---------------------------------------------------
! ABSN2O  : REAL    
! CH4REF  : REAL
! ETAREF  : REAL
! FRACREFA: REAL    
! FRACREFB: REAL
! H2OREF  : REAL
! N2OREF  : REAL
! KA      : REAL     
! KB      : REAL     
! SELFREF : REAL     
!     -----------------------------------------------------------------
END MODULE YOERRTO9
