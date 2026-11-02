MODULE YOERRTO1

USE PARKIND1  ,ONLY : JPIM     ,JPRB, JPRD

IMPLICIT NONE

PUBLIC

SAVE

!     -----------------------------------------------------------------
!*    ** *YOERRTO1* - RRTM ORIGINAL COEFFICIENTS FOR INTERVAL 1
!     BAND 1:  10-250 cm-1 (low - H2O; high - H2O)
!     ABozzo may 2013 update to rrtmg v4.85
!     band 1:  10-350 cm-1
!     -----------------------------------------------------------------

INTEGER(KIND=JPIM), PARAMETER :: NO1  = 16

REAL(KIND=JPRB) :: FRACREFAO(NO1)  , FRACREFBO(NO1)
REAL(KIND=JPRB) :: KAO(NO1,5,13)
REAL(KIND=JPRB) :: KBO(NO1,5,13:59)
REAL(KIND=JPRD) :: KAO_D(NO1,5,13)
REAL(KIND=JPRD) :: KBO_D(NO1,5,13:59)
REAL(KIND=JPRB) :: KAO_MN2(NO1,19) , KBO_MN2(NO1,19)
REAL(KIND=JPRB) :: SELFREFO(NO1,10), FORREFO(NO1,4)


!     -----------------------------------------------------------------
!        * E.C.M.W.F. PHYSICS PACKAGE ** RRTM LW RADIATION **

!     J.-J. MORCRETTE       E.C.M.W.F.      98/07/14

!  NAME     TYPE     PURPOSE
!  ----   : ----   : ---------------------------------------------------
!FRACREFAO: REAL    
!FRACREFBO: REAL
! FORREFO : REAL
! KAO     : REAL     
! KBO     : REAL     
! SELFREFO: REAL     
!     -----------------------------------------------------------------
END MODULE YOERRTO1
