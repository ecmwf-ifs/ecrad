! This is taken from yomlun_ifsaux in the IFS

MODULE YOMLUN_IFSAUX

USE PARKIND1  ,ONLY : JPIM

IMPLICIT NONE

SAVE

!     ------------------------------------------------------------------

!*    Logical units used by code

!     NULOUT :   output unit
!     NULERR :   unit number for comparison with reference run

INTEGER(KIND=JPIM) :: NULOUT = 6
INTEGER(KIND=JPIM) :: NULERR = 0

!     ------------------------------------------------------------------
END MODULE YOMLUN_IFSAUX
