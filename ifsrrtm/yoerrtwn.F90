MODULE YOERRTWN

USE PARKIND1  ,ONLY : JPIM     ,JPRB

IMPLICIT NONE

SAVE

!    -------------------------------------------------------------------

INTEGER(KIND=JPIM) , DIMENSION(16) :: NG
INTEGER(KIND=JPIM) , DIMENSION(16) :: NSPA
INTEGER(KIND=JPIM) , DIMENSION(16) :: NSPB

REAL(KIND=JPRB) , DIMENSION(16) :: DELWAVE

REAL(KIND=JPRB) , DIMENSION(181,16) :: TOTPLNK

!     -----------------------------------------------------------------
!        * E.C.M.W.F. PHYSICS PACKAGE ** RRTM LW RADIATION **

!     J.-J. MORCRETTE       E.C.M.W.F.      98/01/15

!  NAME     TYPE     PURPOSE
!  ----   : ----    : -------
!  NG     : INTEGER : Number of k-coefficients in spectral intervals
!  NSPA   : INTEGER :
!  NSPB   : INTEGER :
! WAVENUM1: REAL    : Lower wavenumber spectral limit
! WAVENUM2: REAL    : Higher wavenumber spectral limit
! DELWAVE : REAL    : Spectral interval width
! TOTPLNK : REAL    :
!     -----------------------------------------------------------------
END MODULE YOERRTWN
