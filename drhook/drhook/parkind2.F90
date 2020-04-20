MODULE PARKIND2
!
!     *** Define huge kinds for strong typing ***
!
IMPLICIT NONE
SAVE
!
!     Integer Kinds
!     -------------
!
INTEGER, PARAMETER :: JPIH = SELECTED_INT_KIND(18)
!
!     Real Kinds
!     ----------
!
#ifdef REALHUGE
INTEGER, PARAMETER :: JPRH = SELECTED_REAL_KIND(31,291)
#else
INTEGER, PARAMETER :: JPRH = SELECTED_REAL_KIND(13,300)
#endif
!
END MODULE PARKIND2
