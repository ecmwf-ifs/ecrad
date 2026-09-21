MODULE YOESRTA16

USE PARKIND1  ,ONLY : JPIM     ,JPRB,JPRD

IMPLICIT NONE

PUBLIC

SAVE

!     -----------------------------------------------------------------
!*    ** *YOESRTA16* - SRTM COEFFICIENTS FOR INTERVAL 16
!     BAND 16:  2600-3250 cm-1 (low - H2O,CH4; high - CH4)
!     -----------------------------------------------------------------

INTEGER(KIND=JPIM), PARAMETER :: JPG=16, NG16 = 16, NGS15 = 0

REAL(KIND=JPRB) :: KA(JPG,9,5,13) 
REAL(KIND=JPRB) :: KB(JPG,5,13:59)
REAL(KIND=JPRD) :: KA_D(JPG,9,5,13) 
REAL(KIND=JPRD) :: KB_D(JPG,5,13:59)
REAL(KIND=JPRB) :: SELFREF(JPG,10),FORREF(JPG,3)
REAL(KIND=JPRB) :: SFLUXREF(JPG)
REAL(KIND=JPRB) :: RAYL            ,STRRAT1
INTEGER(KIND=JPIM) :: LAYREFFR

REAL(KIND=JPRB) :: KAC(NG16,9,5,13),ABSA(NG16,585)
REAL(KIND=JPRB) :: KBC(NG16,5,13:59),ABSB(NG16,235)
REAL(KIND=JPRB) :: SELFREFC(NG16,10),FORREFC(NG16,3)
REAL(KIND=JPRB) :: SFLUXREFC(NG16)

!EQUIVALENCE (KA(1,1,1,1),ABSA(1,1)), (KB(1,13,1),ABSB(1,1))
EQUIVALENCE (KAC(1,1,1,1),ABSA(1,1)), (KBC(1,1,13),ABSB(1,1))

!     -----------------------------------------------------------------
!        * E.C.M.W.F. PHYSICS PACKAGE ** RRTM SW RADIATION **

!     J.-J. MORCRETTE       E.C.M.W.F.      02/10/29
!     M. J. IACONO          AER             12/09/03

!  NAME     TYPE     PURPOSE
!  ----   : ----   : ---------------------------------------------------
! KA      : REAL     absorption coefficient of major absorber
! KB      : REAL     absorption coefficient of secondary absorber
! SELFREF : REAL     self brodening coefficient for water vapour
! FORREF  : REAL     foreign broadening coefficient for water vapour
! SFLUXREF: REAL     Incident solar radiation in the spectral interval
! RAYL    : REAL     Rayleigh scattering parameter
! STRRAT1 : REAL     weighting factor for the transition between tropospheric 
!                    and stratospheric computations
! LAYREFFR: INTEGER  reference level for the transition
! KAC     : REAL     Reduced g-point array for KA
! KBC     : REAL     Reduced g-point array for KB
! SELFREFC: REAL     Reduced g-point array for SELFREF
! FORREFC : REAL     Reduced g-point array for FORREF
!SFLUXREFC: REAL     Reduced g-point array for SFLUXREF
!     -----------------------------------------------------------------
END MODULE YOESRTA16

