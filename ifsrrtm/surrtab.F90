SUBROUTINE SURRTAB

!     -----------------------------------------------------------------
!        * E.C.M.W.F. PHYSICS PACKAGE ** AER'S RRTM LW RADIATION **

!     J.-J. MORCRETTE       E.C.M.W.F.      98/07/14

!     -----------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE YOMHOOK   ,ONLY : LHOOK, DR_HOOK, JPHOOK

USE YOERRTAB , ONLY : TRANS, BPADE

IMPLICIT NONE

INTEGER(KIND=JPIM) :: ITR

REAL(KIND=JPRB) :: ZTAU, ZTFN
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

IF (LHOOK) CALL DR_HOOK('SURRTAB',0,ZHOOK_HANDLE)
BPADE=1.0_JPRB/0.278_JPRB
TRANS(0)   =1.0_JPRB
TRANS(5000)=0.0_JPRB
DO ITR=1,4999
  ZTFN=REAL(ITR)/5000._JPRB
  ZTAU=BPADE*ZTFN/(1.0_JPRB-ZTFN)
  TRANS(ITR)=EXP(-ZTAU)
ENDDO

!     -----------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('SURRTAB',1,ZHOOK_HANDLE)
END SUBROUTINE SURRTAB
