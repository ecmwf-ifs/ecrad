!***************************************************************************
SUBROUTINE RRTM_INIT_140GP(CDIRECTORY)
!***************************************************************************
!     Reformatted for F90 by JJMorcrette, ECMWF, 980714

!     JJMorcrette 20110613 flexible number of g-points

! Parameters
USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE YOMHOOK   ,ONLY : LHOOK, DR_HOOK, JPHOOK

USE PARRRTM  , ONLY : JPBAND   ,JPG
USE YOERRTM  , ONLY : JPGPT
USE YOERRTWN , ONLY : NG
USE YOERRTFTR, ONLY : NGC      ,NGN      ,NGM     , WT
! Output
USE YOERRTBG2, ONLY : CORR1    ,CORR2
USE YOERRTRWT, ONLY : FREFA    ,FREFB    ,FREFADF  ,FREFBDF   ,RWGT
!USE YOMLUN   , ONLY : NULOUT

USE YOMMP0_IFSAUX, ONLY : NPROC, MYPROC
USE MPL_MODULE,     ONLY : MPL_NPROC, MPL_MYRANK

IMPLICIT NONE

CHARACTER(LEN=*), INTENT(IN) :: CDIRECTORY

REAL(KIND=JPRB) :: ZWTSM(JPG)

INTEGER(KIND=JPIM) :: I, IBND, IG, IGC, IGCSM, IND, IPR, IPRSM, IPT

REAL(KIND=JPRB) :: ZFP, ZRTFP, ZWTSUM
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!#include "surrtmcf.intfb.h"
#include "surrtftr.intfb.h"

#include "rrtm_kgb1.intfb.h"
#include "rrtm_kgb10.intfb.h"
#include "rrtm_kgb11.intfb.h"
#include "rrtm_kgb12.intfb.h"
#include "rrtm_kgb13.intfb.h"
#include "rrtm_kgb14.intfb.h"
#include "rrtm_kgb15.intfb.h"
#include "rrtm_kgb16.intfb.h"
#include "rrtm_kgb2.intfb.h"
#include "rrtm_kgb3.intfb.h"
#include "rrtm_kgb4.intfb.h"
#include "rrtm_kgb5.intfb.h"
#include "rrtm_kgb6.intfb.h"
#include "rrtm_kgb7.intfb.h"
#include "rrtm_kgb8.intfb.h"
#include "rrtm_kgb9.intfb.h"

#include "rrtm_cmbgb1.intfb.h"
#include "rrtm_cmbgb10.intfb.h"
#include "rrtm_cmbgb11.intfb.h"
#include "rrtm_cmbgb12.intfb.h"
#include "rrtm_cmbgb13.intfb.h"
#include "rrtm_cmbgb14.intfb.h"
#include "rrtm_cmbgb15.intfb.h"
#include "rrtm_cmbgb16.intfb.h"
#include "rrtm_cmbgb2.intfb.h"
#include "rrtm_cmbgb3.intfb.h"
#include "rrtm_cmbgb4.intfb.h"
#include "rrtm_cmbgb5.intfb.h"
#include "rrtm_cmbgb6.intfb.h"
#include "rrtm_cmbgb7.intfb.h"
#include "rrtm_cmbgb8.intfb.h"
#include "rrtm_cmbgb9.intfb.h"

IF (LHOOK) CALL DR_HOOK('RRTM_INIT_140GP',0,ZHOOK_HANDLE)

! Initialize MPI rank data in YOMMP0_IFSAUX
NPROC=MPL_NPROC()
MYPROC=MPL_MYRANK()

!CALL SURRTMCF
CALL SURRTFTR

! Read the absorption-related coefficients over the 16 x 16 g-points

CALL RRTM_KGB1(CDIRECTORY)
CALL RRTM_KGB2
CALL RRTM_KGB3
CALL RRTM_KGB4
CALL RRTM_KGB5
CALL RRTM_KGB6
CALL RRTM_KGB7
CALL RRTM_KGB8
CALL RRTM_KGB9
CALL RRTM_KGB10
CALL RRTM_KGB11
CALL RRTM_KGB12
CALL RRTM_KGB13
CALL RRTM_KGB14
CALL RRTM_KGB15
CALL RRTM_KGB16

!  Calculate lookup tables for functions needed in routine TAUMOL (TAUGB2)

CORR1(0) = 1.0_JPRB
CORR1(200) = 1.0_JPRB
CORR2(0) = 1.0_JPRB
CORR2(200) = 1.0_JPRB
DO I = 1,199
  ZFP = 0.005_JPRB*REAL(I)
  ZRTFP = SQRT(ZFP)
  CORR1(I) = ZRTFP/ZFP
  CORR2(I) = (1.0_JPRB-ZRTFP)/(1.0_JPRB-ZFP)
ENDDO

!  Perform g-point reduction from 16 per band (256 total points) to
!  a band dependant number (140 total points) for all absorption
!  coefficient input data and Planck fraction input data.
!  Compute relative weighting for new g-point combinations.

IGCSM = 0
DO IBND = 1,JPBAND
  IPRSM = 0
  IF (NGC(IBND) < 16) THEN
    DO IGC = 1,NGC(IBND)
      IGCSM = IGCSM + 1
      ZWTSUM = 0.0_JPRB
      DO IPR = 1, NGN(IGCSM)
        IPRSM = IPRSM + 1
        ZWTSUM = ZWTSUM + WT(IPRSM)
      ENDDO
      ZWTSM(IGC) = ZWTSUM
    ENDDO

    DO IG = 1,NG(IBND)
      IND = (IBND-1)*16 + IG
      RWGT(IND) = WT(IG)/ZWTSM(NGM(IND))
    ENDDO
  ELSE
    DO IG = 1,NG(IBND)
      IGCSM = IGCSM + 1
      IND = (IBND-1)*16 + IG
      RWGT(IND) = 1.0_JPRB
    ENDDO
  ENDIF
ENDDO

!  Initialize arrays for combined Planck fraction data.

DO IPT = 1,13
  DO IPR = 1, JPGPT
    FREFA(IPR,IPT) = 0.0_JPRB
    FREFADF(IPR,IPT) = 0.0_JPRB
  ENDDO
ENDDO
DO IPT = 1,6
  DO IPR = 1, JPGPT
    FREFB(IPR,IPT) = 0.0_JPRB
    FREFBDF(IPR,IPT) = 0.0_JPRB
  ENDDO
ENDDO

!  Reduce g-points for relevant data in each LW spectral band.

CALL RRTM_CMBGB1
CALL RRTM_CMBGB2
CALL RRTM_CMBGB3
CALL RRTM_CMBGB4
CALL RRTM_CMBGB5
CALL RRTM_CMBGB6
CALL RRTM_CMBGB7
CALL RRTM_CMBGB8
CALL RRTM_CMBGB9
CALL RRTM_CMBGB10
CALL RRTM_CMBGB11
CALL RRTM_CMBGB12
CALL RRTM_CMBGB13
CALL RRTM_CMBGB14
CALL RRTM_CMBGB15
CALL RRTM_CMBGB16

IF (LHOOK) CALL DR_HOOK('RRTM_INIT_140GP',1,ZHOOK_HANDLE)
END SUBROUTINE RRTM_INIT_140GP
