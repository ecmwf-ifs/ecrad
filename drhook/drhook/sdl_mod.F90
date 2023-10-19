! (C) Copyright 2014- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
!
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.

MODULE SDL_MOD

!    Interface between user applications and system-dependent intrinsic
!    routines, provided by the computer vendors.

!    All routines which wish to call these routines must contain:
!    USE SDL_MOD

! Author :
! ------
!   11-Apr-2005 R. El Khatib  *METEO-FRANCE*
!   26-Apr-2006 S.T.Saarinen  Dr.Hook trace, calls to EC_RAISE, Intel/ifort traceback

USE PARKIND1  ,ONLY : JPIM  ,JPRB
!USE YOMHOOK   ,ONLY : LHOOK ,DR_HOOK
USE YOMHOOK   ,ONLY : LHOOK 
USE OML_MOD   ,ONLY : OML_MY_THREAD
USE MPL_MPIF  ,ONLY : MPI_COMM_WORLD
IMPLICIT NONE

SAVE

PRIVATE

INTEGER, PARAMETER :: SIGABRT = 6 ! Hardcoded

PUBLIC :: SDL_SRLABORT, SDL_DISABORT, SDL_TRACEBACK

CONTAINS

!-----------------------------------------------------------------------------
SUBROUTINE SDL_TRACEBACK(KTID)

! Purpose :
! -------
!   Traceback

!   KTID : thread 

INTEGER(KIND=JPIM), INTENT(IN), OPTIONAL :: KTID
INTEGER(KIND=JPIM) ITID, IPRINT_OPTION, ILEVEL
CHARACTER(LEN=80) :: CLTRBK
#ifdef NECSX
CHARACTER(LEN=*), PARAMETER :: CLNECMSG = '*** Calling NEC traceback ***'
#endif
IF (PRESENT(KTID)) THEN
  ITID = KTID
ELSE
  ITID = OML_MY_THREAD()
ENDIF
IF (LHOOK) THEN
  IPRINT_OPTION = 2
  ILEVEL = 0
  CALL C_DRHOOK_PRINT(0, ITID, IPRINT_OPTION, ILEVEL) ! from drhook.c
ENDIF
#if defined(VPP)
  CALL ERRTRA
  IF (PRESENT(KTID)) CALL SLEEP(28)
#elif defined(RS6K)
  WRITE(0,*)'SDL_TRACEBACK: Calling XL_TRBK, THRD = ',ITID
  CALL XL__TRBK()
  WRITE(0,*)'SDL_TRACEBACK: Done XL_TRBK, THRD = ',ITID
#elif defined(__INTEL_COMPILER)

  CALL GET_ENVIRONMENT_VARIABLE("EC_LINUX_TRBK",CLTRBK)
  IF (CLTRBK=='1') THEN
    WRITE(0,*)'SDL_TRACEBACK: Calling LINUX_TRBK, THRD = ',ITID
    CALL LINUX_TRBK() ! See ifsaux/utilities/linuxtrbk.c
    WRITE(0,*)'SDL_TRACEBACK: Done LINUX_TRBK, THRD = ',ITID
  ELSE
    WRITE(0,*)'SDL_TRACEBACK: Calling INTEL_TRBK, THRD = ',ITID
    CALL INTEL_TRBK() ! See ifsaux/utilities/gentrbk.F90
    WRITE(0,*)'SDL_TRACEBACK: Done INTEL_TRBK, THRD = ',ITID
  ENDIF
#elif defined(__NEC__)
  ! A traceback using gdb-debugger, if available AND 
  ! activated via 'export GDBDEBUGGER=1'
  WRITE(0,*)'SDL_TRACEBACK: Calling GDB_TRBK, THRD = ',ITID
  CALL GDB_TRBK() ! See ifsaux/utilities/linuxtrbk.c
  WRITE(0,*)'SDL_TRACEBACK: Done GDB_TRBK, THRD = ',ITID
#elif defined(LINUX) || defined(SUN4)
  WRITE(0,*)'SDL_TRACEBACK: Calling LINUX_TRBK, THRD = ',ITID
  CALL LINUX_TRBK() ! See ifsaux/utilities/linuxtrbk.c
  WRITE(0,*)'SDL_TRACEBACK: Done LINUX_TRBK, THRD = ',ITID
#elif defined(NECSX)
! MESPUT writes out onto unit 6
  WRITE(6,*)'SDL_TRACEBACK: Calling NEC/MESPUT, THRD = ',ITID
  CALL NECSX_TRBK(CLNECMSG)
  CALL FLUSH(6)
  WRITE(6,*)'SDL_TRACEBACK: Done NEC/MESPUT, THRD = ',ITID
#else
  WRITE(0,*)'SDL_TRACEBACK: No proper traceback implemented.'
  ! A traceback using dbx-debugger, if available AND 
  ! activated via 'export DBXDEBUGGER=1'
  WRITE(0,*)'SDL_TRACEBACK: Calling DBX_TRBK, THRD = ',ITID
  CALL DBX_TRBK() ! See ifsaux/utilities/linuxtrbk.c
  WRITE(0,*)'SDL_TRACEBACK: Done DBX_TRBK, THRD = ',ITID
  ! A traceback using gdb-debugger, if available AND 
  ! activated via 'export GDBDEBUGGER=1'
  WRITE(0,*)'SDL_TRACEBACK: Calling GDB_TRBK, THRD = ',ITID
  CALL GDB_TRBK() ! See ifsaux/utilities/linuxtrbk.c
  WRITE(0,*)'SDL_TRACEBACK: Done GDB_TRBK, THRD = ',ITID
#endif

END SUBROUTINE SDL_TRACEBACK
!-----------------------------------------------------------------------------
SUBROUTINE SDL_SRLABORT

! Purpose :
! -------
!   To abort in serial environment

CALL EC_RAISE(SIGABRT)
STOP 'SDL_SRLABORT'

END SUBROUTINE SDL_SRLABORT
!-----------------------------------------------------------------------------
SUBROUTINE SDL_DISABORT(KCOMM)

! Purpose :
! -------
!   To abort in distributed environment

!   KCOMM : communicator

INTEGER(KIND=JPIM), INTENT(IN) :: KCOMM

INTEGER(KIND=JPIM) :: IRETURN_CODE,IERROR
CHARACTER(LEN=80) :: CLJOBID
CHARACTER(LEN=80) :: CLTRBK

#ifdef VPP

CALL VPP_ABORT()

#else
#if defined(__INTEL_COMPILER)
! Intel compiler seems to hang in MPI_ABORT -- on all but the failing task(s)
! ... when linux trbk is used. REK
CALL GET_ENVIRONMENT_VARIABLE("EC_LINUX_TRBK",CLTRBK)
IF (CLTRBK=='1') THEN
IF (LHOOK) THEN
  CALL GET_ENVIRONMENT_VARIABLE("SLURM_JOBID",CLJOBID)
  IF (CLJOBID /= ' ') THEN
    CALL SYSTEM("set -x; sleep 10; scancel --signal=TERM "//trim(CLJOBID)//" &")
  ENDIF
ENDIF
ENDIF
#endif

IRETURN_CODE=SIGABRT
!CALL MPI_ABORT(KCOMM,IRETURN_CODE,IERROR)
CALL MPI_ABORT(MPI_COMM_WORLD,IRETURN_CODE,IERROR) ! Tracked by the supervisor/process-damager (manager) -- KCOMM /= MPI_COMM_WORLD may hang as sub-communicator

#endif

CALL EC_RAISE(SIGABRT) ! In case ever ends up here
STOP 'SDL_DISABORT'

END SUBROUTINE SDL_DISABORT
!-----------------------------------------------------------------------------

END MODULE SDL_MOD
