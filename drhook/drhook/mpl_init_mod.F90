MODULE MPL_INIT_MOD

!**** MPL_INIT - Initialises the Message passing environment

!     Purpose.
!     --------
!     Must be called before any other MPL routine.

!**   Interface.
!     ----------
!        CALL MPL_INIT

!        Input required arguments :
!        -------------------------
!           none

!        Input optional arguments :
!        -------------------------
!           KOUTPUT  -  Level of printing for MPL routines
!                       =0: none
!                       =1: intermediate (default)
!                       =2: full trace
!           KUNIT    -  Fortran Unit to receive printed trace
!           LDINFO   -  = .TRUE.  : Print informative msgs from MPL_INIT (default) 
!                       = .FALSE. : Do not print
!           LDENV    -  = .TRUE.  : Propagate environment variables across participating tasks (default)
!                       = .FALSE. : Do not propagate
!

!        Output required arguments :
!        -------------------------
!           none

!        Output optional arguments :
!        -------------------------
!           KERROR   -  return error code.     If not supplied, 
!                       MPL_INIT aborts when an error is detected.
!           KPROCS   -  Number of processes which have been initialised
!                       in the MPI_COMM_WORLD communicator
!     Author.
!     -------
!        D.Dent, M.Hamrud     ECMWF

!     Modifications.
!     --------------
!        Original: 2000-09-01
!        R. El Khatib  14-May-2007 Do not propagate environment if NECSX
!        S. Saarinen 04-Oct-2009 Reduced output & redefined MPL_COMM_OML(1)
!        P. Marguinaud 01-Jan-2011 Add LDENV argument
!        R. El Khatib  24-May-2011 Make MPI2 the default expectation.
!        P. Towers     3-Jul-2014 Add call to ec_cray_meminfo
!     ------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE OML_MOD, ONLY : OML_INIT, OML_MAX_THREADS
USE MPL_MPIF
USE MPL_DATA_MODULE
USE MPL_MESSAGE_MOD
USE MPL_BUFFER_METHOD_MOD
USE MPL_TOUR_TABLE_MOD
USE MPL_LOCOMM_CREATE_MOD
USE MPL_ARG_MOD

IMPLICIT NONE

PUBLIC MPL_INIT

PRIVATE

CONTAINS 

SUBROUTINE MPL_INIT(KOUTPUT,KUNIT,KERROR,KPROCS,LDINFO,LDENV)

#ifdef USE_8_BYTE_WORDS
  USE MPI4TO8, ONLY : &
    MPI_INITIALIZED => MPI_INITIALIZED8, MPI_INIT => MPI_INIT8, &
    MPI_COMM_SIZE => MPI_COMM_SIZE8,     MPI_COMM_RANK => MPI_COMM_RANK8, &
    MPI_BCAST => MPI_BCAST8
#endif



INTEGER(KIND=JPIM),INTENT(IN),OPTIONAL :: KOUTPUT,KUNIT
INTEGER(KIND=JPIM),INTENT(OUT),OPTIONAL :: KERROR,KPROCS
LOGICAL,INTENT(IN),OPTIONAL :: LDINFO,LDENV
INTEGER(KIND=JPIM) :: IERROR,IP,ICOMM,IRANK,JNODE,JROC,ISTA
INTEGER(KIND=JPIM) :: IMAX_THREADS, IRET, IROOT, INUM(2), ICOUNT
INTEGER(KIND=JPIM) :: IREQUIRED,IPROVIDED
INTEGER(KIND=JPIM) :: IWORLD_RANK, IWORLD_SIZE
INTEGER(KIND=JPIM) :: IME
LOGICAL            :: LLABORT=.TRUE., LLINFO
LOGICAL            :: LLINIT
LOGICAL            :: LLENV
CHARACTER(LEN=12)  :: CL_MBX_SIZE
CHARACTER(LEN=12)  :: CL_ARCH
CHARACTER(LEN=12)  :: CL_TASKSPERNODE
CHARACTER(LEN=1024) :: CLENV
CHARACTER(LEN=20)  :: CL_METHOD,CL_HOST
CHARACTER(LEN=1)   :: CL_SET

IF(PRESENT(KOUTPUT)) THEN
  MPL_OUTPUT=MAX(0,KOUTPUT)
ELSE
  MPL_OUTPUT=1
ENDIF

IF(PRESENT(KUNIT)) THEN
  MPL_UNIT=MAX(0,KUNIT)
ELSE
  MPL_UNIT=6
ENDIF

IF(PRESENT(LDINFO)) THEN
  LLINFO = LDINFO
ELSE
  LLINFO = .TRUE.
ENDIF

IF(PRESENT(LDENV)) THEN
  LLENV = LDENV
ELSE
  LLENV = .TRUE.
ENDIF

IF(MPL_NUMPROC /= -1) THEN
!!  We do not want this extra message
!!  CALL MPL_MESSAGE(CDMESSAGE=' MPL_INIT CALLED MULTIPLE TIMES ')
  IF(PRESENT(KERROR)) THEN
    KERROR=0
  ENDIF
  IF(PRESENT(KPROCS)) THEN
    KPROCS=MPL_NUMPROC
  ENDIF
  RETURN
ENDIF

CALL MPI_INITIALIZED(LLINIT, IRET)

IF (.NOT.LLINIT) THEN

  CALL GET_ENVIRONMENT_VARIABLE('ARCH',CL_ARCH)

#ifndef OPS_COMPILE
#ifdef RS6K
  IF(CL_ARCH(1:10)=='ibm_power6')THEN
!   write(0,*)'POWER6: CALLING EC_BIND BEFORE MPI_INIT'
    CALL EC_BIND()
  ENDIF
#endif
#endif


#ifndef MPI1
  IREQUIRED = MPI_THREAD_MULTIPLE
  IPROVIDED = MPI_THREAD_SINGLE
  CALL MPI_INIT_THREAD(IREQUIRED,IPROVIDED,IERROR)
  IF (IERROR /= 0) CALL ABOR1 ('MPL_INIT: MPI_INIT_THREAD FAILED')
  LTHSAFEMPI = (IPROVIDED >= IREQUIRED)
#else
  CALL MPI_INIT(IERROR)
  IF (IERROR /= 0) CALL ABOR1 ('MPL_INIT: MPI_INIT FAILED')
  LTHSAFEMPI = .FALSE.
#endif

CALL MPI_Comm_rank(MPI_COMM_WORLD, IME, IERROR)

! Print out thread safety etc. messages -- must use MPI_Comm_rank since MPL not initialized just yet
IF (IME == 0) THEN
   WRITE(0,'(1X,A,4(1X,I0),1(1X,L1))') &
        & 'MAIN: IREQUIRED, MPI_THREAD_MULTIPLE, MPI_THREAD_SINGLE, IPROVIDED, LTHSAFEMPI =',&
        &        IREQUIRED, MPI_THREAD_MULTIPLE, MPI_THREAD_SINGLE, IPROVIDED, LTHSAFEMPI
ENDIF

#ifndef OPS_COMPILE
#ifdef RS6K
  IF(CL_ARCH(1:10)=='ibm_power4')THEN
!   write(0,*)'POWER5: CALLING EC_BIND AFTER MPI_INIT'
    CALL EC_BIND()
  ENDIF
#endif
#endif

  LINITMPI_VIA_MPL = .TRUE.
!  CALL ec_mpi_atexit() ! ifsaux/support/endian.c: to make sure MPI_FINALIZE gets called

ELSE
  IERROR = 0
ENDIF

IF(PRESENT(KERROR)) THEN
  KERROR=IERROR
ELSE
  IF(IERROR /= 0) THEN
    CALL MPL_MESSAGE(IERROR,CDMESSAGE=' MPL_INIT ERROR ',LDABORT=LLABORT)
  ENDIF
ENDIF

! If LMPLUSERCOMM is not set use MPI_COMM_WORLD
!mps: Sami Saarinen, 29-Nov-2016
! Must be set *AFTER* MPI_INIT*() has ben called (or LLINIT is true)
! Otherwise MPI_COMM_WORLD not defined (at least not in OpenMPI)
IF(LMPLUSERCOMM) THEN
   MPL_COMM = MPLUSERCOMM
ELSE
   MPL_COMM = MPI_COMM_WORLD
ENDIF

CALL MPI_COMM_SIZE(MPL_COMM,MPL_NUMPROC,IERROR)

IF(PRESENT(KPROCS)) THEN
  KPROCS=MPL_NUMPROC
ENDIF

ALLOCATE (MPL_IDS(MPL_NUMPROC))
DO IP=1,MPL_NUMPROC
  MPL_IDS(IP)=IP
ENDDO

CALL MPI_COMM_RANK(MPL_COMM, IRANK, IERROR)
MPL_RANK=IRANK+1

LLINFO = LLINFO .AND. (MPL_RANK <= 1)

IF (LLINFO) THEN
   IF(LMPLUSERCOMM) THEN
      WRITE(MPL_UNIT,'(A)')'MPL_INIT : LMPLUSERCOMM used'
      WRITE(MPL_UNIT,'(A,I0)')'Communicator : ',MPL_COMM
   ELSE
      WRITE(MPL_UNIT,'(A)')'MPL_INIT : LMPLUSERCOMM not used'
      WRITE(MPL_UNIT,'(A,I0)')'Communicator : ',MPL_COMM
   ENDIF
ENDIF

#ifndef NECSX

!-- Propagate environment variables & argument lists
!   Here we have to be careful and use MPI_BCAST directly (not MPL_BROADCAST) since
!   1) MPL_BUFFER_METHOD has not been called
!   2) MPL_COMM_OML has not been initialized since it is possible that only the 
!      master proc knows the # of threads (i.e. OMP_NUM_THREADS may be set only for master)

! Do not propagate on nec machine because the environment variables could be mpi-task-specific.

IF (MPL_NUMPROC > 1 .AND. LLENV) THEN
  IROOT = 0
  !-- Progate environment variables
  INUM(1) = 0 ! The number of environment variables
  INUM(2) = 0 ! Do not (=0) or do (=1) overwrite if particular environment variable already exists (0 = default)
  IF (MPL_RANK == 1) THEN ! Master proc inquires
    CALL EC_NUMENV(INUM(1))        ! ../support/env.c
    CALL EC_OVERWRITE_ENV(INUM(2)) ! ../support/env.c
  ENDIF
  ! The following broadcast does not use "mailbox" nor attached buffer, both potentially yet to be allocated
  CALL MPI_BCAST(INUM(1),2,INT(MPI_INTEGER),IROOT,MPL_COMM,IERROR)
  ICOUNT = LEN(CLENV)
  DO IP=1,INUM(1)
    IF (MPL_RANK == 1) CALL EC_STRENV(IP,CLENV)
    ! The following broadcast does not use "mailbox" nor attached buffer, both potentially yet to be allocated
    CALL MPI_BCAST(CLENV,ICOUNT,INT(MPI_BYTE),IROOT,MPL_COMM,IERROR)
    IF (MPL_RANK > 1) THEN
      IF (INUM(2) == 1) THEN
        CALL EC_PUTENV(CLENV) ! ../support/env.c ; Unconditionally overwrite, even if already exists
      ELSE
        CALL EC_PUTENV_NOOVERWRITE(CLENV) ! ../support/env.c ; Do not overwrite, if exists
      ENDIF
    ENDIF
  ENDDO
  !-- Redo some env. variables (see ../utilities/fnecsx.c)
  CALL EC_ENVREDO()
  !-- Propagate argument list (all under the bonnet using MPL_ARG_MOD-module)
  INUM = MPL_IARGC()
ENDIF

#endif

CALL OML_INIT()
IMAX_THREADS = OML_MAX_THREADS()
ALLOCATE(MPL_COMM_OML(IMAX_THREADS))

IF (LMPLUSERCOMM) THEN
   MPL_COMM_OML(1) = MPLUSERCOMM
   ISTA = 2
ELSE
   ISTA = 1
ENDIF

DO IP=ISTA,IMAX_THREADS
  CALL MPL_LOCOMM_CREATE(MPL_NUMPROC,MPL_COMM_OML(IP))
ENDDO
MPL_COMM = MPL_COMM_OML(1) ! i.e. not necessary MPI_COMM_WORLD anymore

#ifdef VPP
MPL_METHOD=JP_BLOCKING_STANDARD
MPL_MBX_SIZE=4000000
CL_MBX_SIZE=' '
CALL GET_ENVIRONMENT_VARIABLE('VPP_MBX_SIZE',CL_MBX_SIZE)
IF(CL_MBX_SIZE == ' ') THEN
  CALL GET_ENVIRONMENT_VARIABLE('MPL_MBX_SIZE',CL_MBX_SIZE)
ENDIF
IF(CL_MBX_SIZE /= ' ') THEN
  READ(CL_MBX_SIZE,*) MPL_MBX_SIZE
ENDIF
IF (LLINFO) WRITE(MPL_UNIT,'(A)')'MPL_INIT : MPL_METHOD=JP_BLOCKING_STANDARD'
IF (LLINFO) WRITE(MPL_UNIT,'(A,I0)')'MPL_INIT : MAILBOX SIZE=',MPL_MBX_SIZE
LUSEHLMPI = .FALSE.

!#elif defined (LINUX)
!MPL_METHOD=JP_BLOCKING_STANDARD
!MPL_MBX_SIZE=4000000
!CL_MBX_SIZE=' '
!CALL GET_ENVIRONMENT_VARIABLE('VPP_MBX_SIZE',CL_MBX_SIZE)
!IF(CL_MBX_SIZE == ' ') THEN
!  CALL GET_ENVIRONMENT_VARIABLE('MPL_MBX_SIZE',CL_MBX_SIZE)
!ENDIF
!IF(CL_MBX_SIZE /= ' ') THEN
!  READ(CL_MBX_SIZE,*) MPL_MBX_SIZE
!ENDIF
!IF (LLINFO) WRITE(MPL_UNIT,'(A)')'MPL_INIT : MPL_METHOD=JP_BLOCKING_STANDARD'
!IF (LLINFO) WRITE(MPL_UNIT,'(A,I0)')'MPL_INIT : MAILBOX SIZE=',MPL_MBX_SIZE
!LUSEHLMPI = .FALSE.

#else
CL_METHOD=' '
CALL GET_ENVIRONMENT_VARIABLE('MPL_METHOD',CL_METHOD)
IF (CL_METHOD == 'JP_BLOCKING_STANDARD' ) THEN
  MPL_METHOD=JP_BLOCKING_STANDARD
ELSE
  MPL_METHOD=JP_BLOCKING_BUFFERED
ENDIF
MPL_MBX_SIZE=1000000
CL_MBX_SIZE=' '
CALL GET_ENVIRONMENT_VARIABLE('MPL_MBX_SIZE',CL_MBX_SIZE)
IF (CL_MBX_SIZE /= ' ') THEN
  READ(CL_MBX_SIZE,*) MPL_MBX_SIZE
ENDIF
IF (CL_METHOD == 'JP_BLOCKING_STANDARD' ) THEN
  IF (LLINFO) WRITE(MPL_UNIT,'(A)')'MPL_INIT : MPL_METHOD=JP_BLOCKING_STANDARD'
ELSE
  IF (LLINFO) WRITE(MPL_UNIT,'(A)')'MPL_INIT : MPL_METHOD=JP_BLOCKING_BUFFERED'
ENDIF
!IF (LLINFO) WRITE(MPL_UNIT,'(A,I0)')'MPL_INIT : MAILBOX SIZE=',MPL_MBX_SIZE

CALL MPL_BUFFER_METHOD(KMP_TYPE=MPL_METHOD,KMBX_SIZE=MPL_MBX_SIZE,LDINFO=LLINFO)
LUSEHLMPI = .TRUE.
#endif

CALL MPI_COMM_RANK (MPI_COMM_WORLD, IWORLD_RANK, IERROR)
CALL MPI_COMM_SIZE (MPI_COMM_WORLD, IWORLD_SIZE, IERROR)

#ifdef LINUX
CALL LINUX_BIND (IWORLD_RANK, IWORLD_SIZE)
#endif

!-- World-wide tasks
MPL_WORLD_RANK = IWORLD_RANK
MPL_WORLD_SIZE = IWORLD_SIZE

!!!! If you are not at ECMWF this may need changing!!!!
CALL GET_ENVIRONMENT_VARIABLE('EC_TASKS_PER_NODE',CL_TASKSPERNODE)
IF (CL_TASKSPERNODE(1:1) == ' ' ) THEN
   CALL GET_ENVIRONMENT_VARIABLE('HOST',CL_HOST)
   IF(CL_HOST(1:3) == 'cck') THEN ! KNL
      MPL_NCPU_PER_NODE=64
   ELSEIF(CL_HOST(1:3) == 'cct') THEN ! Test-cluster
      MPL_NCPU_PER_NODE=24
   ELSEIF(CL_HOST(1:2) == 'cc') THEN ! cca/ccb
      MPL_NCPU_PER_NODE=36
   ELSEIF(CL_HOST(1:3) == 'lxg') THEN ! GPU-cluster
      MPL_NCPU_PER_NODE=24
   ELSEIF (CL_HOST(1:2) == 'c1') THEN
      MPL_NCPU_PER_NODE=64
   ELSEIF(CL_HOST(1:3) == 'hpc') THEN
      MPL_NCPU_PER_NODE=32
   ELSE
      MPL_NCPU_PER_NODE=1
      IF(LLINFO) WRITE(MPL_UNIT,'(A)')'MPL_INIT CAUTION: MPL_NCPU_PER_NODE=1'
   ENDIF
ELSE
   READ(CL_TASKSPERNODE,*) MPL_NCPU_PER_NODE  
ENDIF
MPL_MAX_TASK_PER_NODE=MAX(1, MPL_NCPU_PER_NODE/IMAX_THREADS)
LFULLNODES=MOD(MPL_NUMPROC,MPL_MAX_TASK_PER_NODE) == 0
MPL_NNODES=(MPL_NUMPROC-1)/MPL_MAX_TASK_PER_NODE+1
ALLOCATE(MPL_TASK_PER_NODE(MPL_NNODES))
ALLOCATE(MPL_NODE(MPL_NUMPROC))
ALLOCATE(MPL_NODE_TASKS(MPL_NNODES,MPL_MAX_TASK_PER_NODE))
MPL_NODE_TASKS(:,:)=-99
ICOUNT=0
DO JNODE=1,MPL_NNODES
  DO JROC=1,MPL_MAX_TASK_PER_NODE
    ICOUNT=ICOUNT+1
    IF (ICOUNT<=MPL_NUMPROC) THEN
       MPL_NODE(ICOUNT)=JNODE
       MPL_TASK_PER_NODE(JNODE) = JROC
       MPL_NODE_TASKS(JNODE,JROC) = ICOUNT
    ENDIF
  ENDDO
ENDDO
MPL_MYNODE=(MPL_RANK-1)/MPL_MAX_TASK_PER_NODE+1
!WRITE(MPL_UNIT,*) 'MPL_INIT : NCPU_PER_NODE,MPL_MAX_TASK_PER_NODE,MPL_NNODES,MPL_MYNODE ',&
! & MPL_NCPU_PER_NODE,MPL_MAX_TASK_PER_NODE,MPL_NNODES,MPL_MYNODE
!WRITE(MPL_UNIT,*) 'MPL_INIT : MPL_NODE_TASKS(MPL_MYNODE,1:MPL_TASK_PER_NODE(MPL_MYNODE)) ', &
! & MPL_NODE_TASKS(MPL_MYNODE,1:MPL_TASK_PER_NODE(MPL_MYNODE))

ALLOCATE(MPL_OPPONENT(MPL_NUMPROC+1))
CALL MPL_TOUR_TABLE(MPL_OPPONENT)

RETURN
END SUBROUTINE MPL_INIT

END MODULE MPL_INIT_MOD
