MODULE MPL_IOINIT_MOD

!
!     Purpose.  initialise parallel IO environment
!     --------
!
!
!     Interface.
!     ----------
!        call mpl_ioinit(...)
!
!        Explicit arguments :
!        --------------------
!
!        input arguments:
!        kop     -  Style of parallel IO
!        kstrout - Number of output processors
!        output arguments:
!        kerror  - error code
!
!        Implicit arguments :
!        --------------------
!
!     Method.
!     -------
!     MPL supports 4 styles of MPIIO
!
!     kop = 1    -  Blocking, non collective, shared file pointer
!                   using MPI_FILE_WRITE_SHARED,
!                         MPI_FILE_READ_SHARED
!     kop = 2    -  Blocking, collective, ordered, shared file pointer
!                   using MPI_FILE_WRITE_ORDERED,
!                         MPI_FILE_READ_ORDERED
!     kop = 3    -  Non Blocking, non collective, shared file pointer
!                   using MPI_FILE_IWRITE_SHARED,
!                         MPI_FILE_IREAD_SHARED
!                   and MPI_WAIT
!     kop = 4    -  Non Blocking, collective, ordered, shared file pointer 
!                   using MPI_FILE_WRITE_ORDERED_BEGIN/END,
!                         MPI_FILE_READ_ORDERED_BEGIN/END
!
!     Externals.
!     ----------
!
!     Reference.
!     ----------
!        none yet
!
!     Author.
!     -------
!        G.Mozdzynski
!
!     Modifications.
!     --------------
!        Original : 2000-12-08 (Based on MPE_IOINIT)
!
!     -----------------------------------------------------------------
!

USE PARKIND1  ,ONLY : JPIM     ,JPRB

USE MPL_MPIF
USE MPL_DATA_MODULE
USE MPL_MESSAGE_MOD

IMPLICIT NONE


INTEGER(KIND=JPIM) :: MPL_NUMIO
INTEGER(KIND=JPIM) :: MPL_IOP
INTEGER(KIND=JPIM) :: MPL_COMM_IO

PRIVATE
PUBLIC :: MPL_IOINIT, MPL_NUMIO, MPL_IOP, MPL_COMM_IO

CONTAINS

SUBROUTINE MPL_IOINIT(KOP,KSTROUT,KERROR)


#ifdef USE_8_BYTE_WORDS
  USE MPI4TO8, ONLY : &
    MPI_COMM_SPLIT => MPI_COMM_SPLIT8
#endif

INTEGER(KIND=JPIM),INTENT(IN) :: KOP,KSTROUT
INTEGER(KIND=JPIM),INTENT(OUT) :: KERROR
INTEGER(KIND=JPIM) :: COLOR,KEY

!
!     -----------------------------------------------------------------
!
!     1.    Preamble
!           --------

IF(KOP < 1 .OR. KOP > 4) THEN

  KERROR = -1
  RETURN

ENDIF

!
!     -----------------------------------------------------------------
!
!     2. Check Style of Operation and take appropriate action
!     -------------------------------------------------------

MPL_NUMIO = KSTROUT

MPL_IOP = KOP

IF(MPL_RANK <= KSTROUT) THEN
  COLOR = 1
ELSE
  COLOR = MPI_UNDEFINED
ENDIF

KEY = 0

CALL MPI_COMM_SPLIT(MPL_COMM,COLOR,KEY,MPL_COMM_IO,KERROR)

!
!
!     -----------------------------------------------------------------
RETURN
END SUBROUTINE MPL_IOINIT

END MODULE MPL_IOINIT_MOD
