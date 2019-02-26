! radiation_io.F90 - Provides logging and abort functionality
!
! Copyright (C) 2015 ECMWF
!
! Author:  Robin Hogan
! Email:   r.j.hogan@ecmwf.int
! License: see the COPYING file for details
!
!  This file provides an interface to the provision of file units used
!  for logging (nulout and nulerr) and for reading data files
!  (nulrad), as well as an abort routine that should do clean-up
!  appropriate for the environment in which the radiation scheme is
!  embedded.
!
!  Rewrite this file as appropriate if the radiation scheme is to be
!  embedded into a model other than the ECMWF Integrated Forecasting
!  System.

module radiation_io

  ! In the IFS, nulout is equivalent to standard output but only
  ! output from the primary node will be logged, while nulerr is
  ! equivalent to standard error, and text sent to this unit from any
  ! node will be logged. Normally, nulerr should only be used before
  ! calling radiation_abort.
  use yomlun_ifsaux, only : nulout, nulerr

  implicit none

  ! This unit may be used for reading radiation configuration files,
  ! but should be closed as soon as the file is read
  integer :: nulrad = 25

  ! The abor1 subroutine is provided in the IFS and cleans up all MPI
  ! nodes.  Radiation routines should call radiation_abort instead.
  interface
    subroutine abor1(text)
      character(len=*) :: text
    end subroutine abor1
  end interface
  private :: abor1

contains

  ! Abort the program with optional error message. Normally you would
  ! log details of the error to nulerr before calling this subroutine.
  subroutine radiation_abort(text)
    character(len=*), intent(in), optional :: text
    if (present(text)) then
      call abor1(text)
    else
      call abor1('Error in radiation scheme')
    end if
  end subroutine radiation_abort

end module radiation_io
