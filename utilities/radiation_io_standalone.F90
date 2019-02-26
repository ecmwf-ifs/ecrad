! radiation_io.F90 - Provides logging and abort functionality
!
! Copyright (C) 2015, 2019 ECMWF
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
!  This is a standalone version that removes the dependence on
!  yomlun_ifsaux.F90 and abor1.F90, and is suitable for the offline
!  ecRad scheme. If ecRad is to be embedded in another model, rewrite
!  this file according to the prefered file unit numbers and preferred
!  abort behaviour.

module radiation_io

  ! In the IFS, nulout is equivalent to standard output but only
  ! output from the primary node will be logged, while nulerr is
  ! equivalent to standard error, and text sent to this unit from any
  ! node will be logged. Normally, nulerr should only be used before
  ! calling radiation_abort.
  integer :: nulout = 6
  integer :: nulerr = 0

  ! This unit may be used for reading radiation configuration files,
  ! but should be closed as soon as the file is read
  integer :: nulrad = 25

contains

  ! Abort the program with optional error message. Normally you would
  ! log details of the error to nulerr before calling this subroutine.
  subroutine radiation_abort(text)

    implicit none

    character(len=*), intent(in), optional :: text

    if (present(text)) then
      write(nulerr,'(a)') text
    else
      write(nulerr,'(a)') '*** Error in radiation scheme'
    end if

    if (nulout >= 0) then
      call flush(nulout)
      if (nulout /= 6) then
        close(nulout)
      end if
    end if

    call abort

  end subroutine radiation_abort

end module radiation_io
