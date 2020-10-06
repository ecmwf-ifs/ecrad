! radiation_io.F90 - Provides logging and abort functionality
!
! (C) Copyright 2015- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
!
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!
! Author:  Robin Hogan
! Email:   r.j.hogan@ecmwf.int
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
  public

  ! This unit may be used for reading radiation configuration files,
  ! but should be closed as soon as the file is read
  integer :: nulrad = 25

contains

  ! Abort the program with optional error message. Normally you would
  ! log details of the error to nulerr before calling this subroutine.
  subroutine radiation_abort(text)
    character(len=*), intent(in), optional :: text
    if (present(text)) then
      write(nulerr,'(a)') text
#ifdef __PGI
      stop 1
#else
      error stop 1
#endif
    else
#ifdef __PGI
      stop 'Error in radiation scheme'
#else
      error stop 'Error in radiation scheme'
#endif
    end if
  end subroutine radiation_abort

end module radiation_io
