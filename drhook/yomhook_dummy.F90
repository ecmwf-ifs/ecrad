! (C) Copyright 2014- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
!
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!
! DR_HOOK is a profiling and debugging system for the IFS, and should
! be called at the beginning and end of each subroutine. This is a
! dummy implementation for offline packages.

module yomhook

  public

  logical :: lhook = .true.

  ! allocate two large arrays
  DOUBLE PRECISION, allocatable :: total_time(:)
  DOUBLE PRECISION, allocatable :: tstart(:), tstop(:)
  integer, allocatable :: ncalls(:)
  CHARACTER(len=80), allocatable :: names(:)
  integer, parameter :: hash_size = 1000

contains

  SUBROUTINE char_to_hash(c, a)
    CHARACTER(len=*), INTENT(in) :: c
    INTEGER, INTENT(out) :: a
    INTEGER :: i

    INTEGER :: p = 31
    INTEGER :: m = hash_size - 3
    INTEGER :: p_pow

    a = 0
    p_pow = 1
    DO i=1,len_trim(c)
        a = MOD(a + (ichar(c(i:i)) + 1) * p_pow, m)
        p_pow = MOD(p_pow * p, m)
    END DO

  END SUBROUTINE

  subroutine dr_hook(proc_name, iswitch, proc_key)

    use parkind1, only : jprb

    character(len=*), intent(in)    :: proc_name
    integer,          intent(in)    :: iswitch
    real(jprb),       intent(inout) :: proc_key
    integer :: idx
    double precision, external :: omp_get_wtime

    idx = INT(proc_key)
    if (iswitch == 0) then
      call char_to_hash(proc_name, idx)
      proc_key = real(idx)
      names(idx) = proc_name
      ncalls(idx) = ncalls(idx) + 1
      tstart(idx) = omp_get_wtime()
    else if (iswitch == 1) then
      tstop(idx) = omp_get_wtime()
      total_time(idx) = total_time(idx) + (tstop(idx) - tstart(idx))
    endif

  end subroutine dr_hook

  subroutine initialize_timers()

    allocate(total_time(hash_size))
    allocate(tstart(hash_size))
    allocate(tstop(hash_size))
    allocate(ncalls(hash_size))
    allocate(names(hash_size))

    total_time = 0
    ncalls = 0

  end subroutine initialize_timers

  subroutine finalize_timers()

    integer :: idx

    open(1, file="timing.txt", action="write")

    do idx = 1,1000
      if(total_time(idx) > 0.0) then
        write(1, '(A80,E10.3,I10)'), names(idx), total_time(idx), ncalls(idx)
      end if
    end do

    close(1)

  end subroutine finalize_timers

end module yomhook
