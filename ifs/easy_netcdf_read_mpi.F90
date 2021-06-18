! easy_netcdf_read_mpi.f90 - Read netcdf file on one task and share with other tasks
!
! (C) Copyright 2017- ECMWF.
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

module easy_netcdf_read_mpi

  use easy_netcdf,   only : netcdf_file_raw => netcdf_file
  use parkind1,      only : jpim, jprb
  use radiation_io,  only : nulout, nulerr, my_abort => radiation_abort

  implicit none

  ! MPI tag for radiation and physics communication
  integer(kind=jpim), parameter :: mtagrad = 2800

  !---------------------------------------------------------------------
  ! An object of this type provides convenient read or write access to
  ! a NetCDF file
  type netcdf_file
    type(netcdf_file_raw) :: file
    logical :: is_master_task = .true.
  contains
    procedure :: open => open_netcdf_file
    procedure :: close => close_netcdf_file
    procedure :: get_real_scalar
    procedure :: get_real_vector
    procedure :: get_real_matrix
    procedure :: get_real_array3
    generic   :: get => get_real_scalar, get_real_vector, &
         &              get_real_matrix, get_real_array3
    procedure :: get_global_attribute

    procedure :: set_verbose
    procedure :: transpose_matrices
    procedure :: exists
  end type netcdf_file

contains

  ! --- GENERIC SUBROUTINES ---

  !---------------------------------------------------------------------
  ! Open a NetCDF file with name "file_name", optionally specifying the
  ! verbosity level (0-5)
  subroutine open_netcdf_file(this, file_name, iverbose)

    USE MPL_MODULE, ONLY : MPL_BROADCAST, MPL_MYRANK

    class(netcdf_file)            :: this
    character(len=*), intent(in)  :: file_name
    integer, intent(in), optional :: iverbose

    integer                       :: istatus

    ! Store verbosity level in object
    if (present(iverbose)) then
      this%file%iverbose = iverbose
    else
      ! By default announce files being opened and closed, but not
      ! variables read/written
      this%file%iverbose = 2
    end if

    ! By default we don't transpose 2D arrays on read
    this%file%do_transpose_2d = .false.

    if (MPL_MYRANK() == 1) then
      this%is_master_task = .true.
      call this%file%open(file_name, iverbose)
    else
      this%is_master_task = .false.
    end if

  end subroutine open_netcdf_file


  !---------------------------------------------------------------------
  ! Close the NetCDF file
  subroutine close_netcdf_file(this)
    class(netcdf_file) :: this
    integer            :: istatus

    if (this%is_master_task) then
      call this%file%close()
    end if

  end subroutine close_netcdf_file


  !---------------------------------------------------------------------
  ! Set the verbosity level from 0 to 5, where the codes have the
  ! following meaning: 0=errors only, 1=warning, 2=info, 3=progress,
  ! 4=detailed, 5=debug
  subroutine set_verbose(this, ival)
    class(netcdf_file) :: this
    integer, optional  :: ival

    if (present(ival)) then
      this%file%iverbose = ival
    else
      this%file%iverbose = 2
    end if

  end subroutine set_verbose



  !---------------------------------------------------------------------
  ! Specify whether 2D arrays should be transposed on read
  subroutine transpose_matrices(this, do_transpose)
    class(netcdf_file) :: this
    logical, optional  :: do_transpose

    if (present(do_transpose)) then
      this%file%do_transpose_2d = do_transpose
    else
      this%file%do_transpose_2d = .true.
    end if

  end subroutine transpose_matrices



  ! --- READING SUBROUTINES ---

  !---------------------------------------------------------------------
  ! Return true if the variable is present, false otherwise
  function exists(this, var_name) result(is_present)

    USE MPL_MODULE, ONLY : MPL_BROADCAST, MPL_NPROC

    class(netcdf_file)           :: this
    character(len=*), intent(in) :: var_name

    logical :: is_present

    if (this%is_master_task) then
      is_present = this%file%exists(var_name)
    end if

    if (MPL_NPROC() > 1) then
      CALL MPL_BROADCAST(is_present, mtagrad, 1, &
           &  CDSTRING='EASY_NETCDF_READ_MPI:EXISTS')
    end if

  end function exists


  !---------------------------------------------------------------------
  ! The method "get" will read either a scalar, vector or matrix
  ! depending on the rank of the output argument. This version reads a
  ! scalar.
  subroutine get_real_scalar(this, var_name, scalar)

    USE MPL_MODULE, ONLY : MPL_BROADCAST, MPL_NPROC

    class(netcdf_file)           :: this
    character(len=*), intent(in) :: var_name
    real(jprb), intent(out)      :: scalar

    if (this%is_master_task) then
      call this%file%get(var_name, scalar)
    end if

    if (MPL_NPROC() > 1) then
      CALL MPL_BROADCAST(scalar, mtagrad, 1, &
           &  CDSTRING='EASY_NETCDF_READ_MPI:GET_REAL_SCALAR')
    end if

  end subroutine get_real_scalar


  !---------------------------------------------------------------------
  ! Read a 1D array into "vector", which must be allocatable and will
  ! be reallocated if necessary
  subroutine get_real_vector(this, var_name, vector)

    USE MPL_MODULE, ONLY : MPL_BROADCAST, MPL_NPROC

    class(netcdf_file)           :: this
    character(len=*), intent(in) :: var_name
    real(jprb), allocatable, intent(out) :: vector(:)

    integer                      :: n  ! Length of vector

    n = 0

    if (this%is_master_task) then
      call this%file%get(var_name, vector)
      n = size(vector)
    end if

    if (MPL_NPROC() > 1) then
      CALL MPL_BROADCAST(n, mtagrad, 1, &
           &  CDSTRING='EASY_NETCDF_READ_MPI:GET_REAL_VECTOR:SIZE')

      if (.not. this%is_master_task) then
        allocate(vector(n))
      end if

      CALL MPL_BROADCAST(vector, mtagrad, 1, &
           &  CDSTRING='EASY_NETCDF_READ_MPI:GET_REAL_VECTOR')
    end if

  end subroutine get_real_vector


  !---------------------------------------------------------------------
  ! Read 2D array into "matrix", which must be allocatable and will be
  ! reallocated if necessary.  Whether to transpose is specifed by the
  ! final optional argument, but can also be specified by the
  ! do_transpose_2d class data member.
  subroutine get_real_matrix(this, var_name, matrix, do_transp)

    USE MPL_MODULE, ONLY : MPL_BROADCAST, MPL_NPROC

    class(netcdf_file)           :: this
    character(len=*), intent(in) :: var_name
    real(jprb), allocatable, intent(out) :: matrix(:,:)
    logical, optional, intent(in):: do_transp ! Transpose data?

    integer                      :: n(2)

    n = 0

    if (this%is_master_task) then
      call this%file%get(var_name, matrix, do_transp)
      n = shape(matrix)
    end if

    if (MPL_NPROC() > 1) then
      CALL MPL_BROADCAST(n, mtagrad, 1, &
           &  CDSTRING='EASY_NETCDF_READ_MPI:GET_REAL_MATRIX:SIZE')

      if (.not. this%is_master_task) then
        allocate(matrix(n(1),n(2)))
      end if

      CALL MPL_BROADCAST(matrix, mtagrad, 1, &
           &  CDSTRING='EASY_NETCDF_READ_MPI:GET_REAL_MATRIX')
    end if

  end subroutine get_real_matrix


  !---------------------------------------------------------------------
  ! Read 3D array into "var", which must be allocatable and will be
  ! reallocated if necessary.  Whether to pemute is specifed by the
  ! final optional argument
  subroutine get_real_array3(this, var_name, var, ipermute)

    USE MPL_MODULE, ONLY : MPL_BROADCAST, MPL_NPROC

    class(netcdf_file)                   :: this
    character(len=*), intent(in)         :: var_name
    real(jprb), allocatable, intent(out) :: var(:,:,:)
    integer, optional, intent(in)        :: ipermute(3)

    integer                              :: n(3)

    n = 0

    if (this%is_master_task) then
      call this%file%get(var_name, var, ipermute)
      n = shape(var)
    end if

    if (MPL_NPROC() > 1) then
      CALL MPL_BROADCAST(n, mtagrad, 1, &
           &  CDSTRING='EASY_NETCDF_READ_MPI:GET_REAL_ARRAY3:SIZE')

      if (.not. this%is_master_task) then
        allocate(var(n(1),n(2),n(3)))
      end if

      CALL MPL_BROADCAST(var, mtagrad, 1, &
           &  CDSTRING='EASY_NETCDF_READ_MPI:GET_REAL_ARRAY3')
    end if

  end subroutine get_real_array3


  !---------------------------------------------------------------------
  ! Get a global attribute as a character string
  subroutine get_global_attribute(this, attr_name, attr_str)

    USE MPL_MODULE, ONLY : MPL_BROADCAST, MPL_NPROC

    class(netcdf_file) :: this

    character(len=*), intent(in)    :: attr_name
    character(len=*), intent(inout) :: attr_str

    if (this%is_master_task) then
      call this%file%get_global_attribute(attr_name, attr_str)
    end if

    if (MPL_NPROC() > 1) then
      CALL MPL_BROADCAST(attr_str, mtagrad, 1, &
           &  CDSTRING='EASY_NETCDF_READ_MPI:GET_GLOBAL_ATTRIBUTE')
    end if

  end subroutine get_global_attribute

end module easy_netcdf_read_mpi
