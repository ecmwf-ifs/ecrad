! radiation_aerosol.F90 - Derived type describing aerosol
!
! (C) Copyright 2014- ECMWF.
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
! Modifications
!   2018-04-15  R. Hogan  Add "direct" option
!   2019-01-14  R. Hogan  Added out_of_physical_bounds routine

module radiation_aerosol

  use parkind1, only : jprb
  use radiation_io, only : nulerr, radiation_abort

  implicit none
  public

  !---------------------------------------------------------------------
  ! Type describing the aerosol content in the atmosphere
  type aerosol_type
     ! The mass mixing ratio of config%n_aerosol_types different
     ! aerosol types dimensioned
     ! (ncol,istartlev:iendlev,config%n_aerosol_types), where ncol is
     ! the number of columns, istartlev:iendlev is the range of model
     ! levels where aerosols are present
     real(jprb), allocatable, dimension(:,:,:) :: &
          &  mixing_ratio  ! mass mixing ratio (kg/kg)

     ! Alternatively, if is_direct=true, the optical properties are
     ! provided directly and are dimensioned
     ! (nband,istartlev:iendlev,ncol)
     real(jprb), allocatable, dimension(:,:,:) :: &
          &  od_sw, ssa_sw, g_sw, & ! Shortwave optical properties
          &  od_lw, ssa_lw, g_lw    ! Longwave optical properties

     ! Range of levels in which the aerosol properties are provided
     integer :: istartlev, iendlev

     ! Are the optical properties going to be provided directly by the
     ! user?
     logical :: is_direct = .false.

   contains
      procedure :: allocate        => allocate_aerosol_arrays
      procedure :: allocate_direct => allocate_aerosol_arrays_direct
      procedure :: deallocate      => deallocate_aerosol_arrays
      procedure :: out_of_physical_bounds
#ifdef _OPENACC
      procedure :: create_device
      procedure :: update_host
      procedure :: update_device
      procedure :: delete_device
#endif
  end type aerosol_type

contains

  !---------------------------------------------------------------------
  ! Allocate array for describing aerosols, although in the offline
  ! code these are allocated when they are read from the NetCDF file
  subroutine allocate_aerosol_arrays(this, ncol, istartlev, iendlev, ntype)

    use yomhook,     only : lhook, dr_hook, jphook

    class(aerosol_type), intent(inout) :: this
    integer, intent(in)                :: ncol  ! Number of columns
    integer, intent(in)                :: istartlev, iendlev ! Level range
    integer, intent(in)                :: ntype ! Number of aerosol types
    real(jphook) :: hook_handle

    if (lhook) call dr_hook('radiation_aerosol:allocate',0,hook_handle)

    allocate(this%mixing_ratio(ncol,istartlev:iendlev,ntype))
    ! !$ACC ENTER DATA CREATE(this%mixing_ratio)
    this%is_direct = .false.
    this%istartlev = istartlev
    this%iendlev   = iendlev
    ! !$ACC UPDATE DEVICE(this%is_direct, this%istartlev, this%iendlev)


    if (lhook) call dr_hook('radiation_aerosol:allocate',1,hook_handle)

  end subroutine allocate_aerosol_arrays


  !---------------------------------------------------------------------
  ! Allocate arrays for describing aerosol optical properties
  subroutine allocate_aerosol_arrays_direct(this, config, &
       &                                    ncol, istartlev, iendlev)

    use yomhook,          only : lhook, dr_hook, jphook
    use radiation_config, only : config_type

    class(aerosol_type), intent(inout) :: this
    type(config_type),   intent(in)    :: config
    integer, intent(in)                :: ncol  ! Number of columns
    integer, intent(in)                :: istartlev, iendlev ! Level range
    integer                            :: jband, jlev, jcol

    real(jphook) :: hook_handle

    if (lhook) call dr_hook('radiation_aerosol:allocate_direct',0,hook_handle)

    this%is_direct = .true.
    this%istartlev = istartlev
    this%iendlev   = iendlev

    if (config%do_sw) then
      allocate(this%od_sw (config%n_bands_sw,istartlev:iendlev,ncol))
      allocate(this%ssa_sw(config%n_bands_sw,istartlev:iendlev,ncol))
      allocate(this%g_sw  (config%n_bands_sw,istartlev:iendlev,ncol))
      ! !$ACC ENTER DATA CREATE(this%od_sw, this%ssa_sw, this%g_sw) ASYNC(1)
    end if

    if (config%do_lw) then
      allocate(this%od_lw (config%n_bands_lw,istartlev:iendlev,ncol))
      allocate(this%ssa_lw(config%n_bands_lw,istartlev:iendlev,ncol))
      allocate(this%g_lw  (config%n_bands_lw,istartlev:iendlev,ncol))
      ! !$ACC ENTER DATA CREATE(this%od_lw, this%ssa_lw, this%g_lw) ASYNC(1)
      ! If longwave scattering by aerosol is not to be represented,
      ! then the user may wish to just provide absorption optical
      ! depth in od_lw, in which case we must set the following two
      ! variables to zero

      ! !$ACC WAIT ! ACCWA (nvhpc 22.7) crashes otherwise

      ! !$ACC PARALLEL DEFAULT(NONE) PRESENT(this, config) ASYNC(1)
      ! !$ACC LOOP GANG VECTOR COLLAPSE(3)
      do jcol = 1,ncol
        do jlev = istartlev,iendlev
          do jband = 1,config%n_bands_lw
            this%ssa_lw(jband,jlev,jcol) = 0.0_jprb
            this%g_lw(jband,jlev,jcol) = 0.0_jprb
          end do
        end do
      end do
      ! !$ACC END PARALLEL
    end if

    if (lhook) call dr_hook('radiation_aerosol:allocate_direct',1,hook_handle)

  end subroutine allocate_aerosol_arrays_direct


  !---------------------------------------------------------------------
  ! Deallocate array
  subroutine deallocate_aerosol_arrays(this)

    use yomhook,     only : lhook, dr_hook, jphook

    class(aerosol_type), intent(inout) :: this

    real(jphook) :: hook_handle

    if (lhook) call dr_hook('radiation_aerosol:deallocate',0,hook_handle)

    ! !$ACC EXIT DATA DELETE(this%mixing_ratio) ASYNC(1) IF(allocated(this%mixing_ratio))
    ! !$ACC EXIT DATA DELETE(this%od_sw) ASYNC(1) IF(allocated(this%od_sw))
    ! !$ACC EXIT DATA DELETE(this%ssa_sw) ASYNC(1) IF(allocated(this%ssa_sw))
    ! !$ACC EXIT DATA DELETE(this%g_sw) ASYNC(1) IF(allocated(this%g_sw))
    ! !$ACC EXIT DATA DELETE(this%od_lw) ASYNC(1) IF(allocated(this%od_lw))
    ! !$ACC EXIT DATA DELETE(this%ssa_lw) ASYNC(1) IF(allocated(this%ssa_lw))
    ! !$ACC EXIT DATA DELETE(this%g_lw) ASYNC(1) IF(allocated(this%g_lw))
    ! !$ACC WAIT
    if (allocated(this%mixing_ratio)) deallocate(this%mixing_ratio)
    if (allocated(this%od_sw))        deallocate(this%od_sw)
    if (allocated(this%ssa_sw))       deallocate(this%ssa_sw)
    if (allocated(this%g_sw))         deallocate(this%g_sw)
    if (allocated(this%od_lw))        deallocate(this%od_lw)
    if (allocated(this%ssa_lw))       deallocate(this%ssa_lw)
    if (allocated(this%g_lw))         deallocate(this%g_lw)

    if (lhook) call dr_hook('radiation_aerosol:deallocate',1,hook_handle)

  end subroutine deallocate_aerosol_arrays


  !---------------------------------------------------------------------
  ! Return .true. if variables are out of a physically sensible range,
  ! optionally only considering columns between istartcol and iendcol
  function out_of_physical_bounds(this, istartcol, iendcol, do_fix) result(is_bad)

    use yomhook,          only : lhook, dr_hook, jphook
    use radiation_check,  only : out_of_bounds_3d

    class(aerosol_type),   intent(inout) :: this
    integer,      optional,intent(in) :: istartcol, iendcol
    logical,      optional,intent(in) :: do_fix
    logical                           :: is_bad

    logical    :: do_fix_local

    real(jphook) :: hook_handle

    if (lhook) call dr_hook('radiation_aerosol:out_of_physical_bounds',0,hook_handle)

    if (present(do_fix)) then
      do_fix_local = do_fix
    else
      do_fix_local = .false.
    end if

    is_bad =    out_of_bounds_3d(this%mixing_ratio, 'aerosol%mixing_ratio', &
         &                       0.0_jprb, 1.0_jprb, do_fix_local, i1=istartcol, i2=iendcol) &
         & .or. out_of_bounds_3d(this%od_sw, 'aerosol%od_sw', &
         &                       0.0_jprb, 100.0_jprb, do_fix_local, k1=istartcol, k2=iendcol) &
         & .or. out_of_bounds_3d(this%od_lw, 'aerosol%od_lw', &
         &                       0.0_jprb, 100.0_jprb, do_fix_local, k1=istartcol, k2=iendcol) &
         & .or. out_of_bounds_3d(this%ssa_sw, 'aerosol%ssa_sw', &
         &                       0.0_jprb, 1.0_jprb, do_fix_local, k1=istartcol, k2=iendcol) &
         & .or. out_of_bounds_3d(this%ssa_lw, 'aerosol%ssa_lw', &
         &                       0.0_jprb, 1.0_jprb, do_fix_local, k1=istartcol, k2=iendcol) &
         & .or. out_of_bounds_3d(this%g_sw, 'aerosol%g_sw', &
         &                       0.0_jprb, 1.0_jprb, do_fix_local, k1=istartcol, k2=iendcol) &
         & .or. out_of_bounds_3d(this%g_lw, 'aerosol%g_lw', &
         &                       0.0_jprb, 1.0_jprb, do_fix_local, k1=istartcol, k2=iendcol)

    if (lhook) call dr_hook('radiation_aerosol:out_of_physical_bounds',1,hook_handle)

  end function out_of_physical_bounds

#ifdef _OPENACC
  !---------------------------------------------------------------------
  ! Creates fields on device
  subroutine create_device(this)

    class(aerosol_type), intent(inout) :: this

    !$ACC ENTER DATA CREATE(this%mixing_ratio) IF(allocated(this%mixing_ratio))
    !$ACC ENTER DATA CREATE(this%od_sw) IF(allocated(this%od_sw))
    !$ACC ENTER DATA CREATE(this%ssa_sw) IF(allocated(this%ssa_sw))
    !$ACC ENTER DATA CREATE(this%g_sw) IF(allocated(this%g_sw))
    !$ACC ENTER DATA CREATE(this%od_lw) IF(allocated(this%od_lw))
    !$ACC ENTER DATA CREATE(this%ssa_lw) IF(allocated(this%ssa_lw))
    !$ACC ENTER DATA CREATE(this%g_lw) IF(allocated(this%g_lw))

  end subroutine create_device

  !---------------------------------------------------------------------
  ! updates fields on host
  subroutine update_host(this)

    class(aerosol_type), intent(inout) :: this

    !$ACC UPDATE HOST(this%mixing_ratio) IF(allocated(this%mixing_ratio))
    !$ACC UPDATE HOST(this%od_sw) IF(allocated(this%od_sw))
    !$ACC UPDATE HOST(this%ssa_sw) IF(allocated(this%ssa_sw))
    !$ACC UPDATE HOST(this%g_sw) IF(allocated(this%g_sw))
    !$ACC UPDATE HOST(this%od_lw) IF(allocated(this%od_lw))
    !$ACC UPDATE HOST(this%ssa_lw) IF(allocated(this%ssa_lw))
    !$ACC UPDATE HOST(this%g_lw) IF(allocated(this%g_lw))

  end subroutine update_host

  !---------------------------------------------------------------------
  ! updates fields on device
  subroutine update_device(this)

    class(aerosol_type), intent(inout) :: this

    !$ACC UPDATE DEVICE(this%mixing_ratio) IF(allocated(this%mixing_ratio))
    !$ACC UPDATE DEVICE(this%od_sw) IF(allocated(this%od_sw))
    !$ACC UPDATE DEVICE(this%ssa_sw) IF(allocated(this%ssa_sw))
    !$ACC UPDATE DEVICE(this%g_sw) IF(allocated(this%g_sw))
    !$ACC UPDATE DEVICE(this%od_lw) IF(allocated(this%od_lw))
    !$ACC UPDATE DEVICE(this%ssa_lw) IF(allocated(this%ssa_lw))
    !$ACC UPDATE DEVICE(this%g_lw) IF(allocated(this%g_lw))

  end subroutine update_device

  !---------------------------------------------------------------------
  ! Deletes fields on device
  subroutine delete_device(this)

    class(aerosol_type), intent(inout) :: this

    !$ACC EXIT DATA DELETE(this%mixing_ratio) IF(allocated(this%mixing_ratio))
    !$ACC EXIT DATA DELETE(this%od_sw) IF(allocated(this%od_sw))
    !$ACC EXIT DATA DELETE(this%ssa_sw) IF(allocated(this%ssa_sw))
    !$ACC EXIT DATA DELETE(this%g_sw) IF(allocated(this%g_sw))
    !$ACC EXIT DATA DELETE(this%od_lw) IF(allocated(this%od_lw))
    !$ACC EXIT DATA DELETE(this%ssa_lw) IF(allocated(this%ssa_lw))
    !$ACC EXIT DATA DELETE(this%g_lw) IF(allocated(this%g_lw))

  end subroutine delete_device
#endif

end module radiation_aerosol
