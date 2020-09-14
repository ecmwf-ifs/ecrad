! radiation_spectral_definition.F90 - Derived type to describe a spectral definition
!
! Copyright (C) 2020 ECMWF
!
! Author:  Robin Hogan
! Email:   r.j.hogan@ecmwf.int
! License: see the COPYING file for details
!

module radiation_spectral_definition

  use parkind1,    only : jprb

  implicit none

  public

  type spectral_definition_type
    
    ! Spectral mapping of g points

    ! Number of wavenumber intervals
    integer :: nwav = 0
    ! Number of k terms / g points
    integer :: ng   = 0
    ! Start and end wavenumber (cm-1), dimensioned (nwav)
    real(jprb), allocatable :: wavenumber1(:)
    real(jprb), allocatable :: wavenumber2(:)
    ! Fraction of each g point in each wavenumber interval,
    ! dimensioned (nwav, ng)
    real(jprb), allocatable :: gpoint_fraction(:,:)

    ! Band information

    ! Number of bands
    integer :: nband = 0
    ! Lower and upper bounds of wavenumber bands (cm-1), dimensioned
    ! (nband)
    real(jprb), allocatable :: wavenumber1_band(:)
    real(jprb), allocatable :: wavenumber2_band(:)
    ! Band (one based) to which each g point belongs
    integer,    allocatable :: i_band_number(:)

  contains
    procedure :: read => read_spectral_definition

  end type spectral_definition_type

contains

  !---------------------------------------------------------------------
  ! Read the description of a spectral definition from a NetCDF
  ! file of the type used to describe an ecCKD model
  subroutine read_spectral_definition(this, file, iverbose)

    use easy_netcdf, only : netcdf_file
    use yomhook,     only : lhook, dr_hook

    class(spectral_definition_type), intent(inout) :: this
    type(netcdf_file),                   intent(inout) :: file
    integer,                   optional, intent(in)    :: iverbose

    real(jprb) :: hook_handle

    if (lhook) call dr_hook('radiation_spectral_definition:read',0,hook_handle)

    ! Read spectral mapping of g points
    call file%get('wavenumber1', this%wavenumber1)
    call file%get('wavenumber2', this%wavenumber2)
    call file%get('gpoint_fraction', this%gpoint_fraction)

    ! Read band information
    call file%get('wavenumber1_band', this%wavenumber1_band)
    call file%get('wavenumber2_band', this%wavenumber2_band)
    call file%get('band_number', this%i_band_number)

    this%nwav  = size(this%wavenumber1)
    this%ng    = size(this%gpoint_fraction, 2);
    this%nband = size(this%wavenumber1_band)

    if (lhook) call dr_hook('radiation_spectral_definition:read',1,hook_handle)

  end subroutine read_spectral_definition


end module radiation_spectral_definition
