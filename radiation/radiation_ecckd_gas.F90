! radiation_ecckd_gas.F90 - type representing a single ecCKD gas
!
! Copyright (C) 2020 ECMWF
!
! Author:  Robin Hogan
! Email:   r.j.hogan@ecmwf.int
! License: see the COPYING file for details
!

module radiation_ecckd_gas

  use parkind1, only : jprb
  use easy_netcdf
  use radiation_gas_constants

  implicit none

  public

  ! Concentration dependence of individual gases
  enum, bind(c)
    enumerator :: IConcDependenceNone = 0, &
         &        IConcDependenceLinear, &
         &        IConcDependenceLUT, &
         &        IConcDependenceRelativeLinear
  end enum

  !---------------------------------------------------------------------
  ! This derived type describes a correlated k-distribution
  ! representation of an individual gas (including composite gases)
  type ckd_gas_type

    ! Code identifying the gas, from the codes in the
    ! radiation_gas_constants module
    integer :: i_gas_code

    ! One of the IConcDependence* enumerators
    integer :: i_conc_dependence

    ! Molar absorption coefficient in m2 mol-1. If
    ! i_conc_dependence==IConcDependenceNone then it is the absorption
    ! cross section per mole of dry air.  If
    ! conc_dependence==IConcDependenceLinear|IConcDependenceRelativeLinear,
    ! it is the absorption cross section per mole of the gas in
    ! question. It is dimensioned (g_point,pressure,temperature).
    real(jprb), allocatable :: molar_abs(:,:,:)
    
    ! If i_conc_dependence==IConcDependenceLUT then we have an
    ! additional dimension for concentration. It is dimensioned
    ! (g_point,pressure,temperature,conc)
    real(jprb), allocatable :: molar_abs_conc(:,:,:,:)

    ! If i_conc_dependence==IConcDependenceRelativeLinear then the
    ! following reference concentration is subtracted from the actual
    ! concentration before the result is multiplied by the mass
    ! absorption coefficient
    real(jprb) :: reference_mole_fraction = 0.0_jprb

    ! Mole fraction coordinate variable if
    ! i_conc_dependence==IConcDependenceLUT
    real(jprb), allocatable :: mole_fraction(:)

  contains

    procedure :: read => read_ckd_gas
!    procedure :: calc_optical_depth => calc_optical_depth_ckd_gas
!    procedure :: deallocate => deallocate_ckd_gas

  end type ckd_gas_type

contains

  !---------------------------------------------------------------------
  ! Read information about the representation of a single gas from a
  ! NetCDF file, identifying it with code i_gas_code
  subroutine read_ckd_gas(this, file, gas_name, i_gas_code)

    class(ckd_gas_type), intent(inout) :: this
    type(netcdf_file),   intent(inout) :: file
    character(len=*),    intent(in)    :: gas_name
    integer,             intent(in)    :: i_gas_code
    
    this%i_gas_code = i_gas_code

    call file%get(gas_name // "_conc_dependence", this%i_conc_dependence)
    if (this%i_conc_dependence == IConcDependenceLut) then
      call file%get(gas_name // "_molar_absorption_coeff", &
           &        this%molar_abs_conc)
      call file%get(gas_name // "_mole_fraction", this%mole_fraction)
    else
      call file%get(gas_name // "_molar_absorption_coeff", &
           &        this%molar_abs)
    end if

    if (this%i_conc_dependence == IConcDependenceRelativeLinear) then
      call file%get(gas_name // "_reference_mole_fraction", &
           &        this%reference_mole_fraction)
    end if

  end subroutine read_ckd_gas

end module radiation_ecckd_gas
