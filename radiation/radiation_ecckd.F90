! radiation_ecckd.F90 - ecCKD generalized gas optics model
!
! Copyright (C) 2020 ECMWF
!
! Author:  Robin Hogan
! Email:   r.j.hogan@ecmwf.int
! License: see the COPYING file for details
!

module radiation_ecckd

  use parkind1, only : jprb
  use radiation_gas

  implicit none

  ! Concentration dependence of individual gases
  enum, bind(c)
    enumerator IConcDependenceNone = 0, &
         &     IConcDependenceLinear, &
         &     IConcDependenceLUT, &
         &     IConcDependenceRelativeLinear
  end enum

  !---------------------------------------------------------------------
  ! This derived type describes a correlated k-distribution
  ! representation of an individual gas (including composite gases)
  type ckd_gas_type

    integer :: i_conc_dependence

    ! Molar absorption coefficient in m2 mol-1. If
    ! i_conc_dependence==IConcDependenceNone then it is the absorption
    ! cross section per mole of dry air.  If
    ! conc_dependence==IConcDependenceLinear|IConcDependenceRelativeLinear,
    ! it is the absorption cross section per mole of the gas in
    ! question. It is dimensioned (g_point,pressure,temperature).
    real(jprb), dimension(:,:,:)   :: molar_abs
    
    ! If i_conc_dependence==IConcDependenceLUT then we have an
    ! additional dimension for concentration. It is dimensioned
    ! (g_point,pressure,temperature,conc)
    real(jprb), dimension(:,:,:,:) :: molar_abs_conc

    ! If i_conc_dependence==IConcDependenceRelativeLinear then the
    ! following reference concentration is subtracted from the actual
    ! concentration before the result is multiplied by the mass
    ! absorption coefficient
    real(jprb) :: reference_mole_fraction

    ! Mole fraction coordinate variable if
    ! i_conc_dependence==IConcDependenceLUT
    real(jprb), dimension(:) :: mole_fraction

  end type ckd_gas_type


  type ckd_model_type

    q


  end type ckd_model_type


contains




end module radiation_ecckd
