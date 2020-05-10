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
  use easy_netcdf
  use radiation_gas

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

    ! Code identifying the gas, from the codes in the radiation_gas
    ! module
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


  type ckd_model_type

    ! Gas information

    ! Number of gases
    integer :: ngas = 0
    ! Array of individual-gas data objects
    type(ckd_gas_type), allocatable :: single_gas(:)
    ! Mapping from the "gas codes" in the radiation_gas module to an
    ! index to the single_gas array, where zero means gas not present
    ! (or part of a composite gas)
    integer :: i_gas_mapping(0:NMaxGases)

    ! Coordinates of main look-up table for absorption coeffts

    ! Number of pressure and temperature points
    integer :: npress = 0
    integer :: ntemp  = 0
    ! Natural logarithm of pressure (Pa), dimensioned (npress)
    real(jprb), allocatable :: log_pressure(:)
    ! First temperature profile (K), dimensioned (npress)
    real(jprb), allocatable :: temperature1(:)
    ! Temperature increment (K)
    real(jprb) :: d_temperature

    ! Look-up table for Planck function

    ! Number of entries
    integer :: nplanck = 0
    ! Temperature of first element of look-up table and increment (K)
    real(jprb), allocatable :: temperature1_planck
    real(jprb), allocatable :: d_temperature_planck
    ! Planck function (black body flux into a horizontal plane) in W
    ! m-2, dimensioned (nplanck)
    real(jprb), allocatable :: planck_function(:)

    ! Normalized solar irradiance in each g point dimensioned (ng)
    real(jprb), allocatable :: norm_solar_irradiance(:)

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
    integer,    allocatable :: i_band_number

    ! Shortwave: true, longwave: false
    logical :: is_sw

  contains

    procedure :: read => read_ckd_model
!    procedure :: calc_optical_depth => calc_optical_depth_ckd_model
!    procedure :: deallocate => deallocate_ckd_model

  end type ckd_model_type


contains

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


  subroutine read_ckd_model(this, filename, iverbose)

    use easy_netcdf, only : netcdf_file
    use yomhook,     only : lhook, dr_hook

    class(ckd_model_type), intent(inout) :: this
    character(len=*),      intent(in)    :: filename
    integer, optional,     intent(in)    :: iverbose

    type(netcdf_file) :: file

    real(jprb), allocatable :: temperature_full(:,:)
    real(jprb), allocatable :: temperature_planck(:)

    character(len=512) :: constituent_id

    integer :: iverbose_local

    ! Loop counters
    integer :: jgas, jjgas

    integer :: istart, inext, nchar, i_gas_code

    real(jprb)         :: hook_handle

    if (lhook) call dr_hook('radiation_ecckd:read_ckd_model',0,hook_handle)

    if (present(iverbose)) then
      iverbose_local = iverbose
    else
      iverbose_local = 3
    end if

    call file%open(trim(filename), iverbose=iverbose_local)

    ! Read temperature and pressure coordinate variables
    call file%get('pressure', this%log_pressure)
    this%log_pressure = log(this%log_pressure)
    this%npress = size(this%log_pressure)
    call file%get('temperature', temperature_full)
    allocate(this%temperature1(this%npress));
    this%temperature1 = temperature_full(1,:)
    this%d_temperature = temperature_full(2,1)-temperature_full(1,1)
    deallocate(temperature_full)
    
    ! Read Planck function or solar irradiance
    if (file%exists('solar_irradiance')) then
      this%is_sw = .true.
      call file%get('solar_irradiance', this%norm_solar_irradiance)
      this%norm_solar_irradiance = this%norm_solar_irradiance &
           &  / sum(this%norm_solar_irradiance)
    else
      this%is_sw = .false.
      call file%get('temperature_planck', temperature_planck)
      this%nplanck = size(temperature_planck)
      this%temperature1_planck = temperature_planck(1)
      this%d_temperature_planck = temperature_planck(2) - temperature_planck(1)
      deallocate(temperature_planck)
      call file%get('planck_function', this%planck_function)
    end if

    ! Read spectral mapping of g points
    call file%get('wavenumber1', this%wavenumber1)
    call file%get('wavenumber2', this%wavenumber2)
    this%nwav = size(this%wavenumber1)
    call file%get('gpoint_fraction', this%gpoint_fraction)

    ! Read band information
    call file%get('wavenumber1_band', this%wavenumber1_band)
    call file%get('wavenumber2_band', this%wavenumber2_band)
    this%nband = size(this%wavenumber1_band)
    call file%get('band_number', this%i_band_number)

    ! Read gases
    call file%get('n_gases', this%ngas)
    allocate(this%single_gas(this%ngas))
    call file%get_global_attribute('constituent_id', constituent_id)
    nchar = len(trim(constituent_id))
    istart = 1
    this%i_gas_mapping = 0
    do jgas = 1, this%ngas
      if (jgas < this%ngas) then
        inext = jgas + scan(constituent_id(istart:nchar), ' ')
      else
        inext = nchar+2
      end if
      ! Find gas code
      i_gas_code = 0
      do jjgas = 1, NMaxGases
        if (constituent_id(istart:inext-2) == GasName(jjgas)) then
          i_gas_code = jjgas
          exit
        end if
      end do
      this%i_gas_mapping(i_gas_code) = jgas;
      call this%single_gas(jgas)%read(file, constituent_id(istart:inext-2), i_gas_code)
      istart = inext
    end do
    
    if (lhook) call dr_hook('radiation_ecckd:read_ckd_model',1,hook_handle)

  end subroutine read_ckd_model


end module radiation_ecckd
