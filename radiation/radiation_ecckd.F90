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
  use radiation_gas_constants
  use radiation_ecckd_gas

  implicit none

  public

  type ckd_model_type

    ! Gas information

    ! Number of gases
    integer :: ngas = 0
    ! Array of individual-gas data objects
    type(ckd_gas_type), allocatable :: single_gas(:)
    ! Mapping from the "gas codes" in the radiation_gas_constants
    ! module to an index to the single_gas array, where zero means gas
    ! not present (or part of a composite gas)
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

  !---------------------------------------------------------------------
  ! Read a complete ecCKD gas optics model from a NetCDF file
  ! "filename"
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


  !---------------------------------------------------------------------
  ! Compute layerwise optical depth for each g point for ncol columns
  ! at nlev layers
  subroutine calc_optical_depth_ckd_model(this, ncol, nlev, nmaxgas, &
       &  pressure_hl, temperature_fl, mole_fraction_fl, &
       &  optical_depth_fl, rayleigh_od_fl)

    ! Input variables

    class(ckd_model_type), intent(in)  :: this
    ! Number of columns, levels and input gases
    integer,               intent(in)  :: ncol, nlev, nmaxgas
    ! Pressure at half levels (Pa)
    real(jprb),            intent(in)  :: pressure_hl(nlev+1,ncol)
    ! Temperature at full levels (K)
    real(jprb),            intent(in)  :: temperature_fl(nlev,ncol)
    ! Gas mole fractions at full levels (mol mol-1)
    real(jprb),            intent(in)  :: mole_fraction_fl(nlev,ncol,nmaxgas)
    
    ! Output variables

    ! Layer absorption optical depth for each g point
    real(jprb),            intent(out) :: optical_depth_fl(this%ng,nlev,ncol)
    ! In the shortwave only, the Rayleigh scattering optical depth
    real(jprb),  optional, intent(out) :: rayleigh_od_fl(this%ng,nlev,ncol)

    ! Local variables

    ! Natural logarithm of pressure at full levels
    real(jprb) :: log_pressure_fl(nlev)

    ! Indices and weights in temperature and pressure interpolation
    real(jprb) :: pindex1, tindex1
    real(jprb) :: pweight1, pweight2
    real(jprb) :: tweight1, tweight2
    integer    :: ip1, it1

    integer :: jcol, jlev

    real(jprb)         :: hook_handle

    if (lhook) call dr_hook('radiation_ecckd:calc_optical_depth',0,hook_handle)

    optical_depth_fl = 0.0_jprb

    global_multiplier = 1.0 / (AccelDueToGravity * 0.001_jprb * AirMolarMass)

    do jcol = 1,ncol

      log_pressure_fl = log(0.5 * (pressure_hl(1:nlev,jcol)+pressure_hl(2:nlev+1,jcol)))

      do jlev = 1,nlev
        ! Find interpolation points in pressure
        pindex1 = (log_pressure_fl(jlev)-this%log_pressure1) &
             &    / this%d_log_pressure
        pindex1 = 1.0_jprb + max(0.0_jprb, min(pindex1, this%npress-1.0001))
        ip1 = aint(pindex1)
        pweight2 = pindex1 - ip1
        pweight1 = 1.0_jprb - pweight2

        ! Find interpolation points in temperature
        temperature1 = pweight1*this%temperature1(ip1) &
             &       + pweight2*this%temperature1(ip1+1)
        tindex1 = (temperature_fl(jlev,jcol) - temperature1) &
             &    / this%d_temperature
        tindex1 = 1.0_jprb + max(0.0_jprb, min(tindex1, this%ntemp-1.0001))
        it1 = aint(tindex1)
        tweight2 = tindex1 - it1
        tweight1 = 1.0_jprb - tweight2

        ! Concentration multiplier
        simple_multiplier = global_multiplier &
             &  * (pressure_hl(jlev+1,jcol) - pressure_hl(jlev,jcol))
        
        do jgas = 1,this%ngas

          

        end do

      end do

    end do

    if (lhook) call dr_hook('radiation_ecckd:calc_optical_depth',1,hook_handle)

  end subroutine calc_optical_depth_ckd_model

end module radiation_ecckd
