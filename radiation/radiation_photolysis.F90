! radiation_photolysis.F90 - Derived type containing data to compute photolysis rates
!
! (C) Copyright 2026- ECMWF.
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

module radiation_photolysis

  use parkind1, only : jprb
  
  implicit none
  public

  integer, parameter :: NMaxGasNameLen = 20
  
  !---------------------------------------------------------------------
  ! This derived type contains all the data needed to calculate
  ! photolysis rates from spectral fluxes output from ecRad
  type photolysis_type

    ! Look-up table versus temperature of scaled cross-sections, so
    ! that when the matrix for a particular temperature is extracted,
    ! it can be matrix-vector multiplied by the vector of actinic
    ! fluxes in a vector of ecCKD g-points and will return a vector of
    ! photolysis rates for a set of gases. Dimensioned (ngas,ng,ntemp).
    real(jprb), allocatable :: cross_section_lut(:,:,:)

    ! Store the gas names as an allocatable array of fixed-length strings
    character(len=NMaxGasNameLen), allocatable :: gas_names(:)

    ! Number of temperatures in look-up table
    integer :: ntemperature = 0

    ! Number of g-points in look-up table, which may be fewer than in
    ! the full ecCKD gas-optics model if only a subset are relevant
    ! for photolysis
    integer :: ng = 0

    ! Number of gases whose photolysis rates are to be computed. If a
    ! particular gas has more than one pathway that needs to be
    ! computed separately, it needs an additional entry here
    integer :: ngas = 0

    ! Indices to the first and last g-point of the ecCKD model that is
    ! relevant for photolysis
    integer :: istartg = 0, iendg = 0

    ! First temperature in look-up table and spacing, in Kelvin
    real(jprb) :: temperature1, dtemperature

  contains
    procedure :: configure
    procedure :: calculate
    procedure :: save
    
  end type photolysis_type
  
  
contains

  !---------------------------------------------------------------------
  ! Configure the photolysis_type structure from a netCDF file and a
  ! list of gases to consider
  subroutine configure(this, config, file_name, gases, iverbose)

    use easy_netcdf,     only : netcdf_file
    use radiation_config,only : config_type
    use radiation_io,    only : nulout, radiation_abort
    use radiation_constants, only : PlanckConstant, SpeedOfLight
    use interpolation,   only : interpolate
    use yomhook,         only : lhook, dr_hook, jphook

    class(photolysis_type), intent(inout) :: this
    ! Name of file containing photolysis cross-sections
    type(config_type),      intent(in)    :: config
    ! Configuration structure, only used for the data directory
    character(len=*),       intent(in)    :: file_name
    ! Array of strings containing gases whose photolysis rates are required
    character(len=*),       intent(in)    :: gases(:)
    ! Verbosity level from 1 (least) to 5 (most verbose)
    integer, optional,      intent(in)    :: iverbose
    
    type(netcdf_file) :: file

    ! Temperature for look-up table
    real(jprb), allocatable :: temperature(:) ! K

    ! Photolysis data for one gas as read from file
    real(jprb), allocatable :: wavelength_nm(:)
    real(jprb), allocatable :: cross_section_cm2(:,:)
    real(jprb), allocatable :: quantum_yield(:)
    real(jprb), allocatable :: wavenumber_cm1(:) ! Reverse order

    ! Photolysis data for one gas interpolated to ecCKD wavenumber
    ! grid
    real(jprb), allocatable :: wavenumber_int_cm1(:) ! cm-1
    real(jprb), allocatable :: cross_section_int_cm2(:)
    real(jprb), allocatable :: quantum_yield_int(:)
    real(jprb), allocatable :: solar_photon_flux(:) ! s-1
    real(jprb), allocatable :: photolysis_multiplier(:)
    
    ! Renormalized gpoint_fraction from ecCKD spectral definition
    ! file, dimensioned (nwavn,ng). The original gpoint_fraction sums
    ! to 1 along the wavenumber direction, whereas the new one sums to
    ! one along the g-point direction, thereby stating what fraction
    ! of the solar energy at a particular wavenumber is dealt with by
    ! each g-point.
    real(jprb), allocatable :: gpoint_fraction_renorm(:,:)
    
    ! Number of wavenumbers in ecCKD spectral description
    integer :: nwavn

    ! Number of wavelengths for a single gas in photolysis file
    integer :: nwavl

    ! Loop indices for gases and temperatures
    integer :: jgas, jt

    ! Is the photolysis of the current gas temperature dependent?
    logical :: is_temperature_dependent
    
    integer :: iverbose_local
    
    real(jphook) :: hook_handle
    
    if (lhook) call dr_hook('radiation_photolysis:configure',0,hook_handle)

    if (present(iverbose)) then
      iverbose_local = iverbose
    else
      iverbose_local = 3
    end if

    associate(ckd_model => config%gas_optics_sw)
    
    if (file_name(1:1) == '/' .or. file_name(1:1) == '.') then
      ! Treat file_name as an absolute path
      call file%open(trim(file_name), iverbose=iverbose_local)
    else
      ! Assume the file is in the ecRad data directory
      call file%open(trim(config%directory_name) // '/' // trim(file_name), &
           &         iverbose=iverbose_local)
    end if

    ! Load temperature and store number, start and difference,
    ! assuming the values to be evenly spaced
    call file%get('temperature', temperature)
    this%ntemperature = size(temperature)
    this%temperature1 = temperature(1)
    this%dtemperature = (temperature(this%ntemperature)-temperature(1)) &
         &            / real(this%ntemperature,jprb)

    this%ngas    = size(gases)
    allocate(this%gas_names(this%ngas))
    
    ! Initially assume all g points relevant for photolysis
    this%istartg = 1
    this%iendg   = ckd_model%ng
    this%ng      = ckd_model%ng

    if (iverbose_local >= 2) then
      write(nulout,'(a,i0,a,i0,a,i0,a)') &
           &  'Setting up photolysis calculation for ', &
           &  this%ngas, ' gases and ', this%ng, &
           &  ' spectral g-points as a look-up table with ',  &
           &  this%ntemperature, ' temperatures'
    end if

    allocate(this%cross_section_lut(this%ngas,this%ng,this%ntemperature))

    ! Allocate variables on wavenumber grid
    nwavn = ckd_model%spectral_def%nwav
    allocate(wavenumber_int_cm1(nwavn))
    allocate(quantum_yield_int(nwavn))
    allocate(cross_section_int_cm2(nwavn))
    allocate(solar_photon_flux(nwavn))
    allocate(photolysis_multiplier(nwavn))
    wavenumber_int_cm1 = 0.5_jprb * (ckd_model%spectral_def%wavenumber1 &
         &                          +ckd_model%spectral_def%wavenumber2)
    wavenumber_int_cm1(1)     = ckd_model%spectral_def%wavenumber1(1)
    wavenumber_int_cm1(nwavn) = ckd_model%spectral_def%wavenumber1(nwavn)
    ! 100 converts wavenumber from units of cm-1 to m-1
    solar_photon_flux = ckd_model%spectral_def%solar_spectral_irradiance &
         &   / (PlanckConstant * SpeedOfLight * wavenumber_int_cm1 * 100.0)

    ! Renormalize gpoint_fraction
    allocate(gpoint_fraction_renorm(nwavn,this%ng))
    gpoint_fraction_renorm = ckd_model%spectral_def%gpoint_fraction &
         &  * spread(ckd_model%spectral_def%solar_spectral_irradiance, 2, this%ng)
    gpoint_fraction_renorm = gpoint_fraction_renorm &
         &  * spread(ckd_model%spectral_def%solar_irradiance &
         &           /sum(gpoint_fraction_renorm, 1), 1, nwavn)
    gpoint_fraction_renorm = gpoint_fraction_renorm &
         &  / spread(sum(gpoint_fraction_renorm, 2), 2, this%ng)
    
    ! Loop over requested gases
    do jgas = 1,this%ngas
      this%gas_names(jgas) = trim(gases(jgas))
      ! Load photolysis data for this gas
      call file%get(trim(gases(jgas)) // "_wavelength", wavelength_nm)
      call file%get(trim(gases(jgas)) // "_quantum_yield", quantum_yield)
      call file%get(trim(gases(jgas)) // "_cross_section", cross_section_cm2)

      ! Interpolate on to ecCKD wavenumber grid
      nwavl = size(wavelength_nm)

      allocate(wavenumber_cm1(nwavl))
      wavenumber_cm1 = 1.0e7_jprb/wavelength_nm(nwavl:1:-1)

      if (size(cross_section_cm2,2) == this%ntemperature) then
        is_temperature_dependent = .true.
      else
        is_temperature_dependent = .false.
      end if
      
      if (iverbose_local >= 2) then
        if (is_temperature_dependent) then
          write(nulout,'(a,a,a,f0.1,a,f0.1,a)') '  Temperature dependent photolysis of gas "', &
               &  trim(gases(jgas)), &
               &  '" sensitive to wavenumbers ', wavenumber_cm1(1), '-', &
               &  wavenumber_cm1(nwavl), ' cm-1'
        else
          write(nulout,'(a,a,a,f0.1,a,f0.1,a)') '  Temperature independent photolysis of gas "', &
               &  trim(gases(jgas)), &
               &  '" sensitive to wavenumbers ', wavenumber_cm1(1), '-', &
               &  wavenumber_cm1(nwavl), ' cm-1'
        end if
      end if
      
      if (wavenumber_int_cm1(1) > wavenumber_cm1(1) &
           &  .or. wavenumber_int_cm1(nwavn) < wavenumber_cm1(nwavl)) then
        if (iverbose >= 1) then
          write(nulout,'(a,a,a,f0.1,a,f0.1,a)') '    Warning: photolysis of "', &
               &  trim(gases(jgas)), &
               &  '" sensitive to wavenumbers out of available range ', &
               &  wavenumber_int_cm1(1), '-', wavenumber_int_cm1(nwavn), ' cm-1'
        end if
      end if
      
      ! Interpolate to wavenumber grid
      call interpolate(wavenumber_cm1, quantum_yield(nwavl:1:-1), &
           &           wavenumber_int_cm1, quantum_yield_int, 0.0_jprb)

      if (is_temperature_dependent) then
        do jt = 1,this%ntemperature
          call interpolate(wavenumber_cm1, cross_section_cm2(nwavl:1:-1,jt), &
               &           wavenumber_int_cm1, cross_section_int_cm2, 0.0_jprb)
          ! 1e-4 converts cross section from cm2 to m2
          photolysis_multiplier = 1.0e-4 * cross_section_int_cm2 * solar_photon_flux;
          this%cross_section_lut(jgas,:,jt) &
               &  = matmul(photolysis_multiplier, gpoint_fraction_renorm) &
               &  / ckd_model%spectral_def%solar_irradiance
        end do
      else
        call interpolate(wavenumber_cm1, cross_section_cm2(nwavl:1:-1,1), &
             &           wavenumber_int_cm1, cross_section_int_cm2, 0.0_jprb)
        photolysis_multiplier = 1.0e-4 * cross_section_int_cm2 * solar_photon_flux;
        this%cross_section_lut(jgas,:,1) &
             &  = matmul(photolysis_multiplier, gpoint_fraction_renorm) &
             &  / ckd_model%spectral_def%solar_irradiance
        this%cross_section_lut(jgas,:,2:this%ntemperature) &
             &  = spread(this%cross_section_lut(jgas,:,1),2,this%ntemperature-1)
      end if
      
      deallocate(wavelength_nm)
      deallocate(wavenumber_cm1)
      deallocate(quantum_yield)
      deallocate(cross_section_cm2)
    end do
    
    call file%close()

    end associate
    
    if (lhook) call dr_hook('radiation_photolysis:configure',1,hook_handle)

  end subroutine configure

  
  !---------------------------------------------------------------------
  ! Calculate photolysis rates from spectral fluxes
  subroutine calculate(this, icol, mu0, temperature_hl, flux, rates, ilay1, ilay2)

    use radiation_flux,  only : flux_type
    use radiation_io,    only : nulerr, radiation_abort
    use yomhook,         only : lhook, dr_hook, jphook

    class(photolysis_type), intent(in)  :: this
    ! Column number
    integer,                     intent(in)  :: icol
    ! Cosine of solar zenith angle
    real(jprb),                  intent(in)  :: mu0 
    ! Half-level temperature (K)
    real(jprb),                  intent(in)  :: temperature_hl(:) ! (nlay+1)
    ! Structure containing spectral fluxes from ecRad
    type(flux_type),             intent(in)  :: flux
    ! Output photodissociation rates (s-1)
    real(jprb),                  intent(out) :: rates(:,:) ! (ngas,nlay)
    ! Optional range of layers to process    
    integer, optional,           intent(in)  :: ilay1, ilay2

    ! Actinic flux at a half-level
    real(jprb) :: actinic_flux(this%ng)

    ! Local cross-sections
    real(jprb) :: cross_section(this%ngas,this%ng)

    ! Rates at half-levels
    real(jprb), allocatable :: rates_hl(:,:)
    
    integer :: il1, il2

    integer :: nlay

    integer :: jlev, jlay, jgas

    ! Index
    integer :: itemp
    real(jprb) :: wtemp
    
    real(jphook) :: hook_handle
    
    if (lhook) call dr_hook('radiation_photolysis:calculate',0,hook_handle)

    ! Checks
    if (.not. allocated(flux%sw_dn_direct_band)) then
      write(nulerr, '(a)') '*** Error: spectral shortwave radiative fluxes not output by ecRad for photolysis'
      call radiation_abort('Radiation configuration error')
    else if (size(flux%sw_dn_direct_band,1) /= this%ng) then
      write(nulerr,'(a,i0,a,i0)') '*** Error: photolysis expects spectral fluxes at ', this%ng, &
           &  ' g-points, got ', size(flux%sw_dn_direct_band,1)
      call radiation_abort('Radiation configuration error')
    end if
      
    ! Check the sun is above the horizon
    if (mu0 > 0.0_jprb) then
    
      if (present(ilay1)) then
        il1 = ilay1
      else
        il1 = 1
      end if

      if (present(ilay2)) then
        il2 = ilay2
      else
        il2 = size(flux%sw_dn_band,3)-1
      end if

      nlay = il2-il1+1

      allocate(rates_hl(this%ngas,nlay+1))
      
      ! Loop over half-levels
      do jlev = il1,il2+1
        ! Assume the diffuse fluxes are isotropic in each hemisphere so
        ! the diffuse fluxes into a horizontal plane are multiplied by 2
        ! and the direct flux into a horizontal plane is scaled by 1/mu0
        actinic_flux = flux%sw_dn_direct_band(this%istartg:this%iendg,icol,jlev) &
             &            * (1.0_jprb/mu0 - 2.0_jprb) &
             &  + 2.0_jprb * (flux%sw_dn_band(this%istartg:this%iendg,icol,jlev) &
             &               +flux%sw_up_band(this%istartg:this%iendg,icol,jlev))
        ! Interpolation points and weights
        wtemp = 1.0_jprb + (temperature_hl(jlev) - this%temperature1) / this%dtemperature
        if (wtemp < 1.0_jprb) then
          itemp = 1
          wtemp = 0.0_jprb
        elseif (wtemp >= this%ntemperature) then
          itemp = this%ntemperature-1
          wtemp = 1.0_jprb
        else
          itemp = int(wtemp)
          wtemp = wtemp - itemp
        end if
        ! Interpolate cross sections
        cross_section = (1.0_jprb - wtemp) * this%cross_section_lut(:,:,itemp) &
             &        +             wtemp  * this%cross_section_lut(:,:,itemp+1)
        
        rates_hl(:,jlev) = matmul(cross_section, actinic_flux)
      end do

      ! Loop over full levels
      do jlay = il1,il2
        ! Assume the photolysis rates for each gas vary exponentially within each layer
        do jgas = 1,this%ngas
          if (rates_hl(jgas,jlay+1) > 0.99_jprb * rates_hl(jgas,jlay)) then
            ! Optically thin layer: very small vertical variation of
            ! rates so take average of values at top and bottom
            ! of layer
            rates(jgas,jlay) = 0.5_jprb * (rates_hl(jgas,jlay) + rates_hl(jgas,jlay+1))
          elseif (rates_hl(jgas,jlay) <= 0.0_jprb) then
            ! No flux
            rates(jgas,jlay) = 0.0_jprb
          else
            ! Assume an exponential variation of actinic flux through
            ! the layer and calculate the layer-mean value
            rates(jgas,jlay) = (rates_hl(jgas,jlay+1) - rates_hl(jgas,jlay)) &
                 &      / log(max(rates_hl(jgas,jlay+1)/rates_hl(jgas,jlay),tiny(1.0_jprb)))
          end if
        end do
      end do
      
    else
      ! Sun below horizon
      rates(:,:) = 0.0_jprb
    end if
    
    if (lhook) call dr_hook('radiation_photolysis:calculate',1,hook_handle)
    
  end subroutine calculate

  
  !---------------------------------------------------------------------
  ! Save computed photolysis rates to a netCDF file
  subroutine save(this, file_name, rates, iverbose)

    use easy_netcdf,     only : netcdf_file
    !use radiation_io,    only : nulout, nulerr, radiation_abort
    use yomhook,         only : lhook, dr_hook, jphook

    class(photolysis_type), intent(inout) :: this
    ! Name of file containing photolysis cross-sections
    character(len=*),            intent(in)    :: file_name
    ! Photolysis rates (s-1) dimensioned (ngas,nlay,ncol)
    real(kind=jprb), allocatable :: rates(:,:,:)
    ! Verbosity level from 1 (least) to 5 (most verbose)
    integer, optional,     intent(in)    :: iverbose

    ! Object for output NetCDF file
    type(netcdf_file) :: out_file

    integer :: nlev, ncol
    integer :: jgas
    integer :: i_local_verbose

    real(jphook) :: hook_handle

    if (lhook) call dr_hook('radiation_photolysis:save',0,hook_handle)

    if (present(iverbose)) then
      i_local_verbose = iverbose
    else
      i_local_verbose = 3
    end if

    ! Open the file
    call out_file%create(trim(file_name), iverbose=i_local_verbose)

    ! Define dimensions
    ncol = size(rates,3)
    nlev = size(rates,2)
    
    call out_file%define_dimension("column", ncol)
    call out_file%define_dimension("level",  nlev)

    ! Put global attributes
    call out_file%put_global_attributes( &
         &   title_str="Photolysis rates computed from the ecRad offline radiation model", &
         &   references_str="Hogan, R. J., and A. Bozzo, 2018: A flexible and efficient radiation " &
         &   //"scheme for the ECMWF model. J. Adv. Modeling Earth Sys., 10, 1990–2008", &
         &   source_str="ecRad offline radiation model")

    do jgas = 1,this%ngas
      call out_file%define_variable(trim(this%gas_names(jgas)) // "_photolysis_rate", &
           &  long_name=trim(this%gas_names(jgas)) // " photolysis rate", units_str="s-1", &
           &  dim2_name="column", dim1_name="level")
    end do

    do jgas = 1,this%ngas
      call out_file%put(trim(this%gas_names(jgas)) // "_photolysis_rate", rates(jgas,:,:))
    end do

    call out_file%close()

    if (lhook) call dr_hook('radiation_photolysis:save',1,hook_handle)

  end subroutine save
    
end module radiation_photolysis

