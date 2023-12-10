! radiation_general_cloud_optics_data.F90 - Type to store generalized cloud optical properties
!
! (C) Copyright 2019- ECMWF.
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
! License: see the COPYING file for details
!

module radiation_general_cloud_optics_data

  use parkind1, only : jprb

  implicit none

  public

#ifdef FLOTSAM
#include <flotsam.inc>
#endif

  ! Enumerators for defining how the phase function shape is to be
  ! characterized: phase function only, Legendre polynomials or
  ! FLOTSAM decomposition
  enum, bind(c) 
     enumerator IPhaseFuncAsymmetry, IPhaseFuncLegendre, IPhaseFuncFLOTSAM
  end enum

  !---------------------------------------------------------------------
  ! This type holds the configuration information to compute optical
  ! properties for a particular type of cloud or hydrometeor in one of
  ! the shortwave or longwave
  type general_cloud_optics_type
    ! Band-specific (or g-point-specific) values as a look-up table
    ! versus effective radius dimensioned (nband,n_effective_radius)

    ! Extinction coefficient per unit mass (m2 kg-1)
    real(jprb), allocatable, dimension(:,:) :: mass_ext
 
    ! Single-scattering albedo (dimensionless)
    real(jprb), allocatable, dimension(:,:) :: ssa

    ! Phase function description for DISORT, dimensioned
    ! (nband,n_effective_radius,ncomponent), or just asymmetry factor
    real(jprb), allocatable, dimension(:,:,:) :: pf
    
    ! Support for RTTOV-style look-up table using water content and
    ! temperature, rather than effective radius; here the look-up
    ! tables are dimensioned (nband, n_water_content, n_temperature)
    real(jprb), allocatable, dimension(:,:,:) &
         &  :: mass_ext_wct, ssa_wct, asymmetry_wct
    
#ifdef FLOTSAM
    ! Phase function description for FLOTSAM
    
    ! Scattering angle (radians)
    real(jprb), allocatable, dimension(:)     :: scattering_angle
    ! Phase function, dimensioned (nangle,nband,n_effective_radius) 
    real(jprb), allocatable, dimension(:,:,:) :: phase_function
    ! Phase function components, dimensioned (ncomp,nband,n_effective_radius) 
    real(jprb), allocatable, dimension(:,:,:) :: phase_function_components
#endif

    ! Number of effective radius coefficients, start value and
    ! interval (m) in look-up table
    integer    :: n_effective_radius = 0
    real(jprb) :: effective_radius_0, d_effective_radius

    ! Number of water content coefficients, natural log of start value
    ! (kg m-3) and interval in look-up table
    integer    :: n_water_content = 0
    real(jprb) :: ln_water_content_0, d_ln_water_content

    ! Number of temperature coefficients, start value and interval (K)
    ! in look-up table
    integer    :: n_temperature = 0
    real(jprb) :: temperature_0, d_temperature
    
    ! Name of cloud/precip type and scattering model
    ! (e.g. "mie_droplet", "fu-muskatel_ice"). These are used to
    ! generate the name of the data file from which the coefficients
    ! are read.
    character(len=511) :: type_name

    ! Phase function mode (asymmetry factor, Legendre or FLOTSAM)
    integer :: pf_mode = IPhaseFuncAsymmetry
    
    ! Number of phase-function components
    integer :: n_pf_components = 1
    
    ! Do we use bands or g-points?
    logical :: use_bands = .false.

    ! Use effective radius, or water content and temperature?
    logical :: use_effective_radius = .true.

   contains
     procedure :: setup => setup_general_cloud_optics
     procedure :: setup_wct => setup_general_cloud_optics_wct
     procedure :: add_optical_properties
#ifdef FLOTSAM
     procedure :: add_optical_properties_flotsam
#endif
     procedure :: save => save_general_cloud_optics_data

  end type general_cloud_optics_type

contains

  ! Provides elemental function "delta_eddington"
#include "radiation_delta_eddington.h"

  !---------------------------------------------------------------------
  ! Setup cloud optics coefficients by reading them from a file
  subroutine setup_general_cloud_optics(this, file_name, specdef, &
       &                                use_bands, use_thick_averaging, &
       &                                weighting_temperature, &
       &                                iverbose, pf_mode, n_pf_components)

    use yomhook,                       only : lhook, dr_hook, jphook
    use easy_netcdf,                   only : netcdf_file
    use radiation_constants,           only : Pi
    use radiation_io,                  only : nulout, nulerr, radiation_abort
    use radiation_constants,           only : Pi
    use radiation_spectral_definition, only : spectral_definition_type

    class(general_cloud_optics_type), intent(inout), target :: this
    character(len=*), intent(in)               :: file_name
    type(spectral_definition_type), intent(in) :: specdef
    logical, intent(in), optional              :: use_bands, use_thick_averaging
    real(jprb), intent(in), optional           :: weighting_temperature ! K
    integer, intent(in), optional              :: iverbose
    integer, intent(in), optional              :: pf_mode
    integer, intent(in), optional              :: n_pf_components

    ! Spectral properties read from file, dimensioned (wavenumber,
    ! n_effective_radius)
    real(jprb), dimension(:,:), allocatable :: mass_ext, & ! m2 kg-1
         &                                     ssa, asymmetry

    ! Asymmetry factor after averaging, dimensioned (band,
    ! n_effective_radius)
    real(jprb), dimension(:,:), allocatable :: asymmetry_band
    
    ! Reflectance of an infinitely thick cloud, needed for thick
    ! averaging
    real(jprb), dimension(:,:), allocatable :: ref_inf

    ! Coordinate variables from file
    real(jprb), dimension(:), allocatable :: wavenumber       ! cm-1
    real(jprb), dimension(:), allocatable :: effective_radius ! m

    ! Matrix mapping optical properties in the file to values per
    ! g-point or band, such that in the thin-averaging case,
    ! this%mass_ext=matmul(mapping,file%mass_ext), so mapping is
    ! dimensioned (ngpoint,nwav)
    real(jprb), dimension(:,:), allocatable :: mapping

    ! Scattering phase function and scattering angles
    real(jprb), dimension(:,:,:), allocatable :: phase_function, phase_function_band
    real(jprb), dimension(:), allocatable :: scattering_angle

    ! The NetCDF file containing the coefficients
    type(netcdf_file)  :: file

    real(jprb) :: diff_spread
    integer    :: iverb
    integer    :: nre  ! Number of effective radii
    integer    :: nwav ! Number of wavenumbers describing cloud
    integer    :: nang ! Number of angles describing phase function
    integer    :: nband! Number of bands (or g points)
    integer    :: n_pf_components_flotsam ! Number of phase-function components for FLOTSAM
    integer    :: istatus ! Return code from flotsam calls
    integer    :: jband, jwav, jang, jre ! Loop indices

    logical    :: use_bands_local, use_thick_averaging_local

    real(jphook) :: hook_handle

    if (lhook) call dr_hook('radiation_general_cloud_optics_data:setup',0,hook_handle)

    ! Set local values of optional inputs
    if (present(iverbose)) then
      iverb = iverbose
    else
      iverb = 2
    end if

    if (present(use_bands)) then
      use_bands_local = use_bands
    else
      use_bands_local = .false.
    end if

    if (present(use_thick_averaging)) then
      use_thick_averaging_local = use_thick_averaging
    else
      use_thick_averaging_local = .false.
    end if

    ! Open the scattering file and configure the way it is read
    call file%open(trim(file_name), iverbose=iverb)
    !call file%transpose_matrices()

    if (file%exists('water_content')) then
      ! File assumed to provide optical properties as a function of
      ! water content and temperature; call alternative routine, close
      ! the file and return
      call this%setup_wct(file_name, file, specdef, use_bands_local, &
           &              use_thick_averaging_local, weighting_temperature, iverb)
      call file%close()
      return
    end if

    ! Read coordinate variables
    call file%get('wavenumber', wavenumber)
    call file%get('effective_radius', effective_radius)

    ! Read the band-specific coefficients
    call file%get('mass_extinction_coefficient', mass_ext)
    call file%get('single_scattering_albedo', ssa)

    if (present(pf_mode)) then
      this%pf_mode = pf_mode
      this%n_pf_components = n_pf_components
      if (pf_mode /= IPhaseFuncAsymmetry) then
        if (file%exists('phase_function')) then
          call file%get('scattering_angle', scattering_angle)
          scattering_angle = scattering_angle * (Pi/180.0_jprb)
          nang  = size(scattering_angle)
          call file%get('phase_function',   phase_function)
        else
          write(nulerr, '(a,a)') '*** Error: phase_function not present in ', &
               &  trim(file_name)
          call radiation_abort('Radiation configuration error')
        end if
      end if
    end if

    call file%get('asymmetry_factor', asymmetry)

#ifdef FLOTSAM
    if (file%exists('phase_function')) then
      call file%get('scattering_angle', this%scattering_angle)
      this%scattering_angle = this%scattering_angle * (Pi/180.0_jprb)
      call file%get('phase_function',   phase_function)      
    end if
#endif

    ! Close scattering file
    call file%close()

    ! Check effective radius is evenly spaced
    nre = size(effective_radius)
    ! Fractional range of differences, should be near zero for evenly
    ! spaced data
    diff_spread = (maxval(effective_radius(2:nre)-effective_radius(1:nre-1))  &
         &        -minval(effective_radius(2:nre)-effective_radius(1:nre-1))) &
         &      /  minval(abs(effective_radius(2:nre)-effective_radius(1:nre-1)))
    if (diff_spread > 0.01_jprb) then
      write(nulerr, '(a,a,a)') '*** Error: effective_radius in ', &
           &  trim(file_name), ', is not evenly spaced to 1%'
      call radiation_abort('Radiation configuration error')
    end if

    ! Set up effective radius coordinate variable
    this%n_effective_radius = nre
    this%effective_radius_0 = effective_radius(1)
    this%d_effective_radius = effective_radius(2) - effective_radius(1)

    nwav = size(wavenumber)

    ! Define the mapping matrix
    call specdef%calc_mapping(wavenumber, mapping, &
         &  weighting_temperature=weighting_temperature, use_bands=use_bands)

    ! Thick averaging should be performed on delta-Eddington scaled
    ! quantities (it makes no difference to thin averaging)
    call delta_eddington(mass_ext, ssa, asymmetry)

    ! Thin averaging
    this%mass_ext  = matmul(mapping, mass_ext)
    this%ssa       = matmul(mapping, mass_ext*ssa) / this%mass_ext
    nband = size(this%mass_ext,1)
    allocate(this%pf(nband,nre,this%n_pf_components))
    
    asymmetry_band = matmul(mapping, mass_ext*ssa*asymmetry) / (this%mass_ext*this%ssa)
    
    if (this%pf_mode /= IPhaseFuncAsymmetry) then
      ! Spectral averaging of phase function, weighted by scattering
      ! coefficient
      allocate(phase_function_band(nang,nband,this%n_effective_radius))
      phase_function_band = 0.0_jprb
      do jre = 1,nre
        do jwav = 1,nwav
          do jband = 1,nband
            phase_function_band(:,jband,jre) = phase_function_band(:,jband,jre) &
                 &  + mapping(jband,jwav) * mass_ext(jwav,jre) * ssa(jwav,jre) &
                 &    * phase_function(:,jwav,jre)
          end do
        end do
      end do
      do jang = 1,nang
        phase_function_band(jang,:,:) = phase_function_band(jang,:,:) / (this%mass_ext*this%ssa)
      end do

!#define OVERRIDE_PHASE_FUNCTION 1
#ifndef OVERRIDE_PHASE_FUNCTION
      call legendre_decomposition(nang, scattering_angle, phase_function_band, &
           &                      n_pf_components, this%pf)
#else
      ! Provide instead a simple delta-Eddington two-stream phase
      ! function
      !this%pf = 0.0_jprb
      !this%pf(:,:,1) = asymmetry_band ! Already delta-Eddington'ed

      ! Remove Delta-Eddington, then fix the subsequent terms to be
      ! equal to g^2 so that all numbers of streams result in
      ! delta-Eddington
      this%pf(:,:,1) = asymmetry_band/(1.0_jprb-asymmetry_band)
      this%pf(:,:,2) = this%pf(:,:,1)*this%pf(:,:,1)
      do jang = 2,n_pf_components
         this%pf(:,:,jang) = this%pf(:,:,2)
      end do

      ! Remove Delta-Eddington
      !this%pf(:,:,1) = asymmetry_band/(1.0_jprb-asymmetry_band)
      !do jang = 2,n_pf_components
        ! Henyey-Greenstein
      !  this%pf(:,:,jang) = this%pf(:,:,jang-1)*this%pf(:,:,1)
      !end do
#endif
      if (use_thick_averaging_local) then
        write(nulout,'(a)') 'Warning: thick averaging not performed when using Legendre decomposition' 
      end if
    else if (use_thick_averaging_local) then
      ! Thick averaging as described by Edwards and Slingo (1996),
      ! modifying only the single-scattering albedo
      allocate(ref_inf(nwav, nre))

      ! Eqs. 18 and 17 of Edwards & Slingo (1996)
      ref_inf = sqrt((1.0_jprb - ssa) / (1.0_jprb - ssa*asymmetry))
      ref_inf = (1.0_jprb - ref_inf) / (1.0_jprb + ref_inf)
      ! Here the left-hand side is actually the averaged ref_inf
      this%ssa = matmul(mapping, ref_inf)
      ! Eq. 19 of Edwards and Slingo (1996)
      this%ssa = 4.0_jprb * this%ssa / ((1.0_jprb + this%ssa)**2 &
           &  - asymmetry_band * (1.0_jprb - this%ssa)**2)

      deallocate(ref_inf)
    end if

#ifdef FLOTSAM
    nang  = size(this%scattering_angle)
    nband = size(this%mass_ext,1)
    if (allocated(phase_function)) then

! Optionally extract phase-function components for all wavenumbers in
! file
!#define ANALYSE_ALL_PHASE_FUNCTIONS 1
#ifdef ANALYSE_ALL_PHASE_FUNCTIONS
      allocate(this%phase_function(nang,nwav,this%n_effective_radius))
      n_pf_components_flotsam = flotsam_n_phase_function_components()
      allocate(this%phase_function_components(n_pf_components_flotsam,nwav,this%n_effective_radius))
      istatus = flotsam_analyse_phase_functions(nwav*this%n_effective_radius, nang, &
           &  this%scattering_angle, phase_function, 4.0_jprb*Pi, &
           &  this%phase_function, this%phase_function_components)
      do jband = 1,this%n_effective_radius
        do jwav = 1,nwav
          write(103,*) effective_radius(jband), wavenumber(jwav), ssa(jwav,jband), &
               &  asymmetry(jwav,jband), this%phase_function_components(:,jwav,jband)
        end do
      end do
      deallocate(this%phase_function)
      deallocate(this%phase_function_components)
#endif

      allocate(phase_function_band(nang,nband,this%n_effective_radius))
      phase_function_band = 0.0_jprb
      do jwav = 1,nwav
        do jband = 1,nband
          phase_function_band(:,jband,:) = phase_function_band(:,jband,:) &
               &  + mapping(jband,jwav) * phase_function(:,jwav,:)
        end do
      end do
      allocate(this%phase_function(nang,nband,this%n_effective_radius))
      n_pf_components_flotsam = flotsam_n_phase_function_components()
      allocate(this%phase_function_components(n_pf_components_flotsam,nband,this%n_effective_radius))
      istatus = flotsam_analyse_phase_functions(nband*this%n_effective_radius, nang, &
           &  this%scattering_angle, phase_function_band, 4.0_jprb*Pi, &
           &  this%phase_function, this%phase_function_components)
    end if
#endif

    deallocate(mapping)

    ! Revert back to unscaled quantities
    call revert_delta_eddington(this%mass_ext, this%ssa, asymmetry_band)

    if (this%pf_mode == IPhaseFuncAsymmetry) then
      this%pf(:,:,1) = asymmetry_band
    end if
    
    deallocate(asymmetry_band)  
    
    if (iverb >= 2) then
      write(nulout,'(a,a)') '  File: ', trim(file_name)
      if (present(weighting_temperature)) then
        write(nulout,'(a,f7.1,a)') '  Weighting temperature: ', weighting_temperature, ' K'
      else
        write(nulout,'(a,f7.1,a)') '  Weighting temperature: ', specdef%reference_temperature, ' K'
      end if
      if (use_thick_averaging_local) then
        write(nulout,'(a)') '  SSA averaging: optically thick limit'
      else
        write(nulout,'(a)') '  SSA averaging: optically thin limit'
      end if
      if (use_bands_local) then
        write(nulout,'(a,i0,a)') '  Spectral discretization: ', specdef%nband, ' bands'
      else
        write(nulout,'(a,i0,a)') '  Spectral discretization: ', specdef%ng, ' g-points'
      end if
      write(nulout,'(a,i0,a,f6.1,a,f6.1,a)') '  Effective radius look-up: ', nre, ' points in range ', &
           &  effective_radius(1)*1.0e6_jprb, '-', effective_radius(nre)*1.0e6_jprb, ' um'
      write(nulout,'(a,i0,a,i0,a)') '  Wavenumber range: ', int(specdef%min_wavenumber()), '-', &
           &  int(specdef%max_wavenumber()), ' cm-1'
#ifdef FLOTSAM
      if (allocated(this%phase_function)) then
        write(nulout,'(a,i0)') '  Phase function angles for FLOTSAM: ', nang
      end if
#endif
      if (this%pf_mode == IPhaseFuncAsymmetry) then
        write(nulout,'(a)') '  Phase function represented by asymmetry factor'
      else if (this%pf_mode == IPhaseFuncLegendre) then
        write(nulout,'(a,i0,a)') '  Phase function represented by ', this%n_pf_components, ' Legendre terms'
      end if
    end if

    if (lhook) call dr_hook('radiation_general_cloud_optics_data:setup',1,hook_handle)

  end subroutine setup_general_cloud_optics


  !---------------------------------------------------------------------
  ! Setup cloud optics coefficients as a function of water content and
  ! temperature
  subroutine setup_general_cloud_optics_wct(this, file_name, file, specdef, &
       &                                    use_bands, use_thick_averaging, &
       &                                    weighting_temperature, &
       &                                    iverbose)

    use yomhook,                       only : lhook, dr_hook, jphook
    use easy_netcdf,                   only : netcdf_file
    use radiation_spectral_definition, only : spectral_definition_type
    use radiation_io,                  only : nulout, nulerr, radiation_abort
    use radiation_constants,           only : Pi

    class(general_cloud_optics_type), intent(inout)    :: this
    ! The NetCDF file containing the coefficients
    character(len=*), intent(in)               :: file_name
    type(netcdf_file), intent(inout)           :: file
    type(spectral_definition_type), intent(in) :: specdef
    logical, intent(in)                        :: use_bands, use_thick_averaging
    real(jprb), intent(in)                     :: weighting_temperature ! K
    integer, intent(in)                        :: iverbose
    
    ! Spectral properties read from file, dimensioned (wavenumber,
    ! water_content, temperature)
    real(jprb), dimension(:,:,:), allocatable :: mass_ext, & ! m2 kg-1
         &                                       ssa, asymmetry

    ! Reflectance of an infinitely thick cloud, needed for thick
    ! averaging
    real(jprb), dimension(:,:,:), allocatable :: ref_inf

    ! Coordinate variables from file
    real(jprb), dimension(:), allocatable :: wavenumber       ! cm-1
    real(jprb), dimension(:), allocatable :: water_content    ! m
    real(jprb), dimension(:), allocatable :: temperature      ! m
    
    ! Matrix mapping optical properties in the file to values per
    ! g-point or band, such that in the thin-averaging case,
    ! this%mass_ext=matmul(mapping,file%mass_ext), so mapping is
    ! dimensioned (ngpoint,nwav)
    real(jprb), dimension(:,:), allocatable :: mapping

    real(jprb) :: diff_spread
    integer    :: iverb
    integer    :: nwc  ! Number of water contents
    integer    :: ntemp! Number of temperatures
    integer    :: nwav ! Number of wavenumbers describing cloud
    integer    :: nband! Number of bands (or g points)

    integer    :: jband, jwav, jt ! Loop indices

    real(jphook) :: hook_handle

    if (lhook) call dr_hook('radiation_general_cloud_optics_data:setup_wct',0,hook_handle)

    ! Read coordinate variables
    call file%get('wavenumber', wavenumber)
    call file%get('water_content', water_content)
    call file%get('temperature', temperature)
    
    ! Read the band-specific coefficients
    call file%get('mass_extinction_coefficient', mass_ext)
    call file%get('single_scattering_albedo', ssa)
    call file%get('asymmetry_factor', asymmetry)

    this%use_effective_radius = .false.

    ! Check temperature is evenly spaced
    ntemp = size(temperature)
    ! Fractional range of differences, should be near zero for evenly
    ! spaced data
    diff_spread = (maxval(temperature(2:ntemp)-temperature(1:ntemp-1))  &
         &        -minval(temperature(2:ntemp)-temperature(1:ntemp-1))) &
         &      /  minval(abs(temperature(2:ntemp)-temperature(1:ntemp-1)))
    if (diff_spread > 0.01_jprb) then
      write(nulerr, '(a,a,a)') '*** Error: temperature in ', &
           &  trim(file_name), ', is not evenly spaced to 1%'
      call radiation_abort('Radiation configuration error')
    end if

    this%n_temperature = ntemp
    this%temperature_0 = temperature(1)
    this%d_temperature = temperature(2) - temperature(1)
    
    ! Check water content is logarithmically spaced
    nwc = size(water_content)
    ! Fractional range of differences, should be near zero for evenly
    ! spaced data
    diff_spread = (maxval(log(water_content(2:nwc))-log(water_content(1:nwc-1))) &
         &        -minval(log(water_content(2:nwc))-log(water_content(1:nwc-1)))) &
         &      /  minval(abs(log(water_content(2:nwc))-log(water_content(1:nwc-1))))
    if (diff_spread > 0.01_jprb) then
      write(nulerr, '(a,a,a)') '*** Error: log(water_content) in ', &
           &  trim(file_name), ', is not evenly spaced to 1%'
      call radiation_abort('Radiation configuration error')
    end if

    this%n_water_content = nwc
    this%ln_water_content_0 = log(water_content(1))
    this%d_ln_water_content = log(water_content(2)/water_content(1))

    nwav = size(wavenumber)

    ! Define the mapping matrix
    call specdef%calc_mapping(wavenumber, mapping, &
         &  weighting_temperature=weighting_temperature, use_bands=use_bands)

    ! Thick averaging should be performed on delta-Eddington scaled
    ! quantities (it makes no difference to thin averaging)
    call delta_eddington(mass_ext, ssa, asymmetry)

    ! Thin averaging
    nband = size(mapping,1)
    allocate(this%mass_ext_wct(nband,nwc,ntemp))
    allocate(this%ssa_wct(nband,nwc,ntemp))
    allocate(this%asymmetry_wct(nband,nwc,ntemp))
    do jt = 1,ntemp
      this%mass_ext_wct(:,:,jt)  = matmul(mapping, mass_ext(:,:,jt))
      this%ssa_wct(:,:,jt)       = matmul(mapping, ssa(:,:,jt))
      this%asymmetry_wct(:,:,jt) = matmul(mapping, asymmetry(:,:,jt))
    end do
    
    if (use_thick_averaging) then
      ! Thick averaging as described by Edwards and Slingo (1996),
      ! modifying only the single-scattering albedo
      allocate(ref_inf(nwav, nwc, ntemp))

      ! Eqs. 18 and 17 of Edwards & Slingo (1996)
      ref_inf = sqrt((1.0_jprb - ssa) / (1.0_jprb - ssa*asymmetry))
      ref_inf = (1.0_jprb - ref_inf) / (1.0_jprb + ref_inf)
      ! Here the left-hand side is actually the averaged ref_inf
      
      do jt = 1,ntemp
        this%ssa_wct(:,:,jt) = matmul(mapping, ref_inf(:,:,jt))
      end do

      ! Eq. 19 of Edwards and Slingo (1996)
      this%ssa_wct = 4.0_jprb * this%ssa_wct / ((1.0_jprb + this%ssa_wct)**2 &
           &  - this%asymmetry_wct * (1.0_jprb - this%ssa_wct)**2)

      deallocate(ref_inf)
    end if

    deallocate(mapping)

    ! Revert back to unscaled quantities
    call revert_delta_eddington(this%mass_ext_wct, this%ssa_wct, this%asymmetry_wct)

    if (iverb >= 2) then
      write(nulout,'(a,a)') '  File: ', trim(file_name)
      write(nulout,'(a,f7.1,a)') '  Weighting temperature: ', weighting_temperature, ' K'
      if (use_thick_averaging) then
        write(nulout,'(a)') '  SSA averaging: optically thick limit'
      else
        write(nulout,'(a)') '  SSA averaging: optically thin limit'
      end if
      if (use_bands) then
        write(nulout,'(a,i0,a)') '  Spectral discretization: ', specdef%nband, ' bands'
      else
        write(nulout,'(a,i0,a)') '  Spectral discretization: ', specdef%ng, ' g-points'
      end if
      write(nulout,'(a,i0,a,f6.1,a,f6.1,a)') '  Water-content look-up: ', nwc, ' log-spaced points in range ', &
           &  water_content(1), '-', water_content(nwc), ' kg m-3'
      write(nulout,'(a,i0,a,f6.1,a,f6.1,a)') '  Temperature look-up: ', ntemp, ' points in range ', &
           &  temperature(1), '-', temperature(ntemp), ' K'
      
      write(nulout,'(a,i0,a,i0,a)') '  Wavenumber range: ', int(specdef%min_wavenumber()), '-', &
           &  int(specdef%max_wavenumber()), ' cm-1'
    end if
    
    if (lhook) call dr_hook('radiation_general_cloud_optics_data:setup_wct',1,hook_handle)

  end subroutine setup_general_cloud_optics_wct
  
    
  !---------------------------------------------------------------------
  ! Add the optical properties of a particular cloud type to the
  ! accumulated optical properties of all cloud types
  subroutine add_optical_properties(this, ng, nlev, ncol, &
       &                            cloud_fraction, &
       &                            water_path, effective_radius, &
       &                            od, scat_od, scat_pf, &
       &                            layer_depth, temperature_fl)

    use yomhook, only : lhook, dr_hook, jphook

    class(general_cloud_optics_type), intent(in) :: this

    ! Number of g points, levels and columns
    integer, intent(in) :: ng, nlev, ncol

    ! Properties of present cloud type, dimensioned (ncol,nlev)
    real(jprb), intent(in) :: cloud_fraction(:,:)
    real(jprb), intent(in) :: water_path(:,:)       ! kg m-2
    real(jprb), intent(in) :: effective_radius(:,:) ! m

    ! Optical properties which are additive per cloud type,
    ! dimensioned (ng,nlev,ncol)
    real(jprb), intent(inout), dimension(ng,nlev,ncol) &
         &  :: od             ! Optical depth of layer
    real(jprb), intent(inout), dimension(ng,nlev,ncol), optional &
         &  :: scat_od        ! Scattering optical depth of layer
    real(jprb), intent(inout), dimension(ng,nlev,ncol,this%n_pf_components), optional &
         &  :: scat_pf        ! Scattering optical depth x phase function component

    ! Layer thickness (m), full-level temperature (K)
    real(jprb), intent(in), optional :: layer_depth(ncol,nlev)
    real(jprb), intent(in), optional :: temperature_fl(ncol,nlev)

    ! Local optical depth
    real(jprb) :: od_local(ng)

    ! Natural logarithm of water content (kg m-3)
    real(jprb) :: ln_water_content
    
    ! Effective radius interpolation index and weights
    real(jprb) :: re_index, weight1, weight2

    ! Water content interpolation index and weights
    real(jprb) :: wc_index, weight_wc1, weight_wc2

    ! Temperature interpolation index and weights
    real(jprb) :: t_index, weight_t1, weight_t2

    ! Interpolation indices to first of the two points
    integer :: ire, iwc, it

    integer :: jcol, jlev, jg, jcomp

    real(jphook) :: hook_handle

    if (lhook) call dr_hook('radiation_general_cloud_optics_data:add_optical_properties',0,hook_handle)

    if (this%use_effective_radius) then

      ! Look-up versus effective radius
      if (present(scat_od)) then
        do jcol = 1,ncol
          do jlev = 1,nlev
            if (cloud_fraction(jcol, jlev) > 0.0_jprb) then
              re_index = max(1.0_jprb, min(1.0_jprb + (effective_radius(jcol,jlev)-this%effective_radius_0) &
                   &              / this%d_effective_radius, this%n_effective_radius-0.0001_jprb))
              ire = int(re_index)
              weight2 = re_index - ire
              weight1 = 1.0_jprb - weight2
              od_local = water_path(jcol, jlev) * (weight1*this%mass_ext(:,ire) &
                   &                              +weight2*this%mass_ext(:,ire+1))
              od(:,jlev,jcol) = od(:,jlev,jcol) + od_local
              od_local = od_local * (weight1*this%ssa(:,ire) &
                   &                +weight2*this%ssa(:,ire+1))
              scat_od(:,jlev,jcol) = scat_od(:,jlev,jcol) + od_local
              ! Loop over components describing the phase function
              ! (might be just asymmetry factor)
              do jcomp = 1,this%n_pf_components
                scat_pf(:,jlev,jcol,jcomp) = scat_pf(:,jlev,jcol,jcomp) &
                     & + od_local * (weight1*this%pf(:,ire,jcomp) &
                     &              +weight2*this%pf(:,ire+1,jcomp))
              end do
            end if
          end do
        end do
      else
        ! No scattering: return the absorption optical depth
        do jcol = 1,ncol
          do jlev = 1,nlev
            if (water_path(jcol, jlev) > 0.0_jprb) then
              re_index = max(1.0, min(1.0_jprb + (effective_radius(jcol,jlev)-this%effective_radius_0) &
                   &              / this%d_effective_radius, this%n_effective_radius-0.0001_jprb))
              ire = int(re_index)
              weight2 = re_index - ire
              weight1 = 1.0_jprb - weight2
              od(:,jlev,jcol) = od(:,jlev,jcol) &
                   &  + water_path(jcol, jlev) * (weight1*this%mass_ext(:,ire) &
                   &                             +weight2*this%mass_ext(:,ire+1)) &
                   &  * (1.0_jprb - (weight1*this%ssa(:,ire)+weight2*this%ssa(:,ire+1)))
            end if
          end do
        end do
      end if

    else

      ! Look-up versus water content and temperature
      if (present(scat_od)) then
        do jcol = 1,ncol
          do jlev = 1,nlev
            if (cloud_fraction(jcol, jlev) > 0.0_jprb) then
              ln_water_content = log(max(1.0e-12, water_path(jcol,jlev)/layer_depth(jcol,jlev)))
              t_index = max(1.0_jprb, min(1.0_jprb + (temperature_fl(jcol,jlev)-this%temperature_0) &
                   &              / this%d_temperature, this%n_temperature-0.0001_jprb))
              wc_index = max(1.0_jprb, min(1.0_jprb + (ln_water_content-this%ln_water_content_0) &
                   &              / this%d_ln_water_content, this%n_water_content-0.0001_jprb))
              it = int(t_index)
              weight_t2 = t_index - it
              weight_t1 = 1.0_jprb - weight_t2
              
              iwc = int(wc_index)
              weight_wc2 = wc_index - iwc
              weight_wc1 = 1.0_jprb - weight_wc2
              od_local = water_path(jcol, jlev) * (weight_t1*(weight_wc1*this%mass_ext_wct(:,iwc,it) &
                   &                                         +weight_wc2*this%mass_ext_wct(:,iwc+1,it)) &
                   &                              +weight_t2*(weight_wc1*this%mass_ext_wct(:,iwc,it+1) &
                   &                                         +weight_wc2*this%mass_ext_wct(:,iwc+1,it+1)))
              od(:,jlev,jcol) = od(:,jlev,jcol) + od_local
              od_local = od_local * (weight_t1*(weight_wc1*this%ssa_wct(:,iwc,it) &
                   &                           +weight_wc2*this%ssa_wct(:,iwc+1,it)) &
                   &                +weight_t2*(weight_wc1*this%ssa_wct(:,iwc,it+1) &
                   &                           +weight_wc2*this%ssa_wct(:,iwc+1,it+1)))
              scat_od(:,jlev,jcol) = scat_od(:,jlev,jcol) + od_local
              scat_pf(:,jlev,jcol,1) = scat_pf(:,jlev,jcol,1) &
                   & + od_local * (weight_t1*(weight_wc1*this%asymmetry_wct(:,iwc,it) &
                   &                         +weight_wc2*this%asymmetry_wct(:,iwc+1,it)) &
                   &              +weight_t2*(weight_wc1*this%asymmetry_wct(:,iwc,it+1) &
                   &                         +weight_wc2*this%asymmetry_wct(:,iwc+1,it+1)))
            end if
          end do
        end do
      else
        ! No scattering: return the absorption optical depth
        do jcol = 1,ncol
          do jlev = 1,nlev
            if (water_path(jcol, jlev) > 0.0_jprb) then
              re_index = max(1.0_jprb, min(1.0_jprb + (effective_radius(jcol,jlev)-this%effective_radius_0) &
                   &              / this%d_effective_radius, this%n_effective_radius-0.0001_jprb))
              ire = int(re_index)
              weight2 = re_index - ire
              weight1 = 1.0_jprb - weight2
              od(:,jlev,jcol) = od(:,jlev,jcol) &
                   &  + water_path(jcol, jlev) * (weight1*this%mass_ext(:,ire) &
                   &                             +weight2*this%mass_ext(:,ire+1)) &
                   &  * (1.0_jprb - (weight1*this%ssa(:,ire)+weight2*this%ssa(:,ire+1)))
            end if
          end do
        end do
      end if
        
    end if

    if (lhook) call dr_hook('radiation_general_cloud_optics_data:add_optical_properties',1,hook_handle)

  end subroutine add_optical_properties

#ifdef FLOTSAM
  !---------------------------------------------------------------------
  ! Add the optical properties of a particular cloud type to the
  ! accumulated optical properties of all cloud types, for FLOTSAM
  subroutine add_optical_properties_flotsam(this, ng, nlev, ncol, npf, &
       &                            cloud_fraction, &
       &                            water_path, effective_radius, &
       &                            scattering_angle, &
       &                            od, scat_od, pf)

    use yomhook, only : lhook, dr_hook, jphook

    class(general_cloud_optics_type), intent(in) :: this

    ! Number of g points, levels, columns and phase-function elements
    integer, intent(in) :: ng, nlev, ncol, npf

    ! Properties of present cloud type, dimensioned (ncol,nlev)
    real(jprb), intent(in) :: cloud_fraction(:,:)
    real(jprb), intent(in) :: water_path(:,:)       ! kg m-2
    real(jprb), intent(in) :: effective_radius(:,:) ! m

    ! Scattering angle (radians) between sun and sensor, dimensioned
    ! (ncol)
    real(jprb), intent(in) :: scattering_angle(ncol)

    ! Optical properties which are additive per cloud type,
    ! dimensioned (ng,nlev,ncol)
    real(jprb), intent(inout), dimension(ng,nlev,ncol) &
         &  :: od, &        ! Optical depth of layer
         &     scat_od      ! Scattering optical depth of layer
    ! Phase function elements
    real(jprb), intent(inout), dimension(ng,nlev,ncol,npf) :: pf

    real(jprb) :: od_local(ng)

    real(jprb) :: re_index, weight1, weight2
    integer :: ire

    integer :: iang
    real(jprb) :: w1, w2

    integer :: jcol, jlev, jcomp

    real(jphook) :: hook_handle

    if (lhook) call dr_hook('radiation_general_cloud_optics_data:add_optical_properties_flotsam',0,hook_handle)

    do jcol = 1,ncol
      ! Find indices and weights to bracketing angles
      iang = 1
      do while (this%scattering_angle(iang+1) < scattering_angle(jcol) &
           &    .and. iang < size(this%scattering_angle)-1)
        iang = iang+1
      end do
      w1 = (this%scattering_angle(iang+1)-scattering_angle(jcol)) &
           &  / (this%scattering_angle(iang+1)-this%scattering_angle(iang))
      w2 = (scattering_angle(jcol)-this%scattering_angle(iang)) &
           &  / (this%scattering_angle(iang+1)-this%scattering_angle(iang))

      do jlev = 1,nlev
        if (cloud_fraction(jcol, jlev) > 0.0_jprb) then
          re_index = max(1.0_jprb, min(1.0_jprb + (effective_radius(jcol,jlev)-this%effective_radius_0) &
               &              / this%d_effective_radius, this%n_effective_radius-0.0001_jprb))
          ire = int(re_index)
          weight2 = re_index - ire
          weight1 = 1.0_jprb - weight2
          ! Local extinction optical depth
          od_local = water_path(jcol, jlev) * (weight1*this%mass_ext(:,ire) &
               &                              +weight2*this%mass_ext(:,ire+1))
          od(:,jlev,jcol) = od(:,jlev,jcol) + od_local
          ! Convert "od_local" to be the scattering optical depth
          od_local = od_local * (weight1*this%ssa(:,ire) &
               &                +weight2*this%ssa(:,ire+1))
          scat_od(:,jlev,jcol) = scat_od(:,jlev,jcol) + od_local
          pf(:,jlev,jcol,1) = pf(:,jlev,jcol,1) + od_local &
               &  * (weight1 * (w1*this%phase_function(iang,:,ire) &
               &               +w2*this%phase_function(iang+1,:,ire)) &
               &    +weight2 * (w1*this%phase_function(iang,:,ire+1) &
               &               +w2*this%phase_function(iang+1,:,ire+1)))
          do jcomp = 1,size(this%phase_function_components,1)
            pf(:,jlev,jcol,jcomp+1) = pf(:,jlev,jcol,jcomp+1) + od_local &
                 &  * (weight1 * this%phase_function_components(jcomp,:,ire) &
                 &    +weight2 * this%phase_function_components(jcomp,:,ire+1))
          end do
        end if
      end do
    end do

    if (lhook) call dr_hook('radiation_general_cloud_optics_data:add_optical_properties_flotsam',1,hook_handle)

  end subroutine add_optical_properties_flotsam
#endif

  !---------------------------------------------------------------------
  ! Return the Planck function (in W m-2 (cm-1)-1) for a given
  ! wavenumber (cm-1) and temperature (K), ensuring double precision
  ! for internal calculation
  elemental function calc_planck_function_wavenumber(wavenumber, temperature)

    use parkind1,            only : jprb, jprd
    use radiation_constants, only : SpeedOfLight, BoltzmannConstant, PlanckConstant

    real(jprb), intent(in) :: wavenumber  ! cm-1
    real(jprb), intent(in) :: temperature ! K
    real(jprb) :: calc_planck_function_wavenumber

    real(jprd) :: freq ! Hz
    real(jprd) :: planck_fn_freq ! W m-2 Hz-1

    freq = 100.0_jprd * real(SpeedOfLight,jprd) * real(wavenumber,jprd)
    planck_fn_freq = 2.0_jprd * real(PlanckConstant,jprd) * freq**3 &
         &  / (real(SpeedOfLight,jprd)**2 * (exp(real(PlanckConstant,jprd)*freq &
         &     /(real(BoltzmannConstant,jprd)*real(temperature,jprd))) - 1.0_jprd))
    calc_planck_function_wavenumber = real(planck_fn_freq * 100.0_jprd * real(SpeedOfLight,jprd), jprb)

  end function calc_planck_function_wavenumber

  !---------------------------------------------------------------------
  ! Perform a Legendre decomposition of the scattering phase function
  subroutine legendre_decomposition(nang, scattering_angle, pf_in, n_pf_components, pf)
    integer,    intent(in)  :: nang
    real(jprb), intent(in)  :: scattering_angle(:) ! radians
    real(jprb), intent(in)  :: pf_in(:,:,:)
    integer,    intent(in)  :: n_pf_components
    real(jprb), intent(out) :: pf(:,:,:)

    real(jprb) :: cos_ang(nang)
    real(jprb) :: weight(nang-1)
    real(jprb) :: p0_legendre(nang)
    real(jprb) :: p1_legendre(nang)
    real(jprb) :: p_legendre(nang)
    real(jprb) :: normalization_factor
    
    integer :: jc, ja

    ! There are two ways of expressing a phase function with Legendre
    ! polynomials Pn(mu): the unnormalized way:
    !   p(mu) = Sum[ c_n*Pn(mu) ]
    ! or the DISORT way:
    !   p(mu) = Sum[ g_n*(2*n+1)*Pn(mu) ]
    ! In the DISORT way the weights g_n lie between -1 and 1, and g_1
    ! is equal to the asymmetry factor. In the unnormalized way the
    ! decomposition is a little simpler but coefficient c_1 is three
    ! times the asymmety factor.
    logical, parameter :: use_disort_normalization = .true.
    
    cos_ang = cos(scattering_angle)
    weight  = cos_ang(1:nang-1)-cos_ang(2:nang)
    
    p0_legendre = 1.0_jprb
    p1_legendre = cos_ang;

    pf = 0.0_jprb
    
    do jc = 1,n_pf_components
      if (jc == 1) then
        p_legendre = cos_ang
      else
        p_legendre = ((2.0_jprb*jc + 1.0_jprb)*cos_ang*p1_legendre - jc*p0_legendre) / (jc+1.0_jprb)
        p0_legendre = p1_legendre
        p1_legendre = p_legendre
      end if
      if (use_disort_normalization) then
        normalization_factor = 1.0_jprb
      else
        normalization_factor = 2.0_jprb*jc + 1.0_jprb
      end if
      do ja = 1,nang-1
        ! Note that the factor of 0.25 here is because we are
        ! averaging over two points (a factor of a half) and we are
        ! integrating over mu in the range -1 to 1 (another factor of
        ! a half required)
        pf(:,:,jc) = pf(:,:,jc) + weight(ja)*normalization_factor*0.25_jprb &
             &  * (p_legendre(ja)*pf_in(ja,:,:) + p_legendre(ja+1)*pf_in(ja+1,:,:))
      end do
    end do
    
  end subroutine legendre_decomposition

  !---------------------------------------------------------------------
  ! Save cloud optical properties in the named file
  subroutine save_general_cloud_optics_data(this, file_name, iverbose)

    use yomhook,     only : lhook, dr_hook, jphook
    use easy_netcdf, only : netcdf_file

    class(general_cloud_optics_type), intent(in) :: this
    character(len=*),                 intent(in) :: file_name
    integer,                optional, intent(in) :: iverbose

    ! Object for output NetCDF file
    type(netcdf_file) :: out_file

    real(jprb) :: effective_radius(this%n_effective_radius)
    integer :: ire
    
    real(jphook) :: hook_handle

    if (lhook) call dr_hook('radiation_general_cloud_optics_data:save',0,hook_handle)

    ! Create the file
    call out_file%create(trim(file_name), iverbose=iverbose)

    ! Define dimensions
    call out_file%define_dimension("band", size(this%mass_ext,1))
    call out_file%define_dimension("effective_radius", this%n_effective_radius)

    ! Put global attributes
    call out_file%put_global_attributes( &
         &   title_str="Optical properties of "//trim(this%type_name) &
         &   //" hydrometeors using the spectral intervals of ecRad", &
         &   source_str="ecRad offline radiation model")

    ! Define variables
    call out_file%define_variable("effective_radius", units_str="m", &
         &  long_name="Effective radius", dim1_name="effective_radius")
    call out_file%define_variable("mass_extinction_coefficient", units_str="m2 kg-1", &
         &  long_name="Mass-extinction coefficient", &
         &  dim2_name="effective_radius", dim1_name="band")
    call out_file%define_variable("single_scattering_albedo", units_str="1", &
         &  long_name="Single scattering albedo", &
         &  dim2_name="effective_radius", dim1_name="band")
    if (this%pf_mode == IPhaseFuncAsymmetry) then
      call out_file%define_variable("asymmetry_factor", units_str="1", &
         &  long_name="Asymmetry factor", &
         &  dim2_name="effective_radius", dim1_name="band")
    end if

    ! Define effective radius
    do ire = 1,this%n_effective_radius
      effective_radius(ire) = this%effective_radius_0 + this%d_effective_radius*(ire-1)
    end do
    
    ! Write variables
    call out_file%put("effective_radius", effective_radius)
    call out_file%put("mass_extinction_coefficient", this%mass_ext)
    call out_file%put("single_scattering_albedo", this%ssa)
    if (this%pf_mode == IPhaseFuncAsymmetry) then
      call out_file%put("asymmetry_factor", this%pf(:,:,1))
    end if
    
    call out_file%close()

    if (lhook) call dr_hook('radiation_general_cloud_optics_data:save',1,hook_handle)

  end subroutine save_general_cloud_optics_data
  
end module radiation_general_cloud_optics_data
