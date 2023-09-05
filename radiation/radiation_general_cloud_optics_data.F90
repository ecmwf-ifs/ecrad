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

    ! Phase function description, dimensioned (nband,n_effective_radius,ncomponent)
    real(jprb), allocatable, dimension(:,:,:) :: pf

    ! Number of effective radius coefficients, start value and
    ! interval in look-up table
    integer    :: n_effective_radius = 0
    real(jprb) :: effective_radius_0, d_effective_radius

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

   contains
     procedure :: setup => setup_general_cloud_optics
     procedure :: add_optical_properties
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
         weighting_temperature=weighting_temperature, use_bands=use_bands)

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
      !print *, 'Asymmetry factor from Legendre: ', this%pf(:,nre,1)/3.0_jprb
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

    ! Revert back to unscaled quantities
    call revert_delta_eddington(this%mass_ext, this%ssa, asymmetry_band)

    if (this%pf_mode == IPhaseFuncAsymmetry) then
      this%pf(:,:,1) = asymmetry_band
    !else
    !  print *, 'Asymmetry factor:              ', asymmetry_band(:,nre)
    end if
    
    deallocate(asymmetry_band)  
    deallocate(mapping)
    
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
      if (this%pf_mode == IPhaseFuncAsymmetry) then
        write(nulout,'(a)') '  Phase function represented by asymmetry factor'
      else if (this%pf_mode == IPhaseFuncLegendre) then
        write(nulout,'(a,i0,a)') '  Phase function represented by ', this%n_pf_components, ' Legendre terms'
      end if
    end if

    if (lhook) call dr_hook('radiation_general_cloud_optics_data:setup',1,hook_handle)

  end subroutine setup_general_cloud_optics


  !---------------------------------------------------------------------
  ! Add the optical properties of a particular cloud type to the
  ! accumulated optical properties of all cloud types
  subroutine add_optical_properties(this, ng, nlev, ncol, &
       &                            cloud_fraction, &
       &                            water_path, effective_radius, &
       &                            od, scat_od, scat_pf)

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

    real(jprb) :: od_local(ng)

    real(jprb) :: re_index, weight1, weight2
    integer :: ire

    integer :: jcol, jlev, jcomp

    real(jphook) :: hook_handle

    if (lhook) call dr_hook('radiation_general_cloud_optics_data:add_optical_properties',0,hook_handle)

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

    if (lhook) call dr_hook('radiation_general_cloud_optics_data:add_optical_properties',1,hook_handle)

  end subroutine add_optical_properties


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
