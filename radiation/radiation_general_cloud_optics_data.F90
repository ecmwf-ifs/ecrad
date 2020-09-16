! radiation_general_cloud_optics_data.F90 - Type to store generalized cloud optical properties
!
! Copyright (C) 2019 ECMWF
!
! Author:  Robin Hogan
! Email:   r.j.hogan@ecmwf.int
! License: see the COPYING file for details
!

module radiation_general_cloud_optics_data

  use parkind1, only : jprb

  implicit none

  !---------------------------------------------------------------------
  ! This type holds the configuration information to compute optical
  ! properties for a particular type of cloud or hydrometeor in one of
  ! the shortwave or longwave
  type general_cloud_optics_type
    ! Band-specific (or g-point-specific) values as a look-up table
    ! versus effective radius dimensioned (nband,n_effective_radius)
    
    ! Extinction coefficient per unit mass (m2 kg-1)
    real(jprb), allocatable, dimension(:,:) :: &
         &  mass_ext
    
    ! Single-scattering albedo and asymmetry factor (dimensionless)
    real(jprb), allocatable, dimension(:,:) :: &
         &  ssa, asymmetry

    ! Number of effective radius coefficients, start value and
    ! interval in look-up table
    integer    :: n_effective_radius = 0
    real(jprb) :: effective_radius_0, d_effective_radius

    ! Name of cloud/precip type (e.g. "liquid", "ice", "rain", "snow")
    ! and the name of the optics scheme.  These two are used to
    ! generate the name of the data file from which the coefficients
    ! are read.
    character(len=511) :: type_name, scheme_name
    
    ! Do we use bands or g-points?
    logical :: use_bands = .false.

   contains
     procedure :: setup => setup_general_cloud_optics

  end type general_cloud_optics_type

contains

  !---------------------------------------------------------------------
  ! Setup cloud optics coefficients by reading them from a file
  subroutine setup_general_cloud_optics(this, file_name, specdef, &
       &                                use_bands, use_thick_averaging, &
       &                                weighting_temperature, &
       &                                iverbose)

    use yomhook,                       only : lhook, dr_hook
    use easy_netcdf,                   only : netcdf_file
    use radiation_spectral_definition, only : spectral_definition_type
    use radiation_io,                  only : nulerr, radiation_abort

    class(general_cloud_optics_type), intent(inout)    :: this
    character(len=*), intent(in)               :: file_name
    type(spectral_definition_type), intent(in) :: specdef
    logical, intent(in), optional              :: use_bands, use_thick_averaging
    real(jprb), intent(in), optional           :: weighting_temperature ! K
    integer, intent(in), optional              :: iverbose
    
    ! Spectral properties read from file, dimensioned (wavenumber,
    ! n_effective_radius)
    real(jprb), dimension(:,:), allocatable :: mass_ext, & ! m2 kg-1
         &                                     ssa, asymmetry

    ! Reflectance of an infinitely thick cloud, needed for thick
    ! averaging
    real(jprb), dimension(:,:), allocatable :: ref_inf

    ! Coordinate variables from file
    real(jprb), dimension(:), allocatable :: wavenumber       ! cm-1
    real(jprb), dimension(:), allocatable :: effective_radius ! m

    ! Spectral weights to apply, same length as wavenumber above
    real(jprb), dimension(:), allocatable :: weight

    ! Matrix mapping optical properties in the file to values per
    ! g-point or band, such that in the thin-averaging case,
    ! this%mass_ext=matmul(mapping,file%mass_ext), so mapping is
    ! dimensioned (ngpoint,nwav)
    real(jprb), dimension(:,:), allocatable :: mapping

    ! The NetCDF file containing the coefficients
    type(netcdf_file)  :: file

    real(jprb) :: diff_spread
    integer    :: iverb
    integer    :: nre  ! Number of effective radii
    integer    :: nwav ! Number of wavenumbers describing cloud

    ! Wavenumbers (cm-1) marking triangle of influence of a cloud
    ! spectral point
    real(jprb) :: wavenum0, wavenum1, wavenum2

    ! Indices to wavenumber intervals in spectral definition structure
    integer    :: isd0, isd1, isd2, isd

    ! Loop indices
    integer    :: jg, jwav

    logical    :: use_bands_local, use_thick_averaging_local

    real(jprb) :: hook_handle

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

    ! Set up weighting
    if (.not. present(weighting_temperature)) then
      write(nulerr, '(a)') '*** Error: weighting_temperature not provided'
      call radiation_abort('Radiation configuration error')
    end if

    nwav = size(wavenumber)

    ! Define the mapping matrix
    if (use_bands_local) then
      ! Cloud properties per band
      allocate(mapping(specdef%nband, nwav))
      allocate(this%mass_ext(specdef%nband, nre))
      allocate(this%ssa(specdef%nband, nre))
      allocate(this%asymmetry(specdef%nband, nre))

    else
      ! Cloud properties separately per g point
      allocate(mapping(specdef%ng, nwav))
      allocate(this%mass_ext(specdef%ng, nre))
      allocate(this%ssa(specdef%ng, nre))
      allocate(this%asymmetry(specdef%ng, nre))

      allocate(weight(specdef%nwav))

      mapping = 0.0_jprb
      ! Loop over wavenumbers representing cloud
      do jwav = 1,nwav
        ! Clear the weights
        weight = 0.0_jprb

        ! Cloud properties are linearly interpolated between each of
        ! the nwav cloud points; therefore, the influence of a
        ! particular cloud point extends as a triangle between
        ! wavenum0 and wavenum2, peaking at wavenum1
        wavenum1 = wavenumber(jwav)
        isd1 = specdef%find(wavenum1)
        if (jwav > 1) then
          wavenum0 = wavenumber(jwav-1)

          ! Map triangle under (wavenum0,0) to (wavenum1,1) to the
          ! wavenumbers in specdef%gpoint_fraction
          isd0 = specdef%find(wavenum0)
          if (isd0 == isd1) then
            ! Triangle completely within the range
            ! specdef%wavenumber1(isd0)-specdef%wavenumber2(isd0)
            weight(isd0) = 0.5_jprb*(wavenum1-wavenum0) &
                 &       / (specdef%wavenumber2(isd0)-specdef%wavenumber1(isd0))
          else
            ! Left part of triangle
            weight(isd0) = 0.5_jprb * (specdef%wavenumber2(isd0)-wavenum0)**2 &
                 &       / ((specdef%wavenumber2(isd0)-specdef%wavenumber1(isd0)) &
                 &         *(wavenum1-wavenum0))
            ! Right part of triangle (trapezium)
            weight(isd1) = 0.5_jprb * (wavenum1-specdef%wavenumber1(isd1)) &
                 &       * (wavenum1 + specdef%wavenumber1(isd1) - 2.0_jprb*wavenum0) &
                 &       / (wavenum1-wavenum0)
            do isd = isd0+1,isd1-1
              ! Intermediate trapezia
              weight(isd) = 0.5_jprb * (specdef%wavenumber1(isd)+specdef%wavenumber2(isd) &
                   &                    - 2.0_jprb*wavenum0) &
                   &      / (wavenum1-wavenum0)
            end do
          end if

        else
          ! First cloud wavenumber: all wavenumbers in the spectral
          ! definition below this will use the first one
          weight(1:isd1-1) = 1.0_jprb
          weight(isd1) = (wavenum1-specdef%wavenumber1(isd1)) &
               &       / (specdef%wavenumber2(isd1)-specdef%wavenumber1(isd1))
        end if

        if (jwav < nwav) then
          wavenum2 = wavenumber(jwav+1)

          ! Map triangle under (wavenum1,1) to (wavenum2,0) to the
          ! wavenumbers in specdef%gpoint_fraction
          isd2 = specdef%find(wavenum2)

          if (isd1 == isd2) then
            ! Triangle completely within the range
            ! specdef%wavenumber1(isd1)-specdef%wavenumber2(isd1)
            weight(isd1) = weight(isd1) + 0.5_jprb*(wavenum2-wavenum1) &
                 &       / (specdef%wavenumber2(isd1)-specdef%wavenumber1(isd1))
          else
            ! Right part of triangle
            weight(isd2) = weight(isd2) + 0.5_jprb * (wavenum2-specdef%wavenumber1(isd2))**2 &
                 &       / ((specdef%wavenumber2(isd2)-specdef%wavenumber1(isd2)) &
                 &         *(wavenum2-wavenum1))
            ! Right part of triangle (trapezium)
            weight(isd1) = weight(isd1) + 0.5_jprb * (specdef%wavenumber2(isd1)-wavenum1) &
                 &       * (wavenum1 + specdef%wavenumber2(isd1) - 2.0_jprb*wavenum2) &
                 &       / (wavenum2-wavenum1)
            do isd = isd1+1,isd2-1
              ! Intermediate trapezia
              weight(isd) = weight(isd) + 0.5_jprb * (specdef%wavenumber1(isd)+specdef%wavenumber2(isd) &
                   &                    - 2.0_jprb*wavenum2) &
                   &      / (wavenum2-wavenum1)
            end do
          end if

        else
          ! Last cloud wavenumber: all wavenumbers in the spectral
          ! definition above this will use the last one
          weight(isd1+1) = 1.0_jprb
          weight(isd1) = (specdef%wavenumber2(isd1)-wavenum1) &
               &       / (specdef%wavenumber2(isd1)-specdef%wavenumber1(isd1))
        end if

        do jg = 1,specdef%ng
          mapping(jg, jwav) = sum(weight * specdef%gpoint_fraction(:,jg))
        end do

      end do

      deallocate(weight)

    end if

    ! Thin averaging
    this%mass_ext  = matmul(mapping, mass_ext)
    this%ssa       = matmul(mapping, mass_ext*ssa) / this%mass_ext
    this%asymmetry = matmul(mapping, mass_ext*ssa*asymmetry) / (this%mass_ext*this%ssa)
    
    if (use_thick_averaging_local) then
      ! Thick averaging as described by Edwards and Slingo (1996),
      ! modifying only the single-scattering albedo
      if (use_bands_local) then
        allocate(ref_inf(specdef%nband, nre))
      else
        allocate(ref_inf(specdef%ng, nre))
      end if

      ! Eqs. 18 and 17 of Edwards & Slingo (1996)
      ref_inf = sqrt((1.0_jprb - ssa) / (1.0_jprb - ssa*asymmetry))
      ref_inf = (1.0_jprb - ref_inf) / (1.0_jprb + ref_inf)
      ! Here the left-hand side is actually the averaged ref_inf
      this%ssa = matmul(mapping, ref_inf)
      ! Eq. 19 of Edwards and Slingo (1996)
      this%ssa = 4.0_jprb * this%ssa / ((1.0_jprb + this%ssa)**2 &
           &  - this%asymmetry * (1.0_jprb - this%ssa)**2)

      deallocate(ref_inf)
    end if

    deallocate(mapping)

    if (lhook) call dr_hook('radiation_general_cloud_optics_data:setup',1,hook_handle)

  end subroutine setup_general_cloud_optics

end module radiation_general_cloud_optics_data
