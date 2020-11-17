! radiation_spectral_definition.F90 - Derived type to describe a spectral definition
!
! (C) Copyright 2020- ECMWF.
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

module radiation_spectral_definition

  use parkind1,    only : jprb

  implicit none

  public

  real(jprb), parameter :: SolarReferenceTemperature = 5777.0_jprb ! K
  real(jprb), parameter :: TerrestrialReferenceTemperature = 273.15_jprb ! K

  !---------------------------------------------------------------------
  ! A derived type describing the contribution of the g points of a
  ! correlated k-distribution gas-optics model from each part of the
  ! spectrum. This is used primarily to map the cloud and aerosol
  ! optical properties on to the gas g points.
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
    procedure :: find => find_wavenumber
    procedure :: calc_mapping

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


  !---------------------------------------------------------------------
  ! Find the index to the highest wavenumber in the spectral
  ! definition that is lower than or equal to "wavenumber", used for
  ! implementing look-up tables
  pure function find_wavenumber(this, wavenumber)
    class(spectral_definition_type), intent(in) :: this
    real(jprb),                      intent(in) :: wavenumber ! cm-1
    integer                                     :: find_wavenumber

    integer :: iwav

    if (wavenumber < this%wavenumber1(1) .or. wavenumber > this%wavenumber2(this%nwav)) then
      ! Wavenumber not present
      find_wavenumber = 0
    else
      find_wavenumber = 1
      do while (wavenumber > this%wavenumber2(find_wavenumber) &
           &    .and. find_wavenumber < this%nwav)
        find_wavenumber = find_wavenumber + 1
      end do
    end if
  end function find_wavenumber


  !---------------------------------------------------------------------
  ! Compute a mapping matrix "mapping" that can be used in an
  ! expression y=matmul(mapping,x) where x is a variable containing
  ! optical properties at each input "wavenumber", and y is this
  ! variable mapped on to the spectral intervals in the spectral
  ! definition "this". Temperature (K) is used to generate a Planck
  ! function to weight each wavenumber appropriately.
  subroutine calc_mapping(this, temperature, wavenumber, mapping, use_bands)

    use yomhook,     only : lhook, dr_hook

    class(spectral_definition_type), intent(in)    :: this
    real(jprb),                      intent(in)    :: temperature   ! K
    real(jprb),                      intent(in)    :: wavenumber(:) ! cm-1
    real(jprb), allocatable,         intent(inout) :: mapping(:,:)
    logical,    optional,            intent(in)    :: use_bands

    ! Spectral weights to apply, same length as wavenumber above
    real(jprb), dimension(:), allocatable :: weight, planck_weight

    ! Wavenumbers (cm-1) marking triangle of influence of a cloud
    ! spectral point
    real(jprb) :: wavenum0, wavenum1, wavenum2

    integer    :: nwav ! Number of wavenumbers describing cloud

    ! Indices to wavenumber intervals in spectral definition structure
    integer    :: isd, isd0, isd1, isd2

    ! Loop indices
    integer    :: jg, jwav

    logical    :: use_bands_local

    real(jprb) :: hook_handle

    if (lhook) call dr_hook('radiation_spectral_definition:calc_mapping',0,hook_handle)

    if (present(use_bands)) then
      use_bands_local = use_bands
    else
      use_bands_local = .false.
    end if

    nwav = size(wavenumber)

    if (allocated(mapping)) then
      deallocate(mapping)
    end if
    
    ! Define the mapping matrix
    if (use_bands_local) then
      ! Cloud properties per band
      allocate(mapping(this%nband, nwav))
    else
      allocate(mapping(this%ng, nwav))
      allocate(weight(this%nwav))
      allocate(planck_weight(this%nwav))

      planck_weight = calc_planck_function_wavenumber(0.5_jprb &
           &  * (this%wavenumber1 + this%wavenumber2), &
           &  temperature)

      mapping = 0.0_jprb
      ! Loop over wavenumbers representing cloud
      do jwav = 1,nwav
        ! Clear the weights. The weight says for one wavenumber in the
        ! cloud file, what is its fractional contribution to each of
        ! the spectral-definition intervals
        weight = 0.0_jprb

        ! Cloud properties are linearly interpolated between each of
        ! the nwav cloud points; therefore, the influence of a
        ! particular cloud point extends as a triangle between
        ! wavenum0 and wavenum2, peaking at wavenum1
        wavenum1 = wavenumber(jwav)
        isd1 = this%find(wavenum1)
        if (isd1 < 1) then
          cycle
        end if
        if (jwav > 1) then
          wavenum0 = wavenumber(jwav-1)

          ! Map triangle under (wavenum0,0) to (wavenum1,1) to the
          ! wavenumbers in this%gpoint_fraction
          isd0 = this%find(wavenum0)
          if (isd0 == isd1) then
            ! Triangle completely within the range
            ! this%wavenumber1(isd0)-this%wavenumber2(isd0)
            weight(isd0) = 0.5_jprb*(wavenum1-wavenum0) &
                 &       / (this%wavenumber2(isd0)-this%wavenumber1(isd0))
          else
            if (isd0 >= 1) then
              ! Left part of triangle
              weight(isd0) = 0.5_jprb * (this%wavenumber2(isd0)-wavenum0)**2 &
                   &       / ((this%wavenumber2(isd0)-this%wavenumber1(isd0)) &
                   &         *(wavenum1-wavenum0))
            end if
            ! Right part of triangle (trapezium)
!            weight(isd1) = 0.5_jprb * (wavenum1-this%wavenumber1(isd1)) &
!                 &       * (wavenum1 + this%wavenumber1(isd1) - 2.0_jprb*wavenum0) &
!                 &       / (wavenum1-wavenum0)
            weight(isd1) = 0.5_jprb * (1.0_jprb &
                 &  + (this%wavenumber1(isd1)-wavenum1)/(wavenum1-wavenum0)) &
                 &  * (wavenum1-this%wavenumber1(isd1)) &
                 &  / (this%wavenumber2(isd1)-this%wavenumber1(isd1))
            if (isd1-isd0 > 1) then
              do isd = isd0+1,isd1-1
                ! Intermediate trapezia
                weight(isd) = 0.5_jprb * (this%wavenumber1(isd)+this%wavenumber2(isd) &
                     &                    - 2.0_jprb*wavenum0) &
                     &      / (wavenum1-wavenum0)
              end do
            end if
          end if

        else
          ! First cloud wavenumber: all wavenumbers in the spectral
          ! definition below this will use the first one
          if (isd1 >= 1) then
            weight(1:isd1-1) = 1.0_jprb
            weight(isd1) = (wavenum1-this%wavenumber1(isd1)) &
                 &       / (this%wavenumber2(isd1)-this%wavenumber1(isd1))
          end if
        end if

        if (jwav < nwav) then
          wavenum2 = wavenumber(jwav+1)

          ! Map triangle under (wavenum1,1) to (wavenum2,0) to the
          ! wavenumbers in this%gpoint_fraction
          isd2 = this%find(wavenum2)

          if (isd1 == isd2) then
            ! Triangle completely within the range
            ! this%wavenumber1(isd1)-this%wavenumber2(isd1)
            weight(isd1) = weight(isd1) + 0.5_jprb*(wavenum2-wavenum1) &
                 &       / (this%wavenumber2(isd1)-this%wavenumber1(isd1))
          else
            if (isd2 >= 1 .and. isd2 <= this%nwav) then
              ! Right part of triangle
              weight(isd2) = weight(isd2) + 0.5_jprb * (wavenum2-this%wavenumber1(isd2))**2 &
                   &       / ((this%wavenumber2(isd2)-this%wavenumber1(isd2)) &
                   &         *(wavenum2-wavenum1))
            end if
            ! Left part of triangle (trapezium)
!            weight(isd1) = weight(isd1) + 0.5_jprb * (this%wavenumber2(isd1)-wavenum1) &
!                 &       * (wavenum1 + this%wavenumber2(isd1) - 2.0_jprb*wavenum2) &
!                 &       / (wavenum2-wavenum1)
            weight(isd1) = weight(isd1) + 0.5_jprb * (1.0_jprb &
                 &  + (wavenum2-this%wavenumber2(isd1)) / (wavenum2-wavenum1)) &
                 &  * (this%wavenumber2(isd1)-wavenum1) &
                 &  / (this%wavenumber2(isd1)-this%wavenumber1(isd1))
            if (isd2-isd1 > 1) then
              do isd = isd1+1,isd2-1
                ! Intermediate trapezia
                weight(isd) = weight(isd) + 0.5_jprb * (2.0_jprb*wavenum2 &
                     & - this%wavenumber1(isd) - this%wavenumber2(isd)) &
                     &      / (wavenum2-wavenum1)
              end do
            end if
          end if

        else
          ! Last cloud wavenumber: all wavenumbers in the spectral
          ! definition above this will use the last one
          if (isd1 <= this%nwav) then
            weight(isd1+1:this%nwav) = 1.0_jprb
            weight(isd1) = (this%wavenumber2(isd1)-wavenum1) &
                 &       / (this%wavenumber2(isd1)-this%wavenumber1(isd1))
          end if
        end if

        weight = weight * planck_weight

        do jg = 1,this%ng
          mapping(jg, jwav) = sum(weight * this%gpoint_fraction(:,jg))
        end do

      end do

      deallocate(weight)
      deallocate(planck_weight)

    end if

    ! Normalize mapping matrix
    do jg = 1,this%ng
      mapping(jg,:) = mapping(jg,:) * (1.0_jprb/sum(mapping(jg,:)))
    end do

    if (lhook) call dr_hook('radiation_spectral_definition:calc_mapping',1,hook_handle)

  end subroutine calc_mapping


  !---------------------------------------------------------------------
  ! Return the Planck function (in W m-2 (cm-1)-1) for a given
  ! wavenumber (cm-1) and temperature (K)
  elemental function calc_planck_function_wavenumber(wavenumber, temperature)

    use radiation_constants, only : SpeedOfLight, BoltzmannConstant, PlanckConstant

    real(jprb), intent(in) :: wavenumber  ! cm-1
    real(jprb), intent(in) :: temperature ! K
    real(jprb) :: calc_planck_function_wavenumber

    real(jprb) :: freq
    freq = wavenumber * 100.0_jprb * SpeedOfLight
    calc_planck_function_wavenumber = (2.0_jprb * PlanckConstant * freq**3 * SpeedOfLight**-2) &
         &  / (exp(PlanckConstant*freq/(BoltzmannConstant*temperature)) - 1.0_jprb)
    calc_planck_function_wavenumber = calc_planck_function_wavenumber * 100.0_jprb * SpeedOfLight

  end function calc_planck_function_wavenumber

end module radiation_spectral_definition
