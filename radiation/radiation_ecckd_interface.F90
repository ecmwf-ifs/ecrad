! radiation_ecckd_interface.F90 - Interface to ecCKD gas optics model
!
! Copyright (C) 2020 ECMWF
!
! Author:  Robin Hogan
! Email:   r.j.hogan@ecmwf.int
! License: see the COPYING file for details
!

module radiation_ecckd_interface

  public  :: setup_gas_optics, gas_optics, planck_function

contains

  !---------------------------------------------------------------------
  ! Setup the ecCKD generalized gas optics model
  subroutine setup_gas_optics(config)

    use parkind1, only : jprb
    use radiation_config

    implicit none

    type(config_type), intent(inout) :: config

    integer :: jj
    
    if (config%do_sw) then

      ! Read shortwave ecCKD gas optics NetCDF file
      call config%gas_optics_sw%read(config%gas_optics_sw_file_name, &
           &                         config%iverbosesetup)

      ! Copy over relevant properties
      config%n_g_sw     = config%gas_optics_sw%ng

      if (config%do_cloud_aerosol_per_sw_g_point) then
        ! Bands and g points are the same
        config%n_bands_sw = config%n_g_sw
      else
        ! Bands are groups of g points and span a continuous region of
        ! wavenumber space
        config%n_bands_sw = config%gas_optics_sw%nband
        allocate(config%wavenumber1_sw(config%n_bands_sw))
        allocate(config%wavenumber2_sw(config%n_bands_sw))
      end if

      allocate(config%i_band_from_g_sw          (config%n_g_sw))
      allocate(config%i_band_from_reordered_g_sw(config%n_g_sw))
      allocate(config%i_g_from_reordered_g_sw   (config%n_g_sw))
        
      if (config%do_cloud_aerosol_per_sw_g_point) then
        config%i_band_from_g_sw           = [ (jj, jj = 1,config%n_g_sw) ]
        config%i_band_from_reordered_g_sw = [ (jj, jj = 1,config%n_g_sw) ]
      else
        config%i_band_from_g_sw           = config%gas_optics_sw%i_band_number
        config%i_band_from_reordered_g_sw = config%gas_optics_sw%i_band_number
      end if
      config%i_g_from_reordered_g_sw      = [ (jj, jj = 1,config%n_g_sw) ]

    end if

    if (config%do_lw) then

      ! Read longwave ecCKD gas optics NetCDF file
      call config%gas_optics_lw%read(config%gas_optics_lw_file_name, &
           &                         config%iverbosesetup)

      ! Copy over relevant properties
      config%n_g_lw     = config%gas_optics_lw%ng

      if (config%do_cloud_aerosol_per_lw_g_point) then
        ! Bands and g points are the same
        config%n_bands_lw = config%n_g_lw
      else
        ! Bands are groups of g points and span a continuous region of
        ! wavenumber space
        config%n_bands_lw = config%gas_optics_lw%nband
        allocate(config%wavenumber1_lw(config%n_bands_lw))
        allocate(config%wavenumber2_lw(config%n_bands_lw))
      end if

      allocate(config%i_band_from_g_lw          (config%n_g_lw))
      allocate(config%i_band_from_reordered_g_lw(config%n_g_lw))
      allocate(config%i_g_from_reordered_g_lw   (config%n_g_lw))
        
      if (config%do_cloud_aerosol_per_lw_g_point) then
        config%i_band_from_g_lw           = [ (jj, jj = 1,config%n_g_lw) ]
        config%i_band_from_reordered_g_lw = [ (jj, jj = 1,config%n_g_lw) ]
      else
        config%i_band_from_g_lw           = config%gas_optics_lw%i_band_number
        config%i_band_from_reordered_g_lw = config%gas_optics_lw%i_band_number
      end if
      config%i_g_from_reordered_g_lw      = [ (jj, jj = 1,config%n_g_lw) ]

    end if

end subroutine setup_gas_optics


  !---------------------------------------------------------------------
  ! Compute gas optical depths, shortwave scattering, Planck function
  ! and incoming shortwave radiation at top-of-atmosphere
  subroutine gas_optics(ncol,nlev,istartcol,iendcol, &
       &  config, single_level, thermodynamics, gas, & 
       &  od_lw, od_sw, ssa_sw, lw_albedo, planck_hl, lw_emission, &
       &  incoming_sw)

    use parkind1, only : jprb
    use radiation_config,         only : config_type
    use radiation_thermodynamics, only : thermodynamics_type
    use radiation_single_level,   only : single_level_type
    use radiation_gas

    implicit none

    integer,                  intent(in) :: ncol               ! number of columns
    integer,                  intent(in) :: nlev               ! number of levels
    integer,                  intent(in) :: istartcol, iendcol ! range of cols to process
    type(config_type),        intent(in) :: config
    type(single_level_type),  intent(in) :: single_level
    type(thermodynamics_type),intent(in) :: thermodynamics
    type(gas_type),           intent(in) :: gas

    ! Longwave albedo of the surface
    real(jprb), dimension(config%n_g_lw,istartcol:iendcol), &
         &  intent(in), optional :: lw_albedo

    ! Gaseous layer optical depth in longwave and shortwave, and
    ! shortwave single scattering albedo (i.e. fraction of extinction
    ! due to Rayleigh scattering) at each g-point
    real(jprb), dimension(config%n_g_lw,nlev,istartcol:iendcol), intent(out) :: &
         &   od_lw
    real(jprb), dimension(config%n_g_sw,nlev,istartcol:iendcol), intent(out) :: &
         &   od_sw, ssa_sw

    ! The Planck function (emitted flux from a black body) at half
    ! levels at each longwave g-point
    real(jprb), dimension(config%n_g_lw,nlev+1,istartcol:iendcol), &
         &   intent(out), optional :: planck_hl
    ! Planck function for the surface (W m-2)
    real(jprb), dimension(config%n_g_lw,istartcol:iendcol), &
         &   intent(out), optional :: lw_emission

    ! The incoming shortwave flux into a plane perpendicular to the
    ! incoming radiation at top-of-atmosphere in each of the shortwave
    ! g-points
    real(jprb), dimension(config%n_g_sw,istartcol:iendcol), &
         &   intent(out), optional :: incoming_sw


  end subroutine gas_optics

  !---------------------------------------------------------------------
  ! Externally facing function for computing the Planck function
  ! without reference to any gas profile; typically this would be used
  ! for computing the emission by a surface.
  subroutine planck_function(config, temperature, planck_surf)

    use parkind1,                 only : jprb
    use radiation_config,         only : config_type

    type(config_type), intent(in) :: config
    real(jprb),        intent(in) :: temperature

    ! Planck function of the surface (W m-2)
    real(jprb), dimension(config%n_g_lw), intent(out) :: planck_surf

  end subroutine planck_function

end module radiation_ecckd_interface
