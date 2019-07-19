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
  ! properties for a particular type of cloud or hydrometeor
  type general_cloud_optics_type
    ! Band-specific coefficients are provided separately in the
    ! shortwave and longwave, and are dimensioned (nband,ncoeff),
    ! where ncoeff depends on the nature of the parameterization but
    ! is usually a look-up table as a function of effective radius
    
    ! Extinction coefficient per unit mass (m2 kg-1)
    real(jprb), allocatable, dimension(:,:) :: &
         &  mass_ext_coeff_lw, mass_ext_coeff_sw
    
    ! Single-scattering albedo and asymmetry factor (dimensionless)
    real(jprb), allocatable, dimension(:,:) :: &
         &  ssa_lw, ssa_sw, g_lw, g_sw

    ! Number of effective radius coefficients, start value and
    ! interval in look-up table
    integer    :: n_effective_radius
    real(jprb) :: effective_radius_0, d_effective_radius

    ! Name of cloud/precip type (e.g. "liquid", "ice", "rain", "snow")
    ! and the name of the optics scheme.  These two are used to
    ! generate the name of the data file from which the coefficients
    ! are read.
    character(len=511) :: type_name, scheme_name
    
   contains
     procedure :: setup => setup_general_cloud_optics

  end type general_cloud_optics_type

contains

  !---------------------------------------------------------------------
  ! Setup cloud optics coefficients by reading them from a file
  subroutine setup_general_cloud_optics(this, liq_file_name, ice_file_name, &
       &                        iverbose)

    use yomhook,  only : lhook, dr_hook
    use easy_netcdf, only : netcdf_file

    class(cloud_optics_type), intent(inout) :: this
    character(len=*), intent(in)            :: liq_file_name, ice_file_name
    integer, intent(in), optional           :: iverbose

    ! The NetCDF file containing the coefficients
    type(netcdf_file)  :: file
    integer            :: iverb
    real(jprb)         :: hook_handle

    if (lhook) call dr_hook('radiation_general_cloud_optics_data:setup',0,hook_handle)

    if (present(iverbose)) then
      iverb = iverbose
    else
      iverb = 2
    end if

    ! Open the droplet scattering file and configure the way it is
    ! read
    call file%open(trim(liq_file_name), iverbose=iverb)
    call file%transpose_matrices()

    ! Read the band-specific coefficients
    call file%get('coeff_lw',this%liq_coeff_lw)
    call file%get('coeff_sw',this%liq_coeff_sw)

    ! Read the general  coefficients
    if (file%exists('coeff_gen')) then
      call file%get('coeff_gen',this%liq_coeff_gen)
    end if

    ! Close droplet scattering file
    call file%close()

    ! Open the ice scattering file and configure the way it is read
    call file%open(trim(ice_file_name), iverbose=iverb)
    call file%transpose_matrices()

    ! Read the band-specific  coefficients
    call file%get('coeff_lw',this%ice_coeff_lw)
    call file%get('coeff_sw',this%ice_coeff_sw)

    ! Read the general  coefficients
    if (file%exists('coeff_gen')) then
      call file%get('coeff_gen',this%ice_coeff_gen)
    end if

    ! Close ice scattering file
    call file%close()

    if (lhook) call dr_hook('radiation_general_cloud_optics_data:setup',1,hook_handle)

  end subroutine setup_general_cloud_optics

end module radiation_general_cloud_optics_data
