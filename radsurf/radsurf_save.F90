! radsurf_save.f90 - Save surface data to NetCDF files
!
! Copyright (C) 2017 ECMWF
!
! Author:  Robin Hogan
! Email:   r.j.hogan@ecmwf.int
! License: see the COPYING file for details
!

module radsurf_save

  public :: save_surface_fluxes

contains

  !------------------------------------------------------------------
  ! Save fluxes in "surface_flux" to NetCDF file_name
  subroutine save_surface_fluxes(file_name, config, surface_flux, iverbose)

    use yomhook,                  only : lhook, dr_hook

    use easy_netcdf

    use radiation_io,             only : nulout
    use radiation_config,         only : config_type, IGasModelMonochromatic
    use radsurf_flux,             only : surface_flux_type

    character(len=*),        intent(in)  :: file_name
    type(config_type),       intent(in)  :: config
    type(surface_flux_type), intent(in)  :: surface_flux

    type(netcdf_file)                    :: out_file
    integer                              :: ncol

    character(5), parameter :: default_lw_units_str = 'W m-2'
    character(5)            :: lw_units_str

    integer, optional, intent(in) :: iverbose
    integer                       :: i_local_verbose

    integer :: nfacet, ntile

    real(jprb) :: hook_handle

    if (lhook) call dr_hook('radsurf_save:save_surface_fluxes',0,hook_handle)

    if (present(iverbose)) then
      i_local_verbose = iverbose
    else
      i_local_verbose = config%iverbose
    end if

    ncol   = 0
    ntile  = 0
    nfacet = 0
    if (allocated(surface_flux%lw_dn_facet)) then
      ncol   = size(surface_flux%lw_dn_facet,1)
      nfacet = size(surface_flux%lw_dn_facet,2)
      ntile  = size(surface_flux%lw_abs_canopy,2)
    else if (allocated(surface_flux%sw_dn_facet)) then
      ncol   = size(surface_flux%sw_dn_facet,1)
      nfacet = size(surface_flux%sw_dn_facet,2)
      ntile  = size(surface_flux%sw_abs_canopy,2)
    end if

    if (ncol == 0) then
      write(nulout,'(a)') 'Warning: surface-flux structure contains no data to write'
    else

      if (config%i_gas_model == IGasModelMonochromatic &
           .and. config%mono_lw_wavelength > 0.0_jprb) then
        lw_units_str = 'W m-3'
      else
        lw_units_str = default_lw_units_str
      end if

      ! Open the file
      call out_file%create(trim(file_name), iverbose=i_local_verbose)

      ! Variables stored internally with column varying fastest, but in
      ! output file column varies most slowly so need to transpose
      call out_file%transpose_matrices(.true.)

      ! Define dimensions
      call out_file%define_dimension("column", ncol)
      if (ntile > 0) then
        call out_file%define_dimension("tile",  ntile)
        call out_file%define_dimension("facet", nfacet)
      end if

      if (allocated(surface_flux%lw_dn_facet)) then
        call out_file%define_variable("flux_dn_lw_facet", &
             &  dim2_name="column", dim1_name="facet", units_str=lw_units_str, &
             &  long_name="Longwave flux into facet of surface")
        call out_file%define_variable("flux_up_lw_facet", &
             &  dim2_name="column", dim1_name="facet", units_str=lw_units_str, &
             &  long_name="Longwave flux out of facet of surface")
        call out_file%define_variable("absorption_lw_canopy", &
             &  dim2_name="column", dim1_name="tile", units_str=lw_units_str, &
             &  long_name="Longwave absorption by tile canopy")
      end if
      if (allocated(surface_flux%sw_dn_facet)) then
        call out_file%define_variable("flux_dn_sw_facet", &
             &  dim2_name="column", dim1_name="facet", units_str="W m-2", &
             &  long_name="Shortwave flux into facet of surface")
        call out_file%define_variable("flux_dn_direct_sw_facet", &
             &  dim2_name="column", dim1_name="facet", units_str="W m-2", &
             &  long_name="Shortwave direct flux into facet of surface")
        call out_file%define_variable("flux_up_sw_facet", &
             &  dim2_name="column", dim1_name="facet", units_str="W m-2", &
             &  long_name="Shortwave flux out of facet of surface")
        call out_file%define_variable("absorption_sw_canopy", &
             &  dim2_name="column", dim1_name="tile", units_str="W m-2", &
             &  long_name="Shortwave absorption by tile canopy")
      end if

      if (allocated(surface_flux%lw_dn_facet)) then
        call out_file%put("flux_dn_lw_facet", surface_flux%lw_dn_facet)
        call out_file%put("flux_up_lw_facet", surface_flux%lw_up_facet)
        call out_file%put("absorption_lw_canopy", surface_flux%lw_abs_canopy)
      end if
      if (allocated(surface_flux%sw_dn_facet)) then
        call out_file%put("flux_dn_sw_facet", surface_flux%sw_dn_facet)
        call out_file%put("flux_dn_direct_sw_facet", surface_flux%sw_dn_direct_facet)
        call out_file%put("flux_up_sw_facet", surface_flux%sw_up_facet)
        call out_file%put("absorption_sw_canopy", surface_flux%sw_abs_canopy)
      end if
      
      ! Close file
      call out_file%close()

    end if
    if (lhook) call dr_hook('radsurf_save:save_surface_fluxes',1,hook_handle)
      
    end subroutine save_surface_fluxes
    
end module radsurf_save
