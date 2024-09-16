! ecrad3d_save.F90 - Save data to NetCDF files preserving 3D structure
!
! (C) Copyright 2014- ECMWF.
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

module ecrad3d_save
  
  use parkind1, only : jprb

  implicit none

  public :: save_fluxes

  contains


  !---------------------------------------------------------------------
  ! Save fluxes in "flux" to NetCDF file_name, plus pressure from the
  ! thermodynamics object
  subroutine save_fluxes(file_name, config, geometry, thermodynamics, flux, &
       &                 iverbose, is_hdf5_file, experiment_name, &
       &                 is_double_precision)

    use yomhook,                  only : lhook, dr_hook, jphook

    use easy_netcdf

    use radiation_io,             only : nulout
    use radiation_config,         only : config_type, IGasModelMonochromatic
    use ecrad3d_geometry,         only : geometry_type
    use radiation_thermodynamics, only : thermodynamics_type
    use radiation_flux,           only : flux_type
    use radiation_save,           only : save_fluxes_2d => save_fluxes 

    character(len=*),           intent(in) :: file_name
    type(config_type),          intent(in) :: config
    type(geometry_type),        intent(in) :: geometry
    type(thermodynamics_type),  intent(in) :: thermodynamics
    type(flux_type),            intent(in) :: flux
    integer,          optional, intent(in) :: iverbose
    logical,          optional, intent(in) :: is_hdf5_file
    logical,          optional, intent(in) :: is_double_precision
    character(len=*), optional, intent(in) :: experiment_name

    type(netcdf_file)                      :: out_file
    integer                                :: ncol, n_lev_plus1
    character(5), parameter                :: default_lw_units_str = 'W m-2'
    character(5)                           :: lw_units_str
    integer                                :: i_local_verbose
    
    ! Shape of output arrays
    integer, dimension(2) :: shape_2d
    integer, dimension(3) :: shape_3d, shape_lw_3d, shape_sw_3d, shape_canopy_lw_3d, shape_canopy_sw_3d
    integer, dimension(4) :: shape_lw_4d, shape_sw_4d

    ! Either "x" and "y" or "lon" and "lat"
    character(3)                           :: xname, yname
    
    real(jphook) :: hook_handle

    if (lhook) call dr_hook('ecrad3d_save:save_fluxes',0,hook_handle)

    if (.not. geometry%is_3d) then
      call save_fluxes_2d(file_name, config, thermodynamics, flux, &
       &                 iverbose, is_hdf5_file, experiment_name, &
       &                 is_double_precision)
      return
    end if
    
    if (present(iverbose)) then
      i_local_verbose = iverbose
    else
      i_local_verbose = config%iverbose
    end if

    ! Work out array dimensions
    if (config%do_sw) then
      ncol = size(flux%sw_up,1)
      n_lev_plus1 = size(flux%sw_up,2)
    elseif (config%do_lw) then
      ncol = size(flux%lw_up,1)
      n_lev_plus1 = size(flux%lw_up,2)
    else
      if (i_local_verbose >= 1) then
        write(nulout,'(a,a,a)') 'Warning: neither longwave nor shortwave computed so ', &
             &                  trim(file_name),' not written'
      end if
      return
    end if

    if (config%i_gas_model_lw == IGasModelMonochromatic &
         .and. config%mono_lw_wavelength > 0.0_jprb) then
      lw_units_str = 'W m-3'
    else
      lw_units_str = default_lw_units_str
    end if

    ! Open the file
    call out_file%create(trim(file_name), iverbose=i_local_verbose, is_hdf5_file=is_hdf5_file)

    ! Use compression if is an hdf5 file, for all variables with 3 or
    ! more dimensions using deflate-level 3 and shuffling
    if (is_hdf5_file) then
      call out_file%deflate_policy(3, 3, .true.)
    end if
    
    ! Variables stored internally with column varying fastest, but in
    ! output file column varies most slowly so need to transpose
    !call out_file%transpose_matrices(.true.)

    ! Set default precision for file, if specified
    if (present(is_double_precision)) then
      call out_file%double_precision(is_double_precision)
    end if

    ! Spectral fluxes in memory are dimensioned (nband,ncol,nlev), but
    ! are reoriented in the output file to be (nband,nlev,ncol), where
    ! the convention here is first dimension varying fastest
    !call out_file%permute_3d_arrays( (/ 1, 3, 2 /) )

    if (geometry%is_lat_lon) then
      xname = 'lon'
      yname = 'lat'
    else
      xname = 'x'
      yname = 'y'
    end if
    
    ! Define dimensions
    call out_file%define_dimension("half_level", n_lev_plus1)
    call out_file%define_dimension(trim(yname), geometry%ny)
    call out_file%define_dimension(trim(xname), geometry%nx)

    shape_2d = [geometry%nx, geometry%ny]
    shape_3d = [geometry%nx, geometry%ny, n_lev_plus1]
    
    if (config%do_save_spectral_flux .or. config%do_toa_spectral_flux) then
      if (config%do_lw) then
        call out_file%define_dimension("band_lw", config%n_spec_lw)
        shape_sw_3d = [config%n_spec_lw, geometry%nx, geometry%ny]
        shape_sw_4d = [config%n_spec_lw, geometry%nx, geometry%ny, n_lev_plus1]
      end if
      if (config%do_sw) then
        call out_file%define_dimension("band_sw", config%n_spec_sw)
        shape_sw_3d = [config%n_spec_sw, geometry%nx, geometry%ny]
        shape_sw_4d = [config%n_spec_sw, geometry%nx, geometry%ny, n_lev_plus1]
     end if
    else if (config%do_surface_sw_spectral_flux) then
      if (config%do_sw) then
        call out_file%define_dimension("band_sw", config%n_bands_sw)
        shape_sw_3d = [config%n_bands_sw, geometry%nx, geometry%ny]
        shape_sw_4d = [config%n_bands_sw, geometry%nx, geometry%ny, n_lev_plus1]
      end if
    end if

    if (config%do_lw .and. config%do_canopy_fluxes_lw) then
      call out_file%define_dimension("canopy_band_lw", &
           &  size(flux%lw_dn_surf_canopy, 1))
      shape_canopy_lw_3d = [size(flux%lw_dn_surf_canopy,1), geometry%nx, geometry%ny]
    end if
    if (config%do_sw .and. config%do_canopy_fluxes_sw) then
      call out_file%define_dimension("canopy_band_sw", &
           &  size(flux%sw_dn_diffuse_surf_canopy, 1))
      shape_canopy_sw_3d = [size(flux%sw_dn_diffuse_surf_canopy,1), geometry%nx, geometry%ny]
    end if

    ! Put global attributes
    call out_file%put_global_attributes( &
         &   title_str="Radiative flux profiles from the ecRad offline radiation model", &
         &   references_str="Hogan, R. J., and A. Bozzo, 2018: A flexible and efficient radiation " &
         &   //"scheme for the ECMWF model. J. Adv. Modeling Earth Sys., 10, 1990–2008", &
         &   source_str="ecRad offline radiation model")

    ! Save "experiment" global attribute if present and not empty
    if (present(experiment_name)) then
      if (experiment_name /= " ") then
        call out_file%put_global_attribute("experiment", experiment_name)
      end if
    end if
    
    ! Define variables
    call out_file%define_variable("pressure_hl", &
         &   dim3_name="half_level", dim2_name=trim(yname), dim1_name=trim(xname), &
         &   units_str="Pa", long_name="Pressure", &
         &   standard_name="air_pressure")

    if (config%do_lw) then
      call out_file%define_variable("flux_up_lw", &
           &   dim3_name="half_level", dim2_name=trim(yname), dim1_name=trim(xname), &
           &   units_str=lw_units_str, long_name="Upwelling longwave flux", &
           &   standard_name="upwelling_longwave_flux_in_air")
      call out_file%define_variable("flux_dn_lw", &
           &   dim3_name="half_level", dim2_name=trim(yname), dim1_name=trim(xname), &
           &   units_str=lw_units_str, long_name="Downwelling longwave flux", &
           &   standard_name="downwelling_longwave_flux_in_air")
      if (config%do_clear) then
        call out_file%define_variable("flux_up_lw_clear", &
             &   dim3_name="half_level", dim2_name=trim(yname), dim1_name=trim(xname), &
             &   units_str=lw_units_str, &
             &   long_name="Upwelling clear-sky longwave flux")
        call out_file%define_variable("flux_dn_lw_clear", &
             &   dim3_name="half_level", dim2_name=trim(yname), dim1_name=trim(xname), &
             &   units_str=lw_units_str, &
             &   long_name="Downwelling clear-sky longwave flux")
      end if

      if (config%do_lw_derivatives) then
        call out_file%define_variable("lw_derivative", &
             &  dim3_name="half_level", dim2_name=trim(yname), dim1_name=trim(xname), &
             &  units_str="1", &
             &  long_name="Derivative of upwelling LW flux w.r.t. surface value")
      end if

      if (config%do_save_spectral_flux) then
        call out_file%define_variable("spectral_flux_up_lw", &
             &   dim4_name="half_level", dim3_name=trim(yname), dim2_name=trim(xname), &
             &   dim1_name="band_lw", units_str=lw_units_str, &
             &   long_name="Spectral upwelling longwave flux")
        call out_file%define_variable("spectral_flux_dn_lw", &
             &   dim4_name="half_level", dim3_name=trim(yname), dim2_name=trim(xname), &
             &   dim1_name="band_lw", units_str=lw_units_str, &
             &   long_name="Spectral downwelling longwave flux")
        if (config%do_clear) then
          call out_file%define_variable("spectral_flux_up_lw_clear", &
               &   dim4_name="half_level", dim3_name=trim(yname), dim2_name=trim(xname), &
               &   dim1_name="band_lw", units_str=lw_units_str, &
               &   long_name="Spectral upwelling clear-sky longwave flux")
          call out_file%define_variable("spectral_flux_dn_lw_clear", &
               &   dim4_name="half_level", dim3_name=trim(yname), dim2_name=trim(xname), &
               &   dim1_name="band_lw", units_str=lw_units_str, &
               &   long_name="Spectral downwelling clear-sky longwave flux")
        end if
      end if
   
      if (config%do_toa_spectral_flux) then
        call out_file%define_variable("spectral_flux_up_lw_toa", &
             &   dim3_name=trim(yname), dim2_name=trim(xname), dim1_name="band_lw", units_str="W m-2", &
             &   long_name="Spectral upwelling longwave flux at top-of-atmosphere")
        if (config%do_clear) then
          call out_file%define_variable("spectral_flux_up_lw_toa_clear", &
               &   dim3_name=trim(yname), dim2_name=trim(xname), dim1_name="band_lw", units_str="W m-2", &
               &   long_name="Spectral upwelling clear-sky longwave flux at top-of-atmosphere")
        end if
      end if
   
      if (config%do_canopy_fluxes_lw) then
        call out_file%define_variable("canopy_flux_dn_lw_surf", &
             &   dim3_name=trim(yname), dim2_name=trim(xname), dim1_name="canopy_band_lw", units_str=lw_units_str, &
             &   long_name="Surface downwelling longwave flux in canopy bands")
      end if

    end if

    if (config%do_sw) then
      call out_file%define_variable("flux_up_sw", &
           &   dim3_name="half_level", dim2_name=trim(yname), dim1_name=trim(xname),  &
           &   units_str="W m-2", long_name="Upwelling shortwave flux", &
           &   standard_name="upwelling_shortwave_flux_in_air")
      call out_file%define_variable("flux_dn_sw", &
           &   dim3_name="half_level", dim2_name=trim(yname), dim1_name=trim(xname),  &
           &   units_str="W m-2", long_name="Downwelling shortwave flux", &
           &   standard_name="downwelling_shortwave_flux_in_air")
      if (config%do_sw_direct) then
        call out_file%define_variable("flux_dn_direct_sw", &
             &   dim3_name="half_level", dim2_name=trim(yname), dim1_name=trim(xname),  &
             &   units_str="W m-2", &
             &   long_name="Downwelling direct shortwave flux")
      end if
      if (config%do_clear) then
        call out_file%define_variable("flux_up_sw_clear", &
             &   dim3_name="half_level", dim2_name=trim(yname), dim1_name=trim(xname),  &
             &   units_str="W m-2", &
             &   long_name="Upwelling clear-sky shortwave flux")
        call out_file%define_variable("flux_dn_sw_clear", &
             &   dim3_name="half_level", dim2_name=trim(yname), dim1_name=trim(xname),  &
             &   units_str="W m-2", &
             &   long_name="Downwelling clear-sky shortwave flux")
        if (config%do_sw_direct) then
          call out_file%define_variable("flux_dn_direct_sw_clear", &
               &   dim3_name="half_level", dim2_name=trim(yname), dim1_name=trim(xname),  &
               &   units_str="W m-2", &
               &   long_name="Downwelling clear-sky direct shortwave flux")
        end if
      end if

      if (config%do_save_spectral_flux) then
        call out_file%define_variable("spectral_flux_up_sw", &
             &   dim4_name="half_level", dim3_name=trim(yname), dim2_name=trim(xname), &
             &   dim1_name="band_sw", units_str="W m-2", &
             &   long_name="Spectral upwelling shortwave flux")
        call out_file%define_variable("spectral_flux_dn_sw", &
             &   dim4_name="half_level", dim3_name=trim(yname), dim2_name=trim(xname), &
             &   dim1_name="band_sw", units_str="W m-2", &
             &   long_name="Spectral downwelling shortwave flux")
        if (config%do_sw_direct) then
          call out_file%define_variable("spectral_flux_dn_direct_sw", &
               &   dim4_name="half_level", dim3_name=trim(yname), dim2_name=trim(xname), &
               &   dim1_name="band_sw", units_str="W m-2", &
               &   long_name="Spectral downwelling direct shortwave flux")
        end if
        if (config%do_clear) then
          call out_file%define_variable("spectral_flux_up_sw_clear", &
               &   dim4_name="half_level", dim3_name=trim(yname), dim2_name=trim(xname), &
               &   dim1_name="band_sw", units_str="W m-2", &
               &   long_name="Spectral upwelling clear-sky shortwave flux")
          call out_file%define_variable("spectral_flux_dn_sw_clear", &
               &   dim4_name="half_level", dim3_name=trim(yname), dim2_name=trim(xname), &
               &   dim1_name="band_sw", units_str="W m-2", &
               &   long_name="Spectral downwelling clear-sky shortwave flux")
          if (config%do_sw_direct) then
            call out_file%define_variable("spectral_flux_dn_direct_sw_clear", &
                 &   dim4_name="half_level", dim3_name=trim(yname), dim2_name=trim(xname), &
                 &   dim1_name="band_sw", units_str="W m-2", &
                 &   long_name="Spectral downwelling clear-sky direct shortwave flux")
          end if
        end if
      else if (config%do_surface_sw_spectral_flux) then
        call out_file%define_variable("spectral_flux_dn_sw_surf", &
             &   dim3_name=trim(yname), dim2_name=trim(xname), dim1_name="band_sw", units_str="W m-2", &
             &   long_name="Spectral downwelling shortwave flux at surface")
        call out_file%define_variable("spectral_flux_dn_direct_sw_surf", &
             &   dim3_name=trim(yname), dim2_name=trim(xname), dim1_name="band_sw", units_str="W m-2", &
             &   long_name="Spectral downwelling direct shortwave flux at surface")
        if (config%do_clear) then
          call out_file%define_variable("spectral_flux_dn_sw_surf_clear", &
               &   dim3_name=trim(yname), dim2_name=trim(xname), dim1_name="band_sw", units_str="W m-2", &
               &   long_name="Spectral downwelling clear-sky shortwave flux at surface")
          call out_file%define_variable("spectral_flux_dn_direct_sw_surf_clear", &
               &   dim3_name=trim(yname), dim2_name=trim(xname), dim1_name="band_sw", units_str="W m-2", &
               &   long_name="Spectral downwelling clear-sky direct shortwave flux at surface")
        end if
      end if

      if (config%do_toa_spectral_flux) then
        call out_file%define_variable("spectral_flux_dn_sw_toa", &
             &   dim3_name=trim(yname), dim2_name=trim(xname), dim1_name="band_sw", units_str="W m-2", &
             &   long_name="Spectral downwelling shortwave flux at top-of-atmosphere")
        call out_file%define_variable("spectral_flux_up_sw_toa", &
             &   dim3_name=trim(yname), dim2_name=trim(xname), dim1_name="band_sw", units_str="W m-2", &
             &   long_name="Spectral upwelling shortwave flux at top-of-atmosphere")
        if (config%do_clear) then
          call out_file%define_variable("spectral_flux_up_sw_toa_clear", &
               &   dim3_name=trim(yname), dim2_name=trim(xname), dim1_name="band_sw", units_str="W m-2", &
               &   long_name="Spectral upwelling clear-sky shortwave flux at top-of-atmosphere")
        end if
      end if
   
      if (config%do_canopy_fluxes_sw) then
        call out_file%define_variable("canopy_flux_dn_diffuse_sw_surf", &
             &   dim3_name=trim(yname), dim2_name=trim(xname), dim1_name="canopy_band_sw", units_str="W m-2", &
             &   long_name="Surface downwelling diffuse shortwave flux in canopy bands")
        call out_file%define_variable("canopy_flux_dn_direct_sw_surf", &
             &   dim3_name=trim(yname), dim2_name=trim(xname), dim1_name="canopy_band_sw", units_str="W m-2", &
             &   long_name="Surface downwelling direct shortwave flux in canopy bands")
      end if

    end if
   
    if (config%do_lw .and. config%do_clouds) then
      call out_file%define_variable("cloud_cover_lw", &
           &  dim2_name=trim(yname), dim1_name=trim(xname), units_str="1", &
           &  long_name="Total cloud cover diagnosed by longwave solver", &
           &  standard_name="cloud_area_fraction")
    end if
    if (config%do_sw .and. config%do_clouds) then
      call out_file%define_variable("cloud_cover_sw", &
           &  dim2_name=trim(yname), dim1_name=trim(xname), units_str="1", &
           &  long_name="Total cloud cover diagnosed by shortwave solver", &
           &  standard_name="cloud_area_fraction")
    end if

    ! Write variables

    call out_file%put("pressure_hl", reshape(thermodynamics%pressure_hl, shape_3d))

    if (config%do_lw) then
      call out_file%put("flux_up_lw", reshape(flux%lw_up, shape_3d))
      call out_file%put("flux_dn_lw", reshape(flux%lw_dn, shape_3d))
      if (config%do_clear) then
        call out_file%put("flux_up_lw_clear", reshape(flux%lw_up_clear, shape_3d))
        call out_file%put("flux_dn_lw_clear", reshape(flux%lw_dn_clear, shape_3d))
      end if

      if (config%do_lw_derivatives) then
        call out_file%put("lw_derivative", reshape(flux%lw_derivatives, shape_3d))
      end if

      if (config%do_save_spectral_flux) then
        call out_file%put("spectral_flux_up_lw", reshape(flux%lw_up_band, shape_lw_4d))
        call out_file%put("spectral_flux_dn_lw", reshape(flux%lw_dn_band, shape_lw_4d))
        if (config%do_clear) then
          call out_file%put("spectral_flux_up_lw_clear", reshape(flux%lw_up_clear_band, shape_lw_4d))
          call out_file%put("spectral_flux_dn_lw_clear", reshape(flux%lw_dn_clear_band, shape_lw_4d))
        end if
      end if

      if (config%do_toa_spectral_flux) then
        call out_file%put("spectral_flux_up_lw_toa", reshape(flux%lw_up_toa_band, shape_lw_3d))
        if (config%do_clear) then
          call out_file%put("spectral_flux_up_lw_toa_clear", reshape(flux%lw_up_toa_clear_band, shape_lw_3d))
        end if
      end if
      
      if (config%do_canopy_fluxes_lw) then
        call out_file%put("canopy_flux_dn_lw_surf", reshape(flux%lw_dn_surf_canopy, shape_canopy_lw_3d))
      end if

    end if

    if (config%do_sw) then
      call out_file%put("flux_up_sw", reshape(flux%sw_up, shape_3d))
      call out_file%put("flux_dn_sw", reshape(flux%sw_dn, shape_3d))
      if (config%do_sw_direct) then
        call out_file%put("flux_dn_direct_sw", reshape(flux%sw_dn_direct, shape_3d))
      end if
      if (config%do_clear) then
        call out_file%put("flux_up_sw_clear", reshape(flux%sw_up_clear, shape_3d))
        call out_file%put("flux_dn_sw_clear", reshape(flux%sw_dn_clear, shape_3d))
        if (config%do_sw_direct) then
          call out_file%put("flux_dn_direct_sw_clear", reshape(flux%sw_dn_direct_clear, shape_3d))
        end if
      end if

      if (config%do_save_spectral_flux) then
        call out_file%put("spectral_flux_up_sw", reshape(flux%sw_up_band, shape_sw_4d))
        call out_file%put("spectral_flux_dn_sw", reshape(flux%sw_dn_band, shape_sw_4d))
        if (config%do_sw_direct) then
          call out_file%put("spectral_flux_dn_direct_sw", &
               &   reshape(flux%sw_dn_direct_band, shape_sw_4d))
        end if
        if (config%do_clear) then
          call out_file%put("spectral_flux_up_sw_clear", reshape(flux%sw_up_clear_band, shape_sw_4d))
          call out_file%put("spectral_flux_dn_sw_clear", reshape(flux%sw_dn_clear_band, shape_sw_4d))
          if (config%do_sw_direct) then
            call out_file%put("spectral_flux_dn_direct_sw_clear", &
                 &   reshape(flux%sw_dn_direct_clear_band, shape_sw_4d))
          end if
        end if
      else if (config%do_surface_sw_spectral_flux) then
        call out_file%put("spectral_flux_dn_sw_surf", reshape(flux%sw_dn_surf_band, shape_sw_3d))
        call out_file%put("spectral_flux_dn_direct_sw_surf", reshape(flux%sw_dn_direct_surf_band, shape_sw_3d))
        if (config%do_clear) then
          call out_file%put("spectral_flux_dn_sw_surf_clear", reshape(flux%sw_dn_surf_clear_band, shape_sw_3d))
          call out_file%put("spectral_flux_dn_direct_sw_surf_clear", &
               &            reshape(flux%sw_dn_direct_surf_clear_band, shape_sw_3d))
        end if
      end if

      if (config%do_toa_spectral_flux) then
        call out_file%put("spectral_flux_dn_sw_toa", reshape(flux%sw_dn_toa_band, shape_sw_3d))
        call out_file%put("spectral_flux_up_sw_toa", reshape(flux%sw_up_toa_band, shape_sw_3d))
        if (config%do_clear) then
          call out_file%put("spectral_flux_up_sw_toa_clear", reshape(flux%sw_up_toa_clear_band, shape_sw_3d))
        end if
      end if
      
      if (config%do_canopy_fluxes_sw) then
        call out_file%put("canopy_flux_dn_diffuse_sw_surf", &
             &            reshape(flux%sw_dn_diffuse_surf_canopy, shape_canopy_sw_3d))
        call out_file%put("canopy_flux_dn_direct_sw_surf", &
             &            reshape(flux%sw_dn_direct_surf_canopy, shape_canopy_sw_3d))
      end if

    end if

    if (config%do_lw .and. config%do_clouds) then
      call out_file%put("cloud_cover_lw", reshape(flux%cloud_cover_lw, shape_2d))
    end if
    if (config%do_sw .and. config%do_clouds) then
      call out_file%put("cloud_cover_sw", reshape(flux%cloud_cover_sw, shape_2d))
    end if

    ! Close file
    call out_file%close()

    if (lhook) call dr_hook('ecrad3d_save:save_fluxes',1,hook_handle)

  end subroutine save_fluxes


  !---------------------------------------------------------------------
  ! Save radiances in "flux" to NetCDF file_name, plus
  ! cos_sensor_zenith_angle from single_level object
  subroutine save_radiances(file_name, config, geometry, single_level, flux, &
       &                    iverbose, is_hdf5_file, experiment_name)

    use yomhook,                  only : lhook, dr_hook, jphook

    use easy_netcdf

    use radiation_io,             only : nulout
    use radiation_config,         only : config_type, IGasModelMonochromatic
    use ecrad3d_geometry,         only : geometry_type
    use radiation_single_level,   only : single_level_type
    use radiation_flux,           only : flux_type

    character(len=*),           intent(in) :: file_name
    type(config_type),          intent(in) :: config
    type(geometry_type),        intent(in) :: geometry
    type(single_level_type),    intent(in) :: single_level
    type(flux_type),            intent(in) :: flux
    integer,          optional, intent(in) :: iverbose
    logical,          optional, intent(in) :: is_hdf5_file
    character(len=*), optional, intent(in) :: experiment_name

    type(netcdf_file)                      :: out_file
    integer                                :: ncol, n_bands_lw, n_bands_sw
    character(10), parameter               :: default_lw_units_str = 'W m-2 sr-1'
    character(10), parameter               :: default_sw_units_str = 'W m-2 sr-1'
    character(10)                          :: lw_units_str, sw_units_str
    integer                                :: i_local_verbose

    ! Shape of output arrays
    integer, dimension(2) :: shape_2d
    integer, dimension(3) :: shape_lw_3d, shape_sw_3d

    ! Either "x" and "y" or "lon" and "lat"
    character(3)                           :: xname, yname
    
    real(jphook) :: hook_handle

    if (lhook) call dr_hook('ecrad3d_save:save_radiances',0,hook_handle)
    
    if (present(iverbose)) then
      i_local_verbose = iverbose
    else
      i_local_verbose = config%iverbose
    end if

    lw_units_str = default_lw_units_str
    sw_units_str = default_sw_units_str

    ncol = size(single_level%cos_sensor_zenith_angle)

    ! Open the file
    call out_file%create(trim(file_name), iverbose=i_local_verbose, is_hdf5_file=is_hdf5_file)

    ! Use compression if is an hdf5 file, for all variables with 3 or
    ! more dimensions using deflate-level 3 and shuffling
    if (is_hdf5_file) then
      call out_file%deflate_policy(3, 3, .true.)
    end if

    if (geometry%is_lat_lon) then
      xname = 'lon'
      yname = 'lat'
    else
      xname = 'x'
      yname = 'y'
    end if
    
    ! Define dimensions
    call out_file%define_dimension(trim(yname), geometry%ny)
    call out_file%define_dimension(trim(xname), geometry%nx)

    shape_2d = [geometry%nx, geometry%ny]

    if (config%do_lw) then
      n_bands_lw = size(flux%lw_radiance_band,1)
      call out_file%define_dimension("band_lw", n_bands_lw)
      shape_lw_3d = [n_bands_lw, geometry%nx, geometry%ny]
    end if
    if (config%do_sw) then
      n_bands_sw = size(flux%sw_radiance_band,1)
      call out_file%define_dimension("band_sw", n_bands_sw)
      shape_sw_3d = [n_bands_sw, geometry%nx, geometry%ny]
    end if

   ! Put global attributes
    call out_file%put_global_attributes( &
         &   title_str="Radiances from the ecRad offline radiation model", &
         &   references_str="Hogan, R. J., and A. Bozzo, 2018: A flexible and efficient radiation " &
         &   //"scheme for the ECMWF model. J. Adv. Modeling Earth Sys., 10, 1990–2008", &
         &   source_str="ecRad offline radiation model")

    ! Save "experiment" global attribute if present and not empty
    if (present(experiment_name)) then
      if (experiment_name /= " ") then
        call out_file%put_global_attribute("experiment", experiment_name)
      end if
    end if

    ! Define variables
    call out_file%define_variable("cos_sensor_zenith_angle", &
         &  dim2_name=trim(yname), dim1_name=trim(xname), units_str="1", &
         &  long_name="Cosine of the sensor zenith angle")

    if (config%do_lw) then

      if (flux%is_brightness_temperature) then
        call out_file%define_variable("brightness_temperature_lw_band", &
             &  dim3_name=trim(yname), dim2_name=trim(xname), dim1_name="band_lw", &
             &  units_str="K", long_name="Brightness temperature")
      else
        call out_file%define_variable("radiance_lw_band", &
             &  dim3_name=trim(yname), dim2_name=trim(xname), dim1_name="band_lw", &
             &  units_str=trim(lw_units_str), long_name="Thermal radiance")
      end if
      call out_file%define_variable("cloud_cover_lw", &
           &  dim2_name=trim(yname), dim1_name=trim(xname), units_str="1", &
           &  long_name="Total cloud cover diagnosed by longwave solver", &
           &  standard_name="cloud_area_fraction")

    end if

    if (config%do_sw) then

      if (.not. config%do_save_spectral_flux &
           &  .and. size(config%gas_optics_sw%spectral_def%wavenumber1_band) == n_bands_sw) then
        call out_file%define_variable("wavenumber1_sw", &
             &  dim1_name="band_sw", units_str="cm-1", &
             &  long_name="Lower wavenumber of shortwave band")
        call out_file%define_variable("wavenumber2_sw", &
             &  dim1_name="band_sw", units_str="cm-1", &
             &  long_name="Upper wavenumber of shortwave band")
      end if

      call out_file%define_variable("radiance_sw_band", &
           &  dim3_name=trim(yname), dim2_name=trim(xname), dim1_name="band_sw",&
           &  units_str=trim(sw_units_str), long_name="Solar radiance")

      if (config%do_clear) then
        call out_file%define_variable("radiance_sw_clear_band", &
             &  dim3_name=trim(yname), dim2_name=trim(xname), dim1_name="band_sw", &
             &  units_str=trim(sw_units_str), long_name="Clear-sky solar radiance")
      end if

    end if

    ! Write variables
    call out_file%put("cos_sensor_zenith_angle", reshape(single_level%cos_sensor_zenith_angle, shape_2d))
    
    if (config%do_lw) then
      !if (.not. config%do_save_spectral_flux) then
      !  call out_file%put("wavenumber1_lw", config%gas_optics_lw%spectral_def%wavenumber1_band)
      !  call out_file%put("wavenumber2_lw", config%gas_optics_lw%spectral_def%wavenumber2_band)
      !end if

      if (flux%is_brightness_temperature) then
        call out_file%put("brightness_temperature_lw_band", reshape(flux%lw_radiance_band, shape_lw_3d))
      else
        call out_file%put("radiance_lw_band", reshape(flux%lw_radiance_band, shape_lw_3d))
      end if
      call out_file%put("cloud_cover_lw", reshape(flux%cloud_cover_lw, shape_2d))
    end if

    if (config%do_sw) then
      if (.not. config%do_save_spectral_flux &
           &  .and. size(config%gas_optics_sw%spectral_def%wavenumber1_band) == n_bands_sw) then
        call out_file%put("wavenumber1_sw", config%gas_optics_sw%spectral_def%wavenumber1_band)
        call out_file%put("wavenumber2_sw", config%gas_optics_sw%spectral_def%wavenumber2_band)
      end if

      call out_file%put("radiance_sw_band", reshape(flux%sw_radiance_band, shape_sw_3d))
      if (config%do_clear) then
        call out_file%put("radiance_sw_clear_band", reshape(flux%sw_radiance_clear_band, shape_sw_3d))
      end if
    end if

    ! Close file
    call out_file%close()

    if (lhook) call dr_hook('ecrad3d_save:save_radiances',1,hook_handle)

  end subroutine save_radiances
 
    
end module ecrad3d_save
