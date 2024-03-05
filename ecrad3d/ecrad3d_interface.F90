! ecrad3d_interface.F90



module ecrad3d_interface

  
contains

  subroutine ecrad3d(config, geometry, &
       &  single_level, thermodynamics, gas, cloud, aerosol, flux)

    use parkind1,                 only : jprb
    use yomhook,                  only : lhook, dr_hook, jphook

    use radiation_io,             only : nulout, nulerr, radiation_abort
    use radiation_config,         only : config_type, &
         &  IGasModelECCKD, IGasModelIFSRRTMG, IGasModelMonochromatic
    use radiation_single_level,   only : single_level_type
    use radiation_thermodynamics, only : thermodynamics_type
    use radiation_gas,            only : gas_type
    use radiation_cloud,          only : cloud_type
    use radiation_aerosol,        only : aerosol_type
    use ecrad3d_geometry,         only : geometry_type
    use ecrad3d_solver_interface, only : solver_interface_sw
    use radiation_flux,           only : flux_type
    use radiation_monochromatic,  only : &
         &   gas_optics_mono     => gas_optics, &
         &   cloud_optics_mono   => cloud_optics
    use radiation_ifs_rrtm,       only : gas_optics_rrtmg => gas_optics
    use radiation_ecckd_interface,only : gas_optics_ecckd => gas_optics
    use radiation_cloud_optics,   only : cloud_optics
    use radiation_general_cloud_optics, only : general_cloud_optics
    use radiation_aerosol_optics, only : add_aerosol_optics

    implicit none
    
    type(config_type),        intent(in)   :: config
    type(geometry_type),      intent(in)   :: geometry
    type(single_level_type),  intent(in)   :: single_level
    type(thermodynamics_type),intent(in)   :: thermodynamics
    type(gas_type),           intent(in)   :: gas
    type(cloud_type),         intent(inout):: cloud
    type(aerosol_type),       intent(in)   :: aerosol
    ! Output
    type(flux_type),          intent(inout):: flux


    ! Local variables

    ! Layer optical depth, single scattering albedo and asymmetry factor of
    ! gases and aerosols at each longwave g-point, where the latter
    ! two variables are only defined if aerosol longwave scattering is
    ! enabled (otherwise both are treated as zero).
    real(jprb), dimension(config%n_g_lw,geometry%nz,geometry%ncol) :: od_lw
    real(jprb), dimension(config%n_g_lw_if_scattering,geometry%nz,geometry%ncol) :: &
         &  ssa_lw, g_lw

    ! Layer in-cloud optical depth, single scattering albedo and
    ! asymmetry factor (or phase function components) of hydrometeors
    ! in each longwave band, where the latter two variables are only
    ! defined if hydrometeor longwave scattering is enabled (otherwise
    ! both are treated as zero).
    real(jprb), dimension(config%n_bands_lw,geometry%nz,geometry%ncol) :: od_lw_cloud
    real(jprb), dimension(config%n_bands_lw_if_scattering,geometry%nz,geometry%ncol) :: &
         &  ssa_lw_cloud
    real(jprb), dimension(config%n_bands_lw_if_scattering,geometry%nz,geometry%ncol,config%n_pf_lw) &
         &  :: pf_lw_cloud

    ! Layer optical depth, single scattering albedo and asymmetry factor of
    ! gases and aerosols at each shortwave g-point
    real(jprb), dimension(config%n_g_sw,geometry%nz,geometry%ncol) :: od_sw, ssa_sw, g_sw

    ! Layer in-cloud optical depth and single scattering albedo of
    ! hydrometeors in each shortwave band
    real(jprb), dimension(config%n_bands_sw,geometry%nz,geometry%ncol)   :: &
         &  od_sw_cloud, ssa_sw_cloud

    ! Variables describing the shortwave phase function; can be just
    ! asymmetry factor or extra variables if FLOTSAM/DISORT is being used
    real(jprb), dimension(config%n_bands_sw,geometry%nz,geometry%ncol,config%n_pf_sw) &
         &  :: pf_sw_cloud

    ! The Planck function (emitted flux from a black body) at half
    ! levels
    real(jprb), dimension(config%n_g_lw,geometry%nz+1,geometry%ncol) :: planck_hl

    ! The longwave emission from and albedo of the surface in each
    ! longwave g-point; note that these are weighted averages of the
    ! values from individual tiles
    real(jprb), dimension(config%n_g_lw,geometry%ncol) :: lw_emission
    real(jprb), dimension(config%n_g_lw,geometry%ncol) :: lw_albedo

    ! Direct and diffuse shortwave surface albedo in each shortwave
    ! g-point; note that these are weighted averages of the values
    ! from individual tiles
    real(jprb), dimension(config%n_g_sw,geometry%ncol) :: sw_albedo_direct
    real(jprb), dimension(config%n_g_sw,geometry%ncol) :: sw_albedo_diffuse

    ! The incoming shortwave flux into a plane perpendicular to the
    ! incoming radiation at top-of-atmosphere in each of the shortwave
    ! g-points
    real(jprb), dimension(config%n_g_sw,geometry%ncol) :: incoming_sw

    integer :: istartcol, iendcol
    
    real(jphook) :: hook_handle

    if (lhook) call dr_hook('ecrad3d_interface:ecrad3d',0,hook_handle)

    istartcol = 1
    iendcol   = geometry%ncol

    ! Extract surface albedos at each gridpoint
    call single_level%get_albedos(istartcol, iendcol, config, &
         &                        sw_albedo_direct, sw_albedo_diffuse, &
         &                        lw_albedo)


    ! Compute gas absorption optical depth in shortwave and
    ! longwave, shortwave single scattering albedo (i.e. fraction of
    ! extinction due to Rayleigh scattering), Planck functions and
    ! incoming shortwave flux at each g-point, for the specified
    ! range of atmospheric columns
    if (config%i_gas_model_sw == IGasModelMonochromatic) then
      call gas_optics_mono(geometry%ncol, geometry%nz, istartcol, iendcol, config, &
           &  single_level, thermodynamics, gas, lw_albedo, &
           &  od_lw, od_sw, ssa_sw, &
           &  planck_hl, lw_emission, incoming_sw)
    else
      if (config%i_gas_model_sw == IGasModelIFSRRTMG &
           &   .or. config%i_gas_model_lw == IGasModelIFSRRTMG) then
        call gas_optics_rrtmg(geometry%ncol, geometry%nz, istartcol, iendcol, config, &
             &  single_level, thermodynamics, gas, &
             &  od_lw, od_sw, ssa_sw, lw_albedo=lw_albedo, &
             &  planck_hl=planck_hl, lw_emission=lw_emission, &
             &  incoming_sw=incoming_sw)
      end if
      if (config%i_gas_model_sw == IGasModelECCKD &
           &   .or. config%i_gas_model_lw == IGasModelECCKD) then
        call gas_optics_ecckd(geometry%ncol, geometry%nz, istartcol, iendcol, config, &
             &  single_level, thermodynamics, gas, &
             &  od_lw, od_sw, ssa_sw, lw_albedo=lw_albedo, &
             &  planck_hl=planck_hl, lw_emission=lw_emission, &
             &  incoming_sw=incoming_sw)
      end if
    end if

    if (config%do_clouds) then
      ! Compute hydrometeor absorption/scattering properties in each
      ! shortwave and longwave band
      if (config%i_gas_model_sw == IGasModelMonochromatic) then
        call cloud_optics_mono(geometry%nz, istartcol, iendcol, &
             &  config, thermodynamics, cloud, &
             &  od_lw_cloud, ssa_lw_cloud, pf_lw_cloud, &
             &  od_sw_cloud, ssa_sw_cloud, pf_sw_cloud)
      elseif (config%use_general_cloud_optics) then
        call general_cloud_optics(geometry%nz, istartcol, iendcol, &
             &  config, single_level, thermodynamics, cloud, & 
             &  od_lw_cloud, ssa_lw_cloud, pf_lw_cloud, &
             &  od_sw_cloud, ssa_sw_cloud, pf_sw_cloud)
      else
        call cloud_optics(geometry%nz, istartcol, iendcol, &
             &  config, thermodynamics, cloud, & 
             &  od_lw_cloud, ssa_lw_cloud, pf_lw_cloud, &
             &  od_sw_cloud, ssa_sw_cloud, pf_sw_cloud)
      end if
    end if
    
    if (config%use_aerosols) then
      if (config%i_gas_model_sw == IGasModelMonochromatic) then
!          call add_aerosol_optics_mono(nlev,istartcol,iendcol, &
!               &  config, thermodynamics, gas, aerosol, & 
!               &  od_lw, ssa_lw, g_lw, od_sw, ssa_sw, g_sw)
      else
        call add_aerosol_optics(geometry%nz, istartcol, iendcol, &
             &  config, thermodynamics, gas, aerosol, & 
             &  od_lw, ssa_lw, g_lw, od_sw, ssa_sw, g_sw)
      end if
    else
      g_sw(:,:,istartcol:iendcol) = 0.0_jprb
      if (config%do_lw_aerosol_scattering) then
        ssa_lw(:,:,istartcol:iendcol) = 0.0_jprb
        g_lw(:,:,istartcol:iendcol)   = 0.0_jprb
      end if
    end if

    if (config%do_lw) then
      if (config%iverbose >= 2) then
        write(nulout,'(a)') 'Computing longwave fluxes'
      end if
      !call solver_lw(config, geometry, single_level, od_lw, ssa_lw, g_lw, &
      !     &  od_lw_cloud, ssa_lw_cloud, pf_lw_cloud, &
      !     &  planck_hl, lw_emission, lw_albedo, flux)
    end if

    if (config%do_sw) then
      if (config%iverbose >= 2) then
        write(nulout,'(a)') 'Computing shortwave fluxes'
      end if
      call solver_interface_sw(config, geometry, single_level, & 
           &  od_sw, ssa_sw, g_sw, od_sw_cloud, ssa_sw_cloud, &
           &  pf_sw_cloud, sw_albedo_direct, sw_albedo_diffuse, &
           &  incoming_sw, flux)
    end if
        
    if (lhook) call dr_hook('ecrad3d_interface:ecrad3d',1,hook_handle)
    
  end subroutine ecrad3d
  
end module ecrad3d_interface
