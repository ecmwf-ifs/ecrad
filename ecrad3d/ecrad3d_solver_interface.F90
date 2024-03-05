module ecrad3d_solver_interface

  !procedure :: solver_sw
  
contains

  subroutine solver_interface_sw(config, geometry, single_level, & 
           &  od_clear, ssa_clear, g_clear, od_cloud, ssa_cloud, &
           &  g_cloud, sw_albedo_direct, sw_albedo_diffuse, &
           &  incoming_sw, flux)
    
    use parkind1,                 only : jprb
    use yomhook,                  only : lhook, dr_hook, jphook

    use radiation_io,             only : nulout, nulerr, radiation_abort
    use radiation_config,         only : config_type
    use ecrad3d_config,           only : config3d_type => config_type
    use ecrad3d_geometry,         only : geometry_type
    use radiation_single_level,   only : single_level_type
    use radiation_flux,           only : flux_type
    use ecrad3d_solver_sw,        only : solver_sw
    
    implicit none
    
    type(config_type),        intent(in)    :: config
    type(geometry_type),      intent(in)    :: geometry
    type(single_level_type),  intent(in)    :: single_level
    type(flux_type),          intent(inout) :: flux
    
    ! Layer optical depth, single scattering albedo and asymmetry factor of
    ! gases and aerosols at each shortwave g-point
    real(jprb), intent(in), dimension(config%n_g_sw,geometry%nz,geometry%ncol) &
         &  :: od_clear, ssa_clear, g_clear

    ! Layer in-cloud optical depth and single scattering albedo of
    ! hydrometeors in each shortwave band
    real(jprb), intent(in), dimension(config%n_bands_sw,geometry%nz,geometry%ncol) &
         &  :: od_cloud, ssa_cloud, g_cloud

    ! Direct and diffuse shortwave surface albedo in each shortwave
    ! g-point; note that these are weighted averages of the values
    ! from individual tiles
    real(jprb), intent(in), dimension(config%n_g_sw,geometry%ncol) :: sw_albedo_direct
    real(jprb), intent(in), dimension(config%n_g_sw,geometry%ncol) :: sw_albedo_diffuse

    ! The incoming shortwave flux into a plane perpendicular to the
    ! incoming radiation at top-of-atmosphere in each of the shortwave
    ! g-points
    real(jprb), intent(in), dimension(config%n_g_sw,geometry%ncol) :: incoming_sw

    ! Local variables

    type(config3d_type) :: config3d
    
    ! Layer optical depth, single scattering albedo and asymmetry
    ! factor
    real(jprb), dimension(geometry%ncol,config%n_g_sw,geometry%nz) &
         &  :: od, ssa, asymmetry

    ! Up and down fluxes (the latter split into direct and diffuse) at
    ! the top of each layer, plus the surface value at the end. Note
    ! that the downward direct flux is into a plane perpendicular to
    ! the solar direction.
    real(jprb), dimension(geometry%ncol,config%n_g_sw,geometry%nz+1) &
         &  :: flux_dn_dir_top, flux_dn_diff_top, flux_up_top

    real(jprb) :: scat_od, scat_od_cloud
    
    integer :: jcol, jlev, jg

    integer :: iband
    
    real(jphook) :: hook_handle

    if (lhook) call dr_hook('ecrad3d_solver_interface:solver_interface_sw',0,hook_handle)

    config3d%do_3d = .true.

    ! Clear-sky calculation

    if (config%do_clear) then
      ! Permute arrays
      do jcol = 1,geometry%ncol
        do jlev = 1,geometry%nz
          do jg = 1,config%n_g_sw
            od(jcol,jg,jlev)        = od_clear(jg,jlev,jcol)
            ssa(jcol,jg,jlev)       = ssa_clear(jg,jlev,jcol)
            asymmetry(jcol,jg,jlev) = g_clear(jg,jlev,jcol)
          end do
        end do
      end do
      ! Run solver on clear-sky optical properties
      call solver_sw(config3d, geometry, geometry%ncol, geometry%nz, config%n_g_sw, &
           &  config%n_g_sw, single_level%cos_sza, single_level%solar_azimuth_angle, &
           &  incoming_sw(:,1), transpose(sw_albedo_direct), transpose(sw_albedo_diffuse), &
           &  od, ssa, asymmetry, flux_dn_dir_top, flux_dn_diff_top, flux_up_top)
      flux%sw_up_clear        = sum(flux_up_top,2)
      flux%sw_dn_direct_clear = sum(flux_dn_dir_top,2)
      flux%sw_dn_clear        = sum(flux_dn_diff_top,2) + flux%sw_dn_direct_clear
    end if

    ! All-sky calculation
    
    ! Combine gas/aerosol and cloud optical properties
    do jcol = 1,geometry%ncol
      do jlev = 1,geometry%nz
        if (od_cloud(1,jlev,jcol) > 0.0_jprb) then
          ! Cloudy layer
          do jg = 1,config%n_g_sw
            iband = config%i_band_from_reordered_g_sw(jg)
            scat_od = od_clear(jg,jlev,jcol)*ssa_clear(jg,jlev,jcol)
            scat_od_cloud = od_cloud(iband,jlev,jcol) &
                 &  * ssa_cloud(iband,jlev,jcol)
            ! Add scaled cloud optical depth to clear-sky value
            od(jcol,jg,jlev) = od_clear(jg,jlev,jcol) &
                 &  + od_cloud(iband,jlev,jcol)
            ! Compute single-scattering albedo and asymmetry
            ! factor of gas-cloud combination
            ssa(jcol,jg,jlev) = (scat_od + scat_od_cloud) &
                 &  / od(jcol,jg,jlev)
            asymmetry(jcol,jg,jlev) = (scat_od*g_clear(jg,jlev,jcol) &
                 &         + scat_od_cloud * g_cloud(iband,jlev,jcol)) &
                 &      / (scat_od + scat_od_cloud)
          end do
        else if (.not. config%do_clear) then
          ! Clear layer
          do jg = 1,config%n_g_sw
            od(jcol,jg,jlev)        = od_clear(jg,jlev,jcol)
            ssa(jcol,jg,jlev)       = ssa_clear(jg,jlev,jcol)
            asymmetry(jcol,jg,jlev) = g_clear(jg,jlev,jcol)
          end do
        end if
      end do
      ! Cloud cover is one or zero
      if (allocated(flux%cloud_cover_sw)) then
        flux%cloud_cover_sw(jcol) = merge(1.0_jprb, 0.0_jprb, any(od(jcol,1,:)>0.0_jprb))
      end if
    end do

    call solver_sw(config3d, geometry, geometry%ncol, geometry%nz, config%n_g_sw, &
         &  config%n_g_sw, single_level%cos_sza, single_level%solar_azimuth_angle, &
         &  incoming_sw(:,1), transpose(sw_albedo_direct), transpose(sw_albedo_diffuse), &
         &  od, ssa, asymmetry, flux_dn_dir_top, flux_dn_diff_top, flux_up_top)

    flux%sw_up        = sum(flux_up_top,2)
    flux%sw_dn_direct = sum(flux_dn_dir_top,2)
    flux%sw_dn        = sum(flux_dn_diff_top,2) + flux%sw_dn_direct
    
    if (lhook) call dr_hook('ecrad3d_solver_interface:solver_interface_sw',1,hook_handle)

  end subroutine solver_interface_sw
  
end module ecrad3d_solver_interface
