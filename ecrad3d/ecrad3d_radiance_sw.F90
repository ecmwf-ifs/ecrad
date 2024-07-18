module ecrad3d_radiance_sw

contains

  subroutine radiance_sw(config, geometry, ncol, nlay, nmaxspec, nspec, &
       &                 cos_sza, solar_azimuth, cos_zenith_angle, azimuth_angle, &
       &                 od, ssa, asymmetry, &
       &                 flux_dn_dir_top, flux_dn_diff_top, flux_up_diff_base, &
       &                 radiance, first_direct_3d_layer)
    
    use parkind1,                only : jprb
    use ecrad3d_config,          only : config_type
    use ecrad3d_geometry,        only : geometry_type
    use ecrad3d_layer_solutions, only : calc_radiance_coeffs_sw

    implicit none

    real(jprb), parameter :: ONE_OVER_PI  = 1.0_jprb / acos(-1.0_jprb)
    
    ! INPUTS
    
    type(config_type),   intent(in) :: config
    type(geometry_type), intent(in) :: geometry

    ! Number of columns, layers and spectral intervals defining the
    ! array dimensions
    integer, intent(in) :: ncol, nlay, nmaxspec


    ! Number of spectral intervals to process, allowing for the
    ! calling routine to process the spectrum in blocks and the final
    ! block may have fewer intervals
    integer, intent(in) :: nspec

    ! Cosine of solar zenith angle and solar azimuth angle (radians)
    real(jprb), intent(in), dimension(ncol) :: cos_sza, solar_azimuth
    
    ! Cosine of sensor zenith angle and sensor azimuth angle (radians)
    real(jprb), intent(in), dimension(ncol) :: cos_zenith_angle, azimuth_angle
    
    ! Layer optical depth, single scattering albedo and asymmetry
    ! factor
    real(jprb), intent(in), dimension(ncol,nmaxspec,nlay) &
         &  :: od, ssa, asymmetry

    ! Downward fluxes, direct and diffuse, at the top of each layer,
    ! plus the surface value at the end. Note that the downward direct
    ! flux is into a plane perpendicular to the solar direction.
    real(jprb), intent(in), dimension(ncol,nmaxspec,nlay+1) &
         &  :: flux_dn_dir_top, flux_dn_diff_top

    ! Flux up at the base of each layer, excluding upwelling at TOA
    real(jprb), intent(in), dimension(ncol,nmaxspec,nlay) :: flux_up_diff_base

    ! OPTIONAL INPUTS
     
    ! Index of first layer at which "advection" of the direct beam is
    ! to be performed, typically the first layer containing cloud
    integer, intent(in), optional :: first_direct_3d_layer
    
    ! OUTPUT
    
    ! Upward radiances
    real(jprb), intent(out), dimension(ncol,nmaxspec) :: radiance

    ! LOCALS

    real(jprb), dimension(:,:), allocatable &
         &  :: dir_dn_to_rad, diff_dn_to_rad, diff_up_to_rad, trans_rad, radiance_tmp

    real(jprb) :: azim_diff(ncol)

    integer :: local_first_direct_3d_layer
    
    integer :: jl, js

    ! Check if sensors are all at nadir, avoiding the requirement for
    ! 3D transport
    logical :: all_nadir
    
    if (present(first_direct_3d_layer)) then
      local_first_direct_3d_layer = first_direct_3d_layer
    else
      local_first_direct_3d_layer = 1
    end if

    if (any(cos_zenith_angle < 0.999_jprb)) then
      all_nadir = .false.
    else
      all_nadir = .true.
    end if
    
    ! Surface
    radiance(:,:) = ONE_OVER_PI * flux_up_diff_base(:,:,nlay)

    allocate(dir_dn_to_rad(ncol,nspec))
    allocate(diff_dn_to_rad(ncol,nspec))
    allocate(diff_up_to_rad(ncol,nspec))
    allocate(trans_rad(ncol,nspec))
    allocate(radiance_tmp(ncol,nspec))

    azim_diff(:) = solar_azimuth - azimuth_angle
    
    do jl = nlay,1,-1

      call calc_radiance_coeffs_sw(ncol, nmaxspec, nspec, cos_sza, cos_zenith_angle, azim_diff, &
           &                       od(:,:,jl), ssa(:,:,jl), asymmetry(:,:,jl), &
           &                       dir_dn_to_rad, diff_dn_to_rad, diff_up_to_rad, trans_rad)

      radiance_tmp(:,:) = radiance(:,:)            * trans_rad(:,:) &
           &            + flux_dn_dir_top(:,:,jl)  * dir_dn_to_rad(:,:) &
           &            + flux_dn_diff_top(:,:,jl) * diff_dn_to_rad(:,:) &
           &            + flux_up_diff_base(:,:,jl)* diff_up_to_rad(:,:)

      !write(*,*) 'Layer ', jl, ' dir_dn_to_rad  =', dir_dn_to_rad(1,:)
      !write(*,*) 'Layer ', jl, ' diff_dn_to_rad =', diff_dn_to_rad(1,:)
      !write(*,*) 'Layer ', jl, ' diff_up_to_rad =', diff_up_to_rad(1,:)
      !write(*,*) 'Layer ', jl, ' trans_rad      =', trans_rad(1,:)

      !write(*,*) 'Layer ', jl, ' radiance = ', radiance_tmp(1,:)
      if (config%do_3d .and. .not. all_nadir .and. jl >= local_first_direct_3d_layer) then
        !print *,'Advecting'
        call geometry%advect(ncol, nspec, jl, cos_zenith_angle, azimuth_angle, &
             &               radiance_tmp, radiance)
      else
        do js = 1,nspec
          radiance(:,js) = radiance_tmp(:,js)
        end do
      end if
    end do

    deallocate(dir_dn_to_rad, diff_dn_to_rad, diff_up_to_rad, trans_rad, radiance_tmp)

  end subroutine radiance_sw
    
end module ecrad3d_radiance_sw
