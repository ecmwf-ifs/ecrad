module ecrad3d_solver_sw

  !procedure :: solver_sw
  
contains

  subroutine solver_sw(config, geometry, ncol, nlay, nmaxspec, nspec, &
       &               cos_sza, solar_azim, incoming, &
       &               surf_albedo_dir, surf_albedo_diff, &
       &               od, ssa, asymmetry, &
       &               flux_dn_dir_top, flux_dn_diff_top, flux_up_top, &
       &               flux_dn_dir_base, flux_dn_diff_base, flux_up_base, &
       &               flux_div, &
       &               first_direct_3d_layer)

    use parkind1,                only : jprb
    use ecrad3d_config,          only : config_type
    use ecrad3d_geometry,        only : geometry_type
    use ecrad3d_layer_solutions, only : calc_ref_trans_sw

    implicit none
    
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
    real(jprb), intent(in), dimension(ncol) :: cos_sza, solar_azim

    ! Incoming solar flux (W m-2) per spectral interval into a plane
    ! perpendicular to the sun's rays
    real(jprb), intent(in), dimension(nspec) :: incoming

    ! Surface albedo to direct and diffuse radiation
    real(jprb), intent(in), dimension(ncol,nspec) &
         &  :: surf_albedo_dir, surf_albedo_diff
    
    ! Layer optical depth, single scattering albedo and asymmetry
    ! factor
    real(jprb), intent(in), dimension(ncol,nmaxspec,nlay) &
         &  :: od, ssa, asymmetry

    ! OUTPUTS

    ! Up and down fluxes (the latter split into direct and diffuse) at
    ! the top of each layer, plus the surface value at the end. Note
    ! that the downward direct flux is into a plane perpendicular to
    ! the solar direction.
    real(jprb), intent(out), dimension(ncol,nmaxspec,nlay+1) &
         &  :: flux_dn_dir_top, flux_dn_diff_top, flux_up_top

    ! OPTIONAL OUTPUTS
    
    ! The up and down fluxes at the base of each layer, excluding the
    ! top-of-atmosphere value.
    real(jprb), intent(out), dimension(ncol,nmaxspec,nlay), optional &
         &  :: flux_dn_dir_base, flux_dn_diff_base, flux_up_base

    ! Flux divergence, intended to be used to compute heating rates
    real(jprb), intent(out), dimension(ncol,nmaxspec,nlay), optional :: flux_div

    ! OPTIONAL INPUTS
     
    ! Index of first layer at which "advection" of the direct beam is
    ! to be performed, typically the first layer containing cloud
    integer, intent(in), optional :: first_direct_3d_layer
    
    ! LOCALS

    ! Fluxes at base of one layer, since not all callers of this
    ! routine need them to be output
    real(jprb), dimension(ncol,nspec) &
         &  :: flux_dn_dir_base_layer, flux_dn_diff_base_layer, flux_up_base_layer

    ! Reflectance and transmittance to direct and diffuse radiation,
    ! where trans_dir_diff and trans_dir_dir are the fractions of
    ! incident direct radiation at the top of a layer that emerge from
    ! the base as diffuse or direct radiation, respectively.
    real(jprb), dimension(ncol,nspec,nlay) &
         &  :: ref_diff, trans_diff, ref_dir, trans_dir_diff, trans_dir_dir

    ! Albedo of entire earth-atmosphere below a half-level, and the
    ! upward radiation due to scattering of the direct beam below a
    ! half-level, both defined at layer base
    real(jprb), dimension(ncol,nspec,nlay) :: albedo_base, source_base

    ! As above, but at layer top. Note that these are not required in
    ! the second pass of the adding method so do not need to be stored
    ! at every layer.
    real(jprb), dimension(ncol,nspec) :: albedo_top, source_top, inv_denom
    
    ! Local version of optional first_direct_3d_layer
    integer :: local_first_direct_3d_layer
    
    ! Loop indices    
    integer :: jc, jl, js

    if (present(first_direct_3d_layer)) then
      local_first_direct_3d_layer = first_direct_3d_layer
    else
      local_first_direct_3d_layer = 1
    end if

    ! Compute two-stream layerwise reflectances and transmittances
    call calc_ref_trans_sw(ncol, nlay, nmaxspec, nspec, cos_sza, &
         &                 od, ssa, asymmetry, ref_diff, trans_diff, &
         &                 ref_dir, trans_dir_diff, trans_dir_dir)
    
    ! DIRECT BEAM

    ! Incoming radiation at top of atmosphere
    do js = 1,nspec
      flux_dn_dir_top(:,js,1) = incoming(js)
    end do
    
    ! Loop down through the atmosphere 
    do jl = 1,nlay

      ! Direct transmittance through a layer
      !$OMP PARALLEL DO PRIVATE(jc)
      do js = 1,nspec
        do jc = 1,ncol
          flux_dn_dir_base_layer(jc,js) &
               &  = flux_dn_dir_top(jc,js,jl) * trans_dir_dir(jc,js,jl)
          ! * exp(-od(jc,js,jl)/cos_sza(jc))
        end do

        ! Optionally store flux at base of layer
        if (present(flux_dn_dir_base)) then
          flux_dn_dir_base(jc,js,jl+1) = flux_dn_dir_base_layer(jc,js)
        end if

        ! Optionally store flux divergence (to be augmented by diffuse
        ! fluxes later)
        if (present(flux_div)) then
          do jc = 1,ncol
            flux_div(jc,js,jl) = cos_sza(jc) * (flux_dn_dir_base_layer(jc,js) &
                 &                                   - flux_dn_dir_top(jc,js,jl))
          end do
        end if
      end do
      !$OMP END PARALLEL DO

      ! Horizontal "advection"
      if (jl < nlay) then
        if (config%do_3d .and. jl >= local_first_direct_3d_layer) then
          call geometry%advect(ncol, nspec, jl,flux_dn_dir_base_layer,flux_dn_dir_top(:,:,jl+1))
        else
          do js = 1,nspec
            flux_dn_dir_top(:,js,jl+1) = flux_dn_dir_base_layer(:,js)
          end do
        end if
      end if

    end do

    ! Store surface direct flux
    flux_dn_dir_top(:,:,nlay+1) = flux_dn_dir_base_layer
    
    ! Assign values at the surface
    do js = 1,nspec
      do jc = 1,ncol
        source_base(jc,js,nlay) = flux_dn_dir_base(jc,js,nlay) &
             &  * surf_albedo_dir(jc,js) * cos_sza(jc)
        albedo_base(jc,js,nlay) = surf_albedo_diff(jc,js)
      end do
    end do

    ! DIFFUSE FLUXES

    ! Loop up through the atmosphere
    do jl = nlay,1,-1
      ! First step of adding method
      !$OMP PARALLEL DO PRIVATE(jc)
      do js = 1,nspec
        do jc = 1,ncol
          ! Lacis and Hansen (1974) Eq 33, Shonk & Hogan (2008) Eq 10:
          inv_denom(jc,js) = 1.0_jprb / (1.0_jprb-albedo_base(jc,js,jl)*ref_diff(jc,js,jl))
          ! Shonk & Hogan (2008) Eq 9, Petty (2006) Eq 13.81:
          albedo_top(jc,js) = ref_diff(jc,js,jl) + trans_diff(jc,js,jl) &
               &  * trans_diff(jc,js,jl) * albedo_base(jc,js,jl) * inv_denom(jc,js)
          ! Shonk & Hogan (2008) Eq 11:
          source_top(jc,js) = ref_dir(jc,js,jl) * flux_dn_dir_top(jc,js,jl) &
               &  + trans_diff(jc,js,jl) * (source_base(jc,js,jl) &
               &  + albedo_base(jc,js,jl)*trans_dir_diff(jc,js,jl)*flux_dn_dir_top(jc,js,jl)) &
               &  * inv_denom(jc,js)
        end do
      end do
      !$OMP END PARALLEL DO

      ! Horizontal "diffusion"
      if (jl > 1) then
        if (config%do_3d) then
          call geometry%diffuse(ncol, nspec, jl, albedo_top, albedo_base(:,:,jl-1))
          call geometry%diffuse(ncol, nspec, jl, source_top, source_base(:,:,jl-1))
        else
          do js = 1,nspec
            albedo_base(:,js,jl-1) = albedo_top(:,js)
            source_base(:,js,jl-1) = source_top(:,js)
          end do
        end if
      end if
    end do

    do js = 1,nspec
      ! At top-of-atmosphere there is no diffuse downwelling radiation
      flux_dn_diff_top(:,js,1) = 0.0_jprb
      ! ...and all upwelling radiation due to scattering by direct beam
      ! below that level
      flux_up_top(:,js,1)      = source_top(:,js)
    end do
      
    ! Loop down through the atmosphere
    do jl = 1,nlay
      ! Second step of adding method
      !$OMP PARALLEL DO PRIVATE(jc)
      do js = 1,nspec
        do jc = 1,ncol
          ! Lacis and Hansen (1974) Eq 33, Shonk & Hogan (2008) Eq 10:
          inv_denom(jc,js) = 1.0_jprb / (1.0_jprb-albedo_base(jc,js,jl)*ref_diff(jc,js,jl))
          ! Shonk & Hogan (2008) Eq 14 (after simplification):
          flux_dn_diff_base_layer(jc,js) &
               &  = (trans_diff(jc,js,jl)*flux_dn_diff_top(jc,js,jl) &
               &     + ref_diff(jc,js,jl)*source_base(jc,js,jl) &
               &     + trans_dir_diff(jc,js,jl)*flux_dn_dir_top(jc,js,jl)) &
               &   * inv_denom(jc,js)
          ! Shonk & Hogan (2008) Eq 12:
          flux_up_base_layer(jc,js) = albedo_base(jc,js,jl)*flux_dn_diff_base_layer(jc,js) &
               &                    + source_base(jc,js,jl)
          flux_up_top(jc,js,jl) = flux_up_base_layer(jc,js)  * trans_diff(jc,js,jl) &
               &                + flux_dn_diff_top(jc,js,jl) * ref_diff(jc,js,jl) &
               &                + flux_dn_dir_top(jc,js,jl)  * ref_dir(jc,js,jl)
        end do

        ! Optionally store additional outputs
        if (present(flux_up_base)) then
          flux_up_base(:,js,jl) = flux_up_base_layer(:,js)
        end if

        if (present(flux_dn_diff_base)) then
          flux_dn_diff_base(:,js,jl) = flux_dn_diff_base_layer(:,js)
        end if

        if (present(flux_div)) then
          do jc = 1,ncol
            flux_div(jc,js,jl) = flux_div(jc,js,jl) &
                 &  - flux_dn_diff_top(jc,js,jl) + flux_dn_diff_base_layer(jc,js) &
                 &  + flux_up_top(jc,js,jl)      - flux_up_base_layer(jc,js)
          end do
        end if
        
      end do
      !$OMP END PARALLEL DO
      
      ! Horizontal "diffusion"
      if (jl < nlay) then
        if (config%do_3d) then
          call geometry%diffuse(ncol, nspec, jl, flux_dn_diff_base, flux_dn_diff_top(:,:,jl+1))
        else
          do js = 1,nspec
            flux_dn_diff_top(:,js,jl+1) = flux_dn_diff_base_layer(:,js)
          end do
        end if
      end if
    end do

  end subroutine solver_sw
    
end module ecrad3d_solver_sw
