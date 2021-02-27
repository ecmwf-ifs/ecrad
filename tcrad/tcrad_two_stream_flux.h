! tcrad_two_stream_flux.h - Fluxes in multi-region atmospheric profile -*- f90 -*-
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
! This file is included in modules specifying the NREGION parameter
! (typically 2 or 3) which makes this routine use a "doubleclouds" or
! "tripleclouds" assumption.
!

!---------------------------------------------------------------------
! Compute fluxes at top and base of each layer and region from
! precomputed layer transmittance, reflectance and upward/downward
! layer sources, and precomputed upward and downward cloud overlap
! matrices, using the Tripleclouds two-stream method of Shonk and
! Hogan (2008).
subroutine calc_two_stream_flux(nspec, nlev, surf_emission, surf_albedo, &
     &  reflectance, transmittance, source_up, source_dn, &
     &  is_cloud_free_layer, u_overlap, v_overlap, &
     &  flux_up_base, flux_dn_base, flux_up_top, flux_dn_top)

  use parkind1, only           : jpim, jprb
  use yomhook,  only           : lhook, dr_hook
  
  implicit none
  
  ! Inputs

  ! Number of spectral intervals and levels
  integer(jpim), intent(in) :: nspec, nlev

  ! Surface upwards emission, in W m-2 (i.e. emissivity multiplied by
  ! Planck function at the surface skin temperature) integrated across
  ! each spectral interval, and albedo in the same intervals
  real(jprb), intent(in),  dimension(nspec) :: surf_emission, surf_albedo
  
  ! Reflectance and transmittance of each layer and region
  real(jprb), intent(in),  dimension(nspec,NREGION,nlev) :: reflectance
  real(jprb), intent(in),  dimension(nspec,NREGION,nlev) :: transmittance

  ! Rate of emission upward from the top of a layer and downwards from
  ! its base, due to emission within the layer, in Watts of power per
  ! square metre of the entire gridbox, so the energy is scaled by the
  ! size of each region
  real(jprb), intent(in),  dimension(nspec,NREGION,nlev) :: source_up
  real(jprb), intent(in),  dimension(nspec,NREGION,nlev) :: source_dn

  ! Where are the cloud free layers?  Includes dummy layers at 0 and
  ! nlev+1 which are also cloud free
  logical,    intent(in) :: is_cloud_free_layer(0:nlev+1)

  ! U and V overlap matrices, defined by Shonk and Hogan (2008)
  real(jprb), intent(in),  dimension(NREGION,NREGION,nlev+1) :: u_overlap, v_overlap
  
  ! Outputs

  ! Upwelling and downwelling fluxes at the top and base of each layer
  ! in the regions of that layer, in Watts of power per square metre
  ! of the entire gridbox, so the energy is scaled by the size of each
  ! region
  real(jprb), intent(out), dimension(nspec,NREGION,nlev) :: flux_up_base, flux_dn_base
  real(jprb), intent(out), dimension(nspec,NREGION,nlev) :: flux_up_top, flux_dn_top
  
  ! Local variables

  ! Total albedo of the atmosphere/surface just above a layer
  ! interface with respect to downwelling diffuse radiation at that
  ! interface, where level index = 1 corresponds to the
  ! top-of-atmosphere
  real(jprb), dimension(nspec, NREGION, nlev+1) :: total_albedo

  ! Upwelling radiation just above a layer interface due to emission
  ! below that interface, where level index = 1 corresponds to the
  ! top-of-atmosphere
  real(jprb), dimension(nspec, NREGION, nlev+1) :: total_source

  ! Total albedo and source of the atmosphere just below a layer interface
  real(jprb), dimension(nspec, NREGION) :: total_albedo_below, total_source_below

  ! Term in the adding method to replace divisions by multiplications
  real(jprb), dimension(nspec, NREGION) :: inv_denom

  ! Index of highest layer containing cloud
  integer(jpim) :: icloudtop

  ! Loop indices for spectral interval, level and region
  integer(jpim) :: jspec, jlev, jreg, jreg2

  real(jprb) :: hook_handle

  if (lhook) call dr_hook('tcrad:calc_two_stream_flux',0,hook_handle)

  ! --------------------------------------------------------
  ! Section 1: Prepare variables and arrays
  ! --------------------------------------------------------

  ! Find cloud top
  icloudtop = nlev
  do jlev = 1,nlev
    if (.not. is_cloud_free_layer(jlev)) then
      icloudtop = jlev
      exit
    end if
  end do

  ! --------------------------------------------------------
  ! Section 2: Clear-sky downwelling fluxes from TOA to cloud top
  ! --------------------------------------------------------
  ! Initialize outputs
  flux_up_base = 0.0_jprb
  flux_dn_base = 0.0_jprb
  flux_up_top  = 0.0_jprb
  flux_dn_top  = 0.0_jprb

  ! Downwelling fluxes from TOA to cloud top, neglecting clear-air
  ! scattering. Initial flux_dn_top in the clear region is already
  ! zero.
  do jlev = 1,icloudtop
    if (jlev > 1) then
      flux_dn_top(:,1,jlev) = flux_dn_base(:,1,jlev-1)
    end if
    flux_dn_base(:,1,jlev) = source_dn(:,1,jlev) &
         &                 + transmittance(:,1,jlev)*flux_dn_top(:,1,jlev)
  end do

  ! --------------------------------------------------------
  ! Section 3: Compute total sources and albedos up to cloud top
  ! --------------------------------------------------------
  total_albedo = 0.0_jprb
  total_source = 0.0_jprb
  ! Calculate the upwelling radiation emitted by the surface, and copy
  ! the surface albedo into total_albedo
  do jreg = 1,NREGION
    do jspec = 1,nspec
      total_source(jspec,jreg,nlev+1) = u_overlap(jreg,1,nlev+1)*surf_emission(jspec)
      total_albedo(jspec,jreg,nlev+1) = surf_albedo(jspec)
    end do
  end do

  ! Work up from the surface computing the total albedo of the
  ! atmosphere and the total upwelling due to emission below each
  ! level below using the Adding Method
  do jlev = nlev,icloudtop,-1
    total_albedo_below = 0.0_jprb

    if (is_cloud_free_layer(jlev)) then
      total_albedo_below = 0.0_jprb
      total_source_below = 0.0_jprb
      do jspec = 1,nspec
        inv_denom(jspec,1) = 1.0_jprb &
             &  / (1.0_jprb - total_albedo(jspec,1,jlev+1)*reflectance(jspec,1,jlev))
        total_albedo_below(jspec,1) = reflectance(jspec,1,jlev) &
             &  + transmittance(jspec,1,jlev)*transmittance(jspec,1,jlev)*total_albedo(jspec,1,jlev+1) &
             &  * inv_denom(jspec,1)
        total_source_below(jspec,1) = source_up(jspec,1,jlev) &
             &  + transmittance(jspec,1,jlev)*(total_source(jspec,1,jlev+1) &
             &  + total_albedo(jspec,1,jlev+1)*source_dn(jspec,1,jlev)) &
             &  * inv_denom(jspec,1)
      end do
    else
      inv_denom = 1.0_jprb / (1.0_jprb - total_albedo(:,:,jlev+1)*reflectance(:,:,jlev))
      total_albedo_below = reflectance(:,:,jlev) &
           &  + transmittance(:,:,jlev)*transmittance(:,:,jlev)*total_albedo(:,:,jlev+1) &
           &  * inv_denom
      total_source_below = source_up(:,:,jlev) &
           &  + transmittance(:,:,jlev)*(total_source(:,:,jlev+1) &
           &  + total_albedo(:,:,jlev+1)*source_dn(:,:,jlev)) &
           &  * inv_denom
    end if
    
    ! Account for cloud overlap when converting albedo below a layer
    ! interface to the equivalent values just above
    if (is_cloud_free_layer(jlev) .and. is_cloud_free_layer(jlev-1)) then
      total_albedo(:,:,jlev) = total_albedo_below(:,:)
      total_source(:,:,jlev) = total_source_below(:,:)
    else
      total_source(:,:,jlev) = singlemat_x_vec(nspec,&
           &  u_overlap(:,:,jlev), total_source_below)
      ! Use overlap matrix and exclude "anomalous" horizontal
      ! transport described by Shonk & Hogan (2008).  Therefore, the
      ! operation we perform is essentially diag(total_albedo) =
      ! matmul(transpose(v_overlap), diag(total_albedo_below)).
      do jreg = 1,NREGION
        do jreg2 = 1,NREGION
          total_albedo(:,jreg,jlev) &
               &  = total_albedo(:,jreg,jlev) &
               &  + total_albedo_below(:,jreg2) * v_overlap(jreg2,jreg,jlev)
        end do
      end do
      
    end if
        
  end do ! Reverse loop over levels

  ! --------------------------------------------------------
  ! Section 4: Compute fluxes up to top-of-atmosphere
  ! --------------------------------------------------------

  if (icloudtop > 1) then
    ! Compute the fluxes in the layer just above the highest cloud
    flux_up_base(:,1,icloudtop-1) = total_source(:,1,icloudtop) & 
         &  + total_albedo(:,1,icloudtop)*flux_dn_base(:,1,icloudtop-1)
    flux_up_top(:,1,icloudtop-1) = source_up(:,1,icloudtop-1) &
         &  + transmittance(:,1,icloudtop-1)*flux_up_base(:,1,icloudtop-1)
    ! Compute fluxes in remaining layers up to TOA
    do jlev = icloudtop-2,1,-1
      flux_up_base(:,1,jlev) = flux_up_top(:,1,jlev+1)
      flux_up_top(:,1,jlev) = source_up(:,1,jlev) &
           &                + transmittance(:,1,jlev)*flux_up_base(:,1,jlev)
    end do
  end if

  ! --------------------------------------------------------
  ! Section 5: Compute fluxes from cloud top down to surface
  ! --------------------------------------------------------

  ! Copy over downwelling spectral fluxes at top of first scattering
  ! layer, using overlap matrix to translate to the regions of the
  ! first layer of cloud
  if (icloudtop > 1) then
    do jreg = 1,NREGION
      flux_dn_top(:,jreg,icloudtop) &
           &  = v_overlap(jreg,1,icloudtop)*flux_dn_base(:,1,icloudtop-1)
    end do
    ! else the highest layer is cloudy, in which case
    ! flux_dn_top(:,:,jlev=1) is already set to zero
  end if


  ! Final loop back down through the atmosphere to compute fluxes
  do jlev = icloudtop,nlev

    if (is_cloud_free_layer(jlev)) then
      do jspec = 1,nspec
        flux_dn_base(jspec,1,jlev) = (transmittance(jspec,1,jlev)*flux_dn_top(jspec,1,jlev) &
             &  + reflectance(jspec,1,jlev)*total_source(jspec,1,jlev+1) &
             &  + source_dn(jspec,1,jlev) ) &
             &  / (1.0_jprb - reflectance(jspec,1,jlev)*total_albedo(jspec,1,jlev+1))
        flux_up_base(jspec,1,jlev) = total_source(jspec,1,jlev+1) &
             &  + flux_dn_base(jspec,1,jlev)*total_albedo(jspec,1,jlev+1)
        flux_up_top(jspec,1,jlev) = flux_up_base(jspec,1,jlev)*transmittance(jspec,1,jlev) &
             &  + source_up(jspec,1,jlev) + flux_dn_top(jspec,1,jlev)*reflectance(jspec,1,jlev)
      end do
    else
      flux_dn_base(:,:,jlev) = (transmittance(:,:,jlev)*flux_dn_top(:,:,jlev) &
           &     + reflectance(:,:,jlev)*total_source(:,:,jlev+1) + source_dn(:,:,jlev) ) &
           &  / (1.0_jprb - reflectance(:,:,jlev)*total_albedo(:,:,jlev+1))
      flux_up_base(:,:,jlev) = total_source(:,:,jlev+1) &
           &  + flux_dn_base(:,:,jlev)*total_albedo(:,:,jlev+1)
      flux_up_top(:,:,jlev) = flux_up_base(:,:,jlev)*transmittance(:,:,jlev) &
           &  + source_up(:,:,jlev) + flux_dn_top(:,:,jlev)*reflectance(:,:,jlev)
    end if
    
    if (jlev < nlev) then
      if (.not. (is_cloud_free_layer(jlev) .and. is_cloud_free_layer(jlev+1))) then
        ! Account for overlap rules in translating fluxes just above a
        ! layer interface to the values just below
        flux_dn_top(:,:,jlev+1) = singlemat_x_vec(nspec, &
             &  v_overlap(:,:,jlev+1), flux_dn_base(:,:,jlev))
      else 
        ! Two clear-sky layers: copy the fluxes
        flux_dn_top(:,1,jlev+1) = flux_dn_base(:,1,jlev)
      end if
    end if

  end do

  if (lhook) call dr_hook('tcrad:calc_two_stream_flux',1,hook_handle)

end subroutine calc_two_stream_flux


