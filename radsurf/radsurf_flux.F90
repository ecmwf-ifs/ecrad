! radsurf_flux.f90 - Derived type to store fluxes into facets of the surface
!
! Copyright (C) 2017 ECMWF
!
! Author:  Robin Hogan
! Email:   r.j.hogan@ecmwf.int
! License: see the COPYING file for details

module radsurf_flux

  use parkind1, only : jprb

  implicit none

  !---------------------------------------------------------------------
  ! This derived type contains the output from the surface radiation
  ! calculation.  
  type surface_flux_type
     ! If there are tiles then the broadband fluxes are reported at
     ! each facet of the surface, with dimension (ncol,nfacet)
     real(jprb), allocatable, dimension(:,:) :: &
          &  lw_dn_facet, lw_up_facet, &
          &  sw_dn_facet, sw_dn_direct_facet, sw_up_facet
     ! Likewise the absorption by a vegetation canopy or the air in a
     ! street canopy is dimensioned (ncol,ntile)
     real(jprb), allocatable, dimension(:,:) :: &
          &  lw_abs_canopy, sw_abs_canopy
     ! If there are tiles and
     ! config%do_surface_sw_spectral_flux==.true. then the
     ! facet/canopy quantities in the shortwave are output per band,
     ! dimensioned (nband,ncol,nfacet) or (nband,ncol,ntile)
     real(jprb), allocatable, dimension(:,:,:) :: &
          &  sw_dn_facet_band, sw_dn_direct_facet_band, sw_up_facet_band
     real(jprb), allocatable, dimension(:,:,:) :: &
          &  sw_abs_canopy_band

   contains
     procedure :: allocate => allocate_surface_flux_type
     procedure :: deallocate => deallocate_surface_flux_type

  end type surface_flux_type

contains

  !---------------------------------------------------------------------
  ! Allocate arrays for surface fluxes at facets. The arrays are
  ! dimensioned for columns between istartcol and iendcol
  subroutine allocate_surface_flux_type(this, config, istartcol, iendcol, i_tile_representation)

    use yomhook,            only : lhook, dr_hook
!    use radiation_io,       only : nulerr, radiation_abort
    use radiation_config,   only : config_type
    use radsurf_properties, only : NTileFacets

    integer, intent(in)                     :: istartcol, iendcol
    class(surface_flux_type), intent(inout) :: this
    type(config_type), intent(in)           :: config
    integer, intent(in)                     :: i_tile_representation(:)

    ! Surface description
    integer :: ntile, nfacet

    real(jprb)                      :: hook_handle

    if (lhook) call dr_hook('radsurf_flux:allocate',0,hook_handle)

    nfacet    = sum(NTileFacets (i_tile_representation))
    ntile     = size(i_tile_representation)

    if (config%do_lw) then
      allocate(this%lw_dn_facet(istartcol:iendcol,nfacet))
      allocate(this%lw_up_facet(istartcol:iendcol,nfacet))
      allocate(this%lw_abs_canopy(istartcol:iendcol,ntile))
    end if
    if (config%do_sw) then
      allocate(this%sw_dn_facet(istartcol:iendcol,nfacet))
      allocate(this%sw_dn_direct_facet(istartcol:iendcol,nfacet))
      allocate(this%sw_up_facet(istartcol:iendcol,nfacet))
      allocate(this%sw_abs_canopy(istartcol:iendcol,ntile))
    end if

    if (lhook) call dr_hook('radsurf_flux:allocate',1,hook_handle)

  end subroutine allocate_surface_flux_type


  subroutine deallocate_surface_flux_type(this)

    use yomhook,          only : lhook, dr_hook

    class(surface_flux_type), intent(inout) :: this
    real(jprb)                              :: hook_handle

    if (lhook) call dr_hook('radsurf_flux:deallocate',0,hook_handle)

    if (allocated(this%lw_dn_facet)) then
      deallocate(this%lw_dn_facet)
      deallocate(this%lw_up_facet)
      deallocate(this%lw_abs_canopy)
    end if
    if (allocated(this%sw_dn_facet)) then
      deallocate(this%sw_dn_facet)
      deallocate(this%sw_dn_direct_facet)
      deallocate(this%sw_up_facet)
      deallocate(this%sw_abs_canopy)
    end if

    if (lhook) call dr_hook('radsurf_flux:deallocate',1,hook_handle)

  end subroutine deallocate_surface_flux_type

end module radsurf_flux
