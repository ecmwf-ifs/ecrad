! radsurf_intermediate.f90 - Derived type for intermediate radiative properties of surface
!
! Copyright (C) 2017-2018 ECMWF
!
! Author:  Robin Hogan
! Email:   r.j.hogan@ecmwf.int
! License: see the COPYING file for details
!

module radsurf_intermediate

  use parkind1, only : jpim, jprb

  implicit none

  !---------------------------------------------------------------------
  ! Derived type storing a description of radiative properties at the
  ! level of individual facets
  type surface_intermediate_type

    ! Surface "facet" properties, dimensioned (ng, nfacet, istartcol:iendcol)

    ! Longwave blackbody emission (W m-2) and emissivity from facets
    real(kind=jprb), allocatable, dimension(:,:,:) :: planck_facet
    real(kind=jprb), allocatable, dimension(:,:,:) :: lw_emissivity
    ! Shortwave direct and diffuse albedo
    real(kind=jprb), allocatable, dimension(:,:,:) :: sw_albedo_direct
    real(kind=jprb), allocatable, dimension(:,:,:) :: sw_albedo_diffuse

    ! Longwave blackbody emission (W m-2) from regions (volumes
    ! e.g. vegetation and urban canopies)
    real(kind=jprb), allocatable, dimension(:,:,:) :: planck_region
  
    ! Volumetric "region" properties of the canopy, e.g. vegetation or
    ! the space between buildings, dimensioned
    ! (ng,nregion,istartcol:iendcol)

    ! Shortwave reflectance and transmittance to diffuse radiation
    real(kind=jprb), allocatable, dimension(:,:,:) :: sw_ref_dif
    real(kind=jprb), allocatable, dimension(:,:,:) :: sw_tra_dif
    ! Shortwave reflectance and transmittance to direct radiation,
    ! where the transmittance has separate direct-to-diffuse
    ! transmittance (including scattering) and direct-to-direct
    ! transmittance (without scattering)
    real(kind=jprb), allocatable, dimension(:,:,:) :: sw_ref_dir
    real(kind=jprb), allocatable, dimension(:,:,:) :: sw_tra_dir_dif
    real(kind=jprb), allocatable, dimension(:,:,:) :: sw_tra_dir_dir
    ! Fraction of direct radiation at canyon top that is absorbed by
    ! the wall and the atmosphere
    real(kind=jprb), allocatable, dimension(:,:,:) :: sw_wall_abs_dir
    real(kind=jprb), allocatable, dimension(:,:,:) :: sw_air_abs_dir
    ! Ratio of diffuse absorption in a street canyon by the wall
    ! (rather than the air or vegetation in the canyon)
    real(kind=jprb), allocatable, dimension(:,:,:) :: sw_wall_abs_frac_dif
    real(kind=jprb), allocatable, dimension(:,:,:) :: lw_wall_abs_frac
    ! Shortwave direct and diffuse albedo at the top of a region
    real(kind=jprb), allocatable, dimension(:,:,:) :: sw_albedo_direct_reg
    real(kind=jprb), allocatable, dimension(:,:,:) :: sw_albedo_diffuse_reg
    ! Longwave reflectance, transmittance and source (note that source
    ! is the same up and down because the canopy temperature is
    ! assumed constant with height)
    real(kind=jprb), allocatable, dimension(:,:,:) :: lw_reflectance
    real(kind=jprb), allocatable, dimension(:,:,:) :: lw_transmittance
    real(kind=jprb), allocatable, dimension(:,:,:) :: lw_source
    ! Longwave emission and emissivity passed to the atmosphere scheme
    ! from the top of a region
    real(kind=jprb), allocatable, dimension(:,:,:) :: lw_emissivity_region
    real(kind=jprb), allocatable, dimension(:,:,:) :: lw_emission_region
    ! Total emission from wall and the air in an urban canopy, needed
    ! at the final partitioning stage to determine net fluxes into
    ! wall and air
    real(kind=jprb), allocatable, dimension(:,:,:) :: lw_total_wall_emission
    real(kind=jprb), allocatable, dimension(:,:,:) :: lw_total_canopy_emission

    ! Number of bands in which calculations are to be performed (can
    ! match either the spectral resolution of the atmosphere or that
    ! of the surface data)
    integer :: nswbands, nlwbands

    ! Do we represent gas radiative transfer in street/vegetation
    ! canopies?
    !logical :: do_canopy_gases_sw = .false.
    !logical :: do_canopy_gases_lw = .false.

    ! Do we use the full spectrum?
    !logical :: use_full_spectrum_sw = .false.
    !logical :: use_full_spectrum_lw = .true.

  contains

    procedure :: allocate    => allocate_surface_intermediate
    procedure :: deallocate  => deallocate_surface_intermediate
    procedure :: calc_boundary_conditions_sw
    procedure :: calc_boundary_conditions_lw
    procedure :: calc_boundary_conditions
    procedure :: partition_fluxes

  end type surface_intermediate_type

contains

  !---------------------------------------------------------------------
  subroutine allocate_surface_intermediate(this, istartcol, iendcol, &
       &                                   config, surface)
    
    use radiation_config,   only : config_type
    use radsurf_properties, only : surface_type

    class(surface_intermediate_type), intent(inout) :: this
    integer(kind=jpim),               intent(in)    :: istartcol, iendcol
    type(config_type),                intent(in)    :: config
    type(surface_type),               intent(in)    :: surface

    call this%deallocate

    !this%use_full_spectrum_sw = config%use_canopy_full_spectrum_sw
    !this%use_full_spectrum_lw = config%use_canopy_full_spectrum_lw
    ! Assume that canopy gases will only be used if we do surface
    ! calculations at the full spectral resolution
    !this%do_canopy_gases_sw   = config%do_canopy_gases_sw
    !this%do_canopy_gases_lw   = config%do_canopy_gases_lw

    if (config%use_canopy_full_spectrum_sw) then
      ! Calculations will be performed at the same spectral resolution
      ! as in the atmosphere
      this%nswbands = config%n_g_sw
    else
      ! Calculations will be performed at the same spectral resolution
      ! as the input surface-property data
      this%nswbands = surface%nalbedobands
    end if

    if (config%use_canopy_full_spectrum_lw) then
      ! Calculations will be performed at the same spectral resolution
      ! as in the atmosphere
      this%nlwbands = config%n_g_lw
    else
      ! Calculations will be performed at the same spectral resolution
      ! as the input surface-property data
      this%nlwbands = surface%nemissbands
    end if

    if (config%do_lw) then
      allocate(this%planck_facet     (this%nlwbands,surface%nfacet,istartcol:iendcol))
      allocate(this%lw_emissivity    (this%nlwbands,surface%nfacet,istartcol:iendcol))
    end if
    if (config%do_sw) then
      allocate(this%sw_albedo_direct (this%nswbands,surface%nfacet,istartcol:iendcol))
      allocate(this%sw_albedo_diffuse(this%nswbands,surface%nfacet,istartcol:iendcol))
    end if

    if (surface%nregion > 0) then
      if (config%do_lw) then
        allocate(this%planck_region   (this%nlwbands,surface%nregion,istartcol:iendcol))
        allocate(this%lw_reflectance  (this%nlwbands,surface%nregion,istartcol:iendcol))
        allocate(this%lw_transmittance(this%nlwbands,surface%nregion,istartcol:iendcol))
        allocate(this%lw_source       (this%nlwbands,surface%nregion,istartcol:iendcol))
        allocate(this%lw_emission_region(this%nlwbands,surface%nregion,istartcol:iendcol))
        allocate(this%lw_emissivity_region(this%nlwbands,surface%nregion,istartcol:iendcol))
        allocate(this%lw_wall_abs_frac(this%nlwbands,surface%nregion,istartcol:iendcol))
        allocate(this%lw_total_wall_emission(this%nlwbands,surface%nregion,istartcol:iendcol))
        allocate(this%lw_total_canopy_emission(this%nlwbands,surface%nregion,istartcol:iendcol))
      end if
      if (config%do_sw) then
        allocate(this%sw_ref_dif    (this%nswbands,surface%nregion,istartcol:iendcol))
        allocate(this%sw_tra_dif    (this%nswbands,surface%nregion,istartcol:iendcol))
        allocate(this%sw_ref_dir    (this%nswbands,surface%nregion,istartcol:iendcol))
        allocate(this%sw_tra_dir_dif(this%nswbands,surface%nregion,istartcol:iendcol))
        allocate(this%sw_tra_dir_dir(this%nswbands,surface%nregion,istartcol:iendcol))
        allocate(this%sw_wall_abs_dir(this%nswbands,surface%nregion,istartcol:iendcol))
        allocate(this%sw_air_abs_dir(this%nswbands,surface%nregion,istartcol:iendcol))
        allocate(this%sw_wall_abs_frac_dif(this%nswbands,surface%nregion,istartcol:iendcol))
        allocate(this%sw_albedo_diffuse_reg(this%nswbands,surface%nregion,istartcol:iendcol))
        allocate(this%sw_albedo_direct_reg(this%nswbands,surface%nregion,istartcol:iendcol))
      end if
    end if

  end subroutine allocate_surface_intermediate

  !---------------------------------------------------------------------
  subroutine deallocate_surface_intermediate(this)

    class(surface_intermediate_type), intent(inout) :: this

    if (allocated(this%planck_facet)) then
      deallocate(this%planck_facet)
    end if
    if (allocated(this%planck_region)) then
      deallocate(this%planck_region)
    end if
    if (allocated(this%lw_emissivity)) then
      deallocate(this%lw_emissivity)
    end if
    if (allocated(this%sw_albedo_direct)) then
      deallocate(this%sw_albedo_direct)
    end if
    if (allocated(this%sw_albedo_diffuse)) then
      deallocate(this%sw_albedo_diffuse)
    end if

    if (allocated(this%lw_reflectance)) then
      deallocate(this%lw_reflectance)
    end if
    if (allocated(this%lw_transmittance)) then
      deallocate(this%lw_transmittance)
    end if
    if (allocated(this%lw_source)) then
      deallocate(this%lw_source)
    end if
    if (allocated(this%lw_emissivity_region)) then
      deallocate(this%lw_emissivity_region)
    end if
    if (allocated(this%lw_emission_region)) then
      deallocate(this%lw_emission_region)
    end if
    if (allocated(this%lw_wall_abs_frac)) then
      deallocate(this%lw_wall_abs_frac)
    end if
    if (allocated(this%lw_total_wall_emission)) then
      deallocate(this%lw_total_wall_emission)
    end if
    if (allocated(this%lw_total_canopy_emission)) then
      deallocate(this%lw_total_canopy_emission)
    end if

    if (allocated(this%sw_ref_dif)) then
      deallocate(this%sw_ref_dif)
    end if
    if (allocated(this%sw_tra_dif)) then
      deallocate(this%sw_tra_dif)
    end if
    if (allocated(this%sw_ref_dir)) then
      deallocate(this%sw_ref_dir)
    end if
    if (allocated(this%sw_tra_dir_dif)) then
      deallocate(this%sw_tra_dir_dif)
    end if
    if (allocated(this%sw_tra_dir_dir)) then
      deallocate(this%sw_tra_dir_dir)
    end if
    if (allocated(this%sw_wall_abs_frac_dif)) then
      deallocate(this%sw_wall_abs_frac_dif)
    end if
    if (allocated(this%sw_wall_abs_dir)) then
      deallocate(this%sw_wall_abs_dir)
    end if
    if (allocated(this%sw_air_abs_dir)) then
      deallocate(this%sw_air_abs_dir)
    end if
    if (allocated(this%sw_albedo_direct_reg)) then
      deallocate(this%sw_albedo_direct_reg)
    end if
    if (allocated(this%sw_albedo_diffuse_reg)) then
      deallocate(this%sw_albedo_diffuse_reg)
    end if

  end subroutine deallocate_surface_intermediate



  !---------------------------------------------------------------------
  ! Use the detailed physical properties of the surface in "surface",
  ! and optionally the atmospheric optical properties of the lowest
  ! atmospheric level, to work out the shortwave direct/diffuse albedo
  ! that is presented to the rest of the radiation scheme, and stored
  ! in "single_level". Also, store the necessary information so that
  ! after the atmospheric radiation scheme has been run, the net
  ! fluxes at each.
  subroutine calc_boundary_conditions_sw(this, istartcol, iendcol, &   ! in
       &  config, surface, &                    ! in
       &  single_level, &                                 ! out
       &  ext_sw_air, ssa_sw_air, g_sw_air)               ! in, optional

    use yomhook,                only : lhook, dr_hook

    use radiation_io,           only : nulerr, radiation_abort
    use radiation_config,       only : config_type
    use radiation_single_level, only : single_level_type
    use radsurf_properties,     only : surface_type, ITileFlat,ITileVegetation,ITileUrban3D
    use radiation_two_stream,   only : calc_two_stream_gammas_sw, &
         &                             calc_reflectance_transmittance_sw, &
         &                             calc_reflectance_transmittance_z_sw 

    use radiation_constants,    only: Pi, GasConstantDryAir, AccelDueToGravity

    class(surface_intermediate_type), intent(inout) :: this
    integer(kind=jpim),               intent(in)    :: istartcol, iendcol
    type(config_type),                intent(in)    :: config
    type(surface_type),               intent(in)    :: surface
    type(single_level_type),          intent(inout) :: single_level

    ! Input properties of the air in the lowest model level:
    ! extinction coefficient (m-1), single-scattering albedo and
    ! asymmetry factor
    real(kind=jprb), intent(in), dimension(config%n_g_sw,istartcol:iendcol), optional &
         &  :: ext_sw_air, ssa_sw_air, g_sw_air

    ! Shortwave region properties
    real(kind=jprb), dimension(this%nswbands) :: od_sw_region, ssa_sw_region, g_sw_region
    ! Extinction coefficient better for urban areas
    real(kind=jprb), dimension(this%nswbands) :: ext_sw_region
    ! Optical depth of air
    real(kind=jprb), dimension(this%nswbands) :: od_sw_air

    real(kind=jprb)    :: tile_fraction, cos_sza

    ! Exchange coefficients (m-1) for direct and diffuse radiation
    ! into building walls
    real(kind=jprb)    :: fdiff, fdir

    ! One minus the building fraction
    real(kind=jprb)    :: canyon_fraction

    ! Tangent of solar zenith angle
    real(kind=jprb)    :: tan_sza

    ! Two-stream coefficients
    real(jprb), dimension(this%nswbands) :: gamma1_sw, gamma2_sw
    real(jprb), dimension(this%nswbands) :: gamma0_sw, gamma3_sw, gamma4_sw

    real(jprb), dimension(this%nswbands) :: inv_denominator_sw

    integer(kind=jpim) :: jcol,jtile,iregion

    ! Mapping from albedo bands to reordered shortwave g points
    integer(kind=jpim) :: i_albedo_from_g(config%n_g_sw)

    ! Indices to different facets
    integer(kind=jpim) :: iground, iwall, iroof

    real(jprb) :: hook_handle

    if (lhook) call dr_hook('radsurf_intermediate:calc_boundary_conditions_sw',0,hook_handle)

    if (present(ext_sw_air)) then
      if (.not. present(ssa_sw_air) .or. .not. present(g_sw_air)) then
        write(nulerr,'(a)') '*** Error: if ext_sw_air is provided then ssa_sw_air and g_sw_air must also be provided'
        call radiation_abort()
      end if
    end if

    if (size(single_level%sw_albedo,2) /= this%nswbands) then
      write(nulerr,'(a,i0,a,i0)') '*** Error: single-level albedo has ', &
           &  size(single_level%sw_albedo,2), ' bands, needs ', this%nswbands
      call radiation_abort()
    end if

    if (config%use_canopy_full_spectrum_sw) then
      ! Put shortwave albedo on g-point grid and permute
      i_albedo_from_g = config%i_albedo_from_band_sw(config%i_band_from_reordered_g_sw)
      do jcol = istartcol,iendcol
        this%sw_albedo_diffuse(:,:,jcol) = surface%sw_albedo(jcol,i_albedo_from_g,:)
        if (allocated(surface%sw_albedo_direct)) then
          this%sw_albedo_direct(:,:,jcol) = surface%sw_albedo_direct(jcol,i_albedo_from_g,:)
        else
          this%sw_albedo_direct = this%sw_albedo_diffuse
        end if
      end do
    else
      do jcol = istartcol,iendcol
        ! No change to spectral resolution: simply permute 
        this%sw_albedo_diffuse(:,:,jcol) = surface%sw_albedo(jcol,:,:)
        if (allocated(surface%sw_albedo_direct)) then
          this%sw_albedo_direct(:,:,jcol) = surface%sw_albedo_direct(jcol,:,:)
        else
          this%sw_albedo_direct = this%sw_albedo_diffuse
        end if
      end do
    end if

    if (surface%is_simple) then
      ! We have a "traditional" representation: one flat tile
      single_level%sw_albedo_direct(istartcol:iendcol,:) &
           &  = transpose(this%sw_albedo_direct(:,1,istartcol:iendcol))
      single_level%sw_albedo       (istartcol:iendcol,:) &
           &  = transpose(this%sw_albedo_diffuse(:,1,istartcol:iendcol))
    else
      ! More complex description of surface

      ! Firstly initialize outputs to zero
      single_level%sw_albedo_direct(istartcol:iendcol,:) = 0.0_jprb
      single_level%sw_albedo       (istartcol:iendcol,:) = 0.0_jprb

      ! Loop over column and tile
      do jcol = istartcol,iendcol
        cos_sza = single_level%cos_sza(jcol)

        do jtile = 1,surface%ntile
          tile_fraction = surface%tile_fraction(jcol,jtile)

          if (tile_fraction > 0.0_jprb) then

            select case (surface%i_representation(jtile))
            case (ITileFlat)
              ! SIMPLE FLAT TILE

              ! Add the contribution from this simple flat tile to the
              ! accumulated values for the column
              iground = surface%i_ground_facet(jcol,jtile)
              single_level%sw_albedo_direct(jcol,:) = single_level%sw_albedo_direct(jcol,:) &
                   &  + tile_fraction * this%sw_albedo_direct(:,iground,jcol)
              single_level%sw_albedo(jcol,:) = single_level%sw_albedo(jcol,:) &
                   &  + tile_fraction * this%sw_albedo_diffuse(:,iground,jcol)

            case (ITileVegetation)
              ! HORIZONTALLY HOMOGENEOUS VEGETATION CANOPY WITH
              ! SELLERS-LIKE FORMULATION

              iground = surface%i_ground_facet(jcol,jtile)
              iregion = surface%i_region_1(jcol,jtile)

              ! Shortwave calculation
              if (present(ext_sw_air)) then
                od_sw_air = surface%canopy_depth(jcol,jtile)*ext_sw_air(:,jcol)
                od_sw_region = od_sw_air + surface%vegetation_optical_depth(jcol,jtile)
                ssa_sw_region= (ssa_sw_air(:,jcol)*od_sw_air &
                     &          + surface%vegetation_optical_depth(jcol,jtile) &
                     &            * surface%vegetation_sw_albedo(jcol,:,jtile)) &
                     &       / od_sw_region
                ! Assume asymmetry factor of vegetation is zero
                g_sw_region = g_sw_air(:,jcol)*ssa_sw_air(:,jcol)*od_sw_air &
                     &  / (ssa_sw_region*od_sw_region)
              else
                od_sw_region = surface%vegetation_optical_depth(jcol,jtile)
                ssa_sw_region= surface%vegetation_sw_albedo(jcol,:,jtile)
                g_sw_region  = 0.0_jprb
              end if
              call calc_two_stream_gammas_sw(this%nswbands, &
                   &  cos_sza, ssa_sw_region, g_sw_region, &
                   &  gamma1_sw, gamma2_sw, gamma3_sw)
              call calc_reflectance_transmittance_sw(this%nswbands, &
                   &  cos_sza, od_sw_region, ssa_sw_region, &
                   &  gamma1_sw, gamma2_sw, gamma3_sw, &
                   &  this%sw_ref_dif(:,iregion,jcol), this%sw_tra_dif(:,iregion,jcol), &
                   &  this%sw_ref_dir(:,iregion,jcol), this%sw_tra_dir_dif(:,iregion,jcol), &
                   &  this%sw_tra_dir_dir(:,iregion,jcol))

              ! Shortwave adding method for a single layer
              inv_denominator_sw(:) = 1.0_jprb &
                   &  / (1.0_jprb - this%sw_albedo_diffuse(:,iground,jcol) & 
                   &               *this%sw_ref_dif(:,iregion,jcol))
              this%sw_albedo_diffuse_reg(:,iregion,jcol) = this%sw_ref_dif(:,iregion,jcol) &
                   &  + this%sw_tra_dif(:,iregion,jcol)**2 &
                   &  * this%sw_albedo_diffuse(:,iground,jcol) * inv_denominator_sw(:)
              this%sw_albedo_direct_reg(:,iregion,jcol) = this%sw_ref_dir(:,iregion,jcol) &
                   &     + (this%sw_tra_dir_dir(:,iregion,jcol)*this%sw_albedo_direct (:,iground,jcol) &
                   &       +this%sw_tra_dir_dif(:,iregion,jcol)*this%sw_albedo_diffuse(:,iground,jcol)) &
                   &     * this%sw_tra_dif(:,iregion,jcol) * inv_denominator_sw(:)
              single_level%sw_albedo(jcol,:) = single_level%sw_albedo(jcol,:) &
                   &  + tile_fraction * this%sw_albedo_diffuse_reg(:,iregion,jcol)
              single_level%sw_albedo_direct(jcol,:) = single_level%sw_albedo_direct(jcol,:) &
                   &  + tile_fraction * this%sw_albedo_direct_reg(:,iregion,jcol)

            case (ITileUrban3D)
              ! URBAN CANOPY WITH NO VEGETATION, TREATED USING
              ! SPARTACUS METHODOLOGY

              iground  = surface%i_ground_facet(jcol,jtile)
              iwall    = surface%i_wall_facet(jcol,jtile)
              iroof    = surface%i_roof_facet(jcol,jtile)
              iregion = surface%i_region_1(jcol,jtile)

              canyon_fraction = 1.0_jprb - surface%building_fraction(jcol,jtile)

              fdiff = 0.5_jprb * surface%building_normalized_perimeter(jcol,jtile) / canyon_fraction

              tan_sza = sqrt(1.0_jprb / (cos_sza*cos_sza) - 1.0_jprb)
              fdir  = surface%building_normalized_perimeter(jcol,jtile) * tan_sza / (Pi * canyon_fraction)

              if (present(ext_sw_air)) then
                ext_sw_region = ext_sw_air(:,jcol)
                ssa_sw_region = ssa_sw_air(:,jcol)
                g_sw_region   = g_sw_air(:,jcol)
              else
                ext_sw_region = 0.0_jprb
                ssa_sw_region = 0.0_jprb
                g_sw_region   = 0.0_jprb
              end if

              ! Get gammas for atmosphere only
              call calc_two_stream_gammas_sw(this%nswbands, &
                   &  cos_sza, ssa_sw_region, g_sw_region, &
                   &  gamma1_sw, gamma2_sw, gamma3_sw)

              ! At this point gamma1_sw-gamma2_sw is the rate of
              ! absorption per unit optical depth of the air in the
              ! canyon
              this%sw_wall_abs_frac_dif(:,iregion,jcol) &
                   &  = fdiff * (1.0_jprb - this%sw_albedo_diffuse(:,iwall,jcol))
              this%sw_wall_abs_frac_dif(:,iregion,jcol) = this%sw_wall_abs_frac_dif(:,iregion,jcol) &
                   &  / max(1.0e-8_jprb,ext_sw_region*(gamma1_sw-gamma2_sw) &
                   &                    + this%sw_wall_abs_frac_dif(:,iregion,jcol))

              gamma4_sw = 1.0_jprb - gamma3_sw
              gamma0_sw = ext_sw_region / cos_sza + fdir
              gamma1_sw = ext_sw_region * gamma1_sw &
                   &    + fdiff * (1.0_jprb - 0.5_jprb*this%sw_albedo_diffuse(:,iwall,jcol))
              gamma2_sw = ext_sw_region * gamma2_sw &
                   &    + fdiff * 0.5_jprb*this%sw_albedo_diffuse(:,iwall,jcol)
              gamma3_sw = ext_sw_region * ssa_sw_region * gamma3_sw &
                   &    + 0.5_jprb * fdir * this%sw_albedo_direct(:,iwall,jcol)
              gamma4_sw = ext_sw_region * ssa_sw_region * gamma4_sw &
                   &    + 0.5_jprb * fdir * this%sw_albedo_direct(:,iwall,jcol)

              call calc_reflectance_transmittance_z_sw(this%nswbands, &
                   &  cos_sza, surface%canopy_depth(jcol,jtile), gamma0_sw, &
                   &  gamma1_sw, gamma2_sw, gamma3_sw, gamma3_sw, &
                   &  this%sw_ref_dif(:,iregion,jcol), this%sw_tra_dif(:,iregion,jcol), &
                   &  this%sw_ref_dir(:,iregion,jcol), this%sw_tra_dir_dif(:,iregion,jcol), &
                   &  this%sw_tra_dir_dir(:,iregion,jcol))

              ! Compute fraction of direct at canyon top that is absorbed by wall and air
              this%sw_wall_abs_dir(:,iregion,jcol) &
                   &  = (1.0_jprb - this%sw_tra_dir_dir(:,iregion,jcol)) &
                   &  * fdir * (1.0_jprb - this%sw_albedo_direct(:,iwall,jcol)) * cos_sza &
                   &  /  max(1.0e-8_jprb, fdir*cos_sza + ext_sw_region)
              this%sw_air_abs_dir(:,iregion,jcol)&
                   &  = (1.0_jprb - this%sw_tra_dir_dir(:,iregion,jcol)) &
                   &  * ext_sw_region * (1.0_jprb - ssa_sw_region) &
                   &  /  max(1.0e-8_jprb, fdir*cos_sza + ext_sw_region)

              ! Add roof component
              single_level%sw_albedo(jcol,:) = single_level%sw_albedo(jcol,:) &
                   &  + tile_fraction * surface%building_fraction(jcol,jtile) &
                   &  * this%sw_albedo_diffuse(:,iroof,jcol)
              single_level%sw_albedo_direct(jcol,:) = single_level%sw_albedo_direct(jcol,:) &
                   &  + tile_fraction * surface%building_fraction(jcol,jtile) &
                   &  * this%sw_albedo_direct(:,iroof,jcol)

              ! Add canyon component using shortwave adding method for a single layer
              inv_denominator_sw(:) = 1.0_jprb &
                   &  / (1.0_jprb - this%sw_albedo_diffuse(:,iground,jcol) & 
                   &               *this%sw_ref_dif(:,iregion,jcol))
              this%sw_albedo_diffuse_reg(:,iregion,jcol) = this%sw_ref_dif(:,iregion,jcol) &
                   &  + this%sw_tra_dif(:,iregion,jcol)**2 &
                   &  * this%sw_albedo_diffuse(:,iground,jcol) * inv_denominator_sw(:)
              this%sw_albedo_direct_reg(:,iregion,jcol) = this%sw_ref_dir(:,iregion,jcol) &
                   &     + (this%sw_tra_dir_dir(:,iregion,jcol)*this%sw_albedo_direct (:,iground,jcol) &
                   &       +this%sw_tra_dir_dif(:,iregion,jcol)*this%sw_albedo_diffuse(:,iground,jcol)) &
                   &     * this%sw_tra_dif(:,iregion,jcol) * inv_denominator_sw(:)
              single_level%sw_albedo(jcol,:) = single_level%sw_albedo(jcol,:) &
                   &  + tile_fraction * canyon_fraction * this%sw_albedo_diffuse_reg(:,iregion,jcol)
              single_level%sw_albedo_direct(jcol,:) = single_level%sw_albedo_direct(jcol,:) &
                   &  + tile_fraction * canyon_fraction * this%sw_albedo_direct_reg(:,iregion,jcol)

            end select
          end if
        end do
      end do
    end if

    if (lhook) call dr_hook('radsurf_intermediate:calc_boundary_conditions_sw',1,hook_handle)

  end subroutine calc_boundary_conditions_sw


  !---------------------------------------------------------------------
  ! As calc_boundary_conditions_sw, but for the longwave.
  subroutine calc_boundary_conditions_lw(this, istartcol, iendcol, &   ! in
       &  config, surface, &                    ! in
       &  single_level, &                                 ! out
       &  ext_lw_air, ssa_lw_air, g_lw_air)               ! in, optional

    use yomhook,                only : lhook, dr_hook

    use radiation_io,           only : nulerr, radiation_abort
    use radiation_config,       only : config_type
    use radiation_single_level, only : single_level_type
    use radsurf_properties,     only : surface_type, ITileFlat,ITileVegetation,ITileUrban3D
    use radiation_two_stream,   only : calc_two_stream_gammas_lw, &
         &                             calc_reflectance_transmittance_isothermal_lw, &
         &                             LwDiffusivity
    use radiation_ifs_rrtm,     only : planck_function
    use radiation_constants,    only: Pi, GasConstantDryAir, AccelDueToGravity, StefanBoltzmann

    class(surface_intermediate_type), intent(inout) :: this
    integer(kind=jpim),               intent(in)    :: istartcol, iendcol
    type(config_type),                intent(in)    :: config
    type(surface_type),               intent(in)    :: surface
    type(single_level_type),          intent(inout) :: single_level

    ! Input properties of the air in the lowest model level:
    ! extinction coefficient (m-1), single-scattering albedo and
    ! asymmetry factor
    real(kind=jprb), intent(in), dimension(config%n_g_lw,istartcol:iendcol), optional &
         &  :: ext_lw_air, ssa_lw_air, g_lw_air

    ! Longwave region properties
    real(kind=jprb), dimension(this%nlwbands) :: od_lw_region, ssa_lw_region, g_lw_region

    ! Optical depth of air
    real(kind=jprb), dimension(this%nlwbands) :: od_lw_air

    ! Effective Planck function of urban canopy as a weighted average
    ! of wall and air emission
    real(kind=jprb), dimension(this%nlwbands) :: planck_canopy

    ! Vegetation emissivity using local spectral representation
    real(kind=jprb), dimension(this%nlwbands) :: vegetation_lw_emissivity

    real(kind=jprb)    :: tile_fraction

    ! Exchange coefficient (m-1) for diffuse radiation into building
    ! walls, multiplied by canyon depth (m) to get a dimensionless
    ! analogue for optical depth
    real(kind=jprb)    :: od_lw_wall

    ! One minus the building fraction
    real(kind=jprb)    :: canyon_fraction

    ! Two-stream coefficients
    real(jprb), dimension(this%nlwbands) :: gamma1_lw, gamma2_lw

    real(jprb), dimension(this%nlwbands) :: inv_denominator_lw

    integer(kind=jpim) :: jcol,jtile,jfacet,iregion

    ! Indices to different facets
    integer(kind=jpim) :: iground, iwall, iroof

    real(jprb) :: hook_handle

    if (lhook) call dr_hook('radsurf_intermediate:calc_boundary_conditions_lw',0,hook_handle)

    if (present(ssa_lw_air)) then
      if (.not. present(g_lw_air)) then
        write(nulerr,'(a)') '*** Error: if ssa_lw_air is provided then g_lw_air must also be provided'
        call radiation_abort()
      end if
    end if

    if (size(single_level%lw_emissivity,2) /= this%nlwbands) then
      write(nulerr,'(a,i0,a,i0)') '*** Error: single-level emissivity has ', &
           &  size(single_level%lw_emissivity,2), ' bands, needs ', this%nlwbands
      call radiation_abort()
    end if

    ! FIX: This section ought to check what tiles are present in each column
    if (config%use_canopy_full_spectrum_lw) then
      ! Put longwave emissivity on g-point grid and permute      
      do jcol = istartcol,iendcol
        this%lw_emissivity(:,:,jcol) = surface%lw_emissivity(jcol,config%i_emiss_from_band_lw( &
             &                           config%i_band_from_reordered_g_lw),:)
      end do
      ! Compute Planck function at each facet
      do jfacet = 1,surface%nfacet
        do jcol = istartcol,iendcol
          ! FIX this function assumes contiguous data for planck_facet
          ! so copies into a temporary
          call planck_function(config, &
               &  surface%skin_temperature(jcol,jfacet), &
               &  this%planck_facet(:,jfacet,jcol))
        end do
      end do
      do jtile = 1,surface%ntile
        do jcol = istartcol,iendcol
          iregion = surface%i_region_1(jcol,jtile)
          if (iregion > 0) then
            call planck_function(config, &
                 &  surface%canopy_temperature(jcol,jtile), &
                 &  this%planck_region(:,iregion,jcol:jcol))
          end if
        end do
      end do
    else
      if (surface%nemissbands /= 1) then
        write(nulerr,'(a)') '*** Error: insufficient information to compute Planck function in emissivity bands'
        call radiation_abort()
      end if
      do jcol = istartcol,iendcol
        ! No change to spectral resolution: simply permute
        this%lw_emissivity(:,:,jcol) = surface%lw_emissivity(jcol,:,:)
      end do
      ! Broadband calculation: Stefan-Boltzmann law
      do jfacet = 1,surface%nfacet
        this%planck_facet(1,jfacet,istartcol:iendcol) &
             &  = StefanBoltzmann*surface%skin_temperature(istartcol:iendcol,jfacet)**4
      end do
      do jtile = 1,surface%ntile
        do jcol = istartcol,iendcol
          iregion = surface%i_region_1(jcol,jtile)
          if (iregion > 0) then
            this%planck_region(1,iregion,jcol) &
                 &  = StefanBoltzmann*surface%canopy_temperature(jcol,jtile)
          end if
        end do
      end do
    end if

    if (surface%is_simple) then
      ! We have a "traditional" representation: one flat tile
      single_level%lw_emissivity(istartcol:iendcol,:) &
           &  = transpose(this%lw_emissivity(:,1,istartcol:iendcol))
      single_level%lw_emission(istartcol:iendcol,:) &
           &  = transpose(this%planck_facet(:,1,istartcol:iendcol) &
           &              * this%lw_emissivity(:,1,istartcol:iendcol))

    else
      ! More complex description of surface

      ! Firstly initialize outputs to zero
      single_level%lw_emissivity(istartcol:iendcol,:) = 0.0_jprb
      single_level%lw_emission  (istartcol:iendcol,:) = 0.0_jprb

      ! Loop over column and tile
      do jcol = istartcol,iendcol

        do jtile = 1,surface%ntile
          tile_fraction = surface%tile_fraction(jcol,jtile)

          if (tile_fraction > 0.0_jprb) then

            select case (surface%i_representation(jtile))
            case (ITileFlat)
              ! SIMPLE FLAT TILE

              ! Add the contribution from this simple flat tile to the
              ! accumulated values for the column
              iground = surface%i_ground_facet(jcol,jtile)
              single_level%lw_emissivity(jcol,:) = single_level%lw_emissivity(jcol,:) &
                   &  + tile_fraction*this%lw_emissivity(:,iground,jcol)
              single_level%lw_emission(jcol,:) = single_level%lw_emission(jcol,:) + tile_fraction &
                   &  * (this%planck_facet(:,iground,jcol)*this%lw_emissivity(:,iground,jcol))

            case (ITileVegetation)
              ! HORIZONTALLY HOMOGENEOUS VEGETATION CANOPY

              iground  = surface%i_ground_facet(jcol,jtile)
              iregion = surface%i_region_1(jcol,jtile)

              ! Convert vegetation emissivity to local spectral representation
              if (config%use_canopy_full_spectrum_lw) then
                vegetation_lw_emissivity = surface%vegetation_lw_emissivity(jcol, &
                     &  config%i_emiss_from_band_lw(config%i_band_from_reordered_g_lw), &
                     &  jtile)
              else
                vegetation_lw_emissivity = surface%vegetation_lw_emissivity(jcol,:,jtile)
              end if

              if (present(ext_lw_air)) then
                od_lw_air = surface%canopy_depth(jcol,jtile)*ext_lw_air(:,jcol)
                od_lw_region = od_lw_air + surface%vegetation_optical_depth(jcol,jtile)
                if (present(ssa_lw_air)) then
                  ssa_lw_region= (ssa_lw_air(:,jcol)*od_lw_air &
                       &          + surface%vegetation_optical_depth(jcol,jtile) &
                       &            * (1.0_jprb-vegetation_lw_emissivity)) &
                       &       / od_lw_region
                  ! Assume asymmetry factor of vegetation is zero
                  g_lw_region = g_lw_air(:,jcol)*ssa_lw_air(:,jcol)*od_lw_air &
                       &  / (ssa_lw_region * od_lw_region)
                else
                  ! No longwave scattering properties for air
                  ssa_lw_region= surface%vegetation_optical_depth(jcol,jtile) &
                       &            * (1.0_jprb-vegetation_lw_emissivity) &
                       &       / od_lw_region
                  ! Assume asymmetry factor of vegetation is zero
                  g_lw_region = 0.0_jprb

                end if
              else
                ! No gases
                od_lw_region = surface%vegetation_optical_depth(jcol,jtile)
                ssa_lw_region= 1.0_jprb-vegetation_lw_emissivity
                g_lw_region  = 0.0_jprb
              end if

              call calc_two_stream_gammas_lw(this%nlwbands, &
                   &  ssa_lw_region, g_lw_region, &
                   &  gamma1_lw, gamma2_lw)
              call calc_reflectance_transmittance_isothermal_lw(this%nlwbands, &
                   &  od_lw_region, gamma1_lw, gamma2_lw, this%planck_region(:,iregion,jcol), &
                   &  this%lw_reflectance(:,iregion,jcol), this%lw_transmittance(:,iregion,jcol), &
                   &  this%lw_source(:,iregion,jcol))
              
              ! Longwave adding method for a single layer
              inv_denominator_lw(:) = 1.0_jprb &
                   &  / (1.0_jprb - (1.0_jprb - this%lw_emissivity(:,iground,jcol)) &
                   &               *this%lw_reflectance(:,iregion,jcol))
              single_level%lw_emissivity(jcol,:) = single_level%lw_emissivity(jcol,:) &
                   &  + tile_fraction * (1.0_jprb - (this%lw_reflectance(:,iregion,jcol) &
                   &  + this%lw_transmittance(:,iregion,jcol)**2 &
                   &  * (1.0_jprb - this%lw_emissivity(:,iground,jcol)) * inv_denominator_lw(:)))
              single_level%lw_emission(jcol,:) = single_level%lw_emission(jcol,:) + tile_fraction &
                   &  * (this%lw_source(:,iregion,jcol) * (1.0_jprb + inv_denominator_lw(:) &
                   &       * (1.0_jprb - this%lw_emissivity(:,iground,jcol)) &
                   &       * this%lw_transmittance(:,iregion,jcol)) &
                   &     +this%planck_facet(:,iground,jcol)*this%lw_emissivity(:,iground,jcol) &
                   &     *this%lw_transmittance(:,iregion,jcol)*inv_denominator_lw(:))

            case (ITileUrban3D)
              ! URBAN CANOPY WITH NO VEGETATION, TREATED USING
              ! SPARTACUS METHODOLOGY

              iground  = surface%i_ground_facet(jcol,jtile)
              iwall    = surface%i_wall_facet(jcol,jtile)
              iroof    = surface%i_roof_facet(jcol,jtile)
              iregion  = surface%i_region_1(jcol,jtile)

              canyon_fraction = 1.0_jprb - surface%building_fraction(jcol,jtile)

              ! fdiff woould be 0.5 * building_normalized_perimeter /
              ! canyon_fraction, but to get equivalent zenith optical
              ! depth we multiply by the building height and divide by
              ! the longwave diffusivity angle
!              od_lw_wall = 0.5_jprb * surface%building_normalized_perimeter(jcol,jtile) &
!                   &     * surface%canopy_depth(jcol,jtile) / (LwDiffusivity * canyon_fraction)

              ! Or first compute H/W
              od_lw_wall = 0.5_jprb * surface%building_normalized_perimeter(jcol,jtile) &
                   &     * surface%canopy_depth(jcol,jtile) / canyon_fraction
              ! And then use the Harman et al. (2004) formula for
              ! street-to-sky transmittance T=sqrt[(H/W)^2+1]-H/W, to
              ! compute equivalent optical depth knowing that the
              ! two-stream scheme treatment will be T=exp(-D*od)
              od_lw_wall = -log(sqrt(od_lw_wall*od_lw_wall + 1) - od_lw_wall) / LwDiffusivity

              if (present(ext_lw_air)) then
                od_lw_air     = ext_lw_air(:,jcol)*surface%canopy_depth(jcol,jtile) !*10.0_jprb
                od_lw_region  = od_lw_air + od_lw_wall
                if (present(ssa_lw_air)) then
                  ssa_lw_region = (od_lw_air * ssa_lw_air(:,jcol) &
                       &    + od_lw_wall * (1.0_jprb - this%lw_emissivity(:,iwall,jcol))) &
                       &               / max(od_lw_region,1.0e-6_jprb)
                  ! We assume that any scattering off the walls is
                  ! equally likely to be up or down, so the effective
                  ! asymmetry factor is zero
                  g_lw_region = (od_lw_air * ssa_lw_air(:,jcol)*g_lw_air(:,jcol)) &
                       &      / max(od_lw_region*ssa_lw_region,1.0e-6_jprb)
                  ! Effective Planck function of the canopy is the
                  ! weighted average of wall and air emission
                  this%lw_total_wall_emission(:,iregion,jcol) = LwDiffusivity &
                       &  * od_lw_wall*this%lw_emissivity(:,iwall,jcol)*this%planck_facet(:,iwall,jcol)
                  this%lw_total_canopy_emission(:,iregion,jcol) = LwDiffusivity &
                       &  * od_lw_air*(1.0_jprb-ssa_lw_air(:,jcol))*this%planck_region(:,iregion,jcol)
                  planck_canopy = (this%lw_total_wall_emission(:,iregion,jcol) &
                       &           +this%lw_total_canopy_emission(:,iregion,jcol)) &
                       &        / max(od_lw_region*(1.0_jprb-ssa_lw_region)*LwDiffusivity,1.0e-6_jprb)
                else
                  ssa_lw_region = od_lw_wall * (1.0_jprb - this%lw_emissivity(:,iwall,jcol)) &
                       &        / max(od_lw_region,1.0e-6_jprb)
                  ! We assume that any scattering off the walls is
                  ! equally likely to be up or down, so the effective
                  ! asymmetry factor is zero
                  g_lw_region = 0.0_jprb
                  ! Effective Planck function of the canopy is the
                  ! weighted average of wall and air emission
                  this%lw_total_wall_emission(:,iregion,jcol) = LwDiffusivity &
                       &  * od_lw_wall*this%lw_emissivity(:,iwall,jcol)*this%planck_facet(:,iwall,jcol)
                  this%lw_total_canopy_emission(:,iregion,jcol) = LwDiffusivity &
                       &  * od_lw_air*this%planck_region(:,iregion,jcol)
                  planck_canopy = (this%lw_total_wall_emission(:,iregion,jcol) &
                       &           +this%lw_total_canopy_emission(:,iregion,jcol)) &
                       &        / max(od_lw_region*(1.0_jprb-ssa_lw_region)*LwDiffusivity,1.0e-6_jprb)
                end if

                ! Compute fraction of canyon absorption by the wall 
                this%lw_wall_abs_frac(:,iregion,jcol) = od_lw_wall * this%lw_emissivity(:,iwall,jcol) &
                     &    / max(od_lw_region*(1.0_jprb-ssa_lw_region),1.0e-6_jprb)


              else
                od_lw_region  = od_lw_wall
                ssa_lw_region = 1.0_jprb - this%lw_emissivity(:,iwall,jcol);
                g_lw_region   = 0.0_jprb

                ! All absorption and emission is by the wall
                this%lw_wall_abs_frac(:,iregion,jcol) = 1.0_jprb
                this%lw_total_wall_emission(:,iregion,jcol) = LwDiffusivity &
                     &  * od_lw_wall*this%lw_emissivity(:,iwall,jcol)*this%planck_facet(:,iwall,jcol)
                this%lw_total_canopy_emission(:,iregion,jcol) = 0.0_jprb
                planck_canopy = this%planck_facet(:,iwall,jcol)
              end if

              call calc_two_stream_gammas_lw(this%nlwbands, &
                   &  ssa_lw_region, g_lw_region, &
                   &  gamma1_lw, gamma2_lw)
              call calc_reflectance_transmittance_isothermal_lw(this%nlwbands, &
                   &  od_lw_region, gamma1_lw, gamma2_lw, planck_canopy, &
                   &  this%lw_reflectance(:,iregion,jcol), this%lw_transmittance(:,iregion,jcol), &
                   &  this%lw_source(:,iregion,jcol))

              ! Add roof component
              single_level%lw_emissivity(jcol,:) = single_level%lw_emissivity(jcol,:) &
                   &  + tile_fraction * surface%building_fraction(jcol,jtile) &
                   &  * this%lw_emissivity(:,iroof,jcol)
              single_level%lw_emission(jcol,:) = single_level%lw_emission(jcol,:) &
                   &  + tile_fraction * surface%building_fraction(jcol,jtile) &
                   &  * this%lw_emissivity(:,iroof,jcol)*this%planck_facet(:,iroof,jcol)
              
              ! Add canyon component: longwave adding method for a single layer
              inv_denominator_lw(:) = 1.0_jprb &
                   &  / (1.0_jprb - (1.0_jprb - this%lw_emissivity(:,iground,jcol)) &
                   &               *this%lw_reflectance(:,iregion,jcol))
              this%lw_emissivity_region(:,iregion,jcol) &
                   &  = (1.0_jprb - (this%lw_reflectance(:,iregion,jcol) &
                   &  + this%lw_transmittance(:,iregion,jcol)**2 &
                   &  * (1.0_jprb - this%lw_emissivity(:,iground,jcol)) * inv_denominator_lw(:)))
              single_level%lw_emissivity(jcol,:) = single_level%lw_emissivity(jcol,:) &
                   &  + tile_fraction * canyon_fraction * this%lw_emissivity_region(:,iregion,jcol) 
              this%lw_emission_region(:,iregion,jcol) &
                   &  = this%lw_source(:,iregion,jcol) * (1.0_jprb + inv_denominator_lw(:) &
                   &       * (1.0_jprb - this%lw_emissivity(:,iground,jcol)) &
                   &       * this%lw_transmittance(:,iregion,jcol)) &
                   &     +this%planck_facet(:,iground,jcol)*this%lw_emissivity(:,iground,jcol) &
                   &     *this%lw_transmittance(:,iregion,jcol)*inv_denominator_lw(:)
              single_level%lw_emission(jcol,:) = single_level%lw_emission(jcol,:) &
                   &  + tile_fraction * canyon_fraction * this%lw_emission_region(:,iregion,jcol)

            end select
          end if
        end do
      end do
    end if

    if (lhook) call dr_hook('radsurf_intermediate:calc_boundary_conditions_lw',1,hook_handle)

  end subroutine calc_boundary_conditions_lw


  !---------------------------------------------------------------------
  ! Compute both shortwave and longwave boundary conditions neglecting
  ! gases within vegetation and urban canopies
  subroutine calc_boundary_conditions_vacuum(this, istartcol, iendcol, &   ! in
       &                              config, surface,          &   ! in
       &                              single_level)                 ! out

    use radiation_config,       only : config_type
    use radiation_single_level, only : single_level_type
    use radsurf_properties,     only : surface_type

    class(surface_intermediate_type), intent(inout) :: this
    integer(kind=jpim),               intent(in)    :: istartcol, iendcol
    type(config_type),                intent(in)    :: config
    type(surface_type),               intent(in)    :: surface
    type(single_level_type),          intent(inout) :: single_level

    call this%calc_boundary_conditions_sw(istartcol, iendcol, &
         &                                config, surface, single_level)
    call this%calc_boundary_conditions_lw(istartcol, iendcol, &
         &                                config, surface, single_level)

  end subroutine calc_boundary_conditions_vacuum


  !---------------------------------------------------------------------
  ! Compute both shortwave and longwave boundary conditions
  subroutine calc_boundary_conditions(this, istartcol, iendcol, &   ! in
       &                              config, surface,          &   ! in
       &                              thermodynamics, gas,      &   ! in
       &                              single_level)                 ! out

    use radiation_config,       only : config_type
    use radiation_thermodynamics,only: thermodynamics_type
    use radiation_gas,          only : gas_type
    use radiation_single_level, only : single_level_type
    use radiation_ifs_rrtm,     only : gas_optics
    use radsurf_properties,     only : surface_type
    use radiation_constants,    only : GasConstantDryAir, &
         &                             AccelDueToGravity

    class(surface_intermediate_type), intent(inout) :: this
    integer(kind=jpim),               intent(in)    :: istartcol, iendcol
    type(config_type),                intent(in)    :: config
    type(surface_type),               intent(in)    :: surface
    type(gas_type),                   intent(in)    :: gas
    type(thermodynamics_type),        intent(in)    :: thermodynamics
    type(single_level_type),          intent(inout) :: single_level

    ! Ratio of gas constant for dry air to acceleration due to gravity
    real(jprb), parameter :: R_over_g = GasConstantDryAir / AccelDueToGravity

    ! Number of leves to request gas optical depths for; we only need
    ! one but gas optics scheme assumes more
    integer, parameter :: NGasLevels = 1

    ! Canopy longwave optical depth
    real(jprb), dimension(config%n_g_lw,NGasLevels,istartcol:iendcol) :: od_lw

    ! Canopy longwave extinction coefficient
    real(jprb), dimension(config%n_g_lw,istartcol:iendcol) :: ext_lw

    ! Canopy shortwave optical depth, single scattering albedo and
    ! asymmetry factor of gases and aerosols at each shortwave g-point
    real(jprb), dimension(config%n_g_sw,NGasLevels,istartcol:iendcol) :: od_sw, ssa_sw, g_sw

    ! Thickness of lowest model level
    real(jprb) :: layer_depth

    ! Index to final level of full model grid
    integer :: iendlev

    ! Column loop counter, number of columns
    integer :: jcol, ncol

    ncol    = ubound(thermodynamics%pressure_hl,1)
    iendlev = ubound(thermodynamics%pressure_hl,2)-1

    call gas_optics(ncol, NGasLevels, istartcol, iendcol, &
         &  config, single_level, thermodynamics, gas, & 
         &  od_lw, od_sw, ssa_sw)

    ! Scale optical depths to extinction
    do jcol = istartcol,iendcol
      layer_depth = R_over_g &
             &  * (thermodynamics%pressure_hl(jcol,iendlev+1) &
             &     - thermodynamics%pressure_hl(jcol,iendlev)) &
             &  * (thermodynamics%temperature_hl(jcol,iendlev) &
             &     + thermodynamics%temperature_hl(jcol,iendlev+1)) &
             &  / (thermodynamics%pressure_hl(jcol,iendlev) &
             &     + thermodynamics%pressure_hl(jcol,iendlev+1))
      !if (config%do_sw) then
      !  ext_sw(jcol) = ext_sw(jcol) / layer_depth
      !end if
      if (config%do_lw) then
        ext_lw(:,jcol) = od_lw(:,NGasLevels,jcol) / layer_depth
      end if
    end do

    call this%calc_boundary_conditions_sw(istartcol, iendcol, &
         &                                config, surface, single_level)
    if (config%do_canopy_gases_lw) then
      call this%calc_boundary_conditions_lw(istartcol, iendcol, &
           &                                config, surface, single_level, &
           &                                ext_lw_air=ext_lw)
    else
      call this%calc_boundary_conditions_lw(istartcol, iendcol, &
           &                                config, surface, single_level)
    end if

  end subroutine calc_boundary_conditions


  !---------------------------------------------------------------------
  subroutine partition_fluxes(this, istartcol, iendcol, config, surface, flux, &
       &                      surface_flux)

    use yomhook,           only : lhook, dr_hook

    use radiation_flux,    only : flux_type
    use radsurf_flux,      only : surface_flux_type
    use radiation_config,  only : config_type
    use radsurf_properties,only : surface_type, ITileFlat,ITileVegetation,ITileUrban3D

    class(surface_intermediate_type), intent(in)    :: this
    integer(kind=jpim),               intent(in)    :: istartcol, iendcol
    type(config_type),                intent(in)    :: config
    type(surface_type),               intent(in)    :: surface
    type(flux_type),                  intent(in)    :: flux
    type(surface_flux_type),          intent(inout) :: surface_flux

    real(kind=jprb), dimension(config%n_g_lw)       :: lw_dn_g, lw_up_g
    real(kind=jprb), dimension(config%n_canopy_bands_sw) &
         &  :: sw_dn_diffuse_g, sw_dn_direct_g, sw_up_g, sw_abs_g
    real(kind=jprb), dimension(config%n_canopy_bands_lw) :: lw_abs_g

    ! Ratio of street planar area to wall frontal area
    real(kind=jprb)    :: wall_scaling

    real(kind=jprb)    :: tile_fraction
    integer(kind=jpim) :: iground, iroof, iwall, iregion
    integer(kind=jpim) :: isurf ! index to lowest flux level (=nlev+1)
    integer(kind=jpim) :: jcol, jtile

    real(kind=jprb) :: hook_handle

    if (lhook) call dr_hook('radsurf_intermediate:partition_fluxes',0,hook_handle)

    if (.not. surface%is_simple) then
      isurf = size(flux%lw_dn,2)

      if (config%do_lw) then
        surface_flux%lw_dn_facet(istartcol:iendcol,:) = 0.0_jprb
        surface_flux%lw_up_facet(istartcol:iendcol,:) = 0.0_jprb
        surface_flux%lw_abs_canopy(istartcol:iendcol,:)=0.0_jprb
      end if
      if (config%do_sw) then
        surface_flux%sw_dn_facet(istartcol:iendcol,:) = 0.0_jprb
        surface_flux%sw_dn_direct_facet(istartcol:iendcol,:) = 0.0_jprb
        surface_flux%sw_up_facet(istartcol:iendcol,:) = 0.0_jprb
        surface_flux%sw_abs_canopy(istartcol:iendcol,:)=0.0_jprb
      end if

      ! Loop over column and tile
      do jcol = istartcol,iendcol
        do jtile = 1,surface%ntile
          tile_fraction = surface%tile_fraction(jcol,jtile)

          if (tile_fraction > 0.0_jprb) then

            select case (surface%i_representation(jtile))
            case (ITileFlat)
              iground = surface%i_ground_facet(jcol,jtile)
              if (config%do_lw) then
                surface_flux%lw_dn_facet(jcol,iground) = flux%lw_dn(jcol,isurf)
                surface_flux%lw_up_facet(jcol,iground) &
                     &  = sum(this%lw_emissivity(:,iground,jcol)*this%planck_facet(:,iground,jcol) &
                     &        + (1.0_jprb - this%lw_emissivity(:,iground,jcol)) &
                     &        * flux%lw_dn_surf_canopy(:,jcol))
              end if
              if (config%do_sw) then
                surface_flux%sw_dn_facet(jcol,iground)        = flux%sw_dn(jcol,isurf)
                surface_flux%sw_dn_direct_facet(jcol,iground) = flux%sw_dn_direct(jcol,isurf)
                surface_flux%sw_up_facet(jcol,iground) &
                     & = sum(this%sw_albedo_diffuse(:,iground,jcol) &
                     &       * flux%sw_dn_diffuse_surf_canopy(:,jcol) &
                     &       + this%sw_albedo_direct (:,iground,jcol) &
                     &       * flux%sw_dn_direct_surf_canopy (:,jcol))
              end if

            case (ITileVegetation)
              iground  = surface%i_ground_facet(jcol,jtile)
              iregion = surface%i_region_1(jcol,jtile)
              ! Adding method in longwave and shortwave
              if (config%do_lw) then
                ! Surface downwelling fluxes at each g point
                lw_dn_g = (this%lw_transmittance(:,iregion,jcol)*flux%lw_dn_surf_canopy(:,jcol) &
                     &    +this%lw_reflectance  (:,iregion,jcol) &
                     &     *this%lw_emissivity(:,iground,jcol)*this%planck_facet(:,iground,jcol) &
                     &    +this%lw_source(:,iregion,jcol)) &
                     &  / (1.0_jprb - (1.0_jprb - this%lw_emissivity(:,iground,jcol)) &
                     &               *this%lw_reflectance(:,iregion,jcol))
                surface_flux%lw_dn_facet(jcol,iground) = sum(lw_dn_g)
                surface_flux%lw_up_facet(jcol,iground) &
                     &  = sum((1.0_jprb-this%lw_emissivity(:,iground,jcol))*lw_dn_g &
                     &       +this%lw_emissivity(:,iground,jcol)*this%planck_facet(:,iground,jcol))
                surface_flux%lw_abs_canopy(jcol,jtile) = flux%lw_dn(jcol,isurf) &
                     &  - flux%lw_up(jcol,isurf) - surface_flux%lw_dn_facet(jcol,iground) &
                     &  + surface_flux%lw_up_facet(jcol,iground)

              end if
              if (config%do_sw) then
                sw_dn_direct_g = this%sw_tra_dir_dir(:,iregion,jcol) &
                     &  * flux%sw_dn_direct_surf_canopy(:,jcol)
                ! Note that the following is initially just the
                ! upwelling due to scattering of the direct beam
                sw_up_g = sw_dn_direct_g * this%sw_albedo_direct(:,iground,jcol)
                sw_dn_diffuse_g &
                     &  = (this%sw_tra_dif(:,iregion,jcol)*flux%sw_dn_diffuse_surf_canopy(:,jcol) &
                     &    +this%sw_ref_dif(:,iregion,jcol)*sw_up_g &
                     &    +this%sw_tra_dir_dif(:,iregion,jcol)*flux%sw_dn_direct_surf_canopy(:,jcol)) &
                     &  / (1.0_jprb - this%sw_albedo_diffuse(:,iground,jcol) & 
                     &               *this%sw_ref_dif(:,iregion,jcol))
                sw_up_g = sw_up_g + sw_dn_diffuse_g * this%sw_albedo_diffuse(:,iground,jcol)
                surface_flux%sw_dn_direct_facet(jcol,iground) = sum(sw_dn_direct_g)
                surface_flux%sw_dn_facet(jcol,iground) = surface_flux%sw_dn_direct_facet(jcol,iground) &
                     &  + sum(sw_dn_diffuse_g)
                surface_flux%sw_up_facet(jcol,iground) = sum(sw_up_g)
                surface_flux%sw_abs_canopy(jcol,jtile) &
                     &  = flux%sw_dn(jcol,isurf) - flux%sw_up(jcol,isurf) &
                     &  - surface_flux%sw_dn_facet(jcol,iground) &
                     &  + surface_flux%sw_up_facet(jcol,iground)
              end if
            case (ITileUrban3D)
              iground = surface%i_ground_facet(jcol,jtile)
              iroof   = surface%i_roof_facet(jcol,jtile)
              iwall   = surface%i_wall_facet(jcol,jtile)
              iregion = surface%i_region_1(jcol,jtile)

              ! We want wall fluxes per unit area of wall, not per
              ! unit area of street
              wall_scaling = (1.0_jprb - surface%building_fraction(jcol,jtile)) &
                   &  / max(1.0e-4_jprb, surface%building_normalized_perimeter(jcol,jtile) &
                   &                     * surface%canopy_depth(jcol,jtile))

              if (config%do_sw) then
                ! Roof fluxes
                surface_flux%sw_dn_facet(jcol,iroof)        = flux%sw_dn(jcol,isurf)
                surface_flux%sw_dn_direct_facet(jcol,iroof) = flux%sw_dn_direct(jcol,isurf)
                surface_flux%sw_up_facet(jcol,iroof) &
                     & = sum(this%sw_albedo_diffuse(:,iroof,jcol) &
                     &       * flux%sw_dn_diffuse_surf_canopy(:,jcol) &
                     &       + this%sw_albedo_direct (:,iroof,jcol) &
                     &       * flux%sw_dn_direct_surf_canopy (:,jcol))

                ! Ground fluxes
                sw_dn_direct_g = this%sw_tra_dir_dir(:,iregion,jcol) &
                     &         * flux%sw_dn_direct_surf_canopy(:,jcol)
                ! Note that the following is initially just the
                ! upwelling due to scattering of the direct beam
                sw_up_g = sw_dn_direct_g * this%sw_albedo_direct(:,iground,jcol)
                sw_dn_diffuse_g &
                     &  = (this%sw_tra_dif(:,iregion,jcol)*flux%sw_dn_diffuse_surf_canopy(:,jcol) &
                     &    +this%sw_ref_dif(:,iregion,jcol)*sw_up_g &
                     &    +this%sw_tra_dir_dif(:,iregion,jcol)*flux%sw_dn_direct_surf_canopy(:,jcol)) &
                     &  / (1.0_jprb - this%sw_albedo_diffuse(:,iground,jcol) & 
                     &               *this%sw_ref_dif(:,iregion,jcol))
                sw_up_g = sw_up_g + sw_dn_diffuse_g * this%sw_albedo_diffuse(:,iground,jcol)
                surface_flux%sw_dn_direct_facet(jcol,iground) = sum(sw_dn_direct_g)
                surface_flux%sw_dn_facet(jcol,iground) = surface_flux%sw_dn_direct_facet(jcol,iground) &
                     &  + sum(sw_dn_diffuse_g)
                surface_flux%sw_up_facet(jcol,iground) = sum(sw_up_g)

                ! Wall fluxes
                sw_abs_g = flux%sw_dn_direct_surf_canopy(:,jcol)*this%sw_wall_abs_dir(:,iregion,jcol)
                surface_flux%sw_dn_direct_facet(jcol,iwall) &
                     &  = wall_scaling * sum(sw_abs_g / (1.0_jprb - this%sw_albedo_direct(:,iwall,jcol)))
                ! Initially just the direct reflection
                surface_flux%sw_up_facet(jcol,iwall) = wall_scaling &
                     &  * sum(sw_abs_g * this%sw_albedo_direct(:,iwall,jcol) &
                     &        / (1.0_jprb - this%sw_albedo_direct(:,iwall,jcol)))
                ! Initially just the direct absorption
                surface_flux%sw_abs_canopy(jcol,jtile) &
                     &  = sum(flux%sw_dn_direct_surf_canopy(:,jcol)*this%sw_air_abs_dir(:,iregion,jcol))

                ! Diffuse absorption within canopy
                sw_abs_g = flux%sw_dn_direct_surf_canopy(:,jcol) &
                     &   * (1.0-this%sw_albedo_direct_reg (:,iregion,jcol)) &
                     &   + flux%sw_dn_diffuse_surf_canopy(:,jcol) &
                     &   * (1.0-this%sw_albedo_diffuse_reg(:,iregion,jcol)) &
                     &   - sw_dn_direct_g - sw_dn_diffuse_g + sw_up_g - sw_abs_g
                ! Add diffuse absorption
                surface_flux%sw_abs_canopy(jcol,jtile) = surface_flux%sw_abs_canopy(jcol,jtile) &
                     &   + sum(sw_abs_g * (1.0_jprb - this%sw_wall_abs_frac_dif(:,iregion,jcol)))
                ! Add diffuse reflection
                surface_flux%sw_up_facet(jcol,iwall) = surface_flux%sw_up_facet(jcol,iwall) &
                     &   +  wall_scaling * sum(sw_abs_g * this%sw_wall_abs_frac_dif(:,iregion,jcol) &
                     &         * this%sw_albedo_diffuse(:,iwall,jcol) &
                     &         / (1.0_jprb - this%sw_albedo_diffuse(:,iwall,jcol)))
                ! Add diffuse into wall
                surface_flux%sw_dn_facet(jcol,iwall) = surface_flux%sw_dn_direct_facet(jcol,iwall) &
                     &   + wall_scaling * sum(sw_abs_g * this%sw_wall_abs_frac_dif(:,iregion,jcol) &
                     &         / (1.0_jprb - this%sw_albedo_diffuse(:,iwall,jcol)))
!!$
!!$                sw_direct_abs_g  = flux%sw_dn_direct_surf_g (:,jcol) - sw_dn_direct_g
!!$                sw_diffuse_abs_g = flux%sw_dn_diffuse_surf_g(:,jcol)*(1.0_jprb - this%sw_albedo_diffuse_reg(:,iregion,jcol)) &
!!$                     &   + flux%sw_dn_direct_surf_g (:,jcol)*(1.0_jprb - this%sw_albedo_direct_reg (:,iregion,jcol)) &
!!$                     &   - sw_dn_direct_g - sw_dn_diffuse_g + sw_up_g
!!$
!!$
!!$                surface_flux%sw_abs_canopy(jcol,jtile) = sum((1.0_jprb - this%sw_wall_abs_frac(:,iregion,jcol))*sw_abs_g)
!!$                surface_flux%sw_dn_facet(jcol,iwall) = sum(this%sw_wall_abs_frac(:,iregion,jcol)*sw_abs_g &
!!$                     &                           / (1.0_jprb - this%sw_albedo_diffuse(:,iwall,jcol)))
!!$                surface_flux%sw_up_facet(jcol,iwall) = sum(this%sw_albedo_diffuse(:,iwall,jcol) &
!!$                     &                             *this%sw_wall_abs_frac(:,iregion,jcol)*sw_abs_g &
!!$                     &                           / (1.0_jprb - this%sw_albedo_diffuse(:,iwall,jcol)))
              end if

              if (config%do_lw) then
                surface_flux%lw_dn_facet(jcol,iroof) = flux%lw_dn(jcol,isurf)
                surface_flux%lw_up_facet(jcol,iroof) &
                     &  = sum(this%lw_emissivity(:,iroof,jcol)*this%planck_facet(:,iroof,jcol) &
                     &  +  (1.0_jprb - this%lw_emissivity(:,iroof,jcol))*flux%lw_dn_surf_canopy(:,jcol))
                lw_dn_g = (this%lw_transmittance(:,iregion,jcol)*flux%lw_dn_surf_canopy(:,jcol) &
                     &    +this%lw_reflectance  (:,iregion,jcol) &
                     &     *this%lw_emissivity(:,iground,jcol)*this%planck_facet(:,iground,jcol) &
                     &    +this%lw_source(:,iregion,jcol)) &
                     &  / (1.0_jprb - (1.0_jprb - this%lw_emissivity(:,iground,jcol)) &
                     &               *this%lw_reflectance(:,iregion,jcol))
                lw_up_g = (1.0_jprb-this%lw_emissivity(:,iground,jcol))*lw_dn_g &
                     &       +this%lw_emissivity(:,iground,jcol)*this%planck_facet(:,iground,jcol)
                surface_flux%lw_dn_facet(jcol,iground) = sum(lw_dn_g)
                surface_flux%lw_up_facet(jcol,iground) = sum(lw_up_g)
!                lw_abs_g = flux%lw_dn(jcol,isurf)-flux%lw_up(jcol,isurf) &
!                     &  - surface_flux%lw_dn_facet(jcol,iground)+surface_flux%lw_up_facet(jcol,iground)
!                lw_abs_g = flux%lw_dn_surf_canopy(:,jcol) &
!                     &  * (1.0_jprb-this%lw_emissivity_region(:,iregion,jcol)) &
!                     &  - this%lw_emission_region(:,iregion,jcol) &
!                     &  - lw_dn_g*(1.0_jprb-this%lw_emissivity(:,iground,jcol)) &
!                     &  + this%lw_emissivity(:,iground,jcol)*this%planck_facet(:,iground,jcol)
!                lw_abs_g = flux%lw_dn_surf_canopy(:,jcol) &
!                     &  * this%lw_emissivity_region(:,iregion,jcol) &
!                     &  - this%lw_emission_region(:,iregion,jcol) &
!                     &  - lw_dn_g*this%lw_emissivity(:,iground,jcol) &
!                     &  + this%lw_emissivity(:,iground,jcol)*this%planck_facet(:,iground,jcol)
                lw_abs_g = (flux%lw_dn_surf_canopy(:,jcol) + lw_up_g) &
                     &   * (1.0_jprb - this%lw_reflectance(:,iregion,jcol) &
                     &               - this%lw_transmittance(:,iregion,jcol)) &
                     &   + this%lw_total_wall_emission(:,iregion,jcol) &
                     &   + this%lw_total_canopy_emission(:,iregion,jcol) &
                     &   - 2.0_jprb * this%lw_source(:,iregion,jcol)

!                surface_flux%lw_up_facet(jcol,iwall) = sum(this%planck_facet(:,iwall,jcol) &
!                     &  + wall_scaling*(1.0_jprb-this%lw_emissivity(:,iwall,jcol)) &
!                     &    *this%lw_wall_abs_frac(:,iregion,jcol)*lw_abs_g &
!                     &    / this%lw_emissivity(:,iwall,jcol))
!                surface_flux%lw_dn_facet(jcol,iwall) = sum((this%planck_facet(:,iwall,jcol) &
!                     &  + this%lw_total_wall_emission(:,iregion,jcol)) / this%lw_emissivity(:,iwall,jcol))

!                surface_flux%lw_dn_facet(jcol,iwall) = wall_scaling*sum((this%lw_total_wall_emission(:,iregion,jcol) &
!                     &  + this%lw_wall_abs_frac(:,iregion,jcol)*lw_abs_g) &
!                     &  / this%lw_emissivity(:,iwall,jcol))
!                surface_flux%lw_up_facet(jcol,iwall) = surface_flux%lw_dn_facet(jcol,iwall) &
!                     &  - wall_scaling*sum(this%lw_wall_abs_frac(:,iregion,jcol)*lw_abs_g)

                surface_flux%lw_dn_facet(jcol,iwall) &
                     &  = wall_scaling*sum(this%lw_wall_abs_frac(:,iregion,jcol)*lw_abs_g &
                     &                     / this%lw_emissivity(:,iwall,jcol))
                surface_flux%lw_up_facet(jcol,iwall) = surface_flux%lw_dn_facet(jcol,iwall) &
                     &  + wall_scaling * sum(this%lw_total_wall_emission(:,iregion,jcol) &
                     &                       - this%lw_wall_abs_frac(:,iregion,jcol)*lw_abs_g)

                surface_flux%lw_abs_canopy(jcol,jtile) &
                     &  = sum(lw_abs_g*(1.0_jprb-this%lw_wall_abs_frac(:,iregion,jcol)) &
                     &        -this%lw_total_canopy_emission(:,iregion,jcol))
              end if
            end select

          end if

        end do
      end do
    end if

    if (lhook) call dr_hook('radsurf_intermediate:partition_fluxes',1,hook_handle)

  end subroutine partition_fluxes

end module radsurf_intermediate
