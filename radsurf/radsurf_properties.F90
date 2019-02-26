! radsurf_properties.f90 - Derived type for surface properties
!
! Copyright (C) 2017-2018 ECMWF
!
! Author:  Robin Hogan
! Email:   r.j.hogan@ecmwf.int
! License: see the COPYING file for details
!

module radsurf_properties

  use parkind1, only : jpim, jprb

  implicit none

  ! Number of tile types
  integer(kind=jpim), parameter :: NTileTypes = 3

  ! Codes for the different type of tile
  enum, bind(c)
    enumerator :: ITileFlat = 1, &
         &        ITileVegetation, &
         &        ITileUrban3D
  end enum

  character(len=*), parameter :: TileRepresentationName(NTileTypes) &
       &  = (/ 'Flat                 ', &
       &       'HomogeneousVegetation', &
       &       'Urban3D              ' /)

  ! Number of surface facets and regions 
  integer(kind=jpim), parameter :: NTileFacets (NTileTypes) = (/ 1, 1, 3 /)
  integer(kind=jpim), parameter :: NTileRegions(NTileTypes) = (/ 0, 1, 1 /)


  !---------------------------------------------------------------------
  ! Derived type storing a physical description of the properties of
  ! the surface tiles needed to compute the boundary conditions for
  ! the radiation scheme.
  type surface_type
    ! Skin temperature (ncol,nfacet)
    real(kind=jprb), allocatable :: skin_temperature(:,:)

    ! Air temperature in canopy (ncol,ntile)
    real(kind=jprb), allocatable :: canopy_temperature(:,:)

    ! Shortwave albedo: if sw_albedo_direct is not allocated then
    ! sw_albedo will be used for both direct and diffuse solar
    ! radiation; otherwise, sw_albedo will be used for diffuse
    ! radiation and sw_albedo_direct for direct
    ! radiation. Dimensioning is (ncol,nalbedobands,nfacet).
    real(kind=jprb), allocatable, dimension(:,:,:) :: &
         &  sw_albedo, sw_albedo_direct

    ! Longwave emissivity: dimensioning is (ncol,nemissbands,nfacet)
    real(kind=jprb), allocatable, dimension(:,:,:) :: &
         &  lw_emissivity

    ! Representation codes (ITileFlat etc) for each tile: dimensioning
    ! is (ntile)
    integer(kind=jpim), allocatable, dimension(:) :: &
         &  i_representation

    ! Indices to the facets: dimensioning is (ncol,ntile)
    integer(kind=jpim), allocatable, dimension(:,:) :: &
         &  i_ground_facet, &
         &  i_roof_facet, &
         &  i_wall_facet

    ! Indices to the regions: dimensioning is (ncol,ntile)
    integer(kind=jpim), allocatable, dimension(:,:) :: &
         &  i_region_1 !, i_region_2

    ! Physical properties dimensioned (ncol,ntile)
    real(kind=jprb), allocatable, dimension(:,:) :: &
         &  tile_fraction, & ! Fraction of column occupied by each tile
         &  canopy_depth     ! Depth of canopy (m)

    ! Building properties in urban tile dimensioned (ncol,ntile)
    real(kind=jprb), allocatable, dimension(:,:) :: &
         &  building_fraction, & ! Fraction of canopy containing buildings
         &  building_normalized_perimeter ! Perimeter length of buildings per unit area (m)

    ! Vegetation physical properties dimensioned (ncol,ntile)
    real(kind=jprb), allocatable, dimension(:,:) :: &
         &  vegetation_optical_depth, & ! Same as LAI?
         &  vegetation_fractional_std   ! Fractional standard deviation

    ! Vegetation optical properties
    real(kind=jprb), allocatable, dimension(:,:,:) :: &
         &  vegetation_sw_albedo, & ! (ncol,nalbedobands,ntile)
         &  vegetation_lw_emissivity! (ncol,nemissbands,ntile)

! When 3D effects in vegetation canopies are represented, these will
! be available
!         &  vegetation_normalized_perimeter
!         &  vegetation_fraction, &

    ! Number of tiles, columns, facets, regions
    integer(kind=jpim) :: ntile, ncol, nfacet, nregion

    ! Number of albedo and emissivity bands
    integer(kind=jpim) :: nalbedobands, nemissbands

    ! Is this a simple surface containing one "flat" tile?
    logical :: is_simple = .true.

  contains
    procedure :: allocate   => allocate_surface
    procedure :: deallocate => deallocate_surface
    procedure :: read => read_from_netcdf
    procedure :: set_facet_indices

  end type surface_type


contains

  !---------------------------------------------------------------------
  subroutine allocate_surface(this, ncol, nalbedobands, nemissbands, &
       &                      i_representation, &
       &                      use_sw_albedo_direct)

    use yomhook,        only : lhook, dr_hook

    class(surface_type), intent(inout) :: this
    integer(kind=jpim),  intent(in)    :: ncol, nalbedobands, nemissbands
    integer(kind=jpim),  intent(in), optional :: i_representation(:)
    logical(kind=jpim),  intent(in), optional :: use_sw_albedo_direct

    ! Number of facets and regions
    integer(kind=jpim) :: nfacet, nregion, ntile

    real(jprb) :: hook_handle

    if (lhook) call dr_hook('radiation_surface:allocate',0,hook_handle)

    call this%deallocate()

    this%ncol   = ncol
    if (present(i_representation)) then
      nfacet    = sum(NTileFacets (i_representation))
      nregion   = sum(NTileRegions(i_representation))
      ntile     = size(i_representation)
      this%is_simple = .false.
    else
      nfacet    = 1
      nregion   = 0
      ntile     = 1
      this%is_simple = .true.
    end if
    this%nfacet = nfacet
    this%nregion= nregion
    this%ntile  = ntile
    this%nalbedobands = nalbedobands
    this%nemissbands  = nemissbands

    allocate(this%skin_temperature(ncol,nfacet))
    allocate(this%sw_albedo(ncol, nalbedobands, nfacet))

    if (present(use_sw_albedo_direct)) then
      if (use_sw_albedo_direct) then
        allocate(this%sw_albedo_direct(ncol, nalbedobands, nfacet))
      end if
    end if

    allocate(this%lw_emissivity(ncol, nemissbands, nfacet))

    if (this%ntile > 1) then
      allocate(this%tile_fraction(ncol,ntile))
    end if

    if (nregion > 0) then
      allocate(this%canopy_depth(ncol,ntile))
      allocate(this%canopy_temperature(ncol,ntile))
      this%canopy_depth = 0.0_jprb
      this%canopy_temperature = -1.0_jprb
    end if

    if (.not. this%is_simple) then
      allocate(this%i_representation(ntile))
      this%i_representation = i_representation

      if (any(i_representation == ITileUrban3D)) then
        allocate(this%building_fraction(ncol,ntile))
        allocate(this%building_normalized_perimeter(ncol,ntile))
        ! Initialize to default values
        this%building_fraction = 0.0_jprb
        this%building_normalized_perimeter = 0.0_jprb
      end if
      
      if (any(i_representation == ITileVegetation)) then
        allocate(this%vegetation_optical_depth(ncol,ntile))
        allocate(this%vegetation_fractional_std(ncol,ntile))
        allocate(this%vegetation_sw_albedo(ncol,nalbedobands,ntile))
        allocate(this%vegetation_lw_emissivity(ncol,nemissbands,ntile))
        ! Initialize to default values
        this%vegetation_optical_depth = 0.0_jprb
        this%vegetation_fractional_std = 0.0_jprb
        this%vegetation_sw_albedo = 0.0_jprb
        this%vegetation_lw_emissivity = 1.0_jprb
      end if
      
      call this%set_facet_indices 

    end if

    if (lhook) call dr_hook('radiation_surface:allocate',1,hook_handle)

  end subroutine allocate_surface


  !---------------------------------------------------------------------
  ! Set the indices to facets
  subroutine set_facet_indices(this)

    use radiation_io,   only : nulerr, radiation_abort

    class(surface_type), intent(inout) :: this

    ! Indices to tiles and facets
    integer(kind=jpim) :: ifacet, iregion, jtile
    if (.not. this%is_simple) then

      allocate(this%i_ground_facet(this%ncol,this%ntile))
      this%i_ground_facet = 0

      if (any(this%i_representation == ITileUrban3D)) then
        allocate(this%i_roof_facet(this%ncol,this%ntile))
        allocate(this%i_wall_facet(this%ncol,this%ntile))
        ! Initialize to default values
        this%i_roof_facet = 0
        this%i_wall_facet = 0
      end if
      
      if (this%nregion > 0) then
        allocate(this%i_region_1(this%ncol,this%ntile))
        this%i_region_1 = 0
      end if

      ! Set the indices to the various facets
      ifacet = 1
      iregion = 1
      do jtile = 1,this%ntile
        this%i_ground_facet(:,jtile) = ifacet
        ifacet = ifacet + 1
        if (this%i_representation(jtile) == ITileVegetation) then
          this%i_region_1(:,jtile) = iregion
          iregion = iregion+1
        else if (this%i_representation(jtile) == ITileUrban3D) then
          this%i_roof_facet(:,jtile) = ifacet
          this%i_wall_facet(:,jtile) = ifacet+1
          ifacet = ifacet+2
          this%i_region_1(:,jtile) = iregion
          iregion = iregion+1
        else if (this%i_representation(jtile) /= ITileFlat) then
          write(nulerr,'(a,i0,a)') '*** Error: tile representation ', &
               &  this%i_representation(jtile), ' not understood'
          call radiation_abort()
        end if
      end do
    end if
  end subroutine set_facet_indices

  !---------------------------------------------------------------------
  subroutine deallocate_surface(this)

    use yomhook, only : lhook, dr_hook

    class(surface_type), intent(inout) :: this

    real(jprb) :: hook_handle

    if (lhook) call dr_hook('radiation_surface:deallocate',0,hook_handle)

    if (allocated(this%skin_temperature)) then
      deallocate(this%skin_temperature)
    end if
    if (allocated(this%sw_albedo)) then
      deallocate(this%sw_albedo)
    end if
    if (allocated(this%sw_albedo_direct)) then
      deallocate(this%sw_albedo_direct)
    end if
    if (allocated(this%lw_emissivity)) then
      deallocate(this%lw_emissivity)
    end if
    if (allocated(this%i_representation)) then
      deallocate(this%i_representation)
    end if
    if (allocated(this%i_ground_facet)) then
      deallocate(this%i_ground_facet)
    end if
    if (allocated(this%i_roof_facet)) then
      deallocate(this%i_roof_facet)
    end if
    if (allocated(this%i_wall_facet)) then
      deallocate(this%i_wall_facet)
    end if
    if (allocated(this%tile_fraction)) then
      deallocate(this%tile_fraction)
    end if
    if (allocated(this%canopy_depth)) then
      deallocate(this%canopy_depth)
    end if
    if (allocated(this%canopy_temperature)) then
      deallocate(this%canopy_temperature)
    end if
    if (allocated(this%building_fraction)) then
      deallocate(this%building_fraction)
    end if
    if (allocated(this%building_normalized_perimeter)) then
      deallocate(this%building_normalized_perimeter)
    end if
    if (allocated(this%vegetation_optical_depth)) then
      deallocate(this%vegetation_optical_depth)
    end if
    if (allocated(this%vegetation_fractional_std)) then
      deallocate(this%vegetation_fractional_std)
    end if
    if (allocated(this%vegetation_sw_albedo)) then
      deallocate(this%vegetation_sw_albedo)
    end if
    if (allocated(this%vegetation_lw_emissivity)) then
      deallocate(this%vegetation_lw_emissivity)
    end if

    this%ntile  = 0
    this%ncol   = 0
    this%nfacet = 0

    if (lhook) call dr_hook('radiation_surface:deallocate',1,hook_handle)

  end subroutine deallocate_surface


  !---------------------------------------------------------------------
  ! Print a description of the surface tile types
  subroutine print_surface_representation(i_representation)

    use radiation_io, only : nulout

    integer(kind=jpim), dimension(:), allocatable, intent(in) :: i_representation

    integer :: ntile, jtile, ifacet, iregion

    write(nulout,'(a)') 'Surface tile representation:'
    if (.not. allocated(i_representation)) then
      write(nulout,'(a)') '  Simple (one flat tile)'
    else
      ntile = size(i_representation,1)
      ifacet  = 1
      iregion = 1;
      do jtile = 1,ntile
        write(nulout,'(a,i0,a,a)') '  Tile ', jtile, ': ', trim(TileRepresentationName(i_representation(jtile)))
        write(nulout,'(a,i0,a)') '    Facet ', ifacet, ': ground'
        ifacet = ifacet+1

        select case (i_representation(jtile))

        case (ITileVegetation)
          write(nulout,'(a,i0,a)') '    Region ', iregion, ': vegetation canopy'
          iregion = iregion+1

        case (ITileUrban3D)
          write(nulout,'(a,i0,a)') '    Facet ', ifacet, ': roof'
          write(nulout,'(a,i0,a)') '    Facet ', ifacet+1, ': wall'
          ifacet = ifacet+1
          write(nulout,'(a,i0,a)') '    Region ', iregion, ': street canyon'
          iregion = iregion+1

        end select

      end do
    end if
    
  end subroutine print_surface_representation


  !---------------------------------------------------------------------
  subroutine read_from_netcdf(this, file)

    use parkind1,           only : jprb, jpim
    use easy_netcdf,        only : netcdf_file
    
    implicit none

    type(netcdf_file),  intent(in)     :: file
    class(surface_type), intent(inout) :: this

    real(kind=jprb), allocatable       :: data_1d(:)

    integer(kind=jpim), parameter :: ipermute(3) = [2,3,1]

    call this%deallocate

    this%is_simple = .false.

    call file%get('skin_temperature', this%skin_temperature, do_transp=.true.)
    call file%get('canopy_temperature', this%canopy_temperature, do_transp=.true.)
    call file%get('sw_albedo', this%sw_albedo, ipermute=ipermute)
    if (file%exists('sw_albedo_direct')) then
      call file%get('sw_albedo', this%sw_albedo_direct, ipermute=ipermute)
    end if
    call file%get('lw_emissivity', this%lw_emissivity, ipermute=ipermute)
    call file%get('tile_representation', data_1d)
    this%ntile = size(data_1d)
    allocate(this%i_representation(this%ntile))
    this%i_representation = int(data_1d)
    call file%get('tile_fraction', this%tile_fraction, do_transp=.true.)
    call file%get('canopy_depth', this%canopy_depth, do_transp=.true.)
    call file%get('building_fraction', this%building_fraction, do_transp=.true.)
    if (file%exists('building_normalized_perimeter')) then
      call file%get('building_normalized_perimeter', &
           &  this%building_normalized_perimeter, do_transp=.true.)
    else
      ! Convert building scale to normalized perimeter
      call file%get('building_scale', this%building_normalized_perimeter, do_transp=.true.)
      this%building_normalized_perimeter = 4.0_jprb * this%building_fraction &
           &  * (1.0_jprb-this%building_fraction) / max(1.0e-8_jprb,this%building_normalized_perimeter)
    end if
    call file%get('vegetation_optical_depth', this%vegetation_optical_depth, do_transp=.true.)
    call file%get('vegetation_sw_albedo', this%vegetation_sw_albedo, ipermute=ipermute)
    call file%get('vegetation_lw_emissivity', this%vegetation_lw_emissivity, ipermute=ipermute)

    this%nregion = sum(NTileRegions(this%i_representation))
    this%ncol    = size(this%skin_temperature,1)
    this%nfacet  = size(this%skin_temperature,2)

    this%nalbedobands = size(this%sw_albedo,2)
    this%nemissbands = size(this%lw_emissivity,2)

    call this%set_facet_indices

  end subroutine read_from_netcdf


end module radsurf_properties
