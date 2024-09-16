module ecrad3d_geometry
  
  use parkind1, only : jprb
  use ecrad3d_config, only : IBoundaryPeriodic
  
  implicit none
  
  public

  type geometry_type

    ! These variables should be allocated and set by an external
    ! routine
    
    ! Horizontal coordinate variables, either distance in metres or
    ! (if is_lat_lon==.true.) longitude and latitude in degrees
    real(jprb), allocatable :: x(:), y(:)

    ! Height of half-levels above mean sea level or a geoid reference,
    ! not above the surface, noting that the first is ignored (height
    ! of top-of-atmosphere may be infinite), dimensioned (ncol,nz+1)
    ! in metres. This is only used if is_3d is true.
    real(jprb), allocatable :: height_hl(:,:)

    ! Horizontal area in m2 of each grid point
    !real(jprb), allocatable :: area(:)

    ! Number of columns
    integer :: ncol
    
    ! Number of x and y grid points (could be lontigude and latitude)
    ! if the horizontal grid is rectangular, must multiply together to
    ! get ncol
    integer :: nx, ny, nz

    integer :: boundary_conditions = IBoundaryPeriodic
    
    ! Is the horizontal arrangement of columns on a rectangular grid
    ! such that a point can be identified by an x and a y index?
    logical :: is_rectangular_grid = .true.

    ! Are the x and y vectors actually longitude and latitude
    ! (respectively) in degrees?
    logical :: is_lat_lon

    logical :: is_3d = .false.

    ! These variables should be computed from those above by the
    ! "consolidate" type-bound routine, rather than being set by an
    ! external routine

    ! Approximate version of 3D radiative transfer, approprimate for
    ! limited domains, assumes that the layer spacing of each column
    ! is the same, and hence only the average layer thickness need be used.
    real(jprb), allocatable :: dz(:) ! (nlev), metres

    ! Size of each box in metres, noting that for a lon-lat grid the x
    ! is not constant with y so needs to be defined separately for
    ! each column
    real(jprb), allocatable :: dx(:) ! (ncol)
    real(jprb), allocatable :: dy(:) ! (ny)

    ! Indices to the points adjacent to a target point, assuming the
    ! x-y grid is east-north aligned.
    integer, allocatable :: inorth(:), isouth(:), ieast(:), iwest(:) ! (ncol)
        
  contains

    procedure :: consolidate
    procedure :: advect
    procedure :: diffuse
    
  end type geometry_type

contains

  subroutine consolidate(self)

    use radiation_constants, only : Pi, PlanetRadius
    
    class(geometry_type), intent(inout) :: self
    integer :: jlev, jx, jy, icol
    real(jprb) :: mean_height_base, mean_height_top

    real(jprb), parameter :: lat_to_m = PlanetRadius * 2.0_jprb * Pi / 360.0_jprb

    real(jprb) :: lon_to_m

    integer :: noffs
    
    if (self%is_3d) then
      
      ! Compute average layer thickness
      allocate(self%dz(self%nz))    
      mean_height_base = sum(self%height_hl(:,self%nz+1))/self%ncol
      do jlev = self%nz,1,-1
        mean_height_top = sum(self%height_hl(:,jlev))/self%ncol
        self%dz(jlev) = mean_height_top-mean_height_base
        mean_height_base = mean_height_top
      end do
      
      allocate(self%dx(self%ncol))
      allocate(self%dy(self%ny))

      allocate(self%inorth(self%ncol))
      allocate(self%isouth(self%ncol))
      allocate(self%ieast (self%ncol))
      allocate(self%iwest (self%ncol))
      
      if (self%is_lat_lon) then

        self%dy(1) = abs(self%y(2)-self%y(1))*lat_to_m
        do jy = 2,self%ny-1
          self%dy(jy) = 0.5_jprb * abs(self%y(jy+1)-self%y(jy-1))*lat_to_m
        end do
        self%dy(self%ny) = abs(self%y(self%ny)-self%y(self%ny-1))*lat_to_m
    
        do jy = 1,self%ny
          lon_to_m = lat_to_m * cos(self%y(jy)*Pi/180.0_jprb)
          icol = (jy-1)*self%nx + 1
          self%dx(icol) = lon_to_m * abs(self%x(2)-self%x(1))
          do jx = 2,self%nx-1
            icol = (jy-1)*self%nx + jx
            self%dx(icol) = lon_to_m * 0.5_jprb * abs(self%x(jx+1)-self%x(jx-1))
          end do
          icol = (jy-1)*self%nx + self%nx
          self%dx(icol) = lon_to_m * abs(self%x(self%nx)-self%x(self%nx-1))
        end do

      else
        
        self%dy(1) = abs(self%y(2)-self%y(1))
        do jy = 2,self%ny-1
          self%dy(jy) = 0.5_jprb * abs(self%y(jy+1)-self%y(jy-1))
        end do
        self%dy(self%ny) = abs(self%y(self%ny)-self%y(self%ny-1))

        ! Do the first row...
        self%dx(1) = abs(self%x(2)-self%x(1))
        do jx = 2,self%nx-1
          self%dx(jx) = 0.5_jprb * abs(self%x(jx+1)-self%x(jx-1))
        end do
        self%dx(self%nx) = abs(self%x(self%nx)-self%x(self%nx-1))
        ! ...and repeat
        do jy = 2,self%ny
          do jx = 1,self%nx
            icol = (jy-1)*self%nx + jx
            self%dx(icol) = self%dx(jx)
          end do
        end do
      end if
      
    end if

    ! Assign indices of points to the immediate east and west of each target point
    do jy = 1,self%ny
      icol = (jy-1)*self%nx + 1
      self%ieast(icol) = icol + 1
      do jx = 2,self%nx-1
        icol = (jy-1)*self%nx + jx
        self%ieast(icol) = icol + 1
        self%iwest(icol) = icol - 1
      end do
      icol = jy*self%nx
      self%iwest(icol) = icol - 1
    end do

    ! ...and likewise north and south
    do jy = 2,self%ny-1
      do jx = 1,self%nx
        icol = (jy-1)*self%nx + jx
        self%inorth(icol) = icol + self%nx
        self%isouth(icol) = icol - self%nx
      end do
      icol = jy*self%nx
      self%iwest(icol) = icol - 1
    end do
    noffs = (self%ny-1) * self%nx
    do jx = 1,self%nx
      self%inorth(jx) = self%nx + jx
      icol = noffs + jx
      self%isouth(icol) = icol - self%nx
    end do
    
    ! ...and final points across the boundary
    if (self%boundary_conditions == IBoundaryPeriodic) then
      do jy = 1,self%ny
        icol = (jy-1)*self%nx + 1
        self%iwest(icol) = icol + self%nx - 1
        icol = jy*self%nx
        self%ieast(icol) = icol - self%nx + 1
      end do
      do jx = 1,self%nx
        self%isouth(jx) = jx + noffs
        icol = noffs + jx
        self%inorth(icol) = icol - noffs
      end do
    else
      error stop 'Boundary conditions not understood'
    end if

  end subroutine consolidate
    
  
  subroutine advect(self, ncol, nspec, ilay, cos_zenith_angle, azimuth_angle, &
       &            field_in, field_out)

    implicit none

    class(geometry_type), intent(in) :: self
    integer, intent(in) :: ncol, nspec
    integer, intent(in) :: ilay
    real(jprb), intent(in)  :: cos_zenith_angle(ncol)
    real(jprb), intent(in)  :: azimuth_angle(ncol) ! radians
    real(jprb), intent(in)  :: field_in(ncol,nspec)
    real(jprb), intent(out) :: field_out(ncol,nspec)

    ! Displacement in gridboxes
    real(jprb) :: xoffs(ncol), yoffs(ncol)

    real(jprb) :: dz_tan_zenith_angle
    
    integer :: icol, icol11, icol12, icol21, icol22
    integer :: jx, jy, jcol, jspec
    
    integer :: ix1(ncol), ix2(ncol), iy1(ncol), iy2(ncol)
    
    do jy = 1,self%ny
      do jx = 1,self%nx
        icol = (jy-1)*self%nx + jx
        if (cos_zenith_angle(icol) > 0.0_jprb) then
          dz_tan_zenith_angle = self%dz(ilay) &
               &   * sqrt(1.0_jprb-cos_zenith_angle(icol))/cos_zenith_angle(icol)
          xoffs(icol) = -(cos(azimuth_angle(icol)) * dz_tan_zenith_angle) / self%dx(icol)
          yoffs(icol) = -(sin(azimuth_angle(icol)) * dz_tan_zenith_angle) / self%dy(jy)
        else
          xoffs(icol) = 0.0_jprb
          yoffs(icol) = 0.0_jprb
        end if

        ix1(icol) = jx + floor(xoffs(icol))
        ix2(icol) = ix1(icol) + 1
        iy1(icol) = jy + floor(yoffs(icol))
        iy2(icol) = iy1(icol) + 1
        xoffs(icol) = xoffs(icol) - floor(xoffs(icol))
        yoffs(icol) = yoffs(icol) - floor(yoffs(icol))
      end do
    end do
    
    if (self%boundary_conditions == IBoundaryPeriodic) then
      ! do jcol = 1,self%ncol
      !   if (ix1(jcol) < 1) then
      !     ix1(jcol) = ix1(jcol)+self%nx
      !   else if (ix1(jcol) > self%nx) then
      !     ix1(jcol) = ix1(jcol)-self%nx
      !   end if
      !   if (ix2(jcol) < 1) then
      !     ix2(jcol) = ix2(jcol)+self%nx
      !   else if (ix2(jcol) > self%nx) then
      !     ix2(jcol) = ix2(jcol)-self%nx
      !   end if
      !   if (iy1(jcol) < 1) then
      !     iy1(jcol) = iy1(jcol)+self%ny
      !   else if (iy1(jcol) > self%ny) then
      !     iy1(jcol) = iy1(jcol)-self%ny
      !   end if
      !   if (iy2(jcol) < 1) then
      !     iy2(jcol) = iy2(jcol)+self%ny
      !   else if (iy2(jcol) > self%ny) then
      !     iy2(jcol) = iy2(jcol)-self%ny
      !   end if
      ! end do
      do jcol = 1,self%ncol
        ix1(jcol) = modulo(ix1(jcol)-1, self%nx)+1
        ix2(jcol) = modulo(ix2(jcol)-1, self%nx)+1
        iy1(jcol) = modulo(iy1(jcol)-1, self%ny)+1
        iy2(jcol) = modulo(iy2(jcol)-1, self%ny)+1
      end do
    else
      error stop 'BC not understood'
    end if

    do jcol = 1,self%ncol
      icol11 = (iy1(jcol)-1)*self%nx + ix1(jcol)
      icol12 = (iy2(jcol)-1)*self%nx + ix1(jcol)
      icol21 = (iy1(jcol)-1)*self%nx + ix2(jcol)
      icol22 = (iy2(jcol)-1)*self%nx + ix2(jcol)
      !write(101,*) jcol, ix1(jcol), ix2(jcol), iy1(jcol), iy2(jcol), icol11, icol12, icol21, icol22, &
      !     &  xoffs(jcol), yoffs(jcol)
      do jspec = 1,nspec
        field_out(jcol,jspec) &
             &  = (1.0_jprb-xoffs(jcol)) * ( (1.0_jprb-yoffs(jcol))*field_in(icol11,jspec) &
             &                              +          yoffs(jcol) *field_in(icol12,jspec) ) &
             &  +           xoffs(jcol)  * ( (1.0_jprb-yoffs(jcol))*field_in(icol21,jspec) &
             &                              +          yoffs(jcol) *field_in(icol22,jspec) )
!        field_out(jcol,jspec) &
!             &  = (1.0_jprb-xoffs(jcol)) * ( (1.0_jprb-yoffs(jcol))*field_in(icol22,jspec) &
!             &                              +          yoffs(jcol) *field_in(icol21,jspec) ) &
!             &  +           xoffs(jcol)  * ( (1.0_jprb-yoffs(jcol))*field_in(icol12,jspec) &
!             &                              +          yoffs(jcol) *field_in(icol11,jspec) )
      end do
    end do
    !field_out = field_in

  end subroutine advect

  subroutine diffuse(self, ncol, nspec, ilay, field_in, field_out)

    implicit none
    
    class(geometry_type), intent(in) :: self
    integer, intent(in) :: ncol, nspec
    integer, intent(in) :: ilay
    real(jprb), intent(in)  :: field_in(ncol,nspec)
    real(jprb), intent(out) :: field_out(ncol,nspec)

    ! 1-sigma smoothing scale in units of horizontal pixels
    real(jprb) :: scale
    real(jprb) :: frac

    real(jprb) :: field_tmp(ncol)

    ! We find that using a diffusion mean-free path of 2*dz/dx results
    ! in the best agreement with Monte-Carlo calculations
    real(jprb), parameter :: DiffusionFactor = 2.0_jprb

    ! We don't do a full Gaussian convolution for the diffusion but
    ! one or two neighbour exchanges, and must limit how much is exchanged
    integer    :: n_exchange
    real(jprb) :: exchange_limit
    real(jprb) :: exchange_factor
    
    real(jprb), parameter :: TwoThirds = 2.0_jprb / 3.0_jprb
    
    integer :: jcol, jspec, iy

    ! Decide whether to do one or two exchanges, and in each case how
    ! much to exchange and the limit on the exchange
    if (any(self%dz(ilay) > self%dx(1:self%ncol)*(sqrt(TwoThirds)/DiffusionFactor))) then
      n_exchange = 2
      exchange_limit = 4.0_jprb / 5.0_jprb
      ! Reduced exchange per step because we will do two
      exchange_factor = DiffusionFactor / sqrt(2.0_jprb)
    else
      n_exchange = 1
      exchange_limit = TwoThirds
      exchange_factor = DiffusionFactor
    end if
    
    do jspec = 1,nspec
      ! East-west diffusion
      do jcol = 1,self%ncol
        scale = exchange_factor * self%dz(ilay) / self%dx(jcol)
        frac  = min(scale*scale, exchange_limit)
        field_tmp(jcol) = (1.0_jprb - frac) * field_in(jcol,jspec) &
             &          +  0.5_jprb * frac * (field_in(self%ieast(jcol),jspec) + field_in(self%iwest(jcol),jspec))
      end do
      ! North-south diffusion
      do jcol = 1,self%ncol
        iy = (jcol-1)/self%nx + 1
        scale = exchange_factor * self%dz(ilay) / self%dy(iy)
        frac  = min(scale*scale, exchange_limit)
        field_out(jcol,jspec) = (1.0_jprb - frac) * field_tmp(jcol) &
             &                +  0.5_jprb * frac * (field_tmp(self%inorth(jcol)) + field_tmp(self%isouth(jcol)))
      end do

      if (n_exchange == 2) then
        ! East-west diffusion
        do jcol = 1,self%ncol
          scale = exchange_factor * self%dz(ilay) / self%dx(jcol)
          frac  = min(scale*scale, exchange_limit)
          field_tmp(jcol) = (1.0_jprb - frac) * field_out(jcol,jspec) &
               &          +  0.5_jprb * frac * (field_out(self%ieast(jcol),jspec) + field_out(self%iwest(jcol),jspec))
        end do
        ! North-south diffusion
        do jcol = 1,self%ncol
          iy = (jcol-1)/self%nx + 1
          scale = exchange_factor * self%dz(ilay) / self%dy(iy)
          frac  = min(scale*scale, exchange_limit)
          field_out(jcol,jspec) = (1.0_jprb - frac) * field_tmp(jcol) &
               &                +  0.5_jprb * frac * (field_tmp(self%inorth(jcol)) + field_tmp(self%isouth(jcol)))
        end do
      end if
      
    end do

  end subroutine diffuse

end module ecrad3d_geometry
