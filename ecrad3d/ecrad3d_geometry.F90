module ecrad3d_geometry
  
  use parkind1, only : jprb

  public

  type geometry_type

    
    real(jprb), allocatable :: x(:), y(:)
    
    real(jprb), allocatable :: dz

    ! Horizontal area in m2 of each grid point
    !real(jprb), allocatable :: area(:)

    ! Number of columns
    integer :: ncol
    
    ! Number of x and y grid points (could be lontigude and latitude)
    ! if the horizontal grid is rectangular, must multiply together to
    ! get ncol
    integer :: nx, ny
    
    ! Is the horizontal arrangement of columns on a rectangular grid
    ! such that a point can be identified by an x and a y index?
    logical :: is_rectangular_grid = .true.

    ! Is t
    logical :: is_lat_lon
    
  contains

    procedure :: advect
    procedure :: diffuse
    
  end type geometry_type

contains

  subroutine advect(self, ncol, nspec, ilay, field_in, field_out)

    implicit none

    class(geometry_type), intent(in) :: self
    integer, intent(in) :: ncol, nspec
    integer, intent(in) :: ilay
    real(jprb), intent(in)  :: field_in(ncol,nspec)
    real(jprb), intent(out) :: field_out(ncol,nspec)

    !field_out = field_in

  end subroutine advect

  subroutine diffuse(self, ncol, nspec, ilay, field_in, field_out)

    implicit none
    
    class(geometry_type), intent(in) :: self
    integer, intent(in) :: ncol, nspec
    integer, intent(in) :: ilay
    real(jprb), intent(in)  :: field_in(ncol,nspec)
    real(jprb), intent(out) :: field_out(ncol,nspec)

    field_out = field_in

  end subroutine diffuse

end module ecrad3d_geometry
