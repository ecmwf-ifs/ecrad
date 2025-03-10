module ecrad_definitions

  use parkind1, only : jpim

  implicit none

  ! By making the number of regions a parameter NREGION, the compiler
  ! can optimize better.  We also provide a preprocessor parameter as
  ! in one or two cases the operations on 2x2 and 3x3 matrices
  ! (e.g. matrix exponentials) are completely different algorithms.
#define NUM_REGIONS 2
  integer(jpim), parameter :: NREGION = NUM_REGIONS

! Allow size of inner dimension (number of g-points) to be known at compile time if NG_LW is defined
#ifdef NG_LW
  integer(jpim), parameter :: ng_lw = NG_LW
#else
#define ng_lw ng_lw_in
#endif


end module ecrad_definitions