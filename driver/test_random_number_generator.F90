program test_random_number_generator

  ! This program tests the radiation random number generator - the
  ! Minstd generator should produce the same results whether or not
  ! working floating-point precision (jprb) is single precision (4) or
  ! double (8)
  
  use parkind1, only : jprb, jpim, jpib
  use radiation_random_numbers, only : rng_type, IRngMinstdVector, IRngNative

  integer(jpim), parameter :: streammax = 512

  ! Choose between two types of random number generator
  !integer(jpim), parameter :: IRngType = IRngNative
  integer(jpim), parameter :: IRngType = IRngMinstdVector
  
  type(rng_type) :: random_number_generator
  real(jprb)     :: vec(streammax)
  integer :: jl

  print *, 'working_fp_precision = ', jprb

  if (IRngType == IRngNative) then
    print *, "rng_type = 'native'"
  else
    print *, "rng_type = 'native'"
  end if
  
  call random_number_generator%initialize(itype=IRngType, iseed=212075152, &
       &                                  nmaxstreams=streammax)
  do jl = 1,1
    print *, 'initial_state = [ ', int(random_number_generator%istate(1:streammax),jpib), ' ]'
    call random_number_generator%uniform_distribution(vec)
    print *, 'uniform_deviates = [ ', vec, ' ]'
    print *, 'final_state = [ ', int(random_number_generator%istate(1:streammax),jpib), ' ]'
  end do

end program test_random_number_generator
