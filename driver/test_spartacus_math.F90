! (C) Copyright 2014- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
!
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.

! This simple program tests the functioning of the matrix operations
! in the radiation_matrix module

#define DO_TWO_REGIONS 1

program test_spartacus_math

  use parkind1
  use radiation_matrix
#ifdef DO_TWO_REGIONS
  use tcrad_2region_solver, only : expm_tridiagonal, inv_tridiagonal
#else
  use tcrad_3region_solver, only : expm_tridiagonal, inv_tridiagonal
#endif

  implicit none

  integer, parameter :: n = 10

  integer, parameter :: m = 2

  real(jprb), dimension(n,m,m) :: A, Asave, B, C, expA
  real(jprb), dimension(n,m) :: v, w

  ! Terms in exchange matrix
  real(jprb), dimension(n) :: aa, bb, cc, dd

  integer :: j

  A = 0.0_jprb
  A(:,1,1) = -2.0_jprb
  A(:,2,2) = -3.0_jprb
  A(:,2,1) = 0.5_jprb
  A(:,1,2) = 0.5_jprb

  B = 2.0_jprb
  B(:,1,1) = -1.4_jprb
  B(:,2,1) = 1.0e4_jprb

  v(:,1) = 2.0_jprb
  v(:,2) = 3.0_jprb
  if (m == 3) then
     v(:,m) = -1.5_jprb
     A(:,m,m) = -6.0
     A(:,2,3) = 0.5
     A(:,3,2) = 0.5
  end if

  Asave = A

  write(*,*) 'A ='
  do j = 1,m
     write(*,*) A(1,j,:)
  end do

  write(*,*) 'B ='
  do j = 1,m
     write(*,*) B(1,j,:)
  end do

  write(*,*) 'v = ', v(1,:)

  C = mat_x_mat(n,n,m,A,B)
  w = mat_x_vec(n,n,m,A,v)

  write(*,*) 'C = A*B ='
  do j = 1,m
     write(*,*) C(1,j,:)
  end do
  write(*,*) 'w = A*v = ', w(1,:)

  w = solve_vec(n,n,m,A,v)
  C = solve_mat(n,n,m,A,B)


  write(*,*) 'C = A\B ='
  do j = 1,m
     write(*,*) C(1,j,:)
  end do
  write(*,*) 'w = A\v = ', w(1,:)

  call expm(n,n,m,A,IMatrixPatternDense)

  write(*,*) 'expm(A) ='
  do j = 1,m
     write(*,*) A(1,j,:)
  end do

  ! Test fast_expm_exchange_3
  aa = 2.0_jprb
  bb = 3.0_jprb
  cc = 5.0_jprb
  dd = 7.0_jprb

  if (m == 2) then

    A=Asave

    call inv_tridiagonal(n,A,expA)

    write(*,*) 'inv(A) = '
    do j = 1,m
      write(*,*) expA(1,j,:)
    end do

    A=Asave

    call expm_tridiagonal(n,A,expA)

    write(*,*) 'expm_tridiagonal(A) = '
    do j = 1,m
      write(*,*) expA(1,j,:)
    end do

    call fast_expm_exchange(n,n,aa,bb,A)
    
    write(*,*) 'fast_expm(A) = '
    do j = 1,m
      write(*,*) A(1,j,:)
    end do

    A = 0.0_jprb
    A(:,1,1) = -aa
    A(:,2,1) = aa
    A(:,1,2) = bb
    A(:,2,2) = -bb
    
    A = Asave

    call expm(n,n,m,A,IMatrixPatternDense)
    
    write(*,*) 'expm(A) = '
    do j = 1,m
      write(*,*) A(1,j,:)
    end do

    ! Test zeros lead to identity matrix
    aa = 0.0_jprb
    call fast_expm_exchange(n,n,aa,aa,A)
    
    write(*,*) 'expm(zeros) = '
    do j = 1,m
      write(*,*) A(1,j,:)
    end do

  else 
    
    call fast_expm_exchange(n,n,aa,bb,cc,dd,A)
    
    write(*,*) 'fast_expm(A) = '
    do j = 1,m
      write(*,*) A(1,j,:)
    end do

    A = 0.0_jprb
    A(:,1,1) = -aa
    A(:,2,1) = aa
    A(:,1,2) = bb
    A(:,2,2) = -bb-cc
    if (m == 3) then
      A(:,3,2) = cc
      A(:,2,3) = dd
      A(:,3,3) = -dd
    end if
    
    call expm(n,n,m,A,IMatrixPatternDense)
    
    write(*,*) 'expm(A) = '
    do j = 1,m
      write(*,*) A(1,j,:)
    end do

    ! Test zeros lead to identity matrix
    aa = 0.0_jprb
    call fast_expm_exchange(n,n,aa,aa,aa,aa,A)
    
    write(*,*) 'expm(zeros) = '
    do j = 1,m
      write(*,*) A(1,j,:)
    end do

    A = Asave

!    A(:,1,1) = -1000_jprb
!    A(:,2,2) = -2000_jprb
!    A(:,3,3) = -3000_jprb

    if (m == 3) then
      A(:,3,1) = 0.0_jprb
      A(:,1,3) = 0.0_jprb
    end if

    write(*,*) 'A = '
    do j = 1,m
      write(*,*) A(1,j,:)
    end do

    call expm_tridiagonal(n,A,expA)

    write(*,*) 'expm3x3tridiag(A) = '
    do j = 1,m
      write(*,*) expA(1,j,:)
    end do

    call inv_tridiagonal(n,A,expA)

    write(*,*) 'inv(A) = '
    do j = 1,m
      write(*,*) expA(1,j,:)
    end do

    call expm(n,n,m,A,IMatrixPatternDense)
    
    write(*,*) 'expm(A) = '
    do j = 1,m
      write(*,*) A(1,j,:)
    end do

  end if
    
end program test_spartacus_math
