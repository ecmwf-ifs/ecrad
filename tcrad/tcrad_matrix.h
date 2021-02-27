! tcrad_matrix.h - Batched matrix operations for TCRAD solver -*- f90 -*-
!
! (C) Copyright 2021- ECMWF.
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
! Treat A as an m-by-m square matrix and b as n NREGION-element vectors
! (with the n dimension varying fastest), and perform n matrix-vector
! multiplications
pure function singlemat_x_vec(n,A,b)

  use parkind1, only : jpim, jprb

  implicit none

  integer,    intent(in)                             :: n
  real(jprb), intent(in), dimension(NREGION,NREGION) :: A
  real(jprb), intent(in), dimension(n,NREGION)       :: b
  real(jprb),             dimension(n,NREGION)       :: singlemat_x_vec
  
  integer    :: j1, j2
  
  ! Array-wise assignment
  singlemat_x_vec = 0.0_jprb
  
  do j1 = 1,NREGION
    do j2 = 1,NREGION
      singlemat_x_vec(:,j1) = singlemat_x_vec(:,j1) + A(j1,j2)*b(:,j2)
    end do
  end do
  
end function singlemat_x_vec

#if NUM_REGIONS == 3

!---------------------------------------------------------------------
! Return X = B A^-1 = (A^-T B)^T optimized for 3x3 matrices, where B
! is a diagonal matrix, using LU factorization and substitution with
! no pivoting.
pure subroutine diag_mat_right_divide_3(n,iend,A,B,X)

  use parkind1, only : jpim, jprb

  integer(jpim), intent(in)  :: n, iend
  real(jprb),    intent(in)  :: A(iend,3,3)
  real(jprb),    intent(in)  :: B(iend,3)
  real(jprb),    intent(out) :: X(n,3,3)
  
  real(jprb), dimension(iend) :: L21, L31, L32
  real(jprb), dimension(iend) :: U22, U23, U33
  real(jprb), dimension(iend) :: y2, y3

  !    associate (U11 => A(:,1,1), U12 => A(:,1,2), U13 => A(1,3))
  ! LU decomposition of the *transpose* of A:
  !       ( 1        )   (U11 U12 U13)
  ! A^T = (L21  1    ) * (    U22 U23)
  !       (L31 L32  1)   (        U33)
  L21 = A(1:iend,1,2) / A(1:iend,1,1)
  L31 = A(1:iend,1,3) / A(1:iend,1,1)
  U22 = A(1:iend,2,2) - L21*A(1:iend,2,1)
  U23 = A(1:iend,3,2) - L21*A(1:iend,3,1)
  L32 =(A(1:iend,2,3) - L31*A(1:iend,2,1)) / U22
  U33 = A(1:iend,3,3) - L31*A(1:iend,3,1) - L32*U23
  
  ! Solve X(1,:) = A^-T ( B(1) )
  !                     (  0   )
  !                     (  0   )
  ! Solve Ly = B(:,:,j) by forward substitution
  ! y1 = B(:,1)
  y2 = - L21*B(1:iend,1)
  y3 = - L31*B(1:iend,1) - L32*y2
  ! Solve UX(:,:,j) = y by back substitution
  X(1:iend,1,3) = y3 / U33
  X(1:iend,1,2) = (y2 - U23*X(1:iend,1,3)) / U22
  X(1:iend,1,1) = (B(1:iend,1) - A(1:iend,2,1)*X(1:iend,1,2) &
       &          - A(1:iend,3,1)*X(1:iend,1,3)) / A(1:iend,1,1)
  
  ! Solve X(2,:) = A^-T (  0   )
  !                     ( B(2) )
  !                     (  0   )
  ! Solve Ly = B(:,:,j) by forward substitution
  ! y1 = 0
  ! y2 = B(1:iend,2)
  y3 = - L32*B(1:iend,2)
  ! Solve UX(:,:,j) = y by back substitution
  X(1:iend,2,3) = y3 / U33
  X(1:iend,2,2) = (B(1:iend,2) - U23*X(1:iend,2,3)) / U22
  X(1:iend,2,1) = (-A(1:iend,2,1)*X(1:iend,2,2) &
       &           -A(1:iend,3,1)*X(1:iend,2,3)) / A(1:iend,1,1)
  
  ! Solve X(3,:) = A^-T (  0   )
  !                     (  0   )
  !                     ( B(3) )
  ! Solve Ly = B(:,:,j) by forward substitution
  ! y1 = 0
  ! y2 = 0
  ! y3 = B(1:iend,3)
  ! Solve UX(:,:,j) = y by back substitution
  X(1:iend,3,3) = B(1:iend,3) / U33
  X(1:iend,3,2) = -U23*X(1:iend,3,3) / U22
  X(1:iend,3,1) = (-A(1:iend,2,1)*X(1:iend,3,2) &
       &          - A(1:iend,3,1)*X(1:iend,3,3)) / A(1:iend,1,1)
  
end subroutine diag_mat_right_divide_3

!---------------------------------------------------------------------
! Matrix exponential of a 3x3 tridiagonal matrix, using Viete's method
! to solve the cubic equation to find the eigenvalues of the matrix,
! then using the diagonalization method with the eigenvalues and
! eigenvectors to perform the matrix exponentiation.
subroutine expm_tridiagonal(n, mat, ans)

  use parkind1, only : jprb

  integer,    intent(in)  :: n
  real(jprb), intent(in)  :: mat(n,3,3)
  real(jprb), intent(out) :: ans(n,3,3)

  real(jprb), parameter :: PI = acos(-1.0_jprb)
  real(jprb), parameter :: TWO_PI_OVER_THREE = 2.0_jprb * PI / 3.0_jprb

  ! Coefficients of a cubic equation x^3+bx^2+cx+d=0
  real(jprb) :: b, c, d
  ! Coefficients of reduced cubic t^3+px+q=0
  real(jprb) :: p, q
  ! Terms in Viete's cubic solution
  real(jprb) :: coeff1, coeff2, b_over_3

  ! Eigenvalues, then exp(eigenvalues)
  real(jprb) :: lambda(n,3)
  ! Eigenvectors
  real(jprb) :: eigenvec(n,3,3)

  ! Result of lambda right-divided by eigenvec
  real(jprb), dimension(n,3,3) :: lambda_rdivide_eigenvec

  integer :: j, jlambda, j1, j2

  ! First compute the three eigenvalues of the matrix by solving a
  ! cubic equation
  do j = 1,n
    ! Coefficients of a cubic equation x^3+bx^2+cx+d=0, where x is an
    ! eigenvalue
    b = -(mat(j,1,1) + mat(j,2,2) + mat(j,3,3))
    c = mat(j,1,1)*mat(j,3,3) + mat(j,1,1)*mat(j,2,2) + mat(j,2,2)*mat(j,3,3) &
         &  - mat(j,1,2)*mat(j,2,1) - mat(j,2,3)*mat(j,3,2)
    d = mat(j,1,1)*mat(j,2,3)*mat(j,3,2) + mat(j,1,2)*mat(j,2,1)*mat(j,3,3) &
         &  - mat(j,1,1)*mat(j,2,2)*mat(j,3,3)
    ! Coefficients of reduced cubic t^3+px+q=0
    b_over_3 = (1.0_jprb/3.0_jprb) * b
    p = c - b_over_3*b
    q = (2.0_jprb/27.0_jprb)*b*b*b - b_over_3*c + d
    ! Viete's solution for the real roots of a cubic equation
    coeff1 = 2.0_jprb*sqrt(-p*(1.0_jprb/3.0_jprb))
    coeff2 = (1.0_jprb/3.0_jprb) * acos(1.5_jprb*sqrt(-3.0_jprb/p)*q/p)
    lambda(j,1) = coeff1*cos(coeff2) - b_over_3
    lambda(j,2) = coeff1*cos(coeff2+TWO_PI_OVER_THREE) - b_over_3
    lambda(j,3) = coeff1*cos(coeff2+TWO_PI_OVER_THREE) - b_over_3
  end do

  do jlambda = 1,3
    eigenvec(:,1,jlambda) = mat(:,1,2)*(lambda(:,jlambda)-mat(:,3,3))
    eigenvec(:,2,jlambda) = (lambda(:,jlambda)-mat(:,1,1))*(lambda(:,jlambda)-mat(:,3,3))
    eigenvec(:,3,jlambda) = mat(:,3,2)*(lambda(:,jlambda)-mat(:,1,1))
  end do

  ! Diagonalization method to exponentiate a matrix

  lambda = exp(lambda)
  
  ! Compute lambda_rdivide_eigenvec = lambda * eigenvec^-1
  call diag_mat_right_divide_3(n,n,eigenvec,lambda,lambda_rdivide_eigenvec)

  ! Compute V * diag_rdivide_V
  do j1 = 1,3
    do j2 = 1,3
      ans(:,j2,j1) = eigenvec(:,j2,1)*lambda_rdivide_eigenvec(:,1,j1) &
           &       + eigenvec(:,j2,2)*lambda_rdivide_eigenvec(:,2,j1) &
           &       + eigenvec(:,j2,3)*lambda_rdivide_eigenvec(:,3,j1)
    end do
  end do

end subroutine expm_tridiagonal


#endif
