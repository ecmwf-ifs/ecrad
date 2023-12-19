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

  use parkind1, only : jprb

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

!---------------------------------------------------------------------
! Treat A as NREGION m-by-m square matrices (with the n dimension varying
! fastest) and b as n m-element vectors, and perform matrix-vector
! multiplications on all pairs
function mat_x_vec(n,A,b)

  use parkind1, only : jprb

  integer,    intent(in)                               :: n
  real(jprb), intent(in), dimension(n,NREGION,NREGION) :: A
  real(jprb), intent(in), dimension(n,NREGION)         :: b
  real(jprb),             dimension(n,NREGION)         :: mat_x_vec
  
  integer :: j1, j2

  ! Array-wise assignment
  mat_x_vec = 0.0_jprb

  do j1 = 1,NREGION
    do j2 = 1,NREGION
      mat_x_vec(:,j1) = mat_x_vec(:,j1) &
           &               + A(:,j1,j2)*b(:,j2)
    end do
  end do
    
end function mat_x_vec


#if NUM_REGIONS == 3

!---------------------------------------------------------------------
! Return X = B A^-1 = (A^-T B)^T optimized for 3x3 matrices, where B
! is a diagonal matrix, using LU factorization and substitution with
! no pivoting.
pure subroutine diag_mat_right_divide_3(n,A,B,X)

  use parkind1, only : jpim, jprb

  integer(jpim), intent(in)  :: n
  real(jprb),    intent(in)  :: A(n,3,3)
  real(jprb),    intent(in)  :: B(n,3)
  real(jprb),    intent(out) :: X(n,3,3)
  
  real(jprb), dimension(n) :: L21, L31, L32
  real(jprb), dimension(n) :: U22, U23, U33
  real(jprb), dimension(n) :: y2, y3

  !    associate (U11 => A(:,1,1), U12 => A(:,1,2), U13 => A(1,3))
  ! LU decomposition of the *transpose* of A:
  !       ( 1        )   (U11 U12 U13)
  ! A^T = (L21  1    ) * (    U22 U23)
  !       (L31 L32  1)   (        U33)
  L21 = A(1:n,1,2) / A(1:n,1,1)
  L31 = A(1:n,1,3) / A(1:n,1,1)
  U22 = A(1:n,2,2) - L21*A(1:n,2,1)
  U23 = A(1:n,3,2) - L21*A(1:n,3,1)
  L32 =(A(1:n,2,3) - L31*A(1:n,2,1)) / U22
  U33 = A(1:n,3,3) - L31*A(1:n,3,1) - L32*U23
  
  ! Solve X(1,:) = A^-T ( B(1) )
  !                     (  0   )
  !                     (  0   )
  ! Solve Ly = B(:,:,j) by forward substitution
  ! y1 = B(:,1)
  y2 = - L21*B(1:n,1)
  y3 = - L31*B(1:n,1) - L32*y2
  ! Solve UX(:,:,j) = y by back substitution
  X(1:n,1,3) = y3 / U33
  X(1:n,1,2) = (y2 - U23*X(1:n,1,3)) / U22
  X(1:n,1,1) = (B(1:n,1) - A(1:n,2,1)*X(1:n,1,2) &
       &      - A(1:n,3,1)*X(1:n,1,3)) / A(1:n,1,1)
  
  ! Solve X(2,:) = A^-T (  0   )
  !                     ( B(2) )
  !                     (  0   )
  ! Solve Ly = B(:,:,j) by forward substitution
  ! y1 = 0
  ! y2 = B(1:n,2)
  y3 = - L32*B(1:n,2)
  ! Solve UX(:,:,j) = y by back substitution
  X(1:n,2,3) = y3 / U33
  X(1:n,2,2) = (B(1:n,2) - U23*X(1:n,2,3)) / U22
  X(1:n,2,1) = (-A(1:n,2,1)*X(1:n,2,2) &
       &        -A(1:n,3,1)*X(1:n,2,3)) / A(1:n,1,1)
  
  ! Solve X(3,:) = A^-T (  0   )
  !                     (  0   )
  !                     ( B(3) )
  ! Solve Ly = B(:,:,j) by forward substitution
  ! y1 = 0
  ! y2 = 0
  ! y3 = B(1:n,3)
  ! Solve UX(:,:,j) = y by back substitution
  X(1:n,3,3) = B(1:n,3) / U33
  X(1:n,3,2) = -U23*X(1:n,3,3) / U22
  X(1:n,3,1) = (-A(1:n,2,1)*X(1:n,3,2) &
                -A(1:n,3,1)*X(1:n,3,3)) / A(1:n,1,1)
  
end subroutine diag_mat_right_divide_3

!---------------------------------------------------------------------
! Return X = B A^-1 = (A^-T B)^T optimized for 3x3 matrices, where B
! is a diagonal matrix, using LU factorization and substitution with
! no pivoting.
pure subroutine diag_mat_right_divide_mask_3(n,A,B,mask,X)

  use parkind1, only : jpim, jprb

  integer(jpim), intent(in)  :: n
  real(jprb),    intent(in)  :: A(n,3,3)
  real(jprb),    intent(in)  :: B(n,3)
  logical,       intent(in)  :: mask(n)
  real(jprb),    intent(out) :: X(n,3,3)
  
  real(jprb) :: L21, L31, L32
  real(jprb) :: U22, U23, U33
  real(jprb) :: y2, y3

  integer :: j

  do j = 1,n
    if (mask(j)) then
      !    associate (U11 => A(:,1,1), U12 => A(:,1,2), U13 => A(1,3))
      ! LU decomposition of the *transpose* of A:
      !       ( 1        )   (U11 U12 U13)
      ! A^T = (L21  1    ) * (    U22 U23)
      !       (L31 L32  1)   (        U33)
      L21 = A(j,1,2) / A(j,1,1)
      L31 = A(j,1,3) / A(j,1,1)
      U22 = A(j,2,2) - L21*A(j,2,1)
      U23 = A(j,3,2) - L21*A(j,3,1)
      L32 =(A(j,2,3) - L31*A(j,2,1)) / U22
      U33 = A(j,3,3) - L31*A(j,3,1) - L32*U23
      
      ! Solve X(1,:) = A^-T ( B(1) )
      !                     (  0   )
      !                     (  0   )
      ! Solve Ly = B(:,:,j) by forward substitution
      ! y1 = B(:,1)
      y2 = - L21*B(j,1)
      y3 = - L31*B(j,1) - L32*y2
      ! Solve UX(:,:,j) = y by back substitution
      X(j,1,3) = y3 / U33
      X(j,1,2) = (y2 - U23*X(j,1,3)) / U22
      X(j,1,1) = (B(j,1) - A(j,2,1)*X(j,1,2) &
           &      - A(j,3,1)*X(j,1,3)) / A(j,1,1)
      
      ! Solve X(2,:) = A^-T (  0   )
      !                     ( B(2) )
      !                     (  0   )
      ! Solve Ly = B(:,:,j) by forward substitution
      ! y1 = 0
      ! y2 = B(j,2)
      y3 = - L32*B(j,2)
      ! Solve UX(:,:,j) = y by back substitution
      X(j,2,3) = y3 / U33
      X(j,2,2) = (B(j,2) - U23*X(j,2,3)) / U22
      X(j,2,1) = (-A(j,2,1)*X(j,2,2) &
           &        -A(j,3,1)*X(j,2,3)) / A(j,1,1)
      
      ! Solve X(3,:) = A^-T (  0   )
      !                     (  0   )
      !                     ( B(3) )
      ! Solve Ly = B(:,:,j) by forward substitution
      ! y1 = 0
      ! y2 = 0
      ! y3 = B(j,3)
      ! Solve UX(:,:,j) = y by back substitution
      X(j,3,3) = B(j,3) / U33
      X(j,3,2) = -U23*X(j,3,3) / U22
      X(j,3,1) = (-A(j,2,1)*X(j,3,2) &
           -A(j,3,1)*X(j,3,3)) / A(j,1,1)
    end if
  end do

end subroutine diag_mat_right_divide_mask_3

!---------------------------------------------------------------------
! Invert a tridiagonal 3x3 matrix using Thomas's algorithm but
! starting at the end to avoid dividing by the mat(:,1,1) element
! which is most likely to be zero
subroutine inv_tridiagonal(n, mat, ans)

  use parkind1, only : jprb

  integer,    intent(in)  :: n
  real(jprb), intent(in)  :: mat(n,3,3)
  real(jprb), intent(out) :: ans(n,3,3)

  real(jprb) :: g1(n), g2(n), r2(n)

  ! First factorize "mat" into a matrix of the form
  ! (1      )
  ! (g2 1   )
  ! (   g1 1)
  g1 = mat(:,3,2) / mat(:,3,3)
  g2 = mat(:,2,1) / (mat(:,2,2)-mat(:,2,3)*g1)

  ! Back-substitution (but actually working forward through the matrix
  ! since we did the first step in the backwards direction) with the
  ! identity matrix on the right-hand side.

  ! First column of answer, for RHS=[1 0 0]^T:
  ans(:,1,1) = 1.0_jprb / (mat(:,1,1) - mat(:,1,2)*g2)
  ans(:,2,1) = -g2 * ans(:,1,1)
  ans(:,3,1) = -g1 * ans(:,2,1)

  ! Second column of answer, for RHS=[0 1 0]^T, in which case after
  ! factorization it becomes [r3 r2 0]^T:
  r2 = 1.0_jprb / (mat(:,2,2)-mat(:,2,3)*g1)
  ans(:,1,2) = -mat(:,1,2)*r2 * ans(:,1,1) ! = r3
  ans(:,2,2) = r2 - g2 * ans(:,1,2)
  ans(:,3,2) = -g1 * ans(:,2,2)

  ! Third column of answer, for RHS=[0 0 1]^T, in which case after
  ! factorization it becomes [r3 r2 r1]^T:
  ! r1 = 1.0/mat(:,3,3)
  r2 = r2 * (-mat(:,2,3)/mat(:,3,3))
  ans(:,1,3) = -mat(:,1,2)*r2 * ans(:,1,1) ! = r3
  ans(:,2,3) = r2 - g2 * ans(:,1,3)
  ans(:,3,3) = 1.0_jprb/mat(:,3,3) - g1 * ans(:,2,3)

end subroutine inv_tridiagonal


!---------------------------------------------------------------------
! Matrix exponential of a 3x3 tridiagonal matrix, using Viete's method
! to solve the cubic equation to find the eigenvalues of the matrix,
! then using the diagonalization method with the eigenvalues and
! eigenvectors to perform the matrix exponentiation. This routine
! makes several assumptions about the properties of the matrix that
! are valid for the way it is used in TCRAD but may not be for other
! applications.
subroutine expm_tridiagonal(n, mat, ans)

  use parkind1, only : jprb

  integer,    intent(in)  :: n
  real(jprb), intent(in)  :: mat(n,3,3)
  real(jprb), intent(out) :: ans(n,3,3)

  real(jprb), parameter :: PI = acos(-1.0_jprb)
  real(jprb), parameter :: TWO_PI_OVER_THREE = 2.0_jprb * PI / 3.0_jprb

  real(jprb), parameter :: MIN_DIAGONAL = -100.0_jprb
  real(jprb), parameter :: MIN_RATIO    = 1.001_jprb

  ! Coefficients of a cubic equation x^3+bx^2+cx+d=0
  real(jprb) :: b, c, d
  ! Coefficients of reduced cubic t^3+px+q=0
  real(jprb) :: p, q
  ! Terms in Viete's cubic solution
  real(jprb) :: coeff1, coeff2, b_over_3, acos_arg

  ! Eigenvalues, then exp(eigenvalues)
  real(jprb) :: lambda(n,3)
  ! Eigenvectors
  real(jprb) :: eigenvec(n,3,3)

  ! Result of lambda right-divided by eigenvec
  real(jprb), dimension(n,3,3) :: lambda_rdivide_eigenvec

  ! If the diagonals are too large or too similar then the eigenvalue
  ! method may not work, in which case we ignore off-diagonals and
  ! just take the exponent of the diagonals
  logical :: eigen_mask(n)

  integer :: j, jlambda, j1, j2

  eigen_mask = .false.

  ! First compute the three eigenvalues of the matrix by solving a
  ! cubic equation
  do j = 1,n
    ! Default values
    lambda(j,1) = mat(j,1,1)
    lambda(j,2) = mat(j,2,2)
    lambda(j,3) = mat(j,3,3)

    if (mat(j,1,1) >= MIN_DIAGONAL & !.and. mat(j,2,2) <= MIN_RATIO*mat(j,1,1) &
         & .and. mat(j,1,2) > 1.0e-8_jprb .and. mat(j,3,2) > 1.0e-8_jprb) then
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
      acos_arg = 1.5_jprb*sqrt(-3.0_jprb/p)*q/p
      if (acos_arg >= -1.0_jprb .and. acos_arg <= 1.0_jprb) then
        coeff2 = (1.0_jprb/3.0_jprb) * acos(acos_arg)
        lambda(j,1) = coeff1*cos(coeff2) - b_over_3
        lambda(j,2) = coeff1*cos(coeff2+TWO_PI_OVER_THREE) - b_over_3
        lambda(j,3) = coeff1*cos(coeff2+2.0_jprb*TWO_PI_OVER_THREE) - b_over_3
        eigen_mask(j) = .true.
      end if
    end if
  end do

  do j = 1,n
    if (eigen_mask(j)) then
      do jlambda = 1,3
        eigenvec(j,1,jlambda) = mat(j,1,2)*(lambda(j,jlambda)-mat(j,3,3))
        eigenvec(j,2,jlambda) = (lambda(j,jlambda)-mat(j,1,1))*(lambda(j,jlambda)-mat(j,3,3))
        eigenvec(j,3,jlambda) = mat(j,3,2)*(lambda(j,jlambda)-mat(j,1,1))
      end do
    end if
  end do

  ! Diagonalization method to exponentiate a matrix

  lambda = exp(lambda)
  
  ! Compute lambda_rdivide_eigenvec = lambda * eigenvec^-1
  call diag_mat_right_divide_mask_3(n,eigenvec,lambda,eigen_mask,lambda_rdivide_eigenvec)

  ! Compute V * diag_rdivide_V
  do j = 1,n
    if (eigen_mask(j)) then
      do j1 = 1,3
        do j2 = 1,3
          ans(j,j2,j1) = eigenvec(j,j2,1)*lambda_rdivide_eigenvec(j,1,j1) &
               &       + eigenvec(j,j2,2)*lambda_rdivide_eigenvec(j,2,j1) &
               &       + eigenvec(j,j2,3)*lambda_rdivide_eigenvec(j,3,j1)
        end do
      end do
    end if
  end do

  if (any(.not. eigen_mask)) then
    do j = 1,n
      if (.not. eigen_mask(j)) then
        ans(j,:,:) = 0.0_jprb
        ans(j,1,1) = lambda(j,1);
        ans(j,2,2) = lambda(j,2);
        ans(j,3,3) = lambda(j,3);
      end if
    end do
  end if

end subroutine expm_tridiagonal


#else


!---------------------------------------------------------------------
! Invert a tridiagonal 2x2 matrix
subroutine inv_tridiagonal(n, mat, ans)

  use parkind1, only : jprb

  integer,    intent(in)  :: n
  real(jprb), intent(in)  :: mat(n,2,2)
  real(jprb), intent(out) :: ans(n,2,2)

  real(jprb) :: inv_det(n)

  inv_det = 1.0_jprb / (mat(:,1,1)*mat(:,2,2) - mat(:,1,2)*mat(:,2,1))
  ans(:,1,1) =  inv_det * mat(:,2,2)
  ans(:,2,1) = -inv_det * mat(:,2,1)
  ans(:,1,2) = -inv_det * mat(:,1,2)
  ans(:,2,2) =  inv_det * mat(:,1,1)

end subroutine inv_tridiagonal

!---------------------------------------------------------------------
! Matrix exponential of a 2x2 matrix
subroutine expm_tridiagonal(n, mat, ans)

  use parkind1, only : jprb

  integer,    intent(in)  :: n
  real(jprb), intent(in)  :: mat(n,2,2)
  real(jprb), intent(out) :: ans(n,2,2)

  real(jprb) :: delta(n), factor(n), sinh_delta(n), cosh_delta(n)

  delta = 0.5_jprb * sqrt((mat(:,1,1)-mat(:,2,2))**2 + 4.0_jprb*mat(:,2,1)*mat(:,1,2))
  factor = exp(0.5_jprb*(mat(:,1,1)+mat(:,2,2))) / delta
  sinh_delta = sinh(delta)
  cosh_delta = cosh(delta)

  ans(:,1,1) = 0.5_jprb * factor * (2.0_jprb*delta * cosh_delta + (mat(:,1,1)-mat(:,2,2))*sinh_delta)
  ans(:,2,1) = mat(:,2,1) * factor * sinh_delta
  ans(:,1,2) = mat(:,1,2) * factor * sinh_delta
  ans(:,2,2) = 0.5_jprb * factor * (2.0_jprb*delta * cosh_delta + (mat(:,2,2)-mat(:,1,1))*sinh_delta)

end subroutine expm_tridiagonal


#endif
