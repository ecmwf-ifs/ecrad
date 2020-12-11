! radiation_matrix.F90 - SPARTACUS matrix operations
!
! (C) Copyright 2014- ECMWF.
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
! Modifications
!   2018-10-15  R. Hogan    Added fast_expm_exchange_[23]
!   2020-12-xx  P. Ukkonen  Added an optimized expm routine for shortwave when nreg=3,
!                           and related kernels
!
! This module provides the neccessary mathematical functions for the
! SPARTACUS radiation scheme: matrix multiplication, matrix solvers
! and matrix exponentiation, but (a) multiple matrices are operated on
! at once with array access indended to facilitate vectorization, and
! (b) optimization for 2x2 and 3x3 matrices.  There is probably
! considerable scope for further optimization. Note that this module
! is not used by the McICA solver.

module radiation_matrix

  use parkind1, only : jprb

  implicit none
  public

  ! Codes to describe sparseness pattern, where the SHORTWAVE
  ! pattern is of the form:
  ! (x x x)
  ! (x x x)
  ! (0 0 x)
  ! where each element may itself be a square matrix.  
  integer, parameter :: IMatrixPatternDense     = 0
  integer, parameter :: IMatrixPatternShortwave = 1

  public  :: mat_x_vec, singlemat_x_vec, mat_x_mat, &
       &     singlemat_x_mat, mat_x_singlemat, &
       &     identity_minus_mat_x_mat, solve_vec, solve_mat, expm, &
       &     fast_expm_exchange_2, fast_expm_exchange_3

  private :: solve_vec_2, solve_vec_3, solve_mat_2, &
       &     solve_mat_3, lu_factorization, lu_substitution, solve_mat_n, &
       &     diag_mat_right_divide_3

  interface fast_expm_exchange
    module procedure fast_expm_exchange_2, fast_expm_exchange_3
  end interface fast_expm_exchange

contains

  ! --- MATRIX-VECTOR MULTIPLICATION ---

  !---------------------------------------------------------------------
  ! Treat A as n m-by-m square matrices (with the n dimension varying
  ! fastest) and b as n m-element vectors, and perform matrix-vector
  ! multiplications on first iend pairs
  function mat_x_vec(n,iend,m,A,b,do_top_left_only_in)

    use yomhook, only : lhook, dr_hook

    integer,    intent(in)                   :: n, m, iend
    real(jprb), intent(in), dimension(:,:,:) :: A
    real(jprb), intent(in), dimension(:,:)   :: b
    logical,    intent(in), optional         :: do_top_left_only_in
    real(jprb),             dimension(iend,m):: mat_x_vec

    integer :: j1, j2
    logical :: do_top_left_only

    real(jprb) :: hook_handle

    if (lhook) call dr_hook('radiation_matrix:mat_x_vec',0,hook_handle)


    if (present(do_top_left_only_in)) then
      do_top_left_only = do_top_left_only_in
    else
      do_top_left_only = .false.
    end if

    ! Array-wise assignment
    mat_x_vec = 0.0_jprb

    if (do_top_left_only) then
      mat_x_vec(1:iend,1) = A(1:iend,1,1)*b(1:iend,1)
    else
      do j1 = 1,m
        do j2 = 1,m
          mat_x_vec(1:iend,j1) = mat_x_vec(1:iend,j1) &
               &               + A(1:iend,j1,j2)*b(1:iend,j2)
        end do
      end do
    end if

    if (lhook) call dr_hook('radiation_matrix:mat_x_vec',1,hook_handle)

  end function mat_x_vec

  pure function mat_x_vec_3(n,A,b)

    integer,    intent(in)                   :: n
    real(jprb), intent(in), dimension(n,3,3) :: A
    real(jprb), intent(in), dimension(n,3)   :: b
    real(jprb), dimension(n,3):: mat_x_vec_3
    integer :: j1

    ! Array-wise assignment
    mat_x_vec_3 = 0.0_jprb

    do j1 = 1,3
        mat_x_vec_3(:,j1) = mat_x_vec_3(:,j1) + A(:,j1,1)*b(:,1)
        mat_x_vec_3(:,j1) = mat_x_vec_3(:,j1) + A(:,j1,2)*b(:,2)
        mat_x_vec_3(:,j1) = mat_x_vec_3(:,j1) + A(:,j1,3)*b(:,3)
    end do

  end function mat_x_vec_3


  !---------------------------------------------------------------------
  ! Treat A as an m-by-m square matrix and b as n m-element vectors
  ! (with the n dimension varying fastest), and perform matrix-vector
  ! multiplications on first iend pairs
  function singlemat_x_vec(n,iend,m,A,b)

    use yomhook, only : lhook, dr_hook

    integer,    intent(in)                    :: n, m, iend
    real(jprb), intent(in), dimension(m,m)    :: A
    real(jprb), intent(in), dimension(:,:)    :: b
    real(jprb),             dimension(iend,m) :: singlemat_x_vec

    integer    :: j1, j2
    real(jprb) :: hook_handle
    
    if (lhook) call dr_hook('radiation_matrix:single_mat_x_vec',0,hook_handle)

    ! Array-wise assignment
    singlemat_x_vec = 0.0_jprb

    do j1 = 1,m
      do j2 = 1,m
        singlemat_x_vec(1:iend,j1) = singlemat_x_vec(1:iend,j1) &
             &                    + A(j1,j2)*b(1:iend,j2)
      end do
    end do

    if (lhook) call dr_hook('radiation_matrix:single_mat_x_vec',1,hook_handle)

  end function singlemat_x_vec


  ! --- SQUARE MATRIX-MATRIX MULTIPLICATION ---

  !---------------------------------------------------------------------
  ! Treat A and B each as n m-by-m square matrices (with the n
  ! dimension varying fastest) and perform matrix multiplications on
  ! all n matrix pairs
  function mat_x_mat(n,iend,m,A,B,i_matrix_pattern)

    use yomhook, only : lhook, dr_hook

    integer,    intent(in)                      :: n, m, iend
    integer,    intent(in), optional            :: i_matrix_pattern
    real(jprb), intent(in), dimension(:,:,:)    :: A, B

    real(jprb),             dimension(iend,m,m) :: mat_x_mat
    integer    :: j1, j2, j3
    integer    :: mblock, m2block
    integer    :: i_actual_matrix_pattern
    real(jprb) :: hook_handle

    if (lhook) call dr_hook('radiation_matrix:mat_x_mat',0,hook_handle)

    if (present(i_matrix_pattern)) then
      i_actual_matrix_pattern = i_matrix_pattern
    else
      i_actual_matrix_pattern = IMatrixPatternDense
    end if

    ! Array-wise assignment
    mat_x_mat = 0.0_jprb

    if (i_actual_matrix_pattern == IMatrixPatternShortwave) then
      ! Matrix has a sparsity pattern
      !     (C D E)
      ! A = (F G H)
      !     (0 0 I)
      mblock = m/3
      m2block = 2*mblock 
      ! Do the top-left (C, D, F, G)
      do j2 = 1,m2block
        do j1 = 1,m2block
          do j3 = 1,m2block
            mat_x_mat(1:iend,j1,j2) = mat_x_mat(1:iend,j1,j2) &
                 &                  + A(1:iend,j1,j3)*B(1:iend,j3,j2)
          end do
        end do
      end do
      do j2 = m2block+1,m
        ! Do the top-right (E & H)
        do j1 = 1,m2block
          do j3 = 1,m
            mat_x_mat(1:iend,j1,j2) = mat_x_mat(1:iend,j1,j2) &
                 &                  + A(1:iend,j1,j3)*B(1:iend,j3,j2)
          end do
        end do
        ! Do the bottom-right (I)
        do j1 = m2block+1,m
          do j3 = m2block+1,m
            mat_x_mat(1:iend,j1,j2) = mat_x_mat(1:iend,j1,j2) &
                 &                  + A(1:iend,j1,j3)*B(1:iend,j3,j2)
          end do
        end do
      end do
    else
      ! Ordinary dense matrix
      do j2 = 1,m
        do j1 = 1,m
          do j3 = 1,m
            mat_x_mat(1:iend,j1,j2) = mat_x_mat(1:iend,j1,j2) &
                 &                  + A(1:iend,j1,j3)*B(1:iend,j3,j2)
          end do
        end do
      end do
    end if

    if (lhook) call dr_hook('radiation_matrix:mat_x_mat',1,hook_handle)

  end function mat_x_mat

  pure function mat_x_mat_dense(n,m,A,B)

    integer,    intent(in)                      :: n, m
    real(jprb), intent(in), dimension(n,m,m)    :: A, B

    real(jprb),             dimension(n,m,m) :: mat_x_mat_dense
    integer    :: j1, j2, j3

    if (m==3) then
      do j2 = 1,m
        do j1 = 1,m
          mat_x_mat_dense(:,j1,j2) = A(:,j1,1)*B(:,1,j2) &
          &  + A(:,j1,2)*B(:,2,j2) + A(:,j1,3)*B(:,3,j2) 
        end do
      end do

    else 
  
      ! Array-wise assignment
      mat_x_mat_dense = 0.0_jprb

      do j2 = 1,m
        do j1 = 1,m
          do j3 = 1,m
            mat_x_mat_dense(:,j1,j2) = mat_x_mat_dense(:,j1,j2) &
                  &                  + A(:,j1,j3)*B(:,j3,j2)
          end do
        end do
      end do

    end if

  end function mat_x_mat_dense
 
  !---------------------------------------------------------------------
  ! Treat A as an m-by-m matrix and B as n m-by-m square matrices
  ! (with the n dimension varying fastest) and perform matrix
  ! multiplications on the first iend matrix pairs
  function singlemat_x_mat(n,iend,m,A,B)

    use yomhook, only : lhook, dr_hook

    integer,    intent(in)                      :: n, m, iend
    real(jprb), intent(in), dimension(m,m)      :: A
    real(jprb), intent(in), dimension(:,:,:)    :: B
    real(jprb),             dimension(iend,m,m) :: singlemat_x_mat

    integer    :: j1, j2, j3
    real(jprb) :: hook_handle

    if (lhook) call dr_hook('radiation_matrix:singlemat_x_mat',0,hook_handle)

    ! Array-wise assignment
    singlemat_x_mat = 0.0_jprb

    do j2 = 1,m
      do j1 = 1,m
        do j3 = 1,m
          singlemat_x_mat(1:iend,j1,j2) = singlemat_x_mat(1:iend,j1,j2) &
               &                        + A(j1,j3)*B(1:iend,j3,j2)
        end do
      end do
    end do

    if (lhook) call dr_hook('radiation_matrix:singlemat_x_mat',1,hook_handle)

  end function singlemat_x_mat


  !---------------------------------------------------------------------
  ! Treat B as an m-by-m matrix and A as n m-by-m square matrices
  ! (with the n dimension varying fastest) and perform matrix
  ! multiplications on the first iend matrix pairs
  function mat_x_singlemat(n,iend,m,A,B)

    use yomhook, only : lhook, dr_hook

    integer,    intent(in)                      :: n, m, iend
    real(jprb), intent(in), dimension(:,:,:)    :: A
    real(jprb), intent(in), dimension(m,m)      :: B

    real(jprb),             dimension(iend,m,m) :: mat_x_singlemat
    integer    :: j1, j2, j3
    real(jprb) :: hook_handle

    if (lhook) call dr_hook('radiation_matrix:mat_x_singlemat',0,hook_handle)

    ! Array-wise assignment
    mat_x_singlemat = 0.0_jprb

    do j2 = 1,m
      do j1 = 1,m
        do j3 = 1,m
          mat_x_singlemat(1:iend,j1,j2) = mat_x_singlemat(1:iend,j1,j2) &
               &                        + A(1:iend,j1,j3)*B(j3,j2)
        end do
      end do
    end do

    if (lhook) call dr_hook('radiation_matrix:mat_x_singlemat',1,hook_handle)

  end function mat_x_singlemat


  !---------------------------------------------------------------------
  ! Compute I-A*B where I is the identity matrix and A & B are n
  ! m-by-m square matrices
  function identity_minus_mat_x_mat(n,iend,m,A,B,i_matrix_pattern)

    use yomhook, only : lhook, dr_hook

    integer,    intent(in)                   :: n, m, iend
    integer,    intent(in), optional         :: i_matrix_pattern
    real(jprb), intent(in), dimension(:,:,:) :: A, B
    real(jprb),             dimension(iend,m,m) :: identity_minus_mat_x_mat

    integer    :: j
    real(jprb) :: hook_handle

    if (lhook) call dr_hook('radiation_matrix:identity_mat_x_mat',0,hook_handle)

    if (present(i_matrix_pattern)) then
      identity_minus_mat_x_mat = mat_x_mat(n,iend,m,A,B,i_matrix_pattern)
    else
      identity_minus_mat_x_mat = mat_x_mat(n,iend,m,A,B)
    end if

    identity_minus_mat_x_mat = - identity_minus_mat_x_mat
    do j = 1,m
      identity_minus_mat_x_mat(1:iend,j,j) &
           &     = 1.0_jprb + identity_minus_mat_x_mat(1:iend,j,j)
    end do

    if (lhook) call dr_hook('radiation_matrix:identity_mat_x_mat',1,hook_handle)

  end function identity_minus_mat_x_mat


  ! --- REPEATEDLY SQUARE A MATRIX ---

  !---------------------------------------------------------------------
  ! Square m-by-m matrix "A" nrepeat times. A will be corrupted by
  ! this function.
  function repeated_square(m,A,nrepeat,i_matrix_pattern)
    integer,    intent(in)           :: m, nrepeat
    real(jprb), intent(inout)        :: A(m,m)
    integer,    intent(in), optional :: i_matrix_pattern
    real(jprb)                       :: repeated_square(m,m)

    integer :: j1, j2, j3, j4
    integer :: mblock, m2block
    integer :: i_actual_matrix_pattern

    if (present(i_matrix_pattern)) then
      i_actual_matrix_pattern = i_matrix_pattern
    else
      i_actual_matrix_pattern = IMatrixPatternDense
    end if

    if (i_actual_matrix_pattern == IMatrixPatternShortwave) then
      ! Matrix has a sparsity pattern
      !     (C D E)
      ! A = (F G H)
      !     (0 0 I)
      mblock = m/3
      m2block = 2*mblock
      do j4 = 1,nrepeat
        repeated_square = 0.0_jprb
        ! Do the top-left (C, D, F & G)
        do j2 = 1,m2block
          do j1 = 1,m2block
            do j3 = 1,m2block
              repeated_square(j1,j2) = repeated_square(j1,j2) &
                   &                 + A(j1,j3)*A(j3,j2)
            end do
          end do
        end do
        do j2 = m2block+1, m
          ! Do the top-right (E & H)
          do j1 = 1,m2block
            do j3 = 1,m
              repeated_square(j1,j2) = repeated_square(j1,j2) &
                   &                 + A(j1,j3)*A(j3,j2)
            end do
          end do
          ! Do the bottom-right (I)
          do j1 = m2block+1, m
            do j3 = m2block+1,m
              repeated_square(j1,j2) = repeated_square(j1,j2) &
                   &                 + A(j1,j3)*A(j3,j2)
            end do
          end do
        end do
        if (j4 < nrepeat) then
          A = repeated_square
        end if
      end do
    else
      ! Ordinary dense matrix
      do j4 = 1,nrepeat
        repeated_square = 0.0_jprb
        do j2 = 1,m
          do j1 = 1,m
            do j3 = 1,m
              repeated_square(j1,j2) = repeated_square(j1,j2) &
                   &                 + A(j1,j3)*A(j3,j2)
            end do
          end do
        end do
        if (j4 < nrepeat) then
          A = repeated_square
        end if
      end do
    end if

  end function repeated_square

  pure subroutine repeated_square_sw_9(nrepeat,A,B)
    integer,    intent(in)            :: nrepeat
    real(jprb), intent(inout)         :: A(9,9)
    real(jprb), intent(out)           :: B(9,9)

    integer :: j1, j2, j3, j4

    do j4 = 1,nrepeat
      B = 0.0_jprb
      ! Do the top-left (C, D, F & G)
      do j2 = 1,6
          do j3 = 1,6
            B(1:6,j2) = B(1:6,j2) + A(1:6,j3)*A(j3,j2)
          end do
      end do
      do j2 = 7, 9
          do j3 = 1,9 ! Do the top-right (E & H)
            B(1:6,j2) = B(1:6,j2) + A(1:6,j3)*A(j3,j2)
          end do
          do j3 = 7,9  ! Do the bottom-right (I)
            B(7:9,j2) = B(7:9,j2) + A(7:9,j3)*A(j3,j2)
          end do
      end do
      if (j4 < nrepeat) then
          A = B
      end if
    end do
  end subroutine repeated_square_sw_9


  ! --- SOLVE LINEAR EQUATIONS ---

  !---------------------------------------------------------------------
  ! Solve Ax=b to obtain x.  Version optimized for 2x2 matrices using
  ! Cramer's method: "A" contains n 2x2 matrices and "b" contains n
  ! 2-element vectors; returns A^-1 b.
  pure subroutine solve_vec_2(n,iend,A,b,x)

    integer,    intent(in)  :: n, iend
    real(jprb), intent(in)  :: A(:,:,:)
    real(jprb), intent(in)  :: b(:,:)
    real(jprb), intent(out) :: x(:,:)

    real(jprb) :: inv_det(iend)

    inv_det = 1.0_jprb / (  A(1:iend,1,1)*A(1:iend,2,2) &
         &                - A(1:iend,1,2)*A(1:iend,2,1))

    x(1:iend,1) = inv_det*(A(1:iend,2,2)*b(1:iend,1)-A(1:iend,1,2)*b(1:iend,2))
    x(1:iend,2) = inv_det*(A(1:iend,1,1)*b(1:iend,2)-A(1:iend,2,1)*b(1:iend,1))

  end subroutine solve_vec_2


  !---------------------------------------------------------------------
  ! Solve AX=B to obtain X, i.e. the matrix right-hand-side version of
  ! solve_vec_2, with A, X and B all containing n 2x2 matrices;
  ! returns A^-1 B using Cramer's method.
  pure subroutine solve_mat_2(n,iend,A,B,X)
    integer,    intent(in)  :: n, iend
    real(jprb), intent(in)  :: A(:,:,:)
    real(jprb), intent(in)  :: B(:,:,:)
    real(jprb), intent(out) :: X(:,:,:)

    real(jprb) :: inv_det(iend)

    inv_det = 1.0_jprb / (  A(1:iend,1,1)*A(1:iend,2,2) &
         &                - A(1:iend,1,2)*A(1:iend,2,1))

    X(1:iend,1,1) = inv_det*( A(1:iend,2,2)*B(1:iend,1,1) &
         &                   -A(1:iend,1,2)*B(1:iend,2,1))
    X(1:iend,2,1) = inv_det*( A(1:iend,1,1)*B(1:iend,2,1) &
         &                   -A(1:iend,2,1)*B(1:iend,1,1))
    X(1:iend,1,2) = inv_det*( A(1:iend,2,2)*B(1:iend,1,2) &
         &                   -A(1:iend,1,2)*B(1:iend,2,2))
    X(1:iend,2,2) = inv_det*( A(1:iend,1,1)*B(1:iend,2,2) &
         &                   -A(1:iend,2,1)*B(1:iend,1,2))

  end subroutine solve_mat_2


  !---------------------------------------------------------------------
  ! Solve Ax=b optimized for 3x3 matrices, using LU
  ! factorization and substitution without pivoting.
  pure subroutine solve_vec_3(n,iend,A,b,x)
    integer,    intent(in)  :: n, iend
    real(jprb), intent(in)  :: A(:,:,:)
    real(jprb), intent(in)  :: b(:,:)
    real(jprb), intent(out) :: x(:,:)

    real(jprb), dimension(iend) :: L21, L31, L32
    real(jprb), dimension(iend) :: U22, U23, U33
    real(jprb), dimension(iend) :: y2, y3

    ! Some compilers unfortunately don't support assocate
    !    associate (U11 => A(:,1,1), U12 => A(:,1,2), U13 => A(1,3), &
    !         y1 => b(:,1), x1 => solve_vec3(:,1), &
    !         x2 => solve_vec3(:,2), x3 => solve_vec3(:,3))

    ! LU decomposition:
    !     ( 1        )   (U11 U12 U13)
    ! A = (L21  1    ) * (    U22 U23)
    !     (L31 L32  1)   (        U33)
    L21 = A(1:iend,2,1) / A(1:iend,1,1)
    L31 = A(1:iend,3,1) / A(1:iend,1,1)
    U22 = A(1:iend,2,2) - L21*A(1:iend,1,2)
    U23 = A(1:iend,2,3) - L21*A(1:iend,1,3)
    L32 =(A(1:iend,3,2) - L31*A(1:iend,1,2)) / U22
    U33 = A(1:iend,3,3) - L31*A(1:iend,1,3) - L32*U23

    ! Solve Ly = b by forward substitution
    y2 = b(1:iend,2) - L21*b(1:iend,1)
    y3 = b(1:iend,3) - L31*b(1:iend,1) - L32*y2

    ! Solve Ux = y by back substitution
    x(1:iend,3) = y3/U33
    x(1:iend,2) = (y2 - U23*x(1:iend,3)) / U22
    x(1:iend,1) = (b(1:iend,1) - A(1:iend,1,2)*x(1:iend,2) &
         &         - A(1:iend,1,3)*x(1:iend,3)) / A(1:iend,1,1)
    !    end associate

  end subroutine solve_vec_3


  !---------------------------------------------------------------------
  ! Solve AX=B optimized for 3x3 matrices, using LU factorization and
  ! substitution with no pivoting.
  pure subroutine solve_mat_3(n,iend,A,B,X)
    integer,    intent(in)  :: n, iend
    real(jprb), intent(in)  :: A(:,:,:)
    real(jprb), intent(in)  :: B(:,:,:)
    real(jprb), intent(out) :: X(:,:,:)

    real(jprb), dimension(iend) :: L21, L31, L32
    real(jprb), dimension(iend) :: U22, U23, U33
    real(jprb), dimension(iend) :: y2, y3

    integer :: j

    !    associate (U11 => A(:,1,1), U12 => A(:,1,2), U13 => A(1,3))
    ! LU decomposition:
    !     ( 1        )   (U11 U12 U13)
    ! A = (L21  1    ) * (    U22 U23)
    !     (L31 L32  1)   (        U33)
    L21 = A(1:iend,2,1) / A(1:iend,1,1)
    L31 = A(1:iend,3,1) / A(1:iend,1,1)
    U22 = A(1:iend,2,2) - L21*A(1:iend,1,2)
    U23 = A(1:iend,2,3) - L21*A(1:iend,1,3)
    L32 =(A(1:iend,3,2) - L31*A(1:iend,1,2)) / U22
    U33 = A(1:iend,3,3) - L31*A(1:iend,1,3) - L32*U23

    do j = 1,3
      ! Solve Ly = B(:,:,j) by forward substitution
      ! y1 = B(:,1,j)
      y2 = B(1:iend,2,j) - L21*B(1:iend,1,j)
      y3 = B(1:iend,3,j) - L31*B(1:iend,1,j) - L32*y2
      ! Solve UX(:,:,j) = y by back substitution
      X(1:iend,3,j) = y3 / U33
      X(1:iend,2,j) = (y2 - U23*X(1:iend,3,j)) / U22
      X(1:iend,1,j) = (B(1:iend,1,j) - A(1:iend,1,2)*X(1:iend,2,j) &
           &          - A(1:iend,1,3)*X(1:iend,3,j)) / A(1:iend,1,1)
    end do

  end subroutine solve_mat_3


  !---------------------------------------------------------------------
  ! Return X = B A^-1 = (A^-T B)^T optimized for 3x3 matrices, where B
  ! is a diagonal matrix, using LU factorization and substitution with
  ! no pivoting.
  pure subroutine diag_mat_right_divide_3(n,A,B,X)
    integer,    intent(in)  :: n
    real(jprb), intent(in)  :: A(n,3,3)
    real(jprb), intent(in)  :: B(n,3)
    real(jprb), intent(out) :: X(n,3,3)

    real(jprb), dimension(n) :: L21, L31, L32
    real(jprb), dimension(n) :: U22, U23, U33
    real(jprb), dimension(n) :: y2, y3

    integer :: j

    !    associate (U11 => A(:,1,1), U12 => A(:,1,2), U13 => A(1,3))
    ! LU decomposition of the *transpose* of A:
    !       ( 1        )   (U11 U12 U13)
    ! A^T = (L21  1    ) * (    U22 U23)
    !       (L31 L32  1)   (        U33)
    L21 = A(:,1,2) / A(:,1,1)
    L31 = A(:,1,3) / A(:,1,1)
    U22 = A(:,2,2) - L21*A(:,2,1)
    U23 = A(:,3,2) - L21*A(:,3,1)
    L32 =(A(:,2,3) - L31*A(:,2,1)) / U22
    U33 = A(:,3,3) - L31*A(:,3,1) - L32*U23

    ! Solve X(1,:) = A^-T ( B(1) )
    !                     (  0   )
    !                     (  0   )
    ! Solve Ly = B(:,:,j) by forward substitution
    ! y1 = B(:,1)
    y2 = - L21*B(:,1)
    y3 = - L31*B(:,1) - L32*y2
    ! Solve UX(:,:,j) = y by back substitution
    X(:,1,3) = y3 / U33
    X(:,1,2) = (y2 - U23*X(:,1,3)) / U22
    X(:,1,1) = (B(:,1) - A(:,2,1)*X(:,1,2) &
         &          - A(:,3,1)*X(:,1,3)) / A(:,1,1)

    ! Solve X(2,:) = A^-T (  0   )
    !                     ( B(2) )
    !                     (  0   )
    ! Solve Ly = B(:,:,j) by forward substitution
    ! y1 = 0
    ! y2 = B(:,2)
    y3 = - L32*B(:,2)
    ! Solve UX(:,:,j) = y by back substitution
    X(:,2,3) = y3 / U33
    X(:,2,2) = (B(:,2) - U23*X(:,2,3)) / U22
    X(:,2,1) = (-A(:,2,1)*X(:,2,2) &
         &           -A(:,3,1)*X(:,2,3)) / A(:,1,1)

    ! Solve X(3,:) = A^-T (  0   )
    !                     (  0   )
    !                     ( B(3) )
    ! Solve Ly = B(:,:,j) by forward substitution
    ! y1 = 0
    ! y2 = 0
    ! y3 = B(:,3)
    ! Solve UX(:,:,j) = y by back substitution
    X(:,3,3) = B(:,3) / U33
    X(:,3,2) = -U23*X(:,3,3) / U22
    X(:,3,1) = (-A(:,2,1)*X(:,3,2) &
         &          - A(:,3,1)*X(:,3,3)) / A(:,1,1)

  end subroutine diag_mat_right_divide_3


  !---------------------------------------------------------------------
  ! Treat A as n m-by-m matrices and return the LU factorization of A
  ! compressed into a single matrice (with L below the diagonal and U
  ! on and above the diagonal; the diagonal elements of L are 1). No
  ! pivoting is performed.
  pure subroutine lu_factorization(n, iend, m, A, LU)
    integer,    intent(in)  :: n, m, iend
    real(jprb), intent(in)  :: A(:,:,:)
    real(jprb), intent(out) :: LU(iend,m,m)

    real(jprb) :: s(iend)
    integer    :: j1, j2, j3

    ! This routine is adapted from an in-place one, so we first copy
    ! the input into the output.
    LU(1:iend,1:m,1:m) = A(1:iend,1:m,1:m)

    do j2 = 1, m
      do j1 = 1, j2-1
        s = LU(1:iend,j1,j2)
        do j3 = 1, j1-1
          s = s - LU(1:iend,j1,j3) * LU(1:iend,j3,j2)
        end do
        LU(1:iend,j1,j2) = s
      end do
      do j1 = j2, m
        s = LU(1:iend,j1,j2)
        do j3 = 1, j2-1
          s = s - LU(1:iend,j1,j3) * LU(1:iend,j3,j2)
        end do
        LU(1:iend,j1,j2) = s
      end do
      if (j2 /= m) then
        s = 1.0_jprb / LU(1:iend,j2,j2)
        do j1 = j2+1, m
          LU(1:iend,j1,j2) = LU(1:iend,j1,j2) * s
        end do
      end if
    end do

  end subroutine lu_factorization


  !---------------------------------------------------------------------
  ! Treat LU as an LU-factorization of an original matrix A, and
  ! return x where Ax=b. LU consists of n m-by-m matrices and b as n
  ! m-element vectors.
  pure subroutine lu_substitution(n,iend,m,LU,b,x)
    ! CHECK: dimensions should be ":"?
    integer,    intent(in) :: n, m, iend
    real(jprb), intent(in) :: LU(iend,m,m)
    real(jprb), intent(in) :: b(:,:)
    real(jprb), intent(out):: x(iend,m)

    integer :: j1, j2

    x(1:iend,1:m) = b(1:iend,1:m)

    ! First solve Ly=b
    do j2 = 2, m
      do j1 = 1, j2-1
        x(1:iend,j2) = x(1:iend,j2) - x(1:iend,j1)*LU(1:iend,j2,j1)
      end do
    end do
    ! Now solve Ux=y
    do j2 = m, 1, -1
      do j1 = j2+1, m
        x(1:iend,j2) = x(1:iend,j2) - x(1:iend,j1)*LU(1:iend,j2,j1)
      end do
      x(1:iend,j2) = x(1:iend,j2) / LU(1:iend,j2,j2)
    end do

  end subroutine lu_substitution


  !---------------------------------------------------------------------
  ! Return matrix X where AX=B. LU, A, X, B all consist of n m-by-m
  ! matrices.
  pure subroutine solve_mat_n(n,iend,m,A,B,X)
    integer,    intent(in) :: n, m, iend
    real(jprb), intent(in) :: A(:,:,:)
    real(jprb), intent(in) :: B(:,:,:)
    real(jprb), intent(out):: X(iend,m,m)

    real(jprb) :: LU(iend,m,m)

    integer :: j

    call lu_factorization(n,iend,m,A,LU)

    do j = 1, m
      call lu_substitution(n,iend,m,LU,B(1:,1:m,j),X(1:iend,1:m,j))
!      call lu_substitution(n,iend,m,LU,B(1:n,1:m,j),X(1:iend,1:m,j))
    end do

  end subroutine solve_mat_n


  !---------------------------------------------------------------------
  ! Solve Ax=b, where A consists of n m-by-m matrices and x and b
  ! consist of n m-element vectors. For m=2 or m=3, this function
  ! calls optimized versions, otherwise it uses general LU
  ! decomposition without pivoting.
  function solve_vec(n,iend,m,A,b)

    use yomhook, only : lhook, dr_hook

    integer,    intent(in) :: n, m, iend
    real(jprb), intent(in) :: A(:,:,:)
    real(jprb), intent(in) :: b(:,:)

    real(jprb)             :: solve_vec(iend,m)
    real(jprb)             :: LU(iend,m,m)
    real(jprb)             :: hook_handle

    if (lhook) call dr_hook('radiation_matrix:solve_vec',0,hook_handle)

    if (m == 2) then
      call solve_vec_2(n,iend,A,b,solve_vec)
    elseif (m == 3) then
      call solve_vec_3(n,iend,A,b,solve_vec)
    else
      call lu_factorization(n,iend,m,A,LU)
      call lu_substitution(n,iend,m,LU,b,solve_vec)
    end if

    if (lhook) call dr_hook('radiation_matrix:solve_vec',1,hook_handle)

  end function solve_vec


  !---------------------------------------------------------------------
  ! Solve AX=B, where A, X and B consist of n m-by-m matrices. For m=2
  ! or m=3, this function calls optimized versions, otherwise it uses
  ! general LU decomposition without pivoting.
  function solve_mat(n,iend,m,A,B)

    use yomhook, only : lhook, dr_hook

    integer,    intent(in)  :: n, m, iend
    real(jprb), intent(in)  :: A(:,:,:)
    real(jprb), intent(in)  :: B(:,:,:)

    real(jprb)              :: solve_mat(iend,m,m)
    real(jprb)              :: hook_handle

    if (lhook) call dr_hook('radiation_matrix:solve_mat',0,hook_handle)

    if (m == 2) then
      call solve_mat_2(n,iend,A,B,solve_mat)
    elseif (m == 3) then
      call solve_mat_3(n,iend,A,B,solve_mat)
    else
      call solve_mat_n(n,iend,m,A,B,solve_mat)
    end if

    if (lhook) call dr_hook('radiation_matrix:solve_mat',1,hook_handle)

  end function solve_mat

  pure subroutine mat_square_sw(m,iend,A,C)

    integer,    intent(in)                        :: m,iend
    real(jprb), intent(in),  dimension(iend,m,m)  :: A
    real(jprb), intent(out), dimension(iend,m,m)  :: C
    integer    :: j1, j2, j3
    integer    :: mblock, m2block

    ! Array-wise assignment
    ! C = 0.0_jprb

    ! Matrix has a sparsity pattern
    !     (C D E)
    ! A = (F G H)
    !     (0 0 I)

    mblock = m/3       ! 
    m2block = 2*mblock ! 

    ! Do the top-left (C, D, F, G)
    do j2 = 1,m2block  !    1,6 
      C(:,:,j2) = 0.0_jprb
      do j1 = 1,m2block !   1,6
        do j3 = 1,m2block ! 1,6
          C(:,j1,j2) = C(:,j1,j2) + A(:,j1,j3)*A(:,j3,j2)
        end do
        ! using sum was faster on GCC+AMD Zen platform but not on Intel
        ! C(:,j1,j2) = sum(A(:,j1,1:m2block)*A(:,1:m2block,j2),2)
      end do
    end do

    do j2 = m2block+1,m  ! 7,9
      C(:,:,j2) = 0.0_jprb
      ! Do the top-right (E & H)
      do j1 = 1,m2block  ! 1,6
        do j3 = 1,m      ! 1,9
          C(:,j1,j2) = C(:,j1,j2) + A(:,j1,j3)*A(:,j3,j2)
        end do
      end do

      ! Do the bottom-right (I)
      do j1 = m2block+1,m   ! 7,9
        do j3 = m2block+1,m ! 7,9
          C(:,j1,j2) =  C(:,j1,j2) + A(:,j1,j3)*A(:,j3,j2)
        end do
      end do
    end do

  end subroutine mat_square_sw

  pure subroutine mat_square_sw_9(iend,A,C)

    integer,    intent(in)                        :: iend
    real(jprb), intent(in),  dimension(iend,9,9)  :: A
    real(jprb), intent(out), dimension(iend,9,9)  :: C
    integer    :: j1, j2, j3, jg

    ! First input matrix has pattern

    !     (C    D     E)
    ! A = (F=-D G=-C  H)
    !     (0    0     I)

    ! As a result, output and subsequent input matrices have pattern
    
    !     (C    D    E)
    ! A = (F=D  G=C  H)
    !     (0    0    I)

    do j2 = 1,3  !    1,3 
      do j1 = 1,6 !   1,6   C, F
        ! do j3 = 1,6 ! 1,6
        !   C(:,j1,j2) = C(:,j1,j2) + A(:,j1,j3)*A(:,j3,j2)
        ! end do
        ! Further speedup: flatten last loop, only one write SIMD instruction
        C(:,j1,j2) = A(:,j1,1)*A(:,1,j2) + A(:,j1,2)*A(:,2,j2) &
        &          + A(:,j1,3)*A(:,3,j2) + A(:,j1,4)*A(:,4,j2) &
        &          + A(:,j1,5)*A(:,5,j2) + A(:,j1,6)*A(:,6,j2)
      end do
    end do
    ! D,G
    C(:,1:3,4:6) = C(:,4:6,1:3) ! D = F
    C(:,4:6,4:6) = C(:,1:3,1:3) ! G = C

    ! Lower left corner with zeros
    C(:,7:9,1:6) = 0.0_jprb

    do j2 = 7,9  ! 7,9
      ! Do the top-right (E & H)
      do j1 = 1,6  ! 1,6
        C(:,j1,j2) = A(:,j1,1)*A(:,1,j2) + A(:,j1,2)*A(:,2,j2) &
        &          + A(:,j1,3)*A(:,3,j2) + A(:,j1,4)*A(:,4,j2) &
        &          + A(:,j1,5)*A(:,5,j2) + A(:,j1,6)*A(:,6,j2) &
        &          + A(:,j1,7)*A(:,7,j2) + A(:,j1,8)*A(:,8,j2) &
        &          + A(:,j1,9)*A(:,9,j2)
      end do
      ! Do the bottom-right (I)
      do j1 = 7,9   ! 7,9
        C(:,j1,j2) = A(:,j1,7)*A(:,7,j2) + A(:,j1,8)*A(:,8,j2) + A(:,j1,9)*A(:,9,j2)
      end do
    end do

  end subroutine mat_square_sw_9

  pure subroutine mat_x_mat_sw_9(iend,A,B,C)

    integer,    intent(in)                      :: iend
    real(jprb), intent(in), dimension(iend,9,9) :: A, B
    real(jprb), intent(out),dimension(iend,9,9) :: C
    integer    :: j1, j2, j3

    do j2 = 1,3  !    1,3 
      do j1 = 1,6 !   1,6   C, F
        C(:,j1,j2) = A(:,j1,1)*B(:,1,j2) + A(:,j1,2)*B(:,2,j2) &
        &          + A(:,j1,3)*B(:,3,j2) + A(:,j1,4)*B(:,4,j2) &
        &          + A(:,j1,5)*B(:,5,j2) + A(:,j1,6)*B(:,6,j2)
      end do
    end do
    ! D,G
    C(:,1:3,4:6) = C(:,4:6,1:3) ! D = F
    C(:,4:6,4:6) = C(:,1:3,1:3) ! G = C

    ! Lower left corner with zeros
    C(:,7:9,1:6) = 0.0_jprb

    do j2 = 7,9  ! 7,9
      ! Do the top-right (E & H)
      do j1 = 1,6  ! 1,6
        C(:,j1,j2) = A(:,j1,1)*B(:,1,j2) + A(:,j1,2)*B(:,2,j2) &
        &          + A(:,j1,3)*B(:,3,j2) + A(:,j1,4)*B(:,4,j2) &
        &          + A(:,j1,5)*B(:,5,j2) + A(:,j1,6)*B(:,6,j2) &
        &          + A(:,j1,7)*B(:,7,j2) + A(:,j1,8)*B(:,8,j2) &
        &          + A(:,j1,9)*B(:,9,j2)
      end do

      ! Do the bottom-right (I)
      do j1 = 7,9   ! 7,9
        C(:,j1,j2) = A(:,j1,7)*B(:,7,j2) + A(:,j1,8)*B(:,8,j2) + A(:,j1,9)*B(:,9,j2)
      end do
    end do

  end subroutine mat_x_mat_sw_9

  !---------------------------------------------------------------------
  ! Solve AX=B, where A, X and B consist of iend m-by-m matrices
  ! Overwrite B with X. A is corrupted
  pure subroutine solve_mat_sw_9(iend,A,B)
    integer,    intent(in) :: iend
    real(jprb), intent(inout) :: A(iend,9,9) ! A=LU is corrupted
    real(jprb), intent(inout) :: B(iend,9,9) ! X = B, both input and output
    ! real(jprb), intent(out):: X(iend,m,m)

    integer :: j,j1, j2, j3, mblock, m2block,m
    real(jprb) :: s(iend)

    !     (C   D  E)
    ! A = (-D -C  H)
    !     (0   0  I)

    m = 9
    mblock = 3 !m/3       
    m2block = 6 !2*mblock 

    ! factorization of A into LU

    ! First do columns 1-6, for which only rows 1-6 have non-negative entries
    do j2 = 1, m2block
      do j1 = 1, j2-1
        do j3 = 1, j1-1
          A(:,j1,j2) = A(:,j1,j2)- A(:,j1,j3) * A(:,j3,j2)
        end do
      end do
      do j1 = j2, m2block
        do j3 = 1, j2-1
          A(:,j1,j2) = A(:,j1,j2) - A(:,j1,j3) * A(:,j3,j2)
        end do
      end do
      s = 1.0_jprb / A(:,j2,j2)
      do j1 = j2+1, m2block
        A(:,j1,j2) = A(:,j1,j2) * s
      end do
    end do

    ! Remaining columns
    do j2 = m2block+1, m
      do j1 = 1, j2-1
        do j3 = 1, j1-1
          A(:,j1,j2) = A(:,j1,j2) - A(:,j1,j3) * A(:,j3,j2)
        end do
      end do
      do j1 = j2, m
        do j3 = 1, j2-1
          A(:,j1,j2)= A(:,j1,j2) - A(:,j1,j3) * A(:,j3,j2)
        end do
      end do
      if (j2 /= m) then
        s = 1.0_jprb / A(:,j2,j2)
        do j1 = j2+1, m
          A(:,j1,j2) = A(:,j1,j2) * s
        end do
      end if
    end do

    !---------------------------------------------------------------------
    ! Treat LU as an LU-factorization of an original matrix A, and
    ! return x where Ax=b. LU consists of n m-by-m matrices and b as n
    ! m-element vectors.
    ! Here B is both input b and output x, and A has been LU factorized, combining L and U
    ! into one matrix where the diagonal is the diagonal of U (L has ones in the diagonal)

    ! A and B both have following structure:
    !     (C   D  E)
    !     (F   G  H)
    !     (0   0  I)

    ! Separate j3 (columns) into two regions to avoid redundant operations with zero
    do j3 = 1,m2block ! in this region B(:,7:9),A(:,7:9) are 0
      ! First solve Ly=b
      do j2 = 2, m2block
        do j1 = 1, j2-1
          B(:,j2,j3) = B(:,j2,j3) - B(:,j1,j3)*A(:,j2,j1)
        end do
        ! No division because diagonal of L is unity
      end do
      ! Now solve Ux=y
      do j2 = m2block, 1, -1
        do j1 = j2+1, m2block
          B(:,j2,j3) = ( B(:,j2,j3) - B(:,j1,j3)*A(:,j2,j1) )
        end do
        B(:,j2,j3) = B(:,j2,j3) / A(:,j2,j2) ! Divide by diagonal of A=U
      end do
    end do

    do j3 = m2block+1,m ! columns 7-9: here B has nonzero values for all rows, but A doesn't 
      ! First solve Ly=b
      ! do j2 = 2, m
      do j2 = 2, m2block 
        do j1 = 1, j2-1
          B(:,j2,j3) = B(:,j2,j3) - B(:,j1,j3)*A(:,j2,j1)
        end do
        ! No division because diagonal of L is unity
      end do
      ! When j2 = 7, the A terms are all 0, because A(7,1:6)=0
      ! when j2 = 8, only the last j1 has nonzero A
      ! When j2 = 9, two last j1 are nonzero
      B(:,8,j3) = B(:,8,j3) - B(:,7,j3)*A(:,8,7)   ! j2 = 8
      B(:,9,j3) = B(:,9,j3) - B(:,8,j3)*A(:,9,8) - B(:,7,j3)*A(:,9,7) ! j2 = 9

      ! Now solve Ux=y
      do j2 = m, 1, -1
        do j1 = j2+1, m
          B(:,j2,j3) = B(:,j2,j3) - B(:,j1,j3)*A(:,j2,j1)
        end do
        B(:,j2,j3) = B(:,j2,j3) / A(:,j2,j2) ! Divide by diagonal of A=U
      end do
    end do

  end subroutine solve_mat_sw_9


  ! --- MATRIX EXPONENTIATION ---
  !---------------------------------------------------------------------
  ! Perform matrix exponential of n m-by-m matrices stored in A (where
  ! index n varies fastest) using the Higham scaling and squaring
  ! method. The result is placed in A. This routine is intended for
  ! speed so is accurate only to single precision.  For simplicity and
  ! to aid vectorization, the Pade approximant of order 7 is used for
  ! all input matrices, perhaps leading to a few too many
  ! multiplications for matrices with a small norm.
  subroutine expm(n,iend,m,A,i_matrix_pattern)

    use yomhook, only : lhook, dr_hook

    integer,    intent(in)      :: n, m, iend
    real(jprb), intent(inout)   :: A(n,m,m)
    integer,    intent(in)      :: i_matrix_pattern

    real(jprb), parameter :: theta(3) = (/4.258730016922831e-01_jprb, &
         &                                1.880152677804762e+00_jprb, &
         &                                3.925724783138660e+00_jprb/) 
    real(jprb), parameter :: c(8) = (/17297280.0_jprb, 8648640.0_jprb, &
         &                1995840.0_jprb, 277200.0_jprb, 25200.0_jprb, &
         &                1512.0_jprb, 56.0_jprb, 1.0_jprb/)

    real(jprb), dimension(iend,m,m) :: A2, A4, A6
    real(jprb), dimension(iend,m,m) :: U, V

    real(jprb) :: normA(iend), sum_column(iend)

    integer    :: j1, j2, j3
    real(jprb) :: frac(iend)
    integer    :: expo(iend)
    real(jprb) :: scaling(iend)

    real(jprb) :: hook_handle

    if (lhook) call dr_hook('radiation_matrix:expm',0,hook_handle)

    normA = 0.0_jprb

    ! Compute the 1-norms of A
    do j3 = 1,m
      sum_column(:) = 0.0_jprb
      do j2 = 1,m
        do j1 = 1,iend
          sum_column(j1) = sum_column(j1) + abs(A(j1,j2,j3))
        end do
      end do
      do j1 = 1,iend
        if (sum_column(j1) > normA(j1)) then
          normA(j1) = sum_column(j1)
        end if
      end do
    end do

    frac = fraction(normA/theta(3))
    expo = exponent(normA/theta(3))
    where (frac == 0.5_jprb)
      expo = expo - 1
    end where

    where (expo < 0)
      expo = 0
    end where

    ! Scale the input matrices by a power of 2
    scaling = 2.0_jprb**(-expo)
    do j3 = 1,m
      do j2 = 1,m
        A(1:iend,j2,j3) = A(1:iend,j2,j3) * scaling
      end do
    end do
    ! Pade approximant of degree 7
    A2 = mat_x_mat(n,iend,m,A, A, i_matrix_pattern)
    A4 = mat_x_mat(n,iend,m,A2,A2,i_matrix_pattern)
    A6 = mat_x_mat(n,iend,m,A2,A4,i_matrix_pattern)

    V = c(8)*A6 + c(6)*A4 + c(4)*A2
    do j3 = 1,m
      V(:,j3,j3) = V(:,j3,j3) + c(2)
    end do
    U = mat_x_mat(n,iend,m,A,V,i_matrix_pattern)
    V = c(7)*A6 + c(5)*A4 + c(3)*A2
    ! Add a multiple of the identity matrix
    do j3 = 1,m
      V(:,j3,j3) = V(:,j3,j3) + c(1)
    end do

    V = V-U
    U = 2.0_jprb*U
    A(1:iend,1:m,1:m) = solve_mat(n,iend,m,V,U)

    ! Add the identity matrix
    do j3 = 1,m
      A(1:iend,j3,j3) = A(1:iend,j3,j3) + 1.0_jprb
    end do

    ! Loop through the matrices
    do j1 = 1,iend
      if (expo(j1) > 0) then
        ! Square matrix j1 expo(j1) times          
        A(j1,:,:) = repeated_square(m,A(j1,:,:),expo(j1),i_matrix_pattern)
      end if
    end do

    if (lhook) call dr_hook('radiation_matrix:expm',1,hook_handle)

  end subroutine expm

  !---------------------------------------------------------------------
  ! Like expm, but optimized for the shortwave, which has
  ! a special matrix structure with zeros and repeated elements. 
  ! Further assumes nreg = 3  =>  m = 9
  subroutine expm_opt(iend,A)

    use yomhook, only : lhook, dr_hook

    integer,    intent(in)      :: iend
    real(jprb), intent(inout)   :: A(iend,9,9)

    real(jprb), parameter :: theta(3) = (/4.258730016922831e-01_jprb, &
         &                                1.880152677804762e+00_jprb, &
         &                                3.925724783138660e+00_jprb/) 
    real(jprb), parameter :: c(8) = (/17297280.0_jprb, 8648640.0_jprb, &
         &                1995840.0_jprb, 277200.0_jprb, 25200.0_jprb, &
         &                1512.0_jprb, 56.0_jprb, 1.0_jprb/)

    real(jprb), dimension(iend,9,9) :: A2, A4, A6
    real(jprb), dimension(iend,9,9) :: U, V
    real(jprb), dimension(9,9)      :: temp_in, temp_out
    real(jprb) :: normA(iend), sum_column(iend)

    integer    :: j1, j2, j3, j4,minexpo, nrepeat
    real(jprb) :: frac(iend)
    integer    :: expo(iend)
    real(jprb) :: scaling(iend)
    real(jprb) :: hook_handle

    if (lhook) call dr_hook('radiation_matrix:expm_opt',0,hook_handle)

    normA = 0.0_jprb

    ! Compute the 1-norms of A
    do j3 = 1,9
      sum_column(:) = 0.0_jprb
      do j2 = 1,9
        do j1 = 1,iend
          sum_column(j1) = sum_column(j1) + abs(A(j1,j2,j3))
        end do
      end do
      do j1 = 1,iend
        if (sum_column(j1) > normA(j1)) then
          normA(j1) = sum_column(j1)
        end if
      end do
    end do

    frac = fraction(normA/theta(3))
    expo = exponent(normA/theta(3))
    where (frac == 0.5_jprb)
      expo = expo - 1
    end where

    where (expo < 0)
      expo = 0
    end where
    
    minexpo = minval(expo)
    ! Scale the input matrices by a power of 2
    scaling = 2.0_jprb**(-expo)
    do j3 = 1,9
      do j2 = 1,9
        A(:,j2,j3) = A(:,j2,j3) * scaling
      end do
    end do
    ! Pade approximant of degree 7
    call mat_square_sw_9(iend,A,A2)    ! These matrices have zeroes in the lower left corner AND repeated elements
    call mat_square_sw_9(iend,A2,A4)

    call mat_x_mat_sw_9(iend,A2,A4,A6) ! These matrices have zeroes in the lower left corner AND repeated elements

    V = c(8)*A6 + c(6)*A4 + c(4)*A2
    do j3 = 1,9
      V(:,j3,j3) = V(:,j3,j3) + c(2)
    end do

    U = mat_x_mat(iend,iend,9,A,V,IMatrixPatternShortwave)

    V = c(7)*A6 + c(5)*A4 + c(3)*A2
    ! Add a multiple of the identity matrix
    do j3 = 1,9
      V(:,j3,j3) = V(:,j3,j3) + c(1)
    end do

    V = V-U
    A = 2.0_jprb*U

    ! A = solve_mat(n,iend,m,V,U)
    call solve_mat_sw_9(iend,V,A)

    ! Add the identity matrix
    do j3 = 1,9
      A(:,j3,j3) = A(:,j3,j3) + 1.0_jprb
    end do

    ! Loop through the matrices    
    ! To improve efficiency, square all matrices with the minimum expo first, and then square individual matrices as needed
    do j1 = 1,minexpo
      call mat_square_sw(9,iend,A,A2)
      A = A2
    end do

    do j1 = 1,iend
      if (expo(j1) > minexpo) then
        nrepeat = expo(j1)-minexpo
        !Square matrix nrepeat times 
        temp_in =  A(j1,:,:)
        call repeated_square_sw_9(nrepeat,temp_in,temp_out)
        A(j1,:,:) = temp_out
      end if
    end do
 
    if (lhook) call dr_hook('radiation_matrix:expm_opt',1,hook_handle)

  end subroutine expm_opt


  !---------------------------------------------------------------------
  ! Return the matrix exponential of n 2x2 matrices representing
  ! conservative exchange between SPARTACUS regions, where the
  ! matrices have the structure
  !   (-a   b)
  !   ( a  -b)
  ! and a and b are assumed to be positive or zero.  The solution uses
  ! Putzer's algorithm - see the appendix of Hogan et al. (GMD 2018)
  subroutine fast_expm_exchange_2(n,a,b,R)

    use yomhook, only : lhook, dr_hook

    integer,                      intent(in)  :: n
    real(jprb), dimension(n),     intent(in)  :: a, b
    real(jprb), dimension(n,2,2), intent(out) :: R

    real(jprb), dimension(n) :: factor

    real(jprb) :: hook_handle

    if (lhook) call dr_hook('radiation_matrix:fast_expm_exchange_2',0,hook_handle)

    ! Security to ensure that if a==b==0 then the identity matrix is returned
    factor = (1.0_jprb - exp(-(a(:)+b(:))))/max(1.0e-12_jprb,a(:)+b(:))

    R(:,1,1) = 1.0_jprb - factor*a(:)
    R(:,2,1) = factor*a(:)
    R(:,1,2) = factor*b(:)
    R(:,2,2) = 1.0_jprb - factor*b(:)

    if (lhook) call dr_hook('radiation_matrix:fast_expm_exchange_2',1,hook_handle)

  end subroutine fast_expm_exchange_2


  !---------------------------------------------------------------------
  ! Return the matrix exponential of n 3x3 matrices representing
  ! conservative exchange between SPARTACUS regions, where the
  ! matrices have the structure
  !   (-a   b   0)
  !   ( a -b-c  d)
  !   ( 0   c  -d)
  ! and a-d are assumed to be positive or zero.  The solution uses the
  ! diagonalization method and is a slight generalization of the
  ! solution provided in the appendix of Hogan et al. (GMD 2018),
  ! which assumed c==d.
  subroutine fast_expm_exchange_3(n,a,b,c,d,R)

    use yomhook, only : lhook, dr_hook

    real(jprb), parameter :: my_epsilon = 1.0e-12_jprb

    integer,                      intent(in)  :: n
    real(jprb), dimension(n),     intent(in)  :: a, b, c, d
    real(jprb), dimension(n,3,3), intent(out) :: R

    ! Eigenvectors
    real(jprb), dimension(n,3,3) :: V

    ! Non-zero Eigenvalues
    real(jprb), dimension(n) :: lambda1, lambda2

    ! Diagonal matrix of the exponential of the eigenvalues
    real(jprb), dimension(n,3) :: diag

    ! Result of diag right-divided by V
    real(jprb), dimension(n,3,3) :: X

    ! Intermediate arrays
    real(jprb), dimension(n) :: y2, y3

    real(jprb), dimension(n) :: L21, L31, L32
    real(jprb), dimension(n) :: U22, U23, U33

    integer :: j1, j2

    real(jprb) :: hook_handle

    if (lhook) call dr_hook('radiation_matrix:fast_expm_exchange_3',0,hook_handle)

    ! Eigenvalues
    y2 = 0.5_jprb * (a(:)+b(:)+c(:)+d(:))
    y3 = sqrt(y2*y2 - (a(:)*c(:) + a(:)*d(:) + b(:)*d(:)))
    lambda1 = -y2 + y3
    lambda2 = -y2 - y3

    ! Eigenvectors, with securities such taht if a--d are all zero
    ! then V is non-singular and the identity matrix is returned in R;
    ! note that lambdaX is typically negative so we need a
    ! sign-preserving security
    V(:,1,1) = max(my_epsilon, b(:)) &
         &  / sign(max(my_epsilon, abs(a(:) + lambda1)), a(:) + lambda1)
    V(:,1,2) = b(:) &
         &  / sign(max(my_epsilon, abs(a(:) + lambda2)), a(:) + lambda2)
    V(:,1,3) = b(:) / max(my_epsilon, a(:))
    V(:,2,:) = 1.0_jprb
    V(:,3,1) = c(:) &
         &  / sign(max(my_epsilon, abs(d(:) + lambda1)), d(:) + lambda1)
    V(:,3,2) = c(:) &
         &  / sign(max(my_epsilon, abs(d(:) + lambda2)), d(:) + lambda2)
    V(:,3,3) = max(my_epsilon, c(:)) / max(my_epsilon, d(:))
    
    diag(:,1) = exp(lambda1)
    diag(:,2) = exp(lambda2)
    diag(:,3) = 1.0_jprb

    ! ------ Compute X = diag * V^-1 ---------
    !  call diag_mat_right_divide_3(n,V,diag,X)

      !    associate (U11 => V(:,1,1), U12 => V(:,1,2), U13 => V(1,3))
    ! LU decomposition of the *transpose* of V:
    !       ( 1        )   (U11 U12 U13)
    ! V^T = (L21  1    ) * (    U22 U23)
    !       (L31 L32  1)   (        U33)
    L21 = V(:,1,2) / V(:,1,1)
    L31 = V(:,1,3) / V(:,1,1)
    U22 = V(:,2,2) - L21*V(:,2,1)
    U23 = V(:,3,2) - L21*V(:,3,1)
    L32 =(V(:,2,3) - L31*V(:,2,1)) / U22
    U33 = V(:,3,3) - L31*V(:,3,1) - L32*U23

    ! Solve X(1,:) = V^-T ( diag(1) )
    !                     (  0   )
    !                     (  0   )
    ! Solve Ly = diag(:,:,j) by forward substitution
    ! y1 = diag(:,1)
    y2 = - L21*diag(:,1)
    y3 = - L31*diag(:,1) - L32*y2
    ! Solve UX(:,:,j) = y by back substitution
    X(:,1,3) = y3 / U33
    X(:,1,2) = (y2 - U23*X(:,1,3)) / U22
    X(:,1,1) = (diag(:,1) - V(:,2,1)*X(:,1,2) &
         &          - V(:,3,1)*X(:,1,3)) / V(:,1,1)

    ! Solve X(2,:) = V^-T (  0   )
    !                     ( diag(2) )
    !                     (  0   )
    ! Solve Ly = diag(:,:,j) by forward substitution
    ! y1 = 0
    ! y2 = diag(:,2)
    y3 = - L32*diag(:,2)
    ! Solve UX(:,:,j) = y by back substitution
    X(:,2,3) = y3 / U33
    X(:,2,2) = (diag(:,2) - U23*X(:,2,3)) / U22
    X(:,2,1) = (-V(:,2,1)*X(:,2,2) &
         &           -V(:,3,1)*X(:,2,3)) / V(:,1,1)

    ! Solve X(3,:) = V^-T (  0   )
    !                     (  0   )
    !                     ( diag(3) )
    ! Solve Ly = diag(:,:,j) by forward substitution
    ! y1 = 0
    ! y2 = 0
    ! y3 = diag(:,3)
    ! Solve UX(:,:,j) = y by back substitution
    X(:,3,3) = diag(:,3) / U33
    X(:,3,2) = -U23*X(:,3,3) / U22
    X(:,3,1) = (-V(:,2,1)*X(:,3,2) &
         &          - V(:,3,1)*X(:,3,3)) / V(:,1,1)

    ! Compute V * X
    do j1 = 1,3
      do j2 = 1,3
        R(:,j2,j1) = V(:,j2,1)*X(:,1,j1) &
             &          + V(:,j2,2)*X(:,2,j1) &
             &          + V(:,j2,3)*X(:,3,j1)
      end do
    end do

    if (lhook) call dr_hook('radiation_matrix:fast_expm_exchange_3',1,hook_handle)

  end subroutine fast_expm_exchange_3

end module radiation_matrix
