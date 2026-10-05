module print_matrix_mod
contains
  subroutine print_vector(vec, name, unit)
    use parkind1, only : jprb
    real(jprb),   intent(in) :: vec(:)
    character(*), intent(in), optional :: name
    integer,      intent(in), optional :: unit
    integer    :: i
    integer    :: unit_local
    if (present(unit)) then
      unit_local = unit
    else
      unit_local = 6
    end if
    if (present(name)) then
      write(unit_local,'(a,a,$)') name, '=['
    end if
    do i = 1,size(vec,1)
       write(unit_local,'(f16.8,$)') vec(i)
    end do
    if (present(name)) then
      write(unit_local,'(a)') ']'
    else
      write(unit_local,'(x)')
    end if

  end subroutine print_vector

  subroutine print_matrix(mat, name, unit)
    use parkind1, only : jprb
    real(jprb),   intent(in) :: mat(:,:)
    character(*), intent(in), optional :: name
    integer,      intent(in), optional :: unit
    integer    :: i, j
    integer    :: unit_local
    if (present(unit)) then
      unit_local = unit
    else
      unit_local = 6
    end if

    if (present(name)) then
      write(unit_local,'(a,a,$)') name, '=['
    end if
    do i = 1,size(mat,1)
       do j = 1,size(mat,2)
          write(unit_local,'(f16.8,$)') mat(i,j)
       end do
       if (present(name) .and. i == size(mat,1)) then
         write(unit_local,'(a)') ']'
       else
         write(unit_local,'(x)')
       end if
    end do
  end subroutine print_matrix
  
  subroutine print_array3(name, mat)
    use parkind1, only : jprb
    character(*), intent(in) :: name
    real(jprb),   intent(in) :: mat(:,:,:)
    integer    :: i, j, k
    write(6,'(a,a)') name, '='
    do k = 1,size(mat,3)
      do i = 1,size(mat,1)
        do j = 1,size(mat,2)
             write(6,'(f16.8,$)') mat(i,j,k)
          end do
          write(6,'(x)')
       end do
       write(6,'(x)')
    end do
  end subroutine print_array3
  

end module print_matrix_mod
