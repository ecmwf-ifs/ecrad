program test_gaussian_quadrature

  use parkind1, only : jprb
  use easy_netcdf
  use radiation_gaussian_quadrature

  integer, parameter :: norder = 14
  integer, parameter :: max_order = 64
  
  integer :: orders(norder)
  real(jprb) :: mu(max_order,norder), wt(max_order,norder)
  
  integer :: jord, jquad
  integer :: jacobi_power = 5
  integer :: quadrature_types(3) = [IQuadratureLegendre, IQuadratureLaguerre, IQuadratureJacobi]

  type(netcdf_file) :: file
  character(len=512) :: quadrature_name, quadrature_name_uc
  
  orders = [1, 2, 3, 4, 5, 6, 7, 8, 10, 12, 14, 16, 32, max_order]

  mu = -1.0_jprb
  wt = 0.0_jprb
  
  do jquad = 1,3
    do jord = 1,norder
      call gaussian_quadrature(quadrature_types(jquad), orders(jord), &
           &  mu(1:orders(jord),jord), wt(1:orders(jord),jord), jacobi_power)
    end do

    if (.not. quadrature_types(jquad) == IQuadratureJacobi) then
      write(quadrature_name,'(a,a)') 'gauss-', trim(QuadratureNameLC(quadrature_types(jquad)))
      write(quadrature_name_uc,'(a,a)') 'Gauss-', trim(QuadratureName(quadrature_types(jquad)))
    else
      write(quadrature_name,'(a,a,a,i0)') 'gauss-', &
           &  trim(QuadratureNameLC(quadrature_types(jquad))), '-', jacobi_power
      write(quadrature_name_uc,'(a,a,a,i0)') 'Gauss-', &
           &  trim(QuadratureName(quadrature_types(jquad))), '-', jacobi_power
    end if
    
    call file%create('quadrature_' // trim(quadrature_name) // '.nc', iverbose=3)
    call file%double_precision(.true.)
    call file%define_dimension('order',norder)
    call file%define_dimension('node',max_order)
    call file%define_variable('order',dim1_name='order',long_name='Order of quadrature', &
         &  data_type_name='short')
    call file%define_variable('mu',dim2_name='order',dim1_name='node', &
         &  long_name='Cosine of zenith angle',fill_value=-1.0_jprb)
    call file%define_variable('weight',dim2_name='order',dim1_name='node', &
         &  long_name='Weight', fill_value=0.0_jprb)
    call file%put_global_attributes( &
         &   title_str=trim(quadrature_name_uc)//' quadrature angles for longwave radiative transfer', &
         &   source_str="ecRad test_gaussian_quadrature")
    call file%put_global_attribute('quadrature_method_id',trim(quadrature_name))
    call file%put_global_attribute('quadrature_method',trim(quadrature_name_uc)//' quadrature')
    call file%put('order',orders)
    call file%put('mu',mu)
    call file%put('weight',wt)
    call file%close()

  end do

end program test_gaussian_quadrature
