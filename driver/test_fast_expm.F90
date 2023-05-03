program test_fast_expm
  ! Test the fast matrix exponential exchange routine used for
  ! SPARTACUS shortwave entrapment calculations, but which can fail in
  ! single precision for very specific inputs leading to the matrix
  ! having two repeated eigenvalues
  
  use parkind1, only : jprb

  use radiation_matrix, only : fast_expm_exchange_3

  real(jprb) :: a(1), b(1), c(1), d(1), r(1,3,3)

  a = 9.0408579E-02_jprb
  b = 9.2716664E-07_jprb
  c = 2.2503915E-03_jprb
  d = 8.8152386E-02_jprb

  call fast_expm_exchange_3(1,1,a(1:1),b(1:1),c(1:1),d(1:1),r)

  print *, r
  print *, 0.0000001
  print *, epsilon(1.0_jprb)

end program test_fast_expm
