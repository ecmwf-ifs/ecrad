program test_ecrad3d_radiance_coeffs
  use parkind1, only : jprb
  use ecrad3d_layer_solutions, only : calc_radiance_coeffs_sw

  integer, parameter :: ncol=1, nmaxspec=1, nspec=1
  real(kind=jprb) :: mu0(ncol), mu_sens(ncol), azim_diff(ncol)
  real(kind=jprb) :: od(ncol,nmaxspec), ssa(ncol,nmaxspec), asymmetry(ncol,nmaxspec)
  real(kind=jprb) :: dir_dn_to_rad(ncol,nmaxspec), diff_dn_to_rad(ncol,nmaxspec), diff_up_to_rad(ncol,nmaxspec)
  real(kind=jprb) :: trans_rad(ncol,nmaxspec)

  ssa = 0.5
  od = 0.5
  asymmetry = 0.5
  mu0 = 0.5
  mu_sens = 1.0
  azim_diff = 0.0

  call calc_radiance_coeffs_sw(ncol, nmaxspec, nspec, &
       &      mu0, mu_sens, azim_diff, od, ssa, asymmetry, &
       &      dir_dn_to_rad, diff_dn_to_rad, diff_up_to_rad, trans_rad)

  print *, dir_dn_to_rad, diff_dn_to_rad, diff_up_to_rad, trans_rad  
  
end program test_ecrad3d_radiance_coeffs
