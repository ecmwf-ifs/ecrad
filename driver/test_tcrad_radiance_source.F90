program test_tcrad_radiance_source

  use parkind1, only : jprb
  use tcrad_layer_solutions, only : calc_radiance_trans_source_exact

  integer, parameter :: NSPEC = 9
  integer, parameter :: NLEV = 2
  integer, parameter :: NREG = 3
  
  real(jprb) :: mu, region_fracs(NREG,NLEV)
  real(jprb) :: planck_hl(NSPEC,NLEV+1)
  real(jprb) :: od(NSPEC,NREG,NLEV)
  real(jprb) :: ssa(NSPEC,2:NREG,NLEV)
  real(jprb) :: asymmetry(NSPEC,NLEV)
  real(jprb) :: flux_up_base(NSPEC,NREG,NLEV),flux_dn_top(NSPEC,NREG,NLEV)
  real(jprb) :: transmittance(NSPEC,NREG,NLEV)
  real(jprb) :: source_up(NSPEC,NREG,NLEV)

  integer :: jspec, jlev, jreg
  
  mu = 0.5_jprb
  region_fracs(1,:) = 1.0_jprb/3.0_jprb
  region_fracs(2,:) = 1.0_jprb/3.0_jprb
  region_fracs(3,:) = 1.0_jprb/3.0_jprb
  ssa = 0.000001_jprb
  asymmetry = 0.3_jprb
  planck_hl(:,1)=1.0_jprb
  planck_hl(:,2)=2.0_jprb
  planck_hl(:,3)=4.0_jprb

  do jlev = 1,NLEV
    do jreg = 1,NREG
      od(:,jreg,jlev) = [0.001,0.01,0.1,0.3,1.0,3.0,10.0,100.0,1000.0]
    end do
  end do
  flux_up_base = 50.0_jprb
  flux_dn_top  = 30.0_jprb
      
  call calc_radiance_trans_source_exact(NSPEC, NLEV, NREG, mu, &
       &  region_fracs, planck_hl, od, ssa, asymmetry, &
       &  flux_up_base, flux_dn_top, &
       &  transmittance, source_up=source_up)

  
  print *, 'TRANSMITTANCE'
  do jlev = 1,NLEV
    do jreg = 1,NREG
      write(*,'(9E14.7)') transmittance(:,jreg,jlev)
    end do
    print *,''
  end do

  print *, 'SOURCE UP'
  do jlev = 1,NLEV
    do jreg = 1,NREG
      write(*,'(9E14.7)') source_up(:,jreg,jlev)
    end do
    print *,''
  end do
  
end program test_tcrad_radiance_source
