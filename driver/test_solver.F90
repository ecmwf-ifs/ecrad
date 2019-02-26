program test_solver

  use parkind1,                 only : jprb

  use radiation_config,         only : config_type, IGasModelMonochromatic
  use radiation_single_level,   only : single_level_type
  use radiation_thermodynamics, only : thermodynamics_type
  use radiation_cloud,          only : cloud_type
  use radiation_flux,           only : flux_type
  use radiation_monochromatic,  only : setup_gas_optics
  use radiation_homogeneous_sw, only : solver_homogeneous_sw
  use radiation_cloud_optics,   only : delta_eddington

  implicit none

  integer, parameter :: nlev   = 1
  integer, parameter :: ncol   = 10
  integer, parameter :: n_g_sw = 1

  integer :: istartcol, iendcol, jcol, jlev, jod
  
  type(config_type)         :: config
  type(single_level_type)   :: single_level
  type(thermodynamics_type) :: thermodynamics
  type(cloud_type)          :: cloud
  type(flux_type)           :: flux

  real(jprb), dimension(n_g_sw,nlev,ncol) :: od_sw, ssa_sw, g_sw
  real(jprb), dimension(n_g_sw,nlev,ncol) :: od_sw_cloud, ssa_sw_cloud, &
       &  scat_sw_cloud, g_sw_cloud
  real(jprb), dimension(n_g_sw,ncol) :: incoming_sw

  istartcol = 1
  iendcol   = ncol

  config%i_gas_model = IGasModelMonochromatic
  call config%consolidate()

  call setup_gas_optics(config, trim(config%directory_name))

  call config%define_sw_albedo_intervals(1)

  call flux%allocate(config, istartcol, iendcol, nlev)
  call cloud%allocate(ncol, nlev)
  call thermodynamics%allocate(ncol, nlev)
  call single_level%allocate(ncol, 1, 1, .false.)

  cloud%fraction = 1.0

  do jcol = 1,ncol
    single_level%cos_sza(jcol) = cos((jcol-1)*acos(-1.0)/(2.0_jprb*(ncol-1)))
    incoming_sw(1,jcol) = 1.0_jprb
  end do
  write(*,*) 'cos_sza = ', single_level%cos_sza
  single_level%sw_albedo = 0.08_jprb

  do jlev = 1,nlev+1
    thermodynamics%pressure_hl(:,jlev) = 100000.0_jprb * (jlev-1) / nlev
  end do

  od_sw = 0.0_jprb
  ssa_sw = 1.0_jprb
  g_sw = 0.0_jprb

  do jod = 1,8
    ssa_sw_cloud = 0.999_jprb
    g_sw_cloud = 0.85_jprb
!    g_sw_cloud = 0.0_jprb

    if (jod == 1) then
      od_sw_cloud = 0.0_jprb
    else
!      od_sw_cloud(1,2:,:) = 10.0**(0.5 * (jod-4)) / (nlev-1)
      od_sw_cloud(1,:,:) = 10.0**(0.5 * (jod-4)) / nlev

      scat_sw_cloud = od_sw_cloud * ssa_sw_cloud
      call delta_eddington(od_sw_cloud, scat_sw_cloud, g_sw_cloud)
      where (od_sw_cloud > 0.0) 
        ssa_sw_cloud = scat_sw_cloud / od_sw_cloud
      end where
    end if

    write(*,*) 'Optical depth = ', sum(od_sw_cloud(1,:,1)), ' g=', g_sw_cloud(1,:,1)


    ! Compute fluxes using the homogeneous solver
    call solver_homogeneous_sw(nlev,istartcol,iendcol, &
         &  config, single_level, thermodynamics, cloud, & 
         &  od_sw, ssa_sw, g_sw, od_sw_cloud, ssa_sw_cloud, &
         &  g_sw_cloud, incoming_sw, flux)
    
    write(0,*) flux%sw_up(:,1)
  end do


end program test_solver
