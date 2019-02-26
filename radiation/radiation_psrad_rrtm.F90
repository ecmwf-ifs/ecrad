! radiation_psrad_rrtm.F90 - Interface to PS-Rad implementation of RRTM-G
!
! Copyright (C) 2014-2016 ECMWF
!
! Author:  Robin Hogan
! Email:   r.j.hogan@ecmwf.int
! License: see the COPYING file for details
!

module radiation_psrad_rrtm

  implicit none

  public  :: setup_gas_optics, gas_optics, setup_cloud_optics, cloud_optics
  private :: planckFunction

contains
  !---------------------------------------------------------------------
  subroutine setup_gas_optics(config, directory)

    use mo_lrtm_setup, only : setup_lrtm
    use mo_srtm_setup, only : setup_srtm, ssi_default
    use mo_rrtm_params, only: nbndsw, ngptsw, nbndlw, ngptlw
    use mo_lrtm_setup, only : ngb_lw => ngb
    use mo_srtm_setup, only : ngb_sw => ngb

    use radiation_config, only : config_type

    type(config_type), intent(inout), target :: config
    character(len=*), intent(in)     :: directory

    ! From the mid-latitude summer atmosphere, the original g-points
    ! are listed from least optically thick to most, as judged by the
    ! number of matrix multiplications needed in an expm computation
!!$    integer, parameter :: RRTM_GPOINT_REORDERING_LW(140) = &
!!$         (/89, 77, 90, 78, 139, 79, 80, 137, 131, 69, 97, 91, 81, 70, 71, &
!!$         53, 82, 72, 123, 54, 98, 92, 55, 83, 132, 124, 73, 56, 99, 57, &
!!$         84, 23, 125, 100, 24, 74, 85, 93, 58, 25, 86, 126, 75, 26, 11, &
!!$         101, 133, 59, 27, 87, 76, 140, 12, 102, 94, 28, 127, 13, 39, 60, &
!!$         88, 103, 109, 14, 29, 115, 40, 95, 61, 15, 41, 110, 104, 1, 42, &
!!$         116, 30, 134, 128, 138, 96, 62, 16, 43, 117, 63, 111, 44, 2, 64, &
!!$         31, 65, 105, 17, 45, 66, 118, 32, 3, 33, 67, 18, 129, 135, 46, &
!!$         112, 34, 106, 68, 35, 4, 119, 36, 47, 107, 19, 37, 38, 113, 48, &
!!$         130, 5, 120, 49, 108, 20, 50, 51, 114, 21, 121, 52, 136, 122, 6, &
!!$         22, 7, 8, 9, 10/)
!!$    integer, parameter :: RRTM_GPOINT_REORDERING_SW(112) = &
!!$         (/45, 35, 19, 27, 36, 57, 20, 46, 58, 28, 21, 67, 68, 55, 69, &
!!$         37, 1, 22, 78, 29, 59, 79, 101, 77, 70, 76, 47, 75, 30, 80, &
!!$         81, 60, 102, 82, 23, 2, 83, 84, 85, 86, 61, 103, 87, 31, 56, &
!!$         38, 88, 71, 48, 89, 3, 90, 62, 24, 7, 91, 49, 32, 104, 92, &
!!$         72, 63, 95, 39, 4, 93, 8, 96, 50, 94, 64, 40, 97, 33, 25, 51, &
!!$         73, 65, 9, 98, 41, 99, 100, 105, 52, 5, 10, 42, 66, 11, 74, &
!!$         34, 53, 26, 6, 106, 12, 43, 13, 54, 44, 107, 14, 108, 15, 16, &
!!$         109, 17, 18, 110, 111, 112/)
    integer, parameter :: RRTM_GPOINT_REORDERING_LW(140) = (/ &
          &    89, 90, 139, 77, 137, 69, 131, 97, 91, 70, 78, 71, 53, 72, 123, 54, 79, 98,  &
          &   92, 55, 80, 132, 124, 81, 73, 56, 99, 82, 57, 23, 125, 100, 24, 74, 93, 58, 25,  &
          &   83, 126, 75, 26, 11, 101, 133, 59, 27, 76, 140, 12, 84, 102, 94, 28, 127, 85,  &
          &   13, 39, 60, 86, 103, 87, 109, 14, 29, 115, 40, 95, 15, 61, 88, 41, 110, 104, 1,  &
          &   116, 42, 30, 134, 128, 138, 96, 62, 16, 43, 117, 63, 111, 44, 2, 64, 31, 65,  &
          &   105, 17, 45, 66, 118, 32, 3, 33, 67, 18, 129, 135, 46, 112, 34, 106, 68, 35, 4,  &
          &   119, 36, 47, 107, 19, 37, 38, 113, 48, 130, 5, 120, 49, 108, 20, 50, 51, 114,  &
          &   21, 121, 52, 136, 122, 6, 22, 7, 8, 9, 10 &
          & /)
    integer, parameter :: RRTM_GPOINT_REORDERING_SW(112) = (/ &
          &    35, 45, 19, 27, 36, 57, 20, 46, 58, 21, 28, 67, 55, 68, 37, 1, 69, 22, 29, 59,  &
          &   78, 101, 79, 77, 70, 76, 47, 75, 30, 81, 60, 102, 80, 82, 23, 2, 83, 84, 85,  &
          &   86, 103, 61, 31, 87, 56, 38, 71, 48, 88, 3, 62, 89, 24, 7, 49, 32, 104, 72, 90,  &
          &   63, 39, 4, 8, 50, 91, 64, 40, 33, 25, 51, 95, 96, 73, 65, 9, 41, 97, 92, 105,  &
          &   52, 5, 98, 10, 42, 99, 100, 66, 11, 74, 34, 53, 26, 6, 106, 12, 43, 13, 54, 93,  &
          &   44, 107, 94, 14, 108, 15, 16, 109, 17, 18, 110, 111, 112 &
          & /)

    call setup_lrtm(directory)
    call setup_srtm(directory)

    config%n_g_sw = ngptsw
    config%n_g_lw = ngptlw
    config%n_bands_sw = nbndsw
    config%n_bands_lw = nbndlw

    ! Wavenumber ranges of each band may be needed so that the user
    ! can compute UV and photosynthetically active radiation for a
    ! particular wavelength range
    allocate(config%wavenumber1_sw(config%n_bands_sw))
    allocate(config%wavenumber2_sw(config%n_bands_sw))
    allocate(config%wavenumber1_lw(config%n_bands_lw))
    allocate(config%wavenumber2_lw(config%n_bands_lw))
    config%wavenumber1_lw = (/ 10, 350, 500, 630, 700, 820, 980, 1080, 1180, 1390, 1480, &
         &  1800, 2080, 2250, 2380, 2600 /)
    config%wavenumber2_lw = (/ 350, 500, 630, 700, 820, 980, 1080, 1180, 1390, 1480, 1800, &
         &  2080, 2250, 2380, 2600, 3250 /)
    config%wavenumber1_sw = (/ 2600, 3250, 4000, 4650, 5150, 6150, 7700, 8050, 12850, &
         &  16000 , 22650, 29000, 38000, 820 /)
    config%wavenumber2_sw = (/ 3250, 4000, 4650, 5150, 6150, 7700, 8050, 12850, 16000, &
         &  22650, 29000, 38000, 50000, 2600 /)

    allocate(config%i_band_from_g_sw          (config%n_g_sw))
    allocate(config%i_band_from_g_lw          (config%n_g_lw))
    allocate(config%i_band_from_reordered_g_sw(config%n_g_sw))
    allocate(config%i_band_from_reordered_g_lw(config%n_g_lw))
    allocate(config%i_g_from_reordered_g_sw(config%n_g_sw))
    allocate(config%i_g_from_reordered_g_lw(config%n_g_lw))

    ! Shortwave starts at 16: need to start at 1
    config%i_band_from_g_sw = ngb_sw - ngb_sw(1)+1
    config%i_band_from_g_lw = ngb_lw

    ! Implied-do for no reordering
    !config%i_g_from_reordered_g_sw = (/ (i, i=1,ngptsw) /)
    config%i_g_from_reordered_g_sw = RRTM_GPOINT_REORDERING_SW

    !config%i_g_from_reordered_g_lw = (/ (i, i=1,ngptlw) /)
    config%i_g_from_reordered_g_lw = RRTM_GPOINT_REORDERING_LW

    config%i_band_from_reordered_g_sw &
         = config%i_band_from_g_sw(config%i_g_from_reordered_g_sw)

    config%i_band_from_reordered_g_lw &
         = config%i_band_from_g_lw(config%i_g_from_reordered_g_lw)

    ! The i_spec_* variables are used solely for storing spectral
    ! data, and this can either be by band or by g-point
    if (config%do_save_spectral_flux) then
      if (config%do_save_gpoint_flux) then
        config%n_spec_sw = config%n_g_sw
        config%n_spec_lw = config%n_g_lw
        config%i_spec_from_reordered_g_sw => config%i_g_from_reordered_g_sw
        config%i_spec_from_reordered_g_lw => config%i_g_from_reordered_g_lw
      else
        config%n_spec_sw = config%n_bands_sw
        config%n_spec_lw = config%n_bands_lw
        config%i_spec_from_reordered_g_sw => config%i_band_from_reordered_g_sw
        config%i_spec_from_reordered_g_lw => config%i_band_from_reordered_g_lw
      end if
    else
      config%n_spec_sw = 0
      config%n_spec_lw = 0
      nullify(config%i_spec_from_reordered_g_sw)
      nullify(config%i_spec_from_reordered_g_lw)
    end if

  end subroutine setup_gas_optics


  !---------------------------------------------------------------------
  ! Scale gas mixing ratios according to required units
  subroutine set_gas_units(gas)

    use radiation_gas,           only : gas_type, IVolumeMixingRatio
    type(gas_type),    intent(inout) :: gas

    call gas%set_units(IVolumeMixingRatio)

  end subroutine set_gas_units


  !---------------------------------------------------------------------
  ! Dummy setup routine for cloud optics: in fact, no setup is
  ! required for monochromatic case
  subroutine setup_cloud_optics(config)

    use radiation_config

    type(config_type), intent(inout) :: config

  end subroutine setup_cloud_optics


  !---------------------------------------------------------------------
  ! Compute gas optical depths, shortwave scattering, Planck function
  ! and incoming shortwave radiation at top-of-atmosphere
  subroutine gas_optics(ncol,nlev,istartcol,iendcol, &
       config, single_level, thermodynamics, gas, & 
       od_lw, od_sw, ssa_sw, planck_hl, planck_surf, &
       incoming_sw)

    use mo_physical_constants, only  : amd, amw, amo3, avo, grav
    use mo_rrtm_params, only         : maxxsec
    use mo_srtm_setup, only          : ssi_default
    use mo_lrtm_gas_optics, only     : gas_optics_lw
    use mo_rrtm_coeffs, only         : lrtm_coeffs, srtm_coeffs
    use mo_math_constants, only      : pi
    use mo_srtm_gas_optics, only     : gpt_taumol

    use parkind1, only               : jprb

    use radiation_io,             only : nulout
    use radiation_config,         only : config_type
    use radiation_thermodynamics, only : thermodynamics_type
    !    use radiation_single_level, only   : single_level_type
    use radiation_single_level
    use radiation_gas

    integer, intent(in) :: ncol               ! number of columns
    integer, intent(in) :: nlev               ! number of model levels
    integer, intent(in) :: istartcol, iendcol ! range of columns to process
    type(config_type), intent(in) :: config
    type(single_level_type),  intent(in) :: single_level
    type(thermodynamics_type),intent(in) :: thermodynamics
    type(gas_type),           intent(in) :: gas
    real(jprb), dimension(istartcol:iendcol,nlev,config%n_g_lw), intent(out) :: od_lw
    real(jprb), dimension(istartcol:iendcol,nlev,config%n_g_sw), intent(out) :: od_sw, ssa_sw
    real(jprb), dimension(istartcol:iendcol,nlev+1,config%n_g_lw), intent(out) :: planck_hl
    real(jprb), dimension(istartcol:iendcol,config%n_g_lw), intent(out) :: planck_surf
    real(jprb), dimension(istartcol:iendcol,config%n_g_sw), intent(out) :: incoming_sw

    real(jprb) :: incoming_sw_unscaled

    real(jprb), parameter :: cfc_scaling = 1.0e-20_jprb
    !    real(jprb), parameter :: cfc_scaling = 1.0e-40_jprb
    real(jprb), parameter :: flux_factor_lw = 1.0e4_jprb * Pi

    integer, parameter :: IH2O_sw = 1
    integer, parameter :: ICO2_sw = 2
    integer, parameter :: ICH4_sw = 3
    integer, parameter :: IO2_sw = 4
    integer, parameter :: IO3_sw = 5
    integer, parameter :: NMaxGases_sw = 5

    real(jprb), dimension(istartcol:iendcol, nlev+1, config%n_bands_lw) :: planck_hl_store
    real(jprb), dimension(istartcol:iendcol, config%n_bands_lw)         :: planck_surf_store

    real(jprb) :: flux_factor_sw

    real(jprb) :: &
                                ! Note that previously the maximum possible number of gases
                                ! here was "maxinpx" from mo_rrtm_params, but that was 38 so
                                ! too large for our needs
         gas_mol_cm2(ncol,NMaxGases,nlev), & !< number of molecules/cm2 of
         cfc_mol_cm2(ncol,maxxsec,nlev)      !< number of molecules/cm2 of

    real(jprb) :: &
         pressure_fl_hPa(ncol,nlev), &  ! Note units
         temperature_fl(ncol,nlev), &
         col_dry_vmr(ncol,nlev)  ! Column dry amount (WHAT ARE THE UNITS!)

    real(jprb) :: &
         od_abs_sw(nlev), &
         od_rayleigh(nlev)

    ! -----------------
    ! Variables for gas optics calculations
    integer ::                    &
         & laytrop (ncol     ), & !< tropopause layer index
         & jp      (ncol,nlev), & !< lookup table index 
         & jt      (ncol,nlev), & !< lookup table index 
         & jt1     (ncol,nlev), & !< lookup table index 
         & indself (ncol,nlev), &
         & indfor  (ncol,nlev), &
         & indminor(ncol,nlev) 

    real(jprb) ::                       &
         & col_gas     (ncol,nlev,NMaxGases_sw), & !< column amount (SW gases)
         & coln2o      (ncol,nlev), & !< column amount (n2o)
         & colco       (ncol,nlev), & !< column amount (co)
         & colbrd      (ncol,nlev), & !< column amount (broadening gases)
         & colmol      (ncol,nlev), & !< column amount
         & selffac     (ncol,nlev), &
         & selffrac    (ncol,nlev), &
         & forfac      (ncol,nlev), &
         & forfrac     (ncol,nlev), &           
         & minorfrac   (ncol,nlev), &
         & scaleminor  (ncol,nlev), &
         & scaleminorn2(ncol,nlev), & 
         & wbrodl      (ncol,nlev) 

    real(jprb) ::      & 
         fac00(ncol,nlev),        &
         fac01(ncol,nlev),        &
         fac10(ncol,nlev),        &
         fac11(ncol,nlev)

    real(jprb) ::                     & 
         rat_h2oco2  (ncol,nlev), &
         rat_h2oco2_1(ncol,nlev), &
         rat_h2oo3   (ncol,nlev), & 
         rat_h2oo3_1 (ncol,nlev), &
         rat_h2on2o  (ncol,nlev), &
         rat_h2on2o_1(ncol,nlev), &
         rat_h2och4  (ncol,nlev), &
         rat_h2och4_1(ncol,nlev), &
         rat_n2oco2  (ncol,nlev), & 
         rat_n2oco2_1(ncol,nlev), &
         rat_o3co2   (ncol,nlev), &
         rat_o3co2_1 (ncol,nlev)
    ! -----------------
    real(jprb) fracs(ncol,nlev) !< Planck fraction per g-point

    real(jprb) :: dpressure ! Pressure thickness of layer (Pa)
    real(jprb) :: amm       ! Molecular weight of moist air
    integer :: jlev, jcol, jgas, jg ! Loop indices for level, column, gas and g-point
    integer :: jgreorder ! Indices for reordered g-points
    integer :: icode, iband, jband

    logical :: lpressure_increasing = .true.

!!$    od_lw = 1.0e-5_jprb
!!$    od_sw = 1.0e-5_jprb
!!$    ssa_sw = 0.999999_jprb
!!$    planck_hl = 1.0_jprb
!!$    planck_surf = 1.0_jprb
!!$    incoming_sw = 1.0_jprb
!!$
!!$    return

    if (config%iverbose >= 2) then
      write(nulout,'(a)') 'Computing gas absorption/scattering properties'
    end if

    ! Check we have gas mixing ratios in the right units
    call gas%assert_units(IVolumeMixingRatio)

    ! Note that these variables are initially set to volume mixing
    ! ratio (m3/m3) and later scaled to convert to mol cm-2.
    gas_mol_cm2 = 0.0_jprb
    cfc_mol_cm2 = 0.0_jprb
    gas_mol_cm2(:,IH2O,:) = gas%mixing_ratio(:,:,IH2O)

    do jlev = 1,nlev
      do jcol = istartcol,iendcol
        pressure_fl_hPa(jcol,jlev) = 0.01*0.5_jprb &
             *(thermodynamics%pressure_hl(jcol,jlev) &
             + thermodynamics%pressure_hl(jcol,jlev+1))
        temperature_fl(jcol,jlev) &
             = 0.5_jprb*(thermodynamics%temperature_hl(jcol,jlev) &
             + thermodynamics%temperature_hl(jcol,jlev+1))
        dpressure = thermodynamics%pressure_hl(jcol,jlev+1) &
             -thermodynamics%pressure_hl(jcol,jlev)
        amm = amd*(1.0_jprb + gas_mol_cm2(jcol,IH2O,jlev)*(amw/amd - 1.0_jprb))
        col_dry_vmr(jcol,jlev) = (0.01_jprb*dpressure)*10.0_jprb*avo/grav/amm &
             / (1.0_jprb+gas_mol_cm2(jcol,IH2O,jlev))
      end do
    end do

    do jgas = 1,gas%ntype
      icode = gas%icode(jgas)
      if (icode == IH2O .or. icode == 0) then
        cycle
      end if
      if (icode == ICFC11) then
        cfc_mol_cm2(:,2,:) = gas%mixing_ratio(:,:,ICFC11)*col_dry_vmr(:,:)*cfc_scaling
      else if (icode == ICFC12) then
        cfc_mol_cm2(:,3,:) = gas%mixing_ratio(:,:,ICFC12)*col_dry_vmr(:,:)*cfc_scaling
      elseif (icode <= IH2O .or. icode > NMaxGases) then
        write(*,*) 'Gas code ', icode, ' out of range'
        stop
      else
        gas_mol_cm2(:,icode,:) = gas%mixing_ratio(:,:,icode)*col_dry_vmr(:,:)
      end if
    end do

    gas_mol_cm2(:,IH2O,:) = gas_mol_cm2(:,IH2O,:)*col_dry_vmr(:,:)

    wbrodl(1:ncol,1:nlev) = col_dry_vmr(1:ncol,1:nlev) &
         - sum(gas_mol_cm2(1:ncol,2:,1:nlev), dim=2)

    if (config%do_lw) then
      call lrtm_coeffs(    ncol         ,ncol         ,nlev         ,&
           & pressure_fl_hPa,temperature_fl,col_dry_vmr,gas_mol_cm2 ,&
           & wbrodl       ,laytrop      ,jp           ,jt           ,&
           & jt1          ,                                          &
           & col_gas(:,:,IH2O_sw)       ,col_gas(:,:,ico2_sw)       ,&
           & col_gas(:,:,IO3_sw)        ,coln2o       ,colco        ,&
           & col_gas(:,:,ICH4_sw)       ,col_gas(:,:,io2_sw)        ,&
           & colbrd       ,&
           & fac00        ,fac01        ,fac10        ,fac11        ,&
           & rat_h2oco2   ,rat_h2oco2_1 ,rat_h2oo3    ,rat_h2oo3_1  ,&
           & rat_h2on2o   ,rat_h2on2o_1 ,rat_h2och4   ,rat_h2och4_1 ,& 
           & rat_n2oco2   ,rat_n2oco2_1 ,rat_o3co2    ,rat_o3co2_1  ,&
           & selffac      ,selffrac     ,indself      ,forfac       ,&
           & forfrac      ,indfor       ,minorfrac    ,scaleminor   ,&
           & scaleminorn2 ,indminor     )

      !       write(*,*) 'pressure_fl_hPa=',pressure_fl_hPa
      !       write(*,*) 'laytrop=',laytrop
      !       do jlev = 1, nlev
      !          write(35,'(i4, 7e14.4)') jlev, gas_mol_cm2(1,1:7,jlev)
      !       end do

      ! Compute and store Planck functions for each band
      do jband = 1,config%n_bands_lw
        planck_hl_store(:,:,jband) = planckFunction(thermodynamics%temperature_hl(:,:), jband)
        planck_surf_store(:,jband) = planckFunction(single_level%skin_temperature(:), jband)
      end do

      do jgreorder = 1, config%n_g_lw
        iband = config%i_band_from_reordered_g_lw(jgreorder) 
        jg = config%i_g_from_reordered_g_lw(jgreorder)
        do jcol = 1, ncol
          !
          ! Gas concentrations in colxx variables are normalized by 1.e-20_jprb in lrtm_coeffs
          !   CFC gas concentrations (wx) need the same normalization
          !   Per Eli Mlawer the k values used in gas optics tables have been multiplied by 1e20 
          ! ALREADY DONE???
          !             wx_loc(:,:) = 1.e-20_jprb * wx(jl,:,:)
          call gas_optics_lw  (                          &
               nlev                ,jg             ,pressure_fl_hPa(jcol,:),&
               cfc_mol_cm2 (jcol,:,:),col_dry_vmr (jcol,:),laytrop     (jcol)  ,&
               jp          (jcol,:),jt          (jcol,:),jt1         (jcol,:),&
               col_gas     (jcol,:,IH2O_sw)             ,col_gas(jcol,:,ICO2_sw),&
               col_gas     (jcol,:,IO3_sw)              ,coln2o      (jcol,:),&
               colco       (jcol,:),col_gas     (jcol,:,ich4_sw)             ,&
               col_gas     (jcol,:,io2_sw)              ,colbrd      (jcol,:),&
               fac00       (jcol,:),&
               fac01       (jcol,:),fac10       (jcol,:),fac11       (jcol,:),&
               rat_h2oco2  (jcol,:),rat_h2oco2_1(jcol,:),rat_h2oo3   (jcol,:),&
               rat_h2oo3_1 (jcol,:),rat_h2on2o  (jcol,:),rat_h2on2o_1(jcol,:),&
               rat_h2och4  (jcol,:),rat_h2och4_1(jcol,:),rat_n2oco2  (jcol,:),&
               rat_n2oco2_1(jcol,:),rat_o3co2   (jcol,:),rat_o3co2_1 (jcol,:),&
               selffac     (jcol,:),selffrac    (jcol,:),indself     (jcol,:),&
               forfac      (jcol,:),forfrac     (jcol,:),indfor      (jcol,:),&
               minorfrac   (jcol,:),scaleminor  (jcol,:),scaleminorn2(jcol,:),&
               indminor    (jcol,:),fracs       (jcol,:),od_lw(jcol,:,jgreorder))
          !             write(37,'(i4,137e11.3)') jg, od_lw(jcol,nlev:1:-1,jg)
        end do  ! Loop over columns

        !          write(36,*) 'G point ', jg, ': ', fracs(1,1), fracs(1,nlev), (fracs(1,1)-fracs(1,nlev))/fracs(1,1)
        planck_hl(:,1:nlev,jgreorder) = planck_hl_store(:,1:nlev,iband)*fracs(:,:)*flux_factor_lw
        planck_hl(:,nlev+1,jgreorder) = planck_hl_store(:,nlev+1,iband)*fracs(:,nlev)*flux_factor_lw
        planck_surf(:,jgreorder) = planck_surf_store(:,iband)*fracs(:,nlev)*flux_factor_lw

      end do ! Loop over g point samples - done with gas optical depth calculations 

      if (config%iverbose >= 5) then
        write(nulout, '(a)') 'Longwave column 1: original g point, optical depth'
        do jgreorder = 1, config%n_g_lw
          write(nulout, '(a,i5,f15.5)') 'LW  ', config%i_g_from_reordered_g_lw(jgreorder), &
               &  sum(od_lw(1,:,jgreorder))
        end do
      end if

      ! Water vapor continuum broadening factors are used differently in LW and SW? 
      if (config%do_sw) then
        forfac(:,:)  =  forfac(:,:)  * col_gas(:,:,ih2o_sw)
        selffac(:,:) =  selffac(:,:) * col_gas(:,:,ih2o_sw)
      end if

    end if

    if (config%do_sw) then
      call srtm_coeffs(     ncol          ,ncol          ,nlev          , &
           & pressure_fl_hPa,temperature_fl,col_dry_vmr  ,gas_mol_cm2   , &
           & laytrop       ,jp            ,jt            ,jt1           , &
           & col_gas(:,:,ICH4_sw)         ,col_gas(:,:,ICO2_sw)         , &
           & col_gas(:,:,IH2O_sw)         ,colmol        ,coln2o        , &
           & col_gas(:,:,IO2_sw)          ,col_gas(:,:,IO3_sw)          , &
           & fac00         ,fac01         ,fac10         , &
           & fac11         ,selffac       ,selffrac      ,indself       , &
           & forfac        ,forfrac       ,indfor)
    end if

    flux_factor_sw = 1.0_jprb
    if (.not. config%use_spectral_solar_scaling) then
      flux_factor_sw = single_level%solar_irradiance/sum(ssi_default)
    end if

    if (config%do_sw) then
      !       do jgas = 1, 5
      !          write(33,'(i4,137e11.3)') jgas, col_gas(1,nlev:1:-1,jgas)
      !       end do
      !
      !  -- 2.2.2 Loop over g-points calculating gas optical properties. 
      !
      ! --------------------------------
      !IBM* ASSERT(NODEPS)
      do jgreorder = 1, config%n_g_sw
        iband = config%i_band_from_reordered_g_sw(jgreorder)
        jg = config%i_g_from_reordered_g_sw(jgreorder)
        if (config%use_spectral_solar_scaling) then
          flux_factor_sw = single_level%solar_irradiance &
               &  * single_level%spectral_solar_scaling(iband) / ssi_default(iband)
        end if

        do jcol = 1, ncol
          call gpt_taumol(nlev  ,jg       ,                              &
               & jp      (jcol,:) ,jt      (jcol,:),jt1   (jcol,:)  ,laytrop(jcol)  , &
               & indself (jcol,:) ,indfor  (jcol,:),lpressure_increasing, &
               & col_gas(jcol,:,:),colmol(jcol,:), &
               & fac00   (jcol,:) ,fac01   (jcol,:),fac10 (jcol,:)  ,fac11  (jcol,:), &
               & selffac (jcol,:) ,selffrac(jcol,:),forfac(jcol,:)  ,forfrac(jcol,:), &
               & incoming_sw_unscaled,od_abs_sw       ,od_rayleigh)
          !             write(34,'(i4,137e11.3)') jg, od_abs_sw(nlev:1:-1)
          od_sw(jcol,:,jgreorder)  = od_abs_sw + od_rayleigh
          ssa_sw(jcol,:,jgreorder) = od_rayleigh / (od_abs_sw + od_rayleigh)
          incoming_sw(jcol,jgreorder) = incoming_sw_unscaled * flux_factor_sw 
          !* single_level%cos_sza(jcol)
          !              write(40,'(i4,137e11.3)') jg, od_abs_sw(nlev:1:-1)
          !              write(41,'(i4,137e11.3)') jg, od_rayleigh(nlev:1:-1)
        end do
      end do

      if (config%iverbose >= 5) then
        write(nulout, '(a)') 'Shortwave column 1: original g point, optical depth'
        do jgreorder = 1, config%n_g_sw
          write(nulout, '(a,i5,f15.5)') 'SW  ', config%i_g_from_reordered_g_sw(jgreorder), &
               &  sum(od_sw(1,:,jgreorder))
        end do
      end if

    end if

  end subroutine gas_optics


  ! Compute cloud optical depth, single-scattering albedo and
  ! asymmetry factor in the longwave and shortwave
  subroutine cloud_optics(ncol,nlev,istartcol,iendcol, &
       config, single_level, thermodynamics, cloud, & 
       od_lw_cloud, ssa_lw_cloud, g_lw_cloud, &
       od_sw_cloud, ssa_sw_cloud, g_sw_cloud)

    use mo_cloud_optics, only          : rrtm_cloud_optics => cloud_optics
    use mo_physical_constants, only    : grav

    use parkind1, only           : jprb
    use radiation_config, only  : config_type
    use radiation_thermodynamics, only : thermodynamics_type
    use radiation_single_level, only   : single_level_type
    use radiation_cloud, only    : cloud_type

    integer, intent(in) :: ncol               ! number of columns
    integer, intent(in) :: nlev               ! number of model levels
    integer, intent(in) :: istartcol, iendcol ! range of columns to process
    type(config_type), intent(in) :: config
    type(single_level_type),  intent(in) :: single_level
    type(thermodynamics_type),intent(in) :: thermodynamics
    type(cloud_type),   intent(in) :: cloud
    real(jprb), dimension(istartcol:iendcol,nlev,config%n_bands_lw), intent(out) :: &
         od_lw_cloud
    real(jprb), dimension(istartcol:iendcol,nlev,config%n_bands_lw_if_scattering), &
         intent(out) :: ssa_lw_cloud, g_lw_cloud
    real(jprb), dimension(istartcol:iendcol,nlev,config%n_bands_sw), intent(out) :: &
         od_sw_cloud, ssa_sw_cloud, g_sw_cloud

    ! In-cloud liquid and ice water path in a layer, in g m-2
    real(jprb), dimension(istartcol:iendcol,nlev) :: LWP_g_m2, IWP_g_m2

    ! Effective radius in microns
    real(jprb), dimension(istartcol:iendcol,nlev) :: re_liq_um, re_ice_um

    integer  :: jlev, jcol
    real(jprb) :: factor

    ! Convert cloud mixing ratio into liquid and ice water path
    ! in each layer
    do jlev = 1, nlev
      do jcol = 1, iendcol
        ! Factor to convert from gridbox-mean mass mixing ratio to
        ! in-cloud water path involves 1000 (to convert from kg to
        ! g), the pressure difference in Pa, acceleration due to
        ! gravity and cloud fraction adjusted to avoid division by
        ! zero.
        factor = 1.0e3_jprb &
             * ( thermodynamics%pressure_hl(jcol,jlev+1)    &
             -thermodynamics%pressure_hl(jcol,jlev  )  ) &
             / (grav &
             * max(epsilon(1.0_jprb), cloud%fraction(jcol,jlev)))
        LWP_g_m2(jcol,jlev) = factor * cloud%q_liq(jcol,jlev)
        IWP_g_m2(jcol,jlev) = factor * cloud%q_ice   (jcol,jlev)
      end do
    end do

    ! Array-wise operations
    re_liq_um = 1.0e6_jprb * cloud%re_liq
    re_ice_um = 1.0e6_jprb * cloud%re_ice

    call rrtm_cloud_optics(iendcol, ncol, nlev, LWP_g_m2, IWP_g_m2, &
         re_liq_um, re_ice_um, &
         od_sw_cloud, ssa_sw_cloud, g_sw_cloud, &
         od_lw_cloud)

    if (config%do_lw_cloud_scattering) then
      ! RRTM cloud optics does not provide single-scattering albedo
      ! and asymmetry factor in the longwave, so we need to set these
      ! values to zero; the following are array-wise operations.
      ssa_lw_cloud = 0.0_jprb
      g_lw_cloud = 0.0_jprb
    end if

  end subroutine cloud_optics


  ELEMENTAL FUNCTION planckFunction(temp, band)
    !
    ! Compute the blackbody emission in a given band as a function of
    ! temperature
    !

    use parkind1, only : jprb
    use rrlw_planck, only    : totplanck
    use mo_lrtm_setup, only  : delwave

    REAL(jprb), INTENT(IN) :: temp
    INTEGER,  INTENT(IN)   :: band 
    REAL(jprb)             :: planckFunction

    INTEGER    :: ind
    REAL(jprb) :: frac

    ind = MIN(MAX(1, INT(temp - 159._jprb)),180)
    frac = temp - 159._jprb - float(ind)

    planckFunction = totplanck(ind, band) &
         + frac * (totplanck(ind+1, band) - totplanck(ind, band))
    planckFunction = planckFunction * delwave(band)
  END FUNCTION planckFunction

end module radiation_psrad_rrtm
