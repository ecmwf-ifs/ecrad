! ifs_blocking.F90 - Reshuffle ecRad data into an NPROMA-blocked data structure
!
! (C) Copyright 2022- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
!
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!
! Author:  Balthasar Reuter
! Email:   balthasar.reuter@ecmwf.int
!

module ifs_blocking

  use parkind1,                 only : jprb, jpim ! Working precision, integer type

  implicit none

  public

  type :: ifs_config_type
    ! Offsets in ZRGP
    integer :: igi, imu0, iamu0, iemiss, its, islm, iccnl,    &
        &     ibas, itop, igelam, igemu, iclon, islon, iald, ialp, iti, ipr, iqs, iwv, iclc, ilwa,    &
        &     iiwa, iswa, irwa, irra, idp, ioz, iecpo3, ihpr, iaprs, ihti, iaero, ifrsod, icdir,      &
        &     ifrted, ifrsodc, ifrtedc, iemit, isudu, iuvdf, iparf, iparcf, itincf, ifdir, ifdif,     &
        &     ilwderivative, iswdirectband, iswdiffuseband, ifrso, iswfc, ifrth, ilwfc, iaer,         &
        &     iich4, iin2o, ino2, ic11, ic12, igix, iico2, iccno, ic22, icl4
    integer :: ire_liq, ire_ice, ioverlap
    integer :: ifldstot
  end type ifs_config_type

contains

integer(kind=jpim) function indrad(knext,kflds,lduse)

  integer(kind=jpim), intent(inout) :: knext
  integer(kind=jpim), intent(in) :: kflds
  logical, intent(in) :: lduse

  if( lduse ) then
    indrad=knext
    knext=knext+kflds
  else
    indrad=-99999999
  endif

end function indrad

subroutine ifs_setup_indices (driver_config, ifs_config, yradiation, nlev)

  use radiation_io,             only : nulout
  use ecrad_driver_config,      only : driver_config_type
  use radiation_setup,          only : tradiation

  ! Configuration specific to this driver
  type(driver_config_type), intent(in)     :: driver_config
  type(ifs_config_type), intent(inout)     :: ifs_config

  ! Configuration for the radiation scheme, IFS style
  type(tradiation), intent(inout)          :: yradiation

  integer, intent(inout) :: nlev

  integer :: ifldsin, ifldsout, inext, iinbeg, iinend, ioutbeg, ioutend
  logical :: llactaero
  logical :: lldebug

  ! Extract some config values
  lldebug=(driver_config%iverbose>4)     ! debug
  llactaero = .false.
  if(yradiation%rad_config%n_aerosol_types > 0 .and.&
    & yradiation%rad_config%n_aerosol_types <= 21 .and. yradiation%yrerad%naermacc == 0) then
    llactaero = .true.
  endif

  !
  ! RADINTG
  !

  !  INITIALISE INDICES FOR VARIABLE

  ! INDRAD is a CONTAIN'd function (now a module function)

  inext  =1
  iinbeg =1                        ! start of input variables
  ifs_config%igi    =indrad(inext,1,lldebug)
  ifs_config%imu0   =indrad(inext,1,.true.)
  ifs_config%iamu0  =indrad(inext,1,.true.)
  ifs_config%iemiss =indrad(inext,yradiation%yrerad%nlwemiss,.true.)
  ifs_config%its    =indrad(inext,1,.true.)
  ifs_config%islm   =indrad(inext,1,.true.)
  ifs_config%iccnl  =indrad(inext,1,.true.)
  ifs_config%iccno  =indrad(inext,1,.true.)
  ifs_config%ibas   =indrad(inext,1,.true.)
  ifs_config%itop   =indrad(inext,1,.true.)
  ifs_config%igelam =indrad(inext,1,.true.)
  ifs_config%igemu  =indrad(inext,1,.true.)
  ifs_config%iclon  =indrad(inext,1,.true.)
  ifs_config%islon  =indrad(inext,1,.true.)
  ifs_config%iald   =indrad(inext,yradiation%yrerad%nsw,.true.)
  ifs_config%ialp   =indrad(inext,yradiation%yrerad%nsw,.true.)
  ifs_config%iti    =indrad(inext,nlev,.true.)
  ifs_config%ipr    =indrad(inext,nlev,.true.)
  ifs_config%iqs    =indrad(inext,nlev,.true.)
  ifs_config%iwv    =indrad(inext,nlev,.true.)
  ifs_config%iclc   =indrad(inext,nlev,.true.)
  ifs_config%ilwa   =indrad(inext,nlev,.true.)
  ifs_config%iiwa   =indrad(inext,nlev,.true.)
  ifs_config%iswa   =indrad(inext,nlev,.true.)
  ifs_config%irwa   =indrad(inext,nlev,.true.)
  ifs_config%irra   =indrad(inext,nlev,.true.)
  ifs_config%idp    =indrad(inext,nlev,.true.)
  ifs_config%ioz    =indrad(inext,nlev,.false.)
  ifs_config%iecpo3 =indrad(inext,nlev ,.false.)
  ifs_config%ihpr   =indrad(inext,nlev+1,.true.) ! not used in ecrad
  ifs_config%iaprs  =indrad(inext,nlev+1,.true.)
  ifs_config%ihti   =indrad(inext,nlev+1,.true.)
  ifs_config%iaero  =indrad(inext,yradiation%rad_config%n_aerosol_types*nlev,&
                          & llactaero .and. yradiation%yrerad%naermacc==0)

  iinend =inext-1                  ! end of input variables

  ioutbeg=inext                    ! start of output variables
  if (yradiation%yrerad%naermacc == 1) then
    ifs_config%iaero = indrad(inext,yradiation%rad_config%n_aerosol_types*nlev,&
                            & yradiation%yrerad%ldiagforcing)
  endif
  ifs_config%ifrsod =indrad(inext,1,.true.)
  ifs_config%ifrted =indrad(inext,yradiation%yrerad%nlwout,.true.)
  ifs_config%ifrsodc=indrad(inext,1,.true.)
  ifs_config%ifrtedc=indrad(inext,1,.true.)
  ifs_config%iemit  =indrad(inext,1,.true.)
  ifs_config%isudu  =indrad(inext,1,.true.)
  ifs_config%iuvdf  =indrad(inext,1,.true.)
  ifs_config%iparf  =indrad(inext,1,.true.)
  ifs_config%iparcf =indrad(inext,1,.true.)
  ifs_config%itincf =indrad(inext,1,.true.)
  ifs_config%ifdir  =indrad(inext,1,.true.)
  ifs_config%ifdif  =indrad(inext,1,.true.)
  ifs_config%icdir  =indrad(inext,1,.true.)
  ifs_config%ilwderivative =indrad(inext,nlev+1, yradiation%yrerad%lapproxlwupdate)
  ifs_config%iswdirectband =indrad(inext,yradiation%yrerad%nsw,yradiation%yrerad%lapproxswupdate)
  ifs_config%iswdiffuseband=indrad(inext,yradiation%yrerad%nsw,yradiation%yrerad%lapproxswupdate)
  ifs_config%ifrso  =indrad(inext,nlev+1,.true.)
  ifs_config%iswfc  =indrad(inext,nlev+1,.true.)
  ifs_config%ifrth  =indrad(inext,nlev+1,.true.)
  ifs_config%ilwfc  =indrad(inext,nlev+1,.true.)
  ifs_config%iaer   =indrad(inext,6*nlev,yradiation%yrerad%ldiagforcing)
  ifs_config%ioz    =indrad(inext,nlev,yradiation%yrerad%ldiagforcing)
  ifs_config%iico2  =indrad(inext,nlev,yradiation%yrerad%ldiagforcing)
  ifs_config%iich4  =indrad(inext,nlev,yradiation%yrerad%ldiagforcing)
  ifs_config%iin2o  =indrad(inext,nlev,yradiation%yrerad%ldiagforcing)
  ifs_config%ino2   =indrad(inext,nlev,yradiation%yrerad%ldiagforcing)
  ifs_config%ic11   =indrad(inext,nlev,yradiation%yrerad%ldiagforcing)
  ifs_config%ic12   =indrad(inext,nlev,yradiation%yrerad%ldiagforcing)
  ifs_config%ic22   =indrad(inext,nlev,yradiation%yrerad%ldiagforcing)
  ifs_config%icl4   =indrad(inext,nlev,yradiation%yrerad%ldiagforcing)
  ifs_config%igix   =indrad(inext,1,lldebug)

  ioutend=inext-1                  ! end of output variables

                                ! start of local variables
  if(.not.yradiation%yrerad%ldiagforcing) then
    if (yradiation%rad_config%n_aerosol_types == 0 .or. yradiation%yrerad%naermacc == 1) then
      ifs_config%iaero = indrad(inext,yradiation%rad_config%n_aerosol_types*nlev,.true.)
    endif
    ifs_config%iaer   =indrad(inext,nlev*6,.true.)
    ifs_config%ioz    =indrad(inext,nlev,.true.)
    ifs_config%iico2  =indrad(inext,nlev,.true.)
    ifs_config%iich4  =indrad(inext,nlev,.true.)
    ifs_config%iin2o  =indrad(inext,nlev,.true.)
    ifs_config%ino2   =indrad(inext,nlev,.true.)
    ifs_config%ic11   =indrad(inext,nlev,.true.)
    ifs_config%ic12   =indrad(inext,nlev,.true.)
    ifs_config%ic22   =indrad(inext,nlev,.true.)
    ifs_config%icl4   =indrad(inext,nlev,.true.)
  endif
                                ! end of local variables

                                  ! start of standalone inputs workaround variables
  ifs_config%ire_liq =indrad(inext,nlev,.true.)
  ifs_config%ire_ice =indrad(inext,nlev,.true.)
  ifs_config%ioverlap =indrad(inext,nlev-1,.true.)
                                  ! end of standalone inputs workaround variables

  ifldsin = iinend - iinbeg +1
  ifldsout= ioutend-ioutbeg +1
  ifs_config%ifldstot= inext  - 1

  if( lldebug )then
    write(nulout,'("imu0   =",i0)')ifs_config%imu0
    write(nulout,'("iamu0  =",i0)')ifs_config%iamu0
    write(nulout,'("iemiss =",i0)')ifs_config%iemiss
    write(nulout,'("its    =",i0)')ifs_config%its
    write(nulout,'("islm   =",i0)')ifs_config%islm
    write(nulout,'("iccnl  =",i0)')ifs_config%iccnl
    write(nulout,'("iccno  =",i0)')ifs_config%iccno
    write(nulout,'("ibas   =",i0)')ifs_config%ibas
    write(nulout,'("itop   =",i0)')ifs_config%itop
    write(nulout,'("igelam =",i0)')ifs_config%igelam
    write(nulout,'("igemu  =",i0)')ifs_config%igemu
    write(nulout,'("iclon  =",i0)')ifs_config%iclon
    write(nulout,'("islon  =",i0)')ifs_config%islon
    write(nulout,'("iald   =",i0)')ifs_config%iald
    write(nulout,'("ialp   =",i0)')ifs_config%ialp
    write(nulout,'("iti    =",i0)')ifs_config%iti
    write(nulout,'("ipr    =",i0)')ifs_config%ipr
    write(nulout,'("iqs    =",i0)')ifs_config%iqs
    write(nulout,'("iwv    =",i0)')ifs_config%iwv
    write(nulout,'("iclc   =",i0)')ifs_config%iclc
    write(nulout,'("ilwa   =",i0)')ifs_config%ilwa
    write(nulout,'("iiwa   =",i0)')ifs_config%iiwa
    write(nulout,'("iswa   =",i0)')ifs_config%iswa
    write(nulout,'("irwa   =",i0)')ifs_config%irwa
    write(nulout,'("irra   =",i0)')ifs_config%irra
    write(nulout,'("idp    =",i0)')ifs_config%idp
    write(nulout,'("ioz    =",i0)')ifs_config%ioz
    write(nulout,'("iecpo3 =",i0)')ifs_config%iecpo3
    write(nulout,'("ihpr   =",i0)')ifs_config%ihpr
    write(nulout,'("iaprs  =",i0)')ifs_config%iaprs
    write(nulout,'("ihti   =",i0)')ifs_config%ihti
    write(nulout,'("ifrsod =",i0)')ifs_config%ifrsod
    write(nulout,'("ifrted =",i0)')ifs_config%ifrted
    write(nulout,'("ifrsodc=",i0)')ifs_config%ifrsodc
    write(nulout,'("ifrtedc=",i0)')ifs_config%ifrtedc
    write(nulout,'("iemit  =",i0)')ifs_config%iemit
    write(nulout,'("isudu  =",i0)')ifs_config%isudu
    write(nulout,'("iuvdf  =",i0)')ifs_config%iuvdf
    write(nulout,'("iparf  =",i0)')ifs_config%iparf
    write(nulout,'("iparcf =",i0)')ifs_config%iparcf
    write(nulout,'("itincf =",i0)')ifs_config%itincf
    write(nulout,'("ifdir  =",i0)')ifs_config%ifdir
    write(nulout,'("ifdif  =",i0)')ifs_config%ifdif
    write(nulout,'("icdir  =",i0)')ifs_config%icdir
    write(nulout,'("ilwderivative  =",i0)')ifs_config%ilwderivative
    write(nulout,'("iswdirectband  =",i0)')ifs_config%iswdirectband
    write(nulout,'("iswdiffuseband =",i0)')ifs_config%iswdiffuseband
    write(nulout,'("ifrso  =",i0)')ifs_config%ifrso
    write(nulout,'("iswfc  =",i0)')ifs_config%iswfc
    write(nulout,'("ifrth  =",i0)')ifs_config%ifrth
    write(nulout,'("ilwfc  =",i0)')ifs_config%ilwfc
    write(nulout,'("igi    =",i0)')ifs_config%igi
    write(nulout,'("iaer   =",i0)')ifs_config%iaer
    write(nulout,'("iaero  =",i0)')ifs_config%iaero
    write(nulout,'("iico2  =",i0)')ifs_config%iico2
    write(nulout,'("iich4  =",i0)')ifs_config%iich4
    write(nulout,'("iin2o  =",i0)')ifs_config%iin2o
    write(nulout,'("ino2   =",i0)')ifs_config%ino2
    write(nulout,'("ic11   =",i0)')ifs_config%ic11
    write(nulout,'("ic12   =",i0)')ifs_config%ic12
    write(nulout,'("ic22   =",i0)')ifs_config%ic22
    write(nulout,'("icl4   =",i0)')ifs_config%icl4
    write(nulout,'("ire_liq=",i0)')ifs_config%ire_liq
    write(nulout,'("ire_ice=",i0)')ifs_config%ire_ice
    write(nulout,'("ioverlap=",i0)')ifs_config%ioverlap
    write(nulout,'("ifldsin =",i0)')ifldsin
    write(nulout,'("ifldsout=",i0)')ifldsout
    write(nulout,'("ifldstot=",i0)')ifs_config%ifldstot
  endif

end subroutine ifs_setup_indices

subroutine ifs_copy_inputs_to_blocked ( &
  & driver_config, ifs_config, yradiation, ncol, nlev, &
  & single_level, thermodynamics, gas, cloud, aerosol, &
  & sin_latitude, longitude_rad, land_frac, pressure_fl, temperature_fl, &
  & zrgp, thermodynamics_out, iseed)

  use radiation_single_level,   only : single_level_type
  use radiation_thermodynamics, only : thermodynamics_type
  use radiation_gas,            only : gas_type, IMassMixingRatio, &
        &   IH2O, ICO2, IO3, IN2O, ICH4, ICFC11, ICFC12, IHCFC22, ICCL4
  use radiation_cloud,          only : cloud_type
  use radiation_aerosol,        only : aerosol_type
  use ecrad_driver_config,      only : driver_config_type
  use radiation_setup,          only : tradiation

  implicit none

  ! Configuration specific to this driver
  type(driver_config_type), intent(in)     :: driver_config

  type(ifs_config_type), intent(in)     :: ifs_config

  ! Configuration for the radiation scheme, IFS style
  type(tradiation), intent(in)          :: yradiation

  integer, intent(in) :: ncol, nlev         ! Number of columns and levels

  ! Derived types for the inputs to the radiation scheme
  type(single_level_type), intent(in)   :: single_level
  type(thermodynamics_type), intent(in) :: thermodynamics
  type(gas_type), intent(in)            :: gas
  type(cloud_type), intent(in)          :: cloud
  type(aerosol_type), intent(in)        :: aerosol

  ! Additional input data, required for effective radii calculation
  real(jprb), dimension(:), intent(in)   :: sin_latitude, longitude_rad, land_frac
  real(jprb), dimension(:,:), intent(in) :: pressure_fl, temperature_fl

  ! monolithic IFS data structure to pass to radiation scheme
  real(kind=jprb), intent(out), allocatable :: zrgp(:,:,:)

  ! Empty thermodynamics type to store pressure_hl for output at the end
  type(thermodynamics_type), intent(inout), optional  :: thermodynamics_out

  ! Seed for random number generator
  integer, intent(out), allocatable, optional :: iseed(:,:)

  ! number of column blocks, block size
  integer :: ngpblks, nproma

  integer :: jrl, ibeg, iend, il, ib, ifld, jemiss, jalb, jlev, joff, jaer

  ! Extract some config values
  nproma=driver_config%nblocksize        ! nproma size
  ngpblks=(ncol-1)/nproma+1              ! number of column blocks

  ! Allocate blocked data structure
  allocate(zrgp(nproma,ifs_config%ifldstot,ngpblks))
  if(present(thermodynamics_out)) allocate(thermodynamics_out%pressure_hl(ncol,nlev+1))
  if(present(iseed)) allocate(iseed(nproma,ngpblks))

  ! First touch
  !$OMP PARALLEL DO SCHEDULE(RUNTIME)&
  !$OMP&PRIVATE(IB,IFLD)
  do ib=1,ngpblks
    do ifld=1,ifs_config%ifldstot
      zrgp(:,ifld,ib) = 0._jprb
    enddo
    if(present(iseed)) iseed(:,ib) = 0
  enddo
  !$OMP END PARALLEL DO

  associate(yderad=>yradiation%yrerad, rad_config=>yradiation%rad_config)

    ! REPLACED ich4 with iich4 due to clash
    ! REPLACED in2o with iin2o due to clash
    ! REPLACED ico2 with iico2 due to clash

    !  -------------------------------------------------------
    !
    !  INPUT LOOP
    !
    !  -------------------------------------------------------

    !$OMP PARALLEL DO SCHEDULE(RUNTIME)&
    !$OMP&PRIVATE(JRL,IBEG,IEND,IL,IB,JAER,JOFF,JLEV,JALB)
    do jrl=1,ncol,nproma

      ibeg=jrl
      iend=min(ibeg+nproma-1,ncol)
      il=iend-ibeg+1
      ib=(jrl-1)/nproma+1

      !* RADINTG:  3.      PREPARE INPUT ARRAYS

      ! zrgp(1:il,imu0,ib)  = ???
      zrgp(1:il,ifs_config%iamu0,ib)  =  single_level%cos_sza(ibeg:iend)   ! cosine of solar zenith ang (mu0)

      do jemiss=1,yderad%nlwemiss
        zrgp(1:il,ifs_config%iemiss+jemiss-1,ib)  =  single_level%lw_emissivity(ibeg:iend,jemiss)
      enddo

      zrgp(1:il,ifs_config%its,ib)      = single_level%skin_temperature(ibeg:iend)  ! skin temperature
      zrgp(1:il,ifs_config%islm,ib)     = land_frac(ibeg:iend) ! land-sea mask
      zrgp(1:il,ifs_config%iccnl,ib)    = yderad%rccnlnd ! CCN over land
      zrgp(1:il,ifs_config%iccno,ib)    = yderad%rccnsea ! CCN over sea
      ! zrgp(1:il,ibas,ib)     = ???
      ! zrgp(1:il,itop,ib)     = ???
      zrgp(1:il,ifs_config%igelam,ib)   = longitude_rad(ibeg:iend) ! longitude
      zrgp(1:il,ifs_config%igemu,ib)    = sin_latitude(ibeg:iend) ! sine of latitude
      ! zrgp(1:il,iclon,ib)    = ???
      ! zrgp(1:il,islon,ib)    = ???

      do jalb=1,yderad%nsw
        zrgp(1:il,ifs_config%iald+jalb-1,ib)  =  single_level%sw_albedo(ibeg:iend,jalb)
      enddo

      if (allocated(single_level%sw_albedo_direct)) then
        do jalb=1,yderad%nsw
          zrgp(1:il,ifs_config%ialp+jalb-1,ib)  =  single_level%sw_albedo_direct(ibeg:iend,jalb)
        end do
      else
        do jalb=1,yderad%nsw
          zrgp(1:il,ifs_config%ialp+jalb-1,ib)  =  single_level%sw_albedo(ibeg:iend,jalb)
        end do
      end if
      
      do jlev=1,nlev
        zrgp(1:il,ifs_config%iti+jlev-1,ib)   = temperature_fl(ibeg:iend,jlev) ! full level temperature
        zrgp(1:il,ifs_config%ipr+jlev-1,ib)   = pressure_fl(ibeg:iend,jlev) ! full level pressure
        ! zrgp(1:il,iqs+jlev-1,ib)   = ???
      enddo

      do jlev=1,nlev
        zrgp(1:il,ifs_config%iwv+jlev-1,ib)   = gas%mixing_ratio(ibeg:iend,jlev,IH2O) ! this is already in MassMixingRatio units
        if (rad_config%do_clouds) then
          zrgp(1:il,ifs_config%iclc+jlev-1,ib)  = cloud%fraction(ibeg:iend,jlev)
          zrgp(1:il,ifs_config%ilwa+jlev-1,ib)  = cloud%q_liq(ibeg:iend,jlev)
          zrgp(1:il,ifs_config%iiwa+jlev-1,ib)  = cloud%q_ice(ibeg:iend,jlev)
        else
          zrgp(1:il,ifs_config%iclc+jlev-1,ib)  = 0._jprb
          zrgp(1:il,ifs_config%ilwa+jlev-1,ib)  = 0._jprb
          zrgp(1:il,ifs_config%iiwa+jlev-1,ib)  = 0._jprb
        endif
        zrgp(1:il,ifs_config%iswa+jlev-1,ib)  = 0._jprb  ! snow
        zrgp(1:il,ifs_config%irwa+jlev-1,ib)  = 0._jprb  ! rain

        ! zrgp(1:il,irra+jlev-1,ib)  = ???
        ! zrgp(1:il,idp+jlev-1,ib)   = ???
        ! zrgp(1:il,ifsd+jlev-1,ib)   = ???
        ! zrgp(1:il,iecpo3+jlev-1,ib) = ???
      enddo

      zrgp(1:il,ifs_config%iaer:ifs_config%iaer+nlev,ib)  =  0._jprb ! old aerosol, not used
      if (yderad%naermacc == 1) then
        joff=ifs_config%iaero
        do jaer=1,rad_config%n_aerosol_types
          do jlev=1,nlev
            zrgp(1:il,joff,ib) = aerosol%mixing_ratio(ibeg:iend,jlev,jaer)
            joff=joff+1
          enddo
        enddo
      endif

      do jlev=1,nlev+1
        ! zrgp(1:il,ihpr+jlev-1,ib)  = ???
        zrgp(1:il,ifs_config%iaprs+jlev-1,ib) = thermodynamics%pressure_hl(ibeg:iend,jlev)
        zrgp(1:il,ifs_config%ihti+jlev-1,ib)  = thermodynamics%temperature_hl(ibeg:iend,jlev)
      enddo

      ! -- by default, globally averaged concentrations (mmr)
      call gas%get(ICO2, IMassMixingRatio, zrgp(1:il,ifs_config%iico2:ifs_config%iico2+nlev-1,ib), istartcol=ibeg)
      call gas%get(ICH4, IMassMixingRatio, zrgp(1:il,ifs_config%iich4:ifs_config%iich4+nlev-1,ib), istartcol=ibeg)
      call gas%get(IN2O, IMassMixingRatio, zrgp(1:il,ifs_config%iin2o:ifs_config%iin2o+nlev-1,ib), istartcol=ibeg)
      call gas%get(ICFC11, IMassMixingRatio, zrgp(1:il,ifs_config%ic11:ifs_config%ic11+nlev-1,ib), istartcol=ibeg)
      call gas%get(ICFC12, IMassMixingRatio, zrgp(1:il,ifs_config%ic12:ifs_config%ic12+nlev-1,ib), istartcol=ibeg)
      call gas%get(IHCFC22,IMassMixingRatio, zrgp(1:il,ifs_config%ic22:ifs_config%ic22+nlev-1,ib), istartcol=ibeg)
      call gas%get(ICCL4,  IMassMixingRatio, zrgp(1:il,ifs_config%icl4:ifs_config%icl4+nlev-1,ib), istartcol=ibeg)
      call gas%get(IO3, IMassMixingRatio, zrgp(1:il,ifs_config%ioz:ifs_config%ioz+nlev-1,ib), istartcol=ibeg)
      ! convert ozone kg/kg to Pa*kg/kg
      ! do jlev=1,nlev
      !   zrgp(1:il,ifs_config%ioz+jlev-1,ib)  = zrgp(1:il,ifs_config%ioz+jlev-1,ib) &
      !         &                       * (thermodynamics%pressure_hl(ibeg:iend,jlev+1) &
      !         &                         - thermodynamics%pressure_hl(ibeg:iend,jlev))
      ! enddo

      ! local workaround variables for standalone input files
      if (rad_config%do_clouds) then
        do jlev=1,nlev
          ! missing full-level temperature and pressure as well as land-sea-mask
          zrgp(1:il,ifs_config%ire_liq+jlev-1,ib) = cloud%re_liq(ibeg:iend,jlev)
          zrgp(1:il,ifs_config%ire_ice+jlev-1,ib) = cloud%re_ice(ibeg:iend,jlev)
        enddo
        do jlev=1,nlev-1
          ! for the love of it, I can't figure this one out. Probably to do with
          ! my crude approach of setting PGEMU?
          zrgp(1:il,ifs_config%ioverlap+jlev-1,ib) = cloud%overlap_param(ibeg:iend,jlev)
        enddo
        if(present(iseed)) iseed(1:il,ib) = single_level%iseed(ibeg:iend)
      else
        do jlev=1,nlev
          ! missing full-level temperature and pressure as well as land-sea-mask
          zrgp(1:il,ifs_config%ire_liq+jlev-1,ib) = 0._jprb
          zrgp(1:il,ifs_config%ire_ice+jlev-1,ib) = 0._jprb
        enddo
        do jlev=1,nlev-1
          zrgp(1:il,ifs_config%ioverlap+jlev-1,ib) = 0._jprb
        enddo
        if(present(iseed)) iseed(1:il,ib) = 0
      endif ! do_clouds
    enddo
    !$OMP END PARALLEL DO

    ! Store pressure for output
    if(present(thermodynamics_out)) thermodynamics_out%pressure_hl(:,:) = thermodynamics%pressure_hl(:,:)

  end associate

end subroutine ifs_copy_inputs_to_blocked

subroutine ifs_copy_fluxes_from_blocked(&
    & driver_config, ifs_config, yradiation, ncol, nlev,&
    & zrgp, flux, flux_sw_direct_normal, flux_uv, flux_par, flux_par_clear,&
    & emissivity_out, flux_diffuse_band, flux_direct_band)
  use ecrad_driver_config,      only : driver_config_type
  use radiation_setup,          only : tradiation
  use radiation_flux,           only : flux_type

  ! Configuration specific to this driver
  type(driver_config_type), intent(in)     :: driver_config

  type(ifs_config_type), intent(in)     :: ifs_config

  ! Configuration for the radiation scheme, IFS style
  type(tradiation), intent(in)          :: yradiation

  integer, intent(in) :: ncol, nlev         ! Number of columns and levels

  ! monolithic IFS data structure passed to radiation scheme
  real(kind=jprb), intent(inout), allocatable :: zrgp(:,:,:)

  ! Derived type containing outputs from the radiation scheme
  type(flux_type), intent(inout)              :: flux

  ! Additional output fluxes as arrays
  real(jprb), dimension(:), intent(inout)     :: flux_sw_direct_normal, flux_uv, flux_par,&
                                                 & flux_par_clear, emissivity_out
  real(jprb), dimension(:,:), intent(inout) :: flux_diffuse_band, flux_direct_band

  ! number of column blocks, block size
  integer :: ngpblks, nproma

  integer :: jrl, ibeg, iend, il, ib, jlev, jg

  ! Extract some config values
  nproma=driver_config%nblocksize        ! nproma size
  ngpblks=(ncol-1)/nproma+1              ! number of column blocks

    !  -------------------------------------------------------
    !
    !  OUTPUT LOOP
    !
    !  -------------------------------------------------------

    !$OMP PARALLEL DO SCHEDULE(RUNTIME)&
    !$OMP&PRIVATE(JRL,IBEG,IEND,IL,IB,JLEV,JG)
    do jrl=1,ncol,nproma
      ibeg=jrl
      iend=min(ibeg+nproma-1,ncol)
      il=iend-ibeg+1
      ib=(jrl-1)/nproma+1

      do jlev=1,nlev+1
        flux%sw_up(ibeg:iend,jlev) = zrgp(1:il,ifs_config%ifrso+jlev-1,ib)
        flux%lw_up(ibeg:iend,jlev) = zrgp(1:il,ifs_config%ifrth+jlev-1,ib)
        flux%sw_up_clear(ibeg:iend,jlev) = zrgp(1:il,ifs_config%iswfc+jlev-1,ib)
        flux%lw_up_clear(ibeg:iend,jlev) = zrgp(1:il,ifs_config%ilwfc+jlev-1,ib)
        if (yradiation%yrerad%lapproxlwupdate) then
          flux%lw_derivatives(ibeg:iend,jlev) = zrgp(1:il,ifs_config%ilwderivative+jlev-1,ib)
        else
          flux%lw_derivatives(ibeg:iend,jlev) = 0.0_jprb
        endif
      end do
      flux%sw_dn(ibeg:iend,nlev+1) = zrgp(1:il,ifs_config%ifrsod,ib)
      flux%lw_dn(ibeg:iend,nlev+1) = zrgp(1:il,ifs_config%ifrted,ib)
      flux%sw_dn_clear(ibeg:iend,nlev+1) = zrgp(1:il,ifs_config%ifrsodc,ib)
      flux%lw_dn_clear(ibeg:iend,nlev+1) = zrgp(1:il,ifs_config%ifrtedc,ib)
      flux%sw_dn_direct(ibeg:iend,nlev+1) = zrgp(1:il,ifs_config%ifdir,ib)
      flux%sw_dn_direct_clear(ibeg:iend,nlev+1) = zrgp(1:il,ifs_config%icdir,ib)
      flux_sw_direct_normal(ibeg:iend) = zrgp(1:il,ifs_config%isudu,ib)
      flux_uv(ibeg:iend) = zrgp(1:il,ifs_config%iuvdf,ib)
      flux_par(ibeg:iend) = zrgp(1:il,ifs_config%iparf,ib)
      flux_par_clear(ibeg:iend) = zrgp(1:il,ifs_config%iparcf,ib)
      flux%sw_dn(ibeg:iend,1) = zrgp(1:il,ifs_config%itincf,ib)
      emissivity_out(ibeg:iend) = zrgp(1:il,ifs_config%iemit,ib)
      if (yradiation%yrerad%lapproxswupdate) then
        do jg=1,yradiation%yrerad%nsw
          flux_diffuse_band(ibeg:iend,jg) = zrgp(1:il,ifs_config%iswdiffuseband+jg-1,ib)
          flux_direct_band(ibeg:iend,jg) = zrgp(1:il,ifs_config%iswdirectband+jg-1,ib)
        end do
      else
        flux_diffuse_band(ibeg:iend,:) = 0.0_jprb
        flux_direct_band(ibeg:iend,:) = 0.0_jprb
      endif
    end do

    deallocate(zrgp)

end subroutine ifs_copy_fluxes_from_blocked

end module ifs_blocking
