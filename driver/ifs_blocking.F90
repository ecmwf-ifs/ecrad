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
#ifdef BITIDENTITY_TESTING
    integer :: ire_liq, ire_ice, ioverlap
#endif
    integer :: ifldstot, iinbeg, ioutbeg, iinend, ioutend
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
#ifdef BITIDENTITY_TESTING
  ! To validate results against standalone ecrad, we overwrite effective
  ! radii, cloud overlap and seed with input values
  ifs_config%ire_liq =indrad(inext,nlev,.true.)
  ifs_config%ire_ice =indrad(inext,nlev,.true.)
  ifs_config%ioverlap =indrad(inext,nlev-1,.true.)
#endif
                                  ! end of standalone inputs workaround variables

  ifldsin = iinend - iinbeg +1
  ifldsout= ioutend-ioutbeg +1
  ifs_config%ifldstot= inext  - 1
  ifs_config%iinbeg = iinbeg
  ifs_config%iinend = iinend
  ifs_config%ioutbeg = ioutbeg
  ifs_config%ioutend = ioutend

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
#ifdef BITIDENTITY_TESTING
    write(nulout,'("ire_liq=",i0)')ifs_config%ire_liq
    write(nulout,'("ire_ice=",i0)')ifs_config%ire_ice
    write(nulout,'("ioverlap=",i0)')ifs_config%ioverlap
#endif
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
  use radiation_io,             only : nulout

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
      call gas%get(gas, ICO2, IMassMixingRatio, zrgp(1:il,ifs_config%iico2:ifs_config%iico2+nlev-1,ib), istartcol=ibeg)
      call gas%get(gas, ICH4, IMassMixingRatio, zrgp(1:il,ifs_config%iich4:ifs_config%iich4+nlev-1,ib), istartcol=ibeg)
      call gas%get(gas, IN2O, IMassMixingRatio, zrgp(1:il,ifs_config%iin2o:ifs_config%iin2o+nlev-1,ib), istartcol=ibeg)
      call gas%get(gas, ICFC11, IMassMixingRatio, zrgp(1:il,ifs_config%ic11:ifs_config%ic11+nlev-1,ib), istartcol=ibeg)
      call gas%get(gas, ICFC12, IMassMixingRatio, zrgp(1:il,ifs_config%ic12:ifs_config%ic12+nlev-1,ib), istartcol=ibeg)
      call gas%get(gas, IHCFC22,IMassMixingRatio, zrgp(1:il,ifs_config%ic22:ifs_config%ic22+nlev-1,ib), istartcol=ibeg)
      call gas%get(gas, ICCL4,  IMassMixingRatio, zrgp(1:il,ifs_config%icl4:ifs_config%icl4+nlev-1,ib), istartcol=ibeg)
      call gas%get(gas, IO3, IMassMixingRatio, zrgp(1:il,ifs_config%ioz:ifs_config%ioz+nlev-1,ib), istartcol=ibeg)
      ! convert ozone kg/kg to Pa*kg/kg
      ! do jlev=1,nlev
      !   zrgp(1:il,ifs_config%ioz+jlev-1,ib)  = zrgp(1:il,ifs_config%ioz+jlev-1,ib) &
      !         &                       * (thermodynamics%pressure_hl(ibeg:iend,jlev+1) &
      !         &                         - thermodynamics%pressure_hl(ibeg:iend,jlev))
      ! enddo

      ! local workaround variables for standalone input files
#ifdef BITIDENTITY_TESTING
      ! To validate results against standalone ecrad, we overwrite effective
      ! radii, cloud overlap and seed with input values
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
#endif
    enddo
    !$OMP END PARALLEL DO

    ! Store pressure for output
    if(present(thermodynamics_out)) thermodynamics_out%pressure_hl(:,:) = thermodynamics%pressure_hl(:,:)

  end associate

end subroutine ifs_copy_inputs_to_blocked


subroutine ifs_copy_inputs_to_blocked_new ( driver_config, ifs_config, &
  &  yradiation, single_level, ncol, nlev, jrl, zrgp, &
  &  PMU0, PTEMPERATURE_SKIN, PALBEDO_DIF, PALBEDO_DIR, &
  &  PSPECTRALEMISS, &
  &  PCCN_LAND, PCCN_SEA, &
  &  PGELAM, PGEMU, PLAND_SEA_MASK, &
  &  PPRESSURE, PTEMPERATURE, &
  &  PPRESSURE_H, PTEMPERATURE_H, &
  &  PQ, PCO2, PCH4, PN2O, PNO2, PCFC11, PCFC12, PHCFC22, PCCL4, PO3, &
  &  PCLOUD_FRAC, PQ_LIQUID, PQ_ICE, PQ_RAIN, PQ_SNOW, &
  &  PAEROSOL_OLD, PAEROSOL, &
  &  PFLUX_SW, PFLUX_LW, PFLUX_SW_CLEAR, PFLUX_LW_CLEAR, &
  &  PFLUX_SW_DN, PFLUX_LW_DN, PFLUX_SW_DN_CLEAR, PFLUX_LW_DN_CLEAR, &
  &  PFLUX_DIR, PFLUX_DIR_CLEAR, PFLUX_DIR_INTO_SUN, &
  &  PFLUX_UV, PFLUX_PAR, PFLUX_PAR_CLEAR, &
  &  PFLUX_SW_DN_TOA, PEMIS_OUT, PLWDERIVATIVE, &
  &  PSWDIFFUSEBAND, PSWDIRECTBAND, &
  ! OPTIONAL ARGUMENTS for bit-identical results in tests
  &  iseed, PRE_LIQ, PRE_ICE, PCLOUD_OVERLAP)

  use radiation_single_level,   only : single_level_type
  use ecrad_driver_config,      only : driver_config_type
  use radiation_setup,          only : tradiation
  use radiation_io,             only : nulout

  implicit none

  ! Configuration specific to this driver
  type(driver_config_type), intent(in)     :: driver_config

  type(ifs_config_type), intent(in)     :: ifs_config

  ! Configuration for the radiation scheme, IFS style
  type(tradiation), intent(in)          :: yradiation

  integer, intent(in) :: ncol, nlev         ! Number of columns and levels

  ! Derived types for the inputs to the radiation scheme
  type(single_level_type), intent(in)   :: single_level

  ! monolithic IFS data structure to pass to radiation scheme
  real(kind=jprb), intent(in), allocatable :: zrgp(:,:,:)

  ! Seed for random number generator
  integer, intent(inout), allocatable, optional :: iseed(:)

  integer, intent(in) :: jrl

  real(kind=jprb), intent(inout), allocatable :: PMU0(:)                        !(KLON, ) ! Cosine of solar zenith ang
  real(kind=jprb), intent(inout), allocatable :: PTEMPERATURE_SKIN(:)           !(KLON) ! (K)
  real(kind=jprb), intent(inout), allocatable :: PALBEDO_DIF(:,:)               !(KLON,YRADIATION%YRERAD%NSW)
  real(kind=jprb), intent(inout), allocatable :: PALBEDO_DIR(:,:)               !(KLON,YRADIATION%YRERAD%NSW)
  real(kind=jprb), intent(inout), allocatable :: PSPECTRALEMISS(:,:)            !(KLON,YRADIATION%YRERAD%NLWEMISS)
  real(kind=jprb), intent(inout), allocatable :: PCCN_LAND(:)!(KLON)
  real(kind=jprb), intent(inout), allocatable :: PCCN_SEA(:)!(KLON)
  real(kind=jprb), intent(inout), allocatable :: PGELAM(:)                      !(KLON)
  real(kind=jprb), intent(inout), allocatable :: PGEMU(:)                       !(KLON)
  real(kind=jprb), intent(inout), allocatable :: PLAND_SEA_MASK(:)              !((KLON)
  real(kind=jprb), intent(inout), allocatable :: PPRESSURE(:,:)                 !((KLON,KLEV)    ! (Pa)
  real(kind=jprb), intent(inout), allocatable :: PTEMPERATURE(:,:)              !(KLON,KLEV) ! (K)
  real(kind=jprb), intent(inout), allocatable :: PPRESSURE_H(:,:)               !(KLON,KLEV+1)    ! (Pa)
  real(kind=jprb), intent(inout), allocatable :: PTEMPERATURE_H(:,:)            !(KLON,KLEV+1) ! (K)
  real(kind=jprb), intent(inout), allocatable :: PQ(:,:)!(KLON,KLEV)
  real(kind=jprb), intent(inout), allocatable :: PCO2(:,:)!(KLON,KLEV)
  real(kind=jprb), intent(inout), allocatable :: PCH4(:,:)!(KLON,KLEV)
  real(kind=jprb), intent(inout), allocatable :: PN2O(:,:)!(KLON,KLEV)
  real(kind=jprb), intent(inout), allocatable :: PNO2(:,:)!(KLON,KLEV)
  real(kind=jprb), intent(inout), allocatable :: PCFC11(:,:)!(KLON,KLEV)
  real(kind=jprb), intent(inout), allocatable :: PCFC12(:,:)!(KLON,KLEV)
  real(kind=jprb), intent(inout), allocatable :: PHCFC22(:,:)!(KLON,KLEV)
  real(kind=jprb), intent(inout), allocatable :: PCCL4(:,:)!(KLON,KLEV)
  real(kind=jprb), intent(inout), allocatable :: PO3(:,:)!(KLON,KLEV)
  real(kind=jprb), intent(inout), allocatable :: PCLOUD_FRAC(:,:)!(KLON,KLEV)
  real(kind=jprb), intent(inout), allocatable :: PQ_LIQUID(:,:)!(KLON,KLEV)
  real(kind=jprb), intent(inout), allocatable :: PQ_ICE(:,:)!(KLON,KLEV)
  real(kind=jprb), intent(inout), allocatable :: PQ_RAIN(:,:)!(KLON,KLEV)
  real(kind=jprb), intent(inout), allocatable :: PQ_SNOW(:,:)!(KLON,KLEV)
  real(kind=jprb), intent(inout), allocatable :: PAEROSOL_OLD(:,:,:)!(KLON,6,KLEV)
  real(kind=jprb), intent(inout), allocatable :: PAEROSOL(:,:,:)!(KLON,KLEV,KAEROSOL)
  ! OUT
  real(kind=jprb), intent(inout), allocatable :: PFLUX_SW(:,:)!(KLON,KLEV+1)
  real(kind=jprb), intent(inout), allocatable :: PFLUX_LW(:,:)!(KLON,KLEV+1)
  real(kind=jprb), intent(inout), allocatable :: PFLUX_SW_CLEAR(:,:)!(KLON,KLEV+1)
  real(kind=jprb), intent(inout), allocatable :: PFLUX_LW_CLEAR(:,:)!(KLON,KLEV+1)
  real(kind=jprb), intent(inout), allocatable :: PFLUX_SW_DN(:)!(KLON)
  real(kind=jprb), intent(inout), allocatable :: PFLUX_LW_DN(:)!(KLON)
  real(kind=jprb), intent(inout), allocatable :: PFLUX_SW_DN_CLEAR(:)!(KLON)
  real(kind=jprb), intent(inout), allocatable :: PFLUX_LW_DN_CLEAR(:)!(KLON)
  real(kind=jprb), intent(inout), allocatable :: PFLUX_DIR(:)!(KLON)
  real(kind=jprb), intent(inout), allocatable :: PFLUX_DIR_CLEAR(:)!(KLON)
  real(kind=jprb), intent(inout), allocatable :: PFLUX_DIR_INTO_SUN(:)!(KLON)
  real(kind=jprb), intent(inout), allocatable :: PFLUX_UV(:)!(KLON)
  real(kind=jprb), intent(inout), allocatable :: PFLUX_PAR(:)!(KLON)
  real(kind=jprb), intent(inout), allocatable :: PFLUX_PAR_CLEAR(:)!(KLON)
  real(kind=jprb), intent(inout), allocatable :: PFLUX_SW_DN_TOA(:)!(KLON)
  real(kind=jprb), intent(inout), allocatable :: PEMIS_OUT(:)!(KLON)
  real(kind=jprb), intent(inout), allocatable :: PLWDERIVATIVE(:,:)!(KLON,KLEV+1)
  real(kind=jprb), intent(inout), allocatable :: PSWDIFFUSEBAND(:,:)!(KLON,YRADIATION%YRERAD%NSW)
  real(kind=jprb), intent(inout), allocatable :: PSWDIRECTBAND(:,:)!(KLON,YRADIATION%YRERAD%NSW)
  real(kind=jprb), optional, intent(inout), allocatable :: PRE_LIQ(:,:)!(KLON, KLEV)
  real(kind=jprb), optional, intent(inout), allocatable :: PRE_ICE(:,:)!(KLON, KLEV)
  real(kind=jprb), optional, intent(inout), allocatable :: PCLOUD_OVERLAP(:,:)!(KLON, KLEV-1)

  ! number of column blocks, block size
  integer :: ngpblks, nproma
  integer :: ibeg, iend, il, ib, ifld, jemiss, jalb, jlev, joff, jaer

  ! Extract some config values
  nproma=driver_config%nblocksize        ! nproma size
  ngpblks=(ncol-1)/nproma+1              ! number of column blocks

  associate(yderad=>yradiation%yrerad, rad_config=>yradiation%rad_config)

  PMU0(:)=0._jprb
  PTEMPERATURE_SKIN(:)=0._jprb
  PALBEDO_DIF(:,:)=0._jprb
  PALBEDO_DIR(:,:)=0._jprb
  PSPECTRALEMISS(:,:)=0._jprb
  PGELAM(:)=0._jprb
  PGEMU(:)=0._jprb
  PLAND_SEA_MASK(:)=0._jprb
  PPRESSURE(:,:)=0._jprb
  PTEMPERATURE(:,:)=0._jprb
  PPRESSURE_H(:,:)=0._jprb
  PTEMPERATURE_H(:,:)=0._jprb
  PQ(:,:)=0._jprb
  PCO2(:,:)=0._jprb
  PCH4(:,:)=0._jprb
  PN2O(:,:)=0._jprb
  PNO2(:,:)=0._jprb
  PCFC11(:,:)=0._jprb
  PCFC12(:,:)=0._jprb
  PHCFC22(:,:)=0._jprb
  PCCL4(:,:)=0._jprb
  PO3(:,:)=0._jprb
  PCLOUD_FRAC(:,:)=0._jprb
  PQ_LIQUID(:,:)=0._jprb
  PQ_ICE(:,:)=0._jprb
  PQ_RAIN(:,:)=0._jprb
  PQ_SNOW(:,:)=0._jprb
  PAEROSOL_OLD(:,:,:)=0._jprb
  PAEROSOL(:,:,:)=0._jprb
  PCCN_LAND(:)=0._jprb
  PCCN_SEA(:)=0._jprb

  PNO2(:,:)=0._jprb
  PFLUX_SW(:,:)=0._jprb
  PFLUX_LW(:,:)=0._jprb
  PFLUX_SW_CLEAR(:,:)=0._jprb
  PFLUX_LW_CLEAR(:,:)=0._jprb
  PFLUX_SW_DN(:)=0._jprb
  PFLUX_LW_DN(:)=0._jprb
  PFLUX_SW_DN_CLEAR(:)=0._jprb
  PFLUX_LW_DN_CLEAR(:)=0._jprb
  PFLUX_DIR(:)=0._jprb
  PFLUX_DIR_CLEAR(:)=0._jprb
  PFLUX_DIR_INTO_SUN(:)=0._jprb
  PFLUX_UV(:)=0._jprb
  PFLUX_PAR(:)=0._jprb
  PFLUX_PAR_CLEAR(:)=0._jprb
  PFLUX_SW_DN_TOA(:)=0._jprb
  PEMIS_OUT(:)=0._jprb
  PLWDERIVATIVE(:,:)=0._jprb
  PSWDIFFUSEBAND(:,:)=0._jprb
  PSWDIRECTBAND(:,:)=0._jprb
#ifdef BITIDENTITY_TESTING
  PRE_LIQ(:,:)=0._jprb
  PRE_ICE(:,:)=0._jprb
  PCLOUD_OVERLAP(:,:)=0._jprb
#endif


  ibeg=jrl
  iend=min(ibeg+nproma-1,ncol)
  il=iend-ibeg+1
  ib=(jrl-1)/nproma+1
  PMU0(1:il)                    =  zrgp(1:il,ifs_config%iamu0,ib)   ! cosine of solar zenith ang (mu0)

  do jemiss=1,yderad%nlwemiss
     PSPECTRALEMISS(1:il,jemiss)             =  zrgp(1:il,ifs_config%iemiss+jemiss-1,ib)
  enddo

  PTEMPERATURE_SKIN(1:il)         = zrgp(1:il,ifs_config%its,ib)
  PLAND_SEA_MASK(1:il)            = zrgp(1:il,ifs_config%islm,ib)
  PCCN_LAND(1:il)                 = zrgp(1:il,ifs_config%iccnl,ib)
  PCCN_SEA(1:il)                  = zrgp(1:il,ifs_config%iccno,ib)
  PGELAM(1:il)                    = zrgp(1:il,ifs_config%igelam,ib)
  PGEMU(1:il)                     = zrgp(1:il,ifs_config%igemu,ib)
  do jalb=1,yderad%nsw
     PALBEDO_DIF(1:il,jalb)             =  zrgp(1:il,ifs_config%iald+jalb-1,ib)
  enddo

  if (allocated(single_level%sw_albedo_direct)) then
     do jalb=1,yderad%nsw
        PALBEDO_DIR(1:il,jalb)             =  zrgp(1:il,ifs_config%ialp+jalb-1,ib)
     end do
  else
     do jalb=1,yderad%nsw
        PALBEDO_DIR(1:il,jalb)             =  zrgp(1:il,ifs_config%ialp+jalb-1,ib)
     end do
  end if
  
  do jlev=1,nlev
     PTEMPERATURE(1:il,jlev)            = zrgp(1:il,ifs_config%iti+jlev-1,ib)
     PPRESSURE(1:il,jlev)               = zrgp(1:il,ifs_config%ipr+jlev-1,ib)
  enddo

  do jlev=1,nlev
     PQ(1:il,jlev)                          = zrgp(1:il,ifs_config%iwv+jlev-1,ib)
     if (rad_config%do_clouds) then
        PCLOUD_FRAC(1:il,jlev)              = zrgp(1:il,ifs_config%iclc+jlev-1,ib)
        PQ_LIQUID(1:il,jlev)                = zrgp(1:il,ifs_config%ilwa+jlev-1,ib)
        PQ_ICE(1:il,jlev)                   = zrgp(1:il,ifs_config%iiwa+jlev-1,ib)
     else
        PCLOUD_FRAC(1:il,jlev)             = 0._jprb
        PQ_LIQUID(1:il,jlev)             = 0._jprb
        PQ_ICE(1:il,jlev)             = 0._jprb
     endif
     PQ_SNOW(1:il,jlev)                  = 0._jprb  ! snow
     PQ_RAIN(1:il,jlev)                  = 0._jprb  ! rain
  enddo

  PAEROSOL_OLD(:,:,:) = 0._jprb
  if (yderad%naermacc == 1) then
     joff=ifs_config%iaero
     do jaer=1,rad_config%n_aerosol_types
        do jlev=1,nlev
           PAEROSOL(1:il,jlev,jaer) = zrgp(1:il,joff,ib)
           joff=joff+1
        enddo
     enddo
  endif

  do jlev=1,nlev+1
     PPRESSURE_H(1:il,jlev)             = zrgp(1:il,ifs_config%iaprs+jlev-1,ib)
     PTEMPERATURE_H(1:il,jlev)          = zrgp(1:il,ifs_config%ihti+jlev-1,ib)
  enddo

  PCO2(1:il,1:nlev) = zrgp(1:il,ifs_config%iico2:ifs_config%iico2+nlev-1,ib)
  PCH4(1:il,1:nlev) = zrgp(1:il,ifs_config%iich4:ifs_config%iich4+nlev-1,ib)
  PN2O(1:il,1:nlev) = zrgp(1:il,ifs_config%iin2o:ifs_config%iin2o+nlev-1,ib)
  PCFC11(1:il,1:nlev) = zrgp(1:il,ifs_config%ic11:ifs_config%ic11+nlev-1,ib)
  PCFC12(1:il,1:nlev) = zrgp(1:il,ifs_config%ic12:ifs_config%ic12+nlev-1,ib)
  PHCFC22(1:il,1:nlev) = zrgp(1:il,ifs_config%ic22:ifs_config%ic22+nlev-1,ib)
  PCCL4(1:il,1:nlev) = zrgp(1:il,ifs_config%icl4:ifs_config%icl4+nlev-1,ib)
  PO3(1:il,1:nlev) = zrgp(1:il,ifs_config%ioz:ifs_config%ioz+nlev-1,ib)

  ! local workaround variables for standalone input files
#ifdef BITIDENTITY_TESTING
  ! To validate results against standalone ecrad, we overwrite effective
  ! radii, cloud overlap and seed with input values
  if (rad_config%do_clouds) then
     do jlev=1,nlev
        ! missing full-level temperature and pressure as well as land-sea-mask
        PRE_LIQ(1:il,jlev) = zrgp(1:il,ifs_config%ire_liq+jlev-1,ib)
        PRE_ICE(1:il,jlev) = zrgp(1:il,ifs_config%ire_ice+jlev-1,ib)
     enddo
     do jlev=1,nlev-1
        ! for the love of it, I can't figure this one out. Probably to do with
        ! my crude approach of setting PGEMU?
        PCLOUD_OVERLAP(1:il,jlev) = zrgp(1:il,ifs_config%ioverlap+jlev-1,ib)
     enddo
     if(present(iseed)) iseed(1:il) = single_level%iseed(ibeg:iend)
          
  else
     do jlev=1,nlev
        ! missing full-level temperature and pressure as well as land-sea-mask
        PRE_LIQ(1:il,jlev) = 0._jprb
        PRE_ICE(1:il,jlev) = 0._jprb
     enddo
     do jlev=1,nlev-1
        PCLOUD_OVERLAP(1:il,jlev) = 0._jprb
     enddo
     if(present(iseed)) iseed(1:il) = 0
  endif ! do_clouds
#endif

end associate

end subroutine ifs_copy_inputs_to_blocked_new

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

subroutine ifs_copy_fluxes_from_blocked_new(&
     & driver_config, ifs_config, yradiation, ncol, nlev,jrl,&
     & PFLUX_SW, PFLUX_LW, PFLUX_SW_CLEAR, PFLUX_LW_CLEAR, &
     & PFLUX_SW_DN, PFLUX_LW_DN, PFLUX_SW_DN_CLEAR, PFLUX_LW_DN_CLEAR, &
     & PFLUX_DIR, PFLUX_DIR_CLEAR, PFLUX_DIR_INTO_SUN, &
     & PFLUX_UV, PFLUX_PAR, PFLUX_PAR_CLEAR, &
     & PFLUX_SW_DN_TOA, PEMIS_OUT, PLWDERIVATIVE, &
     & PSWDIFFUSEBAND, PSWDIRECTBAND, zrgp)
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
  !real(kind=jprb), intent(inout), allocatable :: zrgp(:,:,:)

  real(kind=jprb), intent(inout), allocatable :: PFLUX_SW(:,:)!(KLON,KLEV+1)
  real(kind=jprb), intent(inout), allocatable :: PFLUX_LW(:,:)!(KLON,KLEV+1)
  real(kind=jprb), intent(inout), allocatable :: PFLUX_SW_CLEAR(:,:)!(KLON,KLEV+1)
  real(kind=jprb), intent(inout), allocatable :: PFLUX_LW_CLEAR(:,:)!(KLON,KLEV+1)
  real(kind=jprb), intent(inout), allocatable :: PFLUX_SW_DN(:)!(KLON)
  real(kind=jprb), intent(inout), allocatable :: PFLUX_LW_DN(:)!(KLON)
  real(kind=jprb), intent(inout), allocatable :: PFLUX_SW_DN_CLEAR(:)!(KLON)
  real(kind=jprb), intent(inout), allocatable :: PFLUX_LW_DN_CLEAR(:)!(KLON)
  real(kind=jprb), intent(inout), allocatable :: PFLUX_DIR(:)!(KLON)
  real(kind=jprb), intent(inout), allocatable :: PFLUX_DIR_CLEAR(:)!(KLON)
  real(kind=jprb), intent(inout), allocatable :: PFLUX_DIR_INTO_SUN(:)!(KLON)
  real(kind=jprb), intent(inout), allocatable :: PFLUX_UV(:)!(KLON)
  real(kind=jprb), intent(inout), allocatable :: PFLUX_PAR(:)!(KLON)
  real(kind=jprb), intent(inout), allocatable :: PFLUX_PAR_CLEAR(:)!(KLON)
  real(kind=jprb), intent(inout), allocatable :: PFLUX_SW_DN_TOA(:)!(KLON)
  real(kind=jprb), intent(inout), allocatable :: PEMIS_OUT(:)!(KLON)
  real(kind=jprb), intent(inout), allocatable :: PLWDERIVATIVE(:,:)!(KLON,KLEV+1)
  real(kind=jprb), intent(inout), allocatable :: PSWDIFFUSEBAND(:,:)!(KLON,YRADIATION%YRERAD%NSW)
  real(kind=jprb), intent(inout), allocatable :: PSWDIRECTBAND(:,:)!(KLON,YRADIATION%YRERAD%NSW)

  ! monolithic IFS data structure passed to radiation scheme
  real(kind=jprb), intent(inout), allocatable :: zrgp(:,:,:)
  integer, intent(in) :: jrl

  ! number of column blocks, block size
  integer :: ngpblks, nproma
  integer :: ibeg, iend, il, ib, jlev, jg

  ! Extract some config values
  nproma=driver_config%nblocksize        ! nproma size
  ngpblks=(ncol-1)/nproma+1              ! number of column blocks

  !  -------------------------------------------------------
  !
  !  OUTPUT LOOP
  !
  !  -------------------------------------------------------
  
  ibeg=jrl
  iend=min(ibeg+nproma-1,ncol)
  il=iend-ibeg+1
  ib=(jrl-1)/nproma+1

  do jlev=1,nlev+1
     zrgp(1:il,ifs_config%ifrso+jlev-1,ib) = PFLUX_SW(1:il,jlev)
     zrgp(1:il,ifs_config%ifrth+jlev-1,ib) = PFLUX_LW(1:il,jlev)
     zrgp(1:il,ifs_config%iswfc+jlev-1,ib) = PFLUX_SW_CLEAR(1:il,jlev)
     zrgp(1:il,ifs_config%ilwfc+jlev-1,ib) = PFLUX_LW_CLEAR(1:il,jlev)
     if (yradiation%yrerad%lapproxlwupdate) then
        zrgp(1:il,ifs_config%ilwderivative+jlev-1,ib) = PLWDERIVATIVE(1:il,jlev)
     else
        zrgp(1:il,ifs_config%ilwderivative+jlev-1,ib) = 0.0_jprb
     endif
  end do
  
  zrgp(1:il,ifs_config%ifrsod,ib) = PFLUX_SW_DN(1:il)
  zrgp(1:il,ifs_config%ifrted,ib) = PFLUX_LW_DN(1:il)
  zrgp(1:il,ifs_config%ifrsodc,ib) = PFLUX_SW_DN_CLEAR(1:il)
  zrgp(1:il,ifs_config%ifrtedc,ib) = PFLUX_LW_DN_CLEAR(1:il)
  zrgp(1:il,ifs_config%ifdir,ib) = PFLUX_DIR(1:il)
  zrgp(1:il,ifs_config%icdir,ib) = PFLUX_DIR_CLEAR(1:il)
  zrgp(1:il,ifs_config%isudu,ib) = PFLUX_DIR_INTO_SUN(1:il)
  zrgp(1:il,ifs_config%iuvdf,ib) = PFLUX_UV(1:il)
  zrgp(1:il,ifs_config%iparf,ib) = PFLUX_PAR(1:il)
  zrgp(1:il,ifs_config%iparcf,ib) = PFLUX_PAR_CLEAR(1:il)
  zrgp(1:il,ifs_config%itincf,ib) = PFLUX_SW_DN_TOA(1:il)
  zrgp(1:il,ifs_config%iemit,ib) = PEMIS_OUT(1:il)
  if (yradiation%yrerad%lapproxswupdate) then
     do jg=1,yradiation%yrerad%nsw
        zrgp(1:il,ifs_config%iswdiffuseband+jg-1,ib) = PSWDIFFUSEBAND(1:il,jg)
        zrgp(1:il,ifs_config%iswdirectband+jg-1,ib) = PSWDIRECTBAND(1:il,jg)
     end do
  endif

end subroutine ifs_copy_fluxes_from_blocked_new

end module ifs_blocking
