! radintg_zrgp_mod.fypp - Wrap the block-allocated radiation fields from RADINTG in a FIELD API stack
!
! (C) Copyright 2022- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
!
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.


MODULE RADINTG_ZRGP_MOD

USE PARKIND1     , ONLY : JPRB, JPIM
USE YOMHOOK      , ONLY : LHOOK, DR_HOOK, JPHOOK

IMPLICIT NONE

PRIVATE
PUBLIC :: RADINTG_ZRGP_TYPE

#ifdef BITIDENTITY_TESTING
  LOGICAL, PARAMETER :: LBITIDENTITY = .TRUE.
#else
  LOGICAL, PARAMETER :: LBITIDENTITY = .FALSE.
#endif

TYPE RADINTG_ZRGP_TYPE
  ! Field counts and offset indices for ZRGP
  INTEGER(KIND=JPIM) :: IFLDSIN, IFLDSOUT, IFLDSTOT
  INTEGER(KIND=JPIM) :: IINBEG, IINEND, IOUTBEG, IOUTEND
  INTEGER(KIND=JPIM) :: igi
  INTEGER(KIND=JPIM) :: imu0
  INTEGER(KIND=JPIM) :: iamu0
  INTEGER(KIND=JPIM) :: iemiss
  INTEGER(KIND=JPIM) :: its
  INTEGER(KIND=JPIM) :: islm
  INTEGER(KIND=JPIM) :: iccnl
  INTEGER(KIND=JPIM) :: iccno
  INTEGER(KIND=JPIM) :: ibas
  INTEGER(KIND=JPIM) :: itop
  INTEGER(KIND=JPIM) :: igelam
  INTEGER(KIND=JPIM) :: igemu
  INTEGER(KIND=JPIM) :: iclon
  INTEGER(KIND=JPIM) :: islon
  INTEGER(KIND=JPIM) :: iald
  INTEGER(KIND=JPIM) :: ialp
  INTEGER(KIND=JPIM) :: iti
  INTEGER(KIND=JPIM) :: ipr
  INTEGER(KIND=JPIM) :: iqs
  INTEGER(KIND=JPIM) :: iwv
  INTEGER(KIND=JPIM) :: iclc
  INTEGER(KIND=JPIM) :: ilwa
  INTEGER(KIND=JPIM) :: iiwa
  INTEGER(KIND=JPIM) :: iswa
  INTEGER(KIND=JPIM) :: irwa
  INTEGER(KIND=JPIM) :: irra
  INTEGER(KIND=JPIM) :: idp
  INTEGER(KIND=JPIM) :: ioz
  INTEGER(KIND=JPIM) :: iecpo3
  INTEGER(KIND=JPIM) :: ihpr
  INTEGER(KIND=JPIM) :: iaprs
  INTEGER(KIND=JPIM) :: ihti
  INTEGER(KIND=JPIM) :: ire_liq
  INTEGER(KIND=JPIM) :: ire_ice
  INTEGER(KIND=JPIM) :: ioverlap
  INTEGER(KIND=JPIM) :: iaero
  INTEGER(KIND=JPIM) :: ifrsod
  INTEGER(KIND=JPIM) :: ifrted
  INTEGER(KIND=JPIM) :: ifrsodc
  INTEGER(KIND=JPIM) :: ifrtedc
  INTEGER(KIND=JPIM) :: iemit
  INTEGER(KIND=JPIM) :: isudu
  INTEGER(KIND=JPIM) :: iuvdf
  INTEGER(KIND=JPIM) :: iparf
  INTEGER(KIND=JPIM) :: iparcf
  INTEGER(KIND=JPIM) :: itincf
  INTEGER(KIND=JPIM) :: ifdir
  INTEGER(KIND=JPIM) :: ifdif
  INTEGER(KIND=JPIM) :: icdir
  INTEGER(KIND=JPIM) :: ilwderivative
  INTEGER(KIND=JPIM) :: iswdirectband
  INTEGER(KIND=JPIM) :: iswdiffuseband
  INTEGER(KIND=JPIM) :: ifrso
  INTEGER(KIND=JPIM) :: iswfc
  INTEGER(KIND=JPIM) :: ifrth
  INTEGER(KIND=JPIM) :: ilwfc
  INTEGER(KIND=JPIM) :: iaer
  INTEGER(KIND=JPIM) :: iico2
  INTEGER(KIND=JPIM) :: iich4
  INTEGER(KIND=JPIM) :: iin2o
  INTEGER(KIND=JPIM) :: ino2
  INTEGER(KIND=JPIM) :: ic11
  INTEGER(KIND=JPIM) :: ic12
  INTEGER(KIND=JPIM) :: ic22
  INTEGER(KIND=JPIM) :: icl4
  INTEGER(KIND=JPIM) :: igix

CONTAINS
  PROCEDURE :: SETUP => RADINTG_ZRGP_SETUP

END TYPE RADINTG_ZRGP_TYPE

CONTAINS

INTEGER(KIND=JPIM) FUNCTION INDRAD(KNEXT,KFLDS,LDUSE)
INTEGER(KIND=JPIM),INTENT(INOUT) :: KNEXT
INTEGER(KIND=JPIM),INTENT(IN) :: KFLDS
LOGICAL,INTENT(IN) :: LDUSE
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('RADINTG:INDRAD',0,ZHOOK_HANDLE)

IF( LDUSE )THEN
  INDRAD=KNEXT
  KNEXT=KNEXT+KFLDS
ELSE
  INDRAD=-99999999
ENDIF

IF (LHOOK) CALL DR_HOOK('RADINTG:INDRAD',1,ZHOOK_HANDLE)

END FUNCTION INDRAD


SUBROUTINE RADINTG_ZRGP_SETUP( &
  & SELF, NLEV, NLWEMISS, &
  & NLWOUT, NSW, NFSD, NRFTOTAL_RADGRID, &
  & NPROGAER, NRADAER, &
  & LDEBUG, LSPPRAD, LRAYFM, &
  & LAPPROXLWUPDATE, LAPPROXSWUPDATE, &
  & LEPO3RA, LDIAGFORCING)

  USE YOMLUN       , ONLY : NULOUT

  IMPLICIT NONE

  CLASS(RADINTG_ZRGP_TYPE), INTENT(INOUT):: SELF
  INTEGER, INTENT(IN)                    :: NLEV
  INTEGER, INTENT(IN)                    :: NLWEMISS, NLWOUT, NSW
  INTEGER, INTENT(IN)                    :: NFSD, NRFTOTAL_RADGRID
  INTEGER, INTENT(IN)                    :: NPROGAER, NRADAER
  LOGICAL, INTENT(IN)                    :: LDEBUG, LSPPRAD, LRAYFM
  LOGICAL, INTENT(IN)                    :: LAPPROXLWUPDATE, LAPPROXSWUPDATE
  LOGICAL, INTENT(IN)                    :: LEPO3RA, LDIAGFORCING

  INTEGER(KIND=JPIM), ALLOCATABLE     :: MEMBER_MAP(:)
  INTEGER(KIND=JPIM)                  :: INEXT
  REAL(KIND=JPHOOK)                   :: ZHOOK_HANDLE

  IF (LHOOK) CALL DR_HOOK('RADINTG_ZRGP_SETUP',0,ZHOOK_HANDLE)

  INEXT = 1
  SELF%IINBEG=1
  SELF%igi = INDRAD( INEXT, 1, ldebug)
  SELF%imu0 = INDRAD( INEXT, 1, .true.)
  SELF%iamu0 = INDRAD( INEXT, 1, .true.)
  SELF%iemiss = INDRAD( INEXT, nlwemiss, .true.)
  SELF%its = INDRAD( INEXT, 1, .true.)
  SELF%islm = INDRAD( INEXT, 1, .true.)
  SELF%iccnl = INDRAD( INEXT, 1, .true.)
  SELF%iccno = INDRAD( INEXT, 1, .true.)
  SELF%ibas = INDRAD( INEXT, 1, .true.)
  SELF%itop = INDRAD( INEXT, 1, .true.)
  SELF%igelam = INDRAD( INEXT, 1, .true.)
  SELF%igemu = INDRAD( INEXT, 1, .true.)
  SELF%iclon = INDRAD( INEXT, 1, .true.)
  SELF%islon = INDRAD( INEXT, 1, .true.)
  SELF%iald = INDRAD( INEXT, nsw, .true.)
  SELF%ialp = INDRAD( INEXT, nsw, .true.)
  SELF%iti = INDRAD( INEXT, nlev, .true.)
  SELF%ipr = INDRAD( INEXT, nlev, .true.)
  SELF%iqs = INDRAD( INEXT, nlev, .true.)
  SELF%iwv = INDRAD( INEXT, nlev, .true.)
  SELF%iclc = INDRAD( INEXT, nlev, .true.)
  SELF%ilwa = INDRAD( INEXT, nlev, .true.)
  SELF%iiwa = INDRAD( INEXT, nlev, .true.)
  SELF%iswa = INDRAD( INEXT, nlev, .true.)
  SELF%irwa = INDRAD( INEXT, nlev, .true.)
  SELF%irra = INDRAD( INEXT, nlev, .true.)
  SELF%idp = INDRAD( INEXT, nlev, .true.)
  SELF%ioz = INDRAD( INEXT, nlev, lrayfm)
  SELF%iecpo3 = INDRAD( INEXT, nlev, .not.lrayfm.and.lepo3ra)
  SELF%ihpr = INDRAD( INEXT, nlev+1, .true.)
  SELF%iaprs = INDRAD( INEXT, nlev+1, .true.)
  SELF%ihti = INDRAD( INEXT, nlev+1, .true.)
  SELF%ire_liq = INDRAD( INEXT, nlev, lbitidentity)
  SELF%ire_ice = INDRAD( INEXT, nlev, lbitidentity)
  SELF%ioverlap = INDRAD( INEXT, nlev-1, lbitidentity)
  SELF%IINEND = INEXT-1
  SELF%IOUTBEG = INEXT
  SELF%iaero = INDRAD( INEXT, nradaer*nlev, .true.)
  SELF%ifrsod = INDRAD( INEXT, 1, .true.)
  SELF%ifrted = INDRAD( INEXT, nlwout, .true.)
  SELF%ifrsodc = INDRAD( INEXT, 1, .true.)
  SELF%ifrtedc = INDRAD( INEXT, 1, .true.)
  SELF%iemit = INDRAD( INEXT, 1, .true.)
  SELF%isudu = INDRAD( INEXT, 1, .true.)
  SELF%iuvdf = INDRAD( INEXT, 1, .true.)
  SELF%iparf = INDRAD( INEXT, 1, .true.)
  SELF%iparcf = INDRAD( INEXT, 1, .true.)
  SELF%itincf = INDRAD( INEXT, 1, .true.)
  SELF%ifdir = INDRAD( INEXT, 1, .true.)
  SELF%ifdif = INDRAD( INEXT, 1, .true.)
  SELF%icdir = INDRAD( INEXT, 1, .true.)
  SELF%ilwderivative = INDRAD( INEXT, nlev+1, lapproxlwupdate)
  SELF%iswdirectband = INDRAD( INEXT, nsw, lapproxswupdate)
  SELF%iswdiffuseband = INDRAD( INEXT, nsw, lapproxswupdate)
  SELF%ifrso = INDRAD( INEXT, nlev+1, .true.)
  SELF%iswfc = INDRAD( INEXT, nlev+1, .true.)
  SELF%ifrth = INDRAD( INEXT, nlev+1, .true.)
  SELF%ilwfc = INDRAD( INEXT, nlev+1, .true.)
  SELF%iaer = INDRAD( INEXT, 6*nlev, ldiagforcing)
  SELF%ioz = INDRAD( INEXT, nlev, ldiagforcing)
  SELF%iico2 = INDRAD( INEXT, nlev, ldiagforcing)
  SELF%iich4 = INDRAD( INEXT, nlev, ldiagforcing)
  SELF%iin2o = INDRAD( INEXT, nlev, ldiagforcing)
  SELF%ino2 = INDRAD( INEXT, nlev, ldiagforcing)
  SELF%ic11 = INDRAD( INEXT, nlev, ldiagforcing)
  SELF%ic12 = INDRAD( INEXT, nlev, ldiagforcing)
  SELF%ic22 = INDRAD( INEXT, nlev, ldiagforcing)
  SELF%icl4 = INDRAD( INEXT, nlev, ldiagforcing)
  SELF%igix = INDRAD( INEXT, 1, ldebug)
  SELF%IOUTEND = INEXT-1
  SELF%iaer = INDRAD( INEXT, 6*nlev, .not.ldiagforcing)
  SELF%ioz = INDRAD( INEXT, nlev, .not.(ldiagforcing.or.lrayfm))
  SELF%iico2 = INDRAD( INEXT, nlev, .not.ldiagforcing)
  SELF%iich4 = INDRAD( INEXT, nlev, .not.ldiagforcing)
  SELF%iin2o = INDRAD( INEXT, nlev, .not.ldiagforcing)
  SELF%ino2 = INDRAD( INEXT, nlev, .not.ldiagforcing)
  SELF%ic11 = INDRAD( INEXT, nlev, .not.ldiagforcing)
  SELF%ic12 = INDRAD( INEXT, nlev, .not.ldiagforcing)
  SELF%ic22 = INDRAD( INEXT, nlev, .not.ldiagforcing)
  SELF%icl4 = INDRAD( INEXT, nlev, .not.ldiagforcing)

  IF (LDEBUG) THEN
    WRITE(NULOUT,'("RADINTG: ",A16,"=",I12)') 'IGI',SELF%igi
    WRITE(NULOUT,'("RADINTG: ",A16,"=",I12)') 'IMU0',SELF%imu0
    WRITE(NULOUT,'("RADINTG: ",A16,"=",I12)') 'IAMU0',SELF%iamu0
    WRITE(NULOUT,'("RADINTG: ",A16,"=",I12)') 'IEMISS',SELF%iemiss
    WRITE(NULOUT,'("RADINTG: ",A16,"=",I12)') 'ITS',SELF%its
    WRITE(NULOUT,'("RADINTG: ",A16,"=",I12)') 'ISLM',SELF%islm
    WRITE(NULOUT,'("RADINTG: ",A16,"=",I12)') 'ICCNL',SELF%iccnl
    WRITE(NULOUT,'("RADINTG: ",A16,"=",I12)') 'ICCNO',SELF%iccno
    WRITE(NULOUT,'("RADINTG: ",A16,"=",I12)') 'IBAS',SELF%ibas
    WRITE(NULOUT,'("RADINTG: ",A16,"=",I12)') 'ITOP',SELF%itop
    WRITE(NULOUT,'("RADINTG: ",A16,"=",I12)') 'IGELAM',SELF%igelam
    WRITE(NULOUT,'("RADINTG: ",A16,"=",I12)') 'IGEMU',SELF%igemu
    WRITE(NULOUT,'("RADINTG: ",A16,"=",I12)') 'ICLON',SELF%iclon
    WRITE(NULOUT,'("RADINTG: ",A16,"=",I12)') 'ISLON',SELF%islon
    WRITE(NULOUT,'("RADINTG: ",A16,"=",I12)') 'IALD',SELF%iald
    WRITE(NULOUT,'("RADINTG: ",A16,"=",I12)') 'IALP',SELF%ialp
    WRITE(NULOUT,'("RADINTG: ",A16,"=",I12)') 'ITI',SELF%iti
    WRITE(NULOUT,'("RADINTG: ",A16,"=",I12)') 'IPR',SELF%ipr
    WRITE(NULOUT,'("RADINTG: ",A16,"=",I12)') 'IQS',SELF%iqs
    WRITE(NULOUT,'("RADINTG: ",A16,"=",I12)') 'IWV',SELF%iwv
    WRITE(NULOUT,'("RADINTG: ",A16,"=",I12)') 'ICLC',SELF%iclc
    WRITE(NULOUT,'("RADINTG: ",A16,"=",I12)') 'ILWA',SELF%ilwa
    WRITE(NULOUT,'("RADINTG: ",A16,"=",I12)') 'IIWA',SELF%iiwa
    WRITE(NULOUT,'("RADINTG: ",A16,"=",I12)') 'ISWA',SELF%iswa
    WRITE(NULOUT,'("RADINTG: ",A16,"=",I12)') 'IRWA',SELF%irwa
    WRITE(NULOUT,'("RADINTG: ",A16,"=",I12)') 'IRRA',SELF%irra
    WRITE(NULOUT,'("RADINTG: ",A16,"=",I12)') 'IDP',SELF%idp
    WRITE(NULOUT,'("RADINTG: ",A16,"=",I12)') 'IOZ',SELF%ioz
    WRITE(NULOUT,'("RADINTG: ",A16,"=",I12)') 'IECPO3',SELF%iecpo3
    WRITE(NULOUT,'("RADINTG: ",A16,"=",I12)') 'IHPR',SELF%ihpr
    WRITE(NULOUT,'("RADINTG: ",A16,"=",I12)') 'IAPRS',SELF%iaprs
    WRITE(NULOUT,'("RADINTG: ",A16,"=",I12)') 'IHTI',SELF%ihti
    WRITE(NULOUT,'("RADINTG: ",A16,"=",I12)') 'IRE_LIQ',SELF%ire_liq
    WRITE(NULOUT,'("RADINTG: ",A16,"=",I12)') 'IRE_ICE',SELF%ire_ice
    WRITE(NULOUT,'("RADINTG: ",A16,"=",I12)') 'IOVERLAP',SELF%ioverlap
    WRITE(NULOUT,'("RADINTG: ",A16,"=",I12)') 'IAERO',SELF%iaero
    WRITE(NULOUT,'("RADINTG: ",A16,"=",I12)') 'IFRSOD',SELF%ifrsod
    WRITE(NULOUT,'("RADINTG: ",A16,"=",I12)') 'IFRTED',SELF%ifrted
    WRITE(NULOUT,'("RADINTG: ",A16,"=",I12)') 'IFRSODC',SELF%ifrsodc
    WRITE(NULOUT,'("RADINTG: ",A16,"=",I12)') 'IFRTEDC',SELF%ifrtedc
    WRITE(NULOUT,'("RADINTG: ",A16,"=",I12)') 'IEMIT',SELF%iemit
    WRITE(NULOUT,'("RADINTG: ",A16,"=",I12)') 'ISUDU',SELF%isudu
    WRITE(NULOUT,'("RADINTG: ",A16,"=",I12)') 'IUVDF',SELF%iuvdf
    WRITE(NULOUT,'("RADINTG: ",A16,"=",I12)') 'IPARF',SELF%iparf
    WRITE(NULOUT,'("RADINTG: ",A16,"=",I12)') 'IPARCF',SELF%iparcf
    WRITE(NULOUT,'("RADINTG: ",A16,"=",I12)') 'ITINCF',SELF%itincf
    WRITE(NULOUT,'("RADINTG: ",A16,"=",I12)') 'IFDIR',SELF%ifdir
    WRITE(NULOUT,'("RADINTG: ",A16,"=",I12)') 'IFDIF',SELF%ifdif
    WRITE(NULOUT,'("RADINTG: ",A16,"=",I12)') 'ICDIR',SELF%icdir
    WRITE(NULOUT,'("RADINTG: ",A16,"=",I12)') 'ILWDERIVATIVE',SELF%ilwderivative
    WRITE(NULOUT,'("RADINTG: ",A16,"=",I12)') 'ISWDIRECTBAND',SELF%iswdirectband
    WRITE(NULOUT,'("RADINTG: ",A16,"=",I12)') 'ISWDIFFUSEBAND',SELF%iswdiffuseband
    WRITE(NULOUT,'("RADINTG: ",A16,"=",I12)') 'IFRSO',SELF%ifrso
    WRITE(NULOUT,'("RADINTG: ",A16,"=",I12)') 'ISWFC',SELF%iswfc
    WRITE(NULOUT,'("RADINTG: ",A16,"=",I12)') 'IFRTH',SELF%ifrth
    WRITE(NULOUT,'("RADINTG: ",A16,"=",I12)') 'ILWFC',SELF%ilwfc
    WRITE(NULOUT,'("RADINTG: ",A16,"=",I12)') 'IAER',SELF%iaer
    WRITE(NULOUT,'("RADINTG: ",A16,"=",I12)') 'IICO2',SELF%iico2
    WRITE(NULOUT,'("RADINTG: ",A16,"=",I12)') 'IICH4',SELF%iich4
    WRITE(NULOUT,'("RADINTG: ",A16,"=",I12)') 'IIN2O',SELF%iin2o
    WRITE(NULOUT,'("RADINTG: ",A16,"=",I12)') 'INO2',SELF%ino2
    WRITE(NULOUT,'("RADINTG: ",A16,"=",I12)') 'IC11',SELF%ic11
    WRITE(NULOUT,'("RADINTG: ",A16,"=",I12)') 'IC12',SELF%ic12
    WRITE(NULOUT,'("RADINTG: ",A16,"=",I12)') 'IC22',SELF%ic22
    WRITE(NULOUT,'("RADINTG: ",A16,"=",I12)') 'ICL4',SELF%icl4
    WRITE(NULOUT,'("RADINTG: ",A16,"=",I12)') 'IGIX',SELF%igix
  ENDIF

  SELF%IFLDSIN = SELF%IINEND - SELF%IINBEG + 1
  SELF%IFLDSOUT = SELF%IOUTEND - SELF%IOUTBEG + 1
  SELF%IFLDSTOT = INEXT - 1

  WRITE(NULOUT,'("RADINTG: IFLDSIN   =",I12)')SELF%IFLDSIN
  WRITE(NULOUT,'("RADINTG: IFLDSOUT  =",I12)')SELF%IFLDSOUT
  WRITE(NULOUT,'("RADINTG: IFLDSTOT  =",I12)')SELF%IFLDSTOT

  IF (LHOOK) CALL DR_HOOK('RADINTG_ZRGP_SETUP',1,ZHOOK_HANDLE)

END SUBROUTINE RADINTG_ZRGP_SETUP

END MODULE RADINTG_ZRGP_MOD
