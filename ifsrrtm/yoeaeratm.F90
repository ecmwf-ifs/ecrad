MODULE YOEAERATM

USE PARKIND1  ,ONLY : JPIM     ,JPRB

IMPLICIT NONE

SAVE

!     ------------------------------------------------------------------
!*    ** *YOEAERATM* - CONTROL PARAMETERS FOR AEROSOLS IN THE ATMOSPHERE
!     ------------------------------------------------------------------

INTEGER(KIND=JPIM) :: NAERCONF
INTEGER(KIND=JPIM) :: NINIDAY   
INTEGER(KIND=JPIM) :: NXT3DAER
INTEGER(KIND=JPIM) :: NDD1, NSS1
INTEGER(KIND=JPIM) :: NBCOPTP, NDDOPTP, NOMOPTP, NSSOPTP, NSUOPTP
INTEGER(KIND=JPIM) :: NVISWL
INTEGER(KIND=JPIM) :: NDRYDEP
INTEGER(KIND=JPIM) :: NWHTDDEP, NTDDEP, NINDDDEP(15)
INTEGER(KIND=JPIM) :: NWHTSCAV, NTSCAV, NINDSCAV(15)
INTEGER(KIND=JPIM) :: NWHTSEDM, NTSEDM, NINDSEDM(15)
INTEGER(KIND=JPIM) :: NWHTWDEP, NTWDEP, NINDWDEP(15)

REAL(KIND=JPRB) :: RGRATE
REAL(KIND=JPRB) :: RMASSE(15)

REAL(KIND=JPRB) :: REPSCAER

LOGICAL :: LAERCLIMG, LAERCLIMZ, LAERCLIST, LAERDRYDP, LAERHYGRO, LAERLISI 
LOGICAL :: LAERNGAT , LAERSCAV , LAERSEDIM, LAERSURF , LAERELVS , LAER6SDIA
LOGICAL :: LAERGTOP , LAERRAD  , LAERCCN  , LAEROPT(8),LAERINIT , LAERVOL
LOGICAL :: LAERCSTR , LAERDIAG1, LAERDIAG2, LAERRRTM , LAERUVP  , LUVINDX   
LOGICAL :: LAEREXTR , LAERGBUD , LAERPRNT , LAERCALIP

!     ------------------------------------------------------------------
! NDD1       : location of first bin for desert dust
! NSS1       : location of first bin for sea-salt
! NBCOPTP    : index for choosing the black carbon LW and SW optic.prop. (1: Boucher, 2: Bond, Bergstrom, 3: Stier et al.)
! NDDOPTP    : index for choosing the dust LW and SW optic.prop. (1: Boucher, 2: Highwood, 3: Woodward)
! NOMOPTP    : index for choosing the organic carbon optic.prop.
! NSSOPTP    : index for choosing the sea salt optic.prop.
! NSUOPTP    : index for choosing the sulphate optic.prop.
! NVISWL     : index of wavelength for visibility computations
! RMFMIN     : minimum mass flux for convective aerosol transport
! RGRATE     : transformation rate from hygrophopic to hygrophilic for BC and OM aerosols
! RMASSE     : Molar mass: N.B.: either g/mol or Avogadro number

! REPSCAER   : security on aerosol concentration: always >= 1.E-15

! LAERCLIMG  : .T. to start prognostic aerosols with geographical monthly 
!                  mean climatology
! LAERCLIMZ  : .T. to start prognostic aerosols with zonal annual mean 
!                  climatology
! LAERCLIST  : .T. to start prognostic aerosols with geographical monthly 
!                  mean climatology for background stratospheric only
! LAERDRYDP  : .T. dry deposition is active
! LAERHYDRO  : .T. hygroscopic effects on BC and OM aerosols
! LAERNGAT   : .T. prevents negative aerosol concentrations
! LAERSCAV   : .T. in-cloud and below cloud scavenging is active
! LAERSEDIM  : .T. sedimentation is active
! LAERSURF   : .T. if surface emissions
! LAERELVS   : .T. if "elevated" source
! LAER6SDIA  : .T. if radiance diagnostics with 6S
! LAERGTOP   : .T. if gas-to-particle conversion for SO2/SO4
! LAERRAD    : .T. if there is any prognostic aerosols used for RT
! LAEROPT(.) : .T. if a given aerosol type is radiatively interactive
! LAERCCN    : .T. if prognostic aerosols are used to define the Re of liq.wat.clds
! LAERUVP    : .T. if prognostic aerosols are used in UV-processor
! LAERCSTR   : .T. if climatological stratospheric aerosols are used in radiation 
!                  schemes with the prognostic tropospheric aerosols.
! LAERRRTM   : .T. if RRTM schemes get the information from the prognostic aerosols
! LAERINIT   : .T. if analysed prognostic aerosols are ONLY used as "climatological" aerosols
! LAERVOL    : .T. if volcanic aerosol is considered
! NDRYDEP    : dry deposition 1: a la GEMS; 2: "exp" (a la D_GRG_4.6)
! NWHTDDEP   : =0 only SS+DU;  =1 SS+DU+VOL;  =2 SS+DU+OM+BC+SU+VOL
! NTDDEP     : total number of aerosol species to be in- and below-cloud scavenged
! NINDDDEP   : actual set of aerosol indices to be dry deposited
! NWHTSEDM   : =0 only SS+DU;  =1 SS+DU+VOL;  =2 SS+DU+OM+BC+SU+VOL
! NTSEDM     : total number of aerosol species to be in- and below-cloud scavenged
! NINDSEDM   : actual set of aerosol indices to be sedimented
! NWHTWDEP   : =0 only SS+DU;  =1 SS+DU+VOL;  =2 SS+DU+OM+BC+SU+VOL
! NTWDEP     : total number of aerosol species to be in- and below-cloud scavenged
! NINDWDEP   : actual set of aerosol indices to be scavenged
! NWHTSCAV   : =0 only SS+DU;  =1 SS+DU+VOL;  =2 SS+DU+OM+BC+SU+VOL
! NTSCAV     : total number of aerosol species to be in- and below-cloud scavenged
! NINDSCAV   : actual set of aerosol indices to be scavenged
! LAERCALIP  : .T. works with LAERLISI=t, store only the CALIOP-type profile at 532 nm
!     ------------------------------------------------------------------
END MODULE YOEAERATM

