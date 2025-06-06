SUBROUTINE RRTM_KGB1(CDIRECTORY)

!     Originally by Eli J. Mlawer, Atmospheric & Environmental Research.
!     BAND 1:  10-250 cm-1 (low - H2O; high - H2O)
!     Reformatted for F90 by JJMorcrette, ECMWF
!     R. Elkhatib 12-10-2005 Split for faster and more robust compilation.
!     G.Mozdzynski March 2011 read constants from files
!     ABozzo May 2013 update to RRTMG v4.85
!     band 1:  10-350 cm-1
!     T. Wilhelmsson and K. Yessad (Oct 2013) Geometry and setup refactoring.
!     F. Vana  05-Mar-2015  Support for single precision
!     ------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK
USE YOMLUN    ,ONLY : NULRAD, NULOUT
USE MPL_MODULE,ONLY : MPL_BROADCAST
USE YOMTAG    ,ONLY : MTAGRAD
USE YOMMP0_IFSAUX    ,ONLY : NPROC, MYPROC

USE YOERRTO1 , ONLY : KAO     ,KBO     ,SELFREFO   ,FRACREFAO ,&
 & FRACREFBO  ,FORREFO, KAO_MN2, KBO_MN2, KAO_D, KBO_D

!     ------------------------------------------------------------------

IMPLICIT NONE

CHARACTER(LEN=*), INTENT(IN) :: CDIRECTORY

CHARACTER(LEN=512) :: CLF1

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

#include "abor1.intfb.h"

IF (LHOOK) CALL DR_HOOK('RRTM_KGB1',0,ZHOOK_HANDLE)

IF( MYPROC==1 )THEN
  CLF1 = TRIM(CDIRECTORY) // "/RADRRTM"
  WRITE(NULOUT,'(a,a)') 'Reading RRTMG longwave data file ', TRIM(CLF1)
  OPEN(NULRAD,FILE=TRIM(CLF1),FORM="UNFORMATTED",ACTION="READ",ERR=1000,CONVERT='BIG_ENDIAN')

  READ(NULRAD,ERR=1001) KAO_D,KBO_D
  ! Convert the data into model actual precision.
  KAO = REAL(KAO_D,JPRB)
  KBO = REAL(KBO_D,JPRB)
ENDIF
IF( NPROC>1 )THEN
  CALL MPL_BROADCAST (KAO,MTAGRAD,1,CDSTRING='RRTM_KGB1:')
  CALL MPL_BROADCAST (KBO,MTAGRAD,1,CDSTRING='RRTM_KGB1:')
ENDIF

! Planck fraction mapping level: P = 212.7250 mbar, T = 223.06 K
FRACREFAO(:) = (/ &
 & 2.1227E-01_JPRB,1.8897E-01_JPRB,1.3934E-01_JPRB,1.1557E-01_JPRB,9.5282E-02_JPRB,8.3359E-02_JPRB, &
 & 6.5333E-02_JPRB,5.2016E-02_JPRB,3.4272E-02_JPRB,4.0257E-03_JPRB,3.1857E-03_JPRB,2.6014E-03_JPRB, &
 & 1.9141E-03_JPRB,1.2612E-03_JPRB,5.3169E-04_JPRB,7.6476E-05_JPRB/)

! Planck fraction mapping level: P = 212.7250 mbar, T = 223.06 K
! These Planck fractions were calculated using lower atmosphere
! parameters.
FRACREFBO(:) = (/ &
 & 2.1227E-01_JPRB,1.8897E-01_JPRB,1.3934E-01_JPRB,1.1557E-01_JPRB,9.5282E-02_JPRB,8.3359E-02_JPRB, &
 & 6.5333E-02_JPRB,5.2016E-02_JPRB,3.4272E-02_JPRB,4.0257E-03_JPRB,3.1857E-03_JPRB,2.6014E-03_JPRB, &
 & 1.9141E-03_JPRB,1.2612E-03_JPRB,5.3169E-04_JPRB,7.6476E-05_JPRB/)

!     The array FORREFO contains the coefficient of the water vapor
!     foreign-continuum (including the energy term).  The first
!     index refers to reference temperature (296,260,224,260) and
!     pressure (970,475,219,3 mbar) levels.  The second index
!     runs over the g-channel (1 to 16).

      FORREFO(1,:) = (/ &
     & 3.6742E-02_JPRB,1.0664E-01_JPRB,2.6132E-01_JPRB,2.7906E-01_JPRB,2.8151E-01_JPRB,2.7465E-01_JPRB, &
     & 2.8530E-01_JPRB,2.9123E-01_JPRB,3.0697E-01_JPRB,3.1801E-01_JPRB,3.2444E-01_JPRB,2.7746E-01_JPRB, &
     & 3.1994E-01_JPRB,2.9750E-01_JPRB,2.1226E-01_JPRB,1.2847E-01_JPRB/)
      FORREFO(2,:) = (/ &
     & 4.0450E-02_JPRB,1.1085E-01_JPRB,2.9205E-01_JPRB,3.1934E-01_JPRB,3.1739E-01_JPRB,3.1450E-01_JPRB, &
     & 3.2797E-01_JPRB,3.2223E-01_JPRB,3.3099E-01_JPRB,3.4800E-01_JPRB,3.4046E-01_JPRB,3.5700E-01_JPRB, &
     & 3.8264E-01_JPRB,3.6679E-01_JPRB,3.3481E-01_JPRB,3.2113E-01_JPRB/)
      FORREFO(3,:) = (/ &
     & 4.6952E-02_JPRB,1.1999E-01_JPRB,3.1473E-01_JPRB,3.7015E-01_JPRB,3.6913E-01_JPRB,3.6352E-01_JPRB, &
     & 3.7754E-01_JPRB,3.7402E-01_JPRB,3.7113E-01_JPRB,3.7720E-01_JPRB,3.8365E-01_JPRB,4.0876E-01_JPRB, &
     & 4.2968E-01_JPRB,4.4186E-01_JPRB,4.3468E-01_JPRB,4.7083E-01_JPRB/)
      FORREFO(4,:) = (/ &
     & 7.0645E-02_JPRB,1.6618E-01_JPRB,2.8516E-01_JPRB,3.1819E-01_JPRB,3.0131E-01_JPRB,2.9552E-01_JPRB, &
     & 2.8972E-01_JPRB,2.9348E-01_JPRB,2.8668E-01_JPRB,2.8483E-01_JPRB,2.8130E-01_JPRB,2.7757E-01_JPRB, &
     & 2.9735E-01_JPRB,3.1684E-01_JPRB,3.0681E-01_JPRB,3.6778E-01_JPRB/)


!     ------------------------------------------------------------------

!     The array KAO contains absorption coefs at the 16 chosen g-values
!     for a range of pressure levels > ~100mb and temperatures.  The first
!     index in the array, JT, which runs from 1 to 5, corresponds to
!     different temperatures.  More specifically, JT = 3 means that the
!     data are for the corresponding TREF for this  pressure level,
!     JT = 2 refers to the temperatureTREF-15, JT = 1 is for TREF-30,
!     JT = 4 is for TREF+15, and JT = 5 is for TREF+30.  The second
!     index, JP, runs from 1 to 13 and refers to the corresponding
!     pressure level in PREF (e.g. JP = 1 is for a pressure of 1053.63 mb).
!     The third index, IG, goes from 1 to 16, and tells us which
!     g-interval the absorption coefficients are for.



!     The array KBO contains absorption coefs at the 16 chosen g-values
!     for a range of pressure levels < ~100mb and temperatures. The first
!     index in the array, JT, which runs from 1 to 5, corresponds to
!     different temperatures.  More specifically, JT = 3 means that the
!     data are for the reference temperature TREF for this pressure
!     level, JT = 2 refers to the temperature TREF-15, JT = 1 is for
!     TREF-30, JT = 4 is for TREF+15, and JT = 5 is for TREF+30.
!     The second index, JP, runs from 13 to 59 and refers to the JPth
!     reference pressure level (see taumol.f for the value of these
!     pressure levels in mb).  The third index, IG, goes from 1 to 16,
!     and tells us which g-interval the absorption coefficients are for.



     KAO_MN2(:, 1) = (/ &
     & 5.12042E-08_JPRB, 5.51239E-08_JPRB, 5.93436E-08_JPRB, 6.38863E-08_JPRB, 6.87767E-08_JPRB, &
     & 7.40415E-08_JPRB, 7.97093E-08_JPRB, 8.58110E-08_JPRB, 9.23797E-08_JPRB, 9.94513E-08_JPRB, &
     & 1.07064E-07_JPRB, 1.15260E-07_JPRB, 1.24083E-07_JPRB, 1.33581E-07_JPRB, 1.43807E-07_JPRB, &
     & 1.54815E-07_JPRB, 1.66666E-07_JPRB, 1.79424E-07_JPRB, 1.93159E-07_JPRB/)
      KAO_MN2(:, 2) = (/ &
     & 2.30938E-07_JPRB, 2.41696E-07_JPRB, 2.52955E-07_JPRB, 2.64738E-07_JPRB, 2.77071E-07_JPRB, &
     & 2.89978E-07_JPRB, 3.03486E-07_JPRB, 3.17623E-07_JPRB, 3.32419E-07_JPRB, 3.47904E-07_JPRB, &
     & 3.64111E-07_JPRB, 3.81072E-07_JPRB, 3.98824E-07_JPRB, 4.17402E-07_JPRB, 4.36846E-07_JPRB, &
     & 4.57196E-07_JPRB, 4.78494E-07_JPRB, 5.00784E-07_JPRB, 5.24112E-07_JPRB/)
      KAO_MN2(:, 3) = (/ &
     & 6.70458E-07_JPRB, 7.04274E-07_JPRB, 7.39795E-07_JPRB, 7.77109E-07_JPRB, 8.16304E-07_JPRB, &
     & 8.57476E-07_JPRB, 9.00724E-07_JPRB, 9.46154E-07_JPRB, 9.93876E-07_JPRB, 1.04400E-06_JPRB, &
     & 1.09666E-06_JPRB, 1.15197E-06_JPRB, 1.21008E-06_JPRB, 1.27111E-06_JPRB, 1.33522E-06_JPRB, &
     & 1.40256E-06_JPRB, 1.47331E-06_JPRB, 1.54761E-06_JPRB, 1.62567E-06_JPRB/)
      KAO_MN2(:, 4) = (/ &
     & 1.84182E-06_JPRB, 1.89203E-06_JPRB, 1.94360E-06_JPRB, 1.99658E-06_JPRB, 2.05101E-06_JPRB, &
     & 2.10692E-06_JPRB, 2.16435E-06_JPRB, 2.22335E-06_JPRB, 2.28396E-06_JPRB, 2.34622E-06_JPRB, &
     & 2.41017E-06_JPRB, 2.47587E-06_JPRB, 2.54337E-06_JPRB, 2.61270E-06_JPRB, 2.68392E-06_JPRB, &
     & 2.75708E-06_JPRB, 2.83224E-06_JPRB, 2.90944E-06_JPRB, 2.98875E-06_JPRB/)
      KAO_MN2(:, 5) = (/ &
     & 3.41996E-06_JPRB, 3.32758E-06_JPRB, 3.23770E-06_JPRB, 3.15024E-06_JPRB, 3.06515E-06_JPRB, &
     & 2.98235E-06_JPRB, 2.90180E-06_JPRB, 2.82341E-06_JPRB, 2.74715E-06_JPRB, 2.67294E-06_JPRB, &
     & 2.60074E-06_JPRB, 2.53049E-06_JPRB, 2.46214E-06_JPRB, 2.39563E-06_JPRB, 2.33092E-06_JPRB, &
     & 2.26796E-06_JPRB, 2.20670E-06_JPRB, 2.14709E-06_JPRB, 2.08910E-06_JPRB/)
      KAO_MN2(:, 6) = (/ &
     & 3.38746E-06_JPRB, 3.25966E-06_JPRB, 3.13669E-06_JPRB, 3.01836E-06_JPRB, 2.90449E-06_JPRB, &
     & 2.79491E-06_JPRB, 2.68947E-06_JPRB, 2.58801E-06_JPRB, 2.49037E-06_JPRB, 2.39642E-06_JPRB, &
     & 2.30601E-06_JPRB, 2.21902E-06_JPRB, 2.13530E-06_JPRB, 2.05475E-06_JPRB, 1.97723E-06_JPRB, &
     & 1.90264E-06_JPRB, 1.83086E-06_JPRB, 1.76179E-06_JPRB, 1.69532E-06_JPRB/)
      KAO_MN2(:, 7) = (/ &
     & 3.17530E-06_JPRB, 3.07196E-06_JPRB, 2.97199E-06_JPRB, 2.87527E-06_JPRB, 2.78170E-06_JPRB, &
     & 2.69118E-06_JPRB, 2.60360E-06_JPRB, 2.51887E-06_JPRB, 2.43690E-06_JPRB, 2.35759E-06_JPRB, &
     & 2.28087E-06_JPRB, 2.20664E-06_JPRB, 2.13483E-06_JPRB, 2.06536E-06_JPRB, 1.99814E-06_JPRB, &
     & 1.93312E-06_JPRB, 1.87021E-06_JPRB, 1.80934E-06_JPRB, 1.75046E-06_JPRB/)
      KAO_MN2(:, 8) = (/ &
     & 2.84701E-06_JPRB, 2.77007E-06_JPRB, 2.69521E-06_JPRB, 2.62237E-06_JPRB, 2.55150E-06_JPRB, &
     & 2.48254E-06_JPRB, 2.41545E-06_JPRB, 2.35017E-06_JPRB, 2.28666E-06_JPRB, 2.22486E-06_JPRB, &
     & 2.16473E-06_JPRB, 2.10623E-06_JPRB, 2.04930E-06_JPRB, 1.99392E-06_JPRB, 1.94003E-06_JPRB, &
     & 1.88760E-06_JPRB, 1.83659E-06_JPRB, 1.78695E-06_JPRB, 1.73866E-06_JPRB/)
      KAO_MN2(:, 9) = (/ &
     & 2.79917E-06_JPRB, 2.73207E-06_JPRB, 2.66658E-06_JPRB, 2.60266E-06_JPRB, 2.54027E-06_JPRB, &
     & 2.47937E-06_JPRB, 2.41994E-06_JPRB, 2.36192E-06_JPRB, 2.30530E-06_JPRB, 2.25004E-06_JPRB, &
     & 2.19610E-06_JPRB, 2.14346E-06_JPRB, 2.09208E-06_JPRB, 2.04193E-06_JPRB, 1.99298E-06_JPRB, &
     & 1.94520E-06_JPRB, 1.89857E-06_JPRB, 1.85306E-06_JPRB, 1.80864E-06_JPRB/)
      KAO_MN2(:,10) = (/ &
     & 2.74910E-06_JPRB, 2.64462E-06_JPRB, 2.54412E-06_JPRB, 2.44743E-06_JPRB, 2.35442E-06_JPRB, &
     & 2.26495E-06_JPRB, 2.17887E-06_JPRB, 2.09606E-06_JPRB, 2.01641E-06_JPRB, 1.93978E-06_JPRB, &
     & 1.86606E-06_JPRB, 1.79514E-06_JPRB, 1.72692E-06_JPRB, 1.66129E-06_JPRB, 1.59815E-06_JPRB, &
     & 1.53742E-06_JPRB, 1.47899E-06_JPRB, 1.42278E-06_JPRB, 1.36871E-06_JPRB/)
      KAO_MN2(:,11) = (/ &
     & 2.63952E-06_JPRB, 2.60263E-06_JPRB, 2.56626E-06_JPRB, 2.53039E-06_JPRB, 2.49503E-06_JPRB, &
     & 2.46016E-06_JPRB, 2.42578E-06_JPRB, 2.39188E-06_JPRB, 2.35845E-06_JPRB, 2.32549E-06_JPRB, &
     & 2.29299E-06_JPRB, 2.26094E-06_JPRB, 2.22934E-06_JPRB, 2.19819E-06_JPRB, 2.16747E-06_JPRB, &
     & 2.13717E-06_JPRB, 2.10731E-06_JPRB, 2.07786E-06_JPRB, 2.04882E-06_JPRB/)
      KAO_MN2(:,12) = (/ &
     & 2.94106E-06_JPRB, 2.82819E-06_JPRB, 2.71966E-06_JPRB, 2.61528E-06_JPRB, 2.51492E-06_JPRB, &
     & 2.41841E-06_JPRB, 2.32560E-06_JPRB, 2.23635E-06_JPRB, 2.15053E-06_JPRB, 2.06800E-06_JPRB, &
     & 1.98863E-06_JPRB, 1.91232E-06_JPRB, 1.83893E-06_JPRB, 1.76836E-06_JPRB, 1.70049E-06_JPRB, &
     & 1.63524E-06_JPRB, 1.57248E-06_JPRB, 1.51214E-06_JPRB, 1.45411E-06_JPRB/)
      KAO_MN2(:,13) = (/ &
     & 2.94607E-06_JPRB, 2.87369E-06_JPRB, 2.80309E-06_JPRB, 2.73422E-06_JPRB, 2.66705E-06_JPRB, &
     & 2.60152E-06_JPRB, 2.53760E-06_JPRB, 2.47526E-06_JPRB, 2.41445E-06_JPRB, 2.35513E-06_JPRB, &
     & 2.29726E-06_JPRB, 2.24082E-06_JPRB, 2.18577E-06_JPRB, 2.13207E-06_JPRB, 2.07969E-06_JPRB, &
     & 2.02859E-06_JPRB, 1.97875E-06_JPRB, 1.93014E-06_JPRB, 1.88272E-06_JPRB/)
      KAO_MN2(:,14) = (/ &
     & 2.58051E-06_JPRB, 2.48749E-06_JPRB, 2.39782E-06_JPRB, 2.31139E-06_JPRB, 2.22807E-06_JPRB, &
     & 2.14775E-06_JPRB, 2.07033E-06_JPRB, 1.99570E-06_JPRB, 1.92376E-06_JPRB, 1.85441E-06_JPRB, &
     & 1.78756E-06_JPRB, 1.72313E-06_JPRB, 1.66101E-06_JPRB, 1.60114E-06_JPRB, 1.54342E-06_JPRB, &
     & 1.48778E-06_JPRB, 1.43415E-06_JPRB, 1.38245E-06_JPRB, 1.33262E-06_JPRB/)
      KAO_MN2(:,15) = (/ &
     & 3.03447E-06_JPRB, 2.88559E-06_JPRB, 2.74401E-06_JPRB, 2.60938E-06_JPRB, 2.48135E-06_JPRB, &
     & 2.35961E-06_JPRB, 2.24384E-06_JPRB, 2.13375E-06_JPRB, 2.02906E-06_JPRB, 1.92951E-06_JPRB, &
     & 1.83484E-06_JPRB, 1.74481E-06_JPRB, 1.65921E-06_JPRB, 1.57780E-06_JPRB, 1.50039E-06_JPRB, &
     & 1.42677E-06_JPRB, 1.35677E-06_JPRB, 1.29020E-06_JPRB, 1.22690E-06_JPRB/)
      KAO_MN2(:,16) = (/ &
     & 1.48655E-06_JPRB, 1.48283E-06_JPRB, 1.47913E-06_JPRB, 1.47543E-06_JPRB, 1.47174E-06_JPRB, &
     & 1.46806E-06_JPRB, 1.46439E-06_JPRB, 1.46072E-06_JPRB, 1.45707E-06_JPRB, 1.45343E-06_JPRB, &
     & 1.44979E-06_JPRB, 1.44617E-06_JPRB, 1.44255E-06_JPRB, 1.43894E-06_JPRB, 1.43534E-06_JPRB, &
     & 1.43176E-06_JPRB, 1.42817E-06_JPRB, 1.42460E-06_JPRB, 1.42104E-06_JPRB/)
      KBO_MN2(:, 1) = (/ &
     & 5.12042E-08_JPRB, 5.51239E-08_JPRB, 5.93436E-08_JPRB, 6.38863E-08_JPRB, 6.87767E-08_JPRB, &
     & 7.40415E-08_JPRB, 7.97093E-08_JPRB, 8.58110E-08_JPRB, 9.23797E-08_JPRB, 9.94513E-08_JPRB, &
     & 1.07064E-07_JPRB, 1.15260E-07_JPRB, 1.24083E-07_JPRB, 1.33581E-07_JPRB, 1.43807E-07_JPRB, &
     & 1.54815E-07_JPRB, 1.66666E-07_JPRB, 1.79424E-07_JPRB, 1.93159E-07_JPRB/)
      KBO_MN2(:, 2) = (/ &
     & 2.30938E-07_JPRB, 2.41696E-07_JPRB, 2.52955E-07_JPRB, 2.64738E-07_JPRB, 2.77071E-07_JPRB, &
     & 2.89978E-07_JPRB, 3.03486E-07_JPRB, 3.17623E-07_JPRB, 3.32419E-07_JPRB, 3.47904E-07_JPRB, &
     & 3.64111E-07_JPRB, 3.81072E-07_JPRB, 3.98824E-07_JPRB, 4.17402E-07_JPRB, 4.36846E-07_JPRB, &
     & 4.57196E-07_JPRB, 4.78494E-07_JPRB, 5.00784E-07_JPRB, 5.24112E-07_JPRB/)
      KBO_MN2(:, 3) = (/ &
     & 6.70458E-07_JPRB, 7.04274E-07_JPRB, 7.39795E-07_JPRB, 7.77109E-07_JPRB, 8.16304E-07_JPRB, &
     & 8.57476E-07_JPRB, 9.00724E-07_JPRB, 9.46154E-07_JPRB, 9.93876E-07_JPRB, 1.04400E-06_JPRB, &
     & 1.09666E-06_JPRB, 1.15197E-06_JPRB, 1.21008E-06_JPRB, 1.27111E-06_JPRB, 1.33522E-06_JPRB, &
     & 1.40256E-06_JPRB, 1.47331E-06_JPRB, 1.54761E-06_JPRB, 1.62567E-06_JPRB/)
      KBO_MN2(:, 4) = (/ &
     & 1.84182E-06_JPRB, 1.89203E-06_JPRB, 1.94360E-06_JPRB, 1.99658E-06_JPRB, 2.05101E-06_JPRB, &
     & 2.10692E-06_JPRB, 2.16435E-06_JPRB, 2.22335E-06_JPRB, 2.28396E-06_JPRB, 2.34622E-06_JPRB, &
     & 2.41017E-06_JPRB, 2.47587E-06_JPRB, 2.54337E-06_JPRB, 2.61270E-06_JPRB, 2.68392E-06_JPRB, &
     & 2.75708E-06_JPRB, 2.83224E-06_JPRB, 2.90944E-06_JPRB, 2.98875E-06_JPRB/)
      KBO_MN2(:, 5) = (/ &
     & 3.41996E-06_JPRB, 3.32758E-06_JPRB, 3.23770E-06_JPRB, 3.15024E-06_JPRB, 3.06515E-06_JPRB, &
     & 2.98235E-06_JPRB, 2.90180E-06_JPRB, 2.82341E-06_JPRB, 2.74715E-06_JPRB, 2.67294E-06_JPRB, &
     & 2.60074E-06_JPRB, 2.53049E-06_JPRB, 2.46214E-06_JPRB, 2.39563E-06_JPRB, 2.33092E-06_JPRB, &
     & 2.26796E-06_JPRB, 2.20670E-06_JPRB, 2.14709E-06_JPRB, 2.08910E-06_JPRB/)
      KBO_MN2(:, 6) = (/ &
     & 3.38746E-06_JPRB, 3.25966E-06_JPRB, 3.13669E-06_JPRB, 3.01836E-06_JPRB, 2.90449E-06_JPRB, &
     & 2.79491E-06_JPRB, 2.68947E-06_JPRB, 2.58801E-06_JPRB, 2.49037E-06_JPRB, 2.39642E-06_JPRB, &
     & 2.30601E-06_JPRB, 2.21902E-06_JPRB, 2.13530E-06_JPRB, 2.05475E-06_JPRB, 1.97723E-06_JPRB, &
     & 1.90264E-06_JPRB, 1.83086E-06_JPRB, 1.76179E-06_JPRB, 1.69532E-06_JPRB/)
      KBO_MN2(:, 7) = (/ &
     & 3.17530E-06_JPRB, 3.07196E-06_JPRB, 2.97199E-06_JPRB, 2.87527E-06_JPRB, 2.78170E-06_JPRB, &
     & 2.69118E-06_JPRB, 2.60360E-06_JPRB, 2.51887E-06_JPRB, 2.43690E-06_JPRB, 2.35759E-06_JPRB, &
     & 2.28087E-06_JPRB, 2.20664E-06_JPRB, 2.13483E-06_JPRB, 2.06536E-06_JPRB, 1.99814E-06_JPRB, &
     & 1.93312E-06_JPRB, 1.87021E-06_JPRB, 1.80934E-06_JPRB, 1.75046E-06_JPRB/)
      KBO_MN2(:, 8) = (/ &
     & 2.84701E-06_JPRB, 2.77007E-06_JPRB, 2.69521E-06_JPRB, 2.62237E-06_JPRB, 2.55150E-06_JPRB, &
     & 2.48254E-06_JPRB, 2.41545E-06_JPRB, 2.35017E-06_JPRB, 2.28666E-06_JPRB, 2.22486E-06_JPRB, &
     & 2.16473E-06_JPRB, 2.10623E-06_JPRB, 2.04930E-06_JPRB, 1.99392E-06_JPRB, 1.94003E-06_JPRB, &
     & 1.88760E-06_JPRB, 1.83659E-06_JPRB, 1.78695E-06_JPRB, 1.73866E-06_JPRB/)
      KBO_MN2(:, 9) = (/ &
     & 2.79917E-06_JPRB, 2.73207E-06_JPRB, 2.66658E-06_JPRB, 2.60266E-06_JPRB, 2.54027E-06_JPRB, &
     & 2.47937E-06_JPRB, 2.41994E-06_JPRB, 2.36192E-06_JPRB, 2.30530E-06_JPRB, 2.25004E-06_JPRB, &
     & 2.19610E-06_JPRB, 2.14346E-06_JPRB, 2.09208E-06_JPRB, 2.04193E-06_JPRB, 1.99298E-06_JPRB, &
     & 1.94520E-06_JPRB, 1.89857E-06_JPRB, 1.85306E-06_JPRB, 1.80864E-06_JPRB/)
      KBO_MN2(:,10) = (/ &
     & 2.74910E-06_JPRB, 2.64462E-06_JPRB, 2.54412E-06_JPRB, 2.44743E-06_JPRB, 2.35442E-06_JPRB, &
     & 2.26495E-06_JPRB, 2.17887E-06_JPRB, 2.09606E-06_JPRB, 2.01641E-06_JPRB, 1.93978E-06_JPRB, &
     & 1.86606E-06_JPRB, 1.79514E-06_JPRB, 1.72692E-06_JPRB, 1.66129E-06_JPRB, 1.59815E-06_JPRB, &
     & 1.53742E-06_JPRB, 1.47899E-06_JPRB, 1.42278E-06_JPRB, 1.36871E-06_JPRB/)
      KBO_MN2(:,11) = (/ &
     & 2.63952E-06_JPRB, 2.60263E-06_JPRB, 2.56626E-06_JPRB, 2.53039E-06_JPRB, 2.49503E-06_JPRB, &
     & 2.46016E-06_JPRB, 2.42578E-06_JPRB, 2.39188E-06_JPRB, 2.35845E-06_JPRB, 2.32549E-06_JPRB, &
     & 2.29299E-06_JPRB, 2.26094E-06_JPRB, 2.22934E-06_JPRB, 2.19819E-06_JPRB, 2.16747E-06_JPRB, &
     & 2.13717E-06_JPRB, 2.10731E-06_JPRB, 2.07786E-06_JPRB, 2.04882E-06_JPRB/)
      KBO_MN2(:,12) = (/ &
     & 2.94106E-06_JPRB, 2.82819E-06_JPRB, 2.71966E-06_JPRB, 2.61528E-06_JPRB, 2.51492E-06_JPRB, &
     & 2.41841E-06_JPRB, 2.32560E-06_JPRB, 2.23635E-06_JPRB, 2.15053E-06_JPRB, 2.06800E-06_JPRB, &
     & 1.98863E-06_JPRB, 1.91232E-06_JPRB, 1.83893E-06_JPRB, 1.76836E-06_JPRB, 1.70049E-06_JPRB, &
     & 1.63524E-06_JPRB, 1.57248E-06_JPRB, 1.51214E-06_JPRB, 1.45411E-06_JPRB/)
      KBO_MN2(:,13) = (/ &
     & 2.94607E-06_JPRB, 2.87369E-06_JPRB, 2.80309E-06_JPRB, 2.73422E-06_JPRB, 2.66705E-06_JPRB, &
     & 2.60152E-06_JPRB, 2.53760E-06_JPRB, 2.47526E-06_JPRB, 2.41445E-06_JPRB, 2.35513E-06_JPRB, &
     & 2.29726E-06_JPRB, 2.24082E-06_JPRB, 2.18577E-06_JPRB, 2.13207E-06_JPRB, 2.07969E-06_JPRB, &
     & 2.02859E-06_JPRB, 1.97875E-06_JPRB, 1.93014E-06_JPRB, 1.88272E-06_JPRB/)
      KBO_MN2(:,14) = (/ &
     & 2.58051E-06_JPRB, 2.48749E-06_JPRB, 2.39782E-06_JPRB, 2.31139E-06_JPRB, 2.22807E-06_JPRB, &
     & 2.14775E-06_JPRB, 2.07033E-06_JPRB, 1.99570E-06_JPRB, 1.92376E-06_JPRB, 1.85441E-06_JPRB, &
     & 1.78756E-06_JPRB, 1.72313E-06_JPRB, 1.66101E-06_JPRB, 1.60114E-06_JPRB, 1.54342E-06_JPRB, &
     & 1.48778E-06_JPRB, 1.43415E-06_JPRB, 1.38245E-06_JPRB, 1.33262E-06_JPRB/)
      KBO_MN2(:,15) = (/ &
     & 3.03447E-06_JPRB, 2.88559E-06_JPRB, 2.74401E-06_JPRB, 2.60938E-06_JPRB, 2.48135E-06_JPRB, &
     & 2.35961E-06_JPRB, 2.24384E-06_JPRB, 2.13375E-06_JPRB, 2.02906E-06_JPRB, 1.92951E-06_JPRB, &
     & 1.83484E-06_JPRB, 1.74481E-06_JPRB, 1.65921E-06_JPRB, 1.57780E-06_JPRB, 1.50039E-06_JPRB, &
     & 1.42677E-06_JPRB, 1.35677E-06_JPRB, 1.29020E-06_JPRB, 1.22690E-06_JPRB/)
      KBO_MN2(:,16) = (/ &
     & 1.48655E-06_JPRB, 1.48283E-06_JPRB, 1.47913E-06_JPRB, 1.47543E-06_JPRB, 1.47174E-06_JPRB, &
     & 1.46806E-06_JPRB, 1.46439E-06_JPRB, 1.46072E-06_JPRB, 1.45707E-06_JPRB, 1.45343E-06_JPRB, &
     & 1.44979E-06_JPRB, 1.44617E-06_JPRB, 1.44255E-06_JPRB, 1.43894E-06_JPRB, 1.43534E-06_JPRB, &
     & 1.43176E-06_JPRB, 1.42817E-06_JPRB, 1.42460E-06_JPRB, 1.42104E-06_JPRB/)


!     The array SELFREFO contains the coefficient of the water vapor
!     self-continuum (including the energy term).  The first index
!     refers to temperature in 7.2 degree increments.  For instance,
!     JT = 1 refers to a temperature of 245.6, JT = 2 refers to 252.8,
!     etc.  The second index runs over the g-channel (1 to 16).

      SELFREFO(:, 1) = (/ &
     & 2.16803E+00_JPRB, 1.98236E+00_JPRB, 1.81260E+00_JPRB, 1.65737E+00_JPRB, 1.51544E+00_JPRB, &
     & 1.38567E+00_JPRB, 1.26700E+00_JPRB, 1.15850E+00_JPRB, 1.05929E+00_JPRB, 9.68576E-01_JPRB/)
      SELFREFO(:, 2) = (/ &
     & 3.70149E+00_JPRB, 3.43145E+00_JPRB, 3.18110E+00_JPRB, 2.94902E+00_JPRB, 2.73387E+00_JPRB, &
     & 2.53441E+00_JPRB, 2.34951E+00_JPRB, 2.17810E+00_JPRB, 2.01919E+00_JPRB, 1.87188E+00_JPRB/)
      SELFREFO(:, 3) = (/ &
     & 6.17433E+00_JPRB, 5.62207E+00_JPRB, 5.11920E+00_JPRB, 4.66131E+00_JPRB, 4.24438E+00_JPRB, &
     & 3.86474E+00_JPRB, 3.51906E+00_JPRB, 3.20430E+00_JPRB, 2.91769E+00_JPRB, 2.65672E+00_JPRB/)
      SELFREFO(:, 4) = (/ &
     & 6.56459E+00_JPRB, 5.94787E+00_JPRB, 5.38910E+00_JPRB, 4.88282E+00_JPRB, 4.42410E+00_JPRB, &
     & 4.00848E+00_JPRB, 3.63190E+00_JPRB, 3.29070E+00_JPRB, 2.98155E+00_JPRB, 2.70145E+00_JPRB/)
      SELFREFO(:, 5) = (/ &
     & 6.49581E+00_JPRB, 5.91114E+00_JPRB, 5.37910E+00_JPRB, 4.89494E+00_JPRB, 4.45436E+00_JPRB, &
     & 4.05344E+00_JPRB, 3.68860E+00_JPRB, 3.35660E+00_JPRB, 3.05448E+00_JPRB, 2.77956E+00_JPRB/)
      SELFREFO(:, 6) = (/ &
     & 6.50189E+00_JPRB, 5.89381E+00_JPRB, 5.34260E+00_JPRB, 4.84294E+00_JPRB, 4.39001E+00_JPRB, &
     & 3.97944E+00_JPRB, 3.60727E+00_JPRB, 3.26990E+00_JPRB, 2.96409E+00_JPRB, 2.68687E+00_JPRB/)
      SELFREFO(:, 7) = (/ &
     & 6.64768E+00_JPRB, 6.01719E+00_JPRB, 5.44650E+00_JPRB, 4.92993E+00_JPRB, 4.46236E+00_JPRB, &
     & 4.03914E+00_JPRB, 3.65605E+00_JPRB, 3.30930E+00_JPRB, 2.99543E+00_JPRB, 2.71134E+00_JPRB/)
      SELFREFO(:, 8) = (/ &
     & 6.43744E+00_JPRB, 5.87166E+00_JPRB, 5.35560E+00_JPRB, 4.88490E+00_JPRB, 4.45557E+00_JPRB, &
     & 4.06397E+00_JPRB, 3.70679E+00_JPRB, 3.38100E+00_JPRB, 3.08384E+00_JPRB, 2.81281E+00_JPRB/)
      SELFREFO(:, 9) = (/ &
     & 6.55466E+00_JPRB, 5.99777E+00_JPRB, 5.48820E+00_JPRB, 5.02192E+00_JPRB, 4.59525E+00_JPRB, &
     & 4.20484E+00_JPRB, 3.84759E+00_JPRB, 3.52070E+00_JPRB, 3.22158E+00_JPRB, 2.94787E+00_JPRB/)
      SELFREFO(:,10) = (/ &
     & 6.84510E+00_JPRB, 6.26933E+00_JPRB, 5.74200E+00_JPRB, 5.25902E+00_JPRB, 4.81667E+00_JPRB, &
     & 4.41152E+00_JPRB, 4.04046E+00_JPRB, 3.70060E+00_JPRB, 3.38933E+00_JPRB, 3.10424E+00_JPRB/)
      SELFREFO(:,11) = (/ &
     & 6.83128E+00_JPRB, 6.25536E+00_JPRB, 5.72800E+00_JPRB, 5.24510E+00_JPRB, 4.80291E+00_JPRB, &
     & 4.39799E+00_JPRB, 4.02722E+00_JPRB, 3.68770E+00_JPRB, 3.37681E+00_JPRB, 3.09212E+00_JPRB/)
      SELFREFO(:,12) = (/ &
     & 7.35969E+00_JPRB, 6.61719E+00_JPRB, 5.94960E+00_JPRB, 5.34936E+00_JPRB, 4.80968E+00_JPRB, &
     & 4.32445E+00_JPRB, 3.88817E+00_JPRB, 3.49590E+00_JPRB, 3.14321E+00_JPRB, 2.82610E+00_JPRB/)
      SELFREFO(:,13) = (/ &
     & 7.50064E+00_JPRB, 6.80749E+00_JPRB, 6.17840E+00_JPRB, 5.60744E+00_JPRB, 5.08925E+00_JPRB, &
     & 4.61894E+00_JPRB, 4.19210E+00_JPRB, 3.80470E+00_JPRB, 3.45310E+00_JPRB, 3.13399E+00_JPRB/)
      SELFREFO(:,14) = (/ &
     & 7.40801E+00_JPRB, 6.71328E+00_JPRB, 6.08370E+00_JPRB, 5.51316E+00_JPRB, 4.99613E+00_JPRB, &
     & 4.52759E+00_JPRB, 4.10298E+00_JPRB, 3.71820E+00_JPRB, 3.36950E+00_JPRB, 3.05351E+00_JPRB/)
      SELFREFO(:,15) = (/ &
     & 7.51895E+00_JPRB, 6.68846E+00_JPRB, 5.94970E+00_JPRB, 5.29254E+00_JPRB, 4.70796E+00_JPRB, &
     & 4.18795E+00_JPRB, 3.72538E+00_JPRB, 3.31390E+00_JPRB, 2.94787E+00_JPRB, 2.62227E+00_JPRB/)
      SELFREFO(:,16) = (/ &
     & 7.84774E+00_JPRB, 6.80673E+00_JPRB, 5.90380E+00_JPRB, 5.12065E+00_JPRB, 4.44138E+00_JPRB, &
     & 3.85223E+00_JPRB, 3.34122E+00_JPRB, 2.89800E+00_JPRB, 2.51357E+00_JPRB, 2.18014E+00_JPRB/)


IF (LHOOK) CALL DR_HOOK('RRTM_KGB1',1,ZHOOK_HANDLE)
RETURN

1000 CONTINUE
CALL ABOR1("RRTM_KGB1:ERROR OPENING FILE RADRRTM")
1001 CONTINUE
CALL ABOR1("RRTM_KGB1:ERROR READING FILE RADRRTM")

!     -----------------------------------------------------------------
END SUBROUTINE RRTM_KGB1
