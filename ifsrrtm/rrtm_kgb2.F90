SUBROUTINE RRTM_KGB2

!     Originally by Eli J. Mlawer, Atmospheric & Environmental Research.
!     BAND 2:  250-500 cm-1 (low - H2O; high - H2O)
!     Reformatted for F90 by JJMorcrette, ECMWF
!     R. Elkhatib 12-10-2005 Split for faster and more robust compilation.
!     G.Mozdzynski March 2011 read constants from files
!     ABozzo May 2013 update to RRTMG v4.85
!     band 2:  350-500 cm-1 
!     T. Wilhelmsson and K. Yessad (Oct 2013) Geometry and setup refactoring.
!     ------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE YOMLUN    ,ONLY : NULRAD
USE MPL_MODULE,ONLY : MPL_BROADCAST
USE YOMTAG    ,ONLY : MTAGRAD

USE YOERRTO2 , ONLY : KAO     ,KBO     ,KAO_D,KBO_D,SELFREFO   ,FRACREFAO  ,&
 & FRACREFBO  ,FORREFO  
USE YOMMP0    , ONLY : NPROC, MYPROC

!     ------------------------------------------------------------------

IMPLICIT NONE
REAL(KIND=JPRB) :: ZHOOK_HANDLE

#include "abor1.intfb.h"

IF (LHOOK) CALL DR_HOOK('RRTM_KGB2',0,ZHOOK_HANDLE)

IF( MYPROC==1 )THEN
  READ(NULRAD,ERR=1001) KAO_D,KBO_D
 ! Convert the data into model actual precision.
  KAO = REAL(KAO_D,JPRB)
  KBO = REAL(KBO_D,JPRB) 
ENDIF
IF( NPROC>1 )THEN
  CALL MPL_BROADCAST (KAO,MTAGRAD,1,CDSTRING='RRTM_KGB2:')
  CALL MPL_BROADCAST (KBO,MTAGRAD,1,CDSTRING='RRTM_KGB2:')
ENDIF


! Planck fraction mapping level: P = 1053.630 mbar, T = 294.2 K
      FRACREFAO(:) = (/ &
      &  1.6388e-01_JPRB, 1.5241e-01_JPRB, 1.4290e-01_JPRB, 1.2864e-01_JPRB, &
      &  1.1615e-01_JPRB, 1.0047e-01_JPRB, 8.0013e-02_JPRB, 6.0445e-02_JPRB, &
      &  4.0530e-02_JPRB, 4.3879e-03_JPRB, 3.5726e-03_JPRB, 2.7669e-03_JPRB, &
      &  2.0078e-03_JPRB, 1.2864e-03_JPRB, 4.7630e-04_JPRB, 6.9109e-05_JPRB/)

! Planck fraction mapping level: P = 3.206e-2 mb, T = 197.92 K
      FRACREFBO(:) = (/ &
      &  1.4697e-01_JPRB, 1.4826e-01_JPRB, 1.4278e-01_JPRB, 1.3320e-01_JPRB, &
      &  1.1965e-01_JPRB, 1.0297e-01_JPRB, 8.4170e-02_JPRB, 6.3282e-02_JPRB, &
      &  4.2868e-02_JPRB, 4.6644e-03_JPRB, 3.8619e-03_JPRB, 3.0533e-03_JPRB, &
      &  2.2359e-03_JPRB, 1.4226e-03_JPRB, 5.3642e-04_JPRB, 7.6316e-05_JPRB/)

!     The array FORREFO contains the coefficient of the water vapor
!     foreign-continuum (including the energy term).  The first 
!     index refers to reference temperature (296,260,224,260) and 
!     pressure (970,475,219,3 mbar) levels.  The second index 
!     runs over the g-channel (1 to 16).

      FORREFO(1,:) = (/ &
     & 2.8549e-03_JPRB,4.8281e-03_JPRB,6.2570e-03_JPRB,8.2731e-03_JPRB,7.9056e-03_JPRB,7.7840e-03_JPRB, &
     & 1.0115e-02_JPRB,9.6599e-03_JPRB,1.0153e-02_JPRB,1.0921e-02_JPRB,1.2408e-02_JPRB,1.3496e-02_JPRB, &
     & 1.5059e-02_JPRB,1.4636e-02_JPRB,1.6483e-02_JPRB,1.2394e-02_JPRB/)
      FORREFO(2,:) = (/ &
     & 3.0036e-03_JPRB,5.1093e-03_JPRB,5.7317e-03_JPRB,9.2246e-03_JPRB,8.9829e-03_JPRB,8.6477e-03_JPRB, &
     & 1.1448e-02_JPRB,1.0391e-02_JPRB,1.0211e-02_JPRB,1.2921e-02_JPRB,1.2726e-02_JPRB,1.2426e-02_JPRB, &
     & 1.4609e-02_JPRB,1.5783e-02_JPRB,1.6617e-02_JPRB,1.6858e-02_JPRB/)
      FORREFO(3,:) = (/ &
     & 3.0771e-03_JPRB,5.1206e-03_JPRB,5.8426e-03_JPRB,9.5727e-03_JPRB,1.0338e-02_JPRB,9.3737e-03_JPRB, &
     & 1.2805e-02_JPRB,1.1272e-02_JPRB,1.1353e-02_JPRB,1.1837e-02_JPRB,1.1550e-02_JPRB,1.3020e-02_JPRB, &
     & 1.3536e-02_JPRB,1.6226e-02_JPRB,1.6039e-02_JPRB,2.2578e-02_JPRB/)
      FORREFO(4,:) = (/ &
     & 3.3072e-03_JPRB,5.0240e-03_JPRB,6.8474e-03_JPRB,8.2736e-03_JPRB,8.6151e-03_JPRB,8.6762e-03_JPRB, &
     & 1.1476e-02_JPRB,1.0246e-02_JPRB,1.0819e-02_JPRB,1.0640e-02_JPRB,1.0545e-02_JPRB,1.0533e-02_JPRB, &
     & 1.0496e-02_JPRB,1.0142e-02_JPRB,9.7979e-03_JPRB,1.5255e-02_JPRB/)


!     The following are parameters related to the reference water vapor
!     mixing ratios by REFPARAM(I) = REFH2O(I) / (.002+REFH2O(I)).
!     These parameters are used for the Planck function interpolation.
!REFPARAM( :) = (/&
! & 0.903661_JPRB   , 0.859386_JPRB   , 0.746542_JPRB   , 0.580496_JPRB   , 0.412889_JPRB   ,&
! & 0.275283_JPRB   , 0.162745_JPRB   , 7.63929E-02_JPRB, 1.82553E-02_JPRB, 3.72432E-03_JPRB, &
! & 2.14946E-03_JPRB, 1.66320E-03_JPRB, 1.59940E-03_JPRB/)  

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



!     The array SELFREFO contains the coefficient of the water vapor
!     self-continuum (including the energy term).  The first index
!     refers to temperature in 7.2 degree increments.  For instance,
!     JT = 1 refers to a temperature of 245.6, JT = 2 refers to 252.8,
!     etc.  The second index runs over the g-channel (1 to 16).

    SELFREFO(:, 1) = (/ &
     & 7.25695e-01_JPRB, 6.53591e-01_JPRB, 5.88650e-01_JPRB, 5.30162e-01_JPRB, 4.77485e-01_JPRB, &
     & 4.30042e-01_JPRB, 3.87313e-01_JPRB, 3.48830e-01_JPRB, 3.14170e-01_JPRB, 2.82954e-01_JPRB/)
      SELFREFO(:, 2) = (/ &
     & 9.61996e-01_JPRB, 8.77853e-01_JPRB, 8.01070e-01_JPRB, 7.31003e-01_JPRB, 6.67064e-01_JPRB, &
     & 6.08718e-01_JPRB, 5.55476e-01_JPRB, 5.06890e-01_JPRB, 4.62554e-01_JPRB, 4.22096e-01_JPRB/)
      SELFREFO(:, 3) = (/ &
     & 9.72584e-01_JPRB, 9.02658e-01_JPRB, 8.37760e-01_JPRB, 7.77527e-01_JPRB, 7.21626e-01_JPRB, &
     & 6.69743e-01_JPRB, 6.21591e-01_JPRB, 5.76900e-01_JPRB, 5.35423e-01_JPRB, 4.96927e-01_JPRB/)
      SELFREFO(:, 4) = (/ &
     & 1.24790e+00_JPRB, 1.14353e+00_JPRB, 1.04790e+00_JPRB, 9.60263e-01_JPRB, 8.79956e-01_JPRB, &
     & 8.06364e-01_JPRB, 7.38927e-01_JPRB, 6.77130e-01_JPRB, 6.20501e-01_JPRB, 5.68608e-01_JPRB/)
      SELFREFO(:, 5) = (/ &
     & 1.23574e+00_JPRB, 1.12928e+00_JPRB, 1.03200e+00_JPRB, 9.43096e-01_JPRB, 8.61851e-01_JPRB, &
     & 7.87605e-01_JPRB, 7.19755e-01_JPRB, 6.57750e-01_JPRB, 6.01087e-01_JPRB, 5.49305e-01_JPRB/)
      SELFREFO(:, 6) = (/ &
     & 1.20921e+00_JPRB, 1.10660e+00_JPRB, 1.01270e+00_JPRB, 9.26766e-01_JPRB, 8.48124e-01_JPRB, &
     & 7.76155e-01_JPRB, 7.10293e-01_JPRB, 6.50020e-01_JPRB, 5.94861e-01_JPRB, 5.44384e-01_JPRB/)
      SELFREFO(:, 7) = (/ &
     & 1.38112e+00_JPRB, 1.26727e+00_JPRB, 1.16280e+00_JPRB, 1.06694e+00_JPRB, 9.78990e-01_JPRB, &
     & 8.98287e-01_JPRB, 8.24236e-01_JPRB, 7.56290e-01_JPRB, 6.93945e-01_JPRB, 6.36739e-01_JPRB/)
      SELFREFO(:, 8) = (/ &
     & 1.30321e+00_JPRB, 1.20127e+00_JPRB, 1.10730e+00_JPRB, 1.02068e+00_JPRB, 9.40840e-01_JPRB, &
     & 8.67243e-01_JPRB, 7.99403e-01_JPRB, 7.36870e-01_JPRB, 6.79229e-01_JPRB, 6.26096e-01_JPRB/)
      SELFREFO(:, 9) = (/ &
     & 1.26713e+00_JPRB, 1.17927e+00_JPRB, 1.09750e+00_JPRB, 1.02140e+00_JPRB, 9.50575e-01_JPRB, &
     & 8.84662e-01_JPRB, 8.23319e-01_JPRB, 7.66230e-01_JPRB, 7.13099e-01_JPRB, 6.63653e-01_JPRB/)
      SELFREFO(:,10) = (/ &
     & 1.49824e+00_JPRB, 1.37053e+00_JPRB, 1.25370e+00_JPRB, 1.14683e+00_JPRB, 1.04908e+00_JPRB, &
     & 9.59651e-01_JPRB, 8.77849e-01_JPRB, 8.03020e-01_JPRB, 7.34569e-01_JPRB, 6.71954e-01_JPRB/)
      SELFREFO(:,11) = (/ &
     & 1.44786e+00_JPRB, 1.34594e+00_JPRB, 1.25120e+00_JPRB, 1.16313e+00_JPRB, 1.08125e+00_JPRB, &
     & 1.00514e+00_JPRB, 9.34392e-01_JPRB, 8.68620e-01_JPRB, 8.07477e-01_JPRB, 7.50639e-01_JPRB/)
      SELFREFO(:,12) = (/ &
     & 1.38460e+00_JPRB, 1.30437e+00_JPRB, 1.22880e+00_JPRB, 1.15760e+00_JPRB, 1.09053e+00_JPRB, &
     & 1.02735e+00_JPRB, 9.67825e-01_JPRB, 9.11750e-01_JPRB, 8.58924e-01_JPRB, 8.09159e-01_JPRB/)
      SELFREFO(:,13) = (/ &
     & 1.51953e+00_JPRB, 1.42822e+00_JPRB, 1.34240e+00_JPRB, 1.26173e+00_JPRB, 1.18592e+00_JPRB, &
     & 1.11465e+00_JPRB, 1.04768e+00_JPRB, 9.84720e-01_JPRB, 9.25548e-01_JPRB, 8.69932e-01_JPRB/)
      SELFREFO(:,14) = (/ &
     & 1.62608e+00_JPRB, 1.51021e+00_JPRB, 1.40260e+00_JPRB, 1.30266e+00_JPRB, 1.20983e+00_JPRB, &
     & 1.12363e+00_JPRB, 1.04356e+00_JPRB, 9.69200e-01_JPRB, 9.00138e-01_JPRB, 8.35998e-01_JPRB/)
      SELFREFO(:,15) = (/ &
     & 1.65383e+00_JPRB, 1.54808e+00_JPRB, 1.44910e+00_JPRB, 1.35644e+00_JPRB, 1.26971e+00_JPRB, &
     & 1.18853e+00_JPRB, 1.11254e+00_JPRB, 1.04140e+00_JPRB, 9.74813e-01_JPRB, 9.12484e-01_JPRB/)
      SELFREFO(:,16) = (/ &
     & 1.78105e+00_JPRB, 1.61421e+00_JPRB, 1.46300e+00_JPRB, 1.32595e+00_JPRB, 1.20174e+00_JPRB, &
     & 1.08917e+00_JPRB, 9.87141e-01_JPRB, 8.94670e-01_JPRB, 8.10861e-01_JPRB, 7.34904e-01_JPRB/)


IF (LHOOK) CALL DR_HOOK('RRTM_KGB2',1,ZHOOK_HANDLE)
RETURN

1001 CONTINUE
CALL ABOR1("RRTM_KGB2:ERROR READING FILE RADRRTM")

END SUBROUTINE RRTM_KGB2
