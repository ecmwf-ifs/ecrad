#if defined ( VPP )
      FUNCTION TIMEF()
      DATA IFIRST/0/
      SAVE IFIRST,ZFIRST
!
!   fujitsu special
!
      CALL GETTOD(ZTIMEF)
      ZTIMEF=ZTIMEF/1000.0
!
      IF(IFIRST.EQ.0) THEN
         IFIRST=1
         ZFIRST=ZTIMEF
         TIMEF=0.0
      ELSE
         TIMEF=(ZTIMEF-ZFIRST)
      ENDIF
!
      RETURN
      ENDFUNCTION TIMEF
#elif defined ( RS6K )
      REAL*8 FUNCTION TIMEF()
      INTEGER*8 IRTC
      EXTERNAL IRTC
      DATA IFIRST/0/
      SAVE IFIRST,ZFIRST
!
      ZTIMEF=IRTC()*1.0E-6
!
      IF(IFIRST.EQ.0) THEN
         IFIRST=1
         ZFIRST=ZTIMEF
         TIMEF=0.0
      ELSE
         TIMEF=(ZTIMEF-ZFIRST)
      ENDIF
!
      RETURN
      ENDFUNCTION TIMEF
#elif defined ( _CRAYFTN )
      REAL*8 FUNCTION TIMEF()
      INTEGER*8 I_TIME1,I_RATE
      REAL*8 ZTIMEF,ZFIRST
      DATA IFIRST/0/
      SAVE IFIRST,ZFIRST

C     THIS IS A REAL TIME CLOCK IN MILLISECONDS.
C     USAGE OF FORTRAN90 INTRINSIC FUNCTION.
      CALL SYSTEM_CLOCK (COUNT=I_TIME1,COUNT_RATE=I_RATE)
      ZTIMEF=DFLOAT(I_TIME1)/DFLOAT(I_RATE)

      IF(IFIRST.EQ.0) THEN
         IFIRST=1
         ZFIRST=ZTIMEF
         TIMEF=0.0
      ELSE
         TIMEF=1000.0*(ZTIMEF-ZFIRST)
      ENDIF

      RETURN
      ENDFUNCTION TIMEF

#else

      FUNCTION TIMEF()
!ss: Uses util_walltime_() from DrHook -- 
!    Should now have higher precision (avoiding negative accum. wall-time in ifs.stat (see opdis.F90))
      USE PARKIND1, ONLY : JPIM, JPRD
      IMPLICIT NONE
      REAL(KIND=JPRD)           :: TIMEF
      REAL(KIND=JPRD)           :: ZTIMEF
      REAL(KIND=JPRD),SAVE      :: ZFIRST
      INTEGER(KIND=JPIM),SAVE   :: IFIRST = 0
      REAL(KIND=JPRD), EXTERNAL :: UTIL_WALLTIME ! ifsaux/support/drhook.c
      ZTIMEF = UTIL_WALLTIME()
      IF(IFIRST.EQ.0) THEN
         IFIRST=1
         ZFIRST=ZTIMEF
         TIMEF=0.0_JPRD
      ELSE
         TIMEF=1000.0_JPRD * (ZTIMEF-ZFIRST)
      ENDIF
      ENDFUNCTION TIMEF

#endif
