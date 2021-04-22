SUBROUTINE SRTM_TAUMOL26 &
 & ( KIDIA   , KFDIA    , KLEV,&
 & P_COLMOL  ,K_LAYTROP,&
 & P_SFLUXZEN, P_TAUG   , P_TAUR    , PRMU0   &
 & )  

!     Written by Eli J. Mlawer, Atmospheric & Environmental Research.

!     BAND 26:  22650-29000 cm-1 (low - nothing; high - nothing)

!      PARAMETER (MG=16, MXLAY=203, NBANDS=14)

! Modifications
!        M.Hamrud      01-Oct-2003 CY28 Cleaning

!     JJMorcrette 2003-02-24 adapted to ECMWF environment
!        D.Salmond  31-Oct-2007 Vector version in the style of RRTM from Meteo France & NEC
!     JJMorcrette 20110610 Flexible configuration for number of g-points

USE PARKIND1 , ONLY : JPIM, JPRB
USE YOMHOOK  , ONLY : LHOOK, DR_HOOK
USE PARSRTM  , ONLY : JPG
USE YOESRTM  , ONLY : NG26
USE YOESRTA26, ONLY : SFLUXREFC, RAYLC

IMPLICIT NONE

!-- Output
INTEGER(KIND=JPIM),INTENT(IN)    :: KIDIA, KFDIA 
INTEGER(KIND=JPIM),INTENT(IN)    :: KLEV 
REAL(KIND=JPRB)   ,INTENT(IN)    :: P_COLMOL(KIDIA:KFDIA,KLEV) 
INTEGER(KIND=JPIM),INTENT(IN)    :: K_LAYTROP(KIDIA:KFDIA) 

REAL(KIND=JPRB)   ,INTENT(OUT)   :: P_SFLUXZEN(KIDIA:KFDIA,JPG) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: P_TAUG(KIDIA:KFDIA,KLEV,JPG) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: P_TAUR(KIDIA:KFDIA,KLEV,JPG) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PRMU0(KIDIA:KFDIA)
!- from AER
!- from INTFAC      
!- from INTIND
!- from PRECISE             
!- from PROFDATA             
!- from SELF             
INTEGER(KIND=JPIM) :: IG, I_LAY, I_LAYSOLFR(KIDIA:KFDIA), I_NLAYERS, IPLON
INTEGER(KIND=JPIM) :: laytrop_min, laytrop_max
REAL(KIND=JPRB) :: ZHOOK_HANDLE

    laytrop_min = MINVAL(k_laytrop(KIDIA:KFDIA))
    laytrop_max = MAXVAL(k_laytrop(KIDIA:KFDIA))

    i_nlayers = klev
    i_laysolfr(KIDIA:KFDIA) = k_laytrop(KIDIA:KFDIA)

    DO i_lay = 1, laytrop_min
       DO iplon = KIDIA, KFDIA
!$NEC unroll(NG26)
         DO ig = 1 , ng26
           IF(i_lay == i_laysolfr(iplon)) p_sfluxzen(iplon,ig)=sfluxrefc(ig)
           p_taug(iplon,i_lay,ig) = 0.0_JPRB
           p_taur(iplon,i_lay,ig) = p_colmol(iplon,i_lay) * raylc(ig)
         ENDDO
       ENDDO
    ENDDO

    DO i_lay = laytrop_min+1, laytrop_max
       DO iplon = KIDIA, KFDIA
          IF (i_lay <= k_laytrop(iplon)) THEN
!$NEC unroll(NG26)
            DO ig = 1 , ng26
              IF(i_lay == i_laysolfr(iplon)) p_sfluxzen(iplon,ig)=sfluxrefc(ig)
              p_taug(iplon,i_lay,ig) = 0.0_JPRB
              p_taur(iplon,i_lay,ig) = p_colmol(iplon,i_lay) * raylc(ig)
            ENDDO
          ELSE
!$NEC unroll(NG26)
            DO ig = 1 , ng26
              p_taug(iplon,i_lay,ig) = 0.0_JPRB
              p_taur(iplon,i_lay,ig) = p_colmol(iplon,i_lay) * raylc(ig)
            ENDDO
          ENDIF
       ENDDO
    ENDDO

    DO ig = 1 , ng26
       DO i_lay = laytrop_max+1, i_nlayers
         DO iplon = KIDIA, KFDIA
           p_taug(iplon,i_lay,ig) = 0.0_JPRB
           p_taur(iplon,i_lay,ig) = p_colmol(iplon,i_lay) * raylc(ig)
         ENDDO
       ENDDO
    ENDDO

END SUBROUTINE SRTM_TAUMOL26
