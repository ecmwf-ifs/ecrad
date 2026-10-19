SUBROUTINE COMPUTE_LAYTROP_MIN_MAX(KIDIA,KFDIA,LAYTROP,laytrop_min,laytrop_max)
  USE PARKIND1  ,ONLY : JPIM
  IMPLICIT NONE
  INTEGER(KIND=JPIM),INTENT(IN)    :: KIDIA
  INTEGER(KIND=JPIM),INTENT(IN)    :: KFDIA
  INTEGER(KIND=JPIM),INTENT(IN)    :: LAYTROP(KIDIA:KFDIA)
  INTEGER(KIND=JPIM),INTENT(INOUT) :: laytrop_min, laytrop_max
  INTEGER(KIND=JPIM) :: iplon
  
#if defined(_OPENACC) || defined(OMPGPU)
  laytrop_min = HUGE(laytrop_min)
  laytrop_max = -HUGE(laytrop_max)
  !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO REDUCTION(min:laytrop_min) REDUCTION(max:laytrop_max) MAP(TOFROM: laytrop_min, laytrop_max)
  !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
  !$ACC LOOP GANG VECTOR REDUCTION(min:laytrop_min) REDUCTION(max:laytrop_max)
  do iplon = KIDIA,KFDIA
     laytrop_min = MIN(laytrop_min, laytrop(iplon))
     laytrop_max = MAX(laytrop_max, laytrop(iplon))
  end do
  !$ACC END PARALLEL
  !$ACC WAIT
  !$OMP END TARGET TEAMS DISTRIBUTE PARALLEL DO
#else
  laytrop_min = MINVAL(laytrop(KIDIA:KFDIA))
  laytrop_max = MAXVAL(laytrop(KIDIA:KFDIA))
#endif

END SUBROUTINE COMPUTE_LAYTROP_MIN_MAX
  
