MODULE roctx_profiling

#if defined(HAVE_ROCTX)
  INTERFACE
     SUBROUTINE roctxStartRange(message) BIND(c, name="roctxRangePushA")
       USE ISO_C_BINDING,   ONLY: C_CHAR
       IMPLICIT NONE
       CHARACTER(C_CHAR) :: message(*)
     END SUBROUTINE roctxStartRange
     SUBROUTINE roctxEndRange() BIND(c, name="roctxRangePop")
       IMPLICIT NONE
     END SUBROUTINE roctxEndRange
  END INTERFACE
#endif

END MODULE roctx_profiling
