#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <signal.h>
#include <unistd.h>

extern void abor1_(const char msg[], int msglen);

#pragma weak abor1_

#if 0
void batch_kill_()
{
#ifdef __INTEL_COMPILER
  {
    // Fixes (?) hangs Intel MPI
    char *env = getenv("SLURM_JOBID");
    if (env) {
      static char cmd[128] = "set -x; sleep 10; scancel --signal=TERM ";
      //static char cmd[128] = "set -x; sleep 10; scancel ";
      strcat(cmd,env);
      system(cmd);
    }
  }
#endif
}
#endif

#if 0
void _brexit(int errcode)
{
  batch_kill_();
  _exit(errcode);
}
#endif

// Forward declarations
void LinuxTraceBack(const char *prefix, const char *timestr, void *sigcontextptr);
void cmpi_abort_(int *rc);

void tabort_()
{
  int ret = -1;
  const int sig = SIGABRT;
  int rc = 128 + sig;
  static volatile sig_atomic_t irecur = 0;
  if (++irecur == 1) {
#if 1
    LinuxTraceBack(NULL,NULL,NULL);
    cmpi_abort_(&rc); // see ifsaux/parallel/cmpl_binding.F90 : calls MPI_ABORT with MPI_COMM_WORLD
#else
    ret = raise(sig); /* We get better DrHook & LinuxTrbk's with this than abort() aka SIGABRT */
    // abort(); -- essentially raise(SIGABRT) but with messier output (and may bypass DrHook)
    if (ret == 0) { // Means raise() was okay and tracebacks etc. DrHooks took place
      exit(128 + sig);
    }
#endif
  }
  // Still here ?? get the hell out of here ... now !!
  _exit(rc);
}

void abort_()
{
  if (abor1_) { // Call only if available
    static volatile sig_atomic_t irecur = 0;
    if (++irecur == 1) {
      const char msg[] = "Fortran ABORT()";
      abor1_(msg,strlen(msg));
    }
  }
  tabort_();
}

void _gfortran_abort()
{
  abort_();
}
