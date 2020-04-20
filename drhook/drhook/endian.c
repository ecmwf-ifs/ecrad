/* endian.c */

/* 
   Please note: the following 2 routines
   cannot be named as "is_little_endian()"
   and "is_big_endian()", since there is a clash
   with the new Magics++ library  (by SS, 21-Mar-2006)

   --> consequently "ec_" prefix was added
*/

#include <stdio.h>
#include <stdlib.h>
#include <signal.h>
#include <stdarg.h>
#include <sys/types.h>
#include <sys/stat.h>

int ec_is_little_endian()
{
  /* Little/big-endian runtime auto-detection */
  const unsigned int ulbtest = 0x12345678;
  const unsigned char *clbtest = (const unsigned char *)&ulbtest;
  
  if (*clbtest == 0x78) { 
    /* We are on a little-endian machine */
    return 1;
  }
  else { 
    /* We are on a big-endian machine */
    return 0;
  }
}

int ec_is_big_endian()
{
  return !ec_is_little_endian();
}

/* Fortran interface */

int ec_is_big_endian_()    { return ec_is_big_endian();    }
int ec_is_little_endian_() { return ec_is_little_endian(); }

/* A routine to be called at the very end in case MPI wasn't finalized */
/* Registered *only* by MPL_INIT */
/* Disable this feature via : export EC_MPI_ATEXIT=0 */

void ec_mpi_atexit_(void)
{
  char *env = getenv("EC_MPI_ATEXIT");
  int do_it = env ? atoi(env) : 1;
  static int callnum = 0;
  ++callnum;
  if (do_it) {
    if (callnum == 1) {
      /* register */
      atexit(ec_mpi_atexit_);
    }
    else if (callnum == 2) {
      /* action : finish MPI via F90 cmpl_end (in cmpl_binding.F90) */
      extern void cmpl_end_(int *);
      int ierr = 0;
      cmpl_end_(&ierr);
    }
  }
}

void ec_mpi_atexit(void)
{
  ec_mpi_atexit_();
}

void ec_set_umask_(void)
{
  char *env = getenv("EC_SET_UMASK");
  if (env) {
    int newmask;
    int n = sscanf(env,"%o",&newmask);
    if (n == 1) {
      int oldmask = umask(newmask);
      fprintf(stderr,
	      "*** EC_SET_UMASK : new/old = %o/%o (oct), %d/%d (dec), %x/%x (hex)\n",
	      newmask,oldmask,
	      newmask,oldmask,
	      newmask,oldmask);
    } /* if (n == 1) */
  } /* if (env) */
}


/* CALL ec_raise(6) == CALL abort() */

void ec_raise_(const int *sig) { raise(*sig); }
void ec_raise(const int *sig) { ec_raise_(sig); }

/* CALL ec_exit(iexit_code) */

void ec_exit_(const int *exit_code) { exit(exit_code ? *exit_code : 0); }
void ec_exit(const int *exit_code) { ec_exit_(exit_code); }

#ifdef __NEC__
void exit_(const int *exit_code) { ec_exit_(exit_code); }
#endif

/* snprintf replacement for VPP/VPP5000's */

#if defined(VPP5000) || defined(VPP)

int snprintf(char *str, size_t size, const char *format, ...)
{
  int rc;
  va_list ap;
  va_start(ap, format);
  rc = vsprintf(str, format, ap);
  va_end(ap);
  return rc;
}

#endif
