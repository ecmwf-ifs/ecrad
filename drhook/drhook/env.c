
/* env.c */

/* Implement Fortran-callable ec_getenv and ec_putenv,
   since not all environments have getenv & putenv,
   but Unix/C library always have them */

/* Author: Sami Saarinen, ECMWF, 15-Mar-2006 */


#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <signal.h>
#include <sys/types.h>
#include <sys/time.h>
#include <unistd.h>
#include <limits.h>
#include <errno.h>
#include <time.h>
#include "privpub.h"

#define EC_HOST_NAME_MAX 512

extern char **environ; /* Global Unix var */
static int numenv = 0;

void
ec_numenv_(int *n)
{ /* Returns the number of environment variables currently active */
  int j=0;
  if (environ) {
    for (; environ[j]; j++) { }
  }
  if (n) *n = j;
  numenv = j; /* Not thread-safe */
}


void
ec_numenv(int *n)
{
  ec_numenv_(n);
}


void
ec_overwrite_env_(int *do_overwrite)
{
  if (do_overwrite) {
    char *env = getenv("EC_OVERWRITE_ENV");
    if (env) {
      *do_overwrite = atoi(env);
    }
    else {
      *do_overwrite = 0;
    }
  }
}


void
ec_overwrite_env(int *do_overwrite)
{
  ec_overwrite_env_(do_overwrite);
}


void
ec_strenv_(const int *i,
	   char *value,
	   /* Hidden arguments */
	   const int valuelen)
{ /* Returns (*i)'th environment number; 
     Note: "Fortran", not "C" range between [1..numenv] */
  int j = (i && environ) ? (*i) : 0;
  memset(value, ' ', valuelen);
  if (j >= 1 && j <= numenv) {
    char *p = environ[--j];
    if (p) {
      int len = strlen(p);
      if (valuelen < len) len = valuelen;
      memcpy(value,p,len);
    }
  }
}


void
ec_strenv(const int *i,
	  char *value,
	  /* Hidden arguments */
	  const int valuelen)
{
  ec_strenv_(i, value, valuelen);
}


void
ec_getenv_(const char *s,
	   char *value,
	   /* Hidden arguments */
	   int slen,
	   const int valuelen)
{
  char *env = NULL;
  char *p = malloc(slen+1);
  if (!p) {
    fprintf(stderr,"ec_getenv_(): Unable to allocate %d bytes of memory\n", slen+1);
    ABOR1("ec_getenv_(): Unable to allocate memory");
  }
  memcpy(p,s,slen);
  p[slen]='\0';
  memset(value, ' ', valuelen);
  env = getenv(p);
  if (env) {
    int len = strlen(env);
    if (valuelen < len) len = valuelen;
    memcpy(value,env,len);
  }
  free(p);
}


void
ec_getenv(const char *s,
	  char *value,
	  /* Hidden arguments */
	  int slen,
	  const int valuelen)
{
  ec_getenv_(s, value, slen, valuelen);
}


#ifdef __NEC__
void
getenv_(const char *s,
	char *value,
	/* Hidden arguments */
	int slen,
	const int valuelen)
{
  ec_getenv_(s, value, slen, valuelen);
}
#endif

void
ec_putenv_(const char *s,
	   /* Hidden argument */
	   int slen)
{
  const char *x = &s[slen-1];
  /* strip trailing blanks first */
  while (slen > 0 && *x == ' ') { --slen; --x; }
  /* now go ahead */
  if (slen > 0) {
    char *p = malloc(slen+1);
    if (!p) {
      fprintf(stderr,"ec_putenv_(): Unable to allocate %d bytes of memory\n", slen+1);
      ABOR1("ec_putenv_(): Unable to allocate memory");
    }
    memcpy(p,s,slen);
    p[slen]='\0';
    putenv(p);
    /* Cannot free(p); , since putenv() uses this memory area for good ;-( */
  }
}


void
ec_putenv(const char *s,
	  /* Hidden argument */
	  int slen)
{
  ec_putenv_(s,slen);
}


void
ec_putenv_nooverwrite_(const char *s,
		       /* Hidden argument */
		       int slen)
{
  const char *x = &s[slen-1];
  /* strip trailing blanks first */
  while (slen > 0 && *x == ' ') { --slen; --x; }
  /* now go ahead */
  if (slen > 0) {
    char *eq = NULL;
    char *p = malloc(slen+1);
    if (!p) {
      fprintf(stderr,"ec_putenv_nooverwrite_(): Unable to allocate %d bytes of memory\n", slen+1);
      ABOR1("ec_putenv_nooverwrite_(): Unable to allocate memory");
    }
    memcpy(p,s,slen);
    p[slen]='\0';
    eq = strchr(p,'=');
    if (eq) {
      char *env = NULL;
      *eq = '\0';
      env = getenv(p);
      if (env) {
	/* Already found ==> do not overwrite */
	free(p);
	return;
      }
      else {
	/* Reset '=' back and continue with putenv() */
	*eq = '=';
      }
    }
    putenv(p);
    /* Cannot free(p); , since putenv() uses this memory area for good ;-( */
  }
}


void
ec_putenv_nooverwrite(const char *s,
		      /* Hidden argument */
		      int slen)
{
  ec_putenv_nooverwrite_(s,slen);
}


/*--- sleep_by_spinning ---*/

static int sleep_by_spinning(long secs, long nanosecs) { /* see also drhook.c */
  /* This does not call sleep() at all i.e. is not SIGALRM driven */
  int rc;
  struct timespec req, rem;
  req.tv_sec = secs;
  req.tv_nsec = nanosecs;
  rc = nanosleep(&req, &rem);
  if (rc == -1) {
    if (errno == EINTR) {
      rc = rem.tv_sec;
    }
    else
      rc = 0; /* Can't do much more about this */
  }
  return rc;
}

unsigned int
ec_sleep_(const int *nsec)
{
  //return sleep((nsec && *nsec > 0) ? *nsec : 0);
  return sleep_by_spinning((nsec && *nsec > 0) ? *nsec : 0, 0);
}


unsigned int
ec_sleep(const int *nsec)
{
  return ec_sleep_(nsec);
}

#ifdef __NEC__
void sleep_(const int *nsec) { (void)ec_sleep_(nsec); }

void flush_(const int *io) { } /* temporary fix */
#endif

/* Microsecond-sleep, by S.Saarinen, 25-jan-2008 */

void  /* Global, C-callable, too */
ec_microsleep(int usecs) {
  if (usecs > 0) {
    struct timeval t;
    t.tv_sec =  usecs/1000000;
    t.tv_usec = usecs%1000000;
    // (void) select(0, NULL, NULL, NULL, &t);
    (void) sleep_by_spinning(t.tv_sec, (long)1000*t.tv_usec);
  }
}


void
ec_usleep_(const int *usecs)
{
  if (usecs && *usecs > 0) ec_microsleep(*usecs);
}


void
ec_usleep(const int *usecs)
{
  ec_usleep_(usecs);
}

/* ec_gethostname, by S.Saarinen, 30-sep-2016 */

void ec_gethostname_(char a[], 
		     /* Hidden argument */
		     int alen)
{
  char s[EC_HOST_NAME_MAX];
  memset(a,' ',alen);
  if (gethostname(s,sizeof(s)) == 0) {
    int len;
    char *pdot = strchr(s,'.');
    if (pdot) *pdot = '\0'; // cut short from "." char e.g. hostname.fmi.fi becomes just "hostname"
    len = strlen(s);
    if (len > alen) len = alen;
    memcpy(a,s,len);
  }
}

void ec_gethostname(char a[], 
		     /* Hidden argument */
		     int alen)
{
  ec_gethostname_(a,alen);
}

#ifdef __NEC__
int hostnm_(char a[], int alen) { ec_gethostname_(a,alen); return 0; }
#endif

/* For checking runtime affinities (not setting them, though) */

#if defined(LINUX) && !defined(DARWIN) && !defined(__NEC__)
#include <sched.h>
int sched_getcpu(void);
#define getcpu() sched_getcpu()
#else
#define getcpu() -1
#endif

void ec_coreid_(int *coreid)
{
  if (coreid) *coreid = getcpu();
}

void ec_coreid(int *coreid)
{
  ec_coreid_(coreid);
}

#ifdef DARSHAN
/* Some issues with Darshan -- better to use our own version of MPI_Wtime (mpi_wtime_ in Fortran) */
double mpi_wtime_()
{
  extern double util_walltime_(); /* from drhook.c */
  return util_walltime_();
}
#endif

#if defined(__GNUC__)

/* pthread_attr_init() interception to reset guard region size 
   between thread stacks, by S.Saarinen, 30-sep-2016 */

#include <pthread.h>
#include <dlfcn.h>
#include <sys/types.h>
#include <sys/syscall.h>

#if defined(RTLD_NEXT)
#define PTR_LIBC RTLD_NEXT
#else
#define PTR_LIBC ((void*) -1L)
#endif

#ifndef SYS_gettid
#define SYS_gettid __NR_gettid
#endif

static pid_t gettid() {
#if defined(DARWIN)
  uint64_t tid64;
  pthread_threadid_np(NULL, &tid64);
  pid_t tid = (pid_t)tid64;
#else
  pid_t tid = syscall(SYS_gettid);
#endif
  return tid;
}

int getpid_() { /* GNU Fortran did not recognize this ? Here it comes then */
  return (int)getpid();
}

static int GetMe()
{
  int me = -1; /* MPI task id >= 0 && <= NPES - 1 */
  /* Trying to figure out MPI task id since are potentially doing this *before* MPI_Init*() */
  char *env_procid = getenv("ALPS_APP_PE");
  if (!env_procid) env_procid = getenv("EC_FARM_ID");
  if (!env_procid) env_procid = getenv("PMI_RANK");
  if (!env_procid) env_procid = getenv("OMPI_COMM_WORLD_RANK");
  if (env_procid) me = atoi(env_procid);
  return me;
}

static int (*ptr_pthread_attr_init)(pthread_attr_t *attr) = NULL;
int pthread_attr_init(pthread_attr_t *attr)
{
  int rc;
  static int done = 0;
  FILE *fp = NULL;
  int me = GetMe();
  pid_t pid = getpid();
  pid_t tid = gettid();
  int master = (pid == tid) ? 1 : 0;
  if (!ptr_pthread_attr_init) {
    ptr_pthread_attr_init = (int (*)(pthread_attr_t *a))dlsym(PTR_LIBC, "pthread_attr_init");
    if (!ptr_pthread_attr_init) {
      fprintf(stderr,"***Error: Dynamic linking to pthread_attr_init() failed : errno = %d\n",errno);
      abort();
    }
    /* We intend to output only from MPI-task 0, master thread */
    if (!done && me == 0 && master) fp = stderr;
    done = 1;
  }
  rc = ptr_pthread_attr_init(attr);
  {
    char *env_gs = getenv("THREAD_GUARDSIZE");
    if (env_gs) {
      int pgsize = getpagesize();
      size_t guardsize = atoll(env_gs);
      if (strchr(env_gs,'G')) guardsize *= 1073741824; /* hence, in GiB */
      else if (strchr(env_gs,'M')) guardsize *= 1048576; /* hence, in MiB */
      else if (strchr(env_gs,'K')) guardsize *= 1024; /* hence, in KiB */
      guardsize = RNDUP(guardsize,pgsize);
      if (guardsize > pgsize) { /* Now we *do* bother */
	char *env_omp = getenv("OMP_STACKSIZE");
	size_t omp_stacksize = env_omp ? atoll(env_omp) : 0;
	size_t stacksize = 0;
	int iret = pthread_attr_getstacksize(attr,&stacksize);
#if 0
	if (fp) fprintf(fp,
			"[%s@%s:%d] [pid=%ld:tid=%ld]: Requesting guard region size "
			"between thread stacks : %lld bytes (%s PAGESIZE = %d)\n",
			__FUNCTION__,__FILE__,__LINE__,
			(long int)pid,(long int)tid,
			(long long int)guardsize,
			(guardsize > pgsize) ? ">" : "<=",
			pgsize);
#endif
	if (env_omp) {
	  if (strchr(env_omp,'G')) omp_stacksize *= 1073741824; /* hence, in GiB */
	  else if (strchr(env_omp,'M')) omp_stacksize *= 1048576; /* hence, in MiB */
	  else if (strchr(env_omp,'K')) omp_stacksize *= 1024; /* hence, in KiB */
	}
	if (fp) fprintf(fp,
			"[%s@%s:%d] [pid=%ld:tid=%ld]: Stack size(s) : %lld bytes (def), %lld bytes (OMP) : [iret=%d]\n",
			__FUNCTION__,__FILE__,__LINE__,
			(long int)pid,(long int)tid,
			(long long int)stacksize,
			(long long int)omp_stacksize,
			iret);
	if (iret == 0 && omp_stacksize > guardsize) {
	  iret = pthread_attr_setguardsize(attr,guardsize);
	  (void) pthread_attr_getguardsize(attr,&guardsize);
	  if (fp) fprintf(fp,
			  "[%s@%s:%d] [pid=%ld:tid=%ld]: Guard region size now : %lld bytes : [iret=%d]\n",
			  __FUNCTION__,__FILE__,__LINE__,
			  (long int)pid,(long int)tid,
			  (long long int)guardsize,iret);
	}
      }
    }
  }
  if (fp) fflush(fp);
  return rc;
}

#if 0
/* Opting out for now */
static void MemInfoBeforeMain() __attribute__((constructor));
static void MemInfoBeforeMain()
{
  static int done = 0;
  int me = GetMe();
  if (!done && me == 0) {
    extern void meminfo_(const int *, const int *);
    const int kout = 0;
    const int kstep = -1;
    pid_t pid = getpid();
    pid_t tid = gettid();
    int master = (pid == tid) ? 1 : 0;
    if (me == 0 && master) meminfo_(&kout, &kstep); /* utilities/ec_meminfo.F90 */
    done = 1;
  }
}
#endif

#endif /* defined(__GNUC__) */
