#define _DRHOOK_C_   1

#define _GNU_SOURCE

/* 
   drhook.c

   Author: Sami Saarinen, ECMWF, 14..24-Nov-2003

   Thanks to Bob Walkup & John Hague for IBM Power4 version
   Thanks to Bob Carruthers for Cray X1 (SV2), XD1 and XT3 versions,
   as well as David Tanqueray for the flop routines

   Also thanks to Roland Richter for suggesting the use
   of "call tracebackqq()" function.
   In our environment this is accomplished by calling fortran
   routine intel_trbk() from ifsaux/utilities/gentrbk.F90.

*/

/*
If intending to run on IBM P4+ or newer systems the following definition
should be activated to use pm_initialize() instead of pm_init() of PMAPI-lib ($LIBHPM)
#define PMAPI_POST_P4
*/

/*
If *ALSO* intending to run on IBM P5+ systems, then set also BOTH
#define PMAPI_POST_P4
#define PMAPI_P5_PLUS
*/

/* Thanks to John Hague (IBM) 
 If intending to run on IBM p6 systems, then set also BOTH
#define PMAPI_POST_P4
#define PMAPI_P6
 */

#if defined(PMAPI_P7)
#define ENTRY_4 5
#define ENTRY_6 4
#elif defined(PMAPI_P6)
#define ENTRY_4 5
#define ENTRY_6 4
#elif defined(PMAPI_P5_PLUS)
#define ENTRY_4 5
#define ENTRY_6 4
#else
#define ENTRY_4 4
#define ENTRY_6 6
#endif

#if defined(SV2) || defined(XD1) || defined(XT3)
#define DT_FLOP
#define HPM
#define MAX_COUNTERS 6
#endif

#ifdef RS6K
#pragma options opt=3 halt=e
#include <pthread.h>
#endif


#include <unistd.h>
#if defined(DARWIN)
#include <pthread.h>
#endif

#define EC_HOST_NAME_MAX 512

/* === This doesn't handle recursive calls correctly (yet) === */

#include "drhook.h"
#include "cas.h"
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>

#ifdef __USE_GNU
#include <dlfcn.h>
#endif

static void set_timed_kill();
static void process_options();
static char *TimeStr(char *s, int slen);

int drhook_memtrace = 0; /* set to 1, if opt_memprof or opt_timeline ; used in getcurheap.c to lock stuff */

#if !defined(CACHELINESIZE)
#if defined(LEVEL1_DCACHE_LINESIZE)
#define CACHELINESIZE LEVEL1_DCACHE_LINESIZE
#else
/* ***Note: A hardcoded cache line size in bytes !!! */
#ifdef RS6K
#define CACHELINESIZE 128
#else
#define CACHELINESIZE 64
#endif
#endif
#endif

#include "crc.h"
#include <time.h>

static char *start_stamp = NULL;
static char *end_stamp = NULL;

static int numthreads = 0;
static int myproc = 1;
static int nproc = -1;
static int max_threads = 1;

extern int get_thread_id_();

typedef struct drhook_prefix_t {
  char s[3840];
  char timestr[256];
  int nsigs;
} drhook_prefix_t;

static drhook_prefix_t *ec_drhook = NULL;
static int timestr_len = 0;

#define PREFIX(tid) (ec_drhook && tid >= 1 && tid <= numthreads) ? ec_drhook[tid-1].s : ""
#define TIDNSIGS(tid) (ec_drhook && tid >= 1 && tid <= numthreads) ? ec_drhook[tid-1].nsigs : -1
#define TIMESTR(tid) (timestr_len > 0 && ec_drhook && tid >= 1 && tid <= numthreads) ? TimeStr(ec_drhook[tid-1].timestr,timestr_len) : ""
#define FFL __FUNCTION__,__FILE__,__LINE__

static int drhook_trapfpe_master_init = 0;
static int drhook_trapfpe = 1;
static int drhook_trapfpe_invalid = 1;
static int drhook_trapfpe_divbyzero = 1;
static int drhook_trapfpe_overflow = 1;

#if defined(NECSX)
#pragma cdir options -Nv -Csopt
extern void necsx_trbk_(const char *msg, int msglen); /* from ../utilities/gentrbk.F90 */
#endif

#if defined(LINUX) && !defined(XT3) && !defined(XD1) && !defined(CYGWIN)

#if defined(__GNUC__) && !defined(NO_TRAPFPE)

#if defined(CYGWIN)
#include <mingw/fenv.h>
#else
#include <fenv.h>
#endif

extern int feenableexcept(int excepts);
extern int fedisableexcept(int excepts);
extern int fegetexcept(void);

#if defined(DARWIN)
  /*  A temporary fix to link on MacIntosh. Something more clever will be done later -REK. */
int feenableexcept (int excepts) { return 0; }
int fedisableexcept(int excepts) { return 0; }
int fegetexcept(void) { return 0; }
#endif

#if defined(__NEC__)
int fegetexcept(void) { return 0; }
#endif

static void trapfpe(int silent)
{
  /* Enable some exceptions. At startup all exceptions are masked. */
#if 1
  /* New coding -- honours DR_HOOK_TRAPFPE_{INVALID,DIVBYZERO,OVERLOW} set to 1 (or 0) */
  int tid = get_thread_id_();
  int enable = 0;
  int disable = 0;
  int dummy;
  int rc_enable = 0;
  int rc_disable = 0;
  int excepts_before, excepts_after;
  dummy = drhook_trapfpe_invalid ? (enable |= FE_INVALID) : (disable |= FE_INVALID);
  dummy = drhook_trapfpe_divbyzero ? (enable |= FE_DIVBYZERO) : (disable |= FE_DIVBYZERO);
  dummy = drhook_trapfpe_overflow ? (enable |= FE_OVERFLOW) : (disable |= FE_OVERFLOW);
  if (!silent && myproc == 1) {
    excepts_before = fegetexcept();
  }
  if (enable) rc_enable = feenableexcept(enable); // Turn ON these
  if (disable) rc_disable = fedisableexcept(disable); // Turn OFF these
  if (!silent && myproc == 1) {
    char *pfx = PREFIX(tid);
    excepts_after = fegetexcept();
    fprintf(stderr,
	    "%s %s [%s@%s:%d] DR_HOOK trapfpe() : Exceptions before = 0x%x [%d] -- after = 0x%x [%d]\n",
	    pfx,TIMESTR(tid),FFL,
	    excepts_before, excepts_before,
	    excepts_after, excepts_after);
    fprintf(stderr,
    "%s %s [%s@%s:%d] DR_HOOK trapfpe() : with FE_INVALID = 0x%x [%d] -- FE_DIVBYZERO = 0x%x [%d] -- FE_OVERFLOW = 0x%x [%d]\n",
	    pfx,TIMESTR(tid),FFL,
	    (int)FE_INVALID, (int)FE_INVALID,
	    (int)FE_DIVBYZERO, (int)FE_DIVBYZERO,
	    (int)FE_OVERFLOW, (int)FE_OVERFLOW);
    if (enable) {
      fprintf(stderr,
	      "%s %s [%s@%s:%d] DR_HOOK trapfpe() : feenableexcept(0x%x [%d]) returns rc=%d\n",
	      pfx,TIMESTR(tid),FFL,
	      enable,enable,rc_enable);
    }
    if (disable) {
      fprintf(stderr,
	      "%s %s [%s@%s:%d] DR_HOOK trapfpe() : fedisableexcept(0x%x [%d]) returns rc=%d\n",
	      pfx,TIMESTR(tid),FFL,
	      disable,disable,rc_disable);
    }
    if (tid == 1) drhook_trapfpe_master_init = 1; // go-ahead for slave threads in trapfpe_slave_threads()
  }
#else
#if defined(PARKIND1_SINGLE) && !defined(SGEMM)
  /* For now ... we have issues in SGEMM with IEEE-invalid ... especially with LIBSCI from Cray */
  int rc = feenableexcept(FE_DIVBYZERO|FE_OVERFLOW);
#else
  int rc = feenableexcept(FE_INVALID|FE_DIVBYZERO|FE_OVERFLOW);
#endif
#endif
}

static void untrapfpe(int silent)
{
  /* Disable some exceptions. At startup all exceptions are masked. */
  int rc = fedisableexcept(FE_INVALID|FE_DIVBYZERO|FE_OVERFLOW);
}

#endif /* defined(__GNUC__) */

#endif /* defined(LINUX) && !defined(XT3) && !defined(XD1) */

#if (!defined(LINUX) || defined(CYGWIN) || defined(NO_TRAPFPE)) && defined(__GNUC__)
/* For example Solaris with gcc */
#define trapfpe(x)
#define untrapfpe(x)
#endif

#ifndef drhook_harakiri_timeout_default
#define drhook_harakiri_timeout_default 500
#endif

static int drhook_harakiri_timeout = drhook_harakiri_timeout_default;
static int drhook_use_lockfile = 1;

static int atp_enabled = 0; /* Cray ATP specific */
static int atp_max_cores = 20; /* Cray ATP specific */
static int atp_max_analysis_time = 300; /* Cray ATP specific */
static int atp_ignore_sigterm = 0; /* Cray ATP specific */

static int any_memstat = 0;
static int opt_gethwm = 0;
static int opt_getstk = 0;
static int opt_getrss = 0;
static int opt_getpag = 0;
static int opt_walltime = 0;
static int opt_cputime = 0;
static int opt_wallprof = 0;
static int opt_cpuprof = 0;
static int opt_hpmprof = 0;
static int opt_memprof = 0;
static int opt_trim = 0;
static int opt_calls = 0;
static int opt_self = 1; /* 0=exclude drhook altogether, 
			    1=include, but don't print, 
			    2=also print */
static int opt_propagate_signals = 1;
static int opt_sizeinfo = 1;
static int opt_clusterinfo = 0;
static int opt_callpath = 0;
#define callpath_indent_default 2
static int callpath_indent = callpath_indent_default;
#define callpath_depth_default 50
static int callpath_depth = callpath_depth_default;
static int callpath_packed = 0;

static int opt_calltrace = 0;
static int opt_funcenter = 0;
static int opt_funcexit = 0;

static int opt_timeline = 0; /* myproc or -1 [or 0 for --> timeline feature off (default)] */
static int opt_timeline_thread = 1; /* thread-id control : 
				    <= 0 print for all threads
				       1 -> #1 only [but curheap still SUM of all threads] (default), 
				       n -> print for increasing number of threads separately : [1..n] */
static int opt_timeline_format = 1; /* if 1, print only {wall,hwm,rss,curheap} w/o labels "wall=" etc.; else fully expanded fmt */
static int opt_timeline_unitno = 6; /* Fortran unit number : default = 6 i.e. stdout */
static long long int opt_timeline_freq = 1000000; /* How often to print : every n-th call : default = every 10^6 th call or ... */
static double opt_timeline_MB = 1.0; /* ... rss or curheap jumps up/down by more than this many MBytes (default = 1) : unit MBytes */

static volatile sig_atomic_t opt_gencore = 0;
static int opt_gencore_signal = 0;

static int hpm_grp = 0;
static int opt_random_memstat = 0; /* > 0 if to obtain random memory stats (maxhwm, maxstk) for tid=1. Updated when rand() % opt_random_memstat == 0 */

static double opt_trace_stack = 0; /* if > 0, a multiplier for OMP_STACKSIZE to monitor high master thread stack usage --
				      -- implies opt_random_memstat = 1 (regardless of DR_HOOK_RANDOM_MEMSTAT setting) 
				      -- for master MPI task only (for the moment) */
static long long int drhook_omp_stacksize = 0; /* Slave stack size -- 
						  an indicative stack size even master thread should not exceed */
static long long int drhook_stacksize_threshold = 0;
static long long int slave_stacksize();

/* Begin of developer options */
static char *drhook_timed_kill = NULL; /* Timer assisted simulated kill of procs/threads by signal */
static int drhook_dump_maps = 0; /* Print /proc/<tid>/maps from signal handler (before moving to ATP or below) */
static int drhook_dump_smaps = 0; /* Print /proc/<tid>/smaps from signal handler (before moving to ATP or below) */
static int drhook_dump_buddyinfo = 0; /* Print /proc/buddyinfo from signal handler (before moving to ATP or below) */
static int drhook_dump_meminfo = 0; /* Print /proc/meminfo from signal handler (before moving to ATP or below) */
static int drhook_dump_hugepages = 0;
static double drhook_dump_hugepages_freq = 0;
/* End of developer options */

typedef struct drhook_timeline_t {
  unsigned long long int calls[2]; /* 0=drhook_begin , 1=drhook_end */
  double last_curheap_MB;
  double last_rss_MB;
  double last_stack_MB;
  double last_vmpeak_MB;
//#if CACHELINESIZE > (2*sizeof(unsigned long long int) + 4*sizeof(double)) -- disallowed
#if CACHELINESIZE > (2*8 + 4*8)
  char pad[CACHELINESIZE - (2*sizeof(unsigned long long int) + 4*sizeof(double))]; /* padding : e.g. 64 bytes - 6*8 bytes */
#endif
} drhook_timeline_t; /* cachelinesize optimized --> less false sharing when running with OpenMP */

static drhook_timeline_t *timeline = NULL;

/* HPM-specific */

static long long int opt_hpmstop_threshold = -1;
static        double opt_hpmstop_mflops    =  1000000.0; /* Yes, 1 PetaFlop/s !! */


#define DRHOOK_STRBUF 1000

#ifndef SA_SIGINFO
#define SA_SIGINFO 0
#define SIG_EXTRA_ARGS       /* empty */
#define SIG_PASS_EXTRA_ARGS  /* empty */
#else 
#define SIG_EXTRA_ARGS       , siginfo_t *sigcode, void *sigcontextptr
#define SIG_PASS_EXTRA_ARGS  , sigcode, sigcontextptr
#endif

#define NIL "(nil)"

#undef MIN
#define MIN(a,b) ( (a) < (b) ? (a) :  (b) )

#undef MAX
#define MAX(a,b) ( (a) > (b) ? (a) :  (b) )

#undef ABS
#define ABS(x)   ( (x) >= 0  ? (x) : -(x) )

#define strequ(s1,s2)     ((void *)s1 && (void *)s2 && strcmp(s1,s2) == 0)
#define strnequ(s1,s2,n)  ((void *)s1 && (void *)s2 && memcmp(s1,s2,n) == 0)

extern long long int getstk_();
extern long long int getmaxstk_();
extern long long int gethwm_();
extern long long int getmaxhwm_();
extern long long int getrss_();
extern long long int getmaxrss_();
extern long long int getcurheap_();
extern long long int getmaxcurheap_();
extern long long int getcurheap_thread_(const int *tidnum);    /* *tidnum >= 1 && <= max_threads */
extern long long int getmaxcurheap_thread_(const int *tidnum); /* *tidnum >= 1 && <= max_threads */
extern long long int getpag_();
extern long long int getvmpeak_();

extern void ec_set_umask_();

#if defined(DT_FLOP)
extern double flop_();
#endif

extern double util_cputime_();
extern double util_walltime_();

#ifdef RS6K
static long long int irtc_start = 0;
extern long long int irtc();
#define WALLTIME() ((double)(irtc() - irtc_start)*1.0e-9)
#define CPUTIME() util_cputime_()
#elif defined(CRAYXT)
/* Cray XT3/XT4 with catamount microkernel */
#include <catamount/dclock.h>
static double dclock_start = 0;
#define WALLTIME() (dclock() - dclock_start)
#define CPUTIME()  WALLTIME()
#else
#if defined(SV2)
#include <intrinsics.h>
#endif
#if defined(XD1) || defined(XT3)
extern long long int irtc_();             /* integer*8 irtc() */
extern long long int irtc_rate_();        /* integer*8 irtc_rate() */
#endif
#if defined(SV2) || defined(XD1) || defined(XT3)
static long long int irtc_start = 0;
static double my_irtc_rate = 0;
static double my_inv_irtc_rate = 0;
#if defined(SV2)
#define WALLTIME() ((double)(_rtc() - irtc_start)*my_inv_irtc_rate)
#else
#define WALLTIME() ((double)(irtc_() - irtc_start)*my_inv_irtc_rate)
#endif
#define CPUTIME() util_cputime_()
#else
#define WALLTIME() util_walltime_()
#define CPUTIME() util_cputime_()
#endif
#endif

/* #define RAISE(x) { int tmp = x; c_drhook_raise_(&tmp); } */
#include "raise.h"
#include "cargs.h"

extern void LinuxTraceBack(const char *prefix, const char *timestr, void *sigcontextptr);

/*** typedefs ***/

typedef union {
  struct drhook_key_t *keyptr;
  double d;
  unsigned long long int ull;
} equivalence_t;

typedef struct drhook_key_t {
  char *name;
  unsigned short name_len;
  const equivalence_t *callpath; /* parent's tree down to callpath_depth */
  int callpath_len;
  unsigned int callpath_fullhash;
  unsigned short status; /* 0=inactive, >1 active */
  unsigned long long int calls;
  long long int hwm, maxrss, rssnow, stack, maxstack, paging;
  double wall_in, delta_wall_all, delta_wall_child;
  double cpu_in, delta_cpu_all, delta_cpu_child;
#ifdef HPM
  unsigned char hpm_stopped, counter_stopped;
  double this_delta_wall_child;
  double avg_mipsrate, avg_mflops;
  unsigned long long int hpm_calls;
  double mip_count_in, mflop_count_in;
  long long int *counter_in, *counter_sum;
#endif
  char *filename;         /* the filename where the 1st call (on this routine-name) 
			     to dr_hook() occurred */
  long long int sizeinfo; /* # of data elements, bytes, etc. */
  long long int min_sizeinfo, max_sizeinfo; /* min & max of # of data elements, bytes, etc. */
  /* memprof specific */
  long long int mem_seenmax;
  long long int mem_child, mem_curdelta;
  long long int maxmem_selfdelta, maxmem_alldelta;
  long long int mem_maxhwm, mem_maxrss, mem_maxstk, mem_maxpagdelta;
  long long int paging_in;
  unsigned long long int alloc_count, free_count;
  struct drhook_key_t *next;
} drhook_key_t;

typedef struct drhook_calltree_t {
  int active;
  drhook_key_t *keyptr;
  struct drhook_calltree_t *next;
  struct drhook_calltree_t *prev;
} drhook_calltree_t;

typedef struct drhook_sig_t {
  char name[32];
  struct sigaction new;
  struct sigaction old;
  int active;
  int ignore_atexit;
} drhook_sig_t;

typedef union {
  void (*func1args)(int sig);
  void (*func3args)(int sig SIG_EXTRA_ARGS);
} drhook_sigfunc_t;

typedef struct drhook_prof_t {
  double pc;
  double total;
  double self;
  unsigned long long int calls;
  double percall_ms_self;
  double percall_ms_total;
  double mipsrate, mflops, divpc;
  int index;
  int tid;
  int cluster;
  double *maxval;
  unsigned char is_max;
  char *name;
  char *filename;
  long long int sizeinfo;
  long long int min_sizeinfo, max_sizeinfo;
  double sizespeed, sizeavg;
  const equivalence_t *callpath; /* parent's tree down to callpath_depth */
  int callpath_len;
} drhook_prof_t;

typedef struct drhook_memprof_t {
  double pc;
  long long int self;
  long long int children;
  long long int hwm, rss, stk, pag, leaked;
  unsigned long long int calls, alloc_count, free_count;
  int index;
  int tid;
  int cluster;
  long long int *maxval;
  unsigned char is_max;
  char *name;
  char *filename;
  const equivalence_t *callpath; /* parent's tree down to callpath_depth */
  int callpath_len;
} drhook_memprof_t;

#define MAX_WATCH_FIRST_NBYTES 8

typedef struct drhook_watch_t {
  char *name;
  int tid;
  int active;
  int abort_if_changed;
  const char *ptr;
  int nbytes;
  int watch_first_nbytes;
  char first_nbytes[MAX_WATCH_FIRST_NBYTES];
  unsigned int crc32;
  int printkey;
  int nvals;
  struct drhook_watch_t *next;
} drhook_watch_t;

/*** static (local) variables ***/

static o_lock_t DRHOOK_lock = 0;
static pid_t pid = -1;
static drhook_key_t      **keydata  = NULL;
static drhook_calltree_t **calltree = NULL;
static drhook_calltree_t **thiscall = NULL;
static int signals_set = 0;
static volatile sig_atomic_t signal_handler_called = 0;
static volatile sig_atomic_t signal_handler_ignore_atexit = 0;
static volatile sig_atomic_t unlimited_corefile_retcode = 9999;
static volatile unsigned long long int saved_corefile_hardlimit = 0;
static int allow_coredump = -1; /* -1 denotes ALL MPI-tasks, 1..NPES == myproc, 0 = coredump will not be enabled by DrHook at init */
static drhook_sig_t siglist[1+NSIG] = { 0 };
static char *a_out = NULL;
static char *mon_out = NULL;
static int mon_out_procs = -1;
static double percent_limit = -10; /* Lowest percentage accepted into the printouts */
static drhook_key_t **keyself = NULL; /* pointers to itself (per thread) */
static double *overhead; /* Total Dr.Hook-overhead for every thread in either WALL or CPU secs */
static drhook_key_t **curkeyptr = NULL; /* pointers to current keyptr (per thread) */
static drhook_watch_t *watch = NULL;
static drhook_watch_t *last_watch = NULL;
static int watch_count = 0; /* No. of *active* watch points */

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

// Fortran callable : CALL GETTID_C(ITID) where INTEGER(KIND=4) :: ITID 

void gettid_c_(int *tid)
{
  if (tid) *tid = (int)gettid();
}

void gettid_c(int *tid) { gettid_c_(tid); } 

static void set_ec_drhook_label(const char *hostname, int hlen)
{
  int tid = get_thread_id_();
  int j = tid - 1;
  int slen = sizeof(ec_drhook[j].s);
  pid_t unixtid = gettid();
  snprintf(ec_drhook[j].s,slen,"[EC_DRHOOK:%*s:%d:%d:%lld:%lld]",
	   hlen,hostname,myproc,tid,
	   (long long int)pid, (long long int)unixtid);
}

#define SECS(x) ((int)(x))
#define NSECS(x) ((int)(1000000000 * ((x) - SECS(x))))

#ifndef __timer_t_defined
static void set_killer_timer(const int *ntids, const int *target_omptid, 
                             const int *target_sig, const double *start_time,
                             const char *p, int lenp)
{
  // Definition of timer_t, timer_create, timer_set
  //   is a POSIX extention, not available on e.g. Darwin
}
#else
static void set_killer_timer(const int *ntids, const int *target_omptid, 
                             const int *target_sig, const double *start_time,
                             const char *p, int lenp)
{
  static volatile sig_atomic_t TimedKill = 0;
  if (ntids && target_omptid && target_sig && start_time && p) {
    int tid = get_thread_id_();
    if (*target_omptid == -1 || *target_omptid == tid) {
      char *pfx = PREFIX(tid);
      timer_t timerid = { 0 };
      struct itimerspec its = { 0 } ;
      struct sigevent sev = { 0 } ;
      sev.sigev_signo = *target_sig;

#if defined(SIGEV_THREAD_ID)
      sev.sigev_notify = SIGEV_THREAD_ID | SIGEV_SIGNAL;
      /* sev.sigev_notify_thread_id = gettid(); */
      sev._sigev_un._tid = gettid();
#elif defined(SIGEV_THREAD)
      sev.sigev_notify = SIGEV_THREAD | SIGEV_SIGNAL;
#else
      sev.sigev_notify = SIGEV_SIGNAL;
#endif
      sev.sigev_value.sival_ptr = &timerid;
      
      its.it_value.tv_sec = SECS(*start_time);
      its.it_value.tv_nsec = NSECS(*start_time);
      
      its.it_interval.tv_sec = 0;
      its.it_interval.tv_nsec = 0;
      
      timer_create(CLOCK_MONOTONIC, &sev, &timerid);
      /* timer_create(CLOCK_REALTIME, &sev, &timerid); */
      timer_settime(timerid, 0, &its, NULL);
      
      cas_lock(&TimedKill);
      {
	fprintf(stderr,
		"%s %s [%s@%s:%d] Developer timer (%s) expires"
		" after %.3fs through signal#%d (ntids=%d)\n",
		pfx,TIMESTR(tid),FFL,
		p,
		*start_time, *target_sig, *ntids);
	fflush(NULL);
      }
      cas_unlock(&TimedKill);
    } /* if (target_omptid == -1 || target_omptid == tid) */
  }
}
#endif

#if !defined(NCALLSTACK)
#ifdef PARKIND1_SINGLE
/* > 0 : USE call stack approach : needed for single precision version */
#define NCALLSTACK 64
#else
/* == 0 : do NOT use call stack approach : usually for double precision version */
#define NCALLSTACK  0
#endif
#endif
static int cstklen = NCALLSTACK;

#define HASHSIZE(n) ((unsigned int)1<<(n))
#define HASHMASK(n) (HASHSIZE(n)-1)

#define NHASH    16
#define NHASHMAX 24
static int nhash = NHASH;
static unsigned int hashsize = HASHSIZE(NHASH);
static unsigned int hashmask = HASHMASK(NHASH);

#ifdef HPM
/* HPM-specific (static) protos */

static void stopstart_hpm(int tid, drhook_key_t *pstop, drhook_key_t *pstart);
static void stop_only_hpm(int tid, drhook_key_t *pstop);
static void init_hpm(int tid);
static double mflops_hpm(const drhook_key_t *keyptr);
static double mips_hpm(const drhook_key_t *keyptr);
static double divpc_hpm(const drhook_key_t *keyptr);
static double mflop_count(const drhook_key_t *keyptr);
static double mip_count(const drhook_key_t *keyptr);

#else
/* Dummies for HPM as macros that do nothing */

#define stopstart_hpm(tid, pstop, pstart)
#define stop_only_hpm(tid, pstop)
#define init_hpm(tid)
#define mflops_hpm(keyptr)  0
#define mips_hpm(keyptr)    0
#define divpc_hpm(keyptr)   0
#define mflop_count(keyptr) 0
#define mip_count(keyptr)   0

#endif

/*--- spin ---*/

static int nanospin(int secs, int nanosecs) {
  struct timespec req, rem;
  req.tv_sec = secs;
  req.tv_nsec = nanosecs;
  return nanosleep(&req, &rem);
}

static int spin(int secs) {
  return nanospin(secs, 0);
}


/*--- dump_file ---*/

static void dump_file(const char *pfx, int tid, int sig, int nsigs, const char filename[])
{
  /* Developer option: Will this spoil our ATP trace ... ? */
  FILE *fp;
  char in[256];
  char *tst = TIMESTR(tid);
  if (sig > 0 && nsigs >= 1) {
    fprintf(stderr,
	    "%s %s [%s@%s:%d] Content of the file '%s', signal#%d, nsigs = %d\n",
	    pfx,tst,FFL,filename,sig,nsigs);
  }
  else {
    fprintf(stderr,
	    "%s %s [%s@%s:%d] Content of the file '%s'\n",
	    pfx,tst,FFL,filename);
  }
  fp = fopen(filename,"r");
  if (fp) {
    while (fgets(in,sizeof(in),fp) == in) {
       fprintf(stderr,"%s %s [%s@%s:%d] %s",pfx,tst,FFL,in);
       /* fprintf(stderr,"%s",in); */
    }
    fclose(fp);
  }
}

/*--- dump_hugepages ---*/

// Forward declaration of subroutine in ec_meminfo.F90
void ec_meminfo_( const int* ku,
                  const char* cdstring,
                  const int* kcomm,
                  const int* kbarr,
                  const int* kiotask,
                  const int* kcall,
                  int cdstring_strlen );

static void dump_hugepages(int enforce, const char *pfx, int tid, int sig, int nsigs)
{
  if (enforce || drhook_dump_hugepages) {
    if (enforce || tid == 1) { /* OML-thread id >= 1 */
      static double next_scheduled = -1;
      double wt = WALLTIME();
      if (enforce || wt > next_scheduled) {
	const int kcomm = -1;
	const int kbarr = 0;
	const int kiotask = 0;
	const int kcall = -1;
	const int ftnunitno = 0; /* stderr */
	fflush(NULL);
	ec_meminfo_(&ftnunitno,pfx,&kcomm,&kbarr,&kiotask,&kcall,strlen(pfx));
	fflush(NULL);
	if (drhook_dump_buddyinfo) {
	  dump_file(pfx,tid,sig,nsigs,"/proc/buddyinfo");
	}
	if (drhook_dump_meminfo) {
	  dump_file(pfx,tid,sig,nsigs,"/proc/meminfo");
	}
	wt = WALLTIME();
	next_scheduled = wt + drhook_dump_hugepages_freq;
      }
    }
  }
}


/*--- set_default_handler ---*/

static int set_unlimited_corefile(unsigned long long int *hardlimit);

static int set_default_handler(int sig, int unlimited_corefile, int verbose)
{
  int rc = -2;
  if (sig >= 1 && sig <= NSIG) {
    unsigned long long int hardlimit = 0;
    struct sigaction sa = { 0 };
    sa.sa_handler = SIG_DFL;
    sigemptyset(&sa.sa_mask);
    /*
      sigfillset(&sa.sa_mask); -- if we wanted to block all (catchable) signals whilst in subsequent signal handler SIG_DFL
      sigaddset(&sa.sa_mask, some_signal_to_be_blocked); ... just in case
    */
    sigaction(sig, &sa, NULL);
    if (unlimited_corefile) rc = set_unlimited_corefile(&hardlimit); /* unconditionally */
    if (verbose) {
      int tid = get_thread_id_();
      char *pfx = PREFIX(tid);
      char buf[128] = "";
      if (unlimited_corefile && rc == 0) snprintf(buf,sizeof(buf)," -- hardlimit for core file is now %llu (0x%llx)", hardlimit, hardlimit);
      fprintf(stderr,
	      "%s %s [%s@%s:%d] "
	      "Enabled default signal handler (SIG_DFL) for signal#%d%s\n",
	      pfx,TIMESTR(tid),FFL,
	      sig,buf);
    }
  }
  return rc;
}

/*--- malloc_drhook ---*/

static void *
malloc_drhook(size_t size)
{
  size_t size1 = MAX(1,size);
  void *p = malloc(size1);
  if (!p) {
    fprintf(stderr,
	    "***Error in malloc_drhook(): Unable to allocate space for %lld bytes\n", 
	    (long long int)size1);
    RAISE(SIGABRT);
  }
  return p;
}

/*--- calloc_drhook ---*/

static void *
calloc_drhook(size_t nmemb, size_t size)
{
  size_t n = nmemb * size;
  void *p = malloc_drhook(n);
  memset(p,0,n);
  return p;
}

/*--- free_drhook ---*/

#define free_drhook(x) { if (x) { free(x); x = NULL; } }

/*--- callstack ---*/

/* Note: For single precision calls -- small performance penalty */

typedef struct callstack_t {
  drhook_key_t **keyptr;
  unsigned int next;
  unsigned int maxdepth;
} callstack_t;

static callstack_t **cstk = NULL;

static drhook_key_t *callstack(int tid, void *key, drhook_key_t *keyptr)
{
  /* Single routine -- two usages:

     (1) Upon c_drhook_start_() we call:

         (void) callstack(tid, key, u.keyptr);
	 - store keyptr into thread specific call stack
	 - fill *key up to 4-bytes index stating the position in the aforementioned call stack

     (2) Upon c_drhook_end_() we call:

         u.keyptr = callstack(tid, (void *)key, NULL);
	 - pass 4-byte index in
         - obtain keyptr from call stack
	 - decrement call stack

  */

  static const unsigned int inc = 64;
  unsigned int idx, *Index = key;
  callstack_t *c = cstk[tid-1];
  if (keyptr) {
    if (!c) {
      cstk[tid-1] = c = calloc_drhook(1, sizeof(*c));
      c->keyptr = (drhook_key_t **) calloc_drhook(cstklen, sizeof(drhook_key_t *));
      c->next = 0;
      c->maxdepth = cstklen;
    }
    idx = (c->next)++;
    if (idx >= c->maxdepth) {
      drhook_key_t **kptr;
      unsigned int maxdepth = idx + inc;
      char *pfx = PREFIX(tid);
      fprintf(stderr,
	      "%s %s [%s@%s:%d] "
	      "Call stack index %u out of range [0,%u) : extending the range to [0,%u) for this thread\n",
	      pfx,TIMESTR(tid),FFL,
	      idx,c->maxdepth,maxdepth);
      kptr = (drhook_key_t **) calloc_drhook(maxdepth, sizeof(drhook_key_t *));
      memcpy(kptr,c->keyptr,c->maxdepth * sizeof(drhook_key_t *));
      free_drhook(c->keyptr);
      c->keyptr = kptr;
      c->maxdepth = maxdepth;
    }
    if (idx >= c->maxdepth) {
      char *pfx = PREFIX(tid);
      fprintf(stderr,
	      "%s %s [%s@%s:%d] "
	      "Call stack index %u still out of range [0,%u). Aborting ...\n",
	      pfx,TIMESTR(tid),FFL,
	      idx,c->maxdepth);
      RAISE(SIGABRT);
    }
    c->keyptr[idx] = keyptr;
    *Index = idx;
  }
  else {
    idx = --(c->next);
    if (idx != *Index) {
      char *pfx = PREFIX(tid);
      fprintf(stderr,
	      "%s %s [%s@%s:%d] "
	      "Invalid index to call stack %u : out of range [0,%u). Expecting the exact value of %u\n",
	      pfx,TIMESTR(tid),FFL,
	      idx,c->maxdepth,*Index);
      RAISE(SIGABRT);
    }
    keyptr = c->keyptr[idx];
  }
  return keyptr;
}

/*--- strdup_drhook ---*/

static char *
strdup_drhook(const char *s)
{
  int n = strlen(s);
  char *p = malloc_drhook(n+1);
  memcpy(p,s,n);
  p[n] = 0;
  return p;
}

/*--- strdup2_drhook ---*/

static char *
strdup2_drhook(const char *s, int s_len)
{
  int n = s_len;
  char *p = malloc_drhook(n+1);
  memcpy(p,s,n);
  p[n] = 0;
  return p;
}

/*--- timestamp ---*/

static char *
timestamp()
{
  time_t tp;
  const int bufsize = 64;
  char *buf = malloc_drhook(bufsize+1);
  time(&tp);
  strftime(buf, bufsize, "%Y%m%d %H%M%S", localtime(&tp));
  return buf;
}

/*--- TimeStr ---*/

static char *
TimeStr(char *s, int slen)
{
  if (s) {
    time_t tp;
    char buf[64];
    time(&tp);
    strftime(buf, sizeof(buf), "%Y%m%d:%H%M%S", localtime(&tp));
    snprintf(s,slen,"[%s:%lld:%.3f]",buf,(long long int)tp,WALLTIME());
  }
  return s;
}

/* -- These 2 extern's are called primarily from LinuxTrbk() */

const char *drhook_TIMESTR(int tid)
{
  static const char fixed[] = "";
  if (tid <= 0) coml_my_thread_(&tid);
  {
    char *s = TIMESTR(tid);
    return strlen(s) > 0 ? (const char *)s : fixed;
  }
}

const char *drhook_PREFIX(int tid)
{
  static const char fixed[] = "";
  if (tid <= 0) coml_my_thread_(&tid);
  {
    char *s = PREFIX(tid);
    return strlen(s) > 0 ? (const char *)s : fixed;
  }
}

/*--- hashfunc ---*/

unsigned int
hashfunc(const char *s, int s_len)
{
  unsigned int hashval;
  if (opt_trim) {
    for (hashval = 0; s_len>0 ; s++, s_len--) {
      unsigned char c = islower(*s) ? toupper(*s) : *s;
      hashval = (hashval<<4)^(hashval>>28)^(c);
    }
  }
  else {
    for (hashval = s_len; s_len>0 ; s_len--) {
      hashval = (hashval<<4)^(hashval>>28)^(*s++);
    }
  }
  hashval = (hashval ^ (hashval>>10) ^ (hashval>>20)) & hashmask;
  return hashval;
}

/*--- callpath_hashfunc ---*/

unsigned int
callpath_hashfunc(unsigned int inithash, /* from hashfunc() */
		  const equivalence_t *callpath, int callpath_len,
		  unsigned int *fullhash)
{
  unsigned int hashval;
  for (hashval = inithash; callpath_len>0 ; callpath++, callpath_len--) {
    hashval = (hashval<<4)^(hashval>>28)^(callpath->ull);
  }
  if (fullhash) *fullhash = hashval;
  hashval = (hashval ^ (hashval>>10) ^ (hashval>>20)) & hashmask;
  return hashval;
}

/*--- insert_calltree ---*/

static void
insert_calltree(int tid, drhook_key_t *keyptr)
{
  if (tid >= 1 && tid <= numthreads) {
    drhook_calltree_t *treeptr = thiscall[tid-1];
    while (treeptr->active) {
      if (!treeptr->next) {
	treeptr->next = calloc_drhook(1,sizeof(drhook_calltree_t));
	treeptr->next->prev = treeptr;
      }
      treeptr = treeptr->next;
    }
    treeptr->keyptr = keyptr;
    treeptr->active = 1;
    thiscall[tid-1] = treeptr;
#ifdef HPM
    if (opt_hpmprof) {
      drhook_key_t *kptr = treeptr->keyptr;
      if (!kptr->hpm_stopped) {
	stopstart_hpm(tid,
		      treeptr->prev ? treeptr->prev->keyptr : NULL, /* stop current (i.e. my parent) */
		      kptr);                             /* start to gather for me */
	kptr->this_delta_wall_child = 0;
	kptr->mip_count_in = mip_count(kptr);
	kptr->mflop_count_in = mflop_count(kptr);
#ifdef DEBUG
	fprintf(stderr,"insert[%.*s@%d]: this_delta_wall_child=%.15g, mip#%.15g, mflop#%.15g\n",
		kptr->name_len,kptr->name,
		tid,kptr->this_delta_wall_child,
		kptr->mip_count_in,kptr->mflop_count_in);
#endif
      }
      else {
	stop_only_hpm(tid,
		      treeptr->prev ? treeptr->prev->keyptr : NULL /* stop current (i.e. my parent) */);
      } /* if (!kptr->hpm_stopped) else */
    } /* if (opt_hpmprof) */
#endif
  }
}

/*--- remove_calltree ---*/

static void 
remove_calltree(int tid, drhook_key_t *keyptr, 
		const double *delta_wall, const double *delta_cpu)
{
  if (tid >= 1 && tid <= numthreads) {
    drhook_calltree_t *treeptr = thiscall[tid-1];
    if (treeptr->active && treeptr->keyptr == keyptr) {
      treeptr->active = 0;
      if (treeptr->prev) {
	drhook_key_t *parent_keyptr = treeptr->prev->keyptr;
	if (parent_keyptr) { /* extra security */
	  if (opt_walltime) {
	    parent_keyptr->delta_wall_child += (*delta_wall);
#ifdef HPM
	    if (opt_hpmprof) parent_keyptr->this_delta_wall_child += (*delta_wall);
#endif
	  }
	  if (opt_cputime)  {
	    parent_keyptr->delta_cpu_child  += (*delta_cpu);
	  }
	  if (opt_memprof) {
	    /*
	    const long long int size = 0;
	    c_drhook_memcounter_(&tid, &size, NULL);
	    fprintf(stderr,
		    ">parent(%.*s)->mem_child = %lld ; this(%.*s)->alldelta = %lld, mem_child = %lld\n",
		    parent_keyptr->name_len, parent_keyptr->name, parent_keyptr->mem_child,
		    keyptr->name_len, keyptr->name, keyptr->maxmem_alldelta, keyptr->mem_child);
	    */
	    parent_keyptr->mem_child = MAX(parent_keyptr->mem_child, keyptr->maxmem_alldelta);
	    /*
	    fprintf(stderr,
		    "<parent(%.*s)->mem_child = %lld ; this(%.*s)->alldelta = %lld, mem_child = %lld\n",
		    parent_keyptr->name_len, parent_keyptr->name, parent_keyptr->mem_child,
		    keyptr->name_len, keyptr->name, keyptr->maxmem_alldelta, keyptr->mem_child);
	    */
	  }
	} /* if (parent_keyptr) */
	thiscall[tid-1] = treeptr->prev;
      }
      else {
	thiscall[tid-1] = calltree[tid-1];
      }
#ifdef HPM
      if (opt_hpmprof) {
	drhook_key_t *kptr = treeptr->keyptr;
	if (!kptr->hpm_stopped) {
	  double this_delta_wall_self = *delta_wall - kptr->this_delta_wall_child;
	  stopstart_hpm(tid, 
			kptr, 
			thiscall[tid-1]->keyptr); /* stop current, (re-)start previous */
	  /* Calculate moving average of mipsrate & mflops ; divpc we don't bother */
#ifdef DEBUG
	  fprintf(stderr,"remove[%.*s@%d]: this_delta_wall_self=%.15g i.e. %.15g - %.15g",
		  kptr->name_len,kptr->name,
		  tid,this_delta_wall_self,
		  *delta_wall,kptr->this_delta_wall_child);
#endif
	  if (this_delta_wall_self > 0) {
	    long long int hpm_calls = ++kptr->hpm_calls;
	    double mipsrate, mflops;
	    kptr->mip_count_in = mip_count(kptr) - kptr->mip_count_in;
	    kptr->mflop_count_in = mflop_count(kptr) - kptr->mflop_count_in;
	    mipsrate = kptr->mip_count_in/this_delta_wall_self;
	    kptr->avg_mipsrate = ((hpm_calls-1)*kptr->avg_mipsrate + mipsrate)/hpm_calls;
	    mflops = kptr->mflop_count_in/this_delta_wall_self;
	    kptr->avg_mflops = ((hpm_calls-1)*kptr->avg_mflops + mflops)/hpm_calls;
#ifdef DEBUG
	    fprintf(stderr,
		    ", mip#%.15g, mflop#%.15g : mipsrate=%.15g, avg=%.15g; mflops=%.15g, avg=%.15g",
		    kptr->mip_count_in,kptr->mflop_count_in,
		    mipsrate, kptr->avg_mipsrate,
		    mflops, kptr->avg_mflops);
#endif
	  }
#ifdef DEBUG
	  fprintf(stderr,"\n");
#endif
	  if (opt_hpmstop_threshold > 0 && kptr->calls == opt_hpmstop_threshold) {
	    /* check whether hpm should anymore be called for this routine */
	    if (kptr->avg_mflops < opt_hpmstop_mflops) kptr->hpm_stopped = 1;
	  }
	}
	else {
	  stop_only_hpm(tid,kptr);
	} /* if (!kptr->hpm_stopped) else ... */
      } /* if (opt_hpmprof) */
#endif
      curkeyptr[tid-1] = thiscall[tid-1]->keyptr;
    }
    else {
      curkeyptr[tid-1] = NULL;
    } /* if (treeptr->active && treeptr->keyptr == keyptr) else ... */
  }
}

/*--- memstat ---*/

static long long int
slave_stacksize()
{
  char *env_omp = getenv("OMP_STACKSIZE");
  long long int stacksize = env_omp ? atoll(env_omp) : 0;
  if (env_omp) {
    if      (strchr(env_omp,'G')) stacksize *= (long long int)1073741824; /* hence, in GiB */
    else if (strchr(env_omp,'M')) stacksize *= (long long int)1048576; /* hence, in MiB */
    else if (strchr(env_omp,'K')) stacksize *= (long long int)1024; /* hence, in KiB */
  }
  if (stacksize < 0) stacksize = 0;
  return stacksize;
}

static void
memstat(drhook_key_t *keyptr, const int *thread_id, int in_getkey)
{
  if (any_memstat && keyptr) {
    if (opt_gethwm) keyptr->hwm = gethwm_();
    if (opt_getrss) {
      keyptr->maxrss = getrss_();
      keyptr->rssnow = getcurheap_thread_(thread_id);
    }
    if (opt_getstk) {
      long long int stk = getstk_();
      keyptr->stack = stk;
      keyptr->maxstack = MAX(keyptr->maxstack,stk);
    }
    if (opt_getpag) keyptr->paging = getpag_();
    if (opt_memprof) {
      keyptr->mem_seenmax = getmaxcurheap_thread_(thread_id);
      if (in_getkey) { /* Upon enter of a Dr.Hook'ed routine */
	/* A note for "keyptr->mem_curdelta": 
	   1) do not reset to 0
	   2) initially calloc'ed to 0 while initializing the keydata[] ~ alias keyptr
	   3) remember the previous value --> catches memory leaks, too !! */
	/* keyptr->mem_curdelta = 0; */
	/* Nearly the same holds for "keyptr->mem_child"; 
	   we need to capture the maximum/hwm for child */
	/* keyptr->mem_child = 0; */
	keyptr->paging_in = keyptr->paging;
      }
      else { /* Upon exit of a Dr.Hook'ed routine */
	long long int alldelta = keyptr->mem_curdelta + keyptr->mem_child;
	if (alldelta > keyptr->maxmem_alldelta) keyptr->maxmem_alldelta = alldelta;
	if (keyptr->paging - keyptr->paging_in > keyptr->mem_maxpagdelta)
	  keyptr->mem_maxpagdelta = keyptr->paging - keyptr->paging_in;
      }
      if (keyptr->hwm      > keyptr->mem_maxhwm) keyptr->mem_maxhwm = keyptr->hwm;
      if (keyptr->maxrss   > keyptr->mem_maxrss) keyptr->mem_maxrss = keyptr->maxrss;
      if (keyptr->maxstack > keyptr->mem_maxstk) keyptr->mem_maxstk = keyptr->maxstack;
    }
  }
}

/*--- flptrap ---*/

/*
  -----------------------------------------------------------------------
  If we are trapping Floating-Point Error, then set the processor in SYNC
  modes and enable TRP_INVALID, TRP_DIV_BY_ZERO and TRP_OVERFLOW.
  -----------------------------------------------------------------------
*/

#ifdef RS6K
static void
flptrap(int sig, int silent)
{
  if (sig == SIGFPE) {
    /* From John Hague, IBM, UK (--> thanks a lot, John !!)*/
    int ret = fp_trap(FP_TRAP_FASTMODE);
    if ((ret == FP_TRAP_UNIMPL) || (ret == FP_TRAP_ERROR)) {
      char errmsg[4096];
      sprintf(errmsg, 
      "flptrap(): Call to 'fp_trap' in signal_trap failed (return code = %d)\n (line %d in file %s)\n",
      ret, __LINE__, __FILE__);
      perror(errmsg);
      RAISE(SIGABRT);
    }
    fp_enable(TRP_INVALID | TRP_DIV_BY_ZERO | TRP_OVERFLOW);
  }
}
#elif defined(__GNUC__) && !defined(NO_TRAPFPE)
static void
flptrap(int sig, int silent)
{
  if (sig == SIGFPE) {
    /* Adapted from www.twinkle.ws/arnaud/CompilerTricks.html#Glibc_FP */
    trapfpe(silent); /* No need for pgf90's -Ktrap=fp  now ? */
  }
}
#else
static void
flptrap(int sig, int silent)
{
  return; /* A dummy */
}
#endif

static void signal_gencore(int sig SIG_EXTRA_ARGS);
static void signal_harakiri(int sig SIG_EXTRA_ARGS);
static void signal_drhook(int sig SIG_EXTRA_ARGS);
static void trapfpe_treatment(int sig, int silent);

/*--- catch_signals ---*/

#define CATCHSIG(x) {\
  drhook_sig_t *sl = &siglist[x];\
  if (sl->active == 0) {\
    drhook_sigfunc_t u;\
    u.func3args = signal_drhook;\
    sl->active = 1;\
    sigemptyset(&sl->new.sa_mask);\
    sl->new.sa_handler = u.func1args;\
    sl->new.sa_flags = SA_SIGINFO;\
    sigaction(x,&sl->new,&sl->old);\
    trapfpe_treatment(x,silent);   \
    if (!silent && myproc == 1) {\
      int tid = get_thread_id_(); \
      char *pfx = PREFIX(tid); \
      fprintf(stderr,\
	      "%s %s [%s@%s:%d] DR_HOOK also catches signal#%d : New handler '%s' installed at %p (old at %p)\n", \
	      pfx,TIMESTR(tid),FFL,					\
	      x, "signal_drhook", sl->new.sa_handler, sl->old.sa_handler); \
    }\
  }\
}

static void
catch_signals(int silent)
{
  char *env = getenv("DR_HOOK_CATCH_SIGNALS");
  if (!silent && myproc == 1) {
    int tid = get_thread_id_();
    char *pfx = PREFIX(tid);
    fprintf(stderr,
	    "%s %s [%s@%s:%d] DR_HOOK_CATCH_SIGNALS=%s\n",
	    pfx,TIMESTR(tid),FFL,
	    env ? env : "<undef>");
  }
  if (env) {
    const char delim[] = ", \t/";
    char *p, *s = strdup_drhook(env);
    p = strtok(s,delim);
    while (p) {
      int sig = atoi(p);
      if (sig >= 1 && sig <= NSIG) {
	CATCHSIG(sig);
      }
      else if (sig == -1) { /* Makes ALL (catchable) signals available to DR_HOOK */
	int j;
	for (j=1; j<=NSIG; j++) {
	  CATCHSIG(j);
	} /* for (j=1; j<=NSIG; j++) */
	break;
      }
      p = strtok(NULL,delim);
    }
    free_drhook(s);
  }
}

/*--- trapfpe_treatment ---*/

static void
trapfpe_treatment(int sig, int silent)
{
  if (sig == SIGFPE) {
#if defined(__GNUC__) && !defined(NO_TRAPFPE)
    int tid = get_thread_id_();
    char *pfx = PREFIX(tid);
    if (drhook_trapfpe) {
      if (!silent && myproc == 1) {
	fprintf(stderr,
		"%s %s [%s@%s:%d] DR_HOOK enables SIGFPE-related floating point trapping since DRHOOK_TRAPFPE=%d\n",
		pfx,TIMESTR(tid),FFL,
		drhook_trapfpe);
      }
      flptrap(sig,silent); /* Has FLP-trapping on, regardless */
    }
    else {
      if (!silent && myproc == 1) {
	fprintf(stderr,
		"%s %s [%s@%s:%d] DR_HOOK turns SIGFPE-related floating point trapping off since DRHOOK_TRAPFPE=%d\n",
		pfx,TIMESTR(tid),FFL,
		drhook_trapfpe);
      }
      untrapfpe(silent); /* Turns off a possible -Ktrap=fp from pgf90 */
    }
#endif
  }
}

/* Fortran callable : calls trapfpe() for slave threads if drhook_trapfpe indicated so
   Called from DR_HOOK_UTIL_MULTI after DR_HOOK_UTIL (master thread) has been called
   Matters only for slave threads
   If *silent = 0, then more verbose output */

void
trapfpe_slave_threads_(const int *silent)
{
  int tid = get_thread_id_();
  if (tid > 1) { // slave threads
    if (drhook_trapfpe_master_init) trapfpe_treatment(SIGFPE, *silent);
  }
}

void
trapfpe_slave_threads(const int *silent)
{
  trapfpe_slave_threads_(silent);
}


/*--- restore_default_signals ---*/

static void
restore_default_signals(int silent)
{
  char *env = getenv("DR_HOOK_RESTORE_DEFAULT_SIGNALS");
  if (!silent && myproc == 1) {
    int tid = get_thread_id_();
    char *pfx = PREFIX(tid);
    fprintf(stderr,
	    "%s %s [%s@%s:%d] DR_HOOK_RESTORE_DEFAULT_SIGNALS=%s\n",
	    pfx,TIMESTR(tid),FFL,
	    env ? env : "<undef>");
  }
  if (env) {
    int unlim_core = 1;
    const char delim[] = ", \t/";
    char *p, *s = strdup_drhook(env);
    p = strtok(s,delim);
    while (p) {
      int sig = atoi(p);
      if (sig >= 1 && sig <= NSIG) {
	drhook_sig_t *sl = &siglist[sig];
	if (sl->active == 0) { /* Not touched yet by ignore_signals() */
	  set_default_handler(sig,unlim_core,(!silent && myproc == 1));
	  unlim_core = 0;
	  if (sig == SIGFPE) trapfpe_treatment(sig, (!silent && myproc == 1));
	  sl->active = -2;
	}
      }
      else if (sig == -1) { /* Restore default signals for all available/catchable to DR_HOOK */
	int j;
	for (j=1; j<=NSIG; j++) {
	  drhook_sig_t *sl = &siglist[j];
	  if (sl->active == 0) {  /* Not touched yet by ignore_signals() */
	    set_default_handler(j,unlim_core,(!silent && myproc == 1));
	    unlim_core = 0;
	    if (j == SIGFPE) trapfpe_treatment(j, (!silent && myproc == 1));
	    sl->active = -2;
	  }
	} /* for (j=1; j<=NSIG; j++) */
	break;
      }
      p = strtok(NULL,delim);
    }
    free_drhook(s);
  }
}

/*--- ignore_signals ---*/

static void
ignore_signals(int silent)
{
  char *env = getenv("DR_HOOK_IGNORE_SIGNALS");
  if (!silent && myproc == 1) {
    int tid = get_thread_id_();
    char *pfx = PREFIX(tid);
    fprintf(stderr,
	    "%s %s [%s@%s:%d] DR_HOOK_IGNORE_SIGNALS=%s\n",
	    pfx,TIMESTR(tid),FFL,
	    env ? env : "<undef>");
  }
  if (env) {
    int tid = get_thread_id_();
    char *pfx = PREFIX(tid);
    const char delim[] = ", \t/";
    char *p, *s = strdup_drhook(env);
    p = strtok(s,delim);
    while (p) {
      int sig = atoi(p);
      if (sig >= 1 && sig <= NSIG) {
	drhook_sig_t *sl = &siglist[sig];
	if (!silent && myproc == 1) {
	  fprintf(stderr,
		  "%s %s [%s@%s:%d] DR_HOOK ignores signal#%d altogether\n", 
		  pfx,TIMESTR(tid),FFL,
		  sig);
	}
	sl->active = -1;
      }
      else if (sig == -1) { /* Switches off ALL signals from DR_HOOK */
	int j;
	for (j=1; j<=NSIG; j++) {
	  drhook_sig_t *sl = &siglist[j];
	  if (!silent && myproc == 1) {
	    fprintf(stderr,
		    "%s %s [%s@%s:%d] DR_HOOK ignores signal#%d altogether\n", 
		    pfx,TIMESTR(tid),FFL,
		    j);
	  }
	  sl->active = -1;
	} /* for (j=1; j<=NSIG; j++) */
	break;
      }
      p = strtok(NULL,delim);
    }
    free_drhook(s);
  }
}

/*--- gdb__sigdump ---*/

#if (defined(LINUX) || defined(SUN4)) && !defined(XT3) && !defined(XD1) && !defined(_CRAYC)
static void gdb__sigdump(int sig SIG_EXTRA_ARGS)
{
  static int who = 0; /* Current owner of the lock, if > 0 */
  int is_set = 0;
  int it = get_thread_id_(); 
  drhook_sig_t *sl = &siglist[sig];
  char *pfx = PREFIX(it);

  coml_test_lockid_(&is_set, &DRHOOK_lock);
  if (is_set && who == it) {
    fprintf(stderr,"%s %s [%s@%s:%d] Received (another) signal#%d (%s)\n",
	    pfx,TIMESTR(it),FFL,
	    sig,sl->name);
    fprintf(stderr,"%s %s [%s@%s:%d] Recursive calls by the same thread#%d not allowed. Bailing out\n",
	    pfx,TIMESTR(it),FFL,
	    it);
    return;
  }
  if (!is_set) coml_set_lockid_(&DRHOOK_lock);
  who = it;
  fprintf(stderr,"%s %s [%s@%s:%d] Received signal#%d(%s) : sigcontextptr=%p\n",
	  pfx,TIMESTR(it),FFL,
	  sig,sl->name,sigcontextptr);
  LinuxTraceBack(pfx,TIMESTR(it),sigcontextptr);
  /* LinuxTraceBack(pfx,TIMESTR(tid),NULL); */
  who = 0;
  coml_unset_lockid_(&DRHOOK_lock);
}
#endif

/*--- signal_drhook ---*/

#define SETSIG5(x,ignore_flag,handler_name,preserve_old,xstr) {	\
    drhook_sig_t *sl = &siglist[x];				\
    if (sl->active == 0) {		\
      drhook_sigfunc_t u;					\
      u.func3args = handler_name;				\
      sl->active = 1;							\
      strcpy(sl->name,xstr);						\
      sigemptyset(&sl->new.sa_mask);					\
      sl->new.sa_handler = u.func1args;					\
      sl->new.sa_flags = SA_SIGINFO;					\
      sigaction(x,&sl->new,preserve_old ? &sl->old : NULL);		\
      sl->ignore_atexit = ignore_flag;					\
      trapfpe_treatment(x,silent);					\
      if (!silent && myproc == 1) {					\
	int tid = get_thread_id_();					\
	char *pfx = PREFIX(tid);					\
	const char fmt[] = "%s %s [%s@%s:%d] New signal handler '%s' for signal#%d (%s) at %p (old at %p)\n"; \
	fprintf(stderr,fmt,						\
		pfx,TIMESTR(tid),FFL,					\
		#handler_name,						\
		x, sl->name,						\
		sl->new.sa_handler,					\
		preserve_old ? sl->old.sa_handler : NULL);		\
      }							\
    }						\
}

#define SETSIG(x,ignore_flag) SETSIG5(x,ignore_flag,signal_drhook,1,#x)

#define JSETSIG(x,ignore_flag) {					\
    drhook_sig_t *sl = &siglist[x];					\
    drhook_sigfunc_t u;							\
    /* fprintf(stderr,"JSETSIG: sl->active = %d\n",sl->active); */	\
    u.func3args = signal_harakiri;					\
    sl->active = 1;							\
    strcpy(sl->name,#x);						\
    sigemptyset(&sl->new.sa_mask);					\
    sl->new.sa_handler = u.func1args;					\
    sl->new.sa_flags = SA_SIGINFO;					\
    sigaction(x,&sl->new,&sl->old);					\
    sl->ignore_atexit = ignore_flag;					\
    trapfpe_treatment(x,0);						\
  }

#if 0
    {									\
      int tid = get_thread_id_();					\
      char *pfx = PREFIX(tid);						\
      const char fmt[] = "%s %s [%s@%s:%d] Harakiri signal handler '%s' for signal#%d (%s) installed at %p (old at %p)\n"; \
      fprintf(stderr,fmt,						\
	      pfx,TIMESTR(tid),FFL,					\
	      "signal_harakiri",					\
	      x, sl->name,						\
	      sl->new.sa_handler,					\
	      sl->old.sa_handler);					\
    }									\

#endif

#if defined(RS6K) && defined(__64BIT__)
#define DRH_STRUCT_RLIMIT struct rlimit64
#define DRH_GETRLIMIT getrlimit64
#define DRH_SETRLIMIT setrlimit64
#else
#define DRH_STRUCT_RLIMIT struct rlimit
#define DRH_GETRLIMIT getrlimit
#define DRH_SETRLIMIT setrlimit
#endif

static int set_unlimited_corefile(unsigned long long int *hardlimit)
{
  /* 
     Make sure we *only* set soft-limit (not hard-limit) to 0 in our scripts i.e. :
        $ ulimit -S -c 0
     but *not*
        $ ulimit -c 0
     See man ksh or man bash for more
  */
  int rc = -1;
  if (unlimited_corefile_retcode == 9999) { /* Done only once */
    DRH_STRUCT_RLIMIT r;
    if (DRH_GETRLIMIT(RLIMIT_CORE, &r) == 0) {
      r.rlim_cur = r.rlim_max;
      if (DRH_SETRLIMIT(RLIMIT_CORE, &r) == 0) {
	saved_corefile_hardlimit = r.rlim_cur;
	rc = 0;
      }
    }
    unlimited_corefile_retcode = rc;
  }
  if (hardlimit) *hardlimit = saved_corefile_hardlimit;
  rc = unlimited_corefile_retcode;
  return rc;
}

static void 
signal_gencore(int sig SIG_EXTRA_ARGS)
{
  if (opt_gencore > 0) { 
    opt_gencore = 0; /* A tiny chance for a race condition between threads */
    if (sig == opt_gencore_signal && sig >= 1 && sig <= NSIG) {
      signal(sig, SIG_IGN);
      signal(SIGABRT, SIG_DFL);
      { /* Enable unlimited cores (up to hard-limit) and call abort() --> generates core dump */
	if (set_unlimited_corefile(NULL) == 0) {
	  int tid = get_thread_id_();
	  char *pfx = PREFIX(tid);
	  fprintf(stderr,
		  "%s %s [%s@%s:%d] Received signal#%d and now calling abort() ...\n",
		  pfx,TIMESTR(tid),FFL,
		  sig);
	  LinuxTraceBack(pfx,TIMESTR(tid),NULL);
	  abort(); /* Dump core, too */
	}
      }
      /* Should never end up here */
      fflush(NULL);
      _exit(128+ABS(sig));
    } /* if (sig >= 1 && sig <= NSIG && sig == opt_gencore_signal) */
  }
}

static char *safe_llitoa(long long int i, char b[], int blen)
{
  char const digit[] = "0123456789";
  char *p = b;
  long long int shifter;
  if (i < 0) {
    *p++ = '-';
    i *= -1;
  }
  shifter = i;
  do { /* Move to where representation ends */
    ++p;
    shifter = shifter/10;
  } while (shifter);
  *p = '\0';
  do{ /* Move back, inserting digits as u go */
    *--p = digit[i%10];
    i = i/10;
  } while (i);
  return b;
}


static void 
signal_harakiri(int sig SIG_EXTRA_ARGS)
{
  /* A signal handler that will force to exit the current thread immediately for sure */

  /* The following output should be malloc-free */

  time_t tp;
  int idummy;
  int fd = fileno(stderr);
  int tid = get_thread_id_();
  int nsigs = TIDNSIGS(tid);
  char *pfx = PREFIX(tid);
  char buf[128];
  char s[1024];
  strcpy(s,pfx);
  /* [%s@%s:%d] for FFL below */
  strcat(s," [");
  strcat(s,__FUNCTION__);
  strcat(s,"@");
  strcat(s,__FILE__);
  strcat(s,":");
  strcat(s,safe_llitoa(__LINE__,buf,sizeof(buf)));
  strcat(s,"] [epoch=");
  time(&tp);
  strcat(s,safe_llitoa(tp,buf,sizeof(buf)));
  strcat(s,"] Terminating process to avoid hangs due to signal#");
  strcat(s,safe_llitoa(sig,buf,sizeof(buf)));
  strcat(s," by raising signal SIGKILL = ");
  strcat(s,safe_llitoa(SIGKILL,buf,sizeof(buf)));
  strcat(s,", nsigs = ");
  strcat(s,safe_llitoa(nsigs,buf,sizeof(buf)));

  idummy = write(fd,s,strlen(s));

#if 0
  batch_kill_();
#endif
  
  raise(SIGKILL); /* Use raise, not RAISE here */
  _exit(128+ABS(sig)); /* Should never reach here, bu' in case it does, then ... */
}

static void 
signal_drhook(int sig SIG_EXTRA_ARGS)
{
  volatile int nfirst = drhook_use_lockfile ? 0 : 1;
  int nsigs;
  int trace_size;
  int tid;
  pid_t unixtid;
  char *pfx;
  void *trace[GNUC_BTRACE];
  // Let only one ("fastest") thread per task to this error processing
  static volatile sig_atomic_t been_here_already = 0;
  static volatile sig_atomic_t thing = 0;

  if (sig < 1 || sig > NSIG) return; // .. since have seen this, too :-(
  if (been_here_already++ > 0) return; // avoid calling more than once ... since it leads more often than not into troubles
  
  cas_lock(&thing);

  trace_size = backtrace(trace, GNUC_BTRACE);

  unixtid = gettid();

  tid = get_thread_id_();
  pfx = PREFIX(tid);

  if (signals_set && sig >= 1 && sig <= NSIG) {
    drhook_sig_t *sl = &siglist[sig];
    sigset_t newmask, oldmask;

    /* A tiny chance for a race condition between threads */
    // Using compare-and-swap -stuff from the include cas.h (also in ecProf) 

    /* Signal catching */
    {
      nsigs = (++signal_handler_called);
      if (sl->ignore_atexit) signal_handler_ignore_atexit++;
    }

    if (ec_drhook && tid >= 1 && tid <= numthreads) ec_drhook[tid-1].nsigs = nsigs; /* Store for possible signal_harakiri() */
    
    /*------------------------------------------------------------ 
      Strategy:
      - drhook intercepts most interrupts.
      - 1st interupt will 
      - call alarm(10) to try to make sure 2nd interrupt received
      - try to call tracebacks and exit (which includes atexits)
      - 2nd (and subsequent) interupts will 
      - spin for 20 sec (to give 1st interrupt time to complete tracebacks) 
      - and then call _exit (bypassing atexit)
      ------------------------------------------------------------*/
    
    /* if (sig != SIGTERM) signal(SIGTERM, SIG_DFL); */  /* Let the default SIGTERM to occur */
    
    coml_get_max_threads_(&max_threads);
    if (nsigs == 1) {
      /*---- First call to signal handler: call alarm(drhook_harakiri_timeout), tracebacks,  exit ------*/
      
      if (!nfirst) {
	const char drhook_lockfile[] = "drhook_lock";
	if (access(drhook_lockfile,F_OK) == -1) {
	  int fd = open(drhook_lockfile,O_RDONLY);
	  if (fd == -1) { // File did not exist -- create it
	    fd = open(drhook_lockfile, O_CREAT|O_WRONLY|O_TRUNC|O_EXCL, S_IRUSR|S_IWUSR);
	    if (fd >= 0) {
	      int rc_lock = flock(fd, LOCK_EX | LOCK_NB);
	      if (rc_lock == 0) {
		size_t count = sizeof(myproc);
		ssize_t sz = write(fd,&myproc,count);
		if (sz == count) nfirst = 1;
		//rc_lock = flock(fd, LOCK_UN);
	      }
	      close(fd);
	    }
	  }
	  else { // after all the file already existed
	    close(fd);
	  }
	}
      }

      if (nfirst) {
	/* Enjoy some output (only from the first guy that came in) */
	long long int hwm = gethwm_();
	long long int rss = getmaxrss_();
	long long int maxstack = getmaxstk_();
	long long int vmpeak = getvmpeak_();
	long long int pag = getpag_();
	rss /= 1048576;
	hwm /= 1048576;
	maxstack /= 1048576;
	vmpeak /= 1048576;
	fprintf(stderr,
		"%s %s [%s@%s:%d] Received signal#%d (%s) :: %lldMB (heap),"
		" %lldMB (maxrss), %lldMB (maxstack), %lldMB (vmpeak), %lld (paging), nsigs = %d\n",
		pfx,TIMESTR(tid),FFL,
		sig, sl->name, hwm, rss, maxstack, vmpeak, pag, nsigs);
#if 0
	fprintf(stderr,
		"%s %s [%s@%s:%d] Also activating Harakiri-alarm (SIGALRM=%d) to expire after %ds elapsed to prevent hangs, nsigs = %d\n",
		pfx,TIMESTR(tid),FFL,
		SIGALRM,drhook_harakiri_timeout,nsigs);
#endif
      }
      JSETSIG(SIGALRM,1); /* This will now set another signal handler than signal_drhook */
      fflush(NULL);
      alarm(drhook_harakiri_timeout);

#if defined(SA_SIGINFO) && SA_SIGINFO > 0
      if (sigcode) {
	const char *s = NULL;
	void *addr = sigcode->si_addr;
	void *bt = addr;
	ucontext_t *uc = (ucontext_t *)sigcontextptr;
#ifdef __powerpc64__
	bt = uc ? (void *) uc->uc_mcontext.regs->nip : NULL;   // Trick from PAPI_overflow()
#elif defined(__x86_64__) && defined(REG_RIP) // gcc specific
	bt = uc ? (void *) uc->uc_mcontext.gregs[REG_RIP] : NULL; // RIP: x86_64 specific ; only available in 64-bit mode */
#elif defined(__i386__) && defined(REG_EIP) // gcc specific
	bt = uc ? (void *) uc->uc_mcontext.gregs[REG_EIP] : NULL; // EIP: x86 specific ; only available in 32-bit mode */
#endif
	if (!addr) addr = bt;
	if (sig == SIGFPE) {
	  switch (sigcode->si_code) {
	  case FPE_INTDIV: s = "integer divide by zero"; break;
	  case FPE_INTOVF: s = "integer overflow"; break;
	  case FPE_FLTDIV: s = "floating-point divide by zero"; break;
	  case FPE_FLTOVF: s = "floating-point overflow"; break;
	  case FPE_FLTUND: s = "floating-point underflow"; break;
	  case FPE_FLTRES: s = "floating-point inexact result"; break;
	  case FPE_FLTINV: s = "floating-point invalid operation"; break;
	  case FPE_FLTSUB: s = "subscript out of range"; break;
	  default:
	    s = "unrecognized si_code for SIGFPE";  break;
	  }
	}
	else if (sig == SIGILL) {
	  switch (sigcode->si_code) {
	  case ILL_ILLOPC: s = "illegal opcode"; break;
	  case ILL_ILLOPN: s = "illegal operand"; break;
	  case ILL_ILLADR: s = "illegal addressing mode"; break;
	  case ILL_ILLTRP: s = "illegal trap"; break;
	  case ILL_PRVOPC: s = "privileged opcode"; break;
	  case ILL_PRVREG: s = "privileged register"; break;
	  case ILL_COPROC: s = "coprocessor error"; break;
	  case ILL_BADSTK: s = "internal stack error"; break;
	  default:
	    s = "unrecognized si_code for SIGILL";  break;
	  }
	}
	else if (sig == SIGSEGV) {
	  switch (sigcode->si_code) {
	  case SEGV_MAPERR: s = "address not mapped to object"; break;
          case SEGV_ACCERR: s = "invalid permissions for mapped object"; break;
	  default:
	    s = "unrecognized si_code for SIGSEGV";  break;
	  }
	}
	else if (sig == SIGBUS) {
	  switch (sigcode->si_code) {
	  case BUS_ADRALN: s = "invalid address alignment"; break;
          case BUS_ADRERR: s = "nonexistent physical address"; break;
	  case BUS_OBJERR: s = "object-specific hardware error"; break;
	  default:
	    s = "unrecognized si_code for SIGBUS";  break;
	  }
	}
	else {
	  s = "unrecognized si_code";
	}

	if (s) {
#ifdef __USE_GNU
	  int works = 0;
	  Dl_info dlinfo;
	  if (dladdr(bt,&dlinfo) == 0) {
	    dlinfo.dli_fname = NULL;
	    dlinfo.dli_sname = NULL;
	    dlinfo.dli_fbase = 0;
	  }
	  else
	    works = 1;

	  if (sig == SIGFPE) {
	    int excepts = fegetexcept();
	    fprintf(stderr,
		    "%s %s [%s@%s:%d] Signal#%d was caused by %s [memaddr=%p] [excepts=0x%x [%d]] : %p at %s(%s), nsigs = %d\n",
		    pfx,TIMESTR(tid),FFL,
		    sig, s, 
		    addr,
		    excepts, excepts,
		    bt,
		    dlinfo.dli_fname ? dlinfo.dli_fname : "<unknown_object>",
		    dlinfo.dli_sname ? dlinfo.dli_sname : "<unknown_function>",
		    nsigs);
	  }
	  else {
	    fprintf(stderr,
		    "%s %s [%s@%s:%d] Signal#%d was caused by %s [memaddr=%p] : %p at %s(%s), nsigs = %d\n",
		    pfx,TIMESTR(tid),FFL,
		    sig, s, 
		    addr,
		    bt,
		    dlinfo.dli_fname ? dlinfo.dli_fname : "<unknown_object>",
		    dlinfo.dli_sname ? dlinfo.dli_sname : "<unknown_function>",
		    nsigs);
	  }

	  if (works && trace_size > 0) {
	    int ndigits = (trace_size > 0) ? 1 + (int)log10(trace_size) : 0;
	    int jt;
	    for (jt = 0; jt < trace_size; ++jt) {
	      void *pbt = trace[jt];
	      if (dladdr(pbt,&dlinfo) == 0) {
		dlinfo.dli_fname = NULL;
		dlinfo.dli_sname = NULL;
		dlinfo.dli_fbase = 0;
	      }
	      fprintf(stderr,
		      "%s %s [%s@%s:%d] : [%*.*d]: %s %s %p %p # addr2line\n",
		      pfx,TIMESTR(tid),FFL,
		      ndigits, ndigits, jt,
		      dlinfo.dli_sname ? dlinfo.dli_sname : "<unknown_function>",
		      dlinfo.dli_fname ? dlinfo.dli_fname : "<unknown_object>",
		      dlinfo.dli_fbase,
		      pbt);
	    }
	  }
#else
	  fprintf(stderr,
		  "%s %s [%s@%s:%d] Signal#%d was caused by %s [memaddr=%p], nsigs = %d\n",
		  pfx,TIMESTR(tid),FFL,
		  sig, s, 
		  addr, 
		  nsigs);
#endif
	  fflush(NULL);
	}
      }
#endif

    }

    if (nsigs > 1 || !nfirst) {
      /*----- 2nd (and subsequent) calls to signal handler: spin harakiri-timeout + 60 sec,  _exit ---------*/
      int offset = 60;
      int secs = drhook_harakiri_timeout+offset;
      if (!drhook_use_lockfile) { /* Less output if lockfile was used ... */
	fprintf(stderr,
		"%s %s [%s@%s:%d] Calling signal_harakiri upon receipt of signal#%d"
		" after %ds spin, nsigs = %d, nfirst = %d\n",
		pfx,TIMESTR(tid),FFL,
		sig,secs,nsigs,nfirst);
	fflush(NULL);
      }
      spin(secs);
      signal_harakiri(sig SIG_PASS_EXTRA_ARGS);
    }

    /* All below this point should be nsigs == 1 i.e. the first threat arriving signal_drhook() */
    
#ifdef RS6K
    /*-- llcancel attempted but sometimes hangs ---
      {
      char *env = getenv("LOADL_STEP_ID");
      if (env) {
      char *cancel = "delayed_llcancel ";
      char cmd[80];
      sprintf(cmd,"%s %s &",cancel,env);
      fprintf(stderr,"tid#%d issuing command: %s\n",tid,cmd;
      fflush(NULL);
      system(cmd);
      }
      }
      ------------------------------------*/
#endif
    
    /* sigfillset(&newmask); -- dead code since sigprocmask() was not called */
    /*
      sigemptyset(&newmask);
      sigaddset(&newmask, sig);
    */
    
    /* Start critical region (we don't want any signals to interfere while doing this) */
    /* sigprocmask(SIG_BLOCK, &newmask, &oldmask); */
    
    if (nsigs == 1 && nfirst) { 
      /* Print Dr.Hook traceback */
      const int ftnunitno = 0; /* stderr */
      const int print_option = 2; /* calling tree */
      int level = 0;

      fprintf(stderr,
	      "%s %s [%s@%s:%d] Starting DrHook backtrace for signal#%d, nsigs = %d\n",
	      pfx,TIMESTR(tid),FFL,
	      sig,nsigs);

      dump_hugepages(0,pfx,tid,sig,nsigs); /* We don't wanna enforce anymore -- this the first arg == 0 now */

      if (drhook_dump_smaps) {
	char filename[64];
	snprintf(filename,sizeof(filename),"/proc/%ld/smaps",(long)unixtid);
	dump_file(pfx,tid,sig,nsigs,filename);
      }

      if (drhook_dump_maps) {
	char filename[64];
	snprintf(filename,sizeof(filename),"/proc/%ld/maps",(long)unixtid);
	dump_file(pfx,tid,sig,nsigs,filename);
      }

      if (drhook_dump_buddyinfo) {
	dump_file(pfx,tid,sig,nsigs,"/proc/buddyinfo");
      }

      if (drhook_dump_meminfo) {
	dump_file(pfx,tid,sig,nsigs,"/proc/meminfo");
      }

      fflush(NULL);

      c_drhook_print_(&ftnunitno, &tid, &print_option, &level);
      fflush(NULL);

      /* To make it less likely that another thread generates a signal while we are
	 doing a traceback lets wait a while (seems to fix problems of the traceback
	 terminating abnormally. Probably a better way of doing this involving holding
	 off signals but sigprocmask is not safe in multithreaded code -  P Towers Dec 10 2012 
	 This was originally an issue with the Intel compiler but may be of benefit for other
	 compilers. Cannot see it doing harm - P Towers Aug 29 2013 */ 
      spin(MIN(5,tid));

      if (sig != SIGABRT && sig != SIGTERM) {
#ifdef RS6K
	xl__sigdump(sig SIG_PASS_EXTRA_ARGS); /* Can't use xl__trce(...), since it also stops */
#endif
	
#if 1
	/* Active code ? */
#if (defined(LINUX) || defined(SUN4)) && !defined(XT3) && !defined(XD1)
	LinuxTraceBack(pfx,TIMESTR(tid),NULL);
#endif
#else
	/* Dead code ? */
#if (defined(LINUX) || defined(SUN4)) && !defined(XT3) && !defined(XD1) && !defined(_CRAYC)
	gdb__sigdump(sig SIG_PASS_EXTRA_ARGS);
#endif
#endif
	
#ifdef __INTEL_COMPILER
	intel_trbk_(); /* from ../utilities/gentrbk.F90 */
#endif
	
#if defined(NECSX)
	necsx_trbk_("signal_drhook",13); /* from ../utilities/gentrbk.F90 */
#endif
      }

#ifdef VPP
#if defined(SA_SIGINFO) && SA_SIGINFO > 0
      _TraceCalls(sigcontextptr); /* Need VPP's libmp.a by Pierre Lagier */
#endif
#endif
      
      fprintf(stderr, 
	      "%s %s [%s@%s:%d] DrHook backtrace done for signal#%d, nsigs = %d\n", 
	      pfx,TIMESTR(tid),FFL,
	      sig,nsigs);
      fflush(NULL);
    }
    
    /* sigprocmask(SIG_SETMASK, &oldmask, 0); */
    /* End critical region : the original signal state restored */
    
    {
      int restored = 0, tdiff;
      time_t t1, t2;
      drhook_sigfunc_t u;
      u.func3args = signal_drhook;
      if (opt_propagate_signals &&
	  sl->old.sa_handler != SIG_DFL && 
	  sl->old.sa_handler != SIG_IGN && 
	  sl->old.sa_handler != u.func1args) {
	u.func1args = sl->old.sa_handler;

	if (atp_enabled) {
	  /* Restore the default, core-file creating action to these "ATP" recognized signals */
	  switch (sig) {
	  case SIGTERM:
	    if (atp_ignore_sigterm) break; /* SIGSEGV not reset to SIG_DFL as ATP now ignores SIGTERM */
	    /* Fall thru (see man atp on Cray) */
	  case SIGINT: /* Also, see ifssig.c : used as a RESTART signal, confusingly enough */
	  case SIGFPE:
	  case SIGILL:
	  case SIGTRAP:
	  case SIGABRT:
	  case SIGBUS:
	  case SIGSEGV:
	  case SIGSYS:
	  case SIGXCPU:
#if defined(SIGXFSZ)
	  case SIGXFSZ:
#endif
	    fprintf(stderr,
		    "%s %s [%s@%s:%d] Resetting SIGSEGV (%d) to "
		    "default signal handler (SIG_DFL) before calling ATP for signal#%d, nsigs = %d\n",
		    pfx,TIMESTR(tid),FFL,
		    SIGSEGV,sig,nsigs);
	    set_default_handler(SIGSEGV,1,1);
	    restored = 1;
	    break;
	  default: 
	    break;
	  }
	}

	fprintf(stderr,
		"%s %s [%s@%s:%d] Calling previous signal handler at %p for signal#%d, nsigs = %d\n",
		pfx,TIMESTR(tid),FFL,
		u.func1args,sig,nsigs); 

	time(&t1);
	u.func3args(sig SIG_PASS_EXTRA_ARGS); /* This could now be the ATP */
	time(&t2);
	tdiff = (t2 - t1);

	fprintf(stderr,
                "%s %s [%s@%s:%d] Returned from previous signal handler"
		" (at %p, signal#%d, time taken = %ds), nsigs = %d\n",
		pfx,TIMESTR(tid),FFL,
		u.func1args,sig,tdiff,nsigs); 

	if (atp_enabled && restored && atp_max_cores > 0) {
	  /* Assuming it was indeed ATP, then lets spin a bit to allow other cores be dumped */
	  int secs = MIN(drhook_harakiri_timeout,atp_max_analysis_time);
	  int grace = 60;
	  secs = 60 + MIN(tdiff * (atp_max_cores-1),secs);
	  if (secs > 0) {
	    fprintf(stderr,
		    "%s %s [%s@%s:%d] Before aborting (signal#%d) spin %ds (incl. grace %ds)"
		    " to give ATP time to write all #%d core file(s), nsigs = %d\n",
		    pfx,TIMESTR(tid),FFL,
		    sig,secs,grace,atp_max_cores,nsigs);
	    spin(secs);
	  }
	}

	if (sig != SIGABRT && sig != SIGTERM) {
	  if (atp_enabled && atp_max_cores > 0) {
	    fprintf(stderr,
		    "%s %s [%s@%s:%d] DrHook calls abort() and attempts to dump core (signal#%d), nsigs = %d\n",
		    pfx,TIMESTR(tid),FFL,
		    sig,nsigs);
	    set_default_handler(SIGABRT,1,1);
	    abort();
	  }
	}
	/* Now proceed to definitive _exit() */
      }
      else {
	fprintf(stderr,
		"%s %s [%s@%s:%d] Not configured (DR_HOOK_PROPAGATE_SIGNALS=%d) or "
		"can't call previous signal handler (for signal#%d) in the chain at %p, nsigs = %d\n",
		pfx,TIMESTR(tid),FFL,
		opt_propagate_signals,sig,
		sl->old.sa_handler,nsigs);
      }
    }
  }
   
  {
    int errcode = 128 + ABS(sig);
    /* Make sure that the process/thread really exits now -- immediately !! */
    fprintf(stderr, "%s %s [%s@%s:%d] Error _exit(%d) upon receipt of signal#%d, nsigs = %d\n",
	    pfx,TIMESTR(tid),FFL,
	    errcode,sig,nsigs);
    fflush(NULL);
    _exit(errcode);
  }

  cas_unlock(&thing);
}

void
c_drhook_set_mpi_()
{
  dr_hook_procinfo_(&myproc, &nproc);
}

void
c_drhook_not_mpi_()
{
  /* Emulates in a one call : export DR_HOOK_NOT_MPI=1" */
  /* To have a desired effect, call BEFORE the very first call to DR_HOOK */
  static char s[] = "DR_HOOK_NOT_MPI=1"; /* note: must be static */
  putenv(s);
}


/*--- signal_drhook_init ---*/

static void 
signal_drhook_init(int enforce)
{
  char *env = getenv("DR_HOOK_SILENT");
  int silent = env ? atoi(env) : 0;
  int j;
  dr_hook_procinfo_(&myproc, &nproc);
  if (myproc < 1) myproc = 1; /* Just to enable output as if myproc was == 1 */
  /* Signals may not yet been set, since MPI not initialized 
     Only enforce-parameter can enforce to set these => no output on myproc=1 */
  if (!enforce && (myproc < 1 || nproc < 0)) return; 
  if (signals_set) return; /* Extra safety */
  /* To present sumpini.F90 (f.ex.) initializing DrHook-signals in case of 
     DR_HOOK was turned off (=0), then set also export DR_HOOK_INIT_SIGNALS=0 */
  env = getenv("DR_HOOK_INIT_SIGNALS");
  if (env && *env == '0') {
    signals_set = 2; /* Pretend they are set */
    return; /* Never initialize signals via DrHook (dangerous, but sometimes necessary) */
  }
  if (!ec_drhook) {
    int slen;
    char hostname[EC_HOST_NAME_MAX];
    char *pdot;
    int ntids = 1;
    coml_get_max_threads_(&ntids);
    numthreads = ntids;
    ec_drhook = calloc_drhook(ntids, sizeof(*ec_drhook));
    slen = sizeof(ec_drhook[0].s);
    timestr_len = sizeof(ec_drhook[0].timestr);
    if (gethostname(hostname,sizeof(hostname)) != 0) strcpy(hostname,"unknown");
    pdot = strchr(hostname,'.');
    if (pdot) *pdot = '\0'; // cut short from "." char e.g. hostname.fmi.fi becomes just "hostname"
    if (myproc == 1) {
      fprintf(stderr,"[EC_DRHOOK:hostname:myproc:omptid:pid:unixtid] [YYYYMMDD:HHMMSS:epoch:walltime] [function@file:lineno] -- Max OpenMP threads = %d\n",ntids);
    }
#if 1
    {
      extern void run_fortran_omp_parallel_ipfstr_(const int *, 
						   void (*func)(const char *, int),
						   const char *, int);
      run_fortran_omp_parallel_ipfstr_(&ntids,set_ec_drhook_label,hostname,strlen(hostname));
    }
#else
#pragma omp parallel num_threads(ntids)
    {
      set_ec_drhook_label(hostname,strlen(hostname));
    }
#endif
  }
  env = getenv("ATP_ENABLED");
  atp_enabled = (env && *env == '1') ? 1 : 0;
  if (atp_enabled) {
    env = getenv("ATP_MAX_CORES");
    if (env) atp_max_cores = atoi(env);
    env = getenv("ATP_MAX_ANALYSIS_TIME");
    if (env) atp_max_analysis_time = atoi(env);
    env = getenv("ATP_IGNORE_SIGTERM");
    if (env) atp_ignore_sigterm = atoi(env);
    if (!silent && myproc == 1) {
      int tid = get_thread_id_();
      char *pfx = PREFIX(tid);
      fprintf(stderr,"%s %s [%s@%s:%d] ATP_ENABLED=%d\n",pfx,TIMESTR(tid),FFL,atp_enabled);
      fprintf(stderr,"%s %s [%s@%s:%d] ATP_MAX_CORES=%d\n",pfx,TIMESTR(tid),FFL,atp_max_cores);
      fprintf(stderr,"%s %s [%s@%s:%d] ATP_MAX_ANALYSIS_TIME=%d\n",pfx,TIMESTR(tid),FFL,atp_max_analysis_time);
      fprintf(stderr,"%s %s [%s@%s:%d] ATP_IGNORE_SIGTERM=%d\n",pfx,TIMESTR(tid),FFL,atp_ignore_sigterm);
    }
  }
  process_options();
  for (j=1; j<=NSIG; j++) { /* Initialize */
    drhook_sig_t *sl = &siglist[j];
    sprintf(sl->name, "DR_HOOK_SIG#%d", j);
    sl->active = 0;
    sl->ignore_atexit = 0;
  }
  ignore_signals(silent); /* These signals will not be handled by DR_HOOK */
  restore_default_signals(silent); /* These signals will be restored with SIG_DFL status (regardless if to-be-caught with DrHook or ATP or anyhing else) */
  SETSIG(SIGABRT,0); /* Good to be first */
  SETSIG(SIGBUS,0);
  SETSIG(SIGSEGV,0);
#if defined(SIGEMT)
  SETSIG(SIGEMT,0);
#endif
#if defined(SIGSTKFLT)
  SETSIG(SIGSTKFLT,0); /* Stack fault */
#endif
#if !defined(NECSX)
  /* For the moment turn off these on NEC SX ... */
  SETSIG(SIGFPE,0);
  SETSIG(SIGILL,0);
#endif
  SETSIG(SIGTRAP,0); /* Should be switched off when used with debuggers */
  // SETSIG(SIGINT,0);  /* Also, see ifssig.c : used as a RESTART signal, confusingly enough */
  if (atp_enabled) {
    /* We let ATP to catch SIGQUIT (it uses this for non-failed tasks, we think) -- thus commented out */
    /* SETSIG(SIGQUIT,0); */
    /* Unless ATP ignores SIGTERM, we ignore it from DrHook -- thus conditionally commented out */
    if (atp_ignore_sigterm) SETSIG(SIGTERM,0); /* Means: DrHook does NOT ignore SIGTERM -- ATP does */
  }
  else {
    SETSIG(SIGQUIT,0);
    SETSIG(SIGTERM,0);
  }
#if defined(SIGIOT)
  SETSIG(SIGIOT,0);  /* Same as SIGABRT; Used to be a typo SIGIO ;-( */
#endif
  SETSIG(SIGXCPU,1); /* ignore_atexit == 1 i.e. no profile info via atexit() */
#if defined(SIGXFSZ)
  SETSIG(SIGXFSZ,0);
#endif
#if defined(SIGDANGER)
  SETSIG(SIGDANGER,1); /* To catch the place where paging space gets dangerously low */
#endif
  SETSIG(SIGSYS,0);
  /* SETSIG(SIGCHLD); we may not want to catch this either; may interfere parallel processing */
  /* -- not active
  SETSIG(SIGCHLD);
  SETSIG(SIGHUP);
  SETSIG(SIGCONT);
  */
#if defined(SIGCORE)
  SETSIG(SIGCORE,0); /* NEC SX core dumping */
#endif
#if defined(SIGDEAD)
  SETSIG(SIGDEAD,0); /* NEC SX dead lock */
#endif
#if defined(SIGXMEM)
  SETSIG(SIGXMEM,0); /* NEC SX exceeded memory size limit */
#endif
#if defined(SIGXDSZ)
  SETSIG(SIGXDSZ,0); /* NEC SX exceeded data size limit */
#endif
#if defined(SIGMEM32)
  SETSIG(SIGMEM32,0); /* NEC SX exceeded memory size limit of 32KB */
#endif
#if defined(SIGNMEM)
  SETSIG(SIGNMEM,0); /* NEC SX exce error for no memory */
#endif
#if defined(SIGXABT)
  SETSIG(SIGXABT,0); /* NEC SX distributed parallel program aborted */
#endif
  /*
    #if defined(SIG)
    SETSIG(SIG,0);
    #endif
  */
  catch_signals(silent); /* Additional signals to be seen by DR_HOOK */
  if (opt_gencore > 0 && opt_gencore_signal >= 1 && opt_gencore_signal <= NSIG) {
    drhook_sigfunc_t u;
    u.func3args = signal_gencore;
    signal(opt_gencore_signal, u.func1args); /* A facility to dump core */
  }
  signals_set = 1; /* Signals are set now */
}

/*--- get_mon_out ---*/

static char *
get_mon_out(int me)
{
  char *s = mon_out;
  if (mon_out_procs == me || (mon_out_procs == -1 && me >= 1 && me <= nproc)) {
    if (!mon_out) mon_out = strdup_drhook("drhook.prof.%d");
    s = malloc_drhook((strlen(mon_out) + 20) * sizeof(*s));
    sprintf(s,mon_out,me);
  }
  if (!s) s = strdup_drhook("drhook.prof.0");
  return s;
}

/*--- get_memmon_out ---*/

static char *
get_memmon_out(int me)
{
  char *s = NULL;
  char *p = get_mon_out(me);
  if (p) {
    s = malloc_drhook((strlen(p) + 5) * sizeof(*s));
    sprintf(s,"%s-mem",p);
  }
  if (!s) s = strdup_drhook("drhook.prof.0-mem");
  return s;
}

/*--- random_memstat ---*/

static void
random_memstat(int tid, int enforce)
{
  if (tid == 1 && opt_random_memstat > 0 && opt_random_memstat <= RAND_MAX) {
    int random_number = rand();
    if (enforce || random_number % opt_random_memstat == 0) {
      long long int maxhwm = getmaxhwm_();
      long long int maxstk = getmaxstk_();
      if (drhook_stacksize_threshold > 0 && maxstk > drhook_stacksize_threshold) {
	/* Abort hopefully with traceback */	
	char *pfx = PREFIX(tid);
	long long int vmpeak = getvmpeak_() / (long long int) 1048576;
	long long int threshold = drhook_stacksize_threshold / (long long int) 1048576;
	long long int ompstk = drhook_omp_stacksize / (long long int) 1048576;
	maxstk /= (long long int) 1048576;
	maxhwm /= (long long int) 1048576;
	fprintf(stderr,
		"%s %s [%s@%s:%d] Stack usage [MB] very high : %lld > %lld (= %g x OMP_STACKSIZE=%lld ; maxhwm=%lld ; vmpeak=%lld)\n",
		pfx,TIMESTR(tid),FFL,
		maxstk,threshold,
		opt_trace_stack,ompstk,
		maxhwm,vmpeak);
	RAISE(SIGABRT);
      }
    }
  }
}

/*--- process_options ---*/

static void do_prof();

void /* Fortran callable */
c_drhook_process_options_(const int *lhook, const int *Myproc, const int *Nproc)
{
    c_drhook_set_lhook_(lhook);
    if (Myproc) myproc = *Myproc;
    if (Nproc)  nproc  = *Nproc;
    process_options();
}

#define OPTPRINT(fp,...) if (fp) fprintf(fp,__VA_ARGS__)

static void
process_options()
{
  char *pfx = "";
  char *env;
  FILE *fp = NULL;
  int tid, ienv, newline;
  static int processed = 0;
  if (processed) return;

  tid = get_thread_id_();

  env = getenv("DR_HOOK_SHOW_PROCESS_OPTIONS");
  ienv = env ? atoi(env) : 1;
  if (ienv == -1 || ienv == myproc) fp = stderr;
  if (fp) pfx = PREFIX(tid);

  OPTPRINT(fp,"%s %s [%s@%s:%d] fp = %p\n",pfx,TIMESTR(tid),FFL,fp);

  env = getenv("DR_HOOK_ALLOW_COREDUMP");
  if (env) {
    ienv = atoi(env);
    allow_coredump = (ienv == -1 || ienv == myproc) ? ienv : 0;
  }
  OPTPRINT(fp,"%s %s [%s@%s:%d] DR_HOOK_ALLOW_COREDUMP=%d\n",pfx,TIMESTR(tid),FFL,allow_coredump);
  if (allow_coredump) {
    unsigned long long int hardlimit = 0;
    int rc = set_unlimited_corefile(&hardlimit);
    if (rc == 0) {
      OPTPRINT(fp,"%s %s [%s@%s:%d] Hardlimit for core file is now %llu (0x%llx)\n", 
	       pfx,TIMESTR(tid),FFL,hardlimit,hardlimit);
    }
  }

  env = getenv("DR_HOOK_PROFILE");
  if (env) {
    char *s = calloc_drhook(strlen(env) + 15, sizeof(*s));
    strcpy(s,env);
    if (!strchr(env,'%')) strcat(s,".%d");
    mon_out = strdup_drhook(s);
    free_drhook(s);
  }
  if (mon_out) OPTPRINT(fp,"%s %s [%s@%s:%d] DR_HOOK_PROFILE=%s\n",pfx,TIMESTR(tid),FFL,mon_out);

  env = getenv("DR_HOOK_PROFILE_PROC");
  if (env) {
    mon_out_procs = atoi(env);
  }
  OPTPRINT(fp,"%s %s [%s@%s:%d] DR_HOOK_PROFILE_PROC=%d\n",pfx,TIMESTR(tid),FFL,mon_out_procs);

  env = getenv("DR_HOOK_PROFILE_LIMIT");
  if (env) {
    percent_limit = atof(env);
  }
  OPTPRINT(fp,"%s %s [%s@%s:%d] DR_HOOK_PROFILE_LIMIT=%.3f\n",pfx,TIMESTR(tid),FFL,percent_limit);

  env = getenv("DR_HOOK_FUNCENTER");
  if (env) {
    opt_funcenter = atoi(env);
  }
  if (opt_funcenter) OPTPRINT(fp,"%s %s [%s@%s:%d] DR_HOOK_FUNCENTER=%d\n",pfx,TIMESTR(tid),FFL,opt_funcenter);

  env = getenv("DR_HOOK_FUNCEXIT");
  if (env) {
    opt_funcexit = atoi(env);
  }
  if (opt_funcexit) OPTPRINT(fp,"%s %s [%s@%s:%d] DR_HOOK_FUNCEXIT=%d\n",pfx,TIMESTR(tid),FFL,opt_funcexit);

  if (opt_funcenter || opt_funcexit) {
    opt_gethwm = opt_getstk = 1;
  }

  env = getenv("DR_HOOK_TIMELINE");
  if (env) {
    opt_timeline = atoi(env);
  }
  
  if (opt_timeline) {
    OPTPRINT(fp,"%s %s [%s@%s:%d] DR_HOOK_TIMELINE=%d\n",pfx,TIMESTR(tid),FFL,opt_timeline);

    env = getenv("DR_HOOK_TIMELINE_THREAD");
    if (env) {
      opt_timeline_thread = atoi(env);
    }
    OPTPRINT(fp,"%s %s [%s@%s:%d] DR_HOOK_TIMELINE_THREAD=%d\n",pfx,TIMESTR(tid),FFL,opt_timeline_thread);
    
    env = getenv("DR_HOOK_TIMELINE_FORMAT");
    if (env) {
      opt_timeline_format = atoi(env);
    }
    OPTPRINT(fp,"%s %s [%s@%s:%d] DR_HOOK_TIMELINE_FORMAT=%d\n",pfx,TIMESTR(tid),FFL,opt_timeline_format);
    
    env = getenv("DR_HOOK_TIMELINE_UNITNO");
    if (env) {
      opt_timeline_unitno = atoi(env);
    }
    OPTPRINT(fp,"%s %s [%s@%s:%d] DR_HOOK_TIMELINE_UNITNO=%d\n",pfx,TIMESTR(tid),FFL,opt_timeline_unitno);

    env = getenv("DR_HOOK_TIMELINE_FREQ");
    if (env) {
      opt_timeline_freq = atoi(env);
    }
    OPTPRINT(fp,"%s %s [%s@%s:%d] DR_HOOK_TIMELINE_FREQ=%lld\n",pfx,TIMESTR(tid),FFL,opt_timeline_freq);

    env = getenv("DR_HOOK_TIMELINE_MB");
    if (env) {
      opt_timeline_MB = atof(env);
      if (opt_timeline_MB < 0) opt_timeline_MB = 1.0;
    }
    OPTPRINT(fp,"%s %s [%s@%s:%d] DR_HOOK_TIMELINE_MB=%g\n",pfx,TIMESTR(tid),FFL,opt_timeline_MB);
  }

  if (myproc == 1) { /* Only applicable for master MPI task for now */
    env = getenv("DR_HOOK_TRACE_STACK");
    if (env) {
      opt_trace_stack = atof(env);
      if (opt_trace_stack < 0) 
	opt_trace_stack = 0;
      else {
	drhook_omp_stacksize = slave_stacksize();
	if (drhook_omp_stacksize > 0) {
	  drhook_stacksize_threshold = opt_trace_stack * drhook_omp_stacksize;
	  opt_random_memstat = 1;
	  random_memstat(1,1);
	  OPTPRINT(fp,"%s %s [%s@%s:%d] DR_HOOK_TRACE_STACK=%g\n",pfx,TIMESTR(tid),FFL,opt_trace_stack);
	}
	else
	  opt_trace_stack = 0;
      }
    }
  }

  if (!opt_random_memstat) {
    env = getenv("DR_HOOK_RANDOM_MEMSTAT");
    if (env) {
      opt_random_memstat = atoi(env);
      if (opt_random_memstat < 0) opt_random_memstat = 0;
      if (opt_random_memstat > RAND_MAX) opt_random_memstat = RAND_MAX;
      random_memstat(1,1);
    }
  }

  OPTPRINT(fp,"%s %s [%s@%s:%d] DR_HOOK_RANDOM_MEMSTAT=%d  (RAND_MAX=%d)\n",pfx,TIMESTR(tid),FFL,opt_random_memstat,RAND_MAX);
    
  env = getenv("DR_HOOK_HASHBITS");
  if (env) {
    int value = atoi(env);
    if (value < 1) value = 1;
    else if (value > NHASHMAX) value = NHASHMAX;
    nhash = value;
    hashsize = HASHSIZE(nhash);
    hashmask = HASHMASK(nhash);
  }
  OPTPRINT(fp,"%s %s [%s@%s:%d] DR_HOOK_HASHBITS=%d\n",pfx,TIMESTR(tid),FFL,nhash);

  env = getenv("DR_HOOK_NCALLSTACK");
  if (env) {
    int value = atoi(env);
    if (value < 1) value = NCALLSTACK;
    cstklen = value;
  }
  OPTPRINT(fp,"%s %s [%s@%s:%d] DR_HOOK_NCALLSTACK=%d\n",pfx,TIMESTR(tid),FFL,cstklen);

  env = getenv("DR_HOOK_HARAKIRI_TIMEOUT");
  if (env) {
    int value = atoi(env);
    if (value < 1) value = drhook_harakiri_timeout_default;
    drhook_harakiri_timeout = value;
  }
  OPTPRINT(fp,"%s %s [%s@%s:%d] DR_HOOK_HARAKIRI_TIMEOUT=%d\n",pfx,TIMESTR(tid),FFL,drhook_harakiri_timeout);

  env = getenv("DR_HOOK_USE_LOCKFILE");
  if (env) {
    int value = atoi(env);
    drhook_use_lockfile = (value != 0) ? 1 : 0; /* currently accept just 0 or 1 */
  }
  OPTPRINT(fp,"%s %s [%s@%s:%d] DR_HOOK_USE_LOCKFILE=%d\n",pfx,TIMESTR(tid),FFL,drhook_use_lockfile);

  env = getenv("DR_HOOK_TRAPFPE");
  if (env) {
    int value = atoi(env);
    drhook_trapfpe = (value != 0) ? 1 : 0; /* currently accept just 0 or 1 */
  }
  OPTPRINT(fp,"%s %s [%s@%s:%d] DR_HOOK_TRAPFPE=%d\n",pfx,TIMESTR(tid),FFL,drhook_trapfpe);

  env = getenv("DR_HOOK_TRAPFPE_INVALID");
  if (env) {
    int value = atoi(env);
    drhook_trapfpe_invalid = (value != 0) ? 1 : 0; /* currently accept just 0 or 1 */
  }
  OPTPRINT(fp,"%s %s [%s@%s:%d] DR_HOOK_TRAPFPE_INVALID=%d\n",pfx,TIMESTR(tid),FFL,drhook_trapfpe_invalid);

  env = getenv("DR_HOOK_TRAPFPE_DIVBYZERO");
  if (env) {
    int value = atoi(env);
    drhook_trapfpe_divbyzero = (value != 0) ? 1 : 0; /* currently accept just 0 or 1 */
  }
  OPTPRINT(fp,"%s %s [%s@%s:%d] DR_HOOK_TRAPFPE_DIVBYZERO=%d\n",pfx,TIMESTR(tid),FFL,drhook_trapfpe_divbyzero);

  env = getenv("DR_HOOK_TRAPFPE_OVERFLOW");
  if (env) {
    int value = atoi(env);
    drhook_trapfpe_overflow = (value != 0) ? 1 : 0; /* currently accept just 0 or 1 */
  }
  OPTPRINT(fp,"%s %s [%s@%s:%d] DR_HOOK_TRAPFPE_OVERFLOW=%d\n",pfx,TIMESTR(tid),FFL,drhook_trapfpe_overflow);


  env = getenv("DR_HOOK_TIMED_KILL");
  if (env) {
    drhook_timed_kill = strdup_drhook(env);
  }
  if (drhook_timed_kill) OPTPRINT(fp,"%s %s [%s@%s:%d] DR_HOOK_TIMED_KILL=%s\n",pfx,TIMESTR(tid),FFL,drhook_timed_kill);

  env = getenv("DR_HOOK_DUMP_SMAPS");
  if (env) {
    ienv = atoi(env);
    drhook_dump_smaps = (ienv != 0) ? 1 : 0;
  }
  if (drhook_dump_smaps) OPTPRINT(fp,"%s %s [%s@%s:%d] DR_HOOK_DUMP_SMAPS=%d\n",pfx,TIMESTR(tid),FFL,drhook_dump_smaps);

  env = getenv("DR_HOOK_DUMP_MAPS");
  if (env) {
    ienv = atoi(env);
    drhook_dump_maps = (ienv != 0) ? 1 : 0;
  }
  if (drhook_dump_maps) OPTPRINT(fp,"%s %s [%s@%s:%d] DR_HOOK_DUMP_MAPS=%d\n",pfx,TIMESTR(tid),FFL,drhook_dump_maps);

  env = getenv("DR_HOOK_DUMP_BUDDYINFO");
  if (env) {
    ienv = atoi(env);
    drhook_dump_buddyinfo = (ienv != 0) ? 1 : 0;
  }
  if (drhook_dump_buddyinfo) OPTPRINT(fp,"%s %s [%s@%s:%d] DR_HOOK_DUMP_BUDDYINFO=%d\n",pfx,TIMESTR(tid),FFL,drhook_dump_buddyinfo);

  env = getenv("DR_HOOK_DUMP_MEMINFO");
  if (env) {
    ienv = atoi(env);
    drhook_dump_meminfo = (ienv != 0) ? 1 : 0;
  }
  if (drhook_dump_meminfo) OPTPRINT(fp,"%s %s [%s@%s:%d] DR_HOOK_DUMP_MEMINFO=%d\n",pfx,TIMESTR(tid),FFL,drhook_dump_meminfo);

  env = getenv("DR_HOOK_DUMP_HUGEPAGES");
  if (env) {
    double freq;
    int nel = sscanf(env,"%d,%lf",&ienv,&freq);
    if (nel == 2) {
      drhook_dump_hugepages = (freq > 0 && (ienv == -1 || ienv == myproc)) ? ienv : 0;
      if (drhook_dump_hugepages) drhook_dump_hugepages_freq = freq;
    }
  }
  if (drhook_dump_hugepages) OPTPRINT(fp,"%s %s [%s@%s:%d] DR_HOOK_DUMP_HUGEPAGES=%d,%.6f\n",pfx,TIMESTR(tid),FFL,
				      drhook_dump_hugepages,drhook_dump_hugepages_freq);

  env = getenv("DR_HOOK_GENCORE");
  if (env) {
    opt_gencore = atoi(env);
  }

  if (opt_gencore) {
    OPTPRINT(fp,"%s %s [%s@%s:%d] DR_HOOK_GENCORE=%d\n",pfx,TIMESTR(tid),FFL,opt_gencore);
    
    env = getenv("DR_HOOK_GENCORE_SIGNAL");
    if (env) {
      int itmp = atoi(env);
      if (itmp >= 1 && itmp <= NSIG && itmp != SIGABRT) {
	opt_gencore_signal = itmp;
      }
    }
    OPTPRINT(fp,"%s %s [%s@%s:%d] DR_HOOK_GENCORE_SIGNAL=%d\n",pfx,TIMESTR(tid),FFL,opt_gencore_signal);
  }

  env = getenv("DR_HOOK_HPMSTOP");
  if (env) {
    char *s = strdup_drhook(env);
    long long int a;
    double b;
    int n = 0;
    env = s;
    while (*env) {
      if (isspace(*env) || *env == ',') *env = ' ';
      env++;
    }
    n = sscanf(s,"%lld %lf",&a,&b);
    if (n >= 1) opt_hpmstop_threshold = a;
    if (n >= 2) opt_hpmstop_mflops = b;
    OPTPRINT(fp,"%s %s [%s@%s:%d] DR_HOOK_HPMSTOP=%lld,%.15g\n",
		    pfx,TIMESTR(tid),FFL,opt_hpmstop_threshold,opt_hpmstop_mflops);
    free_drhook(s);
  }

  newline = 0;
  env = getenv("DR_HOOK_OPT");
  if (env) {
    const char delim[] = ", \t/";
    char *comma = " DR_HOOK_OPT=\"";
    char *s = strdup_drhook(env);
    char *p = s;
    while (*p) {
      if (islower(*p)) *p = toupper(*p);
      p++;
    } 
    p = strtok(s,delim);
    /* if (p) OPTPRINT(fp,"%s %s [%s@%s:%d] DR_HOOK_OPT=\"",pfx,TIMESTR(tid)); */
    if (p && fp) {
      fprintf(fp,"%s %s [%s@%s:%d]",pfx,TIMESTR(tid),FFL);
      newline = 1;
    }
    while (p) {
      /* Assume that everything is OFF by default */
      if (strequ(p,"ALL")) { /* all except profiler data */
	opt_gethwm = opt_getstk = opt_getrss = opt_getpag = opt_walltime = opt_cputime = 1;
	opt_calls = 1;
	any_memstat++;
	OPTPRINT(fp,"%s%s",comma,"ALL"); comma = ",";
      }
      else if (strequ(p,"MEM") || strequ(p,"MEMORY")) {
	opt_gethwm = opt_getstk = opt_getrss = 1;
	opt_calls = 1;
	any_memstat++;
	OPTPRINT(fp,"%s%s",comma,"MEMORY"); comma = ",";
      }
      else if (strequ(p,"TIME") || strequ(p,"TIMES")) {
	opt_walltime = opt_cputime = 1;
	opt_calls = 1;
	OPTPRINT(fp,"%s%s",comma,"TIMES"); comma = ",";
      }
      else if (strequ(p,"HWM") || strequ(p,"HEAP")) {
	opt_gethwm = 1;
	opt_calls = 1;
	any_memstat++;
	OPTPRINT(fp,"%s%s",comma,"HEAP"); comma = ",";
      }
      else if (strequ(p,"STK") || strequ(p,"STACK")) {
	opt_getstk = 1;
	opt_calls = 1;
	any_memstat++;
	OPTPRINT(fp,"%s%s",comma,"STACK"); comma = ",";
      }
      else if (strequ(p,"RSS")) {
	opt_getrss = 1;
	opt_calls = 1;
	any_memstat++;
	OPTPRINT(fp,"%s%s",comma,"RSS"); comma = ",";
      }
      else if (strequ(p,"PAG") || strequ(p,"PAGING")) {
	opt_getpag = 1;
	opt_calls = 1;
	any_memstat++;
	OPTPRINT(fp,"%s%s",comma,"PAGING"); comma = ",";
      }
      else if (strequ(p,"WALL") || strequ(p,"WALLTIME")) {
	opt_walltime = 1;
	opt_calls = 1;
	OPTPRINT(fp,"%s%s",comma,"WALLTIME"); comma = ",";
      }
      else if (strequ(p,"CPU") || strequ(p,"CPUTIME")) {
	opt_cputime = 1;
	opt_calls = 1;
	OPTPRINT(fp,"%s%s",comma,"CPUTIME"); comma = ",";
      }
      else if (strequ(p,"CALLS") || strequ(p,"COUNT")) {
	opt_calls = 1;
	OPTPRINT(fp,"%s%s",comma,"CALLS"); comma = ",";
      }
      else if (strequ(p,"MEMPROF")) {
	opt_memprof = 1;
	drhook_memtrace = 1;
	opt_gethwm = opt_getstk = opt_getrss = 1;
	opt_getpag = 1;
	opt_calls = 1;
	any_memstat++;
	OPTPRINT(fp,"%s%s",comma,"MEMPROF"); comma = ",";
      }
      else if (strequ(p,"PROF") || strequ(p,"WALLPROF")) {
	opt_wallprof = 1;
	opt_walltime = 1;
	opt_cpuprof = 0; /* Note: Switches cpuprof OFF */
	opt_calls = 1;
	OPTPRINT(fp,"%s%s",comma,"WALLPROF"); comma = ",";
      }
      else if (strequ(p,"CPUPROF")) {
	opt_cpuprof = 1;
	opt_cputime = 1;
	opt_wallprof = 0; /* Note: Switches walprof OFF */
	opt_calls = 1;
	OPTPRINT(fp,"%s%s",comma,"CPUPROF"); comma = ",";
      }
      else if (strequ(p,"HPM") || strequ(p,"HPMPROF") || strequ(p,"MFLOPS")) {
	opt_hpmprof = 1;
	opt_wallprof = 1; /* Note: Implies wallprof (or prof), not cpuprof */
	opt_walltime = 1;
	opt_cpuprof = 0;  /* Note: Switches cpuprof OFF */
	opt_calls = 1;
	OPTPRINT(fp,"%s%s",comma,"HPMPROF"); comma = ",";
      }
      else if (strequ(p,"TRIM")) {
	opt_trim = 1;
	OPTPRINT(fp,"%s%s",comma,"TRIM"); comma = ",";
      }
      else if (strequ(p,"SELF")) {
	opt_self = 2;
	OPTPRINT(fp,"%s%s",comma,"SELF"); comma = ",";
      }
      else if (strequ(p,"NOSELF")) {
	opt_self = 0;
	OPTPRINT(fp,"%s%s",comma,"NOSELF"); comma = ",";
      }
      else if (strequ(p,"NOPROP") || strequ(p,"NOPROPAGATE") ||
	       strequ(p,"NOPROPAGATE_SIGNALS")) {
	opt_propagate_signals = 0;
	OPTPRINT(fp,"%s%s",comma,"NOPROPAGATE_SIGNALS"); comma = ",";
      }
      else if (strequ(p,"NOSIZE") || strequ(p,"NOSIZEINFO")) {
	opt_sizeinfo = 0;
	OPTPRINT(fp,"%s%s",comma,"NOSIZEINFO"); comma = ",";
      }
      else if (strequ(p,"CLUSTER") || strequ(p,"CLUSTERINFO")) {
	opt_clusterinfo = 1;
	OPTPRINT(fp,"%s%s",comma,"CLUSTERINFO"); comma = ",";
      }
      else if (strequ(p,"CALLPATH")) {
	opt_callpath = 1;
	OPTPRINT(fp,"%s%s",comma,"CALLPATH"); comma = ",";
      }
      p = strtok(NULL,delim);
    }
    free_drhook(s);
    if (*comma == ',') {
      OPTPRINT(fp,"\"\n");
      newline = 0;
    }
    if (newline) OPTPRINT(fp,"\n");

    if (opt_callpath) {
      env = getenv("DR_HOOK_CALLPATH_INDENT");
      if (env) {
	callpath_indent = atoi(env);
	if (callpath_indent < 1 || callpath_indent > 8) callpath_indent = callpath_indent_default;
      }
      OPTPRINT(fp,"%s %s [%s@%s:%d] DR_HOOK_CALLPATH_INDENT=%d\n",pfx,TIMESTR(tid),FFL,callpath_indent);
      
      env = getenv("DR_HOOK_CALLPATH_DEPTH");
      if (env) {
	callpath_depth = atoi(env);
	if (callpath_depth < 0) callpath_depth = callpath_depth_default;
      }
      OPTPRINT(fp,"%s %s [%s@%s:%d] DR_HOOK_CALLPATH_DEPTH=%d\n",pfx,TIMESTR(tid),FFL,callpath_depth);
      
      env = getenv("DR_HOOK_CALLPATH_PACKED");
      if (env) {
	callpath_packed = atoi(env);
      }
      OPTPRINT(fp,"%s %s [%s@%s:%d] DR_HOOK_CALLPATH_PACKED=%d\n",pfx,TIMESTR(tid),FFL,callpath_packed);
      
      env = getenv("DR_HOOK_CALLTRACE");
      if (env) {
	opt_calltrace = atoi(env);
      }
      OPTPRINT(fp,"%s %s [%s@%s:%d] DR_HOOK_CALLTRACE=%d\n",pfx,TIMESTR(tid),FFL,opt_calltrace);
    }

    if (opt_wallprof || opt_cpuprof || opt_memprof || opt_timeline) {
      atexit(do_prof);
    }
  }
  else {
    if (opt_timeline) atexit(do_prof);
  } /* if (env) */

  processed = 1;
}

/*--- trim ---*/

static const char *
trim(const char *name, int *n)
{
  const char *from;
  int len;
  int name_len = *n;
  while (*name && isspace(*name) && name_len > 0) {
    /* skip leading blanks */
    name++;
    name_len--;
  }
  len = 0;
  from = name;
  while (*from && !isspace(*from) && name_len > 0) {
    /* find first space point, if any */
    from++;
    len++;
    name_len--;
  }
  *n = len;
  if (!name) {
    /* Never actually called (unless a true fatality) */
    ABOR1("***Fatal error in drhook.c:trim()-function");
  }
  return name;
}

/*--- insertkey ---*/

static drhook_key_t *
insertkey(int tid, const drhook_key_t *keyptr_in)
{
  drhook_key_t *keyptr = NULL;
  if (tid >= 1 && tid <= numthreads) {
    /* no trimming available for this; just raw eval & insert */
    unsigned int hash = hashfunc(keyptr_in->name, keyptr_in->name_len);
    keyptr = &keydata[tid-1][hash];
    for (;;) {
      if (!keyptr->name) { /* A free slot */
	memcpy(keyptr,keyptr_in,sizeof(*keyptr));
	keyptr->next = NULL;
	break;
      }
      else {
	if (!keyptr->next) {
	  keyptr->next = calloc_drhook(1, sizeof(drhook_key_t)); /* chaining */
	}
	keyptr = keyptr->next;
      }  /* if (!keyptr->name) ... else ... */
    } /* for (;;) */
  } /* if (tid >= 1 && tid <= numthreads) */
  return keyptr;
}

/*--- getkey ---*/

static drhook_key_t *
getkey(int tid, const char *name, int name_len,
       const char *filename, int filename_len,
       const double *walltime, const double *cputime,
       const equivalence_t *callpath, int callpath_len,
       int *free_callpath)
{
  drhook_key_t *keyptr = NULL;
  if (tid >= 1 && tid <= numthreads) {
    unsigned int hash, fullhash;
    if (opt_trim) name = trim(name, &name_len);
    hash = hashfunc(name, name_len);
    if (callpath) {
      callpath_hashfunc(hash, callpath, callpath_len, &fullhash);
#ifdef DEBUG
      fprintf(stderr,
	      "getkey: name='%.*s', name_len=%d, callpath_len=%d, fullhash=%u\n",
	      name_len, name, name_len, callpath_len, fullhash);
#endif
    }
    keyptr = &keydata[tid-1][hash];
    for (;;) {
      int found = 0;
      if (!keyptr->name) { /* A free slot */
	keyptr->name = malloc_drhook((name_len+1)*sizeof(*name));
	keyptr->name_len = name_len;
	if (opt_trim) {
	  const char *from = name;
	  char *to = keyptr->name;
	  int len = name_len;
	  for (; len>0; from++, len--) {
	    *to++ = islower(*from) ? toupper(*from) : *from;
	  }
	  *to = 0;
	}
	else {
	  memcpy(keyptr->name, name, name_len);
	  keyptr->name[name_len] = 0;
	}
	if (filename_len > 0 &&
	    filename && 
	    *filename) {
	  char *psave = NULL;
	  char *p = psave = malloc_drhook((filename_len+1)*sizeof(*filename));
	  memcpy(p, filename, filename_len);
	  p[filename_len] = 0;
	  { /* Strip out dirname */
	    char *s = strrchr(p,'/');
	    if (s) p = s+1;
	  }
	  keyptr->filename = strdup_drhook(p);
	  free_drhook(psave);
	}
	if (callpath) {
	  if (free_callpath) *free_callpath = 0;
	  keyptr->callpath = callpath;
	  keyptr->callpath_len = callpath_len;
	  keyptr->callpath_fullhash = fullhash;
	}
	found = 1;
      }
      if (found || 
	  (keyptr->name_len == name_len &&
	   (!callpath || (callpath && keyptr->callpath && 
			  keyptr->callpath_len == callpath_len &&
			  keyptr->callpath_fullhash == fullhash)) &&
	   ((!opt_trim && *keyptr->name == *name && strnequ(keyptr->name, name, name_len)) ||
	    (opt_trim && strncasecmp(keyptr->name, name, name_len) == 0)))) {
	if (opt_walltime) keyptr->wall_in = walltime ? *walltime : WALLTIME();
	if (opt_cputime) keyptr->cpu_in  = cputime ? *cputime : CPUTIME();
	if (any_memstat) memstat(keyptr,&tid,1);
	if (opt_calls) {
	  keyptr->calls++;
	  keyptr->status++;
	}
	insert_calltree(tid, keyptr);
	break; /* for (;;) */
      }
      else {
	if (!keyptr->next) {
	  keyptr->next = calloc_drhook(1, sizeof(drhook_key_t)); /* chaining */
	}
	keyptr = keyptr->next;
      }  /* if (found ...) else ... */
    } /* for (;;) */
    curkeyptr[tid-1] = keyptr;
  } /* if (tid >= 1 && tid <= numthreads) */
  return keyptr;
}

/*--- putkey ---*/

static void
putkey(int tid, drhook_key_t *keyptr, const char *name, int name_len,
       int sizeinfo,
       double *walltime, double *cputime)
{
  const int sig = SIGABRT;
  const char sl_name[] = "SIGABRT";
  drhook_calltree_t *treeptr = (tid >= 1 && tid <= numthreads) ? thiscall[tid-1] : NULL;
  if (!treeptr || !treeptr->active || treeptr->keyptr != keyptr) {
    char *pfx = PREFIX(tid);
    char *s;
    unsigned int hash;
    if (opt_trim) name = trim(name, &name_len);
    hash = hashfunc(name, name_len);
    s = strdup2_drhook(name,name_len);
    if (opt_trim) {
      char *p = s;
      while (*p) {
	if (islower(*p)) *p = toupper(*p);
	p++;
      }
    }
    fprintf(stderr,
	    "%s %s [%s@%s:%d] [signal#%d(%s)]: Dr.Hook has detected an invalid"
	    " key-pointer/handle while leaving the routine '%s' [hash=%u]\n",
	    pfx,TIMESTR(tid),FFL,
	    sig,sl_name,s,hash);

    if (treeptr) {
      equivalence_t u;

      u.keyptr = treeptr->keyptr;
      hash = (u.keyptr && u.keyptr->name) ? hashfunc(u.keyptr->name,u.keyptr->name_len) : 0;
      fprintf(stderr,
	      "%s %s [%s@%s:%d] [signal#%d(%s)]: Expecting the key-pointer=%p"
	      " and treeptr->active-flag = 1\n",
	      pfx,TIMESTR(tid),FFL,
	      sig,sl_name,u.keyptr);
      fprintf(stderr,
	      "%s %s [%s@%s:%d] [signal#%d(%s)]: A probable routine missing the closing"
	      " DR_HOOK-call is '%s' [hash=%u]\n",
	      pfx,TIMESTR(tid),FFL,
	      sig,sl_name,
	      (u.keyptr && u.keyptr->name) ? u.keyptr->name : NIL, hash);

      u.keyptr = keyptr;
      hash = (u.keyptr && u.keyptr->name) ? hashfunc(u.keyptr->name,u.keyptr->name_len) : 0;
      fprintf(stderr,
	      "%s %s [%s@%s:%d] [signal#%d(%s)]: Got a key-pointer=%p"
	      " and treeptr->active-flag = %d\n",
	      pfx,TIMESTR(tid),FFL,
	      sig,sl_name,u.keyptr,treeptr->active);
      fprintf(stderr,
	      "%s %s [%s@%s:%d] [signal#%d(%s)]: This key-pointer maybe associated with"
	      " the routine '%s' [hash=%u]\n",
	      pfx,TIMESTR(tid),FFL,
	      sig,sl_name,
	      (u.keyptr && u.keyptr->name) ? u.keyptr->name : NIL, hash);

      u.keyptr = curkeyptr[tid-1];
      hash = (u.keyptr && u.keyptr->name) ? hashfunc(u.keyptr->name,u.keyptr->name_len) : 0;
      fprintf(stderr,
	      "%s %s [%s@%s:%d] [signal#%d(%s)]: The current key-pointer (=%p) thinks"
	      " it maybe associated with the routine '%s' [hash=%u]\n",
	      pfx,TIMESTR(tid),FFL,
	      sig,sl_name,
	      u.keyptr,
	      (u.keyptr && u.keyptr->name) ? u.keyptr->name : NIL, hash);
    }
    free_drhook(s);
    fprintf(stderr,
	    "%s %s [%s@%s:%d] [signal#%d(%s)]: Aborting...\n",
	    pfx,TIMESTR(tid),FFL,
	    sig,sl_name);
    RAISE(SIGABRT);
  }
  else if (tid >= 1 && tid <= numthreads) {
    double delta_wall = 0;
    double delta_cpu  = 0;
    if (any_memstat) memstat(keyptr,&tid,0);
    if (opt_calls)   keyptr->status--;
    if (opt_sizeinfo && sizeinfo > 0) {
      if (keyptr->sizeinfo == 0) { /* First time */
	keyptr->min_sizeinfo = sizeinfo;
	keyptr->max_sizeinfo = sizeinfo;
      }
      else {
	keyptr->min_sizeinfo = MIN(keyptr->min_sizeinfo, sizeinfo);
	keyptr->max_sizeinfo = MAX(keyptr->max_sizeinfo, sizeinfo);
      }
      keyptr->sizeinfo += sizeinfo;
    }
    if (opt_cputime && cputime) {
      *cputime = CPUTIME();
      delta_cpu = *cputime - keyptr->cpu_in;
    }
    if (opt_walltime && walltime) {
      *walltime = WALLTIME();
      delta_wall = *walltime - keyptr->wall_in;
    }
    if (opt_walltime) keyptr->delta_wall_all += delta_wall;
    if (opt_cputime)  keyptr->delta_cpu_all  += delta_cpu;
    remove_calltree(tid, keyptr, &delta_wall, &delta_cpu);
  }
}
    
/*--- init_drhook ---*/

static void
init_drhook(int ntids)
{
  if (numthreads == 0 || !keydata || !calltree || !keyself || !overhead || !curkeyptr || !cstk) {
    int j;
    if (pid == -1) { /* Ensure that just called once */
      {
	/* Invoke once : timers, memory counters etc. to "wake them up" */
	(void) WALLTIME();
	(void) CPUTIME();
	(void) gethwm_();
	(void) getmaxhwm_();
	(void) getrss_();
	(void) getmaxrss_();
	(void) getstk_();
	(void) getmaxstk_();
	(void) getpag_();
      }
#ifdef RS6K
      irtc_start = irtc();
#endif
#ifdef CRAYXT
      dclock_start = dclock();
#endif
#if defined(SV2) || defined(XD1) || defined(XT3)
#if defined(SV2)
      irtc_start = _rtc();
#else
      irtc_start = irtc_();
#endif
      my_irtc_rate = irtc_rate_();
      my_inv_irtc_rate = 1.0/my_irtc_rate;
#endif
      start_stamp = timestamp();
      {
	char *env = getenv("DR_HOOK_SHOW_LOCK"); /* export DR_HOOK_SHOW_LOCK=1 to show the lock-info */
	int konoff = env ? atoi(env) : 0;
	int kret = 0;
	if (konoff == 1) coml_set_debug_(&konoff, &kret);
	INIT_LOCKID_WITH_NAME(&DRHOOK_lock,"drhook.c:DRHOOK_lock");
	if (kret != 0) {
	  konoff = 0;
	  coml_set_debug_(&konoff, &kret);
	}
      }
#if defined(NECSX)
      { /* If C-programs compiled with -traceback, then NEC/F90 
	   MESPUT-call will also includes C-routines in the traceback if 
	   in addition 'export C_TRACEBACK=YES' */
	char *env = getenv("C_TRACEBACK");
	if (!env) {
	  /* Override only if C_TRACEBACK hadn't already been defined */
	  static char s[] = "C_TRACEBACK=YES"; /* note: must be static */
	  putenv(s);
	}
      }
#endif
      ec_set_umask_();
      pid = getpid();
      signal_drhook_init(1); /* myproc gets set .. if not earlier */
      process_options();
      set_timed_kill();
      drhook_lhook = 1;
    }
    if (!keydata) {
      keydata = malloc_drhook(sizeof(**keydata) * ntids);
      for (j=0; j<ntids; j++) {
	keydata[j] = calloc_drhook(hashsize, sizeof(drhook_key_t));
      }
    }
    if (!cstk) {
      cstk = calloc_drhook(ntids, sizeof(**cstk));
    }
    if (!calltree) {
      calltree = malloc_drhook(sizeof(**calltree) * ntids);
      thiscall = malloc_drhook(sizeof(**thiscall) * ntids);
      for (j=0; j<ntids; j++) {
	thiscall[j] = calltree[j] = calloc_drhook(1,sizeof(drhook_calltree_t));
      }
    }
    if (!keyself && opt_self && (opt_wallprof || opt_cpuprof || opt_hpmprof)) {
      const char *name = "$drhook";
      int name_len = strlen(name);
      keyself = malloc_drhook(sizeof(**keyself) * ntids);
      for (j=0; j<ntids; j++) {
	drhook_key_t *keyptr = keyself[j] = calloc_drhook(1,sizeof(drhook_key_t));
	keyptr->name = strdup_drhook(name);
	keyptr->name_len = name_len;
      }
    }
    if (!overhead) {
      overhead = calloc_drhook(ntids,sizeof(*overhead));
    }
    if (!curkeyptr) {
      curkeyptr = malloc_drhook(sizeof(**curkeyptr) * ntids);
      for (j=0; j<ntids; j++) {
	curkeyptr[j] = NULL;
      }
    }
    numthreads = ntids;
    if (!timeline) {
      if (opt_timeline_unitno >= 0 && opt_timeline_freq >= 1 &&
	  (opt_timeline == myproc || opt_timeline == -1)) {
	timeline = calloc_drhook(ntids, sizeof(*timeline));
      }
      if (timeline) drhook_memtrace = 1;
      if (timeline) {
	/* The first timeline-call */
	const int ftnunitno = opt_timeline_unitno;
	const int master = 1;
	const int print_option = +7;
	int initlev = 0;
	c_drhook_print_(&ftnunitno, &master, &print_option, &initlev);
      }
    }
    init_hpm(1); /* First thread */
  }
}

/*-- overhead-macro --*/

#define OVERHEAD(tid,walltime_in,cputime_in,delta,calc_delta) \
if (overhead && tid >= 1 && tid <= numthreads) { \
  if (calc_delta) { \
    if      (opt_walltime) delta = WALLTIME() - walltime_in; \
    else if (opt_cputime)  delta = CPUTIME()  - cputime_in; \
    else                   delta = 0; \
  } \
  overhead[tid-1] += delta; \
}

/*--- itself ---*/

#define ITSELF_0 \
  double delta = 0; \
  drhook_key_t *keyptr_self = keyself ? itself(NULL,*thread_id,0,NULL,&walltime,&cputime) : NULL;

#define ITSELF_1 \
  if (keyptr_self) { \
    (void) itself(keyptr_self,*thread_id,1,&delta,&walltime,&cputime); \
    if (opt_wallprof) u.keyptr->delta_wall_child += delta; \
    else              u.keyptr->delta_cpu_child  += delta; \
    OVERHEAD(*thread_id,walltime,cputime,delta,0); \
  } \
  else { \
    OVERHEAD(*thread_id,walltime,cputime,delta,1); \
  }

static drhook_key_t *
itself(drhook_key_t *keyptr_self, 
       int tid, int opt, double *delta_time, 
       const double *walltime, const double *cputime) 
{
  drhook_key_t *keyptr = NULL;
  if (keyself) {
    keyptr = keyptr_self ? keyptr_self : keyself[tid-1];
    if (opt == 0) {
      if (opt_wallprof) keyptr->wall_in = walltime ? *walltime : WALLTIME();
      else              keyptr->cpu_in = cputime ? *cputime : CPUTIME();
      keyptr->calls++;
    }
    else if (opt == 1) {
      double delta = 0;
      if (opt_wallprof) {
	delta = walltime ? (*walltime - keyptr->wall_in) : (WALLTIME() - keyptr->wall_in);
	keyptr->delta_wall_all += delta;
      }
      else {
	delta = cputime ? (*cputime - keyptr->cpu_in) : (CPUTIME() - keyptr->cpu_in);
	keyptr->delta_cpu_all += delta;
      }
      if (delta_time) *delta_time = delta;
    }
  }
  return keyptr;
}

/*--- commie -routines : adds "," i.e. comma after each 3 digit, e.g.:
  1234567890 becomes more readable 1,234,567,890 */

static void 
lld_commie(long long int n, char sd[])
{
  const char comma = ',';
  char s[DRHOOK_STRBUF];
  char *p;
  int len, ncommas;
  sprintf(s,"%lld",n);
  len = strlen(s);
  ncommas = (len-1)/3;
  if (ncommas > 0) {
    char *pd = sd + len + ncommas;
    *pd-- = 0;
    p = s + len - 1;
    len = 0;
    while (p-s >= 0) {
      *pd-- = *p--;
      len++;
      if (p-s >= 0 && len%3 == 0) *pd-- = comma;
    }
  }
  else {
    strcpy(sd,s);
  }
}

static void 
dbl_commie(double n, char sd[])
{
  const char comma = ',';
  char s[DRHOOK_STRBUF];
  char *p;
  int len, ncommas;
  sprintf(s,"%.0f",n);
  len = strlen(s);
  ncommas = (len-1)/3;
  if (ncommas > 0) {
    char *pd = sd + len + ncommas;
    *pd-- = 0;
    p = s + len - 1;
    len = 0;
    while (p-s >= 0) {
      *pd-- = *p--;
      len++;
      if (p-s >= 0 && len%3 == 0) *pd-- = comma;
    }
  }
  else {
    strcpy(sd,s);
  }
}

/*--- callpath as a "pathname" ---*/

static void
unroll_callpath(FILE *fp, int len, 
		const equivalence_t *callpath, int callpath_len)
{
  if (fp && callpath && callpath_len > 0) {
    int j;
    for (j=0; j<callpath_len; callpath++, j++) {
      if (callpath && callpath->keyptr && callpath->keyptr->name) {
	const char *name = callpath->keyptr->name;
	int name_len = callpath->keyptr->name_len;
	len -= callpath_indent;
	if (len < 0) len = 0;
	fprintf(fp,"\n%*s%.*s",len," ",name_len,name);
      }
#ifdef DEBUG
      else {
	fprintf(fp,
		"\n????callpath=%p, callpath->keyptr=%p,  callpath->keyptr->name='%s'",
		callpath, callpath ? callpath->keyptr : 0,
		(callpath && callpath->keyptr && callpath->keyptr->name) ?
		callpath->keyptr->name : NIL);
      }
#endif
    }
  } /* if (fp) */
}


static equivalence_t *
get_callpath(int tid, int *callpath_len)
{
  int depth = 0;
  equivalence_t *callpath = NULL;
  if (tid >= 1 && tid <= numthreads) {
    const drhook_calltree_t *treeptr = thiscall[tid-1];
    while (treeptr && treeptr->active && depth < callpath_depth) {
      depth++;
      treeptr = treeptr->prev;
    }
    if (depth > 0) {
      int j = 0;
      callpath = malloc_drhook(sizeof(*callpath) * depth);
      treeptr = thiscall[tid-1];
      while (treeptr && treeptr->active && j < callpath_depth) {
	callpath[j].keyptr = treeptr->keyptr;
	j++;
	treeptr = treeptr->prev;
      }
    } /* if (depth > 0) */
  } /* if (tid >= 1 && tid <= numthreads) */
  if (callpath_len) *callpath_len = depth;
  return callpath;
}

/*--- profiler output ---*/

static int do_prof_off = 0;

static void
do_prof()
{
  /* to avoid recursive signals while atexit() (e.g. SIGXCPU) */
  if (signal_handler_ignore_atexit) return; 

  if (!do_prof_off && (opt_wallprof || opt_cpuprof)) {
    /* CPU, wall-clock and/or MFlop/s profiling */
    const int ftnunitno = 0;
    const int master = 1;
    const int print_option = 3;
    int initlev = 0;
    c_drhook_print_(&ftnunitno, &master, &print_option, &initlev);
  }

  if (!do_prof_off && opt_memprof) {
    /* Memory profiling */
    const int ftnunitno = 0;
    const int master = 1;
    const int print_option = 4;
    int initlev = 0;
    c_drhook_print_(&ftnunitno, &master, &print_option, &initlev);
  }

  if (!do_prof_off && timeline) {
    /* The last timeline-call */
    const int ftnunitno = opt_timeline_unitno;
    const int master = 1;
    const int print_option = -7;
    int initlev = 0;
    c_drhook_print_(&ftnunitno, &master, &print_option, &initlev);
  }
}

void c_drhook_prof_()
{
  if (ec_drhook) {
    do_prof();
    do_prof_off = 1;
  }
}

/*--- Check watch points ---*/

// Forward declarations of subroutines defined in dr_hook_prt.F90

void dr_hook_prt_logical_( const int* kunit, const void* ptr, const int* n );
void dr_hook_prt_char_( const int* kunit, const void* ptr, const int* n );
void dr_hook_prt_i4_( const int* kunit, const void* ptr, const int* n );
void dr_hook_prt_i8_( const int* kunit, const void* ptr, const int* n );
void dr_hook_prt_r4_( const int* kunit, const void* ptr, const int* n );
void dr_hook_prt_r8_( const int* kunit, const void* ptr, const int* n );

typedef enum { /* See dr_hook_watch_mod.F90 */
  KEYNONE =  0,
  KEYLOG  =  1,
  KEYCHAR =  2,
  KEY_I4  =  4,
  KEY_I8  =  8,
  KEY_R4  = 16,
  KEY_R8  = 32 
} PrintWatchKeys_t;

static void print_watch(int ftnunitno, int key, const void *ptr, int n)
{
  if (ptr && key > KEYNONE && n > 0) {
    int nmax = n;
    if (key == KEYLOG) {
      dr_hook_prt_logical_(&ftnunitno, ptr, &nmax);
    }
    else if (key == KEYCHAR) {
      dr_hook_prt_char_(&ftnunitno, ptr, &nmax);
    }
    else if (key == KEY_I4) {
      dr_hook_prt_i4_(&ftnunitno, ptr, &nmax);
    }
    else if (key == KEY_I8) {
      dr_hook_prt_i8_(&ftnunitno, ptr, &nmax);
    }
    else if (key == KEY_R4) {
      dr_hook_prt_r4_(&ftnunitno, ptr, &nmax);
    }
    else if (key == KEY_R8) {
      dr_hook_prt_r8_(&ftnunitno, ptr, &nmax);
    }
  }
}

static void 
check_watch(const char *label,
	    const char *name,
	    int name_len,
	    int allow_abort)
{
  if (watch) {
    int print_traceback = 1;
    drhook_watch_t *p = watch;
    coml_set_lockid_(&DRHOOK_lock);
    while (p) {
      if (p->active) {
	unsigned int crc32 = 0;
	int calc_crc = 0;
	const char *first_nbytes = p->ptr;
	int changed = memcmp(first_nbytes,p->ptr,p->watch_first_nbytes);
	if (!changed) {
	  /* The first nbytes were still the same; checking if crc has changed ... */
	  crc32_(p->ptr, &p->nbytes, &crc32);
	  changed = (crc32 != p->crc32);
	  calc_crc = 1;
	}
	if (changed) {
	  int tid = get_thread_id_();
	  char *pfx = PREFIX(tid);
	  if (!calc_crc) crc32_(p->ptr, &p->nbytes, &crc32);
	  fprintf(stderr,
		  "%s %s [%s@%s:%d] ***%s: Changed watch point '%s' at %p (%d bytes [#%d values])"
		  " -- %s %.*s : new crc32=%u\n",
		  pfx,TIMESTR(tid),FFL,
		  p->abort_if_changed ? "Error" : "Warning",
		  p->name, p->ptr, p->nbytes, p->nvals,
		  label, name_len, name, crc32);
	  print_watch(0, p->printkey, p->ptr, p->nvals);
	  if (print_traceback) {
	    LinuxTraceBack(pfx,TIMESTR(tid),NULL);
	    print_traceback = 0;
	  }
	  if (allow_abort && p->abort_if_changed) {
	    coml_unset_lockid_(&DRHOOK_lock); /* An important unlocking on Linux; otherwise hangs (until time-out) */
	    RAISE(SIGABRT);
	  }
#if 0
	  p->active = 0; /* No more these messages for this array */
	  watch_count--;
#else
	  p->crc32 = crc32;
#endif
	}
      }
      p = p->next;
    } /* while (p) */
    coml_unset_lockid_(&DRHOOK_lock);
  }
}

void
c_drhook_check_watch_(const char *where,
		      const int *allow_abort
		      /* Hidden length */
		      , int where_len)
{
  if (watch && watch_count > 0) check_watch("whilst at", where, where_len, *allow_abort);
}

/*** PUBLIC ***/

#define TIMERS \
  double walltime = opt_walltime ? WALLTIME() : 0;	\
  double cputime  = opt_cputime ? CPUTIME()  : 0;	\
  long long int hwm = opt_gethwm ? gethwm_() : 0;	\
  long long int stk = opt_getstk ? getstk_() : 0


/*=== c_drhook_set_lhook_ ===*/

void
c_drhook_set_lhook_(const int *lhook)
{
  if (lhook) drhook_lhook = *lhook;
}

/*=== c_drhook_getenv_ ===*/

void 
c_drhook_getenv_(const char *s, 
		 char *value,
		 /* Hidden arguments */
		 int slen,
		 const int valuelen) 
{
  char *env = NULL;
  char *p = malloc_drhook(slen+1);
  if (!p) {
    fprintf(stderr,"c_drhook_getenv_(): Unable to allocate %d bytes of memory\n", slen+1);
    RAISE(SIGABRT);
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
  free_drhook(p);
}


/*=== c_drhook_init_ ===*/

void 
c_drhook_init_(const char *progname,
	       const int *num_threads
	       /* Hidden length */
	       ,int progname_len)
{
  init_drhook(*num_threads);
  max_threads = MAX(1,*num_threads);
  if (a_out) free_drhook(a_out);
  progname = trim(progname, &progname_len);  
  if (progname_len > 0) {
    a_out = calloc_drhook(progname_len+1,sizeof(*progname));
    memcpy(a_out, progname, progname_len);
  }
  else {
    /* progname is a blank string;
       this is most likely due to a Fortran-call to getarg
       from program that has a C-main program, thus Fortran getarg
       may return a blank string */

    const char *arg0 = ec_GetArgs(0);
    if (arg0) {
      const char *pc = arg0;
      progname_len = strlen(pc);
      pc = trim(pc, &progname_len);
      a_out = strdup_drhook(pc);
    }
  }
  if (!a_out) {
    a_out = strdup_drhook("a.out"); /* Failed to obtain the name of the executing program */
  }
}


/*=== c_drhook_watch_ ===*/

void
c_drhook_watch_(const int *onoff,
		const char *array_name,
		const void *array_ptr,
		const int *nbytes,
		const int *abort_if_changed,
		const int *printkey,
		const int *nvals,
		const int *print_traceback_when_set
		/* Hidden length */
		,int array_name_len)
{
  int tid = get_thread_id_();
  drhook_watch_t *p = NULL;
  if (!drhook_lhook) return; 

  coml_set_lockid_(&DRHOOK_lock);

  /* check whether this array_ptr is already registered, but maybe inactive */
  p = watch;
  while (p) {
    if (p->ptr == array_ptr) {
      if (p->active) watch_count--;
      free_drhook(p->name);
      break;
    }
    p = p->next;
  }

  if (!p) {
    /* create new branch */
    p = calloc_drhook(1, sizeof(*p)); /* Implies p->next = NULL */
    if (!last_watch) {
      last_watch = watch = p;
    }
    else {
      last_watch->next = p;
      last_watch = p;
    }
  }

  p->name = strdup2_drhook(array_name,array_name_len);
  p->tid = tid;
  p->active = *onoff;
  if (p->active) watch_count++;
  p->abort_if_changed = *abort_if_changed;
  p->ptr = array_ptr;
  p->nbytes = *nbytes;
  p->watch_first_nbytes = MIN(p->nbytes, MAX_WATCH_FIRST_NBYTES);
  memcpy(p->first_nbytes,p->ptr,p->watch_first_nbytes);
  p->crc32 = 0;
  crc32_(p->ptr, &p->nbytes, &p->crc32);
  p->printkey = *printkey;
  p->nvals = *nvals;
  {
    char *pfx = PREFIX(p->tid);
    int ftnunitno = 0;
    int textlen = strlen(pfx) + strlen(p->name) + 256;
    char *text = malloc_drhook(textlen * sizeof(*text));
    snprintf(text,textlen,
	     "%s ***Warning: Set watch point '%s' at %p (%d bytes [%d values]) : crc32=%u",
	     pfx, p->name, p->ptr, p->nbytes, p->nvals, p->crc32);
    dr_hook_prt_(&ftnunitno, text, strlen(text));
    print_watch(ftnunitno, p->printkey, p->ptr, p->nvals);
    free_drhook(text);
    if (*print_traceback_when_set) LinuxTraceBack(pfx,TIMESTR(p->tid),NULL);
  }

  coml_unset_lockid_(&DRHOOK_lock);
}

/*=== c_drhook_start_ ===*/

void 
c_drhook_start_(const char *name, 
		const int *thread_id, 
		double *key,
		const char *filename,
		const int *sizeinfo
		/* Hidden length */
		,int name_len, int filename_len)
{
  TIMERS;
  equivalence_t u;
  ITSELF_0;
  if (!signals_set) signal_drhook_init(1);
  if (name_len > 0 && opt_funcenter == *thread_id) {
    fprintf(stdout,"<e> %d %d %.*s %lld %lld\n",myproc,*thread_id,name_len,name,hwm,stk);
    fflush(stdout);
  }
  if (watch && watch_count > 0) check_watch("when entering routine", name, name_len, 1);
  if (drhook_dump_hugepages) {
    int tid = *thread_id;
    char *pfx = PREFIX(tid);
    dump_hugepages(0,pfx,tid,0,-1);
  }
  if (!opt_callpath) {
    u.keyptr = getkey(*thread_id, name, name_len, 
		      filename, filename_len,
		      &walltime, &cputime,
		      NULL, 0, NULL);
  }
  else { /* (Much) more overhead */
    int free_callpath = 1;
    int callpath_len = 0;
    equivalence_t *callpath = get_callpath(*thread_id, &callpath_len);
    u.keyptr = getkey(*thread_id, name, name_len, 
		      filename, filename_len,
		      &walltime, &cputime,
		      callpath, callpath_len, &free_callpath);
    if (free_callpath) free_drhook(callpath);
  }
  if (cstklen == 0) {
    /* Double precision */
    *key = u.d;
  }
  else {
    /* Single precision : The variable "*key" is treated like max 4-byte entity -- "an index" */
    (void) callstack(*thread_id, key, u.keyptr);
  }
  ITSELF_1;
  if (opt_calltrace) {      
    coml_set_lockid_(&DRHOOK_lock);
    {
      const int ftnunitno = 0; /* stderr */
      const int print_option = 2; /* calling tree */
      int level = 0;
      c_drhook_print_(&ftnunitno, thread_id, &print_option, &level);
      /* fprintf(stderr,"%d#%d> %*.*s [%llu]\n",myproc,*thread_id,name_len,name_len,name,u.ull); */
    }
    coml_unset_lockid_(&DRHOOK_lock);
  }
  if (timeline) {
    int tid = *thread_id;
    if (opt_timeline_thread <= 0 || tid <= opt_timeline_thread) {
      drhook_timeline_t *tl = &timeline[tid-1];
      int bigjump = 1;
      unsigned long long int mod = (tl->calls[0]++)%opt_timeline_freq;
      double rss = (double)(getrss_()/1048576.0); /* in MBytes */
      double curheap = (opt_timeline_thread == 1 && tid == 1) ?
	(double)(getcurheap_()/1048576.0) : (double)(getcurheap_thread_(&tid)/1048576.0); /* in MBytes */
      double stack = (double)(getstk_()/1048576.0); /* in MBytes */
      double vmpeak = (double)(getvmpeak_()/1048576.0); /* in MBytes */
      if (mod != 0) {
	double inc_MB;
	inc_MB = tl->last_rss_MB - rss;
	if (ABS(inc_MB) < opt_timeline_MB) inc_MB = tl->last_curheap_MB - curheap;
	if (ABS(inc_MB) < opt_timeline_MB) inc_MB = tl->last_stack_MB - stack;
	if (ABS(inc_MB) < opt_timeline_MB) inc_MB = tl->last_vmpeak_MB - vmpeak;
	if (ABS(inc_MB) < opt_timeline_MB) bigjump = 0;
      }
      if (mod == 0 || bigjump) {
	coml_set_lockid_(&DRHOOK_lock);
	{
	  int ftnunitno = opt_timeline_unitno;
	  const int print_option = 5; /* calling "tree" with just the current entry */
	  int level = 0;
	  tl->last_rss_MB = rss;
	  tl->last_curheap_MB = curheap;
	  tl->last_stack_MB = stack;
	  tl->last_vmpeak_MB = vmpeak;
	  c_drhook_print_(&ftnunitno, &tid, &print_option, &level);
	}
	coml_unset_lockid_(&DRHOOK_lock);
      }
    } /* if (opt_timeline_thread <= 0 || tid <= opt_timeline_thread) */
  }
  if (opt_random_memstat > 0) random_memstat(*thread_id,0);
}

/*=== c_drhook_end_ ===*/

void 
c_drhook_end_(const char *name,
	      const int *thread_id,
	      const double *key,
	      const char *filename,
	      const int *sizeinfo
	      /* Hidden length */
	      ,int name_len, int filename_len)
{
  TIMERS;
  equivalence_t u;
  ITSELF_0;
  if (cstklen == 0) {
    /* Double precision */
    u.d = *key;
  }
  else {
    /* Single precision : The variable "*key" is treated like max 4-byte entity -- "an index" */
    u.keyptr = callstack(*thread_id, (void *)key, NULL);
  }
  /*
  if (opt_calltrace) {
    coml_set_lockid_(&DRHOOK_lock);
    fprintf(stderr,"%d#%d< %*.*s [%llu]\n",myproc,*thread_id,name_len,name_len,name,u.ull);
    coml_unset_lockid_(&DRHOOK_lock);
  }
  */
  if (name_len > 0 && opt_funcexit == *thread_id) {
    fprintf(stdout,"<x> %d %d %.*s %lld %lld\n",myproc,*thread_id,name_len,name,hwm,stk);
    fflush(stdout);
  }
  if (timeline) {
    int tid = *thread_id;
    if (opt_timeline_thread <= 0 || tid <= opt_timeline_thread) {
      drhook_timeline_t *tl = &timeline[tid-1];
      int bigjump = 1;
      unsigned long long int mod = (tl->calls[1]++)%opt_timeline_freq;
      double rss = (double)(getrss_()/1048576.0); /* in MBytes */
      double curheap = (opt_timeline_thread == 1 && tid == 1) ?
	(double)(getcurheap_()/1048576.0) : (double)(getcurheap_thread_(&tid)/1048576.0); /* in MBytes */
      double stack = (double)(getstk_()/1048576.0); /* in MBytes */
      double vmpeak = (double)(getvmpeak_()/1048576.0); /* in MBytes */
      if (mod != 0) {
	double inc_MB;
	inc_MB = tl->last_rss_MB - rss;
	if (ABS(inc_MB) < opt_timeline_MB) inc_MB = tl->last_curheap_MB - curheap;
	if (ABS(inc_MB) < opt_timeline_MB) inc_MB = tl->last_stack_MB - stack;
	if (ABS(inc_MB) < opt_timeline_MB) inc_MB = tl->last_vmpeak_MB - vmpeak;
	if (ABS(inc_MB) < opt_timeline_MB) bigjump = 0;
      }
      if (mod == 0 || bigjump) {
	coml_set_lockid_(&DRHOOK_lock);
	{
	  int ftnunitno = opt_timeline_unitno;
	  const int print_option = -5; /* calling "tree" with just the current entry */
	  int level = 0;
	  tl->last_rss_MB = rss;
	  tl->last_curheap_MB = curheap;
	  tl->last_stack_MB = stack;
	  tl->last_vmpeak_MB = vmpeak;
	  c_drhook_print_(&ftnunitno, &tid, &print_option, &level);
	}
	coml_unset_lockid_(&DRHOOK_lock);
      }
    } /* if (opt_timeline_thread <= 0 || tid <= opt_timeline_thread) */
  }
  if (watch && watch_count > 0) check_watch("when leaving routine", name, name_len, 1);
  putkey(*thread_id, u.keyptr, name, name_len, 
	 *sizeinfo,
	 &walltime, &cputime);
  ITSELF_1;
}

/*=== c_drhook_memcounter_ ===*/

void
c_drhook_memcounter_(const int *thread_id,
		     const long long int *size,
		     long long int *keyptr_addr)
{
  int tid = (thread_id && (*thread_id >= 1) && (*thread_id <= numthreads))
    ? *thread_id : get_thread_id_();
  int has_timeline = (timeline && size) ? opt_timeline : 0;
  if (has_timeline) {
    if (opt_timeline_thread <= 1 || tid <= opt_timeline_thread) {
      double size_MB = (double)((*size)/1048576.0); /* In MBytes */
      if (ABS(size_MB) < opt_timeline_MB) has_timeline = 0; /* Do not report */
    }
    else {
      has_timeline = 0; /* Do not report */
    }
  } /* if (has_timeline) */
  if (opt_memprof) {
    if (size) {
      union {
	long long int keyptr_addr;
	drhook_key_t *keyptr;
      } u;
      long long int alldelta;
      if (*size > 0) { /* Memory is being allocated */
	if (curkeyptr[tid-1]) {
	  drhook_key_t *keyptr = curkeyptr[tid-1];
	  keyptr->mem_curdelta += *size;
	  alldelta = keyptr->mem_curdelta + keyptr->mem_child;
	  if (alldelta > keyptr->maxmem_alldelta) keyptr->maxmem_alldelta = alldelta;
	  if (keyptr->mem_curdelta > keyptr->maxmem_selfdelta) 
	    keyptr->maxmem_selfdelta = keyptr->mem_curdelta;
	  if (keyptr_addr) {
	    u.keyptr = keyptr;
	    *keyptr_addr = u.keyptr_addr;
	  }
	  keyptr->alloc_count++;
	}
	else {
	  if (keyptr_addr) *keyptr_addr = 0;
	} /* if (curkeyptr[tid-1]) */
	/*
	fprintf(stderr,
		"memcounter: allocated %lld bytes ; *keyptr_addr = %lld\n",
		*size, *keyptr_addr);
	 */
      }
      else { /* Memory is being freed */
	drhook_key_t *keyptr;
	if (keyptr_addr && (*keyptr_addr)) {
	  u.keyptr_addr = *keyptr_addr;
	  keyptr = u.keyptr;
	}
	else 
	  keyptr = curkeyptr[tid-1];
	/*
	fprintf(stderr,
		"memcounter: DE-allocated %lld bytes ; *keyptr_addr = %lld\n",
		*size, *keyptr_addr);
	 */
	if (keyptr) {
	  long long int prev_curdelta = keyptr->mem_curdelta;
	  keyptr->mem_curdelta += *size;
	  alldelta = prev_curdelta + keyptr->mem_child;
	  if (alldelta > keyptr->maxmem_alldelta) keyptr->maxmem_alldelta = alldelta;
	  if (*size < 0) keyptr->free_count++;
	} /* if (keyptr) */
      } /* if (*size > 0) ... else */
    } /* if (size) */
  } /* if (opt_memprof) */
  if (has_timeline) {
    double curheap = (opt_timeline_thread == 1 && tid == 1) ?
      (double)(getcurheap_()/1048576.0) : (double)(getcurheap_thread_(&tid)/1048576.0); /* in MBytes */
    double rss = (double)(getrss_()/1048576.0); /* in MBytes */
    double stack = (double)(getstk_()/1048576.0); /* in MBytes */
    double vmpeak = (double)(getvmpeak_()/1048576.0); /* in MBytes */
    coml_set_lockid_(&DRHOOK_lock);
    {
      int ftnunitno = opt_timeline_unitno;
      double size_MB = (double)((*size)/1048576.0); /* In MBytes */
      int print_option = (size_MB > 0) ? 6 : -6; /* timeline upon c_drhook_memcounter_ & (big) ALLOCATE or DEALLOCATE */
      int level = 0;
      drhook_timeline_t *tl = &timeline[tid-1];
      tl->last_curheap_MB = curheap;
      tl->last_rss_MB = rss;
      tl->last_stack_MB = stack;
      tl->last_vmpeak_MB = vmpeak;
      c_drhook_print_(&ftnunitno, &tid, &print_option, &level);
    }
    coml_unset_lockid_(&DRHOOK_lock);
  } /* if (has_timeline) */
}

/*=== c_drhook_print_ ===*/

#define PRINT_HWM() \
if (opt_gethwm) { sprintf(s,",hwm=%lldK",keyptr->hwm/1024); s += strlen(s); }

#define PRINT_RSS() \
if (opt_getrss) { \
  sprintf(s,",rss/max=%lldK/%lldK",keyptr->rssnow/1024, keyptr->maxrss/1024); \
  s += strlen(s); \
}

#define PRINT_STK() \
if (opt_getstk) { \
  sprintf(s,",stack/max=%lldK/%lldK",keyptr->stack/1024, keyptr->maxstack/1024); \
  s += strlen(s); \
}

#define PRINT_PAG() \
if (opt_getpag) { \
  sprintf(s,",pag=%lld",keyptr->paging); \
  s += strlen(s); \
}

#define PRINT_WALL() \
if (opt_walltime) { \
  double self = keyptr->delta_wall_all-keyptr->delta_wall_child; \
  if (self < 0) self = 0; \
  sprintf(s,",wall=%.3fs/%.3fs", \
	  keyptr->delta_wall_all, self); \
  s += strlen(s); \
}

#define PRINT_CPU() \
if (opt_cputime) { \
  double self = keyptr->delta_cpu_all-keyptr->delta_cpu_child; \
  if (self < 0) self = 0; \
  sprintf(s,",cpu=%.3fs/%.3fs", \
	  keyptr->delta_cpu_all, self); \
  s += strlen(s); \
}

#define PRINT_CALLS() \
if (opt_calls) { \
  sprintf(s,",#%llu,st=%d",keyptr->calls,keyptr->status); \
  s += strlen(s); \
}

static int
prof_name_comp(const void *v1, const void *v2)
{
  const drhook_prof_t *p1 = v1;
  const drhook_prof_t *p2 = v2;
  return strcmp(p1->name,p2->name);
}

static int
memprof_name_comp(const void *v1, const void *v2)
{
  const drhook_memprof_t *p1 = v1;
  const drhook_memprof_t *p2 = v2;
  return strcmp(p1->name,p2->name);
}

static int
prof_pc_comp_desc(const void *v1, const void *v2)
{
  const drhook_prof_t *p1 = v1;
  const drhook_prof_t *p2 = v2;
  if (p1->pc < p2->pc) return 1;
  else if (p1->pc > p2->pc) return -1;
  else return 0;
}

static int
memprof_pc_comp_desc(const void *v1, const void *v2)
{
  const drhook_memprof_t *p1 = v1;
  const drhook_memprof_t *p2 = v2;
  if (p1->pc < p2->pc) return 1;
  else if (p1->pc > p2->pc) return -1;
  else return 0;
}

static const char *
trim_and_adjust_left(const char *p, int *name_len)
{
  int len = strlen(p);
  if (len > 0) {
    const char *back = &p[len-1];
    while (len > 0 && *back-- == ' ') len--;
    while (len > 0 && *p == ' ') { p++; len--; }
  }
  if (name_len) *name_len = len;
  return p;
}

static void print_routine_name0(FILE * fp, const char * p_name, int p_tid, const char * p_filename, int p_cluster, 
                                const equivalence_t * p_callpath, int p_callpath_len, int len, int cluster_size) 
{
  int name_len = 0; 
  const char *name = trim_and_adjust_left(p_name,&name_len); 

  if (callpath_packed) {

    if (p_callpath && p_callpath_len > 0) {
      const equivalence_t * callpath = &p_callpath[p_callpath_len-1];
      int j;
      for (j=0; j<p_callpath_len; callpath--, j++) 
        if (callpath && callpath->keyptr && callpath->keyptr->name) {
          const char *name = callpath->keyptr->name;
          int name_len = callpath->keyptr->name_len;
          fprintf(fp,"%.*s/",name_len,name);
        }
    } 
  } 

  fprintf(fp,"%.*s@%d%s%s", 
          name_len, name, 
          p_tid, 
          p_filename ? ":" : "", 
          p_filename ? p_filename : ""); 
  
  if (opt_clusterinfo) { 
    fprintf(fp," [%d,%d]", 
            p_cluster, ABS(cluster_size)); 
  } 
    
  if (!callpath_packed) 
    unroll_callpath(fp, len, p_callpath, p_callpath_len); 
  

}

#define print_routine_name(fp, p, len, cluster_size) \
  if (fp && p) { \
    print_routine_name0(fp, p->name, p->tid, p->filename, p->cluster, \
                        p->callpath, p->callpath_len, len, cluster_size);\
  } /* if (fp && p) */


static void
DrHookPrint(int ftnunitno, const char *line)
{
  if (line) {
    FILE *fp = NULL;
    if (ftnunitno <= 0) 
      fp = stderr;
    else if (ftnunitno == 6) 
      fp = stdout;
    else
      dr_hook_prt_(&ftnunitno, line, strlen(line));
    OPTPRINT(fp,"%s\n",line);
  }
}

void 
c_drhook_print_(const int *ftnunitno,
		const int *thread_id,
		const int *print_option, /* 
					    1=raw call counts 
					    2=calling tree
					    3=profiling info
					    4=memory profiling
					    5=timeline upon entering the routine
					   -5=timeline upon leaving the routine
					    6=timeline upon c_drhook_memcounter_ & (big) ALLOCATE
					   -6=timeline upon c_drhook_memcounter_ & (big) DEALLOCATE
					    7=timeline : the very first call (upon setup or dr.hook)
					   -7=timeline : the very last call (in atexit())
					 */
		int *level
		)
{
  static int first_time = 0;
  int tid = (thread_id && (*thread_id >= 1) && (*thread_id <= numthreads))
    ? *thread_id : get_thread_id_();
  int mytid = get_thread_id_();
  char *pfx = PREFIX(tid);
  if (ftnunitno && keydata && calltree) {
    char line[4096];
    int abs_print_option = ABS(*print_option);
    int j;

  /* Mod to call traceback and continue if called with level=99 */
      if(*level == 99) {
        *level=0;
      }
      else {
        if(*print_option == 2) {
          if(first_time == 1) return;
          first_time = 1;
        }
      }
  /* end of Mod  */

    if (*print_option == 1) { /* raw call counts */
      for (j=0; j<hashsize; j++) {
	int nestlevel = 0;
	drhook_key_t *keyptr = &keydata[tid-1][j];
	while (keyptr) {
	  if (keyptr->name) {
	    char *s = line;
	    sprintf(s,
		    "%s %s [%s@%s:%d] [hash#%d,nest=%d] '%s'",
		    pfx,TIMESTR(tid),FFL,
		    j,nestlevel,keyptr->name);
	    s += strlen(s);
	    PRINT_CALLS();
	    PRINT_HWM();
	    PRINT_RSS();
	    PRINT_STK();
	    PRINT_PAG();
	    PRINT_WALL();
	    PRINT_CPU();
	    *s = 0;
	    DrHookPrint(*ftnunitno, line);
	  }
	  keyptr = keyptr->next;
	  nestlevel++;
	} /* while (keyptr) */
      } /* for (j=0; j<hashsize; j++) */
    }

    else if (*print_option == 2 || 
	     abs_print_option == 5 || 
	     abs_print_option == 6 || 
	     abs_print_option == 7
	     ) { /* the current calling tree */
      drhook_calltree_t *treeptr = calltree[tid-1];

      if (*print_option == 2) { 
	long long int hwm = getmaxhwm_()/1048576;
	long long int rss = getmaxrss_()/1048576;
	long long int maxstack = getmaxstk_()/1048576;
	long long int vmpeak = getvmpeak_()/1048576;
	snprintf(line,sizeof(line),
		 "%s %s [%s@%s:%d] %lld MB (maxheap), %lld MB (maxrss), %lld MB (maxstack), %lld MB (vmpeak)",
		 pfx,TIMESTR(tid),FFL,
		 hwm,rss,maxstack,vmpeak);
	DrHookPrint(*ftnunitno, line);
      }

      if (tid > 1) {
	if (*print_option == 2) {
	  /* I'm not a master thread, but my master has the beginning of the calltree */
	  int initlev = 0;
	  const int master = 1;
          first_time = 0;
	  c_drhook_print_(ftnunitno, &master, print_option, &initlev);
	  *level += initlev;
	}
	else if (tid > opt_timeline_thread) {
	  return;
	}
      }

      if (abs_print_option == 7) {
	treeptr = NULL;
      }
      else if (abs_print_option == 5 || abs_print_option == 6) {
	treeptr = thiscall[tid-1];
      }
      else {
	treeptr = calltree[tid-1];
      }

      while (abs_print_option == 7 || (treeptr && treeptr->active)) {
	int do_print = (*print_option == 2 || 
			abs_print_option == 7 ||
			abs_print_option == 5 || abs_print_option == 6);
	if (do_print) {
	  drhook_key_t *keyptr = (abs_print_option == 7) ? NULL : treeptr->keyptr;
	  char *s = line;
	  char is_timeline = 1, kind;
	  switch (*print_option) {
	  case -5: kind = '<'; break;
	  case -6: kind = '-'; break;
	  case -7: kind = 'E'; break;
	  case  5: kind = '>'; break;
	  case  6: kind = '+'; break;
	  case  7: kind = 'B'; break;
	  default:
	  case 2: kind = ':'; is_timeline = 0; break;
	  }
	  if (*print_option == 2 || 
	      (is_timeline && tid > 1 && tid <= opt_timeline_thread))  {
	    sprintf(s,"%s %s [%s@%s:%d] %s%c ",
		    pfx,TIMESTR(tid),FFL,
		    is_timeline ? "tl:" : "",
		    kind);
	  }
	  else if (is_timeline && opt_timeline_thread == 1 && tid == 1) {
	    sprintf(s,"%s %s [%s@%s:%d] %s%c ",
		    pfx,TIMESTR(tid),FFL,
		    is_timeline ? "tl:" : "",
		    kind);
	  }
	  s += strlen(s);
	  (*level)++;
	  for (j=0; j<(*level); j++) *s++ = ' ';
	  if (*print_option == 2) {
            if(mytid != tid) { /* We are printing the master call tree as far as >OMP*/
              if(strncmp(">OMP",keyptr->name,4) == 0) {
                (*level)--;
                return;
              }
            }
	    sprintf(s,"%s ",keyptr->name);
	    s += strlen(s);
	  }
	  if (is_timeline) {
	    double wall = WALLTIME();
	    double rss, curheap, stack, vmpeak;
	    drhook_timeline_t *tl = &timeline[tid-1];
	    if (abs_print_option == 5 || abs_print_option == 6) { /* when called via drhook_begin/_end or memcounter */
	      curheap = tl->last_curheap_MB;
	      rss = tl->last_rss_MB;
	      stack = tl->last_stack_MB;
	      vmpeak = tl->last_vmpeak_MB;
	    }
	    else {
	      rss = (double)(getrss_()/1048576.0); /* in MBytes */
	      curheap = (opt_timeline_thread == 1 && tid == 1) ?
		(double)(getcurheap_()/1048576.0) : (double)(getcurheap_thread_(&tid)/1048576.0); /* in MBytes */
	      stack = (double)(getstk_()/1048576.0); /* in MBytes */
	      vmpeak = (double)(getvmpeak_()/1048576.0); /* in MBytes */
	      tl->last_curheap_MB = curheap;
	      tl->last_rss_MB = rss;
	      tl->last_stack_MB = stack;
	      tl->last_vmpeak_MB = vmpeak;
	    }
	    if (opt_timeline_format == 1) {
	      sprintf(s, "%.6f %.4g %.4g %.4g %.4g", wall, rss, curheap, stack, vmpeak);
	    }
	    else {
	      sprintf(s,
		      "wall=%.6f cpu=%.4g hwm=%.4g rss=%.4g curheap=%.4g stack=%.4g vmpeak=%.4g pag=%lld",
		      wall, CPUTIME(),
		      (double)(gethwm_()/1048576.0), rss,
		      curheap,
		      (double)(getstk_()/1048576.0),
		      (double)(getvmpeak_()/1048576.0),
		      getpag_());
	    }
	    s += strlen(s);
	    *s++ = ' ';
	    if (keyptr) {
	      sprintf(s,"'%s'",keyptr->name);
	    }
	    else {
	      sprintf(s,"'#PROGRAM %s'",(*print_option == 7) ? "BEGIN" : "END");
	    }
	    s += strlen(s);
	    {
	      int current_numth = 0;
	      coml_get_num_threads_(&current_numth);
	      sprintf(s,"[#%d]",current_numth);
	      s += strlen(s);
	    }
	  }
	  else {
	    PRINT_CALLS();
	    PRINT_HWM();
	    PRINT_RSS();
	    PRINT_STK();
	    PRINT_PAG();
	    PRINT_WALL();
	    PRINT_CPU();
	  }
	  *s = 0;
	  DrHookPrint(*ftnunitno, line);
	}
	if (abs_print_option == 7 || abs_print_option == 5 || abs_print_option == 6) break;
	if (treeptr) treeptr = treeptr->next;
      } /* while (abs_print_option == 7 || (treeptr && treeptr->active)) */
    }

    else if (*print_option == 3) { /* profiling (CPU, wall-clock and/or MFlop/s) */
      int len;
      int t;
      double cumul;
      double tottime = 0, max_overhead_pc = 0;
      double *tot = NULL;
      int nprof = 0;
      drhook_prof_t *prof = NULL;
      drhook_prof_t *p;
      double flop_tot = 0, instr_tot = 0;
      double *flop = NULL, *instr = NULL;

      if (!opt_wallprof && !opt_cpuprof) return; /* no profiling info available */
      if (tid > 1) return; /* just master thread allowed ; takes care of siblings, too */
      if (numthreads<=0) return;
      if (do_prof_off) return;
      do_prof_off = 1;

      /* Insert "$drhook" */
      if (keyself && opt_self > 1) {
	for (t=0; t<numthreads; t++) (void) insertkey(t+1,keyself[t]);
      }

      flop = calloc_drhook(numthreads, sizeof(*flop));
      instr = calloc_drhook(numthreads, sizeof(*instr));
      tot = calloc_drhook(numthreads, sizeof(*tot));

      for (t=0; t<numthreads; t++) {
	for (j=0; j<hashsize; j++) {
	  drhook_key_t *keyptr = &keydata[t][j];
	  while (keyptr) {
	    if (keyptr->name && (keyptr->status == 0 || signal_handler_called)) {
	      double self;
	      if (opt_wallprof) {
		self = keyptr->delta_wall_all - keyptr->delta_wall_child;
	      }
	      else {
		self = keyptr->delta_cpu_all - keyptr->delta_cpu_child;
	      }
	      /* if (self < 0) self = 0; */
	      tot[t] += self;
#ifdef HPM
	      flop[t] += keyptr->avg_mflops * self; /* mflop_count(keyptr); */
	      instr[t] += keyptr->avg_mipsrate * self; /* mip_count(keyptr); */
#endif
	      nprof++;
	    }
	    keyptr = keyptr->next;
	  } /* while (keyptr && keyptr->status == 0) */
	} /* for (t=0; t<numthreads; t++) */
      } /* for (j=0; j<hashsize; j++) */

      if (opt_wallprof) { /* a bit unreliable; had not taken max. value of threads wall yet; will be recalculated */
	tottime = tot[0] + ((keyself && opt_self > 1) ? keyself[0]->delta_wall_all : 0);
	for (t=1; t<numthreads; t++) {
	  double tmp = tot[t] + ((keyself && opt_self > 1) ? keyself[t]->delta_wall_all : 0);
	  tottime = MAX(tottime,tmp);
	}
      }
      else { /* ok & reliable (for cpuprof) */
	tottime = 0;
	for (t=0; t<numthreads; t++) tottime += (tot[t] + ((keyself && opt_self > 1) ? keyself[t]->delta_cpu_all : 0));
      }

      if (tottime <= 0) tottime = 1e-10;

      p = prof = calloc_drhook(nprof + 1, sizeof(*prof)); /* Make sure there is at least one entry */

      for (t=0; t<numthreads; t++) {
	for (j=0; j<hashsize; j++) {
	  drhook_key_t *keyptr = &keydata[t][j];
	  while (keyptr) {
	    if (keyptr->name && (keyptr->status == 0 || signal_handler_called)) {
	      p->self = opt_wallprof ?
		keyptr->delta_wall_all - keyptr->delta_wall_child :
		keyptr->delta_cpu_all - keyptr->delta_cpu_child;
	      p->total = opt_wallprof ?
		keyptr->delta_wall_all :
		keyptr->delta_cpu_all;
	      p->calls = keyptr->calls;
	      p->name = keyptr->name;
	      p->pc = (p->self/tottime) * 100.0;
	      if (p->calls > 0) {
		p->percall_ms_self = (p->self/p->calls) * 1000.0;
		p->percall_ms_total = (p->total/p->calls) * 1000.0;
	      }
	      p->tid = t+1;
	      p->index = p - prof;
#ifdef HPM
	      if (opt_hpmprof) {
		p->mflops = keyptr->avg_mflops; /* mflops_hpm(keyptr); */
		p->mipsrate = keyptr->avg_mipsrate; /* mips_hpm(keyptr); */
		p->divpc = divpc_hpm(keyptr);
	      }
#endif
	      p->filename = keyptr->filename;
	      p->sizeinfo = keyptr->sizeinfo;
	      p->min_sizeinfo = keyptr->min_sizeinfo;
	      p->max_sizeinfo = keyptr->max_sizeinfo;
	      p->sizespeed = (p->self > 0 && p->sizeinfo > 0) ? p->sizeinfo/p->self : 0;
	      p->sizeavg = (p->calls > 0 && p->sizeinfo > 0) ? p->sizeinfo/p->calls : 0;
	      p->callpath = keyptr->callpath;
	      p->callpath_len = keyptr->callpath_len;
	      p++;
	    }
	    keyptr = keyptr->next;
	  } /* while (keyptr && keyptr->status == 0) */
	} /* for (j=0; j<hashsize; j++) */
      } /* for (t=0; t<numthreads; t++) */

      do {
	double mflop_rate = 0;
	double mip_rate = 0;
	int numroutines = 0;
	int cluster;
	double *maxval = calloc_drhook(nprof+1, sizeof(*maxval)); /* make sure at least 1 element */
	int *clusize = calloc_drhook(nprof+1, sizeof(*clusize)); /* make sure at least 1 element */
	char *prevname = NULL;
	const char *fmt1 = "%5d %8.2f %12.3f %12.3f %12.3f %14llu %11.2f %11.2f   %s";
	const char *fmt2 = "%5d %8.2f %12.3f %12.3f %12.3f %14llu %7.0f %7.0f %7.1f   %s";
	const char *fmt = opt_hpmprof ? fmt2 : fmt1;
	char *filename = get_mon_out(myproc);
	FILE *fp = NULL;

	if (!filename) break;

	if ((myproc == 1 && mon_out_procs == -1) || mon_out_procs == myproc) {
	  fprintf(stderr,
		  "%s %s [%s@%s:%d] Writing profiling information of proc#%d into file '%s'\n",
		  pfx,TIMESTR(tid),FFL,
		  myproc,filename);
	}

	fp = fopen(filename,"w");
	if (!fp) goto finish_3;
	
	/* alphanumerical sorting to find out clusters of the same routine but on different threads */
	/* also find out total wall clock time */
	/* calculate percentage values */
	
	p = prof;
	qsort(p, nprof, sizeof(*p), prof_name_comp);

	cluster = 0;
	maxval[cluster] = p->self;
	p->maxval = &maxval[cluster];
	clusize[cluster] = 1;
	prevname = p->name;
	p++;
	for (j=1; j<nprof; j++) {
	  if (!strequ(prevname,p->name)) {
	    (p-1)->cluster = cluster;
	    (p-1)->maxval = &maxval[cluster];
	    prevname = p->name;
	    cluster++;
	  }
	  if (p->self > maxval[cluster]) maxval[cluster] = p->self;
	  p->cluster = cluster;
	  p->maxval = &maxval[cluster];
	  clusize[cluster]++;
	  p++;
	} /* for (j=1; j<nprof; j++) */

	numroutines = (nprof > 0) ? (cluster + 1) : 0; /* Active no. of routines */

	if (opt_wallprof) tottime = 0;
	p = prof;
	for (j=0; j<nprof; j++) {
	  int use_this = 0;
	  cluster = p->cluster;
	  if (clusize[cluster] > 1) { /* multiple threads <= numthreads indeed called this routine */
	    p->is_max = (p->self == *p->maxval);
	    if (p->is_max) { /* first max found will be used for total time */
	      clusize[cluster] = -clusize[cluster]; /* ensures that max has been found for this cluster */
	      use_this = opt_wallprof;
	    }
	  }
	  else if (clusize[cluster] == 1) {
	    use_this = opt_wallprof;
	  }
	  if (use_this && opt_wallprof) tottime += p->self;
	  p++;
	}

        if (tottime <= 0) tottime = 1e-10;

	if (opt_wallprof) { /* use re-calculated tottime to define percentages */
	  p = prof;
	  for (j=0; j<nprof; j++) {
	    p->pc = (p->self/tottime) * 100.0;
	    p++;
	  }
	}

	/* sorting with respect to percentage value */

	p = prof;
	qsort(p, nprof, sizeof(*p), prof_pc_comp_desc);

	flop_tot = 0;
	instr_tot = 0;
	max_overhead_pc = 0;
	for (t=0; t<numthreads; t++) {
	  flop_tot += flop[t];
	  instr_tot += instr[t];
	  if (overhead) {
	    max_overhead_pc = MAX(max_overhead_pc,overhead[t]);
#ifdef DEBUG
	    fprintf(fp,"tid#%d: overhead = %.15g s\n",t+1,overhead[t]);
#endif
	  }
	}
#ifdef DEBUG
	fprintf(fp,"max overhead = %.15g s, tottime = %.15g s\n",
		max_overhead_pc, tottime);
#endif
	if (tottime - max_overhead_pc > 0) {
	  max_overhead_pc = 100.0*(max_overhead_pc/(tottime - max_overhead_pc));
	}
	else {
	  max_overhead_pc = 100;
	}

	fprintf(fp,
		"Profiling information for program='%s', proc#%d:\n",a_out, myproc);
	fprintf(fp,"\tNo. of instrumented routines called : %d\n", numroutines);
	fprintf(fp,"\tInstrumentation started : %s\n",start_stamp ? start_stamp : "N/A");
	end_stamp = timestamp();
	fprintf(fp,"\tInstrumentation   ended : %s\n",end_stamp ? end_stamp : "N/A");
	fprintf(fp,"\tInstrumentation overhead: %.2f%%\n",max_overhead_pc);
	{
	  long long int hwm = getmaxhwm_()/1048576;
	  long long int rss = getmaxrss_()/1048576;
	  long long int maxstack = getmaxstk_()/1048576;
	  long long int vmpeak = getvmpeak_()/1048576;
	  long long int pag = getpag_();
	  fprintf(fp,
	  "\tMemory usage : %lld MB (heap), %lld MB (rss), %lld MB (stack), %lld MB (vmpeak), %lld (paging)\n",
		  hwm,rss,maxstack,vmpeak,pag);
	}
	if (opt_hpmprof) {
	  mflop_rate = flop_tot / tottime;
	  mip_rate = instr_tot / tottime;
	  fprintf(fp,
		  "\t%s-time is %.2f sec on proc#%d, %.0f MFlops (ops#%.0f*10^6), %.0f MIPS (ops#%.0f*10^6) (%d procs, %d threads)\n",
		  opt_wallprof ? "Wall" : "Total CPU", tottime, myproc,
		  mflop_rate, flop_tot, mip_rate, instr_tot,
		  nproc, numthreads);
	}
	else {
	  fprintf(fp,
		  "\t%s-time is %.2f sec on proc#%d (%d procs, %d threads)\n",
		  opt_wallprof ? "Wall" : "Total CPU", tottime, myproc,
		  nproc, numthreads);
	}

	if (myproc == 1) {
	  fprintf(stderr,
		  "Profiling information for program='%s', proc#%d:\n",a_out, myproc);
	  fprintf(stderr,"\tNo. of instrumented routines called : %d\n", numroutines);
	  fprintf(stderr,"\tInstrumentation started : %s\n",start_stamp ? start_stamp : "N/A");
	  fprintf(stderr,"\tInstrumentation   ended : %s\n",end_stamp ? end_stamp : "N/A");
	  fprintf(stderr,"\tInstrumentation overhead: %.2f%%\n",max_overhead_pc);
	  if (opt_hpmprof) {
	    fprintf(stderr,
		  "\t%s-time is %.2f sec on proc#%d, %.0f MFlops (ops#%.0f*10^6), %.0f MIPS (ops#%.0f*10^6) (%d procs, %d threads)\n",
		  opt_wallprof ? "Wall" : "Total CPU", tottime, myproc,
		  mflop_rate, flop_tot, mip_rate, instr_tot,
		  nproc, numthreads);
	  }
	  else {
	    fprintf(stderr,
		    "\t%s-time is %.2f sec on proc#%d (%d procs, %d threads)\n",
		    opt_wallprof ? "Wall" : "Total CPU", tottime, myproc,
		    nproc, numthreads);
	  }
	} /* if (myproc == 1) */

	free_drhook(end_stamp);

	for (t=0; t<numthreads; t++) {
	  double tmp = 100.0*(tot[t]/tottime);
	  if (opt_hpmprof && tot[t] > 0) {
	    mflop_rate = flop[t]/tot[t];
	    mip_rate = instr[t]/tot[t];
	  }
	  else {
	    mflop_rate = 0;
	    mip_rate = 0;
	  }
	  fprintf(    fp,"\tThread#%d: %11.2f sec (%.2f%%)",t+1,tot[t],tmp);
	  if (opt_hpmprof) fprintf(    fp,", %.0f MFlops (ops#%.0f*10^6), %.0f MIPS (ops#%.0f*10^6)", mflop_rate, flop[t], mip_rate, instr[t]);
	  fprintf(    fp,"\n");
	  if (myproc == 1) {
	    fprintf(stderr,"\tThread#%d: %11.2f sec (%.2f%%)",t+1,tot[t],tmp);
	    if (opt_hpmprof) fprintf(stderr,", %.0f MFlops (ops#%.0f*10^6), %.0f MIPS (ops#%.0f*10^6)", mflop_rate, flop[t], mip_rate, instr[t]);
	    fprintf(stderr,"\n");
	  }
	}

	fprintf(fp,"\n");
	if (opt_hpmprof) {
	  len = 
	    fprintf(fp,"    #  %% Time         Cumul         Self        Total     # of calls    MIPS  MFlops   Div-%%    ");
	}
	else {
	  len = 
	    fprintf(fp,"    #  %% Time         Cumul         Self        Total     # of calls        Self       Total    ");
	}
	fprintf(fp,"Routine@<thread-id>");
	if (opt_clusterinfo) fprintf(fp," [Cluster:(id,size)]");
	fprintf(fp,"\n");
	if (opt_sizeinfo) fprintf(fp,"%*s %s\n",len-20," ","(Size; Size/sec; Size/call; MinSize; MaxSize)");
	if (opt_hpmprof) {
	  fprintf(fp,  "        (self)        (sec)        (sec)        (sec)                                       \n");
	}
	else {
	  fprintf(fp,  "        (self)        (sec)        (sec)        (sec)                    ms/call     ms/call\n");
	}
	fprintf(fp,"\n");

	cumul = 0;
	for (j=0; j<nprof; ) {
	  int cluster_size = clusize[p->cluster];
	  if (p->pc < percent_limit) break;
	  if (opt_cputime) {
	    cumul += p->self;
	  }
	  else {
	    if (p->is_max || cluster_size == 1) cumul += p->self;
	  }
	  if (opt_hpmprof) {
	    fprintf(fp, fmt,
		    ++j, p->pc, cumul, p->self, p->total, p->calls,
		    p->mipsrate, p->mflops, p->divpc,
		    p->is_max ? "*" : " ");
	  }
	  else {
	    fprintf(fp, fmt,
		    ++j, p->pc, cumul, p->self, p->total, p->calls,
		    p->percall_ms_self, p->percall_ms_total, 
		    p->is_max ? "*" : " ");
	  }

	  print_routine_name(fp, p, len, cluster_size);
	    
	  if (opt_sizeinfo && p->sizeinfo > 0) {
	    char s1[DRHOOK_STRBUF], s2[DRHOOK_STRBUF], s3[DRHOOK_STRBUF];
	    char s4[DRHOOK_STRBUF], s5[DRHOOK_STRBUF];
	    lld_commie(p->sizeinfo,s1);
	    dbl_commie(p->sizespeed,s2);
	    dbl_commie(p->sizeavg,s3);
	    lld_commie(p->min_sizeinfo,s4);
	    lld_commie(p->max_sizeinfo,s5);
	    fprintf(fp,"\n%*s (%s; %s; %s; %s; %s)",len-20," ",s1,s2,s3,s4,s5);
	  }
	  fprintf(fp,"\n");
	  p++;
	} /* for (j=0; j<nprof; ) */
	
	fclose(fp);
      finish_3:
	free_drhook(filename);
	free_drhook(maxval);
	free_drhook(clusize);
      } while (0);

      free_drhook(instr);
      free_drhook(flop);
      free_drhook(tot);
      free_drhook(prof);
      do_prof_off = 0;
    }

    else if (*print_option == 4) { /* Memory profiling */
      int t, len;
      int nprof = 0;
      drhook_memprof_t *prof = NULL;
      drhook_memprof_t *p;
      long long int *tot;
      long long int *maxseen_tot;
      double totmaxmem_delta;

      if (!opt_memprof) return; /* no profiling info available */
      if (tid > 1) return; /* just master thread allowed ; takes care of siblings, too */
      if (numthreads<=0) return;
      if (do_prof_off) return;
      do_prof_off = 1;

      tot = calloc_drhook(numthreads, sizeof(*tot));
      maxseen_tot = calloc_drhook(numthreads, sizeof(*maxseen_tot));

      for (t=0; t<numthreads; t++) {
	for (j=0; j<hashsize; j++) {
	  drhook_key_t *keyptr = &keydata[t][j];
	  while (keyptr) {
	    if (keyptr->name && (keyptr->status == 0 || signal_handler_called)) {

	      long long int self;
	      self = keyptr->maxmem_selfdelta;
	      if (self < 0) self = 0;
	      tot[t] += self;
	      maxseen_tot[t] = MAX(maxseen_tot[t], keyptr->mem_seenmax);
	      nprof++;
	    }
	    keyptr = keyptr->next;
	  } /* while (keyptr && keyptr->status == 0) */
	} /* for (t=0; t<numthreads; t++) */
      } /* for (j=0; j<hashsize; j++) */

      totmaxmem_delta = tot[0];
      for (t=1; t<numthreads; t++) {
	long long int tmp = tot[t];
	totmaxmem_delta = MAX(totmaxmem_delta,tmp);
      }

      if (totmaxmem_delta <= 0) totmaxmem_delta = 1e-10; /* To avoid divide-by-zero */

      p = prof = calloc_drhook(nprof + 1, sizeof(*prof)); /* Make sure there is at least one entry */

      for (t=0; t<numthreads; t++) {
	for (j=0; j<hashsize; j++) {
	  drhook_key_t *keyptr = &keydata[t][j];
	  while (keyptr) {
	    if (keyptr->name && (keyptr->status == 0 || signal_handler_called)) {
	      p->self = keyptr->maxmem_selfdelta;
	      p->children = keyptr->mem_child;
	      p->hwm = keyptr->mem_maxhwm;
	      p->rss = keyptr->mem_maxrss;
	      p->stk = keyptr->mem_maxstk;
	      p->pag = keyptr->mem_maxpagdelta;
	      p->leaked = keyptr->mem_curdelta;
	      p->calls = keyptr->calls;
	      p->alloc_count += keyptr->alloc_count;
	      p->free_count += keyptr->free_count;
	      p->name = keyptr->name;
	      p->pc = (p->self/totmaxmem_delta) * 100.0;
	      p->tid = t+1;
	      p->index = p - prof;
	      p->filename = keyptr->filename;
	      p->callpath = keyptr->callpath;
	      p->callpath_len = keyptr->callpath_len;
	      p++;
	    }
	    keyptr = keyptr->next;
	  } /* while (keyptr && keyptr->status == 0) */
	} /* for (t=0; t<numthreads; t++) */
      } /* for (j=0; j<hashsize; j++) */

      do {
	int numroutines = 0;
	int cluster;
	long long int *maxval = calloc_drhook(nprof+1, sizeof(*maxval)); /* make sure at least 1 element */
	int *clusize = calloc_drhook(nprof+1, sizeof(*clusize)); /* make sure at least 1 element */
	char *prevname = NULL;
	const char *fmt1 = "%5d %9.2f  %14lld %14lld %14lld %14lld %14lld %10lld %10llu %10llu%s%10llu   %s";
	const char *fmt = fmt1;
	char *filename = get_memmon_out(myproc);
	FILE *fp = NULL;

	if (!filename) break;

	if ((myproc == 1 && mon_out_procs == -1) || mon_out_procs == myproc) {
	  fprintf(stderr,"Writing memory-profiling information of proc#%d into file '%s'\n",myproc,filename);
	}

	fp = fopen(filename,"w");
	if (!fp) goto finish_4;
	
	/* alphanumerical sorting to find out clusters of the same routine but on different threads */
	
	p = prof;
	qsort(p, nprof, sizeof(*p), memprof_name_comp);

	cluster = 0;
	maxval[cluster] = p->self;
	p->maxval = &maxval[cluster];
	clusize[cluster] = 1;
	prevname = p->name;
	p++;
	for (j=1; j<nprof; j++) {
	  if (!strequ(prevname,p->name)) {
	    (p-1)->cluster = cluster;
	    (p-1)->maxval = &maxval[cluster];
	    prevname = p->name;
	    cluster++;
	  }
	  if (p->self > maxval[cluster]) maxval[cluster] = p->self;
	  p->cluster = cluster;
	  p->maxval = &maxval[cluster];
	  clusize[cluster]++;
	  p++;
	} /* for (j=1; j<nprof; j++) */

	numroutines = (nprof > 0) ? (cluster + 1) : 0; /* Active no. of routines */

	totmaxmem_delta = 0;
	p = prof;
	for (j=0; j<nprof; j++) {
	  int use_this = 0;
	  cluster = p->cluster;
	  if (clusize[cluster] > 1) { /* multiple threads <= numthreads indeed called this routine */
	    p->is_max = (p->self == *p->maxval);
	    if (p->is_max) { /* first max found will be used for total time */
	      clusize[cluster] = -clusize[cluster]; /* ensures that max has been found for this cluster */
	      use_this = 1;
	    }
	  }
	  else if (clusize[cluster] == 1) {
	    use_this = 1;
	  }
	  if (use_this) totmaxmem_delta += p->self;
	  p++;
	}

        if (totmaxmem_delta <= 0) totmaxmem_delta = 1e-10; /* To avoid divide-by-zero */

	/* use re-calculated totmaxmem_delta to define percentages */
	p = prof;
	for (j=0; j<nprof; j++) {
	  p->pc = (p->self/totmaxmem_delta) * 100.0;
	  p++;
	}

	/* sorting with respect to percentage value */

	p = prof;
	qsort(p, nprof, sizeof(*p), memprof_pc_comp_desc);

	fprintf(fp,
		"Memory-profiling information for program='%s', proc#%d:\n",a_out, myproc);
	fprintf(fp,"\tNo. of instrumented routines called : %d\n", numroutines);
	fprintf(fp,"\tInstrumentation started : %s\n",start_stamp ? start_stamp : "N/A");
	end_stamp = timestamp();
	fprintf(fp,"\tInstrumentation   ended : %s\n",end_stamp ? end_stamp : "N/A");
	{
	  long long int hwm = gethwm_()/1048576;
	  long long int rss = getrss_()/1048576;
	  long long int maxstack = getmaxstk_()/1048576;
	  long long int vmpeak = getvmpeak_()/1048576;
	  long long int pag = getpag_();
	  long long int maxseen = 0;
	  long long int leaked = 0;
	  p = prof;
	  for (j=0; j<nprof; j++) {
	    if (p->leaked > 0) leaked += p->leaked;
	    p++;
	  }
	  for (t=0; t<numthreads; t++) {
	    maxseen += maxseen_tot[t];
	  }
	  maxseen /= 1048576;
	  leaked /= 1048576;
	  fprintf(fp,
	  "\tMemory usage : %lld MB (max.seen), %lld MB (leaked), %lld MB (heap), %lld MB (max.rss), %lld MB (max.stack), %lld MB (vmpeak), %lld (paging)\n",
		  maxseen,leaked,hwm,rss,maxstack,vmpeak,pag);
	  fprintf(fp,"\tNo. of procs/threads: %d procs, %d threads\n",nproc,numthreads);
	}

	if (myproc == 1) {
	  fprintf(stderr,
		  "Memory-profiling information for program='%s', proc#%d:\n",a_out, myproc);
	  fprintf(stderr,"\tNo. of instrumented routines called : %d\n", numroutines);
	  fprintf(stderr,"\tInstrumentation started : %s\n",start_stamp ? start_stamp : "N/A");
	  fprintf(stderr,"\tInstrumentation   ended : %s\n",end_stamp ? end_stamp : "N/A");
	} /* if (myproc == 1) */

	free_drhook(end_stamp);

	fprintf(fp,"\n");
	len = 
	  fprintf(fp,"    #  Memory-%%      Self-alloc     + Children    Self-Leaked          Heap       Max.Stack     Paging     #Calls    #Allocs     #Frees   ");
                   /*"12345-1234567899-12345678901234-12345678901234-12345678901234-12345678901234-12345678901234-12345678901234-12345678901234-123456789012-123456789012"*/
	fprintf(fp,"Routine@<thread-id>");
	if (opt_clusterinfo) fprintf(fp," [Cluster:(id,size)]");
	fprintf(fp,"\n");
	fprintf(fp,  "         (self)         (bytes)        (bytes)        (bytes)        (bytes)        (bytes)    (delta)");
                   /*"12345-1234567899-12345678901234-12345678901234-12345678901234-12345678901234-12345678901234-12345678901234-12345678901234-123456789012-123456789012"*/
	fprintf(fp,"\n");

	p = prof;
	for (j=0; j<nprof; ) {
	  int cluster_size = clusize[p->cluster];
	  if (p->pc < percent_limit) break;
	  t = p->tid - 1;
	  if (p->children > maxseen_tot[t]) p->children = maxseen_tot[t]; /* adjust */
	  fprintf(fp, fmt,
		  ++j, p->pc, 
		  p->self, p->children, p->leaked,
		  p->hwm, p->stk, p->pag,
		  p->calls, p->alloc_count, 
		  (p->alloc_count - p->free_count != 0) ? "*" : " ", p->free_count,
		  p->is_max ? "*" : " ");

	  print_routine_name(fp, p, len, cluster_size);

	  fprintf(fp,"\n");
	  p++;
	} /* for (j=0; j<nprof; ) */
	
	fclose(fp);
      finish_4:
	free_drhook(filename);
	free_drhook(maxval);
	free_drhook(clusize);
      } while (0);

      free_drhook(tot);
      free_drhook(maxseen_tot);
      free_drhook(prof);
      do_prof_off = 0;
    }
  }
}

/*=== c_drhook_init_signals_ ===*/

void
c_drhook_init_signals_(const int *enforce)
{
  signal_drhook_init(*enforce);
}

/*=== c_drhook_raise_ ===*/

/* 
   Just a convenience function for Fortran90 which may not have raise()-signal function
   CALL c_drhook_raise(10)  ! Raise signal#10
*/

void 
c_drhook_raise_(const int *sig) 
{ 
  fflush(NULL);
  raise(*sig);
} 

/**** C-interface to Dr.Hook ****/

void
Dr_Hook(const char *name, int option, double *handle, 
	const char *filename, int sizeinfo,
	int name_len, int filename_len)
{
  static int first_time = 1;
  static int value = 1; /* ON by default */
  if (first_time) { /* Not thread safe */
    extern void *cdrhookinit_(int *value); /* from ifsaux/support/cdrhookinit.F90 */
    cdrhookinit_(&value);
    first_time = 0;
  }
  if (value == 0) return; /* Immediate return if OFF */
  if (value != 0) {
    int tid = get_thread_id_();
    if (option == 0) {
      c_drhook_start_(name, &tid, handle, 
		      filename, &sizeinfo,
		      name_len > 0 ? name_len : strlen(name),
		      filename_len > 0 ? filename_len : strlen(filename));
    }
    else if (option == 1) {
      c_drhook_end_(name, &tid, handle, 
		    filename, &sizeinfo,
		    name_len > 0 ? name_len : strlen(name),
		    filename_len > 0 ? filename_len : strlen(filename));
    }
  }
}


/**** Interface to HPM ****/

/*<<< experimental >>>*/

#ifdef HPM

#ifdef RS6K
/**** Interface to HPM (RS6K) ****/

#include <pmapi.h>

static pthread_mutex_t hpm_lock = PTHREAD_MUTEX_INITIALIZER;

static int *hpm_tid_init = NULL;
static double cycles = 1300000000.0; /* 1.3GHz ; changed via pm_cycles() in init_hpm() */

#define MCYCLES (cycles * 1e-6)

#define TEST_PM_ERROR(name, rc) \
  if (rc != 0) { \
    fprintf(stderr,"PM_ERROR(tid#%d, pthread_self()=%d): rc=%d at %s(), line=%d, file=%s\n",\
	    tid,pthread_self(),rc,name,__LINE__,__FILE__); \
    pm_error((char *)name, rc); \
    spin(tid); \
    RAISE(SIGABRT); \
  }

static void
init_hpm(int tid)
{
  const char *name = "init_hpm";
  int rc;

  if (!hpm_tid_init) {
    hpm_tid_init = calloc_drhook(numthreads, sizeof(*hpm_tid_init));
    cycles = pm_cycles();
  }

  if (!hpm_tid_init[tid-1]) {
#ifdef PMAPI_POST_P4
    pm_info2_t pminfo;
#else
    pm_info_t pminfo;
#endif
    pm_groups_info_t pmgroupsinfo;
    
    /*------------------------------------*/
    /* initialize the performance monitor */
    /*------------------------------------*/
#ifdef PMAPI_POST_P4
    rc = pm_initialize(PM_VERIFIED | PM_UNVERIFIED | PM_CAVEAT | PM_GET_GROUPS, 
		 &pminfo, &pmgroupsinfo, PM_CURRENT);
#else
    rc = pm_init(PM_VERIFIED | PM_UNVERIFIED | PM_CAVEAT | PM_GET_GROUPS, 
		 &pminfo, &pmgroupsinfo);
#endif
    TEST_PM_ERROR((char *)name, rc);

    if (myproc <= 1) fprintf(stderr,
			     ">>>pm_init() for ECMWF/OpenMP-tid#%d, pthread_self()=%d\n",
			     tid,pthread_self());
  }

  if (!hpm_tid_init[tid-1]) {
#if defined(PMAPI_P7)
    char *env = getenv("HPM_GROUP");
    hpm_grp = atoi(env);
    int group;
    fprintf(stderr,"hpm_group = %d\n",hpm_grp);
    if (hpm_grp == 150) group = 150; 
    if (hpm_grp == 141) group = 141; 
    /*-- counters --
     case 150:
       strcpy(group_label, "pm_vsu23, VSU Execution");
       strcpy(label[0], "four flops operation (fdiv,fsqrt) Scalar Instructions only (PM_VSU_FSQRT_FDIV)");
       strcpy(label[1], "VSU0 Finished an instruction (PM_VSU_FIN)");
       strcpy(label[2], "two flops operation (fmadd, fnmadd, fmsub, fnmsub) Scalar instructions only (PM_VSU_FMA)");
       strcpy(label[3], "one flop (fadd, fmul, fsub, fcmp, fsel, fabs, fnabs, fres, fsqrte, fneg) operation finished (PM_VSU_1FLOP)");
       strcpy(label[4], "Run instructions completed(PM_RUN_INST_CMPL)");
       strcpy(label[5], "Run cycles (PM_RUN_CYC)");
       strcpy(label[6], "Nothing");
       strcpy(label[7], "Nothing");
    */
    /*-- counters --
     case 141:
       strcpy(group_label, "pm_vsu14, VSU Execution");
       strcpy(label[0], "one flop (fadd, fmul, fsub, fcmp, fsel, fabs, fnabs, fres, fsqrte, fneg) operation finished (PM_VSU_1FLOP)");
       strcpy(label[1], "four flops operation (scalar fdiv, fsqrt; DP vector version of fmadd, fnmadd, fmsub, SP vector versions of single flop instructions) (PM_VSU_4FLOP)");
       strcpy(label[2], "eight flops operation (DP vector versions of fdiv,fsqrt and SP vector versions of fmadd,fnmadd,fmsub,fnmsub) (PM_VSU_8FLOP)");
       strcpy(label[3], "two flops operation (scalar fmadd, fnmadd, fmsub, fnmsub and DP vector versions of single flop instructions) (PM_VSU_2FLOP)");
       strcpy(label[4], "Run instructions completed(PM_RUN_INST_CMPL)");
       strcpy(label[5], "Run cycles (PM_RUN_CYC)");
       strcpy(label[6], "Nothing");
       strcpy(label[7], "Nothing");
    */
#elif defined(PMAPI_P6)
    const int group = 186; /* pm_hpm1 */
    /*-- counters --
     case 186:
       strcpy(group_label, "HPM group");
       strcpy(label[0], "FPU executed one flop instruction (PM_FPU_1FLOP)");
       strcpy(label[1], "FPU executed multiply-add instruction (PM_FPU_FMA)");
       strcpy(label[2], "FPU executed FSQRT or FDIV instruction (PM_FPU_SQRT_FDIV)");
       strcpy(label[3], "Processor Cycles (PM_CYC [shared chip])");
       strcpy(label[4], "Run instructions completed(PM_RUN_INST_CMPL)");
       strcpy(label[5], "Run cycles (PM_RUN_CYC)");
       strcpy(label[6], "Nothing");
       strcpy(label[7], "Nothing");
    */
#elif defined(PMAPI_P5_PLUS)
    /* IBM Power 5+ specific */
    const int group = 150; /* pm_hpmcount2 */
    /*-- counters -- (from John Hague, IBM/UK, 22-Aug-2006 : Thanx!!)
     case 150:
       strcpy(group_label, "pm_flop, Floating point operations");
       strcpy(label[0], "FPU executed FDIV instruction (PM_FPU_FDIV)");
       strcpy(label[1], "FPU executed multiply-add instruction (PM_FPU_FMA)");
       strcpy(label[2], "FPU executed FSQRT instruction (PM_FPU_SQRT)");
       strcpy(label[3], "FPU executed one flop instruction (PM_FPU_1FLOP)");
       strcpy(label[4], "Run instructions completed(PM_RUN_INST_CMPL)");
       strcpy(label[5], "Run cycles (PM_RUN_CYC)");
       strcpy(label[6], "Nothing");
       strcpy(label[7], "Nothing");
    */
#else
    const int group = 60; /* pm_hpmcount2 */
    /*-- counters --
     case 60:
       strcpy(group_label, "pm_hpmcount2, Hpmcount group for computation intensity analysis");
       strcpy(label[0], "FPU executed FDIV instruction (PM_FPU_FDIV)");
       strcpy(label[1], "FPU executed multiply-add instruction (PM_FPU_FMA)");
       strcpy(label[2], "FPU0 produced a result (PM_FPU0_FIN)");
       strcpy(label[3], "FPU1 produced a result (PM_FPU1_FIN)");
       strcpy(label[4], "Processor cycles (PM_CYC)");
       strcpy(label[5], "FPU executed store instruction (PM_FPU_STF)");
       strcpy(label[6], "Instructions completed (PM_INST_CMPL)");
       strcpy(label[7], "LSU executed Floating Point load instruction (PM_LSU_LDF)");
    */
#endif

    if (myproc <= 1) fprintf(stderr,"group = %d\n",group);

    pm_prog_t pmprog;
    pm_data_t pmdata;
    int i;

    /*---------------------*/
    /* set a default group */
    /*---------------------*/
    for (i=0; i<MAX_COUNTERS; i++) {
      pmprog.events[i] = COUNT_NOTHING;
    }
    pmprog.events[0] = group;
    
    /*-------------------------------------------------------------*/
    /* set the mode for user (not kernel) and thread (not process) */
    /*-------------------------------------------------------------*/
    pmprog.mode.w = 0;
    pmprog.mode.b.user = 1;
    pmprog.mode.b.process = 0;
    /* pmprog.mode.b.process = 1; */
    
    /*------------------------------------------*/
    /* for power-4 you have to use event groups */
    /*------------------------------------------*/
    pmprog.mode.b.is_group = 1;
    
    /*---------------------------------------------------*/
    /* set the mode to not to start counting immediately */
    /*---------------------------------------------------*/
    /* pmprog.mode.b.count = 1; */
    pmprog.mode.b.count = 0;
    
    /*-----------------------------------------*/
    /* initialize the group and start counting */
    /*-----------------------------------------*/
    hpm_tid_init[tid-1] = pthread_self(); /* Always > 0 */

    rc = pm_set_program_mythread(&pmprog); 
    TEST_PM_ERROR((char *)name, rc);

    rc = pm_start_mythread();
    TEST_PM_ERROR((char *)name, rc);
  }
}

static void
stop_only_hpm(int tid, drhook_key_t *pstop) 
{
  const char *name = "stop_only_hpm";
  pm_data_t pmdata;
  int i, rc;

  /* if (numthreads > 1) pthread_mutex_lock(&hpm_lock); */

  if (!hpm_tid_init || !hpm_tid_init[tid-1]) init_hpm(tid);

  /*
  rc = pm_stop_mythread();
  TEST_PM_ERROR((char *)name, rc);
  */

  if (pstop && !pstop->counter_stopped) {
    rc = pm_get_data_mythread(&pmdata);
    TEST_PM_ERROR((char *)name, rc);
    
    if (pstop && pstop->counter_in && !pstop->counter_stopped) {
      for (i=0; i<MAX_COUNTERS; i++) {
	pstop->counter_sum[i] += (pmdata.accu[i] - pstop->counter_in[i]);
      }
      pstop->counter_stopped = 1;
    }
  }

  /*
  rc = pm_start_mythread();
  TEST_PM_ERROR((char *)name, rc);
  */

  /* if (numthreads > 1) pthread_mutex_unlock(&hpm_lock); */
}

static void
stopstart_hpm(int tid, drhook_key_t *pstop, drhook_key_t *pstart)
{
  const char *name = "stopstart_hpm";
  pm_data_t pmdata;
  int i, rc;

  /* if (numthreads > 1) pthread_mutex_lock(&hpm_lock); */

  if (!hpm_tid_init || !hpm_tid_init[tid-1]) init_hpm(tid);

  /*
  rc = pm_stop_mythread();
  TEST_PM_ERROR((char *)name, rc);
  */

  rc = pm_get_data_mythread(&pmdata);
  TEST_PM_ERROR((char *)name, rc);

  if (pstop && pstop->counter_in && !pstop->counter_stopped) {
    for (i=0; i<MAX_COUNTERS; i++) {
      pstop->counter_sum[i] += (pmdata.accu[i] - pstop->counter_in[i]);
    }
    pstop->counter_stopped = 1;
  }

  if (pstart) {
    if (!pstart->counter_in ) pstart->counter_in  = calloc_drhook(MAX_COUNTERS, sizeof(*pstart->counter_in ));
    if (!pstart->counter_sum) pstart->counter_sum = calloc_drhook(MAX_COUNTERS, sizeof(*pstart->counter_sum));
     for (i=0; i<MAX_COUNTERS; i++) {
       pstart->counter_in[i] = pmdata.accu[i];
     }
     pstart->counter_stopped = 0;
  }

  /*
  rc = pm_start_mythread();
  TEST_PM_ERROR((char *)name, rc);
  */

  /* if (numthreads > 1) pthread_mutex_unlock(&hpm_lock); */
}

#else

/**** Interface to HPM (CRAY SV2, XD1 and XT3) ****/

static int *hpm_tid_init = NULL;
static double cycles = 0;

#define MCYCLES (cycles * 1e-6)

#define TEST_PM_ERROR(name, rc) \
  if (rc != 0) { \
    fprintf(stderr,"PM_ERROR(tid#%d, pthread_self()=%d): rc=%d at %s(), line=%d, file=%s\n",\
            tid,pthread_self(),rc,name,__LINE__,__FILE__); \
    pm_error((char *)name, rc); \
    spin(tid); \
    RAISE(SIGABRT); \
  }

static void
init_hpm(int tid)
{
  const char *name = "init_hpm";
  int rc;

  cycles = irtc_rate_();
}

static void
stop_only_hpm(int tid, drhook_key_t *pstop)
{
  const char *name = "stop_only_hpm";
  int i, rc;

  if (!hpm_tid_init || !hpm_tid_init[tid-1]) init_hpm(tid);

  if (pstop && !pstop->counter_stopped) {

    if (pstop && pstop->counter_in && !pstop->counter_stopped) {
#if defined(DT_FLOP)
      pstop->counter_sum[0] += ((long long int) flop_() - pstop->counter_in[0]);
#if defined(SV2)
      pstop->counter_sum[ENTRY_4] += (_rtc() - pstop->counter_in[ENTRY_4]);
#else
      pstop->counter_sum[ENTRY_4] += (irtc_() - pstop->counter_in[ENTRY_4]);
#endif
#endif
      pstop->counter_stopped = 1;
    }
  }
}


static void
stopstart_hpm(int tid, drhook_key_t *pstop, drhook_key_t *pstart)
{
  const char *name = "stopstart_hpm";
  int i, rc;

  if (!hpm_tid_init || !hpm_tid_init[tid-1]) init_hpm(tid);

  if (pstop && pstop->counter_in && !pstop->counter_stopped) {
#if defined(DT_FLOP)
      pstop->counter_sum[0] += ((long long int) flop_() - pstop->counter_in[0]);
#if defined(SV2)
      pstop->counter_sum[ENTRY_4] += (_rtc() - pstop->counter_in[ENTRY_4]);
#else
      pstop->counter_sum[ENTRY_4] += (irtc_() - pstop->counter_in[ENTRY_4]);
#endif
#endif
    pstop->counter_stopped = 1;
  }

  if (pstart) {
    if (!pstart->counter_in ) pstart->counter_in  = calloc_drhook(MAX_COUNTERS, sizeof(*pstart->counter_in ));
    if (!pstart->counter_sum) pstart->counter_sum = calloc_drhook(MAX_COUNTERS, sizeof(*pstart->counter_sum));
#if defined(DT_FLOP)
      pstart->counter_in[0] = (long long int) flop_();
#if defined(SV2)
      pstart->counter_in[ENTRY_4] = _rtc();
#else
      pstart->counter_in[ENTRY_4] = irtc_();
#endif
#endif
     pstart->counter_stopped = 0;
  }
}

#endif /*Interface to RS6K and SV2, XD1, XT3 */

static double
mflops_hpm(const drhook_key_t *keyptr)
{
  double mflops = 0;
  if (keyptr && keyptr->counter_sum && keyptr->counter_sum[ENTRY_4] > 0) {
    long long int sum = 0;
#if defined(DT_FLOP)
    sum = keyptr->counter_sum[0];
#elif defined(PMAPI_P7)
    /* IBM Power 7 specific */
    if(hpm_grp == 150) { 
      sum = 2 * keyptr->counter_sum[2] + keyptr->counter_sum[3];
    }
    if(hpm_grp == 141) { 
      sum = 2 * keyptr->counter_sum[0] + 4 * keyptr->counter_sum[1] + 2 * keyptr->counter_sum[3];
    }
#elif defined(PMAPI_P6)
    /* IBM Power 6 specific */
    sum = keyptr->counter_sum[0] + 2 * keyptr->counter_sum[1];
#elif defined(PMAPI_P5_PLUS)
    /* IBM Power 5+ specific */
    sum = 2 * keyptr->counter_sum[1] + keyptr->counter_sum[3];
#else
    sum = keyptr->counter_sum[1] + keyptr->counter_sum[2] + keyptr->counter_sum[3] - keyptr->counter_sum[5];
#endif
    if (sum > 0)
      mflops = (sum * MCYCLES)/keyptr->counter_sum[ENTRY_4];
  }
  return mflops;
}

static double
mips_hpm(const drhook_key_t *keyptr)
{
  double mipsrate = 0;
#if defined(DT_FLOP)
  mipsrate = 0;
#else
  if (keyptr && keyptr->counter_sum && keyptr->counter_sum[ENTRY_4] > 0) {
    mipsrate = (keyptr->counter_sum[ENTRY_6] * MCYCLES)/keyptr->counter_sum[ENTRY_4];
  }
#endif
  return mipsrate;
}

static double
divpc_hpm(const drhook_key_t *keyptr)
{
  double divpc = 0;
#if defined(DT_FLOP)
  divpc = 0;
#else
  if (keyptr && keyptr->counter_sum) {
    long long int sum = 0;
#if defined(PMAPI_P7)
    /* IBM Power 7 specific */
    if(hpm_grp == 150) { 
      sum = 2 * keyptr->counter_sum[2] + keyptr->counter_sum[3];
      if (sum > 0) divpc = (keyptr->counter_sum[0]*100.0)/sum;
    }
    if(hpm_grp == 141) { 
      sum = 2 * keyptr->counter_sum[0] + 4 * keyptr->counter_sum[1] + 2 * keyptr->counter_sum[3];
      if (sum > 0) divpc = (keyptr->counter_sum[1]*100.0)/sum;
    }
#elif defined(PMAPI_P6)
    /* IBM Power 6 specific */
    sum = keyptr->counter_sum[0] + 2 * keyptr->counter_sum[1];
    if (sum > 0) divpc = (keyptr->counter_sum[2]*100.0)/sum;
#elif defined(PMAPI_P5_PLUS)
    /* IBM Power 5+ specific */
    sum = 2 * keyptr->counter_sum[1] + keyptr->counter_sum[3];
    if (sum > 0) divpc = (keyptr->counter_sum[0]*100.0)/sum;
#else
    sum = keyptr->counter_sum[1] + keyptr->counter_sum[2] + keyptr->counter_sum[3] - keyptr->counter_sum[5];
    if (sum > 0) divpc = (keyptr->counter_sum[0]*100.0)/sum;
#endif
  }
#endif
  return divpc;
}

static double
mflop_count(const drhook_key_t *keyptr)
{
  double sum = 0;
  if (keyptr && keyptr->counter_sum && keyptr->counter_sum[ENTRY_4] > 0) {
#if defined(DT_FLOP)
    sum = (keyptr->counter_sum[0]) * 1e-6;
#elif defined(PMAPI_P7)
    /* IBM Power 7 specific */
    if(hpm_grp == 150) { 
      sum = (2 * keyptr->counter_sum[2] + keyptr->counter_sum[3]) * 1e-6;
    }
    if(hpm_grp == 141) { 
      sum = (2 * keyptr->counter_sum[0] + 4 * keyptr->counter_sum[1] + 2 * keyptr->counter_sum[3]) * 1e-6;
    }
#elif defined(PMAPI_P6)
    /* IBM Power 6 specific */
    sum = (keyptr->counter_sum[0] + 2 * keyptr->counter_sum[1]) * 1e-6;
#elif defined(PMAPI_P5_PLUS)
    /* IBM Power 5+ specific */
    sum = (2 * keyptr->counter_sum[1] + keyptr->counter_sum[3]) * 1e-6;
#else
    sum = (keyptr->counter_sum[1] + keyptr->counter_sum[2] + keyptr->counter_sum[3] - keyptr->counter_sum[5]) * 1e-6;
#endif
    if (sum < 0) sum = 0;
  }
  return sum;
}

static double
mip_count(const drhook_key_t *keyptr)
{
  double sum = 0;
#if defined(DT_FLOP)
  sum = 0;
#else
  if (keyptr && keyptr->counter_sum && keyptr->counter_sum[ENTRY_4] > 0) {
    sum = keyptr->counter_sum[ENTRY_6] * 1e-6;
  }
#endif
  return sum;
}

#endif /* HPM */


/* 
   this is result of moving some code from libodb.a
   (odb/aux/util_ccode.c) for use by libifsaux.a
   directly ; simplifies linking sequences.
*/

#include <stdio.h>
#include <string.h>
/* #include <malloc.h> */
#include <stdlib.h>
#include <signal.h>

#define FORTRAN_CALL

#if defined(CRAY) && !defined(SV2)
#define util_cputime_  UTIL_CPUTIME
#define util_walltime_ UTIL_WALLTIME
#endif

/* Portable CPU-timer (User + Sys) ; also WALL CLOCK-timer */

#include <unistd.h>
#include <sys/types.h>
#include <sys/times.h>
#undef MIN
#undef MAX
#include <sys/param.h>

#include <sys/time.h>

#if !defined(VPP)

FORTRAN_CALL
double util_walltime_()
{
  static double time_init = -1;
  double time_in_secs;
#if !defined(CRAYXT)
  struct timeval tbuf;
  if (gettimeofday(&tbuf,NULL) == -1) perror("UTIL_WALLTIME");

  if (time_init == -1) time_init = 
    (double) tbuf.tv_sec + (tbuf.tv_usec / 1000000.0);

  time_in_secs = 
  (double) tbuf.tv_sec + (tbuf.tv_usec / 1000000.0) - time_init;
#else
  if (time_init == -1) time_init = dclock();
  time_in_secs = dclock() - time_init;
#endif
  return time_in_secs;
}

#if defined(CRAYXT)
/* Cray XT3/XT4 with catamount microkernel */

FORTRAN_CALL
double util_cputime_()
{
  return util_walltime_(); /* In absence of anything better */
}

#else

extern clock_t times (struct tms *buffer);

FORTRAN_CALL
double util_cputime_()
{
  struct tms tbuf;
  static int first_time = 1;
  static double clock_ticks = 0;

  (void) times(&tbuf);

  if (first_time) {
    clock_ticks = (double) sysconf(_SC_CLK_TCK);
    first_time = 0;
  }

  return (tbuf.tms_utime + tbuf.tms_stime +
          tbuf.tms_cutime + tbuf.tms_cstime) / clock_ticks; 
}
#endif

#else 
/* VPP */
FORTRAN_CALL
double util_walltime_() 
{
  double w, time_in_secs;
  static double wallref = 0;
  extern FORTRAN_CALL gettod_(double *);
  if (wallref == 0) gettod_(&wallref);
  gettod_(&w);
  time_in_secs = (w - wallref) * 0.000001;
  return time_in_secs;
}
#endif

#ifdef VPP

#include <sys/types.h>
#include <sys/param.h>
#include <sys/signal.h>
#include <sys/fault.h>
#include <sys/syscall.h>
#include <sys/procfs.h>
#include <sys/proc.h>
#include <fcntl.h>

static int fujitsu_getrusage(int who, struct rusage *rusage)
{
  int rc = -1;

  if (rusage) rusage->ru_maxrss = 0;

  if (who == RUSAGE_SELF && rusage) {
    static int maxrss =  0;
    static int oldpid = -1;
    static char procfile[20] = "";
    static char *pf = NULL;
    /* static prpsinfo_t ps; */
    static proc_t proc;
    int pid = getpid();
    static int fildes = -1;
    unsigned int size;

    if (oldpid != pid) {
      oldpid = pid;
      maxrss = 0;
      pf = NULL;
    }

    if (!pf) {
      sprintf(procfile,"/proc/%d",pid);
      pf = procfile;
      fildes = open(procfile, O_RDONLY);
    }

    if (fildes == -1) return rc;

    /*
    if (ioctl(fildes, PIOCPSINFO, &ps) == -1) {
      perror("ioctl@fujitsu_getrusage(PIOCPSINFO)");
      return rc;
    }
    */

    if (ioctl(fildes, PIOCGETPR, &proc) == -1) {
      perror("ioctl@fujitsu_getrusage(PIOCGETPR)");
      return rc;
    }

    size  = /* ps.pr_usevpmem + */ proc.p_brksize + proc.p_stksize;
    if (size > maxrss) maxrss = size;
    rusage->ru_maxrss = maxrss;

    /* close(fildes); */
    rc = 0;
  }
  return rc;
}
#endif /* VPP */

FORTRAN_CALL
int util_ihpstat_(int *option)
{
  int ret_value = 0;

#if defined(SGI) || defined(VPP)
  if (*option == 1) {
    struct rusage rusage;
#ifdef SGI
    int pagesize = 1024;
    getrusage(0, &rusage);
#endif
#ifdef VPP
    int pagesize = 1; /* getpagesize() */
    fujitsu_getrusage(0, &rusage);
#endif
#if defined(SV2)
    int pagesize = getpagesize();
    getrusage(0, &rusage);
#endif
#if defined(XT3)
    int pagesize = getpagesize();
    getrusage(0, &rusage);
#endif
#if defined(XD1)
    int pagesize = getpagesize();
    getrusage(0, &rusage);
#endif
    ret_value = (rusage.ru_maxrss * pagesize + 7) / 8; /* In 8 byte words */
  }
#endif /* SGI or VPP */

  return ret_value;
}


#ifndef __timer_t_defined
static void set_timed_kill()
{
  // Definition of timer_t, timer_create, timer_set
  //   is a POSIX extention, not available on e.g. Darwin
}
#else
static void set_timed_kill()
{
#if !defined MACOSX
  if (drhook_timed_kill) {
    const char delim[] = ", \t/";
    char *p, *s = strdup_drhook(drhook_timed_kill);
    p = strtok(s,delim);
    while (p) {
      int target_myproc, target_omptid, target_sig;
      double start_time;
      int nelems = sscanf(p,"%d:%d:%d:%lf",
			  &target_myproc, &target_omptid, &target_sig, &start_time);
      int ntids = 1;
      coml_get_max_threads_(&ntids);
      if (nelems == 4 && 
	  (target_myproc == myproc || target_myproc == -1) &&
	  (target_omptid == -1 || (target_omptid >= 1 && target_omptid <= ntids)) &&
	  (target_sig >= 1 && target_sig <= NSIG) &&
	  start_time > 0) {
#if 1
	{
	  extern void run_fortran_omp_parallel_ipfipipipdpstr_(const int *, 
		 void (*func)(const int *, const int *, const int *, const double *, const char *, int len),
			      const int *, const int *, const int *, const double *, const char *, int len);
	  run_fortran_omp_parallel_ipfipipipdpstr_(&ntids,set_killer_timer,
						   &ntids,&target_omptid,&target_sig,&start_time,p,strlen(p));
	}
#else
#pragma omp parallel num_threads(ntids)
	{
	  set_killer_timer(&ntids,&target_omptid,&target_sig,&start_time,p,strlen(p));
	}
#endif
      }
      p = strtok(NULL,delim);
    }
    free_drhook(s);
  }
#endif
}
#endif
