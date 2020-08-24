#if defined(RS6K) && defined(__64BIT__)
#include <pthread.h>
#endif
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <signal.h>
#include <unistd.h>
#include "intercept_alloc.h"
#include "drhook.h"
extern int drhook_memtrace;
#include "raise.h"

typedef  long long int  ll_t;

#define WORDLEN   ((ll_t)sizeof(ll_t))
#define RNDUP(i,n) (( ( (i) + (n) - 1 ) / (n) ) * (n))
#define TRUE_BYTES(x) ((x) + 4*WORDLEN) /* size, keyptr at start & padding at end */

#if defined(CRAY) && !defined(SV2)
#define getcurheap           GETCURHEAP
#define getmaxcurheap        GETMAXCURHEAP
#define getcurheap_thread    GETCURHEAP_THREAD
#define getmaxcurheap_thread GETMAXCURHEAP_THREAD
#define getmaxloc            GETMAXLOC
#define resetmaxloc          RESETMAXLOC
#define profile_heap_get     PROFILE_HEAP_GET
#else
#define getcurheap           getcurheap_
#define getmaxcurheap        getmaxcurheap_
#define getcurheap_thread    getcurheap_thread_
#define getmaxcurheap_thread getmaxcurheap_thread_
#define getmaxloc            getmaxloc_
#define resetmaxloc          resetmaxloc_
#define profile_heap_get     profile_heap_get_
#endif

#if !defined(NTHRDS)
#define NTHRDS 64   /* ***Note: A hardcoded max number of threads !!! */
#endif

#if !defined(CACHELINESIZE)
#define CACHELINESIZE 128 /* ***Note: A hardcoded cache line size in bytes !!! */
#endif

#define TIMES  0    /* 0=Normal, 1=timing locks v threads */

static ll_t maxloc = 0;    /* For stackcheck */
static ll_t begloc = 0;    /* For stackcheck */
static int heapcheck = 0;  /* Fro heapcheck */

extern int get_thread_id_(void);   /* ***Note: Returns YOMOML-value [1..max_threads] */
extern int get_proc_id_(void);  
extern int get_max_threads_(void);
extern ll_t getstk_();
extern ll_t gethwm_();
#ifdef RS6K
extern void xl__trbk_();
#elif defined(NECSX)
extern void necsx_trbk_(const char *msg, int msglen); /* from ../utilities/gentrbk.F90 */
extern void necsx_trbk_fl_(const char *msg, const char *filename, int *lineno,
			   int msglen, int filenamelen); /* from ../utilities/gentrbk.F90 */
#define xl__trbk_() { int lineno = __LINE__; necsx_trbk_fl_("Error", __FILE__, &lineno, 5, sizeof(__FILE__)-1); }
#else
#define xl__trbk_() /* do nothing */
#endif

#if defined(CLSZ_OPT)
#undef CLSZ_OPT
#endif

#if defined(INTERCEPT_ALLOC)

/* Intercepts F90 ALLOCATE/DEALLOCATE */

#if defined(RS6K) && defined(__64BIT__)
/* Assume AIX >= 5.1 with 64-bit addressing */

pthread_mutex_t getcurheap_lock = PTHREAD_MUTEX_INITIALIZER;
#define CLSZ_OPT 1

#define EC_free     __free
#define EC_malloc   __malloc
#define EC_calloc   __calloc
#define EC_realloc  __realloc
#define EC_strdup   __strdup

int EC_malloc_will_abort = 1;

#elif defined(NECSX)

#include <pthread.h>
pthread_mutex_t getcurheap_lock = PTHREAD_MUTEX_INITIALIZER;
#define CLSZ_OPT 1

extern int EC_malloc_will_abort; /* see fnecsx.c */

#else

#undef INTERCEPT_ALLOC

#endif

#endif /* defined(INTERCEPT_ALLOC) */


static ll_t maxcurheap = 0;

#if defined(CLSZ_OPT)

static int profile_heap = -1; /* Profiling:  -1 = UNDEF (initvalue), 0 = OFF (default), 1 = ON */
#define NPROFILE 9 /* byte ranges: 10**1 .. 10^9 */
static ll_t malloc_hits[NPROFILE+1]; /* +1 for ranges >= 10^9 */
static ll_t free_hits[NPROFILE+1];
static ll_t alloc_amount[NPROFILE+1];
static ll_t malloc_hits_thrd[NTHRDS][NPROFILE+1]; 
static ll_t free_hits_thrd[NTHRDS][NPROFILE+1];
static ll_t alloc_amount_thrd[NTHRDS][NPROFILE+1];

static ll_t curalloc = 0;

static struct {
  ll_t curalloca;
  ll_t maxcurheapa;
  char pad[CACHELINESIZE - 2*WORDLEN]; /* padding : e.g. 128 bytes - 2*8 bytes */
} clsz_opt[NTHRDS]; /* cachelinesize optimized --> less false sharing when running with OpenMP */

#define NANS_FILL ((ll_t)0x7FF7FFFF7FF7FFFF)
static int nans_fill = -1; /* NaNS fill:  -1 = First, 1 = ON */
static int free_error     = 0;
static int max_free_error = 10;

#ifdef RS6K

#if 1
/* Turn NODISCLAIM on i.e. do NOT call disclaim() [before free()], since things get very, very expensive */
#if !defined(NODISCLAIM)
#define NODISCLAIM 1
#endif
#endif

#else

#if defined(NODISCLAIM)
#undef NODISCLAIM
#endif
#define NODISCLAIM 1

#endif 

#ifdef RS6K

#if !defined(NODISCLAIM)
#include <sys/shm.h> /* for use by disclaim() function */
/*
   http://publib.boulder.ibm.com/infocenter/pseries/v5r3/index.jsp?topic=/com.ibm.aix.genprogc/doc/genprogc/sys_mem_alloc.htm
   extern int disclaim (char *Address,  unsigned int Length,  unsigned int Flag);

   Note: We could also use 'export MALLOCOPTIONS=disclaim', but a short test showed
   that wall clock time went up by factor of 4, and cpu-time by factor 10 on a multithreaded MPI-code !!

   Now use 'export EC_DISCLAIM_THRESHOLD=50' to disclaim space for (say) allocations >= 50 MBytes
*/
/* Try to disclaim only if size < 2^31 && size >= disclaim_threshold (see below) */
static double disclaim_threshold = 2048; /* In MBytes -- for convenience : max however < 2^31/1024/1024 */
static ll_t disclaim_threshold_test = 0;
static ll_t disclaim_threshold_limit = 2147483647;
#endif /* !defined(NODISCLAIM) */

#endif /* RS6K */

static void
Check_curalloc() /* Normally not called */
{
  const ll_t big = (ll_t) 1000000000000L; /* 1,000,000 million bytes */
  if (curalloc < 0 || curalloc > big) {
    fprintf(stderr,
	    "Check_curalloc(): curalloc has probably gone crazy : Attempt to allocate ==> %lld bytes\n",curalloc);
    xl__trbk_();
    RAISE(SIGABRT);
    _exit(1);  /* Just in case, but shouldn't end up here at all */
  }
}

static void
Profile_heap_init()
{
  if (profile_heap == -1) { /* First time */
    char *env = getenv("EC_PROFILE_HEAP");
    if (env) {
      profile_heap = atoi(env);
      if (profile_heap != 0) profile_heap = 1; /* OFF by export EC_PROFILE_HEAP=0 */
    }
    if (profile_heap != 1) profile_heap = 0;
#if !defined(NODISCLAIM)
    if (disclaim_threshold_test == 0) {
      env = getenv("EC_DISCLAIM_THRESHOLD"); /* In MBytes (can have decimals) */
      if (env) {
	double tmp = atof(env);
	if (tmp > 0 || tmp < 2048) disclaim_threshold = tmp;
      }
      {
	double tmp = disclaim_threshold * 1024 * 1024;
	disclaim_threshold_test = (ll_t) tmp;
      }
    }
#endif
  }
}

static void
Profile_heap_put(ll_t size, int is_malloc)
{
  if (profile_heap == -1) Profile_heap_init();
  if (profile_heap == 1) {
    int j;
    ll_t n = 1; /* initial byte range */
    ll_t *p = is_malloc ? malloc_hits : free_hits;
    for (j=0; j<NPROFILE; j++) { 
      n *= 10; /* increment byte range by 10X */
      /* BTW: Don't want log10() overhead here !! */
      if (size < n) { /* i.e. size < pow(10,j+1) */
	alloc_amount[j] += is_malloc ? size : -size;
	p[j]++;
	return;
      }
    }
    j = NPROFILE;
    alloc_amount[j] += is_malloc ? size : -size;
    p[j]++;
  } /* if (profile_heap == 1) */
}

static void
Profile_heap_put_thrd(ll_t size, int is_malloc, int it)
{
  if (it == 0 && profile_heap == -1) Profile_heap_init(); /* 1st thread only checks this */
  if (profile_heap == 1) {
    int j;
    ll_t n = 1; /* initial byte range */
    ll_t *p = is_malloc ? malloc_hits_thrd[it] : free_hits_thrd[it];
    for (j=0; j<NPROFILE; j++) { 
      n *= 10; /* increment byte range by 10X */
      /* BTW: Don't want log10() overhead here !! */
      if (size < n) { /* i.e. size < pow(10,j+1) */
	alloc_amount_thrd[it][j] += is_malloc ? size : -size;
	p[j]++;
	return;
      }
    }
    j = NPROFILE;
    alloc_amount_thrd[it][j] += is_malloc ? size : -size;
    p[j]++;
  } /* if (profile_heap == 1) */
}

#endif /* defined(CLSZ_OPT) */

void
profile_heap_get(ll_t val[], 
		 const int *Nval, 
		 const int *Icase,
		 int *nret)
     /* Fortran callable */
{
#if defined(INTERCEPT_ALLOC)
  int nval = *Nval;
  int icase = *Icase;
  int j, it, nt;
  if (nval < 0) nval = 0;
  if (nval > NPROFILE+1) nval = NPROFILE+1;
  nt = get_max_threads_();
  for (j=0; j<nval; j++) {
    free_hits[j] = 0;
    malloc_hits[j] = 0;
    alloc_amount[j] = 0;
  }
  for (it=0; it<nt; it++) {
    for (j=0; j<nval; j++) {
      free_hits[j] += free_hits_thrd[it][j];
      malloc_hits[j] += malloc_hits_thrd[it][j];
      alloc_amount[j] += alloc_amount_thrd[it][j];
    }
  }
  if (icase == 0) { /* free() hits */
    for (j=0; j<nval; j++) val[j] = free_hits[j];
  }
  else if (icase == 1) { /* malloc() hits */
    for (j=0; j<nval; j++) val[j] = malloc_hits[j];
  }
  else if (icase == 2) { /* outstanding allocs (malloc minus free) */
    for (j=0; j<nval; j++) val[j] = malloc_hits[j] - free_hits[j];
  }
  else if (icase == 3) { /* allocation amount left per range; in bytes */
    for (j=0; j<nval; j++) val[j] = alloc_amount[j];
  }
  else if (icase == 4) { /* average allocation chunk left per range; in bytes */
    for (j=0; j<nval; j++) {
      ll_t tmp = malloc_hits[j] - free_hits[j];
      val[j] = (tmp > 0) ? RNDUP(alloc_amount[j],tmp)/tmp : 0;
    }
  }
  else {
    nval = 0;
  }
  *nret = nval;
#else
  *nret = 0;
#endif /* defined(CLSZ_OPT) */
}


#if defined(CLSZ_OPT)

void EC_free(void *vptr)
{
  if (vptr) {
    ll_t *p = vptr;
    ll_t dummy = *--p; 
    ll_t adjsize = *--p;
    ll_t keyptr = *--p;
    ll_t true_bytes;
    int it;
    if (nans_fill == 1) {
      ll_t *q = vptr;
      ll_t nans  = NANS_FILL;
      ll_t j = adjsize/WORDLEN;
      if (q[j] != nans) {
        fprintf(stderr,"WARNING: NaNS at end of array overwritten with %e\n",q[j]);
        xl__trbk_();
        free_error++;
        if (free_error > max_free_error) {
	  fprintf(stderr,"ERROR: Too many NaNS overwrites at end of arrays\n");
	  xl__trbk_(); /* Oops !! */
	  RAISE(SIGABRT);
	  _exit(1); /* Just in case, but shouldn't end up here at all */
	}
      }
    }
    it=get_thread_id_()-1;
    true_bytes = -TRUE_BYTES(adjsize);
#if !defined(NODISCLAIM)
    if (it == 0) { /* Do on OpenMP thread#1 only */
      ll_t *addr = (ll_t *)p;
      /* ll_t tb = -true_bytes; */
      ll_t tb = addr[-1]; /* Since addr[-1] contains the true length in bytes */ 
      if (tb >= disclaim_threshold_test &&
	  tb <= disclaim_threshold_limit) {
	unsigned int len = (unsigned int)tb;
	/* Important : Call disclaim(p, ...) BEFORE free(p) */
	int rc = disclaim(p, len, ZERO_MEM);
	/* printf("disclaim: rc=%d, len=%u, addr=%x\n",rc,len,addr); */
      }
    }
#endif
    p[1] = 0; /* Bytes reserved reset */ 
    free(p);
    if (drhook_memtrace) pthread_mutex_lock(&getcurheap_lock);
    clsz_opt[it].curalloca += true_bytes; /* += since true_bytes is negative */
    if (drhook_memtrace) pthread_mutex_unlock(&getcurheap_lock);
    if (drhook_memtrace) { ++it; c_drhook_memcounter_(&it, &true_bytes, &keyptr); --it; }
    if (profile_heap != 0) Profile_heap_put_thrd(true_bytes, 0, it);
  }
}

void *EC_malloc(ll_t size)
{
  double *d = NULL;
  void *vptr = NULL;
  ll_t keyptr = 0;
  ll_t adjsize = size;
  ll_t thesize = size;
  ll_t true_bytes;
  int it;
  if (nans_fill == -1) { /* First time */
    char *env = getenv("EC_FILL_NANS");
    nans_fill=0;
    if (env) {
      if (strcmp(env,"true") == 0 || 
	  strcmp(env,"TRUE") == 0 ||
	  strcmp(env,"1"   ) == 0) nans_fill = 1;
      /* fprintf(stderr,"EC_FILL_NANS ==> env,env,nans_fill= 0x%x %s %d\n",env,env,nans_fill); */
    }
  }
  if (adjsize < 0) adjsize = 0;
  adjsize = RNDUP(adjsize,WORDLEN);
#if !defined(NECSX)
  thesize = adjsize;
#endif
  true_bytes = TRUE_BYTES(adjsize);
  it=get_thread_id_()-1;
  if (TIMES == 1) {
    DRHOOK_START(dummy);
    DRHOOK_END(0);
  }
  if (TIMES == 1) {
    DRHOOK_START(malloc);
    d = (double *)malloc(true_bytes); 
    DRHOOK_END(0);
  }
  else {
    d = (double *)malloc(true_bytes);
  }
  vptr = d;
  if (vptr) {
    if (TIMES == 1) {
      {
        ll_t *p = vptr;
        DRHOOK_START(lock);
        (void) getstk_(); /* to gather near up to date stack statistics */
	*p++ = keyptr;
        *p++ = thesize;
        *p++ = 0;
        pthread_mutex_lock(&getcurheap_lock);
        curalloc += true_bytes;
        if (curalloc > maxcurheap) maxcurheap = curalloc;
        /* Check_curalloc(); */
        Profile_heap_put(true_bytes, 1);
        pthread_mutex_unlock(&getcurheap_lock);
        vptr = p;
        DRHOOK_END(0);
      }
      {
        ll_t *p = vptr;
        DRHOOK_START(thread);
        (void) getstk_(); /* to gather near up to date stack statistics */
        if (drhook_memtrace) pthread_mutex_lock(&getcurheap_lock);
        clsz_opt[it].curalloca += true_bytes;
        if (clsz_opt[it].maxcurheapa < clsz_opt[it].curalloca) 
	    clsz_opt[it].maxcurheapa = clsz_opt[it].curalloca;
	if (drhook_memtrace) { ++it; c_drhook_memcounter_(&it, &true_bytes, &keyptr); --it; }
        if (drhook_memtrace) pthread_mutex_unlock(&getcurheap_lock);
        if (profile_heap != 0) Profile_heap_put_thrd(true_bytes, 1, it);
        DRHOOK_END(0);
      }
    }
    else {
      int ip;
      ll_t *p = vptr;
      ll_t q;
      (void) getstk_(); /* to gather near up to date stack statistics */
      *p++ = keyptr;
      *p++ = thesize;
      *p++ = 0;
      if (drhook_memtrace) pthread_mutex_lock(&getcurheap_lock);
      clsz_opt[it].curalloca += true_bytes;
      if (clsz_opt[it].maxcurheapa < clsz_opt[it].curalloca) 
	  clsz_opt[it].maxcurheapa = clsz_opt[it].curalloca;
      if (drhook_memtrace) { ++it; c_drhook_memcounter_(&it, &true_bytes, &keyptr); --it; }
      if (drhook_memtrace) pthread_mutex_unlock(&getcurheap_lock);
      if (profile_heap != 0) Profile_heap_put_thrd(true_bytes, 1, it);
      if (nans_fill == 1) {
        int j;
        ll_t nans  = NANS_FILL;
        for (j=0; j<(WORDLEN+adjsize)/WORDLEN; j++) {
          p[j]=nans;
        }
      }
      if(heapcheck == 1) {
        it=get_thread_id_();
        ip=get_proc_id_();
        if (it == 1 && ip == 1 ) {
          if (begloc == 0) begloc=(ll_t)p;
          q=(ll_t)p+true_bytes;
          if (q > maxloc) {
            maxloc=q;
/*
            fprintf(stderr,"JJJ pntr= %d %ld %ld %ld %ld %ld\n",it,maxloc-begloc,maxloc,begloc,p,true_bytes);
            xl__trbk_();
*/
          }
        }
      }
      vptr=p;
    }
  }
  else {
    pthread_mutex_lock(&getcurheap_lock);
    fprintf(stderr,
	    "EC_malloc(size=%lld => thesize=%lld => adjsize=%lld, true_bytes=%lld bytes) failed in file=%s, line=%d\n",
	    size, thesize, adjsize, true_bytes, __FILE__, __LINE__);
    xl__trbk_(); /* Oops !! */
    if (EC_malloc_will_abort) RAISE(SIGABRT);
    pthread_mutex_unlock(&getcurheap_lock);
    if (EC_malloc_will_abort) _exit(1); /* Just in case, but shouldn't end up here at all */
  }
  return vptr;
}

void *EC_calloc(ll_t nelem, ll_t elsize)
{
  ll_t totbytes = nelem * elsize;
  void *p = EC_malloc(totbytes);
  if (p) memset(p, 0, totbytes);
  return p;
}

void *EC_realloc(void *vptr, ll_t size)
{
  ll_t *pnew = NULL;
  if (vptr) {
    ll_t *p = vptr;
    ll_t oldsize = p[-1];
    if (oldsize < size) {
      pnew = EC_malloc(size);
      if (pnew) {
	memcpy(pnew, p, oldsize);
	EC_free(p);
      }
    }
    else { /* the old allocation size was already sufficient */
      pnew = p;
    }
  }
  else { /* Revert to malloc() */
    pnew = EC_malloc(size);
  }
  return pnew;
}

char *EC_strdup(const char *s)
{
  ll_t totbytes = sizeof(*s) * strlen(s);
  return EC_malloc(totbytes);
}

#endif /* defined(CLSZ_OPT) */

ll_t
getcurheap()
{
#if defined(CLSZ_OPT)
  int it = get_thread_id_();
  ll_t curvalue = clsz_opt[it-1].curalloca;
  if (it == 1) { /* Only thread#1 sums up */
    int i, nt = get_max_threads_();
    if (drhook_memtrace) pthread_mutex_lock(&getcurheap_lock);
    for (i=1; i<nt; i++) {
      curvalue += clsz_opt[i].curalloca;
    }
    if (drhook_memtrace) pthread_mutex_unlock(&getcurheap_lock);
  }
  return curvalue;
  // Cray linker: if you intend to link with -hstd_alloc and use Cray C compiler, then compile this file with -DSTD_ALLOC too
#elif !defined(STD_ALLOC) && (defined(_CRAYC) || defined(USE_TCMALLOC))
  extern size_t get_tcmalloc_current_allocated_bytes_();
  return get_tcmalloc_current_allocated_bytes_();
#else
  ll_t rc = gethwm_();
  if (rc > maxcurheap) maxcurheap = rc;
  return rc;
#endif
}

ll_t
getcurheap_thread(const int *thread_id)
{
#if defined(CLSZ_OPT)
  int it = (thread_id && (*thread_id > 0)) ? *thread_id : get_thread_id_();
  return clsz_opt[--it].curalloca;
#else
  return getcurheap();
#endif
}

/* Maximum (total) current (virtual mem) allocation encountered */

ll_t
getmaxcurheap()
{
#if defined(CLSZ_OPT)
  ll_t maxcurheap_local=0;
  int it, nt = get_max_threads_();
  if (drhook_memtrace) pthread_mutex_lock(&getcurheap_lock);
  for (it=0; it<nt; it++) {
    maxcurheap_local += clsz_opt[it].maxcurheapa;
  }
  if (maxcurheap < maxcurheap_local) maxcurheap = maxcurheap_local;
  if (drhook_memtrace) pthread_mutex_unlock(&getcurheap_lock);
  return maxcurheap_local;
#else
  return 0;
#endif
}

/* Maximum (total) current (virtual mem) allocation encountered per thread */

ll_t
getmaxcurheap_thread(const int *thread_id) /* ***Note: YOMOML thread id */
{
#if defined(CLSZ_OPT)
  int it = (thread_id && (*thread_id > 0)) ? *thread_id : get_thread_id_();
  return  clsz_opt[--it].maxcurheapa;
#else
  return 0;
#endif
}

#if 1
ll_t
getmaxloc()
{
  ll_t z=maxloc-begloc;
  return z;
}

void
resetmaxloc()
{
  maxloc=0;
}

void

setheapcheck_()
{
  heapcheck=1;
}

#endif

