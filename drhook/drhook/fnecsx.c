/* fnecsx.c */

/* 

   This source file contains reverse-engineered functions that are called
   when ALLOCATE & DEALLOCATE are issued in Fortran90 codes on NEC SX.

   This facilitates memory monitoring and allocations in a similar way
   as we can do on IBM (see also getcurheap.c) and 
   provides unified interface to Dr.Hook with memory tracing/profiling.

   To have any effect, compile with sxcc :

   -DINTERCEPT_ALLOC
   -DNECSX


   Author: Sami Saarinen, ECMWF, 02-Feb-2007

*/

#include "intercept_alloc.h"

#if defined(INTERCEPT_ALLOC)
#if defined(NECSX)

#include "raise.h"
#include <stdlib.h>

typedef long long int ll_t;
typedef unsigned long long int u_ll_t;

/* Maximum no of dimensions allowed in any given Fortran90 array */

#define FNECSX_MAXDIM 7

/*
typedef struct {
  ll_t lo;
  ll_t hi;
  ll_t stride;
} dim_t ;
*/

#define FNECSX_NKEYS  3
#define FNECSX_KEY_LO 0
#define FNECSX_KEY_HI 1
#define FNECSX_KEY_ST 2

typedef struct {
  u_ll_t *p;
  int status;
  int ndims;
  ll_t s[FNECSX_MAXDIM][FNECSX_NKEYS];
} desc_t;

#if 0
/* Logical, f.ex.: if (ALLOCATED(allocatable_array)) ... */
extern ll_t fy_sallocd(desc_t *d);
#endif

/* For sxf90 -dw i.e. 4-byte int for istat in [DE]ALLOCATE(..., stat=istat) is required */
extern void f_alloc(ll_t arg1, desc_t *d, int *stat, ll_t elsize);
extern void f_deallc(ll_t arg1, desc_t *d, int *stat);

/* For sxf90 -ew i.e. 8-byte int for istat in [DE]ALLOCATE(..., stat=istat) is required */
extern void f_allocl(ll_t arg1, desc_t *d, ll_t *stat, ll_t elsize);
extern void f_deallcl(ll_t arg1, desc_t *d, ll_t *stat);

int EC_malloc_will_abort = 0; /* affects getcurheap.c when memory allocation fails in EC_malloc() */

extern void necsx_trbk_fl_(const char *msg, const char *filename, int *lineno,
			   int msglen, int filenamelen); /* from ../utilities/gentrbk.F90 */
#define ERROR_MSG(msg) { \
  int lineno = __LINE__; necsx_trbk_fl_(msg, __FILE__, &lineno, strlen(msg), sizeof(__FILE__)-1); }

#else
#undef INTERCEPT_ALLOC
#endif
#endif /* defined(INTERCEPT_ALLOC) */

#if !defined(INTERCEPT_ALLOC)

/* Other than NEC SX machines or when -DINTERCEPT_ALLOC was NOT supplied */

void ec_envredo_() { }

#else /* is indeed defined(INTERCEPT_ALLOC) */

/* NEC SX with -DINTERCEPT_ALLOC */

static int ec_prtdesc = -1;
static int ec_initheap = 0;
static unsigned int  ec_initval4 = 0;
static u_ll_t        ec_initval8 = 0;
static int ec_malloc = 1;

static
void init4(unsigned int p4[], ll_t n)
{
  u_ll_t tmp = ec_initval4;
  ll_t j;
  for (j=0; j<n; j++) p4[j] = tmp;
}

static
void init8(u_ll_t p8[], ll_t n)
{
  u_ll_t tmp = ec_initval8;
  ll_t j;
  for (j=0; j<n; j++) p8[j] = tmp;
}

static
void envinit()
{
  char *env;

  /* Print array descriptor info */
  env = getenv("EC_PRTDESC");
  if (env) {
    ec_prtdesc = atoi(env);
    if (ec_prtdesc != 0) ec_prtdesc = 1;
  }
  else {
    ec_prtdesc = 0;
  }

  /* Simulating the effect of "-init heap={zero|nan|0xXXXX}" */
  env = getenv("EC_INITHEAP");
  if (env) {
    int len = strlen(env);
    if (strcasecmp(env,"zero") == 0 || strcasecmp(env,"0") == 0) {
      ec_initheap = 1; /* 1-byte long */
      if (ec_prtdesc) fprintf(stderr,"EC_INITHEAP='%s' => ec_initheap = %d\n",env,ec_initheap);
    }
    else if (strcasecmp(env,"nan") == 0) {
      ec_initheap = 4; /* 4-bytes long */
      ec_initval4 = 0x7fffffff;
      if (ec_prtdesc) fprintf(stderr,"EC_INITHEAP='%s' => ec_initheap = %d : value = 0x%x (%u,%d)\n",
			      env,ec_initheap,ec_initval4,ec_initval4,ec_initval4);
    }
    else if (len >= 2 && env[0] == '0' && (env[1] == 'x' || env[1] == 'X')) {
      if (len <= 10) { /* 0x12345678 */
	ec_initheap = 4;
	sscanf(env,"0x%x",&ec_initval4);
	if (ec_prtdesc) fprintf(stderr,"EC_INITHEAP='%s' => ec_initheap = %d : value = 0x%x (%u,%d)\n",
				env,ec_initheap,ec_initval4,ec_initval4,ec_initval4);
      }
      else if (len > 10 && len <= 18) { /* 0x1234567890abcdef */
	ec_initheap = 8;
	sscanf(env,"0x%llx",&ec_initval8);
	if (ec_prtdesc) fprintf(stderr,"EC_INITHEAP='%s' => ec_initheap = %d : value = 0x%llx (%llu,%lld)\n",
				env,ec_initheap,ec_initval8,ec_initval8,ec_initval8);
      }
      else {
	ec_initheap = 0;
      }
    }
    else {
      ec_initheap = 0;
    }
  }
  else {
    ec_initheap = 0;
  }

  /* Use EC_malloc/EC_calloc/EC_free (default) or malloc/calloc/free */
  env = getenv("EC_MALLOC");
  if (env) {
    ec_malloc = atoi(env);
    if (ec_malloc != 0) ec_malloc = 1;
  }
  else {
    ec_malloc = 1;
  }
}

void ec_envredo_() { envinit(); };

static 
void prtdesc(ll_t arg1, FILE *fp, const char *s, const desc_t *d, 
	     const int *stat, ll_t elsize_in)
{
  if (ec_prtdesc == -1) {
    envinit();
    if (ec_prtdesc == 0) return;
  }
  if (fp && s && d) {
    int j;
    int ndims = d->ndims;
    ll_t elsize = d->s[0][FNECSX_KEY_ST];
    ll_t ntot = 1;
    ll_t total_bytes;
    ll_t nsave[FNECSX_MAXDIM];
#pragma cdir altcode,loopcnt=FNECSX_MAXDIM
    for (j=0; j<ndims; j++) {
      ll_t n = d->s[j][FNECSX_KEY_HI] - d->s[j][FNECSX_KEY_LO] + 1;
      ntot *= n;
      nsave[j] = n;
    }
    if (((ntot/2)*2) != ntot) ntot++; /* mod(ntot,2) not 0 i.e. ntot is odd --> add 1 */
    total_bytes = ntot * elsize + 1;
    fprintf(fp, "=== %s : desc = 0x%llx [%lld] : arg1 = %lld, stat (addr) = 0x%llx\n",s,d,d,arg1,stat);
    fprintf(fp,"p = 0x%llx [%lld] : status = %d, ndims = %d, elsize's = (%lld, %lld), total = %lld\n",
	    d->p, d->p, d->status, ndims, elsize, elsize_in, total_bytes);
    for (j=0; j<ndims; j++) {
      ll_t n = nsave[j];
      fprintf(fp,"[dim#%d] : lo = %lld, hi = %lld (n=%lld), stride = %lld\n",
	      j+1, d->s[j][FNECSX_KEY_LO], d->s[j][FNECSX_KEY_HI], n, 
	      d->s[j][FNECSX_KEY_ST]);
    }
    if (d->p) {
      u_ll_t dummy = d->p[-1];
      u_ll_t adjsize = d->p[-2];
      u_ll_t keyptr = d->p[-3];
      fprintf(fp,"\tdummy = %llu\n",dummy);
      fprintf(fp,"\tadjsize (bytes) = %llu\n",adjsize);
      fprintf(fp,"\tDr.Hook keyptr = 0x%llx [%llu]\n",keyptr,keyptr);
    }
  }
}

#if 0
ll_t fy_sallocd(desc_t *d)
{
  if (ec_prtdesc) prtdesc(-1, stderr, "fy_sallocd", d, NULL, -1);
  return (d && d->p && (d->status == 3 || d->status == 1)) ? 1 : 0;
}
#endif

void f_alloc(ll_t arg1, desc_t *d, int *stat, ll_t elsize)
{
  const int istat_offset = 193;
  int istat = 0;
  const char *errmsg[2] = {
    "Value of allocate-object must not be currently allocated array in ALLOCATE.", /* 193 */
    "Could not allocate in ALLOCATE." /* 194 */
  };
  if (stat && stat == (int *)0x1) stat = NULL;
  if (ec_prtdesc) prtdesc(arg1, stderr, "f_alloc>", d, stat, elsize);
  if (d && (arg1 == 3 || arg1 == 2)) { 
    /* arg1 == 3 : ALLOCATABLE array --> d->status becomes = 3
       arg1 == 2 : POINTER array     --> d->status becomes = 1 */
    int j;
    int ndims = d->ndims;
    u_ll_t *p = NULL;
    void *vptr = NULL;
    ll_t ntot = 1;
    ll_t total_bytes;
    ll_t nsave[FNECSX_MAXDIM];
#pragma cdir altcode,loopcnt=FNECSX_MAXDIM
    for (j=0; j<ndims; j++) {
      ll_t n = d->s[j][FNECSX_KEY_HI] - d->s[j][FNECSX_KEY_LO] + 1;
      ntot *= n;
      nsave[j] = n;
    }
    if (((ntot/2)*2) != ntot) ntot++; /* mod(ntot,2) not 0 i.e. ntot is odd --> add 1 */
    total_bytes = ntot * elsize + 1;
    if (ec_initheap == 1) {
      /* EC_INITHEAP=zero */
      if (ec_malloc) {
	vptr = EC_calloc(total_bytes,1);
      }
      else {
	vptr = calloc(total_bytes,1);
      }
    }
    else {
      if (ec_malloc) {
	vptr = EC_malloc(total_bytes);
      }
      else {
	vptr = malloc(total_bytes);
      }
      if (ec_initheap != 0) {
	if (ec_initheap == 4) {
	  /* EC_INITHEAP=nan or max an 8 digit hexadecimal number */
	  ll_t n = total_bytes/4;
	  init4(vptr, n);
	}
	else if (ec_initheap == 8) {
	  /* EC_INITHEAP=nan or a 9-16 digit hexadecimal number */
	  ll_t n = total_bytes/8;
	  init8(vptr, n);
	}
      } /* if (ec_initheap) */
    }
    p = vptr;
    if (p) {
      d->p = p;
      d->status = (arg1 == 3) ? 3 : 1;
      d->s[0][FNECSX_KEY_ST] = elsize;
      if (ndims > 1) {
#pragma cdir altcode,loopcnt=FNECSX_MAXDIM
	for (j=1; j<ndims; j++) {
	  d->s[j][FNECSX_KEY_ST] = nsave[j] * d->s[j-1][FNECSX_KEY_ST];
	}
      } /* if (ndims > 1) */
    }
    else {
      istat = 194;
    }
  }
  else {
    istat = 193;
  }
  if (ec_prtdesc) prtdesc(arg1, stderr, "f_alloc<", d, stat, elsize);
  if (istat != 0 && !stat) {
    const char *msg = errmsg[istat-istat_offset];
    fprintf(stderr,"***Error#%d: %s\n",istat,msg);
    prtdesc(arg1, stderr, "f_alloc", d, stat, elsize);
    ERROR_MSG(msg);
    RAISE(SIGABRT);
    _exit(1); /* Just in case, but shouldn't end up here at all */
  }
  if (stat) *stat = istat;
}


void f_allocl(ll_t arg1, desc_t *d, ll_t *stat, ll_t elsize)
{
  int istat = 0;
  if (stat && stat == (ll_t *)0x1) stat = NULL;
  f_alloc(arg1, d, stat ? &istat : NULL, elsize);
  if (stat) *stat = istat;
}


void f_deallc(ll_t arg1, desc_t *d, int *stat)
{
  const int istat_offset = 195;
  int istat = 0;
  const char *errmsg[2] = {
    "Value of allocate-object must not be disassociated pointer/not allocated array in DEALLOCATE.", /* 195 */
    "Illegal value of allocate-object in DEALLOCATE." /* 196 */
  };
  if (stat && stat == (int *)0x1) stat = NULL;
  if (ec_prtdesc) prtdesc(arg1, stderr, "f_deallc>", d, stat, -1);
  if (d && ((arg1 | d->status) == 3)) {
    int j;
    int ndims = d->ndims;
    u_ll_t *p = d->p;
    if (p) {
      if (ec_malloc) {
	EC_free(p);
      }
      else {
	free(p);
      }
      d->p = NULL;
      d->status = 0;
#pragma cdir altcode,loopcnt=FNECSX_MAXDIM
      for (j=0; j<ndims; j++) {
	d->s[j][FNECSX_KEY_ST] = 0;
      }
    }
    else {
      istat = 196;
    }
  }
  else {
    istat = 195;
  }
  if (ec_prtdesc) prtdesc(arg1, stderr, "f_deallc<", d, stat, -1);
  if (istat != 0 && !stat) {
    const char *msg = errmsg[istat-istat_offset];
    fprintf(stderr,"***Error#%d: %s\n",istat,msg);
    prtdesc(arg1, stderr, "f_deallc", d, stat, -1);
    ERROR_MSG(msg);
    RAISE(SIGABRT);
    _exit(1); /* Just in case, but shouldn't end up here at all */
  }
  if (stat) *stat = istat;
}


void f_deallcl(ll_t arg1, desc_t *d, ll_t *stat)
{
  int istat = 0;
  if (stat && stat == (ll_t *)0x1) stat = NULL;
  f_deallc(arg1, d, stat ? &istat : NULL);
  if (stat) *stat = istat;
}

#endif /* #if !defined(INTERCEPT_ALLOC) ... #else ... */

void ec_envredo() { ec_envredo_(); }
