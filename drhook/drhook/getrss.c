#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include "getstatm.h"

typedef  long long int  ll_t;

#if defined(CRAY) && !defined(SV2)
#define getrss GETRSS
#else
#define getrss getrss_
#endif

#if defined(RS6K) || defined(SGI) || defined(NECSX)
#include <sys/resource.h>

#if defined(RS6K) && defined(__64BIT__)
#if defined(USE_GETPROCS) 
/* From http://www-941.ibm.com/collaboration/wiki/display/WikiPtype/ryo
   Also : Thanks to Oliver Treiber, ECMWF, 16-Jan-2007 for useful discussions !! */
#include <procinfo.h>
#include <sys/time.h>
extern int getprocs64 (struct procentry64 *ProcessBuffer, 
		       int ProcessSize, struct fdsinfo64 *FileBuffer, int FileSize, 
		       pid_t *IndexPointer, int Count);
/* By the way : getprocs() for 32-bit addressing may also be worth trying */
#endif
#endif

ll_t
getrss()
{
  const ll_t scaler = 1024; /* in kilobytes */
  ll_t rc = 0;
#if defined(__64BIT__)
#if defined(USE_GETPROCS)
  /* IBM RS6K with -DUSE_GETPROCS */
  struct procentry64 procs;
  pid_t mypid = getpid();
  rc = getprocs64(&procs, sizeof(procs), NULL, 0, &mypid, 1); 
  rc = (rc == 1) ? (ll_t) procs.pi_drss*4*scaler : 0;
#else
  struct rusage64 r;
  rc = getrusage64(RUSAGE_SELF, &r);
  rc = (rc == 0) ? (ll_t) r.ru_maxrss * scaler : 0;
#endif
#else
  struct rusage r;
  rc = getrusage(RUSAGE_SELF, &r);
  rc = (rc == 0) ? (ll_t) r.ru_maxrss * scaler : 0;
#endif
  return rc;
}

#else

#if defined(LINUX)
static ll_t basesize = -1;
static size_t pagesize = 4096;
ll_t getrss()
{
  struct statm sm;
  ll_t rc = 0;
  if (getstatm(&sm) == 0) {
    if (basesize < 0) { /* the very first time */
      basesize = sm.resident;
      pagesize = getpagesize();
      if (pagesize <= 0) pagesize = 4096;
    }
    rc = (sm.resident - basesize) * pagesize;
  }
  return rc;
}
#else
ll_t getrss()
{
  ll_t rc = (ll_t)((char *)sbrk(0) - (char *)0);
  return rc;
}
#endif

#endif

