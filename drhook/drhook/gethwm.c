#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include "getstatm.h"

typedef  long long int  ll_t;

static ll_t maxhwm = 0;

#if defined(CRAY) && !defined(SV2)
#define gethwm GETHWM
#else
#define gethwm gethwm_
#endif

#ifdef RS6K

#if defined(__64BIT__)
/* Assume AIX >= 5.1 with 64-bit addressing */
#include <fcntl.h>
#include <sys/procfs.h>
ll_t
gethwm()
{
  static int fd = -9999;
  static char *heapbase = NULL;
  ll_t heapsize = 0;

  if (fd == -9999) {
    pstatus_t pstatus;
    char procfile[80];
    int pid = getpid();
    sprintf(procfile,"/proc/%d/status",pid);
    fd = open(procfile, O_RDONLY);
    if (read(fd, &pstatus, sizeof(pstatus)) == sizeof(pstatus)) {
      heapbase = (char *)pstatus.pr_brkbase;
      close(fd);
      fd = 0;
    }
  }

  if (fd == 0 && heapbase != NULL) {
    heapsize = (ll_t)((char *)sbrk(0) - heapbase);
  }

  return heapsize;
}

#else

ll_t
gethwm() 
{ 
  extern ll_t getrss_();
  return getrss_();
}

#endif /* defined(__64BIT__) */

#else  /* non-RS6K */

// Cray linker: if you intend to link with -hstd_alloc and use Cray C compiler, then compile this file with -DSTD_ALLOC too
#if !defined(STD_ALLOC) && (defined(_CRAYC) || defined(USE_TCMALLOC))
ll_t
gethwm()
{
  extern size_t get_tcmalloc_heap_size_();
  return get_tcmalloc_heap_size_();
}
#elif defined(LINUX)
static ll_t basesize = -1;
static size_t pagesize = 4096;
ll_t gethwm()
{
  struct statm sm;
  ll_t rc = 0;
  if (getstatm(&sm) == 0) {
    if (basesize < 0) { /* the very first time */
      basesize = sm.size;
      pagesize = getpagesize();
      if (pagesize <= 0) pagesize = 4096;
    }
    rc = (sm.size - basesize) * pagesize;
    if (rc > maxhwm) maxhwm = rc;
  }
  return rc;
}

#elif defined(NECSX)

ll_t
gethwm() 
{ 
  extern ll_t getrss_();
  return getrss_();
}

#else
ll_t gethwm()
{
  ll_t rc = (ll_t)((char *)sbrk(0) - (char *)0);
  return rc;
}
#endif

#endif

#if defined(SV2)
int getpid_()
{
  return getpid();
}

unsigned int sleep_(unsigned int seconds)
{
  return sleep(seconds);
}

#endif

ll_t getmaxhwm_()
{
  ll_t rc = gethwm_();
  if (rc > maxhwm) maxhwm = rc;
  return maxhwm;
}
