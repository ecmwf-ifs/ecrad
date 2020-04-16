typedef long long int ll_t;
typedef unsigned long long int ull_t;

#ifdef VPP
ll_t getstackusage_dummy_() { return 0L; }
#else

#if defined(LINUX) && defined(USE_MEMORY_MONITOR)
#include <stdio.h>
#include <stdlib.h>

ll_t getstackusage_()
{
  ll_t rc = 0;
  static int dont_bother = 0;
  if (dont_bother) {
    rc = -2;
  }
  else {
    FILE *statfile = fopen ("/proc/self/stat", "r");
    if (!statfile) {
      dont_bother = 1;
      rc = -1;
    }
    else {
      char dm[80];
      ull_t startstack, kstkesp; /* stack start & ESP, the 28th and 29th columns, respectively */
      /* Maybe not the brightest coding, but has to suffice for now (SS) */
      int nelem = fscanf(statfile, 
		     "%s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %llu %llu ",
		     dm,dm,dm,dm,dm,dm,dm,dm,dm,dm,dm,dm,dm,dm,dm,dm,dm,dm,dm,dm,dm,dm,dm,dm,dm,dm,dm,&startstack,&kstkesp);
      if (nelem != 29) {
	dont_bother = 1;
	rc = -3;
      }
      else {
	rc = (ll_t)(startstack - kstkesp);
      }
      fclose(statfile);
    }
  }
  return rc;
}

#else

ll_t getstackusage_() { return 0L; }

#endif
#endif