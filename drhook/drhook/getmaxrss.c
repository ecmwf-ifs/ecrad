typedef  long long int  ll_t;
#include <sys/resource.h>
#include <sys/time.h>
ll_t
getmaxrss_()
{
  const ll_t scaler = 1024; /* in kilobytes */
  ll_t rc = 0;
  struct rusage r;
  rc = getrusage(RUSAGE_SELF, &r);
  rc = (rc == 0) ? (ll_t) r.ru_maxrss * scaler : 0;
  return rc;
}
