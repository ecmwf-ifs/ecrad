/**
 * (C) Copyright 2014- ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 *
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

#if defined(CRAY) && !defined(SV2)
#define getpag GETPAG
#else
#define getpag getpag_
#endif

#if defined(RS6K) || defined(SV2)
#include <sys/resource.h>
long long int
getpag()
{
#if defined(__64BIT__)
  struct rusage64 r;
  long long int rc = getrusage64(RUSAGE_SELF, &r);
#else
  struct rusage r;
  long long int rc = getrusage(RUSAGE_SELF, &r);
#endif
  rc = (rc == 0) ? (long long int) r.ru_majflt : 0;
  return rc;
}

#else


long long int getpag()
{
  return 0L;
}

#endif

