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

/* Address difference in bytes */

void 
addrdiff_(const char *p1, const char *p2, int *diff)
{
  *diff = (p2 - p1);
}

/* loc()-function */

unsigned long long int
loc_addr_(const char *p)
{
  return (unsigned long long int)(p - (const char *)0);
}
