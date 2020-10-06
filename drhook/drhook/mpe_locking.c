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

#ifdef VPP

void mpe_lock_(int *lockid, int *status)
{
  int LockID = (*lockid) + 1;
  *status = VPP_SemWait(0,LockID);
}

void mpe_unlock_(int *lockid,int *status)
{
  int LockID = (*lockid) + 1;
  *status = VPP_SemPost(0,LockID);
}

#else

void mpe_lock_(int *lockid, int *status)
{
  *status = 0;
}

void mpe_unlock_(int *lockid,int *status)
{
  *status = 0;
}
#endif

