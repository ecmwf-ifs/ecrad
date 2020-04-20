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

