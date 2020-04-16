#ifndef _RAISE_H_
#define _RAISE_H_

/* raise.h */

#include <stdio.h>
#include <string.h>
#include <signal.h>
#include <unistd.h>

extern void abor1fl_(const char *filename, const int *linenum, 
		     const char *s, 
		     int filenamelen, int slen);
extern void abor1_(const char *s, int slen);

#define ABOR1(txt) { const char *t = (txt); t ? abor1_(t, strlen(t)) : abor1_("",0); }

#define ABOR1FL(txt) { \
  const char *t = (txt); \
  int linenum=__LINE__; \
  t ? abor1fl_(__FILE__, &linenum,  t, sizeof(__FILE__)-1, strlen(t)) \
    : abor1fl_(__FILE__, &linenum, "", sizeof(__FILE__)-1, 0); \
  _exit(1); /* Should never end up here */ }

#define RAISE(x) { \
  if ((x) == SIGABRT) { \
    ABOR1FL("*** Fatal error; aborting (SIGABRT) ..."); \
    _exit(1); /* Should never end up here */ \
  } \
  else raise(x); \
}

#endif /* _RAISE_H_ */
