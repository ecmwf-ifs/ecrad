#include <stdio.h>
#include <stdlib.h>
#include <cxxabi.h>

extern "C" 
char *cxxdemangle(const char *mangled_name, int *status)
{
  int istat = 0;
#ifdef __NEC__
  // Where is libstdc++ on NEC Aurora ??
  char *demangled_name = NULL;
#else
  char *demangled_name = abi::__cxa_demangle(mangled_name, NULL, NULL, &istat);
#endif
  if (status) *status = istat;
  return demangled_name; // this must be free()'d by the user
}
