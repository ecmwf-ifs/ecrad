/* 
   Wrappers to obtain information from tcmalloc 
   Peter Towers - June 2014
   Made it working with any compiler providing libtcmalloc_minimal.so was linked in
   Sami Saarinen - November 2017
*/

#include "malloc_extension_c.h"
#include <stdio.h>

#pragma weak MallocExtension_GetNumericProperty

size_t get_tcmalloc_heap_size_()
{
  size_t value=0;
  if (MallocExtension_GetNumericProperty) {
    MallocExtension_GetNumericProperty("generic.heap_size", &value);
  }
  return value;
}

size_t get_tcmalloc_current_allocated_bytes_()
{
  size_t value=0;
  if (MallocExtension_GetNumericProperty) {
    MallocExtension_GetNumericProperty("generic.current_allocated_bytes", &value);
  }
  return value;
}

size_t get_tcmalloc_pageheap_free_bytes_()
{
  size_t value=0;
  if (MallocExtension_GetNumericProperty) {
    MallocExtension_GetNumericProperty("tcmalloc.pageheap_free_bytes", &value);
  }
  return value;
}

size_t get_tcmalloc_pageheap_unmapped_bytes_()
{
  size_t value=0;
  if (MallocExtension_GetNumericProperty) {
    MallocExtension_GetNumericProperty("tcmalloc.pageheap_unmapped_bytes", &value);
  }
  return value;
}
