#ifndef _INTERCEPT_ALLOC_H_
#define _INTERCEPT_ALLOC_H_

/* intercept_alloc.h */

#if defined(INTERCEPT_ALLOC)

#if defined(RS6K) && defined(__64BIT__)

#define EC_free     __free
#define EC_malloc   __malloc
#define EC_calloc   __calloc
#define EC_realloc  __realloc
#define EC_strdup   __strdup

#elif defined(NECSX)
/* Do nothing */
#else
/* Illegal to have -DINTERCEPT_ALLOC */
#undef INTERCEPT_ALLOC
#endif

#endif

#if defined(INTERCEPT_ALLOC)

/* For reference, see also ifsaux/utilities/getcurheap.c */

#define THEmalloc  EC_malloc
extern void *EC_malloc(long long int size);
#define THEcalloc  EC_calloc
extern void *EC_calloc(long long int nelem, long long int elsize);
#define THErealloc EC_realloc
extern void *EC_realloc(void *p, long long int size);
#define THEstrdup  EC_strdup
extern char *EC_strdup(const char *s);
#define THEfree    EC_free
extern void EC_free(void *p);

#else /* i.e. !defined(INTERCEPT_ALLOC) */

#define THEmalloc  malloc
#define THEcalloc  calloc
#define THErealloc realloc
#define THEstrdup  strdup
#define THEfree    free

#endif

#endif /* _INTERCEPT_ALLOC_H_ */
