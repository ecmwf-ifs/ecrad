/* cargs.h */

/* prototypes for ifsaux/support/cargs.c */

/* Author: Sami Saarinen, ECMWF, 27-Apr-2006 */

#if defined(__cplusplus)
extern "C" {
#endif

void ec_PutArgs(int argc, char *argv[]);
const char *ec_GetArgs(int argno);
int ec_NumArgs(void);

/* The following two as in C-main : "int main(int argc, char *argv[])" */

int ec_argc(void);
char **ec_argv(void);

/* Fortran interface */

int iargc_c_(void);
int iargc_c (void);

void getarg_c_(const int *argno, char *arg
	       /* Hidden argument */
	       , const int arg_len);

void getarg_c (const int *argno, char *arg
	       /* Hidden argument */
	       , const int arg_len);

void putarg_c_(const int *argno, const char *arg
	       /* Hidden argument */
	       , int arg_len);

void putarg_c (const int *argno, const char *arg
	       /* Hidden argument */
	       , int arg_len);

void putarg_info_(const int *argc, const char *cterm
		  /* Hidden argument */
		  , int cterm_len);

void putarg_info (const int *argc, const char *cterm
		  /* Hidden argument */
		  , int cterm_len);

/* From ifsaux/support/cmpl_binding.F90 */

void cmpl_getarg_(const int *argno, char *arg
		  /* Hidden argument */
		  , const int arg_len);

int cmpl_iargc_(); 

#if defined(__cplusplus)
}
#endif
