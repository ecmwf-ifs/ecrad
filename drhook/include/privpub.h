
/* privpub.h */

#ifndef _PRIVPUB_H_
#define _PRIVPUB_H_

#include <stdio.h>
#include <stdbool.h>
#include <string.h>
#include <errno.h>
#include <signal.h>
#include <unistd.h>
#include <sys/types.h>

#if !defined(IO_BUFSIZE_DEFAULT)
/* The default I/O-buffer size = 4MBytes (unless overridden through ioassign) */
#define IO_BUFSIZE_DEFAULT 4194304
#endif

#define NIL "(nil)"

#define PERROR(msg) \
fprintf(stderr,"%s: %s (in %s:%d) : errno#%d\n",msg?msg:NIL,strerror(errno),__FILE__, __LINE__,errno)

#if defined(VPP5000) || defined(VPP)
/* Fujitsu VPP doesn't have snprintf, so we created one in file ifsaux/support/endian.c */
/* For ODB/SQL compiler the same code sits in odb/compiler/odb98.c */
extern int snprintf(char *str, size_t size, const char *format, ...);
#endif

#if defined(NO_TRUNC) || defined(VPP) || defined(SV2)
/* For systems without trunc() -function [an extension of ANSI-C, but usually available] */
#define trunc(x) ((x) - fmod((x),1))
#else
#if !defined(__cplusplus)
  extern double trunc(double d);
#endif
#endif

#define PRIVATE static
#define PUBLIC

/* fprintf even if "fp" was NULL */
extern int ODB_fprintf(FILE *fp, const char *format, ...);

/* Fortran90 compatible NINT-function */

#define F90nint(x) ( ((x) > 0) ? (int)((x) + 0.5) : (int)((x) - 0.5) )

/* Pi related stuff */

#define pi ((double)3.14159265358979323844e0)
#define half_pi           (pi / 2)
#define two_pi            (2 * pi)
#define four_pi           (4 * pi)
#define pi_over_180       (pi/180)
#define recip_pi_over_180 (180/pi)

/* Radians to degrees & degrees to radians (formerly in odb.h) */

#define R2D(x) ( (180/pi) * ( ((x) >  pi) ? ((x) - 2*pi) : (x) ) )
#define D2R(x) ( (pi/180) * ( ((x) > 180) ? ((x) -  360) : (x) ) )

/* Radius of the Earth both in meters and kilometres (formerly in odb.h) */

#define R_Earth    ((double)6371000e0)
#define R_Earth_km ((double)6371e0)

/* Some common things ... */

#undef MIN
#define MIN(a,b) ( ((a) < (b)) ? (a) :  (b) )

#undef MAX
#define MAX(a,b) ( ((a) > (b)) ? (a) :  (b) )

#undef ABS
#define ABS(x)   ( ((x) >= 0)  ? (x) : -(x) )

#undef RNDUP_DIV
#define RNDUP_DIV(i,n) (( (i) + (n) - 1 ) / (n))

#undef RNDUP
#define RNDUP(i,n) ( RNDUP_DIV(i,n) * (n))

#define STRLEN(x)  (((const char *)(x)) ? (int)strlen(x) : (int)0)

#define _POOLNO  "$pool#"
#define _NPOOLS  "$npools#"
#define _NTABLES "$ntables#"
#define _UNIQNUM "$uniq#"
#define _ROWNUM  "$row#"
#define _COLNUM  "$col#"
#define _NROWS   "$nrows#"
#define _NCOLS   "$ncols#"

#define IS_BSNUM(s)  ((s) && *(s) == '\\' && STRLEN(s) >= 2 && isdigit((s)[1]))
#define IS_HASH(s)   ((s) && *(s) == '#')
#define IS_DOLLAR(s) ((s) && *(s) == '$')
#define IS_TABLE(s) ((s) && (s)[0] == '@')

#define IS_POOLNO(s) ((s) && (((s)[0] == '$' && (s)[1] == '#')|| strequ(s,_POOLNO)))
#define IS_(sym,s) strequ(s,_##sym)
#define IS_USDHASH(s) ((s) && (s)[0] == '$' && (s)[STRLEN(s)-1] == '#')
#define IS_USDDOTHASH(s) (IS_USDHASH(s) && strchr(s,'.') && (strchr(s,'.') == strrchr(s,'.')))

static const char ODB_OFFSET_CHAR[] = "#";
#define GET_OFFSET(s) ((!(s) || IS_DOLLAR(s) || IS_HASH(s) || IS_USDHASH(s)) ? NULL : strrchr(s,*ODB_OFFSET_CHAR))

/* Timers */

extern double util_cputime_();
extern double util_walltime_();

#ifdef RS6K
extern void xl__trbk_();
#else
#define xl__trbk_()
#endif

extern void abor1fl_(const char *filename, const int *linenum, 
		     const char *s, 
		     int filenamelen, int slen);
extern void abor1_(const char *s, int slen);

#define ABOR1(txt) { \
  const char *t = (txt); \
  xl__trbk_(); t ? abor1_(t, strlen(t)) : abor1_("",0); }

#define ABOR1FL(txt) { \
  const char *t = (txt); \
  int linenum=__LINE__; \
  t ? abor1fl_(__FILE__, &linenum,  t, sizeof(__FILE__)-1, strlen(t)) \
    : abor1fl_(__FILE__, &linenum, "", sizeof(__FILE__)-1, 0); \
  _exit(1); /* Should never end up here */ }

#define RAISE(x) { \
  xl__trbk_(); \
  if ((x) == SIGABRT) { \
    const char *txt = "*** Fatal error; aborting (SIGABRT) ..."; \
    ABOR1FL(txt); \
    _exit(1); /* Should never end up here */ \
  } \
  else raise(x); \
}

/* No. of bits for "int" */
#define MAXBITS 32

/* No. of bits per 1-byte */
#define BITS_PER_BYTE 8

#ifndef O_LOCK_DONE
#define O_LOCK_DONE

/* OpenMP/ODB lock type */
/* Keep consistent with "ifsaux/include/drhook.h" */
/* Be ALSO consistent with OML_LOCK_KIND in ifsaux/module/oml_mod.F90 */

typedef long long int o_lock_t; /* i.e. 64-bit integer */

#define INIT_LOCKID_WITH_NAME(mylock, lockname) \
  coml_init_lockid_with_name_(mylock, lockname, strlen(lockname))

extern void coml_set_debug_(const int *konoff, int *kret);
extern void coml_init_lock_();
extern void coml_init_lockid_(o_lock_t *mylock);
extern void coml_init_lockid_with_name_(o_lock_t *mylock, const char *name, int name_len);
extern void coml_set_lock_();
extern void coml_set_lockid_(o_lock_t *mylock);
extern void coml_unset_lock_();
extern void coml_unset_lockid_(o_lock_t *mylock);
extern void coml_test_lock_(int *is_set);
extern void coml_test_lockid_(int *is_set, o_lock_t *mylock);
extern void coml_in_parallel_(int *is_parallel_region);

#endif /* #ifndef O_LOCK_DONE */

/* OpenMP-stuff */

/* My thread-id [1 .. $OMP_NUM_THREADS] as in ifsaux/module/oml_mod.F90 */
extern int
get_thread_id_(); /* From "ifsaux/utilities/get_thread_id.F90" */

/* Static number of threads allocated i.e. $OMP_NUM_THREADS */
extern int
get_max_threads_(); /* From "ifsaux/utilities/get_max_threads.F90" */

/* The current number of threads allocated <= $OMP_NUM_THREADS */
extern int
get_num_threads_(); /* From "ifsaux/utilities/get_num_threads.F90" */

#define DEF_IT    int it = get_thread_id_()
#define IT        (it-1)
#define DEF_INUMT int inumt = get_max_threads_()

/* 0 (zero) Celsius in Kelvin */

#define ZERO_POINT ((double)273.15e0)

/* Aggregate function flags */

#define ODB_AGGR_NONE            0x0000000

#define ODB_AGGRMASK_DISTINCT    0x0100000
#define ODB_AGGRMASK_BOOLEAN     0x1000000

#define ODB_AGGR_COUNT           0x0000001
#define ODB_AGGR_COUNT_DISTINCT  (ODB_AGGRMASK_DISTINCT | ODB_AGGR_COUNT)

#define ODB_AGGR_BCOUNT          (ODB_AGGRMASK_BOOLEAN | ODB_AGGR_COUNT)
#define ODB_AGGR_BCOUNT_DISTINCT (ODB_AGGRMASK_BOOLEAN | ODB_AGGR_COUNT_DISTINCT)

#define ODB_AGGR_MIN             0x0000002
#define ODB_AGGR_MAX             0x0000004

#define ODB_AGGR_SUM             0x0000008
#define ODB_AGGR_SUM_DISTINCT    (ODB_AGGRMASK_DISTINCT | ODB_AGGR_SUM)

#define ODB_AGGR_AVG             0x0000010
#define ODB_AGGR_AVG_DISTINCT    (ODB_AGGRMASK_DISTINCT | ODB_AGGR_AVG)

#define ODB_AGGR_STDEV           0x0000020
#define ODB_AGGR_STDEV_DISTINCT  (ODB_AGGRMASK_DISTINCT | ODB_AGGR_STDEV)

#define ODB_AGGR_RMS             0x0000040
#define ODB_AGGR_RMS_DISTINCT    (ODB_AGGRMASK_DISTINCT | ODB_AGGR_RMS)

#define ODB_AGGR_DOTP            0x0000080
#define ODB_AGGR_DOTP_DISTINCT   (ODB_AGGRMASK_DISTINCT | ODB_AGGR_DOTP)

#define ODB_AGGR_NORM            0x0000100
#define ODB_AGGR_NORM_DISTINCT   (ODB_AGGRMASK_DISTINCT | ODB_AGGR_NORM)

#define ODB_AGGR_VAR             0x0000200
#define ODB_AGGR_VAR_DISTINCT    (ODB_AGGRMASK_DISTINCT | ODB_AGGR_VAR)

#define ODB_AGGR_COVAR           0x0000400
#define ODB_AGGR_CORR            0x0000800

#define ODB_AGGR_LINREGR_A       0x0001000
#define ODB_AGGR_LINREGR_B       0x0002000

#define ODB_AGGR_MEDIAN          0x0010000
#define ODB_AGGR_MEDIAN_DISTINCT (ODB_AGGRMASK_DISTINCT | ODB_AGGR_MEDIAN)

#define ODB_AGGR_MINLOC          0x0020000
#define ODB_AGGR_MAXLOC          0x0040000

#define ODB_AGGR_DENSITY         0x0080000

#define ODB_tag_delim ";"


#ifdef __cplusplus
    typedef bool Bool;
#else
/*  typedef enum { false=0, true=1 } Bool; */
    typedef _Bool Bool; /* Std on C99 */
#endif

/* char : Normally occupies 1 byte */

typedef unsigned int  boolean;

typedef boolean Boolean;

typedef unsigned char  byte;
typedef byte           byte1;
#if !defined(RS6K)
typedef unsigned char uchar;
#endif
typedef char            integer1;
typedef unsigned char  uinteger1;

/* short : Normally occupies 2 bytes */

typedef short           integer2;
#if !defined(RS6K) && !defined(LINUX) && !defined(HPPA) && !defined(VPP) && !defined(SUN4) && !defined(NECSX) && !defined(SV2)
typedef unsigned short int ushort;
#endif
typedef unsigned short uinteger2;

/* int : Normally occupies 4 bytes */

/* typedef unsigned int uint; (defined via system include files) */
typedef int             integer4;
typedef unsigned int   uinteger4;

/* typedef unsigned long ulong; (defined via system include files) */

/* long long: A non-standard way to access 8 byte integer(s) */

typedef          long long int  longlong;
typedef          long long int  ll_t;
typedef unsigned long long int ulonglong;
typedef unsigned long long int u_ll_t;
typedef          long long int  integer8;
typedef unsigned long long int uinteger8;

/* float : 32-bit flp (normally IEEE-754), 4 bytes */

typedef float           real4;

/* double : 64-bit flp (normally IEEE-754), 8 bytes */

typedef double          real8;

/* Formula alias double : For use by SELECT-expressions (29-Jun-2006/SS) */

typedef double          Formula;

/* Supported data types */
/* Moved here from odb.h on 7-Jan-2004/SS */

typedef int number;
typedef int pk1int;
typedef int pk2int;
typedef int pk3int;
typedef int pk4int;
typedef int pk5int;
typedef int pk9int;

typedef pk1int linkoffset_t; /* also follows packing#1 */
typedef pk1int linklen_t;    /* also follows packing#1 */

/* 
   typedef unsigned int Bitfield;
   typedef unsigned int hex4;

   We can't have (at least the above) unsigned ints since the following 
   substitution for unsigned int u and double d gets incorrect on RS6K,
   when d < 0 [... the IBM-s.o.b. RS6K sets the result to zero !!]

   u = d ;

   what appears to be working is, but very weird is:
   u = (int)d;

   but we don't want to use this cast-operator, since its introduction 
   would add even more confusion

   ==> lets make Bitfield and hex4 signed ints (29-May-2002/SS)

*/

typedef int Bitfield;
typedef int hex4;

typedef hex4 bufr;
typedef hex4 grib;

typedef integer4 integer;

/* These inherit corresponding packing method */

typedef pk1int yyyymmdd;
typedef pk1int hhmmss;

typedef double pk2real;
typedef double pk3real;
typedef double pk4real;
typedef double pk5real;
typedef double pk9real;
typedef pk3real string;
/* typedef double string; */

typedef double pk11real;
typedef double pk12real;
typedef double pk13real;
typedef double pk14real;
typedef double pk15real;
typedef double pk16real;
typedef double pk17real;
typedef double pk18real;
typedef double pk19real;

typedef double pk21real;
typedef double pk22real;
typedef double pk23real;
typedef double pk24real;
typedef double pk25real;
typedef double pk26real;
typedef double pk27real;
typedef double pk28real;
typedef double pk29real;

typedef double pk31real;
typedef double pk32real;
typedef double pk33real;
typedef double pk34real;
typedef double pk35real;
typedef double pk36real;
typedef double pk37real;
typedef double pk38real;
typedef double pk39real;

typedef real8 real;

/* S2D-business */

#define S2D     "S2D_"
#define S2DLEN  (sizeof(S2D)-1) /* i.e. strlen(S2D) */

#define S2Dlc     "s2d_"
#define S2DlcLEN  (sizeof(S2Dlc)-1) /* i.e. strlen(S2Dlc) */

#define S2D_all_blanks (0x2020202020202020ull)

typedef union {
  double dval;
  u_ll_t llu;
  uint ival[2];
  char str[9];
  char *saddr;
} S2D_Union;

typedef struct {
  uint ival[2];
} S2D_Type;

#define MASK_0           0U  /*                                 0 */
#define MASK_1           1U  /*                                 1 */
#define MASK_2           3U  /*                                11 */
#define MASK_3           7U  /*                               111 */
#define MASK_4          15U  /*                              1111 */
#define MASK_5          31U  /*                             11111 */
#define MASK_6          63U  /*                            111111 */
#define MASK_7         127U  /*                           1111111 */
#define MASK_8         255U  /*                          11111111 */
#define MASK_9         511U  /*                         111111111 */
#define MASK_10       1023U  /*                        1111111111 */
#define MASK_11       2047U  /*                       11111111111 */
#define MASK_12       4095U  /*                      111111111111 */
#define MASK_13       8191U  /*                     1111111111111 */
#define MASK_14      16383U  /*                    11111111111111 */
#define MASK_15      32767U  /*                   111111111111111 */
#define MASK_16      65535U  /*                  1111111111111111 */
#define MASK_17     131071U  /*                 11111111111111111 */
#define MASK_18     262143U  /*                111111111111111111 */
#define MASK_19     524287U  /*               1111111111111111111 */
#define MASK_20    1048575U  /*              11111111111111111111 */
#define MASK_21    2097151U  /*             111111111111111111111 */
#define MASK_22    4194303U  /*            1111111111111111111111 */
#define MASK_23    8388607U  /*           11111111111111111111111 */
#define MASK_24   16777215U  /*          111111111111111111111111 */
#define MASK_25   33554431U  /*         1111111111111111111111111 */
#define MASK_26   67108863U  /*        11111111111111111111111111 */
#define MASK_27  134217727U  /*       111111111111111111111111111 */
#define MASK_28  268435455U  /*      1111111111111111111111111111 */
#define MASK_29  536870911U  /*     11111111111111111111111111111 */
#define MASK_30 1073741823U  /*    111111111111111111111111111111 */
#define MASK_31 2147483647U  /*   1111111111111111111111111111111 */
#define MASK_32 4294967295U  /*  11111111111111111111111111111111 */

#define MASK(n) MASK_##n

#define IOR(x,y)   ((x) | (y))
#define IAND(x,y)  ((x) & (y))

#define ISHFTL(x,n) ((x) << (n))
#define ISHFTR(x,n) ((x) >> (n))

#define GET_BITS(x, pos, len)      IAND(ISHFTR((int)(x), pos), MASK(len))

#define PUT_BITS(to, x, pos, len) \
  to = IOR(IAND(to, ~(ISHFTL(MASK(len), pos))), \
           IAND((int)(x), ISHFTL(MASK(len), pos)))

#define BYTESIZE(x) (x##len * sizeof(*x))
#define BYTESIZE2(x,elemsize) (x##len * (elemsize))

typedef struct {
#ifdef LITTLE
  unsigned int pmethod:8; /* Value 0 .. 255 ; 
			     0 means not packed;
			     1,2,3,5,9,11-19,21-29 methods implemented
			     255 means not really packed, but PCMA-header exists for MDI's */

  unsigned int signbit:1;  /* signed = 1, unsigned = 0 */

  unsigned int byte_swappable:1; /* byte swapping applicable = 1, N/A i.e. raw data = 0 */

  unsigned int precision_bits:6; /* power of 2 for length of data; range = [0..63] bits */

  unsigned int base_type:2; /* other = 0, int = 1, flp = 2, complex flp = 3 
			       notes: 
			       - if other, look for other_type
			       - precion_bits used for length as 2^precision_bits
			       - if int, check also signed/unsigned part
			    */
  
  unsigned int other_type:14; /* used when base_type == 0 ; currently defined
				 undefined   = 00 000 000 000 000 =     0 = 0x0
				 Bitfield    = 00 000 000 000 001 =     1 = 0x1
				 (note: also precision_bits must be 5 i.e. 2^5 == 32; 
				 preferably unsigned)
				 string      = 00 000 000 000 010 =     2 = 0x2 
				 (note: alias to char(8); 
				 also precision_bits should be 3 i.e. 2^3 == 8 bytes;
				 byte_swappable must be == 0 i.e. raw)
				 yyyymmdd    = 00 000 000 000 100 =     4 = 0x4
				 (note: preferably unsigned)
				 hhmmss      = 00 000 000 001 000 =     8 = 0x8
				 (note: preferably unsigned)
				 linkoffset_t = 00 000 000 010 000 =    16 = 0x10
				 (note: preferably unsigned)
				 linklen_t    = 00 000 000 100 000 =    32 = 0x20
				 (note: preferably unsigned)
				 bufr-data   = 00 000 001 000 000 =    64 = 0x40
				 (note: also byte_swappable must be == 0 i.e. raw;
				 also precision_bits must be 5 i.e. 2^5 == 32)
				 grib-data   = 00 000 010 000 000 =   128 = 0x80
				 (note: also byte_swappable must be == 0 i.e. raw;
				 also precision_bits must be 5 i.e. 2^5 == 32)
				 blob        = 00 000 100 000 000 =   256 = 0x100
				 (note: precision_bits used for length; up to 2^16 bytes; enough ?;
				 byte_swappable must be == 0 i.e. raw)
				 long blob   = 00 001 000 000 000 =  512 = 0x200
				 (note: precision_bits used for length; up to 2^31 bytes; enough ?;
				 byte_swappable must be == 0 i.e. raw)
				 char        = 00 010 000 000 000 =  1024 = 0x400
				 (note: precision_bits used for length;
				 byte_swappable must be == 0 i.e. raw)
				 varchar     = 00 100 000 000 000 =  2048 = 0x800
				 (note: precision_bits used for maximum length;
				 byte_swappable must be == 0 i.e. raw)
				 reserved_1  = 01 000 000 000 000 =  4096 = 0x1000
				 reserved_2  = 10 000 000 000 000 =  8192 = 0x2000
			      */ 
#else /* Big endian */
  unsigned int other_type:14;
  unsigned int base_type:2;
  unsigned int precision_bits:6;
  unsigned int byte_swappable:1;
  unsigned int signbit:1;
  unsigned int pmethod:8;
#endif
} odb_types_t;

PUBLIC uint get_dtnum(const char *dt);
PUBLIC  int get_dtsize(const char *dt);

#define EXTRACT_PMETHOD(x)   GET_BITS(x,0,8)
#define EXTRACT_DATATYPE(x)  GET_BITS(x,8,24)
#define EXTRACT_SIGNBIT(x)   GET_BITS(x,8,1)
#define EXTRACT_SWAPPABLE(x) GET_BITS(x,9,1)
#define EXTRACT_PRECISION(x) GET_BITS(x,10,6)
#define EXTRACT_BASETYPE(x)  GET_BITS(x,16,2)
#define EXTRACT_OTHERTYPE(x) GET_BITS(x,18,14)

/* The following is a product of ../tools/typebits.c -program output */

/* signbit=0, byte_swappable=0, precision_bits=0, base_type=0, other_type=0x0 (0) */
#define DATATYPE_UNDEF                0x0          /* (0 dec) undef */

/* signbit=0, byte_swappable=0, precision_bits=0, base_type=1, other_type=0x0 (0) */
#define DATATYPE_BIT                  0x100        /* (256 dec) bit : 1 bits, 0 bytes */

/* signbit=1, byte_swappable=1, precision_bits=3, base_type=1, other_type=0x0 (0) */
#define DATATYPE_INT1                 0x10f        /* (271 dec) char : 8 bits, 1 bytes */

/* signbit=1, byte_swappable=1, precision_bits=4, base_type=1, other_type=0x0 (0) */
#define DATATYPE_INT2                 0x113        /* (275 dec) short : 16 bits, 2 bytes */

/* signbit=1, byte_swappable=1, precision_bits=5, base_type=1, other_type=0x0 (0) */
#define DATATYPE_INT4                 0x117        /* (279 dec) int : 32 bits, 4 bytes */

/* signbit=1, byte_swappable=1, precision_bits=6, base_type=1, other_type=0x0 (0) */
#define DATATYPE_INT8                 0x11b        /* (283 dec) long long : 64 bits, 8 bytes */

/* signbit=0, byte_swappable=1, precision_bits=3, base_type=1, other_type=0x0 (0) */
#define DATATYPE_UINT1                0x10e        /* (270 dec) uchar : 8 bits, 1 bytes */

/* signbit=0, byte_swappable=1, precision_bits=4, base_type=1, other_type=0x0 (0) */
#define DATATYPE_UINT2                0x112        /* (274 dec) ushort : 16 bits, 2 bytes */

/* signbit=0, byte_swappable=1, precision_bits=5, base_type=1, other_type=0x0 (0) */
#define DATATYPE_UINT4                0x116        /* (278 dec) uint : 32 bits, 4 bytes */

/* signbit=0, byte_swappable=1, precision_bits=6, base_type=1, other_type=0x0 (0) */
#define DATATYPE_UINT8                0x11a        /* (282 dec) ulonglong : 64 bits, 8 bytes */

/* signbit=0, byte_swappable=1, precision_bits=5, base_type=2, other_type=0x0 (0) */
#define DATATYPE_REAL4                0x216        /* (534 dec) float : 32 bits, 4 bytes */

/* signbit=0, byte_swappable=1, precision_bits=6, base_type=2, other_type=0x0 (0) */
#define DATATYPE_REAL8                0x21a        /* (538 dec) double : 64 bits, 8 bytes */

/* signbit=0, byte_swappable=1, precision_bits=7, base_type=2, other_type=0x0 (0) */
#define DATATYPE_REAL16               0x21e        /* (542 dec) long double : 128 bits, 16 bytes */

/* signbit=0, byte_swappable=1, precision_bits=6, base_type=3, other_type=0x0 (0) */
#define DATATYPE_CMPLX4               0x31a        /* (794 dec) complex4 : 64 bits, 8 bytes */

/* signbit=0, byte_swappable=1, precision_bits=7, base_type=3, other_type=0x0 (0) */
#define DATATYPE_CMPLX8               0x31e        /* (798 dec) complex8 : 128 bits, 16 bytes */

/* signbit=0, byte_swappable=1, precision_bits=8, base_type=3, other_type=0x0 (0) */
#define DATATYPE_CMPLX16              0x322        /* (802 dec) complex16 : 256 bits, 32 bytes */

/* signbit=0, byte_swappable=1, precision_bits=5, base_type=0, other_type=0x1 (1) */
#define DATATYPE_BITFIELD             0x416        /* (1046 dec) Bitfield : 32 bits, 4 bytes */

/* signbit=0, byte_swappable=0, precision_bits=6, base_type=0, other_type=0x2 (2) */
#define DATATYPE_STRING               0x818        /* (2072 dec) string : 64 bits, 8 bytes */

/* signbit=0, byte_swappable=1, precision_bits=5, base_type=0, other_type=0x4 (4) */
#define DATATYPE_YYYYMMDD             0x1016       /* (4118 dec) yyyymmdd : 32 bits, 4 bytes */

/* signbit=0, byte_swappable=1, precision_bits=5, base_type=0, other_type=0x8 (8) */
#define DATATYPE_HHMMSS               0x2016       /* (8214 dec) hhmmss : 32 bits, 4 bytes */

/* signbit=1, byte_swappable=1, precision_bits=5, base_type=0, other_type=0x10 (16) */
#define DATATYPE_LINKOFFSET           0x4017       /* (16407 dec) linkoffset_t : 32 bits, 4 bytes */

/* signbit=1, byte_swappable=1, precision_bits=5, base_type=0, other_type=0x20 (32) */
#define DATATYPE_LINKLEN              0x8017       /* (32791 dec) linklen_t : 32 bits, 4 bytes */

/* signbit=0, byte_swappable=1, precision_bits=5, base_type=0, other_type=0x40 (64) */
#define DATATYPE_BUFR                 0x10016      /* (65558 dec) bufr : 32 bits, 4 bytes */

/* signbit=0, byte_swappable=1, precision_bits=5, base_type=0, other_type=0x80 (128) */
#define DATATYPE_GRIB                 0x20016      /* (131094 dec) grib : 32 bits, 4 bytes */

/* signbit=0, byte_swappable=0, precision_bits=19, base_type=0, other_type=0x100 (256) */
#define DATATYPE_BLOB                 0x4004c      /* (262220 dec) blob64kB : 65536 bytes */

/* signbit=0, byte_swappable=0, precision_bits=34, base_type=0, other_type=0x200 (512) */
#define DATATYPE_LONGBLOB             0x80088      /* (524424 dec) blob2GB : 2147483648 bytes */

/* signbit=0, byte_swappable=0, precision_bits=11, base_type=0, other_type=0x400 (1024) */
#define DATATYPE_CHAR                 0x10002c     /* (1048620 dec) char(1:255) : 256 bytes */

/* signbit=0, byte_swappable=0, precision_bits=11, base_type=0, other_type=0x800 (2048) */
#define DATATYPE_VARCHAR              0x20002c     /* (2097196 dec) varchar(1:255) : 256 bytes */

#endif
