#ifdef NAG
!-- The NAG-compiler has its own proper F90_UNIX_* modules
!   We need this for non-NAG compilers to satisfy some
!   make/build-environments only
#define F90_UNIX_IO NAGDUMMY_F90_UNIX_IO
#endif
MODULE F90_UNIX_IO
LOGICAL L_NAGDUMMY
END MODULE F90_UNIX_IO
#ifdef NAG
#undef F90_UNIX_IO
#endif
