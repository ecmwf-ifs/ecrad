ifndef PROFILE
$(info No "PROFILE" variable provided, assuming "gfortran"...)
PROFILE = gfortran
endif

# Include a platform-specific makefile that defines FC, FCFLAGS and
# LIBS
include	Makefile_include.$(PROFILE)
all: build

help:
	@echo "Usage:"
	@echo "  make PROFILE=<prof>"
	@echo "where <prof> is one of gfortran, pgi etc."


# Add single-precision flag if SINGLE_PRECISION=1 was given on the
# "make" command line
ifdef SINGLE_PRECISION
CPPFLAGS += -DSINGLE_PRECISION
endif

# If PRINT_ENTRAPMENT_DATA=1 was given on the "make" command line
# then the SPARTACUS shortwave solver will write data to fort.101 and
# fort.102
ifdef PRINT_ENTRAPMENT_DATA
CPPFLAGS += -DPRINT_ENTRAPMENT_DATA 
endif
# For backwards compatibility we allow the following as well
ifdef PRINT_ENCROACHMENT_DATA
CPPFLAGS += -DPRINT_ENTRAPMENT_DATA 
endif

# Consolidate flags
export FC
export FCFLAGS = $(WARNFLAGS) $(BASICFLAGS) $(CPPFLAGS) -I../include \
	$(OPTFLAGS) $(DEBUGFLAGS) $(NETCDF_INCLUDE) $(OMPFLAG)
export LIBS    = $(LDFLAGS) -L../lib -lradsurf -lradiation -lutilities \
	-lifsrrtm -lifsaux $(FCLIBS) $(NETCDF_LIB) $(OMPFLAG)

ifdef USE_PSRAD
build: libifsaux libutilities libpsradrrtm libifsrrtm libradiation libradsurf driver
else
build: libifsaux libutilities libifsrrtm libradiation libradsurf driver
endif

deps: clean-deps
	cd ifsaux && $(MAKE) deps
	cd ifsrrtm && $(MAKE) deps

clean-deps:
	rm -f include/*.intfb.h

libifsaux:
	cd ifsaux && $(MAKE)

libutilities:
	cd utilities && $(MAKE)

libpsradrrtm:
	cd psradrrtm && $(MAKE)

libifsrrtm:
	cd ifsrrtm && $(MAKE)

libradiation:
	cd radiation && $(MAKE)

libradsurf:
	cd radsurf && $(MAKE)

driver:
	cd driver && $(MAKE)

test: test_ifs test_i3rc

test_ifs:
	cd test/ifs && $(MAKE) test

test_i3rc:
	cd test/i3rc && $(MAKE) test

test_surface:
	cd test/surface && $(MAKE) test

clean: clean-tests clean-toplevel clean-utilities clean-mods

clean-tests:
	cd test/ifs && $(MAKE) clean
	cd test/i3rc && $(MAKE) clean
	cd test/surface && $(MAKE) clean

clean-toplevel:
	cd radiation && $(MAKE) clean
	cd radsurf && $(MAKE) clean
	cd driver && $(MAKE) clean

clean-utilities:
	cd ifsaux && $(MAKE) clean
	cd utilities && $(MAKE) clean
	cd ifsrrtm && $(MAKE) clean

clean-mods:
	rm -f mod/*.mod

clean-autosaves:
	rm -f *~ */*~ */*/*~

.PHONY: libifsaux libpsradrrtm libifsrrtm libradiation libradsurf driver clean clean-toplevel test
