# ecRad Makefile - read the README file before editing

#############################
### --- CONFIGURATION --- ###
#############################

# Use the nf-config utility, if available, to set the NETCDF_INCLUDE
# and NETCDF_LIB flags
HAVE_NFCONFIG := $(shell nf-config --version 2> /dev/null)
ifdef HAVE_NFCONFIG
$(info *** Using nf-config to obtain NetCDF flags)
NETCDF_INCLUDE = $(shell nf-config --fflags)
NETCDF_LIB     = $(shell nf-config --flibs)
ifeq ($(shell nf-config --has-nc4),yes)
NETCDF4        = 1
endif
else
$(info *** nf-config not found)
endif

# make can be invoked using "make PROFILE=<prof>" in which case your
# local configuration parameters will be obtained from
# Makefile_include.<prof>
ifndef PROFILE
$(info *** No "PROFILE" variable provided, assuming "gfortran")
PROFILE = gfortran
endif

# Include a platform-specific makefile that defines FC, FCFLAGS and
# LIBS
include	Makefile_include.$(PROFILE)

# Check for presence of the NETCDF_INCLUDE and NETCDF_LIB flags
ifndef NETCDF_INCLUDE
$(info *** You may need to set NETCDF_INCLUDE manually)
endif
ifndef NETCDF_LIB
$(info *** You may need to set NETCDF_LIB manually)
endif

# Add single-precision flag if SINGLE_PRECISION=1 was given on the
# "make" command line
ifdef SINGLE_PRECISION
CPPFLAGS += -DPARKIND1_SINGLE
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
# Allow the capability to write NetCDF4/HDF5 files, provided the code
# is compiled against the NetCDF4 library
ifdef NETCDF4
$(info *** Building with NetCDF4/HDF5 support)
CPPFLAGS += -DNC_NETCDF4
endif

# Consolidate flags
export FC
export FCFLAGS = $(WARNFLAGS) $(BASICFLAGS) $(CPPFLAGS) -I../include \
	$(OPTFLAGS) $(DEBUGFLAGS) $(NETCDF_INCLUDE) $(OMPFLAG)
export LIBS    = $(LDFLAGS) -L../lib -lradiation -lutilities \
	-lifsrrtm -lifsaux $(FCLIBS) $(NETCDF_LIB) $(OMPFLAG)

# Do we include Dr Hook from ECMWF's fiat library?
ifdef FIATDIR
# Prepend location of yomhook.mod module file from fiat library, so
# that it is found in preference to the dummy one in ecRad
FCFLAGS := -I$(FIATDIR)/module/fiat $(FCFLAGS)
# Append fiat library (usually shared: libfiat.so)
LIBS += -L$(FIATDIR)/lib -Wl,-rpath,$(FIATDIR)/lib -lfiat
else
# Dummy Dr Hook library
LIBS += -ldrhook
endif


#############################
### --- BUILD TARGETS --- ###
#############################

all: build

help:
	@echo "Usage:"
	@echo "  make PROFILE=<prof>"
	@echo "where <prof> is one of gfortran, pgi, intel or cray (see Makefile_include.<prof>)"
	@echo "Other possible arguments are:"
	@echo "  DEBUG=1              Compile with debug settings on and optimizations off"
	@echo "  SINGLE_PRECISION=1   Compile with single precision"
	@echo "  FIATDIR=/my/path     Compile with Dr Hook, specifying the directory containing lib/libfiat.so and module/fiat/yomhook.mod"
	@echo "  test                 Run test cases in test directory"
	@echo "  clean                Remove all compiled files"

ifndef FIATDIR
build: directories libifsaux libdummydrhook libutilities libifsrrtm \
	libradiation driver ifsdriver symlinks
libradiation libutilities: libdummydrhook
else
# Note that if we are using Dr Hook from the fiat library we don't
# want to create mod/yomhook.mod as this can sometimes be found before
# the one in the fiat directory leading to an error at link stage
build: directories libifsaux libutilities libifsrrtm libradiation \
	driver ifsdriver symlinks
endif

# git cannot store empty directories so they may need to be created 
directories: mod lib
mod:
	mkdir -p mod
lib:
	mkdir -p lib

deps: clean-deps
	cd ifsaux && $(MAKE) deps
	cd ifsrrtm && $(MAKE) deps
	cd ifs && $(MAKE) deps

clean-deps:
	rm -f include/*.intfb.h

libifs: libradiation
	cd ifs && $(MAKE)

libifsaux:
	cd ifsaux && $(MAKE)

libdummydrhook: libifsaux
	cd drhook && $(MAKE) dummy

libutilities: libifsaux
	cd utilities && $(MAKE)

libifsrrtm: libifsaux
	cd ifsrrtm && $(MAKE)

libradiation: libifsrrtm libutilities libifsaux
	cd radiation && $(MAKE)

driver: libifsaux libifsrrtm libutilities libradiation
	cd driver && $(MAKE) driver

ifsdriver: libifsaux libifsrrtm libutilities libradiation libifs
	cd driver && $(MAKE) ifs_driver

test_programs: driver
	cd driver && $(MAKE) test_programs

symlinks: clean-symlinks
	cd practical && ln -s ../bin/ecrad
	cd practical && ln -s ../data

test: test_ifs test_i3rc test_ckdmip

test_ifs: driver
	cd test/ifs && $(MAKE) test

test_i3rc: driver
	cd test/i3rc && $(MAKE) test

test_ckdmip:
	cd test/ckdmip && $(MAKE) test

clean: clean-tests clean-toplevel clean-utilities clean-mods clean-symlinks

clean-tests:
	cd test/ifs && $(MAKE) clean
	cd test/i3rc && $(MAKE) clean
	cd test/ckdmip && $(MAKE) clean

clean-toplevel:
	cd radiation && $(MAKE) clean
	cd driver && $(MAKE) clean

clean-utilities:
	cd ifsaux && $(MAKE) clean
	cd utilities && $(MAKE) clean
	cd ifsrrtm && $(MAKE) clean
	cd drhook && $(MAKE) clean
	cd ifs && $(MAKE) clean

clean-mods:
	rm -f mod/*.mod

clean-symlinks:
	rm -f practical/ecrad practical/data

clean-autosaves:
	rm -f *~ .gitignore~ */*~ */*/*~

.PHONY: all build help deps clean-deps libifsaux libdummydrhook libutilities libifsrrtm \
	libradiation driver symlinks clean clean-toplevel test test_ifs ifsdriver \
	test_i3rc clean-tests clean-utilities clean-mods clean-symlinks
