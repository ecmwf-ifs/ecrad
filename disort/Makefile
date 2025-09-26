SOURCES = BDREF.f DISOBRDF.f DISORT.f DISOTESTAUX.f ERRPACK.f LAPACK.f LINPAK.f RDI1MACH.f

OBJECTS := $(SOURCES:.f=.o)
OBJECTS += cgqf.o
LIBDISORT = ../lib/libdisort.a

ifdef SINGLE_PRECISION
  # DISORT takes single-precision arguments by default	
  LOCALFLAGS = 
else
  # gfortran-specific flags to use double precision
  LOCALFLAGS = $(FORCEDOUBLEFLAGS)
endif

all: $(LIBDISORT)

$(LIBDISORT): $(OBJECTS)
	ar r $(LIBDISORT) $(OBJECTS)

%.o: %.f
	$(FC) $(FCFLAGS) $(LOCALFLAGS) -c $<

cgqf.o: cgqf.F90
	$(FC) $(FCFLAGS) $(LOCALFLAGS) -c $<

clean:
	rm -f *.o $(LIBDISORT)

