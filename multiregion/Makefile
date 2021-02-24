SOURCES = tripleclouds_solver.F90

OBJECTS := $(SOURCES:.F90=.o)
LIBMR = ../lib/libmultiregion.a

all: $(LIBMR)

$(LIBMR): $(OBJECTS)
	ar r $(LIBMR) $(OBJECTS)

%.o: %.F90
	$(FC) $(FCFLAGS) -c $<

clean:
	rm -f *.o $(LIBMR)

tripleclouds_solver.o: matrix_functions.h calc_multiregion_flux.h
