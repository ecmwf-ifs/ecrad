SOURCES = tcrad_3region_solver.F90 tcrad_two_stream.F90

OBJECTS := $(SOURCES:.F90=.o)
LIBMR = ../lib/libtcrad.a

all: $(LIBMR)

$(LIBMR): $(OBJECTS)
	ar r $(LIBMR) $(OBJECTS)

%.o: %.F90
	$(FC) $(FCFLAGS) -c $<

clean:
	rm -f *.o $(LIBMR)

tcrad_3region_solver.o: tcrad_two_stream.o

tcrad_3region_solver.o: tcrad_matrix.h tcrad_flux.h tcrad_radiance.h
