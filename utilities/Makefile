SOURCES = easy_netcdf.F90 radiation_io.F90

OBJECTS := $(SOURCES:.F90=.o)
LIBUTILITIES = ../lib/libutilities.a

all: $(LIBUTILITIES)

$(LIBUTILITIES): $(OBJECTS)
	ar r $(LIBUTILITIES) $(OBJECTS)

%.o: %.F90
	$(FC) $(FCFLAGS) -c $<

clean:
	rm -f *.o $(LIBUTILITIES)

easy_netcdf.o: radiation_io.o
