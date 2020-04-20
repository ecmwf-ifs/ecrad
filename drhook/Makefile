MAKE_INCLUDES = ../bin/make_intfbl.1.pl
INCLUDE_DIR = include

OBJSMPISER := $(patsubst %.F,%.o,$(wildcard mpi_serial/*.F))

OBJSHOOK    = $(patsubst %.F90,%.o,$(wildcard drhook/*.F90))
OBJSHOOK   += $(patsubst %.F,%.o,  $(wildcard drhook/*.F))
OBJSHOOK   += $(patsubst %.c,%.o,  $(wildcard drhook/*.c))
OBJSHOOK   += $(patsubst %.cc,%.o, $(wildcard drhook/*.cc))

LIBDRHOOK = ../lib/libdrhook.a

all: $(LIBDRHOOK)

deps: includes

$(LIBDRHOOK): yomhook.o $(OBJSHOOK) $(OBJSMPISER)
	ar -r $(LIBDRHOOK) yomhook.o $(OBJSHOOK) $(OBJSMPISER)

%.o: %.F90
	$(FC) $(FCFLAGS) -I./mpi_serial -I$(INCLUDE_DIR) -c $< -o $@

%.o: %.F
	$(FC) $(FCFLAGS) -I./mpi_serial -I$(INCLUDE_DIR) -c $< -o $@

%.o: %.c
	$(CC) $(CFLAGS)  -I$(INCLUDE_DIR) -c $< -o $@

%.o: %.cc 
	$(CXX) $(CFLAGS) -I$(INCLUDE_DIR) -c $< -o $@


include drhook/.depend

includes:
	LOC_INTFBDIR=$(INCLUDE_DIR) INTFBDIR=$(INCLUDE_DIR) $(MAKE_INCLUDES) $(SOURCES)

# Dummy version of Dr Hook creates the YOMHOOK module but the DR_HOOK
# routine does nothing
dummy: yomhook_dummy.o
	ar -r $(LIBDRHOOK) yomhook_dummy.o

clean:
	rm -f *.o drhook/*.o mpi_serial/*.o $(LIBDRHOOK)

.PHONY: deps includes
