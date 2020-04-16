SOURCES  = parkind1.F90 yomlun_ifsaux.F90 yomcst.F90 abor1.F90 yomhook.F90 \
	yomtag.F90 mpl_module.F90 yommp0.F90 yomdyncore.F90 yomlun.F90 \
	random_numbers_mix.F90



MAKE_INCLUDES = ../bin/make_intfbl.1.pl
INCLUDE_DIR = ../include

OBJSAUX := $(SOURCES:.F90=.o)

OBJSMPISER := $(patsubst %.F,%.o,$(wildcard mpi_serial/source/*.F))

OBJSHOOK    = $(patsubst %.F90,%.o,$(wildcard drhook/*.F90))
OBJSHOOK   += $(patsubst %.F,%.o,  $(wildcard drhook/*.F))
OBJSHOOK   += $(patsubst %.c,%.o,  $(wildcard drhook/*.c))
OBJSHOOK   += $(patsubst %.cc,%.o, $(wildcard drhook/*.cc))


LIBIFSAUX = ../lib/libifsaux.a

all: $(LIBIFSAUX)

deps: includes

$(LIBIFSAUX): $(OBJSAUX) $(OBJSHOOK) $(OBJSMPISER)
	ar -r $(LIBIFSAUX) $(OBJSAUX) $(OBJSHOOK) $(OBJSMPISER)

%.o: %.F90
	$(FC) $(FCFLAGS) -I./mpi_serial/source -I$(INCLUDE_DIR)/drhook -c $< -o $@

%.o: %.F
	$(FC) $(FCFLAGS) -I./mpi_serial/source -I$(INCLUDE_DIR)/drhook -c $< -o $@

%.o: %.c
	$(CC) $(CFLAGS)  -I$(INCLUDE_DIR)/drhook -c $< -o $@

%.o: %.cc 
	$(CXX) $(CFLAGS) -I$(INCLUDE_DIR)/drhook -c $< -o $@


include drhook/.depend

includes:
	LOC_INTFBDIR=$(INCLUDE_DIR) INTFBDIR=$(INCLUDE_DIR) $(MAKE_INCLUDES) $(SOURCES)

clean:
	rm -f *.o drhook/*.o mpi_serial/source/*.o $(LIBIFSAUX)

yomlun_ifsaux.o, yomhook.o, random_numbers_mix.o: parkind1.o
yomlun.o : yomlun_ifsaux.o

.PHONY: deps includes
