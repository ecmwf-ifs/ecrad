SOURCES = parkind1.F90 yomlun_ifsaux.F90 yomcst.F90 abor1.F90 \
	yomtag.F90 mpl_module.F90 yommp0.F90 yomdyncore.F90 yomlun.F90
MAKE_INCLUDES = ../bin/make_intfbl.1.pl
INCLUDE_DIR = ../include

OBJECTS := $(SOURCES:.F90=.o)
LIBIFSAUX = ../lib/libifsaux.a

all: $(LIBIFSAUX)

deps: includes

$(LIBIFSAUX): $(OBJECTS)
	ar r $(LIBIFSAUX) $(OBJECTS)

%.o: %.F90
	$(FC) $(FCFLAGS) -c $<

includes:
	LOC_INTFBDIR=$(INCLUDE_DIR) INTFBDIR=$(INCLUDE_DIR) $(MAKE_INCLUDES) $(SOURCES)

clean:
	rm -f *.o $(LIBIFSAUX)

include Makefile_deps

.PHONY: deps includes
