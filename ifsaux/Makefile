SOURCES = parkind1.F90 yomlun_ifsaux.F90 yomcst.F90 abor1.F90 \
	yomtag.F90 mpl_module.F90 yommp0.F90 yomdyncore.F90 yomlun.F90

OBJECTS := $(SOURCES:.F90=.o)
LIBIFSAUX = ../lib/libifsaux.a

all: $(LIBIFSAUX)

deps: includes

$(LIBIFSAUX): $(OBJECTS)
	ar r $(LIBIFSAUX) $(OBJECTS)

%.o: %.F90
	$(FC) $(FCFLAGS) -c $<

includes:
	cd .. && fcm make --config-file=bin/fcm-make-interfaces.cfg interfaces.ns-incl=ifsaux interfaces.source=$(shell dirname $(realpath $(firstword $(MAKEFILE_LIST))))/..

clean:
	rm -f *.o $(LIBIFSAUX)

yomlun_ifsaux.o yomhook.o random_numbers_mix.o yomcst.o mpl_module.o yomtag.o yommp0.o yomdyncore.o: parkind1.o
abor1.o: yomlun_ifsaux.o parkind1.o
yomlun.o: yomlun_ifsaux.o

.PHONY: deps includes
