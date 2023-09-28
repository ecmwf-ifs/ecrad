OBJECTS = parrrtm.o parsrtm.o rrtm_prepare_gases.o	\
rrtm_gas_optical_depth.o rrtm_setcoef_140gp.o	\
rrtm_taumol1.o rrtm_taumol10.o rrtm_taumol11.o rrtm_taumol12.o		\
rrtm_taumol13.o rrtm_taumol14.o rrtm_taumol15.o rrtm_taumol16.o		\
rrtm_taumol2.o rrtm_taumol3.o rrtm_taumol4.o rrtm_taumol5.o		\
rrtm_taumol6.o rrtm_taumol7.o rrtm_taumol8.o rrtm_taumol9.o		\
srtm_setcoef.o surrtrf.o srtm_gas_optical_depth.o 			\
srtm_taumol16.o srtm_taumol17.o srtm_taumol18.o srtm_taumol19.o		\
srtm_taumol20.o srtm_taumol21.o srtm_taumol22.o srtm_taumol23.o		\
srtm_taumol24.o srtm_taumol25.o srtm_taumol26.o srtm_taumol27.o		\
srtm_taumol28.o srtm_taumol29.o yoerrta1.o yoerrta10.o yoerrta11.o	\
yoerrta12.o yoerrta13.o yoerrta14.o yoerrta15.o yoerrta16.o		\
yoerrta2.o yoerrta3.o yoerrta4.o yoerrta5.o yoerrta6.o yoerrta7.o	\
yoerrta8.o yoerrta9.o yoerrtab.o yoesrta16.o yoesrta17.o yoesrta18.o	\
yoesrta19.o yoesrta20.o yoesrta21.o yoesrta22.o yoesrta23.o		\
yoesrta24.o yoesrta25.o yoesrta26.o yoesrta27.o yoesrta28.o		\
yoesrta29.o yoesrtab.o yoesrtaer.o yoesw.o 				\
yom_ygfl.o yomdimv.o yoerrtm.o yoerrtftr.o yoerrtrf.o			\
yoerrtwn.o yoesrtwn.o yoesrtm.o yoerdi.o rrtm_init_140gp.o		\
yoerrtbg2.o yoerrtrwt.o surdi.o surrtab.o surrtpk.o surrtftr.o		\
srtm_init.o rrtm_cmbgb1.o rrtm_cmbgb10.o rrtm_cmbgb11.o			\
rrtm_cmbgb12.o rrtm_cmbgb13.o rrtm_cmbgb14.o rrtm_cmbgb15.o		\
rrtm_cmbgb16.o rrtm_cmbgb2.o rrtm_cmbgb3.o rrtm_cmbgb4.o		\
rrtm_cmbgb5.o rrtm_cmbgb6.o rrtm_cmbgb7.o rrtm_cmbgb8.o rrtm_cmbgb9.o	\
rrtm_kgb1.o rrtm_kgb10.o rrtm_kgb11.o rrtm_kgb12.o rrtm_kgb13.o		\
rrtm_kgb14.o rrtm_kgb15.o rrtm_kgb16.o rrtm_kgb2.o rrtm_kgb3.o		\
rrtm_kgb4.o rrtm_kgb5.o rrtm_kgb6.o rrtm_kgb7.o rrtm_kgb8.o		\
rrtm_kgb9.o srtm_cmbgb16.o srtm_cmbgb17.o srtm_cmbgb18.o		\
srtm_cmbgb19.o srtm_cmbgb20.o srtm_cmbgb21.o srtm_cmbgb22.o		\
srtm_cmbgb23.o srtm_cmbgb24.o srtm_cmbgb25.o srtm_cmbgb26.o		\
srtm_cmbgb27.o srtm_cmbgb28.o srtm_cmbgb29.o srtm_kgb16.o		\
srtm_kgb17.o srtm_kgb18.o srtm_kgb19.o srtm_kgb20.o srtm_kgb21.o	\
srtm_kgb22.o srtm_kgb23.o srtm_kgb24.o srtm_kgb25.o srtm_kgb26.o	\
srtm_kgb27.o srtm_kgb28.o srtm_kgb29.o yoerrto1.o yoerrto10.o		\
yoerrto11.o yoerrto12.o yoerrto13.o yoerrto14.o yoerrto15.o		\
yoerrto16.o yoerrto2.o yoerrto3.o yoerrto4.o yoerrto5.o yoerrto6.o	\
yoerrto7.o yoerrto8.o yoerrto9.o susrtm.o modify_wv_continuum.o

SOURCES := $(OBJECTS:.o=.F90)

MAKE_DEPS = ../bin/make_deps.sh
MAKE_DUMMY_INCLUDES = ../bin/make_dummy_includes.sh
MAKE_INCLUDES = ../bin/make_intfbl.1.pl
INCLUDE_DIR = ../include

DEPS_FILE = Makefile_deps

LIBIFSRRTM = ../lib/libifsrrtm.a

all: $(LIBIFSRRTM)

#deps: clean-deps $(DEPS_FILE) dummy_includes
deps: clean-deps $(DEPS_FILE) includes

$(LIBIFSRRTM): $(OBJECTS)
	ar r $(LIBIFSRRTM) $(OBJECTS)

%.o: %.F90
	$(FC) $(FCFLAGS) -c $<

$(DEPS_FILE): $(SOURCES) 
	$(MAKE_DEPS) $(SOURCES) > $(DEPS_FILE)

dummy_includes:
	$(MAKE_DUMMY_INCLUDES) $(SOURCES)

includes:
	LOC_INTFBDIR=$(INCLUDE_DIR) INTFBDIR=$(INCLUDE_DIR) $(MAKE_INCLUDES) $(SOURCES)

count:
	echo $(OBJECTS) | wc

clean:
	rm -f *.o $(LIBIFSRRTM)

dist-clean: clean-deps

clean-deps:
	rm -f Makefile_deps

.PHONY: deps dummy_includes includes clean-deps

include $(DEPS_FILE)

