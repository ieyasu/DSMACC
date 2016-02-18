include Makefile.defs
export KPP_HOME=$(PWD)/kpp

DSMACC_SRC = driver.f90 dsmacc.kpp \
	global.inc.f90 rate.inc.f90 util.inc.f90 photolysis.inc.f90 \
	inorganic.kpp organic.kpp depos.kpp

.PHONY: all dsmacc kpp tuv depend check clean distclean

all: kpp dsmacc

dsmacc: bin/dsmacc
bin/dsmacc: src/dsmacc_Main.f90 tuv src/depend.mk Makefile.defs
	cd src && $(MAKE)

src/dsmacc_Main.f90: $(DSMACC_SRC)
	cd stage && ./reprocess.sh

depos.kpp: organic.kpp inorganic.kpp
	idl -e '.run deposition.pro'

kpp: kpp/bin/kpp
kpp/bin/kpp: Makefile.defs
	cd kpp && $(MAKE)

tuv: tuv/libtuv.a
tuv/libtuv.a: tuv/depend.mk
	cd tuv && $(MAKE)

depend: tuv/depend.mk src/depend.mk

src/depend.mk: src/dsmacc_Main.f90
	rm -f src/depend.mk
	cd src && $(MAKE) depend

tuv/depend.mk: $(wildcard tuv/*.f)
	rm -f tuv/depend.mk
	cd tuv && $(MAKE) depend

check: kpp/bin/kpp
	cd test && $(MAKE)

clean:
	cd src && $(MAKE) distclean
	cd test && $(MAKE) clean
	cd tuv && $(MAKE) clean
	cd kpp && $(MAKE) clean

distclean: clean
	@rm -f bin/dsmacc
	cd kpp && $(MAKE) distclean
	@rm -f autom4te.cache config.status config.log Makefile.deps

