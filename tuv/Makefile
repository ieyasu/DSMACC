include ../Makefile.defs

SRCS = functs.f odo3.f rdetfl.f rxn.f setcld.f setso2.f swphys.f	\
	vptmp.f grids.f odrl.f rdinp.f savout.f setno2.f sphers.f TUV.f		\
	wshift.f la_srb.f orbit.f rdxs.f setaer.f seto2.f swbiol.f vpair.f	\
	numer.f qys.f rtrans.f setalb.f setsnw.f swchem.f vpo3.f
OBJS = functs.o odo3.o rdetfl.o rxn.o setcld.o setso2.o swphys.o	\
	vptmp.o grids.o odrl.o rdinp.o savout.o setno2.o sphers.o TUV.o		\
	wshift.o la_srb.o orbit.o rdxs.o setaer.o seto2.o swbiol.o vpair.o	\
	numer.o qys.o rtrans.o setalb.o setsnw.o swchem.o vpo3.o

.SUFFIXES: .o .f
.DUMMY: all depend clean distclean

.f.o:
	$(FC) $(TUVF90FLAGS) -o $*.o -c $*.f

all: libtuv.a

# XXX get ar command from autoconf
libtuv.a: $(OBJS)
	ar -r $@ $(OBJS)

depend: depend.mk
depend.mk: $(SRCS)
	../bin/sfmakedepend --file=$@ *.f

clean:
	rm -f *.o *.mod lib*.a
	rm -f *.log *~ depend.mk.old

distclean: clean
	rm -f depend.mk* tuvlog.txt

# list of dependencies (via USE statements)
include depend.mk
