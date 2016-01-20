include ../Makefile.defs

PROG = ../bin/dsmacc
SRCS1 = $(wildcard dsmacc_*.f90)
SRCS2 = constants.f90 
SRCS3 = $(wildcard ../tuv/*.f)
OBJS = $(SRCS1:.f90=.o) $(SRCS2:.f90=.o) 

.SUFFIXES: .o .f90
.DUMMY: all depend clean distclean

# all object files *.o depend on their source files *.f90
# the object files are created with the "-c" compiler option
%.o: %.f90
	$(FC) -o $*.o -c $*.f90 $(F90FLAGS)

all: $(PROG)

# the dependencies depend on the link
# the executable depends on depend and also on all objects
# the executable is created by linking all objects
$(PROG): depend.mk $(OBJS) ../tuv/libtuv.a
	mkdir -p ../bin
	$(FC)  -o $@ $(F90FLAGS) $(OBJS) ../tuv/libtuv.a

../tuv/libtuv.a: $(SRCS3)
	cd ../tuv && $(MAKE)

# update file dependencies
# sfmakedpend comes from http://www.arsc.edu/~kate/Perl
depend: depend.mk
depend.mk: $(SRCS1) $(SRCS3) $(SRCS2)
	../sfmakedepend --file=depend.mk $(SRCS1) $(SRCS3) $(SRCS2)

clean:
	rm -f $(OBJS) *.mod *.log *~ depend.mk.old
	rm -f $(SRCS1)

distclean: clean
	rm -f $(PROG)
	rm -f depend.mk* 
	rm -f *.nc
	rm -f *.dat
	rm -f Makefile_dsmacc dsmacc.map tuvlog.txt ../bin/dsmacc
	rm -f fflags

# list of dependencies (via USE statements)
include depend.mk
