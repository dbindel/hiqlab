# HiQLab
# Copyright (c): Regents of the University of California

include ../../../make.inc

LIBQLAB=
#../../corelib/libqlab.a

MEXFILES=metismex.$(MEXEXT)

WRAPHFILES=\
	metis.mw

QINCL=

all: mex
mex: $(MEXFILES)
mfiles:

metismex.cc: $(WRAPHFILES)
	$(MWRAP) $(WRAPHFILES) -cppcomplex -mb -mex metismex -c metismex.cc

metismex.$(MEXEXT): metismex.cc metisstubs.cc $(LIBQLAB)

%.$(MEXEXT): %.cc 
	$(MEX) -g $(QINCL) $(METIS_INCLUDE) $(INCLUDE) $^ $(LIBQLAB) $(METIS_LIBS) $(MEXLIBS)

clean:
	rm -f *.$(MEXEXT) metismex.cc *.m *~ Makefile.bak
	rm -f fixup.o fixup.c mex.def

depend: metismex.cc metisstubs.cc 
	makedepend $(QINCL) -Y ${OBJS:.o=.cc} $^


# -- Windows-specific -- build for mingw32

meshmex.dll: metismex.cc metisstubs.cc
	$(CXX) -c $(WMEXFLAGS) $(QINCL) $(INCLUDE) metismex.cc
	$(CXX) -c $(WMEXFLAGS) $(QINCL) $(INCLUDE) metisstubs.cc
	$(CXX) -shared $(TOOLS)/mex-mingw32/mexlib/mex.def -o metismex.dll -s \
		metismex.o metisstubs.o \
		$(MEXLIBS) $(WMEXLIBS)

# DO NOT DELETE

metismex.o: metisstubs.h
metisstubs.o: metisstubs.h
