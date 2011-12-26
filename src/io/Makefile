# HiQLab
# Copyright (c): Regents of the University of California

include ../../make.inc

OBJS=	dxfile.o

QINCL= -I../corelib

all: lib
lib: $(LIBQLAB)

$(LIBQLAB): $(OBJS)
	$(AR) rv $(LIBQLAB) $?
	$(RANLIB) $(LIBQLAB)

.cc.o: 
	$(CXX) $(CXXFLAGS) $(QINCL) -c $<

dxfile.o: dxfile.cc dxfile.h
	$(CXX) $(CXXFLAGS) $(QINCL) -DHIQINT32=$(HIQINT32) -c $<

clean:
	rm -f *.a *.o *~ Makefile.bak

depend:
	makedepend -Y ${OBJS:.o=.cc}

# DO NOT DELETE

dxfile.o: ../../config.h dxfile.h
