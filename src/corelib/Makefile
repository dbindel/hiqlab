# HiQLab
# Copyright (c): Regents of the University of California

include ../../make.inc

OBJS=	qassembly.o	\
	gaussquad.o	\
	mesh.o		\
	element.o	\
	fieldeval.o	\
	shapes.o	\
	luasupport.o

all: lib
lib: $(LIBQLAB)

$(LIBQLAB): $(OBJS)
	mkdir -p ../../lib
	$(AR) rv $(LIBQLAB) $?
	$(RANLIB) $(LIBQLAB)

.cc.o: 
	$(CXX) $(CXXFLAGS) -c $<

clean:
	rm -f *.o *~ Makefile.bak

depend:
	makedepend -Y ${OBJS:.o=.cc}

# DO NOT DELETE

qassembly.o: qassembly.h qcomplex.h qmatrix.h mesh.h element.h fieldeval.h
gaussquad.o: gaussquad.h
mesh.o: luasupport.h mesh.h qcomplex.h element.h qassembly.h qmatrix.h
mesh.o: fieldeval.h coordsorter.h
element.o: qmatrix.h element.h qassembly.h qcomplex.h fieldeval.h mesh.h
fieldeval.o: fieldeval.h mesh.h qcomplex.h element.h qassembly.h qmatrix.h
shapes.o: qmatrix.h gaussquad.h shapes.h qcomplex.h mesh.h element.h
shapes.o: qassembly.h fieldeval.h
luasupport.o: luasupport.h
