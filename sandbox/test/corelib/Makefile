# HiQLab
# Copyright (c): Regents of the University of California

include ../../make.inc

#OBJS=	qassembly.o	\
	gaussquad.o	\
	mesh.o		\
	element.o	\
	fieldeval.o	\
	shapes.o	\
	luasupport.o

gaussquad.x: gaussquad.o
	$(CXX) -o gaussquad.x gaussquad.o $(LIBQLAB)

gaussquad.cc:
	octave --eval "gen_gauss_tester('gaussquad.cc')"

.cc.o: 
	$(CXX) $(CXXFLAGS) -I$(SRCDIR)/corelib -c $<

clean:
	rm -f *.o *.x *~ Makefile.bak

depend:
	makedepend -I $(SRCDIR)/corelib -Y ${OBJS:.o=.cc}

