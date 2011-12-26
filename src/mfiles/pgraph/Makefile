include ../../../make.inc

all: mex

.cc.o:
	$(MEX) -c $*.cc

mex: pgraphmex.cc
	$(MEX) pgraphmex.cc

oct: pgraphmex.cc
	mkoctfile --mex pgraphmex.cc

pgraphmex.cc: pgraph.mw 
	$(MWRAP) -cppcomplex -mex pgraphmex -c pgraphmex.cc -mb pgraph.mw

clean:
	rm -f *.o *.mex* *~
	rm -f pgraphmex.cc
	(mw=$(MWRAP); FILES=`$$mw -cppcomplex -list pgraph.mw`; rm -f $$FILES)

