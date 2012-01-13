.PHONY: mex

all: serial mex

serial:
	(cd corelib;   make all)
	(cd numeric;   make all)
	(cd element;   make all)
	(cd io;        make all)
	#(cd corelibp;  make all)
	(cd lua;       make all)

mex: serial
	(cd mex;           make all)
	(cd mfiles/pgraph; make mex)

oct: serial
	(cd mex;           make oct)
	(cd mfiles/pgraph; make oct)

trilinos: serial
	(cd trilinos/corelib;	make all)
	(cd trilinos/lua;	make all)

petsc: serial
	(cd petsc/corelib;	make all)
	(cd petsc/lua;		make all)

hybrid: trilinos petsc
	(cd hybrid/corelib;	make all)
	(cd hybrid/lua;		make all)

clean:  clean-serial clean-mex 
cleanp: clean-trilinos clean-petsc clean-hybrid 

clean-serial:
	(cd corelib;   make clean)
	(cd numeric;   make clean)
	(cd element;   make clean)
	(cd io;        make clean)
	(cd lua;       make clean)
	(cd corelibp;  make clean)

clean-mex:
	(cd mex;            make clean)
	(cd mfiles/pgraph;  make clean)

clean-trilinos:
	(cd trilinos/corelib;	make clean)
	(cd trilinos/lua;	make clean)

clean-petsc:
	(cd petsc/corelib;	make clean)
	(cd petsc/lua;		make clean)

clean-hybrid:
	(cd hybrid/corelib;	make clean)
	(cd hybrid/lua;		make clean)

depend:
	(cd corelib;   make depend)
	(cd numeric;   make depend)
	(cd element;   make depend)
	(cd io;        make depend)
	(cd lua;       make depend)
	(cd mex;       make depend)

wc:
	wc -l \
		Makefile \
		corelib/Makefile corelib/*.h corelib/*.cc \
		lua/Makefile lua/*.h lua/*.cc lua/*.pkg \
		mex/Makefile mex/*.h mex/*.cc mex/*.pkg \
		mfiles/*
