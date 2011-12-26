.PHONY: serial mex oct clean clean-serial depend

serial:
	(cd corelib;   make all)
	(cd numeric;   make all)
	(cd element;   make all)
	(cd io;        make all)
	(cd lua;       make all)

mex: serial
	(cd mex;           make all)
	(cd mfiles/pgraph; make mex)

oct: serial
	(cd mex;           make oct)
	(cd mfiles/pgraph; make oct)

clean:  clean-serial clean-mex 

clean-serial:
	(cd corelib;   make clean)
	(cd numeric;   make clean)
	(cd element;   make clean)
	(cd io;        make clean)
	(cd lua;       make clean)

clean-mex:
	(cd mex;            make clean)
	(cd mfiles/pgraph;  make clean)

depend:
	(cd corelib;   make depend)
	(cd numeric;   make depend)
	(cd element;   make depend)
	(cd io;        make depend)
	(cd lua;       make depend)
	(cd mex;       make depend)

