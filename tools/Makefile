include ../make.inc

all: $(MESHVFLTK_BUILD)
	(cd lua; make)
	(cd tolua++; make)
	(cd mwrap; make)
	(cd matexpr; make)
	(cd dsbweb; make)
	(cd ARPACK; make all)
	(cd UMFPACK; make)

meshvfltk.build:
	(cd meshvfltk; make)

clean:
	(cd lua;   make clean)
	(cd tolua++; make clean)
	(cd mwrap; make clean)
	(cd matexpr; make clean)
	(cd dsbweb; make clean)
	(cd ARPACK; make clean; rm -f libarpack.a)
	(cd AMD; make purge)
	(cd UMFPACK; make purge)
	(cd meshvfltk; make clean)

configp:
	cp  config/trilinos-conf      trilinos/conf_hiqlab
	cp  config/petsc-conf         petsc/conf_hiqlab
	cp  config/parmetis-make.inc  tools/ParMetis-3.1/make.inc
	cp  config/superlu-make.inc   tools/SuperLU_3.0/make.inc

allp:
	(cd ParMetis-3.1; make)
	(cd SuperLU_3.0; make)
	(cd trilinos; make)
	(cd petsc; make)

cleanp:
	(cd ParMetis-3.1; make clean)
	(cd SuperLU_3.0; make clean)
	(cd trilinos; make clean)
	(cd petsc; make clean)
