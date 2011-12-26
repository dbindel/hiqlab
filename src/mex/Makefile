# HiQLab
# Copyright (c): Regents of the University of California

include ../../make.inc

WRAPHFILES=\
	main.mw \
	luamatlab.mw \
	mesh.mw	\
	sparse.mw \
	tedlinear.mw \
	pzlinear.mw

# These may be used again some day, but currently no code accesses the
# MATLAB mesh building commands.
#	element.mw 
#	meshbuild.mw

QINCL=-I../corelib -I../numeric -I../element -I../io -I../lua 

all: mex

mex: meshmex.$(MEXEXT) Mesh_load.m
oct: meshmex.mex Mesh_load.m

mfiles: Mesh_load.m

Mesh_load.m:
	$(MWRAP) $(WRAPHFILES) -cppcomplex -mb -mex meshmex -c meshmex.cc

meshmex.cc: $(WRAPHFILES)
	$(MWRAP) $(WRAPHFILES) -cppcomplex -mb -mex meshmex -c meshmex.cc

luamex.cc: luamex.pkg
	$(TOLUA) -o luamex.cc -H luamex.h luamex.pkg

meshmex.$(MEXEXT): meshmex.cc mexutil.cc meshstubs.cc luamatlab.cc luamex.cc \
		$(LIBQLAB)
	$(MEX) -g $(QINCL) $(INCLUDE) $^ $(MEXLIBS)

meshmex.mex: meshmex.cc mexutil.cc meshstubs.cc luamatlab.cc luamex.cc \
		$(LIBQLAB)
	echo $(MEXOCT) -g $(QINCL) $(INCLUDE) $^ $(MEXLIBS) > tmp.sh
	`perl octquote.pl tmp.sh`
	rm tmp.sh

clean:
	rm -f *.mex* *.dll *.o meshmex.cc luamex.cc luamex.h *.m *~ Makefile.bak
	rm -f fixup.o fixup.c mex.def

depend: meshmex.cc mexutil.cc meshstubs.cc luamatlab.cc luamex.cc
	makedepend $(QINCL) -Y ${OBJS:.o=.cc} $^


# -- Windows-specific -- build for mingw32

meshmex.dll: meshmex.cc mexutil.cc meshstubs.cc \
	luamatlab.cc luamex.cc $(LIBQLAB) 
	$(CXX) -c $(WMEXFLAGS) $(QINCL) $(INCLUDE) meshmex.cc
	$(CXX) -c $(WMEXFLAGS) $(QINCL) $(INCLUDE) mexutil.cc
	$(CXX) -c $(WMEXFLAGS) $(QINCL) $(INCLUDE) meshstubs.cc
	$(CXX) -c $(WMEXFLAGS) $(QINCL) $(INCLUDE) luamatlab.cc
	$(CXX) -c $(WMEXFLAGS) $(QINCL) $(INCLUDE) luamex.cc
	$(CXX) -shared $(TOOLS)/mex-mingw32/mexlib/mex.def -o meshmex.dll -s \
		meshmex.o mexutil.o meshstubs.o luamatlab.o luamex.o \
		$(LIBQLAB) $(MEXLIBS) $(WMEXLIBS)

# DO NOT DELETE

meshmex.o: ../corelib/hiqlab.h luamatlab.h ../lua/qarray.h
meshmex.o: ../corelib/qassembly.h ../corelib/qcomplex.h ../corelib/qmatrix.h
meshmex.o: ../corelib/element.h ../corelib/fieldeval.h ../corelib/mesh.h
meshmex.o: ../element/material_model.h mexutil.h ../numeric/cscmatrix.h
meshmex.o: meshstubs.h ../numeric/pzlinear.h ../numeric/pmlmode.h
meshmex.o: ../numeric/tedlinear.h ../numeric/coordmatrix.h
mexutil.o: ../numeric/cscmatrix.h ../corelib/qcomplex.h
meshstubs.o: mexutil.h ../numeric/cscmatrix.h ../corelib/qcomplex.h
meshstubs.o: ../corelib/mesh.h ../corelib/element.h ../corelib/qassembly.h
meshstubs.o: ../corelib/qmatrix.h ../corelib/fieldeval.h
meshstubs.o: ../numeric/cscindexer.h ../numeric/mesh_csc_dR.h
meshstubs.o: ../numeric/coordmatrix.h ../numeric/cscassembly.h
meshstubs.o: ../numeric/pzlinear.h
luamatlab.o: luamex.h ../lua/qhelperslua.h ../lua/qassemblylua.h
luamatlab.o: ../lua/qarraylua.h ../lua/cscmatrixlua.h ../lua/meshlua.h
luamatlab.o: ../lua/elementlua.h ../lua/shapeslua.h
luamatlab.o: ../lua/material_modellua.h ../lua/tedlinearlua.h
luamatlab.o: ../lua/pzlinearlua.h ../lua/qarray.h ../corelib/qassembly.h
luamatlab.o: ../corelib/qcomplex.h ../corelib/qmatrix.h
luamatlab.o: ../numeric/cscmatrix.h ../corelib/element.h
luamatlab.o: ../corelib/fieldeval.h ../corelib/mesh.h ../corelib/luasupport.h
luamatlab.o: luamatlab.h ../element/material_model.h
luamex.o: luamatlab.h ../lua/qarray.h ../corelib/qassembly.h
luamex.o: ../corelib/qcomplex.h ../corelib/qmatrix.h ../corelib/element.h
luamex.o: ../corelib/fieldeval.h ../corelib/mesh.h
luamex.o: ../element/material_model.h
