# HiQLab
# Copyright (c): Regents of the University of California

include ../../make.inc

OBJS=	qhelperslua.o	\
	cscmatrixlua.o	\
	meshlua.o	\
	elementlua.o	\
	dxfilelua.o	\
	qassemblylua.o  \
	qarray.o	\
	qarraylua.o	\
	qiarray.o	\
	qiarraylua.o	\
	qlapack.o	\
	qlapacklua.o	\
	shapeslua.o	\
	leigslua.o      \
	pmlmodelua.o    \
	tedlinearlua.o  \
        material_modellua.o\
	dirstuff.o	\
	leigs.o         \
        pzlinearlua.o	\
        testlua.o \
	mesh_partitionlua.o \
	mesh_partitionerlua.o \
	mesh_add_blocklua.o \
	mesh_partitioner_luastubs.o \
	mesh_managerlua.o

QINCL=-I../corelib -I../numeric -I../element -I../io -I../corelibp

all: lib exe
lib: $(LIBQLAB)
exe: $(HIQ)

$(LIBQLAB): $(OBJS)
	$(AR) rv $(LIBQLAB) $?
	$(RANLIB) $(LIBQLAB)

$(HIQ): hiqlab.o $(LIBQLAB)
	$(CXX) -o $(HIQ) hiqlab.o $(LIBQLAB) $(SERIAL_METIS_LIBS) \
		$(READLINE_LIBS) $(MEXLIBS)

hiqlab.o: hiqlab.cc
	$(CXX) $(CXXFLAGS) $(FLTK_CXXFLAGS) $(QINCL) -c hiqlab.cc

.cc.o: 
	$(CXX) $(CXXFLAGS) $(QINCL) -c $<

%lua.cc %lua.h: %.pkg
	$(TOLUA) -o $*lua.cc -H $*lua.h $<

meshlua.cc: mesh_m4.pkg meshbc.lua meshscale.lua meshbuilders.lua
	m4 mesh_m4.pkg > mesh.pkg
	$(TOLUA) -o meshlua.cc -H meshlua.h mesh.pkg
	rm mesh.pkg

mesh_partitionlua.cc: mesh_partition_m4.pkg
	m4 mesh_partition_m4.pkg > mesh_partition.pkg
	$(TOLUA) -o mesh_partitionlua.cc -H mesh_partitionlua.h mesh_partition.pkg
	rm mesh_partition.pkg
clean:
	rm -f *lua.* *.o *~ Makefile.bak

depend: ${OBJS:.o=.cc}
	makedepend $(QINCL) -Y ${OBJS:.o=.cc}

# DO NOT DELETE

cscmatrixlua.o: ../numeric/cscmatrix.h ../corelib/qcomplex.h
cscmatrixlua.o: ../numeric/umfmatrix.h ../numeric/coordmatrix.h
cscmatrixlua.o: ../corelib/qassembly.h ../corelib/qmatrix.h qarray.h
meshlua.o: ../corelib/mesh.h ../corelib/qcomplex.h ../corelib/element.h
meshlua.o: ../corelib/qassembly.h ../corelib/qmatrix.h ../corelib/fieldeval.h
meshlua.o: ../numeric/mesh_csc_dR.h ../numeric/coordmatrix.h
meshlua.o: ../numeric/cscmatrix.h qarray.h
elementlua.o: ../corelib/element.h ../corelib/qassembly.h
elementlua.o: ../corelib/qcomplex.h ../corelib/qmatrix.h
elementlua.o: ../corelib/fieldeval.h ../corelib/mesh.h
elementlua.o: ../element/couple_em2d.h ../element/pmlelement.h
elementlua.o: ../element/pmlfunc.h ../element/elastic2d.h
elementlua.o: ../element/pmlscalar1d.h ../element/pmlscalar2d.h
elementlua.o: ../element/pmlscalar3d.h ../element/pmlscalarax.h
elementlua.o: ../element/pmlelastic2d.h ../element/pmlelastic3d.h
elementlua.o: ../element/pmlelasticax.h ../element/pmlelastictax.h
elementlua.o: ../element/pmlelastic2d_te.h ../element/pmlelasticax_te.h
elementlua.o: ../element/pmlelastic3d_te.h ../element/pmlelastic2d_pz.h
elementlua.o: ../element/pmlelastic2hd_pz.h ../element/pmlelastic3d_pz.h
elementlua.o: ../element/tie_field.h ../element/tie_field2.h
elementlua.o: ../element/model_circuit.h
dxfilelua.o: ../io/dxfile.h ../corelib/mesh.h ../corelib/qcomplex.h
dxfilelua.o: ../corelib/element.h ../corelib/qassembly.h ../corelib/qmatrix.h
dxfilelua.o: ../corelib/fieldeval.h
qassemblylua.o: ../corelib/qassembly.h ../corelib/qcomplex.h
qassemblylua.o: ../corelib/qmatrix.h
qarray.o: qarray.h
qarraylua.o: qarray.h
qiarray.o: qiarray.h
qiarraylua.o: qiarray.h
qlapack.o: qlapack.h qarray.h qiarray.h
qlapacklua.o: qlapack.h qarray.h qiarray.h
shapeslua.o: ../corelib/shapes.h ../corelib/qcomplex.h
leigslua.o: leigs.h qarray.h ../numeric/areigs.h ../corelib/qcomplex.h
pmlmodelua.o: ../corelib/mesh.h ../corelib/qcomplex.h ../corelib/element.h
pmlmodelua.o: ../corelib/qassembly.h ../corelib/qmatrix.h
pmlmodelua.o: ../corelib/fieldeval.h ../numeric/pmlmode.h
pmlmodelua.o: ../numeric/cscmatrix.h qarray.h
tedlinearlua.o: ../corelib/mesh.h ../corelib/qcomplex.h ../corelib/element.h
tedlinearlua.o: ../corelib/qassembly.h ../corelib/qmatrix.h
tedlinearlua.o: ../corelib/fieldeval.h ../numeric/tedlinear.h
tedlinearlua.o: ../numeric/coordmatrix.h ../numeric/cscmatrix.h qarray.h
material_modellua.o: ../element/material_model.h
dirstuff.o: dirstuff.h
leigs.o: ../corelib/qcomplex.h leigs.h qarray.h ../numeric/areigs.h
leigs.o: ../corelib/luasupport.h
pzlinearlua.o: ../numeric/pzlinear.h ../corelib/qcomplex.h ../corelib/mesh.h
pzlinearlua.o: ../corelib/element.h ../corelib/qassembly.h
pzlinearlua.o: ../corelib/qmatrix.h ../corelib/fieldeval.h
testlua.o: ../element/pmlfunc.h
mesh_partitionlua.o: ../corelib/mesh.h ../corelib/qcomplex.h
mesh_partitionlua.o: ../corelib/element.h ../corelib/qassembly.h
mesh_partitionlua.o: ../corelib/qmatrix.h ../corelib/fieldeval.h qarray.h
mesh_partitionlua.o: ../corelibp/mesh_partition.h
mesh_partitionlua.o: ../corelibp/mesh_partitioner.h
mesh_partitionlua.o: ../corelibp/mesh_partition_util.h
mesh_partitionlua.o: ../corelibp/mesh_add_block.h
mesh_partitionerlua.o: ../corelib/mesh.h ../corelib/qcomplex.h
mesh_partitionerlua.o: ../corelib/element.h ../corelib/qassembly.h
mesh_partitionerlua.o: ../corelib/qmatrix.h ../corelib/fieldeval.h
mesh_partitionerlua.o: ../corelibp/adjstructure.h
mesh_partitionerlua.o: ../corelibp/mesh_partition.h ../numeric/mesh_manager.h
mesh_partitionerlua.o: ../numeric/cscmatrix.h ../corelibp/mesh_add_block.h
mesh_partitionerlua.o: ../corelibp/mesh_partitioner.h qarray.h
mesh_partitionerlua.o: ../corelibp/mesh_partitioner_metis.h
mesh_partitionerlua.o: ../corelibp/mesh_partitioner_conform.h
mesh_partitionerlua.o: mesh_partitioner_luastubs.h
mesh_add_blocklua.o: ../corelibp/mesh_add_block.h ../corelib/mesh.h
mesh_add_blocklua.o: ../corelib/qcomplex.h ../corelib/element.h
mesh_add_blocklua.o: ../corelib/qassembly.h ../corelib/qmatrix.h
mesh_add_blocklua.o: ../corelib/fieldeval.h ../numeric/mesh_manager.h
mesh_add_blocklua.o: ../corelibp/mesh_partition.h ../numeric/cscmatrix.h
mesh_add_blocklua.o: ../corelibp/mesh_partitioner.h
mesh_partitioner_luastubs.o: mesh_partitioner_luastubs.h
mesh_partitioner_luastubs.o: ../corelibp/mesh_add_block.h ../corelib/mesh.h
mesh_partitioner_luastubs.o: ../corelib/qcomplex.h ../corelib/element.h
mesh_partitioner_luastubs.o: ../corelib/qassembly.h ../corelib/qmatrix.h
mesh_partitioner_luastubs.o: ../corelib/fieldeval.h
mesh_partitioner_luastubs.o: ../corelibp/mesh_partitioner_metis.h
mesh_partitioner_luastubs.o: ../corelibp/mesh_partitioner.h
mesh_partitioner_luastubs.o: ../corelibp/mesh_partitioner_conform.h
mesh_managerlua.o: qarray.h ../numeric/mesh_manager.h
mesh_managerlua.o: ../corelibp/mesh_partition.h ../corelib/mesh.h
mesh_managerlua.o: ../corelib/qcomplex.h ../corelib/element.h
mesh_managerlua.o: ../corelib/qassembly.h ../corelib/qmatrix.h
mesh_managerlua.o: ../corelib/fieldeval.h ../numeric/cscmatrix.h
mesh_managerlua.o: ../corelibp/mesh_add_block.h
mesh_managerlua.o: ../corelibp/mesh_partitioner.h
