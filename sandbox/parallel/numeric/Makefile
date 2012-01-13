# HiQLab
# Copyright (c): Regents of the University of California

include ../../make.inc

OBJS=	cscmatrix.o	\
	umfmatrix.o	\
	cscindexer.o	\
	cscassembly.o	\
	coordmatrix.o	\
	bcoordmatrix.o	\
	areigs.o 	\
	pmlmode.o	\
        pzlinear.o 	\
	tedlinear.o	\
	mesh_csc_dR.o	\
	meshlabeler.o \
	metispart.o \
	mesh_manager.o \
	csrindexer.o \
	qlapack_util.o

QINCL=-I../corelib -I../corelibp

all: lib
lib: $(LIBQLAB)

$(LIBQLAB): $(OBJS)
	$(AR) rv $(LIBQLAB) $?
	$(RANLIB) $(LIBQLAB)

.cc.o: 
	$(CXX) $(CXXFLAGS) $(QINCL) -c $<

clean:
	rm -f *.a *.o *~ Makefile.bak

depend:
	makedepend $(QINCL) -Y ${OBJS:.o=.cc}

# DO NOT DELETE

cscmatrix.o: cscmatrix.h ../corelib/qcomplex.h
umfmatrix.o: cscmatrix.h ../corelib/qcomplex.h umfmatrix.h
cscindexer.o: ../corelib/qassembly.h ../corelib/qcomplex.h
cscindexer.o: ../corelib/qmatrix.h ../corelib/mesh.h ../corelib/element.h
cscindexer.o: ../corelib/fieldeval.h cscindexer.h cscmatrix.h
cscassembly.o: ../corelib/qcomplex.h ../corelib/mesh.h ../corelib/element.h
cscassembly.o: ../corelib/qassembly.h ../corelib/qmatrix.h
cscassembly.o: ../corelib/fieldeval.h cscmatrix.h cscassembly.h
coordmatrix.o: coordmatrix.h ../corelib/qcomplex.h cscmatrix.h
coordmatrix.o: ../corelib/qassembly.h ../corelib/qmatrix.h
bcoordmatrix.o: bcoordmatrix.h ../corelib/qcomplex.h coordmatrix.h
bcoordmatrix.o: cscmatrix.h ../corelib/qassembly.h ../corelib/qmatrix.h
areigs.o: ../corelib/qcomplex.h areigs.h arinfo.cc
pmlmode.o: ../corelib/qcomplex.h cscmatrix.h umfmatrix.h mesh_csc_dR.h
pmlmode.o: ../corelib/mesh.h ../corelib/element.h ../corelib/qassembly.h
pmlmode.o: ../corelib/qmatrix.h ../corelib/fieldeval.h coordmatrix.h
pmlmode.o: pmlmode.h areigs.h
pzlinear.o: pmlmode.h ../corelib/qcomplex.h cscmatrix.h ../corelib/mesh.h
pzlinear.o: ../corelib/element.h ../corelib/qassembly.h ../corelib/qmatrix.h
pzlinear.o: ../corelib/fieldeval.h pzlinear.h mesh_csc_dR.h coordmatrix.h
pzlinear.o: meshlabeler.h bcoordmatrix.h
tedlinear.o: pmlmode.h ../corelib/qcomplex.h cscmatrix.h ../corelib/mesh.h
tedlinear.o: ../corelib/element.h ../corelib/qassembly.h ../corelib/qmatrix.h
tedlinear.o: ../corelib/fieldeval.h tedlinear.h coordmatrix.h meshlabeler.h
tedlinear.o: mesh_csc_dR.h umfmatrix.h bcoordmatrix.h
mesh_csc_dR.o: mesh_csc_dR.h ../corelib/mesh.h ../corelib/qcomplex.h
mesh_csc_dR.o: ../corelib/element.h ../corelib/qassembly.h
mesh_csc_dR.o: ../corelib/qmatrix.h ../corelib/fieldeval.h coordmatrix.h
mesh_csc_dR.o: cscmatrix.h
meshlabeler.o: ../corelib/mesh.h ../corelib/qcomplex.h ../corelib/element.h
meshlabeler.o: ../corelib/qassembly.h ../corelib/qmatrix.h
meshlabeler.o: ../corelib/fieldeval.h meshlabeler.h
metispart.o: metispart.h
mesh_manager.o: mesh_manager.h ../corelibp/mesh_partition.h ../corelib/mesh.h
mesh_manager.o: ../corelib/qcomplex.h ../corelib/element.h
mesh_manager.o: ../corelib/qassembly.h ../corelib/qmatrix.h
mesh_manager.o: ../corelib/fieldeval.h cscmatrix.h
mesh_manager.o: ../corelibp/mesh_add_block.h ../corelibp/mesh_partitioner.h
mesh_manager.o: cscindexer.h cscassembly.h
csrindexer.o: ../corelib/qassembly.h ../corelib/qcomplex.h
csrindexer.o: ../corelib/qmatrix.h ../corelib/mesh.h ../corelib/element.h
csrindexer.o: ../corelib/fieldeval.h csrindexer.h cscmatrix.h cscindexer.h
qlapack_util.o: qlapack_util.h
