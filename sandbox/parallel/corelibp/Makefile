# HiQLab
# Copyright (c): Regents of the University of California

include ../../make.inc

OBJS=	mesh_partition.o \
	qassemblyp.o \
	adjstructure.o \
	mesh_partitioner_metis.o \
	mesh_partitioner_conform.o \
	mesh_partition_util.o \
	mesh_add_block.o

QINCL=-I../corelib -I../numeric

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

mesh_partition.o: mesh_partition.h ../corelib/mesh.h ../corelib/qcomplex.h
mesh_partition.o: ../corelib/element.h ../corelib/qassembly.h
mesh_partition.o: ../corelib/qmatrix.h ../corelib/fieldeval.h
mesh_partition.o: ../corelib/luasupport.h qassemblyp.h
qassemblyp.o: qassemblyp.h ../corelib/qcomplex.h ../corelib/qmatrix.h
qassemblyp.o: ../corelib/qassembly.h mesh_partition.h ../corelib/mesh.h
qassemblyp.o: ../corelib/element.h ../corelib/fieldeval.h mesh_add_block.h
adjstructure.o: adjstructure.h ../corelib/mesh.h ../corelib/qcomplex.h
adjstructure.o: ../corelib/element.h ../corelib/qassembly.h
adjstructure.o: ../corelib/qmatrix.h ../corelib/fieldeval.h
adjstructure.o: ../numeric/cscindexer.h ../numeric/cscmatrix.h
adjstructure.o: ../numeric/metispart.h
mesh_partitioner_metis.o: mesh_partitioner_metis.h mesh_partitioner.h
mesh_partitioner_metis.o: ../corelib/mesh.h ../corelib/qcomplex.h
mesh_partitioner_metis.o: ../corelib/element.h ../corelib/qassembly.h
mesh_partitioner_metis.o: ../corelib/qmatrix.h ../corelib/fieldeval.h
mesh_partitioner_metis.o: ../numeric/cscindexer.h ../numeric/cscmatrix.h
mesh_partitioner_metis.o: ../numeric/metispart.h
mesh_partitioner_conform.o: mesh_partitioner_conform.h mesh_partitioner.h
mesh_partitioner_conform.o: ../corelib/mesh.h ../corelib/qcomplex.h
mesh_partitioner_conform.o: ../corelib/element.h ../corelib/qassembly.h
mesh_partitioner_conform.o: ../corelib/qmatrix.h ../corelib/fieldeval.h
mesh_partition_util.o: mesh_partition_util.h mesh_partitioner.h
mesh_partition_util.o: mesh_partition.h ../corelib/mesh.h
mesh_partition_util.o: ../corelib/qcomplex.h ../corelib/element.h
mesh_partition_util.o: ../corelib/qassembly.h ../corelib/qmatrix.h
mesh_partition_util.o: ../corelib/fieldeval.h
mesh_add_block.o: mesh_add_block.h ../corelib/mesh.h ../corelib/qcomplex.h
mesh_add_block.o: ../corelib/element.h ../corelib/qassembly.h
mesh_add_block.o: ../corelib/qmatrix.h ../corelib/fieldeval.h
mesh_add_block.o: ../corelib/shapes.h qassemblyp.h
