include ../../../make.inc
include $(PETSC_DIR)/bmake/$(PETSC_ARCH)/petscconf
include ../../../make.inc
include ../make.inc

LIBQLAB_PETSC= $(LIBDIR)/libqlab_petsc.a

OBJS=	qpassembly_petsc.o \
	$(USE_SLEPC_OBJS) \
	qpetsc_pc.o \
	petsc_mesh.o \
	qpassembly_part_petsc.o \
	kaczmarz.o \
	projectedpc.o \
	mgfuncq.o \
	matshellab.o \
	petsc_jdqz_structure.o \
	petsc_jdqz_util.o \
	petsc_jdqz.o

SLEPC_OBJS=qslepc.o

QINCL=-I../../corelib -I../../numeric -I../../corelibp

all: lib
lib: $(LIBQLAB_PETSC)

$(LIBQLAB_PETSC): $(OBJS)
	$(AR) rv $(LIBQLAB_PETSC) $?
	$(RANLIB) $(LIBQLAB_PETSC)

.c.o :
	$(CXX) -c $(CXXFLAGS) $(QINCL) $(PETSC_INCLUDE) $(SLEPC_INCLUDE) $<

.cc.o :
	$(CXX) -c $(CXXFLAGS) $(QINCL) $(PETSC_INCLUDE) $(SLEPC_INCLUDE) $< 

depend:
	makedepend $(QINCL) -Y ${OBJS:.o=.cc} qslepc.cc 

clean:
	${RM} *.o libqlab_petsc.a Makefile.bak

# DO NOT DELETE

qpassembly_petsc.o: ../../corelib/mesh.h ../../corelib/qcomplex.h
qpassembly_petsc.o: ../../corelib/element.h ../../corelib/qassembly.h
qpassembly_petsc.o: ../../corelib/qmatrix.h ../../corelib/fieldeval.h
qpassembly_petsc.o: ../../corelibp/mesh_add_block.h qpassembly_petsc.h
qpassembly_petsc.o: ../../numeric/csrindexer.h ../../numeric/cscmatrix.h
qpassembly_petsc.o: ../../numeric/cscindexer.h
qpetsc_pc.o: qpetsc_pc.h ../../corelib/mesh.h ../../corelib/qcomplex.h
qpetsc_pc.o: ../../corelib/element.h ../../corelib/qassembly.h
qpetsc_pc.o: ../../corelib/qmatrix.h ../../corelib/fieldeval.h
petsc_mesh.o: petsc_mesh.h ../../corelib/mesh.h ../../corelib/qcomplex.h
petsc_mesh.o: ../../corelib/element.h ../../corelib/qassembly.h
petsc_mesh.o: ../../corelib/qmatrix.h ../../corelib/fieldeval.h
petsc_mesh.o: qpassembly_petsc.h
qpassembly_part_petsc.o: ../../corelibp/mesh_partition.h ../../corelib/mesh.h
qpassembly_part_petsc.o: ../../corelib/qcomplex.h ../../corelib/element.h
qpassembly_part_petsc.o: ../../corelib/qassembly.h ../../corelib/qmatrix.h
qpassembly_part_petsc.o: ../../corelib/fieldeval.h
qpassembly_part_petsc.o: ../../corelibp/mesh_add_block.h
qpassembly_part_petsc.o: qpassembly_part_petsc.h qpassembly_petsc.h
qpassembly_part_petsc.o: ../../numeric/csrindexer.h ../../numeric/cscmatrix.h
qpassembly_part_petsc.o: ../../numeric/cscindexer.h
qslepc.o: ../../corelib/mesh.h ../../corelib/qcomplex.h
qslepc.o: ../../corelib/element.h ../../corelib/qassembly.h
qslepc.o: ../../corelib/qmatrix.h ../../corelib/fieldeval.h
qslepc.o: qpassembly_petsc.h qslepc.h qpetsc_pc.h
