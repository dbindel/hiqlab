include ../../../make.inc
include $(PETSC_DIR)/bmake/$(PETSC_ARCH)/petscconf
include ../make.inc

LIBQLAB_PETSC= $(LIBDIR)/libqlab_petsc.a

OBJS=	qpetsc.o \
	qpetsclua.o \
	$(USE_SLEPC_OBJS) \
	qpassembly_petsclua.o \
	qpetsc_pclua.o \
	qpetsc_meshlua.o \
	qpetsc_pc_luastubs.o \
	qpetsc_jdqzlua.o

SLEPC_OBJS=qslepclua.o

QINCL= -I../../corelib -I../../lua -I../corelib -I../../corelibp

all: lib exe
lib: $(LIBQLAB_PETSC)
exe: hiqlab

$(LIBQLAB_PETSC): $(OBJS)
	$(AR) rv $(LIBQLAB_PETSC) $?
	$(RANLIB) $(LIBQLAB_PETSC)

hiqlab: hiqlab.o $(LIBQLAB_PETSC) $(LIBQLAB)
	$(CXX) -o hiqlab hiqlab.o $(LIBQLAB_PETSC) $(LIBQLAB) \
		${SLEPC_LIBS} ${PETSC_LIBS} ${PETSC_KSP_LIB} \
		${HYPRE_LIBS} \
                ${SPOOLES_LIB} \
                ${MUMPS_LIB}\
                ${SCALAPACK_LIB} ${BLACS_LIB} \
		${PROMETHEUS_LIBS} ${METIS_LIBS} \
		${SUPERLU_LIB} \
		${SUPERLU_DIST_LIB} \
		${X_LIBS} \
		$(READLINE_LIBS) $(FLTK_LIBS) $(LIBS) 

hiqlab.o: hiqlab.cc $(OBJS)

.cc.o :
	$(CXX) -c $(CXXFLAGS) $(QINCL) ${PETSC_INCLUDE} ${SLEPC_INCLUDE} $< 

qpetsclua.cc: qpetsc.pkg
	$(TOLUA) -o qpetsclua.cc -H qpetsclua.h qpetsc.pkg

qpassembly_petsclua.cc: qpassembly_petsc.pkg
	$(TOLUA) -o qpassembly_petsclua.cc -H qpassembly_petsclua.h \
		qpassembly_petsc.pkg

qslepclua.cc: qslepc.pkg
	$(TOLUA) -o qslepclua.cc -H qslepclua.h qslepc.pkg

qpetsc_pclua.cc: qpetsc_pc.pkg
	$(TOLUA) -o qpetsc_pclua.cc -H qpetsc_pclua.h qpetsc_pc.pkg

qpetsc_meshlua.cc: qpetsc_mesh.pkg
	$(TOLUA) -o qpetsc_meshlua.cc -H qpetsc_meshlua.h qpetsc_mesh.pkg

qpetsc_jdqzlua.cc: qpetsc_jdqz.pkg
	$(TOLUA) -o qpetsc_jdqzlua.cc -H qpetsc_jdqzlua.h qpetsc_jdqz.pkg

depend: ${OBJS:.o=.cc} qslepclua.cc
	makedepend $(QINCL) -Y ${OBJS:.o=.cc} qslepclua.cc

clean:
	${RM} -f hiqlab qpetsclua* qpetsclua* \
		qpassembly_petsclua* qslepclua* qpetsc_pclua* \
		qpetsc_jdqzlua* \
		qpetsc_meshlua* *.o Makefile.bak

# DO NOT DELETE

qpetsc.o: ../../corelib/mesh.h ../../corelib/qcomplex.h
qpetsc.o: ../../corelib/element.h ../../corelib/qassembly.h
qpetsc.o: ../../corelib/qmatrix.h ../../corelib/fieldeval.h qpetsc.h
qpetsclua.o: ../../corelib/mesh.h ../../corelib/qcomplex.h
qpetsclua.o: ../../corelib/element.h ../../corelib/qassembly.h
qpetsclua.o: ../../corelib/qmatrix.h ../../corelib/fieldeval.h qpetsc.h
qpassembly_petsclua.o: ../../corelib/mesh.h ../../corelib/qcomplex.h
qpassembly_petsclua.o: ../../corelib/element.h ../../corelib/qassembly.h
qpassembly_petsclua.o: ../../corelib/qmatrix.h ../../corelib/fieldeval.h
qpassembly_petsclua.o: ../corelib/qpassembly_petsc.h
qpassembly_petsclua.o: ../corelib/qpassembly_part_petsc.h
qpassembly_petsclua.o: ../../corelibp/mesh_partition_util.h
qpassembly_petsclua.o: ../../corelibp/mesh_partitioner.h
qpassembly_petsclua.o: ../../corelibp/mesh_partition.h
qpetsc_pclua.o: ../corelib/qpetsc_pc.h ../../corelib/mesh.h
qpetsc_pclua.o: ../../corelib/qcomplex.h ../../corelib/element.h
qpetsc_pclua.o: ../../corelib/qassembly.h ../../corelib/qmatrix.h
qpetsc_pclua.o: ../../corelib/fieldeval.h qpetsc_pc_luastubs.h
qpetsc_pclua.o: ../corelib/qpetscmg.h
qpetsc_meshlua.o: ../corelib/petsc_mesh.h ../../corelib/mesh.h
qpetsc_meshlua.o: ../../corelib/qcomplex.h ../../corelib/element.h
qpetsc_meshlua.o: ../../corelib/qassembly.h ../../corelib/qmatrix.h
qpetsc_meshlua.o: ../../corelib/fieldeval.h
qpetsc_pc_luastubs.o: qpetsc_pc_luastubs.h ../corelib/qpetscmg.h
qpetsc_jdqzlua.o: ../corelib/petsc_jdqz.h ../corelib/petsc_jdqz_structure.h
qpetsc_jdqzlua.o: ../corelib/matshellab.h
qslepclua.o: ../../corelib/mesh.h ../../corelib/qcomplex.h
qslepclua.o: ../../corelib/element.h ../../corelib/qassembly.h
qslepclua.o: ../../corelib/qmatrix.h ../../corelib/fieldeval.h
qslepclua.o: ../corelib/qslepc.h
