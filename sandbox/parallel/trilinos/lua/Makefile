include ../../../make.inc
include ../corelib/make_trilinos.inc

LIBQLAB=../../corelib/libqlab.a
LIBQLAB_TL = ../corelib/libqlab_trilinos.a

OBJS=	qpassembly_trilinoslua.o	\
	qtrilinos_epetra.o	\
	qtrilinos_epetralua.o	\
	qtrilinos_amesos.o \
	qtrilinos_amesoslua.o \
	qtrilinos_aztecoo.o \
	qtrilinos_aztecoolua.o \
	qtrilinos_teuchos.o \
	qtrilinos_teuchoslua.o \
	qtrilinos_epetraextlua.o \
	qtrilinos_anasazilua.o \
	qtrilinos_arpacklua.o \
	qtrilinos_ml.o \
	qtrilinos_mllua.o \
	qtrilinos_ifpack.o \
	qtrilinos_ifpacklua.o \
	qtrilinos_komplexlua.o \
	qtrilinos_meshlua.o 

all: lib exe
lib: $(LIBQLAB_TL)
exe: hiqlab

$(LIBQLAB_TL): $(OBJS)
	$(AR) rv $(LIBQLAB_TL) $?
	$(RANLIB) $(LIBQLAB_TL)

hiqlab: hiqlab.o $(LIBQLAB_TL) $(LIBQLAB)
	$(CXX) -o hiqlab hiqlab.o $(LIBQLAB_TL) $(LIBQLAB) $(LIB_TL) \
		$(LIBS) $(F90_LIBS) $(READLINE_LIBS) $(FLTK_LIBS)

hiqlab.o: hiqlab.cc $(OBJS)
	$(CXX) $(CXXFLAGS) $(FLTK_CXXFLAGS) \
		-I../../corelib -I../../lua -I../corelib \
		$(INCLUDE_TL) -c hiqlab.cc

.cc.o: 
	$(CXX) -DHAVE_CONFIG_H $(CXXFLAGS) 
		-I../../corelib -I../corelib $(INCLUDE_TL) -c $<

testlua.cc: test.pkg
	$(TOLUA) -o testlua.cc -H testlua.h test.pkg

qpassembly_trilinoslua.cc: qpassembly_trilinos.pkg
	$(TOLUA) -o qpassembly_trilinoslua.cc -H qpassembly_trilinoslua.h \
		qpassembly_trilinos.pkg

qtrilinos_epetralua.cc: qtrilinos_epetra.pkg
	$(TOLUA) -o qtrilinos_epetralua.cc -H qtrilinos_epetralua.h \
		qtrilinos_epetra.pkg

qtrilinos_amesoslua.cc: qtrilinos_amesos.pkg
	$(TOLUA) -o qtrilinos_amesoslua.cc -H qtrilinos_amesoslua.h \
		qtrilinos_amesos.pkg

qtrilinos_aztecoolua.cc: qtrilinos_aztecoo.pkg
	$(TOLUA) -o qtrilinos_aztecoolua.cc -H qtrilinos_aztecoolua.h \
		qtrilinos_aztecoo.pkg

qtrilinos_teuchoslua.cc: qtrilinos_teuchos_m4.pkg
	m4 qtrilinos_teuchos_m4.pkg > qtrilinos_teuchos.pkg
	$(TOLUA) -o qtrilinos_teuchoslua.cc -H qtrilinos_teuchoslua.h \
		qtrilinos_teuchos.pkg
	rm qtrilinos_teuchos.pkg

qtrilinos_epetraextlua.cc: qtrilinos_epetraext.pkg
	$(TOLUA) -o qtrilinos_epetraextlua.cc -H qtrilinos_epetraextlua.h \
		qtrilinos_epetraext.pkg

qtrilinos_anasazilua.cc: qtrilinos_anasazi.pkg
	$(TOLUA) -o qtrilinos_anasazilua.cc -H qtrilinos_anasazilua.h \
		qtrilinos_anasazi.pkg

qtrilinos_arpacklua.cc: qtrilinos_arpack.pkg
	$(TOLUA) -o qtrilinos_arpacklua.cc -H qtrilinos_arpacklua.h \
		qtrilinos_arpack.pkg

qtrilinos_mllua.cc: qtrilinos_ml.pkg
	$(TOLUA) -o qtrilinos_mllua.cc -H qtrilinos_mllua.h qtrilinos_ml.pkg

qtrilinos_ifpacklua.cc: qtrilinos_ifpack.pkg
	$(TOLUA) -o qtrilinos_ifpacklua.cc -H qtrilinos_ifpacklua.h \
		qtrilinos_ifpack.pkg

qtrilinos_komplexlua.cc: qtrilinos_komplex.pkg
	$(TOLUA) -o qtrilinos_komplexlua.cc -H qtrilinos_komplexlua.h \
		qtrilinos_komplex.pkg

qtrilinos_meshlua.cc: qtrilinos_mesh.pkg
	$(TOLUA) -o qtrilinos_meshlua.cc -H qtrilinos_meshlua.h \
		qtrilinos_mesh.pkg

depend: ${OBJS:.o=.cc}
	makedepend -DHAVE_CONFIG_H \
		-I../../corelib -I../corelib -Y ${OBJS:.o=.cc}

clean:
	rm -f *lua.* *.o hiqlab *~ Makefile.bak

# DO NOT DELETE

qpassembly_trilinoslua.o: ../../corelib/mesh.h ../../corelib/qcomplex.h
qpassembly_trilinoslua.o: ../../corelib/element.h ../../corelib/qassembly.h
qpassembly_trilinoslua.o: ../../corelib/qmatrix.h ../../corelib/coordmatrix.h
qpassembly_trilinoslua.o: ../../corelib/cscmatrix.h
qpassembly_trilinoslua.o: ../../corelib/coordmatrix_g.h
qpassembly_trilinoslua.o: ../corelib/qpassembly_trilinos.h
qpassembly_trilinoslua.o: ../corelib/trilinos_epetra_vector.h
qpassembly_trilinoslua.o: ../corelib/trilinos_epetra_matrix.h
qpassembly_trilinoslua.o: ../corelib/trilinos_epetra_operator.h
qtrilinos_epetra.o: qtrilinos_epetra.h ../corelib/trilinos_epetra_vector.h
qtrilinos_epetra.o: ../../corelib/qcomplex.h
qtrilinos_epetralua.o: qtrilinos_epetra.h ../corelib/trilinos_epetra_vector.h
qtrilinos_epetralua.o: ../../corelib/qcomplex.h
qtrilinos_epetralua.o: ../corelib/trilinos_epetra_linearproblem.h
qtrilinos_epetralua.o: ../corelib/trilinos_epetra_matrix.h
qtrilinos_epetralua.o: ../corelib/trilinos_epetra_operator.h
qtrilinos_amesos.o: qtrilinos_amesos.h
qtrilinos_amesos.o: ../corelib/trilinos_epetra_linearproblem.h
qtrilinos_amesos.o: ../corelib/trilinos_epetra_matrix.h
qtrilinos_amesos.o: ../corelib/trilinos_epetra_vector.h
qtrilinos_amesos.o: ../../corelib/qcomplex.h
qtrilinos_amesos.o: ../corelib/trilinos_epetra_operator.h
qtrilinos_amesos.o: ../corelib/trilinos_amesos_base.h
qtrilinos_amesoslua.o: ../corelib/Amesos_Operator_Complex.h
qtrilinos_amesoslua.o: ../corelib/trilinos_epetra_vector.h
qtrilinos_amesoslua.o: ../../corelib/qcomplex.h
qtrilinos_amesoslua.o: ../corelib/trilinos_epetra_matrix.h
qtrilinos_amesoslua.o: ../corelib/trilinos_epetra_operator.h
qtrilinos_amesoslua.o: ../corelib/trilinos_epetra_linearproblem.h
qtrilinos_amesoslua.o: ../corelib/trilinos_amesos_base.h
qtrilinos_amesoslua.o: ../corelib/Amesos_Operator.h qtrilinos_amesos.h
qtrilinos_amesoslua.o: qtrilinos_epetra.h
qtrilinos_aztecoo.o: qtrilinos_aztecoo.h
qtrilinos_aztecoolua.o: ../corelib/AztecOO_Operator_Komplex.h
qtrilinos_aztecoolua.o: ../corelib/trilinos_epetra_vector.h
qtrilinos_aztecoolua.o: ../../corelib/qcomplex.h
qtrilinos_aztecoolua.o: ../corelib/trilinos_epetra_matrix.h
qtrilinos_aztecoolua.o: ../corelib/trilinos_epetra_operator.h
qtrilinos_aztecoolua.o: ../corelib/trilinos_epetra_linearproblem.h
qtrilinos_aztecoolua.o: ../corelib/trilinos_super_matrix.h
qtrilinos_aztecoolua.o: ../corelib/trilinos_indexmap.h
qtrilinos_aztecoolua.o: ../corelib/trilinos_super_vector.h
qtrilinos_aztecoolua.o: qtrilinos_aztecoo.h
qtrilinos_teuchos.o: qtrilinos_teuchos.h
qtrilinos_teuchoslua.o: qtrilinos_teuchos.h
qtrilinos_anasazilua.o: ../../corelib/mesh.h ../../corelib/qcomplex.h
qtrilinos_anasazilua.o: ../../corelib/element.h ../../corelib/qassembly.h
qtrilinos_anasazilua.o: ../../corelib/qmatrix.h ../../corelib/coordmatrix.h
qtrilinos_anasazilua.o: ../../corelib/cscmatrix.h
qtrilinos_anasazilua.o: ../../corelib/coordmatrix_g.h qtrilinos_anasazi.h
qtrilinos_anasazilua.o: ../corelib/trilinos_anasazi.h
qtrilinos_anasazilua.o: ../corelib/trilinos_epetra_matrix.h
qtrilinos_anasazilua.o: ../corelib/trilinos_epetra_vector.h
qtrilinos_anasazilua.o: ../corelib/trilinos_epetra_operator.h
qtrilinos_anasazilua.o: qtrilinos_epetra.h
qtrilinos_arpacklua.o: qtrilinos_epetra.h ../corelib/trilinos_epetra_vector.h
qtrilinos_arpacklua.o: ../../corelib/qcomplex.h ../corelib/trilinos_arpack.h
qtrilinos_arpacklua.o: ../../corelib/cscmatrix.h ../../corelib/mesh.h
qtrilinos_arpacklua.o: ../../corelib/element.h ../../corelib/qassembly.h
qtrilinos_arpacklua.o: ../../corelib/qmatrix.h ../../corelib/coordmatrix.h
qtrilinos_arpacklua.o: ../../corelib/coordmatrix_g.h
qtrilinos_arpacklua.o: ../corelib/trilinos_epetra_matrix.h
qtrilinos_arpacklua.o: ../corelib/trilinos_epetra_operator.h
qtrilinos_ml.o: qtrilinos_ml.h
qtrilinos_mllua.o: qtrilinos_ml.h
qtrilinos_ifpack.o: qtrilinos_ifpack.h
qtrilinos_ifpacklua.o: qtrilinos_ifpack.h
qtrilinos_komplexlua.o: ../corelib/trilinos_komplex.h
qtrilinos_komplexlua.o: ../corelib/trilinos_epetra_matrix.h
qtrilinos_komplexlua.o: ../corelib/trilinos_epetra_vector.h
qtrilinos_komplexlua.o: ../../corelib/qcomplex.h
qtrilinos_komplexlua.o: ../corelib/trilinos_epetra_operator.h
qtrilinos_komplexlua.o: ../corelib/trilinos_super_matrix.h
qtrilinos_komplexlua.o: ../corelib/trilinos_indexmap.h
qtrilinos_komplexlua.o: ../corelib/trilinos_super_vector.h
qtrilinos_meshlua.o: ../corelib/trilinos_mesh.h ../../corelib/mesh.h
qtrilinos_meshlua.o: ../../corelib/qcomplex.h ../../corelib/element.h
qtrilinos_meshlua.o: ../../corelib/qassembly.h ../../corelib/qmatrix.h
qtrilinos_meshlua.o: ../../corelib/coordmatrix.h ../../corelib/cscmatrix.h
qtrilinos_meshlua.o: ../../corelib/coordmatrix_g.h
