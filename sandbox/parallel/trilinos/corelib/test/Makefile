include ../../../make.inc

LIBQLAB_TL= $(LIBDIR)/libqlab_trilinos.a

OBJS= ../AnasaziEpetraMultiVectorComplexAdapter.o \
      ../AnasaziEpetraOperatorComplexAdapter.o

all: BKSComplexTester7 MVOPTesterEx7_complex MVOPTesterEx7_complex2

.cc.o: 
	$(CXX) -DHAVE_CONFIG_H $(CXXFLAGS) -I../../corelib -I. -I.. -c $<

MVOPTesterEx7_complex.o: MVOPTesterEx7_complex.cc
MVOPTesterEx7_complex2.o: MVOPTesterEx7_complex2.cc
BKSComplexTester7.o: BKSComplexTester7.cc

BKSComplexTester7: BKSComplexTester7.o 
	$(CXX) -DHAVE_CONFIG_H -o BKSComplexTester7 BKSComplexTester7.o \
                $(OBJS) $(LIBQLAB_TL) $(LIBQLAB)\
		$(READLINE_LIBS) $(FLTK_LIBS) $(LIBS)

MVOPTesterEx7_complex: MVOPTesterEx7_complex.o 
	$(CXX) -DHAVE_CONFIG_H -o MVOPTesterEx7_complex MVOPTesterEx7_complex.o \
                $(OBJS) $(LIBQLAB_TL) $(LIBQLAB)\
		$(READLINE_LIBS) $(FLTK_LIBS) $(LIBS)

MVOPTesterEx7_complex2: MVOPTesterEx7_complex2.o 
	$(CXX) -DHAVE_CONFIG_H -o MVOPTesterEx7_complex2 MVOPTesterEx7_complex2.o \
                $(OBJS) $(LIBQLAB_TL) $(LIBQLAB)\
		$(READLINE_LIBS) $(FLTK_LIBS) $(LIBS)

clean:
	rm -f *lua.* *.o *~ MVOPTesterEx7_complex MVOPTesterEx7_complex2 BKSComplexTester7
