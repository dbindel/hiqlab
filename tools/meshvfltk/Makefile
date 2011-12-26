include ../../make.inc

meshv: meshv.o
	$(CXX) -o meshv meshv.o $(FLTK_LIBS) $(GD_LIBS)
	$(FLTK_CONFIG) --post meshv

meshv.o: meshv.cc
	$(CXX) -g -c $(FLTK_CXXFLAGS) $(GD_CFLAGS) meshv.cc

clean:
	rm -f *.o meshv *~
