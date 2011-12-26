include ../make.inc

all: testx
	./testx

testx: testx.cc
	$(CXX) -o testx testx.cc

testx.cc: test-in.cc ../matexpr
	../matexpr test-in.cc > testx.cc

../matexpr:
	(cd ..; make)

clean:
	rm -f testx.cc testx *.o *~
