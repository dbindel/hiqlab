include ../make.inc

# === Test targets ===

test1x: test1-in.cc
	../matexpr test1-in.cc > test1.cc
	$(CXX) -o test1x test1.cc

# === Clean-up targets ===

clean:
	rm -f *.o *~
	rm -f test1.cc test1x

