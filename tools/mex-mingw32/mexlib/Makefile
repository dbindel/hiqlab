include ../make.inc

mexlibs:
	$(CC) -c fixup.c
	$(DLLTOOL) --def $(WMEXINC)/libmx.def  --output-lib mexlib1.lib
	$(DLLTOOL) --def $(WMEXINC)/libmex.def --output-lib mexlib2.lib
	$(DLLTOOL) --def $(WMEXINC)/libmat.def --output-lib mexlib3.lib

clean:
	rm -f fixup.o *.lib
