include ../../make.inc

all: dsbweb dsbweb.pdf

dsbweb: dsbweb.c
	$(CC) -o dsbweb dsbweb.c

dsbweb.pdf: dsbweb.c
	./dsbweb dsbweb.c -d dsbweb.tex -p ""
	pdflatex dsbweb.tex

clean:
	rm -rf *.aux *.log *.pdf dsbweb dsbweb.tex
