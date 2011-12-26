# HiQLab
# Copyright (c): Regents of the University of California
# $Id: Makefile,v 1.27 2006/08/26 18:45:47 dbindel Exp $

include make.inc

DATE=`/bin/date +%F`
VER=$(DATE)
SRC=\
	Makefile		\
	init.m			\
	init.lua		\
	MANIFEST		\
	init.lua.in		\
	config.h.in		\
	hiqlab.h.in		\
	make.inc.in		\
	configure.in		\
	aclocal.m4		\
	configure		\
	config/config.sub	\
	config/config.guess	\
	config/depcomp		\
	config/install-sh	\
	config/ltmain.sh	\
	config/missing		\
	config/mkinstalldirs	\
	`find models -name "*.lua" -print` \
	`find models -name "*.m"   -print` \
	src/Makefile		\
	src/corelib/Makefile 	\
	src/corelib/*.h		\
	src/corelib/*.cc 	\
	src/element/Makefile 	\
	src/element/*.h		\
	src/element/*.cc 	\
	src/doc/Makefile	\
	src/doc/*.tex		\
	src/doc/*.pdf		\
	src/\@dxfile/*.m 	\
	src/lua/Makefile	\
	src/lua/*.h		\
	src/lua/*.cc		\
	src/lua/*.pkg		\
	src/lua/*.lua		\
	src/mex/Makefile	\
	src/mex/*.h		\
	src/mex/*.cc		\
	src/mex/*.pkg		\
	src/mfiles/*.m		\
	test/*.m		\
	test/*.lua

GEN= \
	src/mex/*.mex*	\
	src/mex/*.dll	\
	src/mex/*.m	\
	src/lua/hiqlab.*

all: serial

serial: qtools
	(cd src;   make serial)

mex: qtools serial
	(cd src;   make mex)

oct: qtools serial
	(cd src;   make oct)

trilinos: qtoolsp
	(cd src;   make trilinos)

petsc: qtoolsp
	(cd src;   make petsc)

hybrid: qtoolsp
	(cd src;   make hybrid)

qtools:
	(cd tools; make)

qtoolsp: qtools
	(cd tools; make allp)

tests:
	(cd test_cluster; make test_trilinos; make test_petsc;)

doc:
	(cd src/doc; make; make clean)

tgz: clean
	( ls $(SRC) | sed s:^:hiqlab-$(VER)/: > MANIFEST ; \
	  cat tools/MANIFEST | sed s:^:hiqlab-$(VER)/tools/: > MANIFEST ; \
	  cd .. ; \
	  ln -s hiqlab hiqlab-$(VER) ; \
	  tar -h -czf hiqlab-$(VER).tar.gz --files-from hiqlab/MANIFEST ; \
	  rm hiqlab-$(VER) ; \
	  mv hiqlab-$(VER).tar.gz hiqlab )

bin-tgz: 
	( ls $(SRC) | sed s:^:hiqlab-$(VER)/: > MANIFEST ; \
	  cat tools/MANIFEST | sed s:^:hiqlab-$(VER)/tools/: > MANIFEST ; \
	  ls $(GEN) | sed s:^:hiqlab-$(VER)/: >> MANIFEST ; \
	  cd .. ; \
	  ln -s hiqlab hiqlab-$(VER) ; \
	  tar -czf hiqlab-$(VER)-bin.tar.gz --files-from hiqlab/MANIFEST ; \
	  rm hiqlab-$(VER) ; \
	  mv hiqlab-$(VER)-bin.tar.gz hiqlab )

tag:
	( cvs tag rel-$(VER) . )

clean:
	(cd src;     make clean)
	(cd src/doc; make clean)
	(cd tools;   make clean)

cleanp:
	(cd src;     make cleanp)
	(cd tools;   make cleanp)
