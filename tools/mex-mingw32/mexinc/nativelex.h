/*
 * @(#)nativelex.h    generated by: makeheader 4.21  Tue Jul 26 23:14:04 2005
 *
 *		built from:	../../../src/include/copyright.h
 *				../../../src/include/pragma_interface.h
 *				./mtok.c
 *				./nativelex.cpp
 */

#ifndef nativelex_h
#define nativelex_h


/*
 * Copyright 1984-2003 The MathWorks, Inc.
 * All Rights Reserved.
 */



/* Copyright 2003-2004 The MathWorks, Inc. */

/*
 * Prevent g++ from making copies of vtable and typeinfo data
 * in every compilation unit.  By allowing for only one, we can
 * save space and prevent some situations where the linker fails
 * to coalesce them properly into a single entry.
 *
 * References:
 *    http://gcc.gnu.org/onlinedocs/gcc/Vague-Linkage.html#Vague%20Linkage
 *    http://gcc.gnu.org/onlinedocs/gcc/C---Interface.html
 */

#ifdef __cplusplus
#  ifdef __linux__
#    pragma interface
#  endif
#endif


#ifdef __cplusplus
    extern "C" {
#endif

#ifdef __cplusplus
    }	/* extern "C" */
#endif

#endif /* nativelex_h */