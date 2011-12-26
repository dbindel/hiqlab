/*
 * @(#)mwservices.h    generated by: makeheader 4.21  Tue Aug  2 13:24:04 2005
 *
 *		built from:	../../src/include/copyright.h
 *				../../src/include/pragma_interface.h
 *				config/license.cpp
 *				config/prefs.cpp
 *				distcomp/distcomp.cpp
 *				eventqueue/UserEventQueue.cpp
 *				io/capturebuf.cpp
 *				io/cdroot.c
 *				io/console.cpp
 *				io/efopen.cpp
 *				io/filedir.cpp
 *				io/filestat.cpp
 *				io/history_date.cpp
 *				io/iocbk.cpp
 *				io/iofun.cpp
 *				io/isdos.cpp
 *				io/keybrd.cpp
 *				io/pathset.cpp
 *				io/printca.cpp
 *				io/printmat.cpp
 *				io/printnd.cpp
 *				io/prtopaq.cpp
 *				io/prtstruc.cpp
 *				io/strserv.c
 *				io/uctrans.cpp
 *				io/xlate.cpp
 *				lmgr/dongle.c
 *				lmgr/lmgr.cpp
 *				lmgr/lmgrbusyreg.cpp
 *				lmgr/lmisbusyobj.cpp
 *				lmgr/lmisbusyvector.cpp
 *				mat_thread_req/request_queue.cpp
 *				state/feature.cpp
 *				state/javasup.cpp
 *				state/memcounters.cpp
 *				state/mlcharset.cpp
 *				state/mllang.cpp
 *				state/mlmode.cpp
 *				state/mlruntime.cpp
 *				state/mlterm.cpp
 *				state/mlver.cpp
 *				state/stcbk.cpp
 *				svx/complex_scalarw.cpp
 *				svx/complex_vector.cpp
 *				svx/double_scalarw.cpp
 *				svx/double_vector.cpp
 *				svx/int_math.cpp
 *				svx/int_mmx.cpp
 *				svx/log_scalar.cpp
 *				svx/log_vector.cpp
 *				svx/svMath.cpp
 */

#ifndef mwservices_h
#define mwservices_h


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

#ifdef __cplusplus
    extern "C" {
#endif

#ifdef __cplusplus
    }	/* extern "C" */
#endif

#ifdef __cplusplus
    extern "C" {
#endif

#ifdef __cplusplus
    }	/* extern "C" */
#endif


#include "mwutil.h"


#ifdef __cplusplus
extern "C" {
#endif

#define svDoubleScalarPlusW(a,b,c) utDoubleScalarPlus(a,b)
#define svDoubleScalarMinusW(a,b,c) utDoubleScalarMinus(a,b)
#define svDoubleScalarTimesW(a,b,c) utDoubleScalarTimes(a,b)
#define svDoubleScalarUminusW(a,b) utDoubleScalarUminus(a)
#define svDoubleScalarUplusW(a,b) utDoubleScalarUplus(a)
#define svDoubleScalarAnd(a,b) svDoubleScalarAndW(a,b)
#define svDoubleScalarOr(a,b) svDoubleScalarOrW(a,b)
#define svDoubleScalarXor(a,b) svDoubleScalarXorW(a,b)
#define svDoubleScalarNot(a) svDoubleScalarNotW(a)
#define svDoubleScalarLtW(a,b,c) utDoubleScalarLt(a,b)
#define svDoubleScalarLeW(a,b,c) utDoubleScalarLe(a,b)
#define svDoubleScalarGtW(a,b,c) utDoubleScalarGt(a,b)
#define svDoubleScalarGeW(a,b,c) utDoubleScalarGe(a,b)
#define svDoubleScalarEqW(a,b,c) utDoubleScalarEq(a,b)
#define svDoubleScalarNeW(a,b,c) utDoubleScalarNe(a,b)
#define svDoubleScalarNe(a,b) svDoubleScalarNeW(a,b,NULL)
#define svDoubleScalarEq(a,b) svDoubleScalarEqW(a,b,NULL)
#define svDoubleScalarGe(a,b) svDoubleScalarGeW(a,b,NULL)
#define svDoubleScalarGt(a,b) svDoubleScalarGtW(a,b,NULL)
#define svDoubleScalarLe(a,b) svDoubleScalarLeW(a,b,NULL)
#define svDoubleScalarLt(a,b) svDoubleScalarLtW(a,b,NULL)
#define svDoubleScalarUplus(a) svDoubleScalarUplusW(a,NULL)
#define svDoubleScalarUminus(a) svDoubleScalarUminusW(a,NULL)
#define svDoubleScalarTimes(a,b) svDoubleScalarTimesW(a,b,NULL)
#define svDoubleScalarMinus(a,b) svDoubleScalarMinusW(a,b,NULL)
#define svDoubleScalarPlus(a,b) svDoubleScalarPlusW(a,b,NULL)
#define svDoubleScalarRem(a,b) svDoubleScalarRemW(a,b,NULL)
#define svDoubleScalarLog(a) svDoubleScalarLogW(a,NULL)
#define svDoubleScalarLog10(a) svDoubleScalarLog10W(a,NULL)
#define svDoubleScalarLog2(a) svDoubleScalarLog2W(a,NULL)
#define svDoubleScalarMod(a,b) utDoubleScalarMod(a,b)
#define svDoubleScalarMax(a,b) utDoubleScalarMax(a,b)
#define svDoubleScalarMin(a,b) utDoubleScalarMin(a,b)
#define svDoubleScalarCtranspose(a) utDoubleScalarCtranspose(a)
#define svDoubleScalarTranspose(a) utDoubleScalarTranspose(a)
#define svDoubleScalarAbs(a) utDoubleScalarAbs(a)
#define svDoubleScalarAcos(a) utDoubleScalarAcos(a)
#define svDoubleScalarAcosh(a) utDoubleScalarAcosh(a)
#define svDoubleScalarAcot(a) utDoubleScalarAcot(a)
#define svDoubleScalarAcoth(a) utDoubleScalarAcoth(a)
#define svDoubleScalarAcsc(a) utDoubleScalarAcsc(a)
#define svDoubleScalarAcsch(a) utDoubleScalarAcsch(a)
#define svDoubleScalarAll(a) utDoubleScalarAll(a)
#define svDoubleScalarAny(a) utDoubleScalarAny(a)
#define svDoubleScalarAngle(a) utDoubleScalarAngle(a)
#define svDoubleScalarAsec(a) utDoubleScalarAsec(a)
#define svDoubleScalarAsech(a) utDoubleScalarAsech(a)
#define svDoubleScalarAsin(a) utDoubleScalarAsin(a)
#define svDoubleScalarAsinh(a) utDoubleScalarAsinh(a)
#define svDoubleScalarAtan(a) utDoubleScalarAtan(a)
#define svDoubleScalarAtan2(a,b) utDoubleScalarAtan2(a,b)
#define svDoubleScalarAtanh(a) utDoubleScalarAtanh(a)
#define svDoubleScalarCeil(a) utDoubleScalarCeil(a)
#define svDoubleScalarCos(a) utDoubleScalarCos(a)
#define svDoubleScalarCosh(a) utDoubleScalarCosh(a)
#define svDoubleScalarDot(a) utDoubleScalarDot(a)
#define svDoubleScalarDouble(a) utDoubleScalarDouble(a)
#define svDoubleScalarExp(a) utDoubleScalarExp(a)
#define svDoubleScalarFloor(a) utDoubleScalarFloor(a)
#define svDoubleScalarFull(a) utDoubleScalarFull(a)
#define svDoubleScalarInt8(a) utDoubleScalarInt8(a)
#define svDoubleScalarInt16(a) utDoubleScalarInt16(a)
#define svDoubleScalarInt32(a) utDoubleScalarInt32(a)
#define svDoubleScalarIsempty(a) utDoubleScalarIsempty(a)
#define svDoubleScalarIsfinite(a) utDoubleScalarIsfinite(a)
#define svDoubleScalarIslogical(a) utDoubleScalarIslogical(a)
#define svDoubleScalarIsmember(a,b) utDoubleScalarIsmember(a,b)
#define svDoubleScalarIsreal(a) utDoubleScalarIsreal(a)
#define svDoubleScalarIssparse(a) utDoubleScalarIssparse(a)
#define svDoubleScalarLength(a) utDoubleScalarLength(a)
#define svDoubleScalarNdims(a) utDoubleScalarNdims(a)
#define svDoubleScalarPower(a,b) utDoubleScalarPower(a,b)
#define svDoubleScalarPow2(a) utDoubleScalarPow2(a)
#define svDoubleScalarReal(a) utDoubleScalarReal(a)
#define svDoubleScalarRound(a) utDoubleScalarRound(a)
#define svDoubleScalarSec(a) utDoubleScalarSec(a)
#define svDoubleScalarSech(a) utDoubleScalarSech(a)
#define svDoubleScalarSign(a) utDoubleScalarSign(a)
#define svDoubleScalarSin(a) utDoubleScalarSin(a)
#define svDoubleScalarSingle(a) utDoubleScalarSingle(a)
#define svDoubleScalarSinh(a) utDoubleScalarSinh(a)
#define svDoubleScalarSqrt(a) utDoubleScalarSqrt(a)
#define svDoubleScalarTan(a) utDoubleScalarTan(a)
#define svDoubleScalarTanh(a) utDoubleScalarTanh(a)
#define svDoubleScalarUint8(a) utDoubleScalarUint8(a)
#define svDoubleScalarUint16(a) utDoubleScalarUint16(a)
#define svDoubleScalarUint32(a) utDoubleScalarUint32(a)

#ifdef __cplusplus 
}
#endif
/*
 * NOTE: The accelerator is relying on the fact that the following logical operators are
 *       implemented as macros for performance.
 */
/* a | b */
#define svDoubleScalarOrW(a,b) (utIsNaN(a) || utIsNaN(b)        \
         ? (svNaNConversionError(),                             \
             utGetNaN())                                        \
         : utDoubleScalarOr(a,b))
							
/* ~a */
#define svDoubleScalarNotW(a) (utIsNaN(a)       \
         ? (svNaNConversionError(),             \
            utGetNaN())                         \
         : utDoubleScalarNot(a))

/* a & b */
#define svDoubleScalarAndW(a,b) (utIsNaN(a) || utIsNaN(b)       \
         ? (svNaNConversionError(),                             \
            utGetNaN())                                         \
         : utDoubleScalarAnd(a,b))

/* xor(a,b) */	
#define svDoubleScalarXorW(a,b) ( utIsNaN(b) || utIsNaN(a)      \
        ? (svNaNConversionError(),                              \
           utGetNaN())                                          \
        : utDoubleScalarXor(a,b))
													   
/* cot(a) */									
#define svDoubleScalarCotW(a,b)  ( utEQZero(a)                          \
         ? ( !*b                                                        \
                ? (svDivideByZeroWarning(), *b = true, utGetInf())      \
                : utDoubleScalarCot(a))                                 \
         : utDoubleScalarCot(a))

#define svDoubleScalarCot(a) ( utEQZero(a)              \
         ? (svDivideByZeroWarning(), utGetInf())        \
         : utDoubleScalarCot(a))	

/* coth(a) */	
#define svDoubleScalarCothW(a,b)  ( utEQZero(a)                         \
         ? ( !*b                                                        \
                ? (svDivideByZeroWarning(), *b = true, utGetInf())      \
                : utDoubleScalarCoth(a))                                \
         : utDoubleScalarCoth(a))	

#define svDoubleScalarCoth(a) ( utEQZero(a)             \
         ? (svDivideByZeroWarning(), utGetInf())        \
         : utDoubleScalarCoth(a))		

/* csc(a) */									
#define svDoubleScalarCscW(a,b)  ( utEQZero(a)                          \
         ? ( !*b                                                        \
                ? (svDivideByZeroWarning(), *b = true, utGetInf())      \
                : utDoubleScalarCsc(a))                                 \
         : utDoubleScalarCsc(a))	

#define svDoubleScalarCsc(a)  ( utEQZero(a)             \
         ? (svDivideByZeroWarning(), utGetInf())        \
         : utDoubleScalarCsc(a))		

/* csch(a) */
#define svDoubleScalarCschW(a,b)  ( utEQZero(a)                         \
         ? ( !*b                                                        \
                ? (svDivideByZeroWarning(), *b = true, utGetInf())      \
                : utDoubleScalarCsch(a))                                \
         : utDoubleScalarCsch(a))	

#define svDoubleScalarCsch(a)  ( utEQZero(a)            \
         ? (svDivideByZeroWarning(), utGetInf())        \
         : utDoubleScalarCsch(a))
	

#define svDivZero(x) (utGTZero(x) ? utGetInf() :                \
                        utLTZero(x) ? -utGetInf() : utGetNaN())

/*
 * NOTE: The accelerator is relying on the fact that rdivide and ldivide are
 *       implemented as macros for performance.
 */

/* a/b and a./b */

#define svDoubleScalarRdivide_Inline(a, b)                              \
    ( utNEZero(b)                                                       \
        ? (a)/(b)                                                       \
        : (svDivideByZeroWarning() , (a)/(b)))

#define svDoubleScalarRdivideW_Inline(a, b, c)                                  \
    ( utNEZero(b)                                                               \
        ? (a)/(b)                                                               \
        : ( !*c ? (svDivideByZeroWarning() , *c = true, (a)/(b))                            \
        : (a)/(b)))


/* a\b and a.\b */  
#define svDoubleScalarLdivide_Inline(a, b) svDoubleScalarRdivide_Inline(b, a)
#define svDoubleScalarLdivideW_Inline(a, b, c) svDoubleScalarRdivideW_Inline(b, a, c)

#ifdef __cplusplus
extern "C" {
#endif
extern void svDivideByZeroWarning(void);
#ifdef __cplusplus
}
#endif

#ifdef __cplusplus
extern "C" {
#endif
extern void svNaNConversionError(void);
#ifdef __cplusplus
}
#endif

#ifdef __cplusplus
extern "C" {
#endif
extern double svDoubleScalarRdivideW(double a, double b, bool *c);
#ifdef __cplusplus
}
#endif

#ifdef __cplusplus
extern "C" {
#endif
extern double svDoubleScalarRdivide(double a, double b);
#ifdef __cplusplus
}
#endif

#ifdef __cplusplus
extern "C" {
#endif
extern double svDoubleScalarLdivideW(double a, double b, bool *c);
#ifdef __cplusplus
}
#endif

#ifdef __cplusplus
extern "C" {
#endif
extern double svDoubleScalarLdivide(double a, double b);
#ifdef __cplusplus
}
#endif

#ifdef __cplusplus
extern "C" {
#endif
/* remW(a,b,c) */
/*
	Remainder after division function
        c is a boolean pointer that activates (false) or 
        not the zero warning
*/
extern double svDoubleScalarRemW( double a, double b, bool *c);
#ifdef __cplusplus
}
#endif

#ifdef __cplusplus
extern "C" {
#endif
/* log(a,b) */
/* Natural logarithm of a 
   b is a boolean pointer that activates (false) or not the warning
*/
extern double  svDoubleScalarLogW(double a, bool *b);
#ifdef __cplusplus
}
#endif

#ifdef __cplusplus
extern "C" {
#endif
/* log2(a,b) */
/* Base 2 logarithm of a 
   b is a boolean pointer that activates (false) or not the warning
*/
extern double  svDoubleScalarLog2W(double a, bool *b);
#ifdef __cplusplus
}
#endif

#ifdef __cplusplus
extern "C" {
#endif
/* log10(a,b) */
/* Common (base 10) logarithm of a 
   b is a boolean pointer that activates (false) or not the warning
*/
extern double  svDoubleScalarLog10W(double a, bool *b);
#ifdef __cplusplus
}
#endif

#endif /* mwservices_h */
