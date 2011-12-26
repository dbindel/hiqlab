dnl =====================================================================
dnl Available from the GNU Autoconf Macro Archive at:
dnl http://www.gnu.org/software/ac-archive/htmldoc/acx_blas.html
dnl
AC_DEFUN([ACX_BLAS], [
AC_PREREQ(2.50)
AC_REQUIRE([AC_F77_LIBRARY_LDFLAGS])
acx_blas_ok=no

AC_ARG_WITH(blas,
	[AC_HELP_STRING([--with-blas=<lib>], [use BLAS library <lib>])])
case $with_blas in
	yes | "") ;;
	no) acx_blas_ok=disable ;;
	-* | */* | *.a | *.so | *.so.* | *.o) BLAS_LIBS="$with_blas" ;;
	*) BLAS_LIBS="-l$with_blas" ;;
esac

# Get fortran linker names of BLAS functions to check for.
AC_F77_FUNC(sgemm)
AC_F77_FUNC(dgemm)
AC_F77_FUNC(xerbla)

acx_blas_save_LIBS="$LIBS"
LIBS="$LIBS $FLIBS"

# First, check BLAS_LIBS environment variable
if test $acx_blas_ok = no; then
if test "x$BLAS_LIBS" != x; then
	save_LIBS="$LIBS"; LIBS="$BLAS_LIBS $LIBS"
	AC_MSG_CHECKING([for $sgemm in $BLAS_LIBS])
	AC_TRY_LINK_FUNC($sgemm, [acx_blas_ok=yes], [BLAS_LIBS=""])
	AC_MSG_RESULT($acx_blas_ok)
	LIBS="$save_LIBS"
fi
fi

# BLAS linked to by default?  (happens on some supercomputers)
if test $acx_blas_ok = no; then
	save_LIBS="$LIBS"; LIBS="$LIBS"
	AC_CHECK_FUNC($sgemm, [acx_blas_ok=yes])
	LIBS="$save_LIBS"
fi

# BLAS in libgoto library? 
if test $acx_blas_ok = no; then
	AC_CHECK_LIB(xerbla, $xerbla,
		[AC_CHECK_LIB(goto, $dgemm,
			[acx_blas_ok=yes
			 BLAS_LIBS="-lgoto -lxerbla"],
			[], [-lxerbla])])
fi

# BLAS in ATLAS library? (http://math-atlas.sourceforge.net/)
if test $acx_blas_ok = no; then
	AC_CHECK_LIB(atlas, ATL_xerbla,
		[AC_CHECK_LIB(f77blas, $sgemm,
		[AC_CHECK_LIB(cblas, cblas_dgemm,
			[acx_blas_ok=yes
			 BLAS_LIBS="-lcblas -lf77blas -latlas"],
			[], [-lf77blas -latlas])],
			[], [-latlas])])
fi

# BLAS in PhiPACK libraries? (requires generic BLAS lib, too)
if test $acx_blas_ok = no; then
	AC_CHECK_LIB(blas, $sgemm,
		[AC_CHECK_LIB(dgemm, $dgemm,
		[AC_CHECK_LIB(sgemm, $sgemm,
			[acx_blas_ok=yes; BLAS_LIBS="-lsgemm -ldgemm -lblas"],
			[], [-lblas])],
			[], [-lblas])])
fi

# BLAS in Alpha CXML library?
if test $acx_blas_ok = no; then
	AC_CHECK_LIB(cxml, $sgemm, [acx_blas_ok=yes;BLAS_LIBS="-lcxml"])
fi

# BLAS in Alpha DXML library? (now called CXML, see above)
if test $acx_blas_ok = no; then
	AC_CHECK_LIB(dxml, $sgemm, [acx_blas_ok=yes;BLAS_LIBS="-ldxml"])
fi

# BLAS in Sun Performance library?
if test $acx_blas_ok = no; then
	if test "x$GCC" != xyes; then # only works with Sun CC
		AC_CHECK_LIB(sunmath, acosp,
			[AC_CHECK_LIB(sunperf, $sgemm,
        			[BLAS_LIBS="-xlic_lib=sunperf -lsunmath"
                                 acx_blas_ok=yes],[],[-lsunmath])])
	fi
fi

# BLAS in SCSL library?  (SGI/Cray Scientific Library)
if test $acx_blas_ok = no; then
	AC_CHECK_LIB(scs, $sgemm, [acx_blas_ok=yes; BLAS_LIBS="-lscs"])
fi

# BLAS in SGIMATH library?
if test $acx_blas_ok = no; then
	AC_CHECK_LIB(complib.sgimath, $sgemm,
		     [acx_blas_ok=yes; BLAS_LIBS="-lcomplib.sgimath"])
fi

# BLAS in IBM ESSL library? (requires generic BLAS lib, too)
if test $acx_blas_ok = no; then
	AC_CHECK_LIB(blas, $sgemm,
		[AC_CHECK_LIB(essl, $sgemm,
			[acx_blas_ok=yes; BLAS_LIBS="-lessl -lblas"],
			[], [-lblas $FLIBS])])
fi

# Generic BLAS library?
if test $acx_blas_ok = no; then
	AC_CHECK_LIB(blas, $sgemm, [acx_blas_ok=yes; BLAS_LIBS="-lblas"])
fi

AC_SUBST(BLAS_LIBS)

LIBS="$acx_blas_save_LIBS"

# Finally, execute ACTION-IF-FOUND/ACTION-IF-NOT-FOUND:
if test x"$acx_blas_ok" = xyes; then
        ifelse([$1],,AC_DEFINE(HAVE_BLAS,1,[Define if you have a BLAS library.]),[$1])
        :
else
        acx_blas_ok=no
        $2
fi
])dnl ACX_BLAS


dnl =====================================================================
dnl Available from the GNU Autoconf Macro Archive at:
dnl http://www.gnu.org/software/ac-archive/htmldoc/acx_lapack.html
dnl
AC_DEFUN([ACX_LAPACK], [
AC_REQUIRE([ACX_BLAS])
acx_lapack_ok=no

AC_ARG_WITH(lapack,
        [AC_HELP_STRING([--with-lapack=<lib>], [use LAPACK library <lib>])])
case $with_lapack in
        yes | "") ;;
        no) acx_lapack_ok=disable ;;
        -* | */* | *.a | *.so | *.so.* | *.o) LAPACK_LIBS="$with_lapack" ;;
        *) LAPACK_LIBS="-l$with_lapack" ;;
esac

# Get fortran linker name of LAPACK function to check for.
AC_F77_FUNC(cheev)

# We cannot use LAPACK if BLAS is not found
if test "x$acx_blas_ok" != xyes; then
        acx_lapack_ok=noblas
fi

# First, check LAPACK_LIBS environment variable
if test "x$LAPACK_LIBS" != x; then
        save_LIBS="$LIBS"; LIBS="$LAPACK_LIBS $BLAS_LIBS $LIBS $FLIBS"
        AC_MSG_CHECKING([for $cheev in $LAPACK_LIBS])
        AC_TRY_LINK_FUNC($cheev, [acx_lapack_ok=yes], [LAPACK_LIBS=""])
        AC_MSG_RESULT($acx_lapack_ok)
        LIBS="$save_LIBS"
        if test acx_lapack_ok = no; then
                LAPACK_LIBS=""
        fi
fi

# LAPACK linked to by default?  (is sometimes included in BLAS lib)
if test $acx_lapack_ok = no; then
        save_LIBS="$LIBS"; LIBS="$LIBS $BLAS_LIBS $FLIBS"
        AC_CHECK_FUNC($cheev, [acx_lapack_ok=yes])
        LIBS="$save_LIBS"
fi

# Generic LAPACK library?
for lapack in lapack lapack_rs6k; do
        if test $acx_lapack_ok = no; then
                save_LIBS="$LIBS"; LIBS="$BLAS_LIBS $LIBS"
                AC_CHECK_LIB($lapack, $cheev,
                    [acx_lapack_ok=yes; LAPACK_LIBS="-l$lapack"], [], [$FLIBS])
                LIBS="$save_LIBS"
        fi
done

AC_SUBST(LAPACK_LIBS)

# Finally, execute ACTION-IF-FOUND/ACTION-IF-NOT-FOUND:
if test x"$acx_lapack_ok" = xyes; then
        ifelse([$1],,AC_DEFINE(HAVE_LAPACK,1,[Define if you have LAPACK library.]),[$1])
        :
else
        acx_lapack_ok=no
        $2
fi
])dnl ACX_LAPACK


dnl =====================================================================
dnl Available from the GNU Autoconf Macro Archive at:
dnl http://www.gnu.org/software/ac-archive/htmldoc/vl_lib_readline.html
dnl
AC_DEFUN([ACX_LIB_READLINE], [
  AC_CACHE_CHECK([for a readline compatible library],
                 vl_cv_lib_readline, [
    save_LIBS="$LIBS"
    for readline_lib in readline edit editline; do
      for termcap_lib in "" termcap curses ncurses; do
        if test -z "$termcap_lib"; then
          READLINE_LIBS="-l$readline_lib"
        else
          READLINE_LIBS="-l$readline_lib -l$termcap_lib"
        fi
        LIBS="$save_LIBS $READLINE_LIBS"
        AC_TRY_LINK_FUNC(readline, vl_cv_lib_readline="$READLINE_LIBS")
        LIBS="$save_LIBS"
        if test -n "$vl_cv_lib_readline"; then
          break
        fi
      done
      if test -n "$vl_cv_lib_readline"; then
        break
      fi
    done
    if test -z "$vl_cv_lib_readline"; then
      vl_cv_lib_readline="no"
      READLINE_LIBS=""
    fi
  ])

  if test "$vl_cv_lib_readline" != "no"; then
    AC_DEFINE(HAVE_READLINE, 1,
              [Define if you have a readline compatible library])
    AC_CHECK_HEADERS(readline.h readline/readline.h)
    AC_CACHE_CHECK([whether readline supports history],
                   vl_cv_lib_readline_history, [
      vl_cv_lib_readline_history="no"
      save_LIBS="$LIBS"
      LIBS="$save_LIBS -lhistory"
      AC_TRY_LINK_FUNC(add_history, vl_cv_lib_readline_history="yes")
      LIBS="$save_LIBS"
    ])
    if test "$vl_cv_lib_readline_history" = "yes"; then
      AC_DEFINE(HAVE_READLINE_HISTORY, 1,
                [Define if your readline library has \`add_history'])
      AC_CHECK_HEADERS(history.h readline/history.h)
      save_READLINE_LIBS="$READLINE_LIBS"
      READLINE_LIBS="$save_READLINE_LIBS -lhistory"
    fi
  fi

  AC_SUBST(READLINE_LIBS)
])dnl


dnl =====================================================================
dnl Optional ARPACK link
dnl
AC_DEFUN([AC_WITH_ARPACK],
 [ac_have_arpack="no"
  AC_ARG_WITH(arpack,  
   AC_HELP_STRING([--with-arpack=ARPACK_LIB], 
                  [location of ARPACK installation dir]),
   [ac_have_arpack="yes"; ARPACK_LIB="$with_arpack"],
   [ac_have_arpack="no" ; ARPACK_LIB=""])

 if test x"$ac_have_arpack" != x"no"; then
   ifelse([$1],,
          AC_DEFINE(HAVE_ARPACK, 1,
                    [Define if you have ARPACK]),
          [$1])
   :
 else
   ac_have_arpack=no
   $2
 fi

 AC_SUBST(ARPACK_LIB)
])dnl AC_WITH_ARPACK


dnl =====================================================================
dnl Optional SuperLU link taken from Netsolve 1.4.1 configuration
dnl
AC_DEFUN([AC_WITH_SUPERLU],
 [ac_have_superlu="no"
  AC_ARG_WITH(superlu,  
   AC_HELP_STRING([--with-superlu=DIR], 
                  [location of SuperLU installation dir]),
   [ac_have_superlu="yes";
    SUPERLU_DIR="$with_superlu";
    SUPERLU_INCLUDE="-I$with_superlu/SRC";
    SUPERLU_LIB="$with_superlu/libsuperlu.a"],
   [ac_have_superlu="no";
    SUPERLU_DIR="";
    SUPERLU_INCLUDE="";
    SUPERLU_LIB=""])

  AC_ARG_WITH(superluinclude,
   AC_HELP_STRING([--with-superluinclude=DIR],
                  [SuperLU include directory]),
   [SUPERLU_INCLUDE="-I$with_superluinclude"])

  AC_ARG_WITH(superlulib,
             AC_HELP_STRING([--with-superlulib=<lib>],
                            [SuperLU link line]),
             [SUPERLU_LIB="$with_superlulib"])

 if test x"$ac_have_superlu" != x"no"; then
   ifelse([$1],,
          AC_DEFINE(HAVE_HIQLAB_SUPERLU, 1,
                    [Define if you have SuperLU]),
          [$1])
   :
 else
   ac_have_superlu=no
   $2
 fi

 AC_SUBST(SUPERLU_DIR)
 AC_SUBST(SUPERLU_INCLUDE)
 AC_SUBST(SUPERLU_LIB)
])dnl AC_WITH_SUPERLU


dnl =====================================================================
dnl Optional SuperLU_Dist link copied from SuperLU Netsolve 1.4.1 configuration
dnl
AC_DEFUN([AC_WITH_SUPERLUDIST],
 [ac_have_superludist="no"
  AC_ARG_WITH(superludist,  
   AC_HELP_STRING([--with-superludist=DIR], 
                  [location of SuperLU_Dist installation dir]),
   [ac_have_superludist="yes";
    SUPERLUDIST_DIR="$with_superludist";
    SUPERLUDIST_INCLUDE="-I$with_superludist/SRC";
    SUPERLUDIST_LIB="$with_superludist/libsuperlu_dist.a";
    HAVE_HIQLAB_SUPERLUDIST_TRUE=""],
   [ac_have_superludist="no";
    SUPERLUDIST_DIR="";
    SUPERLUDIST_INCLUDE="";
    SUPERLUDIST_LIB="";
    HAVE_HIQLAB_SUPERLUDIST_TRUE="#"])

  AC_ARG_WITH(superludistinclude,
   AC_HELP_STRING([--with-superludistinclude=DIR],
                  [SuperLU Dist include directory]),
   [SUPERLUDIST_LIB="-I$with_superludistinclude"])

  AC_ARG_WITH(superludistlib,
             AC_HELP_STRING([--with-superludistlib=<lib>],
                            [SuperLU_Dist link line]),
             [SUPERLUDIST_LIB="$with_superludistlib"])

 if test x"$ac_have_superludist" != x"no"; then
   ifelse([$1],,
          AC_DEFINE(HAVE_HIQLAB_SUPERLUDIST, 1,
                    [Define if you have SuperLU Dist]),
          [$1])
   :
 else
   ac_have_superludist=no
   $2
 fi

 AC_SUBST(SUPERLUDIST_DIR)
 AC_SUBST(SUPERLUDIST_INCLUDE)
 AC_SUBST(SUPERLUDIST_LIB)
 AC_SUBST(HAVE_HIQLAB_SUPERLUDIST_TRUE)
])dnl AC_WITH_SUPERLUDIST


dnl =====================================================================
dnl Optional BLACS link copied from SuperLU Netsolve 1.4.1 configuration
dnl
AC_DEFUN([AC_WITH_BLACS],
 [ac_have_blacs="no"
  AC_ARG_WITH(blacs,  
   AC_HELP_STRING([--with-blacs=<lib>], 
                  [use BLACS library <lib>]),
   [ac_have_blacs="yes";
    BLACS_LIB="$with_blacs"],
   [ac_have_blacs="no";
    BLACS_LIB=""])

 if test x"$ac_have_blacs" != x"no"; then
   ifelse([$1],,
          AC_DEFINE(HAVE_BLACS, 1,
                    [Define if you have BLACS]),
          [$1])
   :
 else
   ac_have_blacs=no
   $2
 fi

 AC_SUBST(BLACS_LIB)
])dnl AC_WITH_BLACS


dnl =====================================================================
dnl Optional Scalapack link copied from SuperLU Netsolve 1.4.1 configuration
dnl
AC_DEFUN([AC_WITH_SCALAPACK],
 [ac_have_scalapack="no"
  AC_ARG_WITH(scalapack,  
   AC_HELP_STRING([--with-scalapack=<lib>], 
                  [use SCALAPACK library <lib> (Requires BLACS)]),
   [ac_have_scalapack="yes";
    SCALAPACK_LIB="$with_scalapack"],
   [ac_have_scalapack="no";
    SCALAPACK_LIB=""])

 if test x"$ac_have_scalapack" != x"no"; then
   ifelse([$1],,
          AC_DEFINE(HAVE_HIQLAB_SCALAPACK, 1,
                    [Define if you have SCALAPACK]),
          [$1])
   :
 else
   ac_have_scalapack=no
   $2
 fi

 AC_SUBST(SCALAPACK_LIB)
])dnl AC_WITH_SCALAPACK


dnl =====================================================================
dnl Optional FORTRAN90 link copied from SuperLU Netsolve 1.4.1 configuration
dnl FIXME: Would like to replace if can find AC_WITH_F90. Same goes for 
dnl following library statement.
dnl 
AC_DEFUN([AC_WITH_F90_HIQLAB],
 [ac_have_f90hiqlab="no"
  AC_ARG_WITH(f90hiqlab,  
   AC_HELP_STRING([--with-f90hiqlab=F90], 
                  [use F90 compiler and linker ]),
   [ac_have_f90hiqlab="yes";
    F90="$with_f90hiqlab"],
   [ac_have_f90hiqlab="no";
    F90=""])

 if test x"$ac_have_f90hiqlab" != x"no"; then
   ifelse([$1],,
          AC_DEFINE(HAVE_F90_HIQLAB, 1,
                    [Define if you have F90 compiler and linker]),
          [$1])
   :
 else
   ac_have_f90hiqlab=no
   $2
 fi

 AC_SUBST(F90)
])dnl AC_WITH_F90_HIQLAB


dnl =====================================================================
dnl Optional FORTRAN90 link copied from SuperLU Netsolve 1.4.1 configuration
dnl
AC_DEFUN([AC_WITH_F90LIBS_HIQLAB],
 [ac_have_f90libshiqlab="no"
  AC_ARG_WITH(f90libshiqlab,  
   AC_HELP_STRING([--with-f90libshiqlab=<lib>], 
                  [use F90 library <lib>]),
   [ac_have_f90libshiqlab="yes";
    F90_LIBS="$with_f90libshiqlab"],
   [ac_have_f90libshiqlab="no";
    F90_LIBS=""])

 if test x"$ac_have_f90libshiqlab" != x"no"; then
   ifelse([$1],,
          AC_DEFINE(HAVE_F90LIBS_HIQLAB, 1,
                    [Define if you have F90LIBS]),
          [$1])
   :
 else
   ac_have_f90libshiqlab=no
   $2
 fi

 AC_SUBST(F90_LIBS)
])dnl AC_WITH_F90LIBS_HIQLAB


dnl =====================================================================
dnl Optional Mumps link copied from SuperLU Netsolve 1.4.1 configuration
dnl
AC_DEFUN([AC_WITH_MUMPS],
 [ac_have_mumps="no"
  AC_ARG_WITH(mumps,  
   AC_HELP_STRING([--with-mumps=DIR], 
                  [location of Mumps(double,complexdouble) and Pord installation directory(Requires BLACS and Scalapack)]),
   [ac_have_mumps="yes";
    MUMPS_DIR="$with_mumps";
    MUMPS_INCLUDE="-I$with_mumps/include";
    MUMPS_LIB="$with_mumps/lib/libdmumps.a $with_mumps/lib/libzmumps.a $with_mumps/lib/libpord.a";
    HAVE_HIQLAB_MUMPS_TRUE=""],
   [ac_have_mumps="no";
    MUMPS_DIR="";
    MUMPS_INCLUDE="";
    MUMPS_LIB="";
    HAVE_HIQLAB_MUMPS_TRUE="#"])

  AC_ARG_WITH(mumpsinclude,
   AC_HELP_STRING([--with-mumpsinclude=DIR],
                  [Mumps(double,complexdouble) include directory]),
   [MUMPS_LIB="-I$with_mumpsinclude"])

  AC_ARG_WITH(mumpslib,
   AC_HELP_STRING([--with-mumpslib=<lib>],
                  [Mumps(double,complexdouble) and Pord link line]),
   [MUMPS_LIB="$with_mumpslib"])

 if test x"$ac_have_mumps" != x"no"; then
   ifelse([$1],,
          AC_DEFINE(HAVE_HIQLAB_MUMPS, 1,
                    [Define if you have Mumps(double,doublecomplex)]),
          [$1])
   :
 else
   ac_have_mumps=no
   $2
 fi

 AC_SUBST(MUMPS_DIR)
 AC_SUBST(MUMPS_INCLUDE)
 AC_SUBST(MUMPS_LIB)
 AC_SUBST(HAVE_HIQLAB_MUMPS_TRUE)
])dnl AC_WITH_MUMPS


dnl =====================================================================
dnl Optional METIS link copied from SuperLU Netsolve 1.4.1 configuration
dnl
AC_DEFUN([AC_WITH_METIS],
 [ac_have_metis="no"
  AC_ARG_WITH(metis,  
   AC_HELP_STRING([--with-metis=<lib>], 
                  [use METIS library <lib>]),
   [ac_have_metis="yes";
    METIS_LIB="$with_metis"],
   [ac_have_metis="no";
    METIS_LIB=""])

 if test x"$ac_have_metis" != x"no"; then
   ifelse([$1],,
          AC_DEFINE(HAVE_METIS, 1,
                    [Define if you have METIS]),
          [$1])
   :
 else
   ac_have_metis=no
   $2
 fi

 AC_SUBST(METIS_LIB)
])dnl AC_WITH_METIS


dnl =====================================================================
dnl Optional PARMETIS link copied from SuperLU Netsolve 1.4.1 configuration
dnl
AC_DEFUN([AC_WITH_PARMETIS],
 [ac_have_parmetis="no"
  AC_ARG_WITH(parmetis,  
   AC_HELP_STRING([--with-parmetis=<lib>], 
                  [use PARMETIS library <lib>]),
   [ac_have_parmetis="yes";
    PARMETIS_LIB="$with_parmetis"],
   [ac_have_parmetis="no";
    PARMETIS_LIB=""])

 if test x"$ac_have_parmetis" != x"no"; then
   ifelse([$1],,
          AC_DEFINE(HAVE_PARMETIS, 1,
                    [Define if you have PARMETIS]),
          [$1])
   :
 else
   ac_have_parmetis=no
   $2
 fi

 AC_SUBST(PARMETIS_LIB)
])dnl AC_WITH_PARMETIS


dnl =====================================================================
dnl Optional UMFPACK link taken from Netsolve 1.4.1 configuration
dnl
AC_DEFUN([AC_WITH_UMFPACK], [
AC_REQUIRE([ACX_BLAS])
ac_umfpack_ok="no"

dnl Check for UMFPACK directories
AC_ARG_WITH(umfpack,
  AC_HELP_STRING([--with-umfpack=<dir>],
                 [Use UMFPACK directory <dir>]),
  [UMF_INCLUDE="-I$with_umfpack/UMFPACK/Include -I$with_umfpack/AMD/Include";
   UMF_LIBDIR="-L$with_umfpack/UMFPACK/Lib -L$with_umfpack/AMD/Lib"],
  [UMF_INCLUDE="";
   UMF_LIBDIR="" ])

dnl Don't use UMFPACK if no BLAS
if test "x$acx_blas_ok" != xyes; then
  ac_umfpack_ok=noblas
fi

dnl Check that libamd, libumfpack are valid and umfpack.h is present
ac_umfpack_CFLAGS="$CFLAGS"
ac_umfpack_LIBS="$LIBS"
CFLAGS="$UMF_INCLUDE $UMF_LIBDIR $CFLAGS"
LIBS="$LIBS $BLAS_LIBS $FLIBS"

if test "x$ac_umfpack_ok" == xno; then
  AC_CHECK_LIB(amd, amd_order,
    [AC_CHECK_LIB(umfpack, umfpack_di_solve,
      [AC_CHECK_HEADERS([umfpack.h], [ac_umfpack_ok="yes"])],
      [], ["-lamd"])])
fi

CFLAGS="$ac_umfpack_CFLAGS"
LIBS="$ac_umfpack_LIBS"

UMF_CFLAGS="$UMF_INCLUDE"
UMF_LIBS="$UMF_LIBDIR -lumfpack -lamd"

AC_SUBST(UMF_CFLAGS)
AC_SUBST(UMF_LIBS)

dnl Finally, execute ACTION-IF-FOUND/ACTION-IF-NOT-FOUND:
if test x"$ac_umfpack_ok" = xyes; then
  ifelse([$1],,
         AC_DEFINE(HAVE_UMFPACK,1,[Define if you have UMFPACK.]),
         [$1])
  :
else
  ac_umfpack_ok="no"
  $2
fi

])dnl AC_WITH_UMFPACK


dnl =====================================================================
dnl Check for Matlab MEX program.  Hacked together using macros from
dnl an old version of HQP/Omuses and my own inventions.
dnl
AC_DEFUN([AC_PROG_MEX], [
AC_ARG_WITH(mex,
  AC_HELP_STRING([--with-mex=MEX_PROGRAM], [location of MEX script]),
  [if test x"$with_mex" != "xno"; then
     MEX="$with_mex";
   fi],
  [AC_PATH_PROGS([MEX], [mex mex.bat])
])

MEXEXT=""
target=`uname`
case "$target" in
  Linux*)
    cpu=`uname -p`
    case "$cpu" in
      x86_64*)
        MEXEXT="mexa64"
        MEX="$MEX -largeArrayDims"
        ;;
      *)
        MEXEXT="mexglx"
        ;;
    esac
    ;;
  Darwin*)
    MEXEXT="mexmac"
    ;;
  SunOS*)
    MEXEXT="mexsol"
    ;;
  OSF1*)
    MEXEXT="mexaxp"
    ;;
  IRIX*)
    MEXEXT="mexsg"
    ;;
  HP-UX*)
    MEXEXT="mexhpux"
    ;;
  CYGWIN*|MSYS*|MINGW*)
    MEXEXT="dll"
    ;;
  CL)
    MEXEXT="dll"
    ;;
esac

AC_SUBST(MEX)
AC_SUBST(MEXEXT)

dnl Finally, execute ACTION-IF-FOUND/ACTION-IF-NOT-FOUND:
if test x"$with_mex" = xyes; then
  ifelse([$1],,
         AC_DEFINE(HAVE_MEX,1,[Define if you have MEX program.]),
         [$1])
  :
else
  with_mex=no
  $2
fi

])dnl ACX_MEX


dnl =====================================================================
dnl Available from the GNU Autoconf Macro Archive at:
dnl http://www.gnu.org/software/ac-archive/htmldoc/ac_prog_perl_version.html
dnl
AC_DEFUN([AC_PROG_PERL_VERSION],[dnl
# Make sure we have perl
if test -z "$PERL"; then
AC_CHECK_PROG(PERL,perl,perl)
fi

# Check if version of Perl is sufficient
ac_perl_version="$1"

if test "x$PERL" != "x"; then
  AC_MSG_CHECKING(for perl version greater than or equal to $ac_perl_version)
  # NB: It would be nice to log the error if there is one, but we cannot rely
  # on autoconf internals
  $PERL -e "use $ac_perl_version;" > /dev/null 2>&1
  if test $? -ne 0; then
    AC_MSG_RESULT(no);
    $3
  else
    AC_MSG_RESULT(ok);
    $2
  fi
else
  AC_MSG_WARN(could not find perl)
fi
])dnl

dnl =====================================================================
dnl Optional Trilinos link copied from SuperLU Netsolve 1.4.1 configuration
dnl
AC_DEFUN([AC_WITH_TRILINOS],
 [ac_have_trilinos="no"
  AC_ARG_WITH(trilinos,  
   AC_HELP_STRING([--with-trilinos=DIR], 
                  [location of Trilinos root directory]),
   [ac_have_trilinos="yes";
    TRILINOS_DIR="$with_trilinos";
    TRILINOS_LIB_DIR="$with_trilinos"],
   [ac_have_trilinos="no";
    TRILINOS_DIR="";
    TRILINOS_LIB_DIR=""])

  AC_ARG_WITH(trilinoslibdir,
   AC_HELP_STRING([--with-trilinoslibdir=DIR],
                  [location of Trilinos implementation root directory]),
   [TRILINOS_LIB_DIR="$with_trilinoslibdir"])

 if test x"$ac_have_trilinos" != x"no"; then
   ifelse([$1],,
          AC_DEFINE(HAVE_TRILINOS, 1,
                    [Define if you have Trilinos]),
          [$1])
   :
 else
   ac_have_trilinos=no
   $2
 fi

 AC_SUBST(TRILINOS_DIR)
 AC_SUBST(TRILINOS_LIB_DIR)
])dnl AC_WITH_TRILINOS

dnl =====================================================================
dnl Optional Petsc link copied from SuperLU Netsolve 1.4.1 configuration
dnl
AC_DEFUN([AC_WITH_PETSC],
 [ac_have_petsc="no"
  AC_ARG_WITH(petsc,  
   AC_HELP_STRING([--with-petsc=DIR], 
                  [location of Petsc root directory(PETSC_DIR)]),
   [ac_have_petsc="yes";
    PETSC_DIR="$with_petsc";
    PETSC_ARCH=""],
   [ac_have_petsc="no";
    PETSC_DIR="";
    PETSC_ARCH=""])

  AC_ARG_WITH(petscarch,
   AC_HELP_STRING([--with-petscarch=NAME],
                  [name of Petsc architecture]),
   [PETSC_ARCH="$with_petscarch"])

 if test x"$ac_have_petsc" != x"no"; then
   ifelse([$1],,
          AC_DEFINE(HAVE_PETSC, 1,
                    [Define if you have Petsc]),
          [$1])
   :
 else
   ac_have_petsc=no
   $2
 fi

 AC_SUBST(PETSC_DIR)
 AC_SUBST(PETSC_ARCH)
])dnl AC_WITH_PETSC


dnl =====================================================================
dnl Optional Slepc link copied from SuperLU Netsolve 1.4.1 configuration
dnl
AC_DEFUN([AC_WITH_SLEPC],
 [ac_have_SLEPC="no"
  AC_ARG_WITH(slepc,  
   AC_HELP_STRING([--with-slepc=DIR], 
                  [location of Slepc root directory(SLEPC_DIR)]),
   [ac_have_slepc="yes";
    SLEPC_DIR="$with_slepc";
    HAVE_HIQLAB_SLEPC_TRUE=""],
   [ac_have_slepc="no";
    SLEPC_DIR="";
    HAVE_HIQLAB_SLEPC_TRUE="#"])

 if test x"$ac_have_slepc" != x"no"; then
   ifelse([$1],,
          AC_DEFINE(HAVE_HIQLAB_SLEPC, 1,
                    [Define if you have Slepc]),
          [$1])
   :
 else
   ac_have_slepc=no
   $2
 fi

 AC_SUBST(SLEPC_DIR)
 AC_SUBST(HAVE_HIQLAB_SLEPC_TRUE)
])dnl AC_WITH_SLEPC

