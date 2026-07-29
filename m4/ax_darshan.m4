AC_DEFUN([AX_DARSHAN],[
	AC_ARG_WITH([darshan-libdir],
	AS_HELP_STRING([--with-darshan-libdir=DIR],
	[Specify the library directory of Darshan to link with the profiler.]),
	[	   
	if test "x${withval}" = xno; then
		with_darshan=no
	elif test -d "$withval"; then
		ac_darshan_libdir="$withval";
		DARSHAN_LDFLAGS="-L$ac_darshan_libdir -Wl,-rpath=$ac_darshan_libdir"
		with_darshan=yes
	else
		AC_MSG_ERROR([Darshan library directory '$withval' was not found])
	fi
	],[with_darshan=no])

    ax_darshan_ok=no
    if test "x${with_darshan}" != xno; then
        if test -d "$ac_darshan_libdir"; then
            DARSHAN_LIBS="-ldarshan"
            darshan_save_LIBS="$LIBS"
            darshan_save_LDFLAGS="$LDFLAGS"
            LIBS="$DARSHAN_LIBS $LIBS"
            LDFLAGS="$DARSHAN_LDFLAGS $LDFLAGS"
            AC_MSG_CHECKING([for Darshan in $ac_darshan_libdir])
            AC_LINK_IFELSE([AC_LANG_PROGRAM([], [])],
                           [ax_darshan_ok=yes], [ax_darshan_ok=no])
            AC_MSG_RESULT([$ax_darshan_ok])
            LIBS="$darshan_save_LIBS"
            LDFLAGS="$darshan_save_LDFLAGS"
            if test "x${ax_darshan_ok}" != xyes; then
                AC_MSG_ERROR([Darshan requested but libdarshan was not found])
            fi
            LDFLAGS="$DARSHAN_LDFLAGS $LDFLAGS"
        fi
    fi
    AC_SUBST(DARSHAN_LIBS)
])
