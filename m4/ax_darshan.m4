AC_DEFUN([AX_DARSHAN],[
	AC_ARG_WITH([darshan-libdir],
	AS_HELP_STRING([--with-darshan-libdir=DIR],
	[Specify the library directory of Darshan to link with the profiler.]),
	[	   
	if test -d "$withval"; then
		ac_darshan_libdir="$withval";
		DARSHAN_LDFLAGS="-L$ac_darshan_libdir -Wl,-rpath=$ac_darshan_libdir -ldarshan"
	fi
	],[with_darshan=no])

    if test "x${with_darshan}" != xno; then
        if test -d "$ac_darshan_libdir"; then
            DARSHAN_LIBS="-ldarshan"
            AC_SUBST(DARSHAN_LIBS)
            LDFLAGS="$DARSHAN_LDFLAGS $LDFLAGS"
            export LDFLAGS
        fi
    fi
])
