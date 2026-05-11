AC_DEFUN([AX_HYPRE],[
    AC_ARG_WITH([hypre],
    AS_HELP_STRING([--with-hypre=DIR],
    [Compile with support for hypre]),
    [
    if test -d "$withval"; then
        ac_hypre_path="$withval";
        AS_IF([test -d "$ac_hypre_path/lib64"],[suffix="64"],[suffix=""])
        HYPRE_LDFLAGS="-L$ac_hypre_path/lib$suffix -L$ac_hypre_path"
        HYPRE_CPPFLAGS="-I${ac_hypre_path}/include"
        HYPRE_CFLAGS=${HYPRE_CPPFLAGS}
    fi
    ],[with_hypre=no])

    if test "x${with_hypre}" != xno; then
        HYPRE_LIBS="-lHYPRE"
        if test "x${have_cuda}" = xyes; then
          HYPRE_LIBS="${HYPRE_LIBS} -lm -lcudart -lcublas -lcusparse -lcurand -lcusolver -lstdc++ -lumpire -lcamp"
          HYPRE_LDFLAGS="${HYPRE_LDFLAGS} -L$ac_hypre_path/lib"
        fi
        if test "x${have_hip}" = xyes; then
          HYPRE_LIBS="${HYPRE_LIBS} -lm -lrocblas -lrocsparse -lrocsolver -lrocrand -lumpire -lcamp"
          HYPRE_LDFLAGS="${HYPRE_LDFLAGS} -L$ac_hypre_path/lib"
        fi

        CPPFLAGES_SAVED="$CPPFLAGS"
        CPPFLAGS="${HYPRE_CPPFLAGS} ${CPPFLAGS}"
        LDFLAGS_SAVED="$LDFLAGS_SAVED"
        LDFLAGS="${HYPRE_LDFLAGS} ${LDFLAGS}"
        LIBS_SAVED="$LIBS"
        LIBS="${HYPRE_LIBS} ${LIBS}"

        AC_MSG_NOTICE([ ${LDFLAGS} ])
        AC_MSG_NOTICE([ ${CPPFLAGS} ])
        AC_CHECK_HEADER([HYPRE.h],[have_hypre_h=yes],
           [have_hypre_h=no;AC_MSG_WARN([HYPRE header HYPRE.h not found])])

        if test "x$have_hypre_h" = xyes; then
            AC_CHECK_LIB([HYPRE],[HYPRE_IJMatrixCreate],
            [have_hypre=yes],
            [
              AC_MSG_WARN([Could not link against HYPRE (HYPRE_IJMatrixCreate not found)])
              have_hypre=no
            ])
        else
            have_hypre=no
        fi

        if test "x$have_hypre" = xyes; then
            AC_DEFINE(HAVE_HYPRE,1,[Define if you have the HYPRE library.])
        fi

        CPPFLAGS="$CPPFLAGS_SAVED"
        LDFLAGS="$LDFLAGS_SAVED"
        LIBS="$LIBS_SAVED"
    fi
])# AX_HYPRE
