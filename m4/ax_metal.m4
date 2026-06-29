#
# ----------------------------------------------------------------------------
# "THE BEER-WARE LICENSE" (Revision 42):
# <njansson@kth.se> wrote this file. As long as you retain this notice you
# can do whatever you want with this stuff. If we meet some day, and you think
# this stuff is worth it, you can buy me a beer in return Niclas Jansson
# ----------------------------------------------------------------------------
#

AC_DEFUN([AX_METAL],[
	AC_ARG_WITH([metal],
		    AS_HELP_STRING([--with-metal],
		    [Compile with Apple Metal backend (macOS only)]),
		    [], [with_metal=no])
	metal_bcknd="0"
	have_metal="no"
	if test "x${with_metal}" != xno; then
		case $host_os in
		  darwin*)
		    ;;
		  *)
		    AC_MSG_ERROR([Metal backend is only supported on macOS])
		    ;;
		esac

		dnl Apple GPUs lack fp64 support, all Metal kernels are
		dnl single precision
		if test "x${enable_real}" != xsp; then
		   AC_MSG_ERROR([Metal backend requires --enable-real=sp])
		fi

		AC_MSG_CHECKING([for Metal framework])
		AC_LANG_PUSH([C])
		AC_LANG_ASSERT([C])
		LIBS_SAVED="$LIBS"
		METAL_LIBS="-framework Metal -framework CoreGraphics -framework Foundation"
		LIBS="$METAL_LIBS $LIBS"

		AC_LINK_IFELSE([AC_LANG_SOURCE([
		#include <TargetConditionals.h>
		#if !TARGET_OS_MAC
		#error Not macOS
		#endif
		int main(void) { return 0; }
		])],
		[have_metal=yes],[have_metal=no])

		if test x"${have_metal}" = xyes; then
		   metal_bcknd="1"
		   AC_DEFINE(HAVE_METAL,1,[Define if you have Apple Metal.])
		   AC_MSG_RESULT([yes])
		   METAL_OBJCFLAGS="-fobjc-arc"
		   AC_SUBST(METAL_LIBS)
		   AC_SUBST(METAL_OBJCFLAGS)
		else
		   LIBS="$LIBS_SAVED"
		   AC_MSG_ERROR([Metal framework not found])
		fi
		AC_LANG_POP([C])

		dnl Objective-C compiler for the Metal device layer
		dnl (CC is often an MPI wrapper, which might not be able to
		dnl  compile Objective-C, default to clang)
		AC_CHECK_PROGS(OBJC, [clang cc], [clang])
		if test -z "${OBJCFLAGS}"; then
		   OBJCFLAGS="-g -O2"
		fi
		AC_MSG_CHECKING([whether ${OBJC} accepts Objective-C])
		cat > conftest.m << ACEOF
#import <Foundation/Foundation.h>
int main(void) { @autoreleasepool { } return 0; }
ACEOF
		if ${OBJC} ${OBJCFLAGS} ${METAL_OBJCFLAGS} -c conftest.m \
		     -o conftest.o >/dev/null 2>&1; then
		   AC_MSG_RESULT([yes])
		else
		   AC_MSG_RESULT([no])
		   rm -f conftest.m conftest.o
		   AC_MSG_ERROR([${OBJC} can not compile Objective-C])
		fi
		rm -f conftest.m conftest.o
		AC_SUBST(OBJC)
		AC_SUBST(OBJCFLAGS)

		dnl Check for the Metal shader compiler, required to build
		dnl the embedded shader library (requires Xcode; since
		dnl Xcode 26 the Metal toolchain is a separate download)
		AC_PATH_PROG(XCRUN, xcrun, "no")
		if test "x${XCRUN}" = xno; then
		   AC_MSG_ERROR([xcrun not found (Xcode is required for the Metal backend)])
		fi
		AC_MSG_CHECKING([for Metal shader compiler])
		if ${XCRUN} -sdk macosx metal --version >/dev/null 2>&1; then
		   AC_MSG_RESULT([yes])
		else
		   AC_MSG_RESULT([no])
		   AC_MSG_ERROR([Metal shader compiler not found (install Xcode, or run xcodebuild -downloadComponent MetalToolchain)])
		fi
	fi
	AC_SUBST(metal_bcknd)
])
