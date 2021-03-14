#
# ----------------------------------------------------------------------------
# "THE BEER-WARE LICENSE" (Revision 42):
# <njansson@kth.se> wrote this file. As long as you retain this notice you 
# can do whatever you want with this stuff. If we meet some day, and you think
# this stuff is worth it, you can buy me a beer in return Niclas Jansson
# ----------------------------------------------------------------------------
#

AC_DEFUN([AX_SX],[
	AC_MSG_CHECKING([for a NEC SX-Aurora system])
	AC_LANG_PUSH([C])
	AC_EGREP_CPP(yes,
	[#if defined(__ve__))
	  yes
	 #endif
	],
	[AC_MSG_RESULT([yes])
	is_sx="1"],
	[is_sx="0"
	AC_MSG_RESULT([no])])
	AC_LANG_POP([C])
	AC_SUBST(is_sx)])



