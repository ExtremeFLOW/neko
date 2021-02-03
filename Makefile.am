pkgconfigdir = $(libdir)/pkgconfig
pkgconfig_DATA = neko.pc

SUBDIRS = src contrib

if ENABLE_PFUNIT
SUBDIRS += tests
endif

if ENABLE_DOXYGEN
SUBDIRS += doc
endif


dist_data_DATA = regen.sh
bin_SCRIPTS = makenek
EXTRA_DIST = makenek.in
