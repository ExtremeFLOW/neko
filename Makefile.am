pkgconfigdir = $(libdir)/pkgconfig
pkgconfig_DATA = neko.pc

SUBDIRS = src contrib

if ENABLE_PFUNIT
SUBDIRS += tests
endif

if ENABLE_DOXYGEN
SUBDIRS += doc
endif


bin_SCRIPTS = makeneko
EXTRA_DIST = \
	makeneko.in\
	regen.sh\
	flinter_rc.yml\
	examples/hemi/hemi.case\
	examples/hemi/hemi.nmsh\
	examples/tgv/tgv.case\
	examples/tgv/tgv.f90\
	examples/tgv/512.nmsh\
	examples/poisson/data/512.nmsh\
	examples/poisson/data/1024.nmsh\
	examples/poisson/data/2048.nmsh\
	examples/poisson/data/4096.nmsh\
	examples/poisson/data/8192.nmsh\
	examples/poisson/data/16384.nmsh\
	examples/poisson/data/32768.nmsh\
	examples/poisson/Makefile\
	examples/poisson/ax_poisson.f90\
	examples/poisson/driver.f90\
	examples/poisson/setup.f90

if ENABLE_FLINT
lint:
	$(FLINT) score -d 0 src/ -r flinter_rc.yml

check-local: lint

endif
