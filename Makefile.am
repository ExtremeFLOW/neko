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
	bench/ax/setup.f90\
	bench/ax/Makefile\
	bench/ax/axkernel.f90\
	bench/ax/driver.f90\
	bench/gs/Makefile\
	bench/gs/driver.f90\
	bench/tgv32/tgv.case\
	bench/tgv32/README\
	bench/tgv32/extract.sh\
	bench/tgv32/tgv.f90\
	examples/hemi/hemi.case\
	examples/hemi/hemi.nmsh\
	examples/tgv/tgv.case\
	examples/tgv/tgv.f90\
	examples/tgv/512.nmsh\
	examples/tgv/32768.nmsh\
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
	examples/poisson/setup.f90\
	examples/cyl_boundary_layer/cyl_bl_rot.case\
	examples/cyl_boundary_layer/cyl_bl.f90\
	examples/cyl_boundary_layer/README.md\
	examples/cyl_boundary_layer/cyl_bl_user.case\
	examples/cyl_boundary_layer/cyl_bl_basic.case\
	examples/cyl_boundary_layer/cyl.nmsh\
	reframe/checks.py\
	reframe/settings.py\
	reframe/README.md\
	reframe/src/get-tgv-dns.sh\
	reframe/src/minihemi.case.template\
	reframe/src/rocm_select_gpu_device\
	reframe/src/tgv.case.template\
	reframe/src/tgv.f90\
	examples/rayleigh-benard/rayleigh.case\
	examples/rayleigh-benard/rayleigh.f90\
	examples/rayleigh-benard/box.nmsh

if ENABLE_FLINT
lint:
	$(FLINT) score -d 0 src/ -r flinter_rc.yml

check-local: lint

endif
