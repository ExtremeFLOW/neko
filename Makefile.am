ACLOCAL_AMFLAGS = -I m4 --install
pkgconfigdir = $(libdir)/pkgconfig
pkgconfig_DATA = neko.pc

SUBDIRS = src

if !ENABLE_NVSHMEM
if ENABLE_CONTRIB
SUBDIRS += contrib
endif

if ENABLE_PFUNIT
SUBDIRS += tests/unit
endif
endif

if ENABLE_DOXYGEN
SUBDIRS += doc
endif

CLEANFILES = *.mod\
	     *.smod


bin_SCRIPTS = makeneko
EXTRA_DIST = \
	makeneko.in\
	regen.sh\
	flinter_rc.yml\
	tests/bench/ax/setup.f90\
	tests/bench/ax/Makefile\
	tests/bench/ax/axkernel.f90\
	tests/bench/ax/driver.f90\
	tests/bench/gs/Makefile\
	tests/bench/gs/driver.f90\
	tests/bench/tgv32/tgv.case\
	tests/bench/tgv32/README.md\
	tests/bench/tgv32/extract.sh\
	tests/bench/tgv32/tgv.f90\
	examples/turb_channel/turb_channel.case\
	examples/turb_channel/README.md\
	examples/turb_channel/turb_channel.f90\
	examples/turb_channel/box.nmsh\
	examples/turb_pipe/turb_pipe.case\
	examples/turb_pipe/turb_pipe_source.case\
	examples/turb_pipe/turb_pipe.f90\
	examples/turb_pipe/turb_pipe.nmsh\
	examples/turb_pipe/turb_pipe.geo\
	examples/turb_pipe/README.md\
	examples/hemi/hemi.nmsh\
	examples/hemi/README.md\
	examples/hemi/hemi.case\
	examples/lid/lid.f90\
	examples/lid/lid2d.nmsh\
	examples/lid/lid.nmsh\
	examples/lid/README.md\
	examples/lid/lid2d.case\
	examples/lid/lid.case\
	examples/lid/elem_dist.csv\
	examples/lid/plot_stat.py\
	examples/poisson/ax_poisson.f90\
	examples/poisson/setup.f90\
	examples/poisson/Makefile\
	examples/poisson/README.md\
	examples/poisson/data/4096.nmsh\
	examples/poisson/data/16384.nmsh\
	examples/poisson/data/256.nmsh\
	examples/poisson/data/2048.nmsh\
	examples/poisson/data/128.nmsh\
	examples/poisson/data/32768.nmsh\
	examples/poisson/data/8192.nmsh\
	examples/poisson/data/512.nmsh\
	examples/poisson/data/1024.nmsh\
	examples/poisson/driver.f90\
	examples/2d_cylinder/2d_cylinder.nmsh\
	examples/2d_cylinder/2d_cylinder.case\
	examples/2d_cylinder/2d_cylinder.f90\
	examples/2d_cylinder/2d_cylinder.geo\
	examples/2d_cylinder/README.md\
	examples/cyl_boundary_layer/cyl_bl_rot.case\
	examples/cyl_boundary_layer/cyl_bl.f90\
	examples/cyl_boundary_layer/README.md\
	examples/cyl_boundary_layer/cyl_bl_user.case\
	examples/cyl_boundary_layer/cyl_bl_basic.case\
	examples/cyl_boundary_layer/cyl.nmsh\
	examples/rayleigh_benard/rayleigh.f90\
	examples/rayleigh_benard/README.md\
	examples/rayleigh_benard/box.nmsh\
	examples/rayleigh_benard/rayleigh.case\
	examples/tgv/262144.nmsh\
	examples/tgv/tgv.case\
	examples/tgv/README.md\
	examples/tgv/32768.nmsh\
	examples/tgv/tgv.f90\
	examples/tgv/512.nmsh\
	examples/rayleigh_benard_cylinder/probes.csv\
	examples/rayleigh_benard_cylinder/rayleigh.f90\
	examples/rayleigh_benard_cylinder/README.md\
	examples/rayleigh_benard_cylinder/rayleigh.nmsh\
	examples/rayleigh_benard_cylinder/rayleigh.geo\
	examples/rayleigh_benard_cylinder/rayleigh.case\
	examples/rayleigh_benard_cylinder/postprocess_probes.py\
	examples/api/pyneko/README.md\
	examples/api/pyneko/cylinder.case\
	examples/api/pyneko/cylinder.py\
	examples/api/pyneko/cylinder.nmsh\
	examples/api/c/README.md\
	examples/api/c/cylinder.case\
	examples/api/c/cylinder.c\
	examples/api/c/cylinder.nmsh\
	tests/reframe/checks.py\
	tests/reframe/settings.py\
	tests/reframe/README.md\
	tests/reframe/src/get-tgv-dns.sh\
	tests/reframe/src/minihemi.case.template\
	tests/reframe/src/rayleigh.case.template\
	tests/reframe/src/rayleigh.f90\
	tests/reframe/src/rocm_select_gpu_device.sh\
	tests/reframe/src/tgv.case.template\
	tests/reframe/src/tgv.f90\
	tests/integration/default_cpu_launcher.sh\
	tests/integration/default_cuda_launcher.sh\
	tests/integration/conftest.py\
	tests/integration/testlib.py\
	tests/integration/pytest.ini\
	tests/integration/meshes/small_test_cyl.nmsh\
	tests/integration/tests/test_demo/test_suite.py\
	tests/integration/tests/test_examples/test_examples.py\
	tests/integration/tests/test_cylinder/test_cylinder.py\
	tests/integration/tests/test_cylinder/cylinder_part1.case\
	tests/integration/tests/test_cylinder/reflog1_dp.csv\
	tests/integration/tests/test_cylinder/reflog1_sp.csv\
	tests/integration/tests/test_cylinder/reflog2_dp.csv\
	tests/integration/tests/test_cylinder/reflog2_sp.csv\
	patches/nvhpc_bge.patch\
	patches/json-fortran71.patch

if ENABLE_FLINT
lint:
	$(FLINT) score -d 0 src/ -r flinter_rc.yml

check-local: lint

endif
