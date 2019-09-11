lib_LIBRARIES = libneko.a
libneko_a_SOURCES = \
	num_types.f90\
	nmsh.f90\
	re2.f90\
	map.f90\
	comm.f90\
	mpi_types.f90 \
	datadist.f90 \
	utils.f90\
	mxm_wrapper.f90\
	mxm_std.f90\
	speclib.f90 \
	math.f90\
	space.f90\
	entity.f90\
	point.f90\
	element.f90\
	quad.f90\
	hex.f90\
	htable.f90\
	mesh.f90\
	mesh_field.f90\
	field.f90\
	rea.f90\
	generic_file.f90\
	map_file.f90\
	re2_file.f90\
	rea_file.f90\
	vtk_file.f90\
	nmsh_file.f90\
	file.f90\
	neko.f90


pkginclude_HEADERS =  $(wildcard *.mod)

CLEANFILES =*.mod

if ENABLE_MAKEDEPF90
depend:
	$(MAKEDEPF90) $(libneko_a_SOURCES) > .depends
endif

include	.depends
