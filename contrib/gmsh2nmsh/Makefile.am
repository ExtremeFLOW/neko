# gmsh2nmsh has its own configure script, which Neko's configure runs. It is
# built with the plain Fortran compiler FC and none of Neko's dependencies.
# It can also be compiled on its own: gfortran -O2 gmsh2nmsh.F90 -o gmsh2nmsh
bin_PROGRAMS = gmsh2nmsh
gmsh2nmsh_SOURCES = gmsh2nmsh.F90

man_MANS = gmsh2nmsh.1

EXTRA_DIST = README.md

CLEANFILES = gmsh2nmsh_util.mod gmsh2nmsh_reader.mod
