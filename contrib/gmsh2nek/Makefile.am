bin_PROGRAMS = gmsh2nek

gmsh2nek_SOURCES =  mod_SIZE.f90\
                    gmsh2nek.f90\
		    byte.c\
	            speclib.f90\
		    mxm.f90

gmsh2nek_LDFLAGS = $(LDFLAGS) $(LIBS) 

AM_FCFLAGS = -I@top_builddir@/src
