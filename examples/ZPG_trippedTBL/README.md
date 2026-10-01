# Flow
This example is a zero pressure gradient (ZPG) turbulent boundary layer (TBL)
that starts from a laminar Blasius profile which is then tripped about 10*delta99
downstream of the inlet.

The inspiration for this example comes from Schlatter & Örlü, JFM, 2010
The mesh is a shorter version of the computational domain used in that work.

Compile using:

```
makeneko trip.f90 TBL_ZPG.f90
```

The case file is setup to include all the important parameters, but is not
ready for a production run, and neither is the mesh.

The current mesh was generated using genmeshbox with the following command:

```
genmeshbox 0.0 162.5 0.0 65.0 0.0 78.0 96 64 112 .false. .false. .true.
```

The resolutions should be okay for a reasonable simulation 
(wall-normal resolution is even too fine)
but the streamwise extent needs to be extended if the goal is to reach Re_theta=4000
this can be done with the following genmeshbox command:

```
genmeshbox 0.0 1950.0 0.0 65.0 0.0 78.0 1152 64 112 .false. .false. .true.
```

IMPORTANT NOTE: the z averaging in statistics may not work with FP32 (--enable-real=sp)


# Tripping
Introduces a volume forcing along a line parallel to the y or z direction.

`trip.f90` implements the `trip` module to be used in the user file. See instructions
in the header of the file.

Make sure to enable the user source term in the case file as well.
