Generate the mesh with

```bash
genmeshbox 0 1 0 1 -0.1 1 10 10 1 .false. .false. .true.
```

Run with

```bash
makeneko lid.f90
mpirun -np 4 ./neko lid2d.case
```