# Neko
![CI](https://github.com/ExtremeFLOW/neko/workflows/CI/badge.svg)
## About
Neko is a portable framework for high-order spectral element flow simulations. Written in modern Fortran, Neko adopts an object-oriented approach, allowing multi-tier abstractions of the solver stack and facilitating various hardware backends ranging from general-purpose processors, CUDA and HIP enabled accelerators to SX-Aurora vector processors. Neko has its roots in the spectral element code Nek5000 from UChicago/ANL.

Neko is currently maintained and developed at KTH Royal Institute of Technology.

## Cloning the project

```bash
git clone https://github.com/ExtremeFLOW/neko
```

## Building the project
To build the project you will need: A Fortran compiler supporting the Fortran-08 standard, a working MPI installation and BLAS/lapack.
We use automake to build the project. These instructions should work in general, but as the project is quickly developing, things might change.

```bash
cd neko
./regen.sh
./configure --prefix=/path/to/neko_install --with-pfunit=/path/to/pFUnit/installed/PFUNIT-VERSION
make install
```
## Running examples
After the project has been built

```bash
cd examples/tgv
/path/to/neko_install/bin/makeneko tgv.f90
mpirun -np 4 ./neko tgv.case
```

## Testing the Code
Assuming you configured with pFUnit you should be able to test the code with
```bash
make check
```

## Documentation
To generate the documentation, you need to have both doxygen and dot installed (they will be picked up by configure). Once installed, you should be able to generate the documentation with
```bash
make html
```

## Acknowledgments
The development of Neko was supported by the European Commission Horizon 2020 project grant *EPiGRAM-HS: Exascale Programming Models for Heterogeneous Systems* (grant reference 801039), the Swedish Research Council project grant *Efficient Algorithms for Exascale Computational Fluid Dynamics* (grant reference 2019-04723) and the SeRC Exascale Simulation Software Initiative (SESSI).
