# ![Neko](https://user-images.githubusercontent.com/750135/169531665-313c3471-50d1-4c44-964a-fee7312d6459.png)
![CI](https://github.com/ExtremeFLOW/neko/workflows/CI/badge.svg)
## About
Neko is a portable framework for high-order spectral element flow simulations. Written in modern Fortran, Neko adopts an object-oriented approach, allowing multi-tier abstractions of the solver stack and facilitating various hardware backends ranging from general-purpose processors, CUDA and HIP enabled accelerators to SX-Aurora vector processors. Neko has its roots in the spectral element code Nek5000 from UChicago/ANL, from where many of the namings, code structure and numerical methods are adopted.

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
Documentation for Neko is available at https://extremeflow.github.io/neko/.

To generate the documentation, you need to have both doxygen and dot installed (they will be picked up by configure). Once installed, you should be able to generate the documentation with
```bash
make html
```
## Publications using Neko
* Jansson, N., Karp, M., Podobas, A., Markidis, S. and Schlatter, P., 2021. *Neko: A modern, portable, and scalable framework for high-fidelity computational fluid dynamics*. arXiv preprint arXiv:2107.01243.
* Jansson, N., 2021. *Spectral Element Simulations on the NEC SX-Aurora TSUBASA*. In proc. HPCAsia 2021.
* Karp, M., Podobas, A., Kenter, T., Jansson, N., Plessl, C., Schlatter, P. and Markidis, S., 2022. *A high-fidelity flow solver for unstructured meshes on field-programmable gate arrays: Design, evaluation, and future challenges*. In proc. HPCAsia 2022.
* Karp, M., Jansson, N., Podobas, A., Schlatter, P., and Markidis, S., 2022. *Reducing Communication in the Conjugate Gradient Method: A Case Study on High-Order Finite Elements*. In proc. PASC 2022.
* Karp, M., Massaro, D., Jansson, N., Hart, A., Wahlgren, J., Schlatter, P., and Markidis, S., 2022. *Large-Scale Direct Numerical Simulations of Turbulence Using GPUs and Modern Fortran*. arXiv preprint arXiv:2207:07098.

## Acknowledgments
The development of Neko was supported by the European Commission Horizon 2020 project grant *EPiGRAM-HS: Exascale Programming Models for Heterogeneous Systems* (grant reference 801039), the Swedish Research Council project grant *Efficient Algorithms for Exascale Computational Fluid Dynamics* (grant reference 2019-04723) and the SeRC Exascale Simulation Software Initiative (SESSI). The Neko logo was designed by Robert Hansen Jagrelius.
