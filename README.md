# ![Neko](https://user-images.githubusercontent.com/750135/169531665-313c3471-50d1-4c44-964a-fee7312d6459.png)
![CI](https://github.com/ExtremeFLOW/neko/workflows/CI/badge.svg) ![develop](https://github.com/ExtremeFLOW/neko/workflows/develop/badge.svg) [![DOI](https://zenodo.org/badge/338607716.svg)](https://zenodo.org/doi/10.5281/zenodo.6631055)
## About
Neko is a portable framework for high-order spectral element flow simulations. Written in modern Fortran, Neko adopts an object-oriented approach, allowing multi-tier abstractions of the solver stack and facilitating various hardware backends ranging from general-purpose processors, CUDA and HIP enabled accelerators to SX-Aurora vector processors. Neko has its roots in the spectral element code Nek5000 from UChicago/ANL, from where many of the namings, code structure and numerical methods are adopted.


## Cloning the project

```bash
git clone https://github.com/ExtremeFLOW/neko
```

## Documentation
Documentation for Neko is available at https://neko.cfd and most things related to the code, cases, and different features are described there. The documentation is always improving, in large part due to our active users and if something is missing or hard to understand, don't be afraid to create an issue or create a PR. It is a great way to help us improve and also to start getting involved in the project.

## Building the project
To build the project you will need: A Fortran compiler supporting the Fortran-08 standard, a working MPI installation, JSON-Fortran, and BLAS/lapack. Optional dependencies are gslib and ParMETIS. We use automake to build the project. These instructions should work in general, but as the project is quickly developing, things might change. While we assume MPI and BLAS are installed, if JSON-Fortran is not already available it can be cloned, installed, and the correct paths set with the following commands (Skip this step if you already have an installation of JSON-Fortran).

```bash
export JSON_INSTALL=/path/to/json-fortran_install # Where you want to install json-fortran
```
```bash
git clone --depth 1 https://github.com/ExtremeFLOW/json-fortran/
cd json-fortran
mkdir build && cd build
cmake -DCMAKE_INSTALL_PREFIX=${JSON_INSTALL} -DUSE_GNU_INSTALL_CONVENTION=ON ..
make -j4 && make install && cd ../../
export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${JSON_INSTALL}/lib/ #On some systems lib should be replaced with lib64
export PKG_CONFIG_PATH=${PKG_CONFIG_PATH}:${JSON_INSTALL}/lib/pkgconfig 

```

A basic CPU version of Neko can then be installed according to the following
```bash
cd neko
./regen.sh
./configure --prefix=/path/to/neko_install # Where you want to install neko
make install
```
More detailed installation instructions and all the different options (such as how to install Neko for GPUs) can be found in the documentation available at https://neko.cfd. 

## Running examples
After the project has been built

```bash
cd examples/tgv
/path/to/neko_install/bin/makeneko tgv.f90
mpirun -np 4 ./neko tgv.case
```
If there is not a .f90 (user) file in the example, the standard executable `neko` can also be used, for example:
```bash
cd examples/hemi
mpirun -np 4 /path/to/neko_install/bin/neko hemi.case
```
only uses built-in functions and does not need a compiled user file. Whether you will need a user file or not depends on what functionality you want and this is also documented in the documentation.

## Publications using Neko
* Jansson, N., 2021. *Spectral Element Simulations on the NEC SX-Aurora TSUBASA*. In proc. HPCAsia 2021.
* Karp, M., Podobas, A., Kenter, T., Jansson, N., Plessl, C., Schlatter, P. and Markidis, S., 2022. *A high-fidelity flow solver for unstructured meshes on field-programmable gate arrays: Design, evaluation, and future challenges*. In proc. HPCAsia 2022.
* Karp, M., Jansson, N., Podobas, A., Schlatter, P., and Markidis, S., 2022. *Reducing Communication in the Conjugate Gradient Method: A Case Study on High-Order Finite Elements*. In proc. PASC 2022.
* Karp, M., Massaro, D., Jansson, N., Hart, A., Wahlgren, J., Schlatter, P., and Markidis, S., 2023. *Large-Scale Direct Numerical Simulations of Turbulence Using GPUs and Modern Fortran*. The International Journal of High Performance Computing Applications, 37, 5.
* Jansson, N., Karp, M., Perez, A., Mukha, T., Ju, Y., Liu, J., Páll, S., Laure, E., Weinkauf, T., Schumacher, J., Schlatter, P., Markidis, S., 2023. *Exploring the Ultimate Regime of Turbulent Rayleigh–Bénard Convection Through Unprecedented Spectral-Element Simulations*. SC '23: Proceedings of the International Conference for High Performance Computing, Networking, Storage and Analysis.
* Jansson, N., Karp, M., Podobas, A., Markidis, S. and Schlatter, P., 2024. *Neko: A modern, portable, and scalable framework for high-fidelity computational fluid dynamics*. Computer & Fluids, 275.

## Acknowledgments
The development of Neko was supported by the European Commission Horizon 2020 project grant *EPiGRAM-HS: Exascale Programming Models for Heterogeneous Systems* (grant reference 801039), the Swedish Research Council project grant *Efficient Algorithms for Exascale Computational Fluid Dynamics* (grant reference 2019-04723) and the SeRC Exascale Simulation Software Initiative (SESSI). The Neko logo was designed by Robert Hansen Jagrelius.


[<img src="https://raw.githubusercontent.com/zulip/zulip/143baa42432cde9f288bd202336ef2b11172f6e4/static/images/logo/zulip-icon-128x128.png" width="32"/>](https://zulip.com) Sponsored by Zulip, an open-source modern team chat app designed to keep both live and asynchronous conversations organized.
