# Installing Neko {#installation}

Neko can be installed in various ways, either building directly from source, manually compiling all dependencies and Neko or via Spack. Pre-built Docker images are also provided for each release of Neko.

## Building from source

To build Neko, you will need a Fortran compiler supporting the Fortran-08 standard, autotools, a working MPI installation supporting the Fortran 2008 bindings (`mpi_f08`), BLAS/LAPACK and JSON-Fortran. Optional dependencies are gslib and ParMETIS. 

Follow the steps below to install the less common dependencies (e.g. JSON-Fortran).

### Dependencies

#### Building JSON Fortran 

Download and compile, at least version 0.7.1 of JSON Fortran from the main repository.
@note Neko requires JSON Fortran to be configured with `USE_GNU_INSTALL_CONVENTION`.

``` shell
git clone --depth=1 https://github.com/jacobwilliams/json-fortran.git
cd json-fortran && mkdir b && cd b
cmake -DCMAKE_INSTALL_PREFIX=/path/to/installation -DUSE_GNU_INSTALL_CONVENTION=ON ..
make install
```
Now ad the installation path to `PKG_CONFIG_PATH` (and if needed `LD_LIBRARY_PATH`).
@note On certain systems `lib` should be substituted with `lib64`

``` bash
export PKG_CONFIG_PATH=/path/to/installation/lib/pkgconfig:$PKG_CONFIG_PATH
export LD_LIBRARY_PATH=/path/to/installation/lib:$LD_LIBRARY_PATH
```

#### Building gslib (optional)

If you have a Nek5000 installation, use:

``` bash
export GSLIB=/path/to/Nek5000/3rd_party/gslib/gslib/src
```

If not, you should download and compile `gslib`. In a folder outside of Neko:

``` shell
git clone https://github.com/Nek5000/gslib.git
cd gslib
make
```

Check that `libgs.a` has been created:

``` shell
$ ls build/lib
libgs.a 
```

Now add the path to gslib to an environment variable `GSLIB`

``` shell
export GSLIB=$(pwd)/build
```

Later, when configuring Neko, add the following option to enable gslib

``` shell
 --with-gslib=${GSLIB}
```

Make sure you see the following message during the configuration:

``` shell
checking for fgslib_gs_setup in -lgs... yes
```

#### Building ParMETIS (optional)

The following steps is an example on how to build and install ParMETIS

``` shell
$ wget http://glaros.dtc.umn.edu/gkhome/fetch/sw/parmetis/parmetis-4.0.3.tar.gz
$ tar xzf parmetis-4.0.3.tar.gz && cd parmetis-4.0.3 && make config prefix=/parmetis_install_path
$ make -j$(nproc) && make install
```
@note ParMETIS might not install `metis.h`, check if it is found in `/parmetis_install_path/include`.  If this is not the case, repeat the same `make config prefix=/parmetis_install_path` and `make install` commands from the `metis` subfolder.


#### Bulding PFunit (optional)

To build the PFunit testing framework, please refers to the \subpage testing page

### Building Neko
Neko uses autotools as its build system. The first step is to run the `configure` script, located in the top directory.

``` shell
<path-to-neko>/configure FC=<Fortran compiler> CC=<C compiler> \
                         MPIFC=<MPI Fortran compiler> MPICC=<MPI C ompiler> \
                         FCFLAGS=<Fortran compiler flags> CFLAGS=<C compiler flags> \
                         --prefix=<installation path> [options]
```

In the above command, `[options]` refers to either optional features or packages. 

Features are enabled and disabled by passing either `--enable-FEATURE[=arg]` or `--disable-FEATURE` to `configure`. A list of currently supported features are given in the table below.

Name                  |  Description
----                  | ------------
`--enable-real=Xp`    | Specify working precision of REAL types:<br>`sp` -- `REAL(kind=4)`<br>`dp` -- `REAL(kind=8)` (default)<br>`qp` -- `REAL(kind=16)`<br>
`--enable-contrib`    | Compile various tools
`--enable-device-mpi` | Enable device aware MPI

Optional packages are controlled by passing either `--with-PACKAGE[=ARG]` or `--without-PACKAGE` to `configure`. A list of all supported optional packages are given in the table below.

Name                            | Description
----                            | -----------
`--with-blas=<lib>`             | Use BLAS library `<lib>`
`--with-lapack=<lib>`           | Use LAPACK library `<lib>`
`--with-metis=DIR`              | Directory for metis
`--with-metis-libdir=LIBDIR`    | Directory for metis library (if different)
`--with-parmetis=DIR`           | Compile with support for parmetis library
`--with-parmetis-libdir=LIBDIR` | Directory for parmetis library (if different)
`--with-adios2=DIR`             | Compile with support for ADIOS2
`--with-gslib=DIR`              | Compile with support for gslib
`--with-libxsmm`                | Compile with support for libxsmm
`--with-hip=DIR`                | Compile with HIP backend
`--with-cuda=DIR`               | Compile with CUDA backend
`--with-opencl=DIR`             | Compile with OpenCL backend
`--with-nvtx=DIR`               | Compile with support for NVTX
`--with-roctx=DIR`              | Compile with support for ROCTX
`--with-pfunit=DIR`             | Directory for pFUnit (see \subpage testing)   

@note Accelerators backends are not enabled as a feature in Neko, but rather via optional packages.

Once configured, to compile and install Neko issue `make` followed by `make install`

#### Compiling Neko for CPU or SX-Aurora

For a standard CPU or SX-Aurora build of Neko, simply run the `configure` script as given above, using appropriate compilers and compiler flags, e.g:

```shell
$ ./configure FC=gfortran FCFLAGS="-O2 -pedantic -std=f2008" --prefix=/opt/pkg/neko 
$ make && make install
```

#### Compiling Neko for NVIDIA GPUs
To compile Neko for NVIDIA GPUs
* Make sure you have the CUDA Toolkit installed (e.g. `nvidia-cuda-toolkit`)
* Configure Neko to use CUDA using the `--with-cuda=/path/to/cuda` argument to `configure`, e.g.:
```shell
$ ./configure  --with-cuda=/usr/local/cuda
```
* CUDA compiler flags and options can be passed using `CUDA_CFLAGS`, `CUDA_ARCH` and `NVCC` respectively, e.g:
```shell
$ ./configure  --with-cuda=/usr/local/cuda CUDA_CFLAGS=-O3  CUDA_ARCH=-arch=sm_80 NVCC=/usr/local/cuda/bin/nvcc
```
* Build using `make && make install`    

#### Compiling Neko for AMD GPUs
To compile Neko for AMD GPUs
* Make sure you have the ROCm Toolkit installed
* Configure Neko to use HIP using the `--with-hip=/path/to/hip` argument to `configure`, e.g.:
```shell
$ ./configure  --with-hip=/opt/rocm/hip
```
* HIP compiler flags and options can be passed using HIP_HIPCC_FLAGS and HIPCC, respectively, e.g.:
```shell
$ ./configure  --with-hip=/opt/rocm/hip HIP_HIPCC_FLAGS=-O3  HIPCC=/opt/rocm/hip/bin/hipcc
```

@note More examples, and instructions for specific machines can be found on Neko's [user discussions](https://github.com/ExtremeFLOW/neko/discussions) pages.

## Installing via Spack
Neko is distributed as part of the package manager Spack as `neko`. The package can install releases of Neko as well as the latest commit to the `develop` branch, for most of Neko's supported backends. For a list of all supported variants, see `spack info neko`

### Quick start guide with Spack

To install a CPU build of Neko using Spack, follow the steps below:

``` shell 
$ git clone https://github.com/spack/spack.git
$ cd spack
$ . share/spack/setup-env.sh
$ spack install neko
```
For a GPU build using e.g. CUDA, change the last line to :

``` shell
$ spack install neko+cuda
```

For a more detailed guide on getting started with Spack, please refer to the offical documentation: 
https://spack.readthedocs.io/en/latest/getting_started.html

## Using a Docker container
Perhaps the easiest way to quickly give Neko a try is using a Docker container.
Below we assume that you have Docker up and running on your system. The released
container images can be found here:
https://gitlab.com/ExtremeFLOW/neko/container_registry. Select the image with
the release you want to use. Here we will use version 0.6.1, but in most cases
you will want to simply pick up the latest version available. To the left of
every image there is a button with three dots. Click on it to get the full path
to the image in the correct format for Docker. For our release this is
`registry.gitlab.com/extremeflow/neko/release-0.6.1-ubunut20.04-x86_64-gcc-12.3`.
To get the image on you machine use `docker pull`:

```
docker pull registry.gitlab.com/extremeflow/neko/release-0.6.1-ubunut20.04-x86_64-gcc-12.3
```

Now, let's verify that the image has been added using `docker image ls`. There
should be a row of the following kind in the output.

```
registry.gitlab.com/extremeflow/neko/release-0.6.1-ubunut20.04-x86_64-gcc-12.3   latest    6a9febfaa645   3 months ago    2.71GB
```

The third column contains the ID of the image. We will need that to run Neko in
the container. The typical scenario is that you want to run a case stored on
your computer inside the container. For that we will need to mount the directory
with the case to the container file system. This is done using the the `-v` flag
to the `docker run` command. For example, we will consider that the case resides
in `/home/user/case` and we will mount it to `/case` inside the container. The
full command to execute is the following:

```
docker run --rm -v /home/user/case:/case 6a9febfaa645 bin/bash -c "cd /case && mpirun -n 2 neko case_file.case"
```

The `--rm` flag tells Docker to remove the container after the run is finished.
Note that we use the image ID from before as the third argument. As the run
command we simply use `bash`, followed by a sequence of commands to actually
execute the case. The commands are chained using `&&`, so one can easily add
additional steps, for example, running `makeneko`.
