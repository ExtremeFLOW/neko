# Installing Neko {#installation}

\tableofcontents

Neko can be installed in various ways, either building directly from source,
manually compiling all dependencies and Neko or via tools like Spack, pixi, and
Docker. What to use comes down to personal preference, but here are some rough
guidelines.

- Use pixi to quickly obtain a CPU build with all optional dependencies in an
  isolated environment.
- Spack is a more advanced package manager, targeting HPC environments. It gives
  you a lot of control yet the possibility to easily get neko for several
  compute backends.
- Use Docker if you are a fan of containers.
- Build from source if you want full control over the build parameters and
  environment and don't want to use a package manager.

## Building from source

To build Neko, you will need a Fortran compiler supporting the Fortran-08 standard, autotools, libtool, pkg-config, a working MPI installation supporting the Fortran 2008 bindings (`mpi_f08`), BLAS/LAPACK and JSON-Fortran. Optional dependencies are PFunit, HDF5 and ParMETIS.

Follow the steps below to install the less common dependencies (e.g. JSON-Fortran).

### Dependencies

#### Building JSON Fortran

Download and compile, at least version 0.7.1 of JSON Fortran from the main repository.
@note Neko requires JSON Fortran to be configured with `USE_GNU_INSTALL_CONVENTION`.

```shell
git clone --depth=1 https://github.com/jacobwilliams/json-fortran.git
cd json-fortran && mkdir b && cd b
cmake -DCMAKE_INSTALL_PREFIX=/path/to/installation -DUSE_GNU_INSTALL_CONVENTION=ON ..
make install
```
Now ad the installation path to `PKG_CONFIG_PATH` (and if needed `LD_LIBRARY_PATH`).
@note On certain systems `lib` should be substituted with `lib64`

```bash
export PKG_CONFIG_PATH=/path/to/installation/lib/pkgconfig:$PKG_CONFIG_PATH
export LD_LIBRARY_PATH=/path/to/installation/lib:$LD_LIBRARY_PATH
```

#### Building HDF5 (optional, but highly recommended)

Download and compile at least version 1.14 of HDF5, with Fortran and MPI
support. Note that you may need to adjust the compilers depending on you machine
(for example, on Cray supercomputers `mpifort` is usually replaced with `ftn`,
and so on).

```shell
wget https://github.com/HDFGroup/hdf5/archive/refs/tags/hdf5_1.14.6.tar.gz
tar xvf hdf5_1.14.6.tar.gz
cd hdf5-hdf5_1.14.6/
cmake -B build -S ./ --install-prefix /path/to/installation \
    -DCMAKE_C_COMPILER=mpicc -DCMAKE_CXX_COMPILER=mpicxx \
    -DCMAKE_Fortran_COMPILER=mpifort -DHDF5_ENABLE_PARALLEL=ON \
    -DHDF5_BUILD_FORTRAN=ON -DCMAKE_BUILD_TYPE=Release
cmake --build build/ --parallel
cmake --install build/
```

It can be a good idea to double-check that you have files starting with
`libhdf5_fortran` in the `lib` directory in your install path. If not, try to
inspecting the log in the `build` directory with `cat CMakeCache.txt | grep
FORTRAN` to identify issues.

Similar to `json-fortran`, populate relevant environmental variables.
```bash
export PATH=:/path/to/installation/bin:$PATH
export PKG_CONFIG_PATH=/path/to/installation/lib/pkgconfig:$PKG_CONFIG_PATH
export LD_LIBRARY_PATH=/path/to/installation/lib:$LD_LIBRARY_PATH
```

#### Building ParMETIS (optional)

The following steps is an example on how to build and install ParMETIS

``` shell
$ export PARMETIS_INSTALL=/usr/local/parmetis
$ wget https://github.com/mfem/tpls/raw/refs/heads/gh-pages/parmetis-4.0.3.tar.gz
$ tar xzf parmetis-4.0.3.tar.gz && cd parmetis-4.0.3 && make config prefix=${PARMETIS_INSTALL} && make install
$ cd metis && make config prefix=${PARMETIS_INSTALL} && make install && cd ../..
```

#### Bulding PFunit (optional)

To build the PFunit testing framework, please refers to the \subpage testing page

#### Python (optional) {#deps-python}
Neko itself does not depend on Python, however, it is necessary for the following:

- Running unit tests; only depends on an interpreter.
- Running integration tests; requires `pytest`.
- Using the `lint.sh` and `format.sh` tools under `contrib`; requires `findent`
  and `flinter==0.4.0`.
- Running Python post-processing scripts for selected examples. A dependency of
  note is [`pySEMtools`](https://github.com/ExtremeFLOW/pySEMTools), which is
  the package we recommend for post-processing Neko simulation results.

### Building Neko {#building-neko}
Neko uses autotools as its build system. The first step is to run the `configure` script, located in the top directory.

``` shell
<path-to-neko>/configure FC=<Fortran compiler> CC=<C compiler> \
                         MPIFC=<MPI Fortran compiler> MPICC=<MPI C ompiler> \
                         FCFLAGS=<Fortran compiler flags> CFLAGS=<C compiler flags> \
                         --prefix=<installation path> [options]
```

In the above command, `[options]` refers to either optional features or packages.

Features are enabled and disabled by passing either `--enable-FEATURE[=arg]` or `--disable-FEATURE` to `configure`. A list of currently supported features are given in the table below.

| Name                  | Description                                                                                                       |
| --------------------- | ----------------------------------------------------------------------------------------------------------------- |
| `--enable-real=Xp`    | Specify working precision of REAL types:<br>`sp` -- `REAL(kind=REAL32)` <br>`dp` -- `REAL(kind=REAL64)` (default) |
| `--enable-contrib`    | Compile various tools                                                                                             |
| `--enable-device-mpi` | Enable device aware MPI                                                                                           |
| `--enable-openmp`     | Enable OpenMP                                                                                                     |
| `--enable-shared`     | Build shared libraries (default: no)                                                                              |
| `--enable-static`     | Build static libraries (default: yes)                                                                             |

When configuring Neko with `sp` precision some variables are still stored in
double precision, for example when gathering data from multiple processes.
This can be avoided by specifying `ssp` instead, however, this is not actively maintained.
Optional packages are controlled by passing either `--with-PACKAGE[=ARG]` or `--without-PACKAGE` to `configure`. A list of all supported optional packages are given in the table below.

| Name                            | Description                                   |
| ------------------------------- | --------------------------------------------- |
| `--with-blas=<lib>`             | Use BLAS library `<lib>`                      |
| `--with-lapack=<lib>`           | Use LAPACK library `<lib>`                    |
| `--with-metis=DIR`              | Directory for metis                           |
| `--with-metis-libdir=LIBDIR`    | Directory for metis library (if different)    |
| `--with-parmetis=DIR`           | Compile with support for parmetis library     |
| `--with-parmetis-libdir=LIBDIR` | Directory for parmetis library (if different) |
| `--with-adios2=DIR`             | Compile with support for ADIOS2               |
| `--with-adios2-fortran=DIR`     | Compile with support for ADIOS2 with Fortran  |              
| `--with-libxsmm`                | Compile with support for libxsmm              |
| `--with-hip=DIR`                | Compile with HIP backend                      |
| `--with-cuda=DIR`               | Compile with CUDA backend                     |
| `--with-opencl=DIR`             | Compile with OpenCL backend                   |
| `--with-metal`                  | Compile with Metal backend (macOS only)       |
| `--with-nvtx=DIR`               | Compile with support for NVTX                 |
| `--with-roctx=DIR`              | Compile with support for ROCTX                |
| `--with-nccl=DIR`               | Compiler with support for NCCL                |
| `--with-rccl=DIR`               | Compiler with support for RCCL                |
| `--with-utofu=DIR`              | Compile with support for the native Tofu interconnect (uTofu) |
| `--with-hdf5`                   | Compile with support for HDF5                 |
| `--with-pfunit=DIR`             | Directory for pFUnit (see \subpage testing)   |

@note Accelerators backends are not enabled as a feature in Neko, but rather via optional packages.

Once configured, to compile and install Neko issue `make` followed by `make install`

#### Compute backends and support tiers {#compute-backends}

Neko's compute kernels are provided by several interchangeable backends,
selected at configure time. The CPU backend is always available; the remaining
backends are enabled through their respective optional packages (see the table
above). The backends share a common interface, so a case runs unmodified
regardless of the backend it is built against. Their current level of support
is summarised below.

Support tiers are defined as follows:

- **Tier 0** --- Primary, fully supported backends. Regularly tested and
  expected to support all features.
- **Tier 1** --- Supported and maintained, but tested less extensively and may
  trail Tier 0 in features or performance.
- **Tier 2** --- Experimental or optional. Under active development or provided
  as an optional optimisation; feature coverage may be incomplete.

| Backend   | Target hardware                      | Tier   | Notes                                       |
| --------- | ------------------------------------ | ------ | ------------------------------------------- |
| CPU       | Any CPU (portable Fortran)           | Tier 0 | Always available; reference implementation  |
| CUDA      | NVIDIA GPUs                          | Tier 0 |                                             |
| HIP       | AMD GPUs                            | Tier 0 |                                             |
| OpenCL    | Cross-vendor GPUs and accelerators   | Tier 1 |                                             |
| SX-Aurora | NEC SX-Aurora TSUBASA vector engines | Tier 1 |                                             |
| Metal     | Apple Silicon GPUs (macOS)           | Tier 2 | Single precision only (`--enable-real=sp`)  |
| libxsmm   | x86 CPUs (JIT small matrix multiply) | Tier 2 | Accelerates the CPU backend                 |

#### Communication backends {#communication-backends}

Neko's gather--scatter layer, which performs the halo exchange between elements
and MPI ranks, can use several communication backends. The backend is selected
automatically at runtime (device-aware MPI when available, otherwise host MPI),
but can be overridden through the `NEKO_GS_COMM` environment variable. The
supported backends are listed below.

| Backend          | Mechanism / target                              | Compute backends | Enable (configure)            | `NEKO_GS_COMM` |
| ---------------- | ----------------------------------------------- | ---------------- | ----------------------------- | -------------- |
| MPI              | Standard MPI on host buffers (default)          | All backends     | always available              | `MPI`          |
| Device-aware MPI | MPI directly on device buffers (GPUDirect)      | CUDA, HIP        | `--enable-device-mpi`         | `MPIGPU`       |
| NCCL / RCCL      | Collective comms library on GPUs (NVIDIA / AMD) | CUDA, HIP        | `--with-nccl` / `--with-rccl` | `NCCL`         |
| NVSHMEM          | NVIDIA OpenSHMEM, one-sided on device buffers   | CUDA             | `--with-nvshmem`              | `SHMEM`        |
| OpenSHMEM        | Host OpenSHMEM, one-sided                        | CPU              | `--with-openshmem`           | `SHMEM`        |
| Co-Array Fortran | Fortran coarrays                                | CPU              | coarray-capable compiler      | `CAF`          |
| uTofu            | Native Tofu interconnect, one-sided RDMA        | CPU              | `--with-utofu`                | `UTOFU`        |

@note Every device backend (CUDA, HIP, OpenCL, Metal) can use host `MPI`, which
stages data through host memory before the exchange. The device-resident
backends --- device-aware MPI, NCCL/RCCL and NVSHMEM --- are implemented only
for CUDA and HIP (NVSHMEM is CUDA-only), so OpenCL and Metal builds fall back to
host `MPI`.

@note `NEKO_GS_COMM=SHMEM` selects NVSHMEM on a device (GPU) build and host
OpenSHMEM otherwise.

@note `UTOFU` targets the native Tofu interconnect (Tofu-D, e.g. Fugaku) and
requires the `libtofucom` library; it is a host (CPU) backend. The number of
Tofu network interfaces used for injection is set with `NEKO_GS_UTOFU_NTNI`
(default `1`).

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
* If `--enable-device-mpi` is also given, `CUDA_ARCH` is **required**:
  `configure` will error out rather than guess an architecture, since
  targeting the wrong one silently hurts performance instead of failing to
  build, and Neko is often cross-compiled on a host with no GPU present to
  detect against. Find the right value for the target GPU with e.g.
  `nvidia-smi --query-gpu=compute_cap --format=csv,noheader` on the target
  machine, or the vendor's published compute capability for the target
  cluster's GPUs when cross-compiling.
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

@note On APUs with unified physical memory (e.g. AMD Instinct MI300A) and XNACK enabled (`HSA_XNACK=1`), Neko can map arrays zero-copy: host and device share a single allocation instead of keeping replicated copies, and host-device transfers become no-ops. This roughly halves the memory footprint of mapped data, but is currently slower per step than replicated buffers on MI300A, so it is **off by default** and is intended as an opt-in capacity mode for running larger cases. Enable it at runtime by setting the environment variable `NEKO_HIP_ZEROCOPY=1`; it then activates only on a supported APU with XNACK (discrete GPUs such as MI250X always use replicated buffers, and requesting it there prints a warning and falls back). Since host and device share one allocation under zero-copy, host code (e.g. user routines) must not write to a mapped array while device work touching it may be in flight; synchronize with `device_sync` first.

#### Compiling Neko for Apple Silicon GPUs
To compile Neko for Apple Silicon GPUs (macOS only)
* Make sure you have Xcode with the Metal shader compiler installed (since Xcode 26 the Metal toolchain is a separate download, install via `xcodebuild -downloadComponent MetalToolchain`)
* Configure Neko to use Metal using the `--with-metal` argument to `configure`, e.g.:
```shell
$ ./configure --with-metal --enable-real=sp
```
* The Objective-C compiler and flags can be passed using `OBJC` and `OBJCFLAGS`, respectively (defaults to `clang`)
* Build using `make && make install`

@note Apple GPUs do not support double precision; the Metal backend requires Neko to be configured with single precision (`--enable-real=sp`).

@note On devices with unified memory (Apple Silicon), the Metal backend maps arrays zero-copy: host and device share a single allocation instead of keeping replicated copies, and host-device transfers become no-ops. This roughly halves the memory footprint of mapped data. Zero-copy mapping can be disabled at runtime by setting the environment variable `NEKO_METAL_ZEROCOPY=0`, which restores fully replicated buffers. On GPUs without unified memory (e.g. AMD GPUs in Intel-based Macs), buffers are always replicated.

@note More examples, and instructions for specific machines can be found on Neko's [user discussions](https://github.com/ExtremeFLOW/neko/discussions) pages.

#### Compiling Neko with a collective communications library
To compile Neko to use a collective communcations library on GPUs
* Configure Neko to use NCCL (on NVIDIA GPUs) using the `--with-nccl=/path/to/nccl` argument to `configure`
* Configure Neko to use RCCL (on AMD GPUs) using the `--with-rccl=/path/to/rccl` argument to `configure`

@note This will change all collective communications (reductions etc) during a simulation step to use NCCL/RCCL, while gather-scatter operations still uses MPI.

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

For a more detailed guide on getting started with Spack, please refer to the
offical documentation:
https://spack.readthedocs.io/en/latest/getting_started.html

## Installing using pixi
Pixi is a package managment tool using conda under the hood. It is very easy
to install. For more information see (pixi.sh)[pixi.sh].

```bash
curl -fsSL https://pixi.sh/install.sh | sh
```

Pixi will leverage conda to install all the dependencies, including basic ones
like `gfortran` and `openmpi`. All of these will be installed inside an isolated
environment. So, to install Neko simply clone the repo with git,
and run the following command inside it

```bash
pixi run install-neko-cpu
```

This will give you a double-precision CPU build charged with optional
dependencies: HDF5, pFUnit, and ParMETIS. You can, however choose the precision
of reals, by appending either `sp` or `dp` to the command above. Only CPU builds
are currently supported with pixi.

To use Neko, you need to drop into a shell, where the pixi environment will be
activated. For that run

```bash
pixi shell
```

The `neko` and `makeneko` executables are already be in your `PATH`, so you can
start running cases!

The installed executables, libraries, etc. are all located inside the `install` folder in the repo.

Note that you can use this pixi environment as you like, including manually
`configuring` and building Neko (as per instructions for building from source),
for example, with single precision reals or even with a different backend.

Pixi also provides the option to create a Python environment, which is just an
ordinary `conda` environment underneath. This can be convenient because Neko's
unit and integration tests, some post-processing in examples, and developer
tools like `flinter` require Python.

Running
```bash
pixi shell -e python
```
will drop into a shell, which, in addition to Neko, has a Python interpreter
with `pytest`, `findent`, and `flinter` installed. Additional packages can be
installed manually.

## Installing via FreeBSD ports
Neko is available in the FreeBSD ports collection under `science/neko`. The
port tracks Neko releases and pulls in all required dependencies (MPI,
BLAS/LAPACK, JSON-Fortran, etc.) automatically.

To install from the ports tree:

```shell
cd /usr/ports/science/neko
make install clean
```

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
