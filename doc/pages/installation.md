# Installing Neko {#installation}

Neko can be installed in various ways, either building directly from source, manually compiling all dependencies and Neko or via Spack. Pre-built Docker images are also provided for each release of Neko.

## Building from source

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

@todo Needs to be written

#### Bulding PFunit (optional)

To build the PFunit testing framework, please refers to the \subpage testing page

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
