# Installing Neko {#installation}

## Dependencies

### Building JSON Fortran 

Download and compile, at least  version 0.7.1 of JSON Fortran from the main repository.
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

### Building gslib (optional)

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
export GSLIB=$(pwd)/build/lib
```

Later, when configuring Neko, add the following option to enable gslib

``` shell
 --with-gslib=${GSLIB}
```

Make sure you see the following message during the configuration:

``` shell
checking for fgslib_gs_setup in -lgs... yes
```

### Building ParMETIS (optional)

@todo Needs to be written

### Bulding PFunit (optional)

To build the PFunit testing framework, please refers to the \subpage testing page
