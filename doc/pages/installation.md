# Installing Neko {#installation}

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

### Building ParMETIS (optional)

The following steps is an example on how to build and install ParMETIS

``` shell
$ wget http://glaros.dtc.umn.edu/gkhome/fetch/sw/parmetis/parmetis-4.0.3.tar.gz
$ tar xzf parmetis-4.0.3.tar.gz && cd parmetis-4.0.3 && make config prefix=/parmetis_install_path
$ make -j$(nproc) && make install
```

### Bulding PFunit (optional)

To build the PFunit testing framework, please refers to the \subpage testing page
