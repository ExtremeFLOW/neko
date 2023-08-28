# Install and compile gslib

In a folder outside of Neko:

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

Now add the path to `libgs.a` to an environment variable `GSLIB`

``` shell
export GSLIB=$(pwd)/build/lib
```

# Recompile Neko with gslib

In your Neko installation folder, first make sure you use the correct branch

``` shell
git checkout feature/probes
```

Second, reconfigure Neko with `gslib`

``` shell
./regen.sh
./configure --prefix=/path/to/installation --with-gslib=$GSLIB
```

Make sure you see the following message during the configuration:

``` shell
checking for fgslib_gs_setup in -lgs... yes
```

And then recompile Neko:

``` shell
make install
```
