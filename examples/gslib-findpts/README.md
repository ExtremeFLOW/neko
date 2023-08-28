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

# Running the case

``` shell
makeneko user.f90
./neko user.case
```

## Changing the number of points

At line 44: set `N = ...` with the desired value

## Which errors to expect

The call at line 99 is suspected to mess everything up (`fgslib_findpts`). In the current configuration, it seems to be messing with the `elid` array, which contains the owning element for each point.

Running with one rank should give the following output:

```
Rank #    List of process owners                     List of element owners              error code
0       / 0     0     0     0     0 /         0  962592768     553     554     555 /0 1 0 1 1 
```
In the list of element owners, notice an ugly value at `elid(2)`, whereas `elid(3), elid(4)` and `elid(5)` are showing
correct values since we are running on a 1000 element-mesh.


