# Testing {#testing}

## pFUnit
Neko uses the software pFUnit for unit testing.
To install the software you can use the following commands, in which you should just set the desired installation path
```
git clone https://github.com/Goddard-Fortran-Ecosystem/pFUnit.git -b v4.4.2

cd pFUnit && mkdir b && cd b && cmake -DCMAKE_INSTALL_PREFIX=/pfunit_install_path .. && make -j$(nproc) && make install
```
You will now have a dirctory called PFUNIT-4.4 in your intall path.

## Configuring Neko
To use pFUnit with Neko, we need to specify its location during the `configure` phase of the build.
The following configuration is given here as an example, where again you need to provide the path
```bash
 ./configure FC=${FC} FCFLAGS="-O2 -pedantic -std=f2008" --with-pfunit=/pfunit_install_path/PFUNIT-4.4
```
Make sure that in the output you see it says `pFUnit ... yes` at some point.

## Running the tests
To run the tests you should execute 
```bash
make check
```
This will first compile the software and the tests, and then execute them.
In the end you should get a nice list of tests that were run and a summary, which hopefully says that all tests have passed!

## Adding a new test
The first step is of course to read the documentation for pFUnit to understand how it works.
It is up to you to learn how to write `.pf` files.
Here, we only cover how to incorporate the test into the Neko build system.

The tests are located in the `tests` folder.
Unlike some other software, the structure of the `tests` folder does not follow that of `src`.
Instead each folder loosely corresponds to a single class or sometimes a single abstract class and its children (see the `stack` folder, for example).

We now go through the steps of adding a test.
The instructions will differ somewhat depending on whether you test uses MPI or not.

1. In Add a new folder in `tests`, corresponding to the classes that will be tested.
2. Put the `.pf` files you've written inside that folder.
3. Copy over a `Makefile.in` file from some over test folder.
   If your test does not use MPI, you can copy from e.g. `hex`, otherwise take it from `field`. As usual, the `Makefile.in` will turn into a `Makefile` during the `configure` phase of the build.
4. The `Makefile.in` looks very much like the `Makefile` examples in the documentation of pFUnit.
   1. The `check` target should be changed to the name of your test.
      If yor test does not use MPI, name it the same as the name of the folder  + `_test`. If you need MPI, then don't use `test`, but rather something else, like `_suite`.
   2. The lines defining `_TESTS`, `_OTHER_LIBRARIES`, and the other pFUnit stuff should be prefixed with the name selected above. Similarly with the `eval` statement.
   3. In the `clean` target, also put the same name as in `check`.
5. If you use MPI, you now need to create a runner script.
   Copy `field/field_test` tp start with, which looks like this
   ```bash
   #!/bin/sh
   if which mpirun >/dev/null; then
       mpirun -np 1 ./field/field_suite
   else
       mpiexec -np 1 ./field/field_suite
   fi
   ```
   You need to change `./field/field_suite` to the name of the directory for your test and the name of the `check` target in the `Makefile.in`.
6. Now you have to edit `tests/Makefile.am`. This file contains three lists.
   1. The first one is `SUBDIRS`, to which you should add the directory with your tests.
   2. The second one is `TESTS`, here you should add the file, which ends with `_test`.
      If you followed the conventions above it will be the output from pFUnit for a simple test, or the runner you created in the case of an MPI test.
   3. The last list is `EXTRA_DIST`.
      Here, you should add all your `.pf` files.
      That is it for a simple test, but if you use MPI you should also copy the entry you made in `TESTS` to this list.
      Again, you can use `field_test` as an example.
7. Finally, open the `configure.ac` file in the root folder of `neko`.
   Find `# Config tests` line and the `AC_CONFIG_FILES` list below it.
   Add the path to a `Makefile` (note, no `.in`!!) in your test folder. 
8. Please also add the compilation products of your test to the  `tests/.gitignore` file.
   This helps keeps the version controlled file list unpolluted.