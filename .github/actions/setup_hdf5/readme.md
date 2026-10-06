# Installation and caching of HDF5

This action uses caching to accelerate the building of HDF5. The action
will download the HDF5 source code and build it using the provided
compiler. The action will also cache the build to speed up future builds.

## Inputs

| Name               | Optional | Description                                                                                    | Default                 |
| ------------------ | -------- | ---------------------------------------------------------------------------------------------- | ----------------------- |
| `install-dir`      | Yes      | The directory to install HDF5 to.                                                              | `/home/runner/pkg/hdf5` |
| `working-dir`      | Yes      | The directory to work in.                                                                      | `/home/runner/tmp/hdf5` |
| `os`               | Yes      | The operating system to use for building HDF5. Which should allow the use of matrix workflows. | `runner.os`             |
| `compiler`         | Yes      | The compiler to use for building HDF5. The compiler should be available in the PATH.           | `gfortran`.             |
| `compiler-options` | Yes      | The compiler options to use for building HDF5.                                                 | `-O3`                   |
| `build-options`    | Yes      | The build options to use for building HDF5.                                                    | `--parallel=$(nproc)`.  |
| `version`          | Yes      | The version of HDF5 to build.                                                                  | `8.3.0`.                |

## Outputs

| Name          | Description                                     |
| ------------- | ----------------------------------------------- |
| `install-dir` | The directory where HDF5 was installed. |

## Example usage

The following example uses the HDF5 action to build HDF5 using
the `gfortran` compiler with the `-O3` optimization level and the `--parallel 4`
cmake build option.

Additionally the next step will capture the install location and print where the
HDF5 was installed.

```yaml
name: Build HDF5

on: [push, pull_request]

jobs:
  build:
    runs-on: ubuntu-latest
    steps:

    - name: Setup HDF5
      uses: ./.github/actions/setup_hdf5
      with:
        compiler: 'gfortran'
        compiler-options: '-O3'
        build-options: '--parallel=4'

    - name: Echo install directory
      env: 
        JSON_FORTRAN_DIR: ${{ steps.setup-hdf5.outputs.install-dir }}
      run: echo "HDF5 was installed to $JSON_FORTRAN_DIR"
```