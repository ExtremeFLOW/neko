# Installation and caching of JSON-Fortran

This action uses caching to accelerate the building of JSON-Fortran. The action
will download the JSON-Fortran source code and build it using the provided
compiler. The action will also cache the build to speed up future builds.

## Inputs

| Name               | Optional | Description                                                                                            | Default                         |
| ------------------ | -------- | ------------------------------------------------------------------------------------------------------ | ------------------------------- |
| `install-dir`      | Yes      | The directory to install JSON-Fortran to.                                                              | `/home/runner/pkg/json-fortran` |
| `working-dir`      | Yes      | The directory to work in.                                                                              | `/home/runner/tmp/json-fortran` |
| `os`               | Yes      | The operating system to use for building JSON-Fortran. Which should allow the use of matrix workflows. | `runner.os`                     |
| `compiler`         | Yes      | The compiler to use for building JSON-Fortran. The compiler should be available in the PATH.           | `gfortran`.                     |
| `compiler-options` | Yes      | The compiler options to use for building JSON-Fortran.                                                 | `-O3`                           |
| `build-options`    | Yes      | The build options to use for building JSON-Fortran.                                                    | `--parallel=$(nproc)`.          |
| `version`          | Yes      | The version of JSON-Fortran to build.                                                                  | `8.3.0`.                        |

## Outputs

| Name          | Description                                     |
| ------------- | ----------------------------------------------- |
| `install-dir` | The directory where JSON-Fortran was installed. |

## Example usage

The following example uses the JSON-Fortran action to build JSON-Fortran using
the `gfortran` compiler with the `-O3` optimization level and the `--parallel 4`
cmake build option.

Additionally the next step will capture the install location and print where the
JSON-Fortran was installed.

```yaml
name: Build JSON-Fortran

on: [push, pull_request]

jobs:
  build:
    runs-on: ubuntu-latest
    steps:

    - name: Setup JSON-Fortran
      uses: ./.github/actions/setup_json-fortran
      with:
        compiler: 'gfortran'
        compiler-options: '-O3'
        build-options: '--parallel=4'

    - name: Echo install directory
      env: 
        JSON_FORTRAN_DIR: ${{ steps.setup-json-fortran.outputs.install-dir }}
      run: echo "JSON-Fortran was installed to $JSON_FORTRAN_DIR"
```