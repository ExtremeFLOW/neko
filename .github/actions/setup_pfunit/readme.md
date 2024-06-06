# Installation and caching of pFUnit

This action uses caching to accelerate the building of pFUnit. The action will
download the pFUnit source code and build it using the provided compiler. The
action will also cache the build to speed up future builds.

## Inputs

| Name             | Optional | Description                                  | Default                   |
| ---------------- | -------- | -------------------------------------------- | ------------------------- |
| `install-dir`    | Yes      | The directory to install the pFUnit library. | `/home/runner/pkg/pfunit` |
| `working-dir`    | Yes      | The directory to work in.                    | `/home/runner/tmp/pfunit` |
| `os`             | Yes      | The operating system to use.                 | `${{ runner.os }}`        |
| `compiler`       | Yes      | The compiler to use.                         | `mpif90`                  |
| `compiler-flags` | Yes      | The compiler flag to use.                    | `-O3`                     |
| `build-options`  | Yes      | The build option to use.                     | `--parallel $(nproc)`     |
| `version`        | Yes      | The version of the pFUnit library to use.    | `v4.8.0`                  |

## Outputs

| Name          | Description                               |
| ------------- | ----------------------------------------- |
| `install-dir` | The directory where pFUnit was installed. |

## Example usage

The following example uses the pFUnit action to build pFUnit using the `gfortran`
compiler with the `-O3` optimization level and the `--parallel 4` cmake build
option.

Additionally the next step will capture the install location and print where the
pFUnit was installed.

```yaml
name: Build pFUnit

on: [push, pull_request]

jobs:
  build:
    runs-on: ubuntu-latest
    steps:

    - name: Setup pFUnit
      uses: ./.github/actions/setup_neko
      with:
        compiler: 'gfortran'
        compiler-options: '-O3'
        build-options: '--parallel=4'

    - name: Echo install directory
      env: 
        PFUNIT_DIR: ${{ steps.setup-pfunit.outputs.install-dir }}
      run: echo "pfunit was installed to $PFUNIT_DIR"
```