# Global Interpolation {#global-interpolation}

## Overview

The `global_interpolation` module in Neko provides functionality for interpolating arbitrary points within a computational domain. It is designed to handle distributed parallel environments, enabling interpolation across multiple processes. This module is particularly useful for applications requiring interpolation of fields or data points in spectral element simulations.

---

## Module: `global_interpolation`

### Description

The `global_interpolation` module implements global interpolation for arbitrary points in the domain. It supports distributed parallelism using MPI and provides tools for finding process owners, element candidates, and local coordinates (`rst`) for interpolation.

### Features
- **Distributed Interpolation**: Handles interpolation across multiple processes using MPI.
- **Flexible Point Management**: Supports arbitrary points in the domain.
- **Element and Process Ownership**: Identifies owning elements and processes for each point.
- **Device Support**: Compatible with GPU-based computations via device memory management.
- **Error Handling**: Includes mechanisms to validate points and ensure interpolation accuracy.

---

## Type: `global_interpolation_t`

### Description

The `global_interpolation_t` type encapsulates the functionality for global interpolation. It manages the data structures, processes, and algorithms required for interpolation in a distributed environment.

### Attributes

#### Domain and Space Information
- `gdim`: Geometric dimension of the simulation.
- `nelv`: Number of elements in the local mesh.
- `glb_nelv`: Global number of elements across all processes.
- `Xh`: Spectral space used for interpolation.

#### Point Management
- `n_points`: Total number of points to interpolate.
- `xyz`: Coordinates of the points in the global domain (shape: `3 x n_points`).
- `rst`: Local `rst` coordinates for interpolation (shape: `3 x n_points`).
- `pe_owner`: List of owning processes for each point.
- `el_owner0`: List of owning elements for each point (0-indexed).

#### Local Points
- `n_points_local`: Number of points local to the current process.
- `xyz_local`: Local coordinates of points (shape: `3 x n_points_local`).
- `rst_local`: Local `rst` coordinates for points (shape: `3 x n_points_local`).
- `el_owner0_local`: Owning elements for local points.

#### Interpolation Tools
- `local_interp`: Instance of `local_interpolator_t` for local interpolation.
- `rst_finder`: Object for finding `rst` coordinates (`legendre_rst_finder_t`).
- `gs_comm`: Gather-scatter communication object for distributed interpolation.

#### Parallelism
- `comm`: MPI communicator for distributed interpolation.
- `pe_rank`: Rank of the current process in `comm`.
- `pe_size`: Total number of processes in `comm`.

#### Configuration
- `tol`: Tolerance for Newton iterations to find `rst` coordinates.
- `padding`: Padding for bounding boxes used in element and process finding.

---

### Methods

#### Initialization
- `init_xyz(x, y, z, gdim, nelv, Xh, comm, tol, pad)`: Initializes the global interpolation object using coordinates and mesh information.
- `init_dof(dofmap, comm, tol, pad)`: Initializes the object using a `dofmap_t` instance.

#### Point Management
- `find_points_coords(x, y, z)`: Finds process owners, elements, and `rst` coordinates for given points.
- `find_points_coords1d(x, y, z)`: Similar to `find_points_coords`, but for 1D arrays.
- `find_points_and_redist()`: Finds points and redistributes them to their respective owners, changes the points this process has.

#### Interpolation
- `evaluate(interp_values, field, on_host)`: Evaluates interpolated values at the points using a given field.

#### Validation
- `check_points(x, y, z)`: Validates the points to ensure interpolation accuracy.

#### Memory Management
- `free()`: Frees all allocated resources.
- `free_points()`: Frees resources related to global points.
- `free_points_local()`: Frees resources related to local points.

---

### Example Usage

#### Initialization
```fortran
type(global_interpolation_t) :: glob_interp
type(space_t) :: Xh
type(dofmap_t) :: dof

! Initialize global interpolation object
call glob_interp%init_dof(dof, NEKO_COMM, tol=1e-6, pad=1e-2)
```

#### Finding Points
```fortran
real(kind=rp), allocatable :: x(:), y(:), z(:)
integer :: n_points

! Allocate and set coordinates
n_points = 100
allocate(x(n_points), y(n_points), z(n_points))
x = random_number_array(n_points)
y = random_number_array(n_points)
z = random_number_array(n_points)

! Find points and their owners
call glob_interp%find_points_coords(x, y, z)
```

#### Interpolation
```fortran
real(kind=rp), allocatable :: interp_values(:)
real(kind=rp), allocatable :: field(:)
logical :: on_host

! Allocate field and interpolation values
allocate(field(glob_interp%nelv * Xh%lxyz))
allocate(interp_values(n_points))

! Perform interpolation
on_host = .true.
call glob_interp%evaluate(interp_values, field, on_host)
```

---

### Notes
- **Device Support**: If `NEKO_BCKND_DEVICE` is enabled, the module uses GPU memory for computations.
- **Error Handling**: Ensure `tol` and `pad` are appropriately set for accurate interpolation.
- **Parallelism**: The module relies on MPI for distributed computations. Ensure MPI is properly initialized.

### Environment variables

This module depends on the environemnt variables `NEKO_GLOBAL_INTERP_PE_FINDER` and `NEKO_GLOBAL_INTERP_EL_FINDER`. Unless these are specified to AABB (using Axis aligned bounding boxes to find PEs and elements), the global interpolation object will make use of a structured cartesian grid for this. This grid is distributed among all processes for finding the PEs and local for each PE to find the element candidates. This is in general a lot faster than the AABB trees, but if you run into issues try specifying these environment variables to AABB.

---

### Related Modules
- `local_interpolation`: Provides local interpolation functionality.
- `legendre_rst_finder`: Finds `rst` coordinates for interpolation.
- `aabb_pe_finder` and `cartesian_pe_finder`: Tools for finding process owners.
- `aabb_el_finder` and `cartesian_el_finder`: Tools for finding element candidates on this process.
