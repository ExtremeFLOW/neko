# Neko API {#api}

For some use cases and workflows, more flexibility is required than what can be offered through Neko's user interface or library. To accommodate these usage patterns, Neko provides a C-based API that enables simple integration of Neko into custom workflows and interactive computing environments such as Jupyter notebooks. On top of the API, Neko provides interfaces to different high-level languages, including Python (@ref pyneko) and Julia ([Neko.jl](https://gitlab.com/ExtremeFLOW/Neko.jl)).

Both the API (neko_api) and high-level interfaces provide functionality for defining cases, setting up and controlling simulation workflows, creating callbacks for Neko (e.g., initial and inflow conditions), and methods for direct (zero-copy) access to the underlying data inside Neko.

## Installation

The API is always compiled and installed together with Neko; no special configuration options are necessary to use the C-based API. An example on how to use the C interface is shown in [examples/api/c/cylinder.c](@ref examples/api/c/cylinder.c).

### Python (pyneko)

To install the python interface, configure and install Neko with the
flags `--enable-shared` and `--enable-python`

Finally set `PYTHONPATH` and library paths (e.g. `LD_LIBRARY_PATH`) to
the correct paths under installation prefix. On some systems you might
need to preload (e.g. `LD_PRELOAD`) the MPI library to avoid missing
symbols.

An example on how to use the Python interface is shown in [examples/api/pyneko/cylinder.py](@ref examples/api/pyneko/cylinder.py).

### Julia (Neko.jl)

To interface Neko with Julia, first configure and install Neko with
shared libraries enabled `--enable-shared`. Then install the Neko.jl
packge using:

```julia
(@v1.11) pkg> add https://gitlab.com/ExtremeFLOW/Neko.jl
```

Finally, set library paths (e.g., `LD_LIBRARY_PATH`) to the correct paths under the installation prefix. On some systems, you might need to preload (e.g., `LD_PRELOAD`) the MPI library to avoid missing symbols.

An example on how to use the Julia interface is shown in [examples/api/Neko.jl/cylinder.jl](@ref examples/api/Neko.jl/cylinder.jl).
