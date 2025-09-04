
# Neko.jl -- A Julia interface to Neko

## About
Neko.jl allows you to set up and execute simulation cases in Neko
from Julia. The interface defines a simulation case as a JSON
object enabling easy integration in contemporary computing
environments.

### Instalation

To interface Neko with Julia, first configure and install Neko with
shared libraries enabled `--enable-shared`. Then install the Neko.jl
packge using:

```julia
(@v1.11) pkg> add https://gitlab.com/ExtremeFLOW/Neko.jl
```

Finally, set library paths (e.g., `LD_LIBRARY_PATH`) to the correct
paths under the installation prefix. On some systems, you might need
to preload (e.g., `LD_PRELOAD`) the MPI library to avoid missing
symbols.
