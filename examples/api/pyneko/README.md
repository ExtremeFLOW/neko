 ```
                 _  __ ____ __ __ ____ 
     ___  __ __ / |/ // __// //_// __ \
    / _ \/ // //    // _/ / ,<  / /_/ /
   / .__/\_, //_/|_//___//_/|_| \____/ 
  /_/   /___/                         

  A Python interface to Neko
```  
## About
pyneko allows you to set up and execute simulation cases in Neko
from Python. The interface defines a simulation case as a JSON
object enabling easy integration in contemporary computing
environments.

### Instalation

To install the python interface, configure and install Neko with the
flags `--enable-shared` and `--enable-python`

Finally set `PYTHONPATH` and library paths (e.g. `LD_LIBRARY_PATH`) to
the correct paths under installation prefix. On some systems you might
need to preload (e.g. `LD_PRELOAD`) the MPI library to avoid missing
symbols.

