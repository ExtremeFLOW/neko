 ```
                 _  __ ____ __ __ ____ 
     ___  __ __ / |/ // __// //_// __ \
    / _ \/ // //    // _/ / ,<  / /_/ /
   / .__/\_, //_/|_//___//_/|_| \____/ 
  /_/   /___/                         

  A Python interface to Neko
```  
## About
pyNeko allows you to set up and execute simulation cases in Neko
from Python. The interface defines a simulation case as a JSON
object enabling easy integration in contemporary computing
environments.

### Instalation
To build pyNeko you will need a recent version of Neko, Python
and json-fortran. Most dependencies are found via `pkg-config`. 
Thus installation is as simple as:
```bash
./regen.sh (if configure is not present)
./configure --prefix=/path/to/pyneko_installation
make install
```
Finally set `PYTHONPATH` and library paths (e.g. `LD_LIBRARY_PATH`)
to the correct paths under installation prefix.

