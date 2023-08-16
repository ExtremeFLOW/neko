# UQit-ts
## Time-series analysis module for [UQit](https://github.com/KTH-Nek5000/UQit) (to be released in 2022.)
#
### Saleh Rezaeiravesh, salehr@kth.se

### How to use:
  1. clone the toolbox to your disk
  2. In `~/.bashrc`: `export UQit_ts_path=<path-on-the-disk/>`
  3. `source ~/.bashrc`
  4. In your Python scripts: 
     * `import os`
     * `import sys`
     * `UQit_ts_path=os.getenv("UQit_ts_path")`
     * import a module: e.g.<br/>
       `sys.path.append(UQit_ts_path+'estimators_tUQ/')`<br/>
       `import NOBM`<br/>
#
### Required Libraries:
  - [`numpy`](https://numpy.org/)
  - [`scipy`](https://www.scipy.org/)  
  - [`matplotlib`](https://matplotlib.org/)
  - [`statsmodels`](https://www.statsmodels.org/stable/index.html)
