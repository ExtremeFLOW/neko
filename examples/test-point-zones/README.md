# This will be deleted once the PR will be merged

A small test case to showcase point zones.

```bash
makeneko user.f90
./neko run.case
paraview --state=state.py
```

In the case file we select 2 shapes, a box and a sphere. These zones are then used in the user file, where we set the temperature to 10 to be able to visualize them afterwards.

A state file for paraview is provided.
