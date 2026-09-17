# Immersed Urban STL

The maintained example is [Sodermalm cylinder](sodermalm_cylinder/README.md):
a terrain-following cylindrical hex mesh and embedded buildings represented
by a cached Brinkman mask.

The retained baseline is the successful p9 continuation to t = 500, with a
sharp raw mask, 10 m PDE mask filtering, HPFRT, and an actual 180-degree inlet.
The original mesh must be reused for checkpoint restarts. Geometry generation,
case preparation, execution, and the continuous velocity heatmap/movie
workflow are documented in the cylinder directory.

GIS inputs, generated artifacts, scheduler scripts, and private production
notes are not versioned.
