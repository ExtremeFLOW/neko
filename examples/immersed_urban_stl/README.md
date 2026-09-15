# Immersed Urban STL

The maintained example is [Sodermalm cylinder](sodermalm_cylinder/README.md):
a terrain-following cylindrical hex mesh and embedded buildings represented
by a cached Brinkman mask.

The current configuration is p7, sharp raw mask, 10 m PDE mask filtering,
HPFRT, and an exactly validated 220-degree inlet. Geometry generation,
case preparation, execution, and the continuous velocity heatmap/movie
workflow are documented in the cylinder directory.

GIS inputs, generated artifacts, scheduler scripts, and private production
notes are not versioned.
