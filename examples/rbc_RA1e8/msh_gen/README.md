[gmsh](http://gmsh.info//) script for generating 2D/3D mesh for straight pipe 

* For nomenclature and the mesh parameters, see the attached figure
* Choose either 2D (= circular disc) or 3D (= straight pipe) mesh options in GENERAL SETTING
* Set the mesh controlling parameters in SETTINGS
* To generate 3D mesh: 
     gmsh pipeMesh.geo -3 -order 2
* To generate 2D mesh: 
     gmsh pipeMesh.geo -2 -order 2
* For turbulent pipe flow simulation using Nek5000, see this [technical report](http://urn.kb.se/resolve?urn=urn:nbn:se:kth:diva-265021).
