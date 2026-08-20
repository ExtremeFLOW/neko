# Södermalm Cylinder Wind Smoke Case

This folder contains the local Södermalm immersed-boundary wind example used for
quick CPU tests. The mesh is a cylindrical extruded hex mesh with a northwest
quarter inlet arc, terrain-following bottom, open top, and building geometry
provided as an STL for the Brinkman source term.

## Generate Geometry

```sh
python3 build_sodermalm_cylinder_2x_xy.py
```

The script reads prepared GIS cutouts from `input/prepared`, writes the mesh and
building STL to `generated_2x_xy`, produces diagnostic figures in
`figures_2x_xy`, and runs `mesh_checker` on the generated `.nmsh`.

Current mesh target:

- land element size: `110 m`, giving about `55 m` effective p2 x-y spacing
- water element size: `225 m`
- vertical elements: `5`
- domain top: `350 m` above the lowest terrain/water level
- polynomial order in the case: `2`

## Run Neko

```sh
cd run_sodermalm_wind_2x_xy_t10_np6_soft_buildings
./run_neko.sh
```

Optional overrides:

```sh
NP=4 MPIEXEC=mpirun NEKO_BIN=/path/to/neko ./run_neko.sh
```

The case runs to `t = 10`, starts with a uniform northwest-to-southeast wind,
uses a soft Brinkman building mask, and writes the final field under `fields`.

## Render The Heatmap

```sh
cd run_sodermalm_wind_2x_xy_t10_np6_soft_buildings
python3 render_terrain_plus10_heatmap.py
```

The renderer samples the final velocity field at `terrain + 10 m`, overlays the
Södermalm shoreline, buildings, and inlet arc, and writes
`renders/sodermalm_velocity_terrain_plus10.png`.
