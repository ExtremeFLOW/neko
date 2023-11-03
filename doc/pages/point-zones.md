
# Point zones {#point-zones}
## What are point zones?

Point zones are subsections of the computational domain containing GLL points
which are selected based on a given geometrical criterion.

These zones can then be used in different contexts, e.g. 
localized source terms, probes...

## Predefined geometrical shapes

There are two predefined shapes from which to initialize a point zone in the case
file: boxes and spheres. Each shape is described by its own subtype
`box_point_zone_t` and `sphere_point_zone_t`, extending the abstract class 
`point_zone_t`.

### Box

A box is defined from its `x,y` and `z` boundaries.

~~~~~~~~~~~~~~~{.json}
[
    {
        "name": mybox,
        "geometry": "box",
        "x_bounds": [-1.0, 1.0],
        "y_bounds": [-1.0, 1.0],
        "z_bounds": [-1.0, 1.0]
    }
]
~~~~~~~~~~~~~~~
### Sphere

A sphere is defined by its center and its radius.

~~~~~~~~~~~~~~~{.json}
[
    {
        "name": "mysphere",
        "geometry": "sphere",
        "center": [0.0, 0.0, 0.0],
        "radius": 0.01
    },
]
~~~~~~~~~~~~~~~

## Using point zones

Any point zone defined in the case file will be stored in a point
zone registry, `neko_point_zone_registry`, from which it can be retrieved.
~~~~~~~~~~~~~~{.f90}
class(point_zone_t), pointer :: my_point_zone
my_point_zone => neko_point_zone_registry%get_point_zone("myzone")
~~~~~~~~~~~~~~
