# Point zones {#point-zones}

\tableofcontents

## What are point zones?

Point zones are subsections of the computational domain which are 
selected based on a given geometrical criterion. A point zone is 
defined by the `point_zone_t` abstract type. Each `point_zone_t` object has
a unique `name` attribute, and a `mask` containing a list of linear indices
referring to the GLL points whose coordinates verify the above-mentioned
geometrical criterion. Zones can then be used for different purposes, an example 
being applying a localized source term or probing a particular zone of interest.

## Predefined geometrical shapes

There are three predefined shapes from which to initialize a point zone in the
case file: boxes, spheres and cylinders. Each shape is described by its own
subtype `box_point_zone_t`, `sphere_point_zone_t` and `cylinder_point_zone_t`,
extending the abstract class `point_zone_t`.

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

### Cylinder

A cylinder is defined by its end points and its radius.

~~~~~~~~~~~~~~~{.json}
[
    {
        "name": "mycylinder",
        "geometry": "cylinder",
        "start": [0.0, 0.0, 0.0],
        "end": [0.0, 0.0, 1.0],
        "radius": 0.01
    },
]
~~~~~~~~~~~~~~~

## Operations on point zones

### Inversion

By setting the keyword `invert` to `true`, you can indicate that you want to 
select all the points **but** those inside the point zone.

~~~~~~~~~~~~~~~{.json}
[
    {
        "name": "mycylinder",
        "geometry": "cylinder",
        "start": [0.0, 0.0, 0.0],
        "end": [0.0, 0.0, 1.0],
        "radius": 0.01
        "invert": true,
    },
]
~~~~~~~~~~~~~~~

### Combination

It is possible to combine predefined shapes into a single point zone. This is
achieved by the `combine_point_zone_t` type which can construct a point zone based 
on the geometrical parameters of each sub point zone or refer to already existing 
point zones. The mask of a combined point zone will not contain any duplicate 
points in case the subset point zones overlap.

Specifying which point zones to combine is done through the `subsets` array. 
In that array, you can specify:
- a reference to an already existing point zone, or
- an internal point zone.

@note
The first method is useful for creating a new
point zone out of a combination of already existing point zones. If you
are referencing zones `p1` and `p2`, you will get a new point zone `p3`. All 
three of those will be accessible in the registry. 
The second method is useful if you want to create only one point zone using the
geometric criteria of the internal ones. This means that if you create zone 
`p3` out of internal zones `p1` and `p2`, only `p3` will be created. 
This is useful if you do not explicitly need `p1` and `p2`, which will not
be accessible in the registry. Use this method if you are using very large 
point zones and want to save memory.

In the code snippet below, the combine point zone named 
`"combined_box_and_sphere"` uses a subset of two zones: 
- `"mysphere"` is a reference to the external point zone above, and
- `"mybox"` is built internally. 
This means that `"mysphere"` will be accessible in the registry,
whereas `"mybox"` will not.

The combination will be performed according to a specified operator, which can
take the values `"OR"`, `"AND"` or `"XOR"`. 
For example, if you try to combine two point zones `p1` and `p2`:
- The operator `OR` will select all the points that are both in `p1` 
**or** `p2`,
- The operator `AND` will select all the points that are both in `p1`
**and** `p2`,
- The operator `XOR` (exclusive OR) will select all the points that are **either** in 
`p1` **or** in `p2`, but not both at the same time.

@note You can also invert a `combine` point zone.
@attention You can combine other `combine` point zones as long as they are in the
correct order in the case file, as demonstrated in the code snippet below.

~~~~~~~~~~~~~~~{.json}
[
    {
        "name": "mysphere",
        "geometry": "sphere",
        "center": [0.0, 0.0, 0.0],
        "radius": 0.01
    },

    {
        "name": "combined_box_and_sphere",
        "geometry": "combine",
        "operator": "OR",
        "subsets":
        [
            {
                "name": "mysphere"
            },
            {
                "name": "mybox",
                "geometry": "box",
                "x_bounds": [-1.0, 1.0],
                "y_bounds": [-1.0, 1.0],
                "z_bounds": [-1.0, 1.0]
            },

        ],
    },
    {
        "name": "inverted_combination",
        "geometry": "combine",
        "operator": "OR",
        "invert": true,
        "subsets":
        [
            {"name": "combined_box_and_sphere"},
        ],
    },

]
~~~~~~~~~~~~~~~


## User-defined geometrical shapes

The current version of Neko does not support user-defined shapes from the case 
file. That said, shapes can be defined manually into new types by extending 
`point_zone_t` and implementing the abstract `criterion` interface.


## Using point zones {#point-zones_using-point-zones}

Point zones defined in the case file are stored in a point zone registry, 
`neko_point_zone_registry`. The point zone registry allows for the retrieval of
any `point_zone_t` object when needed. Once a `point_zone_t` object is 
retrieved, it can be used for e.g. applying a source term to a localized zone, as demonstrated below:

```fortan
  subroutine forcing(f,t)
    class(fluid_user_source_term_t), intent(inout) :: f
    real(kind=rp), intent(in) :: t

    integer :: i
    class(point_zone_t), pointer :: my_point_zone
    
    my_point_zone => neko_point_zone_registry%get_point_zone("myzone")

    ! Assign a constant forcing to my_point_zone
    do i = 1, my_point_zone%size
       f%u(my_point_zone%mask(i), 1, 1, 1) = 2.0
    end do

  end subroutine forcing
```
