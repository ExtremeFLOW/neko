# Immersed Zones

This set of examples highlights the use of the Brinkman body force to model
immersed zones in a flow. The Brinkman body force is a penalization method
that is used to model the presence of solid objects in a flow. The method
applies a force to the fluid in the vicinity of the solid object to mimic the
effect of the solid object on the flow. The force is proportional to the fluid 
velocity and is directed opposite to the fluid velocity.

## Examples

All the examples will be conducted in a 3D domain with a size of 4x1x1. The
domain is discretized with a uniform grid of 32x8x8 elements. The number of
elements can be modified by using the `-x`, `-y`, and `-z` command line
arguments. The domain have an inflow boundary condition at the surface at
x=0, and an outflow boundary condition at the surface at x=4. All other
boundaries are treated as walls.

The following examples are included in this directory:

- `block`: A simple example of a flow past a block placed as a half height wall
  half way through the domain.
- `sphere`: A flow past a sphere with radius 0.1 placed in the middle of the
  domain.
- `block_sphere`: A flow past a block and a sphere both with size and placement
  as above.
- `cylinder`: A flow past a cylinder with radius 0.1 placed horizontally in the
  middle of the domain.
  