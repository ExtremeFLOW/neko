# Immersed Zones

This set of examples highlights the use of the Brinkman body force to model
immersed zones in a flow. The Brinkman body force is a penalization method
that is used to model the presence of solid objects in a flow. The method
applies a force to the fluid in the vicinity of the solid object to mimic the
effect of the solid object on the flow. The force is proportional to the fluid 
velocity and is directed opposite to the fluid velocity.

## Examples

The following examples are included in this directory:

- `block`: A simple example of a flow past a block placed as a half height wall
  half way through the domain.
- `sphere`: A flow past a sphere with radius 0.1 placed in the middle of the
  domain.
