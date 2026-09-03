# Euler IDP scalar problems

This collection contains scalar-like transport benchmarks for Neko's
experimental Euler invariant-domain-preserving (IDP) discretization. Each case
embeds a scalar profile in the density field and uses the compressible Euler
solver to exercise the IDP graph viscosity and limiter.

For these proxy problems, the user entropy-pair callback also supplies the
scalar wave-speed field used by the low-order graph scheme. The cases disable
the Euler internal-energy and entropy constraints on the antidiffusive
correction, so density is limited according to the scalar maximum principle.
The three-body rotation case uses the relaxed density bounds from the reference
benchmark; the other cases retain strict bounds.

The cases are numerical diagnostics rather than general scalar-solver examples.
Unless a case states otherwise, they use an affine, fully periodic mesh and the
CPU Euler IDP implementation.

Available cases:

- `burgers`: Two-dimensional Burgers equation with a sonic point, represented
  by the shifted density `rho = u + 1`.
- `kpp`: Nonlinear Kurganov-Petrova-Popov rotating-wave problem with
  flux `(sin(u), cos(u))` embedded in the Euler density equation.
- `three_body_rotation`: Rigid rotation of a slotted disk, cosine hump, and
  cone density profile through one complete revolution.

Additional scalar benchmarks can be added as subdirectories of this
collection.
