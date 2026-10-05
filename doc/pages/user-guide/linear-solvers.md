# Linear solvers and preconditioners {#linear-solvers}

\tableofcontents

What to know when choosing the solver, preconditioner and tolerance of the
pressure, velocity and scalar equations. The keys themselves are described
in the [case file](@ref case-file_linear-solver) page.

## Solvers

| `type` | memory (fields of size `n`) | backends | requires |
|---|---|---|---|
| `gmres` (flexible, space size `m`) | `2m + 2` (62 for the default `m = 30`) | all | nothing; the preconditioner may be nonlinear |
| `cg` | 4 | all | SPD operator; symmetric, positive, fixed, linear preconditioner |
| `fused_cg` | 13 | CUDA, HIP | as `cg` |
| `pipecg` | 19 | not OpenCL or Metal | as `cg`; see below |
| `cacg` | 19 (host) | CPU | as `cg` |
| `bicgstab` | 7 | all | fixed, linear preconditioner; two preconditioner applications per iteration |
| `coupled_cg`, `fused_coupled_cg`, `coupled_bicgstab` | | see the case-file page | as their scalar counterparts; velocity with variable viscosity |

The pressure and velocity operators are symmetric positive (semi-)definite,
so every solver above applies to them.

**GMRES space size.** GMRES(`m`) restarts after `m` iterations. A solve that
converges within `m` iterations is unaffected by `m`, so set
`gmres_space_size` to the iteration count you see in steady operation (with
`phmg`: about 5 at a tolerance of 1e-4, 10–15 at 1e-6, 25–30 at 1e-8). With
a multigrid preconditioner, a space half that size costs about 10 % more
iterations. `gmres_space_size: 5` needs 12 fields instead of 62.

**What the residual measures.** All solvers stop when
\f$ (\sum_i m_i r_i^2 / V)^{1/2} < \f$ `absolute_tolerance`, with
\f$ r = f - Ax \f$ assembled and masked, \f$ m_i \f$ the inverse
multiplicity and \f$ V \f$ the volume. The tolerance is absolute: it is
about 2.8 times looser after each uniform refinement by two, and the smooth
pressure error it allows grows with the square of the domain size. The
residual barely sees the smoothest error modes. GMRES minimises the residual
and removes them last; CG minimises the energy norm of the error and removes
them earlier, so at the same tolerance CG leaves a smaller pressure error,
for 5–50 % more iterations. The velocity sees only the remaining divergence,
about a tenth of the tolerance, whatever the solver. If only the velocity
matters, 1e-4 is enough; if pressure statistics matter, use `cg` at 1e-6 or
tighter.

**PipeCG** reports a recurrence, not the true residual. In tests the two
differed by a factor 0.4 to 2.7, and on a badly distorted mesh PipeCG broke
down. Prefer `cg` or `fused_cg`.

## Preconditioners

| `type` | properties | compatible solvers |
|---|---|---|
| `jacobi` | symmetric, positive, fixed, linear | all |
| `phmg`, default `smoother_cheby_acc: "jacobi"` | symmetric, positive, fixed, linear | all |
| `phmg`, `smoother_cheby_acc: "schwarz"` | **not symmetric**, fixed, linear | `gmres`, `bicgstab` |
| `hsmg`, default CG coarse solve | **not linear** (the coarse solve is a Krylov iteration) | `gmres` |
| `hsmg`, `coarse_grid.solver: "tamg"` | **not symmetric**, fixed, linear | `gmres`, `bicgstab` |
| `ident` | identity | all |

A solver with a preconditioner that lacks what it requires rarely fails
outright. It stalls, swings in iteration count, or converges to a wrong
solution.

### The solver check at start-up

Neko measures these properties at initialisation for the pressure, velocity
and scalar solvers and the PDE filter (the velocity only with component-wise
boundary conditions) and prints, for example,

```
   ----Solver check: Pressure----
   Operator    : symmetric ( 1.84E-17), positive ( 0.679)
   phmg        : symmetric ( 5.34E-17), positive ( 8.46E-02)
                 fixed ( 0.00E+00), linear ( 3.30E-16)
   Pairing     : cg + phmg is compatible
```

The numbers are the relative asymmetry, \f$ \langle x, Mx\rangle /
(\|x\|\,\|Mx\|) \f$ for positivity, and relative deviations for fixedness
and linearity. An incompatible pairing gives a warning and the run
continues. The check uses four scratch fields that the time loop needs
anyway and about a dozen preconditioner applications. A positivity value far
below one (about 1e-3) means a nearly indefinite preconditioner, which
happens on very distorted meshes; expect slow convergence.

### Checking the true residual

`residual_check_interval: N` recomputes \f$ \|f - Ax\| \f$ after the solve
every `N` steps (and always at the first), prints it next to the reported
residual and warns when they differ by more than 20 %; differences below a
tenth of the tolerance are ignored. It costs one operator application per
checked solve and catches a drifted recurrence (PipeCG, or GMRES after very
many iterations).

## The multigrid coarse grid

`phmg` coarsens in polynomial order to `p = 1` and solves the coarsest
problem with a few cycles of a tree-AMG. Keys in the `preconditioner`
object:

| key | default | meaning |
|---|---|---|
| `smoother_iterations` | 3 | Chebyshev iterations per level |
| `smoother_cheby_acc` | `"jacobi"` | `"schwarz"` is stronger but not symmetric |
| `pcoarsening_schedule` | `[3, 1]` | polynomial orders of the coarser levels |
| `coarse_grid.levels` | 3 | tree-AMG levels |
| `coarse_grid.iterations` | 1 | tree-AMG cycles per application |
| `coarse_grid.cheby_degree` | 4 | Chebyshev degree in the tree-AMG |

The coarse grid removes the smooth error modes the residual barely sees.
With the defaults, a solve at 1e-4 left 18 % of the pressure increment as
smooth error on a box mesh and 50–94 % on curved, distorted and
unstructured meshes; four levels and three cycles roughly halved this at the
same iteration counts. On large meshes choose `coarse_grid.levels` so that
the coarsest level has at most a few thousand nodes, with two or three
`coarse_grid.iterations`; the coarse problem is at `p = 1` and costs little.
The start-up log prints the estimated largest eigenvalue of each tree-AMG
level; a level-0 value far above 10–20 means stretched or distorted elements
and more iterations. On very distorted meshes (scaled Jacobian below 0.1)
more coarse-grid cycles can make the preconditioner indefinite, seen as
GMRES stagnating at `max_iterations`; use the defaults there.

## Projection

`projection_space_size` starts a solve from the best combination of
previous solutions and helps when iteration counts are large. With a
multigrid preconditioner at loose tolerances it was observed to make the
smooth pressure error much worse (100 to 3000 times after 150 steps at 1e-4,
with no effect on the velocity): the stored solutions carry the smooth
errors of earlier inexact solves, and the projection is as blind to them as
the residual. Use it with caution for the pressure, and not with `phmg` at
1e-4 to 1e-6.

## Recommended settings

Pressure with `phmg`: `cg` (4 fields) for memory, or `gmres` with
`gmres_space_size` at the steady iteration count for speed; projection off;
1e-4 when only the velocity matters, 1e-6 with `cg` when pressure statistics
matter; on large or non-affine meshes strengthen `coarse_grid` before
tightening the tolerance. Velocity: `cg` with `jacobi`, or a coupled variant
when the viscosity varies.
