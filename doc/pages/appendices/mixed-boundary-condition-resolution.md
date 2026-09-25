# Mixed boundary-condition resolution {#mixed-bc-resolution}

\tableofcontents

Mixed vector boundary conditions constrain components relative to the local
boundary normal. At a node in the interior of a boundary face, that normal is
unambiguous. At an edge or corner, however, several boundary faces can meet and
may impose different kinds of constraint. This appendix defines how those face
constraints are reduced to one homogeneous nodal constraint.

The construction acts on vectors used by a linear solve, such as corrections
or residuals, for which the boundary constraints are homogeneous. Prescribed
non-zero boundary values are outside the scope of this construction.

## Boundary classes and priority {#mixed-bc-priority}

Each incident boundary face has one of four class values:

| Class | Constraint | Admissible homogeneous vectors |
|------:|------------|--------------------------------|
| 0 | All components | Only the zero vector |
| 2 | Normal component | Vectors tangent to the boundary |
| 3 | Tangential components | Vectors parallel to the boundary normal |
| 5 | None | All vectors |

The numerical values define the priority order

$$0 \mathrel{\succ} 2 \mathrel{\succ} 3 \mathrel{\succ} 5,$$

where a smaller value has higher priority. This is a policy for resolving
different boundary classes; it is not an ordering by the dimension of their
constraint spaces.

Let \f$\mathcal{F}(x)\f$ be the set of boundary faces incident on a node
\f$x\f$, and let \f$c_f\f$ denote the class of face \f$f\f$. The nodal class is

$$c(x) = \min_{f \in \mathcal{F}(x)} c_f.$$

Consequently, a fully constrained face takes priority over every other class,
and a normal-component constraint takes priority over a tangential-component
constraint.

## Resolution of the nodal normal {#mixed-bc-normal}

A local normal is required only when \f$c(x)\f$ is class 2 or class 3. Among
the faces incident on \f$x\f$, only those whose class equals the resolved nodal
class contribute. Define

$$
\mathcal{F}_c(x) =
\left\{ f \in \mathcal{F}(x) : c_f = c(x) \right\}.
$$

If \f$\boldsymbol{n}_f(x)\f$ is the unit outward normal of face \f$f\f$ at the
node, the resolved normal is the normalised, equally weighted sum

$$
\boldsymbol{n}(x) =
\frac{\displaystyle\sum_{f \in \mathcal{F}_c(x)}
      \boldsymbol{n}_f(x)}
     {\displaystyle\left\|\sum_{f \in \mathcal{F}_c(x)}
      \boldsymbol{n}_f(x)\right\|_2}.
$$

Thus, a lower-priority face neither changes the nodal class nor contributes to
the resolved normal. Faces of the selected class contribute symmetrically. The
construction requires their normal sum to be non-zero; otherwise, no resolved
normal exists and the boundary configuration is invalid for this rule.

If nodes are identified through a rotated periodic mapping, their face normals
are expressed in a common frame before being summed. The resolved normal is
then expressed in the node's local physical frame.

Choose any two unit vectors \f$\boldsymbol{t}_1(x)\f$ and
\f$\boldsymbol{t}_2(x)\f$ that complete an orthonormal basis

$$
\left\{\boldsymbol{n},\boldsymbol{t}_1,\boldsymbol{t}_2\right\}.
$$

The particular choice of tangential directions does not affect the resulting
constraint because both tangential directions are always treated alike.

## Nodal projectors {#mixed-bc-projectors}

Let \f$I\f$ be the \f$3 \times 3\f$ identity matrix and set

$$N = \boldsymbol{n}\boldsymbol{n}^{\mathsf T}, \qquad
  T = I - \boldsymbol{n}\boldsymbol{n}^{\mathsf T}.$$

For a vector \f$\boldsymbol{u}\f$, the homogeneous nodal constraint is enforced
by \f$\boldsymbol{u}_{c}=P_c\boldsymbol{u}\f$, with

$$
P_0 = 0, \qquad
P_2 = T, \qquad
P_3 = N, \qquad
P_5 = I.
$$

Class 2 therefore removes the normal component and preserves the tangential
plane. Class 3 removes both tangential components and preserves only the normal
component. Each matrix is symmetric and idempotent,

$$P_c^{\mathsf T}=P_c, \qquad P_c^2=P_c,$$

so the nodal operation is an orthogonal projection onto the selected admissible
subspace.

## Interpretation at edges and corners {#mixed-bc-interpretation}

The normal-sum rule selects one representative normal and therefore one nodal
subspace. It does not form the exact intersection of the constraint subspaces
from all incident faces.

For example, suppose two class-2 faces meet with orthogonal normals
\f$\boldsymbol{e}_1\f$ and \f$\boldsymbol{e}_2\f$. The resolved normal is

$$\boldsymbol{n} =
  \frac{\boldsymbol{e}_1+\boldsymbol{e}_2}{\sqrt{2}},$$

and the projected vectors satisfy \f$u_1+u_2=0\f$. Enforcing both face
constraints independently would instead require \f$u_1=u_2=0\f$. The former
is the deliberate compromise made by the normal-sum rule. The same distinction
applies when several class-3 faces meet.

When different mixed classes meet, priority is applied before the normal is
formed. For example, at a class-2/class-3 junction, the node is assigned class
2 and only class-2 face normals enter the sum. The result is therefore always a
single projector associated with the highest-priority incident class.
