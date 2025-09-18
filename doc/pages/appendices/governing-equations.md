# Governing equations {#governing-equations}

\tableofcontents

## Momentum (Fluid)

In general, Neko integrates the Navier-Stokes equations in time formulated as follows:

$$\rho \left( \frac{\partial u_i}{\partial t} +
   u_j \frac{\partial u_i}{\partial x_j} \right) =
  -\frac{\partial p}{\partial x_i} +
  \mu \frac{\partial u_i }{\partial x_j \partial x_j} +
  \rho \sum_j f^u_{i,j} , \quad i=1,2,3,$$

together with the continuity equation

$$ \frac{\partial u_i}{\partial x_i} = 0 \ .$$

Here, \f$ u_i \f$ are the components of the velocity vector, \f$ p \f$ is the
pressure, \f$ \rho \f$ is the density, \f$ \mu \f$ is the dynamic viscosity, and
\f$ f^u_{i,j} \f$ are the components of the source terms. Both \f$ \rho \f$ and
\f$ \mu \f$ are usually treated as constants.

When the viscosity is varying in space (e.g. for LES modelling), the viscous
stress tensor cannot be simplified and the following equations are solved in a
coupled manner (the so-called stress formulation):

$$\rho \left( \frac{\partial u_i}{\partial t} +
   u_j \frac{\partial u_i}{\partial x_j} \right) =
  -\frac{\partial p}{\partial x_i} +
  \frac{\partial}{\partial x_j} \left(\mu_{tot}
  \left(\frac{\partial u_i }{\partial x_j} +
   \frac{\partial u_j }{\partial x_i} \right) \right)  +
  \rho \sum_j f^u_{i,j} , \quad i=1,2,3.$$

Here, \f$ \mu_{tot} \f$ is the total viscosity field, potentially including the
contribution from turbulence modelling.

## Scalar

Optionally, an additional equation for scalar transport can be included. Here we
formulate it as an equation for temperature, but the physical meaning can of course
differ from case to case.

$$\rho c_p \left( \frac{\partial T}{\partial t} +
   u_j \frac{\partial T}{\partial x_j} \right) =
  \frac{\partial}{ \partial x_j}
  \left(\lambda_{tot} \frac{\partial T }{ \partial x_j} \right)+
  \rho c_p \sum_j f_j^s.$$

Here, \f$ T \f$ is the scalar temperature field, \f$ c_p \f$ is the specific
heat capacity, \f$ \lambda_{tot} \f$ is the total thermal conductivity, and \f$ f_j^s
\f$ are the components of the active source terms. As for the momentum,  for constant \f$\lambda\f$ the Laplacian of the temperature may be used.

## Non-dimensionalisation
A non-dimensional form of the Navier-Stokes equations may be found by defining the Reynolds number
$$ Re = \frac{\rho U L}{\mu} \ ,$$
with reference velocity \f$U\f$, reference length \f$L\f$, density \f$\rho\f$ and dynamic viscosity \f$\mu\f$; the latter two may be
combined into the kinematic viscosity \f$\nu=\mu/\rho\f$. In order to get a properly scaled problem, one can define 
the mesh such that the reference length is unity, and set the velocity boundary conditions (or forcing) such that the reference velocity 
is unity. Then specifying the Reynolds number (e.g. via the case file) will set the (now non-dimensional) 
density \f$\rho=1\f$ and the non-dimensional viscosity to \f$\nu = 1/Re\f$. Similarly for a scalar, a Prandtl number \f$Pr\f$ may be specified, which together with the Reynolds number defines the relevant P\'{e}clet number \f$Pe = Pr \cdot Re\f$.
