# Governing equations {#governing-equations}

\tableofcontents

## Fluid

When possible, Neko solves the Navier-Stokes equations formulated as follows.

$$\rho \left( \frac{\partial u_i}{\partial t} + 
   u_j \frac{\partial u_i}{\partial x_j} \right) =
  -\frac{\partial p}{\partial x_i} + 
  \mu \frac{\partial u_i }{\partial x_j \partial x_j} + 
  \rho \sum_j f^u_{i,j} , \quad i=1,2,3.$$

Here, \f$ u_i \f$ are the components of the velocity vector, \f$ p \f$ is the
pressure, \f$ \rho \f$ is the density, \f$ \mu \f$ is the dynamic viscosity, and
\f$ f^u_{i,j} \f$ are the components of the source terms. Both \f$ \rho \f$ and
\f$ \mu \f$ are treated as constants.

When the viscosity is varying in space (e.g. for LES modelling), the viscous
stress tensor cannot be simplified and the following equations are solved in a
coupled manner.

$$\rho \left( \frac{\partial u_i}{\partial t} + 
   u_j \frac{\partial u_i}{\partial x_j} \right) =
  -\frac{\partial p}{\partial x_i} + 
  \frac{\partial}{\partial x_j} \left(\mu_{tot} 
  \left(\frac{\partial u_i }{\partial x_j} +
   \frac{\partial u_j }{\partial x_i} \right) \right)  + 
  \rho \sum_j f^u_{i,j} , \quad i=1,2,3.$$

Here, \f$ \mu_{tot} \f$ is the total viscosity field, potentially including the
contribution of turbulence modelling.

## Scalar

Optionally, and additional equation for scalar transport can be solved. Here we
define it as an equation for temperature, but the physical meaning can of course
differ from case to case.

$$\rho c_p \left( \frac{\partial s}{\partial t} + 
   u_j \frac{\partial s}{\partial x_j} \right) =
  \frac{\partial}{ \partial x_j} 
  \left(\lambda_{tot} \frac{\partial s }{ \partial x_j} \right)+ 
  \rho c_p \sum_j f_j^s.$$

Here, \f$ s \f$ is the scalar field, \f$ c_p \f$ is the specific heat capacity,
\f$ \lambda \f$ is the total thermal conductivity, and \f$ f_j^s \f$ is the
active source terms. 
