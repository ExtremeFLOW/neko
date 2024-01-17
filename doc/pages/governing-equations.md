# Governing equations {#governing-equations}

## Fluid

Neko's solves the Navier-Stokes equations formulated as follows.

$$\rho \left( \frac{\partial u_i}{\partial t} + 
   u_j \frac{\partial u_i}{\partial x_j} \right) =
  -\frac{\partial p}{\partial x_i} + 
  \mu \frac{\partial u_i }{\partial x_j \partial x_j} + 
  \rho \sum_j f^u_{i,j} , \quad i=1,2,3.$$

Here, \f$ u_i \f$ are the components of the velocity vector, \f$ p \f$ is the
pressure, \f$ \rho \f$ is the density, \f$ \mu \f$ is the dynamic viscosity, and
\f$ f^u_{i,j} \f$ are the components of the source terms. Both \f$ \rho \f$ and
\f$ \mu \f$ are treated as constants.

## Scalar

Optionally, and additional equation for scalar transport can be solved. Here we
define it as an equation for temperature, but the physical meaning can of course
differ from case to case.

$$\rho c_p \left( \frac{\partial s}{\partial t} + 
   u_j \frac{\partial s}{\partial x_j} \right) =
  \lambda \frac{\partial^2 s }{\partial x_j \partial x_j} + 
  \rho c_p \sum_j f_j^s.$$

Here, \f$ s \f$ is the scalar field, \f$ c_p \f$ is the specific heat capacity,
\f$ \lambda \f$ is the thermal conductivity, and \f$ f_j^s \f$ is the active
source terms. Both \f$ c_p \f$ and \f$ \lambda \f$ are treated as constants.
