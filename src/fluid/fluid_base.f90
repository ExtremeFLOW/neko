module fluid_base
  use bc, only : bc_t, bc_list_t
  use checkpoint, only : chkp_t
  use coefs, only: coef_t
  use dofmap, only : dofmap_t
  use field, only : field_t
  use field_series, only : field_series_t
  use gather_scatter, only : gs_t
  use json_module, only : json_file
  use num_types, only : rp
  use mesh, only : mesh_t, NEKO_MSH_MAX_ZLBL_LEN
  use space, only : space_t, GLL
  use wall, only : no_slip_wall_t

  !> Base type of all fluid formulations.
  type, abstract :: fluid_base_t
     type(space_t) :: Xh        !< Function space \f$ X_h \f$
     type(dofmap_t) :: dm_Xh    !< Dofmap associated with \f$ X_h \f$
     type(gs_t) :: gs_Xh        !< Gather-scatter associated with \f$ X_h \f$
     type(coef_t) :: c_Xh       !< Coefficients associated with \f$ X_h \f$

     !> The velocity field
     type(field_t), pointer :: u => null()    !< x-component of Velocity
     type(field_t), pointer :: v => null()    !< y-component of Velocity
     type(field_t), pointer :: w => null()    !< z-component of Velocity
     type(field_t), pointer :: p => null()    !< Pressure
     type(field_series_t) :: ulag, vlag, wlag !< fluid field (lag)

     !> The variable density field
     type(field_t) :: rho_field

     !> Boundary conditions
     type(field_t) :: bdry                     !< Boundary markings
     type(no_slip_wall_t) :: bc_wall           !< No-slip wall for velocity
     class(bc_t), allocatable :: bc_inflow     !< Dirichlet inflow for velocity
     type(bc_list_t) :: bclst_vel              !< List of velocity conditions
     type(bc_list_t) :: bclst_vel_neumann      !< List of neumann velocity conditions

     type(json_file), pointer :: params        !< Parameters
     type(mesh_t), pointer :: msh => null()    !< Mesh
     type(chkp_t) :: chkp                      !< Checkpoint

     !> Boundary condition labels (if any)
     character(len=NEKO_MSH_MAX_ZLBL_LEN), allocatable :: bc_labels(:)

     !> Dynamic viscosity
     real(kind=rp) :: mu

     !> The variable mu field
     type(field_t) :: mu_field
     
     
  end type fluid_base_t
end module fluid_base