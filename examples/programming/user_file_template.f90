! A template user file containing the user-defined functions
!
module user
  use neko
  implicit none

contains

  ! Register user-defined functions (see user_intf.f90)
  subroutine user_setup(user)
    type(user_t), intent(inout) :: user

    user%user_startup => user_startup
    user%user_init_modules => user_init_modules
    user%fluid_user_ic => fluid_user_ic
    user%fluid_compressible_user_ic => fluid_compressible_user_ic
    user%scalar_user_ic => scalar_user_ic
    user%user_mesh_setup => user_mesh_setup
    user%user_check => user_check
    user%user_finalize_modules => user_finalize_modules
    user%fluid_user_f => fluid_user_f
    user%fluid_user_f_vector => fluid_user_f_vector
    user%scalar_user_f => scalar_user_f
    user%scalar_user_f_vector => scalar_user_f_vector
    user%scalar_user_bc => scalar_user_bc
    user%user_dirichlet_update => user_dirichlet_update
    user%material_properties => material_properties

  end subroutine user_setup

  subroutine user_startup(params)
    type(json_file), intent(inout) :: params

  end subroutine user_startup

  subroutine user_init_modules(t, u, v, w, p, coef, params)
    real(kind=rp) :: t
    type(field_t), intent(inout) :: u, v, w, p
    type(coef_t), intent(inout) :: coef
    type(json_file), intent(inout) :: params

  end subroutine user_init_modules

  subroutine fluid_user_ic(u, v, w, p, params)
    type(field_t), intent(inout) :: u, v, w, p
    type(json_file), intent(inout) :: params

  end subroutine fluid_user_ic

  subroutine fluid_compressible_user_ic(rho, u, v, w, p, params)
    type(field_t), intent(inout) :: rho, u, v, w, p
    type(json_file), intent(inout) :: params

  end subroutine fluid_compressible_user_ic

  subroutine scalar_user_ic(s, params)
    type(field_t), intent(inout) :: s
    type(json_file), intent(inout) :: params

  end subroutine scalar_user_ic

  subroutine user_mesh_setup(msh)
    type(mesh_t), intent(inout) :: msh

  end subroutine user_mesh_setup

  subroutine user_check(t, tstep, u, v, w, p, coef, params)
    real(kind=rp), intent(in) :: t
    integer, intent(in) :: tstep
    type(field_t), intent(inout) :: u, v, w, p
    type(coef_t), intent(inout) :: coef
    type(json_file), intent(inout) :: params

  end subroutine user_check

  subroutine user_finalize_modules(t, params)
    real(kind=rp) :: t
    type(json_file), intent(inout) :: params

  end subroutine user_finalize_modules

  subroutine fluid_user_f(u, v, w, j, k, l, e, t)
    real(kind=rp), intent(inout) :: u, v, w
    integer, intent(in) :: j, k, l, e
    real(kind=rp), intent(in) :: t

  end subroutine fluid_user_f

  subroutine fluid_user_f_vector(f, t)
    class(fluid_user_source_term_t), intent(inout) :: f
    real(kind=rp), intent(in) :: t

  end subroutine fluid_user_f_vector

  subroutine scalar_user_f(field_name, s, j, k, l, e, t)
    character(len=*), intent(in) :: field_name
    real(kind=rp), intent(inout) :: s
    integer, intent(in) :: j, k, l, e
    real(kind=rp), intent(in) :: t

  end subroutine scalar_user_f

  subroutine scalar_user_f_vector(field_name, f, t)
    character(len=*), intent(in) :: field_name
    class(scalar_user_source_term_t), intent(inout) :: f
    real(kind=rp), intent(in) :: t

  end subroutine scalar_user_f_vector

  subroutine scalar_user_bc(scalar_name, s, x, y, z, nx, ny, nz, ix, iy, iz, ie, t, tstep)
    character(len=*), intent(in) :: scalar_name
    real(kind=rp), intent(inout) :: s
    real(kind=rp), intent(in) :: x, y, z, nx, ny, nz
    integer, intent(in) :: ix, iy, iz, ie
    real(kind=rp), intent(in) :: t
    integer, intent(in) :: tstep

  end subroutine scalar_user_bc

  subroutine user_dirichlet_update(fields, bc, coef, time)
    type(field_list_t), intent(inout) :: fields
    type(field_dirichlet_t), intent(in) :: bc
    type(coef_t), intent(inout) :: coef
    type(time_state_t), intent(in) :: time

  end subroutine user_dirichlet_update

  subroutine material_properties(t, tstep, name, properties)
    real(kind=rp), intent(in) :: t
    integer, intent(in) :: tstep
    character(len=*), intent(in) :: name
    type(field_list_t), intent(inout) :: properties

  end subroutine material_properties

end module user
