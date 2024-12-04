module fluid_scheme_cns
  use dirichlet, only : dirichlet_t
  use field, only : field_t
  use field_registry, only : neko_field_registry
  use fluid_base, only : fluid_base_t
  use json_module, only : json_file
  use logger, only : LOG_SIZE
  use num_types, only : rp
  use mesh, only : mesh_t
  use user_intf, only : user_t


  !> Base type of compressible fluid formulations.
  type, abstract, extends(fluid_base_t) :: fluid_scheme_cns_t
     !> The momentum field
     type(field_t), pointer :: m_x => null()    !< x-component of Momentum
     type(field_t), pointer :: m_y => null()    !< y-component of Momentum
     type(field_t), pointer :: m_z => null()    !< z-component of Momentum
     type(field_t), pointer :: E => null()    !< Total energy

   contains
     !> Constructors
     procedure, pass(this) :: fluid_scheme_cns_init_all
     procedure, pass(this) :: fluid_scheme_cns_init_common
     generic :: scheme_init => fluid_scheme_cns_init_all, fluid_scheme_cns_init_common

     !> Destructor for the base type
     procedure, pass(this) :: scheme_free => fluid_scheme_cns_free

     !> Compute the CFL number
     procedure, pass(this) :: compute_cfl => fluid_scheme_cns_compute_cfl

  end type fluid_scheme_cns_t

contains
  !> Initialize common data for the current scheme
  subroutine fluid_scheme_cns_init_common(this, msh, lx, params, scheme, user, &
                                      kspv_init)
    implicit none
    class(fluid_scheme_cns_t), target, intent(inout) :: this
    type(mesh_t), target, intent(inout) :: msh
    integer, intent(inout) :: lx
    character(len=*), intent(in) :: scheme
    type(json_file), target, intent(inout) :: params
    type(user_t), target, intent(in) :: user
    logical, intent(in) :: kspv_init
    type(dirichlet_t) :: bdry_mask
    character(len=LOG_SIZE) :: log_buf
    real(kind=rp), allocatable :: real_vec(:)
    real(kind=rp) :: real_val, kappa, B, z0
    logical :: logical_val
    integer :: integer_val, ierr
    type(json_file) :: wm_json
    character(len=:), allocatable :: string_val1, string_val2

    ! Fill mu and rho field with the physical value
    call this%mu_field%init(this%dm_Xh, "mu")
    call this%rho_field%init(this%dm_Xh, "rho")
    call field_cfill(this%mu_field, 0.0_rp, this%mu_field%size())
    call field_cfill(this%rho_field, 1.0_rp, this%rho_field%size())

    ! Assign momentum fields
    call neko_field_registry%add_field(this%dm_Xh, "m_x")
    call neko_field_registry%add_field(this%dm_Xh, "m_y")
    call neko_field_registry%add_field(this%dm_Xh, "m_z")
    this%m_x => neko_field_registry%get_field("m_x")
    this%m_y => neko_field_registry%get_field("m_y")
    this%m_z => neko_field_registry%get_field("m_z")
    call this%m_x%init(this%dm_Xh, "m_x")
    call this%m_y%init(this%dm_Xh, "m_y")
    call this%m_z%init(this%dm_Xh, "m_z")

    ! Assign energy field
    call neko_field_registry%add_field(this%dm_Xh, "E")
    this%E => neko_field_registry%get_field("E")
    call this%E%init(this%dm_Xh, "E")

    ! Assign velocity fields
    call neko_field_registry%add_field(this%dm_Xh, "u")
    call neko_field_registry%add_field(this%dm_Xh, "v")
    call neko_field_registry%add_field(this%dm_Xh, "w")
    this%u => neko_field_registry%get_field("u")
    this%v => neko_field_registry%get_field("v")
    this%w => neko_field_registry%get_field("w")
    call this%u%init(this%dm_Xh, "u")
    call this%v%init(this%dm_Xh, "v")
    call this%w%init(this%dm_Xh, "w")

    !! Initialize time-lag fields
    call this%ulag%init(this%u, 1)
    call this%vlag%init(this%v, 1)
    call this%wlag%init(this%w, 1)
  end subroutine fluid_scheme_cns_init_common

  subroutine fluid_scheme_cns_init_all(this, msh, lx, params, kspv_init, &
                                        kspp_init, scheme, user)
    implicit none
    class(fluid_scheme_cns_t), target, intent(inout) :: this
    type(mesh_t), target, intent(inout) :: msh
    integer, intent(inout) :: lx
    type(json_file), target, intent(inout) :: params
    type(user_t), target, intent(in) :: user
    logical :: kspv_init
    logical :: kspp_init
    character(len=*), intent(in) :: scheme
    real(kind=rp) :: abs_tol
    integer :: integer_val, ierr
    logical :: logical_val
    character(len=:), allocatable :: solver_type, precon_type
    character(len=LOG_SIZE) :: log_buf
    real(kind=rp) :: GJP_param_a, GJP_param_b

    call fluid_scheme_cns_init_common(this, msh, lx, params, scheme, user, &
                                  kspv_init)
  end subroutine fluid_scheme_cns_init_all

  subroutine fluid_scheme_cns_free(this)
    class(fluid_scheme_cns_t), intent(inout) :: this
    call this%mu_field%free()
    call this%rho_field%free()
    call this%m_x%free()
    call this%m_y%free()
    call this%m_z%free()
    call this%E%free()

    nullify(this%m_x)
    nullify(this%m_y)
    nullify(this%m_z)

    nullify(this%E)

    nullify(this%u)
    nullify(this%v)
    nullify(this%w)
    nullify(this%p)

    call this%ulag%free()
    call this%vlag%free()
    call this%wlag%free()
  end subroutine fluid_scheme_cns_free

  !> Compute CFL
  function fluid_scheme_cns_compute_cfl(this, dt) result(c)
    class(fluid_scheme_cns_t), intent(in) :: this
    real(kind=rp), intent(in) :: dt
    real(kind=rp) :: c

    c = 1.0_rp

  end function fluid_scheme_cns_compute_cfl
end module fluid_scheme_cns