module fluid_scheme_compressible
  use dirichlet, only : dirichlet_t
  use field, only : field_t
  use field_math, only : field_cfill
  use field_registry, only : neko_field_registry
  use fluid_scheme_base, only : fluid_scheme_base_t
  use json_module, only : json_file
  use logger, only : LOG_SIZE
  use num_types, only : rp
  use mesh, only : mesh_t
  use scratch_registry, only : scratch_registry_t
  use space, only : space_t, GLL
  use user_intf, only : user_t
  use usr_inflow, only : usr_inflow_eval

  !> Base type of compressible fluid formulations.
  type, abstract, extends(fluid_scheme_base_t) :: fluid_scheme_compressible_t
     !> The momentum field
     type(field_t), pointer :: m_x => null()    !< x-component of Momentum
     type(field_t), pointer :: m_y => null()    !< y-component of Momentum
     type(field_t), pointer :: m_z => null()    !< z-component of Momentum
     type(field_t), pointer :: E => null()    !< Total energy

     type(scratch_registry_t) :: scratch       !< Manager for temporary fields

   contains
     !> Constructors
     procedure, pass(this) :: scheme_init => fluid_scheme_compressible_init
     !> Destructor for the base type
     procedure, pass(this) :: scheme_free => fluid_scheme_compressible_free

     !> Validate that all components are properly allocated
     procedure, pass(this) :: validate => fluid_scheme_compressible_validate
     !> Set the user inflow procedure
     procedure, pass(this) :: set_usr_inflow => fluid_scheme_compressible_set_usr_inflow
     !> Compute the CFL number
     procedure, pass(this) :: compute_cfl => fluid_scheme_compressible_compute_cfl
     !> Set rho and mu
     procedure, pass(this) :: update_material_properties => &
          fluid_scheme_compressible_update_material_properties

  end type fluid_scheme_compressible_t

contains
  !> Initialize common data for the current scheme
  subroutine fluid_scheme_compressible_init(this, msh, lx, params, scheme, user)
    implicit none
    class(fluid_scheme_compressible_t), target, intent(inout) :: this
    type(mesh_t), target, intent(inout) :: msh
    integer, intent(inout) :: lx
    character(len=*), intent(in) :: scheme
    type(json_file), target, intent(inout) :: params
    type(user_t), target, intent(in) :: user

    !
    ! SEM simulation fundamentals
    !

    this%msh => msh

    if (msh%gdim .eq. 2) then
       call this%Xh%init(GLL, lx, lx)
    else
       call this%Xh%init(GLL, lx, lx, lx)
    end if

    call this%dm_Xh%init(msh, this%Xh)

    call this%gs_Xh%init(this%dm_Xh)

    call this%c_Xh%init(this%gs_Xh)

    ! Local scratch registry
    this%scratch = scratch_registry_t(this%dm_Xh, 10, 2)

    ! Case parameters
    this%params => params

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

    ! ! Assign velocity fields
    call neko_field_registry%add_field(this%dm_Xh, "u")
    call neko_field_registry%add_field(this%dm_Xh, "v")
    call neko_field_registry%add_field(this%dm_Xh, "w")
    this%u => neko_field_registry%get_field("u")
    this%v => neko_field_registry%get_field("v")
    this%w => neko_field_registry%get_field("w")
    call this%u%init(this%dm_Xh, "u")
    call this%v%init(this%dm_Xh, "v")
    call this%w%init(this%dm_Xh, "w")
    call neko_field_registry%add_field(this%dm_Xh, 'p')
    this%p => neko_field_registry%get_field('p')
    call this%p%init(this%dm_Xh, "p")

    ! !! Initialize time-lag fields
    call this%ulag%init(this%u, 1)
    call this%vlag%init(this%v, 1)
    call this%wlag%init(this%w, 1)
  end subroutine fluid_scheme_compressible_init

  subroutine fluid_scheme_compressible_free(this)
    class(fluid_scheme_compressible_t), intent(inout) :: this
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
  end subroutine fluid_scheme_compressible_free

  !> Validate that all components are properly allocated
  subroutine fluid_scheme_compressible_validate(this)
    class(fluid_scheme_compressible_t), target, intent(inout) :: this
    !> TODO: fill here
  end subroutine fluid_scheme_compressible_validate

  !> Initialize a user defined inflow condition
  subroutine fluid_scheme_compressible_set_usr_inflow(this, usr_eval)
    class(fluid_scheme_compressible_t), intent(inout) :: this
    procedure(usr_inflow_eval) :: usr_eval
    !> TODO: fill here
  end subroutine fluid_scheme_compressible_set_usr_inflow

  !> Compute CFL
  function fluid_scheme_compressible_compute_cfl(this, dt) result(c)
    class(fluid_scheme_compressible_t), intent(in) :: this
    real(kind=rp), intent(in) :: dt
    real(kind=rp) :: c

    c = 1.0_rp
    !> TODO: fill here

  end function fluid_scheme_compressible_compute_cfl

  !> Set rho and mu
  subroutine fluid_scheme_compressible_update_material_properties(this)
    class(fluid_scheme_compressible_t), intent(inout) :: this
    !> TODO: fill here, may not be used?
  end subroutine fluid_scheme_compressible_update_material_properties
end module fluid_scheme_compressible