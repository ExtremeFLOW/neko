!> Neko binary mesh format
module nmsh
  use num_types, only : dp, i8
  implicit none
  private

  !> Neko binary mesh element data
  type, private :: nmsh_t
     integer :: el_idx          !< Element id (global)
  end type nmsh_t

  !> Neko binary mesh vertex data
  type, private :: nmsh_vertex_t
     integer :: v_idx                     !< Vertex id (global)
     real(kind=dp), dimension(3) :: v_xyz ! Vertex coordinates
  end type nmsh_vertex_t

  !> Neko quad element data
  type, public, extends(nmsh_t) :: nmsh_quad_t
     type(nmsh_vertex_t), dimension(4) :: v !< List of vertices
  end type nmsh_quad_t

  !> Neko hex element data
  type, public, extends(nmsh_t) :: nmsh_hex_t
     type(nmsh_vertex_t), dimension(8) :: v !< List of vertices
  end type nmsh_hex_t

  !> Neko zone data
  type, public :: nmsh_zone_t
     integer :: e                        !< Element id (global)
     integer :: f                        !< Facet number
     integer :: p_e                      !< Periodic connection (element)
     integer :: p_f                      !< Periodic connection (facet)
     integer, dimension(4) :: glb_pt_ids !< Global point ids
     integer :: type                     !< Zone type
  end type nmsh_zone_t

  !> Neko curve data
  type, public :: nmsh_curve_el_t
     integer :: e                                !< Element id (global)
     real(kind=dp), dimension(5,12) :: curve_data !< Save 6 values for each edge
     integer, dimension(12) :: type               !< type of curve for each edge
  end type nmsh_curve_el_t

  !> Raw nmsh mesh data need by mesh manager
  type, public :: nmsh_mesh_t
     !> Geometrical/topological dimension; given by mesh manager
     integer :: gdim
     !> Local element number; given by mesh manager
     integer :: nelt
     !> Global element number; given by mesh manager
     integer(i8) :: gnelt
     !> Global element offset; given by mesh manager
     integer(i8) :: offset_el
     !> Local number of zones
     integer :: nzone
     !> Local number of curves
     integer :: ncurve
     !> Element data
     type(nmsh_quad_t), allocatable, dimension(:) :: quad
     type(nmsh_hex_t), allocatable, dimension(:) :: hex
     !> Boundary condition
     type(nmsh_zone_t), allocatable, dimension(:) :: zone
     !> Curvature data
     type(nmsh_curve_el_t), allocatable, dimension(:) :: curve
   contains
     !> Free type structure
     procedure, pass(this) :: free => nmsh_mesh_free
  end type nmsh_mesh_t

contains

  !> Free raw nmsh data type
  subroutine nmsh_mesh_free(this)
    class(nmsh_mesh_t), intent(inout) :: this

    this%gdim = 0
    this%nelt = 0
    this%gnelt = 0
    this%offset_el = 0
    this%nzone = 0
    this%ncurve = 0

    if (allocated(this%quad)) deallocate(this%quad)
    if (allocated(this%hex)) deallocate(this%hex)
    if (allocated(this%zone)) deallocate(this%zone)
    if (allocated(this%curve)) deallocate(this%curve)

  end subroutine nmsh_mesh_free

end module nmsh
