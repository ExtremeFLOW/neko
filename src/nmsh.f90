!> Neko binary mesh format
module nmsh
  use num_types
  implicit none


  !> Neko binary mesh element data
  type, private :: nmsh_t
     integer :: el_idx          !< Element id (global)
  end type nmsh_t

  !> Neko binary mesh vertex data
  type, private :: nmsh_vertex_t
     integer :: v_idx           !< Vertex id (global)
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
     integer :: e               !< Element id (global)
     integer :: f               !< Facet number
     integer :: p_e             !< Periodic connection (element)     
     integer :: p_f             !< Periodic connection (facet)
     integer, dimension(4) :: glb_pt_ids !< Global point ids
     integer :: type            !< Zone type
  end type nmsh_zone_t
   !> Neko curve data
  type, public :: nmsh_curve_el_t
     integer :: e               !< Element id (global)
     real(kind=dp), dimension(6,12) :: curve_data !< Save 6 values for each edge
     integer, dimension(12) :: type !< type of curve for each edge
  end type nmsh_curve_el_t
  

end module nmsh

