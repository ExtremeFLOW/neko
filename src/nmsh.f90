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
  

end module nmsh
