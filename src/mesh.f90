!> Defines a mesh
module mesh
  use num_types
  implicit none

  type mesh_t

     integer :: lelv
     integer :: dim 

     integer :: lx1
     integer :: ly1
     integer :: lz1

  end type mesh_t


end module mesh
