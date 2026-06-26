! Minimal no-op user module for the mesh-sanity smoke test.
module user
  use neko
  implicit none
contains
  subroutine user_setup(u)
    type(user_t), intent(inout) :: u
  end subroutine user_setup
end module user
