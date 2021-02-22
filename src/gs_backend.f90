!> Defines a gather-scatter backend
module gs_backend
  use num_types
  implicit none
    
  type, public, abstract :: gs_backend_t
   contains
     procedure(gs_backend_init), pass, deferred :: init
     procedure(gs_gather), pass(this), deferred :: gather
     procedure(gs_scatter), pass(this), deferred :: scatter
  end type gs_backend_t

  !> Abstract interface for initialising a Gather-Scatter backend
  abstract interface
     subroutine gs_backend_init(this, nlocal, nshared)
       import gs_backend_t
       class(gs_backend_t), intent(inout) :: this
       integer, intent(in) :: nlocal
       integer, intent(in) :: nshared
     end subroutine gs_backend_init
  end interface

  !> Abstract interface for the Gather kernel
  !! \f$ v(dg(i)) = op(v(dg(i)), u(gd(i)) \f$
  abstract interface
     subroutine gs_gather(this, v, m, o, dg, u, n, gd, nb, b, op)
       import gs_backend_t       
       import dp
       integer, intent(inout) :: m
       integer, intent(inout) :: n
       integer, intent(inout) :: nb
       class(gs_backend_t), intent(in) :: this
       real(kind=dp), dimension(m), intent(inout) :: v
       integer, dimension(m), intent(inout) :: dg
       real(kind=dp), dimension(n), intent(inout) :: u
       integer, dimension(m), intent(inout) :: gd
       integer, dimension(nb), intent(inout) :: b
       integer, intent(inout) :: o
       integer :: op
    
     end subroutine gs_gather
  end interface
  
  !> Abstract interface for the Scatter kernel
  !! \f$ u(gd(i) = v(dg(i)) \f$
  abstract interface
     subroutine gs_scatter(this, v, m, dg, u, n, gd, nb, b)
       import gs_backend_t       
       import dp
       integer, intent(in) :: m
       integer, intent(in) :: n
       integer, intent(in) :: nb
       class(gs_backend_t), intent(in) :: this              
       real(kind=dp), dimension(m), intent(inout) :: v
       integer, dimension(m), intent(inout) :: dg
       real(kind=dp), dimension(n), intent(inout) :: u
       integer, dimension(m), intent(inout) :: gd
       integer, dimension(nb), intent(inout) :: b
       integer :: i, j, k, blk_len
       real(kind=dp) :: tmp       
     end subroutine gs_scatter
  end interface

  
end module gs_backend
