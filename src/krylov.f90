!> Krylov solver
module krylov
  use gather_scatter
  use num_types
  use precon
  use mesh
  use field
  use space
  use utils
  implicit none

  !> Defines a scratch vectorfor a Krylov solver
  type, private :: ksp_vector_t
     real(kind=dp), allocatable :: x(:)
  end type ksp_vector_t

  !> Defines a canonical Krylov solver
  type :: ksp_t
     procedure(ksp_method), nopass, pointer :: solve => ksp_nop
     procedure(ksp_ax), nopass, pointer :: Ax
     type(pc_t) :: M
     type(ksp_vector_t), allocatable :: v(:)
  end type ksp_t
  
  !> Abstract interface for a Krylov method, solving \f$ Ax = f \f$
  !!
  !! @param x field to solve for
  !! @param f right hand side 
  !! @param n integer, size of vectors
  !! @param iter iterations necessary to solve system
  abstract interface
     subroutine ksp_method(x, f, n, iter)
       import field_t
       import dp
       implicit none
       type(field_t) :: x
       real(kind=dp), dimension(n), intent(inout) :: f
       integer, intent(inout) :: n
       integer, intent(in) :: iter
     end subroutine ksp_method
  end interface

  !> Abstract interface for computing\f$ Ax \f$ inside a Krylov method
  !!
  !! @param w vector of length @a n
  !! @param z vector of length @a n
  !! @param g geometric factors
  !! @param gs_h gather-scatter handle
  !! @param msh mesh
  !! @param Xh function space \f$ X_h \f$
  !! @param n integer, size of vectors
  abstract interface
     subroutine ksp_ax(w, z, g, gs_h, msh, Xh, n)
       import space_t
       import mesh_t
       import gs_t
       import dp
       implicit none
       type(gs_t), intent(inout) :: gs_h
       type(mesh_t), intent(inout) :: msh
       type(space_t), intent(inout) :: Xh
       integer, intent(inout) :: n
       real(kind=dp), dimension(n), intent(inout) :: w
       real(kind=dp), dimension(n), intent(inout) :: z
       real(kind=dp), intent(inout) :: g(6, Xh%lx, Xh%ly, Xh%lz)
     end subroutine ksp_ax
  end interface

contains

  !> Create a krylov solver with @a nvec of size @a n
  subroutine krylov_init(this, nvec, n)    
    type(ksp_t), intent(inout) :: this
    integer, intent(in) :: nvec !< Number of scratch vectors
    integer, intent(in) :: n    !< Size of each scratch vectors
    integer :: i
    
    call krylov_free(this)

    allocate(this%v(nvec))
    
    do i = 1, nvec
       allocate(this%v(i)%x(n))
       this%v(i)%x = 0d0
    end do
    
  end subroutine krylov_init
  
  !> Deallocate a Krylov solver
  subroutine krylov_free(this)
    type(ksp_t), intent(inout) :: this
    integer :: i

    if (allocated(this%v)) then
       do i = 1, size(this%v)
          if (allocated(this%v(i)%x)) then
             deallocate(this%v(i)%x)
          end if
       end do
       deallocate(this%v)
    end if

    this%solve => ksp_nop
    
  end subroutine krylov_free


  !> Dummy no-op krylov method
  subroutine ksp_nop(x, f, n, iter)
    type(field_t) :: x
    real(kind=dp), dimension(n), intent(inout) :: f
    integer, intent(inout) :: n
    integer, intent(in) :: iter

    call neko_error('No Krylov method defined')
    
  end subroutine ksp_nop
  
end module krylov
