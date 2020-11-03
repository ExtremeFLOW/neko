module krylov
  use num_types
  implicit none

  !> Defines a scratch vectorfor a Krylov solver
  type, private :: ksp_vector_t
     real(kind=dp), allocatable :: x(:)
  end type ksp_vector_t
  
  !> Defines a canonical Krylov solver
  type :: ksp_t
     procedure(ksp_method), nopass, pointer :: solve
     type(ksp_vector_t), allocatable :: v(:)
  end type ksp_t

  !> Abstract interface for a Krylov method, solving \f$ Ax = f \f$
  !!
  !! @param x vector of length @a n
  !! @param f vector of length @a n
  !! @param n integer, size of vectors
  !! @param iter iterations necessary to solve system
  abstract interface
     subroutine ksp_method(x, f, n, iter)
       import dp
       real(kind=dp), dimension(n), intent(inout) :: x
       real(kind=dp), dimension(n), intent(inout) :: f
       integer, intent(in) :: n
     end subroutine ksp_method
  end interface

  !> Abstract interface for computing \f$ Ax \f$ inside a Krylov method
  abstract interface
     subroutine ax(w, z, n)
       import dp
       real(kind=dp), dimension(n), intent(inout) :: w
       real(kind=dp), dimension(n), intent(inout) :: z
       integer, intent(in) :: n
     end subroutine ax
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
    
  end subroutine krylov_free
  
end module krylov
