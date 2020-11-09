!> Defines various Conjugate Gradient methods
module cg
  use krylov
  use math
  use num_types
  implicit none

  !> Standard preconditioned conjugate gradient method
  type, public, extends(ksp_t) :: cg_t
     real(kind=dp), allocatable :: w(:)
     real(kind=dp), allocatable :: c(:)
     real(kind=dp), allocatable :: r(:)
     real(kind=dp), allocatable :: p(:)
     real(kind=dp), allocatable :: z(:)
   contains
     procedure, pass(this) :: init => cg_init
     procedure, pass(this) :: free => cg_free
     procedure, pass(this) :: solve => cg_solve
  end type cg_t

contains

  !> Initialise a standard PCG solver
  subroutine cg_init(this)
    class(cg_t), intent(inout) :: this

    call this%free()
    
    call this%ksp()
       
  end subroutine cg_init

  !> Deallocate a standard PCG solver
  subroutine cg_free(this)
    class(cg_t), intent(inout) :: this

    call this%ksp_free()

    if (allocated(this%w)) then
       deallocate(this%w)
    end if

    if (allocated(this%c)) then
       deallocate(this%c)
    end if

    if (allocated(this%r)) then
       deallocate(this%r)
    end if

    if (allocated(this%p)) then
       deallocate(this%p)
    end if
    
    if (allocated(this%z)) then
       deallocate(this%z)
    end if

  end subroutine cg_free
  
  !> Standard PCG solve
  subroutine cg_solve(this, Ax, x, f, n, niter)
    class(cg_t), intent(inout) :: this
    class(ax_t), intent(inout) :: Ax
    type(field_t), intent(inout) :: x
    real(kind=dp), dimension(n), intent(inout) :: f
    integer, intent(inout) :: n
    integer, optional, intent(in) :: niter
    integer :: iter, max_iter
    real(kind=dp) :: rnorm, rtr, rtr0, rtz2, rtz1
    real(kind=dp) :: beta, pap, alpha, alphm, eps

    if (present(niter)) then
       max_iter = niter
    else
       max_iter = KSP_MAX_ITER
    end if
    
    call rzero(x%x, n)
    call copy(this%r, f, n)
    !> @todo add masking call

    rnorm = sqrt(glsc3(this%r, this%c, this%r, n))
    do iter = 1, niter
       call this%M%solve(this%z, this%r, n)

       rtz2 = rtz1
       rtz1 = glsc3(this%r, this%c, this%z, n)

       beta = rtz1 / rtz2
       if (iter .eq. 1) beta = 0d0
       call add2s1(this%p, this%z, beta, n)

       !       call this%Ax(this%w, this%z, 
       pap = glsc3(this%w, this%c, this%p, n)

       alpha = rtz1 / pap
       alphm = -alpha
       call add2s2(x%x, this%p, alpha, n)
       call add2s2(this%r, this%w, alphm, n)

       rtr = glsc3(this%r, this%c, this%r, n)
       if (iter .eq. 1) rtr0 = rtr
       rnorm = sqrt(rtr)       
    end do
    
  end subroutine cg_solve

end module cg
  

