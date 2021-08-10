!> Source terms
module source
  use dofmap
  use num_types
  implicit none

  !> Defines a source term \f$ f \f$
  type :: source_t
     real(kind=rp), allocatable :: u(:,:,:,:) !< x-component of source term
     real(kind=rp), allocatable :: v(:,:,:,:) !< y-component of source term
     real(kind=rp), allocatable :: w(:,:,:,:) !< w-component of source term
     type(dofmap_t), pointer :: dm            !< dofmap for the given space
     procedure(source_term), pass(f), pointer  :: eval => null()
     procedure(source_term_pw), nopass, pointer  :: eval_pw => null()
  end type source_t

  !> Abstract interface defining how to compute a source term
  abstract interface
     subroutine source_term(f)
       import source_t            
       class(source_t) :: f
     end subroutine source_term
  end interface

  !> Abstract interface defining how to compute a source term pointwise
  abstract interface
     subroutine source_term_pw(u, v, w, j, k, l, e)
       import rp
       real(kind=rp), intent(inout) :: u
       real(kind=rp), intent(inout) :: v
       real(kind=rp), intent(inout) :: w
       integer, intent(inout) :: j
       integer, intent(inout) :: k
       integer, intent(inout) :: l
       integer, intent(inout) :: e
     end subroutine source_term_pw
  end interface
  
contains

  !> Initialize a source term @a f
  subroutine source_init(f, dm)
    type(source_t), intent(inout) :: f
    type(dofmap_t), intent(inout), target :: dm

    call source_free(f)

    f%dm => dm

    allocate(f%u(dm%Xh%lx, dm%Xh%ly, dm%Xh%lz, dm%msh%nelv))
    allocate(f%v(dm%Xh%lx, dm%Xh%ly, dm%Xh%lz, dm%msh%nelv))
    allocate(f%w(dm%Xh%lx, dm%Xh%ly, dm%Xh%lz, dm%msh%nelv))

    f%u = 0d0
    f%v = 0d0
    f%w = 0d0
    
  end subroutine source_init

  !> Deallicate a source term @a f
  subroutine source_free(f)
    type(source_t), intent(inout) :: f

    if (allocated(f%u)) then
       deallocate(f%u)
    end if

    if (allocated(f%v)) then
       deallocate(f%v)
    end if

    if (allocated(f%w)) then
       deallocate(f%w)
    end if

    nullify(f%dm)
    
  end subroutine source_free

  !> Set the eval method for the source term @a f
  subroutine source_set_type(f, f_eval)
    type(source_t), intent(inout) :: f
    procedure(source_term) :: f_eval
    f%eval => f_eval
  end subroutine source_set_type

  !> Set the pointwise eval method for the source term @a f
  subroutine source_set_pw_type(f, f_eval_pw)
    type(source_t), intent(inout) :: f
    procedure(source_term_pw) :: f_eval_pw
    f%eval => source_eval_pw
    f%eval_pw => f_eval_pw    
  end subroutine source_set_pw_type

  !> Eval routine for zero forcing
  !! @note Maybe this should be cache, avoding zeroing at each time-step
  subroutine source_eval_noforce(f)
    class(source_t) :: f
    f%u = 0d0
    f%v = 0d0
    f%w = 0d0    
  end subroutine source_eval_noforce

  !> Driver for all pointwise source term evaluatons
  subroutine source_eval_pw(f)
    class(source_t) :: f
    integer :: j, k, l, e
    integer :: jj,kk,ll,ee

    do e = 1, f%dm%msh%nelv
       ee = e
       do l = 1, f%dm%Xh%lz
          ll = l
          do k = 1, f%dm%Xh%ly
             kk = k
             do j = 1, f%dm%Xh%lx
                jj =j
                call f%eval_pw(f%u(j,k,l,e), f%v(j,k,l,e), f%w(j,k,l,e), &
                     jj, kk, ll, ee)
             end do
          end do
       end do
    end do
    
  end subroutine source_eval_pw
  
end module source
