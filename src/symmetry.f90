!> Mixed Dirichlet-Neumann axis aligned symmetry plane
module symmetry
  use num_types
  use dirichlet
  use coefs
  use math
  use utils
  use stack
  implicit none
  private

  !> Mixed Dirichlet-Neumann symmetry plane condition
  type, public, extends(dirichlet_t) :: symmetry_t
     integer, allocatable :: xaxis_msk(:)
     integer, allocatable :: zaxis_msk(:)
     integer, allocatable :: yaxis_msk(:)
   contains
     procedure, pass(this) :: init_msk => symmetry_init_msk
     procedure, pass(this) :: apply_scalar => symmetry_apply_scalar
     procedure, pass(this) :: apply_vector => symmetry_apply_vector
     final :: symmetry_free
  end type symmetry_t

contains

  !> Initialize symmetry mask for each axis
  subroutine symmetry_init_msk(this, c)
    class(symmetry_t), intent(inout) :: this
    type(coef_t), intent(in) :: c
    type(stack_i4_t) :: xmsk, ymsk, zmsk
    integer :: i, m, k, idx(4), facet
    
    call symmetry_free(this)

    call xmsk%init()
    call ymsk%init()
    call zmsk%init()
    
    associate(nx => c%nx, ny => c%ny, nz => c%nz)
      m = this%msk(0)
      do i = 1, m
         k = this%msk(i)
         facet = this%facet(i)
         idx = nonlinear_index(k, c%Xh%lx, c%Xh%lx, c%Xh%lx)
         select case(facet)
         case(1,2)          
            if (abs(nx(idx(2), idx(3), facet, idx(4)) - 1d0) .lt. 1d-12) then
               call xmsk%push(k)
            end if
            if (abs(ny(idx(2), idx(3), facet, idx(4)) - 1d0) .lt. 1d-12) then
               call ymsk%push(k)
            end if
            if (abs(nz(idx(2), idx(3), facet, idx(4)) - 1d0) .lt. 1d-12) then
               call zmsk%push(k)               
            end if
         case(3,4)
            if (abs(nx(idx(1), idx(3), facet, idx(4)) - 1d0) .lt. 1d-12) then
               call xmsk%push(k)
            end if
            if (abs(ny(idx(1), idx(3), facet, idx(4)) - 1d0) .lt. 1d-12) then
               call ymsk%push(k)
            end if
            if (abs(nz(idx(1), idx(3), facet, idx(4)) - 1d0) .lt. 1d-12) then
               call zmsk%push(k)               
            end if
         case(5,6)
            if (abs(nx(idx(1), idx(2), facet, idx(4)) - 1d0) .lt. 1d-12) then
               call xmsk%push(k)
            end if
            if (abs(ny(idx(1), idx(2), facet, idx(4)) - 1d0) .lt. 1d-12) then
               call ymsk%push(k)
            end if
            if (abs(nz(idx(1), idx(2), facet, idx(4)) - 1d0) .lt. 1d-12) then
               call zmsk%push(k)               
            end if
         end select
      end do
    end associate

    if (xmsk%size() .gt. 0) then
       allocate(this%xaxis_msk(0:xmsk%size()))
       this%xaxis_msk(0) = xmsk%size()
       this%xaxis_msk(1:xmsk%size()) = xmsk%array()
    else
       allocate(this%xaxis_msk(0:1))
       this%xaxis_msk(0) = 0
    end if

    if (ymsk%size() .gt. 0) then
       allocate(this%yaxis_msk(0:ymsk%size()))
       this%yaxis_msk(0) = ymsk%size()
       this%yaxis_msk(1:ymsk%size()) = ymsk%array()
    else
       allocate(this%yaxis_msk(0:1))
       this%yaxis_msk(0) = 0
    end if

    if (zmsk%size() .gt. 0) then
       allocate(this%zaxis_msk(0:zmsk%size()))
       this%zaxis_msk(0) = zmsk%size()
       this%zaxis_msk(1:zmsk%size()) = zmsk%array()
    else
       allocate(this%zaxis_msk(0:1))
       this%zaxis_msk(0) = 0
    end if

    call xmsk%free()
    call ymsk%free()
    call zmsk%free()
    
  end subroutine symmetry_init_msk
  
  subroutine symmetry_free(this)
    type(symmetry_t), intent(inout) :: this

    if (allocated(this%xaxis_msk)) then
       deallocate(this%xaxis_msk)
    end if
    
    if (allocated(this%yaxis_msk)) then
       deallocate(this%yaxis_msk)
    end if

    if (allocated(this%zaxis_msk)) then
       deallocate(this%zaxis_msk)
    end if
    
  end subroutine symmetry_free
  
  !> No-op scalar apply
  subroutine symmetry_apply_scalar(this, x, n)
    class(symmetry_t), intent(inout) :: this
    integer, intent(in) :: n
    real(kind=dp), intent(inout), dimension(n) :: x
  end subroutine symmetry_apply_scalar

  !> Apply symmetry conditions (axis aligned)
  subroutine symmetry_apply_vector(this, x, y, z, n)
    class(symmetry_t), intent(inout) :: this
    integer, intent(in) :: n
    real(kind=dp), intent(inout),  dimension(n) :: x
    real(kind=dp), intent(inout),  dimension(n) :: y
    real(kind=dp), intent(inout),  dimension(n) :: z
    integer :: i, m, k

    m = this%xaxis_msk(0)
    do i = 1, m
       k = this%xaxis_msk(i)
       x(k) = 0d0
    end do

    m = this%yaxis_msk(0)
    do i = 1, m
       k = this%yaxis_msk(i)
       y(k) = 0d0
    end do

    m = this%zaxis_msk(0)
    do i = 1, m
       k = this%zaxis_msk(i)
       z(k) = 0d0
    end do
    
  end subroutine symmetry_apply_vector
      
end module symmetry
