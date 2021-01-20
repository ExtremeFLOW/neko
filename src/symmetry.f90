!> Mixed Dirichlet-Neumann axis aligned symmetry plane
module symmetry
  use num_types
  use dirichlet
  use coefs
  use math
  use utils
  use stack
  use htable
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
    type(htable_i4_t) :: algnf 
    integer :: i, m, j, k, idx(4), facet, ntype, msk_size
    real(kind=dp) :: sx,sy,sz
    real(kind=dp), parameter :: TOL = 1d-3
    
    call symmetry_free(this)

    call xmsk%init()
    call ymsk%init()
    call zmsk%init()
    call algnf%init(128)
    
    associate(nx => c%nx, ny => c%ny, nz => c%nz)
      m = this%msk(0)
      do i = 1, m
         k = this%msk(i)
         facet = this%facet(i)

         if (algnf%get(facet, ntype) .eq. 0) then         
            idx = nonlinear_index(k, c%Xh%lx, c%Xh%lx, c%Xh%lx)
            sx = 0d0
            sy = 0d0
            sz = 0d0
            select case (facet)               
            case(1,2)
               do k = 2, c%Xh%lx - 1
                  do j = 2, c%Xh%lx -1
                     sx = sx + abs(abs(nx(idx(2), idx(3), facet, idx(4))) - 1d0)
                     sy = sy + abs(abs(ny(idx(2), idx(3), facet, idx(4))) - 1d0)
                     sz = sz + abs(abs(nz(idx(2), idx(3), facet, idx(4))) - 1d0)
                  end do
               end do
            case(3,4)
               do k = 1, c%Xh%lx
                  do j = 1, c%Xh%lx
                     sx = sx + abs(abs(nx(idx(1), idx(3), facet, idx(4))) - 1d0)
                     sy = sy + abs(abs(ny(idx(1), idx(3), facet, idx(4))) - 1d0)
                     sz = sz + abs(abs(nz(idx(1), idx(3), facet, idx(4))) - 1d0)
                  end do
               end do
            case(5,6)
               do k = 1, c%Xh%lx
                  do j = 1, c%Xh%lx
                     sx = sx + abs(abs(nx(idx(1), idx(2), facet, idx(4))) - 1d0)
                     sy = sy + abs(abs(ny(idx(1), idx(2), facet, idx(4))) - 1d0)
                     sz = sz + abs(abs(nz(idx(1), idx(2), facet, idx(4))) - 1d0)
                  end do
               end do               
            end select
            sx = sx / (c%Xh%lx - 2)**2
            sy = sy / (c%Xh%lx - 2)**2
            sz = sz / (c%Xh%lx - 2)**2

            ntype = 0
            if (sx .lt. TOL) then
               ntype = iand(ntype, 1)
               call xmsk%push(k)
            end if

            if (sy .lt. TOL) then
               ntype = iand(ntype, 2)
               call ymsk%push(k)
            end if

            if (sz .lt. TOL) then
               ntype = iand(ntype, 4)
               call zmsk%push(k)
            end if

            call algnf%set(facet, ntype)
         else
            if (iand(ntype, 1) .eq. 1) then
               call xmsk%push(k)
            end if
            if (iand(ntype, 2) .eq. 2) then
               call ymsk%push(k)
            end if
            if (iand(ntype, 4) .eq. 4) then
               call zmsk%push(k)
            end if
         end if
      end do
    end associate

    !> @note This is to prevent runtime issues with Cray ftn, gfortran and
    !! msk:size() in the allocate call
    msk_size = xmsk%size()
    if (msk_size .gt. 0) then
       allocate(this%xaxis_msk(0:msk_size))
       this%xaxis_msk(0) = msk_size
       this%xaxis_msk(1:msk_size) = xmsk%array()
    else
       allocate(this%xaxis_msk(0:1))
       this%xaxis_msk(0) = 0
    end if

    msk_size = ymsk%size()
    if (msk_size .gt. 0) then
       allocate(this%yaxis_msk(0:msk_size))
       this%yaxis_msk(0) = msk_size
       this%yaxis_msk(1:msk_size) = ymsk%array()
    else
       allocate(this%yaxis_msk(0:1))
       this%yaxis_msk(0) = 0
    end if

    msk_size = zmsk%size()
    if (msk_size .gt. 0) then
       allocate(this%zaxis_msk(0:msk_size))
       this%zaxis_msk(0) = msk_size
       this%zaxis_msk(1:msk_size) = zmsk%array()
    else
       allocate(this%zaxis_msk(0:1))
       this%zaxis_msk(0) = 0
    end if

    call xmsk%free()
    call ymsk%free()
    call zmsk%free()
    call algnf%free()
    
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
