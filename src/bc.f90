!> Defines a boundary condition
module bc
  use num_types
  use dofmap
  use space
  use mesh
  use stack
  use tuple
  implicit none
  private
  
  !> Base type for a boundary condition
  type, public, abstract :: bc_t
     integer(kind=8), allocatable :: msk(:)
     type(dofmap_t), pointer :: dof
     type(mesh_t), pointer :: msh
     type(space_t), pointer :: Xh
     type(stack_i4t2_t) :: marked_facet
   contains     
     procedure, pass(this) :: init => bc_init
     procedure, pass(this) :: free => bc_free
     procedure, pass(this) :: mark_facet => bc_mark_facet
     procedure, pass(this) :: finalize => bc_finalize
     procedure(bc_apply), pass(this), deferred :: apply
  end type bc_t
  
  abstract interface
     subroutine bc_apply(this, x, n)
       import :: bc_t
       import :: dp
       class(bc_t), intent(in) :: this
       integer, intent(in) :: n
       real(kind=dp), dimension(n) :: x
     end subroutine bc_apply
  end interface

contains

  !> Initialize a boundary condition type
  subroutine bc_init(this, dof)
    class(bc_t), intent(inout) :: this
    type(dofmap_t), target, intent(in) :: dof

    call bc_free(this)

    this%dof => dof
    this%Xh => dof%Xh
    this%msh => dof%msh

    call this%marked_facet%init()

  end subroutine bc_init

  !> Deallocate a boundary condition
  subroutine bc_free(this)
    class(bc_t), intent(inout) :: this

    call this%marked_facet%free()
    
    nullify(this%Xh)
    nullify(this%msh)    
    nullify(this%dof)

    if (allocated(this%msk)) then
       deallocate(this%msk)
    end if
    
  end subroutine bc_free

  !> Mark @a facet on element @a el as part of the boundary condition
  subroutine bc_mark_facet(this, facet, el)
    class(bc_t), intent(inout) :: this
    integer, intent(in) :: facet
    integer, intent(in) :: el
    type(tuple_i4_t) :: t

    t = (/facet, el/)
    call this%marked_facet%push(t)
    
  end subroutine bc_mark_facet

  !> Finalize a boundary condition
  !! @details This will linearize the marked facet's indicies in msk
  subroutine bc_finalize(this)
    class(bc_t), intent(inout) :: this
    type(tuple_i4_t), pointer :: bfp(:)
    type(tuple_i4_t) :: bc_facet
    integer :: facet_size, facet, el
    integer :: i, j, k, l, msk_c
    integer :: lx, ly, lz

    lx = this%Xh%lx
    ly = this%Xh%ly
    lz = this%Xh%lz

    !>@todo add 2D case
    
    ! Note we assume that lx = ly = lz
    facet_size = lx**2
    allocate(this%msk(0:facet_size * this%marked_facet%size()))

    msk_c = 0
    bfp => this%marked_facet%array()
    do i = 1, this%marked_facet%size()
       bc_facet = bfp(i)
       facet = bc_facet%x(1)
       el = bc_facet%x(2)
       select case (facet)
       case (1)
          do l = 1, lz
             do k = 1, ly
                msk_c = msk_c + 1
                this%msk(msk_c) = linear_index(1,k,l,el,lx,ly,lz)
             end do
          end do
       case (2)
          do l = 1, lz
             do k = 1, ly
                msk_c = msk_c + 1
                this%msk(msk_c) = linear_index(lx,k,l,el,lx,ly,lz)
             end do
          end do
       case(3)
          do l = 1, lz
             do j = 1, lx
                msk_c = msk_c + 1
                this%msk(msk_c) = linear_index(j,1,l,el,lx,ly,lz)
             end do
          end do
       case(4)
          do l = 1, lz
             do j = 1, lx
                msk_c = msk_c + 1
                this%msk(msk_c) = linear_index(j,ly,l,el,lx,ly,lz)
             end do
          end do
       case(5)
          do k = 1, ly
             do j = 1, lx
                msk_c = msk_c + 1
                this%msk(msk_c) = linear_index(j,k,1,el,lx,ly,lz)
             end do
          end do
       case(6)
          do k = 1, ly
             do j = 1, lx
                msk_c = msk_c + 1
                this%msk(msk_c) = linear_index(j,k,lz,el,lx,ly,lz)
             end do
          end do
       end select
    end do

    this%msk(0) = msk_c
    
  end subroutine bc_finalize
  
end module bc
