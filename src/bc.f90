!> Defines a boundary condition
module bc
  use num_types
  use dofmap
  use space
  use mesh
  use zone
  use stack
  use tuple
  use utils
  implicit none
  private
  
  !> Base type for a boundary condition
  type, public, abstract :: bc_t
     integer, allocatable :: msk(:)
     type(dofmap_t), pointer :: dof
     type(mesh_t), pointer :: msh
     type(space_t), pointer :: Xh
     type(stack_i4t2_t) :: marked_facet
   contains     
     procedure, pass(this) :: init => bc_init
     procedure, pass(this) :: free => bc_free
     procedure, pass(this) :: mark_facet => bc_mark_facet
     procedure, pass(this) :: mark_facets => bc_mark_facets
     procedure, pass(this) :: mark_zone => bc_mark_zone
     procedure, pass(this) :: finalize => bc_finalize
     procedure(bc_apply), pass(this), deferred :: apply
     procedure(bc_apply_mult), pass(this), deferred :: apply_mult
  end type bc_t

  !> Pointer to boundary condtiion
  type, private :: bcp_t
     class(bc_t), pointer :: bcp
  end type bcp_t
  
  !> A list of boundary conditions
  type, public :: bc_list_t
     type(bcp_t), allocatable :: bc(:)
     integer :: n
     integer :: size
  end type bc_list_t
    
  abstract interface
     subroutine bc_apply(this, x, n)
       import :: bc_t
       import :: dp
       class(bc_t), intent(inout) :: this
       integer, intent(in) :: n
       real(kind=dp), intent(inout), dimension(n) :: x
     end subroutine bc_apply
  end interface

  abstract interface
     subroutine bc_apply_mult(this, x, y, z, n)
       import :: bc_t
       import :: dp
       class(bc_t), intent(inout) :: this
       integer, intent(in) :: n
       real(kind=dp), intent(inout), dimension(n) :: x
       real(kind=dp), intent(inout), dimension(n) :: y
       real(kind=dp), intent(inout), dimension(n) :: z
     end subroutine bc_apply_mult
  end interface

  public :: bc_list_init, bc_list_free, bc_list_add, bc_list_apply
  
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

  !> Mark all facets from a (facet, el) tuple list
  subroutine bc_mark_facets(this, facet_list)
    class(bc_t), intent(inout) :: this
    type(stack_i4t2_t), intent(inout) :: facet_list
    type(tuple_i4_t), pointer :: fp(:)
    integer :: i

    fp => facet_list%array()
    do i = 1, facet_list%size()
       call this%marked_facet%push(fp(i))
    end do
       
  end subroutine bc_mark_facets

  !> Mark all facets from a zone
  subroutine bc_mark_zone(this, bc_zone)
    class(bc_t), intent(inout) :: this
    type(zone_t), intent(inout) :: bc_zone
    integer :: i
    write(*,*) bc_zone%size
    do i = 1, bc_zone%size
       call this%marked_facet%push(bc_zone%facet_el(i))
    end do
  end subroutine bc_mark_zone

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

  !> Initialize a list of boundary conditions
  subroutine bc_list_init(bclst, size)
    type(bc_list_t), intent(inout), target :: bclst
    integer, optional :: size
    integer :: n, i

    call bc_list_free(bclst)

    if (present(size)) then
       n = size
    else
       n = 1
    end if

    allocate(bclst%bc(n))

    do i = 1, n
       bclst%bc(i)%bcp => null()
    end do

    bclst%n = 0
    bclst%size = n
        
  end subroutine bc_list_init

  !> Deallocate a list of boundary conditions
  !! @note This will only nullify all pointers, not deallocate any
  !! conditions pointed to by the list
  subroutine bc_list_free(bclst)
    type(bc_list_t), intent(inout) :: bclst

    if (allocated(bclst%bc)) then
       deallocate(bclst%bc)
    end if

    bclst%n = 0
    bclst%size = 0
    
  end subroutine bc_list_free

  !> Add a condition to a list of boundary conditions
  subroutine bc_list_add(bclst, bc)
    type(bc_list_t), intent(inout) :: bclst
    class(bc_t), intent(inout), target :: bc
    type(bcp_t), allocatable :: tmp(:)
    integer :: i 

    if (bclst%n .ge. bclst%size) then
       allocate(tmp(bclst%size * 2))
       tmp(1:bclst%size) = bclst%bc
       call move_alloc(tmp, bclst%bc)
    end if

    bclst%n = bclst%n + 1
    bclst%bc(bclst%n)%bcp => bc
    
  end subroutine bc_list_add

  !> Apply a list of boundary conditions
  subroutine bc_list_apply(bclst, x, n)
    type(bc_list_t), intent(inout) :: bclst
    integer, intent(in) :: n
    real(kind=dp), intent(inout),  dimension(n) :: x
    integer :: i

    do i = 1, bclst%n
       call bclst%bc(i)%bcp%apply(x, n)
    end do

  end subroutine bc_list_apply
   
  
end module bc
