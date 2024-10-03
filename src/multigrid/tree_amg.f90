module tree_amg
  use num_types
  use utils
  use math
  use coefs, only : coef_t
  use mesh, only : mesh_t
  use space, only : space_t
  use ax_product, only: ax_t
  use bc, only: bc_list_t, bc_list_apply
  use gather_scatter, only : gs_t, GS_OP_ADD
  implicit none
  private

  type, private :: tamg_node_t
    logical :: isleaf = .true.
    integer :: gid = -1
    integer :: lvl = -1
    integer :: ndofs = 0
    integer, allocatable :: dofs(:)
    real(kind=rp) :: xyz(3)
    real(kind=rp), allocatable :: interp_r(:)
    real(kind=rp), allocatable :: interp_p(:)
  end type tamg_node_t

  type, private :: tamg_lvl_t
    integer :: lvl = -1
    integer :: nnodes = 0
    type(tamg_node_t), allocatable :: nodes(:)
    integer :: fine_lvl_dofs = 0
    real(kind=rp), allocatable :: wrk_in(:)
    real(kind=rp), allocatable :: wrk_out(:)
  end type tamg_lvl_t

  type, public :: tamg_hierarchy_t
    integer :: nlvls
    type(tamg_lvl_t), allocatable :: lvl(:)

    !Things needed to do finest level matvec
    class(ax_t), pointer :: ax
    type(mesh_t), pointer :: msh
    type(space_t), pointer :: Xh
    type(coef_t), pointer :: coef
    type(gs_t), pointer :: gs_h
    type(bc_list_t), pointer :: blst

  contains
    procedure, pass(this) :: init => tamg_init
    procedure, pass(this) :: matvec => tamg_matvec
    procedure, pass(this) :: matvec_impl => tamg_matvec_impl
    procedure, pass(this) :: interp_f2c => tamg_restriction_operator
    procedure, pass(this) :: interp_c2f => tamg_prolongation_operator
  end type tamg_hierarchy_t

  public tamg_lvl_init, tamg_node_init

contains

  subroutine tamg_init(this, ax, Xh, coef, msh, gs_h, nlvls, blst)
    class(tamg_hierarchy_t), target, intent(inout) :: this
    class(ax_t), target, intent(in) :: ax
    type(space_t),target, intent(in) :: Xh
    type(coef_t), target, intent(in) :: coef
    type(mesh_t), target, intent(in) :: msh
    type(gs_t), target, intent(in) :: gs_h
    integer, intent(in) :: nlvls
    type(bc_list_t), target, intent(in) :: blst

    this%ax => ax
    this%msh => msh
    this%Xh => Xh
    this%coef => coef
    this%gs_h => gs_h
    this%blst => blst

    if (nlvls .lt. 2) then
      call neko_error("Need to request at least two multigrid levels.")
    end if

    this%nlvls = nlvls
    allocate( this%lvl(this%nlvls) )


  end subroutine tamg_init

  subroutine tamg_lvl_init(tamg_lvl, lvl, nnodes)
    type(tamg_lvl_t), intent(inout) :: tamg_lvl
    integer, intent(in) :: lvl
    integer, intent(in) :: nnodes

    tamg_lvl%lvl = lvl
    tamg_lvl%nnodes = nnodes
    allocate( tamg_lvl%nodes(tamg_lvl%nnodes) )
  end subroutine tamg_lvl_init

  subroutine tamg_node_init(node, gid, ndofs)
    type(tamg_node_t), intent(inout) :: node
    integer, intent(in) :: gid
    integer, intent(in) :: ndofs

    node%gid = gid
    node%ndofs = ndofs
    allocate( node%dofs( node%ndofs) )
    node%dofs = -1
    allocate( node%interp_r( node%ndofs) )
    allocate( node%interp_p( node%ndofs) )
    node%interp_r = 1.0
    node%interp_p = 1.0
  end subroutine tamg_node_init

  recursive subroutine tamg_matvec(this, vec_out, vec_in, lvl_out)
    class(tamg_hierarchy_t), intent(inout) :: this
    real(kind=rp), intent(inout) :: vec_out(:)
    real(kind=rp), intent(inout) :: vec_in(:)
    integer, intent(in) :: lvl_out
    integer :: i, n, e
    call this%matvec_impl(vec_out, vec_in, this%nlvls, lvl_out)
    !call tamg_matvec_flat_impl(this, vec_out, vec_in, this%nlvls, lvl_out)
  end subroutine tamg_matvec

  recursive subroutine tamg_matvec_impl(this, vec_out, vec_in, lvl, lvl_out)
    class(tamg_hierarchy_t), intent(inout) :: this
    real(kind=rp), intent(inout) :: vec_out(:)
    real(kind=rp), intent(inout) :: vec_in(:)
    integer, intent(in) :: lvl
    integer, intent(in) :: lvl_out
    integer :: i, n, e

    vec_out = 0d0

    if (lvl .eq. 0) then !> isleaf true
      n = size(vec_in)
      !> Call local finite element assembly
      call this%gs_h%op(vec_in, n, GS_OP_ADD)
      do i = 1, n
        vec_in(i) = vec_in(i) * this%coef%mult(i,1,1,1)
      end do
      !>
      call this%ax%compute(vec_out, vec_in, this%coef, this%msh, this%Xh)
      !>
      call this%gs_h%op(vec_out, n, GS_OP_ADD)
      call bc_list_apply(this%blst, vec_out, n)
      !>
    else !> pass down through hierarchy
      if (lvl_out .ge. lvl) then
        !> lvl is finer than desired output
        !> project input vector to finer grid
        associate( wrk_in => this%lvl(lvl)%wrk_in, wrk_out => this%lvl(lvl)%wrk_out)
        wrk_in = 0d0
        wrk_out = 0d0
        do n = 1, this%lvl(lvl)%nnodes
          associate (node => this%lvl(lvl)%nodes(n))
          do i = 1, node%ndofs
            wrk_in( node%dofs(i) ) = wrk_in( node%dofs(i) ) + vec_in( node%gid ) * node%interp_p( i )
          end do
          end associate
        end do

        call this%matvec_impl(wrk_out, wrk_in, lvl-1, lvl_out)

        !> restrict to coarser grid
        do n = 1, this%lvl(lvl)%nnodes
          associate (node => this%lvl(lvl)%nodes(n))
          do i = 1, node%ndofs
            vec_out( node%gid ) = vec_out(node%gid ) + wrk_out( node%dofs(i) ) * node%interp_r( i )
          end do
          end associate
        end do
        end associate
      else if (lvl_out .lt. lvl) then
        !> lvl is coarser then desired output. Continue down tree
        call this%matvec_impl(vec_out, vec_in, lvl-1, lvl_out)
      else
        print *, "ERROR"
      end if
    end if
  end subroutine tamg_matvec_impl


  recursive subroutine tamg_matvec_flat_impl(this, vec_out, vec_in, lvl_blah, lvl_out)
    class(tamg_hierarchy_t), intent(inout) :: this
    real(kind=rp), intent(inout) :: vec_out(:)
    real(kind=rp), intent(inout) :: vec_in(:)
    integer, intent(in) :: lvl_blah
    integer, intent(in) :: lvl_out
    integer :: i, n, e, nn, lvl

    vec_out = 0d0
    do lvl = 1, lvl_out+1
      associate( wrk_in => this%lvl(lvl)%wrk_in, wrk_out => this%lvl(lvl)%wrk_out)
        wrk_in = 0d0
        wrk_out = 0d0
      end associate
    end do

    !> copy to make things easier to think about
    do n = 1, this%lvl(lvl_out+1)%nnodes
      do i = 1, this%lvl(lvl_out+1)%nodes(n)%ndofs
        associate( wrk_in => this%lvl(lvl_out+1)%wrk_in, node => this%lvl(lvl_out+1)%nodes(n))
          wrk_in( node%dofs(i) ) = vec_in( node%dofs(i) )
        end associate
      end do!i
    end do!n

    if(lvl_out .gt. 0) then
      do n = 1, this%lvl(lvl_out)%nnodes!> this loop is independent
        do lvl = lvl_out, 1, -1

          associate( wrk_in => this%lvl(lvl)%wrk_in, wrk_out => this%lvl(lvl)%wrk_out)
          do nn = 1, this%lvl(lvl)%nnodes
            associate (node => this%lvl(lvl)%nodes(nn))
            do i = 1, node%ndofs
              wrk_in( node%dofs(i) ) = wrk_in( node%dofs(i) ) + this%lvl(lvl+1)%wrk_in( node%gid ) * node%interp_p( i )
            end do!i
            end associate
          end do!nn
          end associate

        end do!lvl
      end do!n
    end if

    associate( wrk_in => this%lvl(1)%wrk_in, wrk_out => this%lvl(1)%wrk_out)
    !> Do finest level matvec
    n = size(wrk_in)
    !> Call local finite element assembly
    call this%gs_h%op(wrk_in, n, GS_OP_ADD)
    do i = 1, n
      wrk_in(i) = wrk_in(i) * this%coef%mult(i,1,1,1)
    end do
    !>
    call this%ax%compute(wrk_out, wrk_in, this%coef, this%msh, this%Xh)
    !>
    call this%gs_h%op(wrk_out, n, GS_OP_ADD)
    call bc_list_apply(this%blst, wrk_out, n)
    !>
    end associate


    if(lvl_out .gt. 0) then
      do n = 1, this%lvl(lvl_out)%nnodes!> this loop is independent
        do lvl = 2, lvl_out+1

          associate( wrk_in => this%lvl(lvl)%wrk_in, wrk_out => this%lvl(lvl)%wrk_out)
          do nn = 1, this%lvl(lvl)%nnodes
            associate (node => this%lvl(lvl)%nodes(nn))
            do i = 1, node%ndofs
              wrk_out( node%gid ) = wrk_out( node%gid ) + this%lvl(lvl-1)%wrk_out( node%dofs(i) ) * node%interp_r( i )
            end do!i
            end associate
          end do!nn
          end associate

        end do!lvl
      end do!n
    end if

    !> copy here to make things easier to think about
    do n = 1, this%lvl(lvl_out+1)%nnodes
      do i = 1, this%lvl(lvl_out+1)%nodes(n)%ndofs
        associate( wrk_out => this%lvl(lvl_out+1)%wrk_out, node => this%lvl(lvl_out+1)%nodes(n))
          vec_out( node%dofs(i) ) = wrk_out( node%dofs(i) )
        end associate
      end do!i
    end do!n
  end subroutine tamg_matvec_flat_impl


  subroutine tamg_restriction_operator(this, vec_out, vec_in, lvl)
    class(tamg_hierarchy_t), intent(inout) :: this
    real(kind=rp), intent(inout) :: vec_out(:)
    real(kind=rp), intent(inout) :: vec_in(:)
    integer, intent(in) :: lvl
    integer :: i, n

    vec_out = 0d0
    do n = 1, this%lvl(lvl)%nnodes
      associate (node => this%lvl(lvl)%nodes(n))
      do i = 1, node%ndofs
        vec_out( node%gid ) = vec_out( node%gid ) + vec_in( node%dofs(i) ) * node%interp_r( i )
      end do
      end associate
    end do
  end subroutine tamg_restriction_operator

  subroutine tamg_prolongation_operator(this, vec_out, vec_in, lvl)
    class(tamg_hierarchy_t), intent(inout) :: this
    real(kind=rp), intent(inout) :: vec_out(:)
    real(kind=rp), intent(inout) :: vec_in(:)
    integer, intent(in) :: lvl
    integer :: i, n

    vec_out = 0d0
    do n = 1, this%lvl(lvl)%nnodes
      associate (node => this%lvl(lvl)%nodes(n))
      do i = 1, node%ndofs
        vec_out( node%dofs(i) ) = vec_out( node%dofs(i) ) + vec_in( node%gid ) * node%interp_p( i )
      end do
      end associate
    end do
  end subroutine tamg_prolongation_operator

end module tree_amg
