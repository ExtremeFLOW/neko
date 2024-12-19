! Copyright (c) 2024, The Neko Authors
! All rights reserved.
!
! Redistribution and use in source and binary forms, with or without
! modification, are permitted provided that the following conditions
! are met:
!
!   * Redistributions of source code must retain the above copyright
!     notice, this list of conditions and the following disclaimer.
!
!   * Redistributions in binary form must reproduce the above
!     copyright notice, this list of conditions and the following
!     disclaimer in the documentation and/or other materials provided
!     with the distribution.
!
!   * Neither the name of the authors nor the names of its
!     contributors may be used to endorse or promote products derived
!     from this software without specific prior written permission.
!
! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
! "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
! LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
! FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
! COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
! INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
! BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
! LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
! CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
! LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
! ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
! POSSIBILITY OF SUCH DAMAGE.
!
!> Implements the base type for TreeAMG hierarchy structure.
module tree_amg
  use num_types
  use utils
  use math
  use coefs, only : coef_t
  use mesh, only : mesh_t
  use space, only : space_t
  use ax_product, only: ax_t
  use bc_list, only: bc_list_t
  use gather_scatter, only : gs_t, GS_OP_ADD
  implicit none
  private

  !> Type for storing TreeAMG tree node information
  type, private :: tamg_node_t
    logical :: isleaf = .true. !< Is the node a leaf node
    integer :: gid = -1 !< The gid of the node TODO: relative to what? MPI or true global
    integer :: lvl = -1 !< The hierarchy level on which the node lives
    integer :: ndofs = 0 !< The number of dofs on the node
    integer, allocatable :: dofs(:) !< The dofs on the node
    real(kind=rp) :: xyz(3) !< Coordinates of the node
    real(kind=rp), allocatable :: interp_r(:) !< Resriciton factors from dofs to node gid
    real(kind=rp), allocatable :: interp_p(:) !< Prolongation factors from node gid to dofs
  end type tamg_node_t

  !> Type for storing TreeAMG level information
  type, private :: tamg_lvl_t
    integer :: lvl = -1 !< The level id
    integer :: nnodes = 0 !< number of nodes on the level
    type(tamg_node_t), allocatable :: nodes(:) !< TreeAMG tree nodes on the level
    integer :: fine_lvl_dofs = 0 !< Number of dofs on the level(TODO:sum of dofs on each node?)
    real(kind=rp), allocatable :: wrk_in(:) !< Work vector for data coming into the level
    real(kind=rp), allocatable :: wrk_out(:) !< Work vector for data leaving the level
    integer, allocatable :: map_f2c_dof(:)
    integer, allocatable :: map_c2f_dof(:)
    !--!
    integer, allocatable :: nodes_ptr(:)
    integer, allocatable :: nodes_gid(:)
    integer, allocatable :: nodes_dofs(:)
    integer, allocatable :: nodes_gids(:)
    ! could make another array of the same size of nodes_dofs
    ! that stores the parent node gid information
    ! (similar to nodes_gid that stores the gid of each node)
    ! then some loops can be simplified to a single loop
    ! of len(nodes_dofs) instead of going through each node
    ! and looping through nodes_ptr(i) to nodes_ptr(i+1)-1
  end type tamg_lvl_t

  !> Type for a TreeAMG hierarchy
  type, public :: tamg_hierarchy_t
    !> Number of AMG levels in the hierarchy
    integer :: nlvls
    !> Levels of the hierarchy
    type(tamg_lvl_t), allocatable :: lvl(:)

    !> Things needed to do finest level matvec
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

  !> Initialization of TreeAMG hierarchy
  !! @param ax Finest level matvec operator
  !! @param Xh Finest level field
  !! @param coef Finest level coeff thing
  !! @param msh Finest level mesh information
  !! @param gs_h Finest level gather scatter operator
  !! @param nlvls Number of levels for the TreeAMG hierarchy
  !! @param blst Finest level BC list
  subroutine tamg_init(this, ax, Xh, coef, msh, gs_h, nlvls, blst)
    class(tamg_hierarchy_t), target, intent(inout) :: this
    class(ax_t), target, intent(in) :: ax
    type(space_t),target, intent(in) :: Xh
    type(coef_t), target, intent(in) :: coef
    type(mesh_t), target, intent(in) :: msh
    type(gs_t), target, intent(in) :: gs_h
    integer, intent(in) :: nlvls
    type(bc_list_t), target, intent(in) :: blst
    integer :: i, n

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

    do i = 1, nlvls
      allocate( this%lvl(i)%map_f2c_dof( coef%dof%size() ))
      allocate( this%lvl(i)%map_c2f_dof( coef%dof%size() ))
    end do

  end subroutine tamg_init

  !> Initialization of a TreeAMG level
  !! @param tamg_lvl The TreeAMG level
  !! @param lvl The level id
  !! @param nnodes Number of nodes on the level (number of aggregates)
  !! @param ndofs  Number of dofs on the level
  subroutine tamg_lvl_init(tamg_lvl, lvl, nnodes, ndofs)
    type(tamg_lvl_t), intent(inout) :: tamg_lvl
    integer, intent(in) :: lvl
    integer, intent(in) :: nnodes
    integer, intent(in) :: ndofs

    tamg_lvl%lvl = lvl
    tamg_lvl%nnodes = nnodes
    allocate( tamg_lvl%nodes(tamg_lvl%nnodes) )
    allocate( tamg_lvl%nodes_ptr(tamg_lvl%nnodes+1) )
    allocate( tamg_lvl%nodes_gid(tamg_lvl%nnodes) )
    allocate( tamg_lvl%nodes_dofs(ndofs) )
    allocate( tamg_lvl%nodes_gids(ndofs) )

    tamg_lvl%fine_lvl_dofs = ndofs
    allocate( tamg_lvl%wrk_in( ndofs ) )
    allocate( tamg_lvl%wrk_out( ndofs ) )
  end subroutine tamg_lvl_init

  !> Initialization of a TreeAMG tree node
  !! @param node The TreeAMG tree node
  !! @param gid The gid for the node
  !! @param ndofs Number of dofs in the node
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

  !> Wrapper for matrix vector product using the TreeAMG hierarchy structure
  !> b=Ax done as vec_out = A * vec_in
  !! @param vec_out Result of Ax
  !! @param vec_in Vector to be multiplied by linear system. A * vec_in
  !! @param lvl_out Level of the TreeAMG hierarchy on which the matvec is done. This
  !!                also specifies the hieararchy level of the incoming vector
  recursive subroutine tamg_matvec(this, vec_out, vec_in, lvl_out)
    class(tamg_hierarchy_t), intent(inout) :: this
    real(kind=rp), intent(inout) :: vec_out(:)
    real(kind=rp), intent(inout) :: vec_in(:)
    integer, intent(in) :: lvl_out
    integer :: i, n, e
    call this%matvec_impl(vec_out, vec_in, this%nlvls, lvl_out)
    !call tamg_matvec_flat_impl(this, vec_out, vec_in, this%nlvls, lvl_out)
  end subroutine tamg_matvec

  !> Matrix vector product using the TreeAMG hierarchy structure
  !> b=Ax done as vec_out = A * vec_in
  !> This is done on a level by level basis
  !! @param vec_out The vector to be returned by level lvl
  !! @param vec_in The vector pased to level lvl
  !! @param lvl The current level of the matvec (wrt tree traversal)
  !! @param lvl_out Level of the TreeAMG hierarchy on which the matvec is output.
  recursive subroutine tamg_matvec_impl(this, vec_out, vec_in, lvl, lvl_out)
    class(tamg_hierarchy_t), intent(inout) :: this
    real(kind=rp), intent(inout) :: vec_out(:)
    real(kind=rp), intent(inout) :: vec_in(:)
    integer, intent(in) :: lvl
    integer, intent(in) :: lvl_out
    integer :: i, n, e

    vec_out = 0d0

    if (lvl .eq. 0) then !> isleaf true
      !> If on finest level, pass to neko ax_t matvec operator
      n = size(vec_in)
      !> Call local finite element assembly
      call this%gs_h%op(vec_in, n, GS_OP_ADD)
      call col2( vec_in, this%coef%mult(1,1,1,1), n)
      !>
      call this%ax%compute(vec_out, vec_in, this%coef, this%msh, this%Xh)
      !>
      call this%gs_h%op(vec_out, n, GS_OP_ADD)
      call this%blst%apply(vec_out, n)
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
        call neko_error("TAMG: matvec level numbering problem.")
      end if
    end if
  end subroutine tamg_matvec_impl


  !> Ignore this. For piecewise constant, can create index map directly to finest level
  recursive subroutine tamg_matvec_flat_impl(this, vec_out, vec_in, lvl_blah, lvl_out)
    class(tamg_hierarchy_t), intent(inout) :: this
    real(kind=rp), intent(inout) :: vec_out(:)
    real(kind=rp), intent(inout) :: vec_in(:)
    integer, intent(in) :: lvl_blah
    integer, intent(in) :: lvl_out
    integer :: i, n, cdof, lvl

    lvl = lvl_out
    if (lvl .eq. 0) then !> isleaf true
      !> If on finest level, pass to neko ax_t matvec operator
      n = size(vec_in)
      !> Call local finite element assembly
      call this%gs_h%op(vec_in, n, GS_OP_ADD)
      call col2( vec_in, this%coef%mult(1,1,1,1), n)
      !>
      call this%ax%compute(vec_out, vec_in, this%coef, this%msh, this%Xh)
      !>
      call this%gs_h%op(vec_out, n, GS_OP_ADD)
      call this%blst%apply(vec_out, n)
      !>
    else !> pass down through hierarchy

      associate( wrk_in => this%lvl(1)%wrk_in, wrk_out => this%lvl(1)%wrk_out)
      n = size(wrk_in)
      wrk_out = 0d0
      vec_out = 0d0

      !> Map input level to finest level
      do i = 1, n
        cdof = this%lvl(lvl)%map_f2c_dof(i)
        wrk_in(i) = vec_in( cdof )
      end do

      !> Average on overlapping dofs
      call this%gs_h%op(wrk_in, n, GS_OP_ADD)
      call col2( wrk_in, this%coef%mult(1,1,1,1), n)
      !> Finest level matvec (Call local finite element assembly)
      call this%ax%compute(wrk_out, wrk_in, this%coef, this%msh, this%Xh)
      !>
      call this%gs_h%op(wrk_out, n, GS_OP_ADD)
      call this%blst%apply(wrk_out, n)
      call bc_list_apply(this%blst, wrk_out, n)
      !>

      !> Map finest level matvec back to output level
      do i = 1, n
        cdof = this%lvl(lvl)%map_f2c_dof(i)
        vec_out(cdof) = vec_out(cdof) + wrk_out( i )
      end do
      end associate

    end if
  end subroutine tamg_matvec_flat_impl


  !> Restriction operator for TreeAMG. vec_out = R * vec_in
  !! @param vec_out The vector to be returned. On level lvl
  !! @param vec_in The vector pased into operator. On level lvl-1
  !! @param lvl The target level of the returned vector after restrction (wrt tree traversal)
  subroutine tamg_restriction_operator(this, vec_out, vec_in, lvl)
    class(tamg_hierarchy_t), intent(inout) :: this
    real(kind=rp), intent(inout) :: vec_out(:)
    real(kind=rp), intent(inout) :: vec_in(:)
    integer, intent(in) :: lvl
    integer :: i, n, node_start, node_end, node_id

    vec_out = 0d0
    do n = 1, this%lvl(lvl)%nnodes
      associate (node => this%lvl(lvl)%nodes(n))
      do i = 1, node%ndofs
        vec_out( node%gid ) = vec_out( node%gid ) + vec_in( node%dofs(i) ) * node%interp_r( i )
      end do
      end associate
    end do
    !do n = 1, this%lvl(lvl)%nnodes
    !  node_start = this%lvl(lvl)%nodes_ptr(n)
    !  node_end   = this%lvl(lvl)%nodes_ptr(n+1)-1
    !  node_id    = this%lvl(lvl)%nodes_gid(n)
    !  do i = node_start, node_end
    !    vec_out( node_id ) = vec_out( node_id ) + &
    !      vec_in( this%lvl(lvl)%nodes_dofs(i) )
    !  end do
    !end do
  end subroutine tamg_restriction_operator

  !> Prolongation operator for TreeAMG. vec_out = P * vec_in
  !! @param vec_out The vector to be returned. On level lvl
  !! @param vec_in The vector pased into operator. On level lvl-1
  !! @param lvl The target level of the returned vector after prolongation (wrt tree traversal)
  subroutine tamg_prolongation_operator(this, vec_out, vec_in, lvl)
    class(tamg_hierarchy_t), intent(inout) :: this
    real(kind=rp), intent(inout) :: vec_out(:)
    real(kind=rp), intent(inout) :: vec_in(:)
    integer, intent(in) :: lvl
    integer :: i, n, node_start, node_end, node_id

    vec_out = 0d0
    do n = 1, this%lvl(lvl)%nnodes
      associate (node => this%lvl(lvl)%nodes(n))
      do i = 1, node%ndofs
        vec_out( node%dofs(i) ) = vec_out( node%dofs(i) ) + vec_in( node%gid ) * node%interp_p( i )
      end do
      end associate
    end do
    !do n = 1, this%lvl(lvl)%nnodes
    !  node_start = this%lvl(lvl)%nodes_ptr(n)
    !  node_end   = this%lvl(lvl)%nodes_ptr(n+1)-1
    !  node_id    = this%lvl(lvl)%nodes_gid(n)
    !  do i = node_start, node_end
    !    vec_out( this%lvl(lvl)%nodes_dofs(i) ) = vec_out( this%lvl(lvl)%nodes_dofs(i) ) + &
    !      vec_in(node_id)
    !  end do
    !end do
  end subroutine tamg_prolongation_operator

end module tree_amg
