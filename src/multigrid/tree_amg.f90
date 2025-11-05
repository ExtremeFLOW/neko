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
  use num_types, only : rp
  use utils, only : neko_error
  use math, only : rzero, col2
  use device_math , only : device_rzero, device_col2, device_masked_atomic_reduction_0, &
       device_masked_gather_copy_0, device_cfill
  use coefs, only : coef_t
  use mesh, only : mesh_t
  use space, only : space_t
  use ax_product, only: ax_t
  use bc_list, only: bc_list_t
  use gather_scatter, only : gs_t, GS_OP_ADD
  use device, only: device_map, device_free, device_deassociate, &
       device_stream_wait_event, glb_cmd_queue, glb_cmd_event
  use neko_config, only: NEKO_BCKND_DEVICE
  use, intrinsic :: iso_c_binding
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
   contains
     procedure, pass(this) :: free => node_free
  end type tamg_node_t

  !> Type for storing TreeAMG level information
  type, private :: tamg_lvl_t
     integer :: lvl = -1 !< The level id
     integer :: nnodes = 0 !< number of nodes on the level
     type(tamg_node_t), allocatable :: nodes(:) !< TreeAMG tree nodes on the level
     integer :: fine_lvl_dofs = 0 !< Number of dofs on the level(TODO:sum of dofs on each node?)
     real(kind=rp), allocatable :: wrk_in(:) !< Work vector for data coming into the level
     type(c_ptr) :: wrk_in_d = C_NULL_PTR
     real(kind=rp), allocatable :: wrk_out(:) !< Work vector for data leaving the level
     type(c_ptr) :: wrk_out_d = C_NULL_PTR
     integer, allocatable :: map_finest2lvl(:)
     type(c_ptr) :: map_finest2lvl_d = C_NULL_PTR
     !--!
     integer, allocatable :: map_f2c(:)
     type(c_ptr) :: map_f2c_d = C_NULL_PTR
   contains
     procedure, pass(this) :: free => lvl_free
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
     procedure, pass(this) :: free => tamg_free
     procedure, pass(this) :: matvec => tamg_matvec
     procedure, pass(this) :: matvec_impl => tamg_matvec_impl
     procedure, pass(this) :: interp_f2c => tamg_restriction_operator
     procedure, pass(this) :: interp_c2f => tamg_prolongation_operator
     procedure, pass(this) :: interp_f2c_d => tamg_device_restriction_operator
     procedure, pass(this) :: interp_c2f_d => tamg_device_prolongation_operator
     procedure, pass(this) :: device_matvec => tamg_device_matvec_flat_impl
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
       allocate( this%lvl(i)%map_finest2lvl( 0:coef%dof%size() ))
       if (NEKO_BCKND_DEVICE .eq. 1) then
          call device_map(this%lvl(i)%map_finest2lvl, this%lvl(i)%map_finest2lvl_d, coef%dof%size()+1)
       end if
    end do

  end subroutine tamg_init

  !> deallocate tamg hierarchy
  subroutine tamg_free(this)
    class(tamg_hierarchy_t), intent(inout) :: this
    integer :: i
    if (allocated(this%lvl)) then
       ! using size() instead of this%nlvls since
       ! this%nlvls may be less than the allocated number of levels due to early
       ! termination of aggregation
       do i = 1, size(this%lvl)
          call this%lvl(i)%free()
       end do
       deallocate(this%lvl)
    end if
    nullify(this%ax)
    nullify(this%msh)
    nullify(this%Xh)
    nullify(this%coef)
    nullify(this%gs_h)
    nullify(this%blst)
  end subroutine tamg_free

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
    allocate( tamg_lvl%map_f2c(0:ndofs) )
    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_map(tamg_lvl%map_f2c, tamg_lvl%map_f2c_d, ndofs+1)
    end if

    tamg_lvl%fine_lvl_dofs = ndofs
    allocate( tamg_lvl%wrk_in( ndofs ) )
    allocate( tamg_lvl%wrk_out( ndofs ) )
    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_map(tamg_lvl%wrk_in, tamg_lvl%wrk_in_d, ndofs)
       call device_cfill(tamg_lvl%wrk_in_d, 0.0_rp, ndofs)
       call device_map(tamg_lvl%wrk_out, tamg_lvl%wrk_out_d, ndofs)
       call device_cfill(tamg_lvl%wrk_out_d, 0.0_rp, ndofs)
    end if
  end subroutine tamg_lvl_init

  !> deallocate tamg level
  subroutine lvl_free(this)
    class(tamg_lvl_t), intent(inout) :: this
    integer :: i

    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_deassociate(this%wrk_in)
       call device_deassociate(this%wrk_out)
       call device_deassociate(this%map_f2c)
       call device_deassociate(this%map_finest2lvl)
    end if
    if (allocated(this%nodes)) then
       do i = 1, this%nnodes
          call this%nodes(i)%free()
       end do
       deallocate(this%nodes)
    end if
    if (allocated(this%wrk_in)) then
       deallocate(this%wrk_in)
    end if
    if (allocated(this%wrk_out)) then
       deallocate(this%wrk_out)
    end if
    if (allocated(this%map_f2c)) then
       deallocate(this%map_f2c)
    end if
    if (allocated(this%map_finest2lvl)) then
       deallocate(this%map_finest2lvl)
    end if
    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_free(this%wrk_in_d)
       call device_free(this%wrk_out_d)
       call device_free(this%map_f2c_d)
       call device_free(this%map_finest2lvl_d)
    end if
    this%nnodes = 0
    this%lvl = -1
    this%fine_lvl_dofs = 0
  end subroutine lvl_free

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
    node%interp_r = 1.0_rp
    node%interp_p = 1.0_rp
  end subroutine tamg_node_init

  !> deallocate tamg tree node
  subroutine node_free(this)
    class(tamg_node_t), intent(inout) :: this

    if (allocated(this%dofs)) then
       deallocate(this%dofs)
    end if
    if (allocated(this%interp_r)) then
       deallocate(this%interp_r)
    end if
    if (allocated(this%interp_p)) then
       deallocate(this%interp_p)
    end if
  end subroutine node_free

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
    !call this%matvec_impl(vec_out, vec_in, this%nlvls, lvl_out)
    call tamg_matvec_flat_impl(this, vec_out, vec_in, this%nlvls, lvl_out)
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

    if (lvl .eq. 0) then !> isleaf true
       !> If on finest level, pass to neko ax_t matvec operator
       n = size(vec_in)
       !> Call local finite element assembly
       call this%gs_h%op(vec_in, n, GS_OP_ADD)
       call col2( vec_in, this%coef%mult, n)

       call this%ax%compute(vec_out, vec_in, this%coef, this%msh, this%Xh)
       call this%gs_h%op(vec_out, n, GS_OP_ADD)
       call this%blst%apply(vec_out, n)

       if (lvl_out .ne. 0) then
          call col2(vec_out, this%coef%mult, n)
       end if
       !>
    else !> pass down through hierarchy
       if (lvl_out .ge. lvl) then
          !> lvl is finer than desired output
          !> project input vector to finer grid
          associate( wrk_in => this%lvl(lvl)%wrk_in, wrk_out => this%lvl(lvl)%wrk_out)
            n = this%lvl(lvl)%fine_lvl_dofs
            call rzero(wrk_in, n)
            call rzero(vec_out, this%lvl(lvl)%nnodes)
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
    n = this%lvl(1)%fine_lvl_dofs
    if (lvl .eq. 0) then !> isleaf true
       call this%ax%compute(vec_out, vec_in, this%coef, this%msh, this%Xh)
       call this%gs_h%op(vec_out, n, GS_OP_ADD)
       call this%blst%apply(vec_out, n)
    else !> pass down through hierarchy
       associate( wrk_in => this%lvl(1)%wrk_in, wrk_out => this%lvl(1)%wrk_out)
         !> Map input level to finest level
         do i = 1, n
            cdof = this%lvl(lvl)%map_finest2lvl(i)
            wrk_in(i) = vec_in( cdof )
         end do

         !> Average on overlapping dofs
         call this%gs_h%op(wrk_in, n, GS_OP_ADD)
         call col2( wrk_in, this%coef%mult, n)
         call this%blst%apply(wrk_in, n)

         !> Finest level matvec (Call local finite element assembly)
         call this%ax%compute(wrk_out, wrk_in, this%coef, this%msh, this%Xh)
         call this%gs_h%op(wrk_out, n, GS_OP_ADD)
         call this%blst%apply(wrk_out, n)

         call col2(wrk_out, this%coef%mult, n)

         !> Map finest level matvec back to output level
         call rzero(vec_out, this%lvl(lvl)%nnodes)
         do i = 1, n
            cdof = this%lvl(lvl)%map_finest2lvl(i)
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
    if (lvl-1 .eq. 0) then
       call col2(vec_in, this%coef%mult, this%lvl(lvl)%fine_lvl_dofs)
    end if
    do n = 1, this%lvl(lvl)%nnodes
       associate (node => this%lvl(lvl)%nodes(n))
         do i = 1, node%ndofs
            vec_out( node%gid ) = vec_out( node%gid ) + vec_in( node%dofs(i) ) * node%interp_r( i )
         end do
       end associate
    end do
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
    if (lvl-1 .eq. 0) then
       call this%gs_h%op(vec_out, this%lvl(lvl)%fine_lvl_dofs, GS_OP_ADD)
       call col2(vec_out, this%coef%mult, this%lvl(lvl)%fine_lvl_dofs)
       call this%blst%apply(vec_out, n)
    end if
  end subroutine tamg_prolongation_operator


  subroutine tamg_device_matvec_flat_impl(this, vec_out, vec_in, vec_out_d, vec_in_d, lvl_out)
    class(tamg_hierarchy_t), intent(inout) :: this
    real(kind=rp), intent(inout) :: vec_out(:)
    real(kind=rp), intent(inout) :: vec_in(:)
    type(c_ptr) :: vec_out_d
    type(c_ptr) :: vec_in_d
    integer, intent(in) :: lvl_out
    integer :: i, n, cdof, lvl

    lvl = lvl_out
    n = this%lvl(1)%fine_lvl_dofs
    if (lvl .eq. 0) then !> isleaf true
       call this%ax%compute(vec_out, vec_in, this%coef, this%msh, this%Xh)
       call this%gs_h%op(vec_out, n, GS_OP_ADD, glb_cmd_event)
       call device_stream_wait_event(glb_cmd_queue, glb_cmd_event, 0)
       call this%blst%apply(vec_out, n)
    else !> pass down through hierarchy

       associate( wrk_in_d => this%lvl(1)%wrk_in_d, wrk_out_d => this%lvl(1)%wrk_out_d)
         !> Map input level to finest level
         call device_masked_gather_copy_0(wrk_in_d, vec_in_d, this%lvl(lvl)%map_finest2lvl_d, this%lvl(lvl)%nnodes, n)
         !> Average on overlapping dofs
         call this%gs_h%op(this%lvl(1)%wrk_in, n, GS_OP_ADD, glb_cmd_event)
         call device_stream_wait_event(glb_cmd_queue, glb_cmd_event, 0)
         call device_col2( wrk_in_d, this%coef%mult_d, n)
         call this%blst%apply(this%lvl(1)%wrk_in, n)

         !> Finest level matvec (Call local finite element assembly)
         call this%ax%compute(this%lvl(1)%wrk_out, this%lvl(1)%wrk_in, this%coef, this%msh, this%Xh)
         call this%gs_h%op(this%lvl(1)%wrk_out, n, GS_OP_ADD, glb_cmd_event)
         call device_stream_wait_event(glb_cmd_queue, glb_cmd_event, 0)
         call this%blst%apply(this%lvl(1)%wrk_out, n)

         call device_col2( wrk_out_d, this%coef%mult_d, n)

         !> Map finest level matvec back to output level
         call device_rzero(vec_out_d, this%lvl(lvl)%nnodes)
         call device_masked_atomic_reduction_0(vec_out_d, wrk_out_d, this%lvl(lvl)%map_finest2lvl_d, this%lvl(lvl)%nnodes, n)
       end associate

    end if
  end subroutine tamg_device_matvec_flat_impl

  subroutine tamg_device_restriction_operator(this, vec_out_d, vec_in_d, lvl)
    class(tamg_hierarchy_t), intent(inout) :: this
    type(c_ptr) :: vec_out_d
    type(c_ptr) :: vec_in_d
    integer, intent(in) :: lvl
    integer :: i, n, m
    n = this%lvl(lvl)%nnodes
    m = this%lvl(lvl)%fine_lvl_dofs
    if (lvl-1 .eq. 0) then
       call device_col2(vec_in_d, this%coef%mult_d, m)
    end if
    call device_rzero(vec_out_d, n)
    call device_masked_atomic_reduction_0(vec_out_d, vec_in_d, this%lvl(lvl)%map_f2c_d, n, m)
  end subroutine tamg_device_restriction_operator

  subroutine tamg_device_prolongation_operator(this, vec_out_d, vec_in_d, lvl, vec_out)
    class(tamg_hierarchy_t), intent(inout) :: this
    real(kind=rp), intent(inout) :: vec_out(:)
    type(c_ptr) :: vec_out_d
    type(c_ptr) :: vec_in_d
    integer, intent(in) :: lvl
    integer :: i, n, m
    n = this%lvl(lvl)%nnodes
    m = this%lvl(lvl)%fine_lvl_dofs
    call device_masked_gather_copy_0(vec_out_d, vec_in_d, this%lvl(lvl)%map_f2c_d, n, m)
    if (lvl-1 .eq. 0) then
       call this%gs_h%op(vec_out, m, GS_OP_ADD, glb_cmd_event)
       call device_stream_wait_event(glb_cmd_queue, glb_cmd_event, 0)
       call device_col2( vec_out_d, this%coef%mult_d, m)
       call this%blst%apply( vec_out, n)
    end if
  end subroutine tamg_device_prolongation_operator

end module tree_amg
