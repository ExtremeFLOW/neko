! Copyright (c) 2019-2026, The Neko Authors
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
!> Implementation for the nonconforming faces/edges during communication
module gs_interp
  use num_types, only : i4, i8, rp
  use utils, only : neko_error
  use mesh_conn, only : mesh_conn_t
  use amr_interpolate, only : amr_interpolate_t, amr_nchildren
  use field, only : field_t
  use amr_restart_component, only : amr_restart_component_t

  implicit none
  private

  !> Type for face/edge interpolation
  type, public, abstract, extends(amr_restart_component_t) :: gs_interp_t
     !> Polynomial order + 1
     integer :: lx
     !> Number of local elements
     integer :: nel
     !> Mesh connectivity
     type(mesh_conn_t), pointer :: conn
     !> AMR interpolation arrays
     type(amr_interpolate_t) :: interpolate
     !> Is there any element with hanging objects
     logical :: ifhang
     !> Number of elements with hanging objects
     integer :: nhang_el
     !> List of  elements with hanging objects
     ! size of nhang_el
     integer, allocatable, dimension(:) :: hang_el
     !> Number of hanging edges
     integer :: nhang_edg
     !> List of hanging edges (position in the element)
     ! size of nhang_edg
     integer, allocatable, dimension(:) :: hang_edg
     !> List of hanging edges (position on the parent edge)
     ! size of nhang_edg
     integer, allocatable, dimension(:) :: hang_edg_pos
     !> Hanging edges offset
     ! size of nhang_el + 1
     integer, allocatable, dimension(:) :: hang_edg_off
     !> Number of hanging faces
     integer :: nhang_fcs
     !> List of hanging faces (position in the element)
     ! size of nhang_fcs
     integer, allocatable, dimension(:) :: hang_fcs
     !> List of hanging faces (position on the parent face)
     ! size of nhang_fcs
     integer, allocatable, dimension(:) :: hang_fcs_pos
     !> Hanging faces offset
     ! size of nhang_el + 1
     integer, allocatable, dimension(:) :: hang_fcs_off
   contains
     !> Initialise base type
     procedure, pass(this) :: init_base => gs_interp_init_base
     !> Free base type
     procedure, pass(this) :: free_base => gs_interp_free_base
     !> Initialise type
     procedure(gs_interp_init), pass(this), deferred :: init
     !> Initialise multiplicity arrays
     procedure(gs_interp_init_mult), pass(this), deferred :: init_mult
     !> Free type
     procedure(gs_interp_free), pass(this), deferred :: free
     !> Perform face/edge interpolation
     procedure(gs_interp_apply_fld), pass(this), deferred :: apply_j_fld
     procedure(gs_interp_apply_r4), pass(this), deferred :: apply_j_r4
     generic :: apply_j => apply_j_fld, apply_j_r4
     !> Perform inverse face/edge interpolation
     procedure(gs_interp_apply_fld), pass(this), deferred :: apply_ji_fld
     procedure(gs_interp_apply_r4), pass(this), deferred :: apply_ji_r4
     generic :: apply_ji => apply_ji_fld, apply_ji_r4
     !> Perform transposed face/edge interpolation
     procedure(gs_interp_apply_fld), pass(this), deferred :: apply_jt_fld
     procedure(gs_interp_apply_r4), pass(this), deferred :: apply_jt_r4
     generic :: apply_jt => apply_jt_fld, apply_jt_r4
     !> Zero children's nonconforming faces/edges
     procedure(gs_interp_apply_fld), pass(this), deferred :: zero_children_fld
     procedure(gs_interp_apply_r4), pass(this), deferred :: zero_children_r4
     generic :: zero_children => zero_children_fld, zero_children_r4
     !> Set children's nonconforming faces/edges
     procedure(gs_interp_set_fld), pass(this), deferred :: set_children_fld
     procedure(gs_interp_set_r4), pass(this), deferred :: set_children_r4
     generic :: set_children => set_children_fld, set_children_r4
     !> Remove multiplicity for H1
     procedure(gs_interp_apply_fld), pass(this), deferred :: remove_mult_h1_fld
     procedure(gs_interp_apply_r4), pass(this), deferred :: remove_mult_h1_r4
     generic :: remove_mult_h1 => remove_mult_h1_fld, remove_mult_h1_r4
     !> Remove multiplicity for J^T
     procedure(gs_interp_apply_fld), pass(this), deferred :: remove_mult_jt_fld
     procedure(gs_interp_apply_r4), pass(this), deferred :: remove_mult_jt_r4
     generic :: remove_mult_jt => remove_mult_jt_fld, remove_mult_jt_r4
     !> Remove multiplicity for J^-1
     procedure(gs_interp_apply_fld), pass(this), deferred :: remove_mult_ji_fld
     procedure(gs_interp_apply_r4), pass(this), deferred :: remove_mult_ji_r4
     generic :: remove_mult_ji => remove_mult_ji_fld, remove_mult_ji_r4
     !> Add multiplicity for J^T
     procedure(gs_interp_apply_fld), pass(this), deferred :: add_mult_jt_fld
     procedure(gs_interp_apply_r4), pass(this), deferred :: add_mult_jt_r4
     generic :: add_mult_jt => add_mult_jt_fld, add_mult_jt_r4
     !> Add multiplicity for J^-1
     procedure(gs_interp_apply_fld), pass(this), deferred :: add_mult_ji_fld
     procedure(gs_interp_apply_r4), pass(this), deferred :: add_mult_ji_r4
     generic :: add_mult_ji => add_mult_ji_fld
     !> AMR restart of a base type
     procedure, pass(this) :: amr_restart_base => gs_interp_amr_restart_base
  end type gs_interp_t

  abstract interface
     !> Initialise GS interpolation
     subroutine gs_interp_init(this, lx, conn)
       import gs_interp_t, mesh_conn_t
       class(gs_interp_t), intent(inout) :: this
       integer, intent(in) :: lx
       type(mesh_conn_t), target, intent(in) :: conn
     end subroutine gs_interp_init

     !> Initialise multiplicity arrays
     subroutine gs_interp_init_mult(this, mult_h1, mult_jt, mult_ji)
       import gs_interp_t, rp
       class(gs_interp_t), intent(inout) :: this
       real(rp), dimension(:, :, :, :) , intent(in) :: mult_h1, mult_jt, mult_ji
     end subroutine gs_interp_init_mult

     !> Free GS interpolation data
     subroutine gs_interp_free(this)
       import gs_interp_t
       class(gs_interp_t), intent(inout) :: this
     end subroutine gs_interp_free

     !> Children's nonconforming face/edge interpolation/zero
     subroutine gs_interp_apply_fld(this, field)
       import gs_interp_t, field_t
       class(gs_interp_t), intent(inout) :: this
       type(field_t), intent(inout) :: field
     end subroutine gs_interp_apply_fld

     subroutine gs_interp_apply_r4(this, vec)
       import gs_interp_t, rp
       class(gs_interp_t), intent(inout) :: this
       real(rp), dimension(:,  :, :, :), intent(inout) :: vec
     end subroutine gs_interp_apply_r4

     !> Children's nonconforming face/edge filling
     subroutine gs_interp_set_fld(this, field, cnst)
       import gs_interp_t, field_t, rp
       class(gs_interp_t), intent(inout) :: this
       type(field_t), intent(inout) :: field
       real(rp), intent(in) :: cnst
     end subroutine gs_interp_set_fld

     subroutine gs_interp_set_r4(this, vec, cnst)
       import gs_interp_t, rp
       class(gs_interp_t), intent(inout) :: this
       real(rp), dimension(:,  :, :, :), intent(inout) :: vec
       real(rp), intent(in) :: cnst
     end subroutine gs_interp_set_r4
  end interface

contains
  !> Initialise gs interpolation type
  !! @param[in]  lx    polynomial order + 1
  !! @param[in]  conn  mesh connectivity
  subroutine gs_interp_init_base(this, lx, conn)
    class(gs_interp_t), intent(inout) :: this
    integer, intent(in) :: lx
    type(mesh_conn_t), target, intent(in) :: conn

    call this%free_base()

    ! this type depends on the functional space size and connectivity
    ! information
    this%lx = lx
    this%nel = conn%nel
    this%conn => conn

    ! initialise interpolation arrays
    call this%interpolate%init(conn%tdim, lx)

    ! initialise hanging face/edge information
    call gs_interp_init_hang(this)

  end subroutine gs_interp_init_base

  !> Initialise face/edge hanging information
  !! @param[in]  this   GS interpolation type
  subroutine gs_interp_init_hang(this)
    class(gs_interp_t), intent(inout) :: this
    integer :: il, jl, itmp

    if (this%conn%ifhang_set .and. this%conn%ifhang) then
       this%ifhang = .true.
       this%nhang_el = 0
       do il = 1, this%conn%nel
          if (this%conn%hang(il)) this%nhang_el = this%nhang_el + 1
       end do
       if (this%nhang_el .eq. 0) call neko_error('gs_interp_init_hang: no &
            &hanging element on MPI rank')
       allocate(this%hang_el(this%nhang_el), &
            this%hang_edg_off(this%nhang_el + 1), &
            this%hang_fcs_off(this%nhang_el + 1))
       itmp = 0
       do il = 1, this%conn%nel
          if (this%conn%hang(il)) then
             itmp = itmp + 1
             this%hang_el(itmp) = il
          end if
       end do
       ! count hanging edges and mark element boundaries
       itmp = 1
       this%hang_edg_off(1) = itmp
       do il = 1, this%nhang_el
          do jl = 1, this%conn%edg%nobj
             ! face unrelated edges only
             if (this%conn%edg%hang(jl, this%hang_el(il)) .eq. 0 .or. &
                  this%conn%edg%hang(jl, this%hang_el(il)) .eq. 1) &
                  itmp = itmp + 1
          end do
          this%hang_edg_off(il + 1) = itmp
       end do
       this%nhang_edg = itmp - 1
       if (this%nhang_edg .gt. 0) then
          allocate(this%hang_edg(this%nhang_edg), &
               this%hang_edg_pos(this%nhang_edg))
          ! fill arrays
          itmp = 0
          do il = 1, this%nhang_el
             do jl = 1, this%conn%edg%nobj
                ! face unrelated edges only
                if (this%conn%edg%hang(jl, this%hang_el(il)) .eq. 0 .or. &
                     this%conn%edg%hang(jl, this%hang_el(il)) .eq. 1) then
                   itmp = itmp + 1
                   this%hang_edg(itmp) = jl
                   this%hang_edg_pos(itmp) = &
                        this%conn%edg%hang(jl, this%hang_el(il))
                end if
             end do
          end do
       end if

       ! count hanging faces and mark element boundaries
       itmp = 1
       this%hang_fcs_off(1) = itmp
       do il = 1, this%nhang_el
          do jl = 1, this%conn%fcs%nobj
             if (this%conn%fcs%hang(jl, this%hang_el(il)) .ne. -1) &
                  itmp = itmp + 1
          end do
          this%hang_fcs_off(il + 1) = itmp
       end do
       this%nhang_fcs = itmp - 1
       if (this%nhang_fcs .gt. 0) then
          allocate(this%hang_fcs(this%nhang_fcs), &
               this%hang_fcs_pos(this%nhang_fcs))
          ! fill arrays
          itmp = 0
          do il = 1, this%nhang_el
             do jl = 1, this%conn%fcs%nobj
                if (this%conn%fcs%hang(jl, this%hang_el(il)) .ne. -1) then
                   itmp = itmp + 1
                   this%hang_fcs(itmp) = jl
                   this%hang_fcs_pos(itmp) = &
                        this%conn%fcs%hang(jl, this%hang_el(il))
                end if
             end do
          end do
       end if
    end if

  end subroutine gs_interp_init_hang

  !> Free gs interpolation type
  subroutine gs_interp_free_base(this)
    class(gs_interp_t), intent(inout) :: this

    nullify(this%conn)

    this%lx = 0
    this%nel = 0

    call this%interpolate%free()

    call gs_interp_free_hang(this)

  end subroutine gs_interp_free_base

  !> Free gs interpolation type
  subroutine gs_interp_free_hang(this)
    class(gs_interp_t), intent(inout) :: this

    this%ifhang = .false.
    this%nhang_el = 0
    this%nhang_edg = 0
    this%nhang_fcs = 0
    if (allocated(this%hang_el)) deallocate(this%hang_el)
    if (allocated(this%hang_edg)) deallocate(this%hang_edg)
    if (allocated(this%hang_edg_pos)) deallocate(this%hang_edg_pos)
    if (allocated(this%hang_edg_off)) deallocate(this%hang_edg_off)
    if (allocated(this%hang_fcs)) deallocate(this%hang_fcs)
    if (allocated(this%hang_fcs_pos)) deallocate(this%hang_fcs_pos)
    if (allocated(this%hang_fcs_off)) deallocate(this%hang_fcs_off)

  end subroutine gs_interp_free_hang

  !> AMR restart
  !! @param[inout]  reconstruct   data reconstruction type
  !! @param[in]     counter       restart counter
  !! @param[in]     tstep         time step
  subroutine gs_interp_amr_restart_base(this)
    class(gs_interp_t), intent(inout) :: this

    ! update local element number
    this%nel = this%conn%nel
    ! reinitialise hanging face/edge information
    call gs_interp_free_hang(this)
    call gs_interp_init_hang(this)

  end subroutine gs_interp_amr_restart_base

  end module gs_interp
