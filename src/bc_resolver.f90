! Copyright (c) 2026, The Neko Authors
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
!> Scaffolding for future global boundary resolution.
module bc_resolver
  use bc, only : bc_t
  use bc_list, only : bc_list_t
  use mask, only : mask_t
  use coefs, only : coef_t
  use dofmap, only : dofmap_t
  use field, only : field_t
  use scratch_registry, only : neko_scratch_registry
  use math, only : cfill_mask
  use gs_ops, only : GS_OP_MAX, GS_OP_ADD, GS_OP_MIN
  use neko_config, only : NEKO_BCKND_DEVICE
  use num_types, only : rp
  use utils, only : neko_error, nonlinear_index
  use device, only : device_get_ptr, device_memcpy, HOST_TO_DEVICE, &
       DEVICE_TO_HOST
  use device_math, only : device_cfill_mask
  use field_math, only : field_cfill
  use, intrinsic :: iso_c_binding, only : c_ptr, c_null_ptr
  implicit none
  private
  public :: vector_bc_resolver_components

  type, public :: scalar_bc_resolver_t
     type(mask_t) :: dof_mask
   contains
     procedure, pass(this) :: free => scalar_bc_resolver_free
     procedure, pass(this) :: mark_bc => scalar_bc_resolver_mark_bc
     procedure, pass(this) :: mark_bc_list => scalar_bc_resolver_mark_bc_list
     procedure, pass(this) :: apply => scalar_bc_resolver_apply
     generic :: mark => mark_bc, mark_bc_list
  end type scalar_bc_resolver_t

  type, public, abstract :: vector_bc_resolver_t
   contains
     procedure(vector_bc_resolver_free_intrf), pass(this), deferred :: free
     procedure(vector_bc_resolver_finalize_intrf), pass(this), deferred :: &
          finalize
     procedure(vector_bc_resolver_apply_intrf), pass(this), deferred :: apply
     procedure(vector_bc_resolver_mark_bc_intrf), pass(this), deferred :: mark_bc
     procedure(vector_bc_resolver_mark_bc_list_intrf), pass(this), deferred :: &
          mark_bc_list
     generic :: mark => mark_bc, mark_bc_list
  end type vector_bc_resolver_t

  type, public, extends(vector_bc_resolver_t) :: segregated_vector_bc_resolver_t
     type(scalar_bc_resolver_t) :: x
     type(scalar_bc_resolver_t) :: y
     type(scalar_bc_resolver_t) :: z
   contains
     procedure, pass(this) :: free => segregated_vector_bc_resolver_free
     procedure, pass(this) :: finalize => segregated_vector_bc_resolver_finalize
     procedure, pass(this) :: mark_bc => segregated_vector_bc_resolver_mark_bc
     procedure, pass(this) :: mark_bc_list => &
          segregated_vector_bc_resolver_mark_bc_list
     procedure, pass(this) :: apply => segregated_vector_bc_resolver_apply
  end type segregated_vector_bc_resolver_t

  type, public, extends(vector_bc_resolver_t) :: coupled_vector_bc_resolver_t
     type(mask_t) :: dof_mask
     type(bc_list_t), private :: bcs
     type(coef_t), pointer, private :: coef => null()
     type(dofmap_t), pointer, private :: dof => null()
     logical, allocatable :: constraint_n(:)
     logical, allocatable :: constraint_t1(:)
     logical, allocatable :: constraint_t2(:)
     real(kind=rp), allocatable :: n(:,:)
     real(kind=rp), allocatable :: t1(:,:)
     real(kind=rp), allocatable :: t2(:,:)
     type(c_ptr) :: constraint_n_d = c_null_ptr
     type(c_ptr) :: constraint_t1_d = c_null_ptr
     type(c_ptr) :: constraint_t2_d = c_null_ptr
     type(c_ptr) :: n_d = c_null_ptr
     type(c_ptr) :: t1_d = c_null_ptr
     type(c_ptr) :: t2_d = c_null_ptr
   contains
     procedure, pass(this) :: free => coupled_vector_bc_resolver_free
     procedure, pass(this) :: mark_bc => coupled_vector_bc_resolver_mark_bc
     procedure, pass(this) :: mark_bc_list => &
          coupled_vector_bc_resolver_mark_bc_list
     procedure, pass(this) :: finalize => coupled_vector_bc_resolver_finalize
     procedure, pass(this) :: apply => coupled_vector_bc_resolver_apply
  end type coupled_vector_bc_resolver_t

  abstract interface
     subroutine vector_bc_resolver_free_intrf(this)
       import :: vector_bc_resolver_t
       class(vector_bc_resolver_t), intent(inout) :: this
     end subroutine vector_bc_resolver_free_intrf
  end interface

  abstract interface
     subroutine vector_bc_resolver_finalize_intrf(this)
       import :: vector_bc_resolver_t
       class(vector_bc_resolver_t), intent(inout) :: this
     end subroutine vector_bc_resolver_finalize_intrf
  end interface

  abstract interface
     subroutine vector_bc_resolver_apply_intrf(this, x, y, z, n)
       import :: vector_bc_resolver_t, rp
       class(vector_bc_resolver_t), intent(in) :: this
       integer, intent(in) :: n
       real(kind=rp), intent(inout) :: x(n)
       real(kind=rp), intent(inout) :: y(n)
       real(kind=rp), intent(inout) :: z(n)
     end subroutine vector_bc_resolver_apply_intrf
  end interface

  abstract interface
     subroutine vector_bc_resolver_mark_bc_intrf(this, bc, component)
       import :: vector_bc_resolver_t, bc_t
       class(vector_bc_resolver_t), intent(inout) :: this
       class(bc_t), intent(inout), target :: bc
       character(len=1), optional, intent(in) :: component
     end subroutine vector_bc_resolver_mark_bc_intrf
  end interface

  abstract interface
     subroutine vector_bc_resolver_mark_bc_list_intrf(this, bclst, component)
       import :: vector_bc_resolver_t, bc_list_t
       class(vector_bc_resolver_t), intent(inout) :: this
       type(bc_list_t), intent(in) :: bclst
       character(len=1), optional, intent(in) :: component
     end subroutine vector_bc_resolver_mark_bc_list_intrf
  end interface

contains

!
!  ************** scalar_bc_resolver_t TBPs **************
!

  !> Free a scalar boundary resolver.
  subroutine scalar_bc_resolver_free(this)
    class(scalar_bc_resolver_t), intent(inout) :: this

    call this%dof_mask%free()
  end subroutine scalar_bc_resolver_free

  !> Add the constrained dofs from a scalar boundary condition.
  subroutine scalar_bc_resolver_mark_bc(this, bc)
    class(scalar_bc_resolver_t), intent(inout) :: this
    class(bc_t), intent(in) :: bc

    integer, pointer :: current_mask(:)
    integer, allocatable :: merged_mask(:)
    integer :: current_size
    integer :: incoming_size
    integer :: merged_size
    integer :: i

    if (.not. allocated(bc%msk)) then
       call neko_error("Attempting to mark resolver from an unfinalized BC.")
    end if

    incoming_size = bc%msk(0)
    if (incoming_size .eq. 0) return

    if (.not. this%dof_mask%is_set()) then
       call this%dof_mask%init(bc%msk(1:incoming_size), incoming_size)
       return
    end if

    current_size = this%dof_mask%size()
    current_mask => this%dof_mask%get()
    allocate(merged_mask(current_size + incoming_size))

    merged_mask(1:current_size) = current_mask(1:current_size)
    merged_size = current_size

    do i = 1, incoming_size
       if (.not. any(merged_mask(1:merged_size) .eq. bc%msk(i))) then
          merged_size = merged_size + 1
          merged_mask(merged_size) = bc%msk(i)
       end if
    end do

    call this%dof_mask%set(merged_mask(1:merged_size), merged_size)
  end subroutine scalar_bc_resolver_mark_bc

  !> Add the constrained dofs from all boundary conditions in a list.
  subroutine scalar_bc_resolver_mark_bc_list(this, bclst)
    class(scalar_bc_resolver_t), intent(inout) :: this
    type(bc_list_t), intent(in) :: bclst

    integer :: i

    do i = 1, bclst%size()
       call this%mark_bc(bclst%get(i))
    end do
  end subroutine scalar_bc_resolver_mark_bc_list

  !> Apply the scalar boundary constraints by zeroing constrained dofs.
  subroutine scalar_bc_resolver_apply(this, x, n)
    class(scalar_bc_resolver_t), intent(in) :: this
    integer, intent(in) :: n
    real(kind=rp), intent(inout) :: x(n)
    type(c_ptr) :: x_d

    if (.not. this%dof_mask%is_set()) return

    if (NEKO_BCKND_DEVICE .eq. 1) then
       x_d = device_get_ptr(x)
       call device_cfill_mask(x_d, 0.0_rp, n, this%dof_mask%get_d(), &
            this%dof_mask%size())
    else
       call cfill_mask(x, 0.0_rp, n, this%dof_mask%get(), this%dof_mask%size())
    end if
  end subroutine scalar_bc_resolver_apply

!
!  ********** segregated_vector_bc_resolver_t TBPs **********
!

  !> Free a segregated vector boundary resolver.
  subroutine segregated_vector_bc_resolver_free(this)
    class(segregated_vector_bc_resolver_t), intent(inout) :: this

    call this%x%free()
    call this%y%free()
    call this%z%free()
  end subroutine segregated_vector_bc_resolver_free

  !> Finalize a segregated vector boundary resolver.
  subroutine segregated_vector_bc_resolver_finalize(this)
    class(segregated_vector_bc_resolver_t), intent(inout) :: this
  end subroutine segregated_vector_bc_resolver_finalize

  !> Add the constrained dofs from a vector boundary condition.
  subroutine segregated_vector_bc_resolver_mark_bc(this, bc, component)
    class(segregated_vector_bc_resolver_t), intent(inout) :: this
    class(bc_t), intent(inout), target :: bc
    character(len=1), optional, intent(in) :: component

    if (.not. all(bc%constraints)) then
       call neko_error("Segregated vector BC resolver only accepts " // &
            "BCs with constraints = (.true., .true., .true.).")
    end if

    if (.not. present(component)) then
       call neko_error("Segregated vector BC resolver mark requires " // &
            "component='x', 'y', or 'z'.")
    end if

    select case (component)
    case ('x')
       call this%x%mark_bc(bc)
    case ('y')
       call this%y%mark_bc(bc)
    case ('z')
       call this%z%mark_bc(bc)
    case default
       call neko_error("Invalid component for segregated vector BC " // &
            "resolver mark.")
    end select
  end subroutine segregated_vector_bc_resolver_mark_bc

  !> Add the constrained dofs from all vector boundary conditions in a list.
  subroutine segregated_vector_bc_resolver_mark_bc_list(this, bclst, component)
    class(segregated_vector_bc_resolver_t), intent(inout) :: this
    type(bc_list_t), intent(in) :: bclst
    character(len=1), optional, intent(in) :: component
    integer :: i

    do i = 1, bclst%size()
       call this%mark_bc(bclst%get(i), component)
    end do
  end subroutine segregated_vector_bc_resolver_mark_bc_list

  !> Apply the vector boundary constraints component-wise.
  subroutine segregated_vector_bc_resolver_apply(this, x, y, z, n)
    class(segregated_vector_bc_resolver_t), intent(in) :: this
    integer, intent(in) :: n
    real(kind=rp), intent(inout) :: x(n)
    real(kind=rp), intent(inout) :: y(n)
    real(kind=rp), intent(inout) :: z(n)

    call this%x%apply(x, n)
    call this%y%apply(y, n)
    call this%z%apply(z, n)
  end subroutine segregated_vector_bc_resolver_apply

!
!  ************** vector_bc_resolver_t helpers **************
!

  !> Return scalar component resolvers for a segregated vector resolver.
  subroutine vector_bc_resolver_components(this, x, y, z)
    class(vector_bc_resolver_t), target, intent(inout) :: this
    type(scalar_bc_resolver_t), pointer :: x
    type(scalar_bc_resolver_t), pointer :: y
    type(scalar_bc_resolver_t), pointer :: z

    select type (this)
    type is (segregated_vector_bc_resolver_t)
       x => this%x
       y => this%y
       z => this%z
    class default
       call neko_error("Component access is only available for " // &
            "segregated vector BC resolvers.")
    end select
  end subroutine vector_bc_resolver_components

!
!  *********** coupled_vector_bc_resolver_t TBPs ***********
!

  !> Free a coupled vector boundary resolver.
  subroutine coupled_vector_bc_resolver_free(this)
    class(coupled_vector_bc_resolver_t), intent(inout) :: this

    call this%dof_mask%free()
    call this%bcs%free()

    if (allocated(this%constraint_n)) then
       deallocate(this%constraint_n)
    end if

    if (allocated(this%constraint_t1)) then
       deallocate(this%constraint_t1)
    end if

    if (allocated(this%constraint_t2)) then
       deallocate(this%constraint_t2)
    end if

    if (allocated(this%n)) then
       deallocate(this%n)
    end if

    if (allocated(this%t1)) then
       deallocate(this%t1)
    end if

    if (allocated(this%t2)) then
       deallocate(this%t2)
    end if

    this%constraint_n_d = c_null_ptr
    this%constraint_t1_d = c_null_ptr
    this%constraint_t2_d = c_null_ptr
    this%n_d = c_null_ptr
    this%t1_d = c_null_ptr
    this%t2_d = c_null_ptr
    nullify(this%coef)
    nullify(this%dof)
  end subroutine coupled_vector_bc_resolver_free

  !> Add a boundary condition to the coupled resolver.
  subroutine coupled_vector_bc_resolver_mark_bc(this, bc, component)
    class(coupled_vector_bc_resolver_t), intent(inout) :: this
    class(bc_t), intent(inout), target :: bc
    character(len=1), optional, intent(in) :: component

    if (.not. associated(this%coef)) then
       this%coef => bc%coef
       this%dof => bc%dof
       call this%bcs%init()
    end if

    call this%bcs%append(bc)
  end subroutine coupled_vector_bc_resolver_mark_bc

  !> Add all boundary conditions in a list to the coupled resolver.
  subroutine coupled_vector_bc_resolver_mark_bc_list(this, bclst, component)
    class(coupled_vector_bc_resolver_t), intent(inout) :: this
    type(bc_list_t), intent(in) :: bclst
    character(len=1), optional, intent(in) :: component
    integer :: i

    do i = 1, bclst%size()
       call this%mark_bc(bclst%get(i), component)
    end do
  end subroutine coupled_vector_bc_resolver_mark_bc_list

  !> Finalize the coupled resolver by resolving the accumulated BC list.
  subroutine coupled_vector_bc_resolver_finalize(this)
    class(coupled_vector_bc_resolver_t), intent(inout) :: this
    type(field_t), pointer :: work1
    type(field_t), pointer :: work2
    type(field_t), pointer :: work3
    type(field_t), pointer :: work4
    integer :: scratch_idx(4)
    integer, allocatable :: mask_values(:)
    integer :: i, j, n, m
    integer :: idx(4), facet
    logical :: c_n, c_t1, c_t2
    real(kind=rp) :: normal(3), ref(3), t1_vec(3), t2_vec(3), len, prio
    class(bc_t), pointer :: bc

    call this%dof_mask%free()

    if (allocated(this%constraint_n)) deallocate(this%constraint_n)
    if (allocated(this%constraint_t1)) deallocate(this%constraint_t1)
    if (allocated(this%constraint_t2)) deallocate(this%constraint_t2)

    if (.not. associated(this%dof)) return
    if (this%bcs%size() .eq. 0) return

    call neko_scratch_registry%request_field(work1, scratch_idx(1), .true.)
    call neko_scratch_registry%request_field(work2, scratch_idx(2), .true.)
    call neko_scratch_registry%request_field(work3, scratch_idx(3), .true.)
    call neko_scratch_registry%request_field(work4, scratch_idx(4), .true.)

    n = this%dof%size()

    do i = 1, this%bcs%size()
       bc => this%bcs%get(i)

       if (.not. allocated(bc%msk)) then
          call neko_error("Attempting to finalize coupled resolver from an " // &
               "unfinalized BC.")
       end if

       ! Mask all the dofs touched by this BC. Since %msk is propagated to all
       ! local dofs via gather-scatter, work1 will contain all local nodes on
       ! the boundary, including those elements that don't touch it with a face.
       call cfill_mask(work1%x(:,1,1,1), 1.0_rp, n, &
            bc%msk(1:bc%msk(0)), bc%msk(0))
    end do


    ! Set priority values for constraint assignment. Mimics the procedure in
    ! Nek5000 directly.
    ! The values are chosen in a way that a min reduction applies the most
    ! restrictive constraint.
    ! 5 -> unconstrained
    ! 3 -> tangentially constrained
    ! 2 -> normally constrained
    ! 0 -> fully constrained
    ! Fill the field to not mess up gather-scatter reduction later.
    call field_cfill(work2, 5.0_rp)
    do i = 1, this%bcs%size()
       bc => this%bcs%get(i)


       prio = 5.0_rp

       if (all(bc%constraints)) then
          prio = 0.0_rp
       else if ((bc%constraints(1)) .and. (.not. bc%constraints(2)) &
            .and. (.not. bc%constraints(3))) then
          prio = 2.0_rp
       else if ((.not. bc%constraints(1)) .and. bc%constraints(2) &
            .and. bc%constraints(3)) then
          prio = 3.0_rp
       else
          call neko_error("Unsupported constraint combination in " // &
               "vector BC resolver.")
       end if

       ! Note that the face_msk is used, so constraints are only directly
       ! applied to elements that touch the boundary at this point.
       ! The min here ensures that the most restricive constraint is kept 
       ! within a single element.
       do j = 1, bc%facet_msk(0)
          m = bc%facet_msk(j)
          work2%x(m,1,1,1) = min(prio, work2%x(m,1,1,1))
       end do
    end do

    ! Propagate constraints to all local dofs via gather-scatter.
    ! Ensures most restrictive constraint is kept across element boundaries.
    call this%coef%gs_h%op(work2, GS_OP_MIN)

    ! Compute the number of constrained dofs so we can allocate vectors.
    m = 0
    do i = 1, n
       if (work1%x(i,1,1,1) .gt. 0.5_rp) then
          m = m + 1
       end if
    end do

    if (m .eq. 0) then
       call neko_scratch_registry%relinquish_field(scratch_idx)
       return
    end if

    ! Fill in the full mask indices and copy to this%dof_mask.
    allocate(mask_values(m))
    j = 0
    do i = 1, n
       if (work1%x(i,1,1,1) .gt. 0.5_rp) then
          j = j + 1
          mask_values(j) = i
       end if
    end do
    call this%dof_mask%init(mask_values(1:m), m)

    ! Allocate constraints
    allocate(this%constraint_n(m))
    allocate(this%constraint_t1(m))
    allocate(this%constraint_t2(m))

    ! Use the priority value to assign constraints.
    do i = 1, m
       j = mask_values(i)

       if (work2%x(j,1,1,1) .lt. 1.9_rp) then
          this%constraint_n(i) = .true.
          this%constraint_t1(i) = .true.
          this%constraint_t2(i) = .true.
       else if (work2%x(j,1,1,1) .gt. 1.9_rp .and. &
            work2%x(j,1,1,1) .lt. 2.9_rp) then
          this%constraint_n(i) = .true.
          this%constraint_t1(i) = .false.
          this%constraint_t2(i) = .false.
       else if (work2%x(j,1,1,1) .gt. 2.9_rp .and. &
            work2%x(j,1,1,1) .lt. 3.9_rp) then
          this%constraint_n(i) = .false.
          this%constraint_t1(i) = .true.
          this%constraint_t2(i) = .true.
       else if (work2%x(j,1,1,1) .gt. 3.9_rp) then
          this%constraint_n(i) = .false.
          this%constraint_t1(i) = .false.
          this%constraint_t2(i) = .false.
       end if
    end do

    allocate(this%n(3,m))
    allocate(this%t1(3,m))
    allocate(this%t2(3,m))

    associate(nx => work1, ny => work2, nz => work3)
      nx%x = 0.0_rp
      ny%x = 0.0_rp
      nz%x = 0.0_rp

      do i = 1, this%bcs%size()
         bc => this%bcs%get(i)
         do j = 1, bc%facet_msk(0)
            m = bc%facet_msk(j)
            facet = bc%facet(j)
            idx = nonlinear_index(m, this%coef%Xh%lx, this%coef%Xh%ly, &
                 this%coef%Xh%lz)
            normal = this%coef%get_normal(idx(1), idx(2), idx(3), idx(4), facet)
            nx%x(m,1,1,1) = nx%x(m,1,1,1) + normal(1)
            ny%x(m,1,1,1) = ny%x(m,1,1,1) + normal(2)
            nz%x(m,1,1,1) = nz%x(m,1,1,1) + normal(3)
         end do
      end do

      this%n = 0.0_rp
      this%t1 = 0.0_rp
      this%t2 = 0.0_rp

      do i = 1, size(mask_values)
         j = mask_values(i)
         normal(1) = nx%x(j,1,1,1)
         normal(2) = ny%x(j,1,1,1)
         normal(3) = nz%x(j,1,1,1)
         len = sqrt(sum(normal**2))
         if (len .gt. 0.0_rp) then
            this%n(:,i) = normal / len

            if (this%coef%msh%gdim .eq. 2) then
               t1_vec = [ -this%n(2,i), this%n(1,i), 0.0_rp ]
               len = sqrt(sum(t1_vec**2))
               if (len .gt. 0.0_rp) then
                  this%t1(:,i) = t1_vec / len
               end if
               this%t2(:,i) = [ 0.0_rp, 0.0_rp, 1.0_rp ]
            else
               if (abs(this%n(1,i)) .lt. 0.9_rp) then
                  ref = [ 1.0_rp, 0.0_rp, 0.0_rp ]
               else
                  ref = [ 0.0_rp, 1.0_rp, 0.0_rp ]
               end if

               t1_vec = ref - sum(ref * this%n(:,i)) * this%n(:,i)
               len = sqrt(sum(t1_vec**2))
               if (len .gt. 0.0_rp) then
                  this%t1(:,i) = t1_vec / len
               end if

               t2_vec(1) = this%n(2,i) * this%t1(3,i) - &
                    this%n(3,i) * this%t1(2,i)
               t2_vec(2) = this%n(3,i) * this%t1(1,i) - &
                    this%n(1,i) * this%t1(3,i)
               t2_vec(3) = this%n(1,i) * this%t1(2,i) - &
                    this%n(2,i) * this%t1(1,i)
               len = sqrt(sum(t2_vec**2))
               if (len .gt. 0.0_rp) then
                  this%t2(:,i) = t2_vec / len
               end if
            end if
         end if
      end do
    end associate

    call neko_scratch_registry%relinquish_field(scratch_idx)

    if (allocated(mask_values)) then
       deallocate(mask_values)
    end if
  end subroutine coupled_vector_bc_resolver_finalize

  !> Apply the coupled vector boundary constraints in the local basis.
  subroutine coupled_vector_bc_resolver_apply(this, x, y, z, n)
    class(coupled_vector_bc_resolver_t), intent(in) :: this
    integer, intent(in) :: n
    real(kind=rp), intent(inout) :: x(n)
    real(kind=rp), intent(inout) :: y(n)
    real(kind=rp), intent(inout) :: z(n)
    integer, pointer :: msk(:)
    integer :: i, j, m
    real(kind=rp) :: u(3), uloc(3)

    if (.not. this%dof_mask%is_set()) return

    m = this%dof_mask%size()
    msk => this%dof_mask%get()

    do i = 1, m
       j = msk(i)

       u(1) = x(j)
       u(2) = y(j)
       u(3) = z(j)

       uloc(1) = u(1) * this%n(1,i) + u(2) * this%n(2,i) + &
            u(3) * this%n(3,i)
       uloc(2) = u(1) * this%t1(1,i) + u(2) * this%t1(2,i) + &
            u(3) * this%t1(3,i)
       uloc(3) = u(1) * this%t2(1,i) + u(2) * this%t2(2,i) + &
            u(3) * this%t2(3,i)

       if (this%constraint_n(i)) uloc(1) = 0.0_rp
       if (this%constraint_t1(i)) uloc(2) = 0.0_rp
       if (this%constraint_t2(i)) uloc(3) = 0.0_rp

       u = uloc(1) * this%n(:,i) + uloc(2) * this%t1(:,i) + &
            uloc(3) * this%t2(:,i)

       x(j) = u(1)
       y(j) = u(2)
       z(j) = u(3)
    end do
  end subroutine coupled_vector_bc_resolver_apply

end module bc_resolver
