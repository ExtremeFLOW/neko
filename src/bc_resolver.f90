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
  use gs_ops, only : GS_OP_MAX
  use neko_config, only : NEKO_BCKND_DEVICE
  use num_types, only : rp
  use utils, only : neko_error
  use device, only : device_get_ptr
  use device_math, only : device_cfill_mask
  use symmetry, only : symmetry_t
  use non_normal, only : non_normal_t
  use shear_stress, only : shear_stress_t
  use wall_model_bc, only : wall_model_bc_t
  use, intrinsic :: iso_c_binding, only : c_ptr, c_null_ptr
  implicit none
  private

  type, public, abstract :: bc_resolver_t
   contains
     procedure, pass(this) :: init => bc_resolver_init
     procedure, pass(this) :: free => bc_resolver_free
  end type bc_resolver_t

  type, public, extends(bc_resolver_t) :: scalar_bc_resolver_t
     type(mask_t) :: dof_mask
   contains
     procedure, pass(this) :: free => scalar_bc_resolver_free
     procedure, pass(this) :: mark_bc => scalar_bc_resolver_mark_bc
     procedure, pass(this) :: mark_bc_list => scalar_bc_resolver_mark_bc_list
     procedure, pass(this) :: apply => scalar_bc_resolver_apply
     generic :: mark => mark_bc, mark_bc_list
  end type scalar_bc_resolver_t

  type, public, extends(bc_resolver_t) :: vector_bc_resolver_t
     type(scalar_bc_resolver_t) :: x
     type(scalar_bc_resolver_t) :: y
     type(scalar_bc_resolver_t) :: z
   contains
     procedure, pass(this) :: free => vector_bc_resolver_free
     procedure, pass(this) :: mark_bc_x => vector_bc_resolver_mark_bc_x
     procedure, pass(this) :: mark_bc_list_x => vector_bc_resolver_mark_bc_list_x
     procedure, pass(this) :: mark_bc_y => vector_bc_resolver_mark_bc_y
     procedure, pass(this) :: mark_bc_list_y => vector_bc_resolver_mark_bc_list_y
     procedure, pass(this) :: mark_bc_z => vector_bc_resolver_mark_bc_z
     procedure, pass(this) :: mark_bc_list_z => vector_bc_resolver_mark_bc_list_z
     procedure, pass(this) :: apply => vector_bc_resolver_apply
     generic :: mark_x => mark_bc_x, mark_bc_list_x
     generic :: mark_y => mark_bc_y, mark_bc_list_y
     generic :: mark_z => mark_bc_z, mark_bc_list_z
  end type vector_bc_resolver_t

  type, public, extends(bc_resolver_t) :: coupled_vector_bc_resolver_t
     type(mask_t) :: dof_mask
     type(bc_list_t) :: bcs
     type(coef_t), pointer :: coef => null()
     type(dofmap_t), pointer :: dof => null()
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
     procedure, pass(this) :: finalize => coupled_vector_bc_resolver_finalize
   end type coupled_vector_bc_resolver_t

contains

  !> Initialize the boundary resolver.
  subroutine bc_resolver_init(this)
    class(bc_resolver_t), intent(inout) :: this
  end subroutine bc_resolver_init

  !> Free the boundary resolver.
  subroutine bc_resolver_free(this)
    class(bc_resolver_t), intent(inout) :: this
  end subroutine bc_resolver_free

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

  !> Free a vector boundary resolver.
  subroutine vector_bc_resolver_free(this)
    class(vector_bc_resolver_t), intent(inout) :: this

    call this%x%free()
    call this%y%free()
    call this%z%free()
  end subroutine vector_bc_resolver_free

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

  !> Add a physical boundary condition to the coupled resolver.
  subroutine coupled_vector_bc_resolver_mark_bc(this, bc)
    class(coupled_vector_bc_resolver_t), intent(inout) :: this
    class(bc_t), intent(inout), target :: bc

    if (.not. associated(this%coef)) then
       this%coef => bc%coef
       this%dof => bc%dof
       call this%bcs%init()
    end if

    call this%bcs%append(bc)
  end subroutine coupled_vector_bc_resolver_mark_bc

  !> Finalize the coupled resolver by resolving the accumulated BC list.
  subroutine coupled_vector_bc_resolver_finalize(this)
    class(coupled_vector_bc_resolver_t), intent(inout) :: this
    type(field_t), pointer :: boundary_flag
    type(field_t), pointer :: constraint_n_flag
    type(field_t), pointer :: constraint_t1_flag
    type(field_t), pointer :: constraint_t2_flag
    integer :: scratch_idx(4)
    integer, allocatable :: mask_values(:)
    integer :: i, j, n, m
    logical :: c_n, c_t1, c_t2
    class(bc_t), pointer :: bc

    call this%dof_mask%free()

    if (allocated(this%constraint_n)) deallocate(this%constraint_n)
    if (allocated(this%constraint_t1)) deallocate(this%constraint_t1)
    if (allocated(this%constraint_t2)) deallocate(this%constraint_t2)

    if (.not. associated(this%dof)) return
    if (this%bcs%size() .eq. 0) return

    call neko_scratch_registry%request_field(boundary_flag, scratch_idx(1), .true.)
    call neko_scratch_registry%request_field(constraint_n_flag, scratch_idx(2), .true.)
    call neko_scratch_registry%request_field(constraint_t1_flag, scratch_idx(3), .true.)
    call neko_scratch_registry%request_field(constraint_t2_flag, scratch_idx(4), .true.)

    n = this%dof%size()

    do i = 1, this%bcs%size()
       bc => this%bcs%get(i)

       if (.not. allocated(bc%msk)) then
          call neko_error("Attempting to finalize coupled resolver from an unfinalized BC.")
       end if

       c_n = .false.
       c_t1 = .false.
       c_t2 = .false.

       select type (bc)
       type is (symmetry_t)
          c_n = .true.
       type is (shear_stress_t)
          c_n = .true.
       type is (wall_model_bc_t)
          c_n = .true.
       type is (non_normal_t)
          c_t1 = .true.
          c_t2 = .true.
       class default
          if (bc%strong) then
             c_n = .true.
             c_t1 = .true.
             c_t2 = .true.
          end if
       end select

       do j = 1, bc%msk(0)
          m = bc%msk(j)
          boundary_flag%x(m,1,1,1) = 1.0_rp
          if (c_n) constraint_n_flag%x(m,1,1,1) = 1.0_rp
          if (c_t1) constraint_t1_flag%x(m,1,1,1) = 1.0_rp
          if (c_t2) constraint_t2_flag%x(m,1,1,1) = 1.0_rp
       end do
    end do

    call this%coef%gs_h%op(boundary_flag, GS_OP_MAX)
    call this%coef%gs_h%op(constraint_n_flag, GS_OP_MAX)
    call this%coef%gs_h%op(constraint_t1_flag, GS_OP_MAX)
    call this%coef%gs_h%op(constraint_t2_flag, GS_OP_MAX)

    m = 0
    do i = 1, n
      if (boundary_flag%x(i,1,1,1) .gt. 0.5_rp) then
         m = m + 1
      end if
    end do

    if (m .gt. 0) then
       allocate(mask_values(m))
       j = 0
       do i = 1, n
          if (boundary_flag%x(i,1,1,1) .gt. 0.5_rp) then
             j = j + 1
             mask_values(j) = i
          end if
       end do

       call this%dof_mask%init(mask_values(1:m), m)
       allocate(this%constraint_n(m))
       allocate(this%constraint_t1(m))
       allocate(this%constraint_t2(m))

       do i = 1, m
          j = mask_values(i)
          this%constraint_n(i) = constraint_n_flag%x(j,1,1,1) .gt. 0.5_rp
          this%constraint_t1(i) = constraint_t1_flag%x(j,1,1,1) .gt. 0.5_rp
          this%constraint_t2(i) = constraint_t2_flag%x(j,1,1,1) .gt. 0.5_rp
       end do
    end if

    call neko_scratch_registry%relinquish_field(scratch_idx)

    if (allocated(mask_values)) then
       deallocate(mask_values)
    end if
  end subroutine coupled_vector_bc_resolver_finalize

  !> Add the constrained dofs from an x-component boundary condition.
  subroutine vector_bc_resolver_mark_bc_x(this, bc)
    class(vector_bc_resolver_t), intent(inout) :: this
    class(bc_t), intent(in) :: bc

    call this%x%mark_bc(bc)
  end subroutine vector_bc_resolver_mark_bc_x

  !> Add the constrained dofs from all x-component boundary conditions in a list.
  subroutine vector_bc_resolver_mark_bc_list_x(this, bclst)
    class(vector_bc_resolver_t), intent(inout) :: this
    type(bc_list_t), intent(in) :: bclst

    call this%x%mark_bc_list(bclst)
  end subroutine vector_bc_resolver_mark_bc_list_x

  !> Add the constrained dofs from a y-component boundary condition.
  subroutine vector_bc_resolver_mark_bc_y(this, bc)
    class(vector_bc_resolver_t), intent(inout) :: this
    class(bc_t), intent(in) :: bc

    call this%y%mark_bc(bc)
  end subroutine vector_bc_resolver_mark_bc_y

  !> Add the constrained dofs from all y-component boundary conditions in a list.
  subroutine vector_bc_resolver_mark_bc_list_y(this, bclst)
    class(vector_bc_resolver_t), intent(inout) :: this
    type(bc_list_t), intent(in) :: bclst

    call this%y%mark_bc_list(bclst)
  end subroutine vector_bc_resolver_mark_bc_list_y

  !> Add the constrained dofs from a z-component boundary condition.
  subroutine vector_bc_resolver_mark_bc_z(this, bc)
    class(vector_bc_resolver_t), intent(inout) :: this
    class(bc_t), intent(in) :: bc

    call this%z%mark_bc(bc)
  end subroutine vector_bc_resolver_mark_bc_z

  !> Add the constrained dofs from all z-component boundary conditions in a list.
  subroutine vector_bc_resolver_mark_bc_list_z(this, bclst)
    class(vector_bc_resolver_t), intent(inout) :: this
    type(bc_list_t), intent(in) :: bclst

    call this%z%mark_bc_list(bclst)
  end subroutine vector_bc_resolver_mark_bc_list_z

  !> Apply the vector boundary constraints component-wise.
  subroutine vector_bc_resolver_apply(this, x, y, z, n)
    class(vector_bc_resolver_t), intent(in) :: this
    integer, intent(in) :: n
    real(kind=rp), intent(inout) :: x(n)
    real(kind=rp), intent(inout) :: y(n)
    real(kind=rp), intent(inout) :: z(n)

    call this%x%apply(x, n)
    call this%y%apply(y, n)
    call this%z%apply(z, n)
  end subroutine vector_bc_resolver_apply

end module bc_resolver
