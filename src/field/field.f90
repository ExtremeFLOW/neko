! Copyright (c) 2018-2025, The Neko Authors
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
!> Defines a field
module field
  use neko_config, only : NEKO_BCKND_DEVICE
  use device_math, only : device_add2, device_cadd, device_copy
  use num_types, only : rp
  use math, only : add2, copy, cadd
  use mesh, only : mesh_t
  use space, only : space_t, operator(.ne.)
  use dofmap, only : dofmap_t
  use tensor4, only : tensor4_t
  use utils, only : neko_error, NEKO_VARNAME_LEN
  use, intrinsic :: iso_c_binding
  implicit none
  private

  type, public, extends(tensor4_t) :: field_t
     type(space_t), pointer :: Xh !< Function space \f$ X_h \f$
     type(mesh_t), pointer :: msh !< Mesh
     type(dofmap_t), pointer :: dof !< Dofmap

     logical :: internal_dofmap = .false. !< Does the field have an own dofmap
   contains
     procedure, private, pass(this) :: init_common => field_init_common
     procedure, private, pass(this) :: init_external_dof => &
          field_init_external_dof
     procedure, private, pass(this) :: init_internal_dof => &
          field_init_internal_dof
     !> Seals off tensor4_t's raw dimension initialiser: field_t must always
     !! be initialised with a mesh/space or dofmap, never bare dimensions,
     !! since that would leave Xh/msh/dof unassociated. Overriding
     !! `init_dims` here (matching tensor4_t's exact signature and passed-
     !! object argument name) replaces it within the inherited `init`
     !! generic rather than merging alongside it as an ambiguous overload.
     procedure, private, pass(t) :: init_dims => field_init_dims_seal
     procedure, private, pass(this) :: assign_field => field_assign_field
     procedure, private, pass(this) :: add_field => field_add_field
     procedure, private, pass(this) :: add_scalar => field_add_scalar
     procedure, pass(t) :: free => field_free
     !> Seals off tensor4_t's raw tensor assignment so `field = tensor4`
     !! cannot silently bypass Xh/msh/dof; same reasoning as init_dims
     !! above. Not ambiguous with assign_field below since the right-hand
     !! side types differ (tensor4_t vs. field_t).
     procedure, pass(t) :: tensor4_assign_tensor4 => &
          field_assign_tensor4_seal
     !> Scalar assignment is *not* overridden: tensor4_t's own
     !! tensor4_assign_scalar (cfill by total size) is exactly what a
     !! field-level scalar fill needs too, with no Xh/msh/dof involved, so
     !! it is left to merge into `assignment(=)` below unchanged.
     !> Initialise a field
     generic :: init => init_external_dof, init_internal_dof
     !> Assignemnt to current field
     generic :: assignment(=) => assign_field
     !> Add to current field
     !! @note We don't overload operator(+), to avoid
     !! the extra assignemnt operator
     generic :: add => add_field, add_scalar
  end type field_t

  !> field_ptr_t, To easily obtain a pointer to a field
  type, public :: field_ptr_t
     type(field_t), pointer :: ptr => null()
   contains
     !> Constructor. Just assigns the pointer.
     procedure, pass(this) :: init => field_ptr_init
     !> Destructor. Just nullifies the pointer.
     procedure, pass(this) :: free => field_ptr_free
  end type field_ptr_t

  !> field_wrapper_t, used to wrap an allocated field for use in a field list
  type, public :: field_wrapper_t
     type(field_t), pointer :: field => null()
   contains
     !> Constructor. Allocates a field and assigns the pointer.
     generic :: init => init_field, init_internal_dof, init_external_dof
     !> Initialize a field wrapper with an allocated field
     procedure, pass(this) :: init_field => field_wrapper_init_field
     !> Initialize a field wrapper with an internal dofmap
     procedure, pass(this) :: init_internal_dof => &
          field_wrapper_init_internal_dof
     !> Initialize a field wrapper with an external dofmap
     procedure, pass(this) :: init_external_dof => &
          field_wrapper_init_external_dof
     !> Destructor. Frees the field and nullifies the pointer.
     procedure, pass(this) :: free => field_wrapper_free
  end type field_wrapper_t

contains

  !> Initialize a field @a this on the mesh @a msh using an internal dofmap
  subroutine field_init_internal_dof(this, msh, space, fld_name)
    class(field_t), intent(inout) :: this !< Field to be initialized
    type(mesh_t), target, intent(in) :: msh !< underlying mesh of the field
    type(space_t), target, intent(in) :: space !< Function space for the field
    character(len=*), optional :: fld_name !< Name of the field

    call this%free()

    this%Xh => space
    this%msh => msh

    allocate(this%dof)
    call this%dof%init(this%msh, this%Xh)
    this%internal_dofmap = .true.

    if (present(fld_name)) then
       call this%init_common(fld_name)
    else
       call this%init_common()
    end if

  end subroutine field_init_internal_dof

  !> Initialize a field @a this on the mesh @a msh using an internal dofmap
  subroutine field_init_external_dof(this, dof, fld_name)
    class(field_t), intent(inout) :: this !< Field to be initialized
    type(dofmap_t), target, intent(in) :: dof !< External dofmap for the field
    character(len=*), optional :: fld_name !< Name of the field

    call this%free()

    this%dof => dof
    this%Xh => dof%Xh
    this%msh => dof%msh

    if (present(fld_name)) then
       call this%init_common(fld_name)
    else
       call this%init_common()
    end if

  end subroutine field_init_external_dof

  !> Initialize a field @a this
  !! @note Called with @a this%x already deallocated (init_external_dof and
  !! init_internal_dof both call `this%free()` before reaching here), so the
  !! parent init's alloc-always-frees-first semantics are a no-op on the
  !! array and safe to rely on unconditionally.
  subroutine field_init_common(this, fld_name)
    class(field_t), intent(inout) :: this !< Field to be initialized
    character(len=*), optional :: fld_name !< Name of the field

    associate(lx => this%Xh%lx, ly => this%Xh%ly, &
         lz => this%Xh%lz, nelv => this%msh%nelv)

      ! tensor4_t%init (via tensor4_allocate) zeroes the device side
      ! (and synchronizes) before the host-side zero-fill below: under
      ! zero-copy the device then faults the pages first (device first
      ! touch), which gives contiguous physical mappings and thus
      ! better GPU TLB utilisation; rewriting the zeros on the host
      ! afterwards is benign. See tensor4_init/tensor4_allocate.
      if (present(fld_name)) then
         call this%tensor4_t%init(lx, ly, lz, nelv, fld_name)
      else
         call this%tensor4_t%init(lx, ly, lz, nelv, "Field")
      end if
    end associate

  end subroutine field_init_common

  !> Deallocate a field @a f
  subroutine field_free(t)
    class(field_t), intent(inout) :: t

    call t%tensor4_t%free()

    if (t%internal_dofmap) then
       call t%dof%free()
       deallocate(t%dof)
       t%internal_dofmap = .false.
    end if

    nullify(t%msh)
    nullify(t%Xh)
    nullify(t%dof)

  end subroutine field_free

  !> Seals off the inherited raw dimension initialiser (tensor4_t%init_dims)
  !! so it cannot be called on a field_t, which requires Xh/msh/dof to be
  !! set up and must always go through init_external_dof/init_internal_dof.
  subroutine field_init_dims_seal(t, n1, n2, n3, n4, name)
    class(field_t), intent(inout) :: t
    integer, intent(in) :: n1
    integer, intent(in) :: n2
    integer, intent(in) :: n3
    integer, intent(in) :: n4
    character(len=*), intent(in), optional :: name

    call neko_error('field_t must be initialised via init(msh, space, ' // &
         'name) or init(dof, name), not the inherited tensor4 dimension ' // &
         'initialiser')

  end subroutine field_init_dims_seal

  !> Seals off the inherited raw tensor4_t assignment so `field = tensor4`
  !! cannot silently bypass Xh/msh/dof; use assign_field (field = field)
  !! instead.
  subroutine field_assign_tensor4_seal(t, w)
    class(field_t), intent(inout) :: t
    type(tensor4_t), intent(in) :: w

    call neko_error('field_t must be assigned from another field_t, not ' // &
         'a bare tensor4_t')

  end subroutine field_assign_tensor4_seal

  !> Assignment \f$ this = G \f$
  !! @note @a this will be initialized if it has a different size than
  !! @a G or it's not allocated
  subroutine field_assign_field(this, g)
    class(field_t), intent(inout) :: this
    type(field_t), intent(in) :: g

    if (allocated(this%x)) then
       if (.not. associated(this%Xh, g%Xh)) then
          call this%free()
       end if
    end if

    this%Xh => g%Xh
    this%msh => g%msh
    if (len_trim(this%name) == 0) then
       this%name = g%name
    end if

    if (.not. g%internal_dofmap) then
       if (this%internal_dofmap) then
          call this%dof%free()
          deallocate(this%dof)
          this%internal_dofmap = .false.
       end if
       this%dof => g%dof
    else
       if (this%internal_dofmap) then
          call this%dof%free()
       else
          allocate(this%dof)
          this%internal_dofmap = .true.
       end if
       call this%dof%init(this%msh, this%Xh)
    end if

    if (.not. allocated(this%x)) then
       ! Delegates to the parent's alloc+zero-fill+device_map. The name is
       ! copied to a local first: tensor4_t%init resets t%name to "" (via
       ! its internal free()) before applying the name argument, and
       ! passing `this%name` directly would alias that same memory,
       ! clobbering the value set above before init could read it.
       block
         character(len=NEKO_VARNAME_LEN) :: fld_name
         fld_name = this%name
         call this%tensor4_t%init(this%Xh%lx, this%Xh%ly, this%Xh%lz, &
              this%msh%nelv, fld_name)
       end block
    end if

    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_copy(this%x_d, g%x_d, this%size())
    else
       call copy(this%x, g%x, this%dof%size())
    end if

  end subroutine field_assign_field

  !> Add \f$ this(u_1, u_2, ... , u_n) =
  !! this(u_1, u_2, ... , u_n) + G(u_1, u_2, ... , u_n) \f$
  !! @note Component wise
  subroutine field_add_field(this, g)
    class(field_t), intent(inout) :: this
    type(field_t), intent(in) :: g

    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_add2(this%x_d, g%x_d, this%size())
    else
       call add2(this%x, g%x, this%size())
    end if

  end subroutine field_add_field


  !> Add \f$ this(u_1, u_2, ... , u_n) =
  !! this(u_1, u_2, ... , u_n) + a \f$
  subroutine field_add_scalar(this, a)
    class(field_t), intent(inout) :: this
    real(kind=rp), intent(in) :: a

    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_cadd(this%x_d, a, this%size())
    else
       call cadd(this%x, a, this%size())
    end if

  end subroutine field_add_scalar

  ! ========================================================================== !
  ! Field pointer type subroutines

  subroutine field_ptr_init(this, ptr)
    class(field_ptr_t), intent(inout) :: this
    type(field_t), target, intent(in) :: ptr

    call this%free()

    this%ptr => ptr

  end subroutine field_ptr_init

  subroutine field_ptr_free(this)
    class(field_ptr_t), intent(inout) :: this

    if (associated(this%ptr)) then
       nullify(this%ptr)
    end if

  end subroutine field_ptr_free

  ! ========================================================================== !
  ! Field wrapper type subroutines

  subroutine field_wrapper_init_field(this, f)
    class(field_wrapper_t), intent(inout) :: this
    type(field_t), intent(in) :: f

    call this%free()
    allocate(this%field)
    this%field = f

  end subroutine field_wrapper_init_field

  subroutine field_wrapper_init_internal_dof(this, msh, space, fld_name)
    class(field_wrapper_t), intent(inout) :: this
    type(mesh_t), target, intent(in) :: msh
    type(space_t), target, intent(in) :: space
    character(len=*), optional :: fld_name

    call this%free()
    allocate(this%field)
    call this%field%init(msh, space, fld_name)

  end subroutine field_wrapper_init_internal_dof

  subroutine field_wrapper_init_external_dof(this, dof, fld_name)
    class(field_wrapper_t), intent(inout) :: this
    type(dofmap_t), target, intent(in) :: dof
    character(len=*), optional :: fld_name

    call this%free()
    allocate(this%field)
    call this%field%init(dof, fld_name)

  end subroutine field_wrapper_init_external_dof

  subroutine field_wrapper_free(this)
    class(field_wrapper_t), intent(inout) :: this

    if (associated(this%field)) then
       call this%field%free()
       deallocate(this%field)
    end if

  end subroutine field_wrapper_free

end module field
