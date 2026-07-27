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
!     copyright notice, this list of conditions and the following disclaimer
!     in the documentation and/or other materials provided with the
!     distribution.
!
!   * Neither the name of Neko nor the names of its contributors may be used
!     to endorse or promote products derived from this software without
!     specific prior written permission.
!
! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
! AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
! IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
! ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
! LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
! CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
! SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
! INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
! CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
! ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
! POSSIBILITY OF SUCH DAMAGE.
!
!> Format-independent checkpoint payloads.
module checkpoint_payload
  use num_types, only : rp, i8
  use field, only : field_t, field_ptr_t
  use field_series, only : field_series_t, field_series_ptr_t
  use mesh, only : mesh_t
  use space, only : space_t
  use utils, only : neko_error
  use, intrinsic :: iso_c_binding, only : c_ptr, c_null_ptr, c_loc, c_f_pointer
  implicit none
  private

  !> A named real array and its selection in a global checkpoint dataset.
  type, public :: checkpoint_array_t
     !> Dataset name within the payload.
     character(len=:), allocatable :: name
     !> Contiguous local host storage.
     real(kind=rp), pointer :: x(:) => null()
     !> Optional device storage corresponding to `x`.
     type(c_ptr) :: x_d = c_null_ptr
     !> Number of entries in the global dataset.
     integer(kind=i8) :: global_count = 0_i8
     !> Zero-based offset of the local selection in the global dataset.
     integer(kind=i8) :: offset = 0_i8
     !> Whether every rank owns an identical copy of the array.
     logical :: replicated = .false.
  end type checkpoint_array_t

  !> Pointer wrapper around `checkpoint_array_t`.
  type, public :: checkpoint_array_ptr_t
     !> Registered array descriptor.
     type(checkpoint_array_t), pointer :: ptr => null()
  end type checkpoint_array_ptr_t

  !> A named nodal real array distributed over mesh elements.
  type, public :: checkpoint_mesh_array_t
     !> Dataset name within the payload.
     character(len=:), allocatable :: name
     !> Contiguous local host storage.
     real(kind=rp), pointer :: x(:) => null()
     !> Optional device storage corresponding to `x`.
     type(c_ptr) :: x_d = c_null_ptr
     !> Mesh defining the element distribution.
     type(mesh_t), pointer :: msh => null()
     !> Function space defining the nodal layout on each element.
     type(space_t), pointer :: Xh => null()
  end type checkpoint_mesh_array_t

  !> Pointer wrapper around `checkpoint_mesh_array_t`.
  type, public :: checkpoint_mesh_array_ptr_t
     !> Registered mesh-array descriptor.
     type(checkpoint_mesh_array_t), pointer :: ptr => null()
  end type checkpoint_mesh_array_ptr_t

  !> A named collection of live fields to checkpoint together.
  type, public :: checkpoint_payload_t
     !> Payload path below the HDF5 `Payloads` group.
     character(len=:), allocatable :: name
     !> Registered individual fields.
     type(field_ptr_t), allocatable :: fields(:)
     !> Registered field series.
     type(field_series_ptr_t), allocatable :: series(:)
     !> Registered generic arrays.
     type(checkpoint_array_ptr_t), allocatable :: arrays(:)
     !> Registered nodal arrays distributed over mesh elements.
     type(checkpoint_mesh_array_ptr_t), allocatable :: mesh_arrays(:)
   contains
     !> Initialize a payload.
     procedure, pass(this) :: init => checkpoint_payload_init
     !> Release all pointers held by a payload.
     procedure, pass(this) :: free => checkpoint_payload_free
     !> Add a field using its native name as the dataset name.
     procedure, pass(this) :: add_field => checkpoint_payload_add_field
     !> Add a field series using the native names of its lag fields.
     procedure, pass(this) :: add_series => checkpoint_payload_add_series
     !> Add a rank-one real array.
     procedure, pass(this), private :: add_array_1d
     !> Add a rank-two contiguous real array.
     procedure, pass(this), private :: add_array_2d
     !> Add a rank-three contiguous real array.
     procedure, pass(this), private :: add_array_3d
     !> Add a rank-four contiguous real array.
     procedure, pass(this), private :: add_array_4d
     !> Add a contiguous real array of rank one to four.
     generic :: add_array => add_array_1d, add_array_2d, add_array_3d, &
          add_array_4d
     !> Add a four-dimensional nodal mesh array.
     procedure, pass(this) :: add_mesh_array => add_mesh_array_4d
     !> Return the number of fields in a payload.
     procedure, pass(this) :: field_count => checkpoint_payload_field_count
     !> Return the number of field series in a payload.
     procedure, pass(this) :: series_count => checkpoint_payload_series_count
     !> Return the number of arrays in a payload.
     procedure, pass(this) :: array_count => checkpoint_payload_array_count
     !> Return the number of mesh arrays in a payload.
     procedure, pass(this) :: mesh_array_count => &
          checkpoint_payload_mesh_array_count
     !> Find a field by its native name.
     procedure, pass(this) :: find_field => checkpoint_payload_find_field
     !> Find a field series by the native name of its base field.
     procedure, pass(this) :: find_series => checkpoint_payload_find_series
     !> Find an array descriptor by its dataset name.
     procedure, pass(this) :: find_array => checkpoint_payload_find_array
     !> Find a mesh-array descriptor by its dataset name.
     procedure, pass(this) :: find_mesh_array => &
          checkpoint_payload_find_mesh_array
  end type checkpoint_payload_t

  !> Pointer wrapper used to keep payload addresses stable as the list grows.
  type, public :: checkpoint_payload_ptr_t
     !> Registered checkpoint payload.
     type(checkpoint_payload_t), pointer :: ptr => null()
  end type checkpoint_payload_ptr_t

contains

  !> Initialize a payload.
  !! @param name Payload path below the HDF5 `Payloads` group.
  subroutine checkpoint_payload_init(this, name)
    class(checkpoint_payload_t), intent(inout) :: this
    character(len=*), intent(in) :: name

    call this%free()

    if (len_trim(name) .eq. 0) then
       call neko_error("Checkpoint payload name cannot be empty")
    end if
    if (name(1:1) .eq. "/" .or. &
         name(len_trim(name):len_trim(name)) .eq. "/") then
       call neko_error("Checkpoint payload names cannot start or end with '/'")
    end if
    if (index(trim(name), "//") .gt. 0) then
       call neko_error("Checkpoint payload names cannot contain empty groups")
    end if

    this%name = trim(name)

  end subroutine checkpoint_payload_init

  !> Release all pointers held by a payload.
  subroutine checkpoint_payload_free(this)
    class(checkpoint_payload_t), intent(inout) :: this
    integer :: i

    if (allocated(this%fields)) then
       do i = 1, size(this%fields)
          if (associated(this%fields(i)%ptr)) nullify(this%fields(i)%ptr)
       end do
       deallocate(this%fields)
    end if

    if (allocated(this%series)) then
       do i = 1, size(this%series)
          if (associated(this%series(i)%ptr)) nullify(this%series(i)%ptr)
       end do
       deallocate(this%series)
    end if

    if (allocated(this%arrays)) then
       do i = 1, size(this%arrays)
          if (associated(this%arrays(i)%ptr)) then
             if (associated(this%arrays(i)%ptr%x)) then
                nullify(this%arrays(i)%ptr%x)
             end if
             this%arrays(i)%ptr%x_d = c_null_ptr
             deallocate(this%arrays(i)%ptr)
          end if
       end do
       deallocate(this%arrays)
    end if

    if (allocated(this%mesh_arrays)) then
       do i = 1, size(this%mesh_arrays)
          if (associated(this%mesh_arrays(i)%ptr)) then
             if (associated(this%mesh_arrays(i)%ptr%x)) then
                nullify(this%mesh_arrays(i)%ptr%x)
             end if
             this%mesh_arrays(i)%ptr%x_d = c_null_ptr
             nullify(this%mesh_arrays(i)%ptr%msh)
             nullify(this%mesh_arrays(i)%ptr%Xh)
             deallocate(this%mesh_arrays(i)%ptr)
          end if
       end do
       deallocate(this%mesh_arrays)
    end if

    if (allocated(this%name)) deallocate(this%name)

  end subroutine checkpoint_payload_free

  !> Add a field using its native name as the dataset name.
  !! @param fld Field to register.
  subroutine checkpoint_payload_add_field(this, fld)
    class(checkpoint_payload_t), intent(inout) :: this
    type(field_t), target, intent(in) :: fld
    type(field_ptr_t), allocatable :: tmp(:)
    integer :: i, n

    call checkpoint_payload_validate_name(this, fld%name)

    n = this%field_count()
    do i = 1, n
       if (associated(this%fields(i)%ptr, fld)) then
          call neko_error( &
               "Field already registered in checkpoint payload '" // &
               trim(this%name) // "'")
       end if
    end do

    allocate(tmp(n + 1))
    if (n .gt. 0) tmp(1:n) = this%fields
    tmp(n + 1)%ptr => fld
    call move_alloc(tmp, this%fields)

  end subroutine checkpoint_payload_add_field

  !> Add a field series using the native names of its lag fields.
  !! @param fld_series Field series to register.
  subroutine checkpoint_payload_add_series(this, fld_series)
    class(checkpoint_payload_t), intent(inout) :: this
    type(field_series_t), target, intent(in) :: fld_series
    type(field_series_ptr_t), allocatable :: tmp(:)
    integer :: i, j, n

    do j = 1, fld_series%size()
       call checkpoint_payload_validate_name(this, fld_series%lf(j)%name)
       do i = 1, j - 1
          if (trim(fld_series%lf(i)%name) .eq. &
               trim(fld_series%lf(j)%name)) then
             call neko_error("Duplicate field name in checkpoint series")
          end if
       end do
    end do

    n = this%series_count()
    do i = 1, n
       if (associated(this%series(i)%ptr, fld_series)) then
          call neko_error("Field series already registered in checkpoint " // &
               "payload '" // trim(this%name) // "'")
       end if
    end do

    allocate(tmp(n + 1))
    if (n .gt. 0) tmp(1:n) = this%series
    tmp(n + 1)%ptr => fld_series
    call move_alloc(tmp, this%series)

  end subroutine checkpoint_payload_add_series

  !> Add a rank-one real array.
  !! @param name Dataset name within the payload.
  !! @param x Contiguous local array storage.
  !! @param[optional] global_count Number of entries in the global dataset.
  !! @param[optional] offset Zero-based offset of `x` in the global dataset.
  !! @param[optional] device_ptr Device storage corresponding to `x`.
  !! @param[optional] replicated Whether every rank owns a copy of `x`.
  subroutine add_array_1d(this, name, x, global_count, offset, &
       device_ptr, replicated)
    class(checkpoint_payload_t), intent(inout) :: this
    character(len=*), intent(in) :: name
    real(kind=rp), contiguous, target, intent(inout) :: x(:)
    integer(kind=i8), intent(in), optional :: global_count, offset
    type(c_ptr), intent(in), optional :: device_ptr
    logical, intent(in), optional :: replicated
    real(kind=rp), pointer :: flat(:)

    call c_f_pointer(c_loc(x), flat, [size(x)])
    call checkpoint_payload_add_array(this, name, flat, global_count, &
         offset, device_ptr, replicated)

  end subroutine add_array_1d

  !> Add a rank-two contiguous real array.
  !! @param name Dataset name within the payload.
  !! @param x Contiguous local array storage.
  !! @param[optional] global_count Number of entries in the global dataset.
  !! @param[optional] offset Zero-based offset of `x` in the global dataset.
  !! @param[optional] device_ptr Device storage corresponding to `x`.
  !! @param[optional] replicated Whether every rank owns a copy of `x`.
  subroutine add_array_2d(this, name, x, global_count, offset, &
       device_ptr, replicated)
    class(checkpoint_payload_t), intent(inout) :: this
    character(len=*), intent(in) :: name
    real(kind=rp), contiguous, target, intent(inout) :: x(:,:)
    integer(kind=i8), intent(in), optional :: global_count, offset
    type(c_ptr), intent(in), optional :: device_ptr
    logical, intent(in), optional :: replicated
    real(kind=rp), pointer :: flat(:)

    call c_f_pointer(c_loc(x), flat, [size(x)])
    call checkpoint_payload_add_array(this, name, flat, global_count, &
         offset, device_ptr, replicated)

  end subroutine add_array_2d

  !> Add a rank-three contiguous real array.
  !! @param name Dataset name within the payload.
  !! @param x Contiguous local array storage.
  !! @param[optional] global_count Number of entries in the global dataset.
  !! @param[optional] offset Zero-based offset of `x` in the global dataset.
  !! @param[optional] device_ptr Device storage corresponding to `x`.
  !! @param[optional] replicated Whether every rank owns a copy of `x`.
  subroutine add_array_3d(this, name, x, global_count, offset, &
       device_ptr, replicated)
    class(checkpoint_payload_t), intent(inout) :: this
    character(len=*), intent(in) :: name
    real(kind=rp), contiguous, target, intent(inout) :: x(:,:,:)
    integer(kind=i8), intent(in), optional :: global_count, offset
    type(c_ptr), intent(in), optional :: device_ptr
    logical, intent(in), optional :: replicated
    real(kind=rp), pointer :: flat(:)

    call c_f_pointer(c_loc(x), flat, [size(x)])
    call checkpoint_payload_add_array(this, name, flat, global_count, &
         offset, device_ptr, replicated)

  end subroutine add_array_3d

  !> Add a rank-four contiguous real array.
  !! @param name Dataset name within the payload.
  !! @param x Contiguous local array storage.
  !! @param[optional] global_count Number of entries in the global dataset.
  !! @param[optional] offset Zero-based offset of `x` in the global dataset.
  !! @param[optional] device_ptr Device storage corresponding to `x`.
  !! @param[optional] replicated Whether every rank owns a copy of `x`.
  subroutine add_array_4d(this, name, x, global_count, offset, &
       device_ptr, replicated)
    class(checkpoint_payload_t), intent(inout) :: this
    character(len=*), intent(in) :: name
    real(kind=rp), contiguous, target, intent(inout) :: x(:,:,:,:)
    integer(kind=i8), intent(in), optional :: global_count, offset
    type(c_ptr), intent(in), optional :: device_ptr
    logical, intent(in), optional :: replicated
    real(kind=rp), pointer :: flat(:)

    call c_f_pointer(c_loc(x), flat, [size(x)])
    call checkpoint_payload_add_array(this, name, flat, global_count, &
         offset, device_ptr, replicated)

  end subroutine add_array_4d

  !> Add a four-dimensional nodal mesh array.
  !! @param name Dataset name within the payload.
  !! @param x Array shaped as `(lx, ly, lz, nelv)`.
  !! @param msh Mesh defining the element distribution.
  !! @param Xh Function space defining the element-local nodal layout.
  !! @param[optional] device_ptr Device storage corresponding to `x`.
  subroutine add_mesh_array_4d(this, name, x, msh, Xh, device_ptr)
    class(checkpoint_payload_t), intent(inout) :: this
    character(len=*), intent(in) :: name
    real(kind=rp), contiguous, target, intent(inout) :: x(:,:,:,:)
    type(mesh_t), target, intent(in) :: msh
    type(space_t), target, intent(in) :: Xh
    type(c_ptr), intent(in), optional :: device_ptr
    real(kind=rp), pointer :: flat(:)

    if (any(shape(x) .ne. [Xh%lx, Xh%ly, Xh%lz, msh%nelv])) then
       call neko_error("Checkpoint mesh array must have shape " // &
            "(lx, ly, lz, nelv)")
    end if
    call c_f_pointer(c_loc(x), flat, [size(x)])
    call checkpoint_payload_add_mesh_array(this, name, flat, msh, Xh, &
         device_ptr)

  end subroutine add_mesh_array_4d

  !> Store a flattened real array in a payload.
  !! @param name Dataset name within the payload.
  !! @param x Contiguous local array storage.
  !! @param[optional] global_count Number of entries in the global dataset.
  !! @param[optional] offset Zero-based offset of `x` in the global dataset.
  !! @param[optional] device_ptr Device storage corresponding to `x`.
  !! @param[optional] replicated Whether every rank owns a copy of `x`.
  subroutine checkpoint_payload_add_array(this, name, x, global_count, &
       offset, device_ptr, replicated)
    class(checkpoint_payload_t), intent(inout) :: this
    character(len=*), intent(in) :: name
    real(kind=rp), pointer, intent(inout) :: x(:)
    integer(kind=i8), intent(in), optional :: global_count, offset
    type(c_ptr), intent(in), optional :: device_ptr
    logical, intent(in), optional :: replicated
    type(checkpoint_array_ptr_t), allocatable :: tmp(:)
    logical :: is_replicated
    integer :: n

    call checkpoint_payload_validate_name(this, name)

    is_replicated = .false.
    if (present(replicated)) is_replicated = replicated

    if (is_replicated) then
       if (present(global_count) .or. present(offset)) then
          call neko_error("Replicated checkpoint arrays do not take a " // &
               "global count or offset")
       end if
    else if (.not. present(global_count) .or. .not. present(offset)) then
       call neko_error("Distributed checkpoint arrays require a global " // &
            "count and offset")
    end if

    n = this%array_count()
    allocate(tmp(n + 1))
    if (n .gt. 0) tmp(1:n) = this%arrays
    allocate(tmp(n + 1)%ptr)
    tmp(n + 1)%ptr%name = trim(name)
    tmp(n + 1)%ptr%x => x
    tmp(n + 1)%ptr%replicated = is_replicated
    if (is_replicated) then
       tmp(n + 1)%ptr%global_count = int(size(x), i8)
       tmp(n + 1)%ptr%offset = 0_i8
    else
       tmp(n + 1)%ptr%global_count = global_count
       tmp(n + 1)%ptr%offset = offset
       if (tmp(n + 1)%ptr%offset .lt. 0_i8 .or. &
            tmp(n + 1)%ptr%global_count .lt. int(size(x), i8) .or. &
            tmp(n + 1)%ptr%offset + int(size(x), i8) .gt. &
            tmp(n + 1)%ptr%global_count) then
          call neko_error("Invalid checkpoint array selection")
       end if
    end if
    if (present(device_ptr)) tmp(n + 1)%ptr%x_d = device_ptr
    call move_alloc(tmp, this%arrays)

  end subroutine checkpoint_payload_add_array

  !> Store a flattened nodal mesh array in a payload.
  !! @param name Dataset name within the payload.
  !! @param x Contiguous local array storage.
  !! @param msh Mesh defining the element distribution.
  !! @param Xh Function space defining the element-local nodal layout.
  !! @param[optional] device_ptr Device storage corresponding to `x`.
  subroutine checkpoint_payload_add_mesh_array(this, name, x, msh, Xh, &
       device_ptr)
    class(checkpoint_payload_t), intent(inout) :: this
    character(len=*), intent(in) :: name
    real(kind=rp), pointer, intent(inout) :: x(:)
    type(mesh_t), target, intent(in) :: msh
    type(space_t), target, intent(in) :: Xh
    type(c_ptr), intent(in), optional :: device_ptr
    type(checkpoint_mesh_array_ptr_t), allocatable :: tmp(:)
    integer :: n

    call checkpoint_payload_validate_name(this, name)

    if (size(x) .ne. msh%nelv * Xh%lxyz) then
       call neko_error( &
            "Checkpoint mesh array size does not match its mesh " // &
            "and function space")
    end if

    n = this%mesh_array_count()
    allocate(tmp(n + 1))
    if (n .gt. 0) tmp(1:n) = this%mesh_arrays
    allocate(tmp(n + 1)%ptr)
    tmp(n + 1)%ptr%name = trim(name)
    tmp(n + 1)%ptr%x => x
    tmp(n + 1)%ptr%msh => msh
    tmp(n + 1)%ptr%Xh => Xh
    if (present(device_ptr)) tmp(n + 1)%ptr%x_d = device_ptr
    call move_alloc(tmp, this%mesh_arrays)

  end subroutine checkpoint_payload_add_mesh_array

  !> Validate that a dataset name is non-empty and unique in a payload.
  !! @param name Dataset name to validate.
  subroutine checkpoint_payload_validate_name(this, name)
    class(checkpoint_payload_t), intent(in) :: this
    character(len=*), intent(in) :: name
    integer :: i, j

    if (len_trim(name) .eq. 0) then
       call neko_error("Checkpoint field name cannot be empty")
    end if
    if (index(trim(name), "/") .gt. 0) then
       call neko_error("Checkpoint field names cannot contain '/'")
    end if

    do i = 1, this%field_count()
       if (trim(this%fields(i)%ptr%name) .eq. trim(name)) then
          call neko_error("Duplicate field name '" // trim(name) // &
               "' in checkpoint payload '" // trim(this%name) // "'")
       end if
    end do

    do i = 1, this%series_count()
       do j = 1, this%series(i)%ptr%size()
          if (trim(this%series(i)%ptr%lf(j)%name) .eq. trim(name)) then
             call neko_error("Duplicate field name '" // trim(name) // &
                  "' in checkpoint payload '" // trim(this%name) // "'")
          end if
       end do
    end do

    do i = 1, this%array_count()
       if (trim(this%arrays(i)%ptr%name) .eq. trim(name)) then
          call neko_error("Duplicate field or array name '" // trim(name) // &
               "' in checkpoint payload '" // trim(this%name) // "'")
       end if
    end do

    do i = 1, this%mesh_array_count()
       if (trim(this%mesh_arrays(i)%ptr%name) .eq. trim(name)) then
          call neko_error("Duplicate field or array name '" // trim(name) // &
               "' in checkpoint payload '" // trim(this%name) // "'")
       end if
    end do

  end subroutine checkpoint_payload_validate_name

  !> Return the number of fields in a payload.
  !! @return Number of registered fields.
  pure function checkpoint_payload_field_count(this) result(n)
    class(checkpoint_payload_t), intent(in) :: this
    integer :: n

    if (allocated(this%fields)) then
       n = size(this%fields)
    else
       n = 0
    end if

  end function checkpoint_payload_field_count

  !> Return the number of field series in a payload.
  !! @return Number of registered field series.
  pure function checkpoint_payload_series_count(this) result(n)
    class(checkpoint_payload_t), intent(in) :: this
    integer :: n

    if (allocated(this%series)) then
       n = size(this%series)
    else
       n = 0
    end if

  end function checkpoint_payload_series_count

  !> Return the number of arrays in a payload.
  !! @return Number of registered arrays.
  pure function checkpoint_payload_array_count(this) result(n)
    class(checkpoint_payload_t), intent(in) :: this
    integer :: n

    if (allocated(this%arrays)) then
       n = size(this%arrays)
    else
       n = 0
    end if

  end function checkpoint_payload_array_count

  !> Return the number of mesh arrays in a payload.
  !! @return Number of registered mesh arrays.
  pure function checkpoint_payload_mesh_array_count(this) result(n)
    class(checkpoint_payload_t), intent(in) :: this
    integer :: n

    if (allocated(this%mesh_arrays)) then
       n = size(this%mesh_arrays)
    else
       n = 0
    end if

  end function checkpoint_payload_mesh_array_count

  !> Find a field by its native name.
  !! @param name Native field name to find.
  !! @return Pointer to the field, or a disassociated pointer if absent.
  function checkpoint_payload_find_field(this, name) result(fld)
    class(checkpoint_payload_t), intent(in) :: this
    character(len=*), intent(in) :: name
    type(field_t), pointer :: fld
    integer :: i

    nullify(fld)
    do i = 1, this%field_count()
       if (trim(this%fields(i)%ptr%name) .eq. trim(name)) then
          fld => this%fields(i)%ptr
          return
       end if
    end do

  end function checkpoint_payload_find_field

  !> Find a field series by the native name of its base field.
  !! @param name Native base-field name to find.
  !! @return Pointer to the series, or a disassociated pointer if absent.
  function checkpoint_payload_find_series(this, name) result(series)
    class(checkpoint_payload_t), intent(in) :: this
    character(len=*), intent(in) :: name
    type(field_series_t), pointer :: series
    integer :: i

    nullify(series)
    do i = 1, this%series_count()
       if (trim(this%series(i)%ptr%f%name) .eq. trim(name)) then
          series => this%series(i)%ptr
          return
       end if
    end do

  end function checkpoint_payload_find_series

  !> Find an array descriptor by dataset name.
  !! @param name Dataset name to find.
  !! @return Pointer to the descriptor, or a disassociated pointer if absent.
  function checkpoint_payload_find_array(this, name) result(array)
    class(checkpoint_payload_t), intent(in) :: this
    character(len=*), intent(in) :: name
    type(checkpoint_array_t), pointer :: array
    integer :: i

    nullify(array)
    do i = 1, this%array_count()
       if (trim(this%arrays(i)%ptr%name) .eq. trim(name)) then
          array => this%arrays(i)%ptr
          return
       end if
    end do

  end function checkpoint_payload_find_array

  !> Find a mesh-array descriptor by dataset name.
  !! @param name Dataset name to find.
  !! @return Pointer to the descriptor, or a disassociated pointer if absent.
  function checkpoint_payload_find_mesh_array(this, name) result(array)
    class(checkpoint_payload_t), intent(in) :: this
    character(len=*), intent(in) :: name
    type(checkpoint_mesh_array_t), pointer :: array
    integer :: i

    nullify(array)
    do i = 1, this%mesh_array_count()
       if (trim(this%mesh_arrays(i)%ptr%name) .eq. trim(name)) then
          array => this%mesh_arrays(i)%ptr
          return
       end if
    end do

  end function checkpoint_payload_find_mesh_array

end module checkpoint_payload
