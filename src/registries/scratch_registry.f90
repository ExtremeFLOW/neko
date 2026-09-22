! Copyright (c) 2025-2026, The Neko Authors
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
!> Defines a registry for storing and requesting temporary objects
!! This can be used when you have a function that will be called
!! often and you don'ptr want to create temporary objects (work arrays) inside
!! it on each call.
module scratch_registry
  use num_types, only : rp
  use registry_entry, only : registry_entry_t
  use host_array, only : host_array_t
  use device_array, only : device_array_t
  use field, only : field_t
  use vector, only : vector_t
  use matrix, only : matrix_t
  use tensor3, only : tensor3_t
  use tensor4, only : tensor4_t

  use math, only : rzero
  use device_math, only : device_rzero
  use field_math, only : field_rzero
  use vector_math, only : vector_rzero
  use matrix_math, only : matrix_rzero

  use dofmap, only : dofmap_t
  use utils, only : neko_error
  use neko_config, only : NEKO_BCKND_DEVICE
  use, intrinsic :: iso_c_binding, only : c_ptr
  implicit none
  private

  type, public :: scratch_registry_t
     !> List of scratch objects
     type(registry_entry_t), private, allocatable :: entries(:)
     !> Tracks which objects are in use
     logical, private, allocatable :: inuse(:)
     !> Number of allocated objects
     integer, private :: n_entries = 0
     !> The size the objects array is increased by upon reallocation
     integer, private :: expansion_size = 10
     !> Dofmap
     type(dofmap_t), pointer :: dof => null()
   contains
     !> Constructor
     procedure, pass(this) :: init => scratch_registry_init
     !> Destructor
     procedure, pass(this) :: free => scratch_registry_free
     !> Expand the registry by expansion_size
     procedure, private, pass(this) :: expand

     !> Assign a dofmap to the scratch registry
     procedure, pass(this) :: set_dofmap => scratch_registry_set_dofmap
     !> Assign a new expansion size to the scratch registry
     procedure, pass(this) :: set_expansion_size => &
          scratch_registry_set_expansion_size

     !> Return the number of allocated entries in the registry
     procedure, pass(this) :: get_n_entries
     !> Return the number of entries currently in use
     procedure, pass(this) :: get_n_inuse
     !> Return the expansion size of the registry
     procedure, pass(this) :: get_expansion_size
     !> Return the current capacity of the registry (size of the entries array)
     procedure, pass(this) :: get_size
     !> Return the inuse status of a given index
     procedure, pass(this) :: get_inuse

     !> Get a new scratch host array
     procedure, pass(this) :: request_host_array
     procedure, pass(this) :: relinquish_host_array_single
     procedure, pass(this) :: relinquish_host_array_multiple
     !> Free a host array for later reuse
     generic :: relinquish_host_array => relinquish_host_array_single, &
          relinquish_host_array_multiple

     !> Get a new scratch device array
     procedure, pass(this) :: request_device_array
     procedure, pass(this) :: relinquish_device_array_single
     procedure, pass(this) :: relinquish_device_array_multiple
     !> Free a device array for later reuse
     generic :: relinquish_device_array => relinquish_device_array_single, &
          relinquish_device_array_multiple

     !> Get a new scratch vector
     procedure, pass(this) :: request_vector
     procedure, pass(this) :: relinquish_vector_single
     procedure, pass(this) :: relinquish_vector_multiple
     !> Free a vector for later reuse
     generic :: relinquish_vector => relinquish_vector_single, &
          relinquish_vector_multiple

     !> Get a new scratch matrix
     procedure, pass(this) :: request_matrix
     procedure, pass(this) :: relinquish_matrix_single
     procedure, pass(this) :: relinquish_matrix_multiple
     !> Free a matrix for later reuse
     generic :: relinquish_matrix => relinquish_matrix_single, &
          relinquish_matrix_multiple

     !> Get a new scratch tensor3
     procedure, pass(this) :: request_tensor3
     procedure, pass(this) :: relinquish_tensor3_single
     procedure, pass(this) :: relinquish_tensor3_multiple
     !> Free a tensor3 for later reuse
     generic :: relinquish_tensor3 => relinquish_tensor3_single, &
          relinquish_tensor3_multiple

     !> Get a new scratch tensor4
     procedure, pass(this) :: request_tensor4
     procedure, pass(this) :: relinquish_tensor4_single
     procedure, pass(this) :: relinquish_tensor4_multiple
     !> Free a tensor4 for later reuse
     generic :: relinquish_tensor4 => relinquish_tensor4_single, &
          relinquish_tensor4_multiple

     !> Get a new scratch field based on the stored dofmap
     procedure, pass(this) :: request_field_stored_dof
     !> Get a new scratch field based on a provided dofmap
     procedure, pass(this) :: request_field_free_dof
     procedure, pass(this) :: relinquish_field_single
     procedure, pass(this) :: relinquish_field_multiple
     !> Get a new scratch field
     generic :: request_field => request_field_stored_dof, &
          request_field_free_dof
     !> Free a field for later reuse
     generic :: relinquish_field => relinquish_field_single, &
          relinquish_field_multiple

     !> Generic request procedure
     generic :: request => request_host_array, request_device_array, &
          request_vector, request_matrix, request_tensor3, request_tensor4, &
          request_field_stored_dof, request_field_free_dof

     !> Generic relinquish procedure for single objects
     procedure, pass(this) :: relinquish_single
     !> Generic relinquish procedure for multiple objects
     procedure, pass(this) :: relinquish_multiple
     !> Generic relinquish procedure
     generic :: relinquish => relinquish_single, relinquish_multiple
  end type scratch_registry_t

  !> Global scratch registry
  type(scratch_registry_t), public, target :: neko_scratch_registry

contains

  !> Constructor, optionally taking initial registry and expansion
  !! size as argument
  !! @param size Initial size of the registry
  !! @param expansion_size Size to expand the registry by when needed
  !! @param dof Dofmap to associate with the scratch registry
  !!
  !! @note If no DOF map is provided here, it must be set later using
  !!       scratch_registry_t::set_dofmap before requesting fields.
  subroutine scratch_registry_init(this, size, expansion_size, dof)
    class(scratch_registry_t), intent(inout) :: this
    integer, optional, intent(in) :: size
    integer, optional, intent(in) :: expansion_size
    type(dofmap_t), target, intent(in), optional :: dof
    integer :: initial_size

    call this%free()

    ! Assign the dofmap and expansion size if provided
    if (present(dof)) this%dof => dof
    if (present(expansion_size)) this%expansion_size = expansion_size

    ! Allocate the entries and inuse arrays with the initial size
    initial_size = this%expansion_size
    if (present(size)) initial_size = size

    allocate(this%entries(initial_size))
    allocate(this%inuse(initial_size), source = .false.)

  end subroutine scratch_registry_init

  !> Destructor
  subroutine scratch_registry_free(this)
    class(scratch_registry_t), intent(inout):: this
    integer :: i

    if (this%get_n_inuse() .ne. 0) then
       call neko_error("scratch_registry::free: " // &
            "Cannot free scratch registry with entries still in use.")
    end if

    if (allocated(this%entries)) then
       do i = 1, this%n_entries
          call this%entries(i)%free()
       end do
       deallocate(this%entries)
    end if

    if (allocated(this%inuse)) deallocate(this%inuse)
    if (associated(this%dof)) nullify(this%dof)

    ! Reset to default values
    this%n_entries = 0
    this%expansion_size = 10

  end subroutine scratch_registry_free

  !> Assign a dofmap to the scratch registry.
  !! @param dof Dofmap to assign
  !! @note First check if a dofmap is already assigned, and throw an error,
  !!       unless it's the same dofmap.
  subroutine scratch_registry_set_dofmap(this, dof)
    class(scratch_registry_t), intent(inout) :: this
    type(dofmap_t), target, intent(in) :: dof

    if (associated(this%dof, dof)) then
       return
    else if (associated(this%dof)) then
       call neko_error("scratch_registry::set_dofmap: "&
            // "Dofmap is already assigned to scratch registry.")
    end if

    this%dof => dof
  end subroutine scratch_registry_set_dofmap

  !> Get the number of objects stored in the registry
  pure function get_n_entries(this) result(n)
    class(scratch_registry_t), intent(in) :: this
    integer :: n

    n = this%n_entries
  end function get_n_entries

  !> Get the number of objects currently in use
  pure function get_n_inuse(this) result(n)
    class(scratch_registry_t), intent(in) :: this
    integer :: n

    if (allocated(this%inuse)) then
       n = count(this%inuse)
    else
       n = 0
    end if
  end function get_n_inuse

  !> Get the size of the objects array
  pure function get_size(this) result(n)
    class(scratch_registry_t), intent(in) :: this
    integer :: n

    if (allocated(this%entries)) then
       n = size(this%entries)
    else
       n = 0
    end if
  end function get_size

  !> Get the expansion size
  pure function get_expansion_size(this) result(n)
    class(scratch_registry_t), intent(in) :: this
    integer :: n

    n = this%expansion_size
  end function get_expansion_size

  !> Get the inuse status for a given index
  pure logical function get_inuse(this, index)
    class(scratch_registry_t), target, intent(in) :: this
    integer, intent(in) :: index

    get_inuse = this%inuse(index)
  end function get_inuse

  subroutine expand(this)
    class(scratch_registry_t), intent(inout) :: this
    type(registry_entry_t), allocatable :: temp(:)
    logical, allocatable :: temp2(:)
    integer :: i, n

    n = this%get_size()

    if (n .gt. 0) then
       call move_alloc(this%entries, temp)
       call move_alloc(this%inuse, temp2)
    end if

    allocate(this%entries(n + this%expansion_size))
    allocate(this%inuse(n + this%expansion_size), source = .false.)

    if (n .gt. 0) then
       do i = 1, n
          call this%entries(i)%move_from(temp(i))
          this%inuse(i) = temp2(i)
          call temp(i)%free()
       end do
    end if

    if (allocated(temp)) deallocate(temp)
    if (allocated(temp2)) deallocate(temp2)

  end subroutine expand

  !> Get a host array from the registry by assigning it to a pointer.
  !! @param ptr Pointer to the requested host array.
  !! @param index Index of the host array in the registry (for
  !! relinquishing later).
  !! @param n Size of the requested host_array.
  !! @param clear If true, the host_array values are set to zero upon request.
  subroutine request_host_array(this, ptr, index, n, clear)
    class(scratch_registry_t), target, intent(inout) :: this
    real(kind=rp), pointer, dimension(:), intent(inout) :: ptr
    integer, intent(inout) :: index
    integer, intent(in) :: n
    logical, intent(in) :: clear
    type(host_array_t), pointer :: v_scratch

    associate(entries => this%entries, n_entries => this%n_entries)

      do index = 1, this%get_size()
         if (.not. this%inuse(index)) then

            if (.not. entries(index)%is_allocated()) then
               call entries(index)%init_host_array(n)
               n_entries = n_entries + 1
            else if (trim(entries(index)%get_type()) .ne. 'host_array') then
               cycle
            end if

            v_scratch => entries(index)%get_host_array()
            if (v_scratch%size() .ne. n) then
               nullify(v_scratch)
               cycle
            end if

            if (clear) call rzero(v_scratch%x, v_scratch%size())
            this%inuse(index) = .true.
            ptr => v_scratch%x
            nullify(v_scratch)
            return
         end if
      end do

      ! all existing host_arrays in use, we need to expand to add a new one
      index = n_entries + 1
      call this%expand()
      n_entries = n_entries + 1
      this%inuse(n_entries) = .true.
      call this%entries(n_entries)%init_host_array(n)
      v_scratch => this%entries(n_entries)%get_host_array()
      ptr => v_scratch%x
      nullify(v_scratch)

    end associate
  end subroutine request_host_array

  !> Get a device array from the registry by assigning it to a pointer.
  !! @param ptr Pointer to the requested device array.
  !! @param index Index of the device array in the registry (for
  !! relinquishing later).
  !! @param n Size of the requested device array.
  !! @param clear If true, the device array values are set to zero upon request.
  subroutine request_device_array(this, ptr, index, n, clear)
    class(scratch_registry_t), target, intent(inout) :: this
    type(c_ptr), intent(inout) :: ptr
    integer, intent(inout) :: index
    integer, intent(in) :: n
    logical, intent(in) :: clear
    type(device_array_t), pointer :: v_tmp

    associate(entries => this%entries, n_entries => this%n_entries)

      do index = 1, this%get_size()
         if (.not. this%inuse(index)) then

            if (.not. entries(index)%is_allocated()) then
               call entries(index)%init_device_array(n)
               n_entries = n_entries + 1
            else if (trim(entries(index)%get_type()) .ne. 'device_array') then
               cycle
            end if

            v_tmp => entries(index)%get_device_array()
            if (v_tmp%size() .ne. n) then
               nullify(v_tmp)
               cycle
            end if

            ptr = v_tmp%x_d
            if (clear) call device_rzero(ptr, n)
            this%inuse(index) = .true.
            nullify(v_tmp)
            return
         end if
      end do

      ! all existing device_arrays in use, we need to expand to add a new one
      index = n_entries + 1
      call this%expand()
      n_entries = n_entries + 1
      this%inuse(n_entries) = .true.
      call this%entries(n_entries)%init_device_array(n)
      v_tmp => this%entries(n_entries)%get_device_array()
      ptr = v_tmp%x_d
      nullify(v_tmp)

    end associate
  end subroutine request_device_array

  !> Get a vector from the registry by assigning it to a pointer.
  !! @param ptr Pointer to the requested vector.
  !! @param index Index of the vector in the registry (for relinquishing later).
  !! @param n Size of the requested vector.
  !! @param clear If true, the vector values are set to zero upon request.
  subroutine request_vector(this, ptr, index, n, clear)
    class(scratch_registry_t), target, intent(inout) :: this
    type(vector_t), pointer, intent(inout) :: ptr
    integer, intent(inout) :: index
    integer, intent(in) :: n
    logical, intent(in) :: clear

    associate(entries => this%entries, n_entries => this%n_entries)

      do index = 1, this%get_size()
         if (.not. this%inuse(index)) then

            if (.not. entries(index)%is_allocated()) then
               call entries(index)%init_vector(n)
               n_entries = n_entries + 1
            else if (trim(entries(index)%get_type()) .ne. 'vector') then
               cycle
            end if

            ptr => entries(index)%get_vector()
            if (ptr%size() .ne. n) then
               nullify(ptr)
               cycle
            end if

            if (clear) call vector_rzero(ptr)
            this%inuse(index) = .true.
            return
         end if
      end do

      ! all existing vectors in use, we need to expand to add a new one
      index = n_entries + 1
      call this%expand()
      n_entries = n_entries + 1
      this%inuse(n_entries) = .true.
      call this%entries(n_entries)%init_vector(n)
      ptr => this%entries(n_entries)%get_vector()

    end associate
  end subroutine request_vector

  !> Get a matrix from the registry by assigning it to a pointer.
  !! @param ptr Pointer to the requested matrix.
  !! @param index Index of the matrix in the registry (for relinquishing later).
  !! @param nrows Number of rows of the requested matrix.
  !! @param ncols Number of columns of the requested matrix.
  !! @param clear If true, the matrix values are set to zero upon request.
  subroutine request_matrix(this, ptr, index, nrows, ncols, clear)
    class(scratch_registry_t), target, intent(inout) :: this
    type(matrix_t), pointer, intent(inout) :: ptr
    integer, intent(inout) :: index
    integer, intent(in) :: nrows, ncols
    logical, intent(in) :: clear

    associate(entries => this%entries, n_entries => this%n_entries)

      do index = 1, this%get_size()
         if (.not. this%inuse(index)) then

            if (.not. entries(index)%is_allocated()) then
               call entries(index)%init_matrix(nrows, ncols)
               n_entries = n_entries + 1
            else if (trim(entries(index)%get_type()) .ne. 'matrix') then
               cycle
            end if

            ptr => entries(index)%get_matrix()
            if (ptr%get_nrows() .ne. nrows .or. &
                 ptr%get_ncols() .ne. ncols) then
               nullify(ptr)
               cycle
            end if

            if (clear) call matrix_rzero(ptr)
            this%inuse(index) = .true.
            return
         end if
      end do

      ! all existing matrices in use, we need to expand to add a new one
      index = n_entries + 1
      call this%expand()
      n_entries = n_entries + 1
      this%inuse(n_entries) = .true.
      call this%entries(n_entries)%init_matrix(nrows, ncols)
      ptr => this%entries(n_entries)%get_matrix()

    end associate
  end subroutine request_matrix

  !> Get a tensor3 from the registry by assigning it to a pointer.
  !! @param ptr Pointer to the requested tensor3.
  !! @param index Index of the tensor3 in the registry (for relinquishing later).
  !! @param n Number of rows of the requested tensor3.
  !! @param ptr Number of columns of the requested tensor3.
  !! @param l Number of layers of the requested tensor3.
  !! @param clear If true, the tensor3 values are set to zero upon request.
  subroutine request_tensor3(this, ptr, index, n, m, l, clear)
    class(scratch_registry_t), target, intent(inout) :: this
    type(tensor3_t), pointer, intent(inout) :: ptr
    integer, intent(inout) :: index
    integer, intent(in) :: n, m, l
    logical, intent(in) :: clear

    associate(entries => this%entries, n_entries => this%n_entries)

      do index = 1, this%get_size()
         if (.not. this%inuse(index)) then

            if (.not. entries(index)%is_allocated()) then
               call entries(index)%init_tensor3(n, m, l)
               n_entries = n_entries + 1
            else if (trim(entries(index)%get_type()) .ne. 'tensor3') then
               cycle
            end if

            ptr => entries(index)%get_tensor3()
            if (ptr%get_n1() .ne. n .or. &
                 ptr%get_n2() .ne. m .or. &
                 ptr%get_n3() .ne. l) then
               nullify(ptr)
               cycle
            end if

            if (clear .and. NEKO_BCKND_DEVICE .eq. 1) then
               call device_rzero(ptr%x_d, ptr%size())
            else if (clear) then
               call rzero(ptr%x, ptr%size())
            end if
            this%inuse(index) = .true.
            return
         end if
      end do

      ! all existing matrices in use, we need to expand to add a new one
      index = n_entries + 1
      call this%expand()
      n_entries = n_entries + 1
      this%inuse(n_entries) = .true.
      call this%entries(n_entries)%init_tensor3(n, m, l)
      ptr => this%entries(n_entries)%get_tensor3()

    end associate
  end subroutine request_tensor3

  !> Get a tensor4 from the registry by assigning it to a pointer.
  !! @param ptr Pointer to the requested tensor4.
  !! @param index Index of the tensor4 in the registry (for relinquishing later).
  !! @param n Number of rows of the requested tensor4.
  !! @param ptr Number of columns of the requested tensor4.
  !! @param l Number of layers of the requested tensor4.
  !! @param k Number of slices of the requested tensor4.
  !! @param clear If true, the tensor4 values are set to zero upon request.
  subroutine request_tensor4(this, ptr, index, n, m, l, k, clear)
    class(scratch_registry_t), target, intent(inout) :: this
    type(tensor4_t), pointer, intent(inout) :: ptr
    integer, intent(inout) :: index
    integer, intent(in) :: n, m, l, k
    logical, intent(in) :: clear

    associate(entries => this%entries, n_entries => this%n_entries)

      do index = 1, this%get_size()
         if (.not. this%inuse(index)) then

            if (.not. entries(index)%is_allocated()) then
               call entries(index)%init_tensor4(n, m, l, k)
               n_entries = n_entries + 1
            else if (trim(entries(index)%get_type()) .ne. 'tensor4') then
               cycle
            end if

            ptr => entries(index)%get_tensor4()
            if (ptr%get_n1() .ne. n .or. &
                 ptr%get_n2() .ne. m .or. &
                 ptr%get_n3() .ne. l .or. &
                 ptr%get_n4() .ne. k) then
               nullify(ptr)
               cycle
            end if

            if (clear .and. NEKO_BCKND_DEVICE .eq. 1) then
               call device_rzero(ptr%x_d, ptr%size())
            else if (clear) then
               call rzero(ptr%x, ptr%size())
            end if
            this%inuse(index) = .true.
            return
         end if
      end do

      ! all existing matrices in use, we need to expand to add a new one
      index = n_entries + 1
      call this%expand()
      n_entries = n_entries + 1
      this%inuse(n_entries) = .true.
      call this%entries(n_entries)%init_tensor4(n, m, l, k)
      ptr => this%entries(n_entries)%get_tensor4()

    end associate
  end subroutine request_tensor4

  !> Get a field from the registry by assigning it to a pointer
  !! @param ptr Pointer to the requested field.
  !! @param index Index of the field in the registry (for relinquishing later).
  !! @param clear If true, the field values are set to zero upon request.
  subroutine request_field_stored_dof(this, ptr, index, clear)
    class(scratch_registry_t), target, intent(inout) :: this
    type(field_t), pointer, intent(inout) :: ptr
    integer, intent(inout) :: index
    logical, intent(in) :: clear
    character(len=10) :: name

    if (.not. associated(this%dof)) then
       call neko_error("scratch_registry::request_field_stored_dof: "&
            // "No dofmap assigned to scratch registry.")
    end if

    associate(entries => this%entries, n_entries => this%n_entries)

      do index = 1, this%get_size()
         if (.not. this%inuse(index)) then

            if (.not. entries(index)%is_allocated()) then
               write(name, "(A3,I0.3)") "wrk", index
               call entries(index)%init_field(this%dof, trim(name))
               n_entries = n_entries + 1
            else if (entries(index)%get_type() .ne. 'field') then
               cycle
            end if

            ptr => entries(index)%get_field()
            if (clear) call field_rzero(ptr)
            this%inuse(index) = .true.
            return
         end if
      end do

      ! all existing fields in use, we need to expand to add a new one
      index = n_entries + 1
      call this%expand()
      n_entries = n_entries + 1
      this%inuse(n_entries) = .true.
      write (name, "(A3,I0.3)") "wrk", index
      call this%entries(n_entries)%init_field(this%dof, trim(name))
      ptr => this%entries(n_entries)%get_field()

    end associate
  end subroutine request_field_stored_dof

  !> Get a field from the registry by assigning it to a pointer
  !! @param ptr Pointer to the requested field.
  !! @param index Index of the field in the registry (for relinquishing later).
  !! @param dof Dofmap to use for the field.
  !! @param clear If true, the field values are set to zero upon request.
  subroutine request_field_free_dof(this, ptr, index, dof, clear)
    class(scratch_registry_t), target, intent(inout) :: this
    type(field_t), pointer, intent(inout) :: ptr
    integer, intent(inout) :: index
    type(dofmap_t), target, intent(in) :: dof
    logical, intent(in) :: clear
    character(len=10) :: name

    associate(entries => this%entries, n_entries => this%n_entries)

      do index = 1, this%get_size()
         if (.not. this%inuse(index)) then

            if (.not. entries(index)%is_allocated()) then
               write(name, "(A3,I0.3)") "wrk", index
               call entries(index)%init_field(dof, trim(name))
               n_entries = n_entries + 1
            else if (entries(index)%get_type() .ne. 'field') then
               cycle
            end if

            ptr => entries(index)%get_field()
            if (.not. associated(ptr%dof, dof)) then
               nullify(ptr)
               cycle
            end if

            if (clear) call field_rzero(ptr)
            this%inuse(index) = .true.
            return
         end if
      end do

      ! all existing fields in use, we need to expand to add a new one

      index = n_entries + 1
      call this%expand()
      n_entries = n_entries + 1
      this%inuse(n_entries) = .true.
      write (name, "(A3,I0.3)") "wrk", index
      call this%entries(n_entries)%init_field(dof, trim(name))
      ptr => this%entries(n_entries)%get_field()

    end associate
  end subroutine request_field_free_dof

  !> Relinquish the use of a host_array in the registry
  !! @param index The index of the host_array to free
  subroutine relinquish_host_array_single(this, index)
    class(scratch_registry_t), target, intent(inout) :: this
    integer, intent(inout) :: index

    if (trim(this%entries(index)%get_type()) .ne. 'host_array') then
       call neko_error("scratch_registry::relinquish_host_array_single: " &
            // "Register entry is not a host_array.")
    end if

    this%inuse(index) = .false.
  end subroutine relinquish_host_array_single

  !> Relinquish the use of multiple host_arrays in the registry
  !! @param indices The indices of the host_arrays to free
  subroutine relinquish_host_array_multiple(this, indices)
    class(scratch_registry_t), target, intent(inout) :: this
    integer, intent(inout) :: indices(:)
    integer :: i

    do i = 1, size(indices)
       if (trim(this%entries(indices(i))%get_type()) .ne. 'host_array') then
          call neko_error("scratch_registry::relinquish_host_array_single: " &
               // "Register entry is not a host_array.")
       end if

       this%inuse(indices(i)) = .false.
    end do
  end subroutine relinquish_host_array_multiple

  !> Relinquish the use of a device_array in the registry
  !! @param index The index of the device_array to free
  subroutine relinquish_device_array_single(this, index)
    class(scratch_registry_t), target, intent(inout) :: this
    integer, intent(inout) :: index

    if (trim(this%entries(index)%get_type()) .ne. 'device_array') then
       call neko_error("scratch_registry::relinquish_device_array_single: " &
            // "Register entry is not a device_array.")
    end if

    this%inuse(index) = .false.
  end subroutine relinquish_device_array_single

  !> Relinquish the use of multiple device_arrays in the registry
  !! @param indices The indices of the device_arrays to free
  subroutine relinquish_device_array_multiple(this, indices)
    class(scratch_registry_t), target, intent(inout) :: this
    integer, intent(inout) :: indices(:)
    integer :: i

    do i = 1, size(indices)
       if (trim(this%entries(indices(i))%get_type()) .ne. 'device_array') then
          call neko_error("scratch_registry::relinquish_device_array_single: " &
               // "Register entry is not a device_array.")
       end if

       this%inuse(indices(i)) = .false.
    end do
  end subroutine relinquish_device_array_multiple

  !> Relinquish the use of a vector in the registry
  !! @param index The index of the vector to free
  subroutine relinquish_vector_single(this, index)
    class(scratch_registry_t), target, intent(inout) :: this
    integer, intent(inout) :: index

    if (trim(this%entries(index)%get_type()) .ne. 'vector') then
       call neko_error("scratch_registry::relinquish_vector_single: " &
            // "Register entry is not a vector.")
    end if

    this%inuse(index) = .false.
  end subroutine relinquish_vector_single

  !> Relinquish the use of multiple vectors in the registry
  !! @param indices The indices of the vectors to free
  subroutine relinquish_vector_multiple(this, indices)
    class(scratch_registry_t), target, intent(inout) :: this
    integer, intent(inout) :: indices(:)
    integer :: i

    do i = 1, size(indices)
       if (trim(this%entries(indices(i))%get_type()) .ne. 'vector') then
          call neko_error("scratch_registry::relinquish_vector_single: " &
               // "Register entry is not a vector.")
       end if

       this%inuse(indices(i)) = .false.
    end do
  end subroutine relinquish_vector_multiple

  !> Relinquish the use of a matrix in the registry
  !! @param index The index of the matrix to free
  subroutine relinquish_matrix_single(this, index)
    class(scratch_registry_t), target, intent(inout) :: this
    integer, intent(inout) :: index

    if (trim(this%entries(index)%get_type()) .ne. 'matrix') then
       call neko_error("scratch_registry::relinquish_matrix_single: " &
            // "Register entry is not a matrix.")
    end if

    this%inuse(index) = .false.
  end subroutine relinquish_matrix_single

  !> Relinquish the use of multiple matrices in the registry
  !! @param indices The indices of the matrices to free
  subroutine relinquish_matrix_multiple(this, indices)
    class(scratch_registry_t), target, intent(inout) :: this
    integer, intent(inout) :: indices(:)
    integer :: i

    do i = 1, size(indices)
       if (trim(this%entries(indices(i))%get_type()) .ne. 'matrix') then
          call neko_error("scratch_registry::relinquish_matrix_single: " &
               // "Register entry is not a matrix.")
       end if

       this%inuse(indices(i)) = .false.
    end do
  end subroutine relinquish_matrix_multiple

  !> Relinquish the use of a tensor3 in the registry
  !! @param index The index of the tensor3 to free
  subroutine relinquish_tensor3_single(this, index)
    class(scratch_registry_t), target, intent(inout) :: this
    integer, intent(inout) :: index

    if (trim(this%entries(index)%get_type()) .ne. 'tensor3') then
       call neko_error("scratch_registry::relinquish_tensor3_single: " &
            // "Register entry is not a tensor3.")
    end if

    this%inuse(index) = .false.
  end subroutine relinquish_tensor3_single

  !> Relinquish the use of multiple tensor3s in the registry
  !! @param indices The indices of the tensor3s to free
  subroutine relinquish_tensor3_multiple(this, indices)
    class(scratch_registry_t), target, intent(inout) :: this
    integer, intent(inout) :: indices(:)
    integer :: i

    do i = 1, size(indices)
       if (trim(this%entries(indices(i))%get_type()) .ne. 'tensor3') then
          call neko_error("scratch_registry::relinquish_tensor3_single: " &
               // "Register entry is not a tensor3.")
       end if

       this%inuse(indices(i)) = .false.
    end do
  end subroutine relinquish_tensor3_multiple

  !> Relinquish the use of a tensor4 in the registry
  !! @param index The index of the tensor4 to free
  subroutine relinquish_tensor4_single(this, index)
    class(scratch_registry_t), target, intent(inout) :: this
    integer, intent(inout) :: index

    if (trim(this%entries(index)%get_type()) .ne. 'tensor4') then
       call neko_error("scratch_registry::relinquish_tensor4_single: " &
            // "Register entry is not a tensor4.")
    end if

    this%inuse(index) = .false.
  end subroutine relinquish_tensor4_single

  !> Relinquish the use of multiple tensor4s in the registry
  !! @param indices The indices of the tensor4s to free
  subroutine relinquish_tensor4_multiple(this, indices)
    class(scratch_registry_t), target, intent(inout) :: this
    integer, intent(inout) :: indices(:)
    integer :: i

    do i = 1, size(indices)
       if (trim(this%entries(indices(i))%get_type()) .ne. 'tensor4') then
          call neko_error("scratch_registry::relinquish_tensor4_single: " &
               // "Register entry is not a tensor4.")
       end if

       this%inuse(indices(i)) = .false.
    end do
  end subroutine relinquish_tensor4_multiple

  !> Relinquish the use of a field in the registry
  !! @param index The index of the field to free
  subroutine relinquish_field_single(this, index)
    class(scratch_registry_t), target, intent(inout) :: this
    integer, intent(inout) :: index

    if (trim(this%entries(index)%get_type()) .ne. 'field') then
       call neko_error("scratch_registry::relinquish_field_single: " &
            // "Register entry is not a field.")
    end if

    this%inuse(index) = .false.
  end subroutine relinquish_field_single

  !> Relinquish the use of multiple fields in the registry
  !! @param indices The indices of the fields to free
  subroutine relinquish_field_multiple(this, indices)
    class(scratch_registry_t), target, intent(inout) :: this
    integer, intent(inout) :: indices(:)
    integer :: i

    do i = 1, size(indices)
       if (trim(this%entries(indices(i))%get_type()) .ne. 'field') then
          call neko_error("scratch_registry::relinquish_field_single: " &
               // "Register entry is not a field.")
       end if

       this%inuse(indices(i)) = .false.
    end do
  end subroutine relinquish_field_multiple

  !> Relinquish the use of an object in the registry
  !! @param index The index of the object to free
  subroutine relinquish_single(this, index)
    class(scratch_registry_t), target, intent(inout) :: this
    integer, intent(inout) :: index

    this%inuse(index) = .false.
  end subroutine relinquish_single

  !> Relinquish the use of multiple objects in the registry
  !! @param indices The indices of the objects to free
  subroutine relinquish_multiple(this, indices)
    class(scratch_registry_t), target, intent(inout) :: this
    integer, intent(inout) :: indices(:)
    integer :: i

    do i = 1, size(indices)
       this%inuse(indices(i)) = .false.
    end do
  end subroutine relinquish_multiple

end module scratch_registry
