! Copyright (c) 2018-2023, The Neko Authors
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
  use device_math
  use num_types, only : rp
  use device
  use math, only : add2, copy, cadd
  use mesh, only : mesh_t
  use space, only : space_t, operator(.ne.)
  use dofmap, only : dofmap_t
  implicit none
  private

  type, public :: field_t
     real(kind=rp), allocatable :: x(:,:,:,:) !< Field data

     type(space_t), pointer :: Xh   !< Function space \f$ X_h \f$
     type(mesh_t), pointer :: msh   !< Mesh
     type(dofmap_t), pointer :: dof !< Dofmap

     logical :: internal_dofmap = .false. !< Does the field have an own dofmap
     character(len=80) :: name            !< Name of the field
     type(c_ptr) :: x_d = C_NULL_PTR
   contains
     procedure, private, pass(this) :: init_common => field_init_common
     procedure, private, pass(this) :: init_external_dof => &
          field_init_external_dof
     procedure, private, pass(this) :: init_internal_dof => &
          field_init_internal_dof
     procedure, private, pass(this) :: assign_field => field_assign_field
     procedure, private, pass(this) :: assign_scalar => field_assign_scalar
     procedure, private, pass(this) :: add_field => field_add_field
     procedure, private, pass(this) :: add_scalar => field_add_scalar
     procedure, pass(this) :: free => field_free
     !> Return the size of the field.
     procedure, pass(this) :: size => field_size
     !> Initialise a field
     generic :: init => init_external_dof, init_internal_dof
     !> Assignemnt to current field
     generic :: assignment(=) => assign_field, assign_scalar
     !> Add to current field
     !! @note We don't overload operator(+), to avoid
     !! the extra assignemnt operator
     generic :: add => add_field, add_scalar
  end type field_t

  !> field_ptr_t, To easily obtain a pointer to a field
  type, public ::  field_ptr_t
     type(field_t), pointer :: ptr => null()
  end type field_ptr_t

contains

  !> Initialize a field @a this on the mesh @a msh using an internal dofmap
  subroutine field_init_internal_dof(this, msh, space, fld_name)
    class(field_t), intent(inout) :: this      !< Field to be initialized
    type(mesh_t), target, intent(in) :: msh    !< underlying mesh of the field
    type(space_t), target, intent(in) :: space !< Function space for the field
    character(len=*), optional :: fld_name     !< Name of the field

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
    class(field_t), intent(inout) :: this      !< Field to be initialized
    type(dofmap_t), target, intent(in) :: dof  !< External dofmap for the field
    character(len=*), optional :: fld_name     !< Name of the field

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
  subroutine field_init_common(this, fld_name)
    class(field_t), intent(inout) :: this  !< Field to be initialized
    character(len=*), optional :: fld_name !< Name of the field
    integer :: ierr
    integer :: n

    associate(lx => this%Xh%lx, ly => this%Xh%ly, &
         lz => this%Xh%lz, nelv => this%msh%nelv)

      if (.not. allocated(this%x)) then
         allocate(this%x(lx, ly, lz, nelv), stat = ierr)
         this%x = 0d0
      end if

      if (present(fld_name)) then
         this%name = fld_name
      else
         this%name = "Field"
      end if

      if (NEKO_BCKND_DEVICE .eq. 1) then
         n = lx * ly * lz * nelv
         call device_map(this%x, this%x_d, n)
      end if
    end associate

  end subroutine field_init_common

  !> Deallocate a field @a f
  subroutine field_free(this)
    class(field_t), intent(inout) :: this

    if (allocated(this%x)) then
       deallocate(this%x)
    end if

    if (this%internal_dofmap) then
       deallocate(this%dof)
       this%internal_dofmap = .false.
    end if

    nullify(this%msh)
    nullify(this%Xh)
    nullify(this%dof)

    if (c_associated(this%x_d)) then
       call device_free(this%x_d)
    end if

  end subroutine field_free

  !> Assignment \f$ this = G \f$
  !! @note @a this will be initialized if it has a different size than
  !! @a G or it's not allocated
  subroutine field_assign_field(this, g)
    class(field_t), intent(inout) :: this
    type(field_t), intent(in) :: g

    if (allocated(this%x)) then
       if (this%Xh .ne. g%Xh) then
          call this%free()
       end if
    end if

    this%Xh => g%Xh
    this%msh => g%msh
    this%dof => g%dof


    this%Xh%lx = g%Xh%lx
    this%Xh%ly = g%Xh%ly
    this%Xh%lz = g%Xh%lz

    if (.not. allocated(this%x)) then

       allocate(this%x(this%Xh%lx, this%Xh%ly, this%Xh%lz, this%msh%nelv))

       if (NEKO_BCKND_DEVICE .eq. 1) then
          call device_map(this%x, this%x_d, this%size())
       end if

    end if

    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_copy(this%x_d, g%x_d, this%size())
    else
       call copy(this%x, g%x, this%dof%size())
    end if

  end subroutine field_assign_field

  !> Assignment \f$ this = a \f$
  subroutine field_assign_scalar(this, a)
    class(field_t), intent(inout) :: this
    real(kind=rp), intent(in) :: a
    integer :: i, j, k, l

    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_cfill(this%x_d, a, this%size())
    else
       do i = 1, this%msh%nelv
          do l = 1, this%Xh%lz
             do k = 1, this%Xh%ly
                do j = 1, this%Xh%lx
                   this%x(j, k, l, i) = a
                end do
             end do
          end do
       end do
    end if

  end subroutine field_assign_scalar

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

  !> Return the size of the field based on the underlying dofmap.
  pure function field_size(this) result(size)
    class(field_t), intent(in) :: this
    integer :: size

    size = this%dof%size()
  end function field_size

end module field

