! Copyright (c) 2023, The Neko Authors
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
!> Defines a matrix
module matrix
  use neko_config
  use num_types
  use device
  use device_math
  use utils
  use, intrinsic :: iso_c_binding
  implicit none
  private

  type, public ::  matrix_t
     real(kind=rp), allocatable :: x(:,:) !< Matrix entries.
     type(c_ptr) :: x_d = C_NULL_PTR      !< Device pointer.
     integer :: nrows  = 0 !< Number of matrix rows.
     integer :: ncols  = 0 !< Number of matrix columns.
     integer :: n = 0      !< Total size nows*ncols.
   contains
     !> Initialise a matrix of size `nrows*ncols`.
     procedure, pass(m) :: init => matrix_init
     !> Deallocate a matrix.
     procedure, pass(m) :: free => matrix_free
     !> Returns the number of entries in the matrix.
     procedure, pass(m) :: size => matrix_size
     !> Assignment \f$ m = w \f$
     procedure, pass(m) :: matrix_assign_matrix
     !> Assignment \f$ m = s \f$.
     procedure, pass(m) :: matrix_assign_scalar
     generic :: assignment(=) => matrix_assign_matrix, &
          matrix_assign_scalar
  end type matrix_t
  
contains

  !> Initialise a matrix of size `nrows*ncols`.
  !! @param nrows Number of rows.
  !! @param ncols Number of columns.
  subroutine matrix_init(m, nrows, ncols)
    class(matrix_t), intent(inout) :: m
    integer, intent(in) :: nrows
    integer, intent(in) :: ncols

    call m%free()

    allocate(m%x(nrows, ncols))
    m%x = 0.0_rp
    m%nrows = nrows
    m%ncols = ncols
    m%n = nrows*ncols

    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_map(m%x, m%x_d, m%n)
       call device_cfill(m%x_d, 0.0_rp, m%n)
    end if

  end subroutine matrix_init

  !> Deallocate a matrix.
  subroutine matrix_free(m)
    class(matrix_t), intent(inout) :: m

    if (allocated(m%x)) then
       deallocate(m%x)
    end if

    if (c_associated(m%x_d)) then
       call device_free(m%x_d)
    end if

    m%nrows = 0
    m%ncols = 0
    m%n = 0

  end subroutine matrix_free

  !> Returns the number of entries in the matrix.
  function matrix_size(m) result(s)
    class(matrix_t), intent(inout) :: m
    integer :: s
    s = m%n
  end function matrix_size

  !> Assignment \f$ m = w \f$
  subroutine matrix_assign_matrix(m, w)
    class(matrix_t), intent(inout) :: m
    type(matrix_t), intent(in) :: w

    if (allocated(m%x)) then
       call m%free()
    end if

    if (.not. allocated(m%x)) then

       m%nrows = w%nrows
       m%ncols = w%ncols
       m%n = w%n
       allocate(m%x(m%nrows, m%ncols))

       if (NEKO_BCKND_DEVICE .eq. 1) then
          call device_map(m%x, m%x_d, m%n)
       end if

    end if

    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_copy(m%x_d, w%x_d, m%n)
    else
       m%x = w%x
    end if

  end subroutine matrix_assign_matrix

  !> Assignment \f$ m = s \f$.
  subroutine matrix_assign_scalar(m, s)
    class(matrix_t), intent(inout) :: m
    real(kind=rp), intent(in) :: s

    if (.not. allocated(m%x)) then
       call neko_error('matrix not allocated')
    end if

    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_cfill(m%x_d, s, m%n)
    else
       m%x = s
    end if

  end subroutine matrix_assign_scalar


end module matrix
