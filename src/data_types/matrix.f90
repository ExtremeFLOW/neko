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
!> Defines a matrix
module matrix
  use array, only : array_t
  use num_types, only : rp, xp
  use neko_config, only : NEKO_BCKND_DEVICE
  use device, only : HOST_TO_DEVICE, DEVICE_TO_HOST
  use utils, only : neko_warning, neko_error
  implicit none
  private

  type, public, extends(array_t) :: matrix_t
     !> Matrix shaped pointer.
     real(kind=rp), pointer, dimension(:,:) :: x => null()
     !> Number of matrix rows.
     integer, private :: nrows = 0
     !> Number of matrix columns.
     integer, private :: ncols = 0

   contains
     !> Initialise a matrix of size `nrows*ncols`.
     procedure, pass(this) :: init => matrix_init
     !> Deallocate a matrix.
     procedure, pass(this) :: free => matrix_free
     !> Returns the number of rows in the matrix.
     procedure, pass(this) :: get_nrows => matrix_nrows
     !> Returns the number of columns in the matrix.
     procedure, pass(this) :: get_ncols => matrix_ncols

     !> Inverse a matrix.
     procedure, pass(this) :: inverse => matrix_inverse
     procedure, pass(this) :: inverse_on_host => cpu_matrix_inverse

  end type matrix_t

  type, public :: matrix_ptr_t
     type(matrix_t), pointer :: ptr
  end type matrix_ptr_t

contains

  !> Initialise a matrix of size `nrows*ncols`.
  !! @param nrows Number of rows.
  !! @param ncols Number of columns.
  !! @param name Optional name of the object.
  subroutine matrix_init(this, nrows, ncols, name)
    class(matrix_t), intent(inout), target :: this
    integer, intent(in) :: nrows, ncols
    character(len=*), intent(in), optional :: name

    call this%init_base(nrows*ncols, name)
    this%nrows = nrows
    this%ncols = ncols
    this%x(1:nrows, 1:ncols) => this%data(:)

  end subroutine matrix_init

  !> Deallocate a matrix.
  subroutine matrix_free(this)
    class(matrix_t), intent(inout) :: this

    call this%free_base()
    this%nrows = 0
    this%ncols = 0
    nullify(this%x)

  end subroutine matrix_free

  !> Returns the number of rows in the matrix.
  pure function matrix_nrows(this) result(nr)
    class(matrix_t), intent(in) :: this
    integer :: nr
    nr = this%nrows
  end function matrix_nrows

  !> Returns the number of columns in the matrix.
  pure function matrix_ncols(this) result(nc)
    class(matrix_t), intent(in) :: this
    integer :: nc
    nc = this%ncols
  end function matrix_ncols

  !> Compute the inverse of the matrix.
  subroutine matrix_inverse(this, bcknd)
    class(matrix_t), intent(inout) :: this
    integer, optional :: bcknd

    if (NEKO_BCKND_DEVICE .eq. 1 .and. bcknd .eq. NEKO_BCKND_DEVICE) then
       call neko_error("matrix_bcknd_inverse not implemented on accelarators.")
    else
       call this%inverse_on_host()
    end if

  end subroutine matrix_inverse

  subroutine cpu_matrix_inverse(this)
    ! Gauss-Jordan matrix inversion with full pivoting
    ! Num. Rec. p. 30, 2nd Ed., Fortran
    ! this%x     is an square matrix
    ! rmult is this  work array of length nrows = ncols
    class(matrix_t), intent(inout) :: this
    integer :: indr(this%nrows), indc(this%ncols), ipiv(this%ncols)
    real(kind=xp) :: rmult(this%nrows), amx, tmp, piv, eps
    integer :: i, j, k, ir, jc

    if (.not. (this%ncols .eq. this%nrows)) then
       call neko_error("Fatal error: trying to invert this matrix that is not &
       &square")
    end if

    eps = 1e-9_rp
    ipiv = 0

    do k = 1, this%nrows
       amx = 0.0_rp
       do i = 1, this%nrows ! Pivot search
          if (ipiv(i) .ne. 1) then
             do j = 1, this%nrows
                if (ipiv(j) .eq. 0) then
                   if (abs(this%x(i, j)) .ge. amx) then
                      amx = abs(this%x(i, j))
                      ir = i
                      jc = j
                   end if
                else if (ipiv(j) .gt. 1) then
                   return
                end if
             end do
          end if
       end do
       ipiv(jc) = ipiv(jc) + 1

       !  Swap rows
       if (ir .ne. jc) then
          do j = 1, this%ncols
             tmp = this%x(ir, j)
             this%x(ir, j) = this%x(jc, j)
             this%x(jc, j) = tmp
          end do
       end if
       indr(k) = ir
       indc(k) = jc

       if (abs(this%x(jc, jc)) .lt. eps) then
          call neko_error("matrix_inverse error: small Gauss Jordan Piv")
       end if
       piv = 1.0_xp/this%x(jc, jc)
       this%x(jc, jc) = 1.0_xp
       do j = 1, this%ncols
          this%x(jc, j) = this%x(jc, j)*piv
       end do

       do j = 1, this%ncols
          tmp = this%x(jc, j)
          this%x(jc, j) = this%x(1 , j)
          this%x(1 , j) = tmp
       end do
       do i = 2, this%nrows
          rmult(i) = this%x(i, jc)
          this%x(i, jc) = 0.0_rp
       end do

       do j = 1, this%ncols
          do i = 2, this%nrows
             this%x(i, j) = this%x(i, j) - rmult(i)*this%x(1, j)
          end do
       end do

       do j = 1, this%ncols
          tmp = this%x(jc, j)
          this%x(jc, j) = this%x(1 , j)
          this%x(1 , j) = tmp
       end do
    end do

    ! Unscramble matrix
    do j = this%nrows, 1, -1
       if (indr(j) .ne. indc(j)) then
          do i = 1, this%nrows
             tmp = this%x(i, indr(j))
             this%x(i, indr(j)) = this%x(i, indc(j))
             this%x(i, indc(j)) = tmp
          end do
       end if
    end do

  end subroutine cpu_matrix_inverse

end module matrix
