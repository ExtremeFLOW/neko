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
  use neko_config, only: NEKO_BCKND_DEVICE
  use math, only: sub3, chsign, add3, cmult2, cadd2, copy
  use num_types, only: rp
  use device, only: device_map, device_free, c_ptr, C_NULL_PTR
  use device_math, only: device_copy, device_cfill, device_cmult, &
       device_sub3, device_cmult2, device_add3, device_cadd2
  use utils, only: neko_error
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
     !> Matrix-matrix addition \f$ v = m + b \f$.
     procedure, pass(m) :: matrix_add_matrix
     !> Matrix-scalar addition \f$ v = m + c \f$.
     procedure, pass(m) :: matrix_add_scalar_left
     !> Scalar-matrix addition \f$ v = c + m \f$.
     procedure, pass(m) :: matrix_add_scalar_right
     !> Matrix-matrix subtraction \f$ v = m - b \f$.
     procedure, pass(m) :: matrix_sub_matrix
     !> Matrix-scalar subtraction \f$ v = m - c \f$.
     procedure, pass(m) :: matrix_sub_scalar_left
     !> Scalar-matrix subtraction \f$ v = c - m \f$.
     procedure, pass(m) :: matrix_sub_scalar_right
     !> Matrix-scalar multiplication \f$ v = m*c \f$.
     procedure, pass(m) :: matrix_cmult_left
     !> Scalar-matrix multiplication \f$ v = c*m \f$.
     procedure, pass(m) :: matrix_cmult_right
     !> Inverse a matrix.
     procedure, pass(m) :: inverse => matrix_bcknd_inverse

     generic :: assignment(=) => matrix_assign_matrix, &
          matrix_assign_scalar
     generic :: operator(+) => matrix_add_matrix, &
          matrix_add_scalar_left, matrix_add_scalar_right
     generic :: operator(-) => matrix_sub_matrix, &
          matrix_sub_scalar_left, matrix_sub_scalar_right
     generic :: operator(*) => matrix_cmult_left, matrix_cmult_right
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

  !> Matrix-matrix addition \f$ v = m + b \f$.
  function matrix_add_matrix(m, b) result(v)
    class(matrix_t), intent(in) :: m, b
    type(matrix_t) :: v

    if (m%nrows .ne. b%nrows .or. m%ncols .ne. b%ncols) &
         call neko_error("Matrices must be the same size!")

    v%n = m%n
    v%nrows = m%nrows
    v%ncols = m%ncols
    allocate(v%x(v%nrows, v%ncols))

    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_map(v%x, v%x_d, v%n)
    end if

    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_add3(v%x_d, m%x_d, b%x_d, v%n)
    else
       call add3(v%x, m%x, b%x, v%n)
    end if

  end function matrix_add_matrix

  !> Matrix-scalar addition \f$ v = m + c \f$.
  function matrix_add_scalar_left(m, c) result(v)
    class(matrix_t), intent(in) :: m
    real(kind=rp), intent(in) :: c
    type(matrix_t) :: v

    v%n = m%n
    v%nrows = m%nrows
    v%ncols = m%ncols
    allocate(v%x(v%nrows, v%ncols))

    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_map(v%x, v%x_d, v%n)
    end if

    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_cadd2(v%x_d, m%x_d, c, v%n)
    else
       call cadd2(v%x, m%x, c, v%n)
    end if

  end function matrix_add_scalar_left

  !> Scalar-matrix addition \f$ v = c + m \f$.
  function matrix_add_scalar_right(c, m) result(v)
    real(kind=rp), intent(in) :: c
    class(matrix_t), intent(in) :: m
    type(matrix_t) :: v

    v = matrix_add_scalar_left(m, c)

  end function matrix_add_scalar_right

  !> Matrix-matrix subtraction \f$ v = m - b \f$.
  function matrix_sub_matrix(m, b) result(v)
    class(matrix_t), intent(in) :: m, b
    type(matrix_t) :: v

    if (m%nrows .ne. b%nrows .or. m%ncols .ne. b%ncols) &
         call neko_error("Matrices must be the same size!")

    v%n = m%n
    v%nrows = m%nrows
    v%ncols = m%ncols
    allocate(v%x(v%nrows, v%ncols))

    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_map(v%x, v%x_d, v%n)
    end if

    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_sub3(v%x_d, m%x_d, b%x_d, v%n)
    else
       call sub3(v%x, m%x, b%x, v%n)
    end if

  end function matrix_sub_matrix

  !> Matrix-scalar subtraction \f$ v = m - c \f$.
  function matrix_sub_scalar_left(m, c) result(v)
    class(matrix_t), intent(in) :: m
    real(kind=rp), intent(in) :: c
    type(matrix_t) :: v

    v%n = m%n
    v%nrows = m%nrows
    v%ncols = m%ncols
    allocate(v%x(v%nrows, v%ncols))

    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_map(v%x, v%x_d, v%n)
    end if

    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_cadd2(v%x_d, m%x_d, -1.0_rp*c, v%n)
    else
       call cadd2(v%x, m%x, -1.0_rp*c, m%n)
    end if

  end function matrix_sub_scalar_left

  !> Scalar-matrix subtraction \f$ v = c - m \f$.
  function matrix_sub_scalar_right(c, m) result(v)
    real(kind=rp), intent(in) :: c
    class(matrix_t), intent(in) :: m
    type(matrix_t) :: v

    v = matrix_sub_scalar_left(m, c)

    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_cmult(v%x_d, -1.0_rp, v%n)
    else
       call chsign(v%x, v%n)
    end if

  end function matrix_sub_scalar_right

  !> Matrix-scalar multiplication \f$ v = m*c \f$.
  function matrix_cmult_left(m, c) result(v)
    class(matrix_t), intent(in) :: m
    real(kind=rp), intent(in) :: c
    type(matrix_t) :: v

    v%n = m%n
    v%nrows = m%nrows
    v%ncols = m%ncols
    allocate(v%x(v%nrows, v%ncols))

    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_map(v%x, v%x_d, v%n)
    end if

    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_cmult2(v%x_d, m%x_d, c, v%n)
    else
       call cmult2(v%x, m%x, c, v%n)
    end if

  end function matrix_cmult_left

  !> Scalar-matrix multiplication \f$ v = c*m \f$.
  function matrix_cmult_right(c, m) result(v)
    real(kind=rp), intent(in) :: c
    class(matrix_t), intent(in) :: m
    type(matrix_t) :: v

    v = matrix_cmult_left(m, c)

  end function matrix_cmult_right

  subroutine matrix_bcknd_inverse(m)
    class(matrix_t), intent(inout) :: m
    if (NEKO_BCKND_DEVICE .eq. 1) then
       call neko_error("matrix_bcknd_inverse not implemented on accelarators.")
    else
       call cpu_matrix_inverse(m)
    end if
  end subroutine matrix_bcknd_inverse

  subroutine cpu_matrix_inverse(m)
    ! Gauss-Jordan matrix inversion with full pivoting
    ! Num. Rec. p. 30, 2nd Ed., Fortran
    ! m%x     is an sqaure matrix
    ! rmult is m  work array of length nrows = ncols
    class(matrix_t), intent(inout) :: m
    integer :: indr(m%nrows), indc(m%ncols), ipiv(m%ncols)
    real(kind=rp) ::  rmult(m%nrows), amx, tmp, piv, eps
    integer :: i, j, k, ir, jc

    if (.not. (m%ncols .eq. m%nrows)) then
       call neko_error("Fatal error: trying to invert m matrix that is not &
&square")
    end if

    eps = 1e-9_rp
    ipiv = 0

    do k = 1, m%nrows
       amx = 0.0_rp
       do i = 1, m%nrows                    ! Pivot search
          if (ipiv(i) .ne. 1) then
             do j = 1, m%nrows
                if (ipiv(j) .eq. 0) then
                   if (abs(m%x(i, j)) .ge. amx) then
                      amx = abs(m%x(i, j))
                      ir  = i
                      jc  = j
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
          do j = 1, m%ncols
             tmp       = m%x(ir, j)
             m%x(ir, j) = m%x(jc, j)
             m%x(jc, j) = tmp
          end do
       end if
       indr(k) = ir
       indc(k) = jc

       if (abs(m%x(jc, jc)) .lt. eps) then
          call neko_error("matrix_inverse error: small Gauss Jordan Piv")
       end if
       piv = 1.0_rp/m%x(jc, jc)
       m%x(jc, jc) = 1.0_rp
       do j = 1, m%ncols
          m%x(jc, j) = m%x(jc, j)*piv
       end do

       do j = 1, m%ncols
          tmp       = m%x(jc, j)
          m%x(jc, j) = m%x(1 , j)
          m%x(1 , j) = tmp
       end do
       do i = 2, m%nrows
          rmult(i)   = m%x(i, jc)
          m%x(i, jc)  = 0.0_rp
       end do

       do j = 1, m%ncols
          do i = 2, m%nrows
             m%x(i, j) = m%x(i, j) - rmult(i)*m%x(1, j)
          end do
       end do

       do j = 1, m%ncols
          tmp       = m%x(jc, j)
          m%x(jc, j) = m%x(1 , j)
          m%x(1 , j) = tmp
       end do
    end do

    ! Unscramble matrix
    do j= m%nrows, 1, -1
       if (indr(j) .ne. indc(j)) then
          do i = 1, m%nrows
             tmp            = m%x(i, indr(j))
             m%x(i, indr(j)) = m%x(i, indc(j))
             m%x(i, indc(j)) = tmp
          end do
       end if
    end do

    return
  end subroutine cpu_matrix_inverse

end module matrix
