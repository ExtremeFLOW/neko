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
!> Dummy interpolation operators for nonconforming vertices
module subset_interpolation_vertex
  use num_types, only : i4, rp
  use subset_interpolation, only : subset_interpolation_t
  implicit none
  private

  public :: subset_intp_vertex_init, subset_interpolation_vertex_dmy_t

  !> number of various interpolation operators
  integer(i4), public, parameter :: NEKO_INTP_VERTEX_NOPERATION = 0

  !> Dummy interpolation operators
  type, extends(subset_interpolation_t) :: subset_interpolation_vertex_dmy_t
   contains
     !> Is transformation a dummy one
     procedure, pass(this) :: ifdmy => ifdummy_vertex_dmy
     !> Patent-child interpolation
     procedure, nopass :: intp  => transform_vertex_dmy
     !> Transposed interpolation
     procedure, nopass :: intpT => transform_vertex_dmy
     !> Initialise interpolation data
     procedure, pass(this) :: set_jmat  => vertex_set_jmat_dmy
     !> Free interpolation data
     procedure, pass(this) :: free_jmat => vertex_free_jmat_dmy
  end type subset_interpolation_vertex_dmy_t

contains

  !> @brief Allocate a single interpolation operator
  !! @parameter[in]      hng   hanging information
  !! @parameter[inout]   trns  interpolation operator
  subroutine subset_intp_vertex_init(hng, trns)
    integer(i4), intent(in) :: hng
    class(subset_interpolation_t), allocatable, intent(inout) :: trns

    if (allocated(trns)) then
       call trns%free_jmat()
       deallocate(trns)
    end if

    ! there is only a dummy operation for vertex
    allocate(subset_interpolation_vertex_dmy_t :: trns)
    call trns%set_hng(hng)

  end subroutine subset_intp_vertex_init

  !> Function returning dummy operation flag
  !! @return   ifdmy
  pure function ifdummy_vertex_dmy(this) result(ifdmy)
    class(subset_interpolation_vertex_dmy_t), intent(in) :: this
    logical :: ifdmy
    ifdmy = .true.
  end function ifdummy_vertex_dmy

  !> Dummy interpolation operator, do nothing
  !! @parameter[inout]   vec      data vector
  !! @parameter[in]      n1, n2   dimensions
  pure subroutine transform_vertex_dmy(vec, n1, n2)
    integer(i4), intent(in) :: n1, n2
    real(rp), dimension(n1, n2), intent(inout) :: vec
  end subroutine transform_vertex_dmy

  !> Dummy initialisation of the interpolation data
  !! @notice  This routine can be called after space_t gets initialised.
  !! @parameter[in]      lr, ls   array sizes for r and s dimensions
  subroutine vertex_set_jmat_dmy(this, lr, ls)
    class(subset_interpolation_vertex_dmy_t), intent(inout) :: this
    integer(i4), intent(in) :: lr, ls
  end subroutine vertex_set_jmat_dmy

  !> Free the dummy interpolation data
  subroutine vertex_free_jmat_dmy(this)
    class(subset_interpolation_vertex_dmy_t), intent(inout) :: this
  end subroutine vertex_free_jmat_dmy

end module subset_interpolation_vertex
