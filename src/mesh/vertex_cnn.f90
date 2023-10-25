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
!> Connectivity vertex type
module vertex_cnn
  use num_types, only : i2, i4
  use utils, only : neko_error
  use polytope_cnn, only : polytope_cnn_t
  implicit none
  private

  public :: vertex_cnn_t, vertex_cnn_ptr

  ! object information
  integer(i4), public, parameter :: NEKO_VERTEX_DIM = 0
  integer(i4), public, parameter :: NEKO_VERTEX_NFACET = 0
  integer(i4), public, parameter :: NEKO_VERTEX_NRIDGE = 0
  integer(i4), public, parameter :: NEKO_VERTEX_NPEAK = 0

  !> Vertex type for global communication
  !! @details Vertex as the only realisation of zero-dimensional polytope
  !! (monon) and contains unique global id only. Vertex has no alignment.
  type, extends(polytope_cnn_t) :: vertex_cnn_t
   contains
     procedure, pass(this) :: init => vertex_init
     procedure, pass(this) :: equal => vertex_equal
     generic :: operator(.eq.) => equal
  end type vertex_cnn_t

  !> Pointer to a vertex type
  type ::  vertex_cnn_ptr
     type(vertex_cnn_t), pointer :: obj
  end type vertex_cnn_ptr

contains

  !> @brief Initialise vertex with global id
  !! @parameter[in]   id     unique id
  subroutine vertex_init(this, id)
    class(vertex_cnn_t), intent(inout) :: this
    integer(i4), intent(in) :: id

    call this%set_dim(NEKO_VERTEX_DIM)
    call this%set_nelem(NEKO_VERTEX_NFACET, NEKO_VERTEX_NRIDGE,&
         & NEKO_VERTEX_NPEAK)
    call this%set_id(id)

    return
  end subroutine vertex_init

  !> @brief Check if two vertices are the same
  !! @return   equal
  pure function vertex_equal(this, other) result(equal)
    class(vertex_cnn_t), intent(in) :: this
    class(polytope_cnn_t), intent(in) :: other
    logical :: equal

    ! check polygon information
    equal = this%equal_poly(other)
    if (equal) then
       ! check global id
       equal = (this%id() == other%id())
    end if
    return
  end function vertex_equal

end module vertex_cnn
