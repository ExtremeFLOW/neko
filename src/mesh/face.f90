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
!> Connectivity face abstract type
module face
  use num_types, only : i2, i4
  use utils, only : neko_error
  use polytope, only : polytope_t
  use vertex, only : vertex_t, vertex_ptr
  use edge, only : edge_t, edge_ptr, edge_aligned_t, edge_aligned_ptr
  implicit none
  private

  public :: face_t

  ! object information
  integer(i4), parameter :: NEKO_FACE_DIM = 2

  !> Base type for an abstract two-dimensional polytope
  !! @details There are multiple possible realisation, so it is just
  !! an abstract type
  type, extends(polytope_t), abstract :: face_t
     !> Facets are aligned
     type(edge_aligned_t), dimension(:), allocatable :: facet
     !> Ridge pointers
     type(vertex_ptr), dimension(:), allocatable :: ridge
   contains
     !> Initialise face dimension
     procedure, pass(this) :: init_dim => face_init_dim
     !> Free face data
     procedure, pass(this) :: free => face_free
     !> Is edged self-periodic
     procedure, pass(this) :: selfp => face_self_periodic
     !> Get pointers to facets
     procedure, pass(this) :: fct => face_facet
     !> Get pointers to ridges
     procedure, pass(this) :: rdg => face_ridge
  end type face_t

contains

  !> @brief Initialise face dimension
  subroutine face_init_dim(this)
    class(face_t), intent(inout) :: this

    call this%free()

    call this%set_dim(NEKO_FACE_DIM)

    return
  end subroutine face_init_dim

  !> @brief Free face data
  subroutine face_free(this)
    class(face_t), intent(inout) :: this
    !local variables
    integer(i4) :: il

    call this%set_dim(-1)
    if (allocated(this%facet)) then
       do il = 1, this%nfacet
          this%facet(il)%edge%obj => null()
       end do
       deallocate(this%facet)
    end if
    if (allocated(this%ridge)) then
       do il = 1, this%nridge
          this%ridge(il)%obj => null()
       end do
       deallocate(this%ridge)
    end if

    return
  end subroutine face_free

  !> @brief Check if face is self-periodic
  !! @return   selfp
  function face_self_periodic(this) result(selfp)
    class(face_t), intent(in) :: this
    logical :: selfp
    integer(i4) :: il, jl, itmp

    ! count self periodic edges
    itmp = 0
    do il = 1, this%nfacet - 1
       do jl = il + 1, this%nfacet - 1
          selfp = this%facet(il)%edge%obj.eq.this%facet(jl)%edge%obj
          if (selfp) itmp = itmp + 1
       end do
    end do
    ! Do I have to do the same with ridges?
    if (itmp == 0) then
       selfp = .false.
    else
       selfp = .true.
    end if

    return
  end function face_self_periodic

  !> @brief Return pointers to face facets
  !! @parameter[out]  facet   facet pointers array
  subroutine face_facet(this, facet)
    class(face_t), intent(in) :: this
    type(edge_ptr), dimension(:), allocatable, intent(out) :: facet
    integer(i4) :: il

    allocate(facet(this%nfacet))
    do il = 1, this%nfacet
       facet(il)%obj => this%facet(il)%edge%obj
    end do

    return
  end subroutine face_facet

  !> @brief Return pointers to face ridges
  !! @parameter[out]  ridge   ridge pointers array
  subroutine face_ridge(this, ridge)
    class(face_t), intent(in) :: this
    type(vertex_ptr), dimension(:), allocatable, intent(out) :: ridge
    integer(i4) :: il

    allocate(ridge(this%nridge))
    do il = 1, this%nridge
       ridge(il)%obj => this%ridge(il)%obj
    end do

    return
  end subroutine face_ridge

end module face
