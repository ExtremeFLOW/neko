! Copyright (c) 2021, The Neko Authors
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
!> Defines a tetrahedral mesh
!! @details Mesh dervied from an existing hexahedral mesh via bisection
module tet_mesh
  use mesh
  use tet
  use point
  use utils
  implicit none
  private

  integer, public, parameter :: TET_MSH_OTPV = 1, TET_MSH_FVTC = 2, &
       TET_MSH_SVTC = 3

  type, public :: tet_mesh_t
     type(tet_t), allocatable :: el(:) !< Tetrahedron elements
     type(mesh_t), pointer :: msh      !< Hexahedron mesh
     integer :: nelv                   !< Number of Tetrahedrons
   contains
     procedure, pass(this) :: init => tet_mesh_init
     procedure, pass(this) :: free => tet_mesh_free
  end type tet_mesh_t

contains

  !> Initialise a tetrahedral mesh based on a hexahedral mesh @a msh
  subroutine tet_mesh_init(this, msh, mthd)
    class(tet_mesh_t), intent(inout) :: this
    type(mesh_t), intent(in), target :: msh
    integer, intent(in), optional :: mthd
    integer :: bsct_mthd

    call this%free()

    this%msh => msh

    if (present(mthd)) then
       bsct_mthd = mthd
    else
       bsct_mthd = TET_MSH_SVTC
    end if

    if (bsct_mthd .eq. TET_MSH_OTPV) then
       this%nelv = msh%nelv * 8
       allocate(this%el(this%nelv))
       call tet_mesh_bisect_otpv(this)
    else if (bsct_mthd .eq. TET_MSH_FVTC) then
       this%nelv = msh%nelv * 5
       allocate(this%el(this%nelv))
       call tet_mesh_bisect_fvtc(this)
    else if (bsct_mthd .eq. TET_MSH_SVTC) then
       this%nelv = msh%nelv * 6
       allocate(this%el(this%nelv))
       call tet_mesh_bisect_svtc(this)
    else
       call neko_error('Invalid bisection strategy')
    end if

  end subroutine tet_mesh_init

  !> Deallocate a tetrahedral mesh
  subroutine tet_mesh_free(this)
    class(tet_mesh_t), intent(inout) :: this

    if (allocated(this%el)) then
       deallocate(this%el)
    end if

    nullify(this%msh)

  end subroutine tet_mesh_free

  ! Bisect hexahedral mesh into a tetrahedral mesh using
  ! one tet per vertex as described in,
  ! P. D. Bello-Maldonado and P. F. Fischer
  ! SIAM J. Sci. Comput., 41(5), S2–S18, 2019.
  subroutine tet_mesh_bisect_otpv(tet_msh)
    type(tet_mesh_t), intent(inout) :: tet_msh
    integer :: i, j
    type(point_t), pointer :: p1, p2, p3, p4

    j = 0
    do i = 1, tet_msh%msh%nelv

       j = j + 1
       p1 => tet_msh%msh%elements(i)%e%pts(1)%p
       p2 => tet_msh%msh%elements(i)%e%pts(2)%p
       p3 => tet_msh%msh%elements(i)%e%pts(5)%p
       p4 => tet_msh%msh%elements(i)%e%pts(4)%p
       call tet_msh%el(j)%init(j, p1, p2, p3, p4)

       j = j + 1
       p1 => tet_msh%msh%elements(i)%e%pts(1)%p
       p2 => tet_msh%msh%elements(i)%e%pts(2)%p
       p3 => tet_msh%msh%elements(i)%e%pts(3)%p
       p4 => tet_msh%msh%elements(i)%e%pts(6)%p
       call tet_msh%el(j)%init(j, p1, p2, p3, p4)

       j = j + 1
       p1 => tet_msh%msh%elements(i)%e%pts(4)%p
       p2 => tet_msh%msh%elements(i)%e%pts(3)%p
       p3 => tet_msh%msh%elements(i)%e%pts(2)%p
       p4 => tet_msh%msh%elements(i)%e%pts(7)%p
       call tet_msh%el(j)%init(j, p1, p2, p3, p4)

       j = j + 1
       p1 => tet_msh%msh%elements(i)%e%pts(4)%p
       p2 => tet_msh%msh%elements(i)%e%pts(3)%p
       p3 => tet_msh%msh%elements(i)%e%pts(8)%p
       p4 => tet_msh%msh%elements(i)%e%pts(1)%p
       call tet_msh%el(j)%init(j, p1, p2, p3, p4)

       j = j + 1
       p1 => tet_msh%msh%elements(i)%e%pts(1)%p
       p2 => tet_msh%msh%elements(i)%e%pts(5)%p
       p3 => tet_msh%msh%elements(i)%e%pts(6)%p
       p4 => tet_msh%msh%elements(i)%e%pts(8)%p
       call tet_msh%el(j)%init(j, p1, p2, p3, p4)

       j = j + 1
       p1 => tet_msh%msh%elements(i)%e%pts(2)%p
       p2 => tet_msh%msh%elements(i)%e%pts(5)%p
       p3 => tet_msh%msh%elements(i)%e%pts(6)%p
       p4 => tet_msh%msh%elements(i)%e%pts(7)%p
       call tet_msh%el(j)%init(j, p1, p2, p3, p4)

       j = j + 1
       p1 => tet_msh%msh%elements(i)%e%pts(3)%p
       p2 => tet_msh%msh%elements(i)%e%pts(8)%p
       p3 => tet_msh%msh%elements(i)%e%pts(7)%p
       p4 => tet_msh%msh%elements(i)%e%pts(6)%p
       call tet_msh%el(j)%init(j, p1, p2, p3, p4)

       j = j + 1
       p1 => tet_msh%msh%elements(i)%e%pts(4)%p
       p2 => tet_msh%msh%elements(i)%e%pts(8)%p
       p3 => tet_msh%msh%elements(i)%e%pts(7)%p
       p4 => tet_msh%msh%elements(i)%e%pts(5)%p
       call tet_msh%el(j)%init(j, p1, p2, p3, p4)

    end do

  end subroutine tet_mesh_bisect_otpv

  !> Bisect each hexahedron into five tetrahedrons
  subroutine tet_mesh_bisect_fvtc(tet_msh)
    type(tet_mesh_t), intent(inout) :: tet_msh
    integer :: i, j
    type(point_t), pointer :: p1, p2, p3, p4

    j = 0
    do i = 1, tet_msh%msh%nelv

       j = j + 1
       p1 => tet_msh%msh%elements(i)%e%pts(1)%p
       p2 => tet_msh%msh%elements(i)%e%pts(2)%p
       p3 => tet_msh%msh%elements(i)%e%pts(3)%p
       p4 => tet_msh%msh%elements(i)%e%pts(6)%p
       call tet_msh%el(j)%init(j, p1, p2, p3, p4)

       j = j + 1
       p1 => tet_msh%msh%elements(i)%e%pts(1)%p
       p2 => tet_msh%msh%elements(i)%e%pts(4)%p
       p3 => tet_msh%msh%elements(i)%e%pts(3)%p
       p4 => tet_msh%msh%elements(i)%e%pts(8)%p
       call tet_msh%el(j)%init(j, p1, p2, p3, p4)

       j = j + 1
       p1 => tet_msh%msh%elements(i)%e%pts(3)%p
       p2 => tet_msh%msh%elements(i)%e%pts(8)%p
       p3 => tet_msh%msh%elements(i)%e%pts(7)%p
       p4 => tet_msh%msh%elements(i)%e%pts(6)%p
       call tet_msh%el(j)%init(j, p1, p2, p3, p4)

       j = j + 1
       p1 => tet_msh%msh%elements(i)%e%pts(1)%p
       p2 => tet_msh%msh%elements(i)%e%pts(5)%p
       p3 => tet_msh%msh%elements(i)%e%pts(8)%p
       p4 => tet_msh%msh%elements(i)%e%pts(6)%p
       call tet_msh%el(j)%init(j, p1, p2, p3, p4)

       j = j + 1
       p1 => tet_msh%msh%elements(i)%e%pts(8)%p
       p2 => tet_msh%msh%elements(i)%e%pts(3)%p
       p3 => tet_msh%msh%elements(i)%e%pts(1)%p
       p4 => tet_msh%msh%elements(i)%e%pts(6)%p
       call tet_msh%el(j)%init(j, p1, p2, p3, p4)
    end do
  end subroutine tet_mesh_bisect_fvtc

  !> Bisect each hexahedron into six tetrahedrons
  subroutine tet_mesh_bisect_svtc(tet_msh)
    type(tet_mesh_t), intent(inout) :: tet_msh
    integer :: i, j
    type(point_t), pointer :: p1, p2, p3, p4

    j = 0
    do i = 1, tet_msh%msh%nelv

       j = j + 1
       p1 => tet_msh%msh%elements(i)%e%pts(5)%p
       p2 => tet_msh%msh%elements(i)%e%pts(8)%p
       p3 => tet_msh%msh%elements(i)%e%pts(1)%p
       p4 => tet_msh%msh%elements(i)%e%pts(7)%p
       call tet_msh%el(j)%init(j, p1, p2, p3, p4)

       j = j + 1
       p1 => tet_msh%msh%elements(i)%e%pts(5)%p
       p2 => tet_msh%msh%elements(i)%e%pts(1)%p
       p3 => tet_msh%msh%elements(i)%e%pts(7)%p
       p4 => tet_msh%msh%elements(i)%e%pts(6)%p
       call tet_msh%el(j)%init(j, p1, p2, p3, p4)

       j = j + 1
       p1 => tet_msh%msh%elements(i)%e%pts(8)%p
       p2 => tet_msh%msh%elements(i)%e%pts(4)%p
       p3 => tet_msh%msh%elements(i)%e%pts(7)%p
       p4 => tet_msh%msh%elements(i)%e%pts(1)%p
       call tet_msh%el(j)%init(j, p1, p2, p3, p4)

       j = j + 1
       p1 => tet_msh%msh%elements(i)%e%pts(1)%p
       p2 => tet_msh%msh%elements(i)%e%pts(7)%p
       p3 => tet_msh%msh%elements(i)%e%pts(2)%p
       p4 => tet_msh%msh%elements(i)%e%pts(6)%p
       call tet_msh%el(j)%init(j, p1, p2, p3, p4)

       j = j + 1
       p1 => tet_msh%msh%elements(i)%e%pts(1)%p
       p2 => tet_msh%msh%elements(i)%e%pts(2)%p
       p3 => tet_msh%msh%elements(i)%e%pts(7)%p
       p4 => tet_msh%msh%elements(i)%e%pts(3)%p
       call tet_msh%el(j)%init(j, p1, p2, p3, p4)

       j = j + 1
       p1 => tet_msh%msh%elements(i)%e%pts(1)%p
       p2 => tet_msh%msh%elements(i)%e%pts(7)%p
       p3 => tet_msh%msh%elements(i)%e%pts(3)%p
       p4 => tet_msh%msh%elements(i)%e%pts(4)%p
       call tet_msh%el(j)%init(j, p1, p2, p3, p4)

    end do


  end subroutine tet_mesh_bisect_svtc

end module tet_mesh
