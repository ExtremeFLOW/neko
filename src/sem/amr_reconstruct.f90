! Copyright (c) 2019-2025, The Neko Authors
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
!> The vector reconstruction/interpolation routines for AMR
module amr_reconstruct
  use num_types, only : i4, i8, rp, dp
  use mesh_manager_transfer, only : mesh_manager_transfer_t

  implicit none
  private

  !> Type for vector/field reconstruction
  type, public :: amr_reconstruct_t
     !> pointer mesh manager data transfer
     class(mesh_manager_transfer_t), pointer :: transfer
     !> Polynomial order or space?

   contains
     !> Initialise type
     procedure, pass(this) :: init => amr_reconstruct_init
     !> Free type
     procedure, pass(this) :: free => amr_reconstruct_free
  end type amr_reconstruct_t


contains
  !> Initialise type
  !! @param[in]  transfer     mesh manager data transfer type
  subroutine amr_reconstruct_init(this, transfer)
    class(amr_reconstruct_t), intent(inout) :: this
    class(mesh_manager_transfer_t), target, intent(in) :: transfer

    call this%free()

    ! data transfer
    this%transfer => transfer
  end subroutine amr_reconstruct_init

  !> Free type
  subroutine amr_reconstruct_free(this)
    class(amr_reconstruct_t), intent(inout) :: this

    this%transfer => NULL()
  end subroutine amr_reconstruct_free

end module amr_reconstruct
