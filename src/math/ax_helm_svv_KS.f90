! Copyright (c) 2026, The Neko Authors
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
!     copyright notice, this list of conditions and the following disclaimer
!     in the documentation and/or other materials provided with the
!     distribution.
!
!   * Neither the name of the authors nor the names of its contributors may
!     be used to endorse or promote products derived from this software
!     without specific prior written permission.
!
! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
! AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
! IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
! ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
! LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
! CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
! SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
! INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
! CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
! ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
! POSSIBILITY OF SUCH DAMAGE.
!
!> Base type for a Kirby-Sherwin SVV Helmholtz operator.
module ax_helm_svv_KS
  use ax_helm, only : ax_helm_t
  use spectral_vanishing_viscosity, only : svv_t
  implicit none
  private

  !> Helmholtz operator carrying non-owning Kirby-Sherwin SVV object.
  type, public, abstract, extends(ax_helm_t) :: ax_helm_svv_KS_t
     !> SVV coefficients used by the operator.
     type(svv_t), pointer :: svv => null()
   contains
     procedure, pass(this) :: free => ax_helm_svv_KS_free
  end type ax_helm_svv_KS_t

contains

  !> Sever the non-owning link to the SVV object.
  !! @param this SVV Helmholtz operator.
  subroutine ax_helm_svv_KS_free(this)
    class(ax_helm_svv_KS_t), intent(inout) :: this
    nullify(this%svv)
  end subroutine ax_helm_svv_KS_free

end module ax_helm_svv_KS
