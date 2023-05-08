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
!> Compound scheme for the advection and diffusion operators in a transport
!! equation. 
module advection_diffusion_time_scheme
  use neko_config
  use num_types, only : rp
  use bdf_time_scheme, only: bdf_time_scheme_t
  use ext_time_scheme, only: ext_time_scheme_t
  implicit none
  private

  type, public :: advection_diffusion_time_scheme_t
     type(ext_time_scheme_t) :: ext
     type(bdf_time_scheme_t) :: bdf
     
   contains
     procedure, pass(this) :: init => advection_diffusion_scheme_init
     procedure, pass(this) :: free => advection_diffusion_scheme_free
  end type advection_diffusion_time_scheme_t

  contains

  !> Contructor
  !! @param torder Desired order of the scheme: 1, 2, 3.
  subroutine advection_diffusion_scheme_init(this, torder)
     class(advection_diffusion_time_scheme_t) :: this
     integer :: torder 
     if (torder .eq. 1) then 
        call this%ext%init(torder)
     else 
        call this%ext%init(3)
     end if
     call this%bdf%init(3)
  end subroutine advection_diffusion_scheme_init

  !> Destructor
  subroutine advection_diffusion_scheme_free(this)
     class(advection_diffusion_time_scheme_t) :: this
     call this%ext%free()
     call this%bdf%free()
  end subroutine advection_diffusion_scheme_free
  

end module advection_diffusion_time_scheme