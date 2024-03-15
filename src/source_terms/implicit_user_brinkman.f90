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
!> Implements the `implicit_user_brinkman` type.
! Consider a scalar field \chi
! and a 'Brinkman' forcing term f = \chi * u
! 
! applying the forcing in user_f results in a explicit formulation
! f = \chi * u^(n)
! To improve stability, we treat this term implicitly
! f = \chi * u^(n+1)
!
! Applications extend to
!  - Immersed Boundary Methods
!  - Sponge regions
!  - tripping/freestream turbulence etc

module implicit_user_brinkman

	use field, only: field_t 
	use field_registry, only: neko_field_registry
	use coefs, only: coef_t
   use user_intf, only: user_t
   use num_types, only: rp
	private

	public :: implicit_brinkman

	type, public :: implicit_user_brinkman_t
		type(field_t), pointer :: chi
		procedure(implicit_brinkman), nopass, pointer :: compute_user_brinkman => null()		
	contains
		procedure, pass(this) :: init => implicit_user_brinkman_init
		procedure, pass(this) :: compute => implicit_user_brinkman_compute
   end type implicit_user_brinkman_t

  !> Abstract interface for user implicit Brinkman forcing
  abstract interface
     subroutine implicit_brinkman(chi, t)
       import field_t
       import rp
       type(field_t), intent(inout) :: chi
    	 real(kind=rp), intent(in) :: t
     end subroutine implicit_brinkman
  end interface

	contains

	subroutine implicit_user_brinkman_init(this, coef, user_implicit_brinkman)
		class(implicit_user_brinkman_t), intent(inout) :: this
		type(coef_t), intent(in) :: coef
		procedure(implicit_brinkman) :: user_implicit_brinkman

		call neko_field_registry%add_field(coef%dof, "chi")
		this%chi => neko_field_registry%get_field("chi")

		this%compute_user_brinkman => user_implicit_brinkman

	end subroutine implicit_user_brinkman_init
	
	subroutine implicit_user_brinkman_compute(this, t)
		class(implicit_user_brinkman_t), intent(inout) :: this
		real(kind=rp) :: t

		call this%compute_user_brinkman(this%chi, t)

	end subroutine implicit_user_brinkman_compute

end module implicit_user_brinkman
