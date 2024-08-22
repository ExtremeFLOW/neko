! Copyright (c) 2020-2023, The Neko Authors
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
!> Defines an output for a adjoint
module adjoint_output
  use num_types, only : rp
  use adjoint_scheme, only : adjoint_scheme_t
  use scalar_scheme, only : scalar_scheme_t
  use field_list, only : field_list_t
  use neko_config, only : NEKO_BCKND_DEVICE
  use device
  use output, only : output_t
  implicit none
  private

  !> adjoint output
  type, public, extends(output_t) :: adjoint_output_t
     type(field_list_t) :: adjoint
   contains
     procedure, pass(this) :: sample => adjoint_output_sample
  end type adjoint_output_t

  interface adjoint_output_t
     module procedure adjoint_output_init
  end interface adjoint_output_t

contains

  function adjoint_output_init(precision, adjoint, scalar, name, path) result(this)
    integer, intent(inout) :: precision
    class(adjoint_scheme_t), intent(in), target :: adjoint
    class(scalar_scheme_t), intent(in), optional, target :: scalar
    character(len=*), intent(in), optional :: name
    character(len=*), intent(in), optional :: path
    type(adjoint_output_t) :: this
    character(len=1024) :: fname

    if (present(name) .and. present(path)) then
       fname = trim(path) // trim(name) // '.fld'
    else if (present(name)) then
       fname = trim(name) // '.fld'
    else if (present(path)) then
       fname = trim(path) // 'adjoint.fld'
    else
       fname = 'adjoint.fld'
    end if

    call this%init_base(fname, precision)

    if (present(scalar)) then
       call this%adjoint%init(5)
    else
       call this%adjoint%init(4)
    end if

    call this%adjoint%assign(1, adjoint%p_adj)
    call this%adjoint%assign(2, adjoint%u_adj)
    call this%adjoint%assign(3, adjoint%v_adj)
    call this%adjoint%assign(4, adjoint%w_adj)

    if (present(scalar)) then
       call this%adjoint%assign(5, scalar%s)
    end if

  end function adjoint_output_init

  !> Sample a adjoint solution at time @a t
  subroutine adjoint_output_sample(this, t)
    class(adjoint_output_t), intent(inout) :: this
    real(kind=rp), intent(in) :: t
    integer :: i

    if (NEKO_BCKND_DEVICE .eq. 1) then

       associate(fields => this%adjoint%items)
         do i = 1, size(fields)
            call device_memcpy(fields(i)%ptr%x, fields(i)%ptr%x_d, &
                 fields(i)%ptr%dof%size(), DEVICE_TO_HOST, &
                 sync=(i .eq. size(fields))) ! Sync on the last field
         end do
       end associate

    end if

    call this%file_%write(this%adjoint, t)

  end subroutine adjoint_output_sample

end module adjoint_output
