! Copyright (c) 2021-2023, The Neko Authors
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
!> Implements mean_field_t.
!
module mean_field
  use neko_config, only : NEKO_BCKND_DEVICE
  use stats_quant, only : stats_quant_t
  use num_types, only : rp
  use field, only : field_t
  use field_math, only : field_cmult, field_add2s2
  implicit none
  private

  !> Computes the temporal mean of a field.
  type, public, extends(stats_quant_t) :: mean_field_t
     !> Pointer to the averaged field.
     type(field_t), pointer :: f => null()
     !> Stores the mean field.
     type(field_t) :: mf
     !> Total time across which the mean has been computed.
     real(kind=rp) :: time
   contains
     !> Constructor.
     procedure, pass(this) :: init => mean_field_init
     !> Destructor.
     procedure, pass(this) :: free => mean_field_free
     !> Updates the mean value with a new sample.
     procedure, pass(this) :: update => mean_field_update
     !> Resets the mean field.
     procedure, pass(this) :: reset => mean_field_reset
  end type mean_field_t

contains

  !> Constructor.
  !! @param f The field that will be averaged.
  !! @param field_name. Optional name for the mean field. By default the name of
  !! `f` prepended with `mean_` is used.
  subroutine mean_field_init(this, f, field_name)
    class(mean_field_t), intent(inout) :: this
    type(field_t), intent(inout), target :: f
    character(len=*), optional, intent(in) :: field_name
    character(len=80) :: name

    call this%free()

    this%f => f
    this%time = 0.0_rp

    if (present(field_name)) then
       name = field_name
    else
       write(name, '(A,A)') 'mean_', trim(f%name)
    end if

    call this%mf%init(f%dof, name)

  end subroutine mean_field_init

  !> Destructor.
  subroutine mean_field_free(this)
    class(mean_field_t), intent(inout) :: this

    if (associated(this%f)) then
       nullify(this%f)
    end if
    call this%mf%free()

  end subroutine mean_field_free

  !> Resets a the mean field and the averaging time value to zero.
  subroutine mean_field_reset(this)
    class(mean_field_t), intent(inout) :: this

    this%time = 0.0
    this%mf = 0.0_rp
  end subroutine mean_field_reset


  !> Update the mean field with a new sample.
  !! @param k Time since last sample.
  subroutine mean_field_update(this, k)
    class(mean_field_t), intent(inout) :: this
    real(kind=rp), intent(in) :: k

    call field_cmult(this%mf, this%time, size(this%mf%x))
    call field_add2s2(this%mf, this%f, k, size(this%mf%x))
    this%time = this%time + k
    call field_cmult(this%mf, 1.0_rp / this%time, size(this%mf%x))

  end subroutine mean_field_update

end module mean_field

