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
!> Defines a mean field
!
module mean_field
  use neko_config
  use stats_quant
  use device_math
  use num_types
  use field
  use math
  use device_math
  use field_registry
  implicit none
  private
  
  type, public, extends(stats_quant_t) ::  mean_field_t
     type(field_t), pointer :: f => null()
     type(field_t) :: mf
     real(kind=rp) :: time
   contains
     procedure, pass(this) :: init => mean_field_init
     procedure, pass(this) :: free => mean_field_free
     procedure, pass(this) :: update => mean_field_update
     procedure, pass(this) :: reset => mean_field_reset
  end type mean_field_t

contains

  !> Initialize a mean field for a field @a f
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
       write(name, '(A,A)') 'mean_',trim(f%name)
    end if

    call this%mf%init(f%dof, name)

  end subroutine mean_field_init

  !> Deallocates a mean field
  subroutine mean_field_free(this)
    class(mean_field_t), intent(inout) :: this

    if (associated(this%f)) then
       nullify(this%f)
    end if
    call this%mf%free()

  end subroutine mean_field_free

  !> Resets a mean field
  subroutine mean_field_reset(this)
    class(mean_field_t), intent(inout) :: this
    
    this%time = 0.0
    this%mf = 0.0_rp  
  end subroutine mean_field_reset


  !> Update a mean field
  subroutine mean_field_update(this, k)
    class(mean_field_t), intent(inout) :: this
    real(kind=rp), intent(in) :: k !< Time since last sample

    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_cmult(this%mf%x_d, this%time, size(this%mf%x))
       call device_add2s2(this%mf%x_d, this%f%x_d, k, size(this%mf%x))
       this%time = this%time + k
       call device_cmult(this%mf%x_d, 1.0_rp / this%time, size(this%mf%x))
    else
       this%mf%x = this%mf%x * this%time
       call add2s2(this%mf%x, this%f%x, k, this%mf%dof%size())
       this%time = this%time + k
       this%mf%x = this%mf%x / this%time
    end if
       
  end subroutine mean_field_update
  
end module mean_field

