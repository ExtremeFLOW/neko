! Copyright (c) 2022, The Neko Authors
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
!> Defines an output for a mean flow field
module fluid_stats_output
  use fluid_stats, only : fluid_stats_t
  use neko_config, only : NEKO_BCKND_DEVICE
  use num_types, only : rp
  use map_1d
  use device
  use output, only : output_t
  use matrix
  implicit none
  private

  type, public, extends(output_t) :: fluid_stats_output_t
     type(fluid_stats_t), pointer :: stats
     type(map_1d_t) :: map_1d
     real(kind=rp) :: T_begin
   contains
     procedure, pass(this) :: sample => fluid_stats_output_sample
  end type fluid_stats_output_t

  interface fluid_stats_output_t
     module procedure fluid_stats_output_init
  end interface fluid_stats_output_t

contains

  function fluid_stats_output_init(stats, T_begin, hom_dir, name, path) result(this)
    type(fluid_stats_t), intent(inout), target :: stats
    real(kind=rp), intent(in) :: T_begin
    character(len=*), intent(in) :: hom_dir
    character(len=*), intent(in), optional :: name
    character(len=*), intent(in), optional :: path
    type(fluid_stats_output_t) :: this
    character(len=1024) :: fname
    if (trim(hom_dir) .eq. 'none') then
       if (present(name) .and. present(path)) then
          fname = trim(path) // trim(name) // '.fld'
       else if (present(name)) then
          fname = trim(name) // '.fld'
       else if (present(path)) then
          fname = trim(path) // 'fluid_stats.fld'
       else
          fname = 'fluid_stats.fld'
       end if
    else
       if (present(name) .and. present(path)) then
          fname = trim(path) // trim(name) // '.csv'
       else if (present(name)) then
          fname = trim(name) // '.csv'
       else if (present(path)) then
          fname = trim(path) // 'fluid_stats.csv'
       else
          fname = 'fluid_stats.csv'
       end if
       call this%map_1d%init_char(stats%coef, hom_dir, 1e-7_rp)
    end if

    call this%init_base(fname)
    this%stats => stats
    this%T_begin = T_begin
  end function fluid_stats_output_init

  !> Sample a mean flow field at time @a t
  subroutine fluid_stats_output_sample(this, t)
    class(fluid_stats_output_t), intent(inout) :: this
    real(kind=rp), intent(in) :: t
    integer :: i
    type(matrix_t) :: avg_output_1d
    associate (out_fields => this%stats%stat_fields%items)
      if (t .ge. this%T_begin) then
         call this%stats%make_strong_grad()
         if ( NEKO_BCKND_DEVICE .eq. 1) then
            do i = 1, size(out_fields)
               call device_memcpy(out_fields(i)%ptr%x, out_fields(i)%ptr%x_d,&
                  out_fields(i)%ptr%dof%size(), DEVICE_TO_HOST, &
                  sync=(i .eq. size(out_fields))) ! Sync on last field
            end do
         end if
         if (allocated(this%map_1d%pt_lvl)) then
            call this%map_1d%average_planes(avg_output_1d, this%stats%stat_fields)
            call this%file_%write(avg_output_1d, t)
         else
            call this%file_%write(this%stats%stat_fields, t)
         end if
         call this%stats%reset()
      end if
    end associate
  end subroutine fluid_stats_output_sample

end module fluid_stats_output


