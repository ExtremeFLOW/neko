! Copyright (c) 2019-2026, The Neko Authors
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
!> Implementation for the nonconforming faces/edges during communication
module gs_interp_cpu
  use num_types, only : i4, i8, rp
  use mesh_conn, only : mesh_conn_t
  use amr_interpolate, only : amr_interpolate_t, amr_nchildren
  use field, only : field_t
  use gs_interp, only : gs_interp_t
  use amr_reconstruct, only : amr_reconstruct_t

  implicit none
  private

  !> Type for face/edge
  type, public, extends(gs_interp_t) :: gs_interp_cpu_t

   contains
     !> Initialise type
     procedure, pass(this) :: init => gs_interp_cpu_init
     !> Free type
     procedure, pass(this) :: free => gs_interp_cpu_free
     !> AMR restart
     procedure, pass(this) :: amr_restart => gs_interp_cpu_amr_restart
  end type gs_interp_cpu_t

contains
  !> Initialise gs interpolation type
  !! @param[in]  lx    polynomial order + 1
  !! @param[in]  conn  mesh connectivity
  subroutine gs_interp_cpu_init(this, lx, conn)
    class(gs_interp_cpu_t), intent(inout) :: this
    integer, intent(in) :: lx
    type(mesh_conn_t), target, intent(in) :: conn

    
    write(*, *) 'TESTinterpINIT'
     

    call this%free()

    call this%init_base(lx, conn)

  end subroutine gs_interp_cpu_init

  !> Free gs interpolation type
  subroutine gs_interp_cpu_free(this)
    class(gs_interp_cpu_t), intent(inout) :: this

    call this%free_base()

  end subroutine gs_interp_cpu_free

  !> AMR restart
  !! @param[inout]  reconstruct   data reconstruction type
  !! @param[in]     counter       restart counter
  !! @param[in]     tstep         time step
  subroutine gs_interp_cpu_amr_restart(this, reconstruct, counter, tstep)
    class(gs_interp_cpu_t), intent(inout) :: this
    type(amr_reconstruct_t), intent(inout) :: reconstruct
    integer, intent(in) :: counter, tstep

    ! Was this component already restarted?
    if (this%counter .eq. counter) return

    this%counter = counter

    
    write(*, *) 'TESTinterpRESTART'
     

    call this%amr_restart_base()

  end subroutine gs_interp_cpu_amr_restart

end module gs_interp_cpu
