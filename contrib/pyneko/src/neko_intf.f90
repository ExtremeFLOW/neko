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
module neko_intf
  use neko
  use json_case
  use json_module
  use, intrinsic :: iso_c_binding
  implicit none
  private

contains

  subroutine neko_intf_init() bind(c, name='init')
    character(len=LOG_SIZE) :: log_buf
    character(10) :: time
    character(8) :: date

    call neko_init()

    call date_and_time(time=time, date=date)           
    call neko_log%section("Session Information")       
    write(log_buf, '(A,A,A,A,1x,A,1x,A,A,A,A,A)') 'Start time: ',&
         time(1:2),':',time(3:4), '/', date(1:4),'-', date(5:6),'-',date(7:8)
    call neko_log%message(log_buf)
    write(log_buf, '(a)') 'Running on: '
    if (pe_size .lt. 1e1)  then
       write(log_buf(13:), '(i1,a)') pe_size, ' MPI '
       if (pe_size .eq. 1) then
          write(log_buf(19:), '(a)') 'rank'
       else
          write(log_buf(19:), '(a)') 'ranks'
       end if
    else if (pe_size .lt. 1e2) then
       write(log_buf(13:), '(i2,a)') pe_size, ' MPI ranks'
    else if (pe_size .lt. 1e3) then
       write(log_buf(13:), '(i3,a)') pe_size, ' MPI ranks'
    else if (pe_size .lt. 1e4) then
       write(log_buf(13:), '(i4,a)') pe_size, ' MPI ranks'
    else if (pe_size .lt. 1e5) then
       write(log_buf(13:), '(i5,a)') pe_size, ' MPI ranks'
    else
       write(log_buf(13:), '(i6,a)') pe_size, ' MPI ranks'
    end if
    call neko_log%message(log_buf)

    write(log_buf, '(a)') 'Bcknd type: '
    if (NEKO_BCKND_SX .eq. 1) then
       write(log_buf(13:), '(a)') 'SX-Aurora'
    else if (NEKO_BCKND_XSMM .eq. 1) then
       write(log_buf(13:), '(a)') 'CPU (libxsmm)'
    else if (NEKO_BCKND_CUDA .eq. 1) then
       write(log_buf(13:), '(a)') 'Accelerator (CUDA)'
    else if (NEKO_BCKND_HIP .eq. 1) then
       write(log_buf(13:), '(a)') 'Accelerator (HIP)'
    else if (NEKO_BCKND_OPENCL .eq. 1) then
       write(log_buf(13:), '(a)') 'Accelerator (OpenCL)'
    else
       write(log_buf(13:), '(a)') 'CPU'
    end if
    call neko_log%message(log_buf)

    if (NEKO_BCKND_HIP .eq. 1 .or. NEKO_BCKND_CUDA .eq. 1 .or. &
         NEKO_BCKND_OPENCL .eq. 1) then
       write(log_buf, '(a)') 'Dev. name : '
       call device_name(log_buf(13:))
       call neko_log%message(log_buf)
    end if

    write(log_buf, '(a)') 'Real type : '
    select case (rp)
    case (real32)
       write(log_buf(13:), '(a)') 'single precision'
    case (real64)
       write(log_buf(13:), '(a)') 'double precision'
    case (real128)
       write(log_buf(13:), '(a)') 'quad precision'
    end select
    call neko_log%message(log_buf)       
    call neko_log%end()
    call neko_log%newline

  end subroutine neko_intf_init

  subroutine neko_intf_finalize() bind(c, name="finalize")
    call neko_finalize()
  end subroutine neko_intf_finalize

  subroutine neko_intf_solve(pyneko_case, ilen) bind(c, name="solve")
    type(c_ptr) :: pyneko_case
    integer(c_int), value :: ilen
    character(len=:), allocatable :: fpyneko_case
    type(json_file) :: json_case
    type(case_t) :: neko_case
    
    

    if (c_associated(pyneko_case)) then
       block    
         character(kind=c_char,len=ilen+1),pointer :: s
         call c_f_pointer(pyneko_case, s)
         fpyneko_case = s(1:ilen)
         call json_case%load_from_string(fpyneko_case)
         deallocate(fpyneko_case)
         nullify(s)
       end block
    end if
    
    call json_case_create_neko_case(neko_case, json_case)
    call json_case%destroy()

    call neko_solve(neko_case)

    call case_free(neko_case)

  end subroutine neko_intf_solve

end module neko_intf
