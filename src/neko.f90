! Copyright (c) 2019-2021, The Neko Authors
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
!> Master module
module neko
  use num_types
  use parameters    
  use comm
  use utils
  use logger
  use math
  use speclib
  use dofmap
  use space
  use htable
  use uset
  use stack
  use tuple
  use mesh
  use mesh_field
  use map
  use mxm_wrapper
  use file
  use field
  use mpi_types
  use gather_scatter
  use coefs
  use bc
  use wall
  use dirichlet
  use krylov_fctry
  use precon_fctry
  use ax_helm_fctry
  use precon
  use ax_product
  use neko_config
  use case
  use sampler
  use output
  use simulation
  use operators
  use mathops
  use projection
  use user_intf
  use parmetis
  use signal
  use jobctrl
  use device
contains

  subroutine neko_init(C)
    type(case_t), intent(inout), optional :: C
    character(len=NEKO_FNAME_LEN) :: case_file
    character(len=LOG_SIZE) :: log_buf
    character(len=10) :: suffix
    character(10) :: time
    character(8) :: date
    integer :: argc

    call date_and_time(time=time, date=date)           

    call comm_init
    call mpi_types_init
    call jobctrl_init
    call device_init

    call neko_log%init()

    if (pe_rank .eq. 0) then
       write(*,*) ''
       write(*,*) '   _  __  ____  __ __  ____ '
       write(*,*) '  / |/ / / __/ / //_/ / __ \'
       write(*,*) ' /    / / _/  / ,<   / /_/ /'
       write(*,*) '/_/|_/ /___/ /_/|_|  \____/ '
       write(*,*) ''
       write(*,*) '(version: ', trim(NEKO_VERSION),')'
       write(*,*) trim(NEKO_BUILD_INFO)
       write(*,*) ''
    end if

    if (present(C)) then

       argc = command_argument_count()

       if ((argc .lt. 1) .or. (argc .gt. 1)) then
          if (pe_rank .eq. 0) write(*,*) 'Usage: ./neko <case file>'
          stop
       end if

       call get_command_argument(1, case_file)

       call filename_suffix(case_file, suffix)

       if (trim(suffix) .ne. 'case') then
          call neko_error('Invalid case file')
       end if

       !
       ! Job information
       !
       call neko_log%section("Job Information")       
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
          write(log_buf, '(a)') 'Bcknd type: '
          write(log_buf, '(a)') 'Dev. name : '
          call device_name(log_buf(13:))
          call neko_log%message(log_buf)
       end if

       write(log_buf, '(a)') 'Real type : '
       select case (rp)
       case (4)
          write(log_buf(13:), '(a)') 'single precision'
       case (8)
          write(log_buf(13:), '(a)') 'double precision'
       case (16)
          write(log_buf(13:), '(a)') 'quad precision'
       end select
       call neko_log%message(log_buf)

       call neko_log%end()

       !
       ! Create case
       !
       call case_init(C, case_file)
       
    end if
    
  end subroutine neko_init

  subroutine neko_finalize(C)
    type(case_t), intent(inout), optional :: C

    if (present(C)) then
       call case_free(C)
    end if
    
    call device_finalize
    call mpi_types_free
    call comm_free
  end subroutine neko_finalize

end module neko
