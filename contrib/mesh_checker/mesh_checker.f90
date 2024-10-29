! Copyright (c) 2018-2023, The Neko Authors
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
! Mesh diagnostics tool
program mesh_checker
  use neko
  implicit none

  character(len=NEKO_FNAME_LEN) :: inputchar, mesh_fname
  type(file_t) :: mesh_file
  type(mesh_t) :: msh
  integer :: argc, i, n_labeled
  character(len=LOG_SIZE) :: log_buf

  argc = command_argument_count()

  if (argc .lt. 1) then
     if (pe_rank .eq. 0) then
        write(*,*) 'Usage: ./mesh_checker mesh.nmsh'
     end if
     stop
  end if

  call neko_init

  call get_command_argument(1, inputchar)
  read(inputchar, *) mesh_fname

  mesh_file = file_t(trim(mesh_fname))

  call mesh_file%read(msh)


  if (pe_rank .eq. 0) then
      write(*,*) ''
      write(*,*) '--------------Size-------------'
      write(*,*) 'Number of elements: ', msh%glb_nelv
      write(*,*) 'Number of points:   ', msh%glb_mpts
      write(*,*) 'Number of faces:    ', msh%glb_mfcs
      write(*,*) 'Number of edges:    ', msh%glb_meds
      write(*,*) ''
      write(*,*) '--------------Zones------------'
      write(*,*) 'Number of built-in inlet faces:         ', msh%inlet%size
      write(*,*) 'Number of built-in wall faces:          ', msh%wall%size
      write(*,*) 'Number of built-in outlet faces:        ', msh%outlet%size
      write(*,*) 'Number of built-in outlet-normal faces: ', &
           msh%outlet_normal%size
      write(*,*) 'Number of built-in symmetry faces:      ', &
            msh%sympln%size
      write(*,*) 'Number of periodic faces:               ', msh%periodic%size

      write(*,*) 'Labeled zones: '
      do i = 1, size(msh%labeled_zones) 
         if (msh%labeled_zones(i)%size .gt. 0) then
            write(*,'(A,I2,A,I0,A)') '    Zone ', i, ': ', &
                 msh%labeled_zones(i)%size, ' faces'
            
         end if
      end do
  end if

  if (pe_rank .eq. 0) write(*,*) 'Done'
  call neko_finalize

end program mesh_checker

