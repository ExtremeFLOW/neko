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
!> Importation of fields from fld files.
!! @details Import field from fld files. This tool is in a separate file
!! because of multiple import requirements.
module import_field_utils
  use fld_file_data, only: fld_file_data_t
  use file, only: file_t
  use num_types, only: rp
  use field, only: field_t
  use field_list, only: field_list_t
  use utils, only: neko_error, extract_fld_file_index, &
          filename_chsuffix, NEKO_FNAME_LEN
  use logger, only: LOG_SIZE, neko_log
  use device, only: HOST_TO_DEVICE 
  implicit none
  private
  
  public :: import_fields


contains

  !> Reads an fld file and import fields, with/without interpolation.
  subroutine import_fields(fname, mesh_fname, &
                  u, v, w, p, t, s_idx_list, s_tgt_list, interpolate, &
                  tolerance)
    character(len=*), intent(in) :: fname
    character(len=*), intent(in), optional :: mesh_fname
    type(field_t), pointer, intent(inout), optional :: u,v,w,p,t
    type(field_list_t), intent(in), optional :: s_tgt_list
    integer, intent(in), optional :: s_idx_list(:)
    logical, intent(in), optional :: interpolate
    real(kind=rp), intent(in) :: tolerance

    character(len=LOG_SIZE) :: log_buf
    integer :: sample_idx, sample_mesh_idx
    character(len=NEKO_FNAME_LEN) :: fname_, mesh_fname_

    type(file_t) :: f
    type(fld_file_data_t) :: fld_data

    call neko_log%section("Import fields") 
    call neko_log%message("File name     : " // trim(fname))
    write (log_buf, '(A,L1)') "Interpolation : ", interpolate
    call neko_log%message(log_buf)

    !
    ! Handling of file names and reading of data
    ! 

    ! Extract sample index from the file name
    sample_idx = extract_fld_file_index(fname, -1)

    if (sample_idx .eq. -1) &
         call neko_error("Invalid file name for the initial condition. The&
    & file format must be e.g. 'mean0.f00001'")

    ! Change from "field0.f000*" to "field0.fld" for the fld reader
    call filename_chsuffix(fname, fname_, 'fld')

    ! Initialize file object
    call f%init(trim(fname_))

    ! If interpolate, check if we need to read the mesh file
    if (interpolate) then

       write (log_buf, '(A,ES12.6)') "Tolerance     : ", tolerance
       call neko_log%message(log_buf)
       
       ! If no mesh file is specified, use the default file name
       if (mesh_fname .eq. "none") then
          mesh_fname_ = trim(fname_)
          sample_mesh_idx = sample_idx
       else
          mesh_fname_ = trim(mesh_fname)

          ! Extract sample index from the mesh file name
          sample_mesh_idx = extract_fld_file_index(mesh_fname_, -1)

          if (sample_mesh_idx .eq. -1) then
             call neko_error("Invalid file name for the initial condition. &
             &The file format must be e.g. 'mean0.f00001'")
          end if

          write (log_buf, '(A,A)') "Mesh file     : ", &
               trim(mesh_fname_)
          call neko_log%message(log_buf)

       end if ! if mesh_file_name .eq. none

       ! Read the mesh coordinates if they are not in our fld file
       if (sample_mesh_idx .ne. sample_idx) then
          call f%set_counter(sample_mesh_idx)
          call f%read(fld_data)
       end if
       
       ! Sync coordinates to device for the interpolation
       call fld_data%x%copy_from(HOST_TO_DEVICE, .false.)
       call fld_data%y%copy_from(HOST_TO_DEVICE, .false.)
       call fld_data%z%copy_from(HOST_TO_DEVICE, .true.)

    end if

    ! Read the field file containing (u,v,w,p)
    call f%set_counter(sample_idx)
    call f%read(fld_data)

    ! Call the import of fields
    call fld_data%import_fields(u, v, w, p, t, s_idx_list, s_tgt_list, &
            interpolate, tolerance)

    call neko_log%end_section() 

  end subroutine import_fields
end module import_field_utils
