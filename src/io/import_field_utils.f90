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
  use fld_file_data, only : fld_file_data_t
  use fld_file, only : fld_file_t
  use num_types, only : rp
  use field, only : field_t
  use field_list, only : field_list_t
  use utils, only : neko_error, extract_fld_file_index, &
       filename_chsuffix, NEKO_FNAME_LEN
  use logger, only : LOG_SIZE, neko_log
  use device, only : HOST_TO_DEVICE
  implicit none
  private

  public :: import_fields


contains

  !> Imports fields from an fld file, potentially with
  !! interpolation.
  !! @param fname The name of the fld file, e.g. "my_field0.f00019".
  !! @param mesh_fname The name of the fld file containing the spatial
  !! coordinates, if interpolation is enabled and fname does not already
  !! contain them.
  !! @param u The field on which to import the u component of the fld data.
  !! @param v The field on which to import the v component of the fld data.
  !! @param w The field on which to import the w component of the fld data.
  !! @param p The field on which to import the pressure field of the fld data.
  !! @param t The field on which to import the temperature field of the fld
  !! data.
  !! @param s_target_list Field list containing the fields on which to import
  !! the scalar fields of the fld data. Unless a list of target indices is
  !! provided in `s_index_list`, assigns field at position `i` in the list
  !! to scalar `i` in the fld file.
  !! @param s_index_list The list of scalars indices from which to load the
  !! fields provided in `s_target_list`. Must have the same size as
  !! `s_target_list`. For example, s_index_list = (/2,3/) will load scalar #2
  !! in `s_target_list%items(1)` and scalar #3 in `s_target_list%items(2)`.
  !! Index  0 corresponds to temperature by default. Therefore using
  !! `s_index_list = (/0/)` is equivalent to using the argument `t=...`.
  !! @param time Returns the time stamp of the fld file.
  !! @param interpolate Whether or not to interpolate the fld data.
  !! @param tolerance If interpolation is enabled, the tolerance to use for the
  !! point
  !! finding.
  !! @note If interpolation is disabled, space-to-space interpolation is still
  !! performed within each element to allow for seamless change of polynomial
  !! order for the same given mesh.
  !! @note This subroutine also takes care of data movement from host to
  !! to device when necessary, i.e. only the required fields are copied to
  !! device.
  subroutine import_fields(fname, mesh_fname, u, v, w, p, t, s_target_list, &
       s_index_list, time, interpolate, tolerance)
    character(len=*), intent(in) :: fname
    character(len=*), intent(in), optional :: mesh_fname
    type(field_t), pointer, intent(inout), optional :: u,v,w,p,t
    type(field_list_t), intent(inout), optional :: s_target_list
    integer, intent(in), optional :: s_index_list(:)
    real(kind=rp), intent(inout), optional :: time
    logical, intent(in), optional :: interpolate
    real(kind=rp), intent(in), optional :: tolerance

    character(len=LOG_SIZE) :: log_buf
    integer :: sample_idx, sample_mesh_idx, i
    character(len=NEKO_FNAME_LEN) :: fname_, mesh_fname_

    logical :: interpolate_, any_input_present

    type(fld_file_t) :: f
    type(fld_file_data_t) :: fld_data

    ! ---- Default values
    interpolate_ = .false.
    if (present(interpolate)) interpolate_ = interpolate
    mesh_fname_ = "none"
    if (present(mesh_fname)) mesh_fname_ = trim(mesh_fname)
    ! ----

    call neko_log%section("Import fields")
    call neko_log%message("File name     : " // trim(fname))
    write (log_buf, '(A,L1)') "Interpolation : ", interpolate_
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

    ! Select which fields to read in the fld file based on inputs
    f%read_mesh = interpolate_ 
    f%read_velocity = (present(u) .or. present(v) .or. present(w))
    f%read_pressure = present(p)
    f%read_temperature = present(t)
    f%read_scalars = present(s_target_list)

    ! If interpolate, check if we need to read the mesh file
    if (interpolate_) then

       if (present(tolerance)) then
          write (log_buf, '(A,ES12.6)') "Tolerance     : ", tolerance
          call neko_log%message(log_buf)
       end if

       ! If no mesh file is specified, use the default file name
       if (mesh_fname_ .eq. "none") then
          mesh_fname_ = trim(fname_)
          sample_mesh_idx = sample_idx
       else
          mesh_fname_ = trim(mesh_fname)

          ! Extract sample index from the mesh file name
          sample_mesh_idx = extract_fld_file_index(mesh_fname_, -1)

          if (sample_mesh_idx .eq. -1) then
             call neko_error("Invalid file name for the initial condition." // &
                  "The file format must be e.g. 'mean0.f00001'")
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

    end if

    ! Read the field file containing (u,v,w,p)
    call f%set_counter(sample_idx)
    call f%read(fld_data)

    ! Store the time stamp if it is requested
    if (present(time)) time = fld_data%time

    ! Interrupt here if we haven't provided any field as input
    any_input_present = (present(u) .or. present(v) .or. present(w) &
            .or. present(p) .or. present(t) .or. present(s_target_list))
    if (.not. any_input_present) return

    !
    ! Copy all vectors to device (GPU) since everything is read on the CPU
    !
    if (present(u)) call fld_data%u%copy_from(HOST_TO_DEVICE, .true.)
    if (present(v)) call fld_data%v%copy_from(HOST_TO_DEVICE, .true.)
    if (present(w)) call fld_data%w%copy_from(HOST_TO_DEVICE, .true.)
    if (present(p)) call fld_data%p%copy_from(HOST_TO_DEVICE, .true.)
    if (present(t)) call fld_data%t%copy_from(HOST_TO_DEVICE, .true.)
    if (present(s_target_list)) then

       if (present(s_index_list)) then
          if (size(s_index_list) .ne. s_target_list%size()) then
             call neko_error("Scalar lists must have same size!")
          end if

          do i = 1, size(s_index_list)
             ! Take care that if we set i=0 we want temperature
             if (s_index_list(i) .eq. 0) then
                call fld_data%t%copy_from(HOST_TO_DEVICE, .true.)
             else
                ! For scalar fields, require indices in 1:this%n_scalars
                if (s_index_list(i) < 1 .or. &
                     s_index_list(i) > fld_data%n_scalars) then
                   call neko_error("s_index_list entry out of bounds")
                end if
                call fld_data%s(s_index_list(i))%copy_from(HOST_TO_DEVICE, &
                     .true.)
             end if
          end do
       else
          do i = 1, s_target_list%size()
             call fld_data%s(i)%copy_from(HOST_TO_DEVICE, .true.)
          end do
       end if ! if present s_index_list

    end if ! present s_tgt


    if (interpolate_) then
       ! Sync coordinates to device for the interpolation
       call fld_data%x%copy_from(HOST_TO_DEVICE, .false.)
       call fld_data%y%copy_from(HOST_TO_DEVICE, .false.)
       call fld_data%z%copy_from(HOST_TO_DEVICE, .true.)
    end if

    ! Call the import of fields
    call fld_data%import_fields(u, v, w, p, t, s_target_list, s_index_list, &
         interpolate_, tolerance = tolerance)

    call neko_log%end_section()
    call fld_data%free()

  end subroutine import_fields

end module import_field_utils
