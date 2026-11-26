! Copyright (c) 2021-2026, The Neko Authors
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
!> Neko legacy checkpoint file format.
!! @details This module defines the interface for reading and writing `.chkp`
!! files.
module chkp_file
  use generic_file, only : generic_file_t
  use field_series, only : field_series_t
  use checkpoint, only : chkp_t
  use checkpoint_payload, only : checkpoint_payload_t, checkpoint_array_t, &
       checkpoint_mesh_array_t
  use num_types, only : rp, dp, i8
  use field, only : field_t
  use dofmap, only : dofmap_t
  use utils, only : neko_error, filename_suffix_pos
  use space, only : space_t, GLL
  use mesh, only : mesh_t
  use math, only : rzero
  use interpolation, only : interpolator_t
  use neko_mpi_types, only : MPI_REAL_PREC_SIZE, MPI_INTEGER_SIZE, &
       MPI_DOUBLE_PRECISION_SIZE
  use global_interpolation, only : global_interpolation_t
  use logger, only : neko_log, NEKO_LOG_VERBOSE
  use comm, only : NEKO_COMM, pe_rank, pe_size, MPI_REAL_PRECISION
  use mpi_f08, only : MPI_File, MPI_Status, MPI_OFFSET_KIND, MPI_MODE_CREATE, &
       MPI_MODE_RDONLY, MPI_MODE_WRONLY, MPI_INFO_NULL, MPI_INTEGER, &
       MPI_INTEGER8, MPI_CHARACTER, MPI_DOUBLE_PRECISION, MPI_LOGICAL, &
       MPI_SUCCESS, MPI_SUM, MPI_File_open, MPI_File_close, &
       MPI_File_read_all, MPI_File_write_all, MPI_File_read_at_all, &
       MPI_File_write_at_all, MPI_File_get_size, MPI_Bcast, MPI_Exscan, &
       MPI_Allreduce
  implicit none
  private

  !> Bit in the `optional_fields` header word flagging the generic payload
  !! section. Readers that predate it ignore the bit and the trailing data.
  integer, parameter :: CHKP_PAYLOADS_BIT = 64
  !> Marker closing the generic payload section footer ("NEKOPAY1").
  integer(kind=i8), parameter :: CHKP_PAYLOADS_MAGIC = &
       int(z'4E454B4F50415931', i8)
  !> Size in bytes of an 8-byte integer or real in the payload section.
  integer(kind=i8), parameter :: CHKP_WORD_BYTES = 8_i8
  !> Array kinds stored in the generic payload section.
  integer, parameter :: CHKP_ARRAY_REPLICATED = 0
  integer, parameter :: CHKP_ARRAY_RESIZABLE = 1

  !> Temporary view of the payloads required by the positional .chkp format.
  type :: legacy_checkpoint_view_t
     !> Fluid and scalar fields.
     type(field_t), pointer :: u => null(), v => null(), w => null()
     type(field_t), pointer :: p => null(), s => null()
     !> Fluid and scalar extrapolation fields.
     type(field_t), pointer :: abx1 => null(), abx2 => null()
     type(field_t), pointer :: aby1 => null(), aby2 => null()
     type(field_t), pointer :: abz1 => null(), abz2 => null()
     type(field_t), pointer :: abs1 => null(), abs2 => null()
     !> Mesh-velocity fields.
     type(field_t), pointer :: wm_x => null(), wm_y => null(), wm_z => null()
     !> Fluid, scalar, and mesh-velocity field histories.
     type(field_series_t), pointer :: ulag => null(), vlag => null()
     type(field_series_t), pointer :: wlag => null(), slag => null()
     type(field_series_t), pointer :: wm_x_lag => null()
     type(field_series_t), pointer :: wm_y_lag => null()
     type(field_series_t), pointer :: wm_z_lag => null()
     !> Time history.
     real(kind=dp), pointer :: tlag(:) => null(), dtlag(:) => null()
     !> Mesh coordinates and lagged mass matrices.
     real(kind=rp), pointer :: msh_x(:) => null(), msh_y(:) => null()
     real(kind=rp), pointer :: msh_z(:) => null(), Blag(:) => null()
     real(kind=rp), pointer :: Blaglag(:) => null()
     !> Rigid-body state.
     real(kind=rp), pointer :: pivot_pos(:) => null()
     real(kind=rp), pointer :: pivot_vel_lag(:) => null()
     real(kind=rp), pointer :: basis_pos(:) => null()
     real(kind=rp), pointer :: basis_vel_lag(:) => null()
  end type legacy_checkpoint_view_t

  !> Interface for Neko legacy checkpoint files.
  type, public, extends(generic_file_t) :: chkp_file_t
     !> Function space in the loaded checkpoint file.
     type(space_t), pointer :: chkp_Xh
     !> Function space used in the simulation.
     type(space_t), pointer :: sim_Xh
     !> Interpolation used when only the polynomial order changes.
     type(interpolator_t) :: space_interp
     !> Interpolation used when the mesh changes.
     type(global_interpolation_t) :: global_interp
     !> Whether the checkpoint mesh differs from the simulation mesh.
     logical :: mesh2mesh
   contains
     !> Load a Neko legacy checkpoint.
     procedure :: read => chkp_file_read
     !> Read and optionally interpolate one legacy checkpoint field.
     procedure :: read_field => chkp_read_field
     !> Write a Neko legacy checkpoint.
     procedure :: write => chkp_file_write
     !> Get the physical file name generated by the next write.
     procedure :: get_next_output_fname => chkp_file_get_next_output_fname
     !> Set whether checkpoint files may be overwritten.
     procedure :: set_overwrite => chkp_file_set_overwrite
  end type chkp_file_t

contains

  !> Get the physical file name generated by the next write.
  function chkp_file_get_next_output_fname(this) result(fname)
    class(chkp_file_t), intent(in) :: this
    character(len=1024) :: fname
    character(len=5) :: id_str
    integer :: counter, suffix_pos

    if (this%overwrite) then
       fname = this%get_fname()
    else
       counter = this%get_counter()
       if (counter .eq. -1) then
          counter = this%get_start_counter()
       else
          counter = counter + 1
       end if

       fname = this%get_base_fname()
       suffix_pos = filename_suffix_pos(fname)
       write(id_str, '(i5.5)') counter
       fname = trim(fname(1:suffix_pos-1)) // id_str // fname(suffix_pos:)
    end if

  end function chkp_file_get_next_output_fname

  !> Write a Neko legacy checkpoint.
  !! @param data Checkpoint data to write.
  !! @param[optional] t Simulation time.
  subroutine chkp_file_write(this, data, t)
    class(chkp_file_t), intent(inout) :: this
    class(*), target, intent(in) :: data
    real(kind=dp), intent(in), optional :: t
    real(kind=dp) :: time
    character(len=5) :: id_str
    character(len=1024) :: fname
    integer :: ierr, suffix_pos, optional_fields
    type(field_t), pointer :: u, v, w, p, s
    type(field_t), pointer :: wm_x => null()
    type(field_t), pointer :: wm_y => null()
    type(field_t), pointer :: wm_z => null()
    type(field_t), pointer :: abx1, abx2
    type(field_t), pointer :: aby1, aby2
    type(field_t), pointer :: abz1, abz2
    type(field_t), pointer :: abs1, abs2
    type(field_series_t), pointer :: ulag => null()
    type(field_series_t), pointer :: vlag => null()
    type(field_series_t), pointer :: wlag => null()
    type(field_series_t), pointer :: wm_x_lag => null()
    type(field_series_t), pointer :: wm_y_lag => null()
    type(field_series_t), pointer :: wm_z_lag => null()
    type(field_series_t), pointer :: slag => null()
    real(kind=rp), pointer :: msh_x(:) => null()
    real(kind=rp), pointer :: msh_y(:) => null()
    real(kind=rp), pointer :: msh_z(:) => null()
    real(kind=rp), pointer :: Blag(:) => null()
    real(kind=rp), pointer :: Blaglag(:) => null()
    real(kind=rp), pointer :: pivot_pos(:) => null()
    real(kind=rp), pointer :: pivot_vel_lag(:) => null()
    real(kind=rp), pointer :: basis_pos(:) => null()
    real(kind=rp), pointer :: basis_vel_lag(:) => null()
    real(kind=dp), pointer :: dtlag(:), tlag(:)
    type(mesh_t), pointer :: msh
    type(MPI_Status) :: status
    type(MPI_File) :: fh
    integer (kind=MPI_OFFSET_KIND) :: mpi_offset, byte_offset
    integer(kind=i8) :: n_glb_dofs, dof_offset
    logical :: write_lag, write_scalar, write_dtlag
    logical :: write_ale
    logical :: write_scalarlag, write_abvel
    integer :: i
    type(legacy_checkpoint_view_t) :: view
    type(chkp_t), pointer :: chkp
    logical :: write_payloads

    if (present(t)) then
       time = t
    else
       time = 0.0_dp
    end if

    select type (data)
    type is (chkp_t)

       if (data%scalar_payload_count() .gt. 1) then
          call neko_error("The legacy .chkp format supports at most one " // &
               "scalar; use HDF5 checkpointing for multiple scalars")
       end if

       call chkp_build_legacy_view(data, view)
       if ( .not. associated(view%u) .or. &
            .not. associated(view%v) .or. &
            .not. associated(view%w) .or. &
            .not. associated(view%p) ) then
          call neko_error('Checkpoint not initialized')
       end if

       u => view%u
       v => view%v
       w => view%w
       p => view%p
       msh => u%msh

       optional_fields = 0

       if (associated(view%ulag)) then
          ulag => view%ulag
          vlag => view%vlag
          wlag => view%wlag
          write_lag = .true.
          optional_fields = optional_fields + 1
       else
          write_lag = .false.
       end if

       if (associated(view%s)) then
          s => view%s
          write_scalar = .true.
          optional_fields = optional_fields + 2
       else
          write_scalar = .false.
       end if

       if (associated(view%tlag)) then
          tlag => view%tlag
          dtlag => view%dtlag
          write_dtlag = .true.
          optional_fields = optional_fields + 4
       else
          write_dtlag = .false.
       end if

       write_abvel = .false.
       if (associated(view%abx1)) then
          abx1 => view%abx1
          abx2 => view%abx2
          aby1 => view%aby1
          aby2 => view%aby2
          abz1 => view%abz1
          abz2 => view%abz2
          optional_fields = optional_fields + 8
          write_abvel = .true.
       end if

       write_scalarlag = .false.
       if (associated(view%abs1)) then
          slag => view%slag
          abs1 => view%abs1
          abs2 => view%abs2
          optional_fields = optional_fields + 16
          write_scalarlag = .true.
       end if

       write_ale = .false.
       if (associated(view%wm_x)) then
          msh_x => view%msh_x
          msh_y => view%msh_y
          msh_z => view%msh_z
          wm_x => view%wm_x
          wm_y => view%wm_y
          wm_z => view%wm_z
          wm_x_lag => view%wm_x_lag
          wm_y_lag => view%wm_y_lag
          wm_z_lag => view%wm_z_lag
          Blag => view%Blag
          Blaglag => view%Blaglag
          pivot_pos => view%pivot_pos
          pivot_vel_lag => view%pivot_vel_lag
          basis_pos => view%basis_pos
          basis_vel_lag => view%basis_vel_lag
          optional_fields = optional_fields + 32
          write_ale = .true.
       end if

       chkp => data
       write_payloads = chkp_generic_payload_count(chkp) .gt. 0
       if (write_payloads) then
          optional_fields = optional_fields + CHKP_PAYLOADS_BIT
       end if

    class default
       call neko_error('Invalid data')
    end select

    dof_offset = int(msh%offset_el, i8) * int(u%Xh%lx * u%Xh%ly * u%Xh%lz, i8)
    n_glb_dofs = int(u%Xh%lx * u%Xh%ly * u%Xh%lz, i8) * int(msh%glb_nelv, i8)

    ! Retrieve the filename and increment the counter if we are not overwriting
    if (.not. this%overwrite) call this%increment_counter()
    fname = trim(this%get_fname())

    call MPI_File_open(NEKO_COMM, trim(fname), &
         MPI_MODE_WRONLY + MPI_MODE_CREATE, MPI_INFO_NULL, fh, ierr)
    call MPI_File_write_all(fh, msh%glb_nelv, 1, MPI_INTEGER, status, ierr)
    call MPI_File_write_all(fh, msh%gdim, 1, MPI_INTEGER, status, ierr)
    call MPI_File_write_all(fh, u%Xh%lx, 1, MPI_INTEGER, status, ierr)
    call MPI_File_write_all(fh, optional_fields, 1, MPI_INTEGER, status, ierr)
    call MPI_File_write_all(fh, time, 1, MPI_DOUBLE_PRECISION, status, ierr)


    !
    ! Dump mandatory checkpoint data
    !

    byte_offset = 4_i8 * int(MPI_INTEGER_SIZE, i8) + &
         int(MPI_DOUBLE_PRECISION_SIZE, i8)
    byte_offset = byte_offset + &
         dof_offset * int(MPI_REAL_PREC_SIZE, i8)
    call MPI_File_write_at_all(fh, byte_offset,u%x, u%dof%size(), &
         MPI_REAL_PRECISION, status, ierr)
    mpi_offset = 4_i8 * int(MPI_INTEGER_SIZE, i8) + &
         int(MPI_DOUBLE_PRECISION_SIZE, i8)
    mpi_offset = mpi_offset +&
         n_glb_dofs * int(MPI_REAL_PREC_SIZE, i8)

    byte_offset = mpi_offset + &
         dof_offset * int(MPI_REAL_PREC_SIZE, i8)
    call MPI_File_write_at_all(fh, byte_offset, v%x, v%dof%size(), &
         MPI_REAL_PRECISION, status, ierr)
    mpi_offset = mpi_offset + n_glb_dofs * int(MPI_REAL_PREC_SIZE, i8)

    byte_offset = mpi_offset + &
         dof_offset * int(MPI_REAL_PREC_SIZE, i8)
    call MPI_File_write_at_all(fh, byte_offset, w%x, w%dof%size(), &
         MPI_REAL_PRECISION, status, ierr)
    mpi_offset = mpi_offset + n_glb_dofs * int(MPI_REAL_PREC_SIZE, i8)

    byte_offset = mpi_offset + &
         dof_offset * int(MPI_REAL_PREC_SIZE, i8)
    call MPI_File_write_at_all(fh, byte_offset, p%x, p%dof%size(), &
         MPI_REAL_PRECISION, status, ierr)
    mpi_offset = mpi_offset + n_glb_dofs * int(MPI_REAL_PREC_SIZE, i8)

    !
    ! Dump optional payload
    !

    if (write_lag) then

       do i = 1, ulag%size()
          byte_offset = mpi_offset + &
               dof_offset * int(MPI_REAL_PREC_SIZE, i8)
          ! We should not need this extra associate block, and it works
          ! great without it for GNU, Intel, NEC and Cray, but throws an
          ! ICE with NAG.
          associate (x => ulag%lf(i)%x)
            call MPI_File_write_at_all(fh, byte_offset, x, &
                 ulag%lf(i)%dof%size(), MPI_REAL_PRECISION, status, ierr)
          end associate
          mpi_offset = mpi_offset + n_glb_dofs * int(MPI_REAL_PREC_SIZE, i8)
       end do

       do i = 1, vlag%size()
          byte_offset = mpi_offset + &
               dof_offset * int(MPI_REAL_PREC_SIZE, i8)
          ! We should not need this extra associate block, and it works
          ! great without it for GNU, Intel, NEC and Cray, but throws an
          ! ICE with NAG.
          associate (x => vlag%lf(i)%x)
            call MPI_File_write_at_all(fh, byte_offset, x, &
                 vlag%lf(i)%dof%size(), MPI_REAL_PRECISION, status, ierr)
          end associate
          mpi_offset = mpi_offset + n_glb_dofs * int(MPI_REAL_PREC_SIZE, i8)
       end do

       do i = 1, wlag%size()
          byte_offset = mpi_offset + &
               dof_offset * int(MPI_REAL_PREC_SIZE, i8)
          ! We should not need this extra associate block, and it works
          ! great without it for GNU, Intel, NEC and Cray, but throws an
          ! ICE with NAG.
          associate (x => wlag%lf(i)%x)
            call MPI_File_write_at_all(fh, byte_offset, x, &
                 wlag%lf(i)%dof%size(), MPI_REAL_PRECISION, status, ierr)
          end associate
          mpi_offset = mpi_offset + n_glb_dofs * int(MPI_REAL_PREC_SIZE, i8)
       end do

    end if

    if (write_scalar) then
       byte_offset = mpi_offset + &
            dof_offset * int(MPI_REAL_PREC_SIZE, i8)
       call MPI_File_write_at_all(fh, byte_offset, s%x, p%dof%size(), &
            MPI_REAL_PRECISION, status, ierr)
       mpi_offset = mpi_offset + n_glb_dofs * int(MPI_REAL_PREC_SIZE, i8)
    end if

    if (write_dtlag) then
       call MPI_File_write_at_all(fh, mpi_offset, tlag, 10, &
            MPI_DOUBLE_PRECISION, status, ierr)
       mpi_offset = mpi_offset + 10_i8 * int(MPI_DOUBLE_PRECISION_SIZE, i8)
       call MPI_File_write_at_all(fh, mpi_offset, dtlag, 10, &
            MPI_DOUBLE_PRECISION, status, ierr)
       mpi_offset = mpi_offset + 10_i8 * int(MPI_DOUBLE_PRECISION_SIZE, i8)
    end if

    if (write_abvel) then
       byte_offset = mpi_offset + &
            dof_offset * int(MPI_REAL_PREC_SIZE, i8)
       call MPI_File_write_at_all(fh, byte_offset, abx1%x, abx1%dof%size(), &
            MPI_REAL_PRECISION, status, ierr)
       mpi_offset = mpi_offset + n_glb_dofs * int(MPI_REAL_PREC_SIZE, i8)
       byte_offset = mpi_offset + &
            dof_offset * int(MPI_REAL_PREC_SIZE, i8)
       call MPI_File_write_at_all(fh, byte_offset, abx2%x, abx1%dof%size(), &
            MPI_REAL_PRECISION, status, ierr)
       mpi_offset = mpi_offset + n_glb_dofs * int(MPI_REAL_PREC_SIZE, i8)
       byte_offset = mpi_offset + &
            dof_offset * int(MPI_REAL_PREC_SIZE, i8)
       call MPI_File_write_at_all(fh, byte_offset, aby1%x, abx1%dof%size(), &
            MPI_REAL_PRECISION, status, ierr)
       mpi_offset = mpi_offset + n_glb_dofs * int(MPI_REAL_PREC_SIZE, i8)
       byte_offset = mpi_offset + &
            dof_offset * int(MPI_REAL_PREC_SIZE, i8)
       call MPI_File_write_at_all(fh, byte_offset, aby2%x, abx1%dof%size(), &
            MPI_REAL_PRECISION, status, ierr)
       mpi_offset = mpi_offset + n_glb_dofs * int(MPI_REAL_PREC_SIZE, i8)
       byte_offset = mpi_offset + &
            dof_offset * int(MPI_REAL_PREC_SIZE, i8)
       call MPI_File_write_at_all(fh, byte_offset, abz1%x, abx1%dof%size(), &
            MPI_REAL_PRECISION, status, ierr)
       mpi_offset = mpi_offset + n_glb_dofs * int(MPI_REAL_PREC_SIZE, i8)
       byte_offset = mpi_offset + &
            dof_offset * int(MPI_REAL_PREC_SIZE, i8)
       call MPI_File_write_at_all(fh, byte_offset, abz2%x, abx1%dof%size(), &
            MPI_REAL_PRECISION, status, ierr)
       mpi_offset = mpi_offset + n_glb_dofs * int(MPI_REAL_PREC_SIZE, i8)
    end if

    if (write_scalarlag) then
       do i = 1, slag%size()
          byte_offset = mpi_offset + &
               dof_offset * int(MPI_REAL_PREC_SIZE, i8)
          ! We should not need this extra associate block, and it works
          ! great without it for GNU, Intel, NEC and Cray, but throws an
          ! ICE with NAG.
          associate (x => slag%lf(i)%x)
            call MPI_File_write_at_all(fh, byte_offset, x, &
                 slag%lf(i)%dof%size(), MPI_REAL_PRECISION, status, ierr)
          end associate
          mpi_offset = mpi_offset + n_glb_dofs * int(MPI_REAL_PREC_SIZE, i8)
       end do

       byte_offset = mpi_offset + &
            dof_offset * int(MPI_REAL_PREC_SIZE, i8)
       call MPI_File_write_at_all(fh, byte_offset, abs1%x, abx1%dof%size(), &
            MPI_REAL_PRECISION, status, ierr)
       mpi_offset = mpi_offset + n_glb_dofs * int(MPI_REAL_PREC_SIZE, i8)
       byte_offset = mpi_offset + &
            dof_offset * int(MPI_REAL_PREC_SIZE, i8)
       call MPI_File_write_at_all(fh, byte_offset, abs2%x, abx1%dof%size(), &
            MPI_REAL_PRECISION, status, ierr)
       mpi_offset = mpi_offset + n_glb_dofs * int(MPI_REAL_PREC_SIZE, i8)
    end if

    if (write_ale) then
       byte_offset = mpi_offset + &
            dof_offset * int(MPI_REAL_PREC_SIZE, i8)
       call MPI_File_write_at_all(fh, byte_offset, msh_x, size(msh_x), &
            MPI_REAL_PRECISION, status, ierr)
       mpi_offset = mpi_offset + n_glb_dofs * int(MPI_REAL_PREC_SIZE, i8)

       byte_offset = mpi_offset + &
            dof_offset * int(MPI_REAL_PREC_SIZE, i8)
       call MPI_File_write_at_all(fh, byte_offset, msh_y, size(msh_y), &
            MPI_REAL_PRECISION, status, ierr)
       mpi_offset = mpi_offset + n_glb_dofs * int(MPI_REAL_PREC_SIZE, i8)

       byte_offset = mpi_offset + &
            dof_offset * int(MPI_REAL_PREC_SIZE, i8)
       call MPI_File_write_at_all(fh, byte_offset, msh_z, size(msh_z), &
            MPI_REAL_PRECISION, status, ierr)
       mpi_offset = mpi_offset + n_glb_dofs * int(MPI_REAL_PREC_SIZE, i8)

       byte_offset = mpi_offset + &
            dof_offset * int(MPI_REAL_PREC_SIZE, i8)
       call MPI_File_write_at_all(fh, byte_offset, wm_x%x, wm_x%dof%size(), &
            MPI_REAL_PRECISION, status, ierr)
       mpi_offset = mpi_offset + n_glb_dofs * int(MPI_REAL_PREC_SIZE, i8)

       byte_offset = mpi_offset + &
            dof_offset * int(MPI_REAL_PREC_SIZE, i8)
       call MPI_File_write_at_all(fh, byte_offset, wm_y%x, wm_y%dof%size(), &
            MPI_REAL_PRECISION, status, ierr)
       mpi_offset = mpi_offset + n_glb_dofs * int(MPI_REAL_PREC_SIZE, i8)

       byte_offset = mpi_offset + &
            dof_offset * int(MPI_REAL_PREC_SIZE, i8)
       call MPI_File_write_at_all(fh, byte_offset, wm_z%x, wm_z%dof%size(), &
            MPI_REAL_PRECISION, status, ierr)
       mpi_offset = mpi_offset + n_glb_dofs * int(MPI_REAL_PREC_SIZE, i8)

       do i = 1, wm_x_lag%size()
          byte_offset = mpi_offset + &
               dof_offset * int(MPI_REAL_PREC_SIZE, i8)
          ! We should not need this extra associate block, and it works
          ! great without it for GNU, Intel, NEC and Cray, but throws an
          ! ICE with NAG.
          associate (x => wm_x_lag%lf(i)%x)
            call MPI_File_write_at_all(fh, byte_offset, x, &
                 wm_x_lag%lf(i)%dof%size(), MPI_REAL_PRECISION, status, ierr)
          end associate
          mpi_offset = mpi_offset + n_glb_dofs * int(MPI_REAL_PREC_SIZE, i8)
       end do

       do i = 1, wm_y_lag%size()
          byte_offset = mpi_offset + &
               dof_offset * int(MPI_REAL_PREC_SIZE, i8)
          ! We should not need this extra associate block, and it works
          ! great without it for GNU, Intel, NEC and Cray, but throws an
          ! ICE with NAG.
          associate (x => wm_y_lag%lf(i)%x)
            call MPI_File_write_at_all(fh, byte_offset, x, &
                 wm_y_lag%lf(i)%dof%size(), MPI_REAL_PRECISION, status, ierr)
          end associate
          mpi_offset = mpi_offset + n_glb_dofs * int(MPI_REAL_PREC_SIZE, i8)
       end do

       do i = 1, wm_z_lag%size()
          byte_offset = mpi_offset + &
               dof_offset * int(MPI_REAL_PREC_SIZE, i8)
          ! We should not need this extra associate block, and it works
          ! great without it for GNU, Intel, NEC and Cray, but throws an
          ! ICE with NAG.
          associate (x => wm_z_lag%lf(i)%x)
            call MPI_File_write_at_all(fh, byte_offset, x, &
                 wm_z_lag%lf(i)%dof%size(), MPI_REAL_PRECISION, status, ierr)
          end associate
          mpi_offset = mpi_offset + n_glb_dofs * int(MPI_REAL_PREC_SIZE, i8)
       end do

       byte_offset = mpi_offset + &
            dof_offset * int(MPI_REAL_PREC_SIZE, i8)
       call MPI_File_write_at_all(fh, byte_offset, Blag, size(Blag), &
            MPI_REAL_PRECISION, status, ierr)
       mpi_offset = mpi_offset + n_glb_dofs * int(MPI_REAL_PREC_SIZE, i8)

       byte_offset = mpi_offset + &
            dof_offset * int(MPI_REAL_PREC_SIZE, i8)
       call MPI_File_write_at_all(fh, byte_offset, Blaglag, size(Blaglag), &
            MPI_REAL_PRECISION, status, ierr)
       mpi_offset = mpi_offset + n_glb_dofs * int(MPI_REAL_PREC_SIZE, i8)

       call MPI_File_write_at_all(fh, mpi_offset, pivot_pos, size(pivot_pos), &
            MPI_REAL_PRECISION, status, ierr)
       mpi_offset = mpi_offset + &
            int(size(pivot_pos), i8) * int(MPI_REAL_PREC_SIZE, i8)
       call MPI_File_write_at_all(fh, mpi_offset, pivot_vel_lag, &
            size(pivot_vel_lag), MPI_REAL_PRECISION, status, ierr)
       mpi_offset = mpi_offset + &
            int(size(pivot_vel_lag), i8) * int(MPI_REAL_PREC_SIZE, i8)
       call MPI_File_write_at_all(fh, mpi_offset, basis_pos, &
            size(basis_pos), MPI_REAL_PRECISION, status, ierr)
       mpi_offset = mpi_offset + &
            int(size(basis_pos), i8) * int(MPI_REAL_PREC_SIZE, i8)
       call MPI_File_write_at_all(fh, mpi_offset, basis_vel_lag, &
            size(basis_vel_lag), MPI_REAL_PRECISION, status, ierr)
       mpi_offset = mpi_offset + &
            int(size(basis_vel_lag), i8) * int(MPI_REAL_PREC_SIZE, i8)
    end if

    if (write_payloads) then
       call chkp_write_generic_payloads(fh, mpi_offset, chkp)
    end if

    call MPI_File_close(fh, ierr)

    if (ierr .ne. MPI_SUCCESS) then
       call neko_error('Error writing checkpoint file ' // trim(fname))
    end if

  end subroutine chkp_file_write

  !> Load a Neko legacy checkpoint.
  !! @param data Checkpoint data to populate.
  subroutine chkp_file_read(this, data)
    class(chkp_file_t) :: this
    class(*), target, intent(inout) :: data
    type(chkp_t), pointer :: chkp
    character(len=5) :: id_str
    character(len=1024) :: fname
    integer :: ierr, suffix_pos
    type(field_t), pointer :: u, v, w, p, s
    type(field_t), pointer :: wm_x => null()
    type(field_t), pointer :: wm_y => null()
    type(field_t), pointer :: wm_z => null()
    type(field_series_t), pointer :: ulag => null()
    type(field_series_t), pointer :: vlag => null()
    type(field_series_t), pointer :: wlag => null()
    type(field_series_t), pointer :: wm_x_lag => null()
    type(field_series_t), pointer :: wm_y_lag => null()
    type(field_series_t), pointer :: wm_z_lag => null()
    type(field_series_t), pointer :: slag => null()
    type(mesh_t), pointer :: msh
    type(MPI_Status) :: status
    type(MPI_File) :: fh
    type(field_t), pointer :: abx1, abx2
    type(field_t), pointer :: aby1, aby2
    type(field_t), pointer :: abz1, abz2
    type(field_t), pointer :: abs1, abs2
    real(kind=rp), pointer :: Blag(:) => null()
    real(kind=rp), pointer :: Blaglag(:) => null()
    real(kind=rp), pointer :: msh_x(:) => null()
    real(kind=rp), pointer :: msh_y(:) => null()
    real(kind=rp), pointer :: msh_z(:) => null()
    real(kind=rp), pointer :: pivot_pos(:) => null()
    real(kind=rp), pointer :: pivot_vel_lag(:) => null()
    real(kind=rp), pointer :: basis_pos(:) => null()
    real(kind=rp), pointer :: basis_vel_lag(:) => null()
    real(kind=rp), allocatable :: x_coord(:,:,:,:)
    real(kind=rp), allocatable :: y_coord(:,:,:,:)
    real(kind=rp), allocatable :: z_coord(:,:,:,:)
    real(kind=dp), pointer :: dtlag(:), tlag(:)
    integer (kind=MPI_OFFSET_KIND) :: mpi_offset, byte_offset
    integer(kind=i8) :: n_glb_dofs, dof_offset
    integer :: glb_nelv, gdim, lx, have_lag, have_scalar, nel
    integer :: optional_fields, have_dtlag
    integer :: have_abvel, have_scalarlag
    integer :: have_ale, have_payloads
    logical :: read_lag, read_scalar, read_dtlag, read_abvel, read_scalarlag
    logical :: read_ale
    real(kind=dp) :: tol
    real(kind=rp) :: center_x, center_y, center_z
    integer :: i, e
    type(dofmap_t) :: dof
    type(legacy_checkpoint_view_t) :: view

    call this%check_exists()

    select type (data)
    type is (chkp_t)

       if (data%scalar_payload_count() .gt. 1) then
          call neko_error("The legacy .chkp format supports at most one " // &
               "scalar; use HDF5 checkpointing for multiple scalars")
       end if

       call chkp_build_legacy_view(data, view)
       if ( .not. associated(view%u) .or. &
            .not. associated(view%v) .or. &
            .not. associated(view%w) .or. &
            .not. associated(view%p) ) then
          call neko_error('Checkpoint not initialized')
       end if

       u => view%u
       v => view%v
       w => view%w
       p => view%p
       this%chkp_Xh => data%previous_Xh
       !> If checkpoint used another mesh
       if (allocated(data%previous_mesh%elements)) then
          msh => data%previous_mesh
          this%mesh2mesh = .true.
          tol = data%mesh2mesh_tol
       else !< The checkpoint was written on the same mesh
          msh => u%msh
          this%mesh2mesh = .false.
       end if

       if (associated(view%ulag)) then
          ulag => view%ulag
          vlag => view%vlag
          wlag => view%wlag
          read_lag = .true.
       else
          read_lag = .false.
       end if

       if (associated(view%s)) then
          s => view%s
          read_scalar = .true.
       else
          read_scalar = .false.
       end if
       if (associated(view%dtlag)) then
          dtlag => view%dtlag
          tlag => view%tlag
          read_dtlag = .true.
       else
          read_dtlag = .false.
       end if
       read_abvel = .false.
       if (associated(view%abx1)) then
          abx1 => view%abx1
          abx2 => view%abx2
          aby1 => view%aby1
          aby2 => view%aby2
          abz1 => view%abz1
          abz2 => view%abz2
          read_abvel = .true.
       end if
       read_scalarlag = .false.
       if (associated(view%abs1)) then
          slag => view%slag
          abs1 => view%abs1
          abs2 => view%abs2
          read_scalarlag = .true.
       end if

       read_ale = .false.
       if (associated(view%wm_x)) then
          msh_x => view%msh_x
          msh_y => view%msh_y
          msh_z => view%msh_z
          wm_x => view%wm_x
          wm_y => view%wm_y
          wm_z => view%wm_z
          wm_x_lag => view%wm_x_lag
          wm_y_lag => view%wm_y_lag
          wm_z_lag => view%wm_z_lag
          Blag => view%Blag
          Blaglag => view%Blaglag
          pivot_pos => view%pivot_pos
          pivot_vel_lag => view%pivot_vel_lag
          basis_pos => view%basis_pos
          basis_vel_lag => view%basis_vel_lag
          read_ale = .true.
       end if

       chkp => data

    class default
       call neko_error('Invalid data')
    end select

    fname = trim(this%get_fname())
    call neko_log%message("Reading checkpoint from file: " // trim(fname), &
         NEKO_LOG_VERBOSE)
    call MPI_File_open(NEKO_COMM, trim(fname), &
         MPI_MODE_RDONLY, MPI_INFO_NULL, fh, ierr)
    call MPI_File_read_all(fh, glb_nelv, 1, MPI_INTEGER, status, ierr)
    call MPI_File_read_all(fh, gdim, 1, MPI_INTEGER, status, ierr)
    call MPI_File_read_all(fh, lx, 1, MPI_INTEGER, status, ierr)
    call MPI_File_read_all(fh, optional_fields, 1, MPI_INTEGER, status, ierr)
    call MPI_File_read_all(fh, chkp%t, 1, MPI_DOUBLE_PRECISION, status, ierr)

    have_lag = mod(optional_fields,2)/1
    have_scalar = mod(optional_fields,4)/2
    have_dtlag = mod(optional_fields,8)/4
    have_abvel = mod(optional_fields,16)/8
    have_scalarlag = mod(optional_fields,32)/16
    have_ale = mod(optional_fields,64)/32
    have_payloads = mod(optional_fields,128)/64

    if ( ( glb_nelv .ne. msh%glb_nelv ) .or. &
         ( gdim .ne. msh%gdim) .or. &
         ( (have_lag .eq. 0) .and. (read_lag) ) .or. &
         ( (have_scalar .eq. 0) .and. (read_scalar) ) .or. &
         ( (have_ale .eq. 0) .and. (read_ale) ) ) then
       call neko_error('Checkpoint does not match case')
    end if
    nel = msh%nelv
    this%sim_Xh => u%Xh
    if (gdim .eq. 3) then
       call this%chkp_Xh%init(GLL, lx, lx, lx)
    else
       call this%chkp_Xh%init(GLL, lx, lx)
    end if
    if (this%mesh2mesh) then
       if (read_ale) then
          call neko_error('ALE does not yet support mesh2mesh ' // &
               'interpolation for restart!')
       end if
       call dof%init(msh, this%chkp_Xh)
       call this%global_interp%init(dof, NEKO_COMM, tol = tol)
       call this%global_interp%find_points(u%dof%x%x, u%dof%y%x, u%dof%z%x, &
            u%dof%size())
    else
       call this%space_interp%init(this%sim_Xh, this%chkp_Xh)
    end if
    dof_offset = int(msh%offset_el, i8) * int(this%chkp_Xh%lxyz, i8)
    n_glb_dofs = int(this%chkp_Xh%lxyz, i8) * int(msh%glb_nelv, i8)

    !
    ! Read mandatory checkpoint data
    !

    byte_offset = 4_i8 * int(MPI_INTEGER_SIZE, i8) + &
         int(MPI_DOUBLE_PRECISION_SIZE, i8)
    byte_offset = byte_offset + &
         dof_offset * int(MPI_REAL_PREC_SIZE, i8)
    call this%read_field(fh, byte_offset, u%x, nel)
    mpi_offset = 4_i8 * int(MPI_INTEGER_SIZE, i8) + &
         int(MPI_DOUBLE_PRECISION_SIZE, i8)
    mpi_offset = mpi_offset +&
         n_glb_dofs * int(MPI_REAL_PREC_SIZE, i8)

    byte_offset = mpi_offset + &
         dof_offset * int(MPI_REAL_PREC_SIZE, i8)
    call this%read_field(fh, byte_offset, v%x, nel)
    mpi_offset = mpi_offset + n_glb_dofs * int(MPI_REAL_PREC_SIZE, i8)

    byte_offset = mpi_offset + &
         dof_offset * int(MPI_REAL_PREC_SIZE, i8)
    call this%read_field(fh, byte_offset, w%x, nel)
    mpi_offset = mpi_offset + n_glb_dofs * int(MPI_REAL_PREC_SIZE, i8)

    byte_offset = mpi_offset + &
         dof_offset * int(MPI_REAL_PREC_SIZE, i8)
    call this%read_field(fh, byte_offset, p%x, nel)
    mpi_offset = mpi_offset + n_glb_dofs * int(MPI_REAL_PREC_SIZE, i8)

    !
    ! Read optional payload
    !

    if (read_lag) then
       do i = 1, ulag%size()
          byte_offset = mpi_offset + &
               dof_offset * int(MPI_REAL_PREC_SIZE, i8)
          call this%read_field(fh, byte_offset, ulag%lf(i)%x, nel)
          mpi_offset = mpi_offset + n_glb_dofs * int(MPI_REAL_PREC_SIZE, i8)
       end do

       do i = 1, vlag%size()
          byte_offset = mpi_offset + &
               dof_offset * int(MPI_REAL_PREC_SIZE, i8)
          call this%read_field(fh, byte_offset, vlag%lf(i)%x, nel)
          mpi_offset = mpi_offset + n_glb_dofs * int(MPI_REAL_PREC_SIZE, i8)
       end do

       do i = 1, wlag%size()
          byte_offset = mpi_offset + &
               dof_offset * int(MPI_REAL_PREC_SIZE, i8)
          call this%read_field(fh, byte_offset, wlag%lf(i)%x, nel)
          mpi_offset = mpi_offset + n_glb_dofs * int(MPI_REAL_PREC_SIZE, i8)
       end do
    end if

    if (read_scalar) then
       byte_offset = mpi_offset + &
            dof_offset * int(MPI_REAL_PREC_SIZE, i8)
       call this%read_field(fh, byte_offset, s%x, nel)
       mpi_offset = mpi_offset + n_glb_dofs * int(MPI_REAL_PREC_SIZE, i8)
    end if

    if (read_dtlag .and. have_dtlag .eq. 1) then
       call MPI_File_read_at_all(fh, mpi_offset, tlag, 10, &
            MPI_DOUBLE_PRECISION, status, ierr)
       mpi_offset = mpi_offset + 10_i8 * int(MPI_DOUBLE_PRECISION_SIZE, i8)
       call MPI_File_read_at_all(fh, mpi_offset, dtlag, 10, &
            MPI_DOUBLE_PRECISION, status, ierr)
       mpi_offset = mpi_offset + 10_i8 * int(MPI_DOUBLE_PRECISION_SIZE, i8)
    end if

    if (read_abvel .and. have_abvel .eq. 1) then
       byte_offset = mpi_offset + &
            dof_offset * int(MPI_REAL_PREC_SIZE, i8)
       call this%read_field(fh, byte_offset, abx1%x, nel)
       mpi_offset = mpi_offset + n_glb_dofs * int(MPI_REAL_PREC_SIZE, i8)
       byte_offset = mpi_offset + &
            dof_offset * int(MPI_REAL_PREC_SIZE, i8)
       call this%read_field(fh, byte_offset, abx2%x, nel)
       mpi_offset = mpi_offset + n_glb_dofs * int(MPI_REAL_PREC_SIZE, i8)
       byte_offset = mpi_offset + &
            dof_offset * int(MPI_REAL_PREC_SIZE, i8)
       call this%read_field(fh, byte_offset, aby1%x, nel)
       mpi_offset = mpi_offset + n_glb_dofs * int(MPI_REAL_PREC_SIZE, i8)
       byte_offset = mpi_offset + &
            dof_offset * int(MPI_REAL_PREC_SIZE, i8)
       call this%read_field(fh, byte_offset, aby2%x, nel)
       mpi_offset = mpi_offset + n_glb_dofs * int(MPI_REAL_PREC_SIZE, i8)
       byte_offset = mpi_offset + &
            dof_offset * int(MPI_REAL_PREC_SIZE, i8)
       call this%read_field(fh, byte_offset, abz1%x, nel)
       mpi_offset = mpi_offset + n_glb_dofs * int(MPI_REAL_PREC_SIZE, i8)
       byte_offset = mpi_offset + &
            dof_offset * int(MPI_REAL_PREC_SIZE, i8)
       call this%read_field(fh, byte_offset, abz2%x, nel)
       mpi_offset = mpi_offset + n_glb_dofs * int(MPI_REAL_PREC_SIZE, i8)
    end if
    if (read_scalarlag .and. have_scalarlag .eq. 1) then
       do i = 1, slag%size()
          byte_offset = mpi_offset + &
               dof_offset * int(MPI_REAL_PREC_SIZE, i8)
          call this%read_field(fh, byte_offset, slag%lf(i)%x, nel)
          mpi_offset = mpi_offset + n_glb_dofs * int(MPI_REAL_PREC_SIZE, i8)
       end do
       byte_offset = mpi_offset + &
            dof_offset * int(MPI_REAL_PREC_SIZE, i8)
       call this%read_field(fh, byte_offset, abs1%x, nel)
       mpi_offset = mpi_offset + n_glb_dofs * int(MPI_REAL_PREC_SIZE, i8)
       byte_offset = mpi_offset + &
            dof_offset * int(MPI_REAL_PREC_SIZE, i8)
       call this%read_field(fh, byte_offset, abs2%x, nel)
       mpi_offset = mpi_offset + n_glb_dofs * int(MPI_REAL_PREC_SIZE, i8)
    end if

    if (read_ale .and. have_ale .eq. 1) then

       byte_offset = mpi_offset + &
            dof_offset * int(MPI_REAL_PREC_SIZE, i8)
       call this%read_field(fh, byte_offset, msh_x, nel)
       mpi_offset = mpi_offset + n_glb_dofs * int(MPI_REAL_PREC_SIZE, i8)

       byte_offset = mpi_offset + &
            dof_offset * int(MPI_REAL_PREC_SIZE, i8)
       call this%read_field(fh, byte_offset, msh_y, nel)
       mpi_offset = mpi_offset + n_glb_dofs * int(MPI_REAL_PREC_SIZE, i8)

       byte_offset = mpi_offset + &
            dof_offset * int(MPI_REAL_PREC_SIZE, i8)
       call this%read_field(fh, byte_offset, msh_z, nel)
       mpi_offset = mpi_offset + n_glb_dofs * int(MPI_REAL_PREC_SIZE, i8)

       byte_offset = mpi_offset + &
            dof_offset * int(MPI_REAL_PREC_SIZE, i8)
       call this%read_field(fh, byte_offset, wm_x%x, nel)
       mpi_offset = mpi_offset + n_glb_dofs * int(MPI_REAL_PREC_SIZE, i8)

       byte_offset = mpi_offset + &
            dof_offset * int(MPI_REAL_PREC_SIZE, i8)
       call this%read_field(fh, byte_offset, wm_y%x, nel)
       mpi_offset = mpi_offset + n_glb_dofs * int(MPI_REAL_PREC_SIZE, i8)

       byte_offset = mpi_offset + &
            dof_offset * int(MPI_REAL_PREC_SIZE, i8)
       call this%read_field(fh, byte_offset, wm_z%x, nel)
       mpi_offset = mpi_offset + n_glb_dofs * int(MPI_REAL_PREC_SIZE, i8)

       do i = 1, wm_x_lag%size()
          byte_offset = mpi_offset + &
               dof_offset * int(MPI_REAL_PREC_SIZE, i8)
          call this%read_field(fh, byte_offset, wm_x_lag%lf(i)%x, nel)
          mpi_offset = mpi_offset + n_glb_dofs * int(MPI_REAL_PREC_SIZE, i8)
       end do
       do i = 1, wm_y_lag%size()
          byte_offset = mpi_offset + &
               dof_offset * int(MPI_REAL_PREC_SIZE, i8)
          call this%read_field(fh, byte_offset, wm_y_lag%lf(i)%x, nel)
          mpi_offset = mpi_offset + n_glb_dofs * int(MPI_REAL_PREC_SIZE, i8)
       end do
       do i = 1, wm_z_lag%size()
          byte_offset = mpi_offset + &
               dof_offset * int(MPI_REAL_PREC_SIZE, i8)
          call this%read_field(fh, byte_offset, wm_z_lag%lf(i)%x, nel)
          mpi_offset = mpi_offset + n_glb_dofs * int(MPI_REAL_PREC_SIZE, i8)
       end do
       byte_offset = mpi_offset + &
            dof_offset * int(MPI_REAL_PREC_SIZE, i8)
       call this%read_field(fh, byte_offset, Blag, nel)
       mpi_offset = mpi_offset + n_glb_dofs * int(MPI_REAL_PREC_SIZE, i8)

       byte_offset = mpi_offset + &
            dof_offset * int(MPI_REAL_PREC_SIZE, i8)
       call this%read_field(fh, byte_offset, Blaglag, nel)
       mpi_offset = mpi_offset + n_glb_dofs * int(MPI_REAL_PREC_SIZE, i8)

       call MPI_File_read_at_all(fh, mpi_offset, pivot_pos, &
            size(pivot_pos), MPI_REAL_PRECISION, status, ierr)
       mpi_offset = mpi_offset + int(size(pivot_pos), i8) &
            * int(MPI_REAL_PREC_SIZE, i8)

       call MPI_File_read_at_all(fh, mpi_offset, pivot_vel_lag, &
            size(pivot_vel_lag), MPI_REAL_PRECISION, status, ierr)
       mpi_offset = mpi_offset + int(size(pivot_vel_lag), i8) &
            * int(MPI_REAL_PREC_SIZE, i8)

       call MPI_File_read_at_all(fh, mpi_offset, basis_pos, &
            size(basis_pos), MPI_REAL_PRECISION, status, ierr)
       mpi_offset = mpi_offset + &
            int(size(basis_pos), i8) * int(MPI_REAL_PREC_SIZE, i8)

       call MPI_File_read_at_all(fh, mpi_offset, basis_vel_lag, &
            size(basis_vel_lag), MPI_REAL_PRECISION, status, ierr)
       mpi_offset = mpi_offset + &
            int(size(basis_vel_lag), i8) * int(MPI_REAL_PREC_SIZE, i8)
    end if

    call chkp_read_generic_payloads(fh, chkp, have_payloads .eq. 1)

    call MPI_File_close(fh, ierr)

    if (ierr .ne. MPI_SUCCESS) then
       call neko_error('Error reading checkpoint file ' // trim(fname))
    end if

    call this%global_interp%free()
    call this%space_interp%free()
    call dof%free()

  end subroutine chkp_file_read

  !> Read and optionally interpolate one legacy checkpoint field.
  !! @param fh Open MPI file handle.
  !! @param byte_offset Byte offset of the rank-local field data.
  !! @param x Simulation-space field storage to populate.
  !! @param nel Number of local elements.
  subroutine chkp_read_field(this, fh, byte_offset, x, nel)
    class(chkp_file_t) :: this
    type(MPI_File) :: fh
    integer(kind=MPI_OFFSET_KIND) :: byte_offset
    integer, intent(in) :: nel
    real(kind=rp) :: x(this%sim_Xh%lxyz, nel)
    real(kind=rp), allocatable :: read_array(:)
    integer :: nel_stride, frac_space
    type(MPI_Status) :: status
    integer :: ierr, i, n
    logical :: interp_on_host = .true.

    n = this%chkp_xh%lxyz*nel
    allocate(read_array(n))

    call rzero(read_array,n)
    call MPI_File_read_at_all(fh, byte_offset, read_array, &
         n, MPI_REAL_PRECISION, status, ierr)
    if (this%mesh2mesh) then
       x = 0.0_rp
       call this%global_interp%evaluate(x, read_array, interp_on_host)

    else if (this%sim_Xh%lx .ne. this%chkp_Xh%lx) then
       call this%space_interp%map_host(x, read_array, nel, this%sim_Xh)
    else
       do i = 1,n
          x(i,1) = read_array(i)
       end do
    end if
    deallocate(read_array)
  end subroutine chkp_read_field

  !> Build the fixed legacy view from format-independent payloads.
  !! @param chkp Format-independent checkpoint registrations.
  !! @param view Positional view required by the `.chkp` implementation.
  subroutine chkp_build_legacy_view(chkp, view)
    type(chkp_t), intent(in) :: chkp
    type(legacy_checkpoint_view_t), intent(out) :: view
    type(checkpoint_payload_t), pointer :: payload
    type(checkpoint_array_t), pointer :: array
    type(checkpoint_mesh_array_t), pointer :: mesh_array
    character(len=:), allocatable :: scalar_name
    integer :: i

    payload => chkp_find_payload(chkp, "fluid")
    if (.not. associated(payload)) then
       call neko_error("Checkpoint not initialized")
    end if
    view%u => payload%find_field("u")
    view%v => payload%find_field("v")
    view%w => payload%find_field("w")
    view%p => payload%find_field("p")
    view%abx1 => payload%find_field("abx1")
    view%abx2 => payload%find_field("abx2")
    view%aby1 => payload%find_field("aby1")
    view%aby2 => payload%find_field("aby2")
    view%abz1 => payload%find_field("abz1")
    view%abz2 => payload%find_field("abz2")
    view%ulag => payload%find_series("u")
    view%vlag => payload%find_series("v")
    view%wlag => payload%find_series("w")

    call chkp%get_time_history(view%tlag, view%dtlag)

    ! Only one scalar assumed, we grab the first one we find.
    do i = 1, chkp%payload_count()
       if (index(trim(chkp%payloads(i)%ptr%name), "scalars/") .eq. 1) then
          payload => chkp%payloads(i)%ptr
          scalar_name = payload%name(len("scalars/") + 1:)
          view%s => payload%find_field(scalar_name)
          view%slag => payload%find_series(scalar_name)
          view%abs1 => payload%find_field(trim(scalar_name) // "_abx1")
          view%abs2 => payload%find_field(trim(scalar_name) // "_abx2")
          exit
       end if
    end do

    payload => chkp_find_payload(chkp, "ale")
    if (.not. associated(payload)) return

    view%wm_x => payload%find_field("wm_x")
    view%wm_y => payload%find_field("wm_y")
    view%wm_z => payload%find_field("wm_z")
    view%wm_x_lag => payload%find_series("wm_x")
    view%wm_y_lag => payload%find_series("wm_y")
    view%wm_z_lag => payload%find_series("wm_z")

    mesh_array => payload%find_mesh_array("mesh_x")
    if (associated(mesh_array)) view%msh_x => mesh_array%x
    mesh_array => payload%find_mesh_array("mesh_y")
    if (associated(mesh_array)) view%msh_y => mesh_array%x
    mesh_array => payload%find_mesh_array("mesh_z")
    if (associated(mesh_array)) view%msh_z => mesh_array%x
    mesh_array => payload%find_mesh_array("B_lag")
    if (associated(mesh_array)) view%Blag => mesh_array%x
    mesh_array => payload%find_mesh_array("B_laglag")
    if (associated(mesh_array)) view%Blaglag => mesh_array%x
    array => payload%find_array("pivot_position")
    if (associated(array)) view%pivot_pos => array%x
    array => payload%find_array("pivot_velocity_lag")
    if (associated(array)) view%pivot_vel_lag => array%x
    array => payload%find_array("basis_position")
    if (associated(array)) view%basis_pos => array%x
    array => payload%find_array("basis_velocity_lag")
    if (associated(array)) view%basis_vel_lag => array%x

  end subroutine chkp_build_legacy_view

  !> Find a payload without making its absence an error.
  !! @param chkp Checkpoint to search.
  !! @param name Payload path to find.
  !! @return Pointer to the payload, or a disassociated pointer if absent.
  function chkp_find_payload(chkp, name) result(payload)
    type(chkp_t), intent(in) :: chkp
    character(len=*), intent(in) :: name
    type(checkpoint_payload_t), pointer :: payload
    integer :: i

    nullify(payload)
    do i = 1, chkp%payload_count()
       if (trim(chkp%payloads(i)%ptr%name) .eq. trim(name)) then
          payload => chkp%payloads(i)%ptr
          return
       end if
    end do

  end function chkp_find_payload

  !> Whether a payload is stored in the fixed positional `.chkp` sections.
  !! @param name Payload path.
  pure function chkp_is_legacy_payload(name) result(legacy)
    character(len=*), intent(in) :: name
    logical :: legacy

    legacy = trim(name) .eq. "fluid" .or. trim(name) .eq. "time" .or. &
         trim(name) .eq. "ale" .or. index(name, "scalars/") .eq. 1
  end function chkp_is_legacy_payload

  !> Number of payloads stored in the generic payload section.
  !! @param chkp Checkpoint registrations.
  function chkp_generic_payload_count(chkp) result(n)
    type(chkp_t), intent(in) :: chkp
    integer :: n, i

    n = 0
    do i = 1, chkp%payload_count()
       if (.not. chkp_is_legacy_payload(chkp%payloads(i)%ptr%name)) n = n + 1
    end do
  end function chkp_generic_payload_count

  !> Write the payloads that have no fixed `.chkp` section.
  !! @details The section starts at `mpi_offset` and is laid out as
  !! @verbatim
  !! int32 n_payloads
  !! per payload: int32 name_len, char name(name_len), int32 n_arrays
  !!   per array: int32 name_len, char name(name_len), int32 kind,
  !!              int64 n_global, real64 data(n_global)
  !! int64 section_start, int64 magic
  !! @endverbatim
  !! The footer lets readers locate the section from the end of the file,
  !! independently of which legacy sections they consumed. Resizable arrays
  !! are concatenated in rank order, replicated arrays are written by rank 0.
  !! Only replicated and resizable arrays are supported.
  !! @param fh Open MPI file handle.
  !! @param mpi_offset Byte offset of the section start, advanced past it.
  !! @param chkp Checkpoint registrations.
  subroutine chkp_write_generic_payloads(fh, mpi_offset, chkp)
    type(MPI_File), intent(inout) :: fh
    integer(kind=MPI_OFFSET_KIND), intent(inout) :: mpi_offset
    type(chkp_t), intent(in) :: chkp
    type(checkpoint_payload_t), pointer :: payload
    type(checkpoint_array_t), pointer :: array
    type(MPI_Status) :: status
    integer(kind=MPI_OFFSET_KIND) :: section_start
    integer(kind=i8) :: n_local, n_global, offset
    real(kind=dp), allocatable :: tmp(:)
    integer :: i, j, ierr, root_count

    root_count = merge(1, 0, pe_rank .eq. 0)
    section_start = mpi_offset

    call write_int(chkp_generic_payload_count(chkp))
    do i = 1, chkp%payload_count()
       payload => chkp%payloads(i)%ptr
       if (chkp_is_legacy_payload(payload%name)) cycle

       if (payload%field_count() .gt. 0 .or. &
            payload%series_count() .gt. 0 .or. &
            payload%mesh_array_count() .gt. 0) then
          call neko_error("Checkpoint payload '" // trim(payload%name) // &
               "' holds fields, which the .chkp format cannot store; " // &
               "use HDF5 checkpointing")
       end if

       call write_str(payload%name)
       call write_int(payload%array_count())
       do j = 1, payload%array_count()
          array => payload%arrays(j)%ptr
          call write_str(array%name)

          if (array%resizable) then
             if (allocated(array%buf)) then
                tmp = array%buf
             else
                allocate(tmp(0))
             end if
             n_local = int(size(tmp), i8)
             offset = 0_i8
             call MPI_Exscan(n_local, offset, 1, MPI_INTEGER8, MPI_SUM, &
                  NEKO_COMM, ierr)
             if (pe_rank .eq. 0) offset = 0_i8
             call MPI_Allreduce(n_local, n_global, 1, MPI_INTEGER8, &
                  MPI_SUM, NEKO_COMM, ierr)

             call write_int(CHKP_ARRAY_RESIZABLE)
             call write_i8(n_global)
             call MPI_File_write_at_all(fh, &
                  mpi_offset + offset * CHKP_WORD_BYTES, tmp, &
                  int(n_local), MPI_DOUBLE_PRECISION, status, ierr)
             deallocate(tmp)

          else if (array%replicated) then
             if (associated(array%x_dp)) then
                tmp = array%x_dp
             else
                tmp = real(array%x, dp)
             end if
             n_global = int(size(tmp), i8)

             call write_int(CHKP_ARRAY_REPLICATED)
             call write_i8(n_global)
             call MPI_File_write_at_all(fh, mpi_offset, tmp, &
                  root_count * size(tmp), MPI_DOUBLE_PRECISION, status, ierr)
             deallocate(tmp)

          else
             call neko_error("Checkpoint array '" // trim(payload%name) // &
                  "/" // trim(array%name) // "' is distributed with a " // &
                  "fixed layout, which the .chkp format cannot store; " // &
                  "use HDF5 checkpointing")
          end if
          mpi_offset = mpi_offset + n_global * CHKP_WORD_BYTES
       end do
    end do

    call write_i8(int(section_start, i8))
    call write_i8(CHKP_PAYLOADS_MAGIC)

  contains

    !> Write one default integer from rank 0.
    subroutine write_int(val)
      integer, intent(in) :: val
      integer :: buf(1)

      buf(1) = val
      call MPI_File_write_at_all(fh, mpi_offset, buf, root_count, &
           MPI_INTEGER, status, ierr)
      mpi_offset = mpi_offset + int(MPI_INTEGER_SIZE, i8)
    end subroutine write_int

    !> Write one 8-byte integer from rank 0.
    subroutine write_i8(val)
      integer(kind=i8), intent(in) :: val
      integer(kind=i8) :: buf(1)

      buf(1) = val
      call MPI_File_write_at_all(fh, mpi_offset, buf, root_count, &
           MPI_INTEGER8, status, ierr)
      mpi_offset = mpi_offset + CHKP_WORD_BYTES
    end subroutine write_i8

    !> Write a length-prefixed string from rank 0.
    subroutine write_str(str)
      character(len=*), intent(in) :: str
      character(len=len_trim(str)) :: buf

      buf = trim(str)
      call write_int(len(buf))
      call MPI_File_write_at_all(fh, mpi_offset, buf, &
           root_count * len(buf), MPI_CHARACTER, status, ierr)
      mpi_offset = mpi_offset + int(len(buf), i8)
    end subroutine write_str

  end subroutine chkp_write_generic_payloads

  !> Read the generic payload section and mark the payloads found there as
  !! loaded. Payloads in the file that are not registered, and arrays a
  !! registered payload does not declare, are skipped. Resizable arrays are
  !! distributed evenly over the ranks, as for HDF5 checkpoints.
  !! @param fh Open MPI file handle.
  !! @param chkp Checkpoint registrations to populate.
  !! @param have_payloads Whether the header flags a generic payload section.
  subroutine chkp_read_generic_payloads(fh, chkp, have_payloads)
    type(MPI_File), intent(inout) :: fh
    type(chkp_t), intent(inout) :: chkp
    logical, intent(in) :: have_payloads
    type(checkpoint_payload_t), pointer :: payload
    type(checkpoint_array_t), pointer :: array
    type(MPI_Status) :: status
    integer(kind=MPI_OFFSET_KIND) :: pos, file_size
    integer(kind=i8) :: footer(2), n_global, n_local, offset, remainder
    character(len=:), allocatable :: payload_name, array_name
    logical, allocatable :: found(:)
    real(kind=dp), allocatable :: tmp(:)
    integer :: n_payloads, n_arrays, array_kind, i, j, k, ierr

    do i = 1, chkp%payload_count()
       if (.not. chkp_is_legacy_payload(chkp%payloads(i)%ptr%name)) then
          chkp%payloads(i)%ptr%loaded = .false.
       end if
    end do

    if (have_payloads) then
       call MPI_File_get_size(fh, file_size, ierr)
       call MPI_File_read_at_all(fh, file_size - 2_i8 * CHKP_WORD_BYTES, &
            footer, 2, MPI_INTEGER8, status, ierr)
       if (footer(2) .ne. CHKP_PAYLOADS_MAGIC) then
          call neko_error("Corrupt payload section in .chkp file")
       end if
       pos = footer(1)

       call read_int(n_payloads)
       do k = 1, n_payloads
          call read_str(payload_name)
          payload => chkp_find_payload(chkp, payload_name)
          if (associated(payload)) then
             if (allocated(found)) deallocate(found)
             allocate(found(payload%array_count()))
             found = .false.
          end if

          call read_int(n_arrays)
          do j = 1, n_arrays
             call read_str(array_name)
             call read_int(array_kind)
             call read_i8(n_global)

             nullify(array)
             if (associated(payload)) then
                do i = 1, payload%array_count()
                   if (payload%arrays(i)%ptr%name .eq. array_name) then
                      array => payload%arrays(i)%ptr
                      found(i) = .true.
                      exit
                   end if
                end do
             end if

             if (associated(array)) then
                if (array%resizable .neqv. &
                     (array_kind .eq. CHKP_ARRAY_RESIZABLE)) then
                   call neko_error("Checkpoint array '" // payload_name // &
                        "/" // array_name // "' has an incompatible layout")
                end if

                if (array%resizable) then
                   n_local = n_global / int(pe_size, i8)
                   remainder = mod(n_global, int(pe_size, i8))
                   offset = int(pe_rank, i8) * n_local + &
                        min(int(pe_rank, i8), remainder)
                   if (int(pe_rank, i8) .lt. remainder) then
                      n_local = n_local + 1_i8
                   end if
                   array%global_count = n_global
                   array%offset = offset
                   if (allocated(array%buf)) deallocate(array%buf)
                   allocate(array%buf(n_local))
                   call MPI_File_read_at_all(fh, &
                        pos + offset * CHKP_WORD_BYTES, array%buf, &
                        int(n_local), MPI_DOUBLE_PRECISION, status, ierr)
                else
                   if (associated(array%x_dp)) then
                      n_local = int(size(array%x_dp), i8)
                   else
                      n_local = int(size(array%x), i8)
                   end if
                   if (n_local .ne. n_global) then
                      call neko_error("Checkpoint array '" // payload_name &
                           // "/" // array_name // &
                           "' has an incompatible global extent")
                   end if
                   allocate(tmp(n_global))
                   call MPI_File_read_at_all(fh, pos, tmp, int(n_global), &
                        MPI_DOUBLE_PRECISION, status, ierr)
                   if (associated(array%x_dp)) then
                      array%x_dp = tmp
                   else
                      array%x = real(tmp, rp)
                   end if
                   deallocate(tmp)
                end if
             end if

             pos = pos + n_global * CHKP_WORD_BYTES
          end do

          if (associated(payload)) then
             do i = 1, payload%array_count()
                if (.not. found(i)) then
                   call neko_error("Checkpoint array '" // payload_name // &
                        "/" // trim(payload%arrays(i)%ptr%name) // &
                        "' is missing")
                end if
             end do
             payload%loaded = .true.
          end if
       end do
    end if

    do i = 1, chkp%payload_count()
       payload => chkp%payloads(i)%ptr
       if (chkp_is_legacy_payload(payload%name)) cycle
       if (.not. payload%loaded .and. .not. payload%is_optional) then
          call neko_error("Checkpoint payload '" // trim(payload%name) // &
               "' is missing")
       end if
    end do

  contains

    !> Read one default integer.
    subroutine read_int(val)
      integer, intent(out) :: val
      integer :: buf(1)

      call MPI_File_read_at_all(fh, pos, buf, 1, MPI_INTEGER, status, ierr)
      val = buf(1)
      pos = pos + int(MPI_INTEGER_SIZE, i8)
    end subroutine read_int

    !> Read one 8-byte integer.
    subroutine read_i8(val)
      integer(kind=i8), intent(out) :: val
      integer(kind=i8) :: buf(1)

      call MPI_File_read_at_all(fh, pos, buf, 1, MPI_INTEGER8, status, ierr)
      val = buf(1)
      pos = pos + CHKP_WORD_BYTES
    end subroutine read_i8

    !> Read a length-prefixed string.
    subroutine read_str(str)
      character(len=:), allocatable, intent(out) :: str
      integer :: n

      call read_int(n)
      allocate(character(len=n) :: str)
      if (n .gt. 0) then
         call MPI_File_read_at_all(fh, pos, str, n, MPI_CHARACTER, &
              status, ierr)
      end if
      pos = pos + int(n, i8)
    end subroutine read_str

  end subroutine chkp_read_generic_payloads

  !> Set whether checkpoint files may be overwritten.
  !! @param overwrite Whether existing files may be overwritten.
  subroutine chkp_file_set_overwrite(this, overwrite)
    class(chkp_file_t), intent(inout) :: this
    logical, intent(in) :: overwrite
    this%overwrite = overwrite
  end subroutine chkp_file_set_overwrite
end module chkp_file
