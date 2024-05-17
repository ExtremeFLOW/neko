! Copyright (c) 2024, The Neko Authors
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
!> HDF5 file format
module hdf5_file
  use num_types, only : rp
  use generic_file, only : generic_file_t
  use checkpoint, only : chkp_t
  use utils, only : neko_error, filename_suffix_pos
  use mesh, only : mesh_t
  use field, only : field_t, field_ptr_t
  use field_list, only : field_list_t
  use field_series, only : field_series_t, field_series_ptr_t
  use dofmap, only : dofmap_t
  use logger
  use comm
  use mpi_f08
#ifdef HAVE_HDF5
  use hdf5
#endif
  implicit none
  private

  !> Interface for HDF5 files
  type, public, extends(generic_file_t) :: hdf5_file_t
   contains
     procedure :: read => hdf5_file_read
     procedure :: write => hdf5_file_write
  end type hdf5_file_t

contains

#ifdef HAVE_HDF5
  
  !> Write data in HDF5 format
  subroutine hdf5_file_write(this, data, t)
    class(hdf5_file_t), intent(inout) :: this
    class(*), target, intent(in) :: data
    real(kind=rp), intent(in), optional :: t
    type(mesh_t), pointer :: msh
    type(dofmap_t), pointer :: dof
    type(field_ptr_t), allocatable :: fp(:)
    type(field_series_ptr_t), allocatable :: fsp(:)
    real(kind=rp), pointer :: dtlag(:)
    real(kind=rp), pointer :: tlag(:) 
    integer :: ierr, info, drank, i, j
    integer(hid_t) :: plist_id, file_id, dset_id, grp_id, attr_id
    integer(hid_t) :: filespace, memspace
    integer(hsize_t), dimension(1) :: ddim, dcount, doffset
    integer :: suffix_pos
    character(len=5) :: id_str
    character(len=1024) :: fname
 
    call hdf5_file_determine_data(data, msh, dof, fp, fsp, dtlag, tlag)

    suffix_pos = filename_suffix_pos(this%fname)
    write(id_str, '(i5.5)') this%counter
    fname = trim(this%fname(1:suffix_pos-1))//id_str//'.h5'
       
    call h5open_f(ierr)
    call h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, ierr)
    info = MPI_INFO_NULL%mpi_val
    call h5pset_fapl_mpio_f(plist_id, NEKO_COMM%mpi_val, info, ierr)

    call h5fcreate_f(fname, H5F_ACC_TRUNC_F, &
         file_id, ierr, access_prp = plist_id)
    
    call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, ierr)
    call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, ierr)
    
    call h5screate_f(H5S_SCALAR_F, filespace, ierr)
    ddim = 1

    if (present(t)) then
       call h5acreate_f(file_id, "Time", H5T_NATIVE_DOUBLE, filespace, attr_id, &
                        ierr, h5p_default_f, h5p_default_f)
       call h5awrite_f(attr_id, H5T_NATIVE_DOUBLE, t, ddim, ierr)
       call h5aclose_f(attr_id, ierr)
    end if

    if (associated(dof)) then
       call h5acreate_f(file_id, "Lx", H5T_NATIVE_INTEGER, filespace, attr_id, &
                        ierr, h5p_default_f, h5p_default_f)
       call h5awrite_f(attr_id, H5T_NATIVE_INTEGER, dof%Xh%lx, ddim, ierr)
       call h5aclose_f(attr_id, ierr)
    end if
    
    if (associated(msh)) then
       call h5gcreate_f(file_id, "Mesh", grp_id, ierr, lcpl_id=h5p_default_f, &
            gcpl_id=h5p_default_f, gapl_id=h5p_default_f)
    
       call h5acreate_f(grp_id, "Elements", H5T_NATIVE_INTEGER, filespace, attr_id, &
                        ierr, h5p_default_f, h5p_default_f)
       call h5awrite_f(attr_id, H5T_NATIVE_INTEGER, msh%glb_nelv, ddim, ierr)
       call h5aclose_f(attr_id, ierr)

       call h5acreate_f(grp_id, "Dimension", H5T_NATIVE_INTEGER, filespace, attr_id, &
                        ierr, h5p_default_f, h5p_default_f)
       call h5awrite_f(attr_id, H5T_NATIVE_INTEGER, msh%gdim, ddim, ierr)
       call h5aclose_f(attr_id, ierr)

       call h5gclose_f(grp_id, ierr)
    end if

    
    call h5sclose_f(filespace, ierr)
    
    !
    ! Write restart group (tlag, dtlag)
    !
    if (associated(tlag) .and. associated(dtlag)) then    
       call h5gcreate_f(file_id, "Restart", grp_id, ierr, lcpl_id=h5p_default_f, &
            gcpl_id=h5p_default_f, gapl_id=h5p_default_f)

       drank = 1
       ddim = size(tlag)
       doffset(1) = 0
       if (pe_rank .eq. 0) then
          dcount = size(tlag)
       else
          dcount = 0
       end if

       call h5screate_simple_f(drank, ddim, filespace, ierr)

       call h5dcreate_f(grp_id,'tlag', H5T_NATIVE_DOUBLE, &
                        filespace, dset_id, ierr)
       call h5dget_space_f(dset_id, filespace, ierr)
       call h5sselect_hyperslab_f (filespace, H5S_SELECT_SET_F, &
                                   doffset, dcount, ierr)
       call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, tlag, &
                          ddim, ierr, xfer_prp = plist_id)
       call h5dclose_f(dset_id, ierr)

       call h5dcreate_f(grp_id,'dtlag', H5T_NATIVE_DOUBLE, &
                        filespace, dset_id, ierr)
       call h5dget_space_f(dset_id, filespace, ierr)
       call h5sselect_hyperslab_f (filespace, H5S_SELECT_SET_F, &
                                   doffset, dcount, ierr)
       call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, dtlag, &
                          ddim, ierr, xfer_prp = plist_id)
       call h5dclose_f(dset_id, ierr)
       
       call h5sclose_f(filespace, ierr)
       call h5gclose_f(grp_id, ierr)
       
    end if


    !
    ! Write fields group
    ! 
    if (allocated(fp) .or. allocated(fsp)) then
       call h5gcreate_f(file_id, "Fields", grp_id, ierr, lcpl_id=h5p_default_f, &
            gcpl_id=h5p_default_f, gapl_id=h5p_default_f)
    
       dcount(1) = int(dof%size(), 8)
       doffset(1) = int(msh%offset_el, 8) * int((dof%Xh%lx**3),8)
       ddim =  int(dof%size(), 8)
       drank = 1
       call MPI_Allreduce(MPI_IN_PLACE, ddim(1), 1, &
                          MPI_INTEGER8, MPI_SUM, NEKO_COMM, ierr)

       call h5screate_simple_f(drank, ddim, filespace, ierr)
       call h5screate_simple_f(drank, dcount, memspace, ierr)
       

       if (allocated(fp)) then
          do i = 1, size(fp)
             call h5dcreate_f(grp_id, fp(i)%ptr%name, H5T_NATIVE_DOUBLE, &
                              filespace, dset_id, ierr)
             call h5dget_space_f(dset_id, filespace, ierr)
             call h5sselect_hyperslab_f(filespace, H5S_SELECT_SET_F, &
                                        doffset, dcount, ierr)
             call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, &
                             fp(i)%ptr%x(1,1,1,1), &
                             ddim, ierr, file_space_id = filespace, &
                             mem_space_id = memspace, xfer_prp = plist_id)
             call h5dclose_f(dset_id, ierr)
          end do
          deallocate(fp)       
       end if

       if (allocated(fsp)) then
          do i = 1, size(fsp)
             do j = 1, fsp(i)%ptr%size()
                call h5dcreate_f(grp_id, fsp(i)%ptr%lf(j)%name, &
                                 H5T_NATIVE_DOUBLE, filespace, dset_id, ierr)
                call h5dget_space_f(dset_id, filespace, ierr)
                call h5sselect_hyperslab_f(filespace, H5S_SELECT_SET_F, &
                                           doffset, dcount, ierr)
                call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, &
                                fsp(i)%ptr%lf(j)%x(1,1,1,1), &
                                ddim, ierr, file_space_id = filespace, &
                                mem_space_id = memspace, xfer_prp = plist_id)
                call h5dclose_f(dset_id, ierr)
             end do
          end do
          deallocate(fsp)
       end if
    
       call h5gclose_f(grp_id, ierr)
       call h5sclose_f(filespace, ierr)
       call h5sclose_f(memspace, ierr)
    end if
    
    call h5pclose_f(plist_id, ierr)
    call h5fclose_f(file_id, ierr)

    call h5close_f(ierr)

    this%counter = this%counter + 1
    
  end subroutine hdf5_file_write

  !> Read data in HDF5 format
  subroutine hdf5_file_read(this, data)
    class(hdf5_file_t) :: this
    class(*), target, intent(inout) :: data
    integer(hid_t) :: plist_id, file_id, dset_id, grp_id, attr_id
    integer(hid_t) :: filespace, memspace
    integer(hsize_t), dimension(1) :: ddim, dcount, doffset
    integer :: i,j, ierr, info, glb_nelv, gdim, lx, drank
    type(mesh_t), pointer :: msh
    type(dofmap_t), pointer :: dof
    type(field_ptr_t), allocatable :: fp(:)
    type(field_series_ptr_t), allocatable :: fsp(:)
    real(kind=rp), pointer :: dtlag(:)
    real(kind=rp), pointer :: tlag(:) 
    real(kind=rp) :: t

    call hdf5_file_determine_data(data, msh, dof, fp, fsp, dtlag, tlag)

    call h5open_f(ierr)
    call h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, ierr)
    info = MPI_INFO_NULL%mpi_val
    call h5pset_fapl_mpio_f(plist_id, NEKO_COMM%mpi_val, info, ierr)

    call h5fopen_f(trim(this%fname), H5F_ACC_RDONLY_F, &
         file_id, ierr, access_prp = plist_id)

    call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, ierr)
    call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, ierr)

    ddim = 1 
    call h5aopen_name_f(file_id, 'Time', attr_id, ierr)
    call h5aread_f(attr_id, H5T_NATIVE_DOUBLE, t, ddim, ierr)
    call h5aclose_f(attr_id, ierr)

    select type(data)
    type is(chkp_t)
       data%t = t
    end select

    call h5aopen_name_f(file_id, 'Lx', attr_id, ierr)
    call h5aread_f(attr_id, H5T_NATIVE_INTEGER, lx, ddim, ierr)
    call h5aclose_f(attr_id, ierr)

    call h5gopen_f(file_id, 'Mesh', grp_id, ierr, gapl_id=h5p_default_f)

    call h5aopen_name_f(grp_id, 'Elements', attr_id, ierr)
    call h5aread_f(attr_id, H5T_NATIVE_INTEGER, glb_nelv, ddim, ierr)
    call h5aclose_f(attr_id, ierr)

    call h5aopen_name_f(grp_id, 'Dimension', attr_id, ierr)
    call h5aread_f(attr_id, H5T_NATIVE_INTEGER, gdim, ddim, ierr)
    call h5aclose_f(attr_id, ierr)
    call h5gclose_f(grp_id, ierr)


    if (associated(tlag) .and. associated(dtlag)) then
       drank = 1
       ddim = size(tlag)
       doffset(1) = 0
       if (pe_rank .eq. 0) then
          dcount = size(tlag)
       else
          dcount = 0
       end if
       
       call h5gopen_f(file_id, 'Restart', grp_id, ierr, gapl_id=h5p_default_f)
       call h5dopen_f(grp_id, 'tlag', dset_id, ierr)
       call h5dget_space_f(dset_id, filespace, ierr)
       call h5sselect_hyperslab_f (filespace, H5S_SELECT_SET_F, &
                                   doffset, dcount, ierr)
       call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, tlag, ddim, ierr, xfer_prp=plist_id)
       call h5dclose_f(dset_id, ierr)
       call h5sclose_f(filespace, ierr)

       call h5dopen_f(grp_id, 'dtlag', dset_id, ierr)
       call h5dget_space_f(dset_id, filespace, ierr)
       call h5sselect_hyperslab_f (filespace, H5S_SELECT_SET_F, &
                                   doffset, dcount, ierr)
       call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, dtlag, ddim, ierr, xfer_prp=plist_id)
       call h5dclose_f(dset_id, ierr)
       call h5sclose_f(filespace, ierr)
       
       call h5gclose_f(grp_id, ierr)
    end if

   if (allocated(fp) .or. allocated(fsp)) then
      call h5gopen_f(file_id, 'Fields', grp_id, ierr, gapl_id=h5p_default_f)

       dcount(1) = int(dof%size(), 8)
       doffset(1) = int(msh%offset_el, 8) * int((dof%Xh%lx**3),8)
       ddim =  int(dof%size(), 8)
       drank = 1

       dcount(1) = int(dof%size(), 8)
       doffset(1) = int(msh%offset_el, 8) * int((dof%Xh%lx**3),8)
       ddim =  int(dof%size(), 8)
       drank = 1
       call MPI_Allreduce(MPI_IN_PLACE, ddim(1), 1, &
                          MPI_INTEGER8, MPI_SUM, NEKO_COMM, ierr)

       call h5screate_simple_f(drank, dcount, memspace, ierr)
       
      if (allocated(fp)) then
          do i = 1, size(fp)
             call h5dopen_f(grp_id, fp(i)%ptr%name, dset_id, ierr)
             call h5dget_space_f(dset_id, filespace, ierr)
             call h5sselect_hyperslab_f (filespace, H5S_SELECT_SET_F, &
                  doffset, dcount, ierr)
             call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, &
                            fp(i)%ptr%x(1,1,1,1), &
                            ddim, ierr, file_space_id = filespace, &
                            mem_space_id = memspace, xfer_prp=plist_id)
             call h5dclose_f(dset_id, ierr)
             call h5sclose_f(filespace, ierr)
          end do
       end if

       if (allocated(fsp)) then
          do i = 1, size(fsp)
             do j = 1, fsp(i)%ptr%size()
                call h5dopen_f(grp_id, fsp(i)%ptr%lf(j)%name, dset_id, ierr)
                call h5dget_space_f(dset_id, filespace, ierr)
                call h5sselect_hyperslab_f (filespace, H5S_SELECT_SET_F, &
                                            doffset, dcount, ierr)
                call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, &
                               fsp(i)%ptr%lf(j)%x(1,1,1,1), &
                               ddim, ierr, file_space_id = filespace, &
                               mem_space_id = memspace, xfer_prp=plist_id)
                call h5dclose_f(dset_id, ierr)
                call h5sclose_f(filespace, ierr)
             end do
          end do
       end if
       call h5sclose_f(memspace, ierr)
       call h5gclose_f(grp_id, ierr)
    end if
   
    call h5pclose_f(plist_id, ierr)
    call h5fclose_f(file_id, ierr)

    call h5close_f(ierr)

  end subroutine hdf5_file_read


  subroutine hdf5_file_determine_data(data, msh, dof, fp, fsp, dtlag, tlag)
    class(*), target, intent(in) :: data
    type(mesh_t), pointer, intent(inout) :: msh
    type(dofmap_t), pointer, intent(inout) :: dof
    type(field_ptr_t), allocatable, intent(inout) :: fp(:)
    type(field_series_ptr_t), allocatable, intent(inout) :: fsp(:)
    real(kind=rp), pointer, intent(inout) :: dtlag(:)
    real(kind=rp), pointer, intent(inout) :: tlag(:)
    integer :: i, j, fp_size, fp_cur, fsp_size, fsp_cur

    select type(data)
    type is (field_t)
       dof => data%dof
       msh => data%msh
       fp_size = 1
       allocate(fp(fp_size))
       fp(1)%ptr => data

       nullify(dtlag)
       nullify(tlag)
       
    type is (field_list_t)

       if (data%size() .gt. 0) then
          allocate(fp(data%size()))
       
          dof => data%dof(1)
          msh => data%msh(1)

          do i = 1, data%size()
             fp(i)%ptr => data%items(i)%ptr
          end do
       else
          call neko_error('Empty field list')
       end if

       nullify(dtlag)
       nullify(tlag)
       
    type is(chkp_t)

       if ( .not. associated(data%u) .or. &
            .not. associated(data%v) .or. &
            .not. associated(data%w) .or. &
            .not. associated(data%p) ) then
          call neko_error('Checkpoint not initialized')
       end if

       fp_size = 4

       if (associated(data%s)) then
          fp_size = fp_size + 1
       end if

       if (associated(data%abx1)) then
          fp_size = fp_size + 6
       end if

       if (associated(data%abs1)) then
          fp_size = fp_size + 2
       end if
       
       allocate(fp(fp_size))

       fsp_size = 0
       if (associated(data%ulag)) then
          fsp_size = fsp_size + 3
       end if

       if (associated(data%slag)) then
          fsp_size = fsp_size + 1
       end if

       if (fsp_size .gt. 0) then
          allocate(fsp(fsp_size))
          fsp_cur = 1
       end if

       dof => data%u%dof
       msh => data%u%msh
       
       fp(1)%ptr => data%u
       fp(2)%ptr => data%v
       fp(3)%ptr => data%w
       fp(4)%ptr => data%p

       fp_cur = 5
       if (associated(data%s)) then
          fp(fp_cur)%ptr => data%s
          fp_cur = fp_cur + 1
       end if

       if (associated(data%abx1)) then
          fp(fp_cur)%ptr => data%abx1
          fp(fp_cur+1)%ptr => data%abx2
          fp(fp_cur+2)%ptr => data%aby1
          fp(fp_cur+3)%ptr => data%aby2
          fp(fp_cur+4)%ptr => data%abz1
          fp(fp_cur+5)%ptr => data%abz2
          fp_cur = fp_cur + 6
       end if

       if (associated(data%abs1)) then
          fp(fp_cur)%ptr => data%abs1
          fp(fp_cur+1)%ptr => data%abs2
          fp_cur = fp_cur + 2
       end if

       if (associated(data%ulag)) then
          fsp(fsp_cur)%ptr => data%ulag
          fsp(fsp_cur+1)%ptr => data%vlag
          fsp(fsp_cur+2)%ptr => data%wlag
          fsp_cur = fsp_cur + 3
       end if

       if (associated(data%slag)) then
          fsp(fsp_cur)%ptr => data%slag
          fsp_cur = fsp_cur + 1
       end if

       if (associated(data%tlag)) then
          tlag => data%tlag
          dtlag => data%dtlag
       end if       
       
    class default
       call neko_log%error('Invalid data')
    end select
    
  end subroutine hdf5_file_determine_data
  
#else

  !> Write data in HDF5 format
  subroutine hdf5_file_write(this, data, t)
    class(hdf5_file_t), intent(inout) :: this
    class(*), target, intent(in) :: data
    real(kind=rp), intent(in), optional :: t
    call neko_error('Neko needs to be built with HDF5 support')
  end subroutine hdf5_file_write

  !> Read data in HDF5 format
  subroutine hdf5_file_read(this, data)
    class(hdf5_file_t) :: this
    class(*), target, intent(inout) :: data
    call neko_error('Neko needs to be built with HDF5 support')
  end subroutine hdf5_file_read


#endif

end module hdf5_file
