! Copyright (c) 2020-2021, The Neko Authors
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
!> Defines a simulation case
module case
  use num_types
  use fluid_fctry
  use fluid_output
  use scalar_output
  use chkp_output
  use mean_sqr_flow_output
  use mean_flow_output
  use parameters
  use mpi_types
  use mesh_field
  use parmetis
  use redist
  use sampler
  use flow_ic
  use stats
  use file
  use utils
  use mesh
  use comm
  use ext_bdf_scheme
  use logger
  use jobctrl
  use amr
  use user_intf
  use scalar_pnpn ! todo directly load the pnpn? can we have other
  implicit none

  type :: case_t
     type(mesh_t) :: msh
     type(param_t) :: params
     type(ext_bdf_scheme_t) :: ext_bdf
     real(kind=rp), dimension(10) :: tlag
     real(kind=rp), dimension(10) :: dtlag
     type(sampler_t) :: s
     type(fluid_output_t) :: f_out
     type(scalar_output_t) :: s_out
     type(chkp_output_t) :: f_chkp
     type(mean_flow_output_t) :: f_mf
     type(mean_sqr_flow_output_t) :: f_msqrf
     type(stats_t) :: q
     type(user_t) :: usr
     class(fluid_scheme_t), allocatable :: fluid
     type(scalar_pnpn_t), allocatable :: scalar
     type(amr_t) :: amr
  end type case_t

contains

  !> Initialize a case from an input file @a case_file
  subroutine case_init(C, case_file)
    type(case_t), target, intent(inout) :: C
    character(len=*), intent(in) :: case_file

    ! Namelist for case description
    character(len=NEKO_FNAME_LEN) :: mesh_file = ''
    character(len=80) :: fluid_scheme  = ''
    character(len=80) :: source_term = ''
    character(len=80) :: initial_condition = ''
    integer :: lx = 0
    logical :: scalar = .false.
    character(len=80) :: scalar_source_term = ''
    logical :: amr = .false.
    type(param_io_t) :: params
    namelist /NEKO_CASE/ mesh_file, fluid_scheme, lx,  &
         source_term, initial_condition, scalar, scalar_source_term,  &
         amr

    integer :: ierr
    type(file_t) :: msh_file, bdry_file, part_file
    type(mesh_fld_t) :: msh_part
    integer, parameter :: nbytes = NEKO_FNAME_LEN + (4 * 80) + 4 + 4 + 4
    character buffer(nbytes)
    integer :: pack_index
    type(mesh_fld_t) :: parts
   
    call neko_log%section('Case')
    call neko_log%message('Reading case file ' // trim(case_file))
    
    !
    ! Read case description
    !
    if (pe_rank .eq. 0) then
       open(10, file=trim(case_file))
       read(10, nml=NEKO_CASE)
       read(10, *) params
       close(10)
       
       pack_index = 0
       call MPI_Pack(mesh_file, NEKO_FNAME_LEN, MPI_CHARACTER, &
            buffer, nbytes, pack_index, NEKO_COMM, ierr)
       call MPI_Pack(fluid_scheme, 80, MPI_CHARACTER, &
            buffer, nbytes, pack_index, NEKO_COMM, ierr)
       call MPI_Pack(source_term, 80, MPI_CHARACTER, &
            buffer, nbytes, pack_index, NEKO_COMM, ierr)
       call MPI_Pack(initial_condition, 80, MPI_CHARACTER, &
            buffer, nbytes, pack_index, NEKO_COMM, ierr)
       call MPI_Pack(scalar_source_term, 80, MPI_CHARACTER, &
            buffer, nbytes, pack_index, NEKO_COMM, ierr)
       call MPI_Pack(lx, 1, MPI_INTEGER, &
            buffer, nbytes, pack_index, NEKO_COMM, ierr)
       call MPI_Pack(scalar, 1, MPI_LOGICAL, &
            buffer, nbytes, pack_index, NEKO_COMM, ierr)
       call MPI_Pack(amr, 1, MPI_LOGICAL, &
            buffer, nbytes, pack_index, NEKO_COMM, ierr)
       call MPI_Bcast(buffer, nbytes, MPI_PACKED, 0, NEKO_COMM, ierr)
       call MPI_Bcast(params%p, 1, MPI_NEKO_PARAMS, 0, NEKO_COMM, ierr)
    else
       call MPI_Bcast(buffer, nbytes, MPI_PACKED, 0, NEKO_COMM, ierr)
       pack_index = 0

       call MPI_Unpack(buffer, nbytes, pack_index, &
            mesh_file, NEKO_FNAME_LEN, MPI_CHARACTER, NEKO_COMM, ierr)
       call MPI_Unpack(buffer, nbytes, pack_index, &
            fluid_scheme, 80, MPI_CHARACTER, NEKO_COMM, ierr)
       call MPI_Unpack(buffer, nbytes, pack_index, &
            source_term, 80, MPI_CHARACTER, NEKO_COMM, ierr)
       call MPI_Unpack(buffer, nbytes, pack_index, &
            initial_condition, 80, MPI_CHARACTER, NEKO_COMM, ierr)
       call MPI_Unpack(buffer, nbytes, pack_index, &
            scalar_source_term, 80, MPI_CHARACTER, NEKO_COMM, ierr)
       call MPI_Unpack(buffer, nbytes, pack_index, &
            lx, 1, MPI_INTEGER, NEKO_COMM, ierr)
       call MPI_Unpack(buffer, nbytes, pack_index, &
            scalar, 1, MPI_LOGICAL, NEKO_COMM, ierr)
       call MPI_Unpack(buffer, nbytes, pack_index, &
            amr, 1, MPI_LOGICAL, NEKO_COMM, ierr)
       call MPI_Bcast(params%p, 1, MPI_NEKO_PARAMS, 0, NEKO_COMM, ierr)
    end if

    if (amr) then
       ! initialise mesh manager and create mesh
       C%amr = amr_t(mesh_file)
       C%msh = C%amr%get_mesh()
    else
       msh_file = file_t(trim(mesh_file))
       call msh_file%read(C%msh)
    end if

    C%params = params%p

    !
    ! Setup user defined functions
    !
    call C%usr%init()

    ! pre-refine the mesh before starting fluid
    ! this should be done for fresh run only (no restart)
    if (C%amr%ifamr().and.C%params%restart_file=='') then
       prerefine : block
         integer :: il
         ! give possibility to fill all the refinement levels
         ! mark refinment by negative tstep
         do il = 1, C%params%amrlmax
            call C%amr%msh_refine(C%msh, C%params, C%usr, 0.0_rp, -il)
         end do
       end block prerefine
    end if

    !
    ! Load Balancing
    !
    if (pe_size .gt. 1 .and. C%params%loadb) then
       call neko_log%section('Load Balancing')
       call parmetis_partmeshkway(C%msh, parts)
       call redist_mesh(C%msh, parts)
       call neko_log%end_section()       
    end if

    call C%usr%user_mesh_setup(C%msh)
    
    !
    ! Setup fluid scheme
    !
    call fluid_scheme_factory(C%fluid, trim(fluid_scheme))
    call C%fluid%init(C%msh, lx, C%params)

!!$    testing_refine : block
!!$      ! for now just testing
!!$      integer(i4), allocatable, dimension(:) :: ref_mark
!!$      ! initialise restructure data
!!$      call amr_rcn_init(C%msh, C%fluid%Xh)
!!$      write(*,*) 'TEST int init', pe_rank
!!$      ! perform refinement
!!$      allocate(ref_mark(C%msh%nelv))
!!$      ! mark all the elements to be refined
!!$      ref_mark(:) = 1 ! THIS VALUE SHOULD BE TAKEN FORM p4est_wrap.h
!!$      ! refine filed
!!$      call amr_refine(C%msh, C%fluid, C%params, ref_mark)
!!$      deallocate(ref_mark)
!!$    end block testing_refine

!!$    testing_gs : block
!!$      use gs_ops
!!$      use gather_scatter
!!$      integer :: il, jl, kl, ll, ierr, iunit
!!$      character(len=2) :: pe_sid
!!$      call MPI_Barrier(NEKO_COMM, ierr)
!!$      ! testing gs
!!$      C%fluid%u%x(:,:,:,:) = 1.0
!!$      call gs_op(C%fluid%gs_Xh, C%fluid%u, GS_OP_ADD)
!!$      C%fluid%v%x(:,:,:,:) = C%fluid%dm_Xh%dof(:,:,:,:)
!!$      call gs_op(C%fluid%gs_Xh, C%fluid%v, GS_OP_ADD)
!!$      C%fluid%v%x(:,:,:,:) = C%fluid%v%x(:,:,:,:)/C%fluid%u%x(:,:,:,:)
!!$      write(pe_sid, '(i2.2)') pe_rank
!!$      open(newunit=iunit,file='testgs'//pe_sid//'.txt')
!!$      write(iunit,*) 'TEST',pe_rank,C%msh%nelv
!!$      write(iunit,*) '  '
!!$      write(iunit,*) 'DOF test'
!!$      write(iunit,*) '  '
!!$      do il = 1, C%msh%nelv
!!$         do ll = 1, C%fluid%Xh%lz
!!$            do kl = 1, C%fluid%Xh%ly
!!$               write(iunit,'(3i4,4x,8f8.3)') il, ll, kl, &
!!$                    & (C%fluid%v%x(jl,kl,ll,il)-C%fluid%dm_Xh%dof(jl,kl,ll,il),jl=1,C%fluid%Xh%lx)
!!$            end do
!!$            write(iunit,*) '  '
!!$         end do
!!$         write(iunit,*) '======================================='
!!$      end do
!!$      write(iunit,*) '  '
!!$      write(iunit,*) 'COORDINATES test'
!!$      write(iunit,*) '  '
!!$      do il = 1, C%msh%nelv
!!$         do ll = 1, C%fluid%Xh%lz
!!$            do kl = 1, C%fluid%Xh%ly
!!$               write(iunit,'(a,3i4,2x,8f9.5)') 'x',il, ll, kl, &
!!$                    & (C%fluid%dm_Xh%x(jl,kl,ll,il),jl=1,C%fluid%Xh%lx)
!!$               write(iunit,'(a,3i4,2x,8f9.5)') 'y',il, ll, kl, &
!!$                    & (C%fluid%dm_Xh%y(jl,kl,ll,il),jl=1,C%fluid%Xh%lx)
!!$               write(iunit,'(a,3i4,2x,8f9.5)') 'z',il, ll, kl, &
!!$                    & (C%fluid%dm_Xh%z(jl,kl,ll,il),jl=1,C%fluid%Xh%lx)
!!$               write(iunit,*) '  '
!!$            end do
!!$            write(iunit,*) '-----------------------------------------'
!!$         end do
!!$         write(iunit,*) '======================================='
!!$      end do
!!$      close(iunit)
!!$      call MPI_Barrier(NEKO_COMM, ierr)
!!$      call neko_error('This is not error.')
!!$    end block testing_gs

    !
    ! Setup scalar scheme
    !
    ! @todo no scalar factroy for now, probably not needed
    if (scalar) then
       allocate(C%scalar)
       call C%scalar%init(C%msh, C%fluid%c_Xh, C%fluid%gs_Xh, C%params)
       call C%fluid%chkp%add_scalar(C%scalar%s)
    end if
    !
    ! Setup user defined conditions    
    !
    if (trim(C%params%fluid_inflow) .eq. 'user') then
       call C%fluid%set_usr_inflow(C%usr%fluid_user_if)
    end if
    
    !
    ! Setup source term
    ! 
    if (trim(source_term) .eq. 'user') then
       call C%fluid%set_source(trim(source_term), usr_f=C%usr%fluid_user_f)
    else if (trim(source_term) .eq. 'user_vector') then
       call C%fluid%set_source(trim(source_term), &
            usr_f_vec=C%usr%fluid_user_f_vector)
    else
       call C%fluid%set_source(trim(source_term))
    end if

    ! Setup source term for the scalar
    ! @todo should be expanded for user sources etc. Now copies the fluid one
    if (scalar) then
       if (trim(scalar_source_term) .eq. 'user') then
          call C%scalar%set_source(trim(scalar_source_term), &
               usr_f=C%usr%scalar_user_f)
       else if (trim(scalar_source_term) .eq. 'user_vector') then
          call C%scalar%set_source(trim(scalar_source_term), &
               usr_f_vec=C%usr%scalar_user_f_vector)
       else
          call C%scalar%set_source(trim(scalar_source_term))
       end if
    end if

    !
    ! Setup initial conditions
    ! 
    if (len_trim(initial_condition) .gt. 0) then
       if (trim(initial_condition) .ne. 'user') then
          call set_flow_ic(C%fluid%u, C%fluid%v, C%fluid%w, C%fluid%p, &
               C%fluid%c_Xh, C%fluid%gs_Xh, initial_condition, C%params)
       else
          call set_flow_ic(C%fluid%u, C%fluid%v, C%fluid%w, C%fluid%p, &
               C%fluid%c_Xh, C%fluid%gs_Xh, C%usr%fluid_user_ic, C%params)
       end if
    end if

    ! Add initial conditions to BDF scheme (if present)
    select type(f => C%fluid)
    type is(fluid_pnpn_t)
       call f%ulag%set(f%u)
       call f%vlag%set(f%v)
       call f%wlag%set(f%w)
    type is(fluid_plan4_t)
       call f%ulag%set(f%u)
       call f%vlag%set(f%v)
       call f%wlag%set(f%w)
    type is(device_fluid_plan4_t)
       call f%ulag%set(f%u)
       call f%vlag%set(f%v)
       call f%wlag%set(f%w)
    end select


    !
    ! Validate that the case is properly setup for time-stepping
    !
    call C%fluid%validate

    if (scalar) then
       call C%scalar%slag%set(C%scalar%s)
       call C%scalar%validate
    end if

    !
    ! Set order of timestepper
    !
    call C%ext_bdf%set_time_order(C%params%time_order)

    !
    ! Save boundary markings for fluid (if requested)
    ! 
    if (C%params%output_bdry) then
       bdry_file = file_t(trim(C%params%output_dir)//'bdry.fld')
       call bdry_file%write(C%fluid%bdry)
    end if

    !
    ! Save mesh partitions (if requested)
    !
    if (C%params%output_part) then
       call mesh_field_init(msh_part, C%msh, 'MPI_Rank')
       msh_part%data = pe_rank
       part_file = file_t(trim(C%params%output_dir)//'partitions.vtk')
       call part_file%write(msh_part)
       call mesh_field_free(msh_part)
    end if

    !
    ! Setup sampler
    !
    call C%s%init(C%params%nsamples, C%params%T_end)
    C%f_out = fluid_output_t(C%fluid, path=C%params%output_dir)
    call C%s%add(C%f_out)

    if (scalar) then
       C%s_out = scalar_output_t(C%scalar, path=C%params%output_dir)
       call C%s%add(C%s_out)
    end if

    !
    ! Save checkpoints (if requested)
    !
    if (C%params%output_chkp) then
       C%f_chkp = chkp_output_t(C%fluid%chkp, path=C%params%output_dir)
       call C%s%add(C%f_chkp)
    end if

    !
    ! Setup statistics
    !
    call C%q%init(C%params%stats_begin)
    if (C%params%stats_mean_flow) then
       call C%q%add(C%fluid%mean%u)
       call C%q%add(C%fluid%mean%v)
       call C%q%add(C%fluid%mean%w)
       call C%q%add(C%fluid%mean%p)

       if (C%params%output_mean_flow) then
          C%f_mf = mean_flow_output_t(C%fluid%mean, C%params%stats_begin, &
                                      path=C%params%output_dir)
          call C%s%add(C%f_mf)
       end if
    end if

    if (C%params%stats_mean_sqr_flow) then
       call C%q%add(C%fluid%mean_sqr%uu)
       call C%q%add(C%fluid%mean_sqr%vv)
       call C%q%add(C%fluid%mean_sqr%ww)
       call C%q%add(C%fluid%mean_sqr%pp)

       if (C%params%output_mean_sqr_flow) then
          C%f_msqrf = mean_sqr_flow_output_t(C%fluid%mean_sqr, &
                                             C%params%stats_begin, &
                                             path=C%params%output_dir)
          call C%s%add(C%f_msqrf)
       end if
    end if

    !
    ! Setup joblimit
    !
    call jobctrl_set_time_limit(C%params%jlimit)

    call neko_log%end_section()
    
  end subroutine case_init
  
  !> Deallocate a case 
  subroutine case_free(C)
    type(case_t), intent(inout) :: C

    if (C%amr%ifamr()) call C%amr%free()

    if (allocated(C%fluid)) then
       call C%fluid%free()
       deallocate(C%fluid)
    end if

    if (allocated(C%scalar)) then
       call C%scalar%free()
       deallocate(C%scalar)
    end if

    call mesh_free(C%msh)

    call C%s%free()

    call C%q%free()
    
  end subroutine case_free
  
end module case
