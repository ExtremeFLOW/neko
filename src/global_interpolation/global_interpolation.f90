! Copyright (c) 2020-2025, The Neko Authors
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
!> Implements global_interpolation given a dofmap.
!!
module global_interpolation
  use num_types, only: rp, dp, xp
  use neko_config, only : NEKO_BCKND_DEVICE
  use space, only: space_t
  use stack, only: stack_i4_t
  use dofmap, only: dofmap_t
  use logger, only: neko_log, LOG_SIZE
  use utils, only: neko_error
  use local_interpolation, only : local_interpolator_t
  use device, only: device_free, device_map, device_memcpy, &
       device_deassociate, HOST_TO_DEVICE, DEVICE_TO_HOST, &
       device_get_ptr
  use aabb_pe_finder, only: aabb_pe_finder_t
  use aabb_el_finder, only: aabb_el_finder_t
  use cartesian_el_finder, only: cartesian_el_finder_t
  use cartesian_pe_finder, only: cartesian_pe_finder_t
  use legendre_rst_finder, only: legendre_rst_finder_t
  use el_finder, only: el_finder_t
  use pe_finder, only: pe_finder_t
  use comm, only: NEKO_COMM
  use mpi_f08, only: MPI_SUM, MPI_COMM, MPI_Comm_rank, &
       MPI_Comm_size, MPI_Wtime, MPI_Allreduce, MPI_IN_PLACE, MPI_INTEGER, &
       MPI_MIN, MPI_Barrier, MPI_Reduce_Scatter_block, MPI_alltoall, &
       MPI_ISend, MPI_IRecv
  use glb_intrp_comm, only : glb_intrp_comm_t
  use vector, only: vector_t
  use matrix, only: matrix_t
  use math, only: copy, NEKO_EPS
  use structs, only : array_ptr_t
  use, intrinsic :: iso_c_binding, only: c_ptr, C_NULL_PTR, c_associated
  implicit none
  private

  integer, public, parameter :: GLOB_MAP_SIZE = 4096

  !> Implements global interpolation for arbitrary points in the domain.
  type, public :: global_interpolation_t
     !> X coordinates from which to interpolate.
     type(vector_t) :: x
     !> Y coordinates from which to interpolate.
     type(vector_t) :: y
     !> Z coordinates from which to interpolate.
     type(vector_t) :: z
     !> Geometric dimension of the simulation.
     integer :: gdim
     !> Number of elements.
     integer :: nelv
     !> Global number of elements.
     integer :: glb_nelv
     !> Which communicator to do interpolation on
     type(MPI_COMM) :: comm
     !> pe_rank in comm
     integer :: pe_rank
     !> pe_size of comm
     integer :: pe_size
     !> Space.
     type(space_t) :: Xh
     !> Components related to the points we want to evalute
     !> Number of points we want to evaluate
     integer :: n_points
     !> x,y,z coordinates
     real(kind=rp), allocatable :: xyz(:,:)
     !> r,s,t coordinates
     real(kind=rp), allocatable :: rst(:,:)
     !> List of owning processes.
     integer, allocatable :: pe_owner(:)
     !> array of stacks (to avoid expensive reinitialization)
     type(stack_i4_t), allocatable :: points_at_pe(:)
     !> List of owning elements.
     !! Note this is 0 indexed
     integer, allocatable :: el_owner0(:)
     type(c_ptr) :: el_owner0_d = C_NULL_PTR

     !> Local points (points in this ranks domain)
     integer :: n_points_local
     integer, allocatable :: el_owner0_local(:)
     type(c_ptr) :: el_owner0_local_d = C_NULL_PTR
     real(kind=rp), allocatable :: rst_local(:,:)
     real(kind=rp), allocatable :: xyz_local(:,:)
     !> Interpolator for local points.
     type(local_interpolator_t) :: local_interp
     !> If all points are local on this PE.
     !> Turns true if points are redistributed to their respective owners
     logical :: all_points_local = .false.
     !> Tolerance for Newton solve to find the correct rst coordinates.
     real(kind=rp) :: tol = NEKO_EPS*1e3_rp
     !> Padding
     real(kind=rp) :: padding = 1e-2_rp

     !> Mapping of points to ranks.
     !> n_points_pe(pe_rank) = n_points I have at this rank
     !> Sum(n_points_pe) = n_points
     integer, allocatable :: n_points_pe(:)
     integer, allocatable :: n_points_offset_pe(:)
     !> n_points_pe_local(pe_rank) = Number local points to send to this rank.
     !> Sum(n_points_pe_local) = n_points_local
     integer, allocatable :: n_points_pe_local(:)
     integer, allocatable :: n_points_offset_pe_local(:)
     !> Finding elements
     !> Structure to find rank candidates
     class(pe_finder_t), allocatable :: pe_finder
     !> Structure to find element candidates
     class(el_finder_t), allocatable :: el_finder
     !> Object to find rst coordinates
     type(legendre_rst_finder_t) :: rst_finder
     !> Things for communication operation (sending interpolated values back and forth)
     type(glb_intrp_comm_t) :: glb_intrp_comm
     !> Working vectors for global interpolation
     type(vector_t) :: temp_local, temp

   contains
     !> Initialize the global interpolation object based on a set of spectral elements.
     procedure, pass(this) :: init_xyz => global_interpolation_init_xyz
     !> Initialize the global interpolation object based on a dofmap.
     procedure, pass(this) :: init_dof => global_interpolation_init_dof
     !> Destructor
     procedure, pass(this) :: free => global_interpolation_free
     !> Destructor for arrays related to evaluation points
     procedure, pass(this) :: free_points => global_interpolation_free_points
     procedure, pass(this) :: free_points_local => global_interpolation_free_points_local
     procedure, pass(this) :: find_points_and_redist => &
          global_interpolation_find_and_redist
     !> Finds the process owner, global element number,
     !! and local rst coordinates for each point.
     !! Sets up correct values to be able to evalute the points
     procedure, pass(this) :: find_points_coords => &
          global_interpolation_find_coords
     procedure, pass(this) :: find_points_coords1d => &
          global_interpolation_find_coords1d
     !> Subroutine to check if the points are within the tolerance
     procedure, pass(this) :: check_points => &
          global_interpolation_check_points
     procedure, pass(this) :: find_points_xyz => global_interpolation_find_xyz
     generic :: find_points => find_points_xyz, find_points_coords, find_points_coords1d
     !> Evaluate the value of the field in each point.
     procedure, pass(this) :: evaluate => global_interpolation_evaluate

     !> Generic constructor
     generic :: init => init_dof, init_xyz

  end type global_interpolation_t

contains

  !> Initialize the global interpolation object on a dofmap.
  !! @param dof Dofmap on which the interpolation is to be carried out.
  !! @param tol Tolerance for Newton iterations.
  subroutine global_interpolation_init_dof(this, dof, comm, tol, pad)
    class(global_interpolation_t), target, intent(inout) :: this
    type(dofmap_t) :: dof
    type(MPI_COMM), optional, intent(in) :: comm
    real(kind=rp), optional :: tol
    real(kind=rp), optional :: pad
    real(kind=rp) :: padding, tolerance

    if (present(pad)) then
       padding = pad
    else
       padding = 1e-2 ! 1% padding of the bounding boxes
    end if

    if (present(tol)) then
       tolerance = tol
    else
       tolerance = NEKO_EPS*1e3_xp ! 1% padding of the bounding boxes
    end if
    ! NOTE: Passing dof%x(:,1,1,1), etc in init_xyz passes down the entire
    ! dof%x array and not a slice. It is done this way for
    ! to get the right dimension (see global_interpolation_init_xyz).
    call this%init_xyz(dof%x(:,1,1,1), dof%y(:,1,1,1), dof%z(:,1,1,1), &
         dof%msh%gdim, dof%msh%nelv, dof%Xh, comm,tol = tolerance, pad=padding)

  end subroutine global_interpolation_init_dof

  !> Initialize the global interpolation object on a set of coordinates.
  !! @param x x-coordinates.
  !! @param y y-coordinates.
  !! @param z z-coordinates.
  !! @param gdim Geometric dimension.
  !! @param nelv Number of elements of the mesh in which to search for the
  !! points.
  !! @param Xh Space on which to interpolate.
  !! @param tol Tolerance for Newton iterations.
  !! @param pad Padding of the bounding boxes.
  subroutine global_interpolation_init_xyz(this, x, y, z, gdim, nelv, Xh, &
       comm, tol, pad)
    class(global_interpolation_t), target, intent(inout) :: this
    real(kind=rp), intent(in) :: x(:)
    real(kind=rp), intent(in) :: y(:)
    real(kind=rp), intent(in) :: z(:)
    integer, intent(in) :: gdim
    integer, intent(in) :: nelv
    type(MPI_COMM), intent(in), optional :: comm
    type(space_t), intent(in) :: Xh
    real(kind=rp), intent(in), optional :: tol
    real(kind=rp), intent(in), optional :: pad
    integer :: lx, ly, lz, ierr, i, n
    character(len=8000) :: log_buf
    real(kind=dp) :: padding
    real(kind=rp) :: time1, time_start
    character(len=255) :: mode_str
    integer :: boxdim, envvar_len

    call this%free()

    if (present(comm)) then
       this%comm = comm
    else
       this%comm = NEKO_COMM
    end if

    if (present(pad)) then
       padding = pad
       this%padding = pad
    else
       padding = 1e-2 ! 1% padding of the bounding boxes
       this%padding = 1e-2
    end if

    time_start = MPI_Wtime()
    call MPI_Barrier(this%comm)

    call MPI_Comm_rank(this%comm, this%pe_rank, ierr)
    call MPI_Comm_size(this%comm, this%pe_size, ierr)

    this%gdim = gdim
    this%nelv = nelv
    if (present(tol)) then
       this%tol = tol
    else
       this%tol = NEKO_EPS*1e3_xp
    end if

    call MPI_Allreduce(nelv, this%glb_nelv, 1, MPI_INTEGER, &
         MPI_SUM, this%comm, ierr)
    lx = Xh%lx
    ly = Xh%ly
    lz = Xh%lz
    n = nelv * lx*ly*lz
    call this%x%init(n)
    call this%y%init(n)
    call this%z%init(n)
    call copy(this%x%x, x, n)
    call this%x%copy_from(HOST_TO_DEVICE,.false.)
    call copy(this%y%x, y, n)
    call this%y%copy_from(HOST_TO_DEVICE,.false.)
    call copy(this%z%x, z, n)
    call this%z%copy_from(HOST_TO_DEVICE,.false.)
    call this%Xh%init(Xh%t, lx, ly, lz)


    call neko_log%section('Global Interpolation')

    call neko_log%message('Initializing global interpolation')
    call get_environment_variable("NEKO_GLOBAL_INTERP_EL_FINDER", &
         mode_str, envvar_len)

    if (envvar_len .gt. 0) then
       if (mode_str(1:envvar_len) == 'AABB') then
          allocate(aabb_el_finder_t :: this%el_finder)
       end if
    end if
    call get_environment_variable("NEKO_GLOBAL_INTERP_PE_FINDER", &
         mode_str, envvar_len)

    if (envvar_len .gt. 0) then
       if (mode_str(1:envvar_len) == 'AABB') then
          allocate(aabb_pe_finder_t :: this%pe_finder)
       end if
    end if

    if (.not. allocated(this%el_finder)) then
       allocate(cartesian_el_finder_t :: this%el_finder)
    end if
    if (.not. allocated(this%pe_finder)) then
       allocate(cartesian_pe_finder_t :: this%pe_finder)
    end if
    select type(el_find => this%el_finder)
    type is (aabb_el_finder_t)
       call neko_log%message('Using AABB element finder')
       call el_find%init(x, y, z, nelv, Xh, padding)
    type is (cartesian_el_finder_t)
       call neko_log%message('Using Cartesian element finder')
       boxdim = max(lx*int(real(nelv,xp)**(1.0_xp/3.0_xp)),2)
       boxdim = min(boxdim, 300)
       call el_find%init(x, y, z, nelv, Xh, boxdim, padding)
    class default
       call neko_error('Unknown element finder type')
    end select

    select type(pe_find => this%pe_finder)
    type is (aabb_pe_finder_t)
       call neko_log%message('Using AABB PE finder')
       call pe_find%init(this%x%x, this%y%x, this%z%x, &
            nelv, Xh, this%comm, padding)
    type is (cartesian_pe_finder_t)
       call neko_log%message('Using Cartesian PE finder')
       boxdim = lx*int(real(this%glb_nelv,xp)**(1.0_xp/3.0_xp))
       boxdim = max(boxdim,32)
       boxdim = min(boxdim, &
            int(8.0_xp*(30000.0_xp*this%pe_size)**(1.0_xp/3.0_xp)))
       call pe_find%init(this%x%x, this%y%x, this%z%x, &
            nelv, Xh, this%comm, boxdim, padding)
    class default
       call neko_error('Unknown PE finder type')
    end select

    call this%rst_finder%init(this%x%x, this%y%x, this%z%x, nelv, Xh, this%tol)
    if (allocated(this%n_points_pe)) deallocate(this%n_points_pe)
    if (allocated(this%n_points_pe_local)) deallocate(this%n_points_pe_local)
    if (allocated(this%n_points_offset_pe_local)) &
         deallocate(this%n_points_offset_pe_local)
    if (allocated(this%n_points_offset_pe)) deallocate(this%n_points_offset_pe)
    allocate(this%n_points_pe(0:(this%pe_size-1)))
    allocate(this%n_points_offset_pe(0:(this%pe_size-1)))
    allocate(this%n_points_pe_local(0:(this%pe_size-1)))
    allocate(this%n_points_offset_pe_local(0:(this%pe_size-1)))
    allocate(this%points_at_pe(0:(this%pe_size-1)))
    do i = 0, this%pe_size-1
       call this%points_at_pe(i)%init()
    end do
    call MPI_Barrier(this%comm)
    time1 = MPI_Wtime()
    write(log_buf, '(A,E15.7)') &
         'Global interpolation initialized (s):', time1-time_start
    call neko_log%message(log_buf)
    call neko_log%end_section()
  end subroutine global_interpolation_init_xyz


  !> Destructor
  subroutine global_interpolation_free(this)
    class(global_interpolation_t), target, intent(inout) :: this
    integer :: i

    call this%x%free()
    call this%y%free()
    call this%z%free()
    call this%Xh%free()

    this%nelv = 0
    this%gdim = 0

    call this%free_points()
    call this%free_points_local()
    call this%local_interp%free()
    if (allocated(this%el_finder)) then
       call this%el_finder%free()
       deallocate(this%el_finder)
    end if
    if (allocated(this%pe_finder)) then
       call this%pe_finder%free()
       deallocate(this%pe_finder)
    end if
    call this%rst_finder%free()

    call this%temp_local%free()
    call this%temp%free()
    if (allocated(this%points_at_pe)) then
       do i = 0, this%pe_size-1
          call this%points_at_pe(i)%free()
       end do
       deallocate(this%points_at_pe)
    end if
    if (allocated(this%n_points_pe)) deallocate(this%n_points_pe)
    if (allocated(this%n_points_pe_local)) deallocate(this%n_points_pe_local)
    if (allocated(this%n_points_offset_pe_local)) &
         deallocate(this%n_points_offset_pe_local)
    if (allocated(this%n_points_offset_pe)) deallocate(this%n_points_offset_pe)



  end subroutine global_interpolation_free

  !> Destructor for point arrays
  subroutine global_interpolation_free_points(this)
    class(global_interpolation_t), target, intent(inout) :: this

    this%n_points = 0
    this%all_points_local = .false.

    if (allocated(this%xyz)) deallocate(this%xyz)
    if (allocated(this%rst)) deallocate(this%rst)
    if (allocated(this%pe_owner)) deallocate(this%pe_owner)
    if (allocated(this%el_owner0)) deallocate(this%el_owner0)

    if (c_associated(this%el_owner0_d)) then
       call device_free(this%el_owner0_d)
    end if

    call this%glb_intrp_comm%free()


  end subroutine global_interpolation_free_points


  subroutine global_interpolation_free_points_local(this)
    class(global_interpolation_t), target, intent(inout) :: this

    this%n_points_local = 0
    this%all_points_local = .false.

    if (allocated(this%xyz_local)) deallocate(this%xyz_local)
    if (allocated(this%rst_local)) deallocate(this%rst_local)
    if (allocated(this%el_owner0_local)) deallocate(this%el_owner0_local)

    if (c_associated(this%el_owner0_local_d)) then
       call device_free(this%el_owner0_local_d)
    end if

  end subroutine global_interpolation_free_points_local


  !> Common routine for finding the points.
  subroutine global_interpolation_find_common(this)
    class(global_interpolation_t), target, intent(inout) :: this
    character(len=8000) :: log_buf
    type(vector_t) :: x_t
    type(vector_t) :: y_t
    type(vector_t) :: z_t
    type(matrix_t) :: rst_local_cand
    type(vector_t) :: resx
    type(vector_t) :: resy
    type(vector_t) :: resz
    type(c_ptr) :: el_cands_d
    type(matrix_t) :: res
    integer :: i, j, stupid_intent
    type(stack_i4_t), target :: all_el_candidates
    integer, allocatable :: n_el_cands(:)
    integer, contiguous, pointer :: el_cands(:), point_ids(:), send_recv(:)
    real(kind=rp), allocatable :: res_results(:,:)
    real(kind=rp), allocatable :: rst_results(:,:)
    integer, allocatable :: el_owner_results(:)
    integer :: ierr, ii, n_point_cand, n_glb_point_cand, point_id, rank
    real(kind=rp) :: time1, time2, time_start
    !Temp stuff for glb_intrp_comm
    type(stack_i4_t) :: send_pe, recv_pe
    type(glb_intrp_comm_t) :: glb_intrp_find, glb_intrp_find_back
    type(stack_i4_t) :: send_pe_find, recv_pe_find

    el_cands_d = C_NULL_PTR
    call neko_log%section('Global Interpolation')
    call glb_intrp_find%init_dofs(this%pe_size)
    call send_pe_find%init()
    call recv_pe_find%init()
    call MPI_Barrier(this%comm)
    time_start = MPI_Wtime()
    write(log_buf,'(A)') 'Global interpolation, finding points'
    call neko_log%message(log_buf)
    ! Find pe candidates that the points i want may be at
    ! Add number to n_points_pe_local
    !Working arrays
    this%n_points_pe = 0
    call this%pe_finder%find_batch(this%xyz, this%n_points, &
         this%points_at_pe, this%n_points_pe)
    call MPI_Barrier(this%comm)
    time1 = MPI_Wtime()
    write(log_buf, '(A,E15.7)') &
         'Found PE candidates time since start of findpts (s):', time1-time_start
    call neko_log%message(log_buf)

    !Send number of points I want to candidates
    ! n_points_local -> how many points might be at this rank
    ! n_points_pe_local -> how many points local on this rank that other pes might want
    this%n_points_pe_local = 0
    this%n_points_local = 0
    call MPI_Reduce_scatter_block(this%n_points_pe, this%n_points_local, &
         1, MPI_INTEGER, MPI_SUM, this%comm, ierr)
    call MPI_Alltoall(this%n_points_pe, 1, MPI_INTEGER,&
         this%n_points_pe_local, 1, MPI_INTEGER, this%comm, ierr)

    !Set up offset arrays
    this%n_points_offset_pe_local(0) = 0
    this%n_points_offset_pe(0) = 0
    do i = 1, (this%pe_size - 1)
       this%n_points_offset_pe_local(i) = this%n_points_pe_local(i-1)&
            + this%n_points_offset_pe_local(i-1)
       this%n_points_offset_pe(i) = this%n_points_pe(i-1)&
            + this%n_points_offset_pe(i-1)
    end do
    do i = 0, (this%pe_size-1)
       if (this%n_points_pe(i) .gt. 0) then
          call send_pe_find%push(i)
          point_ids => this%points_at_pe(i)%array()
          do j = 1, this%n_points_pe(i)
             call glb_intrp_find%send_dof(i)%push(3*(point_ids(j)-1)+1)
             call glb_intrp_find%send_dof(i)%push(3*(point_ids(j)-1)+2)
             call glb_intrp_find%send_dof(i)%push(3*(point_ids(j)-1)+3)
          end do
       end if
       if (this%n_points_pe_local(i) .gt. 0) then
          call recv_pe_find%push(i)
          do j = 1, this%n_points_pe_local(i)
             call glb_intrp_find%recv_dof(i)%push(3*(j+this%n_points_offset_pe_local(i)-1)+1)
             call glb_intrp_find%recv_dof(i)%push(3*(j+this%n_points_offset_pe_local(i)-1)+2)
             call glb_intrp_find%recv_dof(i)%push(3*(j+this%n_points_offset_pe_local(i)-1)+3)
          end do
       end if
    end do



    call glb_intrp_find%init(send_pe_find, recv_pe_find, this%comm)

    call glb_intrp_find_back%init_dofs(this%pe_size)
    ii = 0
    do i = 0, (this%pe_size-1)
       send_recv => glb_intrp_find%recv_dof(i)%array()
       do j = 1, glb_intrp_find%recv_dof(i)%size()
          call glb_intrp_find_back%send_dof(i)%push(send_recv(j))
       end do
       send_recv => glb_intrp_find%send_dof(i)%array()
       do j = 1, glb_intrp_find%send_dof(i)%size()
          ii = ii + 1
          call glb_intrp_find_back%recv_dof(i)%push(ii)
       end do
    end do

    call glb_intrp_find_back%init(recv_pe_find, send_pe_find, this%comm)


    if (allocated(this%xyz_local)) then
       deallocate(this%xyz_local)
    end if
    allocate(this%xyz_local(3, this%n_points_local))
    call glb_intrp_find%sendrecv(this%xyz, this%xyz_local, this%n_points*3, &
         this%n_points_local*3)

    call MPI_Barrier(this%comm)
    time1 = MPI_Wtime()
    write(log_buf, '(A,E15.7)') &
         'Sent to points to PE candidates, time since start of find_points (s):', time1-time_start
    call neko_log%message(log_buf)
    !Okay, now we need to find the rst...
    call all_el_candidates%init()

    if (allocated(n_el_cands)) then
       deallocate(n_el_cands)
    end if

    allocate(n_el_cands(this%n_points_local))
    !> Find element candidates at this rank
    call this%el_finder%find_batch(this%xyz_local, this%n_points_local, &
         all_el_candidates, n_el_cands)

    n_point_cand = all_el_candidates%size()
    if (n_point_cand .gt. 1e8) then
       print *,'Warning, many point candidates on rank', this%pe_rank,'cands:', n_point_cand, &
            'Consider increasing number of ranks'
    end if
    call x_t%init(n_point_cand)
    call y_t%init(n_point_cand)
    call z_t%init(n_point_cand)
    ii = 0
    !> Copy xyz coords to each element candidate
    do i = 1 , this%n_points_local
       do j = 1, n_el_cands(i)
          ii = ii + 1
          x_t%x(ii) = this%xyz_local(1,i)
          y_t%x(ii) = this%xyz_local(2,i)
          z_t%x(ii) = this%xyz_local(3,i)
       end do
    end do

    call MPI_Barrier(this%comm)
    time1 = MPI_Wtime()
    write(log_buf, '(A,E15.7)') &
         'Element candidates found, now time for finding rst,time since start of find_points (s):', time1-time_start
    call neko_log%message(log_buf)
    call rst_local_cand%init(3,n_point_cand)
    call resx%init(n_point_cand)
    call resy%init(n_point_cand)
    call resz%init(n_point_cand)

    ! Find rst within all element candidates for target xyz (x_t, y_t, z_t)
    call MPI_Barrier(this%comm)
    time1 = MPI_Wtime()
    el_cands => all_el_candidates%array()
    if (NEKO_BCKND_DEVICE .eq. 1) then
       call x_t%copy_from(HOST_TO_DEVICE,.false.)

       call y_t%copy_from(HOST_TO_DEVICE,.false.)

       call z_t%copy_from(HOST_TO_DEVICE,.false.)
       call device_map(el_cands, el_cands_d,n_point_cand)
       call device_memcpy(el_cands, el_cands_d,n_point_cand, &
            HOST_TO_DEVICE, .true.)
    end if

    call this%rst_finder%find(rst_local_cand, &
         x_t, y_t, z_t, &
         el_cands, n_point_cand, &
         resx, resy, resz)
    if (NEKO_BCKND_DEVICE .eq. 1) then
       call rst_local_cand%copy_from(DEVICE_TO_HOST,.false.)
       call resx%copy_from(DEVICE_TO_HOST,.false.)
       call resy%copy_from(DEVICE_TO_HOST,.false.)
       call resz%copy_from(DEVICE_TO_HOST,.true.)
       call device_deassociate(el_cands)
       call device_free(el_cands_d)
    end if
    call MPI_Barrier(this%comm)

    time2 = MPI_Wtime()
    write(log_buf, '(A,E15.7)') &
         'Found rst with Newton iteration, time (s):', time2-time1
    call neko_log%message(log_buf)

    write(log_buf,'(A,E15.7)') &
         'Tolerance: ', this%tol
    call neko_log%message(log_buf)
    write(log_buf,'(A)') &
         'Checking validity of points and choosing best candidates.'
    call neko_log%message(log_buf)
    call MPI_Barrier(this%comm,ierr)

    if (allocated(this%rst_local)) deallocate(this%rst_local)
    if (allocated(this%el_owner0_local)) deallocate(this%el_owner0_local)
    allocate(this%rst_local(3,this%n_points_local))
    allocate(this%el_owner0_local(this%n_points_local))
    ! Choose the best candidate at this rank
    ii = 0
    do i = 1 , this%n_points_local
       this%xyz_local(1,i) = 10.0
       this%xyz_local(2,i) = 10.0
       this%xyz_local(3,i) = 10.0
       this%rst_local(1,i) = 10.0
       this%rst_local(2,i) = 10.0
       this%rst_local(3,i) = 10.0
       this%el_owner0_local(i) = -1
       do j = 1, n_el_cands(i)
          ii = ii + 1
          if (rst_cmp(this%rst_local(:,i), rst_local_cand%x(:,ii),&
               this%xyz_local(:,i), (/resx%x(ii),resy%x(ii),resz%x(ii)/), this%padding)) then
             this%rst_local(1,i) = rst_local_cand%x(1,ii)
             this%rst_local(2,i) = rst_local_cand%x(2,ii)

             this%rst_local(3,i) = rst_local_cand%x(3,ii)
             this%xyz_local(1,i) = resx%x(ii)
             this%xyz_local(2,i) = resy%x(ii)
             this%xyz_local(3,i) = resz%x(ii)
             this%el_owner0_local(i) = el_cands(ii)
          end if
          ! if (this%pe_rank .eq. 0) print *,i,  this%rst_local(:,i), &
          !  this%xyz_local(:,i), this%el_owner0_local(i)
       end do
    end do
    call res%init(3,this%n_points)
    n_glb_point_cand = sum(this%n_points_pe)
    if (allocated(rst_results)) deallocate(rst_results)
    if (allocated(res_results)) deallocate(res_results)
    if (allocated(el_owner_results)) deallocate(el_owner_results)
    allocate(rst_results(3,n_glb_point_cand))
    allocate(res_results(3,n_glb_point_cand))
    allocate(el_owner_results(n_glb_point_cand))
    res = 1e2_rp
    this%rst = 1e2
    this%pe_owner = -1
    this%el_owner0 = -1
    call glb_intrp_find_back%sendrecv(this%xyz_local, res_results, &
         this%n_points_local*3, n_glb_point_cand*3)
    call glb_intrp_find_back%sendrecv(this%rst_local, rst_results, &
         this%n_points_local*3, n_glb_point_cand*3)
    do i = 1, size(glb_intrp_find_back%send_pe)
       rank = glb_intrp_find_back%send_pe(i)
       call MPI_Isend(this%el_owner0_local(this%n_points_offset_pe_local(rank)+1),&
            this%n_points_pe_local(rank), &
            MPI_INTEGER, rank, 0, &
            this%comm, glb_intrp_find_back%send_buf(i)%request, ierr)
       glb_intrp_find_back%send_buf(i)%flag = .false.
    end do
    do i = 1, size(glb_intrp_find_back%recv_pe)
       rank = glb_intrp_find_back%recv_pe(i)
       call MPI_IRecv(el_owner_results(this%n_points_offset_pe(rank)+1),&
            this%n_points_pe(rank), &
            MPI_INTEGER, rank, 0, &
            this%comm, glb_intrp_find_back%recv_buf(i)%request, ierr)
       glb_intrp_find_back%recv_buf(i)%flag = .false.
    end do
    call glb_intrp_find_back%nbwait_no_op()
    ii = 0
    do i = 1, size(glb_intrp_find_back%recv_pe)
       point_ids => this%points_at_pe(glb_intrp_find_back%recv_pe(i))%array()
       do j = 1, this%n_points_pe(glb_intrp_find_back%recv_pe(i))
          point_id = point_ids(j)
          ii = ii + 1
          if (rst_cmp(this%rst(:,point_id), rst_results(:,ii), &
               res%x(:,point_id), res_results(:,ii), this%padding) .or. &
               this%pe_owner(point_ids(j)) .eq. -1 ) then
             this%rst(:,point_ids(j)) = rst_results(:,ii)
             res%x(:,point_ids(j)) = res_results(:,ii)
             this%pe_owner(point_ids(j)) = glb_intrp_find_back%recv_pe(i)
             this%el_owner0(point_ids(j)) = el_owner_results(ii)
          end if
          !  if (this%pe_rank .eq. 0) print *,point_id,  &
          !this%rst(:,point_ids(j)),res%x(:,point_ids(j)), this%el_owner0(point_ids(j))
       end do
    end do

    !OK, now I know the correct rst values of the points I want We now
    !send the correct rsts to the correct rank (so a point only
    !belongs to one rank)
    do i = 0, this%pe_size-1
       call this%points_at_pe(i)%clear()
       this%n_points_pe(i) = 0
    end do

    do i = 1, this%n_points
       stupid_intent = i
       if (this%pe_owner(i) .eq. -1 .or. this%el_owner0(i) .eq. -1) then
          print *, 'No owning rank found for',&
               ' point ', stupid_intent, ' with coords', this%xyz(:,i), &
               ' Interpolation will always yield 0.0. Try increase padding.'
       else
          call this%points_at_pe(this%pe_owner(i))%push(stupid_intent)
          this%n_points_pe(this%pe_owner(i)) = this%n_points_pe(this%pe_owner(i)) + 1
       end if
    end do
    call MPI_Reduce_scatter_block(this%n_points_pe, this%n_points_local, 1, MPI_INTEGER, &
         MPI_SUM, this%comm, ierr)
    call MPI_Alltoall(this%n_points_pe, 1, MPI_INTEGER,&
         this%n_points_pe_local, 1, MPI_INTEGER, this%comm, ierr)
    this%n_points_offset_pe_local(0) = 0
    this%n_points_offset_pe(0) = 0
    do i = 1, (this%pe_size - 1)
       this%n_points_offset_pe_local(i) = this%n_points_pe_local(i-1)&
            + this%n_points_offset_pe_local(i-1)
       this%n_points_offset_pe(i) = this%n_points_pe(i-1)&
            + this%n_points_offset_pe(i-1)
    end do
    call send_pe_find%free()
    call recv_pe_find%free()
    call glb_intrp_find%free()
    call send_pe_find%init()
    call recv_pe_find%init()
    call glb_intrp_find%init_dofs(this%pe_size)
    !setup comm to send xyz and rst to chosen ranks
    do i = 0, (this%pe_size-1)
       if (this%n_points_pe(i) .gt. 0) then
          call send_pe_find%push(i)
          point_ids => this%points_at_pe(i)%array()
          do j = 1, this%n_points_pe(i)
             call glb_intrp_find%send_dof(i)%push(3*(point_ids(j)-1)+1)
             call glb_intrp_find%send_dof(i)%push(3*(point_ids(j)-1)+2)
             call glb_intrp_find%send_dof(i)%push(3*(point_ids(j)-1)+3)
          end do
       end if
       if (this%n_points_pe_local(i) .gt. 0) then
          call recv_pe_find%push(i)
          do j = 1, this%n_points_pe_local(i)
             call glb_intrp_find%recv_dof(i)%push(3*(j+this%n_points_offset_pe_local(i)-1)+1)
             call glb_intrp_find%recv_dof(i)%push(3*(j+this%n_points_offset_pe_local(i)-1)+2)
             call glb_intrp_find%recv_dof(i)%push(3*(j+this%n_points_offset_pe_local(i)-1)+3)
          end do
       end if
    end do


    call glb_intrp_find%init(send_pe_find, recv_pe_find, this%comm)
    call glb_intrp_find%sendrecv(this%xyz, this%xyz_local, this%n_points*3, &
         this%n_points_local*3)
    call glb_intrp_find%sendrecv(this%rst, this%rst_local, this%n_points*3, &
         this%n_points_local*3)
    ii = 0
    do i = 1, size(glb_intrp_find%send_pe)
       rank = glb_intrp_find%send_pe(i)
       point_ids => this%points_at_pe(rank)%array()
       do j = 1, this%n_points_pe(rank)
          ii = ii + 1
          el_owner_results(ii) = this%el_owner0(point_ids(j))
       end do
       call MPI_Isend(el_owner_results(this%n_points_offset_pe(rank)+1),&
            this%n_points_pe(rank), &
            MPI_INTEGER, rank, 0, &
            this%comm, glb_intrp_find%send_buf(i)%request, ierr)
       glb_intrp_find%send_buf(i)%flag = .false.
    end do
    do i = 1, size(glb_intrp_find%recv_pe)
       rank = glb_intrp_find%recv_pe(i)
       call MPI_IRecv(this%el_owner0_local(this%n_points_offset_pe_local(rank)+1), &
            this%n_points_pe_local(rank), &
            MPI_INTEGER, rank, 0, &
            this%comm, glb_intrp_find%recv_buf(i)%request, ierr)
       glb_intrp_find%recv_buf(i)%flag = .false.
    end do
    call glb_intrp_find%nbwait_no_op()

    call glb_intrp_find%free()

    !Set up final way of doing communication
    call send_pe%init()
    call recv_pe%init()
    call this%glb_intrp_comm%init_dofs(this%pe_size)
    do i = 0, (this%pe_size-1)
       if (this%n_points_pe(i) .gt. 0) then
          call recv_pe%push(i)
          point_ids => this%points_at_pe(i)%array()
          do j = 1, this%n_points_pe(i)
             call this%glb_intrp_comm%recv_dof(i)%push(point_ids(j))
          end do
       end if
       if (this%n_points_pe_local(i) .gt. 0) then
          call send_pe%push(i)
          do j = 1, this%n_points_pe_local(i)
             call this%glb_intrp_comm%send_dof(i)%push(j+this%n_points_offset_pe_local(i))
          end do
       end if
    end do
    call this%glb_intrp_comm%init(send_pe, recv_pe,this%comm)

    !Initialize working arrays for evaluation
    call this%temp_local%init(this%n_points_local)
    call this%temp%init(this%n_points)

    !Initialize interpolator for local interpolation
    call this%local_interp%init(this%Xh, this%rst_local, &
         this%n_points_local)


    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_memcpy(this%el_owner0, this%el_owner0_d, &
            this%n_points, HOST_TO_DEVICE, sync = .true.)
       call device_map(this%el_owner0_local, this%el_owner0_local_d, this%n_points_local)
       call device_memcpy(this%el_owner0_local, this%el_owner0_local_d, &
            this%n_points_local, HOST_TO_DEVICE, sync = .true.)
    end if

    call this%check_points(this%x%x,this%y%x, this%z%x)

    !Free stuff
    call send_pe%free()
    call recv_pe%free()
    call glb_intrp_find_back%free()
    call send_pe_find%free()
    call recv_pe_find%free()
    call x_t%free()
    call y_t%free()
    call z_t%free()
    call rst_local_cand%free()
    call resx%free()
    call resy%free()
    call resz%free()
    call res%free()
    call all_el_candidates%free()

    if (allocated(n_el_cands)) deallocate(n_el_cands)
    if (allocated(rst_results)) deallocate(rst_results)
    if (allocated(res_results)) deallocate(res_results)
    if (allocated(el_owner_results)) deallocate(el_owner_results)
    call MPI_Barrier(this%comm, ierr)
    time2 = MPI_Wtime()
    write(log_buf, '(A,E15.7)') 'Global interpolation find points done, time (s):', &
         time2-time_start
    call neko_log%message(log_buf)
    call neko_log%end_section()
    call neko_log%newline()
  end subroutine global_interpolation_find_common

  !> Check the points for validity
  !! This is used to check that the points are valid and that the interpolation
  !! is correct. It checks that the points are within the tolerance.
  subroutine global_interpolation_check_points(this, x, y, z)
    class(global_interpolation_t), target, intent(inout) :: this
    real(kind=rp), intent(inout) :: x(:)
    real(kind=rp), intent(inout) :: y(:)
    real(kind=rp), intent(inout) :: z(:)
    integer :: i, j
    character(len=8000) :: log_buf
    real(kind=rp) :: xdiff, ydiff, zdiff
    logical :: isdiff
    type(vector_t) :: x_check, y_check, z_check

    call x_check%init(this%n_points)
    call y_check%init(this%n_points)
    call z_check%init(this%n_points)
    call this%evaluate(x_check%x, x, on_host=.true.)
    call this%evaluate(y_check%x, y, on_host=.true.)
    call this%evaluate(z_check%x, z, on_host=.true.)
    write(log_buf,'(A)') 'Checking validity of points.'
    call neko_log%message(log_buf)
    j = 0
    do i = 1 , this%n_points
       ! Check validity of points
       isdiff = .false.
       xdiff = x_check%x(i)-this%xyz(1,i)
       ydiff = y_check%x(i)-this%xyz(2,i)
       zdiff = z_check%x(i)-this%xyz(3,i)
       isdiff = norm2(real((/xdiff,ydiff,zdiff/),xp)) > this%tol
       if (isdiff) then
          write(*,*) 'Point ', i,'at rank ', this%pe_rank, 'with coordinates: ', &
               this%xyz(1, i), this%xyz(2, i), this%xyz(3, i), &
               'Differ from interpolated coords: ', &
               x_check%x(i), y_check%x(i), z_check%x(i), &
               'Actual difference: ', &
               xdiff, ydiff, zdiff, norm2(real((/xdiff,ydiff,zdiff/),xp)),&
               'Process, element: ', &
               this%pe_owner(i), this%el_owner0(i)+1, &
               'Calculated rst: ', &
               this%rst(1,i), this%rst(2,i), this%rst(3,i)
          j = j + 1
       end if
    end do
    call x_check%free()
    call y_check%free()
    call z_check%free()

  end subroutine global_interpolation_check_points

  !> Finds the corresponding r,s,t coordinates
  !! in the correct global element as well as which process that owns the point.
  !! After this the values at these points can be evaluated.
  !! @param x The x-coordinates of the points.
  !! @param y The y-coordinates of the points.
  !! @param z The z-coordinates of the points.
  !! @param n_points The number of points.
  subroutine global_interpolation_find_coords(this, x, y, z, n_points)
    class(global_interpolation_t), intent(inout) :: this
    integer :: n_points
    !!Perhaps this should be kind dp
    !!this is to get around that x,y,z often is 4 dimensional...
    !!Should maybe add interface for 1d aswell
    real(kind=rp) :: x(n_points,1,1,1)
    real(kind=rp) :: y(n_points,1,1,1)
    real(kind=rp) :: z(n_points,1,1,1)
    integer :: i

    call this%free_points()
    call this%free_points_local()

    this%n_points = n_points

    call global_interpolation_init_point_arrays(this)

    !Deepcopy of coordinates
    do i = 1, n_points
       this%xyz(1, i) = x(i,1,1,1)
       this%xyz(2, i) = y(i,1,1,1)
       this%xyz(3, i) = z(i,1,1,1)
    end do

    call global_interpolation_find_common(this)

  end subroutine global_interpolation_find_coords
  !> Finds the corresponding r,s,t coordinates
  !! in the correct global element as well as which process that owns the point.
  !! After this the values at these points can be evaluated.
  !! @param x The x-coordinates of the points.
  !! @param y The y-coordinates of the points.
  !! @param z The z-coordinates of the points.
  !! @param n_points The number of points.
  subroutine global_interpolation_find_coords1d(this, x, y, z, n_points)
    class(global_interpolation_t), intent(inout) :: this
    integer :: n_points
    real(kind=rp) :: x(n_points)
    real(kind=rp) :: y(n_points)
    real(kind=rp) :: z(n_points)
    integer :: i

    call this%free_points()
    call this%free_points_local()

    this%n_points = n_points

    call global_interpolation_init_point_arrays(this)

    !Deepcopy of coordinates
    do i = 1, n_points
       this%xyz(1, i) = x(i)
       this%xyz(2, i) = y(i)
       this%xyz(3, i) = z(i)
    end do

    call global_interpolation_find_common(this)

  end subroutine global_interpolation_find_coords1d


  subroutine global_interpolation_init_point_arrays(this)
    class(global_interpolation_t) :: this

    if (allocated(this%xyz)) deallocate(this%xyz)
    if (allocated(this%rst)) deallocate(this%rst)
    if (allocated(this%pe_owner)) deallocate(this%pe_owner)
    if (allocated(this%el_owner0)) deallocate(this%el_owner0)

    allocate(this%pe_owner(this%n_points))
    allocate(this%el_owner0(this%n_points))
    allocate(this%xyz(3, this%n_points))
    allocate(this%rst(3, this%n_points))
    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_map(this%el_owner0, this%el_owner0_d, this%n_points)
    end if

  end subroutine global_interpolation_init_point_arrays

  !> Finds the corresponding r,s,t coordinates
  !! in the correct global element as well as which process that owns the point.
  !! After this the values at these points can be evaluated.
  !! If the locations of the points change this must be called again.
  !! @param xyz The coordinates of the points.
  !! @param n_points The number of points.
  subroutine global_interpolation_find_xyz(this, xyz, n_points)
    class(global_interpolation_t), intent(inout) :: this
    integer, intent(in) :: n_points
    !!Perhaps this should be kind dp
    real(kind=rp), intent(inout) :: xyz(3, n_points)


    call this%free_points()
    call this%free_points_local()

    this%n_points = n_points

    call global_interpolation_init_point_arrays(this)

    !> make deep copy incase xyz goes out of scope or deallocated
    call copy(this%xyz, xyz, 3 * n_points)

    call global_interpolation_find_common(this)

  end subroutine global_interpolation_find_xyz

  !> Finds the corresponding r,s,t coordinates and redistributes the points to
  !! the owning rank in the correct global element as well as which process
  !! that owns the point.
  !! There will be some redundancy as for example xyz and xyz_local
  !! will have the same values
  !! After this the values at these points can be evaluated without MPI.
  !! @param xyz The coordinates of the points.
  !! @param n_points The number of points.
  subroutine global_interpolation_find_and_redist(this, xyz, n_points)
    class(global_interpolation_t), intent(inout) :: this
    integer, intent(inout) :: n_points
    !!Perhaps this should be kind dp
    real(kind=rp), allocatable, intent(inout) :: xyz(:,:)


    call this%free_points()
    call this%free_points_local()


    this%n_points = n_points
    call global_interpolation_init_point_arrays(this)

    !> make deep copy incase xyz goes out of scope or deallocated
    call copy(this%xyz, xyz, 3 * n_points)

    call global_interpolation_find_common(this)
    call this%free_points()
    this%n_points = this%n_points_local
    n_points = this%n_points_local
    call global_interpolation_init_point_arrays(this)
    if (allocated(xyz)) then
       deallocate(xyz)
    end if
    allocate(xyz(3,n_points))

    call copy(xyz, this%xyz_local, 3*n_points)
    call copy(this%rst, this%rst_local, 3*n_points)
    call copy(this%xyz, this%xyz_local, 3*n_points)
    this%pe_owner = this%pe_rank
    this%el_owner0 = this%el_owner0_local
    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_memcpy(this%el_owner0, this%el_owner0_d, &
            this%n_points, HOST_TO_DEVICE, sync = .true.)
    end if
    this%all_points_local = .true.

  end subroutine global_interpolation_find_and_redist

  !> Evalute the interpolated value in the points given a field
  !! @param interp_values Array of values in the given points.
  !! @param field Array of values used for interpolation.
  !! @param on_host If interpolation should be carried out on the host
  subroutine global_interpolation_evaluate(this, interp_values, field, on_host)
    class(global_interpolation_t), target, intent(inout) :: this
    real(kind=rp), intent(inout), target :: interp_values(this%n_points)
    real(kind=rp), intent(inout), target :: field(this%nelv*this%Xh%lxyz)
    logical, intent(in) :: on_host
    type(c_ptr) :: interp_d

    if (.not. this%all_points_local) then
       call this%local_interp%evaluate(this%temp_local%x, this%el_owner0_local, &
            field, this%nelv, on_host)
       if (NEKO_BCKND_DEVICE .eq. 1 .and. .not. on_host) then
          call device_memcpy(this%temp_local%x, this%temp_local%x_d, &
               this%n_points_local, DEVICE_TO_HOST, .true.)
       end if
       interp_values = 0.0_rp
       call this%glb_intrp_comm%sendrecv(this%temp_local%x, interp_values, &
            this%n_points_local, this%n_points)
       if (NEKO_BCKND_DEVICE .eq. 1 .and. .not. on_host) then
          interp_d = device_get_ptr(interp_values)
          call device_memcpy(interp_values, interp_d, &
               this%n_points, HOST_TO_DEVICE, .false.)
       end if
    else
       call this%local_interp%evaluate(interp_values, this%el_owner0_local, &
            field, this%nelv, on_host)
    end if

  end subroutine global_interpolation_evaluate


  !> Compares two sets of rst coordinates and checks whether rst2 is
  !! better than rst1 given a tolerance res1 and res2 are the
  !! distances to the interpolated xyz coordinate and true xyz coord
  !! for point 1 and 2 tol specifies the range for the rst coordinate
  !! to be within: (r,s,t) in (-1+tol,1+tol)^3.
  !! @param rst1 local coordinates for point 1
  !! @param rst2 local coordinates for point 2
  !! @param rst1 distance between xyz(rst1) and the true xyz coordinate
  !! @param rst2 distance between xyz(rst2) and the true xyz coordinate
  !! @param tol Tolerance for how much rst1 and rst2 can over/undershoot [-1,1]
  function rst_cmp(rst1, rst2,res1, res2, tol) result(rst2_better)
    real(kind=rp) :: rst1(3), res1(3)
    real(kind=rp) :: rst2(3), res2(3)
    real(kind=rp) :: tol
    logical :: rst2_better
    !If rst1 is invalid and rst2 is valid, take rst2
    ! If both invalidl, take smallest residual
    if (abs(rst1(1)) .gt. 1.0_xp+tol .or. &
         abs(rst1(2)) .gt. 1.0_xp+tol .or. &
         abs(rst1(3)) .gt. 1.0_xp+tol) then
       if (abs(rst2(1)) .le. 1.0_xp+tol .and. &
            abs(rst2(2)) .le. 1.0_xp+tol .and. &
            abs(rst2(3)) .le. 1.0_xp+tol) then
          rst2_better = .true.
       else
          rst2_better = (norm2(real(res2,xp)) .lt. norm2(real(res1,xp)))
       end if
    else
       !> Else we check rst2 is inside and has a smaller distance
       rst2_better = (norm2(real(res2,xp)) .lt. norm2(real(res1,xp)) .and.&
            abs(rst2(1)) .le. 1.0_xp+tol .and. &
            abs(rst2(2)) .le. 1.0_xp+tol .and. &
            abs(rst2(3)) .le. 1.0_xp+tol)
    end if
  end function rst_cmp

end module global_interpolation
