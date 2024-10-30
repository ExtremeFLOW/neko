! Copyright (c) 2020-2023, The Neko Authors
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
  use num_types, only: rp, dp
  use space, only: space_t
  use stack
  use dofmap, only: dofmap_t
  use mesh, only: mesh_t
  use logger, only: neko_log, LOG_SIZE
  use utils, only: neko_error, neko_warning
  use local_interpolation
  use point
  use comm
  use aabb, only: aabb_t
  use aabb_tree, only: aabb_tree_t
  use vector, only: vector_t
  use math, only: copy
  use neko_mpi_types
  use structs, only: array_ptr_t
  use, intrinsic :: iso_c_binding
  implicit none
  private
  !> Implements global interpolation for arbitrary points in the domain.
  type, public :: global_interpolation_t
     !> X coordinates from which to interpolate.
     type(array_ptr_t) :: x
     !> Y coordinates from which to interpolate.
     type(array_ptr_t) :: y
     !> Z coordinates from which to interpolate.
     type(array_ptr_t) :: z
     !> Geometric dimension of the simulation.
     integer :: gdim
     !> Number of elements.
     integer :: nelv
     !> Space.
     type(space_t), pointer :: Xh
     !> Components related to the points we want to evalute
     !> Number of points we want to evaluate
     integer :: n_points
     !> x,y,z coordinates
     real(kind=rp), allocatable :: xyz(:,:)
     !> r,s,t coordinates
     real(kind=rp), allocatable :: rst(:,:)
     !> List of owning processes.
     integer, allocatable :: pe_owner(:)
     !> List of owning elements.
     !! Note this is 0 indexed
     integer, allocatable :: el_owner0(:)
     type(c_ptr) :: el_owner0_d = c_null_ptr
     !> Map to original indices after communication
     integer, allocatable :: pt_ids(:)

     !> Local points (points in this ranks domain)
     integer :: n_points_local
     integer, allocatable :: el_owner0_local(:)
     type(c_ptr) :: el_owner0_local_d = c_null_ptr
     real(kind=rp), allocatable :: rst_local(:,:)
     real(kind=rp), allocatable :: xyz_local(:,:)
     !> Interpolator for local points.
     type(local_interpolator_t) :: local_interp
     !> If all points are local on this PE.
     !> Turns true if points are redistributed to their respective owners
     logical :: all_points_local = .false.
     !> Tolerance for Newton solve to find the correct rst coordinates.
     real(kind=rp) :: tol = NEKO_EPS*1e3_rp
     
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
     type(aabb_t), allocatable :: global_aabb(:)
     type(aabb_tree_t) :: global_aabb_tree
     !> Structure to find element candidates
     type(aabb_t), allocatable :: local_aabb(:)
     type(aabb_tree_t) :: local_aabb_tree
     !> Things for gather-scatter operation (sending interpolated values back and forth)
     !! Essentially the information in n_points_pe, 
     !! n_points_pe_local organized for mpisend/recv 
     !! send variables are related to the local points that should be sent to other ranks
     !! recv variables are related to the points this rank wants from other ranks
     logical, allocatable :: mpi_send_flag(:), mpi_recv_flag(:)
     integer :: n_send_pe, n_recv_pe
     integer, allocatable :: offset_send_pe(:), n_points_send_pe(:)
     integer, allocatable :: offset_recv_pe(:), n_points_recv_pe(:)
     integer, allocatable :: pe_send_id(:), pe_recv_id(:)
     type(MPI_request), allocatable :: mpi_send_request(:), mpi_recv_request(:)
     !> Working vectors for global interpolation
     type(vector_t) :: temp_local, temp

   contains
     !> Initialize the global interpolation object on a dofmap.
     procedure, pass(this) :: init_xyz => global_interpolation_init_xyz
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
     procedure, pass(this) :: find_points_xyz => global_interpolation_find_xyz
     generic :: find_points => find_points_xyz, find_points_coords
     !> Evaluate the value of the field in each point.
     procedure, pass(this) :: evaluate => global_interpolation_evaluate

     !> Generic constructor
     generic :: init => init_dof, init_xyz

  end type global_interpolation_t

contains

  !> Initialize the global interpolation object on a dofmap.
  !! @param dof Dofmap on which the interpolation is to be carried out.
  !! @param tol Tolerance for Newton iterations.
  subroutine global_interpolation_init_dof(this, dof, tol)
    class(global_interpolation_t), intent(inout) :: this
    type(dofmap_t), target :: dof
    real(kind=rp), optional :: tol

    ! NOTE: Passing dof%x(:,1,1,1), etc in init_xyz passes down the entire
    ! dof%x array and not a slice. It is done this way for
    ! this%x%ptr to point to dof%x (see global_interpolation_init_xyz).
    call this%init_xyz(dof%x(:,1,1,1), dof%y(:,1,1,1), dof%z(:,1,1,1), &
         dof%msh%gdim, dof%msh%nelv, dof%Xh, tol = tol)

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
  subroutine global_interpolation_init_xyz(this, x, y, z, gdim, nelv, Xh, tol)
    class(global_interpolation_t), intent(inout) :: this
    real(kind=rp), intent(in), target :: x(:)
    real(kind=rp), intent(in), target :: y(:)
    real(kind=rp), intent(in), target :: z(:)
    integer, intent(in) :: gdim
    integer, intent(in) :: nelv
    type(space_t), intent(in), target :: Xh
    real(kind=rp), intent(in), optional :: tol
    integer :: lx, ly, lz, max_pts_per_iter, ierr, i, id1, id2, n
    real(kind=dp), allocatable :: rank_xyz_max(:,:), rank_xyz_min(:,:)
    type(stack_i4_t) :: pe_candidates
    real(kind=dp) :: max_xyz(3), min_xyz(3), padding

    this%x%ptr => x
    this%y%ptr => y
    this%z%ptr => z
    this%gdim = gdim
    this%nelv = nelv
    this%Xh => Xh
    if (present(tol)) this%tol = tol

    ! Number of points to iterate on simultaneosuly
    lx = Xh%lx
    ly = Xh%ly
    lz = Xh%lz
    n = nelv * lx*ly*lz
    allocate(rank_xyz_max(3,pe_size))
    allocate(rank_xyz_min(3,pe_size))
    max_xyz = (/maxval(x(1:n)), maxval(y(1:n)), maxval(z(1:n))/)
    min_xyz = (/minval(x(1:n)), minval(y(1:n)), minval(z(1:n))/)
    call MPI_Allgather(max_xyz, 3, MPI_DOUBLE_PRECISION, &
    rank_xyz_max, 3, MPI_DOUBLE_PRECISION, NEKO_COMM, ierr) 
    call MPI_Allgather(min_xyz, 3, MPI_DOUBLE_PRECISION, &
    rank_xyz_min, 3, MPI_DOUBLE_PRECISION, NEKO_COMM, ierr) 
    
    allocate(this%global_aabb(pe_size))
    allocate(this%local_aabb(nelv))
    
    do i = 1, pe_size
       call this%global_aabb(i)%init(rank_xyz_min(:,i), rank_xyz_max(:,i))

    end do
    padding = 1e-2
    call this%global_aabb_tree%init(pe_size+1)
    call this%global_aabb_tree%build(this%global_aabb, padding)

    do i = 1, nelv
       id1 = lx*ly*lz*(i-1)
       id2 = lx*ly*lz*(i)
       max_xyz = (/maxval(this%x%ptr(id1:id2)), &
       maxval(this%y%ptr(id1:id2)), maxval(this%z%ptr(id1:id2))/)
       min_xyz = (/minval(this%x%ptr(id1:id2)), &
       minval(this%y%ptr(id1:id2)), minval(this%z%ptr(id1:id2))/)
       call this%local_aabb(i)%init(min_xyz, max_xyz)
    end do

    call this%local_aabb_tree%init(nelv)
    call this%local_aabb_tree%build(this%local_aabb, padding)

  end subroutine global_interpolation_init_xyz


  !> Destructor
  subroutine global_interpolation_free(this)
    class(global_interpolation_t), intent(inout) :: this

    nullify(this%x%ptr)
    nullify(this%y%ptr)
    nullify(this%z%ptr)
    nullify(this%Xh)

    this%nelv = 0
    this%gdim = 0

    call this%free_points()
    call this%free_points_local()
    call this%local_interp%free()
   
    if (allocated(this%mpi_send_flag)) deallocate(this%mpi_send_flag)
    if (allocated(this%offset_send_pe)) deallocate(this%offset_send_pe)
    if (allocated(this%n_points_send_pe)) deallocate(this%n_points_send_pe)
    if (allocated(this%pe_send_id)) deallocate(this%pe_send_id)
    if (allocated(this%mpi_send_request)) deallocate(this%mpi_send_request)

    if (allocated(this%mpi_recv_flag)) deallocate(this%mpi_recv_flag)
    if (allocated(this%offset_recv_pe)) deallocate(this%offset_recv_pe)
    if (allocated(this%n_points_recv_pe)) deallocate(this%n_points_recv_pe)
    if (allocated(this%pe_recv_id)) deallocate(this%pe_recv_id)
    if (allocated(this%mpi_recv_request)) deallocate(this%mpi_recv_request)


    call this%temp_local%free()
    call this%temp%free()
    

  end subroutine global_interpolation_free

  !> Destructor for point arrays
  subroutine global_interpolation_free_points(this)
    class(global_interpolation_t), intent(inout) :: this

    this%n_points = 0
    this%all_points_local = .false.

    if (allocated(this%xyz)) deallocate(this%xyz)
    if (allocated(this%rst)) deallocate(this%rst)
    if (allocated(this%pe_owner)) deallocate(this%pe_owner)
    if (allocated(this%el_owner0)) deallocate(this%el_owner0)

    if (c_associated(this%el_owner0_d)) then
       call device_free(this%el_owner0_d)
    end if
    if (allocated(this%pt_ids)) deallocate(this%pt_ids)

  end subroutine global_interpolation_free_points

  subroutine global_interpolation_free_points_local(this)
    class(global_interpolation_t), intent(inout) :: this

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
    class(global_interpolation_t), intent(inout) :: this
    !!Perhaps this should be kind dp
    real(kind=xp) :: xdiff, ydiff, zdiff
    character(len=8000) :: log_buf
    type(vector_t) :: x_check, x_vec
    type(vector_t) :: y_check, y_vec
    type(vector_t) :: z_check, z_vec
    real(kind=rp), allocatable :: x_t(:)
    real(kind=rp), allocatable :: y_t(:)
    real(kind=rp), allocatable :: z_t(:),  rst_local_cand(:,:)
    real(kind=rp), allocatable :: resx(:)
    real(kind=rp), allocatable :: resy(:)
    real(kind=rp), allocatable :: resz(:)
    real(kind=rp), allocatable :: rsts(:,:)
    real(kind=rp), allocatable :: res(:,:)
    logical :: isdiff
    real(kind=dp) :: pt_xyz(3), res1
    integer :: i, j, stupid_intent
    type(point_t), allocatable :: my_point(:)
    type(point_t), allocatable :: my_points(:)
    type(stack_i4_t) :: all_el_candidates
    type(stack_i4_t), allocatable :: points_at_pe(:)
    type(stack_i4_t), allocatable :: pe_candidates(:)
    type(stack_i4_t), allocatable :: el_candidates(:)
    integer, allocatable :: n_el_cands(:)
    integer, pointer :: pe_cands(:)
    integer, pointer :: el_cands(:)
    integer, pointer :: point_ids(:)
    real(kind=rp), allocatable :: xyz_send_to_pe(:,:)
    real(kind=rp), allocatable :: rst_send_to_pe(:,:)
    real(kind=rp), allocatable :: rst_recv_from_pe(:,:)
    real(kind=rp), allocatable :: res_recv_from_pe(:,:)
    integer, allocatable :: el_owner0s(:), el_send_to_pe(:)
    integer :: ierr, max_n_points_to_send, ii, n_point_cand, point_id, i_send, i_recv
    real(kind=rp) :: time1, time2, time_start


    call MPI_Barrier(NEKO_COMM)
    time_start = MPI_Wtime()
    write(log_buf,'(A)') 'Setting up global interpolation'
    call neko_log%message(log_buf)  
    ! Find pe candidates that the points i want may be at
    ! Add number to n_points_pe_local
    allocate(this%n_points_pe(0:(pe_size-1)))
    allocate(this%n_points_pe_local(0:(pe_size-1)))
    allocate(this%n_points_offset_pe_local(0:(pe_size-1)))
    allocate(this%n_points_offset_pe(0:(pe_size-1)))
    !Working arrays
    allocate(points_at_pe(0:(pe_size-1)))
    allocate(pe_candidates(this%n_points))
    allocate(my_point(this%n_points))
    this%n_points_pe = 0
    do i = 0, pe_size-1
       call points_at_pe(i)%init() 
    end do
    do i = 1, this%n_points
       pt_xyz = (/ this%xyz(1,i),this%xyz(2,i),this%xyz(3,i) /)
       call pe_candidates(i)%init() 
       call my_point(i)%init(pt_xyz) 
    end do
    do i = 1, this%n_points
       call this%global_aabb_tree%query_overlaps(my_point(i),-1, pe_candidates(i))
       pe_cands => pe_candidates(i)%array()
       do j = 1, pe_candidates(i)%size()
          this%n_points_pe(pe_cands(j)-1) = this%n_points_pe(pe_cands(j)-1) + 1
          stupid_intent = i
          call points_at_pe(pe_cands(j)-1)%push(stupid_intent)
       end do
       if (pe_candidates(i)%size() .lt. 1) then
          write (*,*) 'Warning, point', this%xyz(:,i), 'found to be outside domain'
       end if
    end do
    !Send number of points I want to candidates
    ! n_points_local -> how many points might be at this rank
    ! n_points_pe_local -> how many points I have that pe i might want
    this%n_points_pe_local = 0
    do i = 0, (pe_size - 1)
       call MPI_Reduce(this%n_points_pe(i), this%n_points_local, 1, MPI_INTEGER, &
            MPI_SUM, i, NEKO_COMM, ierr)
       !n_points_pe_local gives the number of points I will receive from every rank
       call MPI_Gather(this%n_points_pe(i), 1, MPI_INTEGER,&
                      this%n_points_pe_local, 1, MPI_INTEGER, i, NEKO_COMM, ierr)
    end do

    this%n_points_offset_pe_local(0) = 0
    this%n_points_offset_pe(0) = 0
    do i = 1, (pe_size - 1)
       this%n_points_offset_pe_local(i) = this%n_points_pe_local(i-1)&
                                 + this%n_points_offset_pe_local(i-1)
       this%n_points_offset_pe(i) = this%n_points_pe(i-1)&
                                 + this%n_points_offset_pe(i-1)
    end do
    
    allocate(this%xyz_local(3, this%n_points_local))
    max_n_points_to_send = max(maxval(this%n_points_pe),1)
    allocate(xyz_send_to_pe(3, max_n_points_to_send))
    do i = 0, (pe_size - 1)
       point_ids => points_at_pe(i)%array()
       do j = 1, this%n_points_pe(i)
          xyz_send_to_pe(:,j) = this%xyz(:,point_ids(j))
       end do
       call MPI_Gatherv(xyz_send_to_pe,3*this%n_points_pe(i), &
                        MPI_REAL_PRECISION, this%xyz_local,3*this%n_points_pe_local, &
                        3*this%n_points_offset_pe_local, &
                        MPI_REAL_PRECISION, i, NEKO_COMM, ierr)
    end do



    !Okay, now we need to find the rst...
    call all_el_candidates%init()
    allocate(el_candidates(this%n_points_local))
    allocate(my_points(this%n_points_local))
    do i = 1, this%n_points_local
       call el_candidates(i)%init()
       pt_xyz = (/ this%xyz_local(1,i),this%xyz_local(2,i),this%xyz_local(3,i) /)
       call my_points(i)%init(pt_xyz) 
    end do
    allocate(n_el_cands(this%n_points_local))
    do i = 1, this%n_points_local
       call this%local_aabb_tree%query_overlaps(my_points(i),-1, el_candidates(i))
       el_cands => el_candidates(i)%array()
       do j = 1, el_candidates(i)%size()
          stupid_intent = el_cands(j) -1
          call all_el_candidates%push(stupid_intent) !< OBS c indexing
       end do
       n_el_cands(i) = el_candidates(i)%size()
    end do

      
    n_point_cand = all_el_candidates%size()
    allocate(x_t(n_point_cand))
    allocate(y_t(n_point_cand))
    allocate(z_t(n_point_cand))
    ii = 0
    do i = 1 , this%n_points_local
       do j = 1, n_el_cands(i) 
          ii = ii + 1
          x_t(ii) = this%xyz_local(1,i)
          y_t(ii) = this%xyz_local(2,i)
          z_t(ii) = this%xyz_local(3,i)
       end do
    end do

    allocate(rst_local_cand(3,n_point_cand))
    allocate(resx(n_point_cand))
    allocate(resy(n_point_cand))
    allocate(resz(n_point_cand))
    allocate(this%rst_local(3,this%n_points_local))
    allocate(this%el_owner0_local(this%n_points_local))
    

    el_cands => all_el_candidates%array() 
    
    call MPI_Barrier(NEKO_COMM)
    time1 = MPI_Wtime()
    call find_rst_legendre(rst_local_cand, x_t, y_t, z_t, this%Xh, &
                           this%x%ptr, this%y%ptr, this%z%ptr, &
                           el_cands, n_point_cand, this%nelv, &
                           resx, resy, resz, this%tol)
 
    call MPI_Barrier(NEKO_COMM)
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

    ii = 0
    do i = 1 , this%n_points_local
       this%xyz_local(1,i) = 10.0
       this%xyz_local(2,i) = 10.0
       this%xyz_local(3,i) = 10.0
       this%rst_local(1,i) = 10.0
       this%rst_local(2,i) = 10.0
       this%rst_local(3,i) = 10.0
       do j = 1, n_el_cands(i) 
          ii = ii + 1
          if (rst_cmp(this%rst_local(:,i), rst_local_cand(:,ii),&
             this%xyz_local(:,i), (/resx(ii),resy(ii),resz(ii)/), this%tol)) then
             this%rst_local(1,i) = rst_local_cand(1,ii)
             this%rst_local(2,i) = rst_local_cand(2,ii)
             this%rst_local(3,i) = rst_local_cand(3,ii)
             this%xyz_local(1,i) = resx(ii)
             this%xyz_local(2,i) = resy(ii)
             this%xyz_local(3,i) = resz(ii)
             this%el_owner0_local(i) = el_cands(ii)
          end if
       end do
    end do
    allocate(rsts(3,this%n_points))
    allocate(res(3,this%n_points))
    allocate(rst_recv_from_pe(3, max_n_points_to_send))
    allocate(res_recv_from_pe(3, max_n_points_to_send))
    allocate(el_owner0s(max_n_points_to_send))
    res = 1e2
    this%rst = 1e2
    this%pe_owner = -1
    !> Send rst and res to rank who want this point  S
    do i = 0, (pe_size - 1)
       call MPI_Scatterv(this%rst_local,3*this%n_points_pe_local, &
                         3*this%n_points_offset_pe_local,&
                         MPI_REAL_PRECISION, rst_recv_from_pe,3*this%n_points_pe(i), &
                         MPI_REAL_PRECISION, i, NEKO_COMM, ierr)
       call MPI_Scatterv(this%xyz_local,3*this%n_points_pe_local, &
                         3*this%n_points_offset_pe_local,&
                         MPI_REAL_PRECISION, res_recv_from_pe,3*this%n_points_pe(i), &
                         MPI_REAL_PRECISION, i, NEKO_COMM, ierr)
       call MPI_Scatterv(this%el_owner0_local,this%n_points_pe_local, &
                         this%n_points_offset_pe_local,&
                         MPI_INTEGER, el_owner0s,this%n_points_pe(i), &
                         MPI_INTEGER, i, NEKO_COMM, ierr)
       point_ids => points_at_pe(i)%array()
       do j = 1, this%n_points_pe(i)
          point_id = point_ids(j)
          if (rst_cmp(this%rst(:,point_id), rst_recv_from_pe(:,j), &
                      res(:,point_id), res_recv_from_pe(:,j), this%tol)) then
             this%rst(:,point_ids(j)) = rst_recv_from_pe(:,j)
             res(:,point_ids(j)) = res_recv_from_pe(:,j)
             this%pe_owner(point_ids(j)) = i
             this%el_owner0(point_ids(j)) = el_owner0s(j)
          end if
       end do
    end do
     
    !OK, now I know the correct rst values 
    !of the points I want
    !Time to send the correct rst values to the correct ranks 
    
    do i = 0, pe_size-1
       call points_at_pe(i)%free() 
       call points_at_pe(i)%init() 
       this%n_points_pe(i) = 0
    end do
    
    do i = 1, this%n_points
       stupid_intent = i
       if (this%pe_owner(i) .eq. -1) print *, 'wth', stupid_intent, this%xyz(:,i)
       call points_at_pe(this%pe_owner(i))%push(stupid_intent)
       
       this%n_points_pe(this%pe_owner(i)) =  this%n_points_pe(this%pe_owner(i)) + 1
    end do

    do i = 0, (pe_size - 1)
       call MPI_Reduce(this%n_points_pe(i), this%n_points_local, 1, MPI_INTEGER, &
            MPI_SUM, i, NEKO_COMM, ierr)
       !n_points_pe_local gives the number of points I will receive from every rank
       call MPI_Gather(this%n_points_pe(i), 1, MPI_INTEGER,&
                      this%n_points_pe_local, 1, MPI_INTEGER, i, NEKO_COMM, ierr)
    end do
    this%n_points_offset_pe_local(0) = 0
    this%n_points_offset_pe(0) = 0
    do i = 1, (pe_size - 1)
       this%n_points_offset_pe_local(i) = this%n_points_pe_local(i-1)&
                                 + this%n_points_offset_pe_local(i-1)
       this%n_points_offset_pe(i) = this%n_points_pe(i-1)&
                                 + this%n_points_offset_pe(i-1)
    end do
    allocate(rst_send_to_pe(3, max_n_points_to_send))
    allocate(el_send_to_pe(max_n_points_to_send))
    allocate(this%pt_ids(this%n_points))
    ii = 0
    do i = 0, (pe_size - 1)
       point_ids => points_at_pe(i)%array()
       do j = 1, this%n_points_pe(i)
          ii = ii + 1
          xyz_send_to_pe(:,j) = this%xyz(:,point_ids(j))
          rst_send_to_pe(:,j) = this%rst(:,point_ids(j))
          el_send_to_pe(j) = this%el_owner0(point_ids(j))
          this%pt_ids(ii) = point_ids(j)
       end do
       call MPI_Gatherv(xyz_send_to_pe,3*this%n_points_pe(i), &
                        MPI_REAL_PRECISION, this%xyz_local, &
                        3*this%n_points_pe_local, &
                        3*this%n_points_offset_pe_local, &
                        MPI_REAL_PRECISION, i, NEKO_COMM, ierr)
       call MPI_Gatherv(rst_send_to_pe,3*this%n_points_pe(i), &
                        MPI_REAL_PRECISION, this%rst_local, &
                        3*this%n_points_pe_local, &
                        3*this%n_points_offset_pe_local, &
                        MPI_REAL_PRECISION, i, NEKO_COMM, ierr)
       call MPI_Gatherv(el_send_to_pe,this%n_points_pe(i), &
                        MPI_INTEGER, this%el_owner0_local,this%n_points_pe_local, &
                        this%n_points_offset_pe_local, &
                        MPI_INTEGER, i, NEKO_COMM, ierr)

    end do
    this%n_send_pe = 0 
    this%n_recv_pe = 0 
    do i = 0, (pe_size-1)
       if (this%n_points_pe(i) .gt. 0) then
          this%n_recv_pe = this%n_recv_pe + 1
       end if
       if (this%n_points_pe_local(i) .gt. 0) then
          this%n_send_pe = this%n_send_pe + 1
       end if
    end do
    allocate(this%offset_send_pe(this%n_send_pe), &
             this%n_points_send_pe(this%n_send_pe), &
             this%pe_send_id(this%n_send_pe), &
             this%mpi_send_flag(this%n_send_pe), &
             this%mpi_send_request(this%n_send_pe))
    allocate(this%offset_recv_pe(this%n_recv_pe), &
             this%n_points_recv_pe(this%n_recv_pe), &
             this%pe_recv_id(this%n_recv_pe), &
             this%mpi_recv_flag(this%n_recv_pe), &
             this%mpi_recv_request(this%n_recv_pe))
    i_send = 0
    i_recv = 0
    do i = 0, (pe_size-1)
       if (this%n_points_pe(i) .gt. 0) then
          i_recv = i_recv + 1
          this%n_points_recv_pe(i_recv) = this%n_points_pe(i)
          this%offset_recv_pe(i_recv) = this%n_points_offset_pe(i)
          this%pe_recv_id(i_recv) = i
       end if
       if (this%n_points_pe_local(i) .gt. 0) then
          i_send = i_send + 1
          this%n_points_send_pe(i_send) = this%n_points_pe_local(i)
          this%offset_send_pe(i_send) = this%n_points_offset_pe_local(i)
          this%pe_send_id(i_send) = i
       end if
    end do
    call this%temp_local%init(this%n_points_local)
    call this%temp%init(this%n_points)

    call x_check%init(this%n_points) 
    call y_check%init(this%n_points) 
    call z_check%init(this%n_points) 
  
    call this%local_interp%init(this%Xh, this%rst_local(1,:),&
                                this%rst_local(2,:), &
                                this%rst_local(3,:), this%n_points_local)
    call this%evaluate(x_check%x, this%x%ptr, .true.)
    call this%evaluate(y_check%x, this%y%ptr, .true.)
    call this%evaluate(z_check%x, this%z%ptr, .true.)
    
    do i = 1 , this%n_points

       ! Check validity of points
       isdiff = .false.
       xdiff = x_check%x(i)-this%xyz(1,i)
       ydiff = y_check%x(i)-this%xyz(2,i)
       zdiff = z_check%x(i)-this%xyz(3,i)
       isdiff = norm2(real((/xdiff,ydiff,zdiff/),xp)) > this%tol
       isdiff = isdiff .or. abs(this%rst(1,i)) > 1.0_xp + this%tol 
       isdiff = isdiff .or. abs(this%rst(2,i)) > 1.0_xp + this%tol 
       isdiff = isdiff .or. abs(this%rst(3,i)) > 1.0_xp + this%tol 
       if (isdiff ) then
          write(*,*) 'Point with coordinates: ', &
                this%xyz(1, i), this%xyz(2, i), this%xyz(3, i), &
                'Differ from interpolated coords: ', &
                x_check%x(i), y_check%x(i), z_check%x(i), &
                'Actual difference: ', &
                xdiff, ydiff, zdiff, norm2(real((/xdiff,ydiff,zdiff/),xp)),&
                'Expected difference: ', &
                res(:,i), norm2(real(res(:,i),xp)),&
                'Process, element: ', &
                this%pe_owner(i), this%el_owner0(i)+1, &
                'rst coords: ', &
                this%rst(:,i), &
                ' radius', sqrt(this%xyz(1,i)**2.0_xp+this%xyz(2,i)**2.0_xp)
       end if
    end do



    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_memcpy(this%el_owner0, this%el_owner0_d, &
            this%n_points, HOST_TO_DEVICE, sync = .true.)
       call device_map(this%el_owner0_local, this%el_owner0_local_d, this%n_points)
       call device_memcpy(this%el_owner0_local, this%el_owner0_local_d, &
            this%n_points, HOST_TO_DEVICE, sync = .true.)
    end if

    call MPI_Barrier(NEKO_COMM)
    time2 = MPI_Wtime()
    write(log_buf, '(A,E15.7)') 'Global interpolation setup done, time (s):', &
                                time2-time_start
    call neko_log%message(log_buf)  

  end subroutine global_interpolation_find_common

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

    this%n_points = n_points

    call global_interpolation_init_point_arrays(this)

    do i = 1, n_points
       this%xyz(1, i) = x(i,1,1,1)
       this%xyz(2, i) = y(i,1,1,1)
       this%xyz(3, i) = z(i,1,1,1)
    end do

    call global_interpolation_find_common(this)

  end subroutine global_interpolation_find_coords

  subroutine global_interpolation_init_point_arrays(this)
    class(global_interpolation_t) :: this

    allocate(this%xyz(3, this%n_points))
    allocate(this%rst(3, this%n_points))
    allocate(this%pe_owner(this%n_points))
    allocate(this%el_owner0(this%n_points))

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

    this%n_points = n_points

    call global_interpolation_init_point_arrays(this)

    !> make deep copy incase xyz goes out of scope or deallocated
    call copy(this%xyz, xyz, 3 * n_points)

    call global_interpolation_find_common(this)

  end subroutine global_interpolation_find_xyz

  !> Finds the corresponding r,s,t coordinates and redistributes the points to
  !! the owning rank in the correct global element as well as which process
  !! that owns the point.
  !! After this the values at these points can be evaluated without MPI.
  !! @param xyz The coordinates of the points.
  !! @param n_points The number of points.
  subroutine global_interpolation_find_and_redist(this, xyz, n_points)
    class(global_interpolation_t), intent(inout) :: this
    integer, intent(inout) :: n_points
    !!Perhaps this should be kind dp
    real(kind=rp), allocatable, intent(inout) :: xyz(:,:)
    integer :: i


    call this%free_points()

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
    this%pe_owner = pe_rank
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
  subroutine global_interpolation_evaluate(this, interp_values, field, on_host)
    class(global_interpolation_t), intent(inout) :: this
    real(kind=rp), intent(inout) :: interp_values(this%n_points)
    real(kind=rp), intent(inout) :: field(this%nelv*this%Xh%lxyz)
    logical, intent(in) :: on_host
    integer :: ierr, i
    real(kind=rp) :: time1, time2
    type(c_ptr) :: field_d, interp_d
    integer :: nreqs
    if (.not. this%all_points_local) then

       call this%local_interp%evaluate(this%temp_local%x, this%el_owner0_local, &
                                          field, this%nelv, on_host)
       if (NEKO_BCKND_DEVICE .eq. 1 .and. .not. on_host) then 
          call device_memcpy(this%temp_local%x, this%temp_local%x_d, &
               this%n_points_local, DEVICE_TO_HOST, .true.)
       end if
       !Post sends
       do i = 1, this%n_send_pe
          call MPI_Isend(this%temp_local%x(this%offset_send_pe(i)+1), &
                         this%n_points_send_pe(i), &
                         MPI_REAL_PRECISION, this%pe_send_id(i), 0, &
                         NEKO_COMM, this%mpi_send_request(i), ierr)  
          this%mpi_send_flag(i) = .false.
       end do
       !Post receives
       do i = 1, this%n_recv_pe
          call MPI_Irecv(this%temp%x(this%offset_recv_pe(i)+1), &
                         this%n_points_recv_pe(i), &
                         MPI_REAL_PRECISION, this%pe_recv_id(i), 0, &
                         NEKO_COMM, this%mpi_recv_request(i), ierr) 
           this%mpi_recv_flag(i) = .false.
       end do
       !> Check if I have my stuff
       nreqs = 0
       do while (nreqs .lt. this%n_recv_pe)
          do i = 1, this%n_recv_pe
             if (.not. this%mpi_recv_flag(i)) then
                call MPI_Test(this%mpi_recv_request(i), &
                              this%mpi_recv_flag(i), &
                              MPI_STATUS_IGNORE, ierr)
                if (this%mpi_recv_flag(i)) nreqs = nreqs + 1
             end if
          end do
       end do
       !Fix so my points go back to their ortiginal id
       do i = 1, this%n_points
          interp_values(this%pt_ids(i)) = this%temp%x(i)
       end do
       if (NEKO_BCKND_DEVICE .eq. 1 .and. .not. on_host) then
          interp_d = device_get_ptr(interp_values)
          call device_memcpy(interp_values, interp_d, &
               this%n_points, HOST_TO_DEVICE, .false.)
       end if  
       !> Lastly, check that my sends are ok
       nreqs = 0
       do while (nreqs .lt. this%n_send_pe)
          do i = 1, this%n_send_pe
             if (.not. this%mpi_send_flag(i)) then
                call MPI_Test(this%mpi_send_request(i), &
                              this%mpi_send_flag(i), &
                              MPI_STATUS_IGNORE, ierr) 
                if (this%mpi_send_flag(i)) nreqs = nreqs + 1
             end if
          end do
       end do
    else 
       call this%local_interp%evaluate(interp_values, this%el_owner0_local, &
                                       field, this%nelv, on_host)
    end if
 
  end subroutine global_interpolation_evaluate

end module global_interpolation
