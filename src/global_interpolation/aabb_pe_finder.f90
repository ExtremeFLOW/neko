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
!> Implements aabb_pe_finder given a dofmap.
!!
module aabb_pe_finder
   use num_types, only: rp, dp, xp
   use neko_config, only : NEKO_BCKND_DEVICE
   use space, only: space_t
   use stack, only: stack_i4_t, stack_i4t2_t
   use utils, only: neko_error, neko_warning
   use tuple, only: tuple_i4_t
   use htable, only: htable_i4_t
   use point, only: point_t
   use comm, only: NEKO_COMM, MPI_REAL_PRECISION, pe_rank, pe_size
   use mpi_f08, only: MPI_SUM, MPI_Reduce, MPI_COMM, MPI_Comm_rank, &
      MPI_Comm_size, MPI_Wtime, MPI_INTEGER, MPI_IN_PLACE, &
      MPI_MIN, MPI_Allgather, MPI_Barrier, MPI_DOUBLE_PRECISION
   use aabb, only: aabb_t
   use aabb_tree, only: aabb_tree_t, aabb_node_t, AABB_NULL_NODE
   use vector, only: vector_t
   use matrix, only: matrix_t
   use tensor, only: tnsr3d
   use math, only: copy, glsum, NEKO_M_LN2, NEKO_EPS
   use structs, only : array_ptr_t
   use, intrinsic :: iso_c_binding, only: c_ptr, c_null_ptr, c_associated, &
      c_sizeof, c_bool, c_loc
   implicit none
   private

   !> Minimum number of total boxes in the aabb tree
   integer, public, parameter :: GLOB_MAP_SIZE = 4096

   !> Implements global interpolation for arbitrary points in the domain.
   type, public :: aabb_pe_finder_t
      !> Which communicator to find things on
      type(MPI_COMM) :: comm
      !> pe_rank in comm
      integer :: pe_rank
      !> pe_size of comm
      integer :: pe_size
      real(kind=dp) :: padding
      !> Structure to find rank candidates
      type(aabb_t), allocatable :: global_aabb(:)
      type(aabb_tree_t) :: global_aabb_tree
      !> Number of boxes per rank
      integer :: pe_box_num
      !> Number of boxes in total
      integer :: glob_map_size

   contains
      procedure, pass(this) :: init => aabb_pe_finder_init
      procedure, pass(this) :: free => aabb_pe_finder_free
      procedure, pass(this) :: find => aabb_pe_finder_find_candidates
      procedure, pass(this) :: find_batch => aabb_pe_finder_find_candidates_batch

   end type aabb_pe_finder_t

contains

   !> Initialize the global interpolation object on a set of coordinates.
   !! @param x x-coordinates.
   !! @param y y-coordinates.
   !! @param z z-coordinates.
   !! @param gdim Geometric dimension.
   !! @param nelv Number of elements of the mesh in which to search for the
   !! points.
   !! @param Xh Space on which to interpolate.
   !! @param tol Tolerance for Newton iterations.
   subroutine aabb_pe_finder_init(this, x, y, z, nelv, Xh, comm, padding)
      class(aabb_pe_finder_t), intent(inout) :: this
      real(kind=rp), intent(in), target :: x(:)
      real(kind=rp), intent(in), target :: y(:)
      real(kind=rp), intent(in), target :: z(:)
      integer, intent(in) :: nelv
      type(MPI_COMM), intent(in), optional :: comm
      type(space_t), intent(in), target :: Xh
      real(kind=dp), intent(in) :: padding
      integer :: lx, ly, lz, max_pts_per_iter, ierr, i, id1, id2, n, j
      real(kind=dp), allocatable :: rank_xyz_max(:,:), rank_xyz_min(:,:), max_xyz(:,:), min_xyz(:,:)
      type(stack_i4_t) :: pe_candidates
      type(stack_i4t2_t) :: traverse_stack
      type(aabb_tree_t) :: local_aabb_tree
      type(aabb_t), allocatable :: local_aabb(:)
      integer :: start, last, n_box_el, lvl
      type(tuple_i4_t) :: id_lvl, temp_id_lvl
      real(kind=rp) :: time1, time_start
      type(aabb_node_t) :: node

      call this%free()

      this%comm = comm
      this%padding = padding

      call MPI_Comm_rank(this%comm, this%pe_rank, ierr)
      call MPI_Comm_size(this%comm, this%pe_size, ierr)


      time_start = MPI_Wtime()
      lx = Xh%lx
      ly = Xh%ly
      lz = Xh%lz
      n = nelv * lx*ly*lz
      !> Create a local tree for each element at this rank
      call local_aabb_tree%init(nelv)
      if (allocated(local_aabb)) deallocate(local_aabb)
      allocate(local_aabb(nelv))

      do i = 1, nelv
         id1 = lx*ly*lz*(i-1)+1
         id2 = lx*ly*lz*(i)
         call local_aabb(i)%init( (/minval(x(id1:id2)), &
            minval(y(id1:id2)), &
            minval(z(id1:id2))/), &
            (/maxval(x(id1:id2)), &
            maxval(y(id1:id2)), &
            maxval(z(id1:id2))/))
         call local_aabb_tree%insert_object(local_aabb(i),i)
      end do




      this%pe_box_num = min(GLOB_MAP_SIZE/this%pe_size,nelv)
      this%pe_box_num = ishft(1, ceiling(log(real(this%pe_box_num, rp)) / NEKO_M_LN2))
      call MPI_Allreduce(MPI_IN_PLACE, this%pe_box_num, 1, MPI_INTEGER, &
         MPI_MIN, this%comm, ierr)
      this%pe_box_num = max(this%pe_box_num,2)
      this%glob_map_size = this%pe_box_num*this%pe_size
      if (pe_rank .eq. 0) print *, this%pe_box_num, this%glob_map_size, 'hey'
      allocate(rank_xyz_max(3,this%glob_map_size))
      allocate(rank_xyz_min(3,this%glob_map_size))
      allocate(min_xyz(3,this%pe_box_num))
      allocate(max_xyz(3,this%pe_box_num))
      i = 1
      id_lvl = (/local_aabb_tree%get_root_index(),0/)
      call traverse_stack%init()
      call traverse_stack%push(id_lvl)
      lvl = 0
      !> Traverse the local tree and find ther top boxes
      do while (traverse_stack%size() > 0)
         id_lvl = traverse_stack%pop()
         lvl = id_lvl%x(2)
         node = local_aabb_tree%get_node(id_lvl%x(1))
         if (2**lvl == this%pe_box_num .or. node%is_leaf()) then
            min_xyz(:,i) = node%aabb%get_min()
            max_xyz(:,i) = node%aabb%get_max()
            i = i + 1
         else if (2**lvl < this%pe_box_num) then
            if (node%get_left_index() .ne. AABB_NULL_NODE) then
               temp_id_lvl = (/node%get_left_index(),lvl+1/)
               call traverse_stack%push(temp_id_lvl)
            end if
            if (node%get_right_index() .ne. AABB_NULL_NODE) then
               temp_id_lvl = (/node%get_right_index(),lvl+1/)
               call traverse_stack%push(temp_id_lvl)
            end if
         end if
      end do
      !Needs to be something in the domain
      !If somehow we dont need all boxes we just put a point here
      do j = i, this%pe_box_num
         min_xyz(:,j) =[x(1), y(1), z(1)]
         max_xyz(:,j) =[x(1), y(1), z(1)]
      end do
      !> Get boxes from all ranks
      call MPI_Allgather(max_xyz, 3*this%pe_box_num, MPI_DOUBLE_PRECISION, &
         rank_xyz_max, 3*this%pe_box_num, MPI_DOUBLE_PRECISION, this%comm, ierr)
      call MPI_Allgather(min_xyz, 3*this%pe_box_num, MPI_DOUBLE_PRECISION, &
         rank_xyz_min, 3*this%pe_box_num, MPI_DOUBLE_PRECISION, this%comm, ierr)
      call MPI_Barrier(this%comm)
      if (allocated(this%global_aabb)) deallocate(this%global_aabb)
      allocate(this%global_aabb(this%glob_map_size))
      !> Create global tree for each rank
      do i = 1, this%glob_map_size
         call this%global_aabb(i)%init(rank_xyz_min(:,i), rank_xyz_max(:,i))
      end do
      call this%global_aabb_tree%build(this%global_aabb,padding)
   end subroutine aabb_pe_finder_init


   !> Destructor
   subroutine aabb_pe_finder_free(this)
      class(aabb_pe_finder_t), intent(inout) :: this

      if (allocated(this%global_aabb)) deallocate(this%global_aabb)

   end subroutine aabb_pe_finder_free

   !> Find pe candidates
   !! @param my_point Point to find candidates for.
   !! @param pe_candidates Candidates for the point.
   subroutine aabb_pe_finder_find_candidates(this, my_point, pe_candidates)
      class(aabb_pe_finder_t), intent(inout) :: this
      type(point_t), intent(inout) :: my_point
      type(stack_i4_t), intent(inout) :: pe_candidates
      integer, pointer :: pe_cands(:) => Null()
      integer :: i

      call this%global_aabb_tree%query_overlaps(my_point, -1, pe_candidates)
      pe_cands => pe_candidates%array()
      do i = 1, pe_candidates%size()
         pe_cands(i) = (pe_cands(i)-1)/this%pe_box_num
      end do

   end subroutine aabb_pe_finder_find_candidates

   subroutine aabb_pe_finder_find_candidates_batch(this, points, n_points, points_at_pe, n_points_pe)
      class(aabb_pe_finder_t), intent(inout) :: this
      integer, intent(in) :: n_points
      real(kind=rp), intent(inout) :: points(3,n_points)
      type(stack_i4_t) :: points_at_pe(0:(this%pe_size-1))
      integer, intent(inout) :: n_points_pe(0:(this%pe_size-1))
      type(stack_i4_t) :: pe_candidates
      type(point_t) :: my_point
      integer :: i, j, temp_intent, pe_id, htable_data
      real(kind=rp) :: pt_xyz(3)
      integer, pointer :: pe_cands(:) => Null()
      type(htable_i4_t) :: marked_rank

      do i = 0, this%pe_size-1
         call points_at_pe(i)%clear()
         n_points_pe(i) = 0
      end do
      
      call marked_rank%init(32,htable_data)
      call pe_candidates%init()

      !> Check which ranks might have this point
      do i = 1, n_points
         !Should probably be a htable instead
         call marked_rank%clear()
         pt_xyz = (/ points(1,i),points(2,i),points(3,i) /)
         call my_point%init(pt_xyz)
         call pe_candidates%clear()
         call this%find(my_point, pe_candidates)
         pe_cands => pe_candidates%array()
         do j = 1, pe_candidates%size()
            pe_id = pe_cands(j)
            ! Check if this rank has already been marked 
            if (marked_rank%get(pe_id,htable_data) .ne. 0) then
               n_points_pe(pe_id) = n_points_pe(pe_id) + 1
               temp_intent = i
               call points_at_pe(pe_id)%push(temp_intent)
               ! Doesnt matter, I use htable as a set
               htable_data = 100
               call marked_rank%set(pe_id, htable_data)
            end if
         end do
         if (pe_candidates%size() .lt. 1) then
            write (*,*) 'WARNING, point', points(:,i), &
               'found to be outside domain, try increasing the padding to find rank candidates'
         end if
      end do
      call marked_rank%free()
      
   end subroutine aabb_pe_finder_find_candidates_batch

   !!> Common routine for finding the points.
   !subroutine aabb_pe_finder_find(xyz, n_points, pe_candidates, n_cands_per_point, point_cand_offset)
   !   class(aabb_pe_finder_t), intent(inout) :: this
   !   integer, intent(in) :: n_points
   !   real(kind=rp), intent(in) :: xyz(3,n_points)
   !   integer, intent(out), allocatable :: n_cands_per_point
   !   integer, intent(out), allocatable :: point_cand_offset

   !   !!Perhaps this should be kind dp
   !   real(kind=xp) :: xdiff, ydiff, zdiff
   !   character(len=8000) :: log_buf
   !   type(vector_t) :: x_check, x_vec
   !   type(vector_t) :: y_check, y_vec
   !   type(vector_t) :: z_check, z_vec
   !   type(vector_t) :: x_t
   !   type(vector_t) :: y_t
   !   type(vector_t) :: z_t
   !   type(matrix_t) :: rst_local_cand
   !   type(vector_t) :: resx
   !   type(vector_t) :: resy
   !   type(vector_t) :: resz
   !   type(vector_t) :: x_hat, y_hat, z_hat
   !   logical(kind=c_bool), allocatable, target :: conv_pts(:)
   !   type(c_ptr) :: conv_pts_d = c_null_ptr
   !   type(c_ptr) :: el_cands_d = c_null_ptr
   !   type(c_ptr) :: null_ptr = c_null_ptr
   !   real(kind=rp), allocatable :: res(:,:)
   !   logical :: isdiff
   !   real(kind=dp) :: pt_xyz(3), res1
   !   integer :: i, j, stupid_intent, iter
   !   integer(kind=8) :: bytes
   !   type(point_t) :: my_point
   !   type(stack_i4_t) :: all_el_candidates
   !   type(stack_i4_t), allocatable :: points_at_pe(:)
   !   type(stack_i4_t) :: pe_candidates, temp_stack
   !   type(stack_i4_t) :: el_candidates
   !   integer, allocatable :: n_el_cands(:)
   !   integer, pointer :: pe_cands(:) => Null()
   !   integer, pointer :: el_cands(:) => Null()
   !   integer, pointer :: point_ids(:) => NUll()
   !   integer, pointer :: send_recv(:) => NUll()
   !   real(kind=rp), allocatable :: xyz_send_to_pe(:,:)
   !   real(kind=rp), allocatable :: rst_send_to_pe(:,:)
   !   real(kind=rp), allocatable :: rst_recv_from_pe(:,:)
   !   real(kind=rp), allocatable :: res_recv_from_pe(:,:)
   !   real(kind=rp), allocatable :: res_results(:,:)
   !   real(kind=rp), allocatable :: rst_results(:,:)
   !   integer, allocatable :: el_owner0s(:), el_send_to_pe(:), el_owner_results(:)
   !   integer :: ierr, max_n_points_to_send, ii, n_point_cand, n_glb_point_cand, point_id, rank
   !   real(kind=rp) :: time1, time2, time_start
   !   logical :: converged
   !   logical, allocatable :: marked_rank(:)
   !   !Temp stuff for gs_comm
   !   type(stack_i4_t) :: send_pe, recv_pe
   !   type(gs_mpi_t) :: gs_find, gs_find_back
   !   type(stack_i4_t) :: send_pe_find, recv_pe_find


   !   call gs_find%init_dofs()
   !   call send_pe_find%init()
   !   call recv_pe_find%init()
   !   call MPI_Barrier(this%comm)
   !   time_start = MPI_Wtime()
   !   write(log_buf,'(A)') 'Setting up global interpolation'
   !   call neko_log%message(log_buf)
   !   ! Find pe candidates that the points i want may be at
   !   ! Add number to n_points_pe_local
   !   if (allocated(this%n_points_pe)) deallocate(this%n_points_pe)
   !   if (allocated(this%n_points_pe_local)) deallocate(this%n_points_pe_local)
   !   if (allocated(this%n_points_offset_pe_local)) deallocate(this%n_points_offset_pe_local)
   !   if (allocated(this%n_points_offset_pe)) deallocate(this%n_points_offset_pe)
   !   allocate(this%n_points_pe(0:(this%pe_size-1)))
   !   allocate(marked_rank(0:(this%pe_size-1)))
   !   allocate(this%n_points_offset_pe(0:(this%pe_size-1)))
   !   allocate(this%n_points_pe_local(0:(this%pe_size-1)))
   !   allocate(this%n_points_offset_pe_local(0:(this%pe_size-1)))
   !   !Working arrays
   !   allocate(points_at_pe(0:(this%pe_size-1)))
   !   this%n_points_pe = 0
   !   do i = 0, this%pe_size-1
   !      call points_at_pe(i)%init()
   !   end do
   !   call pe_candidates%init()
   !   call temp_stack%init()
   !   !> Check which ranks might have this point
   !   do i = 1, this%n_points
   !      marked_rank = .false.
   !      pt_xyz = (/ this%xyz(1,i),this%xyz(2,i),this%xyz(3,i) /)
   !      call my_point%init(pt_xyz)
   !      call pe_candidates%clear()
   !      pe_cands => pe_candidates%array()
   !      do j = 1, pe_candidates%size()
   !         !rank = pe_cands(j)
   !         rank = (pe_cands(j)-1)/this%pe_box_num
   !         if (.not. marked_rank(rank)) then
   !            this%n_points_pe(rank) = this%n_points_pe(rank) + 1
   !            stupid_intent = i
   !            call points_at_pe(rank)%push(stupid_intent)
   !            marked_rank(rank) = .true.
   !         end if
   !      end do
   !      if (pe_candidates%size() .lt. 1) then
   !         write (*,*) 'WARNING, point', this%xyz(:,i), &
   !            'found to be outside domain, something is likely very wrong'
   !      end if
   !   end do
   !   call MPI_Barrier(this%comm)
   !   time1 = MPI_Wtime()
   !   write(log_buf, '(A,E15.7)') &
   !      'GPU Found PE candidates time since start of findpts (s):', time1-time_start
   !   call neko_log%message(log_buf)
   !   !Send number of points I want to candidates
   !   ! n_points_local -> how many points might be at this rank
   !   ! n_points_pe_local -> how many points local on this rank that other pes might want
   !   this%n_points_pe_local = 0
   !   this%n_points_local = 0
   !   call MPI_Reduce_scatter_block(this%n_points_pe, this%n_points_local, 1, MPI_INTEGER, &
   !      MPI_SUM, this%comm, ierr)
   !   call MPI_Alltoall(this%n_points_pe, 1, MPI_INTEGER,&
   !      this%n_points_pe_local, 1, MPI_INTEGER, this%comm, ierr)
   !   !Set up offset arrays
   !   this%n_points_offset_pe_local(0) = 0
   !   this%n_points_offset_pe(0) = 0
   !   do i = 1, (this%pe_size - 1)
   !      this%n_points_offset_pe_local(i) = this%n_points_pe_local(i-1)&
   !         + this%n_points_offset_pe_local(i-1)
   !      this%n_points_offset_pe(i) = this%n_points_pe(i-1)&
   !         + this%n_points_offset_pe(i-1)
   !   end do
   !   do i = 0, (this%pe_size-1)
   !      if (this%n_points_pe(i) .gt. 0) then
   !         call send_pe_find%push(i)
   !         point_ids => points_at_pe(i)%array()
   !         do j = 1, this%n_points_pe(i)
   !            call gs_find%send_dof(i)%push(3*(point_ids(j)-1)+1)
   !            call gs_find%send_dof(i)%push(3*(point_ids(j)-1)+2)
   !            call gs_find%send_dof(i)%push(3*(point_ids(j)-1)+3)
   !         end do
   !      end if
   !      if (this%n_points_pe_local(i) .gt. 0) then
   !         call recv_pe_find%push(i)
   !         do j = 1, this%n_points_pe_local(i)
   !            call gs_find%recv_dof(i)%push(3*(j+this%n_points_offset_pe_local(i)-1)+1)
   !            call gs_find%recv_dof(i)%push(3*(j+this%n_points_offset_pe_local(i)-1)+2)
   !            call gs_find%recv_dof(i)%push(3*(j+this%n_points_offset_pe_local(i)-1)+3)
   !         end do
   !      end if
   !   end do



   !   call gs_find%init(send_pe_find, recv_pe_find, this%comm)

   !   call gs_find_back%init_dofs()
   !   ii = 0
   !   do i = 0, (this%pe_size-1)
   !      send_recv => gs_find%recv_dof(i)%array()
   !      do j = 1, gs_find%recv_dof(i)%size()
   !         call gs_find_back%send_dof(i)%push(send_recv(j))
   !      end do
   !      send_recv => gs_find%send_dof(i)%array()
   !      do j = 1, gs_find%send_dof(i)%size()
   !         ii = ii + 1
   !         call gs_find_back%recv_dof(i)%push(ii)
   !      end do
   !   end do

   !   call gs_find_back%init(recv_pe_find, send_pe_find, this%comm)


   !   if (allocated(this%xyz_local)) deallocate(this%xyz_local)
   !   allocate(this%xyz_local(3, this%n_points_local))
   !   max_n_points_to_send = max(maxval(this%n_points_pe),1)
   !   allocate(xyz_send_to_pe(3, max_n_points_to_send))
   !   this%xyz_local = 0.0_rp
   !   call gs_find%nbrecv()
   !   call gs_find%nbsend(this%xyz, this%n_points*3, null_ptr, null_ptr)
   !   call gs_find%nbwait(this%xyz_local, this%n_points_local*3, GS_OP_ADD, null_ptr)

   !   call MPI_Barrier(this%comm)
   !   time1 = MPI_Wtime()
   !   write(log_buf, '(A,E15.7)') &
   !      'Sent to points to PE candidates, time since start of find_points (s):', time1-time_start
   !   call neko_log%message(log_buf)
   !   !Okay, now we need to find the rst...
   !   call all_el_candidates%init()
   !   call el_candidates%init()
   !   allocate(n_el_cands(this%n_points_local))
   !   !> Find element candidates at this rank
   !   do i = 1, this%n_points_local
   !      pt_xyz = (/ this%xyz_local(1,i),this%xyz_local(2,i),this%xyz_local(3,i) /)
   !      call my_point%init(pt_xyz)
   !      call el_candidates%clear()
   !      call this%local_aabb_tree%query_overlaps(my_point,-1, el_candidates)
   !      el_cands => el_candidates%array()
   !      do j = 1, el_candidates%size()
   !         stupid_intent = el_cands(j) - 1
   !         call all_el_candidates%push(stupid_intent) !< OBS c indexing
   !      end do
   !      n_el_cands(i) = el_candidates%size()
   !   end do


   !   n_point_cand = all_el_candidates%size()
   !   call x_t%init(n_point_cand)
   !   call y_t%init(n_point_cand)
   !   call z_t%init(n_point_cand)
   !   ii = 0
   !   !> Copy xyz coords to each element candidate
   !   do i = 1 , this%n_points_local
   !      do j = 1, n_el_cands(i)
   !         ii = ii + 1
   !         x_t%x(ii) = this%xyz_local(1,i)
   !         y_t%x(ii) = this%xyz_local(2,i)
   !         z_t%x(ii) = this%xyz_local(3,i)
   !      end do
   !   end do

   !   call MPI_Barrier(this%comm)
   !   time1 = MPI_Wtime()
   !   write(log_buf, '(A,E15.7)') &
   !      'Element candidates found, now time for finding rst,time since start of find_points (s):', time1-time_start
   !   call neko_log%message(log_buf)
   !   call rst_local_cand%init(3,n_point_cand)
   !   call resx%init(n_point_cand)
   !   call resy%init(n_point_cand)
   !   call resz%init(n_point_cand)

   !   if (allocated(this%rst_local)) deallocate(this%rst_local)
   !   if (allocated(this%el_owner0_local)) deallocate(this%el_owner0_local)
   !   allocate(this%rst_local(3,this%n_points_local))
   !   allocate(this%el_owner0_local(this%n_points_local))

   !   ! Find rst within all element candidates for target xyz (x_t, y_t, z_t)
   !   call MPI_Barrier(this%comm)
   !   time1 = MPI_Wtime()
   !   el_cands => all_el_candidates%array()
   !   if ( NEKO_BCKND_DEVICE .ne. 1) then
   !      call find_rst_legendre(rst_local_cand%x, x_t%x, y_t%x, z_t%x, this%Xh, &
   !         this%x%ptr, this%y%ptr, this%z%ptr, &
   !         el_cands, n_point_cand, this%nelv, &
   !         resx%x, resy%x, resz%x, this%tol)


   !   end if
   !   if (NEKO_BCKND_DEVICE .eq. 1 .and. n_point_cand .gt. 0) then
   !      ! Initialize working arrays
   !      call x_hat%init(this%nelv*this%Xh%lxyz)
   !      call y_hat%init(this%nelv*this%Xh%lxyz)
   !      call z_hat%init(this%nelv*this%Xh%lxyz)
   !      call x_t%copyto(HOST_TO_DEVICE,.false.)
   !      call y_t%copyto(HOST_TO_DEVICE,.false.)
   !      call z_t%copyto(HOST_TO_DEVICE,.false.)

   !      call tnsr3d(x_hat%x, this%Xh%lx, this%x%ptr, &
   !         this%Xh%lx, this%Xh%vinv, &
   !         this%Xh%vinvt, this%Xh%vinvt, this%nelv)
   !      call tnsr3d(y_hat%x, this%Xh%lx, this%y%ptr, &
   !         this%Xh%lx, this%Xh%vinv, &
   !         this%Xh%vinvt, this%Xh%vinvt, this%nelv)
   !      call tnsr3d(z_hat%x, this%Xh%lx,this%z%ptr, &
   !         this%Xh%lx, this%Xh%vinv, &
   !         this%Xh%vinvt, this%Xh%vinvt, this%nelv)
   !      allocate(conv_pts(n_point_cand))
   !      conv_pts = .false.
   !      bytes = n_point_cand*c_sizeof(conv_pts(1))
   !      call device_alloc(conv_pts_d,bytes)
   !      call device_memcpy(c_loc(conv_pts), conv_pts_d, bytes, HOST_TO_DEVICE, .false., glb_cmd_queue)
   !      call device_map(el_cands, el_cands_d,n_point_cand)
   !      call device_memcpy(el_cands, el_cands_d,n_point_cand,HOST_TO_DEVICE, .true.)
   !      iter = 0
   !      converged = .false.
   !      rst_local_cand = 0.0_rp
   !      !Iterate until found, not heavily optimized
   !      do while (.not. converged)
   !         call device_find_rst_legendre(rst_local_cand%x_d, x_t%x_d, y_t%x_d, z_t%x_d, &
   !            x_hat%x_d, y_hat%x_d, z_hat%x_d, &
   !            resx%x_d, resy%x_d, resz%x_d, &
   !            this%Xh%lx,el_cands_d, n_point_cand, this%tol, &
   !            conv_pts_d)
   !         call device_memcpy(c_loc(conv_pts), conv_pts_d, bytes, DEVICE_TO_HOST, .true., glb_cmd_queue)
   !         converged = .true.
   !         iter = iter + 1
   !         do i = 1, n_point_cand
   !            converged = converged .and. conv_pts(i)
   !         end do
   !         if( iter .ge. 50) converged = .true.
   !      end do
   !      call rst_local_cand%copyto(DEVICE_TO_HOST,.false.)
   !      call resx%copyto(DEVICE_TO_HOST,.false.)
   !      call resy%copyto(DEVICE_TO_HOST,.false.)
   !      call resz%copyto(DEVICE_TO_HOST,.true.)
   !      call device_deassociate(el_cands)
   !      call device_free(el_cands_d)
   !      call device_free(conv_pts_d)
   !   end if

   !   call MPI_Barrier(this%comm)
   !   time2 = MPI_Wtime()
   !   write(log_buf, '(A,E15.7)') &
   !      'GPU Found rst with Newton iteration, time (s):', time2-time1
   !   call neko_log%message(log_buf)

   !   write(log_buf,'(A,E15.7)') &
   !      'Tolerance: ', this%tol
   !   call neko_log%message(log_buf)
   !   write(log_buf,'(A)') &
   !      'Checking validity of points and choosing best candidates.'
   !   call neko_log%message(log_buf)

   !   ! Choose the best candidate at this rank
   !   ii = 0
   !   do i = 1 , this%n_points_local
   !      this%xyz_local(1,i) = 10.0
   !      this%xyz_local(2,i) = 10.0
   !      this%xyz_local(3,i) = 10.0
   !      this%rst_local(1,i) = 10.0
   !      this%rst_local(2,i) = 10.0
   !      this%rst_local(3,i) = 10.0
   !      do j = 1, n_el_cands(i)
   !         ii = ii + 1
   !         if (rst_cmp(this%rst_local(:,i), rst_local_cand%x(:,ii),&
   !            this%xyz_local(:,i), (/resx%x(ii),resy%x(ii),resz%x(ii)/), this%tol)) then
   !            this%rst_local(1,i) = rst_local_cand%x(1,ii)
   !            this%rst_local(2,i) = rst_local_cand%x(2,ii)
   !            this%rst_local(3,i) = rst_local_cand%x(3,ii)
   !            this%xyz_local(1,i) = resx%x(ii)
   !            this%xyz_local(2,i) = resy%x(ii)
   !            this%xyz_local(3,i) = resz%x(ii)
   !            this%el_owner0_local(i) = el_cands(ii)
   !         end if
   !      end do
   !   end do
   !   allocate(res(3,this%n_points))
   !   allocate(rst_recv_from_pe(3, max_n_points_to_send))
   !   allocate(res_recv_from_pe(3, max_n_points_to_send))
   !   allocate(el_owner0s(max_n_points_to_send))
   !   n_glb_point_cand = sum(this%n_points_pe)
   !   allocate(rst_results(3,n_glb_point_cand))
   !   allocate(res_results(3,n_glb_point_cand))
   !   allocate(el_owner_results(n_glb_point_cand))
   !   res = 1e2
   !   this%rst = 1e2
   !   this%pe_owner = -1
   !   res_results = 0.0
   !   rst_results = 0.0
   !   call gs_find_back%nbrecv()
   !   call gs_find_back%nbsend(this%xyz_local, this%n_points_local*3, null_ptr, null_ptr)
   !   call gs_find_back%nbwait(res_results, n_glb_point_cand*3, GS_OP_ADD, null_ptr)
   !   call gs_find_back%nbrecv()
   !   call gs_find_back%nbsend(this%rst_local, this%n_points_local*3, null_ptr, null_ptr)
   !   call gs_find_back%nbwait(rst_results, n_glb_point_cand*3, GS_OP_ADD, null_ptr)
   !   do i = 1, size(gs_find_back%send_pe)
   !      rank = gs_find_back%send_pe(i)
   !      call MPI_Isend(this%el_owner0_local(this%n_points_offset_pe_local(rank)+1),&
   !         this%n_points_pe_local(rank), &
   !         MPI_INTEGER, rank, 0, &
   !         this%comm, gs_find_back%send_buf(i)%request, ierr)
   !      gs_find_back%send_buf(i)%flag = .false.
   !   end do
   !   do i = 1, size(gs_find_back%recv_pe)
   !      rank = gs_find_back%recv_pe(i)
   !      call MPI_IRecv(el_owner_results(this%n_points_offset_pe(rank)+1),&
   !         this%n_points_pe(rank), &
   !         MPI_INTEGER, rank, 0, &
   !         this%comm, gs_find_back%recv_buf(i)%request, ierr)
   !      gs_find_back%recv_buf(i)%flag = .false.
   !   end do
   !   call gs_find_back%nbwait_no_op()
   !   ii = 0
   !   do i = 1, size(gs_find_back%recv_pe)
   !      point_ids => points_at_pe(gs_find_back%recv_pe(i))%array()
   !      do j = 1, this%n_points_pe(gs_find_back%recv_pe(i))
   !         point_id = point_ids(j)
   !         ii = ii + 1
   !         if (rst_cmp(this%rst(:,point_id), rst_results(:,ii), &
   !            res(:,point_id), res_results(:,ii), this%tol)) then
   !            this%rst(:,point_ids(j)) = rst_results(:,ii)
   !            res(:,point_ids(j)) = res_results(:,ii)
   !            this%pe_owner(point_ids(j)) = gs_find_back%recv_pe(i)
   !            this%el_owner0(point_ids(j)) = el_owner_results(ii)
   !         end if
   !      end do
   !   end do

   !   !OK, now I know the correct rst values
   !   !of the points I want
   !   !We now send the correct rsts to the correct rank (so a point only belongs to one rank)
   !   do i = 0, this%pe_size-1
   !      call points_at_pe(i)%free()
   !      call points_at_pe(i)%init()
   !      this%n_points_pe(i) = 0
   !   end do

   !   do i = 1, this%n_points
   !      stupid_intent = i
   !      if (this%pe_owner(i) .eq. -1) print *, 'Something is not right for global interpolation',&
   !         ' rank, point coords', stupid_intent, this%xyz(:,i)
   !      call points_at_pe(this%pe_owner(i))%push(stupid_intent)

   !      this%n_points_pe(this%pe_owner(i)) =  this%n_points_pe(this%pe_owner(i)) + 1
   !   end do
   !   call MPI_Reduce_scatter_block(this%n_points_pe, this%n_points_local, 1, MPI_INTEGER, &
   !      MPI_SUM, this%comm, ierr)
   !   call MPI_Alltoall(this%n_points_pe, 1, MPI_INTEGER,&
   !      this%n_points_pe_local, 1, MPI_INTEGER, this%comm, ierr)
   !   this%n_points_offset_pe_local(0) = 0
   !   this%n_points_offset_pe(0) = 0
   !   do i = 1, (this%pe_size - 1)
   !      this%n_points_offset_pe_local(i) = this%n_points_pe_local(i-1)&
   !         + this%n_points_offset_pe_local(i-1)
   !      this%n_points_offset_pe(i) = this%n_points_pe(i-1)&
   !         + this%n_points_offset_pe(i-1)
   !   end do
   !   call send_pe_find%free()
   !   call recv_pe_find%free()
   !   call gs_find%free()
   !   call send_pe_find%init()
   !   call recv_pe_find%init()
   !   call gs_find%init_dofs()
   !   !setup comm to send xyz and rst to chosen ranks
   !   do i = 0, (this%pe_size-1)
   !      if (this%n_points_pe(i) .gt. 0) then
   !         call send_pe_find%push(i)
   !         point_ids => points_at_pe(i)%array()
   !         do j = 1, this%n_points_pe(i)
   !            call gs_find%send_dof(i)%push(3*(point_ids(j)-1)+1)
   !            call gs_find%send_dof(i)%push(3*(point_ids(j)-1)+2)
   !            call gs_find%send_dof(i)%push(3*(point_ids(j)-1)+3)
   !         end do
   !      end if
   !      if (this%n_points_pe_local(i) .gt. 0) then
   !         call recv_pe_find%push(i)
   !         do j = 1, this%n_points_pe_local(i)
   !            call gs_find%recv_dof(i)%push(3*(j+this%n_points_offset_pe_local(i)-1)+1)
   !            call gs_find%recv_dof(i)%push(3*(j+this%n_points_offset_pe_local(i)-1)+2)
   !            call gs_find%recv_dof(i)%push(3*(j+this%n_points_offset_pe_local(i)-1)+3)
   !         end do
   !      end if
   !   end do


   !   call gs_find%init(send_pe_find, recv_pe_find)
   !   call gs_find%nbrecv()
   !   this%xyz_local = 0.0
   !   call gs_find%nbsend(this%xyz, this%n_points*3, null_ptr, null_ptr)
   !   call gs_find%nbwait(this%xyz_local, this%n_points_local*3, GS_OP_ADD, null_ptr)
   !   call gs_find%nbrecv()
   !   this%rst_local = 0.0
   !   call gs_find%nbsend(this%rst, this%n_points*3, null_ptr, null_ptr)
   !   call gs_find%nbwait(this%rst_local, this%n_points_local*3, GS_OP_ADD, null_ptr)
   !   ii = 0
   !   do i = 1, size(gs_find%send_pe)
   !      rank = gs_find%send_pe(i)
   !      point_ids => points_at_pe(rank)%array()
   !      do j = 1, this%n_points_pe(rank)
   !         ii = ii + 1
   !         el_owner_results(ii) = this%el_owner0(point_ids(j))
   !      end do
   !      call MPI_Isend(el_owner_results(this%n_points_offset_pe(rank)+1),&
   !         this%n_points_pe(rank), &
   !         MPI_INTEGER, rank, 0, &
   !         this%comm, gs_find%send_buf(i)%request, ierr)
   !      gs_find%send_buf(i)%flag = .false.
   !   end do
   !   do i = 1, size(gs_find%recv_pe)
   !      rank = gs_find%recv_pe(i)
   !      call MPI_IRecv(this%el_owner0_local(this%n_points_offset_pe_local(rank)+1),&
   !         this%n_points_pe_local(rank), &
   !         MPI_INTEGER, rank, 0, &
   !         this%comm, gs_find%recv_buf(i)%request, ierr)
   !      gs_find%recv_buf(i)%flag = .false.
   !   end do
   !   call gs_find%nbwait_no_op()

   !   call gs_find%free()

   !   !Set up final way of doing communication
   !   call send_pe%init()
   !   call recv_pe%init()
   !   call this%gs_comm%init_dofs()
   !   do i = 0, (this%pe_size-1)
   !      if (this%n_points_pe(i) .gt. 0) then
   !         call recv_pe%push(i)
   !         point_ids => points_at_pe(i)%array()
   !         do j = 1, this%n_points_pe(i)
   !            call this%gs_comm%recv_dof(i)%push(point_ids(j))
   !         end do
   !      end if
   !      if (this%n_points_pe_local(i) .gt. 0) then
   !         call send_pe%push(i)
   !         do j = 1, this%n_points_pe_local(i)
   !            call this%gs_comm%send_dof(i)%push(j+this%n_points_offset_pe_local(i))
   !         end do
   !      end if
   !   end do
   !   call this%gs_comm%init(send_pe, recv_pe,this%comm)


   !   if (allocated(this%pt_ids)) deallocate(this%pt_ids)
   !   allocate(this%pt_ids(this%n_points))
   !   ii = 0
   !   do i = 0, (this%pe_size - 1)
   !      point_ids => points_at_pe(i)%array()
   !      do j = 1, this%n_points_pe(i)
   !         ii = ii + 1
   !         this%pt_ids(ii) = point_ids(j)
   !      end do
   !   end do
   !   !Initialize working arrays for evaluation
   !   call this%temp_local%init(this%n_points_local)
   !   call this%temp%init(this%n_points)

   !   !Initialize arrays to double check interpolation
   !   call x_check%init(this%n_points)
   !   call y_check%init(this%n_points)
   !   call z_check%init(this%n_points)
   !   call this%local_interp%free()
   !   !Initialize interpolator for local interpolation
   !   call this%local_interp%init(this%Xh, this%rst_local(1,:),&
   !      this%rst_local(2,:), &
   !      this%rst_local(3,:), this%n_points_local)
   !   call this%evaluate(x_check%x, this%x%ptr, .true.)
   !   call this%evaluate(y_check%x, this%y%ptr, .true.)
   !   call this%evaluate(z_check%x, this%z%ptr, .true.)

   !   j = 0
   !   do i = 1 , this%n_points

   !      ! Check validity of points
   !      isdiff = .false.
   !      xdiff = x_check%x(i)-this%xyz(1,i)
   !      ydiff = y_check%x(i)-this%xyz(2,i)
   !      zdiff = z_check%x(i)-this%xyz(3,i)
   !      isdiff = norm2(real((/xdiff,ydiff,zdiff/),xp)) > this%tol
   !      isdiff = isdiff .or. abs(this%rst(1,i)) > 1.0_xp + this%tol
   !      isdiff = isdiff .or. abs(this%rst(2,i)) > 1.0_xp + this%tol
   !      isdiff = isdiff .or. abs(this%rst(3,i)) > 1.0_xp + this%tol
   !      if (isdiff ) then
   !         write(*,*) 'Point with coordinates: ', &
   !            this%xyz(1, i), this%xyz(2, i), this%xyz(3, i), &
   !            'Differ from interpolated coords: ', &
   !            x_check%x(i), y_check%x(i), z_check%x(i), &
   !            'Actual difference: ', &
   !            xdiff, ydiff, zdiff, norm2(real((/xdiff,ydiff,zdiff/),xp)),&
   !            'Expected difference: ', &
   !            res(:,i), norm2(real(res(:,i),xp)),&
   !            'Process, element: ', &
   !            this%pe_owner(i), this%el_owner0(i)+1, &
   !            'rst coords: ', &
   !            this%rst(:,i), &
   !            ' radius', sqrt(this%xyz(1,i)**2.0_xp+this%xyz(2,i)**2.0_xp)
   !         j = j + 1
   !      end if
   !   end do

   !   if (NEKO_BCKND_DEVICE .eq. 1) then
   !      call device_memcpy(this%el_owner0, this%el_owner0_d, &
   !         this%n_points, HOST_TO_DEVICE, sync = .true.)
   !      call device_map(this%el_owner0_local, this%el_owner0_local_d, this%n_points_local)
   !      call device_memcpy(this%el_owner0_local, this%el_owner0_local_d, &
   !         this%n_points_local, HOST_TO_DEVICE, sync = .true.)
   !   end if

   !   !Free stuff
   !   call send_pe%free()
   !   call recv_pe%free()
   !   call gs_find_back%free()
   !   call send_pe_find%free()
   !   call recv_pe_find%free()
   !   call x_check%free()
   !   call x_vec%free()
   !   call y_check%free()
   !   call y_vec%free()
   !   call z_check%free()
   !   call z_vec%free()
   !   call x_t%free()
   !   call y_t%free()
   !   call z_t%free()
   !   call rst_local_cand%free()
   !   call resx%free()
   !   call resy%free()
   !   call resz%free()
   !   call x_hat%free()
   !   call y_hat%free()
   !   call z_hat%free()
   !   if (allocated(conv_pts)) deallocate(conv_pts)
   !   if (allocated(res)) deallocate(res)
   !   call all_el_candidates%free()
   !   if (allocated(points_at_pe)) then
   !      do i = 0, this%pe_size-1
   !         call points_at_pe(i)%free()
   !      end do
   !      deallocate(points_at_pe)
   !   end if
   !   call pe_candidates%free()
   !   call el_candidates%free()
   !   if (associated(pe_cands)) pe_cands => Null()
   !   if (associated(el_cands)) pe_cands => Null()
   !   if (associated(point_ids)) point_ids => Null()
   !   if (allocated(xyz_send_to_pe)) deallocate(xyz_send_to_pe)
   !   if (allocated(rst_send_to_pe)) deallocate(rst_send_to_pe)
   !   if (allocated(rst_recv_from_pe)) deallocate(rst_recv_from_pe)
   !   if (allocated(res_recv_from_pe)) deallocate(res_recv_from_pe)
   !   if (allocated(el_owner0s)) deallocate(el_owner0s)
   !   if (allocated(el_send_to_pe)) deallocate(el_send_to_pe)

   !   call MPI_Barrier(this%comm)
   !   time2 = MPI_Wtime()
   !   write(log_buf, '(A,E15.7)') 'Global interpolation find points done, time (s):', &
   !      time2-time_start
   !   call neko_log%message(log_buf)

   !end subroutine aabb_pe_finder_find


end module aabb_pe_finder

