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
  use num_types, only : rp, dp, xp
  use neko_config, only : NEKO_BCKND_DEVICE
  use space, only : space_t
  use pe_finder, only : pe_finder_t
  use stack, only : stack_i4_t, stack_i4t2_t
  use utils, only : neko_error, neko_warning
  use tuple, only : tuple_i4_t
  use htable, only : htable_i4_t
  use point, only : point_t
  use comm, only : NEKO_COMM, MPI_REAL_PRECISION, pe_rank, pe_size
  use mpi_f08, only : MPI_SUM, MPI_Reduce, MPI_COMM, MPI_Comm_rank, &
       MPI_Comm_size, MPI_Wtime, MPI_INTEGER, MPI_IN_PLACE, &
       MPI_MIN, MPI_Allgather, MPI_Barrier, MPI_DOUBLE_PRECISION, &
       MPI_Allreduce
  use aabb, only : aabb_t
  use aabb_tree, only : aabb_tree_t, aabb_node_t, AABB_NULL_NODE
  use math, only : NEKO_M_LN2, NEKO_EPS
  use, intrinsic :: iso_c_binding, only : c_ptr, c_null_ptr, c_associated, &
       c_sizeof, c_bool, c_loc
  implicit none
  private

  !> Minimum number of total boxes in the aabb tree
  integer, public, parameter :: GLOB_MAP_SIZE = 4096

  !> Implements global interpolation for arbitrary points in the domain.
  type, public, extends(pe_finder_t) :: aabb_pe_finder_t
     !> Which communicator to find things on
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
    real(kind=dp), allocatable :: rank_xyz_max(:,:), rank_xyz_min(:,:)
    real(kind=dp), allocatable :: max_xyz(:,:), min_xyz(:,:)
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
    ! Create a local tree for each element at this rank
    call local_aabb_tree%init(nelv)

    allocate(local_aabb(nelv))

    do i = 1, nelv
       id1 = lx*ly*lz*(i-1)+1
       id2 = lx*ly*lz*(i)
       call local_aabb(i)%init( real((/minval(x(id1:id2)), &
            minval(y(id1:id2)), &
            minval(z(id1:id2))/), dp), &
            real((/maxval(x(id1:id2)), &
            maxval(y(id1:id2)), &
            maxval(z(id1:id2))/), dp))
       call local_aabb_tree%insert_object(local_aabb(i), i)
    end do

    this%pe_box_num = min(GLOB_MAP_SIZE/this%pe_size, nelv)
    this%pe_box_num = &
         max(1, ishft(1, ceiling(log(real(this%pe_box_num, rp)) / NEKO_M_LN2)))

    call MPI_Allreduce(MPI_IN_PLACE, this%pe_box_num, 1, MPI_INTEGER, &
         MPI_MIN, this%comm, ierr)

    this%pe_box_num = max(this%pe_box_num,2) !> At least 2 boxes
    this%glob_map_size = this%pe_box_num*this%pe_size

    if (pe_rank .eq. 0) then
       print *, this%pe_box_num, this%glob_map_size
    end if

    allocate(rank_xyz_max(3,this%glob_map_size))
    allocate(rank_xyz_min(3,this%glob_map_size))
    allocate(min_xyz(3,this%pe_box_num))
    allocate(max_xyz(3,this%pe_box_num))

    i = 1
    id_lvl = (/local_aabb_tree%get_root_index(), 0/)
    call traverse_stack%init()
    call traverse_stack%push(id_lvl)

    lvl = 0
    !> Traverse the local tree and find ther top boxes
    do while (traverse_stack%size() .gt. 0)
       id_lvl = traverse_stack%pop()
       lvl = id_lvl%x(2)
       node = local_aabb_tree%get_node(id_lvl%x(1))
       if (2**lvl .eq. this%pe_box_num .or. node%is_leaf()) then
          min_xyz(:,i) = node%aabb%get_min()
          max_xyz(:,i) = node%aabb%get_max()
          i = i + 1
       else if (2**lvl < this%pe_box_num) then
          if (node%get_left_index() .ne. AABB_NULL_NODE) then
             temp_id_lvl = (/node%get_left_index(), lvl+1/)
             call traverse_stack%push(temp_id_lvl)
          end if
          if (node%get_right_index() .ne. AABB_NULL_NODE) then
             temp_id_lvl = (/node%get_right_index(), lvl+1/)
             call traverse_stack%push(temp_id_lvl)
          end if
       end if
    end do
    call traverse_stack%free()

    !If somehow we dont need all boxes we just put a point here
    if (nelv .eq. 0) then
       !> Set the boxes to be empty
       do j = 1, this%pe_box_num
          min_xyz(:,j) = [NEKO_EPS, NEKO_EPS, NEKO_EPS]
          max_xyz(:,j) = [NEKO_EPS, NEKO_EPS, NEKO_EPS]
       end do
    else
       !Needs to be something in the domain
       do j = i, this%pe_box_num
          min_xyz(:,j) = [x(1), y(1), z(1)]
          max_xyz(:,j) = [x(1), y(1), z(1)]
       end do
    end if
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
    call this%global_aabb_tree%build_from_aabb(this%global_aabb,padding)
    deallocate(local_aabb)
    deallocate(rank_xyz_max)
    deallocate(rank_xyz_min)
    deallocate(min_xyz)
    deallocate(max_xyz)
  end subroutine aabb_pe_finder_init


  !> Destructor
  subroutine aabb_pe_finder_free(this)
    class(aabb_pe_finder_t), intent(inout) :: this

    if (allocated(this%global_aabb)) then
       deallocate(this%global_aabb)
    end if

  end subroutine aabb_pe_finder_free

  !> Find pe candidates
  !! @param my_point Point to find candidates for.
  !! @param pe_candidates Candidates for the point.
  subroutine aabb_pe_finder_find_candidates(this, my_point, pe_candidates)
    class(aabb_pe_finder_t), intent(inout) :: this
    type(point_t), intent(in) :: my_point
    type(stack_i4_t), intent(inout) :: pe_candidates
    integer, pointer :: pe_cands(:)
    integer :: i

    call this%global_aabb_tree%query_overlaps(my_point, -1, pe_candidates)
    pe_cands => pe_candidates%array()
    do i = 1, pe_candidates%size()
       pe_cands(i) = (pe_cands(i)-1) / this%pe_box_num
    end do

  end subroutine aabb_pe_finder_find_candidates

  subroutine aabb_pe_finder_find_candidates_batch(this, points, &
       n_points, points_at_pe, n_points_pe)
    class(aabb_pe_finder_t), intent(inout) :: this
    integer, intent(in) :: n_points
    real(kind=rp), intent(in) :: points(3,n_points)
    type(stack_i4_t), intent(inout) :: points_at_pe(0:(this%pe_size-1))
    integer, intent(inout) :: n_points_pe(0:(this%pe_size-1))
    type(stack_i4_t) :: pe_candidates
    type(point_t) :: my_point
    integer :: i, j, temp_intent, pe_id, htable_data
    real(kind=dp) :: pt_xyz(3)
    integer, pointer :: pe_cands(:)
    type(htable_i4_t) :: marked_rank

    do i = 0, this%pe_size-1
       call points_at_pe(i)%clear()
       n_points_pe(i) = 0
    end do

    call marked_rank%init(32, htable_data)
    call pe_candidates%init()

    !> Check which ranks might have this point
    do i = 1, n_points
       call marked_rank%clear()
       pt_xyz = (/ points(1,i), points(2,i), points(3,i) /)
       call my_point%init(pt_xyz)
       call pe_candidates%clear()
       call this%find(my_point, pe_candidates)

       pe_cands => pe_candidates%array()
       do j = 1, pe_candidates%size()
          pe_id = pe_cands(j)
          ! Check if this rank has already been marked
          if (marked_rank%get(pe_id, htable_data) .ne. 0) then
             n_points_pe(pe_id) = n_points_pe(pe_id) + 1
             temp_intent = i
             call points_at_pe(pe_id)%push(temp_intent)
             ! Doesnt matter, I use htable as a set
             htable_data = 100
             call marked_rank%set(pe_id, htable_data)
          end if
       end do

       if (pe_candidates%size() .lt. 1) then
          write (*,*) 'Point', points(:,i), &
               'found to be outside domain, try increasing the padding to find rank candidates.'
       end if
    end do
    call marked_rank%free()
    call pe_candidates%free()

  end subroutine aabb_pe_finder_find_candidates_batch

end module aabb_pe_finder
