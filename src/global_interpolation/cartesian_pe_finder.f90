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

!> Implements cartesian_pe_finder given a dofmap.
!!
module cartesian_pe_finder
  use num_types, only: rp, dp, xp, i8
  use neko_config, only : NEKO_BCKND_DEVICE
  use space, only: space_t
  use pe_finder, only: pe_finder_t
  use stack, only: stack_i4_t, stack_i8_t
  use utils, only: neko_error, neko_warning, linear_index
  use tuple, only: tuple_i4_t
  use htable, only: htable_i8_t
  use point, only: point_t
  use comm, only: NEKO_COMM, MPI_REAL_PRECISION, pe_rank, pe_size
  use mpi_f08, only: MPI_MAX, MPI_Allreduce, MPI_COMM, MPI_Comm_rank, &
       MPI_Comm_size, MPI_Wtime, MPI_INTEGER, MPI_INTEGER8, &
       MPI_MIN, MPI_SUM, MPI_Irecv, MPI_Isend, MPI_Wait, &
       MPI_Exscan, MPI_Request, MPI_Status, &
       MPI_Alltoall, MPI_IN_PLACE, MPI_Barrier
  implicit none
  private

  !> Minimum number of total boxes in the aabb tree
  integer, public, parameter :: GLOB_MAP_SIZE = 4096

  type, private :: i8_mpi_t
     !> status of this operation
     type(MPI_Status) :: status
     !> unique request for this buffer
     type(MPI_Request) :: request
     !> flag = false when operation is in process
     !! Set to true when operation has succeeded
     logical :: flag
     !> Buffer with data to send/recieve
     real(kind=dp), allocatable :: data(:)
     integer :: size = 0
   contains
     procedure, pass(this) :: free => i8_mpi_free
  end type i8_mpi_t


  !> Implements global interpolation for arbitrary points in the domain.
  type, public, extends(pe_finder_t) :: cartesian_pe_finder_t
     !> Which communicator to find things on
     real(kind=dp) :: padding
     !> global number of boxes
     integer(kind=i8) :: glob_n_boxes
     integer(kind=i8) :: n_boxes
     !Number of local boxes in x direction
     integer :: n_boxes_per_pe
     integer :: nelv
     integer :: offset_el
     real(kind=rp) :: max_x_global, max_y_global, max_z_global
     real(kind=rp) :: min_x_global, min_y_global, min_z_global

     ! Resolution of the boxes
     real(kind=xp) :: res_x, res_y, res_z

     type(i8_mpi_t), allocatable :: recv_buf(:)
     type(i8_mpi_t), allocatable :: send_buf(:)
     type(stack_i4_t), allocatable :: pe_map(:)

   contains
     procedure, pass(this) :: init => cartesian_pe_finder_init
     procedure, pass(this) :: free => cartesian_pe_finder_free
     procedure, pass(this) :: find => cartesian_pe_finder_find
     procedure, pass(this) :: find_batch => cartesian_pe_finder_find_batch

  end type cartesian_pe_finder_t

contains

  !> Destructor for i8_mpi_t.
  subroutine i8_mpi_free(this)
    class(i8_mpi_t), intent(inout) :: this

    if (allocated(this%data)) then
       deallocate(this%data)
    end if

  end subroutine i8_mpi_free

  !> Initialize the global interpolation object on a set of coordinates.
  !! @param x x-coordinates.
  !! @param y y-coordinates.
  !! @param z z-coordinates.
  !! @param gdim Geometric dimension.
  !! @param nelv Number of elements of the mesh in which to search for the
  !! points.
  !! @param Xh Space on which to interpolate.
  !! @param tol Tolerance for Newton iterations.
  subroutine cartesian_pe_finder_init(this, x, y, z, nelv, Xh, comm, n_boxes, padding)
    class(cartesian_pe_finder_t), intent(inout) :: this
    real(kind=rp), intent(in), target :: x(:)
    real(kind=rp), intent(in), target :: y(:)
    real(kind=rp), intent(in), target :: z(:)
    integer, intent(in) :: nelv
    type(MPI_COMM), intent(in), optional :: comm
    type(space_t), intent(in), target :: Xh
    integer, intent(in) :: n_boxes
    real(kind=dp), intent(in) :: padding
    integer :: ierr
    integer :: i, j, k, e, i2, j2, k2
    integer :: pe_id, lin_idx
    integer :: lxyz, lx2
    integer(kind=i8) :: glob_id, loc_id
    integer(kind=i8) :: min_id(3), max_id(3)
    real(kind=rp) :: center_x, center_y, center_z
    type(stack_i8_t), allocatable :: glob_ids(:), recv_ids(:)
    integer(i8), pointer :: glb_ids(:)
    integer(kind=i8) :: htable_data, temp ! We just use it as a set
    integer, allocatable :: n_recv(:), n_send(:)
    type(htable_i8_t) :: marked_box
    real(kind=rp) :: min_bb_x, max_bb_x
    real(kind=rp) :: min_bb_y, max_bb_y
    real(kind=rp) :: min_bb_z, max_bb_z
    real(kind=rp) :: el_x(Xh%lxyz), el_y(Xh%lxyz), el_z(Xh%lxyz)

    call this%free()
    this%comm = comm
    this%padding = padding

    call MPI_Comm_rank(this%comm, this%pe_rank, ierr)
    call MPI_Comm_size(this%comm, this%pe_size, ierr)

    this%glob_n_boxes = int(n_boxes,i8)**3
    this%n_boxes = n_boxes
    this%nelv = nelv
    this%n_boxes_per_pe = (this%glob_n_boxes+int(this%pe_size-1,i8))/int(this%pe_size,i8)

    allocate(this%send_buf(0:this%pe_size-1))
    allocate(this%recv_buf(0:this%pe_size-1))
    allocate(this%pe_map(0:this%n_boxes_per_pe-1))
    do i = 0, this%n_boxes_per_pe-1
       call this%pe_map(i)%init()
    end do
    call MPI_Exscan(this%nelv, this%offset_el, 1, &
         MPI_INTEGER, MPI_SUM, this%comm, ierr)
    if (this%nelv .gt. 0) then
       this%max_x_global = maxval(x(1:nelv*Xh%lxyz))
       this%max_y_global = maxval(y(1:nelv*Xh%lxyz))
       this%max_z_global = maxval(z(1:nelv*Xh%lxyz))
       this%min_x_global = minval(x(1:nelv*Xh%lxyz))
       this%min_y_global = minval(y(1:nelv*Xh%lxyz))
       this%min_z_global = minval(z(1:nelv*Xh%lxyz))
    else
       this%max_x_global = -1e20
       this%max_y_global = -1e20
       this%max_z_global = -1e20
       this%min_x_global = 1e20
       this%min_y_global = 1e20
       this%min_z_global = 1e20
    end if


    call MPI_Allreduce(MPI_IN_PLACE, this%max_x_global, 1, MPI_REAL_PRECISION, &
         MPI_MAX, this%comm, ierr)
    call MPI_Allreduce(MPI_IN_PLACE, this%max_y_global, 1, MPI_REAL_PRECISION, &
         MPI_MAX, this%comm, ierr)
    call MPI_Allreduce(MPI_IN_PLACE, this%max_z_global, 1, MPI_REAL_PRECISION, &
         MPI_MAX, this%comm, ierr)
    call MPI_Allreduce(MPI_IN_PLACE, this%min_x_global, 1, MPI_REAL_PRECISION, &
         MPI_MIN, this%comm, ierr)
    call MPI_Allreduce(MPI_IN_PLACE, this%min_y_global, 1, MPI_REAL_PRECISION, &
         MPI_MIN, this%comm, ierr)
    call MPI_Allreduce(MPI_IN_PLACE, this%min_z_global, 1, MPI_REAL_PRECISION, &
         MPI_MIN, this%comm, ierr)

    center_x = (this%max_x_global + this%min_x_global) / 2.0_dp
    center_y = (this%max_y_global + this%min_y_global) / 2.0_dp
    center_z = (this%max_z_global + this%min_z_global) / 2.0_dp

    this%max_x_global = this%max_x_global - center_x
    this%max_y_global = this%max_y_global - center_y
    this%max_z_global = this%max_z_global - center_z
    this%min_x_global = this%min_x_global - center_x
    this%min_y_global = this%min_y_global - center_y
    this%min_z_global = this%min_z_global - center_z
    this%max_x_global = this%max_x_global*(1.0_xp+this%padding) + center_x
    this%max_y_global = this%max_y_global*(1.0_xp+this%padding) + center_y
    this%max_z_global = this%max_z_global*(1.0_xp+this%padding) + center_z
    this%min_x_global = this%min_x_global*(1.0_xp+this%padding) + center_x
    this%min_y_global = this%min_y_global*(1.0_xp+this%padding) + center_y
    this%min_z_global = this%min_z_global*(1.0_xp+this%padding) + center_z


    this%res_x = (this%max_x_global - this%min_x_global) / real(this%n_boxes, xp)
    this%res_y = (this%max_y_global - this%min_y_global) / real(this%n_boxes, xp)
    this%res_z = (this%max_z_global - this%min_z_global) / real(this%n_boxes, xp)
    if (allocated(recv_ids)) then
       do i = 0, this%pe_size-1
          call recv_ids(i)%free()
       end do
       deallocate(recv_ids)
    end if
    if (allocated(glob_ids)) then
       do i = 0, this%pe_size-1
          call glob_ids(i)%free()
       end do
       deallocate(glob_ids)
    end if
    if (allocated(n_recv)) deallocate(n_recv)
    if (allocated(n_send)) deallocate(n_send)
    allocate(n_recv(0:this%pe_size-1))
    allocate(n_send(0:this%pe_size-1))
    allocate(recv_ids(0:this%pe_size-1))
    allocate(glob_ids(0:this%pe_size-1))
    do i = 0, this%pe_size-1
       call recv_ids(i)%init()
       call glob_ids(i)%init()
       if (allocated(this%recv_buf(i)%data)) then
          deallocate(this%recv_buf(i)%data)
       end if
       if (allocated(this%send_buf(i)%data)) then
          deallocate(this%send_buf(i)%data)
       end if
       this%send_buf(i)%size = 0
       this%recv_buf(i)%size = 0
    end do
    n_recv = 0
    n_send = 0

    lxyz = Xh%lxyz
    call marked_box%init(this%nelv*lxyz,htable_data)

    do e = 1, this%nelv
       !move it to do scaling
       lx2 = Xh%lx/2
       if (mod(Xh%lx,2) .eq. 0) then
          lin_idx = linear_index(lx2, lx2, lx2, e, Xh%lx, Xh%lx, Xh%lx)
          center_x = x(lin_idx)
          center_y = y(lin_idx)
          center_z = z(lin_idx)
       else
          center_x = 0d0
          center_y = 0d0
          center_z = 0d0
          do i = lx2, lx2+1
             do j = lx2, lx2 + 1
                do k = lx2, lx2 + 1
                   lin_idx = linear_index(i, j, k, e, Xh%lx, Xh%lx, Xh%lx)
                   center_x = center_x + x(lin_idx)
                   center_y = center_y + y(lin_idx)
                   center_z = center_z + z(lin_idx)
                end do
             end do
          end do
          center_x = center_x / 8.0_xp
          center_y = center_y / 8.0_xp
          center_z = center_z / 8.0_xp
       end if

       !Maybe this is a stupid way to do it

       el_x = x((e-1)*lxyz+1:e*lxyz) - center_x
       el_y = y((e-1)*lxyz+1:e*lxyz) - center_y
       el_z = z((e-1)*lxyz+1:e*lxyz) - center_z
       el_x = el_x * (1.0_rp+padding) + center_x
       el_y = el_y * (1.0_rp+padding) + center_y
       el_z = el_z * (1.0_rp+padding) + center_z
       !Padded and ready

       !Now we go the bounding boxes of all subboxes in the element
       do i = 1, Xh%lx - 1
          do j = 1, Xh%ly - 1
             do k = 1, Xh%lz - 1
                lin_idx = linear_index(i, j, k, 1, Xh%lx, Xh%lx, Xh%lx)
                max_bb_x = el_x(lin_idx)
                min_bb_x = el_x(lin_idx)
                max_bb_y = el_y(lin_idx)
                min_bb_y = el_y(lin_idx)
                max_bb_z = el_z(lin_idx)
                min_bb_z = el_z(lin_idx)
                do i2 = 0, 1
                   do j2 = 0, 1
                      do k2 = 0, 1
                         lin_idx = linear_index(i+i2,j+j2,k+k2, 1, Xh%lx, Xh%lx, Xh%lx)
                         max_bb_x = max(max_bb_x, el_x(lin_idx))
                         min_bb_x = min(min_bb_x, el_x(lin_idx))
                         max_bb_y = max(max_bb_y, el_y(lin_idx))
                         min_bb_y = min(min_bb_y, el_y(lin_idx))
                         max_bb_z = max(max_bb_z, el_z(lin_idx))
                         min_bb_z = min(min_bb_z, el_z(lin_idx))
                      end do
                   end do
                end do


                min_id = get_global_idx(this, min_bb_x, min_bb_y, min_bb_z)
                max_id = get_global_idx(this, max_bb_x, max_bb_y, max_bb_z)
                do i2 = min_id(1), max_id(1)
                   do j2 = min_id(2), max_id(2)
                      do k2 = min_id(3), max_id(3)
                         if (i2 .ge. 0 .and. i2 .lt. this%n_boxes .and. &
                              j2 .ge. 0 .and. j2 .lt. this%n_boxes .and. &
                              k2 .ge. 0 .and. k2 .lt. this%n_boxes) then
                            glob_id = i2 + j2*this%n_boxes + k2*this%n_boxes**2
                            pe_id = get_pe_idx(this, glob_id)
                            if (marked_box%get(glob_id,htable_data) .ne. 0) then
                               call glob_ids(pe_id)%push(glob_id)
                               n_send(pe_id) = n_send(pe_id) + 1
                               call marked_box%set(glob_id, htable_data)
                            end if
                         end if
                      end do
                   end do
                end do
             end do
          end do
       end do
    end do
    !temp = 0
    !do i = 0, this%pe_size-1
    !   temp = temp + glob_ids(i)%size()
    !end do
    !print *, 'size of globids', temp, this%pe_rank
    call send_recv_data(this, recv_ids, n_recv, glob_ids, n_send)

    !temp = 0
    !do i = 0, this%pe_size-1
    !   temp = temp + recv_ids(i)%size()
    !end do
    !print *, 'size of globids', temp, this%pe_rank
    !Buffers are likely way too large for other calls
    do i = 0, this%pe_size-1
       if (allocated(this%recv_buf(i)%data)) then
          deallocate(this%recv_buf(i)%data)
       end if
       if (allocated(this%send_buf(i)%data)) then
          deallocate(this%send_buf(i)%data)
       end if
       this%send_buf(i)%size = 0
       this%recv_buf(i)%size = 0
    end do


    call MPI_Barrier(this%comm, ierr)
    do i = 0, this%pe_size-1
       if (n_recv(i) .gt. 0) then
          glb_ids => recv_ids(i)%array()
          do j = 1, n_recv(i)
             glob_id = glb_ids(j)
             loc_id = glob_id - int(this%pe_rank,i8) * &
                  int(this%n_boxes_per_pe,i8)
             if (loc_id .ge. this%n_boxes_per_pe .or. loc_id < 0) then
                call neko_warning('loc_id out of bounds')
             end if
             call this%pe_map(int(loc_id))%push(i)
          end do
       end if
    end do

    if (allocated(recv_ids)) then
       do i = 0, this%pe_size-1
          call recv_ids(i)%free()
       end do
       deallocate(recv_ids)
    end if
    if (allocated(glob_ids)) then
       do i = 0, this%pe_size-1
          call glob_ids(i)%free()
       end do
       deallocate(glob_ids)
    end if
    if (allocated(n_recv)) deallocate(n_recv)
    if (allocated(n_send)) deallocate(n_send)
  end subroutine cartesian_pe_finder_init

  subroutine cartesian_pe_finder_find(this, my_point, pe_candidates)
    class(cartesian_pe_finder_t), intent(inout) :: this
    type(point_t), intent(in) :: my_point
    type(stack_i4_t), intent(inout) :: pe_candidates

    call neko_error('cartesian_pe_finder_find not implemented')
  end subroutine cartesian_pe_finder_find

  function get_global_idx(this, x, y, z) result(global_box_id)
    class(cartesian_pe_finder_t), intent(in) :: this
    real(kind=rp), intent(in) :: x, y, z
    integer(kind=i8) :: global_box_id(3)
    integer(kind=i8) :: pe_size, n_boxes

    pe_size = this%pe_size
    n_boxes = this%n_boxes


    global_box_id(1) = int((x - this%min_x_global) / this%res_x, i8)
    global_box_id(2) = int((y - this%min_y_global) / this%res_y, i8)
    global_box_id(3) = int((z - this%min_z_global) / this%res_z, i8)
  end function get_global_idx

  function get_pe_idx(this, global_idx) result(pe_id)
    class(cartesian_pe_finder_t), intent(in) :: this
    integer(kind=i8), intent(in) :: global_idx
    integer :: pe_id
    !Get x id and then divide by the number of x boxes per rank to get the correct pe id
    pe_id = global_idx / int(this%n_boxes_per_pe, i8)
  end function get_pe_idx


  !> Destructor
  subroutine cartesian_pe_finder_free(this)
    class(cartesian_pe_finder_t), intent(inout) :: this
    integer :: i

    if (allocated(this%send_buf)) then
       do i = 0, this%pe_size-1
          call this%send_buf(i)%free()
       end do
       deallocate(this%send_buf)
    end if
    if (allocated(this%recv_buf)) then
       do i = 0, this%pe_size-1
          call this%recv_buf(i)%free()
       end do
       deallocate(this%recv_buf)
    end if
    if (allocated(this%pe_map)) then
       do i = 0, this%n_boxes_per_pe-1
          call this%pe_map(i)%free()
       end do
       deallocate(this%pe_map)
    end if

  end subroutine cartesian_pe_finder_free

  subroutine cartesian_pe_finder_find_batch(this, points, n_points, points_at_pe, n_points_pe)
    class(cartesian_pe_finder_t), intent(inout) :: this
    integer, intent(in) :: n_points
    real(kind=rp), intent(in) :: points(3,n_points)
    type(stack_i4_t), intent(inout) :: points_at_pe(0:(this%pe_size-1))
    integer, intent(inout) :: n_points_pe(0:(this%pe_size-1))
    integer :: i, j, k
    integer(kind=i8) :: glob_id(3)
    integer(kind=i8) :: pe_id
    integer(kind=i8) :: lin_glob_id
    integer(kind=i8) :: loc_id
    type(stack_i8_t), allocatable :: work_pe_ids(:), work_pt_ids(:)
    type(stack_i8_t), allocatable :: temp_pe_ids(:), temp_pt_ids(:)
    integer, allocatable :: n_work_ids(:), n_temp_ids(:)
    integer(i8), pointer :: ids(:)
    integer(i8), pointer :: pt_id(:)
    integer, pointer :: pe_cands(:)
    integer(i8), pointer :: pe_cands8(:)
    integer(i8), pointer :: pt_ids(:)
    integer :: ierr

    allocate(work_pe_ids(0:this%pe_size-1))
    allocate(work_pt_ids(0:this%pe_size-1))
    allocate(temp_pe_ids(0:this%pe_size-1))
    allocate(temp_pt_ids(0:this%pe_size-1))
    allocate(n_temp_ids(0:this%pe_size-1))
    allocate(n_work_ids(0:this%pe_size-1))

    do i = 0, this%pe_size-1
       n_work_ids(i) = 0
       n_temp_ids(i) = 0
       call work_pt_ids(i)%init()
       call work_pe_ids(i)%init()
       call temp_pe_ids(i)%init()
       call temp_pt_ids(i)%init()
    end do

    do i = 0, this%pe_size-1
       call points_at_pe(i)%clear()
       n_points_pe(i) = 0
    end do

    ! Compute global ids for the points
    ! and the pe id for each point
    n_work_ids = 0
    do i = 1, n_points
       glob_id = get_global_idx(this, points(1,i), points(2,i), points(3,i))
       lin_glob_id = glob_id(1) + &
            glob_id(2)*this%n_boxes + &
            glob_id(3)*this%n_boxes**2
       pe_id = get_pe_idx(this, lin_glob_id)
       if (glob_id(1) .ge. 0 .and. glob_id(1) .lt. this%n_boxes .and. &
            glob_id(2) .ge. 0 .and. glob_id(2) .lt. this%n_boxes .and. &
            glob_id(3) .ge. 0 .and. glob_id(3) .lt. this%n_boxes) then
          call work_pe_ids(pe_id)%push(lin_glob_id)
          call work_pt_ids(pe_id)%push(int(i,i8))
          n_work_ids(pe_id) = n_work_ids(pe_id) + 1
       else
          print *, 'Point found outside domain:', points(1,i), &
               points(2,i), points(3,i), &
               'Computed global id:', lin_glob_id, &
               'Global box id (x,y,z):', glob_id(1), glob_id(2), glob_id(3)
       end if
    end do


    ! Send the global ids to the correct pe
    ! and get the global ids from the other pes with points I own
    ! Also send point ids to the other pes
    call send_recv_data(this, temp_pe_ids, n_temp_ids, work_pe_ids, n_work_ids)
    call send_recv_data(this, temp_pt_ids, n_temp_ids, work_pt_ids, n_work_ids)
    call MPI_Barrier(this%comm, ierr)
    ! Get the local ids for the points I own
    ! and the pe candidates for the points I own
    n_work_ids = 0
    do i = 0, this%pe_size-1
       call work_pe_ids(i)%clear()
       call work_pt_ids(i)%clear()
    end do

    do i =0 , this%pe_size-1

       if (n_temp_ids(i) .gt. 0) then
          ids => temp_pe_ids(i)%array()
          pt_ids => temp_pt_ids(i)%array()
          do j = 1, n_temp_ids(i)
             loc_id = ids(j) - int(this%pe_rank,i8) * &
                  int(this%n_boxes_per_pe,i8)
             pe_cands => this%pe_map(int(loc_id))%array()
             do k = 1, this%pe_map(int(loc_id))%size()
                call work_pe_ids(i)%push(int(pe_cands(k),i8))
                call work_pt_ids(i)%push(pt_ids(j))
                n_work_ids(i) = n_work_ids(i) + 1
             end do
             if (this%pe_map(int(loc_id))%size() .lt. 1) then
                print *, 'No PE candidates found for point:', &
                     points(1,pt_ids(j)), &
                     points(2,pt_ids(j)), points(3,pt_ids(j))
             end if
          end do
       end if
    end do
    ! pe candidates found for the points I am responsible for
    ! Now i need to send the data to the other pes
    ! and get the candidates from the other pes for my own points
    call send_recv_data(this, temp_pe_ids, n_temp_ids, work_pe_ids, n_work_ids)
    call send_recv_data(this, temp_pt_ids, n_temp_ids, work_pt_ids, n_work_ids)
    ! Gotten candidates for my points
    ! Now I organize the candidates for the points I own
    n_points_pe = 0
    do i = 0, this%pe_size-1
       if (n_temp_ids(i) .gt. 0) then
          pe_cands8 => temp_pe_ids(i)%array()
          pt_ids => temp_pt_ids(i)%array()
          do j = 1, n_temp_ids(i)
             pe_id = pe_cands8(j)
             call points_at_pe(pe_id)%push(int(pt_ids(j)))
             n_points_pe(pe_id) = n_points_pe(pe_id) + 1
          end do
       end if
    end do

    if (allocated(work_pe_ids)) then
       do i = 0, this%pe_size-1
          call work_pe_ids(i)%free()
       end do
       deallocate(work_pe_ids)
    end if
    if (allocated(temp_pe_ids)) then
       do i = 0, this%pe_size-1
          call temp_pe_ids(i)%free()
       end do
       deallocate(temp_pe_ids)
    end if
    if (allocated(temp_pt_ids)) then
       do i = 0, this%pe_size-1
          call temp_pt_ids(i)%free()
       end do
       deallocate(temp_pt_ids)
    end if
    if (allocated(work_pt_ids)) then
       do i = 0, this%pe_size-1
          call work_pt_ids(i)%free()
       end do
       deallocate(work_pt_ids)
    end if

    if (allocated(n_work_ids)) then
       deallocate(n_work_ids)
    end if

    if (allocated(n_temp_ids)) then
       deallocate(n_temp_ids)
    end if
  end subroutine cartesian_pe_finder_find_batch


  subroutine send_recv_data(this, recv_values, n_recv_values, &
       send_values, n_send_values)
    class(cartesian_pe_finder_t), intent(inout) :: this
    type(stack_i8_t), intent(inout) :: recv_values(0:this%pe_size-1)
    type(stack_i8_t), intent(inout) :: send_values(0:this%pe_size-1)
    integer, intent(inout) :: n_recv_values(0:this%pe_size-1)
    integer, intent(inout) :: n_send_values(0:this%pe_size-1)
    integer :: i, j, ierr
    integer(i8) :: idx
    integer(i8) , pointer :: sp(:)

    call MPI_Alltoall(n_send_values, 1, MPI_INTEGER, &
         n_recv_values, 1, MPI_INTEGER, this%comm, ierr)

    do i = 0, this%pe_size-1
       if (n_recv_values(i) .gt. 0) then
          if (this%recv_buf(i)%size .lt. n_recv_values(i)) then
             if (allocated(this%recv_buf(i)%data)) deallocate(this%recv_buf(i)%data)
             allocate(this%recv_buf(i)%data(n_recv_values(i)))
             this%recv_buf(i)%size = n_recv_values(i)
          end if
          call MPI_Irecv(this%recv_buf(i)%data, n_recv_values(i), MPI_INTEGER8, &
               i, 0, this%comm, this%recv_buf(i)%request, ierr)
       end if
    end do

    do i = 0, this%pe_size-1
       if (n_send_values(i) .gt. 0) then
          if (this%send_buf(i)%size .lt. n_send_values(i)) then
             if (allocated(this%send_buf(i)%data)) deallocate(this%send_buf(i)%data)
             allocate(this%send_buf(i)%data(n_send_values(i)))
             this%send_buf(i)%size = n_send_values(i)
          end if
          ! Copy the data to the send buffer
          sp => send_values(i)%array()
          do j = 1, n_send_values(i)
             this%send_buf(i)%data(j) = sp(j)
          end do
          call MPI_Isend(this%send_buf(i)%data, n_send_values(i), MPI_INTEGER8, &
               i, 0, this%comm, this%send_buf(i)%request, ierr)
       end if
    end do

    do i = 0, this%pe_size-1
       call recv_values(i)%clear()
    end do

    do i = 0, this%pe_size-1
       if (n_recv_values(i) .gt. 0) then
          call MPI_Wait(this%recv_buf(i)%request, this%recv_buf(i)%status, ierr)
          do j = 1, n_recv_values(i)
             idx = this%recv_buf(i)%data(j)
             call recv_values(i)%push(idx)
          end do
       end if
    end do
    do i = 0, this%pe_size-1
       if (n_send_values(i) .gt. 0) then
          call MPI_Wait(this%send_buf(i)%request, this%send_buf(i)%status, ierr)
       end if
    end do


  end subroutine send_recv_data
end module cartesian_pe_finder
