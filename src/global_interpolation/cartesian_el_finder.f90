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

module cartesian_el_finder
  use num_types, only: rp, xp, dp
  use neko_config, only : NEKO_BCKND_DEVICE
  use el_finder, only: el_finder_t
  use space, only: space_t
  use stack, only: stack_i4_t
  use tuple, only: tuple_i4_t
  use point, only: point_t
  use htable, only: htable_i4_t
  use utils, only: linear_index
  use mpi_f08, only: MPI_Wtime
  use tensor_cpu, only: tnsr3d_cpu
  use fast3d, only: setup_intp
  implicit none
  private

  !> Implements global interpolation for arbitrary points in the domain.
  type, public, extends(el_finder_t) :: cartesian_el_finder_t
     !> Structure to find rank candidates
     type(stack_i4_t), allocatable :: el_map(:)
     integer :: n_boxes, nel
     real(kind=xp) :: max_x, min_x
     real(kind=xp) :: max_y, min_y
     real(kind=xp) :: max_z, min_z
     real(kind=xp) :: x_res, y_res, z_res
     real(kind=xp) :: padding

   contains
     procedure, pass(this) :: init => cartesian_el_finder_init
     procedure, pass(this) :: free => cartesian_el_finder_free
     procedure, pass(this) :: find => cartesian_el_finder_find_candidates
     procedure, pass(this) :: find_batch => cartesian_el_finder_find_candidates_batch
     procedure, pass(this) :: compute_idx => cartesian_el_finder_compute_idx
     procedure, pass(this) :: compute_3idx => cartesian_el_finder_compute_xyz_idxs
  end type cartesian_el_finder_t

contains

  subroutine cartesian_el_finder_init(this, x, y, z, nel, Xh, n_boxes, padding)
    class(cartesian_el_finder_t), intent(inout) :: this
    type(space_t), intent(in) :: Xh
    real(kind=rp), intent(in), target :: x(:), y(:), z(:)
    integer, intent(in) :: nel
    integer, intent(in) :: n_boxes
    real(kind=dp), intent(in) :: padding
    type(htable_i4_t) :: marked_box
    integer :: i, j, k, e
    integer :: el_npts, el_x_npts, el_y_npts, el_z_npts
    real(kind=rp) :: el_x_max, el_x_min
    real(kind=rp) :: el_y_max, el_y_min
    real(kind=rp) :: el_z_max, el_z_min
    real(kind=rp) :: el_x(Xh%lxyz), el_y(Xh%lxyz), el_z(Xh%lxyz)
    integer :: lxyz
    real(kind=rp), allocatable :: el_rst(:)
    real(kind=rp), allocatable :: res_x(:), res_y(:), res_z(:)
    real(kind=rp), allocatable :: interp_mat(:,:), interp_matT(:,:)
    integer :: el_idx
    integer :: htable_data
    real(kind=rp) :: center_x, center_y, center_z, r_res
    integer :: current_size
    real(kind=rp) :: time_start, time_end
    integer :: min_id(3), max_id(3)
    integer :: lin_idx
    integer :: i2, j2, k2
    real(kind=rp) :: min_bb_x, max_bb_x
    real(kind=rp) :: min_bb_y, max_bb_y
    real(kind=rp) :: min_bb_z, max_bb_z




    call this%free()
    allocate(this%el_map(n_boxes**3))
    do i = 1, n_boxes**3
       call this%el_map(i)%init()
    end do

    time_start = MPI_Wtime()
    current_size = 0
    this%padding = padding
    this%max_x = maxval(x(1:nel*Xh%lxyz))
    this%min_x = minval(x(1:nel*Xh%lxyz))
    this%max_y = maxval(y(1:nel*Xh%lxyz))
    this%min_y = minval(y(1:nel*Xh%lxyz))
    this%max_z = maxval(z(1:nel*Xh%lxyz))
    this%min_z = minval(z(1:nel*Xh%lxyz))
    
    center_x = (this%max_x + this%min_x) /2.0_xp
    center_y = (this%max_y + this%min_y) /2.0_xp
    center_z = (this%max_z + this%min_z) /2.0_xp

    this%max_x = this%max_x - center_x
    this%max_y = this%max_y - center_y
    this%max_z = this%max_z - center_z
    this%min_x = this%min_x - center_x
    this%min_y = this%min_y - center_y
    this%min_z = this%min_z - center_z
    this%max_x = this%max_x * (1.0_xp+this%padding) + center_x
    this%max_y = this%max_y * (1.0_xp+this%padding) + center_y
    this%max_z = this%max_z * (1.0_xp+this%padding) + center_z
    this%min_x = this%min_x * (1.0_xp+this%padding) + center_x
    this%min_y = this%min_y * (1.0_xp+this%padding) + center_y
    this%min_z = this%min_z * (1.0_xp+this%padding) + center_z

    

    !Resulting resolution
    this%x_res = (this%max_x - this%min_x) / real(n_boxes, xp)
    this%y_res = (this%max_y - this%min_y) / real(n_boxes, xp)
    this%z_res = (this%max_z - this%min_z) / real(n_boxes, xp)

    this%nel = nel
    this%n_boxes = n_boxes
    call marked_box%init(nel)
    lxyz = Xh%lxyz
    !do e = 1, nel
    !   call marked_box%clear()
    !   el_x_max = maxval(x((e-1)*lxyz+1:e*lxyz))
    !   el_x_min = minval(x((e-1)*lxyz+1:e*lxyz))
    !   el_y_max = maxval(y((e-1)*lxyz+1:e*lxyz))
    !   el_y_min = minval(y((e-1)*lxyz+1:e*lxyz))
    !   el_z_max = maxval(z((e-1)*lxyz+1:e*lxyz))
    !   el_z_min = minval(z((e-1)*lxyz+1:e*lxyz))

    !   el_x_npts = 2*max(int((el_x_max - el_x_min) / this%x_res),Xh%lx)
    !   el_y_npts = 2*max(int((el_y_max - el_y_min) / this%y_res),Xh%ly)
    !   el_z_npts = 2*max(int((el_z_max - el_z_min) / this%z_res),Xh%lz)

    !   el_npts = max(el_x_npts, el_y_npts, el_z_npts)
    !   ! If this happens often, we can split the element into multiple rst parts
    !   ! This check is to avoid allocating too much memory on the stack
    !   if (el_npts .gt. 128) then
    !      print *, "Warning, in cartesian_el_finder el_npts: ", el_npts, "set to 128"
    !      el_npts = 128
    !   end if

    !   !move it to close to origo
    !   center_x = (el_x_max + el_x_min) / 2.0_xp
    !   center_y = (el_y_max + el_y_min) / 2.0_xp
    !   center_z = (el_z_max + el_z_min) / 2.0_xp
    !   !Maybe this is a stupid way to do it

    !   el_x = x((e-1)*lxyz+1:e*lxyz) - center_x
    !   el_y = y((e-1)*lxyz+1:e*lxyz) - center_y
    !   el_z = z((e-1)*lxyz+1:e*lxyz) - center_z
    !   el_x = el_x * (1.0+padding) + center_x
    !   el_y = el_y * (1.0+padding) + center_y
    !   el_z = el_z * (1.0+padding) + center_z
    !   r_res = 2.0_xp / real(el_npts-1, xp)
    !   if (el_npts .gt. current_size) then
    !      current_size = el_npts
    !      if (allocated(el_rst)) deallocate(el_rst)
    !      if (allocated(res_x)) deallocate(res_x)
    !      if (allocated(res_y)) deallocate(res_y)
    !      if (allocated(res_z)) deallocate(res_z)
    !      if (allocated(interp_mat)) deallocate(interp_mat)
    !      if (allocated(interp_matT)) deallocate(interp_matT)
    !      allocate(res_x(el_npts**3))
    !      allocate(res_y(el_npts**3))
    !      allocate(res_z(el_npts**3))
    !      allocate(interp_mat(Xh%lx, el_npts), &
    !               interp_matT(el_npts, Xh%lx))
    !      allocate(el_rst(el_npts))
    !   end if
    !   do i = 1, el_npts
    !      el_rst(i) = -1_xp+r_res*(i-1)
    !   end do

    !   call setup_intp(interp_mat, interp_matT, &
    !   el_rst, Xh%zg, el_npts, Xh%lx, 0)

    !   call tnsr3d_cpu(res_x, el_npts, el_x, &
    !               Xh%lx,interp_mat, &
    !               interp_matT, interp_matT, 1)
    !   call tnsr3d_cpu(res_y, el_npts, el_y, &
    !               Xh%lx,interp_mat, &
    !               interp_matT, interp_matT, 1)
    !   call tnsr3d_cpu(res_z, el_npts, el_z, &
    !               Xh%lx,interp_mat, &
    !               interp_matT, interp_matT, 1)


    !   do i = 1, el_npts**3
    !      el_idx = this%compute_idx(real(res_x(i),xp), &
    !               real(res_y(i),xp), real(res_z(i),xp))
    !      if (el_idx .ge. 1 .and. el_idx .lt. this%n_boxes**3) then
    !         if (marked_box%get(el_idx,htable_data) .ne. 0)then
    !            call marked_box%set(el_idx, htable_data)
    !            call this%el_map(el_idx)%push(e)
    !         end if
    !      end if
    !   end do
    !end do
    do e = 1, nel
       call marked_box%clear()
       el_x_max = maxval(x((e-1)*lxyz+1:e*lxyz))
       el_x_min = minval(x((e-1)*lxyz+1:e*lxyz))
       el_y_max = maxval(y((e-1)*lxyz+1:e*lxyz))
       el_y_min = minval(y((e-1)*lxyz+1:e*lxyz))
       el_z_max = maxval(z((e-1)*lxyz+1:e*lxyz))
       el_z_min = minval(z((e-1)*lxyz+1:e*lxyz))

       !move it to close to origo
       !center_x = (el_x_max + el_x_min) / 2.0_xp
       !center_y = (el_y_max + el_y_min) / 2.0_xp
       !center_z = (el_z_max + el_z_min) / 2.0_xp
       center_x = 0d0
       center_y = 0d0
       center_z = 0d0
       do i = 1, Xh%lxyz
          center_x = center_x + x((e-1)*lxyz+i)
          center_y = center_y + y((e-1)*lxyz+i)
          center_z = center_z + z((e-1)*lxyz+i)
       end do
       center_x = center_x / Xh%lxyz
       center_y = center_y / Xh%lxyz
       center_z = center_z / Xh%lxyz

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
                lin_idx = linear_index(i,j,k, 1, Xh%lx, Xh%lx, Xh%lx)
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


                min_id = this%compute_3idx(min_bb_x, min_bb_y, min_bb_z)
                max_id = this%compute_3idx(max_bb_x, max_bb_y, max_bb_z)
                do i2 = min_id(1), max_id(1)
                   do j2 = min_id(2), max_id(2)
                      do k2 = min_id(3), max_id(3)
                         if (i2 .ge. 1 .and. i2 .le. this%n_boxes .and. &
                             j2 .ge. 1 .and. j2 .le. this%n_boxes .and. &
                             k2 .ge. 1 .and. k2 .le. this%n_boxes) then
                            el_idx = linear_index(i2, j2, k2, 1, &
                            this%n_boxes, this%n_boxes, this%n_boxes)
                            if (marked_box%get(el_idx,htable_data) .ne. 0)then
                               call marked_box%set(el_idx, htable_data)
                               call this%el_map(el_idx)%push(e)
                            end if
                         end if
                      end do
                   end do
                end do
             end do
          end do
       end do
    end do
    call marked_box%free()
    !print *, "Time for cartesian_el_finder_init: ", MPI_Wtime() - time_start
  end subroutine cartesian_el_finder_init

  subroutine cartesian_el_finder_free(this)
    class(cartesian_el_finder_t), intent(inout) :: this

  end subroutine cartesian_el_finder_free

  function cartesian_el_finder_compute_idx(this, x, y, z) result(idx)
    class(cartesian_el_finder_t), intent(in) :: this
    real(kind=rp), intent(in) :: x, y, z
    integer :: idx
    integer :: ids(3)

    ids = this%compute_3idx(x, y, z)
    idx = linear_index(ids(1), ids(2), ids(3), 1, &
           this%n_boxes, this%n_boxes, this%n_boxes)

  end function cartesian_el_finder_compute_idx

  function cartesian_el_finder_compute_xyz_idxs(this, x, y, z) result(idxs)
    class(cartesian_el_finder_t), intent(in) :: this
    real(kind=rp), intent(in) :: x, y, z
    integer :: idxs(3)
    integer :: x_id, y_id, z_id

    x_id = int(real(x - this%min_x,xp) / this%x_res)
    y_id = int(real(y - this%min_y,xp) / this%y_res)
    z_id = int(real(z - this%min_z,xp) / this%z_res)
    if (x_id .eq. -1) then
       x_id = 0
    end if
    if (x_id .eq. this%n_boxes) then
       x_id = this%n_boxes - 1
    end if
    if (y_id .eq. -1) then
       y_id = 0
    end if
    if (y_id .eq. this%n_boxes) then
       y_id = this%n_boxes - 1
    end if
    if (z_id .eq. -1) then
       z_id = 0
    end if
    if (z_id .eq. this%n_boxes) then
       z_id = this%n_boxes - 1
    end if
    x_id = x_id + 1
    y_id = y_id + 1
    z_id = z_id + 1

    idxs(1) = x_id
    idxs(2) = y_id
    idxs(3) = z_id
  end function cartesian_el_finder_compute_xyz_idxs

  subroutine cartesian_el_finder_find_candidates(this, my_point, el_candidates)
    class(cartesian_el_finder_t), intent(inout) :: this
    type(point_t), intent(in) :: my_point
    type(stack_i4_t), intent(inout) :: el_candidates
    integer :: idx, i, stupid_intent
    integer, pointer :: el_cands(:)

    call el_candidates%clear()
    idx = this%compute_idx(my_point%x(1),my_point%x(2),my_point%x(3))
    el_cands => this%el_map(idx)%array()
    do i = 1, this%el_map(idx)%size()
       stupid_intent = el_cands(i) - 1
       call el_candidates%push(stupid_intent)
    end do
  end subroutine cartesian_el_finder_find_candidates

  ! In order to get more cache hits
  subroutine cartesian_el_finder_find_candidates_batch(this, points, n_points, &
             all_el_candidates, n_el_cands)
    class(cartesian_el_finder_t), intent(inout) :: this
    integer, intent(in) :: n_points
    real(kind=rp), intent(in) :: points(3,n_points)
    type(stack_i4_t), intent(inout) :: all_el_candidates
    integer, intent(inout) :: n_el_cands(n_points)
    integer :: i, j, stupid_intent
    integer :: idx3(3), idx
    integer, pointer :: el_cands(:)

    call all_el_candidates%clear()
    n_el_cands = 0

    do i = 1, n_points
       idx3 = this%compute_3idx(real(points(1,i),xp),&
             real(points(2,i),xp),real(points(3,i),xp))
       if (idx3(1) .ge. 1 .and. idx3(1) .le. this%n_boxes .and. &
           idx3(2) .ge. 1 .and. idx3(2) .le. this%n_boxes .and. &
           idx3(3) .ge. 1 .and. idx3(3) .le. this%n_boxes) then
       idx = this%compute_idx(real(points(1,i),xp),&
             real(points(2,i),xp),real(points(3,i),xp))
       el_cands => this%el_map(idx)%array()
       do j = 1, this%el_map(idx)%size()
          stupid_intent = el_cands(j) - 1
          call all_el_candidates%push(stupid_intent) !< OBS c indexing
       end do
       n_el_cands(i) = this%el_map(idx)%size()
       end if
    end do

  end subroutine cartesian_el_finder_find_candidates_batch

end module cartesian_el_finder
