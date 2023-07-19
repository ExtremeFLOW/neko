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
!> Implements probes.
!! @note This modules uses functions from `gslib`, namely `findpts_setup`,
!! `findpts`, and `findpts_eval`. A full description of these subroutines can
!! be found at https://github.com/Nek5000/gslib/blob/master/src/findpts.c
module probes
  use num_types, only: rp
  use matrix, only: matrix_t
  use logger, only: neko_log, LOG_SIZE
  use utils, only: neko_error
  use comm
  use mpi_types
  use field, only: field_t
  use coefs, only: coef_t
  implicit none

  character(len=LOG_SIZE) :: log_buf ! For logging status

  !@todo reformat the rst and xyz variables with array of points
  type probes_t

     integer :: handle !< handle to pass in each findpts call

     integer :: n_probes !< Number of probes
     integer :: n_fields !< Number of output fields
     integer, allocatable :: proc_owner(:) !< List of owning processes
     integer, allocatable :: el_owner(:)   !< List of owning elements
     integer, allocatable :: dist2(:) !< Distance squared between original
     ! and interpolated point (in xyz space)
     integer, allocatable :: error_code(:) !< Error code for each point
     real(kind=rp), allocatable :: rst(:) !< r,s,t coordinates, findpts format
     real(kind=rp), allocatable :: xyz(:,:) !< x,y,z coordinates, findpts format
     real(kind=rp), allocatable :: out_fields(:,:) !< interpolated fields

     contains
       procedure, pass(this) :: init => probes_init
       procedure, pass(this) :: free => probes_free
       procedure, pass(this) :: show => probes_show
       procedure, pass(this) :: setup => probes_setup
       procedure, pass(this) :: map => probes_map
       procedure, pass(this) :: interpolate => probes_interpolate

    end type probes_t

contains

  !> Initialize probes object
  !! @param n_probes Number of probes
  !! @param n_fields Number of fields to interpolate
  subroutine probes_init(this, n_probes, n_fields)
    class(probes_t), intent(inout) :: this
    integer, intent(in) :: n_probes
    integer, intent(in) :: n_fields

    this%n_probes = n_probes
    this%n_fields = n_fields

    allocate(this%xyz(3, n_probes))
    allocate(this%proc_owner(n_probes))
    allocate(this%el_owner(n_probes))
    allocate(this%dist2(n_probes))
    allocate(this%error_code(n_probes))
    allocate(this%rst(3*n_probes))
    allocate(this%out_fields(n_fields, n_probes))

  end subroutine probes_init

  !> Destructor
  subroutine probes_free(this)
    class(probes_t), intent(inout) :: this

    if (allocated(this%xyz))        deallocate(this%xyz)
    if (allocated(this%rst))        deallocate(this%rst)
    if (allocated(this%proc_owner)) deallocate(this%proc_owner)
    if (allocated(this%el_owner))   deallocate(this%el_owner)
    if (allocated(this%dist2))       deallocate(this%dist2)
    if (allocated(this%error_code)) deallocate(this%error_code)
    if (allocated(this%out_fields)) deallocate(this%out_fields)

    call fgslib_findpts_free(this%handle)

  end subroutine probes_free

  !> Show probes
  subroutine probes_show(this)
    class(probes_t), intent(in) :: this
    integer :: i

    write(log_buf, '(A,I6)') "Number of probes: ", this%n_probes
    call neko_log%message(log_buf)
    call neko_log%message("xyz-coordinates:")
    do i=1,this%n_probes
       write(log_buf, '("(",F10.6,",",F10.6,",",F10.6,")")') this%xyz(:,i)
       call neko_log%message(log_buf)
    end do

  end subroutine probes_show

  !> Setup for the mapping process.
  subroutine probes_setup(this, coef)
    class(probes_t), intent(inout) :: this
    type(coef_t), intent(in) :: coef

    real(kind=rp) :: tolerance
    integer :: lx, ly, lz, nelv, max_pts_per_iter

    tolerance = 5d-13 ! Tolerance for Newton iterations
    lx = coef%xh%lx
    ly = coef%xh%ly
    lz = coef%xh%lz
    nelv = coef%msh%nelv
    max_pts_per_iter = 128 ! number of points to iterate on simultaneously

    call fgslib_findpts_setup(this%handle, &
         NEKO_COMM, pe_size, &
         coef%msh%gdim, &
         coef%dof%x, coef%dof%y, coef%dof%z, & ! Physical nodal values
         lx, ly, lz, nelv, & ! Mesh dimensions
         2*lx, 2*ly, 2*lz, & ! Mesh size for bounding box computation
         0.01, & ! relative size to expand bounding boxes by
         lx*ly*lz*nelv, lx*ly*lz*nelv, & ! local/global hash mesh sizes
         max_pts_per_iter, tolerance)

  end subroutine probes_setup

  !> Maps `x,y,z` to `r,s,t` coordinates. The result of the mapping for each
  !! point can be found in the following variables:
  !! - `error_code`: returns `0` if point found, `1` if closest point on a border
  !! (check dist2), `2` if not found
  !! - `dist2`: distance squared (used to compare the points found by each
  !! processor)
  subroutine probes_map(this, coef)
    class(probes_t), intent(inout) :: this
    type(coef_t), intent(in) :: coef

    call fgslib_findpts(this%handle, &
         this%error_code, 1, &
         this%proc_owner, 1, &
         this%el_owner, 1, &
         this%rst, coef%msh%gdim, &
         this%dist2, 1, &
         this%xyz(1,1), coef%msh%gdim, &
         this%xyz(2,1), coef%msh%gdim, &
         this%xyz(3,1), coef%msh%gdim, this%n_probes)

  end subroutine probes_map

  !> Interpolate each probe from its `r,s,t` coordinates.
  !! @note At the moment the interpolation is performed on fields
  !! `u,v,w,p`.
  subroutine probes_interpolate(this, coef, u, v, w, p)
    class(probes_t), intent(inout) :: this
    type(coef_t), intent(in) :: coef
    type(field_t), intent(in) :: u
    type(field_t), intent(in) :: v
    type(field_t), intent(in) :: w
    type(field_t), intent(in) :: p

    call fgslib_findpts_eval(this%handle, this%out_fields(1,1), this%n_fields, &
         this%error_code, 1, &
         this%proc_owner, 1, &
         this%el_owner, 1, &
         this%rst, coef%msh%gdim, &
         this%n_probes, coef%dof%x)

    call fgslib_findpts_eval(this%handle, this%out_fields(2,1), this%n_fields, &
         this%error_code, 1, &
         this%proc_owner, 1, &
         this%el_owner, 1, &
         this%rst, coef%msh%gdim, &
         this%n_probes, coef%dof%y)

    call fgslib_findpts_eval(this%handle, this%out_fields(3,1), this%n_fields, &
         this%error_code, 1, &
         this%proc_owner, 1, &
         this%el_owner, 1, &
         this%rst, coef%msh%gdim, &
         this%n_probes, coef%dof%z)

    call fgslib_findpts_eval(this%handle, this%out_fields(this%n_fields,1), 4, &
         this%error_code, 1, &
         this%proc_owner, 1, &
         this%el_owner, 1, &
         this%rst, coef%msh%gdim, &
         this%n_probes, u%x)

  end subroutine probes_interpolate

end module probes
