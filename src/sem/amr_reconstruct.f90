! Copyright (c) 2019-2025, The Neko Authors
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
!> The vector reconstruction/interpolation routines for AMR
module amr_reconstruct
  use neko_config, only : NEKO_BCKND_DEVICE
  use num_types, only : rp
  use utils, only : neko_error, neko_warning
  use amr_interpolate, only : amr_interpolate_t, amr_nchildren
  use device, only : device_sync
  use mxm_wrapper, only: mxm
  use mesh_manager_transfer, only : mesh_manager_transfer_t
  use, intrinsic :: iso_c_binding, only : c_ptr

  implicit none
  private

  !> Type for vector/field reconstruction
  type, public :: amr_reconstruct_t
     !> pointer mesh manager data transfer
     class(mesh_manager_transfer_t), pointer :: transfer
     !> AMR interpolation arrays
     type(amr_interpolate_t) :: interpolate

     !> Work space for single element refinement
     real(rp), allocatable, dimension(:, :, :, :) :: tmp

     ! Arrays sizes
     !> Old number of elements
     integer :: nold
     !> New number of elements
     integer :: nnew
     !> Refinement mapping size
     integer :: nref
     !> Coarsening mapping size
     integer :: ncrs
     !> Refinement mapping (element and child position)
     integer, dimension(:, :), allocatable :: rmap
     !> Coarsening mapping (element and child position)
     integer, dimension(:, :), allocatable :: cmap
     !> Output vector for refinement
     real(rp), allocatable, dimension(:, :, :, :) :: vout
     !> Vector with additional data for coarsening
     real(rp), allocatable, dimension(:, :, :, :, :) :: vcrs
     !> Is local mesh changed
     logical :: ifchange

     !
     ! Device pointers (if present)
     !
     ! should tmp be on device?
   contains
     !> Initialise type
     procedure, pass(this) :: init => amr_reconstruct_init
     !> Free type
     procedure, pass(this) :: free => amr_reconstruct_free
     !> Get refinement/coarsening mapping
     procedure, pass(this) :: map_get => amr_reconstruct_map_get
     !> Free refinement/coarsening mapping
     procedure, pass(this) :: map_free => amr_reconstruct_map_free
     !> Perform refinement/coarsening on a single vector
     procedure, pass(this) :: refine_coarsen => amr_reconstruct_refine_coarsen
     !> Map single coarse to fine element
     procedure, pass(this) :: map_c2f => amr_reconstruct_map_c2f
     !> Map single fine to coarse element
     procedure, pass(this) :: map_f2c => amr_reconstruct_map_f2c
  end type amr_reconstruct_t

contains
  !> Initialise type
  !! @param[in]  transfer     mesh manager data transfer type
  !! @param[in]  gdim         geometrical mesh dimension
  !! @param[in]  lx           polynomial order + 1
  subroutine amr_reconstruct_init(this, transfer, gdim, lx)
    class(amr_reconstruct_t), intent(inout) :: this
    class(mesh_manager_transfer_t), target, intent(in) :: transfer
    integer, intent(in) :: gdim, lx

    call this%free()

    ! mesh manager data transfer
    this%transfer => transfer

    call this%interpolate%init(gdim, lx)

    associate(Xh => this%interpolate%Xh)

      ! work space
      allocate(this%tmp(Xh%lx, Xh%ly, Xh%lz, 3))

      if (NEKO_BCKND_DEVICE .eq. 1) then
         ! should tmp be on device?
      end if

    end associate

    call device_sync()

  end subroutine amr_reconstruct_init

  !> Free type
  subroutine amr_reconstruct_free(this)
    class(amr_reconstruct_t), intent(inout) :: this

    this%transfer => NULL()

    call this%interpolate%free()

    this%nold = 0
    this%nnew = 0
    this%nref = 0
    this%ncrs = 0

    if (allocated(this%tmp)) deallocate(this%tmp)

    if (allocated(this%rmap)) deallocate(this%rmap)
    if (allocated(this%cmap)) deallocate(this%cmap)
    if (allocated(this%vout)) deallocate(this%vout)
    if (allocated(this%vcrs)) deallocate(this%vcrs)

    this%ifchange = .false.

    !
    ! Cleanup the device (if present)
    !
    ! should tmp be on device?

  end subroutine amr_reconstruct_free

  !> Get refinement/coarsening mapping
  subroutine amr_reconstruct_map_get(this)
    class(amr_reconstruct_t), intent(inout) :: this
    integer :: nchildren

    call this%transfer%vector_map(this%nold, this%nnew, this%nref, this%ncrs, &
         this%rmap, this%cmap, nchildren, this%ifchange, &
         this%interpolate%Xh%lx, this%interpolate%Xh%ly, this%interpolate%Xh%lz)

    ! sanity check
    if (amr_nchildren .ne. nchildren) &
         call neko_error('Inconsistent children number for mesh manager and &
         &refinement routines.')

    associate(Xh => this%interpolate%Xh)
    ! get coarsening vector
      if (this%ifchange) allocate(this%vcrs(Xh%lx, Xh%ly, Xh%lz, &
           amr_nchildren, this%ncrs))
    end associate
  end subroutine amr_reconstruct_map_get

  !> Free refinement/coarsening mapping
  subroutine amr_reconstruct_map_free(this)
    class(amr_reconstruct_t), intent(inout) :: this

    call this%transfer%vector_map_free()

    this%nold = 0
    this%nnew = 0
    this%nref = 0
    this%ncrs = 0

    if (allocated(this%rmap)) deallocate(this%rmap)
    if (allocated(this%cmap)) deallocate(this%cmap)
    if (allocated(this%vout)) deallocate(this%vout)
    if (allocated(this%vcrs)) deallocate(this%vcrs)

    this%ifchange = .false.

  end subroutine amr_reconstruct_map_free

  !> Perform refinement/coarsening on a single vector
  !! @param[inout] vec    vector for refinement/coarsening
  !! @param[inout] vec_d  device pointer to vector
  subroutine amr_reconstruct_refine_coarsen(this, vec, vec_d)
    class(amr_reconstruct_t), intent(inout) :: this
    real(rp), dimension(:, :, :, :), allocatable, intent(inout) :: vec
    type(c_ptr), optional, intent(inout) :: vec_d
    integer :: il, jl, itmp
    integer, dimension(3) :: ch_pos

    if (this%ifchange) then
       associate(int => this%interpolate)
         ! check device pointer
         if (present(vec_d) .and. NEKO_BCKND_DEVICE .eq. 1) then
            call neko_error('AMR reconstruction; nothing done for device')
         end if

         ! CPU part
         allocate(this%vout(int%Xh%lx, int%Xh%ly, int%Xh%lz, this%nnew))
         call this%transfer%vector_constr(vec, this%vout, this%vcrs)
         ! refinement
         if (this%nref .gt. 0) then
            do il = 1, this%nref
               ! get child position with respect to the parent
               ! (0...amr_nchildren - 1)
               itmp = this%rmap(2, il)
               ch_pos(3) = itmp / 4 + 1 ! z position
               ch_pos(2) = mod(itmp / 2, 2) + 1 ! y position
               ch_pos(1) = mod(itmp, 2) + 1 ! x position
               call this%map_c2f(int%gdim, ch_pos, &
                    this%vout(:, :, :, this%rmap(1, il)))
            end do
         end if
         ! coarsening
         if (this%ncrs .gt. 0) then
            do il = 1, this%ncrs
               this%vout(:, :, :, this%cmap(1, il)) = 0.0_rp
               do jl = 1, amr_nchildren
                  ! get child position with respect to the parent
                  ! (0...amr_nchildren - 1)
                  itmp = this%cmap(1 + jl, il)
                  ch_pos(3) = itmp / 4 + 1 ! z position
                  ch_pos(2) = mod(itmp / 2, 2) + 1 ! y position
                  ch_pos(1) = mod(itmp, 2) + 1 ! x position
                  call this%map_f2c(int%gdim, ch_pos, &
                       this%vcrs(:, :, :, jl, il), this%tmp(:, :, :, 3))
                  ! accumulate result
                  this%vout(:, :, :, this%cmap(1, il)) = &
                       this%vout(:, :, :, this%cmap(1, il)) + &
                       this%tmp(:, :, :, 3)
               end do
               ! reduce multiplicity
               this%vout(:, :, :, this%cmap(1, il)) = &
                    this%vout(:, :, :, this%cmap(1, il)) * int%el_mult(:, :, :)
            end do
         end if
         call move_alloc(this%vout, vec)

         ! check device pointer
         if (present(vec_d) .and. NEKO_BCKND_DEVICE .eq. 1) then
            call neko_error('AMR reconstruction; nothing done for device; 2')
         end if
       end associate
    end if

  end subroutine amr_reconstruct_refine_coarsen

  !> @brief Map a single coarse element to a fine one
  !! @param[in]    gdim    geometrical dimension
  !! @param[in]    ch_pos  child position
  !! @param[inout] vcf     coarse (input) and fine (output) element vector
  subroutine amr_reconstruct_map_c2f(this, gdim, ch_pos, vcf)
    class(amr_reconstruct_t), intent(inout) :: this
    integer, intent(in) :: gdim
    integer, dimension(3), intent(in) :: ch_pos
    real(rp), dimension(:,:,:), intent(inout) :: vcf
    integer :: iz

    associate(int => this%interpolate)
      if (gdim == 3) then ! 3D
         call mxm(int%x_cr2fn(:,:, ch_pos(1)), int%Xh%lx, vcf, &
              & int%Xh%lx, this%tmp(:,:,:, 1), int%Xh%lyz)
         do iz = 1, int%Xh%lz
            call mxm(this%tmp(:,:, iz, 1), int%Xh%lx,&
                 & int%y_cr2fnT(:,:, ch_pos(2)),&
                 & int%Xh%ly, this%tmp(:,:, iz, 2), int%Xh%ly)
         end do
         call mxm(this%tmp(:,:,:, 2), int%Xh%lxy,&
              & int%z_cr2fnT(:,:, ch_pos(3)),&
              & int%Xh%lz, vcf, int%Xh%lz)
      else ! 2D
         call mxm(int%x_cr2fn(:,:, ch_pos(1)), int%Xh%lx, vcf, int%Xh%lz,&
              & this%tmp(:,:,1,1), int%Xh%lyz)
         call mxm(this%tmp(:,:,1,1), int%Xh%lx, int%y_cr2fnT(:,:, ch_pos(2)),&
              & int%Xh%ly, vcf, int%Xh%ly)
      end if
    end associate

  end subroutine amr_reconstruct_map_c2f

  !> @brief Map a single fine element to a coarse one
  !! @param[in]    gdim    geometrical dimension
  !! @param[in]    ch_pos  child position
  !! @param[in]    vf      fine element vector
  !! @param[out]   vc      coarse element vector
  subroutine amr_reconstruct_map_f2c(this, gdim, ch_pos, vf, vc)
    class(amr_reconstruct_t), intent(inout) :: this
    integer, intent(in) :: gdim
    integer, dimension(3), intent(in) :: ch_pos
    real(rp), dimension(:,:,:), intent(in) :: vf
    real(rp), dimension(:,:,:), intent(out) :: vc
    integer :: iz

    associate(int => this%interpolate)
      if (gdim == 3) then ! 3D
         call mxm(int%x_fn2cr(:,:, ch_pos(1)), int%Xh%lx, vf, &
              & int%Xh%lx, this%tmp(:,:,:, 1), int%Xh%lyz)
         do iz = 1, int%Xh%lz
            call mxm(this%tmp(:,:, iz, 1), int%Xh%lx,&
                 & int%y_fn2crT(:,:, ch_pos(2)),&
                 & int%Xh%ly, this%tmp(:,:, iz, 2), int%Xh%ly)
         end do
         call mxm(this%tmp(:,:,:, 2), int%Xh%lxy,&
              & int%z_fn2crT(:,:, ch_pos(3)),&
              & int%Xh%lz, vc, int%Xh%lz)
      else ! 2D
         call mxm(int%x_fn2cr(:,:, ch_pos(1)), int%Xh%lx, vf, int%Xh%lz,&
              & this%tmp(:,:,1,1), int%Xh%lyz)
         call mxm(this%tmp(:,:,1,1), int%Xh%lx, int%y_fn2crT(:,:, ch_pos(2)),&
              & int%Xh%ly, vc, int%Xh%ly)
      end if
    end associate

  end subroutine amr_reconstruct_map_f2c

end module amr_reconstruct
