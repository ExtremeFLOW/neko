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

  !> Type for field reconstruction for given lx
  type :: amr_reconstruct_lx_t
     !> Initialisation flag
     logical :: ifinit = .false.
     !> Polynomial order + 1
     integer :: lx
     !> 3D element size
     integer :: lxyz
     !> AMR interpolation arrays
     type(amr_interpolate_t) :: interpolate
     !> Work space for single element refinement
     real(rp), allocatable, dimension(:, :, :, :) :: tmp
     !> Output vector for refinement
     real(rp), allocatable, dimension(:, :, :, :) :: vout
     !> Vector with additional data for coarsening
     real(rp), allocatable, dimension(:, :, :, :, :) :: vcrs
     !> Receive buffer
     real(rp), allocatable, dimension(:, :, :, :) :: buff_rcv
     !> Send buffer
     real(rp), allocatable, dimension(:, :, :, :) :: buff_snd
     !
     ! Device pointers (if present)
     !
     ! should some stuff be on device?
   contains
     !> Initialise type
     procedure, pass(this) :: init => amr_reconstruct_lx_init
     !> Free type
     procedure, pass(this) :: free => amr_reconstruct_lx_free
     !> Initialise type
     procedure, pass(this) :: buff_get => amr_reconstruct_lx_buff_get
     !> Free type
     procedure, pass(this) :: buff_free => amr_reconstruct_lx_buff_free
     !> Perform refinement/coarsening on a single vector for given lx
     procedure, pass(this) :: refine_coarsen => &
          amr_reconstruct_lx_refine_coarsen
     !> Map single coarse to fine element
     procedure, pass(this) :: map_c2f => amr_reconstruct_lx_map_c2f
     !> Map single fine to coarse element
     procedure, pass(this) :: map_f2c => amr_reconstruct_lx_map_f2c
  end type amr_reconstruct_lx_t

  !> Type for vector/field reconstruction
  type, public :: amr_reconstruct_t
     !> Primary number of grid points in 1D
     integer :: lx
     !> Primary number of grid points in 2D
     integer :: lxy
     !> Primary number of grid points in 3D
     integer :: lxyz
     !> pointer mesh manager data transfer
     class(mesh_manager_transfer_t), pointer :: transfer
     !> Size of interpolation array
     integer :: int_lx
     !> AMR interpolation operators for various lx
     type(amr_reconstruct_lx_t), dimension(:), allocatable :: recon

     ! Arrays sizes
     !> Old number of elements
     integer :: nold
     !> New number of elements
     integer :: nnew
     !> Refinement mapping size
     integer :: nref
     !> Coarsening mapping size
     integer :: ncrs
     !> Receive buffer size
     integer :: nrcv
     !> Send buffer size
     integer :: nsnd
     !> Refinement mapping (element and child position)
     integer, dimension(:, :), allocatable :: rmap
     !> Coarsening mapping (element and child position)
     integer, dimension(:, :), allocatable :: cmap
     !> Is local mesh changed
     logical :: ifchange

     !
     ! Device pointers (if present)
     !
     ! should some stuff be on device?
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
  end type amr_reconstruct_t

contains
  !> Initialise type
  !! @param[in]  gdim         geometrical mesh dimension
  !! @param[in]  lx           polynomial order + 1
  subroutine amr_reconstruct_lx_init(this, gdim, lx)
    class(amr_reconstruct_lx_t), intent(inout) :: this
    integer, intent(in) :: gdim, lx

    call this%free()

    this%ifinit = .true.
    ! I assume lx = ly = lz.
    this%lx = lx
    this%lxyz = lx * lx * lx

    call this%interpolate%init(gdim, lx)

    ! work space
    allocate(this%tmp(this%lx, this%lx, this%lx, 3))

!    if (NEKO_BCKND_DEVICE .eq. 1) then
!       ! should tmp be on device?
!    end if
!    call device_sync()

  end subroutine amr_reconstruct_lx_init

  !> Free type
  subroutine amr_reconstruct_lx_free(this)
    class(amr_reconstruct_lx_t), intent(inout) :: this

    call this%interpolate%free()

    this%ifinit = .false.
    this%lx = 0
    this%lxyz = 0

    if (allocated(this%tmp)) deallocate(this%tmp)
    call this%buff_free()

    !
    ! Cleanup the device (if present)
    !
    ! should tmp be on device?

  end subroutine amr_reconstruct_lx_free

  !> Allocate buffer arrays
  !! @param[in]    nrcv      Receive buffer size
  !! @param[in]    nsnd      Send buffer size
  !! @param[in]    ncrs      Coarsening mapping size
  subroutine amr_reconstruct_lx_buff_get(this, nrcv, nsnd, ncrs)
    class(amr_reconstruct_lx_t), intent(inout) :: this
    integer, intent(in) :: nrcv, nsnd, ncrs

    if (nrcv .gt. 0) allocate(this%buff_rcv(this%lx, this%lx, this%lx, nrcv))
    if (nsnd .gt. 0) allocate(this%buff_snd(this%lx, this%lx, this%lx, nsnd))
    if (ncrs .gt. 0) allocate(this%vcrs(this%lx, this%lx, this%lx, &
         amr_nchildren, ncrs))

  end subroutine amr_reconstruct_lx_buff_get

  !> Free buffer arrays
  subroutine amr_reconstruct_lx_buff_free(this)
    class(amr_reconstruct_lx_t), intent(inout) :: this

    if (allocated(this%vout)) deallocate(this%vout)
    if (allocated(this%vcrs)) deallocate(this%vcrs)
    if (allocated(this%buff_rcv)) deallocate(this%buff_rcv)
    if (allocated(this%buff_snd)) deallocate(this%buff_snd)

  end subroutine amr_reconstruct_lx_buff_free

  !> Perform refinement/coarsening on a single vector
  !! @param[inout] transfer  data transfer type
  !! @param[in]    nnew      new number of elements
  !! @param[in]    nref      Refinement mapping size
  !! @param[in]    ncrs      Coarsening mapping size
  !! @param[in]    rmap      Refinement mapping
  !! @param[in]    cmap      Coarsening mapping
  !! @param[inout] vec       vector for refinement/coarsening
  !! @param[inout] vec_d     device pointer to vector
  subroutine amr_reconstruct_lx_refine_coarsen(this, transfer, nnew, nref, &
       ncrs, rmap, cmap, vec, vec_d)
    class(amr_reconstruct_lx_t), intent(inout) :: this
    class(mesh_manager_transfer_t), intent(inout) :: transfer
    integer, intent(in) :: nnew, nref, ncrs
    integer, dimension(:, :),intent(in) :: rmap, cmap
    real(rp), dimension(:, :, :, :), allocatable, intent(inout) :: vec
    type(c_ptr), optional, intent(inout) :: vec_d
    integer :: il, jl, itmp
    integer, dimension(3) :: ch_pos

    associate(int => this%interpolate)
      ! check device pointer
      if (present(vec_d) .and. NEKO_BCKND_DEVICE .eq. 1) then
         call neko_error('AMR reconstruction; nothing done for device')
      end if

      ! CPU part
      allocate(this%vout(this%lx, this%lx, this%lx, nnew))
      call transfer%vector_constr(vec, this%vout, this%vcrs, this%buff_rcv, &
           this%buff_snd, this%lxyz)
      ! refinement
      if (nref .gt. 0) then
         do il = 1, nref
            ! get child position with respect to the parent
            ! (0...nchildren - 1)
            itmp = rmap(2, il)
            ch_pos(3) = itmp / 4 + 1 ! z position
            ch_pos(2) = mod(itmp / 2, 2) + 1 ! y position
            ch_pos(1) = mod(itmp, 2) + 1 ! x position
            call this%map_c2f(int%gdim, ch_pos, this%vout(:, :, :, rmap(1, il)))
         end do
      end if
      ! coarsening
      if (ncrs .gt. 0) then
         do il = 1, ncrs
            this%vout(:, :, :, cmap(1, il)) = 0.0_rp
            do jl = 1, amr_nchildren
               ! get child position with respect to the parent
               ! (0...amr_nchildren - 1)
               itmp = cmap(1 + jl, il)
               ch_pos(3) = itmp / 4 + 1 ! z position
               ch_pos(2) = mod(itmp / 2, 2) + 1 ! y position
               ch_pos(1) = mod(itmp, 2) + 1 ! x position
               call this%map_f2c(int%gdim, ch_pos, &
                       this%vcrs(:, :, :, jl, il), this%tmp(:, :, :, 3))
               ! accumulate result
               this%vout(:, :, :, cmap(1, il)) = &
                    this%vout(:, :, :, cmap(1, il)) + this%tmp(:, :, :, 3)
            end do
            ! reduce multiplicity
            this%vout(:, :, :, cmap(1, il)) = &
                 this%vout(:, :, :, cmap(1, il)) * int%el_mult(:, :, :)
         end do
      end if
      call move_alloc(this%vout, vec)

      ! check device pointer
      if (present(vec_d) .and. NEKO_BCKND_DEVICE .eq. 1) then
         call neko_error('AMR reconstruction; nothing done for device; 2')
      end if
    end associate

  end subroutine amr_reconstruct_lx_refine_coarsen

  !> @brief Map a single coarse element to a fine one
  !! @param[in]    gdim    geometrical dimension
  !! @param[in]    ch_pos  child position
  !! @param[inout] vcf     coarse (input) and fine (output) element vector
  subroutine amr_reconstruct_lx_map_c2f(this, gdim, ch_pos, vcf)
    class(amr_reconstruct_lx_t), intent(inout) :: this
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

  end subroutine amr_reconstruct_lx_map_c2f

  !> @brief Map a single fine element to a coarse one
  !! @param[in]    gdim    geometrical dimension
  !! @param[in]    ch_pos  child position
  !! @param[in]    vf      fine element vector
  !! @param[out]   vc      coarse element vector
  subroutine amr_reconstruct_lx_map_f2c(this, gdim, ch_pos, vf, vc)
    class(amr_reconstruct_lx_t), intent(inout) :: this
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

  end subroutine amr_reconstruct_lx_map_f2c

  !> Initialise type
  !! @param[in]  transfer     mesh manager data transfer type
  !! @param[in]  gdim         geometrical mesh dimension
  !! @param[in]  lx           polynomial order + 1
  subroutine amr_reconstruct_init(this, transfer, gdim, lx)
    class(amr_reconstruct_t), intent(inout) :: this
    class(mesh_manager_transfer_t), target, intent(in) :: transfer
    integer, intent(in) :: gdim, lx
    integer :: il

    call this%free()

    ! Primary number of point in 1D, 2D and 3D. Required by solvers and
    ! projections, as they just get a whole vector length at initialisation.
    ! I assume lx = ly = lz.
    this%lx = lx
    this%lxy = this%lx * lx
    this%lxyz = this%lxy * lx

    ! mesh manager data transfer
    this%transfer => transfer

    ! THIS ISN'T SMART, BUT FOR NOW I INITIALISE ALL
    ! AMR reconstruction operators for various lx
    this%int_lx = lx + 2 ! due to Schwarz
    allocate(this%recon(2: this%int_lx))
    do il = 2, this%int_lx
       call this%recon(il)%init(gdim, il)
    end do

!    if (NEKO_BCKND_DEVICE .eq. 1) then
!       ! should tmp be on device?
!    end if
!    call device_sync()

  end subroutine amr_reconstruct_init

  !> Free type
  subroutine amr_reconstruct_free(this)
    class(amr_reconstruct_t), intent(inout) :: this
    integer :: il

    this%transfer => NULL()

    if (allocated(this%recon)) then
       do il = 2, this%int_lx
          call this%recon(il)%free()
       end do
       deallocate(this%recon)
    end if

    this%lx = 0
    this%lxy = 0
    this%lxyz = 0
    this%int_lx = 0
    this%nold = 0
    this%nnew = 0
    this%nref = 0
    this%ncrs = 0
    this%nrcv = 0
    this%nsnd = 0
    this%ifchange = .false.

    if (allocated(this%rmap)) deallocate(this%rmap)
    if (allocated(this%cmap)) deallocate(this%cmap)

    !
    ! Cleanup the device (if present)
    !
    ! should some stuff be on device?

  end subroutine amr_reconstruct_free

  !> Get refinement/coarsening mapping
  subroutine amr_reconstruct_map_get(this)
    class(amr_reconstruct_t), intent(inout) :: this
    integer :: nchildren, il

    call this%transfer%vector_map(this%nold, this%nnew, this%nref, this%ncrs, &
         this%rmap, this%cmap, nchildren, this%ifchange, this%nrcv, this%nsnd)

    ! sanity check
    if (amr_nchildren .ne. nchildren) &
         call neko_error('Inconsistent children number for mesh manager and &
         &refinement routines.')

    ! THIS ISN'T SMART, BUT FOR NOW I ALLOCATE ALL
    ! get work space
    if (this%ifchange) then
       do il = 2, this%int_lx
          call this%recon(il)%buff_get(this%nrcv, this%nsnd, this%ncrs)
       end do
    end if

  end subroutine amr_reconstruct_map_get

  !> Free refinement/coarsening mapping
  subroutine amr_reconstruct_map_free(this)
    class(amr_reconstruct_t), intent(inout) :: this
    integer :: il

    call this%transfer%vector_map_free()

    this%nold = 0
    this%nnew = 0
    this%nref = 0
    this%ncrs = 0
    this%nrcv = 0
    this%nsnd = 0
    this%ifchange = .false.

    if (allocated(this%rmap)) deallocate(this%rmap)
    if (allocated(this%cmap)) deallocate(this%cmap)

    do il = 2, this%int_lx
       call this%recon(il)%buff_free()
    end do

  end subroutine amr_reconstruct_map_free

  !> Perform refinement/coarsening on a single vector
  !! @param[inout] vec    vector for refinement/coarsening
  !! @param[in]    lx     polynomial order + 1
  !! @param[inout] vec_d  device pointer to vector
  subroutine amr_reconstruct_refine_coarsen(this, vec, lx, vec_d)
    class(amr_reconstruct_t), intent(inout) :: this
    real(rp), dimension(:, :, :, :), allocatable, intent(inout) :: vec
    integer, intent(in) :: lx
    type(c_ptr), optional, intent(inout) :: vec_d
    integer :: il, jl, itmp
    integer, dimension(3) :: ch_pos

    if (lx .lt. 2 .or. lx .gt. this%int_lx) &
         call neko_error('Polynomial order outside the assumed range')

    if (this%ifchange) then
       if (present(vec_d) .and. NEKO_BCKND_DEVICE .eq. 1) then
          call this%recon(lx)%refine_coarsen(this%transfer, this%nnew, &
               this%nref, this%ncrs, this%rmap, this%cmap, vec, vec_d)
       else
          call this%recon(lx)%refine_coarsen(this%transfer, this%nnew, &
               this%nref, this%ncrs, this%rmap, this%cmap, vec)
       end if
    end if

  end subroutine amr_reconstruct_refine_coarsen

end module amr_reconstruct
