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
  use neko_config
  use num_types, only : i8, rp
  use utils, only : neko_error, neko_warning
  use speclib, only : igllm
  use space, only : space_t, GLL
  use device
  use mxm_wrapper, only: mxm
  use math, only : invcol1
  use mesh_manager_transfer, only : mesh_manager_transfer_t
  use, intrinsic :: iso_c_binding

  implicit none
  private

  ! number of children; 3D only
  integer, parameter :: amr_nchildren = 8

  !> Type for vector/field reconstruction
  type, public :: amr_reconstruct_t
     !> pointer mesh manager data transfer
     class(mesh_manager_transfer_t), pointer :: transfer
     !> Function space \f$ X_h \f$
     type(space_t) :: Xh
     !> Interpolation operators; coarse to fine
     real(rp), allocatable, dimension(:, :, :) :: x_cr2fn, x_cr2fnT
     real(rp), allocatable, dimension(:, :, :) :: y_cr2fn, y_cr2fnT
     real(rp), allocatable, dimension(:, :, :) :: z_cr2fn, z_cr2fnT
     !> Interpolation operators; fine to coarse
     real(rp), allocatable, dimension(:, :, :) :: x_fn2cr, x_fn2crT
     real(rp), allocatable, dimension(:, :, :) :: y_fn2cr, y_fn2crT
     real(rp), allocatable, dimension(:, :, :) :: z_fn2cr, z_fn2crT
     ! used for coarsening after children face summation
     ! I assume lx=ly=lz, but in general there should be 3 face and edge arrays
     !> Point multiplicity; element
     real(rp), allocatable, dimension(:, :, :) :: el_mult
     !> Point multiplicity; face
     real(rp), allocatable, dimension(:, :) :: fc_mult
     !> Point multiplicity; edge
     real(rp), allocatable, dimension(:) :: ed_mult

     !> Work space for single element refinement
     real(rp), allocatable, dimension(:, :, :, :) :: tmp

     ! Arrays sizes
     !> Old element number
     integer :: nold
     !> New element number
     integer :: nnew
     !> Refinement mapping size
     integer :: nref
     !> Coarsening mapping size
     integer :: ncrs
     !> Refinement mapping (element and child position)
     integer, dimension(:, :), allocatable :: rmap
     !> Coarsening mapping (element position)
     integer, dimension(:), allocatable :: cmap
     !> Output vector for refinement
     real(rp), allocatable, dimension(:, :, :, :) :: vout
     !> Vector with additional data for coarsening
     real(rp), allocatable, dimension(:, :, :, :, :) :: vcrs

     !
     ! Device pointers (if present)
     !
     type(c_ptr) :: x_cr2fn_d = C_NULL_PTR
     type(c_ptr) :: x_cr2fnT_d = C_NULL_PTR
     type(c_ptr) :: y_cr2fn_d = C_NULL_PTR
     type(c_ptr) :: y_cr2fnT_d = C_NULL_PTR
     type(c_ptr) :: z_cr2fn_d = C_NULL_PTR
     type(c_ptr) :: z_cr2fnT_d = C_NULL_PTR
     type(c_ptr) :: x_fn2cr_d = C_NULL_PTR
     type(c_ptr) :: x_fn2crT_d = C_NULL_PTR
     type(c_ptr) :: y_fn2cr_d = C_NULL_PTR
     type(c_ptr) :: y_fn2crT_d = C_NULL_PTR
     type(c_ptr) :: z_fn2cr_d = C_NULL_PTR
     type(c_ptr) :: z_fn2crT_d = C_NULL_PTR
     type(c_ptr) :: el_mult_d = C_NULL_PTR
     type(c_ptr) :: fc_mult_d = C_NULL_PTR
     type(c_ptr) :: ed_mult_d = C_NULL_PTR
     ! should tmp be on device?
   contains
     !> Initialise type
     procedure, pass(this) :: init => amr_reconstruct_init
     !> Free type
     procedure, pass(this) :: free => amr_reconstruct_free
     !> Get refinement/coarsening mapping
     procedure, pass(this) :: map_get => amr_reconstruct_map_get
     !> free refinement/coarsening mapping
     procedure, pass(this) :: map_free => amr_reconstruct_map_free
     !> Single field coarsening operation
     procedure, pass(this) :: coarsen_single => amr_reconstruct_coarsen_single
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
    integer :: il, jl, kl, nt2
    real(rp), allocatable, dimension(:) :: tmpl

    call this%free()

    ! mesh manager data transfer
    this%transfer => transfer

    ! functional space
    if (gdim .ne. 3) &
         call neko_error('AMR reconstruction supports 3D domain only')
    call this%Xh%init(GLL, lx, lx, lx)

    associate(Xh => this%Xh)

      ! work array
      allocate(tmpl(max(Xh%lx, Xh%ly, Xh%lz)))

      ! interpolation operators
      ! x-direction
      allocate(this%x_cr2fn(Xh%lx, Xh%lx, 2), source = 0.0_rp)
      allocate(this%x_cr2fnT, source = this%x_cr2fn)
      allocate(this%x_fn2cr, source = this%x_cr2fn)
      allocate(this%x_fn2crT, source = this%x_cr2fn)

      ! coarse -> fine
      ! negative
      do concurrent (il = 1: Xh%lx)
         tmpl(il) = 0.5_rp * (Xh%zg(il,1) - 1.0_rp)
      end do
      call igllm(this%x_cr2fn, this%x_cr2fnT, Xh%zg(:, 1), tmpl, Xh%lx, Xh%lx, &
           Xh%lx, Xh%lx)
      ! positive; we use symmetry
      do concurrent (jl = 1: Xh%lx, il = 1: Xh%lx)
         this%x_cr2fn(Xh%lx-il+1, Xh%lx-jl+1, 2) = this%x_cr2fn(il, jl, 1)
         this%x_cr2fnT(Xh%lx-il+1, Xh%lx-jl+1, 2) = this%x_cr2fnT(il, jl, 1)
      end do
      ! fine -> coarse
      ! negative
      nt2 = Xh%lx/2 + mod(Xh%lx, 2)
      do concurrent (il=1: nt2)
         tmpl(il) = 2.0_rp * Xh%zg(il, 1) + 1.0_rp
      end do
      call igllm(this%x_fn2cr, this%x_fn2crT, Xh%zg(:, 1), tmpl, Xh%lx, &
           nt2, Xh%lx, Xh%lx)
      ! positive; we use symmetry
      do concurrent (jl = 1: Xh%lx, il = 1: nt2)
         this%x_fn2cr(Xh%lx-il+1, Xh%lx-jl+1, 2) = this%x_fn2cr(il, jl, 1)
         this%x_fn2crT(Xh%lx-il+1, Xh%lx-jl+1, 2) = this%x_fn2crT(il, jl, 1)
      end do

      ! y-direction
      allocate(this%y_cr2fn(Xh%ly, Xh%ly, 2), source = 0.0_rp)
      allocate(this%y_cr2fnT, source = this%y_cr2fn)
      allocate(this%y_fn2cr, source = this%y_cr2fn)
      allocate(this%y_fn2crT, source = this%y_cr2fn)

      ! coarse -> fine
      ! negative
      do concurrent (il = 1: Xh%ly)
         tmpl(il) = 0.5_rp * (Xh%zg(il, 2) - 1.0_rp)
      end do
      call igllm(this%y_cr2fn, this%y_cr2fnT, Xh%zg(:, 2), tmpl, Xh%ly, Xh%ly, &
           Xh%ly, Xh%ly)
      ! positive; we use symmetry
      do concurrent (jl = 1: Xh%ly, il = 1: Xh%ly)
         this%y_cr2fn(Xh%ly-il+1, Xh%ly-jl+1, 2) = this%y_cr2fn(il, jl, 1)
         this%y_cr2fnT(Xh%ly-il+1, Xh%ly-jl+1, 2) = this%y_cr2fnT(il, jl, 1)
      end do
      ! fine -> coarse
      ! negative
      nt2 = Xh%ly/2 + mod(Xh%ly, 2)
      do concurrent (il = 1: nt2)
         tmpl(il) = 2.0_rp * Xh%zg(il, 2) + 1.0_rp
      end do
      call igllm(this%y_fn2cr, this%y_fn2crT, Xh%zg(:, 2), tmpl, Xh%ly, nt2, &
           Xh%ly, Xh%ly)
      ! positive; we use symmetry
      do concurrent (jl = 1: Xh%ly, il = 1: nt2)
         this%y_fn2cr(Xh%ly-il+1, Xh%ly-jl+1, 2) = this%y_fn2cr(il, jl, 1)
         this%y_fn2crT(Xh%ly-il+1, Xh%ly-jl+1, 2) = this%y_fn2crT(il, jl, 1)
      end do

      ! z-direction
      if (gdim == 3) then
         allocate(this%z_cr2fn(Xh%lz, Xh%lz, 2), source = 0.0_rp)
         allocate(this%z_cr2fnT, source = this%z_cr2fn)
         allocate(this%z_fn2cr, source = this%z_cr2fn)
         allocate(this%z_fn2crT, source = this%z_cr2fn)

         ! coarse -> fine
         ! negative
         do concurrent (il = 1: Xh%lz)
            tmpl(il) = 0.5_rp * (Xh%zg(il, 3) - 1.0_rp)
         end do
         call igllm(this%z_cr2fn, this%z_cr2fnT, Xh%zg(:, 3), &
              &tmpl, Xh%lz, Xh%lz, Xh%lz, Xh%lz)
         ! positive; we use symmetry
         do concurrent (jl = 1: Xh%lz, il = 1: Xh%lz)
            this%z_cr2fn(Xh%lz-il+1, Xh%lz-jl+1, 2) = this%z_cr2fn(il, jl, 1)
            this%z_cr2fnT(Xh%lz-il+1, Xh%lz-jl+1, 2) = this%z_cr2fnT(il, jl, 1)
         end do
         ! fine -> coarse
         ! negative
         nt2 = Xh%lz/2 + mod(Xh%lz, 2)
         do concurrent (il = 1: nt2)
            tmpl(il) = 2.0_rp * Xh%zg(il, 3) + 1.0_rp
         end do
         call igllm(this%z_fn2cr, this%z_fn2crT, Xh%zg(:, 3), tmpl, Xh%lz, &
              nt2, Xh%lz, Xh%lz)
         ! positive; we use symmetry
         do concurrent (jl = 1: Xh%lz, il = 1: nt2)
            this%z_fn2cr(Xh%lz-il+1, Xh%lz-jl+1, 2) = this%z_fn2cr(il, jl, 1)
            this%z_fn2crT(Xh%lz-il+1, Xh%lz-jl+1, 2) = this%z_fn2crT(il, jl, 1)
         end do
      end if

      ! multiplicity arrays
      if ((Xh%lx /= Xh%ly) .or. ((Xh%lz /= 1) .and. (Xh%lz /= Xh%lx))) &
           call neko_error('AMR does not support lx /= ly /= lz (for lz /= 1)')
      allocate(this%el_mult(Xh%lx, Xh%ly, Xh%lz), source = 0.0_rp)
      allocate(this%fc_mult(Xh%lx, Xh%lx), source = 0.0_rp)
      allocate(this%ed_mult(Xh%lx), source = 0.0_rp)
      ! X
      if (mod(Xh%lx, 2) == 1) then
         il = Xh%lx/2 + 1
         do concurrent (kl = 1: Xh%lz, jl = 1: Xh%ly)
            this%el_mult(il, jl, kl) = this%el_mult(il, jl, kl) + 1.0_rp
         end do
      end if
      ! Y
      if (mod(Xh%ly, 2) == 1) then
         jl = Xh%ly/2 + 1
         do concurrent (kl = 1: Xh%lz, il = 1: Xh%lx)
            this%el_mult(il, jl, kl) = this%el_mult(il, jl, kl) + 1.0_rp
         end do
         if (mod(Xh%lx, 2) == 1) then
            il = Xh%lx/2 + 1
            do concurrent (kl = 1: Xh%lz)
               this%el_mult(il, jl, kl) = this%el_mult(il, jl, kl) + 1.0_rp
            end do
         end if
      end if
      ! Z
      if (gdim == 3) then
         if (mod(Xh%lz, 2) == 1) then
            kl = Xh%lz/2 + 1
            do concurrent (jl = 1: Xh%ly, il = 1: Xh%lx)
               this%el_mult(il, jl, kl) = this%el_mult(il, jl, kl) + 1.0_rp
            end do
            if (mod(Xh%lx, 2) == 1) then
               il = Xh%lx/2 + 1
               do concurrent (jl = 1: Xh%ly)
                  this%el_mult(il, jl, kl) = this%el_mult(il, jl, kl) + 1.0_rp
               end do
            end if
            if (mod(Xh%ly, 2) == 1) then
               jl = Xh%ly/2 + 1
               do concurrent (il = 1: Xh%lx)
                  this%el_mult(il, jl, kl) = this%el_mult(il, jl, kl) + 1.0_rp
               end do
            end if
            if ((mod(Xh%lx, 2) == 1).and.(mod(Xh%ly,2 ) == 1)) then
               il = Xh%lx/2 + 1
               jl = Xh%ly/2 + 1
               this%el_mult(il, jl, kl) = this%el_mult(il, jl, kl) + 1.0_rp
            end if
         end if
      end if
      ! calculate inverse
      nt2 = Xh%lx * Xh%ly * Xh%lz
      call invcol1(this%el_mult, nt2)

      ! to get proper J^-1 on faces and edges for fast diagonalisation method
      ! I assume here LX1 = LY1 = LZ1, so only one array is needed
      !call ftovecl(this%fc_mult,this%el_mult,1,Xh%lx,Xh%ly,Xh%lz)
      !call etovec(this%ed_mult,1,this%el_mult,Xh%lx,Xh%ly,Xh%lz)
      ! THIS SHOLD BE CHECKED, BUT SHOULD BE FINE
      this%fc_mult(:,:) = this%el_mult(:,:, 1)
      this%ed_mult(:) = this%el_mult(:, 1, 1)

      ! work space
      allocate(this%tmp(Xh%lx, Xh%ly, Xh%lz, 3))

      if (NEKO_BCKND_DEVICE .eq. 1) then
         nt2 = Xh%lxy * 2
         call device_map(this%x_cr2fn, this%x_cr2fn_d, nt2)
         call device_map(this%x_cr2fnT, this%x_cr2fnT_d, nt2)
         call device_map(this%x_fn2cr, this%x_fn2cr_d, nt2)
         call device_map(this%x_fn2crT, this%x_fn2crT_d, nt2)
         call device_map(this%y_cr2fn, this%y_cr2fn_d, nt2)
         call device_map(this%y_cr2fnT, this%y_cr2fnT_d, nt2)
         call device_map(this%y_fn2cr, this%y_fn2cr_d, nt2)
         call device_map(this%y_fn2crT, this%y_fn2crT_d, nt2)
         call device_map(this%z_cr2fn, this%x_cr2fn_d, nt2)
         call device_map(this%z_cr2fnT, this%z_cr2fnT_d, nt2)
         call device_map(this%z_fn2cr, this%z_fn2cr_d, nt2)
         call device_map(this%z_fn2crT, this%z_fn2crT_d, nt2)
         call device_map(this%el_mult, this%el_mult_d, Xh%lxyz)
         call device_map(this%fc_mult, this%fc_mult_d, Xh%lxy)
         call device_map(this%ed_mult, this%ed_mult_d, Xh%lx)
         ! should tmp be on device?
      end if

    end associate

    call device_sync()

    deallocate(tmpl)

  end subroutine amr_reconstruct_init

  !> Free type
  subroutine amr_reconstruct_free(this)
    class(amr_reconstruct_t), intent(inout) :: this

    this%transfer => NULL()

    call this%Xh%free()

    this%nold = 0
    this%nnew = 0
    this%nref = 0
    this%ncrs = 0

    if (allocated(this%x_cr2fn)) deallocate(this%x_cr2fn)
    if (allocated(this%x_cr2fnT)) deallocate(this%x_cr2fnT)
    if (allocated(this%x_fn2cr)) deallocate(this%x_fn2cr)
    if (allocated(this%x_fn2crT)) deallocate(this%x_fn2crT)

    if (allocated(this%y_cr2fn)) deallocate(this%y_cr2fn)
    if (allocated(this%y_cr2fnT)) deallocate(this%y_cr2fnT)
    if (allocated(this%y_fn2cr)) deallocate(this%y_fn2cr)
    if (allocated(this%y_fn2crT)) deallocate(this%y_fn2crT)

    if (allocated(this%z_cr2fn)) deallocate(this%z_cr2fn)
    if (allocated(this%z_cr2fnT)) deallocate(this%z_cr2fnT)
    if (allocated(this%z_fn2cr)) deallocate(this%z_fn2cr)
    if (allocated(this%z_fn2crT)) deallocate(this%z_fn2crT)

    if (allocated(this%el_mult)) deallocate(this%el_mult)
    if (allocated(this%fc_mult)) deallocate(this%fc_mult)
    if (allocated(this%ed_mult)) deallocate(this%ed_mult)

    if (allocated(this%tmp)) deallocate(this%tmp)

    if (allocated(this%rmap)) deallocate(this%rmap)
    if (allocated(this%cmap)) deallocate(this%cmap)
    if (allocated(this%vout)) deallocate(this%vout)
    if (allocated(this%vcrs)) deallocate(this%vcrs)

    !
    ! Cleanup the device (if present)
    !
    if (c_associated(this%x_cr2fn_d)) call device_free(this%x_cr2fn_d)
    if (c_associated(this%x_cr2fnT_d)) call device_free(this%x_cr2fnT_d)
    if (c_associated(this%y_cr2fn_d)) call device_free(this%y_cr2fn_d)
    if (c_associated(this%y_cr2fnT_d)) call device_free(this%y_cr2fnT_d)
    if (c_associated(this%z_cr2fn_d)) call device_free(this%z_cr2fn_d)
    if (c_associated(this%z_cr2fnT_d)) call device_free(this%z_cr2fnT_d)
    if (c_associated(this%x_fn2cr_d)) call device_free(this%x_fn2cr_d)
    if (c_associated(this%x_fn2crT_d)) call device_free(this%x_fn2crT_d)
    if (c_associated(this%y_fn2cr_d)) call device_free(this%y_fn2cr_d)
    if (c_associated(this%y_fn2crT_d)) call device_free(this%y_fn2crT_d)
    if (c_associated(this%z_fn2cr_d)) call device_free(this%z_fn2cr_d)
    if (c_associated(this%z_fn2crT_d)) call device_free(this%z_fn2crT_d)
    if (c_associated(this%el_mult_d)) call device_free(this%el_mult_d)
    if (c_associated(this%fc_mult_d)) call device_free(this%fc_mult_d)
    if (c_associated(this%ed_mult_d)) call device_free(this%ed_mult_d)
    ! should tmp be on device?

  end subroutine amr_reconstruct_free

  !> Get refinement/coarsening mapping
  subroutine amr_reconstruct_map_get(this)
    class(amr_reconstruct_t), intent(inout) :: this

    call this%transfer%vector_map(this%nold, this%nnew, this%nref, this%ncrs, &
         this%rmap, this%cmap)

    ! get coarsening vector
    allocate(this%vcrs(this%Xh%lx, this%Xh%ly, this%Xh%lz, amr_nchildren, &
         this%ncrs))
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

  end subroutine amr_reconstruct_map_free

  !> @brief Perform a single field coarsening operation
  !! @param[inout] vcf     coarsened vector
  subroutine amr_reconstruct_coarsen_single(this, vfc)
    class(amr_reconstruct_t), intent(inout) :: this
    real(rp), dimension(:,:,:,:), intent(inout) :: vfc
    integer :: il, jl, itmp
    integer, dimension(3) :: ch_pos

!!$    ! JUST A PLACEHOLDER FOR NOW
!!$    ! I assume el_lst(1) gives position of the coarse block
!!$    ! and final ch_pos() = 1,1,1
!!$    ! loop over coarsened elements
!!$    do il= 1, this%transfer%coarsen_nr
!!$       ! loop over all the children
!!$       do jl= 1, amr_nchildren
!!$          ! get child position
!!$          itmp = this%transfer%coarsen(2, jl, il) ! new position in the array
!!$          ch_pos(3) = (jl-1)/4 + 1 ! z position
!!$          ch_pos(2) = mod((jl - 1)/2, 2) + 1 ! y position
!!$          ch_pos(1) = mod(jl - 1, 2) +1 ! x position
!!$          ! coarsen; 3D only
!!$          call this%map_f2c(3, ch_pos, vfc(:,:,:,itmp), this%tmp(:,:,:,3))
!!$          ! sum contributions
!!$          if (jl == 1) then
!!$             vfc(:,:,:,itmp) = this%tmp(:,:,:,3)
!!$          else
!!$             vfc(:,:,:,itmp) =  vfc(:,:,:,itmp) + this%tmp(:,:,:,3)
!!$          end if
!!$       end do
!!$       vfc(:,:,:,itmp) =  vfc(:,:,:,itmp)*this%el_mult(:,:,:)
!!$    end do

  end subroutine amr_reconstruct_coarsen_single

  !> @brief Map a single coarse element to a fine one
  !! @param[in]    gdim    geometrical dimension
  !! @param[in]    ch_pos  child position
  !! @param[in]    vc      coarse element vector
  !! @param[out]   vf      fine element vector
  subroutine amr_reconstruct_map_c2f(this, gdim, ch_pos, vc, vf)
    class(amr_reconstruct_t), intent(inout) :: this
    integer, intent(in) :: gdim
    integer, dimension(3), intent(in) :: ch_pos
    real(rp), dimension(:,:,:), intent(in) :: vc
    real(rp), dimension(:,:,:), intent(out) :: vf
    integer :: iz

    if (gdim == 3) then ! 3D
       call mxm(this%x_cr2fn(:,:, ch_pos(1)), this%Xh%lx, vc, &
            & this%Xh%lx, this%tmp(:,:,:, 1), this%Xh%lyz)
       do iz = 1, this%Xh%lz
          call mxm(this%tmp(:,:, iz, 1), this%Xh%lx,&
               & this%y_cr2fnT(:,:, ch_pos(2)),&
               & this%Xh%ly, this%tmp(:,:, iz, 2), this%Xh%ly)
       end do
       call mxm(this%tmp(:,:,:, 2), this%Xh%lxy,&
            & this%z_cr2fnT(:,:, ch_pos(3)),&
            & this%Xh%lz, vf, this%Xh%lz)
    else ! 2D
       call mxm(this%x_cr2fn(:,:, ch_pos(1)), this%Xh%lx, vc, this%Xh%lz,&
            & this%tmp(:,:,1,1), this%Xh%lyz)
       call mxm(this%tmp(:,:,1,1), this%Xh%lx, this%y_cr2fnT(:,:, ch_pos(2)),&
            & this%Xh%ly, vf, this%Xh%ly)
    end if

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

    if (gdim == 3) then ! 3D
       call mxm(this%x_fn2cr(:,:, ch_pos(1)), this%Xh%lx, vf, &
            & this%Xh%lx, this%tmp(:,:,:, 1), this%Xh%lyz)
       do iz = 1, this%Xh%lz
          call mxm(this%tmp(:,:, iz, 1), this%Xh%lx,&
               & this%y_fn2crT(:,:, ch_pos(2)),&
               & this%Xh%ly, this%tmp(:,:, iz, 2), this%Xh%ly)
       end do
       call mxm(this%tmp(:,:,:, 2), this%Xh%lxy,&
            & this%z_fn2crT(:,:, ch_pos(3)),&
            & this%Xh%lz, vc, this%Xh%lz)
    else ! 2D
       call mxm(this%x_fn2cr(:,:, ch_pos(1)), this%Xh%lx, vf, this%Xh%lz,&
            & this%tmp(:,:,1,1), this%Xh%lyz)
       call mxm(this%tmp(:,:,1,1), this%Xh%lx, this%y_fn2crT(:,:, ch_pos(2)),&
            & this%Xh%ly, vc, this%Xh%ly)
    end if

  end subroutine amr_reconstruct_map_f2c

end module amr_reconstruct
