! Copyright (c) 2020-2021, The Neko Authors
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
!> Mixed Dirichlet-Neumann axis aligned symmetry plane
module symmetry
  use device_symmetry
  use neko_config
  use num_types
  use dirichlet
  use device
  use coefs
  use math
  use utils
  use stack
  use, intrinsic :: iso_c_binding
  implicit none
  private

  !> Mixed Dirichlet-Neumann symmetry plane condition
  type, public, extends(dirichlet_t) :: symmetry_t
     integer, allocatable :: xaxis_msk(:)
     integer, allocatable :: zaxis_msk(:)
     integer, allocatable :: yaxis_msk(:)
     type(c_ptr) :: xaxis_msk_d = C_NULL_PTR
     type(c_ptr) :: yaxis_msk_d = C_NULL_PTR
     type(c_ptr) :: zaxis_msk_d = C_NULL_PTR
   contains
     procedure, pass(this) :: init_msk => symmetry_init_msk
     procedure, pass(this) :: apply_scalar => symmetry_apply_scalar
     procedure, pass(this) :: apply_vector => symmetry_apply_vector
     procedure, pass(this) :: apply_scalar_dev => symmetry_apply_scalar_dev
     procedure, pass(this) :: apply_vector_dev => symmetry_apply_vector_dev
     final :: symmetry_free
  end type symmetry_t

contains

  !> Initialize symmetry mask for each axis
  subroutine symmetry_init_msk(this, c)
    class(symmetry_t), intent(inout) :: this
    type(coef_t), intent(in) :: c
    type(stack_i4_t), target :: xmsk, ymsk, zmsk
    integer :: i, m, j, k, l, idx(4), facet, ntype, msk_size
    integer, pointer :: sp(:)        
    real(kind=rp) :: sx,sy,sz
    real(kind=rp), parameter :: TOL = 1d-3
    
    call symmetry_free(this)

    call xmsk%init()
    call ymsk%init()
    call zmsk%init()
    
    associate(nx => c%nx, ny => c%ny, nz => c%nz)
      m = this%msk(0)
      do i = 1, m
         k = this%msk(i)
         facet = this%facet(i)
         idx = nonlinear_index(k, c%Xh%lx, c%Xh%lx, c%Xh%lx)
         sx = 0d0
         sy = 0d0
         sz = 0d0
         select case (facet)               
         case(1,2)
            do l = 2, c%Xh%lx - 1
               do j = 2, c%Xh%lx -1
                  sx = sx + abs(abs(nx(l, j, facet, idx(4))) - 1d0)
                  sy = sy + abs(abs(ny(l, j, facet, idx(4))) - 1d0)
                  sz = sz + abs(abs(nz(l, j, facet, idx(4))) - 1d0)
               end do
            end do
         case(3,4)
            do l = 2, c%Xh%lx - 1
               do j = 2, c%Xh%lx - 1
                  sx = sx + abs(abs(nx(l, j, facet, idx(4))) - 1d0)
                  sy = sy + abs(abs(ny(l, j, facet, idx(4))) - 1d0)
                  sz = sz + abs(abs(nz(l, j, facet, idx(4))) - 1d0)
               end do
            end do
         case(5,6)
            do l = 2, c%Xh%lx - 1
               do j = 2, c%Xh%lx - 1
                  sx = sx + abs(abs(nx(l, j, facet, idx(4))) - 1d0)
                  sy = sy + abs(abs(ny(l, j, facet, idx(4))) - 1d0)
                  sz = sz + abs(abs(nz(l, j, facet, idx(4))) - 1d0)
               end do
            end do               
         end select
         sx = sx / (c%Xh%lx - 2)**2
         sy = sy / (c%Xh%lx - 2)**2
         sz = sz / (c%Xh%lx - 2)**2

         ntype = 0
         if (sx .lt. TOL) then
            ntype = iand(ntype, 1)
            call xmsk%push(k)
         end if

         if (sy .lt. TOL) then
            ntype = iand(ntype, 2)
            call ymsk%push(k)
         end if

         if (sz .lt. TOL) then
            ntype = iand(ntype, 4)
            call zmsk%push(k)
         end if

      end do
    end associate

    !> @note This is to prevent runtime issues with Cray ftn, gfortran and
    !! msk:size() in the allocate call
    msk_size = xmsk%size()
    if (msk_size .gt. 0) then
       allocate(this%xaxis_msk(0:msk_size))
       this%xaxis_msk(0) = msk_size
       sp => xmsk%array()
       do i = 1, msk_size
          this%xaxis_msk(i) = sp(i)
       end do
       if ((NEKO_BCKND_HIP .eq. 1) .or. (NEKO_BCKND_CUDA .eq. 1) .or. &
            (NEKO_BCKND_OPENCL .eq. 1)) then
          call device_map(this%xaxis_msk, this%xaxis_msk_d, msk_size + 1)
          call device_memcpy(this%xaxis_msk, this%xaxis_msk_d, &
               msk_size + 1, HOST_TO_DEVICE)
       end if
    else
       allocate(this%xaxis_msk(0:1))
       this%xaxis_msk(0) = 0
       if ((NEKO_BCKND_HIP .eq. 1) .or. (NEKO_BCKND_CUDA .eq. 1) .or. &
            (NEKO_BCKND_OPENCL .eq. 1)) then
          call device_map(this%xaxis_msk, this%xaxis_msk_d, 2)
          call device_memcpy(this%xaxis_msk, this%xaxis_msk_d, &
               2, HOST_TO_DEVICE)
       end if
    end if

    msk_size = ymsk%size()
    if (msk_size .gt. 0) then
       allocate(this%yaxis_msk(0:msk_size))
       this%yaxis_msk(0) = msk_size
       sp => ymsk%array()
       do i = 1, msk_size
          this%yaxis_msk(i) = sp(i)
       end do
       if ((NEKO_BCKND_HIP .eq. 1) .or. (NEKO_BCKND_CUDA .eq. 1) .or. &
            (NEKO_BCKND_OPENCL .eq. 1)) then
          call device_map(this%yaxis_msk, this%yaxis_msk_d, msk_size + 1)
          call device_memcpy(this%yaxis_msk, this%yaxis_msk_d, &
               msk_size + 1, HOST_TO_DEVICE)
       end if
    else
       allocate(this%yaxis_msk(0:1))
       this%yaxis_msk(0) = 0
       if ((NEKO_BCKND_HIP .eq. 1) .or. (NEKO_BCKND_CUDA .eq. 1) .or. &
            (NEKO_BCKND_OPENCL .eq. 1)) then
          call device_map(this%yaxis_msk, this%yaxis_msk_d, 2)
          call device_memcpy(this%yaxis_msk, this%yaxis_msk_d, &
               2, HOST_TO_DEVICE)
       end if
    end if

    msk_size = zmsk%size()
    if (msk_size .gt. 0) then
       allocate(this%zaxis_msk(0:msk_size))
       this%zaxis_msk(0) = msk_size
       sp => zmsk%array()
       do i = 1, msk_size
          this%zaxis_msk(i) = sp(i)
       end do
       if ((NEKO_BCKND_HIP .eq. 1) .or. (NEKO_BCKND_CUDA .eq. 1) .or. &
            (NEKO_BCKND_OPENCL .eq. 1)) then
          call device_map(this%zaxis_msk, this%zaxis_msk_d, msk_size + 1)
          call device_memcpy(this%zaxis_msk, this%zaxis_msk_d, &
               msk_size + 1, HOST_TO_DEVICE)
       end if
    else
       allocate(this%zaxis_msk(0:1))
       this%zaxis_msk(0) = 0
       if ((NEKO_BCKND_HIP .eq. 1) .or. (NEKO_BCKND_CUDA .eq. 1) .or. &
            (NEKO_BCKND_OPENCL .eq. 1)) then
          call device_map(this%zaxis_msk, this%zaxis_msk_d, 2)
          call device_memcpy(this%zaxis_msk, this%zaxis_msk_d, &
               2, HOST_TO_DEVICE)
       end if
    end if    

    nullify(sp)
    call xmsk%free()
    call ymsk%free()
    call zmsk%free()
    
  end subroutine symmetry_init_msk
  
  subroutine symmetry_free(this)
    type(symmetry_t), intent(inout) :: this

    if (allocated(this%xaxis_msk)) then
       deallocate(this%xaxis_msk)
    end if
    
    if (allocated(this%yaxis_msk)) then
       deallocate(this%yaxis_msk)
    end if

    if (allocated(this%zaxis_msk)) then
       deallocate(this%zaxis_msk)
    end if

    if (c_associated(this%xaxis_msk_d)) then
       call device_free(this%xaxis_msk_d)
       this%xaxis_msk_d = C_NULL_PTR
    end if

    if (c_associated(this%yaxis_msk_d)) then
       call device_free(this%yaxis_msk_d)
       this%yaxis_msk_d = C_NULL_PTR
    end if

    if (c_associated(this%zaxis_msk_d)) then
       call device_free(this%zaxis_msk_d)
       this%zaxis_msk_d = C_NULL_PTR
    end if

  end subroutine symmetry_free
  
  !> No-op scalar apply
  subroutine symmetry_apply_scalar(this, x, n)
    class(symmetry_t), intent(inout) :: this
    integer, intent(in) :: n
    real(kind=rp), intent(inout), dimension(n) :: x
  end subroutine symmetry_apply_scalar

  !> Apply symmetry conditions (axis aligned)
  subroutine symmetry_apply_vector(this, x, y, z, n)
    class(symmetry_t), intent(inout) :: this
    integer, intent(in) :: n
    real(kind=rp), intent(inout),  dimension(n) :: x
    real(kind=rp), intent(inout),  dimension(n) :: y
    real(kind=rp), intent(inout),  dimension(n) :: z
    integer :: i, m, k

    m = this%xaxis_msk(0)
    do i = 1, m
       k = this%xaxis_msk(i)
       x(k) = 0d0
    end do

    m = this%yaxis_msk(0)
    do i = 1, m
       k = this%yaxis_msk(i)
       y(k) = 0d0
    end do

    m = this%zaxis_msk(0)
    do i = 1, m
       k = this%zaxis_msk(i)
       z(k) = 0d0
    end do
    
  end subroutine symmetry_apply_vector

  !> No-op scalar apply (device version)
  subroutine symmetry_apply_scalar_dev(this, x_d)
    class(symmetry_t), intent(inout), target :: this
    type(c_ptr) :: x_d
  end subroutine symmetry_apply_scalar_dev

  !> Apply symmetry conditions (axis aligned) (device version)
  subroutine symmetry_apply_vector_dev(this, x_d, y_d, z_d)
    class(symmetry_t), intent(inout), target :: this
    type(c_ptr) :: x_d
    type(c_ptr) :: y_d
    type(c_ptr) :: z_d

    call device_symmetry_apply_vector(this%xaxis_msk_d, this%yaxis_msk_d, &
                                      this%zaxis_msk_d, x_d, y_d, z_d, &
                                      this%xaxis_msk(0), &
                                      this%yaxis_msk(0), &
                                      this%zaxis_msk(0))


  end subroutine symmetry_apply_vector_dev
      
end module symmetry
