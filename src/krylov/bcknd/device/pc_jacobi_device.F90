! Copyright (c) 2021-2022, The Neko Authors
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
!> Jacobi preconditioner accelerator backend
module device_jacobi
  use precon
  use coefs
  use dofmap
  use num_types
  use device_math
  use device
  use gather_scatter
  implicit none
  private

  !> Defines a jacobi preconditioner
  type, public, extends(pc_t) :: device_jacobi_t
     real(kind=rp), allocatable :: d(:,:,:,:)
     type(c_ptr) :: d_d = C_NULL_PTR
     type(gs_t), pointer :: gs_h
     type(dofmap_t), pointer :: dof
     type(coef_t), pointer :: coef
   contains
     procedure, pass(this) :: init => device_jacobi_init
     procedure, pass(this) :: free => device_jacobi_free
     procedure, pass(this) :: solve => device_jacobi_solve
     procedure, pass(this) :: update => device_jacobi_update
  end type device_jacobi_t

  interface
     subroutine hip_jacobi_update(d_d, dxt_d, dyt_d, dzt_d, &
          G11_d, G22_d, G33_d, G12_d, G13_d, G23_d, nelv, lx) &
          bind(c, name='hip_jacobi_update')
       use, intrinsic :: iso_c_binding
       type(c_ptr), value :: d_d, dxt_d, dyt_d, dzt_d
       type(c_ptr), value :: G11_d, G22_d, G33_d, G12_d, G13_d, G23_d
       integer(c_int) :: nelv, lx
     end subroutine hip_jacobi_update
  end interface

  interface
     subroutine cuda_jacobi_update(d_d, dxt_d, dyt_d, dzt_d, &
          G11_d, G22_d, G33_d, G12_d, G13_d, G23_d, nelv, lx) &
          bind(c, name='cuda_jacobi_update')
       use, intrinsic :: iso_c_binding
       type(c_ptr), value :: d_d, dxt_d, dyt_d, dzt_d
       type(c_ptr), value :: G11_d, G22_d, G33_d, G12_d, G13_d, G23_d
       integer(c_int) :: nelv, lx
     end subroutine cuda_jacobi_update
  end interface

  interface
     subroutine opencl_jacobi_update(d_d, dxt_d, dyt_d, dzt_d, &
          G11_d, G22_d, G33_d, G12_d, G13_d, G23_d, nelv, lx) &
          bind(c, name='opencl_jacobi_update')
       use, intrinsic :: iso_c_binding
       type(c_ptr), value :: d_d, dxt_d, dyt_d, dzt_d
       type(c_ptr), value :: G11_d, G22_d, G33_d, G12_d, G13_d, G23_d
       integer(c_int) :: nelv, lx
     end subroutine opencl_jacobi_update
  end interface

contains

  subroutine device_jacobi_init(this, coef, dof, gs_h)
    class(device_jacobi_t), intent(inout) :: this
    type(coef_t), intent(inout), target :: coef
    type(dofmap_t), intent(inout), target :: dof
    type(gs_t), intent(inout), target :: gs_h

    call this%free()

    this%gs_h => gs_h
    this%dof => dof
    this%coef => coef

    allocate(this%d(dof%Xh%lx,dof%Xh%ly,dof%Xh%lz, dof%msh%nelv))

    call device_map(this%d, this%d_d, size(this%d))

    call device_jacobi_update(this)

  end subroutine device_jacobi_init

  subroutine device_jacobi_free(this)
    class(device_jacobi_t), intent(inout) :: this

    if (c_associated(this%d_d)) then
       call device_free(this%d_d)
       this%d_d = C_NULL_PTR
    end if

    if (allocated(this%d)) then
       deallocate(this%d)
    end if

    nullify(this%dof)
    nullify(this%gs_h)
    nullify(this%coef)
  end subroutine device_jacobi_free

  !> The jacobi preconditioner \f$ J z = r \f$
  !! \f$ z = J^{-1}r\f$ where \f$ J^{-1} ~= 1/diag(A) \f$
  subroutine device_jacobi_solve(this, z, r, n)
    integer, intent(in) :: n
    class(device_jacobi_t), intent(inout) :: this
    real(kind=rp), dimension(n), intent(inout) :: z
    real(kind=rp), dimension(n), intent(inout) :: r
    type(c_ptr) :: z_d, r_d

    z_d = device_get_ptr(z)
    r_d = device_get_ptr(r)

    call device_col3(z_d, r_d, this%d_d, n)

  end subroutine device_jacobi_solve

  subroutine device_jacobi_update(this)
    class(device_jacobi_t), intent(inout) :: this
    integer :: lz, ly, lx
    associate(dof => this%dof, coef => this%coef, Xh => this%dof%Xh, &
         gs_h => this%gs_h, nelv => this%dof%msh%nelv)

      lx = Xh%lx
      ly = Xh%ly
      lz = Xh%lz

#ifdef HAVE_HIP
      call hip_jacobi_update(this%d_d, Xh%dxt_d, Xh%dyt_d, Xh%dzt_d, &
                             coef%G11_d, coef%G22_d, coef%G33_d, &
                             coef%G12_d, coef%G13_d, coef%G23_d, &
                             nelv, lx)
#elif HAVE_CUDA
      call cuda_jacobi_update(this%d_d, Xh%dxt_d, Xh%dyt_d, Xh%dzt_d, &
                             coef%G11_d, coef%G22_d, coef%G33_d, &
                             coef%G12_d, coef%G13_d, coef%G23_d, &
                             nelv, lx)
#elif HAVE_OPENCL
      call opencl_jacobi_update(this%d_d, Xh%dxt_d, Xh%dyt_d, Xh%dzt_d, &
                                coef%G11_d, coef%G22_d, coef%G33_d, &
                                coef%G12_d, coef%G13_d, coef%G23_d, &
                                nelv, lx)
#endif

      call device_col2(this%d_d, coef%h1_d, coef%dof%size())

      if (coef%ifh2) then
         call device_addcol3(this%d_d, coef%h2_d, coef%B_d, coef%dof%size())
      end if

      call gs_h%op(this%d, dof%size(), GS_OP_ADD)

      call device_invcol1(this%d_d, dof%size())
    end associate
  end subroutine device_jacobi_update

end module device_jacobi
