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
!> Source terms
module source
  use neko_config
  use num_types
  use dofmap
  use utils
  use device
  use device_math
  use, intrinsic :: iso_c_binding
  implicit none

  !> Defines a source term \f$ f \f$
  type :: source_t
     real(kind=rp), allocatable :: u(:,:,:,:) !< x-component of source term
     real(kind=rp), allocatable :: v(:,:,:,:) !< y-component of source term
     real(kind=rp), allocatable :: w(:,:,:,:) !< w-component of source term
     type(dofmap_t), pointer :: dm            !< dofmap for the given space
     type(c_ptr) :: u_d = C_NULL_PTR          !< dev. ptr for x-component
     type(c_ptr) :: v_d = C_NULL_PTR          !< dev. ptr for y-component
     type(c_ptr) :: w_d = C_NULL_PTR          !< dev. ptr for z-component
     procedure(source_term), pass(f), pointer  :: eval => null()
     procedure(source_term_pw), nopass, pointer  :: eval_pw => null()
  end type source_t

  !> Abstract interface defining how to compute a source term
  abstract interface
     subroutine source_term(f, t)
       import source_t
       import rp
       class(source_t), intent(inout) :: f
       real(kind=rp), intent(in) :: t
     end subroutine source_term
  end interface

  !> Abstract interface defining how to compute a source term pointwise
  abstract interface
     subroutine source_term_pw(u, v, w, j, k, l, e, t)
       import rp
       real(kind=rp), intent(inout) :: u
       real(kind=rp), intent(inout) :: v
       real(kind=rp), intent(inout) :: w
       integer, intent(in) :: j
       integer, intent(in) :: k
       integer, intent(in) :: l
       integer, intent(in) :: e
       real(kind=rp), intent(in) :: t
     end subroutine source_term_pw
  end interface
  
contains

  !> Initialize a source term @a f
  subroutine source_init(f, dm)
    type(source_t), intent(inout) :: f
    type(dofmap_t), intent(inout), target :: dm

    call source_free(f)

    f%dm => dm

    allocate(f%u(dm%Xh%lx, dm%Xh%ly, dm%Xh%lz, dm%msh%nelv))
    allocate(f%v(dm%Xh%lx, dm%Xh%ly, dm%Xh%lz, dm%msh%nelv))
    allocate(f%w(dm%Xh%lx, dm%Xh%ly, dm%Xh%lz, dm%msh%nelv))

    f%u = 0d0
    f%v = 0d0
    f%w = 0d0

    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_map(f%u, f%u_d, dm%size())
       call device_map(f%v, f%v_d, dm%size())
       call device_map(f%w, f%w_d, dm%size())
    end if
    
  end subroutine source_init

  !> Deallocate a source term @a f
  subroutine source_free(f)
    type(source_t), intent(inout) :: f

    if (allocated(f%u)) then
       deallocate(f%u)
    end if

    if (allocated(f%v)) then
       deallocate(f%v)
    end if

    if (allocated(f%w)) then
       deallocate(f%w)
    end if

    nullify(f%dm)

    if (c_associated(f%u_d)) then
       call device_free(f%u_d)
    end if

    if (c_associated(f%v_d)) then
       call device_free(f%v_d)
    end if

    if (c_associated(f%w_d)) then
       call device_free(f%w_d)
    end if
    
  end subroutine source_free

  !> Set the eval method for the source term @a f
  subroutine source_set_type(f, f_eval)
    type(source_t), intent(inout) :: f
    procedure(source_term) :: f_eval
    f%eval => f_eval
  end subroutine source_set_type

  !> Set the pointwise eval method for the source term @a f
  subroutine source_set_pw_type(f, f_eval_pw)
    type(source_t), intent(inout) :: f
    procedure(source_term_pw) :: f_eval_pw
    if (NEKO_BCKND_DEVICE .eq. 1) then
       call neko_error('Pointwise source terms not supported on accelerators')
    end if
    f%eval => source_eval_pw
    f%eval_pw => f_eval_pw
  end subroutine source_set_pw_type

  !> Eval routine for zero forcing
  !! @note Maybe this should be cache, avoding zeroing at each time-step
  subroutine source_eval_noforce(f, t)
    class(source_t), intent(inout) :: f
    real(kind=rp), intent(in) :: t
    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_rzero(f%u_d, f%dm%size())
       call device_rzero(f%v_d, f%dm%size())
       call device_rzero(f%w_d, f%dm%size())
    else
       f%u = 0d0
       f%v = 0d0
       f%w = 0d0
    end if
  end subroutine source_eval_noforce

  !> Driver for all pointwise source term evaluatons
  subroutine source_eval_pw(f, t)
    class(source_t), intent(inout) :: f
    real(kind=rp), intent(in) :: t
    integer :: j, k, l, e
    integer :: jj,kk,ll,ee

    do e = 1, f%dm%msh%nelv
       ee = e
       do l = 1, f%dm%Xh%lz
          ll = l
          do k = 1, f%dm%Xh%ly
             kk = k
             do j = 1, f%dm%Xh%lx
                jj =j
                call f%eval_pw(f%u(j,k,l,e), f%v(j,k,l,e), f%w(j,k,l,e), &
                     jj, kk, ll, ee, t)
             end do
          end do
       end do
    end do
    
  end subroutine source_eval_pw
  
end module source
