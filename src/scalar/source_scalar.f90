! Copyright (c) 2020-2021, The Neko Authors
! All rights reserved.
!
! Redistribution and use in source_scalar and binary forms, with or without
! modification, are permitted provided that the following conditions
! are met:
!
!   * Redistributions of source_scalar code must retain the above copyright
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
!> Source terms for scalars
module source_scalar
  use neko_config
  use num_types
  use dofmap
  use utils
  use device
  use device_math
  use, intrinsic :: iso_c_binding
  implicit none

  !> Defines a source term for the scalar transport equation term \f$ f \f$
  type :: source_scalar_t
     real(kind=rp), allocatable :: s(:,:,:,:) !< the source values
     type(dofmap_t), pointer :: dm            !< dofmap for the given space
     type(c_ptr) :: s_d = C_NULL_PTR          !< dev. ptr for the forces
     procedure(source_scalar_term), pass(f), pointer  :: eval => null()
     procedure(source_scalar_term_pw), nopass, pointer  :: eval_pw => null()
  end type source_scalar_t

  !> Abstract interface defining how to compute a source_scalar term
  abstract interface
     subroutine source_scalar_term(f, t)
       import source_scalar_t
       import rp
       class(source_scalar_t), intent(inout) :: f
       real(kind=rp), intent(in) :: t
     end subroutine source_scalar_term
  end interface

  !> Abstract interface defining how to compute a source_scalar term pointwise
  abstract interface
     subroutine source_scalar_term_pw(s, j, k, l, e, t)
       import rp
       real(kind=rp), intent(inout) :: s
       integer, intent(in) :: j
       integer, intent(in) :: k
       integer, intent(in) :: l
       integer, intent(in) :: e
       real(kind=rp), intent(in) :: t
     end subroutine source_scalar_term_pw
  end interface
  
contains

  !> Initialize a source_scalar term @a f
  subroutine source_scalar_init(f, dm)
    type(source_scalar_t), intent(inout) :: f
    type(dofmap_t), intent(inout), target :: dm

    call source_scalar_free(f)

    f%dm => dm

    allocate(f%s(dm%Xh%lx, dm%Xh%ly, dm%Xh%lz, dm%msh%nelv))

    f%s = 0d0

    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_map(f%s, f%s_d, dm%size())
    end if
    
  end subroutine source_scalar_init

  !> Deallocate a source_scalar term @a f
  subroutine source_scalar_free(f)
    type(source_scalar_t), intent(inout) :: f

    if (allocated(f%s)) then
       deallocate(f%s)
    end if

    nullify(f%dm)

    if (c_associated(f%s_d)) then
       call device_free(f%s_d)
    end if
    
  end subroutine source_scalar_free

  !> Set the eval method for the source_scalar term @a f
  subroutine source_scalar_set_type(f, f_eval)
    type(source_scalar_t), intent(inout) :: f
    procedure(source_scalar_term) :: f_eval
    f%eval => f_eval
  end subroutine source_scalar_set_type

  !> Set the pointwise eval method for the source_scalar term @a f
  subroutine source_scalar_set_pw_type(f, f_eval_pw)
    type(source_scalar_t), intent(inout) :: f
    procedure(source_scalar_term_pw) :: f_eval_pw
    if (NEKO_BCKND_DEVICE .eq. 1) then
       call neko_error('Pointwise source_scalar terms not supported on accelerators')
    end if
    f%eval => source_scalar_eval_pw
    f%eval_pw => f_eval_pw
  end subroutine source_scalar_set_pw_type

  !> Eval routine for zero forcing
  !! @note Maybe this should be cache, avoding zeroing at each time-step
  subroutine source_scalar_eval_noforce(f, t)
    class(source_scalar_t), intent(inout) :: f
    real(kind=rp), intent(in) :: t
    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_rzero(f%s_d, f%dm%size())
    else
       f%s = 0d0
    end if
  end subroutine source_scalar_eval_noforce

  !> Driver for all pointwise source_scalar term evaluatons
  subroutine source_scalar_eval_pw(f, t)
    class(source_scalar_t), intent(inout) :: f
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
                call f%eval_pw(f%s(j,k,l,e), jj, kk, ll, ee, t)
             end do
          end do
       end do
    end do
    
  end subroutine source_scalar_eval_pw
  
end module source_scalar
