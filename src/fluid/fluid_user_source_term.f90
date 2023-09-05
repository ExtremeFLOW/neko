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
module fluid_user_source_term
  use neko_config, only : NEKO_BCKND_DEVICE
  use num_types, only : rp
  use dofmap, only : dofmap_t
  use utils, only : neko_error
  use device
  use device_math
  use, intrinsic :: iso_c_binding
  implicit none

  !> Defines a source term \f$ f \f$
  type :: fluid_user_source_term_t
     real(kind=rp), allocatable :: u(:,:,:,:) !< x-component of source term
     real(kind=rp), allocatable :: v(:,:,:,:) !< y-component of source term
     real(kind=rp), allocatable :: w(:,:,:,:) !< w-component of source term
     type(dofmap_t), pointer :: dm            !< dofmap for the given space
     type(c_ptr) :: u_d = C_NULL_PTR          !< dev. ptr for x-component
     type(c_ptr) :: v_d = C_NULL_PTR          !< dev. ptr for y-component
     type(c_ptr) :: w_d = C_NULL_PTR          !< dev. ptr for z-component
     procedure(source_term), pass(f), pointer  :: eval => null()
     procedure(source_term_pw), nopass, pointer  :: eval_pw => null()
   contains
     !> Constructor.
     procedure, pass(this) :: init => fluid_user_source_term_init
     !> Destructor.
     procedure, pass(this) :: free => fluid_user_source_term_free
     !> Set the source type (no force, user pointwise or user vector)
     procedure, pass(this) :: set_source_type => &
       fluid_user_source_term_set_source_type
  end type fluid_user_source_term_t

  !> Abstract interface defining how to compute a source term
  abstract interface
     subroutine source_term(f, t)
       import fluid_user_source_term_t
       import rp
       class(fluid_user_source_term_t), intent(inout) :: f
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

  !> Costructor.
  subroutine fluid_user_source_term_init(this,dm)
    class(fluid_user_source_term_t), intent(inout) :: this
    type(dofmap_t), intent(inout), target :: dm

    call this%free()

    this%dm => dm

    allocate(this%u(dm%Xh%lx, dm%Xh%ly, dm%Xh%lz, dm%msh%nelv))
    allocate(this%v(dm%Xh%lx, dm%Xh%ly, dm%Xh%lz, dm%msh%nelv))
    allocate(this%w(dm%Xh%lx, dm%Xh%ly, dm%Xh%lz, dm%msh%nelv))

    this%u = 0d0
    this%v = 0d0
    this%w = 0d0

    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_map(this%u, this%u_d, dm%size())
       call device_map(this%v, this%v_d, dm%size())
       call device_map(this%w, this%w_d, dm%size())
    end if
    
  end subroutine fluid_user_source_term_init

  !> Destructctor.
  subroutine fluid_user_source_term_free(this)
    class(fluid_user_source_term_t), intent(inout) :: this

    if (allocated(this%u)) then
       deallocate(this%u)
    end if

    if (allocated(this%v)) then
       deallocate(this%v)
    end if

    if (allocated(this%w)) then
       deallocate(this%w)
    end if

    nullify(this%dm)

    if (c_associated(this%u_d)) then
       call device_free(this%u_d)
    end if

    if (c_associated(this%v_d)) then
       call device_free(this%v_d)
    end if

    if (c_associated(this%w_d)) then
       call device_free(this%w_d)
    end if
  end subroutine fluid_user_source_term_free

  !> Set the source type (no force, user pointwise or user vector).
  !! @param source_term_type The name of the type of the term: 'noforce', 'user'
  !! or 'user_vector'.
  subroutine fluid_user_source_term_set_source_type(this, source_term_type, &
    user_proc_pw, user_proc_vector)
    class(fluid_user_source_term_t), intent(inout) :: this
    character(len=*) :: source_term_type
    procedure(source_term_pw), optional :: user_proc_pw
    procedure(source_term), optional :: user_proc_vector

    if (trim(source_term_type) .eq. 'noforce') then
       call source_set_type(this, source_eval_noforce)
    else if (trim(source_term_type) .eq. 'user' .and. &
              present(user_proc_pw)) then
       call source_set_pw_type(this, user_proc_pw)
    else if (trim(source_term_type) .eq. 'user_vector' .and. &
             present(user_proc_vector)) then
       call source_set_type(this, user_proc_vector)
    else
       call neko_error('Invalid fluid source term '//source_term_type)
    end if
 
  end subroutine fluid_user_source_term_set_source_type

  !> Set the eval method for the source term @a f
  subroutine source_set_type(f, f_eval)
    type(fluid_user_source_term_t), intent(inout) :: f
    procedure(source_term) :: f_eval
    f%eval => f_eval
  end subroutine source_set_type

  !> Set the pointwise eval method for the source term @a f
  subroutine source_set_pw_type(f, f_eval_pw)
    type(fluid_user_source_term_t), intent(inout) :: f
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
    class(fluid_user_source_term_t), intent(inout) :: f
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
    class(fluid_user_source_term_t), intent(inout) :: f
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
  
end module fluid_user_source_term
