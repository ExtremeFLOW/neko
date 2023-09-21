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
!> Implements data streaming.
module data_streamer

  use num_types, only: rp
  use field, only: field_t
  use coefs, only: coef_t
 
  use device 
  use comm
  use mpi_types
  use neko_config

  use, intrinsic :: iso_c_binding
  implicit none
  private
  
  type, public :: data_streamer_t
       
     !> Define if the execution is asyncrhonous
     integer :: if_asynch

     !> global element numbers
     integer, allocatable :: lglel(:)

     contains
       !> Constructor
       procedure, pass(this) :: init => data_streamer_init
       !> Destructor
       procedure, pass(this) :: free => data_streamer_free
       !> Stream data
       procedure, pass(this) :: stream => data_streamer_stream

    end type data_streamer_t

contains

  !> Constructor
  subroutine data_streamer_init(this, coef, if_asynch)
    class(data_streamer_t), intent(inout) :: this
    type(coef_t), intent(in) :: coef
    integer, intent(in) :: if_asynch
    integer :: nelb, nelb2, nelv, nelgv,npts,e

    !Allocate and initialize the global element number
    allocate(this%lglel(coef%msh%nelv))
    do e = 1, coef%msh%nelv
       this%lglel(e) = e + coef%msh%offset_el
    end do

    !Assign if the streaming is asynchronous
    this%if_asynch = if_asynch

    !Assign the set up parameters
    nelv  = coef%msh%nelv
    npts  = coef%Xh%lx*coef%Xh%ly*coef%Xh%lz
    nelgv = coef%msh%glb_nelv
    nelb  = coef%msh%offset_el
    ! Alternative way to get nelb:
    !nelb = elem_running_sum(nelv)
    !nelb = nelb - nelv
!!#ifdef HAVE_ADIOS
    call adios2_setup(npts,nelv,nelb,nelgv, &
                nelgv,coef%dof%x,coef%dof%y,  &
                coef%dof%z,if_asynch,NEKO_COMM)
!!#else
    call neko_error('NEKO needs to be built with ADIOS2 support')
!!#endif


  end subroutine data_streamer_init

  !> Destructor
  subroutine data_streamer_free(this)
    class(data_streamer_t), intent(inout) :: this

    if (allocated(this%lglel))        deallocate(this%lglel)

!!#ifdef HAVE_ADIOS2
    call adios2_finalize()
!!#else
    call neko_error('NEKO needs to be built with ADIOS2 support')
!!#endif
    


  end subroutine data_streamer_free
  
  !> streamer
  subroutine data_streamer_stream(this,u,v,w,p,coef)
    class(data_streamer_t), intent(inout) :: this
    type(coef_t), intent(inout) :: coef
    type(field_t), intent(inout) :: u
    type(field_t), intent(inout) :: v
    type(field_t), intent(inout) :: w
    type(field_t), intent(inout) :: p
    integer :: nelv, npts
    
    nelv  = coef%msh%nelv
    npts  = coef%Xh%lx*coef%Xh%ly*coef%Xh%lz
    
    if (NEKO_BCKND_DEVICE .eq. 1) then 
       ! Move the data to the CPU to be able to write it
       call device_memcpy(u%x,  u%x_d, nelv*npts,DEVICE_TO_HOST)
       call device_memcpy(v%x,  v%x_d, nelv*npts,DEVICE_TO_HOST)
       call device_memcpy(w%x,  w%x_d, nelv*npts,DEVICE_TO_HOST)
       call device_memcpy(p%x,  p%x_d, nelv*npts,DEVICE_TO_HOST)
    end if

#if HAVE_ADIOS2
    call adios2_stream(this%lglel,p%x, u%x, v%x, w%x, coef%B, u%x)
#else
    call neko_error('NEKO needs to be built with ADIOS2 support')
#endif


  end subroutine data_streamer_stream
  
  
  function elem_running_sum(nelv) result(rbuff)

    integer, intent(in) :: nelv
    integer ::  ierr,xbuff,wbuff,rbuff

    xbuff = nelv  ! running sum
    wbuff = nelv  ! working buff
    rbuff = 0   ! recv buff
     
    call mpi_scan(xbuff,rbuff,1,mpi_integer,mpi_sum,NEKO_COMM,ierr)

  end function elem_running_sum

end module data_streamer
