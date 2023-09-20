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

  use neko

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

    call adios2_setup(npts,nelv,nelb,nelgv, &
                nelgv,coef%dof%x,coef%dof%y,  &
                coef%dof%z,if_asynch,NEKO_COMM)

  end subroutine data_streamer_init

  !> Destructor
  subroutine data_streamer_free(this)
    class(data_streamer_t), intent(inout) :: this

    if (allocated(this%lglel))        deallocate(this%lglel)
    
    call adios2_finalize()

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

    call adios2_stream(this%lglel,p%x, u%x, v%x, w%x, coef%B, u%x)

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

module user
  use neko
  use data_streamer

  implicit none

  !> Variables to store the Rayleigh and Prandlt numbers
  real(kind=rp) :: Ra = 0
  real(kind=rp) :: Re = 0
  real(kind=rp) :: Pr = 0

  type(data_streamer_t) :: dstream

contains
  ! Register user defined functions (see user_intf.f90)
  subroutine user_setup(u)
    type(user_t), intent(inout) :: u
    u%user_init_modules => user_initialize
    u%user_finalize_modules => user_finalize
    u%fluid_user_ic => set_initial_conditions_for_u_and_s
    u%scalar_user_bc => set_scalar_boundary_conditions
    u%fluid_user_f_vector => set_bousinesq_forcing_term
    u%user_check => check
  end subroutine user_setup
 
  subroutine set_scalar_boundary_conditions(s, x, y, z, nx, ny, nz, ix, iy, iz, ie, t, tstep)
    real(kind=rp), intent(inout) :: s
    real(kind=rp), intent(in) :: x
    real(kind=rp), intent(in) :: y
    real(kind=rp), intent(in) :: z
    real(kind=rp), intent(in) :: nx
    real(kind=rp), intent(in) :: ny
    real(kind=rp), intent(in) :: nz
    integer, intent(in) :: ix
    integer, intent(in) :: iy
    integer, intent(in) :: iz
    integer, intent(in) :: ie
    real(kind=rp), intent(in) :: t
    integer, intent(in) :: tstep

    !> Variables for bias
    real(kind=rp) :: arg, bias
    
    ! This will be used on all zones without labels
    ! e.g. the ones hardcoded to 'v', 'w', etcetc
    s = 1.0_rp - z

  end subroutine set_scalar_boundary_conditions

  subroutine set_initial_conditions_for_u_and_s(u, v, w, p, params)
    type(field_t), intent(inout) :: u
    type(field_t), intent(inout) :: v
    type(field_t), intent(inout) :: w
    type(field_t), intent(inout) :: p
    type(json_file), intent(inout) :: params
    type(field_t), pointer :: s
    integer :: i, j, k, e
    real(kind=rp) :: rand, r,z
    s => neko_field_registry%get_field('s')

    !> Initialize with zero velocity
    call rzero(u%x,u%dof%size())
    call rzero(v%x,v%dof%size())
    call rzero(w%x,w%dof%size())

    !> Initialize with rand perturbations on temperature
    call rzero(s%x,w%dof%size())
    do i = 1, s%dof%size()
       s%x(i,1,1,1) = 1-s%dof%z(i,1,1,1) 
    end do
    ! perturb not on element boundaries
    ! Maybe not necessary, but lets be safe
    do e = 1, s%msh%nelv
       do k = 2,s%Xh%lx-1
          do j = 2,s%Xh%lx-1
             do i = 2,s%Xh%lx-1

                !call random_number(rand)
                !Somewhat random
                rand = cos(real(e+s%msh%offset_el,rp)*real(i*j*k,rp))
                r = sqrt(s%dof%x(i,j,k,e)**2+s%dof%y(i,j,k,e)**2)
                z = s%dof%z(i,j,k,e)
                s%x(i,j,k,e) = 1-z + 0.0001*rand*s%dof%x(i,j,k,e)*&
                                                    sin(3*pi*r/0.05_rp)*sin(10*pi*z)
             end do
          end do
       end do
    end do
    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_memcpy(s%x,s%x_d,s%dof%size(),HOST_TO_DEVICE)
    end if

  end subroutine set_initial_conditions_for_u_and_s

  subroutine user_initialize(t, u, v, w, p, coef, params)
    real(kind=rp) :: t
    type(field_t), intent(inout) :: u
    type(field_t), intent(inout) :: v
    type(field_t), intent(inout) :: w
    type(field_t), intent(inout) :: p
    type(coef_t), intent(inout) :: coef
    type(json_file), intent(inout) :: params

    !> Log variable
    character(len=LOG_SIZE) :: log_buf ! For logging status

    !> Support variables for probes 
    integer :: nelb, nelb2, nelv, nelgv,npts, if_asynch

    !> Recalculate the non dimensional parameters
    call json_get(params, 'case.scalar.Pr', Pr)
    call json_get(params, 'case.fluid.Re', Re)
    Ra = (Re**2)*Pr
    write(log_buf,*) 'Rayleigh Number is Ra=', Ra
    call neko_log%message(log_buf)
    
    nelv  = u%msh%nelv
    npts  = u%Xh%lx**3
    nelgv = u%msh%glb_nelv
    if_asynch = 1
    nelb = elem_running_sum(nelv)    
    nelb = nelb - nelv 
    nelb2 = u%msh%offset_el

    
    !write(*,*) 'my sum is', nelb
    !write(*,*) 'my sum_neko is', nelb2
    !write(*,*) 'my number of points is', npts
    !write(*,*) 'my number of elements is', nelv
    !write(*,*) 'total number of elements is', nelgv

    !call adios2_setup(npts,nelv,nelb2,nelgv, &
    !            nelgv,u%dof%x,u%dof%y,  &
    !            u%dof%z,if_asynch,NEKO_COMM)

    call dstream%init(coef, 1)


  end subroutine user_initialize

  subroutine user_finalize(t, param)
    real(kind=rp) :: t
    type(json_file), intent(inout) :: param
         
    !call adios2_finalize()
    call dstream%free()
     
  end subroutine user_finalize

  subroutine set_bousinesq_forcing_term(f, t)
    class(source_t), intent(inout) :: f
    real(kind=rp), intent(in) :: t
    integer :: i
    type(field_t), pointer :: u, v, w, s
    real(kind=rp) :: rapr, ta2pr
    u => neko_field_registry%get_field('u')
    v => neko_field_registry%get_field('v')
    w => neko_field_registry%get_field('w')
    s => neko_field_registry%get_field('s')

    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_rzero(f%u_d,f%dm%size())
       call device_rzero(f%v_d,f%dm%size())
       call device_copy(f%w_d,s%x_d,f%dm%size())
    else
       call rzero(f%u,f%dm%size())
       call rzero(f%v,f%dm%size())
       call copy(f%w,s%x,f%dm%size())
    end if
  end subroutine set_bousinesq_forcing_term

  subroutine check(t, tstep, u, v, w, p, coef, params)
    real(kind=rp), intent(in) :: t
    integer, intent(in) :: tstep
    type(coef_t), intent(inout) :: coef
    type(json_file), intent(inout) :: params
    type(field_t), intent(inout) :: u
    type(field_t), intent(inout) :: v
    type(field_t), intent(inout) :: w
    type(field_t), intent(inout) :: p
    
    
    integer ::  lglel(u%msh%nelv)
    integer ::  n, npts, nelv

    nelv = u%msh%nelv
    npts  = u%Xh%lx**3
    n = u%msh%nelv*npts

    !if ((NEKO_BCKND_HIP .eq. 1) .or. (NEKO_BCKND_CUDA .eq. 1) .or. &
    !       (NEKO_BCKND_OPENCL .eq. 1)) then

    !       call device_memcpy(u%x,  u%x_d, &
    !                   nelv*npts,                         &
    !                   DEVICE_TO_HOST)
    !    write(*,*) "copy the data to the cpu"
    !    !call device_memcpy(u%x, u%x_d, n, DEVICE_TO_HOST,synch = .true.)
    
    !end if

    !!Write the fields into a file
    !call adios2_stream(lglel,u%x, u%x, &
    !                    u%x, u%x, coef%B, u%x)

    call dstream%stream(u,v,w,p,coef)

  end subroutine check
  
  
  function elem_running_sum(nelv) result(rbuff)

    integer, intent(in) :: nelv
    integer ::  ierr,xbuff,wbuff,rbuff

    xbuff = nelv  ! running sum
    wbuff = nelv  ! working buff
    rbuff = 0   ! recv buff
     
    call mpi_scan(xbuff,rbuff,1,mpi_integer,mpi_sum,NEKO_COMM,ierr)

  end function elem_running_sum

end module user
