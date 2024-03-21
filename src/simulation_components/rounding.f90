! Copyright (c) 2023, The Neko Authors
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
!
!> Implements the `rounding_t` type.

module rounding
  use num_types, only : rp, dp, sp
  use json_module, only : json_file
  use simulation_component, only : simulation_component_t
  use field_registry, only : neko_field_registry
  use field, only : field_t
  use neko_config
  use case, only : case_t
  use cpfloat_f
  use pcg_f
  use math
  use gather_scatter
  use device
  use comm
  implicit none
  private

  !> A simulation component that computes the rounding field.
  !! Added to the field registry as `omega_x`, `omega_y``, and `omega_z`.
  type, public, extends(simulation_component_t) :: rounding_t
     !> X velocity component.
     type(field_t), pointer :: u
     !> Y velocity component.
     type(field_t), pointer :: v
     !> Z velocity component.
     type(field_t), pointer :: w
     !> Pressure.
     type(field_t), pointer :: p

     real(kind=rp), allocatable :: mask(:)
     type(optstruct) :: fpopts
     integer(c_int128_t) :: seedinit, seedseq
     real :: noise 
   contains
     !> Constructor from json, wrapping the actual constructor.
     procedure, pass(this) :: init => rounding_init_from_json
     !> Actual constructor.
     procedure, pass(this) :: init_from_attributes => &
        rounding_init_from_attributes
     !> Destructor.
     procedure, pass(this) :: free => rounding_free
     !> Compute the rounding field.
     procedure, pass(this) :: compute_ => rounding_compute
  end type rounding_t

contains

  !> Constructor from json.
  subroutine rounding_init_from_json(this, json, case)
    class(rounding_t), intent(inout) :: this
    type(json_file), intent(inout) :: json
    class(case_t), intent(inout), target :: case
    character(len=:), allocatable :: filename
    character(len=:), allocatable :: precision

    call this%init_base(json, case)

    call rounding_init_from_attributes(this)
  end subroutine rounding_init_from_json

  !> Actual constructor.
  subroutine rounding_init_from_attributes(this, filename, precision)
    class(rounding_t), intent(inout) :: this
    character(len=*), intent(in), optional :: filename
    integer, intent(in), optional :: precision
    type(field_t) :: temp
    integer :: i, n, argc, ierr
    character(len=1024) :: inputchar

    this%u => neko_field_registry%get_field_by_name("u")
    this%v => neko_field_registry%get_field_by_name("v")
    this%w => neko_field_registry%get_field_by_name("w")
    this%p => neko_field_registry%get_field_by_name("p")

    n = this%u%dof%size()
    call temp%init(this%u%dof,'rand')
    allocate(this%mask(this%u%dof%size()))
    ! give each point a unique value based on numbering
    do i = 1, this%u%dof%size()
       temp%x(i,1,1,1) = dble(i)+dble(this%u%msh%offset_el)*dble(this%u%Xh%lx**3)
    end do
    call copy(this%mask,temp%x,this%u%dof%size())
    if (NEKO_BCKND_DEVICE .eq. 1) &
       call device_memcpy(temp%x, temp%x_d, this%u%dof%size(), HOST_TO_DEVICE,sync=.true.)
    call this%case%fluid%c_Xh%gs_h%op(temp%x,this%u%dof%size(),GS_OP_MAX)
    if (NEKO_BCKND_DEVICE .eq. 1) &
       call device_memcpy(temp%x, temp%x_d, this%u%dof%size(), DEVICE_TO_HOST,sync=.true.)

    do i = 1, this%u%dof%size()
       ! If this point is the one with the highest number, set mask to 1.0
       if (abs(temp%x(i,1,1,1) -this%mask(i)) .lt. 1e-8) then
          this%mask(i) = 1.0
       else
          this%mask(i) = 0.0
       end if
    end do
    argc = command_argument_count()
    if (argc .ne. 4) then
       if(pe_rank .eq. 0) print *, 'ERROR: not enough args'
       if(pe_rank .eq. 0) print *, './bla precision emax rounding'
       if(pe_rank .eq. 0) print *, 'FP64: precision=53 emax=1023 '
       if(pe_rank .eq. 0) print *, 'FP32: 24 127'
       if(pe_rank .eq. 0) print *, 'FP16: 11 15'
       if(pe_rank .eq. 0) print *, 'bfloat16: 8 127'
       if(pe_rank .eq. 0) print *, 'FP8-e4m3: 4 7'
       if(pe_rank .eq. 0) print *, 'FP8-e5m2: 3 15'
       if(pe_rank .eq. 0) print *, 'Rounding: 1=round to nearest'
       if(pe_rank .eq. 0) print *, 'Rounding: 5=Stochastic rounding'
       if(pe_rank .eq. 0) print *, 'Rounding: 10=noise'
       call exit()
    end if
    call get_command_argument(2, inputchar)
    read(inputchar, *) this%fpopts%precision
    call get_command_argument(3, inputchar)
    read(inputchar, *) this%fpopts%emax
    call get_command_argument(4, inputchar)
    read(inputchar, *) this%fpopts%round
 
    !fpopts%precision = 4
    !fpopts%emax = 9
    !fpopts%round = 1 !CPFLOAT_RND_TO EVEN;         
    !fpopts%round = 5 !CPFLOAT_RND_prop_stoc;         
    this%fpopts%subnormal = 1!   // Support for subnormals is on.
    this%fpopts%flip = 0 !CPFLOAT_NO_SOFTERR;      // Bit flips are off.
    this%fpopts%p = 0!;                          // Bit flip probability (not used).
    this%fpopts%explim = 1!CPFLOAT_EXPRANGE_TARG; // Limited exponent in target format.

    ierr = cpfloat_validate_optstruct(this%fpopts)
    this%seedinit = pe_rank
    this%seedseq = pe_rank 
    call pcg64_srandom(this%seedinit, this%seedseq)

    if (pe_rank.eq.0) then
       write(*,*) 'cpfloat init', ierr, n
       write(*,*) 'precision', this%fpopts%precision
       write(*,*) 'emax', this%fpopts%emax
       write(*,*) 'subnormal', this%fpopts%subnormal
       write(*,*) 'rounding 5=stochastic', this%fpopts%round
       write(*,*) 'flip', this%fpopts%flip
       write(*,*) 'flip prob', this%fpopts%p
       write(*,*) 'explim', this%fpopts%explim
    end if

  end subroutine rounding_init_from_attributes

  subroutine perturb(u, n, fpopts)
    type(optstruct), intent(in) :: fpopts
    integer(kind=INT64) :: n 
    real(kind=rp), intent(inout) :: u(n)
    integer :: ierr, i
    real(kind=rp) :: noise
    integer(c_int64_t) :: rand
    
    if (fpopts%round .eq. 10) then 
       do i = 1, n
          rand = pcg64_random()
          noise = scale(dble(abs(rand)),-63)
          noise = (noise -0.5)* 2.0**(-fpopts%precision+1)
          if (abs(noise) .gt. 2.0**(-fpopts%precision)) then
             print *,'noise',noise
             call exit()
          end if
          u(i) = (1.0+noise)*u(i)
       end do  
    else 
       ierr = cpfloat(u, u, n, fpopts)
    end if
  end subroutine perturb


  !> Destructor.
  subroutine rounding_free(this)
    class(rounding_t), intent(inout) :: this
  end subroutine rounding_free

  !> Compute the rounding field.
  !! @param t The time value.
  !! @param tstep The current time-step
  subroutine rounding_compute(this, t, tstep)
    class(rounding_t), intent(inout) :: this
    real(kind=rp), intent(in) :: t
    integer, intent(in) :: tstep
    integer(kind=INT64) :: n
    n = this%u%dof%size()

    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_memcpy(this%u%x, this%u%x_d, this%u%dof%size(), DEVICE_TO_HOST,sync=.true.)
       call device_memcpy(this%v%x, this%v%x_d, this%u%dof%size(), DEVICE_TO_HOST,sync=.true.)
       call device_memcpy(this%w%x, this%w%x_d, this%u%dof%size(), DEVICE_TO_HOST,sync=.true.)
       call device_memcpy(this%p%x, this%p%x_d, this%u%dof%size(), DEVICE_TO_HOST,sync=.true.)
    end if

    call perturb(this%u%x, n,this%fpopts) 
    call perturb(this%v%x, n,this%fpopts) 
    call perturb(this%w%x, n,this%fpopts) 
    call perturb(this%p%x, n,this%fpopts) 
    call col2(this%u%x,this%mask,this%u%dof%size())
    call col2(this%v%x,this%mask,this%u%dof%size())
    call col2(this%w%x,this%mask,this%u%dof%size())
    call col2(this%p%x,this%mask,this%u%dof%size())

    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_memcpy(this%u%x, this%u%x_d, this%u%dof%size(), HOST_TO_DEVICE,sync=.false.)
       call device_memcpy(this%v%x, this%v%x_d, this%u%dof%size(), HOST_TO_DEVICE,sync=.false.)
       call device_memcpy(this%w%x, this%w%x_d, this%u%dof%size(), HOST_TO_DEVICE,sync=.false.)
       call device_memcpy(this%p%x, this%p%x_d, this%u%dof%size(), HOST_TO_DEVICE,sync=.true.)
    end if
    call this%case%fluid%c_Xh%gs_h%op(this%u%x,this%u%dof%size(),GS_OP_ADD)
    call this%case%fluid%c_Xh%gs_h%op(this%v%x,this%u%dof%size(),GS_OP_ADD)
    call this%case%fluid%c_Xh%gs_h%op(this%w%x,this%u%dof%size(),GS_OP_ADD)
    call this%case%fluid%c_Xh%gs_h%op(this%p%x,this%u%dof%size(),GS_OP_ADD)

   end subroutine rounding_compute

end module rounding
