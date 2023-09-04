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
!> A simulation component that computes lambda2
!! The values are stored in the field registry under the name 'lambda2'

module lambda2
  use num_types, only : rp
  use neko_config
  use json_module, only : json_file
  use simulation_component, only : simulation_component_t
  use field_registry, only : neko_field_registry
  use field, only : field_t, field_free, field_init, field_ptr_t
  use field_list, only : field_list_t
  use operators, only : opgrad
  use math, only : pi, sort
  use case, only : case_t
  use device
  implicit none
  private
  
  type, public, extends(simulation_component_t) :: lambda2_t
     !> X velocity component.
     type(field_t), pointer :: u
     !> Y velocity component.
     type(field_t), pointer :: v
     !> Z velocity component.
     type(field_t), pointer :: w

     !> X lambda2 component.
     type(field_t), pointer :: lambda2

     !> Work arrays.
     type(field_t) :: temp1
     type(field_t) :: temp2

     contains
      !> Constructor from json.
      procedure, pass(this) :: init => lambda2_init_from_json
      !> Actual constructor.
      procedure, pass(this) :: init_from_attributes => &
        lambda2_init_from_attributes
      !> Destructor.
      procedure, pass(this) :: free => lambda2_free
      !> Compute the lambda2 field
     procedure, pass(this) :: compute => lambda2_compute
  end type lambda2_t
  
  contains
  
  !> Constructor from json.
  subroutine lambda2_init_from_json(this, json, case)
       class(lambda2_t), intent(inout) :: this
       type(json_file), intent(inout) :: json
       class(case_t), intent(inout), target ::case 
       
       call this%init_base(json, case)
       print *, 'hey'

       call lambda2_init_from_attributes(this)
  end subroutine lambda2_init_from_json

  !> Actual constructor.
  subroutine lambda2_init_from_attributes(this)
       class(lambda2_t), intent(inout) :: this

       this%u => neko_field_registry%get_field_by_name("u")
       this%v => neko_field_registry%get_field_by_name("v")
       this%w => neko_field_registry%get_field_by_name("w")

       if (.not. neko_field_registry%field_exists("lambda2")) then
          call neko_field_registry%add_field(this%u%dof, "lambda2")
       end if
       this%lambda2 => neko_field_registry%get_field_by_name("lambda2")

       call field_init(this%temp1, this%u%dof)
       call field_init(this%temp2, this%u%dof)
       call this%case%f_out%fluid%append(this%lambda2)
  end subroutine lambda2_init_from_attributes

  !> Destructor.
  subroutine lambda2_free(this)
       class(lambda2_t), intent(inout) :: this
       call field_free(this%temp1)
       call field_free(this%temp2)
  end subroutine lambda2_free

  !> Compute the lambda2 field.
  !! @param t The time value.
  !! @param tstep The current time-step
  subroutine lambda2_compute(this, t, tstep)
    class(lambda2_t), intent(inout) :: this
    real(kind=rp), intent(in) :: t
    integer, intent(in) :: tstep
    integer :: temp_indices(9), e, i, ind_sort(3)
    real(kind=rp) :: eigen(3), B, C, D, q, r, theta
    real(kind=rp) :: s11, s22, s33, s12, s13, s23, o12, o13, o23
    real(kind=rp) :: a11, a22, a33, a12, a13, a23
    type(field_ptr_t) :: grad(3,3)


    !> Get some fields to work with
    call this%case%fluid%scratch%request_field(grad(1,1)%f, temp_indices(1))
    call this%case%fluid%scratch%request_field(grad(1,2)%f, temp_indices(2))
    call this%case%fluid%scratch%request_field(grad(1,3)%f, temp_indices(3))
    call this%case%fluid%scratch%request_field(grad(2,1)%f, temp_indices(4))
    call this%case%fluid%scratch%request_field(grad(2,2)%f, temp_indices(5))
    call this%case%fluid%scratch%request_field(grad(2,3)%f, temp_indices(6))
    call this%case%fluid%scratch%request_field(grad(3,1)%f, temp_indices(7))
    call this%case%fluid%scratch%request_field(grad(3,2)%f, temp_indices(8))
    call this%case%fluid%scratch%request_field(grad(3,3)%f, temp_indices(9))

    call opgrad(grad(1,1)%f%x,grad(1,2)%f%x,grad(1,2)%f%x,this%u%x,this%case%fluid%c_Xh)
    call opgrad(grad(2,1)%f%x,grad(2,2)%f%x,grad(2,2)%f%x,this%v%x,this%case%fluid%c_Xh)
    call opgrad(grad(3,1)%f%x,grad(3,2)%f%x,grad(3,2)%f%x,this%w%x,this%case%fluid%c_Xh)

    if (NEKO_BCKND_DEVICE .eq. 1) then
        call device_memcpy(grad(1,1)%f%x,grad(1,1)%f%x_d,grad(1,1)%f%dof%size(),DEVICE_TO_HOST)
        call device_memcpy(grad(1,2)%f%x,grad(1,2)%f%x_d,grad(1,1)%f%dof%size(),DEVICE_TO_HOST)
        call device_memcpy(grad(1,3)%f%x,grad(1,3)%f%x_d,grad(1,1)%f%dof%size(),DEVICE_TO_HOST)
        call device_memcpy(grad(2,1)%f%x,grad(2,1)%f%x_d,grad(1,1)%f%dof%size(),DEVICE_TO_HOST)
        call device_memcpy(grad(2,2)%f%x,grad(2,2)%f%x_d,grad(1,1)%f%dof%size(),DEVICE_TO_HOST)
        call device_memcpy(grad(2,3)%f%x,grad(2,3)%f%x_d,grad(1,1)%f%dof%size(),DEVICE_TO_HOST)
        call device_memcpy(grad(3,1)%f%x,grad(3,1)%f%x_d,grad(1,1)%f%dof%size(),DEVICE_TO_HOST)
        call device_memcpy(grad(3,2)%f%x,grad(3,2)%f%x_d,grad(1,1)%f%dof%size(),DEVICE_TO_HOST)
        call device_memcpy(grad(3,3)%f%x,grad(3,3)%f%x_d,grad(1,1)%f%dof%size(),DEVICE_TO_HOST)
    end if
    do e = 1, this%lambda2%msh%nelv

       do i = 1, this%lambda2%Xh%lxyz
          s11 = grad(1,1)%f%x(i,1,1,e)
          s22 = grad(2,2)%f%x(i,1,1,e)
          s33 = grad(3,3)%f%x(i,1,1,e)

          
          s12 = 0.5*(grad(1,2)%f%x(i,1,1,e)+grad(1,2)%f%x(i,1,1,e))
          s13 = 0.5*(grad(1,3)%f%x(i,1,1,e)+grad(3,1)%f%x(i,1,1,e))
          s23 = 0.5*(grad(2,3)%f%x(i,1,1,e)+grad(3,2)%f%x(i,1,1,e))

          o12 = 0.5*(grad(1,2)%f%x(i,1,1,e)-grad(1,2)%f%x(i,1,1,e))
          o13 = 0.5*(grad(1,3)%f%x(i,1,1,e)-grad(3,1)%f%x(i,1,1,e))
          o23 = 0.5*(grad(2,3)%f%x(i,1,1,e)-grad(3,2)%f%x(i,1,1,e))

          a11 = s11*s11 + s12*s12 + s13*s13 - o12*o12 - o13*o13
          a12 = s11 * s12  +  s12 * s22  +  s13 * s23 - o13 * o23
          a13 = s11 * s13  +  s12 * s23  +  s13 * s33 + o12 * o23

          a22 = s12*s12 + s22*s22 + s23*s23 - o12*o12 - o23*o23
          a23 = s12 * s13 + s22 * s23 + s23 * s33 - o12 * o13
          a33 = s13*s13 + s23*s23 + s33*s33 - o13*o13 - o23*o23


          B = -(a11 + a22 + a33)
          C = -(a12*a12 + a13*a13 + a23*a23 - a11 * a22 - a11 * a33 - a22 * a33)
          D = -(2.0 * a12 * a13 * a23 - a11 * a23*a23 - a22 * a13*a13 - a33 * a12*a12  +  a11 * a22 * a33)


          q = (3.0 * C - B*B) / 9.0
          r = (9.0 * C * B - 27.0 * D - 2.0 * B*B*B) / 54.0
          theta = acos( r / sqrt(-q*q*q) )

          eigen(1) = 2.0 * sqrt(-q) * cos(theta / 3.0) - B / 3.0
          eigen(2) = 2.0 * sqrt(-q) * cos((theta + 2.0 * pi) / 3.0) - B / 3.0
          eigen(3) = 2.0 * sqrt(-q) * cos((theta + 4.0 * pi) / 3.0) - B / 3.0

          call sort(eigen,ind_sort,3)
          this%lambda2%x(i,1,1,e) = eigen(2) 
       end do
    end do
    if (NEKO_BCKND_DEVICE .eq. 1) then
        call device_memcpy(this%lambda2%x,this%lambda2%x_d,this%lambda2%dof%size(),HOST_TO_DEVICE)
    end if
    call this%case%fluid%scratch%relinquish_field(temp_indices)

  end subroutine lambda2_compute
  
end module lambda2
