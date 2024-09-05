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
!> Implements the `force_torque_t` type.

module force_torque
  use num_types, only : rp, dp, sp
  use json_module, only : json_file
  use simulation_component, only : simulation_component_t
  use field_registry, only : neko_field_registry
  use field, only : field_t
  use operators, only : curl
  use case, only : case_t
  use fld_file_output, only : fld_file_output_t
  use json_utils, only : json_get, json_get_or_default
  use field_writer, only : field_writer_t
  use coefs, only : coef_t
  use operators, only : strain_rate
  use vector, only : vector_t
  use dirichlet
  use drag_torque
  use comm
  use math, only : masked_red_copy, cadd, glsum
  use device_math, only : device_masked_red_copy, device_cadd, device_glsum
  
  implicit none
  private

  !> A simulation component that computes the force_torque field.
  !! Added to the field registry as `omega_x`, `omega_y``, and `omega_z`.
  type, public, extends(simulation_component_t) :: force_torque_t
     !> X velocity component.
     type(field_t), pointer :: u
     !> Y velocity component.
     type(field_t), pointer :: v
     !> Z velocity component.
     type(field_t), pointer :: w
     !> Pressure.
     type(field_t), pointer :: p

     type(field_t) :: s11, s22, s33, s12, s13, s23
     !Masked working arrays
     type(vector_t) :: n1, n2, n3
     type(vector_t) :: r1, r2, r3
     type(vector_t) :: force1, force2, force3
     type(vector_t) :: force4, force5, force6
     type(vector_t) :: pmsk
     type(vector_t) :: s11msk, s22msk, s33msk, s12msk, s13msk, s23msk
     real(kind=rp) :: center(3) = 0.0_rp
     integer :: zone_id
     character(len=20) :: zone_name 
     type(coef_t), pointer :: coef
     type(dirichlet_t) :: bc

   contains
     !> Constructor from json, wrapping the actual constructor.
     procedure, pass(this) :: init => force_torque_init_from_json
     !> Actual constructor.
     procedure, pass(this) :: init_from_attributes => &
        force_torque_init_from_attributes
     !> Destructor.
     procedure, pass(this) :: free => force_torque_free
     !> Compute the force_torque field.
     procedure, pass(this) :: compute_ => force_torque_compute
  end type force_torque_t

contains

  !> Constructor from json.
  subroutine force_torque_init_from_json(this, json, case)
    class(force_torque_t), intent(inout) :: this
    type(json_file), intent(inout) :: json
    class(case_t), intent(inout), target :: case
    integer :: zone_id
    real(kind=rp), allocatable :: center(:)
    character(len=:), allocatable :: zone_name
    real(kind=rp) :: def_center(3)

    ! Add fields keyword to the json so that the field_writer picks it up.
    ! Will also add fields to the registry.

    call this%init_base(json, case)
    def_center = 0.0
    call json_get(json, 'zone_id',zone_id)
    call json_get(json, 'zone_name',zone_name)
    call json_get(json, 'center', center)

    call force_torque_init_from_attributes(this, zone_id, zone_name, &
                                           center, case%fluid%c_xh)
  end subroutine force_torque_init_from_json

  !> Actual constructor.
  subroutine force_torque_init_from_attributes(this, zone_id, zone_name, center, coef)
    class(force_torque_t), intent(inout) :: this
    real(kind=rp) :: center(3)
    character(len=*):: zone_name
    integer :: zone_id
    type(coef_t), target :: coef
    integer :: n_pts

    this%coef => coef
    this%zone_id = zone_id
    this%center = center
    this%zone_name = zone_name

    this%u => neko_field_registry%get_field_by_name("u")
    this%v => neko_field_registry%get_field_by_name("v")
    this%w => neko_field_registry%get_field_by_name("w")
    this%p => neko_field_registry%get_field_by_name("p")

    call this%s11%init(this%u%dof) 
    call this%s22%init(this%u%dof) 
    call this%s33%init(this%u%dof) 
    call this%s12%init(this%u%dof) 
    call this%s13%init(this%u%dof) 
    call this%s23%init(this%u%dof) 
    
    call this%bc%init_base(this%coef) 
    call this%bc%mark_zone(this%case%msh%labeled_zones(this%zone_id))
    call this%bc%finalize()
    n_pts = this%bc%msk(0)
    call this%n1%init(n_pts)
    call this%n2%init(n_pts)
    call this%n3%init(n_pts)
    call this%r1%init(n_pts)
    call this%r2%init(n_pts)
    call this%r3%init(n_pts)
    call this%force1%init(n_pts)
    call this%force2%init(n_pts)
    call this%force3%init(n_pts)
    call this%force4%init(n_pts)
    call this%force5%init(n_pts)
    call this%force6%init(n_pts)
    call this%s11msk%init(n_pts)
    call this%s22msk%init(n_pts)
    call this%s33msk%init(n_pts)
    call this%s12msk%init(n_pts)
    call this%s13msk%init(n_pts)
    call this%s23msk%init(n_pts)
    call this%pmsk%init(n_pts)
    call setup_normals(this%coef, this%bc%msk, this%bc%facet, &
                       this%n1%x, this%n2%x, this%n3%x, n_pts)
    call masked_red_copy(this%r1%x,this%coef%dof%x,this%bc%msk,this%u%size(),n_pts)
    call masked_red_copy(this%r2%x,this%coef%dof%y,this%bc%msk,this%u%size(),n_pts)
    call masked_red_copy(this%r3%x,this%coef%dof%z,this%bc%msk,this%u%size(),n_pts)
    call cadd(this%r1%x,-center(1),n_pts)
    call cadd(this%r2%x,-center(2),n_pts)
    call cadd(this%r3%x,-center(3),n_pts)

  end subroutine force_torque_init_from_attributes

  !> Destructor.
  subroutine force_torque_free(this)
    class(force_torque_t), intent(inout) :: this
    call this%free_base()

    nullify(this%u)
    nullify(this%v)
    nullify(this%w)
    nullify(this%p)
    nullify(this%coef)
  end subroutine force_torque_free

  !> Compute the force_torque field.
  !! @param t The time value.
  !! @param tstep The current time-step
  subroutine force_torque_compute(this, t, tstep)
    class(force_torque_t), intent(inout) :: this
    real(kind=rp), intent(in) :: t
    integer, intent(in) :: tstep
    real(kind=rp) :: dgtq(12)
    integer :: n_pts
    n_pts = this%bc%msk(0)
    call strain_rate(this%s11%x, this%s22%x, this%s33%x, this%s12%x, this%s13%x,&
                     this%s23%x, this%u, this%v, this%w, this%coef)
    if (NEKO_BCKND_DEVICE .eq. 0) then
       call masked_red_copy(this%s11msk%x,this%s11%x,this%bc%msk,this%u%size(),n_pts)
       call masked_red_copy(this%s22msk%x,this%s22%x,this%bc%msk,this%u%size(),n_pts)
       call masked_red_copy(this%s33msk%x,this%s33%x,this%bc%msk,this%u%size(),n_pts)
       call masked_red_copy(this%s12msk%x,this%s12%x,this%bc%msk,this%u%size(),n_pts)
       call masked_red_copy(this%s13msk%x,this%s13%x,this%bc%msk,this%u%size(),n_pts)
       call masked_red_copy(this%s23msk%x,this%s23%x,this%bc%msk,this%u%size(),n_pts)
       call masked_red_copy(this%pmsk%x,this%p%x,this%bc%msk,this%u%size(),n_pts)
       call calc_force_array(this%force1%x, this%force2%x, this%force3%x,&
                             this%force4%x, this%force5%x, this%force6%x,&
                             this%s11msk%x,&
                             this%s22msk%x,&
                             this%s33msk%x,&
                             this%s12msk%x,&
                             this%s13msk%x,&
                             this%s23msk%x,&
                             this%pmsk%x,&
                             this%n1%x,&
                             this%n2%x,&
                             this%n3%x,&
                             this%case%material_properties%mu,&
                             n_pts)
                                
       dgtq(1) = glsum(this%force1%x,n_pts)
       dgtq(2) = glsum(this%force2%x,n_pts)
       dgtq(3) = glsum(this%force3%x,n_pts)
       dgtq(4) = glsum(this%force4%x,n_pts)
       dgtq(5) = glsum(this%force5%x,n_pts)
       dgtq(6) = glsum(this%force6%x,n_pts)
    else
       call device_masked_red_copy(this%s11msk%x_d,this%s11%x_d,this%bc%msk_d,this%u%size(),n_pts)
       call device_masked_red_copy(this%s22msk%x_d,this%s22%x_d,this%bc%msk_d,this%u%size(),n_pts)
       call device_masked_red_copy(this%s33msk%x_d,this%s33%x_d,this%bc%msk_d,this%u%size(),n_pts)
       call device_masked_red_copy(this%s12msk%x_d,this%s12%x_d,this%bc%msk_d,this%u%size(),n_pts)
       call device_masked_red_copy(this%s13msk%x_d,this%s13%x_d,this%bc%msk_d,this%u%size(),n_pts)
       call device_masked_red_copy(this%s23msk%x_d,this%s23%x_d,this%bc%msk_d,this%u%size(),n_pts)
       call device_masked_red_copy(this%pmsk%x_d,this%p%x_d,this%bc%msk_d,this%u%size(),n_pts)

      call device_calc_force_array(this%force1%x_d, this%force2%x_d, this%force3%x_d,&
                             this%force4%x_d, this%force5%x_d, this%force6%x_d,&
                             this%s11msk%x_d,&
                             this%s22msk%x_d,&
                             this%s33msk%x_d,&
                             this%s12msk%x_d,&
                             this%s13msk%x_d,&
                             this%s23msk%x_d,&
                             this%pmsk%x_d,&
                             this%n1%x_d,&
                             this%n2%x_d,&
                             this%n3%x_d,&
                             this%case%material_properties%mu,&
                             n_pts)    
       dgtq(1) = device_glsum(this%force1%x_d,n_pts)
       dgtq(2) = device_glsum(this%force2%x_d,n_pts)
       dgtq(3) = device_glsum(this%force3%x_d,n_pts)
       dgtq(4) = device_glsum(this%force4%x_d,n_pts)
       dgtq(5) = device_glsum(this%force5%x_d,n_pts)
       dgtq(6) = device_glsum(this%force6%x_d,n_pts)

    end if
    if (pe_rank .eq. 0) then
       write(*,*) 'Zone id', this%zone_id, this%zone_name,'calc forces and torque'
       write(*,*) tstep,dgtq(1)+dgtq(4),dgtq(1),dgtq(4),'dragx'
       write(*,*) tstep,dgtq(2)+dgtq(5),dgtq(2),dgtq(5),'dragy'
       write(*,*) tstep,dgtq(3)+dgtq(6),dgtq(3),dgtq(6),'dragz'
    end if
    call drag_torque_zone(dgtq,tstep, &
                          this%case%msh%labeled_zones(this%zone_id),&
                          this%center,&
                          this%s11%x, this%s22%x, this%s33%x, this%s12%x,&
                          this%s13%x, this%s23%x, this%p, this%coef,&
                          this%case%material_properties%mu)
     if (pe_rank .eq. 0) then
        write(*,*) 'Zone id', this%zone_id, this%zone_name,'calc forces and torque'
        write(*,*) tstep,dgtq(1)+dgtq(4),dgtq(1),dgtq(4),'dragx'
        write(*,*) tstep,dgtq(2)+dgtq(5),dgtq(2),dgtq(5),'dragy'
        write(*,*) tstep,dgtq(3)+dgtq(6),dgtq(3),dgtq(6),'dragz'
     end if
  end subroutine force_torque_compute

end module force_torque
