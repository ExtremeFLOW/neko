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
  use scratch_registry, only : neko_scratch_registry
  use field, only : field_t
  use operators, only : curl
  use case, only : case_t
  use fld_file_output, only : fld_file_output_t
  use json_utils, only : json_get, json_get_or_default
  use field_writer, only : field_writer_t
  use coefs, only : coef_t
  use operators, only : strain_rate
  use vector, only : vector_t
  use dirichlet, only : dirichlet_t
  use drag_torque
  use logger, only : LOG_SIZE, neko_log
  use comm
  use math, only : masked_red_copy, cadd, glsum, vcross
  use device_math, only : device_masked_red_copy, device_cadd, &
                          device_glsum, device_vcross
  use device
  
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

     !Masked working arrays
     type(vector_t) :: n1, n2, n3
     type(vector_t) :: r1, r2, r3
     type(vector_t) :: force1, force2, force3
     type(vector_t) :: force4, force5, force6
     type(vector_t) :: pmsk
     type(vector_t) :: s11msk, s22msk, s33msk, s12msk, s13msk, s23msk
     real(kind=rp) :: center(3) = 0.0_rp
     real(kind=rp) :: scale
     integer :: zone_id
     character(len=20) :: zone_name 
     type(coef_t), pointer :: coef
     type(dirichlet_t) :: bc
     character(len=80) :: print_format

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
    real(kind=rp) :: scale
    logical :: long_print
    
    ! Add fields keyword to the json so that the field_writer picks it up.
    ! Will also add fields to the registry.

    call this%init_base(json, case)

    call json_get(json, 'zone_id', zone_id)
    call json_get_or_default(json, 'zone_name', zone_name, ' ')
    call json_get_or_default(json, 'scale', scale, 1.0_rp)
    call json_get_or_default(json, 'long_print', long_print, .false.)
    call json_get(json, 'center', center)
    call force_torque_init_from_attributes(this, zone_id, zone_name, &
                                           center, scale, case%fluid%c_xh, &
                                           long_print)
  end subroutine force_torque_init_from_json

  !> Actual constructor.
  subroutine force_torque_init_from_attributes(this, zone_id, zone_name, &
                                               center, scale, coef, long_print)
    class(force_torque_t), intent(inout) :: this
    real(kind=rp), intent(in) :: center(3)
    real(kind=rp), intent(in) :: scale
    character(len=*), intent(in) :: zone_name
    integer, intent(in) :: zone_id
    type(coef_t), target, intent(in) :: coef
    logical, intent(in) :: long_print
    integer :: n_pts, glb_n_pts, ierr
    character(len=1000) :: log_buf

    this%coef => coef
    this%zone_id = zone_id
    this%center = center
    this%scale = scale
    this%zone_name = zone_name
    if (long_print ) then
       this%print_format = '(I7,E20.10,E20.10,E20.10,E20.10,A)'    
    else
       this%print_format = '(I7,E13.5,E13.5,E13.5,E13.5,A)'    
    end if

    this%u => neko_field_registry%get_field_by_name("u")
    this%v => neko_field_registry%get_field_by_name("v")
    this%w => neko_field_registry%get_field_by_name("w")
    this%p => neko_field_registry%get_field_by_name("p")
    
    
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
    call masked_red_copy(this%r1%x, this%coef%dof%x, this%bc%msk, &
         this%u%size(), n_pts)
    call masked_red_copy(this%r2%x, this%coef%dof%y, this%bc%msk, &
         this%u%size(), n_pts)
    call masked_red_copy(this%r3%x, this%coef%dof%z, this%bc%msk, & 
         this%u%size(), n_pts)

    call MPI_Allreduce(n_pts, glb_n_pts, 1, &
         MPI_INTEGER, MPI_SUM, NEKO_COMM, ierr)
    ! Print some information
    call neko_log%section('Force/torque calculation')
    write(log_buf, '(A,I4,A,A)') 'Zone ', zone_id, '   ', trim(zone_name)
    call neko_log%message(log_buf)
    write(log_buf, '(A,I6, I6)') 'Global number of GLL points in zone: ', glb_n_pts
    call neko_log%message(log_buf)
    write(log_buf, '(A,E15.7,E15.7,E15.7)') 'Average of zone''s coordinates: ', &
                                            glsum(this%r1%x, n_pts)/glb_n_pts, &
                                            glsum(this%r2%x, n_pts)/glb_n_pts, &
                                            glsum(this%r3%x, n_pts)/glb_n_pts
    call neko_log%message(log_buf)
    write(log_buf, '(A,E15.7,E15.7,E15.7)') 'Center for torque calculation: ', center
    call neko_log%message(log_buf)
    write(log_buf, '(A,E15.7)') 'Scale: ', scale
    call neko_log%message(log_buf)
    call neko_log%end_section()

    call cadd(this%r1%x,-center(1), n_pts)
    call cadd(this%r2%x,-center(2), n_pts)
    call cadd(this%r3%x,-center(3), n_pts)
    if (NEKO_BCKND_DEVICE .eq. 1 .and. n_pts .gt. 0) then
       call device_memcpy(this%n1%x, this%n1%x_d, n_pts, HOST_TO_DEVICE, &
            .false.)   
       call device_memcpy(this%n2%x, this%n2%x_d, n_pts, HOST_TO_DEVICE, &
            .false.)   
       call device_memcpy(this%n3%x, this%n3%x_d, n_pts, HOST_TO_DEVICE, &
            .true.)   
       call device_memcpy(this%r1%x, this%r1%x_d, n_pts, HOST_TO_DEVICE, &
            .false.)   
       call device_memcpy(this%r2%x, this%r2%x_d, n_pts, HOST_TO_DEVICE, &
            .false.)   
       call device_memcpy(this%r3%x, this%r3%x_d, n_pts, HOST_TO_DEVICE, &
            .true.)   
    end if

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
    real(kind=rp) :: dgtq(12) = 0.0_rp
    integer :: n_pts, temp_indices(6)
    type(field_t), pointer :: s11, s22, s33, s12, s13, s23
    character(len=1000) :: log_buf

    n_pts = this%bc%msk(0)

    call neko_scratch_registry%request_field(s11, temp_indices(1))
    call neko_scratch_registry%request_field(s12, temp_indices(2))
    call neko_scratch_registry%request_field(s13, temp_indices(3))
    call neko_scratch_registry%request_field(s22, temp_indices(4))
    call neko_scratch_registry%request_field(s23, temp_indices(5))
    call neko_scratch_registry%request_field(s33, temp_indices(6))

    call strain_rate(s11%x, s22%x, s33%x, s12%x, &
         s13%x, s23%x, this%u, this%v, this%w, this%coef)

    ! On the CPU we can actually just use the original subroutines...
    if (NEKO_BCKND_DEVICE .eq. 0) then
       call masked_red_copy(this%s11msk%x, s11%x, this%bc%msk, &
            this%u%size(), n_pts)
       call masked_red_copy(this%s22msk%x, s22%x, this%bc%msk, &
            this%u%size(), n_pts)
       call masked_red_copy(this%s33msk%x, s33%x, this%bc%msk, &
            this%u%size(), n_pts)
       call masked_red_copy(this%s12msk%x, s12%x, this%bc%msk, &
            this%u%size(), n_pts)
       call masked_red_copy(this%s13msk%x, s13%x, this%bc%msk, &
            this%u%size(), n_pts)
       call masked_red_copy(this%s23msk%x, s23%x, this%bc%msk, &
            this%u%size(), n_pts)
       call masked_red_copy(this%pmsk%x, this%p%x, this%bc%msk, &
            this%u%size(), n_pts)
       call calc_force_array(this%force1%x, this%force2%x, this%force3%x, &
                             this%force4%x, this%force5%x, this%force6%x, &
                             this%s11msk%x, &
                             this%s22msk%x, &
                             this%s33msk%x, &
                             this%s12msk%x, &
                             this%s13msk%x, &
                             this%s23msk%x, &
                             this%pmsk%x, &
                             this%n1%x, &
                             this%n2%x, &
                             this%n3%x, &
                             this%case%fluid%mu, &
                             n_pts)
       dgtq(1) = glsum(this%force1%x, n_pts)
       dgtq(2) = glsum(this%force2%x, n_pts)
       dgtq(3) = glsum(this%force3%x, n_pts)
       dgtq(4) = glsum(this%force4%x, n_pts)
       dgtq(5) = glsum(this%force5%x, n_pts)
       dgtq(6) = glsum(this%force6%x, n_pts)
       !Overwriting masked s11, s22, s33 as they are no longer needed
       this%s11msk = 0.0_rp
       this%s22msk = 0.0_rp
       this%s33msk = 0.0_rp
       call vcross(this%s11msk%x, this%s22msk%x, this%s33msk%x, &
                   this%r1%x, this%r2%x, this%r3%x, &
                   this%force1%x, this%force2%x, this%force3%x, n_pts)
       
       dgtq(7) = glsum(this%s11msk%x, n_pts)
       dgtq(8) = glsum(this%s22msk%x, n_pts)
       dgtq(9) = glsum(this%s33msk%x, n_pts)
       this%s11msk = 0.0_rp
       this%s22msk = 0.0_rp
       this%s33msk = 0.0_rp
       call vcross(this%s11msk%x, this%s22msk%x, this%s33msk%x, &
                    this%r1%x, this%r2%x, this%r3%x, &
                    this%force4%x, this%force5%x, this%force6%x, n_pts)
       dgtq(10) = glsum(this%s11msk%x, n_pts)
       dgtq(11) = glsum(this%s22msk%x, n_pts)
       dgtq(12) = glsum(this%s33msk%x, n_pts)
    else
       if (n_pts .gt. 0) then
          call device_masked_red_copy(this%s11msk%x_d, s11%x_d, &
                                      this%bc%msk_d, this%u%size(), n_pts)
          call device_masked_red_copy(this%s22msk%x_d, s22%x_d, &
                                      this%bc%msk_d, this%u%size(), n_pts)
          call device_masked_red_copy(this%s33msk%x_d, s33%x_d, &
                                      this%bc%msk_d, this%u%size(), n_pts)
          call device_masked_red_copy(this%s12msk%x_d, s12%x_d, &
                                      this%bc%msk_d, this%u%size(), n_pts)
          call device_masked_red_copy(this%s13msk%x_d, s13%x_d, &
                                      this%bc%msk_d, this%u%size(), n_pts)
          call device_masked_red_copy(this%s23msk%x_d, s23%x_d, &
                                      this%bc%msk_d, this%u%size(), n_pts)
          call device_masked_red_copy(this%pmsk%x_d, this%p%x_d,  &
                                      this%bc%msk_d, this%u%size(), n_pts)

          call device_calc_force_array(this%force1%x_d, this%force2%x_d, &
                             this%force3%x_d, &
                             this%force4%x_d, &
                             this%force5%x_d, &
                             this%force6%x_d, &
                             this%s11msk%x_d, &
                             this%s22msk%x_d, &
                             this%s33msk%x_d, &
                             this%s12msk%x_d, &
                             this%s13msk%x_d, &
                             this%s23msk%x_d, &
                             this%pmsk%x_d, &
                             this%n1%x_d, &
                             this%n2%x_d, &
                             this%n3%x_d, &
                             this%case%fluid%mu, &
                             n_pts)    
          !Overwriting masked s11, s22, s33 as they are no longer needed
          call device_vcross(this%s11msk%x_d, this%s22msk%x_d, &
                             this%s33msk%x_d, &
                             this%r1%x_d, this%r2%x_d, this%r3%x_d, &
                             this%force1%x_d, this%force2%x_d, &
                             this%force3%x_d, n_pts)
          call device_vcross(this%s12msk%x_d,this%s13msk%x_d,this%s23msk%x_d, &
                      this%r1%x_d, this%r2%x_d, this%r3%x_d, &
                      this%force4%x_d, this%force5%x_d, this%force6%x_d, n_pts)
       end if   
       dgtq(1) = device_glsum(this%force1%x_d, n_pts)
       dgtq(2) = device_glsum(this%force2%x_d, n_pts)
       dgtq(3) = device_glsum(this%force3%x_d, n_pts)
       dgtq(4) = device_glsum(this%force4%x_d, n_pts)
       dgtq(5) = device_glsum(this%force5%x_d, n_pts)
       dgtq(6) = device_glsum(this%force6%x_d, n_pts)
       dgtq(7) = device_glsum(this%s11msk%x_d, n_pts)
       dgtq(8) = device_glsum(this%s22msk%x_d, n_pts)
       dgtq(9) = device_glsum(this%s33msk%x_d, n_pts)
       dgtq(10) = device_glsum(this%s12msk%x_d, n_pts)
       dgtq(11) = device_glsum(this%s13msk%x_d, n_pts)
       dgtq(12) = device_glsum(this%s23msk%x_d, n_pts)
    end if
    dgtq = this%scale*dgtq
    write(log_buf,'(A, I4, A, A)') 'Force and torque on zone ', &
          this%zone_id,'  ', this%zone_name
    call neko_log%message(log_buf)
    write(log_buf,'(A)') &
          'Time step, time, total force/torque, pressure, viscous, direction'
    call neko_log%message(log_buf)
    write(log_buf, this%print_format) &
          tstep,t,dgtq(1)+dgtq(4),dgtq(1),dgtq(4),', forcex'
    call neko_log%message(log_buf)
    write(log_buf, this%print_format) &
         tstep,t,dgtq(2)+dgtq(5),dgtq(2),dgtq(5),', forcey'
    call neko_log%message(log_buf)
    write(log_buf, this%print_format) &
         tstep,t,dgtq(3)+dgtq(6),dgtq(3),dgtq(6),', forcez'
    call neko_log%message(log_buf)
    write(log_buf, this%print_format) &
         tstep,t,dgtq(7)+dgtq(10),dgtq(7),dgtq(10),', torquex'
    call neko_log%message(log_buf)
    write(log_buf, this%print_format) &
         tstep,t,dgtq(8)+dgtq(11),dgtq(8),dgtq(11),', torquey'
    call neko_log%message(log_buf)
    write(log_buf, this%print_format) &
         tstep,t,dgtq(9)+dgtq(12),dgtq(9),dgtq(12),', torquez'
    call neko_log%message(log_buf)
    call neko_scratch_registry%relinquish_field(temp_indices)

  end subroutine force_torque_compute

end module force_torque
