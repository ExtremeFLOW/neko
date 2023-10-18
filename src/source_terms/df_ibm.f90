
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
!> Implements the `df_ibm_source_term_t`
module df_ibm_source_term
  use num_types, only : rp
  use field_list, only : field_list_t
  use json_module, only : json_file
  use json_utils, only : json_get
  use source_term, only : source_term_t
  use coefs, only : coef_t
  use utils, only : neko_error
  use file
  use tri_mesh, only : tri_mesh_t
  use point, only : point_t
  use point_interpolator, only : point_interpolator_t
  use global_interpolation, only : global_interpolation_t
  implicit none
  private

  !> Direct-forcing immersed boundary (DF-IBM) source term.
  type, public, extends(source_term_t) :: df_ibm_source_term_t
     type(tri_mesh_t) :: model
     type(global_interpolation_t) :: global_interp
     type(point_interpolator_t) :: point_interp
     real(kind=rp), allocatable :: xyz(:,:)
     real(kind=rp), allocatable :: F_ib(:,:)
     real(kind=rp), allocatable :: weights_r(:,:)
     real(kind=rp), allocatable :: weights_s(:,:)
     real(kind=rp), allocatable :: weights_t(:,:)
     integer, allocatable :: local_el_owner(:)
     integer, allocatable :: local_to_global(:)
   contains
     !> The common constructor using a JSON object.
     procedure, pass(this) :: init => df_ibm_source_term_init_from_json
     !> The constructor from type components.
     procedure, pass(this) :: init_from_compenents => & 
          df_ibm_source_term_init_from_components
     !> Destructor.
     procedure, pass(this) :: free => df_ibm_source_term_free
     !> Computes the source term and adds the result to `fields`.
     procedure, pass(this) :: compute => df_ibm_source_term_compute     
  end type df_ibm_source_term_t
 
contains

  !> The common constructor using a JSON object.
  !! @param json The JSON object for the source.
  !! @param fields A list of fields for adding the source values.
  !! @param coef The SEM coeffs.
  subroutine df_ibm_source_term_init_from_json(this, json, fields, coef) 
    class(df_ibm_source_term_t), intent(inout) :: this
    type(json_file), intent(inout) :: json
    type(field_list_t), intent(inout), target :: fields
    type(coef_t), intent(inout) :: coef
    character(len=:), allocatable :: fname

    call json_get(json, "model", fname)
        
    call df_ibm_source_term_init_from_components(this, fields, fname, coef)

  end subroutine df_ibm_source_term_init_from_json

  !> The constructor from type components.
  !! @param fields A list of fields for adding the source values.
  !! @param fname STL filename for immersed boundary
  !! @param coef The SEM coeffs.
  !! @todo Throw an error unless fields is <u,v,w>
  subroutine df_ibm_source_term_init_from_components(this, fields, fname, coef) 
    class(df_ibm_source_term_t), intent(inout) :: this
    class(field_list_t), intent(inout), target :: fields
    character(len=*), intent(in) :: fname
    type(coef_t) :: coef
    type(file_t) :: model_file
    character(len=LOG_SIZE) :: log_buf
    integer :: i

    call this%free()
    call this%init_base(fields, coef)

    ! Load STL of the model representing the immersed boundary
    model_file = file_t(fname)
    call model_file%read(this%model)

    call neko_log%section("Direct-Forcing IBM")
    write(log_buf, '(A, A)') 'Model      : ', trim(fname)
    call neko_log%message(log_buf)
    if (this%model%mpts .lt. 1e1) then
       write(log_buf, '(A, I1)') '`-Points   : ', this%model%mpts
    else if (this%model%mpts .lt. 1e2) then
       write(log_buf, '(A, I2)') '`-Points   : ', this%model%mpts
    else if (this%model%mpts .lt. 1e3) then
       write(log_buf, '(A, I3)') '`-Points   : ', this%model%mpts
    else if (this%model%mpts .lt. 1e4) then
       write(log_buf, '(A, I4)') '`-Points   : ', this%model%mpts
    else if (this%model%mpts .lt. 1e5) then
       write(log_buf, '(A, I5)') '`-Points   : ', this%model%mpts
    else if (this%model%mpts .ge. 1e6) then
       write(log_buf, '(A, I6)') '`-Points   : ', this%model%mpts
    else if (this%model%mpts .ge. 1e7) then
       write(log_buf, '(A, I7)') '`-Points   : ', this%model%mpts
    else if (this%model%mpts .ge. 1e8) then
       write(log_buf, '(A, I8)') '`-Points   : ', this%model%mpts
    end if
    call neko_log%message(log_buf)
    call neko_log%end_section()


    call this%global_interp%init(coef%dof)
    call this%point_interp%init(coef%Xh)

    allocate(this%xyz(3, this%model%mpts))
    allocate(this%weights_r(coef%Xh%lx, this%model%mpts))
    allocate(this%weights_s(coef%Xh%lx, this%model%mpts))
    allocate(this%weights_t(coef%Xh%lx, this%model%mpts))
    allocate(this%F_ib(this%model%mpts, size(fields%fields)))
    
    do i = 1, this%model%mpts
       this%xyz(1, i) = this%model%points(i)%x(1)
       this%xyz(2, i) = this%model%points(i)%x(2)
       this%xyz(3, i) = this%model%points(i)%x(3)
    end do

    call this%global_interp%find_points(this%xyz, this%model%mpts)
    call this%point_interp%compute_weights(real(this%model%points(:)%x(1), rp), &
                                           real(this%model%points(:)%x(2), rp), &
                                           real(this%model%points(:)%x(3), rp), &
                                           this%weights_r, &
                                           this%weights_s, &
                                           this%weights_t)
    
    allocate(this%local_el_owner(this%model%mpts))
    allocate(this%local_to_global(this%model%mpts))

    do i = 1, this%model%mpts
       this%local_to_global(i) = i
       this%local_el_owner(i) = this%global_interp%el_owner(i)
    end do
    
  end subroutine df_ibm_source_term_init_from_components

  !> Destructor.
  subroutine df_ibm_source_term_free(this) 
    class(df_ibm_source_term_t), intent(inout) :: this

    call this%model%free()
    !call this%global_interp%free()
    !call this%point_interp%free()

    if (allocated(this%xyz)) then
       deallocate(this%xyz)
    end if

    if (allocated(this%F_ib)) then
       deallocate(this%F_ib)
    end if

    if (allocated(this%weights_r)) then
       deallocate(this%weights_r)
    end if

    if (allocated(this%weights_s)) then
       deallocate(this%weights_s)
    end if

    if (allocated(this%weights_t)) then
       deallocate(this%weights_t)
    end if

    if (allocated(this%local_el_owner)) then
       deallocate(this%local_el_owner)
    end if

    if (allocated(this%local_to_global)) then
       deallocate(this%local_to_global)
    end if   
    
    call this%free_base()
  end subroutine df_ibm_source_term_free

  !> Computes the source term and adds the result to `fields`.
  !! @param t The time value.
  !! @param tstep The current time-step.
  subroutine df_ibm_source_term_compute(this, t, tstep) 
    class(df_ibm_source_term_t), intent(inout) :: this
    real(kind=rp), intent(in) :: t
    integer, intent(in) :: tstep
    integer :: n_fields, i, n

    this%F_ib = this%point_interp%interpolate(this%model%points, &
                                              this%local_el_owner, &
                                              this%fields, &
                                              this%weights_r, &
                                              this%weights_s, &
                                              this%weights_t)
      end subroutine df_ibm_source_term_compute
  
end module df_ibm_source_term

