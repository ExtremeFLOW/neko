
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
  use json_utils, only : json_get, json_get_or_default
  use source_term, only : source_term_t
  use coefs, only : coef_t
  use field, only : field_t
  use utils, only : neko_error
  use tri_mesh, only : tri_mesh_t
  use point, only : point_t
  use global_interpolation, only : global_interpolation_t
  use field_registry, only : neko_field_registry
  use uset, only : uset_i4_t
  use gather_scatter
  use file
  use stack
  use math
  implicit none
  private

  !> Direct-forcing immersed boundary (DF-IBM) source term.
  type, public, extends(source_term_t) :: df_ibm_source_term_t
     type(tri_mesh_t) :: model
     type(global_interpolation_t) :: global_interp
     real(kind=rp), allocatable :: xyz(:,:)
     real(kind=rp), allocatable :: F_ib(:,:)
     type(stack_i4_t), allocatable :: neigh_el(:)
     type(field_t) :: w
     real(kind=rp), allocatable :: rmax(:)
     real(kind=rp) :: pwr_param
     type(field_list_t) :: sampled_fields
     logical :: stationary = .true.
     logical :: w_computed = .false.
     type(gs_t) :: gs
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
    real(kind=rp) :: rmax_p

    call json_get(json, "model", fname)
    call json_get_or_default(json, "rmax_percentage", rmax_p, 0.25_rp)
    call json_get_or_default(json, "power_parameter", this%pwr_param, 0.5_rp)
        
    call df_ibm_source_term_init_from_components(this, fields, &
                                                 fname, rmax_p, coef)

  end subroutine df_ibm_source_term_init_from_json

  !> The constructor from type components.
  !! @param fields A list of fields for adding the source values.
  !! @param fname STL filename for immersed boundary
  !! @param coef The SEM coeffs.
  !! @todo Throw an error unless fields is <u,v,w>
  subroutine df_ibm_source_term_init_from_components(this, fields, &
                                                     fname, rmax_p, coef) 
    class(df_ibm_source_term_t), intent(inout) :: this
    class(field_list_t), intent(inout), target :: fields
    character(len=*), intent(in) :: fname
    real(kind=rp), intent(in) :: rmax_p
    type(coef_t) :: coef
    type(file_t) :: model_file, vtk_out
    character(len=LOG_SIZE) :: log_buf
    integer :: i, j, neigh, k, pt_lid, el
    integer, pointer :: pt_neigh(:)
    type(point_t) :: pt
    type(uset_i4_t) :: neigh_el

    
    call this%free()
    call this%init_base(fields, coef)

    ! Load STL of the model representing the immersed boundary
    model_file = file_t(fname)
    call model_file%read(this%model)

    call neko_log%section("Direct-Forcing IBM")
    write(log_buf, '(A, A)') 'Model      : ', trim(fname)
    call neko_log%message(log_buf)
    if (this%model%mpts .lt. 1e1) then
       write(log_buf, '(A, I1)') ' `-Points  : ', this%model%mpts
    else if (this%model%mpts .lt. 1e2) then
       write(log_buf, '(A, I2)') ' `-Points  : ', this%model%mpts
    else if (this%model%mpts .lt. 1e3) then
       write(log_buf, '(A, I3)') ' `-Points  : ', this%model%mpts
    else if (this%model%mpts .lt. 1e4) then
       write(log_buf, '(A, I4)') ' `-Points  : ', this%model%mpts
    else if (this%model%mpts .lt. 1e5) then
       write(log_buf, '(A, I5)') ' `-Points  : ', this%model%mpts
    else if (this%model%mpts .ge. 1e6) then
       write(log_buf, '(A, I6)') ' `-Points  : ', this%model%mpts
    else if (this%model%mpts .ge. 1e7) then
       write(log_buf, '(A, I7)') ' `-Points  : ', this%model%mpts
    else if (this%model%mpts .ge. 1e8) then
       write(log_buf, '(A, I8)') ' `-Points  : ', this%model%mpts
    end if
    call neko_log%message(log_buf)

    write(log_buf, '(A,f5.2)')  'R_max %    : ', rmax_p*100.0_rp
    call neko_log%message(log_buf)

    write(log_buf, '(A,f5.2)')  'IDW Pwr    : ', this%pwr_param
    call neko_log%message(log_buf)


    
    call this%global_interp%init(coef%dof)
    call this%gs%init(coef%dof)

    allocate(this%xyz(3, this%model%mpts))
    
    do i = 1, this%model%mpts
       this%xyz(1, i) = this%model%points(i)%x(1)
       this%xyz(2, i) = this%model%points(i)%x(2)
       this%xyz(3, i) = this%model%points(i)%x(3)
    end do

    call this%global_interp%find_points_xyz(this%xyz, this%model%mpts)

    allocate(this%F_ib(this%model%mpts, size(fields%fields)))

    allocate(this%neigh_el(this%model%mpts))

    allocate(this%rmax(this%model%mpts))
    call neigh_el%init(32)
    
    do i = 1, this%model%mpts

       if (this%global_interp%proc_owner(i) .ne. pe_rank) then
          continue
       end if
       
       el = this%global_interp%el_owner(i)
       
       call this%neigh_el(i)%init()
       call this%neigh_el(i)%push(el)

       this%rmax(i) = this%coef%msh%elements(el)%e%diameter() * rmax_p

       call neigh_el%clear

       ! Add point neighbours (corners etc.)
       do j = 1, this%coef%msh%elements(el)%e%npts()
          pt = this%coef%msh%elements(el)%e%p(j)
          pt_lid = this%coef%msh%get_local(pt)
          pt_neigh => this%coef%msh%point_neigh(pt_lid)%array()
          do k = 1, this%coef%msh%point_neigh(pt_lid)%size()

             if (pt_neigh(k) .le. coef%msh%nelv .and. pt_neigh(k) .gt. 0) then
                
                call neigh_el%add(pt_neigh(k))


                this%rmax(i) = min(this%rmax(i), &
                     this%coef%msh%elements(pt_neigh(k))%e%diameter() * rmax_p)
             end if
          end do
       end do

       call neigh_el%iter_init()
       do while(neigh_el%iter_next())
          call this%neigh_el(i)%push(neigh_el%iter_value())
       end do

    end do

    call this%w%init(coef%dof, 'ib_weight')

    allocate(this%sampled_fields%fields(3))

    call neko_log%end_section()
  end subroutine df_ibm_source_term_init_from_components

  !> Destructor.
  subroutine df_ibm_source_term_free(this) 
    class(df_ibm_source_term_t), intent(inout) :: this

    call this%model%free()

    call this%global_interp%free()

    if (allocated(this%xyz)) then
       deallocate(this%xyz)
    end if

    if (allocated(this%F_ib)) then
       deallocate(this%F_ib)
    end if

    if (allocated(this%neigh_el)) then
       deallocate(this%neigh_el)
    end if

    if (allocated(this%rmax)) then
       deallocate(this%rmax)
    end if

    call this%w%free()

    call this%sampled_fields%free()
    
    this%stationary = .true.
    this%w_computed = .false.

    call this%gs%free()
    
    call this%free_base()
  end subroutine df_ibm_source_term_free

  !> Computes the source term and adds the result to `fields`.
  !! @param t The time value.
  !! @param tstep The current time-step.
  subroutine df_ibm_source_term_compute(this, t, tstep) 
    class(df_ibm_source_term_t), intent(inout) :: this
    real(kind=rp), intent(in) :: t
    integer, intent(in) :: tstep
    integer :: n_fields, i, n, e, j, k, l, lx, ei
    real(kind=rp) :: r, idw, rmax
    integer, pointer :: neigh_el(:)
    real(kind=rp) :: dt

    this%sampled_fields%fields(1)%f => neko_field_registry%get_field('u')
    this%sampled_fields%fields(2)%f => neko_field_registry%get_field('v')
    this%sampled_fields%fields(3)%f => neko_field_registry%get_field('w')

    !> @todo Change this once we have variable time-stepping
    dt = t / real(tstep, rp)

    associate( lx => this%coef%Xh%lx, global_interp => this%global_interp, &
         mpts => this%model%mpts, F_ib => this%F_ib, dm_Xh => this%coef%dof, &
         points => this%model%points, fields => this%fields%fields, w => this%w%x, &
         p => this%pwr_param)

      do i = 1, 3
         call global_interp%evaluate(F_ib(:,i), this%sampled_fields%fields(i)%f%x)
      end do

      if (.not. this%stationary .or. .not. this%w_computed) then
         w = 0.0_rp

         do i = 1, mpts

            if (global_interp%proc_owner(i) .ne. pe_rank) then
               continue
            end if
            
            neigh_el => this%neigh_el(i)%array()
            rmax = this%rmax(i)
            do ei = 1, this%neigh_el(i)%size()
               e = neigh_el(ei)
               do l = 1, lx 
                  do k = 1, lx
                     do j = 1, lx
                        r = sqrt((dm_Xh%x(j,k,l,e) - points(i)%x(1))**2 &
                               + (dm_Xh%y(j,k,l,e) - points(i)%x(2))**2 &
                               + (dm_Xh%z(j,k,l,e) - points(i)%x(3))**2)
                        
                        w(j, k, l, e) = w(j, k, l, e) &
                                      + inv_dist_weight(r, rmax, p)
                     end do
                  end do
               end do
            end do
         end do

         this%w_computed = .true.

         call this%gs%op(this%w, GS_OP_ADD)
         
      end if
   
      do i = 1, mpts
         if (global_interp%proc_owner(i) .ne. pe_rank) then
            continue
         end if
         
         neigh_el => this%neigh_el(i)%array()
         rmax = this%rmax(i)
         do ei = 1, this%neigh_el(i)%size()
            e = neigh_el(ei)
            do l = 1, lx 
               do k = 1, lx
                  do j = 1, lx
                     r = sqrt((dm_Xh%x(j,k,l,e) - points(i)%x(1))**2 &
                            + (dm_Xh%y(j,k,l,e) - points(i)%x(2))**2 &
                            + (dm_Xh%z(j,k,l,e) - points(i)%x(3))**2)
                     
                     idw = inv_dist_weight(r, rmax, p)
                     if (w(j,k,l,e) .gt. 0.0_rp) then
                        fields(1)%f%x(j,k,l,e) = fields(1)%f%x(j,k,l,e) &
                             + (-this%F_ib(i, 1) * idw) / (w(j, k, l, e) * dt)
                        fields(2)%f%x(j,k,l,e) = fields(2)%f%x(j,k,l,e) &
                             + (-this%F_ib(i, 2) * idw) / (w(j, k, l, e) * dt)
                        fields(3)%f%x(j,k,l,e) = fields(3)%f%x(j,k,l,e) &
                             + (-this%F_ib(i, 3) * idw) / (w(j, k, l, e) * dt)
                     end if
                  end do
               end do
            end do
         end do
      end do

      !> @todo To be removed once we have decided if we have a global gs op
      !! after all source terms has been computed
      do i = 1, 3
         call this%gs%op(fields(i)%f, GS_OP_ADD)
      end do
      
    end associate
       
  end subroutine df_ibm_source_term_compute

  !> Inverse distance weighting coefficient
  !! @param r Radial distance to Lagrangian point.
  !! @param rmax Radial distance for sphere of influence.
  pure function inv_dist_weight(r, rmax, p) result(idw)    
    real(kind=rp), intent(in) :: r
    real(kind=rp), intent(in) :: rmax
    real(kind=rp), intent(in) :: p
    real(kind=rp) :: idw

    if(r .ge. rmax) then
       idw = 0.0_rp
    else
       idw = ((rmax-r)/(rmax * r + NEKO_EPS))**p
    end if

  end function inv_dist_weight
  
end module df_ibm_source_term

