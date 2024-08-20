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
!> A "do nothing" filter 

module no_filter
! these were copied for lambda2
  use num_types, only : rp
  use json_module, only : json_file
  use field_registry, only : neko_field_registry
  use field, only : field_t
  use operators, only : lambda2op

! these are the ones I want
    use json_utils, only: json_get, json_get_or_default, json_extract_item
  	 use coefs, only: coef_t
    use ax_helm_fctry, only: ax_helm_factory
    use ax_product, only : ax_t
    use krylov, only : ksp_t, ksp_monitor_t
    use krylov_fctry, only : krylov_solver_factory
    use precon, only : pc_t
    use bc
    !use bc, only : bc_list_t, bc_list_init, bc_list_apply_scalar
    use neumann, only : neumann_t
    use profiler, only : profiler_start_region, profiler_end_region
    use gather_scatter, only : gs_t, GS_OP_ADD
    use pnpn_residual, only: pnpn_prs_res_t
    use mesh, only : mesh_t, NEKO_MSH_MAX_ZLBLS, NEKO_MSH_MAX_ZLBL_LEN
    use fld_file_output
    use field_registry, only : neko_field_registry
    use logger, only : neko_log, LOG_SIZE
    use filter, only: filter_t
    use scratch_registry, only: neko_scratch_registry
    use math, only: copy


    implicit none
  private

  type, public, extends(filter_t) :: no_filter_t
   contains
     !> Constructor from json.
     procedure, pass(this) :: init => no_filter_init_from_json
     !> Actual constructor.
     procedure, pass(this) :: init_from_attributes => &
          no_filter_init_from_attributes
     !> Destructor.
     procedure, pass(this) :: free => no_filter_free
     !> Compute the lambda2 field
     procedure, pass(this) :: apply => no_filter_apply
  end type no_filter_t

contains

  !> Constructor from json.
  subroutine no_filter_init_from_json(this, json, coef)
    class(no_filter_t), intent(inout) :: this
    type(json_file), intent(inout) :: json
    type(coef_t), intent(inout) :: coef
    character(len=20) :: fields(1)
    real(kind=rp) :: tmp_real



    call this%init_base(json, coef)
    call no_filter_init_from_attributes(this, coef)
   
  end subroutine no_filter_init_from_json

  !> Actual constructor.
  subroutine no_filter_init_from_attributes(this, coef)
    class(no_filter_t), intent(inout) :: this
    type(coef_t), intent(inout) :: coef
    integer :: n

  end subroutine no_filter_init_from_attributes

  !> Destructor.
  subroutine no_filter_free(this)
    class(no_filter_t), intent(inout) :: this
    call this%free_base()
  end subroutine no_filter_free

  !> Apply the filter
  !! @param F_out filtered field
  !! @param F_in unfiltered field
  subroutine no_filter_apply(this, F_out, F_in)
    class(no_filter_t), intent(inout) :: this
    type(field_t), intent(in) ::  F_in
    type(field_t), intent(inout) ::  F_out
    integer :: n, i
    n = this%coef%dof%size()
	call copy(F_out%x,F_in%x,n)
    ! I think it's a good idea to to trim everything!!!
    do i = 1, n
    	if(F_out%x(i,1,1,1).gt.1) F_out%x(i,1,1,1) = 1.0_rp
    	if(F_out%x(i,1,1,1).lt.0) F_out%x(i,1,1,1) = 0.0_rp
    end do



  end subroutine no_filter_apply

end module no_filter
