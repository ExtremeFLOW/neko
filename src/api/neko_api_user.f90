! Copyright (c) 2025, The Neko Authors
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
!> Neko API user callbacks
submodule(neko_api) neko_api_user
  implicit none

  !> Abstract interface for initial condition callbacks
  abstract interface
     subroutine api_ic_callback(scheme_name, scheme_name_len) bind(c)
       use, intrinsic :: iso_c_binding
       implicit none
       character(kind=c_char), dimension(*) :: scheme_name
       integer(c_int), value :: scheme_name_len
     end subroutine api_ic_callback
  end interface

  !> Abstract interface for boundary condition callbacks
  abstract interface
     subroutine api_bc_callback(msk, msk_size, t, tstep) bind(c)
       use, intrinsic :: iso_c_binding
       import c_rp
       implicit none
       type(c_ptr), value :: msk
       integer(c_int), value :: msk_size
       real(kind=c_rp), value :: t
       integer(c_int), value :: tstep
     end subroutine api_bc_callback
  end interface

  !> Abstract interface for callbacks requiring a field list and time
  !! Used for material properties and source terms
  abstract interface
     subroutine api_ft_callback(scheme_name, scheme_name_len, t, tstep) bind(c)
       use, intrinsic :: iso_c_binding
       import c_rp
       implicit none
       character(kind=c_char), dimension(*) :: scheme_name
       integer(c_int), value :: scheme_name_len
       real(kind=c_rp), value :: t
       integer(c_int), value :: tstep
     end subroutine api_ft_callback
  end interface

  !> Abstract interface for generic callbacks requiring only time
  !! Used for preprocess and compute callbacks
  abstract interface
     subroutine api_gn_callback(t, tstep) bind(c)
       use, intrinsic :: iso_c_binding
       import c_rp
       implicit none
       real(kind=c_rp), value :: t
       integer(c_int), value :: tstep
     end subroutine api_gn_callback
  end interface

  !> Type defining all supported callbacks via the API
  type api_user_cb
     procedure(api_ic_callback), nopass, pointer :: initial
     procedure(api_gn_callback), nopass, pointer :: preprocess
     procedure(api_gn_callback), nopass, pointer :: compute
     procedure(api_bc_callback), nopass, pointer :: dirichlet
     procedure(api_ft_callback), nopass, pointer :: material
     procedure(api_ft_callback), nopass, pointer :: source
  end type api_user_cb

  !> Registered callbacks in the API
  type(api_user_cb), allocatable :: neko_api_user_cb

  !> Pointer to an active field_list_t in a callback
  type(field_list_t), pointer :: neko_api_cb_field_list => null()

contains

  !> Register callbacks
  !! @param user User interface type
  !! @param initial_cb Initial condition callback
  !! @param preprocess_cb Pre timestep callback
  !! @param compute_cb End of timestep callback
  !! @param dirichlet_cb User boundary condition callback
  !! @param material_cb Material properties callback
  !! @param source_cb Source term callback
  module subroutine neko_api_user_cb_register(user, initial_cb, preprocess_cb, &
       compute_cb, dirichlet_cb, material_cb, source_cb)
    type(user_t), intent(inout) :: user
    type(c_funptr), value :: initial_cb, preprocess_cb, compute_cb
    type(c_funptr), value :: dirichlet_cb, material_cb, source_cb

    ! Keeping neko_api_user_cb as an alloctable is a work around for
    ! NAG which throws an incompatbile function pointer warning for
    ! single precision
    if (.not. allocated(neko_api_user_cb)) then
       allocate(neko_api_user_cb)
       neko_api_user_cb%initial => null()
       neko_api_user_cb%preprocess => null()
       neko_api_user_cb%compute => null()
       neko_api_user_cb%dirichlet => null()
       neko_api_user_cb%material => null()
       neko_api_user_cb%source => null()
    end if

    ! We need the block construct in the following if statements to
    ! adhere strictly with the f2008 standard, and support compilers
    ! not implementing TS29133 (mainly GNU Fortran with -std=f2008)
    if (c_associated(initial_cb)) then
       user%initial_conditions => neko_api_user_initial_condition
       block
         procedure(api_ic_callback), pointer :: tmp
         call c_f_procpointer(initial_cb, tmp)
         neko_api_user_cb%initial => tmp
       end block
    end if

    if (c_associated(preprocess_cb)) then
       user%preprocess => neko_api_user_preprocess
       block
         procedure(api_gn_callback), pointer :: tmp
         call c_f_procpointer(preprocess_cb, tmp)
         neko_api_user_cb%preprocess => tmp
       end block
    end if

    if (c_associated(compute_cb)) then
       user%compute => neko_api_user_compute
       block
         procedure(api_gn_callback), pointer :: tmp
         call c_f_procpointer(compute_cb, tmp)
         neko_api_user_cb%compute => tmp
       end block
    end if

    if (c_associated(dirichlet_cb)) then
       user%dirichlet_conditions => neko_api_user_dirichlet_condition
       block
         procedure(api_bc_callback), pointer :: tmp
         call c_f_procpointer(dirichlet_cb, tmp)
         neko_api_user_cb%dirichlet => tmp
       end block
    end if

    if (c_associated(material_cb)) then
       user%material_properties => neko_api_user_material_properties
       block
         procedure(api_ft_callback), pointer :: tmp
         call c_f_procpointer(material_cb, tmp)
         neko_api_user_cb%material => tmp
       end block
    end if

    if (c_associated(source_cb)) then
       user%source_term => neko_api_user_source_term
       block
         procedure(api_ft_callback), pointer :: tmp
         call c_f_procpointer(source_cb, tmp)
         neko_api_user_cb%source => tmp
       end block
    end if

  end subroutine neko_api_user_cb_register

  !> API user initial condition callback caller
  subroutine neko_api_user_initial_condition(scheme_name, fields)
    character(len=*), intent(in) :: scheme_name
    type(field_list_t), intent(inout) :: fields

    if (associated(neko_api_user_cb%initial)) then
       call neko_api_user_cb%initial(trim(scheme_name), len_trim(scheme_name))
    else
       call neko_error("Initial condition callback not defined")
    end if

  end subroutine neko_api_user_initial_condition

  !> API user preprocessing callback caller
  subroutine neko_api_user_preprocess(time)
    type(time_state_t), intent(in) :: time

    if (associated(neko_api_user_cb%preprocess)) then
       call neko_api_user_cb%preprocess(time%t, time%tstep)
    else
       call neko_error("Preprocessing callback not defined")
    end if

  end subroutine neko_api_user_preprocess

  !> API user compute callback caller
  subroutine neko_api_user_compute(time)
    type(time_state_t), intent(in) :: time

    if (associated(neko_api_user_cb%compute)) then
       call neko_api_user_cb%compute(time%t, time%tstep)
    else
       call neko_error("Compute callback not defined")
    end if

  end subroutine neko_api_user_compute

  !> API user dirichlet condition callback caller
  subroutine neko_api_user_dirichlet_condition(fields, bc, time)
    type(field_list_t), intent(inout) :: fields
    type(field_dirichlet_t), intent(in) :: bc
    type(time_state_t), intent(in) :: time
    type(c_ptr) :: bc_msk

    call neko_api_set_cb_field_list(fields)

    bc_msk = neko_api_user_bc_msk_ptr(bc)

    if (associated(neko_api_user_cb%dirichlet)) then
       call neko_api_user_cb%dirichlet(bc_msk, bc%msk(0), time%t, time%tstep)
    else
       call neko_error("Dirichlet condition callback not defined")
    end if
    nullify(neko_api_cb_field_list)

  contains

    !> Helper function to extract a pointer to the mask
    function neko_api_user_bc_msk_ptr(bc) result(bc_ptr)
      type(field_dirichlet_t), intent(in), target :: bc
      type(c_ptr) :: bc_ptr
      bc_ptr = c_loc(bc%msk(1))
    end function neko_api_user_bc_msk_ptr

  end subroutine neko_api_user_dirichlet_condition

  !> API user material properties callback caller
  subroutine neko_api_user_material_properties(scheme_name, properties, time)
    character(len=*), intent(in) :: scheme_name
    type(field_list_t), intent(inout) :: properties
    type(time_state_t), intent(in) :: time

    call neko_api_set_cb_field_list(properties)

    if (associated(neko_api_user_cb%material)) then
       call neko_api_user_cb%material(trim(scheme_name), &
            len_trim(scheme_name),time%t, time%tstep)
    else
       call neko_error("Material properties callback not defined")
    end if

    nullify(neko_api_cb_field_list)
  end subroutine neko_api_user_material_properties

  !> API user source term callback caller
  subroutine neko_api_user_source_term(scheme_name, rhs, time)
    character(len=*), intent(in) :: scheme_name
    type(field_list_t), intent(inout) :: rhs
    type(time_state_t), intent(in) :: time

    call neko_api_set_cb_field_list(rhs)

    if (associated(neko_api_user_cb%source)) then
       call neko_api_user_cb%source(trim(scheme_name), &
            len_trim(scheme_name),time%t, time%tstep)
    else
       call neko_error("Source term callback not defined")
    end if

    nullify(neko_api_cb_field_list)
  end subroutine neko_api_user_source_term

  !> Set the callbacks active field list
  subroutine neko_api_set_cb_field_list(fields)
    type(field_list_t), target, intent(inout) :: fields

    if (associated(neko_api_cb_field_list)) then
       call neko_error("Callback field list already defined")
    end if
    neko_api_cb_field_list => fields
  end subroutine neko_api_set_cb_field_list

  !> Retrive a pointer to a field for the currently active callback
  !! @param field_name Field list entry
  module function neko_api_user_cb_get_field_by_name(field_name) result(f)
    character(len=*), intent(in) :: field_name
    type(field_t), pointer :: f

    if (.not. associated(neko_api_cb_field_list)) then
       call neko_error("Callback field list not defined")
    end if

    f => neko_api_cb_field_list%get(trim(field_name))

  end function neko_api_user_cb_get_field_by_name

  !> Retrive a pointer to a field for the currently active callback
  !! @param field_idx Field index in the field list
  module function neko_api_user_cb_get_field_by_index(field_idx) result(f)
    integer, intent(in) :: field_idx
    type(field_t), pointer :: f

    if (.not. associated(neko_api_cb_field_list)) then
       call neko_error("Callback field list not defined")
    end if

    f => neko_api_cb_field_list%get(field_idx)

  end function neko_api_user_cb_get_field_by_index

end submodule neko_api_user
