! Copyright (c) 2022, The Neko Authors
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
module json_case
  use neko
  use json_module
  implicit none

contains
  
  subroutine json_case_create_neko_case(neko_case, json_case)
    type(case_t), intent(inout) :: neko_case
    type(json_file) :: json_case
    type(json_core) :: jcore
    type(json_value), pointer :: json_param, p
    character(len=:), allocatable :: mesh_file
    character(len=:), allocatable :: fluid_scheme
    character(len=:), allocatable :: source_term
    character(len=:), allocatable :: initial_condition
    integer :: lx
    logical found_key
    type(file_t) :: msh_file



    call neko_log%section('Case')
    call neko_log%message('Reading case from pyNeko json file')

    ! Fill default parameters 
    call param_default(neko_case%params)
    call json_case%get_core(jcore)
    call json_case%get('parameters', json_param, found_key)
    if (found_key) then
       call json_case_create_neko_params(neko_case%params, json_param, jcore)
    else
       call neko_warning('No parameter block found')
    end if

    ! Read case description

    call json_case%get('case.mesh_file', mesh_file, found_key)
    if (.not. found_key) then
       call neko_error('Key mesh_file not found')
    else
       msh_file = file_t(trim(mesh_file))
       call msh_file%read(neko_case%msh)
    end if

    call neko_case%usr%init()
    call neko_case%usr%usr_msh_setup(neko_case%msh)

    call json_case%get('case.lx', lx, found_key)
    if (.not. found_key) then
       call neko_error('lx not defined')
    end if

    call json_case%get('case.fluid_scheme', fluid_scheme, found_key)
    if (.not. found_key) then
       call neko_error('Fluid scheme not defined')
    else       
       call fluid_scheme_factory(neko_case%fluid, trim(fluid_scheme))
       call neko_case%fluid%init(neko_case%msh, lx, neko_case%params)
    end if

    call json_case%get('case.source_term', source_term, found_key)
    if (.not. found_key) then
       call neko_error('Source term not defined')
    else
       call neko_case%fluid%set_source(trim(source_term))
    end if

    call json_case%get('case.initial_condition', initial_condition, found_key)
    if (.not. found_key) then
       call neko_error('Intiail condition not defined')
    else       
       call set_flow_ic(neko_case%fluid%u, neko_case%fluid%v, neko_case%fluid%w,&
       neko_case%fluid%p, neko_case%fluid%c_Xh, neko_case%fluid%gs_Xh,&
       initial_condition, neko_case%params)
    end if


   select type(f => neko_case%fluid)
    type is(fluid_pnpn_t)
       call f%ulag%set(f%u)
       call f%vlag%set(f%v)
       call f%wlag%set(f%w)
    end select

    call neko_case%fluid%validate()

    call neko_case%ab_bdf%set_time_order(neko_case%params%time_order)

    call neko_case%s%init(neko_case%params%nsamples, neko_case%params%T_end)
    neko_case%f_out = fluid_output_t(neko_case%fluid, &
                                     path=neko_case%params%output_dir)
    call neko_case%s%add(neko_case%f_out)

    call neko_case%q%init(neko_case%params%stats_begin)
    if (neko_case%params%stats_mean_flow) then
       call neko_case%q%add(neko_case%fluid%mean%u)
       call neko_case%q%add(neko_case%fluid%mean%v)
       call neko_case%q%add(neko_case%fluid%mean%w)
       call neko_case%q%add(neko_case%fluid%mean%p)

       if (neko_case%params%output_mean_flow) then
          neko_case%f_mf = mean_flow_output_t(neko_case%fluid%mean, &
               neko_case%params%stats_begin, path=neko_case%params%output_dir)
          call neko_case%s%add(neko_case%f_mf)
       end if
    end if

    if (neko_case%params%stats_mean_sqr_flow) then
       call neko_case%q%add(neko_case%fluid%mean_sqr%uu)
       call neko_case%q%add(neko_case%fluid%mean_sqr%vv)
       call neko_case%q%add(neko_case%fluid%mean_sqr%ww)
       call neko_case%q%add(neko_case%fluid%mean_sqr%pp)

       if (neko_case%params%output_mean_sqr_flow) then
          neko_case%f_msqrf = mean_sqr_flow_output_t(neko_case%fluid%mean_sqr, &
                                             neko_case%params%stats_begin, &
                                             path=neko_case%params%output_dir)
          call neko_case%s%add(neko_case%f_msqrf)
       end if
    end if

    call jobctrl_set_time_limit(neko_case%params%jlimit)
    
    call neko_log%end_section()

  end subroutine json_case_create_neko_case

  subroutine json_case_create_neko_params(params, json_param, jcore)
    type(param_t), intent(inout) :: params
    type(json_value), pointer :: json_param
    type(json_core)  :: jcore
    type(json_value), pointer :: p
    logical :: found_key
    real(kind=rp), allocatable :: uinf(:)
    character(len=:), allocatable :: str

    call jcore%get_child(json_param, 'nsamples', p, found_key)
    if (found_key) call jcore%get(p, params%nsamples)

    call jcore%get_child(json_param, 'output_bdry', p, found_key)
    if (found_key) call jcore%get(p, params%output_bdry)
    call jcore%get_child(json_param, 'output_part', p, found_key)
    if (found_key) call jcore%get(p, params%output_part)
    call jcore%get_child(json_param, 'output_chkp', p, found_key)
    if (found_key) call jcore%get(p, params%output_chkp)

    call jcore%get_child(json_param, 'dt', p, found_key)
    if (found_key) call jcore%get(p, params%dt)
    call jcore%get_child(json_param, 'T_end', p, found_key)
    if (found_key) call jcore%get(p, params%T_end)
    call jcore%get_child(json_param, 'rho', p, found_key)
    if (found_key) call jcore%get(p, params%rho)
    call jcore%get_child(json_param, 'mu', p, found_key)
    if (found_key) call jcore%get(p, params%mu)
    call jcore%get_child(json_param, 'Re', p, found_key)
    if (found_key) call jcore%get(p, params%Re)
    call jcore%get_child(json_param, 'uinf', p, found_key)
    if (found_key) then
       call jcore%get(p, uinf)
       params%uinf = uinf
    end if


    call jcore%get_child(json_param, 'ksp_prs.type', p, found_key)
    if (found_key) then
       call jcore%get(p, str)
       params%ksp_prs = str
       deallocate(str)
    end if
    
    call jcore%get_child(json_param, 'ksp_prs.pc', p, found_key)
    if (found_key) then
       call jcore%get(p, str)
       params%pc_prs = str
       deallocate(str)
    end if
    call jcore%get_child(json_param, 'ksp_prs.abs_tol', p, found_key)
    if (found_key) call jcore%get(p, params%abstol_prs)

    call jcore%get_child(json_param, 'ksp_vel.type', p, found_key)
    if (found_key) then
       call jcore%get(p, str)
       params%ksp_vel = str
       deallocate(str)
    end if
    
    call jcore%get_child(json_param, 'ksp_vel.pc', p, found_key)
    if (found_key) then
       call jcore%get(p, str)
       params%pc_vel = str
       deallocate(str)
    end if
    call jcore%get_child(json_param, 'ksp_vel.abs_tol', p, found_key)
    if (found_key) call jcore%get(p, params%abstol_vel)

    
    call jcore%get_child(json_param, 'fluid_inflow', p, found_key)
    if (found_key) then
       call jcore%get(p, str)
       params%fluid_inflow = str
       deallocate(str)
    end if

    call jcore%get_child(json_param, 'vol_flow_dir', p, found_key)
    if (found_key) call jcore%get(p, params%vol_flow_dir)

    call jcore%get_child(json_param, 'proj_prs_dim', p, found_key)
    if (found_key) call jcore%get(p, params%proj_prs_dim)
    call jcore%get_child(json_param, 'proj_vel_dim', p, found_key)
    if (found_key) call jcore%get(p, params%proj_vel_dim)

    call jcore%get_child(json_param, 'stats.begin', p, found_key)
    if (found_key) call jcore%get(p, params%stats_begin)
    call jcore%get_child(json_param, 'stats.mean_flow', p, found_key)
    if (found_key) call jcore%get(p, params%stats_mean_flow)
    call jcore%get_child(json_param, 'stats.output_mean_flow', p, found_key)
    if (found_key) call jcore%get(p, params%output_mean_flow)
    call jcore%get_child(json_param, 'stats.mean_sqr_flow', p, found_key)
    if (found_key) call jcore%get(p, params%stats_mean_flow)
    call jcore%get_child(json_param, 'stats.output_mean_sqr_flow', p, found_key)
    if (found_key) call jcore%get(p, params%output_mean_flow)

    call jcore%get_child(json_param, 'dealias', p, found_key)
    if (found_key) call jcore%get(p, params%dealias)
    call jcore%get_child(json_param, 'lxd', p, found_key)
    if (found_key) call jcore%get(p, params%lxd)
    
  end subroutine json_case_create_neko_params

end module json_case
