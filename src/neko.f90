! Copyright (c) 2019-2024, The Neko Authors
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
!> Master module
module neko
  use num_types, only : rp, sp, dp, qp
  use comm
  use utils
  use logger
  use math, only : abscmp, rzero, izero, row_zero, rone, copy, cmult, cadd, &
       cfill, glsum, glmax, glmin, chsign, vlmax, vlmin, invcol1, invcol3, &
       invers2, vcross, vdot2, vdot3, vlsc3, vlsc2, add2, add3, add4, sub2, &
       sub3, add2s1, add2s2, addsqr2s2, cmult2, invcol2, col2, col3, subcol3, &
       add3s2, subcol4, addcol3, addcol4, ascol5, p_update, x_update, glsc2, &
       glsc3, glsc4, sort, masked_copy, cfill_mask, relcmp, glimax, glimin, &
       swap, reord, flipv, cadd2, pi, absval
  use speclib
  use dofmap, only : dofmap_t
  use space, only : space_t, GL, GLL, GJ
  use htable
  use uset
  use stack
  use tuple
  use mesh, only : mesh_t
  use point, only : point_t
  use mesh_field, only : mesh_fld_t
  use map
  use mxm_wrapper, only : mxm
  use global_interpolation
  use file
  use field, only : field_t, field_ptr_t
  use neko_mpi_types
  use gather_scatter
  use coefs, only : coef_t
  use bc
  use wall, only : no_slip_wall_t
  use dirichlet, only : dirichlet_t
  use ax_product, only : ax_t, ax_helm_factory
  use parmetis, only : parmetis_partgeom, parmetis_partmeshkway
  use neko_config
  use case, only : case_t, case_init, case_free
  use output_controller, only : output_controller_t
  use output, only : output_t
  use simulation, only : neko_solve
  use operators, only : dudxyz, opgrad, ortho, cdtp, conv1, curl, cfl,&
       lambda2op, strain_rate, div, grad
  use mathops, only : opchsign, opcolv, opcolv3c, opadd2cm, opadd2col
  use projection
  use user_intf
  use signal
  use jobctrl, only : jobctrl_init, jobctrl_set_time_limit, &
       jobctrl_time_limit, jobctrl_jobtime
  use device
  use device_math, only : device_copy, device_rzero, device_rone, &
       device_cmult, device_cmult2, device_cadd, device_cfill, device_add2, &
       device_add2s1, device_add2s2, device_addsqr2s2, device_add3s2, &
       device_invcol1, device_invcol2, device_col2, device_col3, &
       device_subcol3, device_sub2, device_sub3, device_addcol3, &
       device_addcol4, device_vdot3, device_vlsc3, device_glsc3, &
       device_glsc3_many, device_add2s2_many, device_glsc2, device_glsum, &
       device_masked_copy, device_cfill_mask, device_add3, device_cadd2, &
       device_absval
  use map_1d, only : map_1d_t
  use map_2d, only : map_2d_t
  use cpr, only : cpr_t, cpr_init, cpr_free
  use fluid_stats, only : fluid_stats_t
  use field_list, only : field_list_t
  use fluid_user_source_term
  use scalar_user_source_term
  use vector, only : vector_t, vector_ptr_t
  use matrix, only : matrix_t
  use tensor
  use simulation_component, only : simulation_component_t, &
       simulation_component_wrapper_t
  use probes, only : probes_t
  use spectral_error_indicator
  use system, only : system_cpu_name, system_cpuid
  use drag_torque, only : drag_torque_zone, drag_torque_facet, drag_torque_pt
  use field_registry, only : neko_field_registry
  use scratch_registry, only : neko_scratch_registry
  use simcomp_executor, only : neko_simcomps
  use data_streamer, only : data_streamer_t
  use time_interpolator, only : time_interpolator_t
  use point_interpolator, only : point_interpolator_t
  use point_zone, only: point_zone_t
  use box_point_zone, only: box_point_zone_t
  use sphere_point_zone, only: sphere_point_zone_t
  use point_zone_registry, only: neko_point_zone_registry
  use field_dirichlet, only : field_dirichlet_t
  use field_dirichlet_vector, only : field_dirichlet_vector_t
  use runtime_stats, only : neko_rt_stats
  use json_module, only : json_file
  use json_utils, only : json_get, json_get_or_default, json_extract_item
  use, intrinsic :: iso_fortran_env
  !$ use omp_lib
  implicit none

contains

  subroutine neko_init(C)
    type(case_t), target, intent(inout), optional :: C
    character(len=NEKO_FNAME_LEN) :: case_file
    character(len=LOG_SIZE) :: log_buf
    character(len=10) :: suffix
    character(10) :: time
    character(8) :: date
    integer :: argc, nthrds, rw, sw

    call date_and_time(time = time, date = date)

    call comm_init
    call neko_mpi_types_init
    call jobctrl_init
    call device_init

    call neko_log%init()
    call neko_field_registry%init()

    call neko_log%header(NEKO_VERSION, NEKO_BUILD_INFO)

    if (present(C)) then

       argc = command_argument_count()

       if ((argc .lt. 1) .or. (argc .gt. 1)) then
          if (pe_rank .eq. 0) write(*,*) 'Usage: ./neko < case file >'
          stop
       end if

       call get_command_argument(1, case_file)

       call filename_suffix(case_file, suffix)

       if (trim(suffix) .ne. 'case') then
          call neko_error('Invalid case file')
       end if

       ! Check the device count against the number of MPI ranks
       if (NEKO_BCKND_DEVICE .eq. 1) then
          if (device_count() .ne. 1) then
             call neko_error('Only one device is supported per MPI rank')
          end if
       end if

       !
       ! Job information
       !
       call neko_log%section("Job Information")
       write(log_buf, '(A,A,A,A,1x,A,1x,A,A,A,A,A)') 'Start time: ',&
            time(1:2), ':', time(3:4), &
            '/', date(1:4), '-', date(5:6), '-', date(7:8)
       call neko_log%message(log_buf, NEKO_LOG_QUIET)
       write(log_buf, '(a)') 'Running on: '
       sw = 10
       if (pe_size .lt. 1e1) then
          write(log_buf(13:), '(i1,a)') pe_size, ' MPI '
          if (pe_size .eq. 1) then
             write(log_buf(19:), '(a)') 'rank'
             sw = 9
          else
             write(log_buf(19:), '(a)') 'ranks'
          end if
          rw = 1
       else if (pe_size .lt. 1e2) then
          write(log_buf(13:), '(i2,a)') pe_size, ' MPI ranks'
          rw = 2
       else if (pe_size .lt. 1e3) then
          write(log_buf(13:), '(i3,a)') pe_size, ' MPI ranks'
          rw = 3
       else if (pe_size .lt. 1e4) then
          write(log_buf(13:), '(i4,a)') pe_size, ' MPI ranks'
          rw = 4
       else if (pe_size .lt. 1e5) then
          write(log_buf(13:), '(i5,a)') pe_size, ' MPI ranks'
          rw = 5
       else
          write(log_buf(13:), '(i6,a)') pe_size, ' MPI ranks'
          rw = 6
       end if

       nthrds = 1
       !$omp parallel
       !$omp master
       !$ nthrds = omp_get_num_threads()
       !$omp end master
       !$omp end parallel

       if (nthrds .gt. 1) then
          if (nthrds .lt. 1e1) then
             write(log_buf(13 + rw + sw:), '(a,i1,a)') ', using ', &
                  nthrds, ' thrds each'
          else if (nthrds .lt. 1e2) then
             write(log_buf(13 + rw + sw:), '(a,i2,a)') ', using ', &
                  nthrds, ' thrds each'
          else if (nthrds .lt. 1e3) then
             write(log_buf(13 + rw + sw:), '(a,i3,a)') ', using ', &
                  nthrds, ' thrds each'
          else if (nthrds .lt. 1e4) then
             write(log_buf(13 + rw + sw:), '(a,i4,a)') ', using ', &
                  nthrds, ' thrds each'
          end if
       end if
       call neko_log%message(log_buf, NEKO_LOG_QUIET)

       write(log_buf, '(a)') 'CPU type  : '
       call system_cpu_name(log_buf(13:))
       call neko_log%message(log_buf, NEKO_LOG_QUIET)

       write(log_buf, '(a)') 'Bcknd type: '
       if (NEKO_BCKND_SX .eq. 1) then
          write(log_buf(13:), '(a)') 'SX-Aurora'
       else if (NEKO_BCKND_XSMM .eq. 1) then
          write(log_buf(13:), '(a)') 'CPU (libxsmm)'
       else if (NEKO_BCKND_CUDA .eq. 1) then
          write(log_buf(13:), '(a)') 'Accelerator (CUDA)'
       else if (NEKO_BCKND_HIP .eq. 1) then
          write(log_buf(13:), '(a)') 'Accelerator (HIP)'
       else if (NEKO_BCKND_OPENCL .eq. 1) then
          write(log_buf(13:), '(a)') 'Accelerator (OpenCL)'
       else
          write(log_buf(13:), '(a)') 'CPU'
       end if
       call neko_log%message(log_buf, NEKO_LOG_QUIET)

       if (NEKO_BCKND_HIP .eq. 1 .or. NEKO_BCKND_CUDA .eq. 1 .or. &
            NEKO_BCKND_OPENCL .eq. 1) then
          write(log_buf, '(a)') 'Dev. name : '
          call device_name(log_buf(13:))
          call neko_log%message(log_buf, NEKO_LOG_QUIET)
       end if

       write(log_buf, '(a)') 'Real type : '
       select case (rp)
         case (real32)
          write(log_buf(13:), '(a)') 'single precision'
         case (real64)
          write(log_buf(13:), '(a)') 'double precision'
         case (real128)
          write(log_buf(13:), '(a)') 'quad precision'
       end select
       call neko_log%message(log_buf, NEKO_LOG_QUIET)

       call neko_log%end()

       !
       ! Create case
       !
       call case_init(C, case_file)

       !
       ! Setup runtime statistics
       !
       call neko_rt_stats%init(C%params)


       !
       ! Create simulation components
       !
       call neko_simcomps%init(C)

    end if

  end subroutine neko_init

  subroutine neko_finalize(C)
    type(case_t), intent(inout), optional :: C

    call neko_rt_stats%report()
    call neko_rt_stats%free()

    if (present(C)) then
       call case_free(C)
    end if

    call neko_field_registry%free()
    call neko_scratch_registry%free()
    call device_finalize
    call neko_mpi_types_free
    call comm_free
  end subroutine neko_finalize

end module neko
