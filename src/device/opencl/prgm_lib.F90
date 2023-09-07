!> OpenCL JIT program library
module opencl_prgm_lib
  use opencl_intf
  use utils
  use, intrinsic :: iso_c_binding
  implicit none

#ifdef HAVE_OPENCL

  !> Device math kernels
  type(c_ptr), bind(c) :: math_program = C_NULL_PTR

  !> Device mathops kernels
  type(c_ptr), bind(c) :: mathops_program = C_NULL_PTR

  !> Device Dirichlet kernels
  type(c_ptr), bind(c) :: dirichlet_program = C_NULL_PTR

  !> Device Inflow kernels
  type(c_ptr), bind(c) :: inflow_program = C_NULL_PTR

  !> Device No-slip wall kernels
  type(c_ptr), bind(c) :: no_slip_wall_program = C_NULL_PTR

  !> Device Symmetry kernels
  type(c_ptr), bind(c) :: symmetry_program = C_NULL_PTR

  !> Device Facet normal kernels
  type(c_ptr), bind(c) :: facet_normal_program = C_NULL_PTR

  !> Device Blasius profile kernel
  type(c_ptr), bind(c) :: inhom_dirichlet_program = C_NULL_PTR

  !> Device Derivative kernels
  type(c_ptr), bind(c) :: dudxyz_program = C_NULL_PTR

  !> Device \f$ D^T X \f$ kernels
  type(c_ptr), bind(c) :: cdtp_program = C_NULL_PTR

  !> Device onvective kernels
  type(c_ptr), bind(c) :: conv1_program = C_NULL_PTR
  
  !> Device CFL kernels
  type(c_ptr), bind(c) :: cfl_program = C_NULL_PTR

  !> Device Velocity gradient kernels
  type(c_ptr), bind(c) :: opgrad_program = C_NULL_PTR

  !> Device Gather-Scatter kernels
  type(c_ptr), bind(c) :: gs_program = C_NULL_PTR

  !> Device Ax helm kernels
  type(c_ptr), bind(c) :: ax_helm_program = C_NULL_PTR

  !> Device jacobi kernels
  type(c_ptr), bind(c) :: jacobi_program = C_NULL_PTR

  !> Device rhs_maker kernels
  type(c_ptr), bind(c) :: rhs_maker_program = C_NULL_PTR

  !> Device pnpn residual kernels
  type(c_ptr), bind(c) :: pnpn_res_program = C_NULL_PTR

  !> Device fdm kernels
  type(c_ptr), bind(c) :: fdm_program = C_NULL_PTR

  !> Device tensor kernels
  type(c_ptr), bind(c) :: tensor_program = C_NULL_PTR

  !> Device schwarz kernels
  type(c_ptr), bind(c) :: schwarz_program = C_NULL_PTR

  !> Device dong kernels
  type(c_ptr), bind(c) :: dong_program = C_NULL_PTR

  !> Device coef kernels
  type(c_ptr), bind(c) :: coef_program = C_NULL_PTR

  !> Device scalar residual kernels
  type(c_ptr), bind(c) :: scalar_residual_program = C_NULL_PTR

  !> Device lambda2 kernels
  type(c_ptr), bind(c) :: lambda2_program = C_NULL_PTR

contains

  subroutine opencl_prgm_lib_release

    if (c_associated(math_program)) then
       if(clReleaseProgram(math_program) .ne. CL_SUCCESS) then
          call neko_error('Failed to release program')
       end if
       math_program = C_NULL_PTR
    end if

    if (c_associated(mathops_program)) then
       if(clReleaseProgram(mathops_program) .ne. CL_SUCCESS) then
          call neko_error('Failed to release program')
       end if
       mathops_program = C_NULL_PTR
    end if

    if (c_associated(dirichlet_program)) then
       if(clReleaseProgram(dirichlet_program) .ne. CL_SUCCESS) then
          call neko_error('Failed to release program')
       end if
       dirichlet_program = C_NULL_PTR
    end if

    if (c_associated(inflow_program)) then
       if(clReleaseProgram(inflow_program) .ne. CL_SUCCESS) then
          call neko_error('Failed to release program')
       end if
       inflow_program = C_NULL_PTR
    end if

    if (c_associated(no_slip_wall_program)) then
       if(clReleaseProgram(no_slip_wall_program) .ne. CL_SUCCESS) then
          call neko_error('Failed to release program')
       end if
       no_slip_wall_program = C_NULL_PTR
    end if

    if (c_associated(symmetry_program)) then
       if(clReleaseProgram(symmetry_program) .ne. CL_SUCCESS) then
          call neko_error('Failed to release program')
       end if
       symmetry_program = C_NULL_PTR
    end if

    if (c_associated(facet_normal_program)) then
       if(clReleaseProgram(facet_normal_program) .ne. CL_SUCCESS) then
          call neko_error('Failed to release program')
       end if
       facet_normal_program = C_NULL_PTR
    end if

    if (c_associated(inhom_dirichlet_program)) then
       if(clReleaseProgram(inhom_dirichlet_program) .ne. CL_SUCCESS) then
          call neko_error('Failed to release program')
       end if
       inhom_dirichlet_program = C_NULL_PTR
    end if

    if (c_associated(dudxyz_program)) then
       if(clReleaseProgram(dudxyz_program) .ne. CL_SUCCESS) then
          call neko_error('Failed to release program')
       end if
       dudxyz_program = C_NULL_PTR
    end if

    if (c_associated(cdtp_program)) then
       if(clReleaseProgram(cdtp_program) .ne. CL_SUCCESS) then
          call neko_error('Failed to release program')
       end if
       cdtp_program = C_NULL_PTR
    end if

    if (c_associated(conv1_program)) then
       if(clReleaseProgram(conv1_program) .ne. CL_SUCCESS) then
          call neko_error('Failed to release program')
       end if
       conv1_program = C_NULL_PTR
    end if

    if (c_associated(cfl_program)) then
       if(clReleaseProgram(cfl_program) .ne. CL_SUCCESS) then
          call neko_error('Failed to release program')
       end if
       cfl_program = C_NULL_PTR
    end if

    if (c_associated(opgrad_program)) then
       if(clReleaseProgram(opgrad_program) .ne. CL_SUCCESS) then
          call neko_error('Failed to release program')
       end if
       opgrad_program = C_NULL_PTR
    end if

    if (c_associated(gs_program)) then
       if(clReleaseProgram(gs_program) .ne. CL_SUCCESS) then
          call neko_error('Failed to release program')
       end if
       gs_program = C_NULL_PTR
    end if

    if (c_associated(ax_helm_program)) then
       if(clReleaseProgram(ax_helm_program) .ne. CL_SUCCESS) then
          call neko_error('Failed to release program')
       end if
       ax_helm_program = C_NULL_PTR
    end if

    if (c_associated(jacobi_program)) then
       if(clReleaseProgram(jacobi_program) .ne. CL_SUCCESS) then
          call neko_error('Failed to release program')
       end if
       jacobi_program = C_NULL_PTR
    end if

    if (c_associated(rhs_maker_program)) then
       if(clReleaseProgram(rhs_maker_program) .ne. CL_SUCCESS) then
          call neko_error('Failed to release program')
       end if
       rhs_maker_program = C_NULL_PTR
    end if

    if (c_associated(pnpn_res_program)) then
       if(clReleaseProgram(pnpn_res_program) .ne. CL_SUCCESS) then
          call neko_error('Failed to release program')
       end if
       pnpn_res_program = C_NULL_PTR
    end if

    if (c_associated(fdm_program)) then
       if(clReleaseProgram(fdm_program) .ne. CL_SUCCESS) then
          call neko_error('Failed to release program')
       end if
       fdm_program = C_NULL_PTR
    end if

    if (c_associated(tensor_program)) then
       if(clReleaseProgram(tensor_program) .ne. CL_SUCCESS) then
          call neko_error('Failed to release program')
       end if
       tensor_program = C_NULL_PTR
    end if

    if (c_associated(schwarz_program)) then
       if(clReleaseProgram(schwarz_program) .ne. CL_SUCCESS) then
          call neko_error('Failed to release program')
       end if
       schwarz_program = C_NULL_PTR
    end if

    if (c_associated(dong_program)) then
       if(clReleaseProgram(dong_program) .ne. CL_SUCCESS) then
          call neko_error('Failed to release program')
       end if
       dong_program = C_NULL_PTR
    end if

    if (c_associated(coef_program)) then
       if(clReleaseProgram(coef_program) .ne. CL_SUCCESS) then
          call neko_error('Failed to release program')
       end if
       coef_program = C_NULL_PTR
    end if

    if (c_associated(scalar_residual_program)) then
       if(clReleaseProgram(scalar_residual_program) .ne. CL_SUCCESS) then
          call neko_error('Failed to release program')
       end if
       scalar_residual_program = C_NULL_PTR
    end if
    
    if (c_associated(lambda2_program)) then
       if(clReleaseProgram(lambda2_program) .ne. CL_SUCCESS) then
          call neko_error('Failed to release program')
       end if
       lambda2_program = C_NULL_PTR
    end if
    
  end subroutine opencl_prgm_lib_release

#endif
  
end module opencl_prgm_lib
