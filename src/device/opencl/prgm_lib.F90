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

  !> Device Derivative kernels
  type(c_ptr), bind(c) :: dudxyz_program = C_NULL_PTR

  !> Device \f$ D^T X \f$ kernels
  type(c_ptr), bind(c) :: cdtp_program = C_NULL_PTR

  !> Device onvective kernels
  type(c_ptr), bind(c) :: conv1_program = C_NULL_PTR

  !> Device Velocity gradient kernels
  type(c_ptr), bind(c) :: opgrad_program = C_NULL_PTR

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

    if (c_associated(opgrad_program)) then
       if(clReleaseProgram(opgrad_program) .ne. CL_SUCCESS) then
          call neko_error('Failed to release program')
       end if
       opgrad_program = C_NULL_PTR
    end if
    
  end subroutine opencl_prgm_lib_release

#endif
  
end module opencl_prgm_lib
