!> OpenCL JIT program library
module opencl_prgm_lib
  use opencl_intf
  use utils
  use, intrinsic :: iso_c_binding
  implicit none

#ifdef HAVE_OPENCL

  !> Device math kernels
  type(c_ptr), bind(c) :: math_program = C_NULL_PTR

contains

  subroutine opencl_prgm_lib_release

    if (c_associated(math_program)) then
       if(clReleaseProgram(math_program) .ne. CL_SUCCESS) then
          call neko_error('Failed to release program')
       end if
       math_program = C_NULL_PTR
    end if
    
  end subroutine opencl_prgm_lib_release

#endif
  
end module opencl_prgm_lib
