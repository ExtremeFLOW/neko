module device_facet_normal
  use num_types
  use utils
  use, intrinsic :: iso_c_binding
  implicit none

#ifdef HAVE_HIP
  interface
     subroutine hip_facet_normal_apply_surfvec(msk, facet, x, y, z, u, v, w, &
                                               nx, ny, nz, area, lx, m) &
          bind(c, name='hip_facet_normal_apply_surfvec')
       use, intrinsic :: iso_c_binding
       implicit none
       integer(c_int) :: m, lx
       type(c_ptr), value  :: msk, facet, x, y, z, u, v, w, nx, ny, nz, area
     end subroutine hip_facet_normal_apply_surfvec
  end interface
#elif HAVE_CUDA
  interface
     subroutine cuda_facet_normal_apply_surfvec(msk, facet, x, y, z, u, v, w, &
                                                nx, ny, nz, area, lx, m) &
          bind(c, name='cuda_facet_normal_apply_surfvec')
       use, intrinsic :: iso_c_binding
       implicit none
       integer(c_int) :: m, lx
       type(c_ptr), value  :: msk, facet, x, y, z, u, v, w, nx, ny, nz, area
     end subroutine cuda_facet_normal_apply_surfvec
  end interface
#endif
  
contains

  subroutine device_facet_normal_apply_surfvec(msk, facet, x, y, z, u, v, w, &
                                               nx, ny, nz, area, lx, m)
    integer, intent(in) :: m, lx
    type(c_ptr) :: msk, facet, x, y, z, u, v, w, nx, ny, nz, area

#ifdef HAVE_HIP
    call hip_facet_normal_apply_surfvec(msk, facet, x, y, z, u, v, w, &
                                        nx, ny, nz, area, lx, m)
#elif HAVE_CUDA
    call cuda_facet_normal_apply_surfvec(msk, facet, x, y, z, u, v, w, &
                                         nx, ny, nz, area, lx, m)
#else
    call neko_error('No device backend configured')
#endif
    
  end subroutine device_facet_normal_apply_surfvec
  
end module device_facet_normal
