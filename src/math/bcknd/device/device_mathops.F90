module device_mathops
  use num_types
  use utils
  use comm
  use, intrinsic :: iso_c_binding
  implicit none
  
#ifdef HAVE_HIP
  interface
     subroutine hip_opchsign(a1_d, a2_d, a3_d, gdim, n) &
          bind(c, name='hip_opchsign')
       use, intrinsic :: iso_c_binding
       type(c_ptr), value :: a1_d, a2_d, a3_d
       integer(c_int) :: gdim, n
     end subroutine hip_opchsign
  end interface

  interface
     subroutine hip_opcolv(a1_d, a2_d, a3_d, c_d, gdim, n) &
          bind(c, name='hip_opcolv')
       use, intrinsic :: iso_c_binding
       type(c_ptr), value :: a1_d, a2_d, a3_d, c_d
       integer(c_int) :: gdim, n
     end subroutine hip_opcolv
  end interface

  interface
     subroutine hip_opcolv3c(a1_d, a2_d, a3_d, b1_d, b2_d, b3_d, c_d, d, gdim, n) &
          bind(c, name='hip_opcolv3c')
       use, intrinsic :: iso_c_binding
       import c_rp
       type(c_ptr), value :: a1_d, a2_d, a3_d, b1_d, b2_d, b3_d, c_d
       real(c_rp) :: d
       integer(c_int) :: gdim, n
     end subroutine hip_opcolv3c
  end interface

  interface
     subroutine hip_opadd2cm(a1_d, a2_d, a3_d, b1_d, b2_d, b3_d, c, gdim, n) &
          bind(c, name='hip_opadd2cm')
       use, intrinsic :: iso_c_binding
       import c_rp
       type(c_ptr), value :: a1_d, a2_d, a3_d, b1_d, b2_d, b3_d
       real(c_rp) :: c
       integer(c_int) :: gdim, n
     end subroutine hip_opadd2cm
  end interface

  interface
     subroutine hip_opadd2col(a1_d, a2_d, a3_d, b1_d, b2_d, b3_d, c_d, gdim, n) &
          bind(c, name='hip_opadd2col')
       use, intrinsic :: iso_c_binding
       type(c_ptr), value :: a1_d, a2_d, a3_d, b1_d, b2_d, b3_d, c_d
       integer(c_int) :: gdim, n
     end subroutine hip_opadd2col
  end interface
#endif

contains

  subroutine device_opchsign(a1_d, a2_d, a3_d, gdim, n)
    type(c_ptr) :: a1_d, a2_d, a3_d
    integer :: n, gdim
#ifdef HAVE_HIP
    call hip_opchsign(a1_d, a2_d, a3_d, gdim, n)
#else
    call neko_error('No device backend configured')
#endif
  end subroutine device_opchsign

  subroutine device_opcolv(a1_d, a2_d, a3_d, c_d, gdim, n)
    type(c_ptr) :: a1_d, a2_d, a3_d, c_d
    integer :: n, gdim
#ifdef HAVE_HIP
    call hip_opcolv(a1_d, a2_d, a3_d, c_d, gdim, n)
#else
    call neko_error('No device backend configured')
#endif
  end subroutine device_opcolv
 
  subroutine device_opcolv3c(a1_d, a2_d, a3_d, &
                             b1_d, b2_d, b3_d, c_d, d, n, gdim)
    type(c_ptr) :: a1_d, a2_d, a3_d, b1_d, b2_d, b3_d, c_d
    real(kind=rp) :: d
    integer :: n, gdim
#ifdef HAVE_HIP
    call hip_opcolv3c(a1_d, a2_d, a3_d, b1_d, b2_d, b3_d, c_d, d, gdim, n)
#else
    call neko_error('No device backend configured')
#endif
  end subroutine device_opcolv3c

  subroutine device_opadd2cm (a1_d, a2_d, a3_d, b1_d, b2_d, b3_d, c, n, gdim)
    type(c_ptr) :: a1_d, a2_d, a3_d, b1_d, b2_d, b3_d
    real(kind=rp) :: c
    integer :: n, gdim
#ifdef HAVE_HIP
    call hip_opadd2cm(a1_d, a2_d, a3_d, b1_d, b2_d, b3_d, c, gdim, n)
#else
    call neko_error('No device backend configured')
#endif
  end subroutine device_opadd2cm

  subroutine device_opadd2col (a1_d, a2_d, a3_d, b1_d, b2_d, b3_d, c_d, n, gdim)
    type(c_ptr) :: a1_d, a2_d, a3_d, b1_d, b2_d, b3_d, c_d
    integer :: n, gdim
#ifdef HAVE_HIP
    call hip_opadd2col(a1_d, a2_d, a3_d, b1_d, b2_d, b3_d, c_d, gdim, n)
#else
    call neko_error('No device backend configured')
#endif
  end subroutine device_opadd2col

end module device_mathops
