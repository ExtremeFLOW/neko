!> Interface to NoaSci
module noa
  use, intrinsic :: iso_c_binding
  implicit none

  enum, bind(c)
     enumerator :: INT = 0
     enumerator :: FLOAT = 1
     enumerator :: DOUBLE = 2
  end enum

  enum, bind(c)
     enumerator :: BINARY = 0
     enumerator :: HDF5 = 1
     enumerator :: VTK  = 2
  end enum

  enum, bind(c)
     enumerator :: POSIX = 0
     enumerator :: MERO = 1
  end enum

  interface
     integer (c_int) function noa_init &
          (mero_config_filename, block_size, socket, tier) &
          bind(c, name='noa_init')
       use, intrinsic :: iso_c_binding
       implicit none
       character(kind=c_char), dimension(*) :: mero_config_filename
       integer(c_size_t) :: block_size
       integer(c_int) :: socket, tier
     end function noa_init
  end interface

  interface
     integer (c_int) function noa_finalize() &
          bind(c, name='noa_finalize')
       use, intrinsic :: iso_c_binding
       implicit none      
     end function noa_finalize
  end interface

end module noa
