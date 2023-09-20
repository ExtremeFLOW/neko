!> Interface to NVTX
module nvtx
  use, intrinsic :: iso_c_binding
  implicit none
  private

  integer, parameter :: NVTX_MAX_LEN = 256

#ifdef HAVE_NVTX
  interface nvtxRangePushA
     subroutine nvtxRangePushA(name) bind(C, name='nvtxRangePushA')
       use iso_c_binding
       character(kind=c_char) :: name(256)
     end subroutine nvtxRangePushA
  end interface nvtxRangePushA
  
  interface nvtxRangePop
     subroutine nvtxRangePop() bind(C, name='nvtxRangePop')
     end subroutine nvtxRangePop
  end interface nvtxRangePop

  public : nvtxStartRange, nvtxRangePushA, nvtxRangePop
  
contains
  
  subroutine nvtxStartRange(name)
    character(kind=c_char,len=*) :: name
    character :: c_name(NVTX_MAX_LEN)
    integer:: i, str_len
    
    str_len = len(trim(name))    
    do i = 1, len(trim(name))
       c_name(i) = name(i:i)
    end do
    c_name(str_len+1) = C_NULL_CHAR
    
    call nvtxRangePushA(c_name)

  end subroutine nvtxStartRange
  
#endif
end module nvtx
