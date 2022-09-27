!> Interface to NVTX
module nvtx
  use iso_c_binding
  implicit none

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
#endif

contains
  
  subroutine nvtxStartRange(name)
  character(kind=c_char,len=*) :: name
  character(kind=c_char,len=256) :: trimmed_name
  integer:: i
  character :: tempName(256)
  
  trimmed_name=trim(name)//c_null_char

  ! move scalar trimmed_name into character array tempName
  do i=1,LEN(trim(name)) + 1
     tempName(i) = trimmed_name(i:i)
  enddo

  call nvtxRangePushA(tempName)

end subroutine nvtxStartRange

end module nvtx
