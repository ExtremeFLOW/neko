!> Interface to NVTX
!! Based on https://github.com/maxcuda/NVTX_example
module nvtx
  use, intrinsic :: iso_c_binding
  implicit none
  private

  integer, parameter :: NVTX_MAX_LEN = 256
  integer, parameter :: color(24) = [int(Z'00A6CEE3'), int(Z'001F78B4'), &
                                     int(Z'00B2DF8A'), int(Z'0033A02C'), &
                                     int(Z'00FB9A99'), int(Z'00E31A1C'), &
                                     int(Z'00FDBF6F'), int(Z'00FF7F00'), &
                                     int(Z'00CAB2D6'), int(Z'006A3D9A'), &
                                     int(Z'00FFFF99'), int(Z'00B15928'), &
                                     int(Z'008DD3C7'), int(Z'00FFFFB3'), &
                                     int(Z'00BEBADA'), int(Z'00FB8072'), &
                                     int(Z'0080B1D3'), int(Z'00FDB462'), &
                                     int(Z'00B3DE69'), int(Z'00FCCDE5'), &
                                     int(Z'00D9D9D9'), int(Z'00BC89BD'), &
                                     int(Z'00CCEBC5'), int(Z'00FFED6F')]

#ifdef HAVE_NVTX

  type, bind(c) :: nvtxEventAttributes
     integer(c_int16_t) :: version = 1
     integer(c_int16_t) :: size = 48
     integer(c_int32_t) :: category = 0
     integer(c_int32_t) :: colortype = 1
     integer(c_int32_t) :: color
     integer(c_int32_t) :: payloadtype = 0
     integer(c_int32_t) :: reserved0
     integer(c_int64_t) :: payload
     integer(c_int) :: messagetype = 1
     type(c_ptr) :: message
  end type nvtxEventAttributes

  interface nvtxRangePushA
     subroutine nvtxRangePushA(name) bind(C, name='nvtxRangePushA')
       use iso_c_binding
       character(kind=c_char) :: name(256)
     end subroutine nvtxRangePushA
  end interface nvtxRangePushA

  interface nvtxRangePushEx
     subroutine nvtxRangePushEx(event) bind(C, name='nvtxRangePushEx')
       use iso_c_binding
       import :: nvtxEventAttributes
       type(nvtxEventAttributes) :: event
     end subroutine nvtxRangePushEx
  end interface nvtxRangePushEx

  interface nvtxRangePop
     subroutine nvtxRangePop() bind(C, name='nvtxRangePop')
     end subroutine nvtxRangePop
  end interface nvtxRangePop

  public :: nvtxStartRange, nvtxRangePushA, nvtxRangePop

contains

  subroutine nvtxStartRange(name, region_id)
    character(kind=c_char,len=*) :: name
    integer, optional :: region_id
    type(nvtxEventAttributes) :: event
    character, target :: c_name(NVTX_MAX_LEN)
    integer :: i, str_len

    str_len = len(trim(name))
    do i = 1, len(trim(name))
       c_name(i) = name(i:i)
    end do
    c_name(str_len+1) = C_NULL_CHAR

    if (present(region_id)) then
       event%color = color(mod(region_id, 24) + 1)
       event%message = c_loc(c_name)
       call nvtxRangePushEx(event)
    else
       call nvtxRangePushA(c_name)
    end if

  end subroutine nvtxStartRange

#endif
end module nvtx
