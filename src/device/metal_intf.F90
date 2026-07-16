! Copyright (c) 2025, The Neko Authors
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
!> Fortran interface to the Metal device layer (see device/metal/metal.m)
module metal_intf
  use utils, only : neko_error
  use, intrinsic :: iso_c_binding
  implicit none

#ifdef HAVE_METAL

  !> Success return code from Metal wrappers
  integer, parameter :: metalSuccess = 0

  interface

     !> Allocate a Metal buffer of @a s bytes
     integer(c_int) function metalAlloc(ptr_d, s) &
          bind(c, name = 'metal_alloc')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr) :: ptr_d
       integer(c_size_t), value :: s
     end function metalAlloc

     !> Map host memory to a Metal buffer of @a s bytes
     !! (zero-copy on unified memory when the host allocation allows it)
     integer(c_int) function metalMap(ptr_d, ptr_h, s) &
          bind(c, name = 'metal_map')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr) :: ptr_d
       type(c_ptr), value :: ptr_h
       integer(c_size_t), value :: s
     end function metalMap

     !> Free a Metal buffer
     integer(c_int) function metalFree(ptr_d) &
          bind(c, name = 'metal_free')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: ptr_d
     end function metalFree

     !> Copy @a s bytes from host to device
     integer(c_int) function metalMemcpyHtoD(dst_d, src, s) &
          bind(c, name = 'metal_memcpy_htod')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: dst_d, src
       integer(c_size_t), value :: s
     end function metalMemcpyHtoD

     !> Copy @a s bytes from device to host
     integer(c_int) function metalMemcpyDtoH(dst, src_d, s) &
          bind(c, name = 'metal_memcpy_dtoh')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: dst, src_d
       integer(c_size_t), value :: s
     end function metalMemcpyDtoH

     !> Copy @a s bytes device-to-device
     integer(c_int) function metalMemcpyDtoD(dst_d, src_d, s) &
          bind(c, name = 'metal_memcpy_dtod')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: dst_d, src_d
       integer(c_size_t), value :: s
     end function metalMemcpyDtoD

     !> Memset on a Metal buffer
     integer(c_int) function metalMemset(buf_d, v, s) &
          bind(c, name = 'metal_memset')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: buf_d
       integer(c_int), value :: v
       integer(c_size_t), value :: s
     end function metalMemset

     !> Synchronize the Metal device (global queue)
     integer(c_int) function metalDeviceSynchronize() &
          bind(c, name = 'metal_device_sync')
       use, intrinsic :: iso_c_binding
       implicit none
     end function metalDeviceSynchronize

     !> Synchronize a specific Metal command queue
     integer(c_int) function metalStreamSynchronize(stream) &
          bind(c, name = 'metal_stream_sync')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: stream
     end function metalStreamSynchronize

     !> Create a Metal command queue
     integer(c_int) function metalStreamCreate(stream) &
          bind(c, name = 'metal_stream_create')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr) :: stream
     end function metalStreamCreate

     !> Destroy a Metal command queue
     integer(c_int) function metalStreamDestroy(stream) &
          bind(c, name = 'metal_stream_destroy')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: stream
     end function metalStreamDestroy

     !> Create a Metal event
     integer(c_int) function metalEventCreate(event) &
          bind(c, name = 'metal_event_create')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr) :: event
     end function metalEventCreate

     !> Destroy a Metal event
     integer(c_int) function metalEventDestroy(event) &
          bind(c, name = 'metal_event_destroy')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: event
     end function metalEventDestroy

     !> Record a Metal event on a command queue
     integer(c_int) function metalEventRecord(event, stream) &
          bind(c, name = 'metal_event_record')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: event, stream
     end function metalEventRecord

     !> Wait for a Metal event
     integer(c_int) function metalEventSynchronize(event) &
          bind(c, name = 'metal_event_sync')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: event
     end function metalEventSynchronize

     !> Make a command queue wait for an event
     integer(c_int) function metalStreamWaitEvent(stream, event) &
          bind(c, name = 'metal_stream_wait_event')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: stream, event
     end function metalStreamWaitEvent

     !> Query the Metal device name (C wrapper)
     subroutine metal_device_name_c(name, len) &
          bind(c, name = 'metal_device_name')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: name
       integer(c_int), value :: len
     end subroutine metal_device_name_c

     !> Return the number of available Metal devices
     integer(c_int) function metal_device_count_c() &
          bind(c, name = 'metal_device_count')
       use, intrinsic :: iso_c_binding
       implicit none
     end function metal_device_count_c

  end interface

contains

  !> Initialise the Metal device and create the global/aux command queues
  subroutine metal_init(glb_cmd_queue, aux_cmd_queue)
    type(c_ptr), intent(inout) :: glb_cmd_queue
    type(c_ptr), intent(inout) :: aux_cmd_queue

    interface
       subroutine metal_init_c(glb_q, aux_q) &
            bind(c, name = 'metal_init')
         use, intrinsic :: iso_c_binding
         implicit none
         type(c_ptr) :: glb_q
         type(c_ptr) :: aux_q
       end subroutine metal_init_c
    end interface

    call metal_init_c(glb_cmd_queue, aux_cmd_queue)
  end subroutine metal_init

  !> Finalise the Metal device and release command queues
  subroutine metal_finalize(glb_cmd_queue, aux_cmd_queue)
    type(c_ptr), intent(inout) :: glb_cmd_queue
    type(c_ptr), intent(inout) :: aux_cmd_queue

    interface
       subroutine metal_finalize_c(glb_q, aux_q) &
            bind(c, name = 'metal_finalize')
         use, intrinsic :: iso_c_binding
         implicit none
         type(c_ptr), value :: glb_q
         type(c_ptr), value :: aux_q
       end subroutine metal_finalize_c
    end interface

    call metal_finalize_c(glb_cmd_queue, aux_cmd_queue)
  end subroutine metal_finalize

  !> Query the Metal device name
  subroutine metal_device_name(name)
    character(len=*), intent(inout) :: name
    character(kind=c_char, len=256), target :: c_name
    integer :: end_pos

    call metal_device_name_c(c_loc(c_name), 256)

    end_pos = scan(c_name, C_NULL_CHAR)
    if (end_pos .ge. 2) then
       name(1:end_pos-1) = c_name(1:end_pos-1)
    end if
  end subroutine metal_device_name

  !> Return the number of available Metal devices
  integer function metal_device_count()
    metal_device_count = metal_device_count_c()
  end function metal_device_count

#endif

end module metal_intf
