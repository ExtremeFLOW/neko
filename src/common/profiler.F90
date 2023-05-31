! Copyright (c) 2022, The Neko Authors
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
!> Profiling interface
module profiler
  use neko_config
  use device
  use utils
  use nvtx
  use roctx
  use craypat
  implicit none

contains

  !> Start profiling
  subroutine profiler_start
    if ((NEKO_BCKND_CUDA .eq. 1)) then
#if defined(HAVE_NVTX) 
       call device_profiler_start
#endif
    else
#ifdef CRAYPAT
       call craypat_record_start
#endif
    end if
  end subroutine profiler_start

  !> Stop profiling
  subroutine profiler_stop
    if ((NEKO_BCKND_CUDA .eq. 1)) then
#if defined(HAVE_NVTX)
       call device_profiler_stop
#endif
    else
#ifdef CRAYPAT
       call craypat_record_stop
#endif
    end if
  end subroutine profiler_stop
  
  !> Started a named (@a name) profiler region
  subroutine profiler_start_region(name)
    character(kind=c_char,len=*) :: name

#ifdef HAVE_NVTX
    call nvtxStartRange(name)
#elif HAVE_ROCTX
    call roctxStartRange(name)
#elif CRAYPAT
 !   call craypat_region_begin(name)
#endif
    
  end subroutine profiler_start_region

  !> End the most recently started profiler region
  subroutine profiler_end_region

#ifdef HAVE_NVTX
    call nvtxRangePop
#elif HAVE_ROCTX
    call roctxRangePop
#elif CRAYPAT
 !   call craypat_region_end
#endif
    
  end subroutine profiler_end_region
  
end module profiler
