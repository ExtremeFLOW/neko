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
  use device, only : device_profiler_start, device_profiler_stop
  use nvtx
  use roctx
  use craypat
  use neko_profiler, only: init_profiler, free_profiler, start_region, end_region
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

    if (NEKO_ENABLE_PROFILING) then
       call init_profiler()
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

    if (NEKO_ENABLE_PROFILING) then
       call free_profiler()
    end if
  end subroutine profiler_stop

  !> Started a named (@a name) profiler region
  subroutine profiler_start_region(name, region_id)
    character(kind=c_char,len=*) :: name
    integer, optional :: region_id

#ifdef HAVE_NVTX
    if (present(region_id)) then
       call nvtxStartRange(name, region_id)
    else
       call nvtxStartRange(name)
    end if
#elif HAVE_ROCTX
    call roctxStartRange(name)
#elif CRAYPAT
    !   call craypat_region_begin(name)
#endif

    if (NEKO_ENABLE_PROFILING) then
       call start_region(name, region_id)
    end if
  end subroutine profiler_start_region

  !> End the most recently started profiler region or the one with the given
  !! name
  subroutine profiler_end_region(name, region_id)
    character(kind=c_char,len=*), optional :: name
    integer, optional :: region_id

#ifdef HAVE_NVTX
    call nvtxRangePop
#elif HAVE_ROCTX
    call roctxRangePop
#elif CRAYPAT
    !   call craypat_region_end
#endif

    if (NEKO_ENABLE_PROFILING) then
       call end_region(name, region_id)
    end if
  end subroutine profiler_end_region

end module profiler
