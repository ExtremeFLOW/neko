! Copyright (c) 2023, The Neko Authors
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
!> Interface to system information routines
module system
  use, intrinsic :: iso_c_binding
  implicit none
  private

  interface
     !> Interface to a C function to retrieve the CPU name (type).
     !! @param name Stores the retrieved name. Should be a `character`
     !! of length `len` and `c_char` kind.
     !! @param len The maximum anticipated length of the retrieved name.
     subroutine system_cpuid(name, len) &
          bind(c, name='system_cpuid')
       use, intrinsic :: iso_c_binding
       type(c_ptr), value :: name
       integer(c_int), value :: len
     end subroutine system_cpuid
  end interface

  public :: system_cpuid, system_cpu_name

 contains

   !> Retrieve the system's CPU name (type)
   !! @param name Stores the retrieved name.
   subroutine system_cpu_name(name)
     character(len=*), intent(inout) :: name
     character(kind=c_char, len=80), target :: c_name
     integer :: end_pos

     call system_cpuid(c_loc(c_name), 80)

     end_pos = scan(c_name, C_NULL_CHAR)
     if(end_pos .ge. 2) then
        name(1:end_pos-1) = c_name(1:end_pos-1)
     end if
   end subroutine system_cpu_name
  
end module system
