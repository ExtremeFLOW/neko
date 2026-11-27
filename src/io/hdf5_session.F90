! Copyright (c) 2026, The Neko Authors
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
!> HDF5 session management
module hdf5_session
  use utils, only : neko_error
#ifdef HAVE_HDF5
  use hdf5, only : h5open_f, h5close_f
#endif
  implicit none
  private

  public :: hdf5_session_init, hdf5_session_finalize

  !> Whether a session currently holds the library open. The file backends
  !! consult this so they do not close the library out from under it.
  integer :: session_active = 0

contains

  !> Initialise the global HDF5 session
  !! HDF5 does not maintain a global reference count, so we must.
  subroutine hdf5_session_init
#ifdef HAVE_HDF5
    integer :: ierr

    if (session_active .eq. 0) then
       call h5open_f(ierr)
       if (ierr .ne. 0) call neko_error('Failed to initialize HDF5')
    end if
    session_active = session_active + 1
#endif
  end subroutine hdf5_session_init

  !> Finalize the global HDF5 session
  !! @note Must be called after all HDF5 file objects have been closed,
  !! and before MPI is finalized, since the library holds duplicated
  !! communicators for parallel I/O
  subroutine hdf5_session_finalize
#ifdef HAVE_HDF5
    integer :: ierr

    session_active = session_active - 1
    if (session_active .eq. 0) then
       call h5close_f(ierr)
       if (ierr .ne. 0) call neko_error('Failed to close HDF5')
    end if
#endif
  end subroutine hdf5_session_finalize

end module hdf5_session
