! Copyright (c) 2008-2020, UCHICAGO ARGONNE, LLC. 
!
! The UChicago Argonne, LLC as Operator of Argonne National
! Laboratory holds copyright in the Software. The copyright holder
! reserves all rights except those expressly granted to licensees,
! and U.S. Government license rights.
!
! Redistribution and use in source and binary forms, with or without
! modification, are permitted provided that the following conditions
! are met:
!
! 1. Redistributions of source code must retain the above copyright
! notice, this list of conditions and the disclaimer below.
!
! 2. Redistributions in binary form must reproduce the above copyright
! notice, this list of conditions and the disclaimer (as noted below)
! in the documentation and/or other materials provided with the
! distribution.
!
! 3. Neither the name of ANL nor the names of its contributors
! may be used to endorse or promote products derived from this software
! without specific prior written permission.
!
! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS 
! "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT 
! LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS 
! FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL 
! UCHICAGO ARGONNE, LLC, THE U.S. DEPARTMENT OF 
! ENERGY OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, 
! SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED 
! TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, 
! DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY 
! THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT 
! (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE 
! OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
!
! Additional BSD Notice
! ---------------------
! 1. This notice is required to be provided under our contract with
! the U.S. Department of Energy (DOE). This work was produced at
! Argonne National Laboratory under Contract 
! No. DE-AC02-06CH11357 with the DOE.
!
! 2. Neither the United States Government nor UCHICAGO ARGONNE, 
! LLC nor any of their employees, makes any warranty, 
! express or implied, or assumes any liability or responsibility for the
! accuracy, completeness, or usefulness of any information, apparatus,
! product, or process disclosed, or represents that its use would not
! infringe privately-owned rights.
!
! 3. Also, reference herein to any specific commercial products, process, 
! or services by trade name, trademark, manufacturer or otherwise does 
! not necessarily constitute or imply its endorsement, recommendation, 
! or favoring by the United States Government or UCHICAGO ARGONNE LLC. 
! The views and opinions of authors expressed 
! herein do not necessarily state or reflect those of the United States 
! Government or UCHICAGO ARGONNE, LLC, and shall 
! not be used for advertising or product endorsement purposes.
!

  !> Collection of vector field operations operating on \f$ a_i \f$ and \f$b_i\f$.
  !! Note that in general the indices \f$i=1 \ldots gdim\f$ and \f$j=1 \ldots n\f$.
  !! \f$gdim\f$ is assumed to be either 2 or 3 only.

module mathops
  use num_types
  implicit none

contains

  !> \f$ a_i(j) = -a_i(j) \f$ for \f$j=1 \ldots n\f$ and \f$i=1 \ldots gdim\f$.
  subroutine opchsign(a1, a2, a3, gdim, n)
    integer, intent(in) :: n, gdim
    real(kind=rp), dimension(n), intent(inout) :: a1, a2, a3
    integer :: i

    if (gdim .eq. 3) then
       do i = 1, n
          a1(i) = -a1(i)
          a2(i) = -a2(i)
          a3(i) = -a3(i)
       end do
    else
       do i = 1, n
          a1(i) = -a1(i)
          a2(i) = -a2(i)
       end do
    end if

  end subroutine opchsign
  
  !> \f$ a_i(j) = a_i(j) * c(j) \f$
  subroutine opcolv(a1, a2, a3, c, gdim, n)
    integer, intent(in) :: n, gdim
    real(kind=rp), dimension(n), intent(inout) :: a1, a2, a3
    real(kind=rp), dimension(n), intent(in) :: c
    integer :: i

    if (gdim .eq. 3) then
       do i = 1, n
          a1(i) = a1(i)*c(i)
          a2(i) = a2(i)*c(i)
          a3(i) = a3(i)*c(i)
       end do
    else
       do i = 1, n
          a1(i) = a1(i)*c(i)
          a2(i) = a2(i)*c(i)
       end do
    end if

  end subroutine opcolv

  !> \f$ a_i(j) = b_i(j) * c(j) * d \f$ 
  subroutine opcolv3c(a1, a2, a3, b1, b2, b3, c, d, n, gdim)
    integer, intent(in) :: n, gdim
    real(kind=rp), dimension(n), intent(inout) :: a1, a2, a3
    real(kind=rp), dimension(n), intent(in) :: b1, b2, b3
    real(kind=rp), intent(in) :: c(n), d
    integer :: i

    if (gdim .eq. 3) then
       do i = 1, n
          a1(i) = b1(i)*c(i)*d
          a2(i) = b2(i)*c(i)*d
          a3(i) = b3(i)*c(i)*d
       end do
    else
       do i = 1, n
          a1(i) =  b1(i)*c(i)*d
          a2(i) =  b2(i)*c(i)*d
       end do
    endif

  end subroutine opcolv3c

  !> \f$ a_i(j) = a_i(j) + b_i(j) * c \f$ 
  subroutine opadd2cm(a1, a2, a3, b1, b2, b3, c, n, gdim)
    integer, intent(in) :: n, gdim
    real(kind=rp), dimension(n), intent(inout) :: a1, a2, a3
    real(kind=rp), dimension(n), intent(in) :: b1, b2, b3
    real(kind=rp), intent(in) :: c
    integer :: i

    if (gdim .eq. 3) then
       do i = 1, n
          a1(i) = a1(i) + b1(i)*c
          a2(i) = a2(i) + b2(i)*c
          a3(i) = a3(i) + b3(i)*c
       end do
    else
       do i = 1, n
          a1(i) = a1(i) + b1(i)*c
          a2(i) = a2(i) + b2(i)*c
       end do
    endif

  end subroutine opadd2cm

  !> \f$ a_i(j) = a_i(j) + b_i(j) * c(j) \f$
  subroutine opadd2col(a1, a2, a3, b1, b2, b3, c, n, gdim)
    integer, intent(in) :: n, gdim
    real(kind=rp), dimension(n), intent(inout) :: a1, a2, a3
    real(kind=rp), dimension(n), intent(in) :: b1, b2, b3
    real(kind=rp), intent(in) :: c(n)
    integer :: i
    
    if (gdim .eq. 3) then
       do i = 1, n
          a1(i) = a1(i) + b1(i)*c(i)
          a2(i) = a2(i) + b2(i)*c(i)
          a3(i) = a3(i) + b3(i)*c(i)
       end do
    else
       do i = 1, n
          a1(i) = a1(i) + b1(i)*c(i)
          a2(i) = a2(i) + b2(i)*c(i)
       end do
    endif
    
  end subroutine opadd2col
  
end module mathops
