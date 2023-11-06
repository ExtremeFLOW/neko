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
c-----------------------------------------------------------------------
      subroutine mxmf2(A,N1,B,N2,C,N3)
c
c     unrolled loop version 
c
      real a(n1,n2),b(n2,n3),c(n1,n3)

      if (n2.le.8) then
         if (n2.eq.1) then
            call mxf1(a,n1,b,n2,c,n3)
         elseif (n2.eq.2) then
            call mxf2(a,n1,b,n2,c,n3)
         elseif (n2.eq.3) then
            call mxf3(a,n1,b,n2,c,n3)
         elseif (n2.eq.4) then
            call mxf4(a,n1,b,n2,c,n3)
         elseif (n2.eq.5) then
            call mxf5(a,n1,b,n2,c,n3)
         elseif (n2.eq.6) then
            call mxf6(a,n1,b,n2,c,n3)
         elseif (n2.eq.7) then
            call mxf7(a,n1,b,n2,c,n3)
         else
            call mxf8(a,n1,b,n2,c,n3)
         endif
      elseif (n2.le.16) then
         if (n2.eq.9) then
            call mxf9(a,n1,b,n2,c,n3)
         elseif (n2.eq.10) then
            call mxf10(a,n1,b,n2,c,n3)
         elseif (n2.eq.11) then
            call mxf11(a,n1,b,n2,c,n3)
         elseif (n2.eq.12) then
            call mxf12(a,n1,b,n2,c,n3)
         elseif (n2.eq.13) then
            call mxf13(a,n1,b,n2,c,n3)
         elseif (n2.eq.14) then
            call mxf14(a,n1,b,n2,c,n3)
         elseif (n2.eq.15) then
            call mxf15(a,n1,b,n2,c,n3)
         else
            call mxf16(a,n1,b,n2,c,n3)
         endif
      elseif (n2.le.24) then
         if (n2.eq.17) then
            call mxf17(a,n1,b,n2,c,n3)
         elseif (n2.eq.18) then
            call mxf18(a,n1,b,n2,c,n3)
         elseif (n2.eq.19) then
            call mxf19(a,n1,b,n2,c,n3)
         elseif (n2.eq.20) then
            call mxf20(a,n1,b,n2,c,n3)
         elseif (n2.eq.21) then
            call mxf21(a,n1,b,n2,c,n3)
         elseif (n2.eq.22) then
            call mxf22(a,n1,b,n2,c,n3)
         elseif (n2.eq.23) then
            call mxf23(a,n1,b,n2,c,n3)
         elseif (n2.eq.24) then
            call mxf24(a,n1,b,n2,c,n3)
         endif
      else
         call mxm44_0(a,n1,b,n2,c,n3)
      endif
c
      return
      end
c-----------------------------------------------------------------------
      subroutine mxf1(a,n1,b,n2,c,n3)
c
      real a(n1,1),b(1,n3),c(n1,n3)
c
      do j=1,n3
         do i=1,n1
            c(i,j) = a(i,1)*b(1,j)
         enddo
      enddo
      return
      end
c-----------------------------------------------------------------------
      subroutine mxf2(a,n1,b,n2,c,n3)
c
      real a(n1,2),b(2,n3),c(n1,n3)
c
      do j=1,n3
         do i=1,n1
            c(i,j) = a(i,1)*b(1,j)
     $             + a(i,2)*b(2,j)
         enddo
      enddo
      return
      end
c-----------------------------------------------------------------------
      subroutine mxf3(a,n1,b,n2,c,n3)
c
      real a(n1,3),b(3,n3),c(n1,n3)
c
      do j=1,n3
         do i=1,n1
            c(i,j) = a(i,1)*b(1,j)
     $             + a(i,2)*b(2,j)
     $             + a(i,3)*b(3,j)
         enddo
      enddo
      return
      end
c-----------------------------------------------------------------------
      subroutine mxf4(a,n1,b,n2,c,n3)
c
      real a(n1,4),b(4,n3),c(n1,n3)
c
      do j=1,n3
         do i=1,n1
            c(i,j) = a(i,1)*b(1,j)
     $             + a(i,2)*b(2,j)
     $             + a(i,3)*b(3,j)
     $             + a(i,4)*b(4,j)
         enddo
      enddo
      return
      end
c-----------------------------------------------------------------------
      subroutine mxf5(a,n1,b,n2,c,n3)
c
      real a(n1,5),b(5,n3),c(n1,n3)
c
      do j=1,n3
         do i=1,n1
            c(i,j) = a(i,1)*b(1,j)
     $             + a(i,2)*b(2,j)
     $             + a(i,3)*b(3,j)
     $             + a(i,4)*b(4,j)
     $             + a(i,5)*b(5,j)
         enddo
      enddo
      return
      end
c-----------------------------------------------------------------------
      subroutine mxf6(a,n1,b,n2,c,n3)
c
      real a(n1,6),b(6,n3),c(n1,n3)
c
      do j=1,n3
         do i=1,n1
            c(i,j) = a(i,1)*b(1,j)
     $             + a(i,2)*b(2,j)
     $             + a(i,3)*b(3,j)
     $             + a(i,4)*b(4,j)
     $             + a(i,5)*b(5,j)
     $             + a(i,6)*b(6,j)
         enddo
      enddo
      return
      end
c-----------------------------------------------------------------------
      subroutine mxf7(a,n1,b,n2,c,n3)
c
      real a(n1,7),b(7,n3),c(n1,n3)
c
      do j=1,n3
         do i=1,n1
            c(i,j) = a(i,1)*b(1,j)
     $             + a(i,2)*b(2,j)
     $             + a(i,3)*b(3,j)
     $             + a(i,4)*b(4,j)
     $             + a(i,5)*b(5,j)
     $             + a(i,6)*b(6,j)
     $             + a(i,7)*b(7,j)
         enddo
      enddo
      return
      end
c-----------------------------------------------------------------------
      subroutine mxf8(a,n1,b,n2,c,n3)
c
      real a(n1,8),b(8,n3),c(n1,n3)
c
      do j=1,n3
         do i=1,n1
            c(i,j) = a(i,1)*b(1,j)
     $             + a(i,2)*b(2,j)
     $             + a(i,3)*b(3,j)
     $             + a(i,4)*b(4,j)
     $             + a(i,5)*b(5,j)
     $             + a(i,6)*b(6,j)
     $             + a(i,7)*b(7,j)
     $             + a(i,8)*b(8,j)
         enddo
      enddo
      return
      end
c-----------------------------------------------------------------------
      subroutine mxf9(a,n1,b,n2,c,n3)
c
      real a(n1,9),b(9,n3),c(n1,n3)
c
      do j=1,n3
         do i=1,n1
            c(i,j) = a(i,1)*b(1,j)
     $             + a(i,2)*b(2,j)
     $             + a(i,3)*b(3,j)
     $             + a(i,4)*b(4,j)
     $             + a(i,5)*b(5,j)
     $             + a(i,6)*b(6,j)
     $             + a(i,7)*b(7,j)
     $             + a(i,8)*b(8,j)
     $             + a(i,9)*b(9,j)
         enddo
      enddo
      return
      end
c-----------------------------------------------------------------------
      subroutine mxf10(a,n1,b,n2,c,n3)
c
      real a(n1,10),b(10,n3),c(n1,n3)
c
      do j=1,n3
         do i=1,n1
            c(i,j) = a(i,1)*b(1,j)
     $             + a(i,2)*b(2,j)
     $             + a(i,3)*b(3,j)
     $             + a(i,4)*b(4,j)
     $             + a(i,5)*b(5,j)
     $             + a(i,6)*b(6,j)
     $             + a(i,7)*b(7,j)
     $             + a(i,8)*b(8,j)
     $             + a(i,9)*b(9,j)
     $             + a(i,10)*b(10,j)
         enddo
      enddo
      return
      end
c-----------------------------------------------------------------------
      subroutine mxf11(a,n1,b,n2,c,n3)
c
      real a(n1,11),b(11,n3),c(n1,n3)
c
      do j=1,n3
         do i=1,n1
            c(i,j) = a(i,1)*b(1,j)
     $             + a(i,2)*b(2,j)
     $             + a(i,3)*b(3,j)
     $             + a(i,4)*b(4,j)
     $             + a(i,5)*b(5,j)
     $             + a(i,6)*b(6,j)
     $             + a(i,7)*b(7,j)
     $             + a(i,8)*b(8,j)
     $             + a(i,9)*b(9,j)
     $             + a(i,10)*b(10,j)
     $             + a(i,11)*b(11,j)
         enddo
      enddo
      return
      end
c-----------------------------------------------------------------------
      subroutine mxf12(a,n1,b,n2,c,n3)
c
      real a(n1,12),b(12,n3),c(n1,n3)
c
      do j=1,n3
         do i=1,n1
            c(i,j) = a(i,1)*b(1,j)
     $             + a(i,2)*b(2,j)
     $             + a(i,3)*b(3,j)
     $             + a(i,4)*b(4,j)
     $             + a(i,5)*b(5,j)
     $             + a(i,6)*b(6,j)
     $             + a(i,7)*b(7,j)
     $             + a(i,8)*b(8,j)
     $             + a(i,9)*b(9,j)
     $             + a(i,10)*b(10,j)
     $             + a(i,11)*b(11,j)
     $             + a(i,12)*b(12,j)
         enddo
      enddo
      return
      end
c-----------------------------------------------------------------------
      subroutine mxf13(a,n1,b,n2,c,n3)
c
      real a(n1,13),b(13,n3),c(n1,n3)
c
      do j=1,n3
         do i=1,n1
            c(i,j) = a(i,1)*b(1,j)
     $             + a(i,2)*b(2,j)
     $             + a(i,3)*b(3,j)
     $             + a(i,4)*b(4,j)
     $             + a(i,5)*b(5,j)
     $             + a(i,6)*b(6,j)
     $             + a(i,7)*b(7,j)
     $             + a(i,8)*b(8,j)
     $             + a(i,9)*b(9,j)
     $             + a(i,10)*b(10,j)
     $             + a(i,11)*b(11,j)
     $             + a(i,12)*b(12,j)
     $             + a(i,13)*b(13,j)
         enddo
      enddo
      return
      end
c-----------------------------------------------------------------------
      subroutine mxf14(a,n1,b,n2,c,n3)
c
      real a(n1,14),b(14,n3),c(n1,n3)
c
      do j=1,n3
         do i=1,n1
            c(i,j) = a(i,1)*b(1,j)
     $             + a(i,2)*b(2,j)
     $             + a(i,3)*b(3,j)
     $             + a(i,4)*b(4,j)
     $             + a(i,5)*b(5,j)
     $             + a(i,6)*b(6,j)
     $             + a(i,7)*b(7,j)
     $             + a(i,8)*b(8,j)
     $             + a(i,9)*b(9,j)
     $             + a(i,10)*b(10,j)
     $             + a(i,11)*b(11,j)
     $             + a(i,12)*b(12,j)
     $             + a(i,13)*b(13,j)
     $             + a(i,14)*b(14,j)
         enddo
      enddo
      return
      end
c-----------------------------------------------------------------------
      subroutine mxf15(a,n1,b,n2,c,n3)
c
      real a(n1,15),b(15,n3),c(n1,n3)
c
      do j=1,n3
         do i=1,n1
            c(i,j) = a(i,1)*b(1,j)
     $             + a(i,2)*b(2,j)
     $             + a(i,3)*b(3,j)
     $             + a(i,4)*b(4,j)
     $             + a(i,5)*b(5,j)
     $             + a(i,6)*b(6,j)
     $             + a(i,7)*b(7,j)
     $             + a(i,8)*b(8,j)
     $             + a(i,9)*b(9,j)
     $             + a(i,10)*b(10,j)
     $             + a(i,11)*b(11,j)
     $             + a(i,12)*b(12,j)
     $             + a(i,13)*b(13,j)
     $             + a(i,14)*b(14,j)
     $             + a(i,15)*b(15,j)
         enddo
      enddo
      return
      end
c-----------------------------------------------------------------------
      subroutine mxf16(a,n1,b,n2,c,n3)
c
      real a(n1,16),b(16,n3),c(n1,n3)
c
      do j=1,n3
         do i=1,n1
            c(i,j) = a(i,1)*b(1,j)
     $             + a(i,2)*b(2,j)
     $             + a(i,3)*b(3,j)
     $             + a(i,4)*b(4,j)
     $             + a(i,5)*b(5,j)
     $             + a(i,6)*b(6,j)
     $             + a(i,7)*b(7,j)
     $             + a(i,8)*b(8,j)
     $             + a(i,9)*b(9,j)
     $             + a(i,10)*b(10,j)
     $             + a(i,11)*b(11,j)
     $             + a(i,12)*b(12,j)
     $             + a(i,13)*b(13,j)
     $             + a(i,14)*b(14,j)
     $             + a(i,15)*b(15,j)
     $             + a(i,16)*b(16,j)
         enddo
      enddo
      return
      end
c-----------------------------------------------------------------------
      subroutine mxf17(a,n1,b,n2,c,n3)
c
      real a(n1,17),b(17,n3),c(n1,n3)
c
      do j=1,n3
         do i=1,n1
            c(i,j) = a(i,1)*b(1,j)
     $             + a(i,2)*b(2,j)
     $             + a(i,3)*b(3,j)
     $             + a(i,4)*b(4,j)
     $             + a(i,5)*b(5,j)
     $             + a(i,6)*b(6,j)
     $             + a(i,7)*b(7,j)
     $             + a(i,8)*b(8,j)
     $             + a(i,9)*b(9,j)
     $             + a(i,10)*b(10,j)
     $             + a(i,11)*b(11,j)
     $             + a(i,12)*b(12,j)
     $             + a(i,13)*b(13,j)
     $             + a(i,14)*b(14,j)
     $             + a(i,15)*b(15,j)
     $             + a(i,16)*b(16,j)
     $             + a(i,17)*b(17,j)
         enddo
      enddo
      return
      end
c-----------------------------------------------------------------------
      subroutine mxf18(a,n1,b,n2,c,n3)
c
      real a(n1,18),b(18,n3),c(n1,n3)
c
      do j=1,n3
         do i=1,n1
            c(i,j) = a(i,1)*b(1,j)
     $             + a(i,2)*b(2,j)
     $             + a(i,3)*b(3,j)
     $             + a(i,4)*b(4,j)
     $             + a(i,5)*b(5,j)
     $             + a(i,6)*b(6,j)
     $             + a(i,7)*b(7,j)
     $             + a(i,8)*b(8,j)
     $             + a(i,9)*b(9,j)
     $             + a(i,10)*b(10,j)
     $             + a(i,11)*b(11,j)
     $             + a(i,12)*b(12,j)
     $             + a(i,13)*b(13,j)
     $             + a(i,14)*b(14,j)
     $             + a(i,15)*b(15,j)
     $             + a(i,16)*b(16,j)
     $             + a(i,17)*b(17,j)
     $             + a(i,18)*b(18,j)
         enddo
      enddo
      return
      end
c-----------------------------------------------------------------------
      subroutine mxf19(a,n1,b,n2,c,n3)
c
      real a(n1,19),b(19,n3),c(n1,n3)
c
      do j=1,n3
         do i=1,n1
            c(i,j) = a(i,1)*b(1,j)
     $             + a(i,2)*b(2,j)
     $             + a(i,3)*b(3,j)
     $             + a(i,4)*b(4,j)
     $             + a(i,5)*b(5,j)
     $             + a(i,6)*b(6,j)
     $             + a(i,7)*b(7,j)
     $             + a(i,8)*b(8,j)
     $             + a(i,9)*b(9,j)
     $             + a(i,10)*b(10,j)
     $             + a(i,11)*b(11,j)
     $             + a(i,12)*b(12,j)
     $             + a(i,13)*b(13,j)
     $             + a(i,14)*b(14,j)
     $             + a(i,15)*b(15,j)
     $             + a(i,16)*b(16,j)
     $             + a(i,17)*b(17,j)
     $             + a(i,18)*b(18,j)
     $             + a(i,19)*b(19,j)
         enddo
      enddo
      return
      end
c-----------------------------------------------------------------------
      subroutine mxf20(a,n1,b,n2,c,n3)
c
      real a(n1,20),b(20,n3),c(n1,n3)
c
      do j=1,n3
         do i=1,n1
            c(i,j) = a(i,1)*b(1,j)
     $             + a(i,2)*b(2,j)
     $             + a(i,3)*b(3,j)
     $             + a(i,4)*b(4,j)
     $             + a(i,5)*b(5,j)
     $             + a(i,6)*b(6,j)
     $             + a(i,7)*b(7,j)
     $             + a(i,8)*b(8,j)
     $             + a(i,9)*b(9,j)
     $             + a(i,10)*b(10,j)
     $             + a(i,11)*b(11,j)
     $             + a(i,12)*b(12,j)
     $             + a(i,13)*b(13,j)
     $             + a(i,14)*b(14,j)
     $             + a(i,15)*b(15,j)
     $             + a(i,16)*b(16,j)
     $             + a(i,17)*b(17,j)
     $             + a(i,18)*b(18,j)
     $             + a(i,19)*b(19,j)
     $             + a(i,20)*b(20,j)
         enddo
      enddo
      return
      end
c-----------------------------------------------------------------------
      subroutine mxf21(a,n1,b,n2,c,n3)
c
      real a(n1,21),b(21,n3),c(n1,n3)
c
      do j=1,n3
         do i=1,n1
            c(i,j) = a(i,1)*b(1,j)
     $             + a(i,2)*b(2,j)
     $             + a(i,3)*b(3,j)
     $             + a(i,4)*b(4,j)
     $             + a(i,5)*b(5,j)
     $             + a(i,6)*b(6,j)
     $             + a(i,7)*b(7,j)
     $             + a(i,8)*b(8,j)
     $             + a(i,9)*b(9,j)
     $             + a(i,10)*b(10,j)
     $             + a(i,11)*b(11,j)
     $             + a(i,12)*b(12,j)
     $             + a(i,13)*b(13,j)
     $             + a(i,14)*b(14,j)
     $             + a(i,15)*b(15,j)
     $             + a(i,16)*b(16,j)
     $             + a(i,17)*b(17,j)
     $             + a(i,18)*b(18,j)
     $             + a(i,19)*b(19,j)
     $             + a(i,20)*b(20,j)
     $             + a(i,21)*b(21,j)
         enddo
      enddo
      return
      end
c-----------------------------------------------------------------------
      subroutine mxf22(a,n1,b,n2,c,n3)
c
      real a(n1,22),b(22,n3),c(n1,n3)
c
      do j=1,n3
         do i=1,n1
            c(i,j) = a(i,1)*b(1,j)
     $             + a(i,2)*b(2,j)
     $             + a(i,3)*b(3,j)
     $             + a(i,4)*b(4,j)
     $             + a(i,5)*b(5,j)
     $             + a(i,6)*b(6,j)
     $             + a(i,7)*b(7,j)
     $             + a(i,8)*b(8,j)
     $             + a(i,9)*b(9,j)
     $             + a(i,10)*b(10,j)
     $             + a(i,11)*b(11,j)
     $             + a(i,12)*b(12,j)
     $             + a(i,13)*b(13,j)
     $             + a(i,14)*b(14,j)
     $             + a(i,15)*b(15,j)
     $             + a(i,16)*b(16,j)
     $             + a(i,17)*b(17,j)
     $             + a(i,18)*b(18,j)
     $             + a(i,19)*b(19,j)
     $             + a(i,20)*b(20,j)
     $             + a(i,21)*b(21,j)
     $             + a(i,22)*b(22,j)
         enddo
      enddo
      return
      end
c-----------------------------------------------------------------------
      subroutine mxf23(a,n1,b,n2,c,n3)
c
      real a(n1,23),b(23,n3),c(n1,n3)
c
      do j=1,n3
         do i=1,n1
            c(i,j) = a(i,1)*b(1,j)
     $             + a(i,2)*b(2,j)
     $             + a(i,3)*b(3,j)
     $             + a(i,4)*b(4,j)
     $             + a(i,5)*b(5,j)
     $             + a(i,6)*b(6,j)
     $             + a(i,7)*b(7,j)
     $             + a(i,8)*b(8,j)
     $             + a(i,9)*b(9,j)
     $             + a(i,10)*b(10,j)
     $             + a(i,11)*b(11,j)
     $             + a(i,12)*b(12,j)
     $             + a(i,13)*b(13,j)
     $             + a(i,14)*b(14,j)
     $             + a(i,15)*b(15,j)
     $             + a(i,16)*b(16,j)
     $             + a(i,17)*b(17,j)
     $             + a(i,18)*b(18,j)
     $             + a(i,19)*b(19,j)
     $             + a(i,20)*b(20,j)
     $             + a(i,21)*b(21,j)
     $             + a(i,22)*b(22,j)
     $             + a(i,23)*b(23,j)
         enddo
      enddo
      return
      end
c-----------------------------------------------------------------------
      subroutine mxf24(a,n1,b,n2,c,n3)
c
      real a(n1,24),b(24,n3),c(n1,n3)
c
      do j=1,n3
         do i=1,n1
            c(i,j) = a(i,1)*b(1,j)
     $             + a(i,2)*b(2,j)
     $             + a(i,3)*b(3,j)
     $             + a(i,4)*b(4,j)
     $             + a(i,5)*b(5,j)
     $             + a(i,6)*b(6,j)
     $             + a(i,7)*b(7,j)
     $             + a(i,8)*b(8,j)
     $             + a(i,9)*b(9,j)
     $             + a(i,10)*b(10,j)
     $             + a(i,11)*b(11,j)
     $             + a(i,12)*b(12,j)
     $             + a(i,13)*b(13,j)
     $             + a(i,14)*b(14,j)
     $             + a(i,15)*b(15,j)
     $             + a(i,16)*b(16,j)
     $             + a(i,17)*b(17,j)
     $             + a(i,18)*b(18,j)
     $             + a(i,19)*b(19,j)
     $             + a(i,20)*b(20,j)
     $             + a(i,21)*b(21,j)
     $             + a(i,22)*b(22,j)
     $             + a(i,23)*b(23,j)
     $             + a(i,24)*b(24,j)
         enddo
      enddo
      return
      end
c-----------------------------------------------------------------------
      subroutine mxm44_0(a, m, b, k, c, n)
c
c matrix multiply with a 4x4 pencil 
c
      real a(m,k), b(k,n), c(m,n)
      real s11, s12, s13, s14, s21, s22, s23, s24
      real s31, s32, s33, s34, s41, s42, s43, s44

      mresid = iand(m,3) 
      nresid = iand(n,3) 
      m1 = m - mresid + 1
      n1 = n - nresid + 1

      do i=1,m-mresid,4
        do j=1,n-nresid,4
          s11 = 0.0d0
          s21 = 0.0d0
          s31 = 0.0d0
          s41 = 0.0d0
          s12 = 0.0d0
          s22 = 0.0d0
          s32 = 0.0d0
          s42 = 0.0d0
          s13 = 0.0d0
          s23 = 0.0d0
          s33 = 0.0d0
          s43 = 0.0d0
          s14 = 0.0d0
          s24 = 0.0d0
          s34 = 0.0d0
          s44 = 0.0d0
          do l=1,k
            s11 = s11 + a(i,l)*b(l,j)
            s12 = s12 + a(i,l)*b(l,j+1)
            s13 = s13 + a(i,l)*b(l,j+2)
            s14 = s14 + a(i,l)*b(l,j+3)

            s21 = s21 + a(i+1,l)*b(l,j)
            s22 = s22 + a(i+1,l)*b(l,j+1)
            s23 = s23 + a(i+1,l)*b(l,j+2)
            s24 = s24 + a(i+1,l)*b(l,j+3)

            s31 = s31 + a(i+2,l)*b(l,j)
            s32 = s32 + a(i+2,l)*b(l,j+1)
            s33 = s33 + a(i+2,l)*b(l,j+2)
            s34 = s34 + a(i+2,l)*b(l,j+3)

            s41 = s41 + a(i+3,l)*b(l,j)
            s42 = s42 + a(i+3,l)*b(l,j+1)
            s43 = s43 + a(i+3,l)*b(l,j+2)
            s44 = s44 + a(i+3,l)*b(l,j+3)
          enddo
          c(i,j)     = s11 
          c(i,j+1)   = s12 
          c(i,j+2)   = s13
          c(i,j+3)   = s14

          c(i+1,j)   = s21 
          c(i+2,j)   = s31 
          c(i+3,j)   = s41 

          c(i+1,j+1) = s22
          c(i+2,j+1) = s32
          c(i+3,j+1) = s42

          c(i+1,j+2) = s23
          c(i+2,j+2) = s33
          c(i+3,j+2) = s43

          c(i+1,j+3) = s24
          c(i+2,j+3) = s34
          c(i+3,j+3) = s44
        enddo
* Residual when n is not multiple of 4
        if (nresid .ne. 0) then
          if (nresid .eq. 1) then
            s11 = 0.0d0
            s21 = 0.0d0
            s31 = 0.0d0
            s41 = 0.0d0
            do l=1,k
              s11 = s11 + a(i,l)*b(l,n)
              s21 = s21 + a(i+1,l)*b(l,n)
              s31 = s31 + a(i+2,l)*b(l,n)
              s41 = s41 + a(i+3,l)*b(l,n)
            enddo
            c(i,n)     = s11 
            c(i+1,n)   = s21 
            c(i+2,n)   = s31 
            c(i+3,n)   = s41 
          elseif (nresid .eq. 2) then
            s11 = 0.0d0
            s21 = 0.0d0
            s31 = 0.0d0
            s41 = 0.0d0
            s12 = 0.0d0
            s22 = 0.0d0
            s32 = 0.0d0
            s42 = 0.0d0
            do l=1,k
              s11 = s11 + a(i,l)*b(l,j)
              s12 = s12 + a(i,l)*b(l,j+1)

              s21 = s21 + a(i+1,l)*b(l,j)
              s22 = s22 + a(i+1,l)*b(l,j+1)

              s31 = s31 + a(i+2,l)*b(l,j)
              s32 = s32 + a(i+2,l)*b(l,j+1)

              s41 = s41 + a(i+3,l)*b(l,j)
              s42 = s42 + a(i+3,l)*b(l,j+1)
            enddo
            c(i,j)     = s11 
            c(i,j+1)   = s12

            c(i+1,j)   = s21 
            c(i+2,j)   = s31 
            c(i+3,j)   = s41 

            c(i+1,j+1) = s22
            c(i+2,j+1) = s32
            c(i+3,j+1) = s42
          else
            s11 = 0.0d0
            s21 = 0.0d0
            s31 = 0.0d0
            s41 = 0.0d0
            s12 = 0.0d0
            s22 = 0.0d0
            s32 = 0.0d0
            s42 = 0.0d0
            s13 = 0.0d0
            s23 = 0.0d0
            s33 = 0.0d0
            s43 = 0.0d0
            do l=1,k
              s11 = s11 + a(i,l)*b(l,j)
              s12 = s12 + a(i,l)*b(l,j+1)
              s13 = s13 + a(i,l)*b(l,j+2)

              s21 = s21 + a(i+1,l)*b(l,j)
              s22 = s22 + a(i+1,l)*b(l,j+1)
              s23 = s23 + a(i+1,l)*b(l,j+2)

              s31 = s31 + a(i+2,l)*b(l,j)
              s32 = s32 + a(i+2,l)*b(l,j+1)
              s33 = s33 + a(i+2,l)*b(l,j+2)

              s41 = s41 + a(i+3,l)*b(l,j)
              s42 = s42 + a(i+3,l)*b(l,j+1)
              s43 = s43 + a(i+3,l)*b(l,j+2)
            enddo
            c(i,j)     = s11 
            c(i+1,j)   = s21 
            c(i+2,j)   = s31 
            c(i+3,j)   = s41 
            c(i,j+1)   = s12 
            c(i+1,j+1) = s22
            c(i+2,j+1) = s32
            c(i+3,j+1) = s42
            c(i,j+2)   = s13
            c(i+1,j+2) = s23
            c(i+2,j+2) = s33
            c(i+3,j+2) = s43
          endif
        endif
      enddo

* Residual when m is not multiple of 4
      if (mresid .eq. 0) then
        return
      elseif (mresid .eq. 1) then
        do j=1,n-nresid,4
          s11 = 0.0d0
          s12 = 0.0d0
          s13 = 0.0d0
          s14 = 0.0d0
          do l=1,k
            s11 = s11 + a(m,l)*b(l,j)
            s12 = s12 + a(m,l)*b(l,j+1)
            s13 = s13 + a(m,l)*b(l,j+2)
            s14 = s14 + a(m,l)*b(l,j+3)
          enddo
          c(m,j)     = s11 
          c(m,j+1)   = s12 
          c(m,j+2)   = s13
          c(m,j+3)   = s14
        enddo
* mresid is 1, check nresid
        if (nresid .eq. 0) then
          return
        elseif (nresid .eq. 1) then
          s11 = 0.0d0
          do l=1,k
            s11 = s11 + a(m,l)*b(l,n)
          enddo
          c(m,n) = s11
          return
        elseif (nresid .eq. 2) then
          s11 = 0.0d0
          s12 = 0.0d0
          do l=1,k
            s11 = s11 + a(m,l)*b(l,n-1)
            s12 = s12 + a(m,l)*b(l,n)
          enddo
          c(m,n-1) = s11
          c(m,n) = s12
          return
        else
          s11 = 0.0d0
          s12 = 0.0d0
          s13 = 0.0d0
          do l=1,k
            s11 = s11 + a(m,l)*b(l,n-2)
            s12 = s12 + a(m,l)*b(l,n-1)
            s13 = s13 + a(m,l)*b(l,n)
          enddo
          c(m,n-2) = s11
          c(m,n-1) = s12
          c(m,n) = s13
          return
        endif          
      elseif (mresid .eq. 2) then
        do j=1,n-nresid,4
          s11 = 0.0d0
          s12 = 0.0d0
          s13 = 0.0d0
          s14 = 0.0d0
          s21 = 0.0d0
          s22 = 0.0d0
          s23 = 0.0d0
          s24 = 0.0d0
          do l=1,k
            s11 = s11 + a(m-1,l)*b(l,j)
            s12 = s12 + a(m-1,l)*b(l,j+1)
            s13 = s13 + a(m-1,l)*b(l,j+2)
            s14 = s14 + a(m-1,l)*b(l,j+3)

            s21 = s21 + a(m,l)*b(l,j)
            s22 = s22 + a(m,l)*b(l,j+1)
            s23 = s23 + a(m,l)*b(l,j+2)
            s24 = s24 + a(m,l)*b(l,j+3)
          enddo
          c(m-1,j)   = s11 
          c(m-1,j+1) = s12 
          c(m-1,j+2) = s13
          c(m-1,j+3) = s14
          c(m,j)     = s21
          c(m,j+1)   = s22 
          c(m,j+2)   = s23
          c(m,j+3)   = s24
        enddo
* mresid is 2, check nresid
        if (nresid .eq. 0) then
          return
        elseif (nresid .eq. 1) then
          s11 = 0.0d0
          s21 = 0.0d0
          do l=1,k
            s11 = s11 + a(m-1,l)*b(l,n)
            s21 = s21 + a(m,l)*b(l,n)
          enddo
          c(m-1,n) = s11
          c(m,n) = s21
          return
        elseif (nresid .eq. 2) then
          s11 = 0.0d0
          s21 = 0.0d0
          s12 = 0.0d0
          s22 = 0.0d0
          do l=1,k
            s11 = s11 + a(m-1,l)*b(l,n-1)
            s12 = s12 + a(m-1,l)*b(l,n)
            s21 = s21 + a(m,l)*b(l,n-1)
            s22 = s22 + a(m,l)*b(l,n)
          enddo
          c(m-1,n-1) = s11
          c(m-1,n)   = s12
          c(m,n-1)   = s21
          c(m,n)     = s22
          return
        else
          s11 = 0.0d0
          s21 = 0.0d0
          s12 = 0.0d0
          s22 = 0.0d0
          s13 = 0.0d0
          s23 = 0.0d0
          do l=1,k
            s11 = s11 + a(m-1,l)*b(l,n-2)
            s12 = s12 + a(m-1,l)*b(l,n-1)
            s13 = s13 + a(m-1,l)*b(l,n)
            s21 = s21 + a(m,l)*b(l,n-2)
            s22 = s22 + a(m,l)*b(l,n-1)
            s23 = s23 + a(m,l)*b(l,n)
          enddo
          c(m-1,n-2) = s11
          c(m-1,n-1) = s12
          c(m-1,n)   = s13
          c(m,n-2)   = s21
          c(m,n-1)   = s22
          c(m,n)     = s23
          return
        endif
      else
* mresid is 3
        do j=1,n-nresid,4
          s11 = 0.0d0
          s21 = 0.0d0
          s31 = 0.0d0

          s12 = 0.0d0
          s22 = 0.0d0
          s32 = 0.0d0

          s13 = 0.0d0
          s23 = 0.0d0
          s33 = 0.0d0

          s14 = 0.0d0
          s24 = 0.0d0
          s34 = 0.0d0

          do l=1,k
            s11 = s11 + a(m-2,l)*b(l,j)
            s12 = s12 + a(m-2,l)*b(l,j+1)
            s13 = s13 + a(m-2,l)*b(l,j+2)
            s14 = s14 + a(m-2,l)*b(l,j+3)

            s21 = s21 + a(m-1,l)*b(l,j)
            s22 = s22 + a(m-1,l)*b(l,j+1)
            s23 = s23 + a(m-1,l)*b(l,j+2)
            s24 = s24 + a(m-1,l)*b(l,j+3)

            s31 = s31 + a(m,l)*b(l,j)
            s32 = s32 + a(m,l)*b(l,j+1)
            s33 = s33 + a(m,l)*b(l,j+2)
            s34 = s34 + a(m,l)*b(l,j+3)
          enddo
          c(m-2,j)   = s11 
          c(m-2,j+1) = s12 
          c(m-2,j+2) = s13
          c(m-2,j+3) = s14

          c(m-1,j)   = s21 
          c(m-1,j+1) = s22
          c(m-1,j+2) = s23
          c(m-1,j+3) = s24

          c(m,j)     = s31 
          c(m,j+1)   = s32
          c(m,j+2)   = s33
          c(m,j+3)   = s34
        enddo
* mresid is 3, check nresid
        if (nresid .eq. 0) then
          return
        elseif (nresid .eq. 1) then
          s11 = 0.0d0
          s21 = 0.0d0
          s31 = 0.0d0
          do l=1,k
            s11 = s11 + a(m-2,l)*b(l,n)
            s21 = s21 + a(m-1,l)*b(l,n)
            s31 = s31 + a(m,l)*b(l,n)
          enddo
          c(m-2,n) = s11
          c(m-1,n) = s21
          c(m,n) = s31
          return
        elseif (nresid .eq. 2) then
          s11 = 0.0d0
          s21 = 0.0d0
          s31 = 0.0d0
          s12 = 0.0d0
          s22 = 0.0d0
          s32 = 0.0d0
          do l=1,k
            s11 = s11 + a(m-2,l)*b(l,n-1)
            s12 = s12 + a(m-2,l)*b(l,n)
            s21 = s21 + a(m-1,l)*b(l,n-1)
            s22 = s22 + a(m-1,l)*b(l,n)
            s31 = s31 + a(m,l)*b(l,n-1)
            s32 = s32 + a(m,l)*b(l,n)
          enddo
          c(m-2,n-1) = s11
          c(m-2,n)   = s12
          c(m-1,n-1) = s21
          c(m-1,n)   = s22
          c(m,n-1)   = s31
          c(m,n)     = s32
          return
        else
          s11 = 0.0d0
          s21 = 0.0d0
          s31 = 0.0d0
          s12 = 0.0d0
          s22 = 0.0d0
          s32 = 0.0d0
          s13 = 0.0d0
          s23 = 0.0d0
          s33 = 0.0d0
          do l=1,k
            s11 = s11 + a(m-2,l)*b(l,n-2)
            s12 = s12 + a(m-2,l)*b(l,n-1)
            s13 = s13 + a(m-2,l)*b(l,n)
            s21 = s21 + a(m-1,l)*b(l,n-2)
            s22 = s22 + a(m-1,l)*b(l,n-1)
            s23 = s23 + a(m-1,l)*b(l,n)
            s31 = s31 + a(m,l)*b(l,n-2)
            s32 = s32 + a(m,l)*b(l,n-1)
            s33 = s33 + a(m,l)*b(l,n)
          enddo
          c(m-2,n-2) = s11
          c(m-2,n-1) = s12
          c(m-2,n)   = s13
          c(m-1,n-2) = s21
          c(m-1,n-1) = s22
          c(m-1,n)   = s23
          c(m,n-2)   = s31
          c(m,n-1)   = s32
          c(m,n)     = s33
          return
        endif
      endif

      return
      end
c-----------------------------------------------------------------------
