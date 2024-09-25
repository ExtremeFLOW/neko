! Copyright (c) 2021-2024, The Neko Authors
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
!> conv1 kernels
submodule (opr_cpu) cpu_conv1
  implicit none

contains

  module subroutine opr_cpu_conv1(du, u, vx, vy, vz, Xh, coef, e_start, e_end)
    type(space_t), intent(in) :: Xh
    type(coef_t), intent(in) :: coef
    integer, intent(in) :: e_start, e_end
    real(kind=rp), intent(inout) ::  du(Xh%lxyz, e_end-e_start+1)
    real(kind=rp), intent(inout) ::  u(Xh%lx, Xh%ly, Xh%lz, e_end-e_start+1)
    real(kind=rp), intent(inout) ::  vx(Xh%lx, Xh%ly, Xh%lz, e_end-e_start+1)
    real(kind=rp), intent(inout) ::  vy(Xh%lx, Xh%ly, Xh%lz, e_end-e_start+1)
    real(kind=rp), intent(inout) ::  vz(Xh%lx, Xh%ly, Xh%lz, e_end-e_start+1)
    integer :: e_len

    e_len = e_end-e_start+1

    if (e_len .eq. 1) then
       call opr_cpu_conv1_single(du, u, vx, vy, vz, Xh, coef, e_start)
    else
       call opr_cpu_conv1_many(du, u, vx, vy, vz, Xh, coef, e_start, e_len)
    end if
  end subroutine opr_cpu_conv1

  subroutine opr_cpu_conv1_many(du, u, vx, vy, vz, Xh, coef, e_start, e_len)
    type(space_t), intent(in) :: Xh
    type(coef_t), intent(in) :: coef
    integer, intent(in) :: e_start, e_len
    real(kind=rp), intent(inout) ::  du(Xh%lxyz, e_len)
    real(kind=rp), intent(inout) ::  u(Xh%lx, Xh%ly, Xh%lz, e_len)
    real(kind=rp), intent(inout) ::  vx(Xh%lx, Xh%ly, Xh%lz, e_len)
    real(kind=rp), intent(inout) ::  vy(Xh%lx, Xh%ly, Xh%lz, e_len)
    real(kind=rp), intent(inout) ::  vz(Xh%lx, Xh%ly, Xh%lz, e_len)
    
    associate(drdx => coef%drdx, drdy => coef%drdy, drdz => coef%drdz, &
         dsdx => coef%dsdx, dsdy => coef%dsdy, dsdz => coef%dsdz, &
         dtdx => coef%dtdx, dtdy => coef%dtdy, dtdz => coef%dtdz, &
         jacinv => coef%jacinv)
      select case (Xh%lx)
      case (14)
         call cpu_conv1_lx14(du, u, vx, vy, vz, Xh%dx, Xh%dy, Xh%dz, &
              drdx(1,1,1, e_start), dsdx(1,1,1, e_start), dtdx(1,1,1, e_start),&
              drdy(1,1,1, e_start), dsdy(1,1,1, e_start), dtdy(1,1,1, e_start),&
              drdz(1,1,1, e_start), dsdz(1,1,1, e_start), dtdz(1,1,1, e_start),&
              jacinv(1,1,1, e_start), e_len)
      case (13)
         call cpu_conv1_lx13(du, u, vx, vy, vz, Xh%dx, Xh%dy, Xh%dz, &
              drdx(1,1,1, e_start), dsdx(1,1,1, e_start), dtdx(1,1,1, e_start),&
              drdy(1,1,1, e_start), dsdy(1,1,1, e_start), dtdy(1,1,1, e_start),&
              drdz(1,1,1, e_start), dsdz(1,1,1, e_start), dtdz(1,1,1, e_start),&
              jacinv(1,1,1, e_start), e_len)
      case (12)
         call cpu_conv1_lx12(du, u, vx, vy, vz, Xh%dx, Xh%dy, Xh%dz, &
              drdx(1,1,1, e_start), dsdx(1,1,1, e_start), dtdx(1,1,1, e_start),&
              drdy(1,1,1, e_start), dsdy(1,1,1, e_start), dtdy(1,1,1, e_start),&
              drdz(1,1,1, e_start), dsdz(1,1,1, e_start), dtdz(1,1,1, e_start),&
              jacinv(1,1,1, e_start), e_len)
      case (11)
         call cpu_conv1_lx11(du, u, vx, vy, vz, Xh%dx, Xh%dy, Xh%dz, &
              drdx(1,1,1, e_start), dsdx(1,1,1, e_start), dtdx(1,1,1, e_start),&
              drdy(1,1,1, e_start), dsdy(1,1,1, e_start), dtdy(1,1,1, e_start),&
              drdz(1,1,1, e_start), dsdz(1,1,1, e_start), dtdz(1,1,1, e_start),&
              jacinv(1,1,1, e_start), e_len)
      case (10)
         call cpu_conv1_lx10(du, u, vx, vy, vz, Xh%dx, Xh%dy, Xh%dz, &
              drdx(1,1,1, e_start), dsdx(1,1,1, e_start), dtdx(1,1,1, e_start),&
              drdy(1,1,1, e_start), dsdy(1,1,1, e_start), dtdy(1,1,1, e_start),&
              drdz(1,1,1, e_start), dsdz(1,1,1, e_start), dtdz(1,1,1, e_start),&
              jacinv(1,1,1, e_start), e_len)
      case (9)
         call cpu_conv1_lx9(du, u, vx, vy, vz, Xh%dx, Xh%dy, Xh%dz, &
              drdx(1,1,1, e_start), dsdx(1,1,1, e_start), dtdx(1,1,1, e_start),&
              drdy(1,1,1, e_start), dsdy(1,1,1, e_start), dtdy(1,1,1, e_start),&
              drdz(1,1,1, e_start), dsdz(1,1,1, e_start), dtdz(1,1,1, e_start),&
              jacinv(1,1,1, e_start), e_len)
      case (8)
         call cpu_conv1_lx8(du, u, vx, vy, vz, Xh%dx, Xh%dy, Xh%dz, &
              drdx(1,1,1, e_start), dsdx(1,1,1, e_start), dtdx(1,1,1, e_start),&
              drdy(1,1,1, e_start), dsdy(1,1,1, e_start), dtdy(1,1,1, e_start),&
              drdz(1,1,1, e_start), dsdz(1,1,1, e_start), dtdz(1,1,1, e_start),&
              jacinv(1,1,1, e_start), e_len)
      case (7)
         call cpu_conv1_lx7(du, u, vx, vy, vz, Xh%dx, Xh%dy, Xh%dz, &
              drdx(1,1,1, e_start), dsdx(1,1,1, e_start), dtdx(1,1,1, e_start),&
              drdy(1,1,1, e_start), dsdy(1,1,1, e_start), dtdy(1,1,1, e_start),&
              drdz(1,1,1, e_start), dsdz(1,1,1, e_start), dtdz(1,1,1, e_start),&
              jacinv(1,1,1, e_start), e_len)
      case (6)
         call cpu_conv1_lx6(du, u, vx, vy, vz, Xh%dx, Xh%dy, Xh%dz, &
              drdx(1,1,1, e_start), dsdx(1,1,1, e_start), dtdx(1,1,1, e_start),&
              drdy(1,1,1, e_start), dsdy(1,1,1, e_start), dtdy(1,1,1, e_start),&
              drdz(1,1,1, e_start), dsdz(1,1,1, e_start), dtdz(1,1,1, e_start),&
              jacinv(1,1,1, e_start), e_len)
      case (5)
         call cpu_conv1_lx5(du, u, vx, vy, vz, Xh%dx, Xh%dy, Xh%dz, &
              drdx(1,1,1, e_start), dsdx(1,1,1, e_start), dtdx(1,1,1, e_start),&
              drdy(1,1,1, e_start), dsdy(1,1,1, e_start), dtdy(1,1,1, e_start),&
              drdz(1,1,1, e_start), dsdz(1,1,1, e_start), dtdz(1,1,1, e_start),&
              jacinv(1,1,1, e_start), e_len)
      case (4)
         call cpu_conv1_lx4(du, u, vx, vy, vz, Xh%dx, Xh%dy, Xh%dz, &
              drdx(1,1,1, e_start), dsdx(1,1,1, e_start), dtdx(1,1,1, e_start),&
              drdy(1,1,1, e_start), dsdy(1,1,1, e_start), dtdy(1,1,1, e_start),&
              drdz(1,1,1, e_start), dsdz(1,1,1, e_start), dtdz(1,1,1, e_start),&
              jacinv(1,1,1, e_start), e_len)
      case (3)
         call cpu_conv1_lx3(du, u, vx, vy, vz, Xh%dx, Xh%dy, Xh%dz, &
              drdx(1,1,1, e_start), dsdx(1,1,1, e_start), dtdx(1,1,1, e_start),&
              drdy(1,1,1, e_start), dsdy(1,1,1, e_start), dtdy(1,1,1, e_start),&
              drdz(1,1,1, e_start), dsdz(1,1,1, e_start), dtdz(1,1,1, e_start),&
              jacinv(1,1,1, e_start), e_len)
      case (2)
         call cpu_conv1_lx2(du, u, vx, vy, vz, Xh%dx, Xh%dy, Xh%dz, &
              drdx(1,1,1, e_start), dsdx(1,1,1, e_start), dtdx(1,1,1, e_start),&
              drdy(1,1,1, e_start), dsdy(1,1,1, e_start), dtdy(1,1,1, e_start),&
              drdz(1,1,1, e_start), dsdz(1,1,1, e_start), dtdz(1,1,1, e_start),&
              jacinv(1,1,1, e_start), e_len)
      case default
         call cpu_conv1_lx(du, u, vx, vy, vz, Xh%dx, Xh%dy, Xh%dz, &
              drdx(1,1,1, e_start), dsdx(1,1,1, e_start), dtdx(1,1,1, e_start),&
              drdy(1,1,1, e_start), dsdy(1,1,1, e_start), dtdy(1,1,1, e_start),&
              drdz(1,1,1, e_start), dsdz(1,1,1, e_start), dtdz(1,1,1, e_start),&
              jacinv(1,1,1, e_start), e_len, Xh%lx)
      end select
    end associate

  end subroutine opr_cpu_conv1_many

  subroutine cpu_conv1_lx(du, u, vx, vy, vz, dx, dy, dz, &
       drdx, dsdx, dtdx, drdy, dsdy, dtdy, drdz, dsdz, dtdz, &
       jacinv, nelv, lx)
    integer, intent(in) :: nelv, lx
    real(kind=rp), dimension(lx, lx, lx, nelv), intent(inout) ::  du
    real(kind=rp), dimension(lx, lx, lx, nelv), intent(in) ::  u, vx, vy, vz
    real(kind=rp), dimension(lx, lx, lx, nelv), intent(in) :: drdx, dsdx, dtdx
    real(kind=rp), dimension(lx, lx, lx, nelv), intent(in) :: drdy, dsdy, dtdy
    real(kind=rp), dimension(lx, lx, lx, nelv), intent(in) :: drdz, dsdz, dtdz
    real(kind=rp), dimension(lx, lx, lx, nelv), intent(in) :: jacinv
    real(kind=rp), dimension(lx, lx), intent(in) :: dx, dy, dz
    real(kind=rp), dimension(lx, lx, lx) ::  dudr
    real(kind=rp), dimension(lx, lx, lx) ::  duds
    real(kind=rp), dimension(lx, lx, lx) ::  dudt
    real(kind=rp) :: tmp
    integer :: e, i, j, k, l

    do e = 1, nelv
       do j = 1, lx * lx
          do i = 1, lx
             tmp = 0.0_rp
             do k = 1, lx
                tmp = tmp + dx(i,k) * u(k,j,1,e)
             end do
             dudr(i,j,1) = tmp
          end do
       end do

       do k = 1, lx
          do j = 1, lx
             do i = 1, lx
                tmp = 0.0_rp
                do l = 1, lx
                   tmp = tmp + dy(j,l) * u(i,l,k,e)
                end do
                duds(i,j,k) = tmp
             end do
          end do
       end do

       do k = 1, lx
          do i = 1, lx*lx
             tmp = 0.0_rp
             do l = 1, lx
                tmp = tmp + dz(k,l) * u(i,1,l,e)
             end do
             dudt(i,1,k) = tmp
          end do
       end do

       do i = 1, lx * lx * lx
          du(i,1,1,e) = jacinv(i,1,1,e) &
                      * ( vx(i,1,1,e) &
                        * ( drdx(i,1,1,e) * dudr(i,1,1) &
                          + dsdx(i,1,1,e) * duds(i,1,1) &
                          + dtdx(i,1,1,e) * dudt(i,1,1) ) &
                        + vy(i,1,1,e) &
                        * ( drdy(i,1,1,e) * dudr(i,1,1) &
                          + dsdy(i,1,1,e) * duds(i,1,1) &
                          + dtdy(i,1,1,e) * dudt(i,1,1) ) &
                        + vz(i,1,1,e) &
                        * ( drdz(i,1,1,e) * dudr(i,1,1) &
                          + dsdz(i,1,1,e) * duds(i,1,1) &
                          + dtdz(i,1,1,e) * dudt(i,1,1) ) )
       end do
    end do

  end subroutine cpu_conv1_lx

  subroutine cpu_conv1_lx14(du, u, vx, vy, vz, dx, dy, dz, &
       drdx, dsdx, dtdx, drdy, dsdy, dtdy, drdz, dsdz, dtdz, &
       jacinv, nelv)
    integer, parameter :: lx = 14
    integer, intent(in) :: nelv
    real(kind=rp), dimension(lx, lx, lx, nelv), intent(inout) ::  du
    real(kind=rp), dimension(lx, lx, lx, nelv), intent(in) ::  u, vx, vy, vz
    real(kind=rp), dimension(lx, lx, lx, nelv), intent(in) :: drdx, dsdx, dtdx
    real(kind=rp), dimension(lx, lx, lx, nelv), intent(in) :: drdy, dsdy, dtdy
    real(kind=rp), dimension(lx, lx, lx, nelv), intent(in) :: drdz, dsdz, dtdz
    real(kind=rp), dimension(lx, lx, lx, nelv), intent(in) :: jacinv
    real(kind=rp), dimension(lx, lx), intent(in) :: dx, dy, dz
    real(kind=rp), dimension(lx, lx, lx) ::  dudr
    real(kind=rp), dimension(lx, lx, lx) ::  duds
    real(kind=rp), dimension(lx, lx, lx) ::  dudt
    integer :: e, i, j, k

    do e = 1, nelv
       do j = 1, lx * lx
          do i = 1, lx
             dudr(i,j,1) = dx(i,1) * u(1,j,1,e) &
                         + dx(i,2) * u(2,j,1,e) &
                         + dx(i,3) * u(3,j,1,e) &
                         + dx(i,4) * u(4,j,1,e) &
                         + dx(i,5) * u(5,j,1,e) &
                         + dx(i,6) * u(6,j,1,e) &
                         + dx(i,7) * u(7,j,1,e) &
                         + dx(i,8) * u(8,j,1,e) &
                         + dx(i,9) * u(9,j,1,e) &
                         + dx(i,10) * u(10,j,1,e) &
                         + dx(i,11) * u(11,j,1,e) &
                         + dx(i,12) * u(12,j,1,e) &
                         + dx(i,13) * u(13,j,1,e) &
                         + dx(i,14) * u(14,j,1,e)
          end do
       end do

       do k = 1, lx
          do j = 1, lx
             do i = 1, lx
                duds(i,j,k) = dy(j,1) * u(i,1,k,e) &
                            + dy(j,2) * u(i,2,k,e) &
                            + dy(j,3) * u(i,3,k,e) &
                            + dy(j,4) * u(i,4,k,e) &
                            + dy(j,5) * u(i,5,k,e) &
                            + dy(j,6) * u(i,6,k,e) &
                            + dy(j,7) * u(i,7,k,e) &
                            + dy(j,8) * u(i,8,k,e) &
                            + dy(j,9) * u(i,9,k,e) &
                            + dy(j,10) * u(i,10,k,e) &
                            + dy(j,11) * u(i,11,k,e) &
                            + dy(j,12) * u(i,12,k,e) &
                            + dy(j,13) * u(i,13,k,e) &
                            + dy(j,14) * u(i,14,k,e)
             end do
          end do
       end do

       do k = 1, lx
          do i = 1, lx*lx
             dudt(i,1,k) = dz(k,1) * u(i,1,1,e) &
                         + dz(k,2) * u(i,1,2,e) &
                         + dz(k,3) * u(i,1,3,e) &
                         + dz(k,4) * u(i,1,4,e) &
                         + dz(k,5) * u(i,1,5,e) &
                         + dz(k,6) * u(i,1,6,e) &
                         + dz(k,7) * u(i,1,7,e) &
                         + dz(k,8) * u(i,1,8,e) &
                         + dz(k,9) * u(i,1,9,e) &
                         + dz(k,10) * u(i,1,10,e) &
                         + dz(k,11) * u(i,1,11,e) &
                         + dz(k,12) * u(i,1,12,e) &
                         + dz(k,13) * u(i,1,13,e) &
                         + dz(k,14) * u(i,1,14,e)
          end do
       end do

       do i = 1, lx * lx * lx
          du(i,1,1,e) = jacinv(i,1,1,e) &
                      * ( vx(i,1,1,e) &
                        * ( drdx(i,1,1,e) * dudr(i,1,1) &
                          + dsdx(i,1,1,e) * duds(i,1,1) &
                          + dtdx(i,1,1,e) * dudt(i,1,1) ) &
                        + vy(i,1,1,e) &
                        * ( drdy(i,1,1,e) * dudr(i,1,1) &
                          + dsdy(i,1,1,e) * duds(i,1,1) &
                          + dtdy(i,1,1,e) * dudt(i,1,1) ) &
                        + vz(i,1,1,e) &
                        * ( drdz(i,1,1,e) * dudr(i,1,1) &
                          + dsdz(i,1,1,e) * duds(i,1,1) &
                          + dtdz(i,1,1,e) * dudt(i,1,1) ) )
       end do
    end do

  end subroutine cpu_conv1_lx14

  subroutine cpu_conv1_lx13(du, u, vx, vy, vz, dx, dy, dz, &
       drdx, dsdx, dtdx, drdy, dsdy, dtdy, drdz, dsdz, dtdz, &
       jacinv, nelv)
    integer, parameter :: lx = 13
    integer, intent(in) :: nelv
    real(kind=rp), dimension(lx, lx, lx, nelv), intent(inout) ::  du
    real(kind=rp), dimension(lx, lx, lx, nelv), intent(in) ::  u, vx, vy, vz
    real(kind=rp), dimension(lx, lx, lx, nelv), intent(in) :: drdx, dsdx, dtdx
    real(kind=rp), dimension(lx, lx, lx, nelv), intent(in) :: drdy, dsdy, dtdy
    real(kind=rp), dimension(lx, lx, lx, nelv), intent(in) :: drdz, dsdz, dtdz
    real(kind=rp), dimension(lx, lx, lx, nelv), intent(in) :: jacinv
    real(kind=rp), dimension(lx, lx), intent(in) :: dx, dy, dz
    real(kind=rp), dimension(lx, lx, lx) ::  dudr
    real(kind=rp), dimension(lx, lx, lx) ::  duds
    real(kind=rp), dimension(lx, lx, lx) ::  dudt
    integer :: e, i, j, k

    do e = 1, nelv
       do j = 1, lx * lx
          do i = 1, lx
             dudr(i,j,1) = dx(i,1) * u(1,j,1,e) &
                         + dx(i,2) * u(2,j,1,e) &
                         + dx(i,3) * u(3,j,1,e) &
                         + dx(i,4) * u(4,j,1,e) &
                         + dx(i,5) * u(5,j,1,e) &
                         + dx(i,6) * u(6,j,1,e) &
                         + dx(i,7) * u(7,j,1,e) &
                         + dx(i,8) * u(8,j,1,e) &
                         + dx(i,9) * u(9,j,1,e) &
                         + dx(i,10) * u(10,j,1,e) &
                         + dx(i,11) * u(11,j,1,e) &
                         + dx(i,12) * u(12,j,1,e) &
                         + dx(i,13) * u(13,j,1,e)
          end do
       end do

       do k = 1, lx
          do j = 1, lx
             do i = 1, lx
                duds(i,j,k) = dy(j,1) * u(i,1,k,e) &
                            + dy(j,2) * u(i,2,k,e) &
                            + dy(j,3) * u(i,3,k,e) &
                            + dy(j,4) * u(i,4,k,e) &
                            + dy(j,5) * u(i,5,k,e) &
                            + dy(j,6) * u(i,6,k,e) &
                            + dy(j,7) * u(i,7,k,e) &
                            + dy(j,8) * u(i,8,k,e) &
                            + dy(j,9) * u(i,9,k,e) &
                            + dy(j,10) * u(i,10,k,e) &
                            + dy(j,11) * u(i,11,k,e) &
                            + dy(j,12) * u(i,12,k,e) &
                            + dy(j,13) * u(i,13,k,e)
             end do
          end do
       end do

       do k = 1, lx
          do i = 1, lx*lx
             dudt(i,1,k) = dz(k,1) * u(i,1,1,e) &
                         + dz(k,2) * u(i,1,2,e) &
                         + dz(k,3) * u(i,1,3,e) &
                         + dz(k,4) * u(i,1,4,e) &
                         + dz(k,5) * u(i,1,5,e) &
                         + dz(k,6) * u(i,1,6,e) &
                         + dz(k,7) * u(i,1,7,e) &
                         + dz(k,8) * u(i,1,8,e) &
                         + dz(k,9) * u(i,1,9,e) &
                         + dz(k,10) * u(i,1,10,e) &
                         + dz(k,11) * u(i,1,11,e) &
                         + dz(k,12) * u(i,1,12,e) &
                         + dz(k,13) * u(i,1,13,e)
          end do
       end do

       do i = 1, lx * lx * lx
          du(i,1,1,e) = jacinv(i,1,1,e) &
                     * ( vx(i,1,1,e) &
                       * ( drdx(i,1,1,e) * dudr(i,1,1) &
                         + dsdx(i,1,1,e) * duds(i,1,1) &
                         + dtdx(i,1,1,e) * dudt(i,1,1) ) &
                       + vy(i,1,1,e) &
                       * ( drdy(i,1,1,e) * dudr(i,1,1) &
                         + dsdy(i,1,1,e) * duds(i,1,1) &
                         + dtdy(i,1,1,e) * dudt(i,1,1) ) &
                       + vz(i,1,1,e) &
                       * ( drdz(i,1,1,e) * dudr(i,1,1) &
                         + dsdz(i,1,1,e) * duds(i,1,1) &
                         + dtdz(i,1,1,e) * dudt(i,1,1) ) )
       end do
    end do

  end subroutine cpu_conv1_lx13

  subroutine cpu_conv1_lx12(du, u, vx, vy, vz, dx, dy, dz, &
       drdx, dsdx, dtdx, drdy, dsdy, dtdy, drdz, dsdz, dtdz, &
       jacinv, nelv)
    integer, parameter :: lx = 12
    integer, intent(in) :: nelv
    real(kind=rp), dimension(lx, lx, lx, nelv), intent(inout) ::  du
    real(kind=rp), dimension(lx, lx, lx, nelv), intent(in) ::  u, vx, vy, vz
    real(kind=rp), dimension(lx, lx, lx, nelv), intent(in) :: drdx, dsdx, dtdx
    real(kind=rp), dimension(lx, lx, lx, nelv), intent(in) :: drdy, dsdy, dtdy
    real(kind=rp), dimension(lx, lx, lx, nelv), intent(in) :: drdz, dsdz, dtdz
    real(kind=rp), dimension(lx, lx, lx, nelv), intent(in) :: jacinv
    real(kind=rp), dimension(lx, lx), intent(in) :: dx, dy, dz
    real(kind=rp), dimension(lx, lx, lx) ::  dudr
    real(kind=rp), dimension(lx, lx, lx) ::  duds
    real(kind=rp), dimension(lx, lx, lx) ::  dudt
    integer :: e, i, j, k

    do e = 1, nelv
       do j = 1, lx * lx
          do i = 1, lx
             dudr(i,j,1) = dx(i,1) * u(1,j,1,e) &
                         + dx(i,2) * u(2,j,1,e) &
                         + dx(i,3) * u(3,j,1,e) &
                         + dx(i,4) * u(4,j,1,e) &
                         + dx(i,5) * u(5,j,1,e) &
                         + dx(i,6) * u(6,j,1,e) &
                         + dx(i,7) * u(7,j,1,e) &
                         + dx(i,8) * u(8,j,1,e) &
                         + dx(i,9) * u(9,j,1,e) &
                         + dx(i,10) * u(10,j,1,e) &
                         + dx(i,11) * u(11,j,1,e) &
                         + dx(i,12) * u(12,j,1,e)
          end do
       end do

       do k = 1, lx
          do j = 1, lx
             do i = 1, lx
                duds(i,j,k) = dy(j,1) * u(i,1,k,e) &
                            + dy(j,2) * u(i,2,k,e) &
                            + dy(j,3) * u(i,3,k,e) &
                            + dy(j,4) * u(i,4,k,e) &
                            + dy(j,5) * u(i,5,k,e) &
                            + dy(j,6) * u(i,6,k,e) &
                            + dy(j,7) * u(i,7,k,e) &
                            + dy(j,8) * u(i,8,k,e) &
                            + dy(j,9) * u(i,9,k,e) &
                            + dy(j,10) * u(i,10,k,e) &
                            + dy(j,11) * u(i,11,k,e) &
                            + dy(j,12) * u(i,12,k,e)
             end do
          end do
       end do

       do k = 1, lx
          do i = 1, lx*lx
             dudt(i,1,k) = dz(k,1) * u(i,1,1,e) &
                         + dz(k,2) * u(i,1,2,e) &
                         + dz(k,3) * u(i,1,3,e) &
                         + dz(k,4) * u(i,1,4,e) &
                         + dz(k,5) * u(i,1,5,e) &
                         + dz(k,6) * u(i,1,6,e) &
                         + dz(k,7) * u(i,1,7,e) &
                         + dz(k,8) * u(i,1,8,e) &
                         + dz(k,9) * u(i,1,9,e) &
                         + dz(k,10) * u(i,1,10,e) &
                         + dz(k,11) * u(i,1,11,e) &
                         + dz(k,12) * u(i,1,12,e)
          end do
       end do

       do i = 1, lx * lx * lx
          du(i,1,1,e) = jacinv(i,1,1,e) &
                       * ( vx(i,1,1,e) &
                         * ( drdx(i,1,1,e) * dudr(i,1,1) &
                           + dsdx(i,1,1,e) * duds(i,1,1) &
                           + dtdx(i,1,1,e) * dudt(i,1,1) ) &
                         + vy(i,1,1,e) &
                         * ( drdy(i,1,1,e) * dudr(i,1,1) &
                           + dsdy(i,1,1,e) * duds(i,1,1) &
                           + dtdy(i,1,1,e) * dudt(i,1,1) ) &
                         + vz(i,1,1,e) &
                         * ( drdz(i,1,1,e) * dudr(i,1,1) &
                           + dsdz(i,1,1,e) * duds(i,1,1) &
                           + dtdz(i,1,1,e) * dudt(i,1,1) ) )
       end do
    end do

  end subroutine cpu_conv1_lx12

  subroutine cpu_conv1_lx11(du, u, vx, vy, vz, dx, dy, dz, &
       drdx, dsdx, dtdx, drdy, dsdy, dtdy, drdz, dsdz, dtdz, &
       jacinv, nelv)
    integer, parameter :: lx = 11
    integer, intent(in) :: nelv
    real(kind=rp), dimension(lx, lx, lx, nelv), intent(inout) ::  du
    real(kind=rp), dimension(lx, lx, lx, nelv), intent(in) ::  u, vx, vy, vz
    real(kind=rp), dimension(lx, lx, lx, nelv), intent(in) :: drdx, dsdx, dtdx
    real(kind=rp), dimension(lx, lx, lx, nelv), intent(in) :: drdy, dsdy, dtdy
    real(kind=rp), dimension(lx, lx, lx, nelv), intent(in) :: drdz, dsdz, dtdz
    real(kind=rp), dimension(lx, lx, lx, nelv), intent(in) :: jacinv
    real(kind=rp), dimension(lx, lx), intent(in) :: dx, dy, dz
    real(kind=rp), dimension(lx, lx, lx) ::  dudr
    real(kind=rp), dimension(lx, lx, lx) ::  duds
    real(kind=rp), dimension(lx, lx, lx) ::  dudt
    integer :: e, i, j, k

    do e = 1, nelv
       do j = 1, lx * lx
          do i = 1, lx
             dudr(i,j,1) = dx(i,1) * u(1,j,1,e) &
                         + dx(i,2) * u(2,j,1,e) &
                         + dx(i,3) * u(3,j,1,e) &
                         + dx(i,4) * u(4,j,1,e) &
                         + dx(i,5) * u(5,j,1,e) &
                         + dx(i,6) * u(6,j,1,e) &
                         + dx(i,7) * u(7,j,1,e) &
                         + dx(i,8) * u(8,j,1,e) &
                         + dx(i,9) * u(9,j,1,e) &
                         + dx(i,10) * u(10,j,1,e) &
                         + dx(i,11) * u(11,j,1,e)
          end do
       end do

       do k = 1, lx
          do j = 1, lx
             do i = 1, lx
                duds(i,j,k) = dy(j,1) * u(i,1,k,e) &
                            + dy(j,2) * u(i,2,k,e) &
                            + dy(j,3) * u(i,3,k,e) &
                            + dy(j,4) * u(i,4,k,e) &
                            + dy(j,5) * u(i,5,k,e) &
                            + dy(j,6) * u(i,6,k,e) &
                            + dy(j,7) * u(i,7,k,e) &
                            + dy(j,8) * u(i,8,k,e) &
                            + dy(j,9) * u(i,9,k,e) &
                            + dy(j,10) * u(i,10,k,e) &
                            + dy(j,11) * u(i,11,k,e)
             end do
          end do
       end do

       do k = 1, lx
          do i = 1, lx*lx
             dudt(i,1,k) = dz(k,1) * u(i,1,1,e) &
                         + dz(k,2) * u(i,1,2,e) &
                         + dz(k,3) * u(i,1,3,e) &
                         + dz(k,4) * u(i,1,4,e) &
                         + dz(k,5) * u(i,1,5,e) &
                         + dz(k,6) * u(i,1,6,e) &
                         + dz(k,7) * u(i,1,7,e) &
                         + dz(k,8) * u(i,1,8,e) &
                         + dz(k,9) * u(i,1,9,e) &
                         + dz(k,10) * u(i,1,10,e) &
                         + dz(k,11) * u(i,1,11,e)
          end do
       end do

       do i = 1, lx * lx * lx
          du(i,1,1,e) = jacinv(i,1,1,e) &
                       * ( vx(i,1,1,e) &
                         * ( drdx(i,1,1,e) * dudr(i,1,1) &
                           + dsdx(i,1,1,e) * duds(i,1,1) &
                           + dtdx(i,1,1,e) * dudt(i,1,1) ) &
                         + vy(i,1,1,e) &
                         * ( drdy(i,1,1,e) * dudr(i,1,1) &
                           + dsdy(i,1,1,e) * duds(i,1,1) &
                           + dtdy(i,1,1,e) * dudt(i,1,1) ) &
                         + vz(i,1,1,e) &
                         * ( drdz(i,1,1,e) * dudr(i,1,1) &
                           + dsdz(i,1,1,e) * duds(i,1,1) &
                           + dtdz(i,1,1,e) * dudt(i,1,1) ) )
       end do
    end do

  end subroutine cpu_conv1_lx11

  subroutine cpu_conv1_lx10(du, u, vx, vy, vz, dx, dy, dz, &
       drdx, dsdx, dtdx, drdy, dsdy, dtdy, drdz, dsdz, dtdz, &
       jacinv, nelv)
    integer, parameter :: lx = 10
    integer, intent(in) :: nelv
    real(kind=rp), dimension(lx, lx, lx, nelv), intent(inout) ::  du
    real(kind=rp), dimension(lx, lx, lx, nelv), intent(in) ::  u, vx, vy, vz
    real(kind=rp), dimension(lx, lx, lx, nelv), intent(in) :: drdx, dsdx, dtdx
    real(kind=rp), dimension(lx, lx, lx, nelv), intent(in) :: drdy, dsdy, dtdy
    real(kind=rp), dimension(lx, lx, lx, nelv), intent(in) :: drdz, dsdz, dtdz
    real(kind=rp), dimension(lx, lx, lx, nelv), intent(in) :: jacinv
    real(kind=rp), dimension(lx, lx), intent(in) :: dx, dy, dz
    real(kind=rp), dimension(lx, lx, lx) ::  dudr
    real(kind=rp), dimension(lx, lx, lx) ::  duds
    real(kind=rp), dimension(lx, lx, lx) ::  dudt
    integer :: e, i, j, k

    do e = 1, nelv
       do j = 1, lx * lx
          do i = 1, lx
             dudr(i,j,1) = dx(i,1) * u(1,j,1,e) &
                         + dx(i,2) * u(2,j,1,e) &
                         + dx(i,3) * u(3,j,1,e) &
                         + dx(i,4) * u(4,j,1,e) &
                         + dx(i,5) * u(5,j,1,e) &
                         + dx(i,6) * u(6,j,1,e) &
                         + dx(i,7) * u(7,j,1,e) &
                         + dx(i,8) * u(8,j,1,e) &
                         + dx(i,9) * u(9,j,1,e) &
                         + dx(i,10) * u(10,j,1,e)
          end do
       end do

       do k = 1, lx
          do j = 1, lx
             do i = 1, lx
                duds(i,j,k) = dy(j,1) * u(i,1,k,e) &
                            + dy(j,2) * u(i,2,k,e) &
                            + dy(j,3) * u(i,3,k,e) &
                            + dy(j,4) * u(i,4,k,e) &
                            + dy(j,5) * u(i,5,k,e) &
                            + dy(j,6) * u(i,6,k,e) &
                            + dy(j,7) * u(i,7,k,e) &
                            + dy(j,8) * u(i,8,k,e) &
                            + dy(j,9) * u(i,9,k,e) &
                            + dy(j,10) * u(i,10,k,e)
             end do
          end do
       end do

       do k = 1, lx
          do i = 1, lx*lx
             dudt(i,1,k) = dz(k,1) * u(i,1,1,e) &
                         + dz(k,2) * u(i,1,2,e) &
                         + dz(k,3) * u(i,1,3,e) &
                         + dz(k,4) * u(i,1,4,e) &
                         + dz(k,5) * u(i,1,5,e) &
                         + dz(k,6) * u(i,1,6,e) &
                         + dz(k,7) * u(i,1,7,e) &
                         + dz(k,8) * u(i,1,8,e) &
                         + dz(k,9) * u(i,1,9,e) &
                         + dz(k,10) * u(i,1,10,e)
          end do
       end do

       do i = 1, lx * lx * lx
          du(i,1,1,e) = jacinv(i,1,1,e) &
                       * ( vx(i,1,1,e) &
                         * ( drdx(i,1,1,e) * dudr(i,1,1) &
                           + dsdx(i,1,1,e) * duds(i,1,1) &
                           + dtdx(i,1,1,e) * dudt(i,1,1) ) &
                         + vy(i,1,1,e) &
                         * ( drdy(i,1,1,e) * dudr(i,1,1) &
                           + dsdy(i,1,1,e) * duds(i,1,1) &
                           + dtdy(i,1,1,e) * dudt(i,1,1) ) &
                         + vz(i,1,1,e) &
                         * ( drdz(i,1,1,e) * dudr(i,1,1) &
                           + dsdz(i,1,1,e) * duds(i,1,1) &
                           + dtdz(i,1,1,e) * dudt(i,1,1) ) )
       end do
    end do

  end subroutine cpu_conv1_lx10

  subroutine cpu_conv1_lx9(du, u, vx, vy, vz, dx, dy, dz, &
       drdx, dsdx, dtdx, drdy, dsdy, dtdy, drdz, dsdz, dtdz, &
       jacinv, nelv)
    integer, parameter :: lx = 9
    integer, intent(in) :: nelv
    real(kind=rp), dimension(lx, lx, lx, nelv), intent(inout) ::  du
    real(kind=rp), dimension(lx, lx, lx, nelv), intent(in) ::  u, vx, vy, vz
    real(kind=rp), dimension(lx, lx, lx, nelv), intent(in) :: drdx, dsdx, dtdx
    real(kind=rp), dimension(lx, lx, lx, nelv), intent(in) :: drdy, dsdy, dtdy
    real(kind=rp), dimension(lx, lx, lx, nelv), intent(in) :: drdz, dsdz, dtdz
    real(kind=rp), dimension(lx, lx, lx, nelv), intent(in) :: jacinv
    real(kind=rp), dimension(lx, lx), intent(in) :: dx, dy, dz
    real(kind=rp), dimension(lx, lx, lx) ::  dudr
    real(kind=rp), dimension(lx, lx, lx) ::  duds
    real(kind=rp), dimension(lx, lx, lx) ::  dudt
    integer :: e, i, j, k

    do e = 1, nelv
       do j = 1, lx * lx
          do i = 1, lx
             dudr(i,j,1) = dx(i,1) * u(1,j,1,e) &
                         + dx(i,2) * u(2,j,1,e) &
                         + dx(i,3) * u(3,j,1,e) &
                         + dx(i,4) * u(4,j,1,e) &
                         + dx(i,5) * u(5,j,1,e) &
                         + dx(i,6) * u(6,j,1,e) &
                         + dx(i,7) * u(7,j,1,e) &
                         + dx(i,8) * u(8,j,1,e) &
                         + dx(i,9) * u(9,j,1,e)
          end do
       end do

       do k = 1, lx
          do j = 1, lx
             do i = 1, lx
                duds(i,j,k) = dy(j,1) * u(i,1,k,e) &
                            + dy(j,2) * u(i,2,k,e) &
                            + dy(j,3) * u(i,3,k,e) &
                            + dy(j,4) * u(i,4,k,e) &
                            + dy(j,5) * u(i,5,k,e) &
                            + dy(j,6) * u(i,6,k,e) &
                            + dy(j,7) * u(i,7,k,e) &
                            + dy(j,8) * u(i,8,k,e) &
                            + dy(j,9) * u(i,9,k,e)
             end do
          end do
       end do

       do k = 1, lx
          do i = 1, lx*lx
             dudt(i,1,k) = dz(k,1) * u(i,1,1,e) &
                         + dz(k,2) * u(i,1,2,e) &
                         + dz(k,3) * u(i,1,3,e) &
                         + dz(k,4) * u(i,1,4,e) &
                         + dz(k,5) * u(i,1,5,e) &
                         + dz(k,6) * u(i,1,6,e) &
                         + dz(k,7) * u(i,1,7,e) &
                         + dz(k,8) * u(i,1,8,e) &
                         + dz(k,9) * u(i,1,9,e)
          end do
       end do

       do i = 1, lx * lx * lx
          du(i,1,1,e) = jacinv(i,1,1,e) &
                       * ( vx(i,1,1,e) &
                         * ( drdx(i,1,1,e) * dudr(i,1,1) &
                           + dsdx(i,1,1,e) * duds(i,1,1) &
                           + dtdx(i,1,1,e) * dudt(i,1,1) ) &
                         + vy(i,1,1,e) &
                         * ( drdy(i,1,1,e) * dudr(i,1,1) &
                           + dsdy(i,1,1,e) * duds(i,1,1) &
                           + dtdy(i,1,1,e) * dudt(i,1,1) ) &
                         + vz(i,1,1,e) &
                         * ( drdz(i,1,1,e) * dudr(i,1,1) &
                           + dsdz(i,1,1,e) * duds(i,1,1) &
                           + dtdz(i,1,1,e) * dudt(i,1,1) ) )
       end do
    end do

  end subroutine cpu_conv1_lx9

  subroutine cpu_conv1_lx8(du, u, vx, vy, vz, dx, dy, dz, &
       drdx, dsdx, dtdx, drdy, dsdy, dtdy, drdz, dsdz, dtdz, &
       jacinv, nelv)
    integer, parameter :: lx = 8
    integer, intent(in) :: nelv
    real(kind=rp), dimension(lx, lx, lx, nelv), intent(inout) ::  du
    real(kind=rp), dimension(lx, lx, lx, nelv), intent(in) ::  u, vx, vy, vz
    real(kind=rp), dimension(lx, lx, lx, nelv), intent(in) :: drdx, dsdx, dtdx
    real(kind=rp), dimension(lx, lx, lx, nelv), intent(in) :: drdy, dsdy, dtdy
    real(kind=rp), dimension(lx, lx, lx, nelv), intent(in) :: drdz, dsdz, dtdz
    real(kind=rp), dimension(lx, lx, lx, nelv), intent(in) :: jacinv
    real(kind=rp), dimension(lx, lx), intent(in) :: dx, dy, dz
    real(kind=rp), dimension(lx, lx, lx) ::  dudr
    real(kind=rp), dimension(lx, lx, lx) ::  duds
    real(kind=rp), dimension(lx, lx, lx) ::  dudt
    integer :: e, i, j, k

    do e = 1, nelv
       do j = 1, lx * lx
          do i = 1, lx
             dudr(i,j,1) = dx(i,1) * u(1,j,1,e) &
                         + dx(i,2) * u(2,j,1,e) &
                         + dx(i,3) * u(3,j,1,e) &
                         + dx(i,4) * u(4,j,1,e) &
                         + dx(i,5) * u(5,j,1,e) &
                         + dx(i,6) * u(6,j,1,e) &
                         + dx(i,7) * u(7,j,1,e) &
                         + dx(i,8) * u(8,j,1,e)
          end do
       end do

       do k = 1, lx
          do j = 1, lx
             do i = 1, lx
                duds(i,j,k) = dy(j,1) * u(i,1,k,e) &
                            + dy(j,2) * u(i,2,k,e) &
                            + dy(j,3) * u(i,3,k,e) &
                            + dy(j,4) * u(i,4,k,e) &
                            + dy(j,5) * u(i,5,k,e) &
                            + dy(j,6) * u(i,6,k,e) &
                            + dy(j,7) * u(i,7,k,e) &
                            + dy(j,8) * u(i,8,k,e)
             end do
          end do
       end do

       do k = 1, lx
          do i = 1, lx*lx
             dudt(i,1,k) = dz(k,1) * u(i,1,1,e) &
                         + dz(k,2) * u(i,1,2,e) &
                         + dz(k,3) * u(i,1,3,e) &
                         + dz(k,4) * u(i,1,4,e) &
                         + dz(k,5) * u(i,1,5,e) &
                         + dz(k,6) * u(i,1,6,e) &
                         + dz(k,7) * u(i,1,7,e) &
                         + dz(k,8) * u(i,1,8,e)
          end do
       end do

       do i = 1, lx * lx * lx
          du(i,1,1,e) = jacinv(i,1,1,e) &
                       * ( vx(i,1,1,e) &
                         * ( drdx(i,1,1,e) * dudr(i,1,1) &
                           + dsdx(i,1,1,e) * duds(i,1,1) &
                           + dtdx(i,1,1,e) * dudt(i,1,1) ) &
                         + vy(i,1,1,e) &
                         * ( drdy(i,1,1,e) * dudr(i,1,1) &
                           + dsdy(i,1,1,e) * duds(i,1,1) &
                           + dtdy(i,1,1,e) * dudt(i,1,1) ) &
                         + vz(i,1,1,e) &
                         * ( drdz(i,1,1,e) * dudr(i,1,1) &
                           + dsdz(i,1,1,e) * duds(i,1,1) &
                           + dtdz(i,1,1,e) * dudt(i,1,1) ) )
       end do
    end do

  end subroutine cpu_conv1_lx8

  subroutine cpu_conv1_lx7(du, u, vx, vy, vz, dx, dy, dz, &
       drdx, dsdx, dtdx, drdy, dsdy, dtdy, drdz, dsdz, dtdz, &
       jacinv, nelv)
    integer, parameter :: lx = 7
    integer, intent(in) :: nelv
    real(kind=rp), dimension(lx, lx, lx, nelv), intent(inout) ::  du
    real(kind=rp), dimension(lx, lx, lx, nelv), intent(in) ::  u, vx, vy, vz
    real(kind=rp), dimension(lx, lx, lx, nelv), intent(in) :: drdx, dsdx, dtdx
    real(kind=rp), dimension(lx, lx, lx, nelv), intent(in) :: drdy, dsdy, dtdy
    real(kind=rp), dimension(lx, lx, lx, nelv), intent(in) :: drdz, dsdz, dtdz
    real(kind=rp), dimension(lx, lx, lx, nelv), intent(in) :: jacinv
    real(kind=rp), dimension(lx, lx), intent(in) :: dx, dy, dz
    real(kind=rp), dimension(lx, lx, lx) ::  dudr
    real(kind=rp), dimension(lx, lx, lx) ::  duds
    real(kind=rp), dimension(lx, lx, lx) ::  dudt
    integer :: e, i, j, k

    do e = 1, nelv
       do j = 1, lx * lx
          do i = 1, lx
             dudr(i,j,1) = dx(i,1) * u(1,j,1,e) &
                         + dx(i,2) * u(2,j,1,e) &
                         + dx(i,3) * u(3,j,1,e) &
                         + dx(i,4) * u(4,j,1,e) &
                         + dx(i,5) * u(5,j,1,e) &
                         + dx(i,6) * u(6,j,1,e) &
                         + dx(i,7) * u(7,j,1,e)
          end do
       end do

       do k = 1, lx
          do j = 1, lx
             do i = 1, lx
                duds(i,j,k) = dy(j,1) * u(i,1,k,e) &
                            + dy(j,2) * u(i,2,k,e) &
                            + dy(j,3) * u(i,3,k,e) &
                            + dy(j,4) * u(i,4,k,e) &
                            + dy(j,5) * u(i,5,k,e) &
                            + dy(j,6) * u(i,6,k,e) &
                            + dy(j,7) * u(i,7,k,e)
             end do
          end do
       end do

       do k = 1, lx
          do i = 1, lx*lx
             dudt(i,1,k) = dz(k,1) * u(i,1,1,e) &
                         + dz(k,2) * u(i,1,2,e) &
                         + dz(k,3) * u(i,1,3,e) &
                         + dz(k,4) * u(i,1,4,e) &
                         + dz(k,5) * u(i,1,5,e) &
                         + dz(k,6) * u(i,1,6,e) &
                         + dz(k,7) * u(i,1,7,e)
          end do
       end do

       do i = 1, lx * lx * lx
          du(i,1,1,e) = jacinv(i,1,1,e) &
                       * ( vx(i,1,1,e) &
                         * ( drdx(i,1,1,e) * dudr(i,1,1) &
                           + dsdx(i,1,1,e) * duds(i,1,1) &
                           + dtdx(i,1,1,e) * dudt(i,1,1) ) &
                         + vy(i,1,1,e) &
                         * ( drdy(i,1,1,e) * dudr(i,1,1) &
                           + dsdy(i,1,1,e) * duds(i,1,1) &
                           + dtdy(i,1,1,e) * dudt(i,1,1) ) &
                         + vz(i,1,1,e) &
                         * ( drdz(i,1,1,e) * dudr(i,1,1) &
                           + dsdz(i,1,1,e) * duds(i,1,1) &
                           + dtdz(i,1,1,e) * dudt(i,1,1) ) )
       end do
    end do

  end subroutine cpu_conv1_lx7

  subroutine cpu_conv1_lx6(du, u, vx, vy, vz, dx, dy, dz, &
       drdx, dsdx, dtdx, drdy, dsdy, dtdy, drdz, dsdz, dtdz, &
       jacinv, nelv)
    integer, parameter :: lx = 6
    integer, intent(in) :: nelv
    real(kind=rp), dimension(lx, lx, lx, nelv), intent(inout) ::  du
    real(kind=rp), dimension(lx, lx, lx, nelv), intent(in) ::  u, vx, vy, vz
    real(kind=rp), dimension(lx, lx, lx, nelv), intent(in) :: drdx, dsdx, dtdx
    real(kind=rp), dimension(lx, lx, lx, nelv), intent(in) :: drdy, dsdy, dtdy
    real(kind=rp), dimension(lx, lx, lx, nelv), intent(in) :: drdz, dsdz, dtdz
    real(kind=rp), dimension(lx, lx, lx, nelv), intent(in) :: jacinv
    real(kind=rp), dimension(lx, lx), intent(in) :: dx, dy, dz
    real(kind=rp), dimension(lx, lx, lx) ::  dudr
    real(kind=rp), dimension(lx, lx, lx) ::  duds
    real(kind=rp), dimension(lx, lx, lx) ::  dudt
    integer :: e, i, j, k

    do e = 1, nelv
       do j = 1, lx * lx
          do i = 1, lx
             dudr(i,j,1) = dx(i,1) * u(1,j,1,e) &
                         + dx(i,2) * u(2,j,1,e) &
                         + dx(i,3) * u(3,j,1,e) &
                         + dx(i,4) * u(4,j,1,e) &
                         + dx(i,5) * u(5,j,1,e) &
                         + dx(i,6) * u(6,j,1,e)
          end do
       end do

       do k = 1, lx
          do j = 1, lx
             do i = 1, lx
                duds(i,j,k) = dy(j,1) * u(i,1,k,e) &
                            + dy(j,2) * u(i,2,k,e) &
                            + dy(j,3) * u(i,3,k,e) &
                            + dy(j,4) * u(i,4,k,e) &
                            + dy(j,5) * u(i,5,k,e) &
                            + dy(j,6) * u(i,6,k,e)
             end do
          end do
       end do

       do k = 1, lx
          do i = 1, lx*lx
             dudt(i,1,k) = dz(k,1) * u(i,1,1,e) &
                         + dz(k,2) * u(i,1,2,e) &
                         + dz(k,3) * u(i,1,3,e) &
                         + dz(k,4) * u(i,1,4,e) &
                         + dz(k,5) * u(i,1,5,e) &
                         + dz(k,6) * u(i,1,6,e)
          end do
       end do

       do i = 1, lx * lx * lx
          du(i,1,1,e) = jacinv(i,1,1,e) &
                       * ( vx(i,1,1,e) &
                         * ( drdx(i,1,1,e) * dudr(i,1,1) &
                           + dsdx(i,1,1,e) * duds(i,1,1) &
                           + dtdx(i,1,1,e) * dudt(i,1,1) ) &
                         + vy(i,1,1,e) &
                         * ( drdy(i,1,1,e) * dudr(i,1,1) &
                           + dsdy(i,1,1,e) * duds(i,1,1) &
                           + dtdy(i,1,1,e) * dudt(i,1,1) ) &
                         + vz(i,1,1,e) &
                         * ( drdz(i,1,1,e) * dudr(i,1,1) &
                           + dsdz(i,1,1,e) * duds(i,1,1) &
                           + dtdz(i,1,1,e) * dudt(i,1,1) ) )
       end do
    end do

  end subroutine cpu_conv1_lx6

  subroutine cpu_conv1_lx5(du, u, vx, vy, vz, dx, dy, dz, &
       drdx, dsdx, dtdx, drdy, dsdy, dtdy, drdz, dsdz, dtdz, &
       jacinv, nelv)
    integer, parameter :: lx = 5
    integer, intent(in) :: nelv
    real(kind=rp), dimension(lx, lx, lx, nelv), intent(inout) ::  du
    real(kind=rp), dimension(lx, lx, lx, nelv), intent(in) ::  u, vx, vy, vz
    real(kind=rp), dimension(lx, lx, lx, nelv), intent(in) :: drdx, dsdx, dtdx
    real(kind=rp), dimension(lx, lx, lx, nelv), intent(in) :: drdy, dsdy, dtdy
    real(kind=rp), dimension(lx, lx, lx, nelv), intent(in) :: drdz, dsdz, dtdz
    real(kind=rp), dimension(lx, lx, lx, nelv), intent(in) :: jacinv
    real(kind=rp), dimension(lx, lx), intent(in) :: dx, dy, dz
    real(kind=rp), dimension(lx, lx, lx) ::  dudr
    real(kind=rp), dimension(lx, lx, lx) ::  duds
    real(kind=rp), dimension(lx, lx, lx) ::  dudt
    integer :: e, i, j, k

    do e = 1, nelv
       do j = 1, lx * lx
          do i = 1, lx
             dudr(i,j,1) = dx(i,1) * u(1,j,1,e) &
                         + dx(i,2) * u(2,j,1,e) &
                         + dx(i,3) * u(3,j,1,e) &
                         + dx(i,4) * u(4,j,1,e) &
                         + dx(i,5) * u(5,j,1,e)
          end do
       end do

       do k = 1, lx
          do j = 1, lx
             do i = 1, lx
                duds(i,j,k) = dy(j,1) * u(i,1,k,e) &
                            + dy(j,2) * u(i,2,k,e) &
                            + dy(j,3) * u(i,3,k,e) &
                            + dy(j,4) * u(i,4,k,e) &
                            + dy(j,5) * u(i,5,k,e)
             end do
          end do
       end do

       do k = 1, lx
          do i = 1, lx*lx
             dudt(i,1,k) = dz(k,1) * u(i,1,1,e) &
                         + dz(k,2) * u(i,1,2,e) &
                         + dz(k,3) * u(i,1,3,e) &
                         + dz(k,4) * u(i,1,4,e) &
                         + dz(k,5) * u(i,1,5,e)
          end do
       end do

       do i = 1, lx * lx * lx
          du(i,1,1,e) = jacinv(i,1,1,e) &
                       * ( vx(i,1,1,e) &
                         * ( drdx(i,1,1,e) * dudr(i,1,1) &
                           + dsdx(i,1,1,e) * duds(i,1,1) &
                           + dtdx(i,1,1,e) * dudt(i,1,1) ) &
                         + vy(i,1,1,e) &
                         * ( drdy(i,1,1,e) * dudr(i,1,1) &
                           + dsdy(i,1,1,e) * duds(i,1,1) &
                           + dtdy(i,1,1,e) * dudt(i,1,1) ) &
                         + vz(i,1,1,e) &
                         * ( drdz(i,1,1,e) * dudr(i,1,1) &
                           + dsdz(i,1,1,e) * duds(i,1,1) &
                           + dtdz(i,1,1,e) * dudt(i,1,1) ) )
       end do
    end do

  end subroutine cpu_conv1_lx5

  subroutine cpu_conv1_lx4(du, u, vx, vy, vz, dx, dy, dz, &
       drdx, dsdx, dtdx, drdy, dsdy, dtdy, drdz, dsdz, dtdz, &
       jacinv, nelv)
    integer, parameter :: lx = 4
    integer, intent(in) :: nelv
    real(kind=rp), dimension(lx, lx, lx, nelv), intent(inout) ::  du
    real(kind=rp), dimension(lx, lx, lx, nelv), intent(in) ::  u, vx, vy, vz
    real(kind=rp), dimension(lx, lx, lx, nelv), intent(in) :: drdx, dsdx, dtdx
    real(kind=rp), dimension(lx, lx, lx, nelv), intent(in) :: drdy, dsdy, dtdy
    real(kind=rp), dimension(lx, lx, lx, nelv), intent(in) :: drdz, dsdz, dtdz
    real(kind=rp), dimension(lx, lx, lx, nelv), intent(in) :: jacinv
    real(kind=rp), dimension(lx, lx), intent(in) :: dx, dy, dz
    real(kind=rp), dimension(lx, lx, lx) ::  dudr
    real(kind=rp), dimension(lx, lx, lx) ::  duds
    real(kind=rp), dimension(lx, lx, lx) ::  dudt
    integer :: e, i, j, k

    do e = 1, nelv
       do j = 1, lx * lx
          do i = 1, lx
             dudr(i,j,1) = dx(i,1) * u(1,j,1,e) &
                         + dx(i,2) * u(2,j,1,e) &
                         + dx(i,3) * u(3,j,1,e) &
                         + dx(i,4) * u(4,j,1,e)
          end do
       end do

       do k = 1, lx
          do j = 1, lx
             do i = 1, lx
                duds(i,j,k) = dy(j,1) * u(i,1,k,e) &
                            + dy(j,2) * u(i,2,k,e) &
                            + dy(j,3) * u(i,3,k,e) &
                            + dy(j,4) * u(i,4,k,e)
             end do
          end do
       end do

       do k = 1, lx
          do i = 1, lx*lx
             dudt(i,1,k) = dz(k,1) * u(i,1,1,e) &
                         + dz(k,2) * u(i,1,2,e) &
                         + dz(k,3) * u(i,1,3,e) &
                         + dz(k,4) * u(i,1,4,e)
          end do
       end do

       do i = 1, lx * lx * lx
          du(i,1,1,e) = jacinv(i,1,1,e) &
                       * ( vx(i,1,1,e) &
                         * ( drdx(i,1,1,e) * dudr(i,1,1) &
                           + dsdx(i,1,1,e) * duds(i,1,1) &
                           + dtdx(i,1,1,e) * dudt(i,1,1) ) &
                         + vy(i,1,1,e) &
                         * ( drdy(i,1,1,e) * dudr(i,1,1) &
                           + dsdy(i,1,1,e) * duds(i,1,1) &
                           + dtdy(i,1,1,e) * dudt(i,1,1) ) &
                         + vz(i,1,1,e) &
                         * ( drdz(i,1,1,e) * dudr(i,1,1) &
                           + dsdz(i,1,1,e) * duds(i,1,1) &
                           + dtdz(i,1,1,e) * dudt(i,1,1) ) )
       end do
    end do

  end subroutine cpu_conv1_lx4

  subroutine cpu_conv1_lx3(du, u, vx, vy, vz, dx, dy, dz, &
       drdx, dsdx, dtdx, drdy, dsdy, dtdy, drdz, dsdz, dtdz, &
       jacinv, nelv)
    integer, parameter :: lx = 3
    integer, intent(in) :: nelv
    real(kind=rp), dimension(lx, lx, lx, nelv), intent(inout) ::  du
    real(kind=rp), dimension(lx, lx, lx, nelv), intent(in) ::  u, vx, vy, vz
    real(kind=rp), dimension(lx, lx, lx, nelv), intent(in) :: drdx, dsdx, dtdx
    real(kind=rp), dimension(lx, lx, lx, nelv), intent(in) :: drdy, dsdy, dtdy
    real(kind=rp), dimension(lx, lx, lx, nelv), intent(in) :: drdz, dsdz, dtdz
    real(kind=rp), dimension(lx, lx, lx, nelv), intent(in) :: jacinv
    real(kind=rp), dimension(lx, lx), intent(in) :: dx, dy, dz
    real(kind=rp), dimension(lx, lx, lx) ::  dudr
    real(kind=rp), dimension(lx, lx, lx) ::  duds
    real(kind=rp), dimension(lx, lx, lx) ::  dudt
    integer :: e, i, j, k

    do e = 1, nelv
       do j = 1, lx * lx
          do i = 1, lx
             dudr(i,j,1) = dx(i,1) * u(1,j,1,e) &
                         + dx(i,2) * u(2,j,1,e) &
                         + dx(i,3) * u(3,j,1,e)
          end do
       end do

       do k = 1, lx
          do j = 1, lx
             do i = 1, lx
                duds(i,j,k) = dy(j,1) * u(i,1,k,e) &
                            + dy(j,2) * u(i,2,k,e) &
                            + dy(j,3) * u(i,3,k,e)
             end do
          end do
       end do

       do k = 1, lx
          do i = 1, lx*lx
             dudt(i,1,k) = dz(k,1) * u(i,1,1,e) &
                         + dz(k,2) * u(i,1,2,e) &
                         + dz(k,3) * u(i,1,3,e)
          end do
       end do

       do i = 1, lx * lx * lx
          du(i,1,1,e) = jacinv(i,1,1,e) &
                       * ( vx(i,1,1,e) &
                         * ( drdx(i,1,1,e) * dudr(i,1,1) &
                           + dsdx(i,1,1,e) * duds(i,1,1) &
                           + dtdx(i,1,1,e) * dudt(i,1,1) ) &
                         + vy(i,1,1,e) &
                         * ( drdy(i,1,1,e) * dudr(i,1,1) &
                           + dsdy(i,1,1,e) * duds(i,1,1) &
                           + dtdy(i,1,1,e) * dudt(i,1,1) ) &
                         + vz(i,1,1,e) &
                         * ( drdz(i,1,1,e) * dudr(i,1,1) &
                           + dsdz(i,1,1,e) * duds(i,1,1) &
                           + dtdz(i,1,1,e) * dudt(i,1,1) ) )
       end do
    end do

  end subroutine cpu_conv1_lx3

  subroutine cpu_conv1_lx2(du, u, vx, vy, vz, dx, dy, dz, &
       drdx, dsdx, dtdx, drdy, dsdy, dtdy, drdz, dsdz, dtdz, &
       jacinv, nelv)
    integer, parameter :: lx = 2
    integer, intent(in) :: nelv
    real(kind=rp), dimension(lx, lx, lx, nelv), intent(inout) ::  du
    real(kind=rp), dimension(lx, lx, lx, nelv), intent(in) ::  u, vx, vy, vz
    real(kind=rp), dimension(lx, lx, lx, nelv), intent(in) :: drdx, dsdx, dtdx
    real(kind=rp), dimension(lx, lx, lx, nelv), intent(in) :: drdy, dsdy, dtdy
    real(kind=rp), dimension(lx, lx, lx, nelv), intent(in) :: drdz, dsdz, dtdz
    real(kind=rp), dimension(lx, lx, lx, nelv), intent(in) :: jacinv
    real(kind=rp), dimension(lx, lx), intent(in) :: dx, dy, dz
    real(kind=rp), dimension(lx, lx, lx) ::  dudr
    real(kind=rp), dimension(lx, lx, lx) ::  duds
    real(kind=rp), dimension(lx, lx, lx) ::  dudt
    integer :: e, i, j, k

    do e = 1, nelv
       do j = 1, lx * lx
          do i = 1, lx
             dudr(i,j,1) = dx(i,1) * u(1,j,1,e) &
                         + dx(i,2) * u(2,j,1,e)
          end do
       end do

       do k = 1, lx
          do j = 1, lx
             do i = 1, lx
                duds(i,j,k) = dy(j,1) * u(i,1,k,e) &
                            + dy(j,2) * u(i,2,k,e)
             end do
          end do
       end do

       do k = 1, lx
          do i = 1, lx*lx
             dudt(i,1,k) = dz(k,1) * u(i,1,1,e) &
                         + dz(k,2) * u(i,1,2,e)
          end do
       end do

       do i = 1, lx * lx * lx
          du(i,1,1,e) = jacinv(i,1,1,e) &
                       * ( vx(i,1,1,e) &
                         * ( drdx(i,1,1,e) * dudr(i,1,1) &
                           + dsdx(i,1,1,e) * duds(i,1,1) &
                           + dtdx(i,1,1,e) * dudt(i,1,1) ) &
                         + vy(i,1,1,e) &
                         * ( drdy(i,1,1,e) * dudr(i,1,1) &
                           + dsdy(i,1,1,e) * duds(i,1,1) &
                           + dtdy(i,1,1,e) * dudt(i,1,1) ) &
                         + vz(i,1,1,e) &
                         * ( drdz(i,1,1,e) * dudr(i,1,1) &
                           + dsdz(i,1,1,e) * duds(i,1,1) &
                           + dtdz(i,1,1,e) * dudt(i,1,1) ) )
       end do
    end do

  end subroutine cpu_conv1_lx2

  subroutine opr_cpu_conv1_single(du, u, vx, vy, vz, Xh, coef, e)
    integer, parameter :: e_len = 1
    type(space_t), intent(in) :: Xh
    type(coef_t), intent(in) :: coef
    integer, intent(in) :: e
    real(kind=rp), intent(inout) ::  du(Xh%lxyz, e_len)
    real(kind=rp), intent(inout) ::  u(Xh%lx, Xh%ly, Xh%lz, e_len)
    real(kind=rp), intent(inout) ::  vx(Xh%lx, Xh%ly, Xh%lz, e_len)
    real(kind=rp), intent(inout) ::  vy(Xh%lx, Xh%ly, Xh%lz, e_len)
    real(kind=rp), intent(inout) ::  vz(Xh%lx, Xh%ly, Xh%lz, e_len)
    
    associate(drdx => coef%drdx, drdy => coef%drdy, drdz => coef%drdz, &
         dsdx => coef%dsdx, dsdy => coef%dsdy, dsdz => coef%dsdz, &
         dtdx => coef%dtdx, dtdy => coef%dtdy, dtdz => coef%dtdz, &
         jacinv => coef%jacinv)
      select case (Xh%lx)
      case (14)
         call cpu_conv1_lx14_single(du, u, vx, vy, vz, Xh%dx, Xh%dy, Xh%dz, &
              drdx(1,1,1,e), dsdx(1,1,1,e), dtdx(1,1,1,e),&
              drdy(1,1,1,e), dsdy(1,1,1,e), dtdy(1,1,1,e),&
              drdz(1,1,1,e), dsdz(1,1,1,e), dtdz(1,1,1,e),&
              jacinv(1,1,1,e))
      case (13)
         call cpu_conv1_lx13_single(du, u, vx, vy, vz, Xh%dx, Xh%dy, Xh%dz, &
              drdx(1,1,1,e), dsdx(1,1,1,e), dtdx(1,1,1,e),&
              drdy(1,1,1,e), dsdy(1,1,1,e), dtdy(1,1,1,e),&
              drdz(1,1,1,e), dsdz(1,1,1,e), dtdz(1,1,1,e),&
              jacinv(1,1,1,e))
      case (12)
         call cpu_conv1_lx12_single(du, u, vx, vy, vz, Xh%dx, Xh%dy, Xh%dz, &
              drdx(1,1,1,e), dsdx(1,1,1,e), dtdx(1,1,1,e),&
              drdy(1,1,1,e), dsdy(1,1,1,e), dtdy(1,1,1,e),&
              drdz(1,1,1,e), dsdz(1,1,1,e), dtdz(1,1,1,e),&
              jacinv(1,1,1,e))
      case (11)
         call cpu_conv1_lx11_single(du, u, vx, vy, vz, Xh%dx, Xh%dy, Xh%dz, &
              drdx(1,1,1,e), dsdx(1,1,1,e), dtdx(1,1,1,e),&
              drdy(1,1,1,e), dsdy(1,1,1,e), dtdy(1,1,1,e),&
              drdz(1,1,1,e), dsdz(1,1,1,e), dtdz(1,1,1,e),&
              jacinv(1,1,1,e))
      case (10)
         call cpu_conv1_lx10_single(du, u, vx, vy, vz, Xh%dx, Xh%dy, Xh%dz, &
              drdx(1,1,1,e), dsdx(1,1,1,e), dtdx(1,1,1,e),&
              drdy(1,1,1,e), dsdy(1,1,1,e), dtdy(1,1,1,e),&
              drdz(1,1,1,e), dsdz(1,1,1,e), dtdz(1,1,1,e),&
              jacinv(1,1,1,e))
      case (9)
         call cpu_conv1_lx9_single(du, u, vx, vy, vz, Xh%dx, Xh%dy, Xh%dz, &
              drdx(1,1,1,e), dsdx(1,1,1,e), dtdx(1,1,1,e),&
              drdy(1,1,1,e), dsdy(1,1,1,e), dtdy(1,1,1,e),&
              drdz(1,1,1,e), dsdz(1,1,1,e), dtdz(1,1,1,e),&
              jacinv(1,1,1,e))
      case (8)
         call cpu_conv1_lx8_single(du, u, vx, vy, vz, Xh%dx, Xh%dy, Xh%dz, &
              drdx(1,1,1,e), dsdx(1,1,1,e), dtdx(1,1,1,e),&
              drdy(1,1,1,e), dsdy(1,1,1,e), dtdy(1,1,1,e),&
              drdz(1,1,1,e), dsdz(1,1,1,e), dtdz(1,1,1,e),&
              jacinv(1,1,1,e))
      case (7)
         call cpu_conv1_lx7_single(du, u, vx, vy, vz, Xh%dx, Xh%dy, Xh%dz, &
              drdx(1,1,1,e), dsdx(1,1,1,e), dtdx(1,1,1,e),&
              drdy(1,1,1,e), dsdy(1,1,1,e), dtdy(1,1,1,e),&
              drdz(1,1,1,e), dsdz(1,1,1,e), dtdz(1,1,1,e),&
              jacinv(1,1,1,e))
      case (6)
         call cpu_conv1_lx6_single(du, u, vx, vy, vz, Xh%dx, Xh%dy, Xh%dz, &
              drdx(1,1,1,e), dsdx(1,1,1,e), dtdx(1,1,1,e),&
              drdy(1,1,1,e), dsdy(1,1,1,e), dtdy(1,1,1,e),&
              drdz(1,1,1,e), dsdz(1,1,1,e), dtdz(1,1,1,e),&
              jacinv(1,1,1,e))
      case (5)
         call cpu_conv1_lx5_single(du, u, vx, vy, vz, Xh%dx, Xh%dy, Xh%dz, &
              drdx(1,1,1,e), dsdx(1,1,1,e), dtdx(1,1,1,e),&
              drdy(1,1,1,e), dsdy(1,1,1,e), dtdy(1,1,1,e),&
              drdz(1,1,1,e), dsdz(1,1,1,e), dtdz(1,1,1,e),&
              jacinv(1,1,1,e))
      case (4)
         call cpu_conv1_lx4_single(du, u, vx, vy, vz, Xh%dx, Xh%dy, Xh%dz, &
              drdx(1,1,1,e), dsdx(1,1,1,e), dtdx(1,1,1,e),&
              drdy(1,1,1,e), dsdy(1,1,1,e), dtdy(1,1,1,e),&
              drdz(1,1,1,e), dsdz(1,1,1,e), dtdz(1,1,1,e),&
              jacinv(1,1,1,e))
      case (3)
         call cpu_conv1_lx3_single(du, u, vx, vy, vz, Xh%dx, Xh%dy, Xh%dz, &
              drdx(1,1,1,e), dsdx(1,1,1,e), dtdx(1,1,1,e),&
              drdy(1,1,1,e), dsdy(1,1,1,e), dtdy(1,1,1,e),&
              drdz(1,1,1,e), dsdz(1,1,1,e), dtdz(1,1,1,e),&
              jacinv(1,1,1,e))
      case (2)
         call cpu_conv1_lx2_single(du, u, vx, vy, vz, Xh%dx, Xh%dy, Xh%dz, &
              drdx(1,1,1,e), dsdx(1,1,1,e), dtdx(1,1,1,e),&
              drdy(1,1,1,e), dsdy(1,1,1,e), dtdy(1,1,1,e),&
              drdz(1,1,1,e), dsdz(1,1,1,e), dtdz(1,1,1,e),&
              jacinv(1,1,1,e))
      case default
         call cpu_conv1_lx_single(du, u, vx, vy, vz, Xh%dx, Xh%dy, Xh%dz, &
              drdx(1,1,1,e), dsdx(1,1,1,e), dtdx(1,1,1,e),&
              drdy(1,1,1,e), dsdy(1,1,1,e), dtdy(1,1,1,e),&
              drdz(1,1,1,e), dsdz(1,1,1,e), dtdz(1,1,1,e),&
              jacinv(1,1,1,e), Xh%lx)
      end select
    end associate

  end subroutine opr_cpu_conv1_single

  subroutine cpu_conv1_lx_single(du, u, vx, vy, vz, dx, dy, dz, &
       drdx, dsdx, dtdx, drdy, dsdy, dtdy, drdz, dsdz, dtdz, &
       jacinv, lx)
    integer, intent(in) :: lx
    real(kind=rp), dimension(lx, lx, lx), intent(inout) ::  du
    real(kind=rp), dimension(lx, lx, lx), intent(in) ::  u, vx, vy, vz
    real(kind=rp), dimension(lx, lx, lx), intent(in) :: drdx, dsdx, dtdx
    real(kind=rp), dimension(lx, lx, lx), intent(in) :: drdy, dsdy, dtdy
    real(kind=rp), dimension(lx, lx, lx), intent(in) :: drdz, dsdz, dtdz
    real(kind=rp), dimension(lx, lx, lx), intent(in) :: jacinv
    real(kind=rp), dimension(lx, lx), intent(in) :: dx, dy, dz
    real(kind=rp), dimension(lx, lx, lx) ::  dudr
    real(kind=rp), dimension(lx, lx, lx) ::  duds
    real(kind=rp), dimension(lx, lx, lx) ::  dudt
    real(kind=rp) :: tmp
    integer :: i, j, k, l

    do j = 1, lx * lx
       do i = 1, lx
          tmp = 0.0_rp
          do k = 1, lx
             tmp = tmp + dx(i,k) * u(k,j,1)
          end do
          dudr(i,j,1) = tmp
       end do
    end do

    do k = 1, lx
       do j = 1, lx
          do i = 1, lx
             tmp = 0.0_rp
             do l = 1, lx
                tmp = tmp + dy(j,l) * u(i,l,k)
             end do
             duds(i,j,k) = tmp
          end do
       end do
    end do
    
    do k = 1, lx
       do i = 1, lx*lx
          tmp = 0.0_rp
          do l = 1, lx
             tmp = tmp + dz(k,l) * u(i,1,l)
          end do
          dudt(i,1,k) = tmp
       end do
    end do
    
    do i = 1, lx * lx * lx
       du(i,1,1) = jacinv(i,1,1) &
                * ( vx(i,1,1) &
                  * ( drdx(i,1,1) * dudr(i,1,1) &
                    + dsdx(i,1,1) * duds(i,1,1) &
                    + dtdx(i,1,1) * dudt(i,1,1) ) &
                  + vy(i,1,1) &
                  * ( drdy(i,1,1) * dudr(i,1,1) &
                    + dsdy(i,1,1) * duds(i,1,1) &
                    + dtdy(i,1,1) * dudt(i,1,1) ) &
                  + vz(i,1,1) &
                  * ( drdz(i,1,1) * dudr(i,1,1) &
                    + dsdz(i,1,1) * duds(i,1,1) &
                    + dtdz(i,1,1) * dudt(i,1,1) ) )
    end do
    
  end subroutine cpu_conv1_lx_single

  subroutine cpu_conv1_lx14_single(du, u, vx, vy, vz, dx, dy, dz, &
       drdx, dsdx, dtdx, drdy, dsdy, dtdy, drdz, dsdz, dtdz, jacinv)
    integer, parameter :: lx = 14
    real(kind=rp), dimension(lx, lx, lx), intent(inout) ::  du
    real(kind=rp), dimension(lx, lx, lx), intent(in) ::  u, vx, vy, vz
    real(kind=rp), dimension(lx, lx, lx), intent(in) :: drdx, dsdx, dtdx
    real(kind=rp), dimension(lx, lx, lx), intent(in) :: drdy, dsdy, dtdy
    real(kind=rp), dimension(lx, lx, lx), intent(in) :: drdz, dsdz, dtdz
    real(kind=rp), dimension(lx, lx, lx), intent(in) :: jacinv
    real(kind=rp), dimension(lx, lx), intent(in) :: dx, dy, dz
    real(kind=rp), dimension(lx, lx, lx) ::  dudr
    real(kind=rp), dimension(lx, lx, lx) ::  duds
    real(kind=rp), dimension(lx, lx, lx) ::  dudt
    integer :: i, j, k

    do j = 1, lx * lx
       do i = 1, lx
          dudr(i,j,1) = dx(i,1) * u(1,j,1) &
                      + dx(i,2) * u(2,j,1) &
                      + dx(i,3) * u(3,j,1) &
                      + dx(i,4) * u(4,j,1) &
                      + dx(i,5) * u(5,j,1) &
                      + dx(i,6) * u(6,j,1) &
                      + dx(i,7) * u(7,j,1) &
                      + dx(i,8) * u(8,j,1) &
                      + dx(i,9) * u(9,j,1) &
                      + dx(i,10) * u(10,j,1) &
                      + dx(i,11) * u(11,j,1) &
                      + dx(i,12) * u(12,j,1) &
                      + dx(i,13) * u(13,j,1) &
                      + dx(i,14) * u(14,j,1)
       end do
    end do

    do k = 1, lx
       do j = 1, lx
          do i = 1, lx
             duds(i,j,k) = dy(j,1) * u(i,1,k) &
                         + dy(j,2) * u(i,2,k) &
                         + dy(j,3) * u(i,3,k) &
                         + dy(j,4) * u(i,4,k) &
                         + dy(j,5) * u(i,5,k) &
                         + dy(j,6) * u(i,6,k) &
                         + dy(j,7) * u(i,7,k) &
                         + dy(j,8) * u(i,8,k) &
                         + dy(j,9) * u(i,9,k) &
                         + dy(j,10) * u(i,10,k) &
                         + dy(j,11) * u(i,11,k) &
                         + dy(j,12) * u(i,12,k) &
                         + dy(j,13) * u(i,13,k) &
                         + dy(j,14) * u(i,14,k)
          end do
       end do
    end do

    do k = 1, lx
       do i = 1, lx*lx
          dudt(i,1,k) = dz(k,1) * u(i,1,1) &
                      + dz(k,2) * u(i,1,2) &
                      + dz(k,3) * u(i,1,3) &
                      + dz(k,4) * u(i,1,4) &
                      + dz(k,5) * u(i,1,5) &
                      + dz(k,6) * u(i,1,6) &
                      + dz(k,7) * u(i,1,7) &
                      + dz(k,8) * u(i,1,8) &
                      + dz(k,9) * u(i,1,9) &
                      + dz(k,10) * u(i,1,10) &
                      + dz(k,11) * u(i,1,11) &
                      + dz(k,12) * u(i,1,12) &
                      + dz(k,13) * u(i,1,13) &
                      + dz(k,14) * u(i,1,14)
       end do
    end do
    
    do i = 1, lx * lx * lx
       du(i,1,1) = jacinv(i,1,1) &
                   * ( vx(i,1,1) &
                     * ( drdx(i,1,1) * dudr(i,1,1) &
                       + dsdx(i,1,1) * duds(i,1,1) &
                       + dtdx(i,1,1) * dudt(i,1,1) ) &
                     + vy(i,1,1) &
                     * ( drdy(i,1,1) * dudr(i,1,1) &
                       + dsdy(i,1,1) * duds(i,1,1) &
                       + dtdy(i,1,1) * dudt(i,1,1) ) &
                     + vz(i,1,1) &
                     * ( drdz(i,1,1) * dudr(i,1,1) &
                       + dsdz(i,1,1) * duds(i,1,1) &
                       + dtdz(i,1,1) * dudt(i,1,1) ) )
    end do
  
  end subroutine cpu_conv1_lx14_single

  subroutine cpu_conv1_lx13_single(du, u, vx, vy, vz, dx, dy, dz, &
       drdx, dsdx, dtdx, drdy, dsdy, dtdy, drdz, dsdz, dtdz, jacinv)
    integer, parameter :: lx = 13
    real(kind=rp), dimension(lx, lx, lx), intent(inout) ::  du
    real(kind=rp), dimension(lx, lx, lx), intent(in) ::  u, vx, vy, vz
    real(kind=rp), dimension(lx, lx, lx), intent(in) :: drdx, dsdx, dtdx
    real(kind=rp), dimension(lx, lx, lx), intent(in) :: drdy, dsdy, dtdy
    real(kind=rp), dimension(lx, lx, lx), intent(in) :: drdz, dsdz, dtdz
    real(kind=rp), dimension(lx, lx, lx), intent(in) :: jacinv
    real(kind=rp), dimension(lx, lx), intent(in) :: dx, dy, dz
    real(kind=rp), dimension(lx, lx, lx) ::  dudr
    real(kind=rp), dimension(lx, lx, lx) ::  duds
    real(kind=rp), dimension(lx, lx, lx) ::  dudt
    integer :: i, j, k

    do j = 1, lx * lx
       do i = 1, lx
          dudr(i,j,1) = dx(i,1) * u(1,j,1) &
                      + dx(i,2) * u(2,j,1) &
                      + dx(i,3) * u(3,j,1) &
                      + dx(i,4) * u(4,j,1) &
                      + dx(i,5) * u(5,j,1) &
                      + dx(i,6) * u(6,j,1) &
                      + dx(i,7) * u(7,j,1) &
                      + dx(i,8) * u(8,j,1) &
                      + dx(i,9) * u(9,j,1) &
                      + dx(i,10) * u(10,j,1) &
                      + dx(i,11) * u(11,j,1) &
                      + dx(i,12) * u(12,j,1) &
                      + dx(i,13) * u(13,j,1)
       end do
    end do
    
    do k = 1, lx
       do j = 1, lx
          do i = 1, lx
             duds(i,j,k) = dy(j,1) * u(i,1,k) &
                         + dy(j,2) * u(i,2,k) &
                         + dy(j,3) * u(i,3,k) &
                         + dy(j,4) * u(i,4,k) &
                         + dy(j,5) * u(i,5,k) &
                         + dy(j,6) * u(i,6,k) &
                         + dy(j,7) * u(i,7,k) &
                         + dy(j,8) * u(i,8,k) &
                         + dy(j,9) * u(i,9,k) &
                         + dy(j,10) * u(i,10,k) &
                         + dy(j,11) * u(i,11,k) &
                         + dy(j,12) * u(i,12,k) &
                         + dy(j,13) * u(i,13,k)
          end do
       end do
    end do
    
    do k = 1, lx
       do i = 1, lx*lx
          dudt(i,1,k) = dz(k,1) * u(i,1,1) &
                      + dz(k,2) * u(i,1,2) &
                      + dz(k,3) * u(i,1,3) &
                      + dz(k,4) * u(i,1,4) &
                      + dz(k,5) * u(i,1,5) &
                      + dz(k,6) * u(i,1,6) &
                      + dz(k,7) * u(i,1,7) &
                      + dz(k,8) * u(i,1,8) &
                      + dz(k,9) * u(i,1,9) &
                      + dz(k,10) * u(i,1,10) &
                      + dz(k,11) * u(i,1,11) &
                      + dz(k,12) * u(i,1,12) &
                      + dz(k,13) * u(i,1,13)
       end do
    end do

    do i = 1, lx * lx * lx
       du(i,1,1) = jacinv(i,1,1) &
                  * ( vx(i,1,1) &
                    * ( drdx(i,1,1) * dudr(i,1,1) &
                      + dsdx(i,1,1) * duds(i,1,1) &
                      + dtdx(i,1,1) * dudt(i,1,1) ) &
                      + vy(i,1,1) &
                    * ( drdy(i,1,1) * dudr(i,1,1) &
                      + dsdy(i,1,1) * duds(i,1,1) &
                      + dtdy(i,1,1) * dudt(i,1,1) ) &
                    + vz(i,1,1) &
                    * ( drdz(i,1,1) * dudr(i,1,1) &
                      + dsdz(i,1,1) * duds(i,1,1) &
                      + dtdz(i,1,1) * dudt(i,1,1) ) )
    end do

  end subroutine cpu_conv1_lx13_single

  
  subroutine cpu_conv1_lx12_single(du, u, vx, vy, vz, dx, dy, dz, &
       drdx, dsdx, dtdx, drdy, dsdy, dtdy, drdz, dsdz, dtdz, jacinv)
    integer, parameter :: lx = 12
    real(kind=rp), dimension(lx, lx, lx), intent(inout) ::  du
    real(kind=rp), dimension(lx, lx, lx), intent(in) ::  u, vx, vy, vz
    real(kind=rp), dimension(lx, lx, lx), intent(in) :: drdx, dsdx, dtdx
    real(kind=rp), dimension(lx, lx, lx), intent(in) :: drdy, dsdy, dtdy
    real(kind=rp), dimension(lx, lx, lx), intent(in) :: drdz, dsdz, dtdz
    real(kind=rp), dimension(lx, lx, lx), intent(in) :: jacinv
    real(kind=rp), dimension(lx, lx), intent(in) :: dx, dy, dz
    real(kind=rp), dimension(lx, lx, lx) ::  dudr
    real(kind=rp), dimension(lx, lx, lx) ::  duds
    real(kind=rp), dimension(lx, lx, lx) ::  dudt
    integer :: i, j, k

    do j = 1, lx * lx
       do i = 1, lx
          dudr(i,j,1) = dx(i,1) * u(1,j,1) &
                      + dx(i,2) * u(2,j,1) &
                      + dx(i,3) * u(3,j,1) &
                      + dx(i,4) * u(4,j,1) &
                      + dx(i,5) * u(5,j,1) &
                      + dx(i,6) * u(6,j,1) &
                      + dx(i,7) * u(7,j,1) &
                      + dx(i,8) * u(8,j,1) &
                      + dx(i,9) * u(9,j,1) &
                      + dx(i,10) * u(10,j,1) &
                      + dx(i,11) * u(11,j,1) &
                      + dx(i,12) * u(12,j,1)
       end do
    end do

    do k = 1, lx
       do j = 1, lx
          do i = 1, lx
             duds(i,j,k) = dy(j,1) * u(i,1,k) &
                         + dy(j,2) * u(i,2,k) &
                         + dy(j,3) * u(i,3,k) &
                         + dy(j,4) * u(i,4,k) &
                         + dy(j,5) * u(i,5,k) &
                         + dy(j,6) * u(i,6,k) &
                         + dy(j,7) * u(i,7,k) &
                         + dy(j,8) * u(i,8,k) &
                         + dy(j,9) * u(i,9,k) &
                         + dy(j,10) * u(i,10,k) &
                         + dy(j,11) * u(i,11,k) &
                         + dy(j,12) * u(i,12,k)
          end do
       end do
    end do

    do k = 1, lx
       do i = 1, lx*lx
          dudt(i,1,k) = dz(k,1) * u(i,1,1) &
                      + dz(k,2) * u(i,1,2) &
                      + dz(k,3) * u(i,1,3) &
                      + dz(k,4) * u(i,1,4) &
                      + dz(k,5) * u(i,1,5) &
                      + dz(k,6) * u(i,1,6) &
                      + dz(k,7) * u(i,1,7) &
                      + dz(k,8) * u(i,1,8) &
                      + dz(k,9) * u(i,1,9) &
                      + dz(k,10) * u(i,1,10) &
                      + dz(k,11) * u(i,1,11) &
                      + dz(k,12) * u(i,1,12)
       end do
    end do
    
    do i = 1, lx * lx * lx
       du(i,1,1) = jacinv(i,1,1) &
                    * ( vx(i,1,1) &
                      * ( drdx(i,1,1) * dudr(i,1,1) &
                        + dsdx(i,1,1) * duds(i,1,1) &
                        + dtdx(i,1,1) * dudt(i,1,1) ) &
                      + vy(i,1,1) &
                      * ( drdy(i,1,1) * dudr(i,1,1) &
                        + dsdy(i,1,1) * duds(i,1,1) &
                        + dtdy(i,1,1) * dudt(i,1,1) ) &
                      + vz(i,1,1) &
                      * ( drdz(i,1,1) * dudr(i,1,1) &
                        + dsdz(i,1,1) * duds(i,1,1) &
                        + dtdz(i,1,1) * dudt(i,1,1) ) )
    end do

  end subroutine cpu_conv1_lx12_single

  subroutine cpu_conv1_lx11_single(du, u, vx, vy, vz, dx, dy, dz, &
       drdx, dsdx, dtdx, drdy, dsdy, dtdy, drdz, dsdz, dtdz, jacinv)
    integer, parameter :: lx = 11
    real(kind=rp), dimension(lx, lx, lx), intent(inout) ::  du
    real(kind=rp), dimension(lx, lx, lx), intent(in) ::  u, vx, vy, vz
    real(kind=rp), dimension(lx, lx, lx), intent(in) :: drdx, dsdx, dtdx
    real(kind=rp), dimension(lx, lx, lx), intent(in) :: drdy, dsdy, dtdy
    real(kind=rp), dimension(lx, lx, lx), intent(in) :: drdz, dsdz, dtdz
    real(kind=rp), dimension(lx, lx, lx), intent(in) :: jacinv
    real(kind=rp), dimension(lx, lx), intent(in) :: dx, dy, dz
    real(kind=rp), dimension(lx, lx, lx) ::  dudr
    real(kind=rp), dimension(lx, lx, lx) ::  duds
    real(kind=rp), dimension(lx, lx, lx) ::  dudt
    integer :: i, j, k

    do j = 1, lx * lx
       do i = 1, lx
          dudr(i,j,1) = dx(i,1) * u(1,j,1) &
                      + dx(i,2) * u(2,j,1) &
                      + dx(i,3) * u(3,j,1) &
                      + dx(i,4) * u(4,j,1) &
                      + dx(i,5) * u(5,j,1) &
                      + dx(i,6) * u(6,j,1) &
                      + dx(i,7) * u(7,j,1) &
                      + dx(i,8) * u(8,j,1) &
                      + dx(i,9) * u(9,j,1) &
                      + dx(i,10) * u(10,j,1) &
                      + dx(i,11) * u(11,j,1)
       end do
    end do

    do k = 1, lx
       do j = 1, lx
          do i = 1, lx
             duds(i,j,k) = dy(j,1) * u(i,1,k) &
                         + dy(j,2) * u(i,2,k) &
                         + dy(j,3) * u(i,3,k) &
                         + dy(j,4) * u(i,4,k) &
                         + dy(j,5) * u(i,5,k) &
                         + dy(j,6) * u(i,6,k) &
                         + dy(j,7) * u(i,7,k) &
                         + dy(j,8) * u(i,8,k) &
                         + dy(j,9) * u(i,9,k) &
                         + dy(j,10) * u(i,10,k) &
                         + dy(j,11) * u(i,11,k)
          end do
       end do
    end do
    
    do k = 1, lx
       do i = 1, lx*lx
          dudt(i,1,k) = dz(k,1) * u(i,1,1) &
                      + dz(k,2) * u(i,1,2) &
                      + dz(k,3) * u(i,1,3) &
                      + dz(k,4) * u(i,1,4) &
                      + dz(k,5) * u(i,1,5) &
                      + dz(k,6) * u(i,1,6) &
                      + dz(k,7) * u(i,1,7) &
                      + dz(k,8) * u(i,1,8) &
                      + dz(k,9) * u(i,1,9) &
                      + dz(k,10) * u(i,1,10) &
                      + dz(k,11) * u(i,1,11)
       end do
    end do

    do i = 1, lx * lx * lx
       du(i,1,1) = jacinv(i,1,1) &
                    * ( vx(i,1,1) &
                      * ( drdx(i,1,1) * dudr(i,1,1) &
                         + dsdx(i,1,1) * duds(i,1,1) &
                         + dtdx(i,1,1) * dudt(i,1,1) ) &
                      + vy(i,1,1) &
                      * ( drdy(i,1,1) * dudr(i,1,1) &
                        + dsdy(i,1,1) * duds(i,1,1) &
                        + dtdy(i,1,1) * dudt(i,1,1) ) &
                      + vz(i,1,1) &
                      * ( drdz(i,1,1) * dudr(i,1,1) &
                        + dsdz(i,1,1) * duds(i,1,1) &
                        + dtdz(i,1,1) * dudt(i,1,1) ) )
    end do

  end subroutine cpu_conv1_lx11_single

  subroutine cpu_conv1_lx10_single(du, u, vx, vy, vz, dx, dy, dz, &
       drdx, dsdx, dtdx, drdy, dsdy, dtdy, drdz, dsdz, dtdz, jacinv)
    integer, parameter :: lx = 10
    real(kind=rp), dimension(lx, lx, lx), intent(inout) ::  du
    real(kind=rp), dimension(lx, lx, lx), intent(in) ::  u, vx, vy, vz
    real(kind=rp), dimension(lx, lx, lx), intent(in) :: drdx, dsdx, dtdx
    real(kind=rp), dimension(lx, lx, lx), intent(in) :: drdy, dsdy, dtdy
    real(kind=rp), dimension(lx, lx, lx), intent(in) :: drdz, dsdz, dtdz
    real(kind=rp), dimension(lx, lx, lx), intent(in) :: jacinv
    real(kind=rp), dimension(lx, lx), intent(in) :: dx, dy, dz
    real(kind=rp), dimension(lx, lx, lx) ::  dudr
    real(kind=rp), dimension(lx, lx, lx) ::  duds
    real(kind=rp), dimension(lx, lx, lx) ::  dudt
    integer :: i, j, k

    do j = 1, lx * lx
       do i = 1, lx
          dudr(i,j,1) = dx(i,1) * u(1,j,1) &
                      + dx(i,2) * u(2,j,1) &
                      + dx(i,3) * u(3,j,1) &
                      + dx(i,4) * u(4,j,1) &
                      + dx(i,5) * u(5,j,1) &
                      + dx(i,6) * u(6,j,1) &
                      + dx(i,7) * u(7,j,1) &
                      + dx(i,8) * u(8,j,1) &
                      + dx(i,9) * u(9,j,1) &
                      + dx(i,10) * u(10,j,1)
       end do
    end do
    
    do k = 1, lx
       do j = 1, lx
          do i = 1, lx
             duds(i,j,k) = dy(j,1) * u(i,1,k) &
                         + dy(j,2) * u(i,2,k) &
                         + dy(j,3) * u(i,3,k) &
                         + dy(j,4) * u(i,4,k) &
                         + dy(j,5) * u(i,5,k) &
                         + dy(j,6) * u(i,6,k) &
                         + dy(j,7) * u(i,7,k) &
                         + dy(j,8) * u(i,8,k) &
                         + dy(j,9) * u(i,9,k) &
                         + dy(j,10) * u(i,10,k)
          end do
       end do
    end do
    
    do k = 1, lx
       do i = 1, lx*lx
          dudt(i,1,k) = dz(k,1) * u(i,1,1) &
                      + dz(k,2) * u(i,1,2) &
                      + dz(k,3) * u(i,1,3) &
                      + dz(k,4) * u(i,1,4) &
                      + dz(k,5) * u(i,1,5) &
                      + dz(k,6) * u(i,1,6) &
                      + dz(k,7) * u(i,1,7) &
                      + dz(k,8) * u(i,1,8) &
                      + dz(k,9) * u(i,1,9) &
                      + dz(k,10) * u(i,1,10)
       end do
    end do
       
    do i = 1, lx * lx * lx
       du(i,1,1) = jacinv(i,1,1) &
                    * ( vx(i,1,1) &
                      * ( drdx(i,1,1) * dudr(i,1,1) &
                        + dsdx(i,1,1) * duds(i,1,1) &
                        + dtdx(i,1,1) * dudt(i,1,1) ) &
                      + vy(i,1,1) &
                      * ( drdy(i,1,1) * dudr(i,1,1) &
                        + dsdy(i,1,1) * duds(i,1,1) &
                        + dtdy(i,1,1) * dudt(i,1,1) ) &
                      + vz(i,1,1) &
                      * ( drdz(i,1,1) * dudr(i,1,1) &
                        + dsdz(i,1,1) * duds(i,1,1) &
                        + dtdz(i,1,1) * dudt(i,1,1) ) )
    end do

  end subroutine cpu_conv1_lx10_single

  subroutine cpu_conv1_lx9_single(du, u, vx, vy, vz, dx, dy, dz, &
       drdx, dsdx, dtdx, drdy, dsdy, dtdy, drdz, dsdz, dtdz, jacinv)
    integer, parameter :: lx = 9
    real(kind=rp), dimension(lx, lx, lx), intent(inout) ::  du
    real(kind=rp), dimension(lx, lx, lx), intent(in) ::  u, vx, vy, vz
    real(kind=rp), dimension(lx, lx, lx), intent(in) :: drdx, dsdx, dtdx
    real(kind=rp), dimension(lx, lx, lx), intent(in) :: drdy, dsdy, dtdy
    real(kind=rp), dimension(lx, lx, lx), intent(in) :: drdz, dsdz, dtdz
    real(kind=rp), dimension(lx, lx, lx), intent(in) :: jacinv
    real(kind=rp), dimension(lx, lx), intent(in) :: dx, dy, dz
    real(kind=rp), dimension(lx, lx, lx) ::  dudr
    real(kind=rp), dimension(lx, lx, lx) ::  duds
    real(kind=rp), dimension(lx, lx, lx) ::  dudt
    integer :: i, j, k

    do j = 1, lx * lx
       do i = 1, lx
          dudr(i,j,1) = dx(i,1) * u(1,j,1) &
                      + dx(i,2) * u(2,j,1) &
                      + dx(i,3) * u(3,j,1) &
                      + dx(i,4) * u(4,j,1) &
                      + dx(i,5) * u(5,j,1) &
                      + dx(i,6) * u(6,j,1) &
                      + dx(i,7) * u(7,j,1) &
                      + dx(i,8) * u(8,j,1) &
                      + dx(i,9) * u(9,j,1)
       end do
    end do
    
    do k = 1, lx
       do j = 1, lx
          do i = 1, lx
             duds(i,j,k) = dy(j,1) * u(i,1,k) &
                         + dy(j,2) * u(i,2,k) &
                         + dy(j,3) * u(i,3,k) &
                         + dy(j,4) * u(i,4,k) &
                         + dy(j,5) * u(i,5,k) &
                         + dy(j,6) * u(i,6,k) &
                         + dy(j,7) * u(i,7,k) &
                         + dy(j,8) * u(i,8,k) &
                         + dy(j,9) * u(i,9,k)
          end do
       end do
    end do
    
    do k = 1, lx
       do i = 1, lx*lx
          dudt(i,1,k) = dz(k,1) * u(i,1,1) &
                      + dz(k,2) * u(i,1,2) &
                      + dz(k,3) * u(i,1,3) &
                      + dz(k,4) * u(i,1,4) &
                      + dz(k,5) * u(i,1,5) &
                      + dz(k,6) * u(i,1,6) &
                      + dz(k,7) * u(i,1,7) &
                      + dz(k,8) * u(i,1,8) &
                      + dz(k,9) * u(i,1,9)
       end do
    end do
    
    do i = 1, lx * lx * lx
       du(i,1,1) = jacinv(i,1,1) &
                    * ( vx(i,1,1) &
                      * ( drdx(i,1,1) * dudr(i,1,1) &
                        + dsdx(i,1,1) * duds(i,1,1) &
                        + dtdx(i,1,1) * dudt(i,1,1) ) &
                      + vy(i,1,1) &
                      * ( drdy(i,1,1) * dudr(i,1,1) &
                        + dsdy(i,1,1) * duds(i,1,1) &
                        + dtdy(i,1,1) * dudt(i,1,1) ) &
                      + vz(i,1,1) &
                      * ( drdz(i,1,1) * dudr(i,1,1) &
                        + dsdz(i,1,1) * duds(i,1,1) &
                        + dtdz(i,1,1) * dudt(i,1,1) ) )
    end do
    
  end subroutine cpu_conv1_lx9_single

  subroutine cpu_conv1_lx8_single(du, u, vx, vy, vz, dx, dy, dz, &
       drdx, dsdx, dtdx, drdy, dsdy, dtdy, drdz, dsdz, dtdz, jacinv)
    integer, parameter :: lx = 8
    real(kind=rp), dimension(lx, lx, lx), intent(inout) ::  du
    real(kind=rp), dimension(lx, lx, lx), intent(in) ::  u, vx, vy, vz
    real(kind=rp), dimension(lx, lx, lx), intent(in) :: drdx, dsdx, dtdx
    real(kind=rp), dimension(lx, lx, lx), intent(in) :: drdy, dsdy, dtdy
    real(kind=rp), dimension(lx, lx, lx), intent(in) :: drdz, dsdz, dtdz
    real(kind=rp), dimension(lx, lx, lx), intent(in) :: jacinv
    real(kind=rp), dimension(lx, lx), intent(in) :: dx, dy, dz
    real(kind=rp), dimension(lx, lx, lx) ::  dudr
    real(kind=rp), dimension(lx, lx, lx) ::  duds
    real(kind=rp), dimension(lx, lx, lx) ::  dudt
    integer :: i, j, k

    do j = 1, lx * lx
       do i = 1, lx
          dudr(i,j,1) = dx(i,1) * u(1,j,1) &
                      + dx(i,2) * u(2,j,1) &
                      + dx(i,3) * u(3,j,1) &
                      + dx(i,4) * u(4,j,1) &
                      + dx(i,5) * u(5,j,1) &
                      + dx(i,6) * u(6,j,1) &
                      + dx(i,7) * u(7,j,1) &
                      + dx(i,8) * u(8,j,1)
       end do
    end do

    do k = 1, lx
       do j = 1, lx
          do i = 1, lx
             duds(i,j,k) = dy(j,1) * u(i,1,k) &
                         + dy(j,2) * u(i,2,k) &
                         + dy(j,3) * u(i,3,k) &
                         + dy(j,4) * u(i,4,k) &
                         + dy(j,5) * u(i,5,k) &
                         + dy(j,6) * u(i,6,k) &
                         + dy(j,7) * u(i,7,k) &
                         + dy(j,8) * u(i,8,k)
          end do
       end do
    end do
       
    do k = 1, lx
       do i = 1, lx*lx
          dudt(i,1,k) = dz(k,1) * u(i,1,1) &
                      + dz(k,2) * u(i,1,2) &
                      + dz(k,3) * u(i,1,3) &
                      + dz(k,4) * u(i,1,4) &
                      + dz(k,5) * u(i,1,5) &
                      + dz(k,6) * u(i,1,6) &
                      + dz(k,7) * u(i,1,7) &
                      + dz(k,8) * u(i,1,8)
       end do
    end do

    do i = 1, lx * lx * lx
       du(i,1,1) = jacinv(i,1,1) &
                    * ( vx(i,1,1) &
                      * ( drdx(i,1,1) * dudr(i,1,1) &
                        + dsdx(i,1,1) * duds(i,1,1) &
                        + dtdx(i,1,1) * dudt(i,1,1) ) &
                      + vy(i,1,1) &
                      * ( drdy(i,1,1) * dudr(i,1,1) &
                        + dsdy(i,1,1) * duds(i,1,1) &
                        + dtdy(i,1,1) * dudt(i,1,1) ) &
                      + vz(i,1,1) &
                      * ( drdz(i,1,1) * dudr(i,1,1) &
                        + dsdz(i,1,1) * duds(i,1,1) &
                        + dtdz(i,1,1) * dudt(i,1,1) ) )
    end do
    
  end subroutine cpu_conv1_lx8_single

  subroutine cpu_conv1_lx7_single(du, u, vx, vy, vz, dx, dy, dz, &
       drdx, dsdx, dtdx, drdy, dsdy, dtdy, drdz, dsdz, dtdz, jacinv)
    integer, parameter :: lx = 7
    real(kind=rp), dimension(lx, lx, lx), intent(inout) ::  du
    real(kind=rp), dimension(lx, lx, lx), intent(in) ::  u, vx, vy, vz
    real(kind=rp), dimension(lx, lx, lx), intent(in) :: drdx, dsdx, dtdx
    real(kind=rp), dimension(lx, lx, lx), intent(in) :: drdy, dsdy, dtdy
    real(kind=rp), dimension(lx, lx, lx), intent(in) :: drdz, dsdz, dtdz
    real(kind=rp), dimension(lx, lx, lx), intent(in) :: jacinv
    real(kind=rp), dimension(lx, lx), intent(in) :: dx, dy, dz
    real(kind=rp), dimension(lx, lx, lx) ::  dudr
    real(kind=rp), dimension(lx, lx, lx) ::  duds
    real(kind=rp), dimension(lx, lx, lx) ::  dudt
    integer :: i, j, k

    do j = 1, lx * lx
       do i = 1, lx
          dudr(i,j,1) = dx(i,1) * u(1,j,1) &
                      + dx(i,2) * u(2,j,1) &
                      + dx(i,3) * u(3,j,1) &
                      + dx(i,4) * u(4,j,1) &
                      + dx(i,5) * u(5,j,1) &
                      + dx(i,6) * u(6,j,1) &
                      + dx(i,7) * u(7,j,1)
       end do
    end do
    
    do k = 1, lx
       do j = 1, lx
          do i = 1, lx
             duds(i,j,k) = dy(j,1) * u(i,1,k) &
                         + dy(j,2) * u(i,2,k) &
                         + dy(j,3) * u(i,3,k) &
                         + dy(j,4) * u(i,4,k) &
                         + dy(j,5) * u(i,5,k) &
                         + dy(j,6) * u(i,6,k) &
                         + dy(j,7) * u(i,7,k)
          end do
       end do
    end do
    
    do k = 1, lx
       do i = 1, lx*lx
          dudt(i,1,k) = dz(k,1) * u(i,1,1) &
                      + dz(k,2) * u(i,1,2) &
                      + dz(k,3) * u(i,1,3) &
                      + dz(k,4) * u(i,1,4) &
                      + dz(k,5) * u(i,1,5) &
                      + dz(k,6) * u(i,1,6) &
                      + dz(k,7) * u(i,1,7)
       end do
    end do
    
    do i = 1, lx * lx * lx
       du(i,1,1) = jacinv(i,1,1) &
                    * ( vx(i,1,1) &
                      * ( drdx(i,1,1) * dudr(i,1,1) &
                        + dsdx(i,1,1) * duds(i,1,1) &
                        + dtdx(i,1,1) * dudt(i,1,1) ) &
                      + vy(i,1,1) &
                      * ( drdy(i,1,1) * dudr(i,1,1) &
                        + dsdy(i,1,1) * duds(i,1,1) &
                        + dtdy(i,1,1) * dudt(i,1,1) ) &
                      + vz(i,1,1) &
                      * ( drdz(i,1,1) * dudr(i,1,1) &
                        + dsdz(i,1,1) * duds(i,1,1) &
                        + dtdz(i,1,1) * dudt(i,1,1) ) )
    end do

  end subroutine cpu_conv1_lx7_single

  subroutine cpu_conv1_lx6_single(du, u, vx, vy, vz, dx, dy, dz, &
       drdx, dsdx, dtdx, drdy, dsdy, dtdy, drdz, dsdz, dtdz, jacinv)
    integer, parameter :: lx = 6
    real(kind=rp), dimension(lx, lx, lx), intent(inout) ::  du
    real(kind=rp), dimension(lx, lx, lx), intent(in) ::  u, vx, vy, vz
    real(kind=rp), dimension(lx, lx, lx), intent(in) :: drdx, dsdx, dtdx
    real(kind=rp), dimension(lx, lx, lx), intent(in) :: drdy, dsdy, dtdy
    real(kind=rp), dimension(lx, lx, lx), intent(in) :: drdz, dsdz, dtdz
    real(kind=rp), dimension(lx, lx, lx), intent(in) :: jacinv
    real(kind=rp), dimension(lx, lx), intent(in) :: dx, dy, dz
    real(kind=rp), dimension(lx, lx, lx) ::  dudr
    real(kind=rp), dimension(lx, lx, lx) ::  duds
    real(kind=rp), dimension(lx, lx, lx) ::  dudt
    integer :: i, j, k

    do j = 1, lx * lx
       do i = 1, lx
          dudr(i,j,1) = dx(i,1) * u(1,j,1) &
                      + dx(i,2) * u(2,j,1) &
                      + dx(i,3) * u(3,j,1) &
                      + dx(i,4) * u(4,j,1) &
                      + dx(i,5) * u(5,j,1) &
                      + dx(i,6) * u(6,j,1)
       end do
    end do
    
    do k = 1, lx
       do j = 1, lx
          do i = 1, lx
             duds(i,j,k) = dy(j,1) * u(i,1,k) &
                         + dy(j,2) * u(i,2,k) &
                         + dy(j,3) * u(i,3,k) &
                         + dy(j,4) * u(i,4,k) &
                         + dy(j,5) * u(i,5,k) &
                         + dy(j,6) * u(i,6,k)
          end do
       end do
    end do
    
    do k = 1, lx
       do i = 1, lx*lx
          dudt(i,1,k) = dz(k,1) * u(i,1,1) &
                      + dz(k,2) * u(i,1,2) &
                      + dz(k,3) * u(i,1,3) &
                      + dz(k,4) * u(i,1,4) &
                      + dz(k,5) * u(i,1,5) &
                      + dz(k,6) * u(i,1,6)
       end do
    end do
    
    do i = 1, lx * lx * lx
       du(i,1,1) = jacinv(i,1,1) &
                    * ( vx(i,1,1) &
                      * ( drdx(i,1,1) * dudr(i,1,1) &
                        + dsdx(i,1,1) * duds(i,1,1) &
                        + dtdx(i,1,1) * dudt(i,1,1) ) &
                      + vy(i,1,1) &
                      * ( drdy(i,1,1) * dudr(i,1,1) &
                        + dsdy(i,1,1) * duds(i,1,1) &
                        + dtdy(i,1,1) * dudt(i,1,1) ) &
                      + vz(i,1,1) &
                      * ( drdz(i,1,1) * dudr(i,1,1) &
                        + dsdz(i,1,1) * duds(i,1,1) &
                        + dtdz(i,1,1) * dudt(i,1,1) ) )
    end do

  end subroutine cpu_conv1_lx6_single

  subroutine cpu_conv1_lx5_single(du, u, vx, vy, vz, dx, dy, dz, &
       drdx, dsdx, dtdx, drdy, dsdy, dtdy, drdz, dsdz, dtdz, jacinv)
    integer, parameter :: lx = 5
    real(kind=rp), dimension(lx, lx, lx), intent(inout) ::  du
    real(kind=rp), dimension(lx, lx, lx), intent(in) ::  u, vx, vy, vz
    real(kind=rp), dimension(lx, lx, lx), intent(in) :: drdx, dsdx, dtdx
    real(kind=rp), dimension(lx, lx, lx), intent(in) :: drdy, dsdy, dtdy
    real(kind=rp), dimension(lx, lx, lx), intent(in) :: drdz, dsdz, dtdz
    real(kind=rp), dimension(lx, lx, lx), intent(in) :: jacinv
    real(kind=rp), dimension(lx, lx), intent(in) :: dx, dy, dz
    real(kind=rp), dimension(lx, lx, lx) ::  dudr
    real(kind=rp), dimension(lx, lx, lx) ::  duds
    real(kind=rp), dimension(lx, lx, lx) ::  dudt
    integer :: i, j, k

    do j = 1, lx * lx
       do i = 1, lx
          dudr(i,j,1) = dx(i,1) * u(1,j,1) &
                      + dx(i,2) * u(2,j,1) &
                      + dx(i,3) * u(3,j,1) &
                      + dx(i,4) * u(4,j,1) &
                      + dx(i,5) * u(5,j,1)
       end do
    end do

    do k = 1, lx
       do j = 1, lx
          do i = 1, lx
             duds(i,j,k) = dy(j,1) * u(i,1,k) &
                         + dy(j,2) * u(i,2,k) &
                         + dy(j,3) * u(i,3,k) &
                         + dy(j,4) * u(i,4,k) &
                         + dy(j,5) * u(i,5,k)
          end do
       end do
    end do
    
    do k = 1, lx
       do i = 1, lx*lx
          dudt(i,1,k) = dz(k,1) * u(i,1,1) &
                      + dz(k,2) * u(i,1,2) &
                      + dz(k,3) * u(i,1,3) &
                      + dz(k,4) * u(i,1,4) &
                      + dz(k,5) * u(i,1,5)
       end do
    end do
    
    do i = 1, lx * lx * lx
       du(i,1,1) = jacinv(i,1,1) &
                    * ( vx(i,1,1) &
                      * ( drdx(i,1,1) * dudr(i,1,1) &
                        + dsdx(i,1,1) * duds(i,1,1) &
                        + dtdx(i,1,1) * dudt(i,1,1) ) &
                      + vy(i,1,1) &
                      * ( drdy(i,1,1) * dudr(i,1,1) &
                        + dsdy(i,1,1) * duds(i,1,1) &
                        + dtdy(i,1,1) * dudt(i,1,1) ) &
                      + vz(i,1,1) &
                      * ( drdz(i,1,1) * dudr(i,1,1) &
                        + dsdz(i,1,1) * duds(i,1,1) &
                        + dtdz(i,1,1) * dudt(i,1,1) ) )
    end do

  end subroutine cpu_conv1_lx5_single

  subroutine cpu_conv1_lx4_single(du, u, vx, vy, vz, dx, dy, dz, &
       drdx, dsdx, dtdx, drdy, dsdy, dtdy, drdz, dsdz, dtdz, jacinv)
    integer, parameter :: lx = 4
    real(kind=rp), dimension(lx, lx, lx), intent(inout) ::  du
    real(kind=rp), dimension(lx, lx, lx), intent(in) ::  u, vx, vy, vz
    real(kind=rp), dimension(lx, lx, lx), intent(in) :: drdx, dsdx, dtdx
    real(kind=rp), dimension(lx, lx, lx), intent(in) :: drdy, dsdy, dtdy
    real(kind=rp), dimension(lx, lx, lx), intent(in) :: drdz, dsdz, dtdz
    real(kind=rp), dimension(lx, lx, lx), intent(in) :: jacinv
    real(kind=rp), dimension(lx, lx), intent(in) :: dx, dy, dz
    real(kind=rp), dimension(lx, lx, lx) ::  dudr
    real(kind=rp), dimension(lx, lx, lx) ::  duds
    real(kind=rp), dimension(lx, lx, lx) ::  dudt
    integer :: i, j, k

    do j = 1, lx * lx
       do i = 1, lx
          dudr(i,j,1) = dx(i,1) * u(1,j,1) &
                      + dx(i,2) * u(2,j,1) &
                      + dx(i,3) * u(3,j,1) &
                      + dx(i,4) * u(4,j,1)
       end do
    end do
    
    do k = 1, lx
       do j = 1, lx
          do i = 1, lx
             duds(i,j,k) = dy(j,1) * u(i,1,k) &
                         + dy(j,2) * u(i,2,k) &
                         + dy(j,3) * u(i,3,k) &
                         + dy(j,4) * u(i,4,k)
          end do
       end do
    end do
    
    do k = 1, lx
       do i = 1, lx*lx
          dudt(i,1,k) = dz(k,1) * u(i,1,1) &
                      + dz(k,2) * u(i,1,2) &
                      + dz(k,3) * u(i,1,3) &
                      + dz(k,4) * u(i,1,4)
       end do
    end do
    
    do i = 1, lx * lx * lx
       du(i,1,1) = jacinv(i,1,1) &
                    * ( vx(i,1,1) &
                      * ( drdx(i,1,1) * dudr(i,1,1) &
                        + dsdx(i,1,1) * duds(i,1,1) &
                        + dtdx(i,1,1) * dudt(i,1,1) ) &
                      + vy(i,1,1) &
                      * ( drdy(i,1,1) * dudr(i,1,1) &
                        + dsdy(i,1,1) * duds(i,1,1) &
                        + dtdy(i,1,1) * dudt(i,1,1) ) &
                      + vz(i,1,1) &
                      * ( drdz(i,1,1) * dudr(i,1,1) &
                        + dsdz(i,1,1) * duds(i,1,1) &
                        + dtdz(i,1,1) * dudt(i,1,1) ) )
    end do

  end subroutine cpu_conv1_lx4_single

  subroutine cpu_conv1_lx3_single(du, u, vx, vy, vz, dx, dy, dz, &
       drdx, dsdx, dtdx, drdy, dsdy, dtdy, drdz, dsdz, dtdz, jacinv)
    integer, parameter :: lx = 3
    real(kind=rp), dimension(lx, lx, lx), intent(inout) ::  du
    real(kind=rp), dimension(lx, lx, lx), intent(in) ::  u, vx, vy, vz
    real(kind=rp), dimension(lx, lx, lx), intent(in) :: drdx, dsdx, dtdx
    real(kind=rp), dimension(lx, lx, lx), intent(in) :: drdy, dsdy, dtdy
    real(kind=rp), dimension(lx, lx, lx), intent(in) :: drdz, dsdz, dtdz
    real(kind=rp), dimension(lx, lx, lx), intent(in) :: jacinv
    real(kind=rp), dimension(lx, lx), intent(in) :: dx, dy, dz
    real(kind=rp), dimension(lx, lx, lx) ::  dudr
    real(kind=rp), dimension(lx, lx, lx) ::  duds
    real(kind=rp), dimension(lx, lx, lx) ::  dudt
    integer :: i, j, k

    do j = 1, lx * lx
       do i = 1, lx
          dudr(i,j,1) = dx(i,1) * u(1,j,1) &
                      + dx(i,2) * u(2,j,1) &
                      + dx(i,3) * u(3,j,1)
       end do
    end do
    
    do k = 1, lx
       do j = 1, lx
          do i = 1, lx
             duds(i,j,k) = dy(j,1) * u(i,1,k) &
                         + dy(j,2) * u(i,2,k) &
                         + dy(j,3) * u(i,3,k)
          end do
       end do
    end do

    do k = 1, lx
       do i = 1, lx*lx
          dudt(i,1,k) = dz(k,1) * u(i,1,1) &
                      + dz(k,2) * u(i,1,2) &
                      + dz(k,3) * u(i,1,3)
       end do
    end do
    
    do i = 1, lx * lx * lx
       du(i,1,1) = jacinv(i,1,1) &
                    * ( vx(i,1,1) &
                      * ( drdx(i,1,1) * dudr(i,1,1) &
                        + dsdx(i,1,1) * duds(i,1,1) &
                        + dtdx(i,1,1) * dudt(i,1,1) ) &
                      + vy(i,1,1) &
                      * ( drdy(i,1,1) * dudr(i,1,1) &
                        + dsdy(i,1,1) * duds(i,1,1) &
                        + dtdy(i,1,1) * dudt(i,1,1) ) &
                      + vz(i,1,1) &
                      * ( drdz(i,1,1) * dudr(i,1,1) &
                        + dsdz(i,1,1) * duds(i,1,1) &
                        + dtdz(i,1,1) * dudt(i,1,1) ) )
    end do

  end subroutine cpu_conv1_lx3_single

  subroutine cpu_conv1_lx2_single(du, u, vx, vy, vz, dx, dy, dz, &
       drdx, dsdx, dtdx, drdy, dsdy, dtdy, drdz, dsdz, dtdz, jacinv)
    integer, parameter :: lx = 2
    real(kind=rp), dimension(lx, lx, lx), intent(inout) ::  du
    real(kind=rp), dimension(lx, lx, lx), intent(in) ::  u, vx, vy, vz
    real(kind=rp), dimension(lx, lx, lx), intent(in) :: drdx, dsdx, dtdx
    real(kind=rp), dimension(lx, lx, lx), intent(in) :: drdy, dsdy, dtdy
    real(kind=rp), dimension(lx, lx, lx), intent(in) :: drdz, dsdz, dtdz
    real(kind=rp), dimension(lx, lx, lx), intent(in) :: jacinv
    real(kind=rp), dimension(lx, lx), intent(in) :: dx, dy, dz
    real(kind=rp), dimension(lx, lx, lx) ::  dudr
    real(kind=rp), dimension(lx, lx, lx) ::  duds
    real(kind=rp), dimension(lx, lx, lx) ::  dudt
    integer :: i, j, k

    do j = 1, lx * lx
       do i = 1, lx
          dudr(i,j,1) = dx(i,1) * u(1,j,1) &
                      + dx(i,2) * u(2,j,1)
       end do
    end do
    
    do k = 1, lx
       do j = 1, lx
          do i = 1, lx
             duds(i,j,k) = dy(j,1) * u(i,1,k) &
                         + dy(j,2) * u(i,2,k)
          end do
       end do
    end do
    
    do k = 1, lx
       do i = 1, lx*lx
          dudt(i,1,k) = dz(k,1) * u(i,1,1) &
                      + dz(k,2) * u(i,1,2)
       end do
    end do
    
    do i = 1, lx * lx * lx
       du(i,1,1) = jacinv(i,1,1) &
                    * ( vx(i,1,1) &
                      * ( drdx(i,1,1) * dudr(i,1,1) &
                        + dsdx(i,1,1) * duds(i,1,1) &
                        + dtdx(i,1,1) * dudt(i,1,1) ) &
                      + vy(i,1,1) &
                      * ( drdy(i,1,1) * dudr(i,1,1) &
                        + dsdy(i,1,1) * duds(i,1,1) &
                        + dtdy(i,1,1) * dudt(i,1,1) ) &
                      + vz(i,1,1) &
                      * ( drdz(i,1,1) * dudr(i,1,1) &
                        + dsdz(i,1,1) * duds(i,1,1) &
                        + dtdz(i,1,1) * dudt(i,1,1) ) )
    end do

  end subroutine cpu_conv1_lx2_single

end submodule cpu_conv1
