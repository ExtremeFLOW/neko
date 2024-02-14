! Copyright (c) 2021, The Neko Authors
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
!> Jacobi preconditioner SX-Aurora backend
module sx_jacobi
  use math
  use precon
  use coefs
  use dofmap
  use num_types
  use gather_scatter
  implicit none
  private

  !> Defines a jacobi preconditioner for SX-Aurora
  type, public, extends(pc_t) :: sx_jacobi_t
     real(kind=rp), allocatable :: d(:,:,:,:)
     type(gs_t), pointer :: gs_h
     type(dofmap_t), pointer :: dof
     type(coef_t), pointer :: coef
   contains
     procedure, pass(this) :: init => sx_jacobi_init
     procedure, pass(this) :: free => sx_jacobi_free
     procedure, pass(this) :: solve => sx_jacobi_solve
     procedure, pass(this) :: update => sx_jacobi_update
  end type sx_jacobi_t

contains

  subroutine sx_jacobi_init(this, coef, dof, gs_h)
    class(sx_jacobi_t), intent(inout) :: this
    type(coef_t), intent(inout), target :: coef
    type(dofmap_t), intent(inout), target :: dof
    type(gs_t), intent(inout), target :: gs_h

    call this%free()
    this%gs_h => gs_h
    this%dof => dof
    this%coef => coef
    allocate(this%d(dof%Xh%lx,dof%Xh%ly,dof%Xh%lz, dof%msh%nelv))
    call sx_jacobi_update(this)

  end subroutine sx_jacobi_init

  subroutine sx_jacobi_free(this)
    class(sx_jacobi_t), intent(inout) :: this
    if (allocated(this%d)) then
       deallocate(this%d)
    end if
    nullify(this%dof)
    nullify(this%gs_h)
    nullify(this%coef)
  end subroutine sx_jacobi_free

  !> The jacobi preconditioner \f$ J z = r \f$
  !! \f$ z = J^{-1}r\f$ where \f$ J^{-1} ~= 1/diag(A) \f$
  subroutine sx_jacobi_solve(this, z, r, n)
    integer, intent(in) :: n
    class(sx_jacobi_t), intent(inout) :: this
    real(kind=rp), dimension(n), intent(inout) :: z
    real(kind=rp), dimension(n), intent(inout) :: r
    call col3(z,r,this%d,n)
  end subroutine sx_jacobi_solve

  subroutine sx_jacobi_update(this)
    class(sx_jacobi_t), intent(inout) :: this
    integer :: lz, ly, lx
    associate(dof => this%dof, coef => this%coef, &
         gs_h => this%gs_h, nelv => this%dof%msh%nelv)

      lx = dof%Xh%lx
      ly = dof%Xh%ly
      lz = dof%Xh%lz

      this%d = 0d0

      select case(lx)
      case (14)
         call sx_update_lx13(this%d, coef%Xh%dxt, coef%Xh%dyt, coef%Xh%dzt, &
              coef%G11, coef%G22, coef%G33, coef%G12, coef%G13, coef%G23, nelv)
      case (13)
         call sx_update_lx13(this%d, coef%Xh%dxt, coef%Xh%dyt, coef%Xh%dzt, &
              coef%G11, coef%G22, coef%G33, coef%G12, coef%G13, coef%G23, nelv)
      case (12)
         call sx_update_lx12(this%d, coef%Xh%dxt, coef%Xh%dyt, coef%Xh%dzt, &
              coef%G11, coef%G22, coef%G33, coef%G12, coef%G13, coef%G23, nelv)
      case (11)
         call sx_update_lx11(this%d, coef%Xh%dxt, coef%Xh%dyt, coef%Xh%dzt, &
              coef%G11, coef%G22, coef%G33, coef%G12, coef%G13, coef%G23, nelv)
      case (10)
         call sx_update_lx10(this%d, coef%Xh%dxt, coef%Xh%dyt, coef%Xh%dzt, &
              coef%G11, coef%G22, coef%G33, coef%G12, coef%G13, coef%G23, nelv)
      case (9)
         call sx_update_lx9(this%d, coef%Xh%dxt, coef%Xh%dyt, coef%Xh%dzt, &
              coef%G11, coef%G22, coef%G33, coef%G12, coef%G13, coef%G23, nelv)
      case (8)
         call sx_update_lx8(this%d, coef%Xh%dxt, coef%Xh%dyt, coef%Xh%dzt, &
              coef%G11, coef%G22, coef%G33, coef%G12, coef%G13, coef%G23, nelv)
      case (7)
         call sx_update_lx7(this%d, coef%Xh%dxt, coef%Xh%dyt, coef%Xh%dzt, &
              coef%G11, coef%G22, coef%G33, coef%G12, coef%G13, coef%G23, nelv)
      case (6)
         call sx_update_lx6(this%d, coef%Xh%dxt, coef%Xh%dyt, coef%Xh%dzt, &
              coef%G11, coef%G22, coef%G33, coef%G12, coef%G13, coef%G23, nelv)
      case (5)
         call sx_update_lx5(this%d, coef%Xh%dxt, coef%Xh%dyt, coef%Xh%dzt, &
              coef%G11, coef%G22, coef%G33, coef%G12, coef%G13, coef%G23, nelv)
      case (4)
         call sx_update_lx4(this%d, coef%Xh%dxt, coef%Xh%dyt, coef%Xh%dzt, &
              coef%G11, coef%G22, coef%G33, coef%G12, coef%G13, coef%G23, nelv)
      case (2)
         call sx_update_lx2(this%d, coef%Xh%dxt, coef%Xh%dyt, coef%Xh%dzt, &
              coef%G11, coef%G22, coef%G33, coef%G12, coef%G13, coef%G23, nelv)
      case default
         call sx_update_lx(this%d, coef%Xh%dxt, coef%Xh%dyt, coef%Xh%dzt, &
              coef%G11, coef%G22, coef%G33, coef%G12, coef%G13, coef%G23, &
              nelv, lx)
      end select

      call col2(this%d, coef%h1, coef%dof%size())
      if (coef%ifh2) call addcol3(this%d, coef%h2, coef%B, coef%dof%size())
      call gs_h%op(this%d, dof%size(), GS_OP_ADD)
      if (.not. coef%ifh2) call col2(this%d, coef%mult, coef%dof%size())
      call invcol1(this%d, dof%size())
    end associate
  end subroutine sx_jacobi_update

  subroutine sx_update_lx(d, dxt, dyt, dzt, G11, G22, G33, G12, G13, G23, n, lx)
    integer, intent(in) :: n, lx
    real(kind=rp), intent(inout) :: d(lx, lx, lx, n)
    real(kind=rp), intent(in) :: dxt(lx,lx)
    real(kind=rp), intent(in) :: dyt(lx,lx)
    real(kind=rp), intent(in) :: dzt(lx,lx)
    real(kind=rp), intent(in) :: G11(lx, lx, lx, n)
    real(kind=rp), intent(in) :: G22(lx, lx, lx, n)
    real(kind=rp), intent(in) :: G33(lx, lx, lx, n)
    real(kind=rp), intent(in) :: G12(lx, lx, lx, n)
    real(kind=rp), intent(in) :: G13(lx, lx, lx, n)
    real(kind=rp), intent(in) :: G23(lx, lx, lx, n)
    integer :: i, j, k, l, e


    do l = 1, lx
       do k = 1, lx
          do j = 1, lx
             do i = 1, lx
                do e = 1, n
                   d(i,j,k,e) = d(i,j,k,e) + &
                        G11(l,j,k,e) * dxt(i,l)**2

                   d(i,j,k,e) = d(i,j,k,e) + &
                        G22(i,l,k,e) * dyt(j,l)**2

                   d(i,j,k,e) = d(i,j,k,e) + &
                        G33(i,j,l,e) * dzt(k,l)**2

                end do
             end do
          end do
       end do
    end do

    do j = 1, lx, lx-1
       do k = 1, lx, lx-1
          do e = 1, n
             d(1,j,k,e) = d(1,j,k,e) &
                  + G12(1,j,k,e) * dxt(1,1)*dyt(j,j) &
                  + G13(1,j,k,e) * dxt(1,1)*dzt(k,k)
             d(lx,j,k,e) = d(lx,j,k,e) &
                  + G12(lx,j,k,e) * dxt(lx,lx)*dyt(j,j) &
                  + G13(lx,j,k,e) * dxt(lx,lx)*dzt(k,k)
          end do
       end do
    end do

    do i = 1, lx, lx-1
       do k = 1, lx, lx-1
          do e = 1, n
             d(i,1,k,e) = d(i,1,k,e) &
                  + G12(i,1,k,e) * dyt(1,1)*dxt(i,i) &
                  + G23(i,1,k,e) * dyt(1,1)*dzt(k,k)
             d(i,lx,k,e) = d(i,lx,k,e) &
                  + G12(i,lx,k,e) * dyt(lx,lx)*dxt(i,i) &
                  + G23(i,lx,k,e) * dyt(lx,lx)*dzt(k,k)
          end do
       end do
    end do

    do i = 1, lx, lx-1
       do j = 1, lx, lx-1
          do e = 1, n
             d(i,j,1,e) = d(i,j,1,e) &
                  + G13(i,j,1,e) * dzt(1,1)*dxt(i,i) &
                  + G23(i,j,1,e) * dzt(1,1)*dyt(j,j)
             d(i,j,lx,e) = d(i,j,lx,e) &
                  + G13(i,j,lx,e) * dzt(lx,lx)*dxt(i,i) &
                  + G23(i,j,lx,e) * dzt(lx,lx)*dyt(j,j)
          end do
       end do
    end do

  end subroutine sx_update_lx

  subroutine sx_update_lx14(d, dxt, dyt, dzt, G11, G22, G33, G12, G13, G23, n)
    integer, parameter :: lx = 14
    integer, parameter :: ly = 14
    integer, parameter :: lz = 14
    integer, intent(in) :: n
    real(kind=rp), intent(inout) :: d(lx, lx, lx, n)
    real(kind=rp), intent(in) :: dxt(lx,lx)
    real(kind=rp), intent(in) :: dyt(lx,lx)
    real(kind=rp), intent(in) :: dzt(lx,lx)
    real(kind=rp), intent(in) :: G11(lx, lx, lx, n)
    real(kind=rp), intent(in) :: G22(lx, lx, lx, n)
    real(kind=rp), intent(in) :: G33(lx, lx, lx, n)
    real(kind=rp), intent(in) :: G12(lx, lx, lx, n)
    real(kind=rp), intent(in) :: G13(lx, lx, lx, n)
    real(kind=rp), intent(in) :: G23(lx, lx, lx, n)
    integer :: i, j, k, l, e


    do l = 1, lx
       do k = 1, lz
          do j = 1, ly
             do i = 1, lx
                do e = 1, n
                   d(i,j,k,e) = d(i,j,k,e) + &
                        G11(l,j,k,e) * dxt(i,l)**2

                   d(i,j,k,e) = d(i,j,k,e) + &
                        G22(i,l,k,e) * dyt(j,l)**2

                   d(i,j,k,e) = d(i,j,k,e) + &
                        G33(i,j,l,e) * dzt(k,l)**2

                end do
             end do
          end do
       end do
    end do

    do j = 1, ly, ly-1
       do k = 1, lz, lz-1
          do e = 1, n
             d(1,j,k,e) = d(1,j,k,e) &
                  + G12(1,j,k,e) * dxt(1,1)*dyt(j,j) &
                  + G13(1,j,k,e) * dxt(1,1)*dzt(k,k)
             d(lx,j,k,e) = d(lx,j,k,e) &
                  + G12(lx,j,k,e) * dxt(lx,lx)*dyt(j,j) &
                  + G13(lx,j,k,e) * dxt(lx,lx)*dzt(k,k)
          end do
       end do
    end do

    do i = 1, lx, lx-1
       do k = 1, lz, lz-1
          do e = 1, n
             d(i,1,k,e) = d(i,1,k,e) &
                  + G12(i,1,k,e) * dyt(1,1)*dxt(i,i) &
                  + G23(i,1,k,e) * dyt(1,1)*dzt(k,k)
             d(i,ly,k,e) = d(i,ly,k,e) &
                  + G12(i,ly,k,e) * dyt(ly,ly)*dxt(i,i) &
                  + G23(i,ly,k,e) * dyt(ly,ly)*dzt(k,k)
          end do
       end do
    end do

    do i = 1, lx, lx-1
       do j = 1, ly, ly-1
          do e = 1, n
             d(i,j,1,e) = d(i,j,1,e) &
                  + G13(i,j,1,e) * dzt(1,1)*dxt(i,i) &
                  + G23(i,j,1,e) * dzt(1,1)*dyt(j,j)
             d(i,j,lz,e) = d(i,j,lz,e) &
                  + G13(i,j,lz,e) * dzt(lz,lz)*dxt(i,i) &
                  + G23(i,j,lz,e) * dzt(lz,lz)*dyt(j,j)
          end do
       end do
    end do

  end subroutine sx_update_lx14

  subroutine sx_update_lx13(d, dxt, dyt, dzt, G11, G22, G33, G12, G13, G23, n)
    integer, parameter :: lx = 13
    integer, parameter :: ly = 13
    integer, parameter :: lz = 13
    integer, intent(in) :: n
    real(kind=rp), intent(inout) :: d(lx, lx, lx, n)
    real(kind=rp), intent(in) :: dxt(lx,lx)
    real(kind=rp), intent(in) :: dyt(lx,lx)
    real(kind=rp), intent(in) :: dzt(lx,lx)
    real(kind=rp), intent(in) :: G11(lx, lx, lx, n)
    real(kind=rp), intent(in) :: G22(lx, lx, lx, n)
    real(kind=rp), intent(in) :: G33(lx, lx, lx, n)
    real(kind=rp), intent(in) :: G12(lx, lx, lx, n)
    real(kind=rp), intent(in) :: G13(lx, lx, lx, n)
    real(kind=rp), intent(in) :: G23(lx, lx, lx, n)
    integer :: i, j, k, l, e


    do l = 1, lx
       do k = 1, lz
          do j = 1, ly
             do i = 1, lx
                do e = 1, n
                   d(i,j,k,e) = d(i,j,k,e) + &
                        G11(l,j,k,e) * dxt(i,l)**2

                   d(i,j,k,e) = d(i,j,k,e) + &
                        G22(i,l,k,e) * dyt(j,l)**2

                   d(i,j,k,e) = d(i,j,k,e) + &
                        G33(i,j,l,e) * dzt(k,l)**2

                end do
             end do
          end do
       end do
    end do

    do j = 1, ly, ly-1
       do k = 1, lz, lz-1
          do e = 1, n
             d(1,j,k,e) = d(1,j,k,e) &
                  + G12(1,j,k,e) * dxt(1,1)*dyt(j,j) &
                  + G13(1,j,k,e) * dxt(1,1)*dzt(k,k)
             d(lx,j,k,e) = d(lx,j,k,e) &
                  + G12(lx,j,k,e) * dxt(lx,lx)*dyt(j,j) &
                  + G13(lx,j,k,e) * dxt(lx,lx)*dzt(k,k)
          end do
       end do
    end do

    do i = 1, lx, lx-1
       do k = 1, lz, lz-1
          do e = 1, n
             d(i,1,k,e) = d(i,1,k,e) &
                  + G12(i,1,k,e) * dyt(1,1)*dxt(i,i) &
                  + G23(i,1,k,e) * dyt(1,1)*dzt(k,k)
             d(i,ly,k,e) = d(i,ly,k,e) &
                  + G12(i,ly,k,e) * dyt(ly,ly)*dxt(i,i) &
                  + G23(i,ly,k,e) * dyt(ly,ly)*dzt(k,k)
          end do
       end do
    end do

    do i = 1, lx, lx-1
       do j = 1, ly, ly-1
          do e = 1, n
             d(i,j,1,e) = d(i,j,1,e) &
                  + G13(i,j,1,e) * dzt(1,1)*dxt(i,i) &
                  + G23(i,j,1,e) * dzt(1,1)*dyt(j,j)
             d(i,j,lz,e) = d(i,j,lz,e) &
                  + G13(i,j,lz,e) * dzt(lz,lz)*dxt(i,i) &
                  + G23(i,j,lz,e) * dzt(lz,lz)*dyt(j,j)
          end do
       end do
    end do

  end subroutine sx_update_lx13

  subroutine sx_update_lx12(d, dxt, dyt, dzt, G11, G22, G33, G12, G13, G23, n)
    integer, parameter :: lx = 12
    integer, parameter :: ly = 12
    integer, parameter :: lz = 12
    integer, intent(in) :: n
    real(kind=rp), intent(inout) :: d(lx, lx, lx, n)
    real(kind=rp), intent(in) :: dxt(lx,lx)
    real(kind=rp), intent(in) :: dyt(lx,lx)
    real(kind=rp), intent(in) :: dzt(lx,lx)
    real(kind=rp), intent(in) :: G11(lx, lx, lx, n)
    real(kind=rp), intent(in) :: G22(lx, lx, lx, n)
    real(kind=rp), intent(in) :: G33(lx, lx, lx, n)
    real(kind=rp), intent(in) :: G12(lx, lx, lx, n)
    real(kind=rp), intent(in) :: G13(lx, lx, lx, n)
    real(kind=rp), intent(in) :: G23(lx, lx, lx, n)
    integer :: i, j, k, l, e


    do l = 1, lx
       do k = 1, lz
          do j = 1, ly
             do i = 1, lx
                do e = 1, n
                   d(i,j,k,e) = d(i,j,k,e) + &
                        G11(l,j,k,e) * dxt(i,l)**2

                   d(i,j,k,e) = d(i,j,k,e) + &
                        G22(i,l,k,e) * dyt(j,l)**2

                   d(i,j,k,e) = d(i,j,k,e) + &
                        G33(i,j,l,e) * dzt(k,l)**2

                end do
             end do
          end do
       end do
    end do

    do j = 1, ly, ly-1
       do k = 1, lz, lz-1
          do e = 1, n
             d(1,j,k,e) = d(1,j,k,e) &
                  + G12(1,j,k,e) * dxt(1,1)*dyt(j,j) &
                  + G13(1,j,k,e) * dxt(1,1)*dzt(k,k)
             d(lx,j,k,e) = d(lx,j,k,e) &
                  + G12(lx,j,k,e) * dxt(lx,lx)*dyt(j,j) &
                  + G13(lx,j,k,e) * dxt(lx,lx)*dzt(k,k)
          end do
       end do
    end do

    do i = 1, lx, lx-1
       do k = 1, lz, lz-1
          do e = 1, n
             d(i,1,k,e) = d(i,1,k,e) &
                  + G12(i,1,k,e) * dyt(1,1)*dxt(i,i) &
                  + G23(i,1,k,e) * dyt(1,1)*dzt(k,k)
             d(i,ly,k,e) = d(i,ly,k,e) &
                  + G12(i,ly,k,e) * dyt(ly,ly)*dxt(i,i) &
                  + G23(i,ly,k,e) * dyt(ly,ly)*dzt(k,k)
          end do
       end do
    end do

    do i = 1, lx, lx-1
       do j = 1, ly, ly-1
          do e = 1, n
             d(i,j,1,e) = d(i,j,1,e) &
                  + G13(i,j,1,e) * dzt(1,1)*dxt(i,i) &
                  + G23(i,j,1,e) * dzt(1,1)*dyt(j,j)
             d(i,j,lz,e) = d(i,j,lz,e) &
                  + G13(i,j,lz,e) * dzt(lz,lz)*dxt(i,i) &
                  + G23(i,j,lz,e) * dzt(lz,lz)*dyt(j,j)
          end do
       end do
    end do

  end subroutine sx_update_lx12

  subroutine sx_update_lx11(d, dxt, dyt, dzt, G11, G22, G33, G12, G13, G23, n)
    integer, parameter :: lx = 11
    integer, parameter :: ly = 11
    integer, parameter :: lz = 11
    integer, intent(in) :: n
    real(kind=rp), intent(inout) :: d(lx, lx, lx, n)
    real(kind=rp), intent(in) :: dxt(lx,lx)
    real(kind=rp), intent(in) :: dyt(lx,lx)
    real(kind=rp), intent(in) :: dzt(lx,lx)
    real(kind=rp), intent(in) :: G11(lx, lx, lx, n)
    real(kind=rp), intent(in) :: G22(lx, lx, lx, n)
    real(kind=rp), intent(in) :: G33(lx, lx, lx, n)
    real(kind=rp), intent(in) :: G12(lx, lx, lx, n)
    real(kind=rp), intent(in) :: G13(lx, lx, lx, n)
    real(kind=rp), intent(in) :: G23(lx, lx, lx, n)
    integer :: i, j, k, l, e


    do l = 1, lx
       do k = 1, lz
          do j = 1, ly
             do i = 1, lx
                do e = 1, n
                   d(i,j,k,e) = d(i,j,k,e) + &
                        G11(l,j,k,e) * dxt(i,l)**2

                   d(i,j,k,e) = d(i,j,k,e) + &
                        G22(i,l,k,e) * dyt(j,l)**2

                   d(i,j,k,e) = d(i,j,k,e) + &
                        G33(i,j,l,e) * dzt(k,l)**2

                end do
             end do
          end do
       end do
    end do

    do j = 1, ly,ly-1
       do k = 1, lz,lz-1
          do e = 1, n
             d(1,j,k,e) = d(1,j,k,e) &
                  + G12(1,j,k,e) * dxt(1,1)*dyt(j,j) &
                  + G13(1,j,k,e) * dxt(1,1)*dzt(k,k)
             d(lx,j,k,e) = d(lx,j,k,e) &
                  + G12(lx,j,k,e) * dxt(lx,lx)*dyt(j,j) &
                  + G13(lx,j,k,e) * dxt(lx,lx)*dzt(k,k)
          end do
       end do
    end do

    do i = 1, lx, lx-1
       do k = 1, lz, lz-1
          do e = 1, n
             d(i,1,k,e) = d(i,1,k,e) &
                  + G12(i,1,k,e) * dyt(1,1)*dxt(i,i) &
                  + G23(i,1,k,e) * dyt(1,1)*dzt(k,k)
             d(i,ly,k,e) = d(i,ly,k,e) &
                  + G12(i,ly,k,e) * dyt(ly,ly)*dxt(i,i) &
                  + G23(i,ly,k,e) * dyt(ly,ly)*dzt(k,k)
          end do
       end do
    end do

    do i = 1, lx, lx-1
       do j = 1, ly, ly-1
          do e = 1, n
             d(i,j,1,e) = d(i,j,1,e) &
                  + G13(i,j,1,e) * dzt(1,1)*dxt(i,i) &
                  + G23(i,j,1,e) * dzt(1,1)*dyt(j,j)
             d(i,j,lz,e) = d(i,j,lz,e) &
                  + G13(i,j,lz,e) * dzt(lz,lz)*dxt(i,i) &
                  + G23(i,j,lz,e) * dzt(lz,lz)*dyt(j,j)
          end do
       end do
    end do

  end subroutine sx_update_lx11

  subroutine sx_update_lx10(d, dxt, dyt, dzt, G11, G22, G33, G12, G13, G23, n)
    integer, parameter :: lx = 10
    integer, parameter :: ly = 10
    integer, parameter :: lz = 10
    integer, intent(in) :: n
    real(kind=rp), intent(inout) :: d(lx, lx, lx, n)
    real(kind=rp), intent(in) :: dxt(lx,lx)
    real(kind=rp), intent(in) :: dyt(lx,lx)
    real(kind=rp), intent(in) :: dzt(lx,lx)
    real(kind=rp), intent(in) :: G11(lx, lx, lx, n)
    real(kind=rp), intent(in) :: G22(lx, lx, lx, n)
    real(kind=rp), intent(in) :: G33(lx, lx, lx, n)
    real(kind=rp), intent(in) :: G12(lx, lx, lx, n)
    real(kind=rp), intent(in) :: G13(lx, lx, lx, n)
    real(kind=rp), intent(in) :: G23(lx, lx, lx, n)
    integer :: i, j, k, l, e


    do l = 1, lx
       do k = 1, lz
          do j = 1, ly
             do i = 1, lx
                do e = 1, n
                   d(i,j,k,e) = d(i,j,k,e) + &
                        G11(l,j,k,e) * dxt(i,l)**2

                   d(i,j,k,e) = d(i,j,k,e) + &
                        G22(i,l,k,e) * dyt(j,l)**2

                   d(i,j,k,e) = d(i,j,k,e) + &
                        G33(i,j,l,e) * dzt(k,l)**2

                end do
             end do
          end do
       end do
    end do

    do j = 1, ly, ly-1
       do k = 1, lz, lz-1
          do e = 1, n
             d(1,j,k,e) = d(1,j,k,e) &
                  + G12(1,j,k,e) * dxt(1,1)*dyt(j,j) &
                  + G13(1,j,k,e) * dxt(1,1)*dzt(k,k)
             d(lx,j,k,e) = d(lx,j,k,e) &
                  + G12(lx,j,k,e) * dxt(lx,lx)*dyt(j,j) &
                  + G13(lx,j,k,e) * dxt(lx,lx)*dzt(k,k)
          end do
       end do
    end do

    do i = 1, lx, lx-1
       do k = 1, lz, lz-1
          do e = 1, n
             d(i,1,k,e) = d(i,1,k,e) &
                  + G12(i,1,k,e) * dyt(1,1)*dxt(i,i) &
                  + G23(i,1,k,e) * dyt(1,1)*dzt(k,k)
             d(i,ly,k,e) = d(i,ly,k,e) &
                  + G12(i,ly,k,e) * dyt(ly,ly)*dxt(i,i) &
                  + G23(i,ly,k,e) * dyt(ly,ly)*dzt(k,k)
          end do
       end do
    end do

    do i = 1, lx, lx-1
       do j = 1, ly, ly-1
          do e = 1, n
             d(i,j,1,e) = d(i,j,1,e) &
                  + G13(i,j,1,e) * dzt(1,1)*dxt(i,i) &
                  + G23(i,j,1,e) * dzt(1,1)*dyt(j,j)
             d(i,j,lz,e) = d(i,j,lz,e) &
                  + G13(i,j,lz,e) * dzt(lz,lz)*dxt(i,i) &
                  + G23(i,j,lz,e) * dzt(lz,lz)*dyt(j,j)
          end do
       end do
    end do

  end subroutine sx_update_lx10

  subroutine sx_update_lx9(d, dxt, dyt, dzt, G11, G22, G33, G12, G13, G23, n)
    integer, parameter :: lx = 9
    integer, parameter :: ly = 9
    integer, parameter :: lz = 9
    integer, intent(in) :: n
    real(kind=rp), intent(inout) :: d(lx, lx, lx, n)
    real(kind=rp), intent(in) :: dxt(lx,lx)
    real(kind=rp), intent(in) :: dyt(lx,lx)
    real(kind=rp), intent(in) :: dzt(lx,lx)
    real(kind=rp), intent(in) :: G11(lx, lx, lx, n)
    real(kind=rp), intent(in) :: G22(lx, lx, lx, n)
    real(kind=rp), intent(in) :: G33(lx, lx, lx, n)
    real(kind=rp), intent(in) :: G12(lx, lx, lx, n)
    real(kind=rp), intent(in) :: G13(lx, lx, lx, n)
    real(kind=rp), intent(in) :: G23(lx, lx, lx, n)
    integer :: i, j, k, l, e


    do l = 1, lx
       do k = 1, lz
          do j = 1, ly
             do i = 1, lx
                do e = 1, n
                   d(i,j,k,e) = d(i,j,k,e) + &
                        G11(l,j,k,e) * dxt(i,l)**2

                   d(i,j,k,e) = d(i,j,k,e) + &
                        G22(i,l,k,e) * dyt(j,l)**2

                   d(i,j,k,e) = d(i,j,k,e) + &
                        G33(i,j,l,e) * dzt(k,l)**2

                end do
             end do
          end do
       end do
    end do

    do j = 1, ly, ly-1
       do k = 1, lz, lz-1
          do e = 1, n
             d(1,j,k,e) = d(1,j,k,e) &
                  + G12(1,j,k,e) * dxt(1,1)*dyt(j,j) &
                  + G13(1,j,k,e) * dxt(1,1)*dzt(k,k)
             d(lx,j,k,e) = d(lx,j,k,e) &
                  + G12(lx,j,k,e) * dxt(lx,lx)*dyt(j,j) &
                  + G13(lx,j,k,e) * dxt(lx,lx)*dzt(k,k)
          end do
       end do
    end do

    do i = 1, lx, lx-1
       do k = 1, lz, lz-1
          do e = 1, n
             d(i,1,k,e) = d(i,1,k,e) &
                  + G12(i,1,k,e) * dyt(1,1)*dxt(i,i) &
                  + G23(i,1,k,e) * dyt(1,1)*dzt(k,k)
             d(i,ly,k,e) = d(i,ly,k,e) &
                  + G12(i,ly,k,e) * dyt(ly,ly)*dxt(i,i) &
                  + G23(i,ly,k,e) * dyt(ly,ly)*dzt(k,k)
          end do
       end do
    end do

    do i = 1, lx, lx-1
       do j = 1, ly, ly-1
          do e = 1, n
             d(i,j,1,e) = d(i,j,1,e) &
                  + G13(i,j,1,e) * dzt(1,1)*dxt(i,i) &
                  + G23(i,j,1,e) * dzt(1,1)*dyt(j,j)
             d(i,j,lz,e) = d(i,j,lz,e) &
                  + G13(i,j,lz,e) * dzt(lz,lz)*dxt(i,i) &
                  + G23(i,j,lz,e) * dzt(lz,lz)*dyt(j,j)
          end do
       end do
    end do

  end subroutine sx_update_lx9

  subroutine sx_update_lx8(d, dxt, dyt, dzt, G11, G22, G33, G12, G13, G23, n)
    integer, parameter :: lx = 8
    integer, parameter :: ly = 8
    integer, parameter :: lz = 8
    integer, intent(in) :: n
    real(kind=rp), intent(inout) :: d(lx, lx, lx, n)
    real(kind=rp), intent(in) :: dxt(lx,lx)
    real(kind=rp), intent(in) :: dyt(lx,lx)
    real(kind=rp), intent(in) :: dzt(lx,lx)
    real(kind=rp), intent(in) :: G11(lx, lx, lx, n)
    real(kind=rp), intent(in) :: G22(lx, lx, lx, n)
    real(kind=rp), intent(in) :: G33(lx, lx, lx, n)
    real(kind=rp), intent(in) :: G12(lx, lx, lx, n)
    real(kind=rp), intent(in) :: G13(lx, lx, lx, n)
    real(kind=rp), intent(in) :: G23(lx, lx, lx, n)
    integer :: i, j, k, l, e


    do l = 1, lx
       do k = 1, lz
          do j = 1, ly
             do i = 1, lx
                do e = 1, n
                   d(i,j,k,e) = d(i,j,k,e) + &
                        G11(l,j,k,e) * dxt(i,l)**2

                   d(i,j,k,e) = d(i,j,k,e) + &
                        G22(i,l,k,e) * dyt(j,l)**2

                   d(i,j,k,e) = d(i,j,k,e) + &
                        G33(i,j,l,e) * dzt(k,l)**2

                end do
             end do
          end do
       end do
    end do

    do j = 1, ly, ly-1
       do k = 1, lz, lz-1
          do e = 1, n
             d(1,j,k,e) = d(1,j,k,e) &
                  + G12(1,j,k,e) * dxt(1,1)*dyt(j,j) &
                  + G13(1,j,k,e) * dxt(1,1)*dzt(k,k)
             d(lx,j,k,e) = d(lx,j,k,e) &
                  + G12(lx,j,k,e) * dxt(lx,lx)*dyt(j,j) &
                  + G13(lx,j,k,e) * dxt(lx,lx)*dzt(k,k)
          end do
       end do
    end do

    do i = 1, lx, lx-1
       do k = 1, lz, lz-1
          do e = 1, n
             d(i,1,k,e) = d(i,1,k,e) &
                  + G12(i,1,k,e) * dyt(1,1)*dxt(i,i) &
                  + G23(i,1,k,e) * dyt(1,1)*dzt(k,k)
             d(i,ly,k,e) = d(i,ly,k,e) &
                  + G12(i,ly,k,e) * dyt(ly,ly)*dxt(i,i) &
                  + G23(i,ly,k,e) * dyt(ly,ly)*dzt(k,k)
          end do
       end do
    end do

    do i = 1, lx, lx-1
       do j = 1, ly, ly-1
          do e = 1, n
             d(i,j,1,e) = d(i,j,1,e) &
                  + G13(i,j,1,e) * dzt(1,1)*dxt(i,i) &
                  + G23(i,j,1,e) * dzt(1,1)*dyt(j,j)
             d(i,j,lz,e) = d(i,j,lz,e) &
                  + G13(i,j,lz,e) * dzt(lz,lz)*dxt(i,i) &
                  + G23(i,j,lz,e) * dzt(lz,lz)*dyt(j,j)
          end do
       end do
    end do

  end subroutine sx_update_lx8

  subroutine sx_update_lx7(d, dxt, dyt, dzt, G11, G22, G33, G12, G13, G23, n)
    integer, parameter :: lx = 7
    integer, parameter :: ly = 7
    integer, parameter :: lz = 7
    integer, intent(in) :: n
    real(kind=rp), intent(inout) :: d(lx, lx, lx, n)
    real(kind=rp), intent(in) :: G11(lx, lx, lx, n)
    real(kind=rp), intent(in) :: G22(lx, lx, lx, n)
    real(kind=rp), intent(in) :: G33(lx, lx, lx, n)
    real(kind=rp), intent(in) :: G12(lx, lx, lx, n)
    real(kind=rp), intent(in) :: G13(lx, lx, lx, n)
    real(kind=rp), intent(in) :: G23(lx, lx, lx, n)
    real(kind=rp), intent(in) :: dxt(lx,lx)
    real(kind=rp), intent(in) :: dyt(lx,lx)
    real(kind=rp), intent(in) :: dzt(lx,lx)
    integer :: i, j, k, l, e

    do l = 1, lx
       do k = 1, lz
          do j = 1, ly
             do i = 1, lx
                do e = 1, n
                   d(i,j,k,e) = d(i,j,k,e) + &
                        G11(l,j,k,e) * dxt(i,l)**2

                   d(i,j,k,e) = d(i,j,k,e) + &
                        G22(i,l,k,e) * dyt(j,l)**2

                   d(i,j,k,e) = d(i,j,k,e) + &
                        G33(i,j,l,e) * dzt(k,l)**2

                end do
             end do
          end do
       end do
    end do

    do j = 1, ly, ly-1
       do k = 1, lz, lz-1
          do e = 1, n
             d(1,j,k,e) = d(1,j,k,e) &
                  + G12(1,j,k,e) * dxt(1,1)*dyt(j,j) &
                  + G13(1,j,k,e) * dxt(1,1)*dzt(k,k)
             d(lx,j,k,e) = d(lx,j,k,e) &
                  + G12(lx,j,k,e) * dxt(lx,lx)*dyt(j,j) &
                  + G13(lx,j,k,e) * dxt(lx,lx)*dzt(k,k)
          end do
       end do
    end do

    do i = 1, lx, lx-1
       do k = 1, lz, lz-1
          do e = 1, n
             d(i,1,k,e) = d(i,1,k,e) &
                  + G12(i,1,k,e) * dyt(1,1)*dxt(i,i) &
                  + G23(i,1,k,e) * dyt(1,1)*dzt(k,k)
             d(i,ly,k,e) = d(i,ly,k,e) &
                  + G12(i,ly,k,e) * dyt(ly,ly)*dxt(i,i) &
                  + G23(i,ly,k,e) * dyt(ly,ly)*dzt(k,k)
          end do
       end do
    end do

    do i = 1, lx, lx-1
       do j = 1, ly, ly-1
          do e = 1, n
             d(i,j,1,e) = d(i,j,1,e) &
                  + G13(i,j,1,e) * dzt(1,1)*dxt(i,i) &
                  + G23(i,j,1,e) * dzt(1,1)*dyt(j,j)
             d(i,j,lz,e) = d(i,j,lz,e) &
                  + G13(i,j,lz,e) * dzt(lz,lz)*dxt(i,i) &
                  + G23(i,j,lz,e) * dzt(lz,lz)*dyt(j,j)
          end do
       end do
    end do


  end subroutine sx_update_lx7

  subroutine sx_update_lx6(d, dxt, dyt, dzt, G11, G22, G33, G12, G13, G23, n)
    integer, parameter :: lx = 6
    integer, parameter :: ly = 6
    integer, parameter :: lz = 6
    integer, intent(in) :: n
    real(kind=rp), intent(inout) :: d(lx, lx, lx, n)
    real(kind=rp), intent(in) :: G11(lx, lx, lx, n)
    real(kind=rp), intent(in) :: G22(lx, lx, lx, n)
    real(kind=rp), intent(in) :: G33(lx, lx, lx, n)
    real(kind=rp), intent(in) :: G12(lx, lx, lx, n)
    real(kind=rp), intent(in) :: G13(lx, lx, lx, n)
    real(kind=rp), intent(in) :: G23(lx, lx, lx, n)
    real(kind=rp), intent(in) :: dxt(lx,lx)
    real(kind=rp), intent(in) :: dyt(lx,lx)
    real(kind=rp), intent(in) :: dzt(lx,lx)
    integer :: i, j, k, l, e

    do l = 1, lx
       do k = 1, lz
          do j = 1, ly
             do i = 1, lx
                do e = 1, n
                   d(i,j,k,e) = d(i,j,k,e) + &
                        G11(l,j,k,e) * dxt(i,l)**2

                   d(i,j,k,e) = d(i,j,k,e) + &
                        G22(i,l,k,e) * dyt(j,l)**2

                   d(i,j,k,e) = d(i,j,k,e) + &
                        G33(i,j,l,e) * dzt(k,l)**2

                end do
             end do
          end do
       end do
    end do

    do j = 1, ly, ly-1
       do k = 1, lz, lz-1
          do e = 1, n
             d(1,j,k,e) = d(1,j,k,e) &
                  + G12(1,j,k,e) * dxt(1,1)*dyt(j,j) &
                  + G13(1,j,k,e) * dxt(1,1)*dzt(k,k)
             d(lx,j,k,e) = d(lx,j,k,e) &
                  + G12(lx,j,k,e) * dxt(lx,lx)*dyt(j,j) &
                  + G13(lx,j,k,e) * dxt(lx,lx)*dzt(k,k)
          end do
       end do
    end do

    do i = 1, lx, lx-1
       do k = 1, lz, lz-1
          do e = 1, n
             d(i,1,k,e) = d(i,1,k,e) &
                  + G12(i,1,k,e) * dyt(1,1)*dxt(i,i) &
                  + G23(i,1,k,e) * dyt(1,1)*dzt(k,k)
             d(i,ly,k,e) = d(i,ly,k,e) &
                  + G12(i,ly,k,e) * dyt(ly,ly)*dxt(i,i) &
                  + G23(i,ly,k,e) * dyt(ly,ly)*dzt(k,k)
          end do
       end do
    end do

    do i = 1, lx, lx-1
       do j = 1, ly, ly-1
          do e = 1, n
             d(i,j,1,e) = d(i,j,1,e) &
                  + G13(i,j,1,e) * dzt(1,1)*dxt(i,i) &
                  + G23(i,j,1,e) * dzt(1,1)*dyt(j,j)
             d(i,j,lz,e) = d(i,j,lz,e) &
                  + G13(i,j,lz,e) * dzt(lz,lz)*dxt(i,i) &
                  + G23(i,j,lz,e) * dzt(lz,lz)*dyt(j,j)
          end do
       end do
    end do


  end subroutine sx_update_lx6

  subroutine sx_update_lx5(d, dxt, dyt, dzt, G11, G22, G33, G12, G13, G23, n)
    integer, parameter :: lx = 5
    integer, parameter :: ly = 5
    integer, parameter :: lz = 5
    integer, intent(in) :: n
    real(kind=rp), intent(inout) :: d(lx, lx, lx, n)
    real(kind=rp), intent(in) :: G11(lx, lx, lx, n)
    real(kind=rp), intent(in) :: G22(lx, lx, lx, n)
    real(kind=rp), intent(in) :: G33(lx, lx, lx, n)
    real(kind=rp), intent(in) :: G12(lx, lx, lx, n)
    real(kind=rp), intent(in) :: G13(lx, lx, lx, n)
    real(kind=rp), intent(in) :: G23(lx, lx, lx, n)
    real(kind=rp), intent(in) :: dxt(lx,lx)
    real(kind=rp), intent(in) :: dyt(lx,lx)
    real(kind=rp), intent(in) :: dzt(lx,lx)
    integer :: i, j, k, l, e

    do l = 1, lx
       do k = 1, lz
          do j = 1, ly
             do i = 1, lx
                do e = 1, n
                   d(i,j,k,e) = d(i,j,k,e) + &
                        G11(l,j,k,e) * dxt(i,l)**2

                   d(i,j,k,e) = d(i,j,k,e) + &
                        G22(i,l,k,e) * dyt(j,l)**2

                   d(i,j,k,e) = d(i,j,k,e) + &
                        G33(i,j,l,e) * dzt(k,l)**2

                end do
             end do
          end do
       end do
    end do

    do j = 1, ly, ly-1
       do k = 1, lz, lz-1
          do e = 1, n
             d(1,j,k,e) = d(1,j,k,e) &
                  + G12(1,j,k,e) * dxt(1,1)*dyt(j,j) &
                  + G13(1,j,k,e) * dxt(1,1)*dzt(k,k)
             d(lx,j,k,e) = d(lx,j,k,e) &
                  + G12(lx,j,k,e) * dxt(lx,lx)*dyt(j,j) &
                  + G13(lx,j,k,e) * dxt(lx,lx)*dzt(k,k)
          end do
       end do
    end do

    do i = 1, lx, lx-1
       do k = 1, lz, lz-1
          do e = 1, n
             d(i,1,k,e) = d(i,1,k,e) &
                  + G12(i,1,k,e) * dyt(1,1)*dxt(i,i) &
                  + G23(i,1,k,e) * dyt(1,1)*dzt(k,k)
             d(i,ly,k,e) = d(i,ly,k,e) &
                  + G12(i,ly,k,e) * dyt(ly,ly)*dxt(i,i) &
                  + G23(i,ly,k,e) * dyt(ly,ly)*dzt(k,k)
          end do
       end do
    end do

    do i = 1, lx,lx-1
       do j = 1, ly,ly-1
          do e = 1, n
             d(i,j,1,e) = d(i,j,1,e) &
                  + G13(i,j,1,e) * dzt(1,1)*dxt(i,i) &
                  + G23(i,j,1,e) * dzt(1,1)*dyt(j,j)
             d(i,j,lz,e) = d(i,j,lz,e) &
                  + G13(i,j,lz,e) * dzt(lz,lz)*dxt(i,i) &
                  + G23(i,j,lz,e) * dzt(lz,lz)*dyt(j,j)
          end do
       end do
    end do


  end subroutine sx_update_lx5

  subroutine sx_update_lx4(d, dxt, dyt, dzt, G11, G22, G33, G12, G13, G23, n)
    integer, parameter :: lx = 4
    integer, parameter :: ly = 4
    integer, parameter :: lz = 4
    integer, intent(in) :: n
    real(kind=rp), intent(inout) :: d(lx, lx, lx, n)
    real(kind=rp), intent(in) :: G11(lx, lx, lx, n)
    real(kind=rp), intent(in) :: G22(lx, lx, lx, n)
    real(kind=rp), intent(in) :: G33(lx, lx, lx, n)
    real(kind=rp), intent(in) :: G12(lx, lx, lx, n)
    real(kind=rp), intent(in) :: G13(lx, lx, lx, n)
    real(kind=rp), intent(in) :: G23(lx, lx, lx, n)
    real(kind=rp), intent(in) :: dxt(lx,lx)
    real(kind=rp), intent(in) :: dyt(lx,lx)
    real(kind=rp), intent(in) :: dzt(lx,lx)
    integer :: i, j, k, l, e

    do l = 1, lx
       do k = 1, lz
          do j = 1, ly
             do i = 1, lx
                do e = 1, n
                   d(i,j,k,e) = d(i,j,k,e) + &
                        G11(l,j,k,e) * dxt(i,l)**2

                   d(i,j,k,e) = d(i,j,k,e) + &
                        G22(i,l,k,e) * dyt(j,l)**2

                   d(i,j,k,e) = d(i,j,k,e) + &
                        G33(i,j,l,e) * dzt(k,l)**2

                end do
             end do
          end do
       end do
    end do

    do j = 1, ly, ly-1
       do k = 1, lz, lz-1
          do e = 1, n
             d(1,j,k,e) = d(1,j,k,e) &
                  + G12(1,j,k,e) * dxt(1,1)*dyt(j,j) &
                  + G13(1,j,k,e) * dxt(1,1)*dzt(k,k)
             d(lx,j,k,e) = d(lx,j,k,e) &
                  + G12(lx,j,k,e) * dxt(lx,lx)*dyt(j,j) &
                  + G13(lx,j,k,e) * dxt(lx,lx)*dzt(k,k)
          end do
       end do
    end do

    do i = 1, lx, lx-1
       do k = 1, lz, lz-1
          do e = 1, n
             d(i,1,k,e) = d(i,1,k,e) &
                  + G12(i,1,k,e) * dyt(1,1)*dxt(i,i) &
                  + G23(i,1,k,e) * dyt(1,1)*dzt(k,k)
             d(i,ly,k,e) = d(i,ly,k,e) &
                  + G12(i,ly,k,e) * dyt(ly,ly)*dxt(i,i) &
                  + G23(i,ly,k,e) * dyt(ly,ly)*dzt(k,k)
          end do
       end do
    end do

    do i = 1, lx, lx-1
       do j = 1, ly, ly-1
          do e = 1, n
             d(i,j,1,e) = d(i,j,1,e) &
                  + G13(i,j,1,e) * dzt(1,1)*dxt(i,i) &
                  + G23(i,j,1,e) * dzt(1,1)*dyt(j,j)
             d(i,j,lz,e) = d(i,j,lz,e) &
                  + G13(i,j,lz,e) * dzt(lz,lz)*dxt(i,i) &
                  + G23(i,j,lz,e) * dzt(lz,lz)*dyt(j,j)
          end do
       end do
    end do


  end subroutine sx_update_lx4

  subroutine sx_update_lx3(d, dxt, dyt, dzt, G11, G22, G33, G12, G13, G23, n)
    integer, parameter :: lx = 3
    integer, parameter :: ly = 3
    integer, parameter :: lz = 3
    integer, intent(in) :: n
    real(kind=rp), intent(inout) :: d(lx, lx, lx, n)
    real(kind=rp), intent(in) :: G11(lx, lx, lx, n)
    real(kind=rp), intent(in) :: G22(lx, lx, lx, n)
    real(kind=rp), intent(in) :: G33(lx, lx, lx, n)
    real(kind=rp), intent(in) :: G12(lx, lx, lx, n)
    real(kind=rp), intent(in) :: G13(lx, lx, lx, n)
    real(kind=rp), intent(in) :: G23(lx, lx, lx, n)
    real(kind=rp), intent(in) :: dxt(lx,lx)
    real(kind=rp), intent(in) :: dyt(lx,lx)
    real(kind=rp), intent(in) :: dzt(lx,lx)
    integer :: i, j, k, l, e


    do l = 1, lx
       do k = 1, lz
          do j = 1, ly
             do i = 1, lx
                do e = 1, n
                   d(i,j,k,e) = d(i,j,k,e) + &
                        G11(l,j,k,e) * dxt(i,l)**2

                   d(i,j,k,e) = d(i,j,k,e) + &
                        G22(i,l,k,e) * dyt(j,l)**2

                   d(i,j,k,e) = d(i,j,k,e) + &
                        G33(i,j,l,e) * dzt(k,l)**2

                end do
             end do
          end do
       end do
    end do

    do j = 1, ly, ly-1
       do k = 1, lz, lz-1
          do e = 1, n
             d(1,j,k,e) = d(1,j,k,e) &
                  + G12(1,j,k,e) * dxt(1,1)*dyt(j,j) &
                  + G13(1,j,k,e) * dxt(1,1)*dzt(k,k)
             d(lx,j,k,e) = d(lx,j,k,e) &
                  + G12(lx,j,k,e) * dxt(lx,lx)*dyt(j,j) &
                  + G13(lx,j,k,e) * dxt(lx,lx)*dzt(k,k)
          end do
       end do
    end do

    do i = 1, lx, lx-1
       do k = 1, lz, lz-1
          do e = 1, n
             d(i,1,k,e) = d(i,1,k,e) &
                  + G12(i,1,k,e) * dyt(1,1)*dxt(i,i) &
                  + G23(i,1,k,e) * dyt(1,1)*dzt(k,k)
             d(i,ly,k,e) = d(i,ly,k,e) &
                  + G12(i,ly,k,e) * dyt(ly,ly)*dxt(i,i) &
                  + G23(i,ly,k,e) * dyt(ly,ly)*dzt(k,k)
          end do
       end do
    end do

    do i = 1, lx, lx-1
       do j = 1, ly, ly-1
          do e = 1, n
             d(i,j,1,e) = d(i,j,1,e) &
                  + G13(i,j,1,e) * dzt(1,1)*dxt(i,i) &
                  + G23(i,j,1,e) * dzt(1,1)*dyt(j,j)
             d(i,j,lz,e) = d(i,j,lz,e) &
                  + G13(i,j,lz,e) * dzt(lz,lz)*dxt(i,i) &
                  + G23(i,j,lz,e) * dzt(lz,lz)*dyt(j,j)
          end do
       end do
    end do

  end subroutine sx_update_lx3

  subroutine sx_update_lx2(d, dxt, dyt, dzt, G11, G22, G33, G12, G13, G23, n)
    integer, parameter :: lx = 2
    integer, parameter :: ly = 2
    integer, parameter :: lz = 2
    integer, intent(in) :: n
    real(kind=rp), intent(inout) :: d(lx, lx, lx, n)
    real(kind=rp), intent(in) :: G11(lx, lx, lx, n)
    real(kind=rp), intent(in) :: G22(lx, lx, lx, n)
    real(kind=rp), intent(in) :: G33(lx, lx, lx, n)
    real(kind=rp), intent(in) :: G12(lx, lx, lx, n)
    real(kind=rp), intent(in) :: G13(lx, lx, lx, n)
    real(kind=rp), intent(in) :: G23(lx, lx, lx, n)
    real(kind=rp), intent(in) :: dxt(lx,lx)
    real(kind=rp), intent(in) :: dyt(lx,lx)
    real(kind=rp), intent(in) :: dzt(lx,lx)
    integer :: i, j, k, l, e


    do l = 1, lx
       do k = 1, lz
          do j = 1, ly
             do i = 1, lx
                do e = 1, n
                   d(i,j,k,e) = d(i,j,k,e) + &
                        G11(l,j,k,e) * dxt(i,l)**2

                   d(i,j,k,e) = d(i,j,k,e) + &
                        G22(i,l,k,e) * dyt(j,l)**2

                   d(i,j,k,e) = d(i,j,k,e) + &
                        G33(i,j,l,e) * dzt(k,l)**2

                end do
             end do
          end do
       end do
    end do

    do j = 1, ly, ly-1
       do k = 1, lz, lz-1
          do e = 1, n
             d(1,j,k,e) = d(1,j,k,e) &
                  + G12(1,j,k,e) * dxt(1,1)*dyt(j,j) &
                  + G13(1,j,k,e) * dxt(1,1)*dzt(k,k)
             d(lx,j,k,e) = d(lx,j,k,e) &
                  + G12(lx,j,k,e) * dxt(lx,lx)*dyt(j,j) &
                  + G13(lx,j,k,e) * dxt(lx,lx)*dzt(k,k)
          end do
       end do
    end do

    do i = 1, lx, lx-1
       do k = 1, lz, lz-1
          do e = 1, n
             d(i,1,k,e) = d(i,1,k,e) &
                  + G12(i,1,k,e) * dyt(1,1)*dxt(i,i) &
                  + G23(i,1,k,e) * dyt(1,1)*dzt(k,k)
             d(i,ly,k,e) = d(i,ly,k,e) &
                  + G12(i,ly,k,e) * dyt(ly,ly)*dxt(i,i) &
                  + G23(i,ly,k,e) * dyt(ly,ly)*dzt(k,k)
          end do
       end do
    end do

    do i = 1, lx, lx-1
       do j = 1, ly, ly-1
          do e = 1, n
             d(i,j,1,e) = d(i,j,1,e) &
                  + G13(i,j,1,e) * dzt(1,1)*dxt(i,i) &
                  + G23(i,j,1,e) * dzt(1,1)*dyt(j,j)
             d(i,j,lz,e) = d(i,j,lz,e) &
                  + G13(i,j,lz,e) * dzt(lz,lz)*dxt(i,i) &
                  + G23(i,j,lz,e) * dzt(lz,lz)*dyt(j,j)
          end do
       end do
    end do

  end subroutine sx_update_lx2

end module sx_jacobi

