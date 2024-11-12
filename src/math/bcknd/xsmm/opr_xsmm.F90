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
!> Operators libxsmm backend
module opr_xsmm
  use num_types, only : rp
  use mxm_wrapper, only : mxm
  use space, only : space_t
  use coefs, only : coef_t
  use math, only : rzero, col2, col3, sub3, add2, addcol3, invcol2, copy
  use mesh, only : mesh_t
  use field, only : field_t
  use interpolation, only : interpolator_t
  use gather_scatter, only : gs_t, GS_OP_ADD
  use mathops, only : opcolv
#ifdef HAVE_LIBXSMM
  use libxsmm, only: libxsmm_mmcall => libxsmm_dmmcall_abc, &
                     libxsmm_dmmfunction, libxsmm_dispatch, &
                     LIBXSMM_PREFETCH_AUTO
#endif
  implicit none
  private

  public :: opr_xsmm_dudxyz, opr_xsmm_opgrad, opr_xsmm_cdtp, opr_xsmm_conv1, &
            opr_xsmm_curl, opr_xsmm_convect_scalar, opr_xsmm_set_convect_rst

#ifdef HAVE_LIBXSMM
  type(libxsmm_dmmfunction), private :: lgrad_xmm1
  type(libxsmm_dmmfunction), private :: lgrad_xmm2
  type(libxsmm_dmmfunction), private :: lgrad_xmm3
#endif

contains

  subroutine opr_xsmm_dudxyz(du, u, dr, ds, dt, coef)
    type(coef_t), intent(in), target :: coef
    real(kind=rp), dimension(coef%Xh%lx, coef%Xh%ly, coef%Xh%lz, &
                             coef%msh%nelv), intent(inout) ::  du
    real(kind=rp), dimension(coef%Xh%lx, coef%Xh%ly, coef%Xh%lz, &
                             coef%msh%nelv), intent(in) ::  u, dr, ds, dt
#ifdef HAVE_LIBXSMM
    real(kind=rp) :: drst(coef%Xh%lx, coef%Xh%ly, coef%Xh%lz)
    type(space_t), pointer :: Xh
    type(mesh_t), pointer :: msh
    integer :: e, k, lxy, lyz, lxyz
    type(libxsmm_dmmfunction), save :: dudxyz_xmm1
    type(libxsmm_dmmfunction), save :: dudxyz_xmm2
    type(libxsmm_dmmfunction), save :: dudxyz_xmm3
    logical, save :: dudxyz_xsmm_init = .false.

    Xh => coef%Xh
    msh => coef%msh
    lxy  = Xh%lx*Xh%ly
    lyz  = Xh%ly*Xh%lz
    lxyz = Xh%lx*Xh%ly*Xh%lz

    if (.not. dudxyz_xsmm_init) then
       call libxsmm_dispatch(dudxyz_xmm1, Xh%lx, Xh%ly*Xh%lz, Xh%lx, &
            alpha = 1d0, beta = 0d0, prefetch = LIBXSMM_PREFETCH_AUTO)
       call libxsmm_dispatch(dudxyz_xmm2, Xh%lx, Xh%ly, Xh%ly, &
            alpha = 1d0, beta = 0d0, prefetch = LIBXSMM_PREFETCH_AUTO)
       call libxsmm_dispatch(dudxyz_xmm3, Xh%lx*Xh%ly, Xh%lz, Xh%lz, &
            alpha = 1d0, beta = 0d0, prefetch = LIBXSMM_PREFETCH_AUTO)
       dudxyz_xsmm_init = .true.
    end if

    do e = 1, msh%nelv
       if (msh%gdim .eq. 2) then
          call mxm(Xh%dx, Xh%lx, u(1,1,1,e), Xh%lx, du(1,1,1,e), lyz)
          call col2(du(1,1,1,e), dr(1,1,1,e), lxyz)
          call mxm(U(1,1,1,e), Xh%lx, Xh%dyt, Xh%ly, drst, Xh%ly)
          call addcol3(du(1,1,1,e), drst, ds(1,1,1,e), lxyz)
       else
          call libxsmm_mmcall(dudxyz_xmm1, Xh%dx, u(1,1,1,e), du(1,1,1,e))
          call col2(du(1,1,1,e), dr(1,1,1,e), lxyz)
          do k = 1, Xh%lz
             call libxsmm_mmcall(dudxyz_xmm2, u(1,1,k,e), Xh%dyt, drst(1,1,k))
          end do
          call addcol3(du(1,1,1,e), drst, ds(1,1,1,e), lxyz)
          call libxsmm_mmcall(dudxyz_xmm3, u(1,1,1,e), Xh%dzt, drst)
          call addcol3(du(1,1,1,e), drst, dt(1,1,1,e), lxyz)
       end if
    end do
    call col2(du, coef%jacinv, coef%dof%n_dofs)

#endif
  end subroutine opr_xsmm_dudxyz

  subroutine opr_xsmm_opgrad(ux, uy, uz, u, coef)
    type(coef_t), intent(in) :: coef
    real(kind=rp), dimension(coef%Xh%lxyz, coef%msh%nelv), intent(inout) :: ux
    real(kind=rp), dimension(coef%Xh%lxyz, coef%msh%nelv), intent(inout) :: uy
    real(kind=rp), dimension(coef%Xh%lxyz, coef%msh%nelv), intent(inout) :: uz
    real(kind=rp), dimension(coef%Xh%lxyz, coef%msh%nelv), intent(in) :: u
#ifdef HAVE_LIBXSMM
    real(kind=rp) :: ur(coef%Xh%lxyz)
    real(kind=rp) :: us(coef%Xh%lxyz)
    real(kind=rp) :: ut(coef%Xh%lxyz)
    logical, save :: lgrad_xsmm_init = .false.
    integer, save :: init_size = 0
    integer :: e, i, N
    N = coef%Xh%lx - 1


    if ((.not. lgrad_xsmm_init) .or. &
         (init_size .gt. 0 .and. init_size .ne. N)) then
       call libxsmm_dispatch(lgrad_xmm1, (N+1), (N+1)**2, (N+1), &
            alpha = 1d0, beta = 0d0, prefetch = LIBXSMM_PREFETCH_AUTO)
       call libxsmm_dispatch(lgrad_xmm2, (N+1), (N+1), (N+1), &
            alpha = 1d0, beta = 0d0, prefetch = LIBXSMM_PREFETCH_AUTO)
       call libxsmm_dispatch(lgrad_xmm3, (N+1)**2, (N+1), (N+1), &
            alpha = 1d0, beta = 0d0, prefetch = LIBXSMM_PREFETCH_AUTO)
       lgrad_xsmm_init = .true.
       init_size = N
    end if

    do e = 1, coef%msh%nelv
       if (coef%msh%gdim .eq. 3) then
          call local_grad3_xsmm(ur, us, ut, u(1,e), N, coef%Xh%dx, coef%Xh%dxt)
          do i = 1, coef%Xh%lxyz
             ux(i,e) = coef%Xh%w3(i,1,1) * (ur(i) * coef%drdx(i,1,1,e) &
                  + us(i) * coef%dsdx(i,1,1,e) &
                  + ut(i) * coef%dtdx(i,1,1,e) )
             uy(i,e) = coef%Xh%w3(i,1,1) * (ur(i) * coef%drdy(i,1,1,e) &
                  + us(i) * coef%dsdy(i,1,1,e) &
                  + ut(i) * coef%dtdy(i,1,1,e) )
             uz(i,e) = coef%Xh%w3(i,1,1) * (ur(i) * coef%drdz(i,1,1,e) &
                  + us(i) * coef%dsdz(i,1,1,e) &
                  + ut(i) * coef%dtdz(i,1,1,e) )
          end do
       else

          call local_grad2(ur, us, u(1,e), N, coef%Xh%dx, coef%Xh%dyt)

          do i = 1, coef%Xh%lxyz
             ux(i,e) = coef%Xh%w3(i,1,1) * (ur(i) * coef%drdx(i,1,1,e) &
                  + us(i) * coef%dsdx(i,1,1,e) )
             uy(i,e) = coef%Xh%w3(i,1,1) * (ur(i) * coef%drdy(i,1,1,e) &
                  + us(i) * coef%dsdy(i,1,1,e) )
          end do
       end if
    end do
#endif
  end subroutine opr_xsmm_opgrad

  subroutine local_grad3_xsmm(ur, us, ut, u, n, D, Dt)
    integer, intent(in) :: n
    real(kind=rp), intent(inout) :: ur(0:n, 0:n, 0:n)
    real(kind=rp), intent(inout) :: us(0:n, 0:n, 0:n)
    real(kind=rp), intent(inout) :: ut(0:n, 0:n, 0:n)
    real(kind=rp), intent(in) :: u(0:n, 0:n, 0:n)
    real(kind=rp), intent(in) :: D(0:n, 0:n)
    real(kind=rp), intent(in) :: Dt(0:n, 0:n)
#ifdef HAVE_LIBXSMM
    integer :: m1, m2, k

    m1 = n + 1
    m2 = m1*m1

    call libxsmm_mmcall(lgrad_xmm1, D, u, ur)
    do k = 0, n
       call libxsmm_mmcall(lgrad_xmm2, u(0,0,k), Dt, us(0,0,k))
    end do
    call libxsmm_mmcall(lgrad_xmm3, u, Dt, ut)
#endif

  end subroutine local_grad3_xsmm

  subroutine local_grad2(ur, us, u, n, D, Dt)
    integer, intent(in) :: n
    real(kind=rp), intent(inout) :: ur(0:n, 0:n)
    real(kind=rp), intent(inout) :: us(0:n, 0:n)
    real(kind=rp), intent(in) :: u(0:n, 0:n)
    real(kind=rp), intent(in) :: D(0:n, 0:n)
    real(kind=rp), intent(in) :: Dt(0:n, 0:n)
    integer :: m1

    m1 = n + 1

    call mxm(D, m1, u, m1, ur, m1)
    call mxm(u, m1, Dt, m1, us, m1)

  end subroutine local_grad2

  subroutine opr_xsmm_cdtp(dtx, x, dr, ds, dt, coef)
    type(coef_t), intent(in) :: coef
    real(kind=rp), dimension(coef%Xh%lxyz, coef%msh%nelv), intent(inout) :: dtx
    real(kind=rp), dimension(coef%Xh%lxyz, coef%msh%nelv), intent(inout) :: x
    real(kind=rp), dimension(coef%Xh%lxyz, coef%msh%nelv), intent(in) :: dr
    real(kind=rp), dimension(coef%Xh%lxyz, coef%msh%nelv), intent(in) :: ds
    real(kind=rp), dimension(coef%Xh%lxyz, coef%msh%nelv), intent(in) :: dt
#ifdef HAVE_LIBXSMM
    real(kind=rp) :: wx(coef%Xh%lxyz)
    real(kind=rp) :: ta1(coef%Xh%lxyz)
    real(kind=rp) :: ta2(coef%Xh%lxyz)
    real(kind=rp) :: ta3(coef%Xh%lxyz)
    integer :: e, i1, i2, n1, n2, iz
    type(space_t), pointer :: Xh

    type(libxsmm_dmmfunction), save :: cdtp_xmm1
    type(libxsmm_dmmfunction), save :: cdtp_xmm2
    type(libxsmm_dmmfunction), save :: cdtp_xmm3
    logical, save :: cdtp_xsmm_init = .false.

    Xh => coef%Xh
    n1 = Xh%lx*Xh%ly
    n2 = Xh%lx*Xh%ly

    if (.not. cdtp_xsmm_init) then
       call libxsmm_dispatch(cdtp_xmm1, Xh%lx, Xh%ly*Xh%lz, Xh%lx, &
            alpha = 1d0, beta = 0d0, prefetch = LIBXSMM_PREFETCH_AUTO)
       call libxsmm_dispatch(cdtp_xmm2, Xh%lx, Xh%ly, Xh%ly, &
            alpha = 1d0, beta = 0d0, prefetch = LIBXSMM_PREFETCH_AUTO)
       call libxsmm_dispatch(cdtp_xmm3, Xh%lx*Xh%ly, Xh%lz, Xh%lz, &
            alpha = 1d0, beta = 0d0, prefetch = LIBXSMM_PREFETCH_AUTO)
       cdtp_xsmm_init = .true.
    end if

    do e = 1, coef%msh%nelv
       call col3(wx, coef%B(1,1,1,e), x(1,e), Xh%lxyz)
       call invcol2(wx, coef%jac(1,1,1,e), Xh%lxyz)
       call col3(ta1, wx, dr(1,e), Xh%lxyz)
       call libxsmm_mmcall(cdtp_xmm1, Xh%dxt, ta1, dtx(1,e))
       call col3 (ta1, wx, ds(1,e), Xh%lxyz)
       i1 = 1
       i2 = 1
       do iz = 1, Xh%lz
          call libxsmm_mmcall(cdtp_xmm2, ta1(i2), Xh%dy, ta2(i1))
          i1 = i1 + n1
          i2 = i2 + n2
       end do
       call add2(dtx(1,e), ta2, Xh%lxyz)
       call col3(ta1, wx, dt(1,e), Xh%lxyz)
       call libxsmm_mmcall(cdtp_xmm3, ta1, Xh%dz, ta2)
       call add2 (dtx(1,e), ta2, Xh%lxyz)
    end do
#endif
  end subroutine opr_xsmm_cdtp

  subroutine opr_xsmm_conv1(du,u, vx, vy, vz, Xh, coef, nelv, gdim)
    type(space_t), intent(in) :: Xh
    type(coef_t), intent(in) :: coef
    integer, intent(in) :: nelv, gdim
    real(kind=rp), intent(inout) ::  du(Xh%lxyz, nelv)
    real(kind=rp), intent(inout), dimension(Xh%lx, Xh%ly, Xh%lz, nelv) ::  u
    real(kind=rp), intent(inout), dimension(Xh%lx, Xh%ly, Xh%lz, nelv) ::  vx
    real(kind=rp), intent(inout), dimension(Xh%lx, Xh%ly, Xh%lz, nelv) ::  vy
    real(kind=rp), intent(inout), dimension(Xh%lx, Xh%ly, Xh%lz, nelv) ::  vz
#ifdef HAVE_LIBXSMM
    !   Store the inverse jacobian to speed this operation up
    real(kind=rp), dimension(Xh%lx, Xh%ly, Xh%lz) :: dudr, duds, dudt
    integer :: ie, iz, i

    type(libxsmm_dmmfunction), save :: conv1_xmm1
    type(libxsmm_dmmfunction), save :: conv1_xmm2
    type(libxsmm_dmmfunction), save :: conv1_xmm3
    logical, save :: conv1_xsmm_init = .false.

    if (.not. conv1_xsmm_init) then
       call libxsmm_dispatch(conv1_xmm1, Xh%lx, Xh%ly*Xh%lx, Xh%lx, &
            alpha = 1d0, beta = 0d0, prefetch = LIBXSMM_PREFETCH_AUTO)
       call libxsmm_dispatch(conv1_xmm2, Xh%lx, Xh%ly, Xh%ly, &
            alpha = 1d0, beta = 0d0, prefetch = LIBXSMM_PREFETCH_AUTO)
       call libxsmm_dispatch(conv1_xmm3, Xh%lx*Xh%ly, Xh%lz, Xh%lz, &
            alpha = 1d0, beta = 0d0, prefetch = LIBXSMM_PREFETCH_AUTO)
       conv1_xsmm_init = .true.
    end if

    !   Compute vel.grad(u)
    do ie = 1, nelv
       if (gdim .eq. 3) then
          call libxsmm_mmcall(conv1_xmm1, Xh%dx, u(1,1,1, ie), dudr)
          do iz = 1, Xh%lz
             call libxsmm_mmcall(conv1_xmm2, u(1,1, iz, ie), Xh%dyt,&
                                 duds(1,1, iz))
          end do
          call libxsmm_mmcall(conv1_xmm3, u(1,1,1, ie), Xh%dzt, dudt)
          do i = 1, Xh%lxyz
             du(i, ie) = coef%jacinv(i,1,1, ie) * ( &
                  vx(i,1,1, ie) * ( &
                  coef%drdx(i,1,1, ie) * dudr(i,1,1) &
                  + coef%dsdx(i,1,1, ie) * duds(i,1,1) &
                  + coef%dtdx(i,1,1, ie) * dudt(i,1,1)) &
                  + vy(i,1,1, ie) * ( &
                  coef%drdy(i,1,1, ie) * dudr(i,1,1) &
                  + coef%dsdy(i,1,1, ie) * duds(i,1,1) &
                  + coef%dtdy(i,1,1, ie) * dudt(i,1,1)) &
                  + vz(i,1,1, ie) * ( &
                  coef%drdz(i,1,1, ie) * dudr(i,1,1) &
                  + coef%dsdz(i,1,1, ie) * duds(i,1,1) &
                  + coef%dtdz(i,1,1, ie) * dudt(i,1,1)))
          end do
       else
          !        2D
          call mxm(Xh%dx, Xh%lx, u(1,1,1, ie), Xh%lx, dudr, Xh%lyz)
          call mxm(u(1,1,1, ie), Xh%lx, Xh%dyt, Xh%ly, duds, Xh%ly)
          do i = 1, Xh%lxyz
             du(i, ie) = coef%jacinv(i,1,1, ie) * ( &
                  vx(i,1,1, ie) * ( &
                  coef%drdx(i,1,1, ie) * dudr(i,1,1) &
                  + coef%dsdx(i,1,1, ie) * duds(i,1,1)) &
                  + vy(i,1,1, ie) * ( &
                  coef%drdy(i,1,1, ie) * dudr(i,1,1) &
                  + coef%dsdy(i,1,1, ie) * duds(i,1,1)))
          end do
       end if
    end do

#endif

  end subroutine opr_xsmm_conv1

  subroutine opr_xsmm_convect_scalar(du, u, c, Xh_GLL, Xh_GL, coef_GLL, &
                                     coef_GL, GLL_to_GL)
    type(space_t), intent(in) :: Xh_GL
    type(space_t), intent(in) :: Xh_GLL
    type(coef_t), intent(in) :: coef_GLL
    type(coef_t), intent(in) :: coef_GL
    type(interpolator_t), intent(inout) :: GLL_to_GL
    real(kind=rp), intent(inout) :: du(Xh_GLL%lx, Xh_GLL%ly, Xh_GLL%lz, &
                                       coef_GL%msh%nelv)
    real(kind=rp), intent(inout) :: u(Xh_GL%lxyz, coef_GL%msh%nelv)
    real(kind=rp), intent(inout) :: c(Xh_GL%lxyz, coef_GL%msh%nelv, 3)
    real(kind=rp) :: ur(Xh_GL%lxyz)
    real(kind=rp) :: us(Xh_GL%lxyz)
    real(kind=rp) :: ut(Xh_GL%lxyz)
    real(kind=rp) :: ud(Xh_GL%lxyz, coef_GL%msh%nelv)
    logical, save :: lgrad_xsmm_init = .false.
    integer, save :: init_size = 0
    integer :: e, i, N, n_GLL
    N = coef_GL%Xh%lx - 1
    n_GLL = coef_GLL%msh%nelv * Xh_GLL%lxyz

#ifdef HAVE_LIBXSMM
    if ((.not. lgrad_xsmm_init) .or. &
         (init_size .gt. 0 .and. init_size .ne. N)) then
       call libxsmm_dispatch(lgrad_xmm1, (N+1), (N+1)**2, (N+1), &
            alpha = 1d0, beta = 0d0, prefetch = LIBXSMM_PREFETCH_AUTO)
       call libxsmm_dispatch(lgrad_xmm2, (N+1), (N+1), (N+1), &
            alpha = 1d0, beta = 0d0, prefetch = LIBXSMM_PREFETCH_AUTO)
       call libxsmm_dispatch(lgrad_xmm3, (N+1)**2, (N+1), (N+1), &
            alpha = 1d0, beta = 0d0, prefetch = LIBXSMM_PREFETCH_AUTO)
       lgrad_xsmm_init = .true.
       init_size = N
    end if

    do e = 1, coef_GLL%msh%nelv
       call local_grad3_xsmm(ur, us, ut, u(1,e), N, Xh_GL%dx, Xh_GL%dxt)
       do i = 1, Xh_GL%lxyz
          ud(i,e) = c(i,e,1) * ur(i) + c(i,e,2) * us(i) + c(i,e,3) * ut(i)
       end do
    end do
#endif
    call GLL_to_GL%map(du, ud, coef_GL%msh%nelv, Xh_GLL)
    call coef_GLL%gs_h%op(du, n_GLL, GS_OP_ADD)
    call col2(du, coef_GLL%Binv, n_GLL)
  end subroutine opr_xsmm_convect_scalar

  subroutine opr_xsmm_curl(w1, w2, w3, u1, u2, u3, work1, work2, c_Xh)
    type(field_t), intent(inout) :: w1
    type(field_t), intent(inout) :: w2
    type(field_t), intent(inout) :: w3
    type(field_t), intent(inout) :: u1
    type(field_t), intent(inout) :: u2
    type(field_t), intent(inout) :: u3
    type(field_t), intent(inout) :: work1
    type(field_t), intent(inout) :: work2
    type(coef_t), intent(in)  :: c_Xh
    integer :: gdim, n

    n = w1%dof%size()
    gdim = c_Xh%msh%gdim

    !     this%work1=dw/dy ; this%work2=dv/dz
    call opr_xsmm_dudxyz(work1%x, u3%x, c_Xh%drdy, c_Xh%dsdy, c_Xh%dtdy, c_Xh)
    if (gdim .eq. 3) then
       call opr_xsmm_dudxyz(work2%x, u2%x, c_Xh%drdz, c_Xh%dsdz, &
                            c_Xh%dtdz, c_Xh)
       call sub3(w1%x, work1%x, work2%x, n)
    else
       call copy(w1%x, work1%x, n)
    end if
    !     this%work1=du/dz ; this%work2=dw/dx
    if (gdim .eq. 3) then
       call opr_xsmm_dudxyz(work1%x, u1%x, c_Xh%drdz, c_Xh%dsdz, &
                            c_Xh%dtdz, c_Xh)
       call opr_xsmm_dudxyz(work2%x, u3%x, c_Xh%drdx, c_Xh%dsdx, &
                            c_Xh%dtdx, c_Xh)
       call sub3(w2%x, work1%x, work2%x, n)
    else
       call rzero (work1%x, n)
       call opr_xsmm_dudxyz(work2%x, u3%x, c_Xh%drdx, c_Xh%dsdx, &
                            c_Xh%dtdx, c_Xh)
       call sub3(w2%x, work1%x, work2%x, n)
    end if
    !     this%work1=dv/dx ; this%work2=du/dy
    call opr_xsmm_dudxyz(work1%x, u2%x, c_Xh%drdx, c_Xh%dsdx, c_Xh%dtdx, c_Xh)
    call opr_xsmm_dudxyz(work2%x, u1%x, c_Xh%drdy, c_Xh%dsdy, c_Xh%dtdy, c_Xh)
    call sub3(w3%x, work1%x, work2%x, n)
    !!    BC dependent, Needs to change if cyclic

    call opcolv(w1%x, w2%x, w3%x, c_Xh%B, gdim, n)
    call c_Xh%gs_h%op(w1, GS_OP_ADD)
    call c_Xh%gs_h%op(w2, GS_OP_ADD)
    call c_Xh%gs_h%op(w3, GS_OP_ADD)
    call opcolv(w1%x, w2%x, w3%x, c_Xh%Binv, gdim, n)

  end subroutine opr_xsmm_curl

  subroutine opr_xsmm_set_convect_rst(cr, cs, ct, cx, cy, cz, Xh, coef)
    type(space_t), intent(inout) :: Xh
    type(coef_t), intent(inout) :: coef
    real(kind=rp), dimension(Xh%lxyz, coef%msh%nelv), &
                   intent(inout) :: cr, cs, ct
    real(kind=rp), dimension(Xh%lxyz, coef%msh%nelv), &
                   intent(in) :: cx, cy, cz
    integer :: e, i, t, nxyz

    associate(drdx => coef%drdx, drdy => coef%drdy, drdz => coef%drdz, &
      dsdx => coef%dsdx, dsdy => coef%dsdy, dsdz => coef%dsdz, &
      dtdx => coef%dtdx, dtdy => coef%dtdy, dtdz => coef%dtdz, &
      nelv => coef%msh%nelv, lx => Xh%lx, w3 => Xh%w3)
      nxyz = lx * lx * lx
      do e = 1, nelv
         do i = 1, nxyz
            cr(i,e) = w3(i,1,1) * (cx(i,e) * drdx(i,1,1,e) &
                        + cy(i,e) * drdy(i,1,1,e) &
                        + cz(i,e) * drdz(i,1,1,e))
            cs(i,e) = w3(i,1,1) * (cx(i,e) * dsdx(i,1,1,e) &
                        + cy(i,e) * dsdy(i,1,1,e) &
                        + cz(i,e) * dsdz(i,1,1,e))
            ct(i,e) = w3(i,1,1) * (cx(i,e) * dtdx(i,1,1,e) &
                        + cy(i,e) * dtdy(i,1,1,e) &
                        + cz(i,e) * dtdz(i,1,1,e))
         end do
      end do
    end associate

  end subroutine opr_xsmm_set_convect_rst

end module opr_xsmm
