! Copyright (c) 2020-2024, The Neko Authors
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
!> Operators
module operators
  use neko_config, only : NEKO_BCKND_SX, NEKO_BCKND_DEVICE, NEKO_BCKND_XSMM,&
                          NEKO_DEVICE_MPI
  use num_types, only : rp
  use opr_cpu, only : opr_cpu_cfl, opr_cpu_curl, opr_cpu_opgrad, opr_cpu_conv1,&
                      opr_cpu_cdtp, opr_cpu_dudxyz
  use opr_sx, only : opr_sx_cfl, opr_sx_curl, opr_sx_dudxyz, opr_sx_opgrad, &
                     opr_sx_cdtp, opr_sx_conv1
  use opr_xsmm, only : opr_xsmm_cdtp, opr_xsmm_conv1, opr_xsmm_curl, &
                       opr_xsmm_dudxyz, opr_xsmm_opgrad
  use opr_device, only : opr_device_cdtp, opr_device_cfl, opr_device_curl, &
                         opr_device_conv1, opr_device_dudxyz, &
                         opr_device_lambda2, opr_device_opgrad
  use space, only : space_t
  use coefs, only : coef_t
  use field, only : field_t
  use math, only : glsum, pi, cmult, add2, cadd
  use device, only : c_ptr, device_get_ptr
  use device_math, only : device_add2, device_cmult
  use comm
  implicit none
  private

  public :: dudxyz, opgrad, ortho, cdtp, conv1, curl, cfl,&
            lambda2op, strain_rate
  
contains

  !> Compute derivative of a scalar field along a single direction.
  !! @param du Holds the resulting derivative values.
  !! @param u The values of the field.
  !! @param dr The derivative of r with respect to the chosen direction.
  !! @param ds The derivative of s with respect to the chosen direction.
  !! @param dt The derivative of t with respect to the chosen direction.
  !! @param coef The SEM coefficients.
  subroutine dudxyz (du, u, dr, ds, dt, coef)
    type(coef_t), intent(in), target :: coef
    real(kind=rp), dimension(coef%Xh%lx,coef%Xh%ly,coef%Xh%lz,coef%msh%nelv), &
         intent(inout) ::  du
    real(kind=rp), dimension(coef%Xh%lx,coef%Xh%ly,coef%Xh%lz,coef%msh%nelv), &
         intent(in) ::  u, dr, ds, dt

    if (NEKO_BCKND_SX .eq. 1) then
       call opr_sx_dudxyz(du, u, dr, ds, dt, coef)
    else if (NEKO_BCKND_XSMM .eq. 1) then
       call opr_xsmm_dudxyz(du, u, dr, ds, dt, coef)
    else if (NEKO_BCKND_DEVICE .eq. 1) then
       call opr_device_dudxyz(du, u, dr, ds, dt, coef)
    else
       call opr_cpu_dudxyz(du, u, dr, ds, dt, coef)
    end if

  end subroutine dudxyz

  !> Compute the gradient of a scalar field.
  !! @details By providing `es` and `ee`, it is possible to compute only for a
  !! range of element indices.
  !! @param ux Will store the x component of the gradient.
  !! @param uy Will store the y component of the gradient.
  !! @param uz Will store the z component of the gradient.
  !! @param u The values of the field.
  !! @param coef The SEM coefficients.
  !! @param es Starting element index, optional, defaults to 1.
  !! @param ee Ending element index, optional, defaults to `nelv`.
  !! @note Equals wgradm1 in Nek5000, the weak form of the gradient.
  subroutine opgrad(ux, uy, uz, u, coef, es, ee)
    type(coef_t), intent(in) :: coef
    real(kind=rp), dimension(coef%Xh%lxyz,coef%msh%nelv), intent(inout) :: ux
    real(kind=rp), dimension(coef%Xh%lxyz,coef%msh%nelv), intent(inout) :: uy
    real(kind=rp), dimension(coef%Xh%lxyz,coef%msh%nelv), intent(inout) :: uz
    real(kind=rp), dimension(coef%Xh%lxyz,coef%msh%nelv), intent(in) :: u
    integer, optional :: es, ee
    integer :: eblk_start, eblk_end

    if (present(es)) then
       eblk_start = es
    else
       eblk_start = 1
    end if

    if (present(ee)) then
       eblk_end = ee
    else
       eblk_end = coef%msh%nelv
    end if

    if (NEKO_BCKND_SX .eq. 1) then
       call opr_sx_opgrad(ux, uy, uz, u, coef)
    else if (NEKO_BCKND_XSMM .eq. 1) then
       call opr_xsmm_opgrad(ux, uy, uz, u, coef)
    else if (NEKO_BCKND_DEVICE .eq. 1) then
       call opr_device_opgrad(ux, uy, uz, u, coef)
    else
       call opr_cpu_opgrad(ux, uy, uz, u, coef, eblk_start, eblk_end)
    end if

  end subroutine opgrad

  !> Othogonalize with regard to vector (1,1,1,1,1,1...,1)^T.
  !! @param x The vector to orthogonolize.
  !! @param n The size of `x`.
  !! @param glb_n The global number of elements of `x` across all MPI ranks. Be careful with overflow!
  subroutine ortho(x, n, glb_n)
    integer, intent(in) :: n
    integer, intent(in) :: glb_n
    real(kind=rp), dimension(n), intent(inout) :: x
    real(kind=rp) :: rlam

    rlam = glsum(x, n)/glb_n
    call cadd(x, -rlam, n)

  end subroutine ortho

  !> Apply D^T to a scalar field, where D is the derivative matrix.
  !! @param dtx Will store the result.
  !! @param x The values of the field.
  !! @param dr The derivative of r with respect to the chosen direction.
  !! @param ds The derivative of s with respect to the chosen direction.
  !! @param dt The derivative of t with respect to the chosen direction.
  !! @param coef The SEM coefficients.
  !> @note This needs to be revised... the loop over n1,n2 is probably
  !! unesccssary
  subroutine cdtp (dtx, x, dr, ds, dt, coef)
    type(coef_t), intent(in) :: coef
    real(kind=rp), dimension(coef%Xh%lxyz,coef%msh%nelv), intent(inout) :: dtx
    real(kind=rp), dimension(coef%Xh%lxyz,coef%msh%nelv), intent(inout) :: x
    real(kind=rp), dimension(coef%Xh%lxyz,coef%msh%nelv), intent(in) :: dr
    real(kind=rp), dimension(coef%Xh%lxyz,coef%msh%nelv), intent(in) :: ds
    real(kind=rp), dimension(coef%Xh%lxyz,coef%msh%nelv), intent(in) :: dt

    if (NEKO_BCKND_SX .eq. 1) then
       call opr_sx_cdtp(dtx, x, dr, ds, dt, coef)
    else if (NEKO_BCKND_XSMM .eq. 1) then
       call opr_xsmm_cdtp(dtx, x, dr, ds, dt, coef)
    else if (NEKO_BCKND_DEVICE .eq. 1) then
       call opr_device_cdtp(dtx, x, dr, ds, dt, coef)
    else
       call opr_cpu_cdtp(dtx, x, dr, ds, dt, coef)
    end if

  end subroutine cdtp

  !> Compute the advection term.
  !! @param du Holds the result.
  !! @param u The advected field.
  !! @param vx The x component of the advecting velocity.
  !! @param vy The y component of the advecting velocity.
  !! @param vz The z component of the advecting velocity.
  !! @param Xh The function space for the fields involved.
  !! @param coef The SEM coefficients.
  !! @param es Starting element index, defaults to 1.
  !! @param ee Last element index, defaults to mesh size.
  subroutine conv1(du, u, vx, vy, vz, Xh, coef, es, ee)
    type(space_t), intent(inout) :: Xh
    type(coef_t), intent(inout) :: coef
    real(kind=rp), intent(inout) :: du(Xh%lxyz,coef%msh%nelv)
    real(kind=rp), intent(inout) :: u(Xh%lx,Xh%ly,Xh%lz,coef%msh%nelv)
    real(kind=rp), intent(inout) :: vx(Xh%lx,Xh%ly,Xh%lz,coef%msh%nelv)
    real(kind=rp), intent(inout) :: vy(Xh%lx,Xh%ly,Xh%lz,coef%msh%nelv)
    real(kind=rp), intent(inout) :: vz(Xh%lx,Xh%ly,Xh%lz,coef%msh%nelv)
    integer, optional :: es, ee
    integer :: eblk_end, eblk_start

    associate(nelv => coef%msh%nelv, gdim => coef%msh%gdim)
      if (present(es)) then
         eblk_start = es
      else
         eblk_start = 1
      end if

      if (present(ee)) then
         eblk_end = ee
      else
         eblk_end = coef%msh%nelv
      end if

      if (NEKO_BCKND_SX .eq. 1) then
         call opr_sx_conv1(du, u, vx, vy, vz, Xh, coef, nelv, gdim)
      else if (NEKO_BCKND_XSMM .eq. 1) then
         call opr_xsmm_conv1(du, u, vx, vy, vz, Xh, coef, nelv, gdim)
      else if (NEKO_BCKND_DEVICE .eq. 1) then
         call opr_device_conv1(du, u, vx, vy, vz, Xh, coef, nelv, gdim)
      else
         call opr_cpu_conv1(du, u, vx, vy, vz, Xh, coef, eblk_start, eblk_end)
      end if
    end associate

  end subroutine conv1

  !! Compute the curl fo a vector field.
  !! @param w1 Will store the x component of the curl.
  !! @param w2 Will store the y component of the curl.
  !! @param w3 Will store the z component of the curl.
  !! @param u1 The x component of the vector field.
  !! @param u2 The y component of the vector field.
  !! @param u3 The z component of the vector field.
  !! @param work1 A temporary array for computations.
  !! @param work2 A temporary array for computations.
  !! @param coef The SEM coefficients.
  subroutine curl(w1, w2, w3, u1, u2, u3, work1, work2, coef)
    type(field_t), intent(inout) :: w1
    type(field_t), intent(inout) :: w2
    type(field_t), intent(inout) :: w3
    type(field_t), intent(inout) :: u1
    type(field_t), intent(inout) :: u2
    type(field_t), intent(inout) :: u3
    type(field_t), intent(inout) :: work1
    type(field_t), intent(inout) :: work2
    type(coef_t), intent(in)  :: coef

    if (NEKO_BCKND_SX .eq. 1) then
       call opr_sx_curl(w1, w2, w3, u1, u2, u3, work1, work2, coef)
    else if (NEKO_BCKND_XSMM .eq. 1) then
       call opr_xsmm_curl(w1, w2, w3, u1, u2, u3, work1, work2, coef)
    else if (NEKO_BCKND_DEVICE .eq. 1) then
       call opr_device_curl(w1, w2, w3, u1, u2, u3, work1, work2, coef)
    else
       call opr_cpu_curl(w1, w2, w3, u1, u2, u3, work1, work2, coef)
    end if

  end subroutine curl

  !! Compute the CFL number
  !! @param dt The timestep.
  !! @param u The x component of velocity.
  !! @param v The y component of velocity.
  !! @param w The z component of velocity.
  !! @param Xh The SEM function space.
  !! @param coef The SEM coefficients.
  !! @param nelv The total number of elements.
  !! @param gdim Number of geometric dimensions.
  function cfl(dt, u, v, w, Xh, coef, nelv, gdim)
    type(space_t), intent(in) :: Xh
    type(coef_t), intent(in) :: coef
    integer, intent(in) :: nelv, gdim
    real(kind=rp), intent(in) :: dt
    real(kind=rp), dimension(Xh%lx,Xh%ly,Xh%lz,nelv), intent(in) ::  u, v, w
    real(kind=rp) :: cfl
    integer :: ierr

    if (NEKO_BCKND_SX .eq. 1) then
       cfl = opr_sx_cfl(dt, u, v, w, Xh, coef, nelv, gdim)
    else if (NEKO_BCKND_DEVICE .eq. 1) then
       cfl = opr_device_cfl(dt, u, v, w, Xh, coef, nelv, gdim)
    else
       cfl = opr_cpu_cfl(dt, u, v, w, Xh, coef, nelv, gdim)
    end if

    if (.not. NEKO_DEVICE_MPI) then
       call MPI_Allreduce(MPI_IN_PLACE, cfl, 1, &
            MPI_REAL_PRECISION, MPI_MAX, NEKO_COMM, ierr)
    end if

  end function cfl

  !> Compute double the strain rate tensor, i.e du_i/dx_j + du_j/dx_i
  !! @param s11 Will hold the 1,1 component of the strain rate tensor.
  !! @param s22 Will hold the 2,2 component of the strain rate tensor.
  !! @param s33 Will hold the 3,3 component of the strain rate tensor.
  !! @param s12 Will hold the 1,2 component of the strain rate tensor.
  !! @param s13 Will hold the 1,3 component of the strain rate tensor.
  !! @param s23 Will hold the 2,3 component of the strain rate tensor.
  !! @param u The x component of velocity.
  !! @param v The y component of velocity.
  !! @param w The z component of velocity.
  !! @param coef The SEM coefficients.
  !! @note Similar to comp_sij in Nek5000.
  subroutine strain_rate(s11, s22, s33, s12, s13, s23, u, v, w, coef)
    type(field_t), intent(in) :: u, v, w !< velocity components
    type(coef_t), intent(in) :: coef
    real(kind=rp), intent(inout) :: s11(u%Xh%lx, u%Xh%ly, u%Xh%lz, u%msh%nelv)
    real(kind=rp), intent(inout) :: s22(u%Xh%lx, u%Xh%ly, u%Xh%lz, u%msh%nelv)
    real(kind=rp), intent(inout) :: s33(u%Xh%lx, u%Xh%ly, u%Xh%lz, u%msh%nelv)
    real(kind=rp), intent(inout) :: s12(u%Xh%lx, u%Xh%ly, u%Xh%lz, u%msh%nelv)
    real(kind=rp), intent(inout) :: s13(u%Xh%lx, u%Xh%ly, u%Xh%lz, u%msh%nelv)
    real(kind=rp), intent(inout) :: s23(u%Xh%lx, u%Xh%ly, u%Xh%lz, u%msh%nelv)

    type(c_ptr) :: s11_d, s22_d, s33_d, s12_d, s23_d, s13_d

    integer :: nelv, lxyz

    if (NEKO_BCKND_DEVICE .eq. 1) then
       s11_d = device_get_ptr(s11)
       s22_d = device_get_ptr(s22)
       s33_d = device_get_ptr(s33)
       s12_d = device_get_ptr(s12)
       s23_d = device_get_ptr(s23)
       s13_d = device_get_ptr(s13)
    endif

    nelv = u%msh%nelv
    lxyz = u%Xh%lxyz

    ! we use s11 as a work array here
    call dudxyz (s12, u%x, coef%drdy, coef%dsdy, coef%dtdy, coef)
    call dudxyz (s11, v%x, coef%drdx, coef%dsdx, coef%dtdx, coef)
    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_add2(s12_d, s11_d, nelv*lxyz)
    else
       call add2(s12, s11, nelv*lxyz)
    endif

    call dudxyz (s13, u%x, coef%drdz, coef%dsdz, coef%dtdz, coef)
    call dudxyz (s11, w%x, coef%drdx, coef%dsdx, coef%dtdx, coef)
    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_add2(s13_d, s11_d, nelv*lxyz)
    else
       call add2(s13, s11, nelv*lxyz)
    endif

    call dudxyz (s23, v%x, coef%drdz, coef%dsdz, coef%dtdz, coef)
    call dudxyz (s11, w%x, coef%drdy, coef%dsdy, coef%dtdy, coef)
    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_add2(s23_d, s11_d, nelv*lxyz)
    else
       call add2(s23, s11, nelv*lxyz)
    endif

    call dudxyz (s11, u%x, coef%drdx, coef%dsdx, coef%dtdx, coef)
    call dudxyz (s22, v%x, coef%drdy, coef%dsdy, coef%dtdy, coef)
    call dudxyz (s33, w%x, coef%drdz, coef%dsdz, coef%dtdz, coef)

    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_cmult(s11_d, 2.0_rp, nelv*lxyz)
       call device_cmult(s22_d, 2.0_rp, nelv*lxyz)
       call device_cmult(s33_d, 2.0_rp, nelv*lxyz)
    else
       call cmult(s11, 2.0_rp, nelv*lxyz)
       call cmult(s22, 2.0_rp, nelv*lxyz)
       call cmult(s33, 2.0_rp, nelv*lxyz)
    endif

  end subroutine strain_rate

  !> Compute the Lambda2 field for a given velocity field.
  !! @param lambda2 Holds the computed Lambda2 field.
  !! @param u The x-velocity.
  !! @param v The y-velocity.
  !! @param w the z-velocity.
  !! @param coef The SEM coefficients.
  subroutine lambda2op(lambda2, u, v, w, coef)
    type(coef_t), intent(in) :: coef
    type(field_t), intent(inout) :: lambda2
    type(field_t), intent(in) :: u, v, w
    real(kind=rp) :: grad(coef%Xh%lxyz,3,3)
    integer :: temp_indices(9), e, i, ind_sort(3)
    real(kind=rp) :: eigen(3), B, C, D, q, r, theta, l2
    real(kind=rp) :: s11, s22, s33, s12, s13, s23, o12, o13, o23
    real(kind=rp) :: a11, a22, a33, a12, a13, a23
    real(kind=rp) :: msk1, msk2, msk3

    if (NEKO_BCKND_DEVICE .eq. 1) then
       call opr_device_lambda2(lambda2, u, v, w, coef)
    else
       do e = 1, coef%msh%nelv
          call opgrad(grad(1,1,1), grad(1,1,2), grad(1,1,3), &
                              u%x(1,1,1,e),coef,e,e)
          call opgrad(grad(1,2,1), grad(1,2,2), grad(1,2,3), &
                              v%x(1,1,1,e),coef,e,e)
          call opgrad(grad(1,3,1), grad(1,3,2), grad(1,3,3), &
                              w%x(1,1,1,e),coef,e,e)

          do i = 1, coef%Xh%lxyz
             s11 = grad(i,1,1)
             s22 = grad(i,2,2)
             s33 = grad(i,3,3)


             s12 = 0.5*(grad(i,1,2) + grad(i,2,1))
             s13 = 0.5*(grad(i,1,3) + grad(i,3,1))
             s23 = 0.5*(grad(i,2,3) + grad(i,3,2))

             o12 = 0.5*(grad(i,1,2) - grad(i,2,1))
             o13 = 0.5*(grad(i,1,3) - grad(i,3,1))
             o23 = 0.5*(grad(i,2,3) - grad(i,3,2))

             a11 = s11*s11 + s12*s12 + s13*s13 - o12*o12 - o13*o13
             a12 = s11 * s12  +  s12 * s22  +  s13 * s23 - o13 * o23
             a13 = s11 * s13  +  s12 * s23  +  s13 * s33 + o12 * o23

             a22 = s12*s12 + s22*s22 + s23*s23 - o12*o12 - o23*o23
             a23 = s12 * s13 + s22 * s23 + s23 * s33 - o12 * o13
             a33 = s13*s13 + s23*s23 + s33*s33 - o13*o13 - o23*o23


             B = -(a11 + a22 + a33)
             C = -(a12*a12 + a13*a13 + a23*a23 &
                  - a11 * a22 - a11 * a33 - a22 * a33)
             D = -(2.0 * a12 * a13 * a23 - a11 * a23*a23 &
                  - a22 * a13*a13 - a33 * a12*a12  +  a11 * a22 * a33)


             q = (3.0 * C - B*B) / 9.0
             r = (9.0 * C * B - 27.0 * D - 2.0 * B*B*B) / 54.0
             theta = acos( r / sqrt(-q*q*q) )

             eigen(1) = 2.0 * sqrt(-q) * cos(theta / 3.0) - B / 3.0
             eigen(2) = 2.0 * sqrt(-q) * cos((theta + 2.0 * pi) / 3.0) - B / 3.0
             eigen(3) = 2.0 * sqrt(-q) * cos((theta + 4.0 * pi) / 3.0) - B / 3.0

             msk1 = merge(1.0_rp, 0.0_rp, eigen(2) .le. eigen(1) &
                          .and. eigen(1) .le. eigen(3) .or.  eigen(3) &
                          .le. eigen(1) .and. eigen(1) .le. eigen(2) )
             msk2 = merge(1.0_rp, 0.0_rp, eigen(1) .le. eigen(2) &
                          .and. eigen(2) .le. eigen(3) .or. eigen(3) &
                          .le. eigen(2) .and. eigen(2) .le. eigen(1))
             msk3 = merge(1.0_rp, 0.0_rp, eigen(1) .le. eigen(3) &
                          .and. eigen(3) .le. eigen(2) .or. eigen(2) &
                          .le. eigen(3) .and. eigen(3) .le. eigen(1))

             l2 = msk1 * eigen(1) + msk2 * eigen(2) + msk3 * eigen(3)

             lambda2%x(i,1,1,e) = l2/(coef%B(i,1,1,e)**2)
          end do
       end do
    end if

  end subroutine lambda2op

end module operators
