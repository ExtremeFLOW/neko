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
  use neko_config, only : NEKO_BCKND_SX, NEKO_BCKND_DEVICE, NEKO_BCKND_XSMM, &
                          NEKO_DEVICE_MPI
  use num_types, only : rp
  use opr_cpu, only : opr_cpu_cfl, opr_cpu_curl, opr_cpu_opgrad, &
                      opr_cpu_conv1, opr_cpu_convect_scalar, opr_cpu_cdtp, &
                      opr_cpu_dudxyz, opr_cpu_lambda2, opr_cpu_set_convect_rst
  use opr_sx, only : opr_sx_cfl, opr_sx_curl, opr_sx_opgrad, &
                     opr_sx_conv1, opr_sx_convect_scalar, opr_sx_cdtp, &
                     opr_sx_dudxyz, opr_sx_lambda2, opr_sx_set_convect_rst
  use opr_xsmm, only : opr_xsmm_cdtp, opr_xsmm_conv1, opr_xsmm_curl, &
                       opr_xsmm_dudxyz, opr_xsmm_opgrad, &
                       opr_xsmm_convect_scalar, opr_xsmm_set_convect_rst
  use opr_device, only : opr_device_cdtp, opr_device_cfl, opr_device_curl, &
                         opr_device_conv1, opr_device_dudxyz, &
                         opr_device_lambda2, opr_device_opgrad
  use space, only : space_t
  use coefs, only : coef_t
  use field, only : field_t
  use interpolation, only : interpolator_t
  use math, only : glsum, cmult, add2, add3s2, cadd, copy, col2, invcol2, &
                   invcol3, rzero
  use device, only : c_ptr, device_get_ptr
  use device_math, only : device_add2, device_cmult, device_copy
  use scratch_registry, only : neko_scratch_registry
  use comm
  implicit none
  private

  public :: dudxyz, opgrad, ortho, cdtp, conv1, curl, cfl, &
            lambda2op, strain_rate, div, grad, set_convect_rst, runge_kutta

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
    real(kind=rp), dimension(coef%Xh%lx, coef%Xh%ly, coef%Xh%lz, &
         coef%msh%nelv), intent(inout) :: du
    real(kind=rp), dimension(coef%Xh%lx, coef%Xh%ly, coef%Xh%lz, &
         coef%msh%nelv), intent(in) :: u, dr, ds, dt

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

  !> Compute the divergence of a vector field.
  !! @param res Holds the resulting divergence values.
  !! @param ux The x component  of the vector field.
  !! @param uy The y component  of the vector field.
  !! @param uz The z component  of the vector field.
  !! @param coef The SEM coefficients.
  subroutine div(res, ux, uy, uz, coef)
    type(coef_t), intent(in), target :: coef
    real(kind=rp), dimension(coef%Xh%lx, coef%Xh%ly, coef%Xh%lz, &
         coef%msh%nelv), intent(inout) :: res
    real(kind=rp), dimension(coef%Xh%lx, coef%Xh%ly, coef%Xh%lz, &
         coef%msh%nelv), intent(in) :: ux, uy, uz
    type(field_t), pointer :: work
    integer :: ind
    type(c_ptr) :: res_d

    if (NEKO_BCKND_DEVICE .eq. 1) then
       res_d = device_get_ptr(res)
    end if

    call neko_scratch_registry%request_field(work, ind)

    ! Get dux / dx
    call dudxyz(res, ux, coef%drdx, coef%dsdx, coef%dtdx, coef)

    ! Get duy / dy
    call dudxyz(work%x, uy, coef%drdy, coef%dsdy, coef%dtdy, coef)
    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_add2(res_d, work%x_d, work%size())
    else
       call add2(res, work%x, work%size())
    end if

    ! Get dux / dz
    call dudxyz(work%x, uz, coef%drdz, coef%dsdz, coef%dtdz, coef)
    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_add2(res_d, work%x_d, work%size())
    else
       call add2(res, work%x, work%size())
    end if

    call neko_scratch_registry%relinquish_field(ind)

  end subroutine div

  !> Compute the gradient of a scalar field.
  !! @param ux Will store the x component of the gradient.
  !! @param uy Will store the y component of the gradient.
  !! @param uz Will store the z component of the gradient.
  !! @param u The values of the field.
  !! @param coef The SEM coefficients.
  subroutine grad(ux, uy, uz, u, coef)
    type(coef_t), intent(in) :: coef
    real(kind=rp), dimension(coef%Xh%lxyz, coef%msh%nelv), intent(inout) :: ux
    real(kind=rp), dimension(coef%Xh%lxyz, coef%msh%nelv), intent(inout) :: uy
    real(kind=rp), dimension(coef%Xh%lxyz, coef%msh%nelv), intent(inout) :: uz
    real(kind=rp), dimension(coef%Xh%lxyz, coef%msh%nelv), intent(in) :: u

    call dudxyz(ux, u, coef%drdx, coef%dsdx, coef%dtdx, coef)
    call dudxyz(uy, u, coef%drdy, coef%dsdy, coef%dtdy, coef)
    call dudxyz(uz, u, coef%drdz, coef%dsdz, coef%dtdz, coef)

  end subroutine grad

  !> Compute the weak gradient of a scalar field, i.e. the gradient multiplied
  !! by the mass matrix.
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
    real(kind=rp), dimension(coef%Xh%lxyz, coef%msh%nelv), intent(inout) :: ux
    real(kind=rp), dimension(coef%Xh%lxyz, coef%msh%nelv), intent(inout) :: uy
    real(kind=rp), dimension(coef%Xh%lxyz, coef%msh%nelv), intent(inout) :: uz
    real(kind=rp), dimension(coef%Xh%lxyz, coef%msh%nelv), intent(in) :: u
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
  !! @param es Starting element index, optional, defaults to 1.
  !! @param ee Ending element index, optional, defaults to `nelv`.
  !> @note This needs to be revised... the loop over n1,n2 is probably
  !! unesccssary
  subroutine cdtp (dtx, x, dr, ds, dt, coef, es, ee)
    type(coef_t), intent(in) :: coef
    real(kind=rp), dimension(coef%Xh%lxyz, coef%msh%nelv), intent(inout) :: dtx
    real(kind=rp), dimension(coef%Xh%lxyz, coef%msh%nelv), intent(inout) :: x
    real(kind=rp), dimension(coef%Xh%lxyz, coef%msh%nelv), intent(in) :: dr
    real(kind=rp), dimension(coef%Xh%lxyz, coef%msh%nelv), intent(in) :: ds
    real(kind=rp), dimension(coef%Xh%lxyz, coef%msh%nelv), intent(in) :: dt
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
       call opr_sx_cdtp(dtx, x, dr, ds, dt, coef)
    else if (NEKO_BCKND_XSMM .eq. 1) then
       call opr_xsmm_cdtp(dtx, x, dr, ds, dt, coef)
    else if (NEKO_BCKND_DEVICE .eq. 1) then
       call opr_device_cdtp(dtx, x, dr, ds, dt, coef)
    else
       call opr_cpu_cdtp(dtx, x, dr, ds, dt, coef, eblk_start, eblk_end)
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
    real(kind=rp), intent(inout) :: du(Xh%lxyz, coef%msh%nelv)
    real(kind=rp), intent(inout) :: u(Xh%lx, Xh%ly, Xh%lz, coef%msh%nelv)
    real(kind=rp), intent(inout) :: vx(Xh%lx, Xh%ly, Xh%lz, coef%msh%nelv)
    real(kind=rp), intent(inout) :: vy(Xh%lx, Xh%ly, Xh%lz, coef%msh%nelv)
    real(kind=rp), intent(inout) :: vz(Xh%lx, Xh%ly, Xh%lz, coef%msh%nelv)
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
         call opr_sx_conv1(du, u, vx, vy, vz, Xh, coef, nelv)
      else if (NEKO_BCKND_XSMM .eq. 1) then
         call opr_xsmm_conv1(du, u, vx, vy, vz, Xh, coef, nelv, gdim)
      else if (NEKO_BCKND_DEVICE .eq. 1) then
         call opr_device_conv1(du, u, vx, vy, vz, Xh, coef, nelv, gdim)
      else
         call opr_cpu_conv1(du, u, vx, vy, vz, Xh, coef, eblk_start, eblk_end)
      end if
    end associate

  end subroutine conv1

  !> Apply the convecting velocity c to the to the scalar field u, used in the OIFS scheme.
  !! @param du Holds the result
  !! @param c The convecting velocity
  !! @param u The convected scalar field
  !! @param Xh_GLL The GLL space used in simulation
  !! @param Xh_GL The GL space used for dealiasing
  !! @param coef The coefficients of the original space in simulation
  !! @param coef_GL The coefficients of the GL space used for dealiasing
  !! @param GLL_to_GL the interpolator between the GLL and GL spaces
  !! @note This subroutine is equal to the convop_fst_3d of the NEK5000.
  !! @note This subroutine is used specifically in the OIFS scheme, calculateing eq(17)
  !! in https://publications.anl.gov/anlpubs/2017/12/140626.pdf. The convecting term is
  !! calculated in the rst format and the GL grid. Then converted back to the GLL grid,
  !! going through an ADD gatter scatter operation at the element boundaries,
  !! before being multiplied by inverse of mass matrix.
  subroutine convect_scalar(du, u, c, Xh_GLL, Xh_GL, coef_GLL, &
                            coef_GL, GLL_to_GL)
    type(space_t), intent(in) :: Xh_GL
    type(space_t), intent(in) :: Xh_GLL
    type(coef_t), intent(in) :: coef_GLL
    type(coef_t), intent(in) :: coef_GL
    type(interpolator_t), intent(inout) :: GLL_to_GL
    real(kind=rp), intent(inout) :: &
                   du(Xh_GLL%lx, Xh_GLL%ly, Xh_GLL%lz, coef_GL%msh%nelv)
    real(kind=rp), intent(inout) :: &
                   u(Xh_GL%lx, Xh_GL%lx, Xh_GL%lx, coef_GL%msh%nelv)
    real(kind=rp), intent(inout) :: c(Xh_GL%lxyz, coef_GL%msh%nelv, 3)

    if (NEKO_BCKND_SX .eq. 1) then
       call opr_sx_convect_scalar(du, u, c, Xh_GLL, Xh_GL, &
                               coef_GLL, coef_GL, GLL_to_GL)
    else if (NEKO_BCKND_XSMM .eq. 1) then
       call opr_xsmm_convect_scalar(du, u, c, Xh_GLL, Xh_GL, &
                                 coef_GLL, coef_GL, GLL_to_GL)
    else
       call opr_cpu_convect_scalar(du, u, c, Xh_GLL, Xh_GL, &
                                coef_GLL, coef_GL, GLL_to_GL)
    end if

  end subroutine convect_scalar

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
    type(coef_t), intent(in) :: coef

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
    real(kind=rp), dimension(Xh%lx, Xh%ly, Xh%lz, nelv), intent(in) :: u, v, w
    real(kind=rp) :: cfl
    integer :: ierr

    if (NEKO_BCKND_SX .eq. 1) then
       cfl = opr_sx_cfl(dt, u, v, w, Xh, coef, nelv)
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

  !> Compute the strain rate tensor, i.e 0.5 * du_i/dx_j + du_j/dx_i
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
    end if

    nelv = u%msh%nelv
    lxyz = u%Xh%lxyz

    ! we use s11 as a work array here
    call dudxyz (s12, u%x, coef%drdy, coef%dsdy, coef%dtdy, coef)
    call dudxyz (s11, v%x, coef%drdx, coef%dsdx, coef%dtdx, coef)
    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_add2(s12_d, s11_d, nelv*lxyz)
    else
       call add2(s12, s11, nelv*lxyz)
    end if

    call dudxyz (s13, u%x, coef%drdz, coef%dsdz, coef%dtdz, coef)
    call dudxyz (s11, w%x, coef%drdx, coef%dsdx, coef%dtdx, coef)
    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_add2(s13_d, s11_d, nelv*lxyz)
    else
       call add2(s13, s11, nelv*lxyz)
    end if

    call dudxyz (s23, v%x, coef%drdz, coef%dsdz, coef%dtdz, coef)
    call dudxyz (s11, w%x, coef%drdy, coef%dsdy, coef%dtdy, coef)
    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_add2(s23_d, s11_d, nelv*lxyz)
    else
       call add2(s23, s11, nelv*lxyz)
    end if

    call dudxyz (s11, u%x, coef%drdx, coef%dsdx, coef%dtdx, coef)
    call dudxyz (s22, v%x, coef%drdy, coef%dsdy, coef%dtdy, coef)
    call dudxyz (s33, w%x, coef%drdz, coef%dsdz, coef%dtdz, coef)

    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_cmult(s12_d, 0.5_rp, nelv*lxyz)
       call device_cmult(s13_d, 0.5_rp, nelv*lxyz)
       call device_cmult(s23_d, 0.5_rp, nelv*lxyz)
    else
       call cmult(s12, 0.5_rp, nelv*lxyz)
       call cmult(s13, 0.5_rp, nelv*lxyz)
       call cmult(s23, 0.5_rp, nelv*lxyz)
    end if

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

    if (NEKO_BCKND_SX .eq. 1) then
       call opr_sx_lambda2(lambda2, u, v, w, coef)
    else if (NEKO_BCKND_DEVICE .eq. 1) then
       call opr_device_lambda2(lambda2, u, v, w, coef)
    else
       call opr_cpu_lambda2(lambda2, u, v, w, coef)
    end if

  end subroutine lambda2op

  !> Transforms the convecting velocity field to the rst form of the GL space
  !! @param cr convecting velocity in r-direction
  !! @param cs convecting velocity in s-direction
  !! @param ct convecting velocity in t-direction
  !! @param cx convecting velocity in x-direction
  !! @param cy convecting velocity in y-direction
  !! @param cz convecting velocity in z-direction
  !! @param Xh The GL space used for dealiasing
  !! @param coef The coeffiecients of the GL space used for dealiasing
  !! @note This subroutine is equal to the set_convect_new subroutine of NEK5000
  subroutine set_convect_rst(cr, cs, ct, cx, cy, cz, Xh, coef)
    type(space_t), intent(inout) :: Xh
    type(coef_t), intent(inout) :: coef
    real(kind=rp), dimension(Xh%lxyz, coef%msh%nelv), &
                   intent(inout) :: cr, cs, ct
    real(kind=rp), dimension(Xh%lxyz, coef%msh%nelv), &
                   intent(in) :: cx, cy, cz

    if (NEKO_BCKND_SX .eq. 1) then
       call opr_sx_set_convect_rst(cr, cs, ct, cx, cy, cz, Xh, coef)
    else if (NEKO_BCKND_XSMM .eq. 1) then
       call opr_xsmm_set_convect_rst(cr, cs, ct, cx, cy, cz, Xh, coef)
    else
       call opr_cpu_set_convect_rst(cr, cs, ct, cx, cy, cz, Xh, coef)
    end if

  end subroutine set_convect_rst

  !> Compute one step of Runge Kutta time interpolation for OIFS scheme
  !! @param phi The iterpolated field
  !! @param c_r1 The covecting velocity for the first stage
  !! @param c_r23 The convecting velocity for the second and third stage
  !! @param c_r4 The convecting velocity for the fourth stage
  !! @param Xh_GLL The GLL space used in simulation
  !! @param Xh_GL The GL space used for dealiasing
  !! @param coef The coefficients of the original space in simulation
  !! @param coef_GL The coefficients of the GL space used for dealiasing
  !! @param GLL_to_GL the interpolator between the GLL and GL spaces
  !! @param tau The the starting time
  !! @param dtau The time step used for the Runge Kutta scheme
  !! @param n size of phi
  !! @param nel Total number of elements
  !! @param n_GL the size in the GL space
  subroutine runge_kutta(phi, c_r1, c_r23, c_r4, Xh_GLL, Xh_GL, coef, &
                         coef_GL, GLL_to_GL, tau, dtau, n, nel, n_GL)
    type(space_t), intent(inout) :: Xh_GLL
    type(space_t), intent(inout) :: Xh_GL
    type(coef_t), intent(inout) :: coef
    type(coef_t), intent(inout) :: coef_GL
    type(interpolator_t) :: GLL_to_GL
    real(kind=rp), intent(inout) :: tau, dtau
    integer, intent(in) :: n, nel, n_GL
    real(kind=rp), dimension(n), intent(inout) :: phi
    real(kind=rp), dimension(3 * n_GL), intent(inout) :: c_r1, c_r23, c_r4
    real(kind=rp) :: c1, c2, c3
    ! Work Arrays
    real(kind=rp), dimension(n) ::  u1, r1, r2, r3, r4
    real(kind=rp), dimension(n_GL) :: u1_GL
    integer :: i, e

    c1 = 1.
    c2 = -dtau/2.
    c3 = -dtau

    ! Stage 1:
    call invcol3 (u1, phi, coef%B, n)
    call GLL_to_GL%map(u1_GL, u1, nel, Xh_GL)
    call convect_scalar(r1, u1_GL, c_r1, Xh_GLL, Xh_GL, coef, &
                        coef_GL, GLL_to_GL)
    call col2(r1, coef%B, n)

    ! Stage 2:
    call add3s2 (u1, phi, r1, c1, c2, n)
    call invcol2 (u1, coef%B, n)
    call GLL_to_GL%map(u1_GL, u1, nel, Xh_GL)
    call convect_scalar(r2, u1_GL, c_r23, Xh_GLL, Xh_GL, coef, &
                        coef_GL, GLL_to_GL)
    call col2(r2, coef%B, n)

    ! Stage 3:
    call add3s2 (u1, phi, r2, c1, c2, n)
    call invcol2 (u1,  coef%B, n)
    call GLL_to_GL%map(u1_GL, u1, nel, Xh_GL)
    call convect_scalar(r3, u1_GL, c_r23, Xh_GLL, Xh_GL, coef, &
                        coef_GL, GLL_to_GL)
    call col2(r3, coef%B, n)

    ! Stage 4:
    call add3s2 (u1, phi, r3, c1, c3, n)
    call invcol2 (u1, coef%B, n)
    call GLL_to_GL%map(u1_GL, u1, nel, Xh_GL)
    call convect_scalar(r4, u1_GL, c_r4, Xh_GLL, Xh_GL, coef, &
                        coef_GL, GLL_to_GL)
    call col2(r4, coef%B, n)

    c1 = -dtau/6.
    c2 = -dtau/3.
    do i = 1, n
       phi(i) = phi(i) + c1 * (r1(i) + r4(i)) + c2 * (r2(i) + r3(i))
    end do

  end subroutine runge_kutta

end module operators
