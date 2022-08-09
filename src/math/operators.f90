! Copyright (c) 2020-2022, The Neko Authors
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
  use neko_config
  use num_types  
  use opr_cpu
  use opr_sx
  use opr_xsmm
  use opr_device
  use space  
  use coefs
  use field
  use math
  implicit none
  private

  public :: dudxyz, opgrad, ortho, cdtp, conv1, curl, cfl
  
contains
  
  !> Compute dU/dx or dU/dy or dU/dz 
  subroutine dudxyz (du,u,dr,ds,dt,coef)
    type(coef_t), intent(in), target :: coef
    real(kind=rp), dimension(coef%Xh%lx,coef%Xh%ly,coef%Xh%lz,coef%msh%nelv), &
         intent(inout) ::  du
    real(kind=rp), dimension(coef%Xh%lx,coef%Xh%ly,coef%Xh%lz,coef%msh%nelv), &
         intent(in) ::  u, dr, ds, dt

    if (NEKO_BCKND_SX .eq. 1) then 
       call opr_sx_dudxyz(du, u, dr, ds, dt, coef)
    else if (NEKO_BCKND_XSMM .eq. 1) then
       call opr_xsmm_dudxyz(du, u, dr, ds, dt, coef)
    else if ((NEKO_BCKND_HIP .eq. 1) .or. (NEKO_BCKND_CUDA .eq. 1) &
         .or. (NEKO_BCKND_OPENCL .eq. 1)) then       
       call opr_device_dudxyz(du, u, dr, ds, dt, coef)
    else
       call opr_cpu_dudxyz(du, u, dr, ds, dt, coef)
    end if
    
  end subroutine dudxyz

  !> Equals wgradm1 in nek5000. Gradient of velocity vectors.
  subroutine opgrad(ux,uy,uz,u,coef, es, ee) ! weak form of grad     
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
    else if ((NEKO_BCKND_HIP .eq. 1) .or. (NEKO_BCKND_CUDA .eq. 1) &
         .or. (NEKO_BCKND_OPENCL .eq. 1)) then       
       call opr_device_opgrad(ux, uy, uz, u, coef)
    else
       call opr_cpu_opgrad(ux, uy, uz, u, coef, eblk_start, eblk_end)
    end if
    
  end subroutine opgrad
  
  !> Othogonalize with regard to vector (1,1,1,1,1,1...,1)^T.
  subroutine ortho(x,n ,glb_n)
    integer, intent(in) :: n
    integer, intent(in) :: glb_n
    real(kind=rp), dimension(n), intent(inout) :: x
    real(kind=rp) :: rlam

    rlam = glsum(x,n)/glb_n
    call cadd(x,-rlam,n)

  end subroutine ortho
  
  !> Compute DT*X (entire field)
  !> This needs to be revised... the loop over n1,n2 is probably unesccssary
  subroutine cdtp (dtx,x,dr,ds,dt, coef)
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
    else if ((NEKO_BCKND_HIP .eq. 1) .or. (NEKO_BCKND_CUDA .eq. 1) &
         .or. (NEKO_BCKND_OPENCL .eq. 1)) then       
       call opr_device_cdtp(dtx, x, dr, ds, dt, coef)
    else
       call opr_cpu_cdtp(dtx, x, dr, ds, dt, coef)
    end if
    
  end subroutine cdtp
   
  subroutine conv1(du,u, vx, vy, vz, Xh, coef, es, ee)
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
      else if ((NEKO_BCKND_HIP .eq. 1) .or. (NEKO_BCKND_CUDA .eq. 1) &
           .or. (NEKO_BCKND_OPENCL .eq. 1)) then
         call opr_device_conv1(du, u, vx, vy, vz, Xh, coef, nelv, gdim)
      else
         call opr_cpu_conv1(du, u, vx, vy, vz, Xh, coef, eblk_start, eblk_end)
      end if
    end associate

  end subroutine conv1

  subroutine curl(w1, w2, w3, u1, u2, u3, work1, work2, c_Xh)
    type(field_t), intent(inout) :: w1
    type(field_t), intent(inout) :: w2
    type(field_t), intent(inout) :: w3
    type(field_t), intent(inout) :: u1
    type(field_t), intent(inout) :: u2
    type(field_t), intent(inout) :: u3
    type(field_t), intent(inout) :: work1
    type(field_t), intent(inout) :: work2
    type(coef_t), intent(in)  :: c_Xh

    if (NEKO_BCKND_SX .eq. 1) then
       call opr_sx_curl(w1, w2, w3, u1, u2, u3, work1, work2, c_Xh)
    else if (NEKO_BCKND_XSMM .eq. 1) then
       call opr_xsmm_curl(w1, w2, w3, u1, u2, u3, work1, work2, c_Xh)
    else if ((NEKO_BCKND_HIP .eq. 1) .or. (NEKO_BCKND_CUDA .eq. 1) &
         .or. (NEKO_BCKND_OPENCL .eq. 1)) then       
       call opr_device_curl(w1, w2, w3, u1, u2, u3, work1, work2, c_Xh)
    else
       call opr_cpu_curl(w1, w2, w3, u1, u2, u3, work1, work2, c_Xh)
    end if

  end subroutine curl

  function cfl(dt, u, v, w, Xh, coef, nelv, gdim)
    type(space_t) :: Xh
    type(coef_t) :: coef
    integer :: nelv, gdim
    real(kind=rp) :: dt
    real(kind=rp), dimension(Xh%lx,Xh%ly,Xh%lz,nelv) ::  u, v, w
    real(kind=rp) :: cfl
    integer :: ierr

    if (NEKO_BCKND_SX .eq. 1) then
       cfl = opr_sx_cfl(dt, u, v, w, Xh, coef, nelv, gdim)
    else if ((NEKO_BCKND_HIP .eq. 1) .or. (NEKO_BCKND_CUDA .eq. 1) &
         .or. (NEKO_BCKND_OPENCL .eq. 1)) then  
       cfl = opr_device_cfl(dt, u, v, w, Xh, coef, nelv, gdim)
    else
       cfl = opr_cpu_cfl(dt, u, v, w, Xh, coef, nelv, gdim)
    end if

    call MPI_Allreduce(MPI_IN_PLACE, cfl, 1, &
         MPI_REAL_PRECISION, MPI_MAX, NEKO_COMM, ierr)
    
  end function cfl
  
end module operators
