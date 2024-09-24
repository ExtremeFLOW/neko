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

module drag_torque
  use field, only :  field_t
  use coefs, only : coef_t
  use mesh, only : mesh_t
  use facet_zone, only : facet_zone_t
  use comm
  use math, only : rzero, masked_red_copy, col3, vdot3, cmult
  use space, only : space_t
  use num_types, only : rp
  use utils, only : nonlinear_index
  use device
  use device_math, only : device_cmult, device_masked_red_copy, &
                          device_col3, device_vdot3, device_rzero
  implicit none
  private
  !> Some functions to calculate the lift/drag and torque
  !! Calculation can be done on a zone, a facet, or a point
  !! Currently everything is CPU only
  public :: drag_torque_zone, drag_torque_facet, drag_torque_pt, setup_normals,&
            calc_force_array, device_calc_force_array

contains
  !> Calculate drag and torque over a zone.
  !! @param dgtq, the computed drag and torque
  !! @param tstep, the time step
  !! @param zone, the zone which we compute the drag and toqure over
  !! @param center, the point around which we calculate the torque
  !! @param s11-s23, the strain rate tensor
  !! @param p, the pressure
  !! @param coef, coefficents
  !! @param visc, the viscosity
  subroutine drag_torque_zone(dgtq, tstep, zone, center, s11, s22, s33, s12, s13, s23,&
                              p, coef, visc)
    integer, intent(in) :: tstep
    type(facet_zone_t) :: zone
    type(coef_t), intent(inout) :: coef
    real(kind=rp), intent(inout) :: s11(coef%Xh%lx,coef%Xh%lx,coef%Xh%lz,coef%msh%nelv)
    real(kind=rp), intent(inout) :: s22(coef%Xh%lx,coef%Xh%lx,coef%Xh%lz,coef%msh%nelv)
    real(kind=rp), intent(inout) :: s33(coef%Xh%lx,coef%Xh%lx,coef%Xh%lz,coef%msh%nelv)
    real(kind=rp), intent(inout) :: s12(coef%Xh%lx,coef%Xh%lx,coef%Xh%lz,coef%msh%nelv)
    real(kind=rp), intent(inout) :: s13(coef%Xh%lx,coef%Xh%lx,coef%Xh%lz,coef%msh%nelv)
    real(kind=rp), intent(inout) :: s23(coef%Xh%lx,coef%Xh%lx,coef%Xh%lz,coef%msh%nelv)
    type(field_t), intent(inout) :: p
    real(kind=rp), intent(in) :: visc, center(3)
    real(kind=rp) :: dgtq(3,4)
    real(kind=rp) :: dragpx = 0.0_rp ! pressure
    real(kind=rp) :: dragpy = 0.0_rp
    real(kind=rp) :: dragpz = 0.0_rp
    real(kind=rp) :: dragvx = 0.0_rp ! viscous
    real(kind=rp) :: dragvy = 0.0_rp
    real(kind=rp) :: dragvz = 0.0_rp
    real(kind=rp) :: torqpx = 0.0_rp ! pressure
    real(kind=rp) :: torqpy = 0.0_rp
    real(kind=rp) :: torqpz = 0.0_rp
    real(kind=rp) :: torqvx = 0.0_rp ! viscous
    real(kind=rp) :: torqvy = 0.0_rp
    real(kind=rp) :: torqvz = 0.0_rp
    real(kind=rp) :: dragx, dragy, dragz
    integer :: ie, ifc, mem, ierr
    dragx = 0.0
    dragy = 0.0
    dragz = 0.0

!
!     Fill up viscous array w/ default
!
      dragpx = 0.0
      dragpy = 0.0
      dragpz = 0.0
      dragvx = 0.0
      dragvy = 0.0
      dragvz = 0.0
      do mem  = 1,zone%size
         ie   = zone%facet_el(mem)%x(2)
         ifc   = zone%facet_el(mem)%x(1)
         call drag_torque_facet(dgtq,coef%dof%x,coef%dof%y,coef%dof%z,&
                                center,&
                                s11, s22, s33, s12, s13, s23,&
                                p%x,visc,ifc,ie, coef, coef%Xh)

         dragpx = dragpx + dgtq(1,1)  ! pressure
         dragpy = dragpy + dgtq(2,1)
         dragpz = dragpz + dgtq(3,1)

         dragvx = dragvx + dgtq(1,2)  ! viscous
         dragvy = dragvy + dgtq(2,2)
         dragvz = dragvz + dgtq(3,2)

         torqpx = torqpx + dgtq(1,3)  ! pressure
         torqpy = torqpy + dgtq(2,3)
         torqpz = torqpz + dgtq(3,3)

         torqvx = torqvx + dgtq(1,4)  ! viscous
         torqvy = torqvy + dgtq(2,4)
         torqvz = torqvz + dgtq(3,4)
      end do
!
!     Sum contributions from all processors
!
      call MPI_Allreduce(MPI_IN_PLACE,dragpx, 1, &
         MPI_REAL_PRECISION, MPI_SUM, NEKO_COMM, ierr)
      call MPI_Allreduce(MPI_IN_PLACE,dragpy, 1, &
         MPI_REAL_PRECISION, MPI_SUM, NEKO_COMM, ierr)
      call MPI_Allreduce(MPI_IN_PLACE,dragpz, 1, &
         MPI_REAL_PRECISION, MPI_SUM, NEKO_COMM, ierr)
      call MPI_Allreduce(MPI_IN_PLACE,dragvx, 1, &
         MPI_REAL_PRECISION, MPI_SUM, NEKO_COMM, ierr)
      call MPI_Allreduce(MPI_IN_PLACE,dragvy, 1, &
         MPI_REAL_PRECISION, MPI_SUM, NEKO_COMM, ierr)
      call MPI_Allreduce(MPI_IN_PLACE,dragvz, 1, &
         MPI_REAL_PRECISION, MPI_SUM, NEKO_COMM, ierr)
      !Torque
      call MPI_Allreduce(MPI_IN_PLACE,torqpx, 1, &
         MPI_REAL_PRECISION, MPI_SUM, NEKO_COMM, ierr)
      call MPI_Allreduce(MPI_IN_PLACE,torqpy, 1, &
         MPI_REAL_PRECISION, MPI_SUM, NEKO_COMM, ierr)
      call MPI_Allreduce(MPI_IN_PLACE,torqpz, 1, &
         MPI_REAL_PRECISION, MPI_SUM, NEKO_COMM, ierr)
      call MPI_Allreduce(MPI_IN_PLACE,torqvx, 1, &
         MPI_REAL_PRECISION, MPI_SUM, NEKO_COMM, ierr)
      call MPI_Allreduce(MPI_IN_PLACE,torqvy, 1, &
         MPI_REAL_PRECISION, MPI_SUM, NEKO_COMM, ierr)
      call MPI_Allreduce(MPI_IN_PLACE,torqvz, 1, &
         MPI_REAL_PRECISION, MPI_SUM, NEKO_COMM, ierr)

      dgtq(1,1) = dragpx  ! pressure
      dgtq(2,1) = dragpy
      dgtq(3,1) = dragpz

      dgtq(1,2) = dragvx  ! viscous
      dgtq(2,2) = dragvy
      dgtq(3,2) = dragvz

      dgtq(1,3) = torqpx  ! pressure
      dgtq(2,3) = torqpy
      dgtq(3,3) = torqpz

      dgtq(1,4) = torqvx  ! viscous
      dgtq(2,4) = torqvy
      dgtq(3,4) = torqvz

  end subroutine drag_torque_zone

  !> Calculate drag and torque over a facet.
  !! @param dgtq, the computed drag and torque
  !! @param tstep, the time step
  !! @param xm0, the x coords
  !! @param ym0, the y coords
  !! @param zm0, the z coords
  !! @param center, the point around which we calculate the torque
  !! @param s11-s23, the strain rate tensor
  !! @param p, the pressure
  !! @param coef, coefficents
  !! @param visc, the viscosity
  subroutine drag_torque_facet(dgtq,xm0,ym0,zm0, center,&
                               s11, s22, s33, s12, s13, s23,&
                               pm1,visc,f,e, coef, Xh)
    type(coef_t), intent(in) :: coef
    type(space_t), intent(in) :: Xh
    real(kind=rp), intent(out) :: dgtq(3,4)
    real(kind=rp), intent(in) :: center(3)
    real(kind=rp), intent(in) :: xm0 (Xh%lx,xh%ly,Xh%lz,coef%msh%nelv)
    real(kind=rp), intent(in) :: ym0 (Xh%lx,xh%ly,Xh%lz,coef%msh%nelv)
    real(kind=rp), intent(in) :: zm0 (Xh%lx,xh%ly,Xh%lz,coef%msh%nelv)
    real(kind=rp), intent(in) :: s11 (Xh%lx,xh%ly,Xh%lz,coef%msh%nelv)
    real(kind=rp), intent(in) :: s22 (Xh%lx,xh%ly,Xh%lz,coef%msh%nelv)
    real(kind=rp), intent(in) :: s33 (Xh%lx,xh%ly,Xh%lz,coef%msh%nelv)
    real(kind=rp), intent(in) :: s12 (Xh%lx,xh%ly,Xh%lz,coef%msh%nelv)
    real(kind=rp), intent(in) :: s13 (Xh%lx,xh%ly,Xh%lz,coef%msh%nelv)
    real(kind=rp), intent(in) :: s23 (Xh%lx,xh%ly,Xh%lz,coef%msh%nelv)
    real(kind=rp), intent(in) :: pm1 (Xh%lx,xh%ly,Xh%lz,coef%msh%nelv)
    real(kind=rp), intent(in) :: visc
    integer, intent(in) :: f,e
    integer :: pf, i, j1, j2
    real(kind=rp) ::    n1,n2,n3, a, v, dgtq_i(3,4)
    integer :: skpdat(6,6), NX, NY, NZ
    integer :: js1
    integer :: jf1
    integer :: jskip1
    integer :: js2
    integer :: jf2
    integer :: jskip2
    real(kind=rp) :: s11_, s12_, s22_, s13_, s23_, s33_


    NX = Xh%lx
    NY = Xh%ly
    NZ = Xh%lz
    SKPDAT(1,1)=1
    SKPDAT(2,1)=NX*(NY-1)+1
    SKPDAT(3,1)=NX
    SKPDAT(4,1)=1
    SKPDAT(5,1)=NY*(NZ-1)+1
    SKPDAT(6,1)=NY

    SKPDAT(1,2)=1             + (NX-1)
    SKPDAT(2,2)=NX*(NY-1)+1   + (NX-1)
    SKPDAT(3,2)=NX
    SKPDAT(4,2)=1
    SKPDAT(5,2)=NY*(NZ-1)+1
    SKPDAT(6,2)=NY

    SKPDAT(1,3)=1
    SKPDAT(2,3)=NX
    SKPDAT(3,3)=1
    SKPDAT(4,3)=1
    SKPDAT(5,3)=NY*(NZ-1)+1
    SKPDAT(6,3)=NY

    SKPDAT(1,4)=1           + NX*(NY-1)
    SKPDAT(2,4)=NX          + NX*(NY-1)
    SKPDAT(3,4)=1
    SKPDAT(4,4)=1
    SKPDAT(5,4)=NY*(NZ-1)+1
    SKPDAT(6,4)=NY

    SKPDAT(1,5)=1
    SKPDAT(2,5)=NX
    SKPDAT(3,5)=1
    SKPDAT(4,5)=1
    SKPDAT(5,5)=NY
    SKPDAT(6,5)=1

    SKPDAT(1,6)=1           + NX*NY*(NZ-1)
    SKPDAT(2,6)=NX          + NX*NY*(NZ-1)
    SKPDAT(3,6)=1
    SKPDAT(4,6)=1
    SKPDAT(5,6)=NY
    SKPDAT(6,6)=1
    pf = f
    js1    = skpdat(1,pf)
    jf1    = skpdat(2,pf)
    jskip1 = skpdat(3,pf)
    js2    = skpdat(4,pf)
    jf2    = skpdat(5,pf)
    jskip2 = skpdat(6,pf)
    call rzero(dgtq,12)
    i = 0
    a = 0
    do j2=js2,jf2,jskip2
       do j1=js1,jf1,jskip1
         i = i+1
         n1 = coef%nx(i,1,f,e)*coef%area(i,1,f,e)
         n2 = coef%ny(i,1,f,e)*coef%area(i,1,f,e)
         n3 = coef%nz(i,1,f,e)*coef%area(i,1,f,e)
         a  = a +          coef%area(i,1,f,e)
         v  = visc
         s11_ = s11(j1,j2,1,e)
         s12_ = s12(j1,j2,1,e)
         s22_ = s22(j1,j2,1,e)
         s13_ = s13(j1,j2,1,e)
         s23_ = s23(j1,j2,1,e)
         s33_ = s33(j1,j2,1,e)
         call drag_torque_pt(dgtq_i,xm0(j1,j2,1,e), ym0(j1,j2,1,e),zm0(j1,j2,1,e), center,&
                             s11_, s22_, s33_, s12_, s13_, s23_,&
                             pm1(j1,j2,1,e), n1, n2, n3, v)
         dgtq = dgtq + dgtq_i
       end do
    end do
  end subroutine drag_torque_facet

  !> Calculate drag and torque from one point
  !! @param dgtq, the computed drag and torque
  !! @param xm0, the x coord
  !! @param ym0, the y coord
  !! @param zm0, the z coord
  !! @param center, the point around which we calculate the torque
  !! @param s11-s23, the strain rate tensor
  !! @param p, the pressure
  !! @param n1, normal vector x
  !! @param n2, normal vector y
  !! @param n3, normal vector z
  !! @param v, the viscosity
  subroutine drag_torque_pt(dgtq, x, y, z, center, s11, s22, s33, s12, s13, s23,&
                            p, n1, n2, n3, v)
    real(kind=rp), intent(inout) :: dgtq(3,4)
    real(kind=rp), intent(in) :: x
    real(kind=rp), intent(in) :: y
    real(kind=rp), intent(in) :: z
    real(kind=rp), intent(in) :: p
    real(kind=rp), intent(in) :: v
    real(kind=rp), intent(in) :: n1, n2, n3, center(3)
    real(kind=rp), intent(in) :: s11, s12, s22, s13, s23, s33
    real(kind=rp) ::  s21, s31, s32, r1, r2, r3
    call rzero(dgtq,12)
    s21 = s12
    s32 = s23
    s31 = s13
    !pressure drag
    dgtq(1,1) = p*n1
    dgtq(2,1) = p*n2
    dgtq(3,1) = p*n3
    ! viscous drag
    dgtq(1,2) = -2*v*(s11*n1 + s12*n2 + s13*n3)
    dgtq(2,2) = -2*v*(s21*n1 + s22*n2 + s23*n3)
    dgtq(3,2) = -2*v*(s31*n1 + s32*n2 + s33*n3)
    r1 = x-center(1)
    r2 = y-center(2)
    r3 = z-center(3)
    !pressure torque
    dgtq(1,3) = (r2*dgtq(3,1)-r3*dgtq(2,1))
    dgtq(2,3) = (r3*dgtq(1,1)-r1*dgtq(3,1))
    dgtq(3,3) = (r1*dgtq(2,1)-r2*dgtq(1,1))
    !viscous torque
    dgtq(1,4) = (r2*dgtq(3,2)-r3*dgtq(2,2))
    dgtq(2,4) = (r3*dgtq(1,2)-r1*dgtq(3,2))
    dgtq(3,4) = (r1*dgtq(2,2)-r2*dgtq(1,2))
  end subroutine drag_torque_pt

  !> Calculate drag and torque from array of points
  !! @param force, the computed force
  !! @param s11-s23, the strain rate tensor
  !! @param p, the pressure
  !! @param n1, normal vector x
  !! @param n2, normal vector y
  !! @param n3, normal vector z
  !! @param v, the viscosity
  !! @param n_pts, the number of points
  subroutine calc_force_array(force1, force2, force3, force4, force5, force6,&
                              s11, s22, s33, s12, s13, s23,&
                              p, n1, n2, n3, v, n_pts)
    integer :: n_pts
    real(kind=rp), intent(inout),dimension(n_pts) :: force1, force2, force3
    real(kind=rp), intent(inout),dimension(n_pts) :: force4, force5, force6
    real(kind=rp), intent(in) :: p(n_pts)
    real(kind=rp), intent(in) :: v
    real(kind=rp), intent(in) :: n1(n_pts), n2(n_pts), n3(n_pts)
    real(kind=rp), intent(in), dimension(n_pts) :: s11, s12, s22, s13, s23, s33
    real(kind=rp) :: v2
    call rzero(force4,n_pts)
    call rzero(force5,n_pts)
    call rzero(force6,n_pts)
    !pressure force
    call col3(force1,p,n1,n_pts)
    call col3(force2,p,n2,n_pts)
    call col3(force3,p,n3,n_pts)
    ! viscous force
    v2 = -2.0_rp*v
    call vdot3(force4,s11,s12,s13,n1,n2,n3,n_pts)
    call vdot3(force5,s12,s22,s23,n1,n2,n3,n_pts)
    call vdot3(force6,s13,s23,s33,n1,n2,n3,n_pts)
    call cmult(force4,v2,n_pts)
    call cmult(force5,v2,n_pts)
    call cmult(force6,v2,n_pts)
  end subroutine calc_force_array
  !> Calculate drag and torque from array of points
  !! @param force, the computed force
  !! @param s11-s23, the strain rate tensor
  !! @param p, the pressure
  !! @param n1, normal vector x
  !! @param n2, normal vector y
  !! @param n3, normal vector z
  !! @param v, the viscosity
  !! @param n_pts, the number of points
  subroutine device_calc_force_array(force1, force2, force3,&
                                     force4, force5, force6,&
                                     s11, s22, s33, s12, s13, s23,&
                                     p, n1, n2, n3, v, n_pts)
    integer :: n_pts
    type(c_ptr) :: force1, force2, force3
    type(c_ptr) :: force4, force5, force6
    type(c_ptr) :: p, n1, n2, n3
    type(c_ptr) :: s11, s12, s22, s13, s23, s33
    real(kind=rp) :: v
    real(kind=rp) :: v2
    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_rzero(force4,n_pts)
       call device_rzero(force5,n_pts)
       call device_rzero(force6,n_pts)
       !pressure force
       call device_col3(force1,p,n1,n_pts)
       call device_col3(force2,p,n2,n_pts)
       call device_col3(force3,p,n3,n_pts)
       ! viscous force
       v2 = -2.0_rp*v
       call device_vdot3(force4,s11,s12,s13,n1,n2,n3,n_pts)
       call device_vdot3(force5,s12,s22,s23,n1,n2,n3,n_pts)
       call device_vdot3(force6,s13,s23,s33,n1,n2,n3,n_pts)
       call device_cmult(force4,v2,n_pts)
       call device_cmult(force5,v2,n_pts)
       call device_cmult(force6,v2,n_pts)
    else
       call neko_error('error in drag_torque, no device bcklnd configured')
    end if
  end subroutine device_calc_force_array


  subroutine setup_normals(coef,mask,facets,n1,n2,n3,n_pts)
    type(coef_t) :: coef
    integer :: n_pts
    real(kind=rp), dimension(n_pts) :: n1, n2, n3
    integer :: mask(0:n_pts), facets(0:n_pts), fid, idx(4)
    real(kind=rp) :: normal(3), area(3)
    integer :: i
    
    do i = 1, n_pts
       fid = facets(i)
       idx = nonlinear_index(mask(i), coef%Xh%lx, coef%Xh%lx,&
                             coef%Xh%lx)
       normal = coef%get_normal(idx(1), idx(2), idx(3), idx(4), fid)
       area = coef%get_area(idx(1), idx(2), idx(3), idx(4), fid)
       n1(i) = normal(1)*area(1)
       n2(i) = normal(2)*area(2)
       n3(i) = normal(3)*area(3)

    end do

  end subroutine setup_normals

end module drag_torque
