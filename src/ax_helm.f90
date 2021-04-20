module ax_helm
  use ax_product
  implicit none

  type, public, extends(ax_t) :: ax_helm_t
  contains
     procedure, nopass :: compute => ax_helm_compute
  end type ax_helm_t

contains 
  
  subroutine ax_helm_compute(w, u, coef, msh, Xh)
    type(mesh_t), intent(inout) :: msh
    type(space_t), intent(inout) :: Xh
    type(coef_t), intent(inout) :: coef
    real(kind=rp), intent(inout) :: w(Xh%lx, Xh%ly, Xh%lz, msh%nelv)
    real(kind=rp), intent(inout) :: u(Xh%lx, Xh%ly, Xh%lz, msh%nelv)
  
    real(kind=rp) :: dudr(Xh%lx,Xh%ly,Xh%lz)
    real(kind=rp) :: duds(Xh%lx,Xh%ly,Xh%lz)
    real(kind=rp) :: dudt(Xh%lx,Xh%ly,Xh%lz)
    real(kind=rp) :: tmp1(Xh%lx,Xh%ly,Xh%lz)
    real(kind=rp) :: tmp2(Xh%lx,Xh%ly,Xh%lz)
    real(kind=rp) :: tmp3(Xh%lx,Xh%ly,Xh%lz)
    real(kind=rp) :: tm1(Xh%lx,Xh%ly,Xh%lz)
    real(kind=rp) :: tm2(Xh%lx,Xh%ly,Xh%lz)
    real(kind=rp) :: tm3(Xh%lx,Xh%ly,Xh%lz)
    integer :: e, k, lxy, lxz, lyz, lxyz

    lxy = Xh%lx*Xh%ly
    lxz = Xh%lx*Xh%lz
    lyz = Xh%ly*Xh%lz
    lxyz = Xh%lx*Xh%ly*Xh%lz

  
    do e = 1, msh%nelv
       if(msh%gdim .eq. 2) then
         call mxm  (Xh%dx,Xh%lx,u(1,1,1,e),Xh%lx,dudr,lyz)
         call mxm  (u(1,1,1,e),Xh%lx,Xh%dyt,Xh%ly,duds,Xh%ly)
         call col3 (tmp1,dudr,coef%G1(1,1,1,e),lxyz)
         call col3 (tmp2,duds,coef%G2(1,1,1,e),lxyz)
         if (msh%dfrmd_el(e)) then
            call addcol3 (tmp1,duds,coef%G4(1,1,1,e),lxyz)
            call addcol3 (tmp2,dudr,coef%G4(1,1,1,e),lxyz)
         end if
         call col2 (tmp1,coef%h1(1,1,1,e),lxyz)
         call col2 (tmp2,coef%h1(1,1,1,e),lxyz)
         call mxm  (Xh%dxt,Xh%lx,tmp1,Xh%lx,tm1,lyz)
         call mxm  (tmp2,Xh%lx,Xh%dy,Xh%ly,tm2,Xh%ly)
         call add3 (w(1,1,1,e),tm1,tm2,lxyz)

       ! 3D evaluation!
       else
         call mxm(Xh%dx,Xh%lx,u(1,1,1,e),Xh%lx,dudr,lyz)
         do k=1,Xh%lz
            call mxm(u(1,1,k,e),Xh%lx,Xh%dyt,Xh%ly,duds(1,1,k),Xh%ly)
         end do
         call mxm     (u(1,1,1,e),lxy,Xh%dzt,Xh%lz,dudt,Xh%lz)
         call col3    (tmp1,dudr,coef%G1(1,1,1,e),lxyz)
         call col3    (tmp2,duds,coef%G2(1,1,1,e),lxyz)
         call col3    (tmp3,dudt,coef%G3(1,1,1,e),lxyz)
         if (msh%dfrmd_el(e)) then
            call addcol3 (tmp1,duds,coef%G4(1,1,1,e),lxyz)
            call addcol3 (tmp1,dudt,coef%G5(1,1,1,e),lxyz)
            call addcol3 (tmp2,dudr,coef%G4(1,1,1,e),lxyz)
            call addcol3 (tmp2,dudt,coef%G6(1,1,1,e),lxyz)
            call addcol3 (tmp3,dudr,coef%G5(1,1,1,e),lxyz)
            call addcol3 (tmp3,duds,coef%G6(1,1,1,e),lxyz)
         end if
         call col2 (tmp1,coef%h1(1,1,1,e),lxyz)
         call col2 (tmp2,coef%h1(1,1,1,e),lxyz)
         call col2 (tmp3,coef%h1(1,1,1,e),lxyz)
         call mxm  (Xh%dxt,Xh%lx,tmp1,Xh%lx,tm1,lyz)
         do k=1,Xh%lz
            call mxm(tmp2(1,1,k),Xh%lx,Xh%dy,Xh%ly,tm2(1,1,k),Xh%ly)
         end do
         call mxm  (tmp3,lxy,Xh%dz,Xh%lz,tm3,Xh%lz)
         call add4(w(1,1,1,e),tm1,tm2,tm3,lxyz)
       end if 
    end do
    
    if (coef%ifh2) call addcol4 (w,coef%h2,coef%B,u,coef%dof%n_dofs)

  end subroutine ax_helm_compute
  
end module ax_helm
