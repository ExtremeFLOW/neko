module ax_helm_xsmm
  use ax_product
  use num_types
#ifdef HAVE_LIBXSMM
  use libxsmm, libxsmm_mmcall => libxsmm_dmmcall_abc
#endif
  implicit none
  private
  
  type, public, extends(ax_t) :: ax_helm_xsmm_t
   contains
     procedure, nopass :: compute => ax_helm_xsmm_compute
  end type ax_helm_xsmm_t

contains 
  
  subroutine ax_helm_xsmm_compute(w, u, coef, msh, Xh)
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
#ifdef HAVE_LIBXSMM
    type(libxsmm_dmmfunction), save :: ax_helm_xmm1
    type(libxsmm_dmmfunction), save :: ax_helm_xmm2
    type(libxsmm_dmmfunction), save :: ax_helm_xmm3
    integer, save :: ax_helm_xsmm_lx = 0
    logical, save :: ax_helm_xsmm_init = .false.

    lxy = Xh%lx*Xh%ly
    lxz = Xh%lx*Xh%lz
    lyz = Xh%ly*Xh%lz
    lxyz = Xh%lx*Xh%ly*Xh%lz

    if (.not. ax_helm_xsmm_init .or. (ax_helm_xsmm_lx .ne. Xh%lx)) then
       call libxsmm_dispatch(ax_helm_xmm1, Xh%lx, Xh%ly*Xh%lz, Xh%lx, &
            alpha=1d0, beta=0d0, prefetch=LIBXSMM_PREFETCH_AUTO)
       call libxsmm_dispatch(ax_helm_xmm2, Xh%lx, Xh%ly, Xh%ly, &
            alpha=1d0, beta=0d0, prefetch=LIBXSMM_PREFETCH_AUTO)
       call libxsmm_dispatch(ax_helm_xmm3, Xh%lx*Xh%ly, Xh%lz, Xh%lz, &
            alpha=1d0, beta=0d0, prefetch=LIBXSMM_PREFETCH_AUTO)
       ax_helm_xsmm_init = .true.
       ax_helm_xsmm_lx = Xh%lx
    end if

  
    do e = 1, msh%nelv
       if(msh%gdim .eq. 2) then
         call mxm(Xh%dx, Xh%lx,u(1,1,1,e), Xh%lx, dudr, lyz)
         call mxm(u(1,1,1,e), Xh%lx, Xh%dyt, Xh%ly, duds, Xh%ly)
         call col3(tmp1, dudr, coef%G11(1,1,1,e), lxyz)
         call col3(tmp2, duds, coef%G22(1,1,1,e), lxyz)
         if (msh%dfrmd_el(e)) then
            call addcol3(tmp1, duds, coef%G12(1,1,1,e), lxyz)
            call addcol3(tmp2, dudr, coef%G12(1,1,1,e), lxyz)
         end if
         call col2(tmp1, coef%h1(1,1,1,e), lxyz)
         call col2(tmp2, coef%h1(1,1,1,e), lxyz)
         call mxm(Xh%dxt, Xh%lx, tmp1, Xh%lx, tm1, lyz)
         call mxm(tmp2, Xh%lx, Xh%dy, Xh%ly, tm2, Xh%ly)
         call add3(w(1,1,1,e), tm1, tm2, lxyz)

       ! 3D evaluation!
      else
         call libxsmm_mmcall(ax_helm_xmm1, Xh%dx, u(1,1,1,e), dudr)
         do k = 1,Xh%lz
            call libxsmm_mmcall(ax_helm_xmm2, u(1,1,k,e), Xh%dyt, duds(1,1,k))
         end do
         call libxsmm_mmcall(ax_helm_xmm3, u(1,1,1,e), Xh%dzt, dudt)
         call col3(tmp1, dudr, coef%G11(1,1,1,e), lxyz)
         call col3(tmp2, duds, coef%G22(1,1,1,e), lxyz)
         call col3(tmp3, dudt, coef%G33(1,1,1,e), lxyz)
         if (msh%dfrmd_el(e)) then
            call addcol3(tmp1, duds, coef%G12(1,1,1,e), lxyz)
            call addcol3(tmp1, dudt, coef%G13(1,1,1,e), lxyz)
            call addcol3(tmp2, dudr, coef%G12(1,1,1,e), lxyz)
            call addcol3(tmp2, dudt, coef%G23(1,1,1,e), lxyz)
            call addcol3(tmp3, dudr, coef%G13(1,1,1,e), lxyz)
            call addcol3(tmp3, duds, coef%G23(1,1,1,e), lxyz)
         end if
         call col2(tmp1, coef%h1(1,1,1,e), lxyz)
         call col2(tmp2, coef%h1(1,1,1,e), lxyz)
         call col2(tmp3, coef%h1(1,1,1,e), lxyz)
         call libxsmm_mmcall(ax_helm_xmm1, Xh%dxt, tmp1, tm1)
         do k = 1,Xh%lz
            call libxsmm_mmcall(ax_helm_xmm2, tmp2(1,1,k), Xh%dy, tm2(1,1,k))
         end do
         call libxsmm_mmcall(ax_helm_xmm3, tmp3, Xh%dz, tm3)
         call add4(w(1,1,1,e), tm1, tm2, tm3, lxyz)
       end if 
    end do
    
    if (coef%ifh2) call addcol4 (w,coef%h2,coef%B,u,coef%dof%n_dofs)
#endif

  end subroutine ax_helm_xsmm_compute
  
end module ax_helm_xsmm
