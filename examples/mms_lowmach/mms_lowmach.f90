!> Manufactured-solution verification for the low-Mach viscous stress in
!> BOTH residuals (lowmach_vel_res_cpu_compute and lowmach_prs_res_cpu_compute),
!> with a spatially varying viscosity.
!>
!> Deviatoric stress:  tau = mu (grad u + grad u^T) - (2/3) mu (div u) I.
!>
!> Manufactured fields on the fully periodic box [0,2pi]^3:
!>   u   = (sin x, sin y, sin z)            (so div u = cos x + cos y + cos z)
!>   Q_T = cos x + cos y + cos z            (= div u, exactly, by construction)
!>   mu  = mu0 ( 1 + a sin x sin y sin z )  (variable viscosity)
!>
!> Because grad u is diagonal, the stress divergence collapses to the clean form
!>   (div tau)_i = -(4/3) mu sin(x_i) + d_i(mu) ( 2 cos(x_i) - (2/3) Q_T ),
!> which reduces to -(4/3) mu0 sin(x_i) for constant mu. This is the analytic
!> target `au,av,aw`.
!>
!> Test A (momentum): the shipped velocity residual with p=0, f=0, bd=0 returns
!>   the weak div(tau); gather-scatter + inverse mass -> nodal div(tau), compared
!>   to (au,av,aw). Exercises ax_helm_full (incl. the variable-mu S^T grad mu via
!>   h1=mu) and the explicit dilatation -(2/3) grad(mu Q_T).
!>
!> Test B (pressure): the shipped pressure residual with p=0, f=0, bd=0, rho=1
!>   returns the weak div(div tau). Comparing it to the discrete divergence of the
!>   analytic stress vector (cdtp on au,av,aw, same post-processing) verifies the
!>   rotational viscous form plus the new low-Mach corrections
!>   -(4/3) mu grad(Q_T) + (2/3) Q_T grad(mu) and the S^T grad mu term.
!>
!> A wrong sign/coefficient in either residual leaves an O(1) error that does not
!> converge under p-refinement.
program mms_lowmach
  use num_types, only : rp
  use mesh, only : mesh_t
  use file, only : file_t
  use space, only : space_t, GLL
  use dofmap, only : dofmap_t
  use gather_scatter, only : gs_t, GS_OP_ADD
  use coefs, only : coef_t
  use field, only : field_t
  use facet_normal, only : facet_normal_t
  use math, only : glsc3
  use operators, only : cdtp
  use scratch_registry, only : neko_scratch_registry
  use ax_product, only : ax_t, ax_helm_factory
  use lowmach_res_cpu, only : lowmach_vel_res_cpu_t, lowmach_prs_res_cpu_t
  use neko, only : neko_init, neko_finalize
  use, intrinsic :: iso_c_binding, only : c_ptr, C_NULL_PTR
  implicit none

  type(mesh_t) :: msh
  type(file_t) :: msh_file
  integer :: lxs(4), il

  call neko_init()

  call msh_file%init("box.nmsh")
  call msh_file%read(msh)

  write(*,'(A)') ''
  write(*,'(A)') '  Low-Mach viscous-stress MMS (variable mu)'
  write(*,'(A)') '   lx    momentum: Linf      relL2        pressure: Linf      relL2'
  write(*,'(A)') '  --------------------------------------------------------------------'

  lxs = [4, 6, 8, 10]
  do il = 1, size(lxs)
     call run_resolution(msh, lxs(il))
  end do
  write(*,'(A)') ''

  call neko_finalize()

contains

  subroutine run_resolution(msh, lx)
    type(mesh_t), intent(inout) :: msh
    integer, intent(in) :: lx
    type(space_t) :: Xh
    type(dofmap_t), target :: dof
    type(gs_t) :: gs_h
    type(coef_t) :: coef
    class(ax_t), allocatable :: ax_vel, ax_prs
    type(lowmach_vel_res_cpu_t) :: vres
    type(lowmach_prs_res_cpu_t) :: pres
    type(facet_normal_t) :: bc_prs, bc_sym
    type(c_ptr) :: event
    type(field_t) :: u, v, w, p, fx, fy, fz, mu, rho, qt
    type(field_t) :: ur, vr, wr, au, av, aw, pr, t1, t2, t3, la
    real(kind=rp), parameter :: mu0 = 0.7_rp
    real(kind=rp), parameter :: a = 0.5_rp
    real(kind=rp), parameter :: c43 = 4.0_rp / 3.0_rp
    real(kind=rp), parameter :: c23 = 2.0_rp / 3.0_rp
    real(kind=rp) :: xx, yy, zz, sx, sy, sz, cx, cy, cz, m, qv
    real(kind=rp) :: dmx, dmy, dmz, lv, l2e, l2a, lp
    integer :: i, n

    event = C_NULL_PTR

    ! SEM machinery for this polynomial order.
    call Xh%init(GLL, lx, lx, lx)
    call dof%init(msh, Xh)
    call gs_h%init(dof)
    call coef%init(gs_h)
    n = dof%size()

    ! Point the global scratch registry at this dofmap.
    call neko_scratch_registry%free()
    call neko_scratch_registry%init(16, 16, dof)

    call u%init(dof, 'u');   call v%init(dof, 'v');   call w%init(dof, 'w')
    call p%init(dof, 'p')
    call fx%init(dof, 'fx'); call fy%init(dof, 'fy'); call fz%init(dof, 'fz')
    call mu%init(dof, 'mu'); call rho%init(dof, 'rho'); call qt%init(dof, 'qt')
    call ur%init(dof, 'ur'); call vr%init(dof, 'vr'); call wr%init(dof, 'wr')
    call au%init(dof, 'au'); call av%init(dof, 'av'); call aw%init(dof, 'aw')
    call pr%init(dof, 'pr'); call la%init(dof, 'la')
    call t1%init(dof, 't1'); call t2%init(dof, 't2'); call t3%init(dof, 't3')

    call ax_helm_factory(ax_vel, full_formulation = .true.)
    call ax_helm_factory(ax_prs, full_formulation = .false.)

    ! Manufactured fields and analytic div(tau).
    do i = 1, n
       xx = dof%x(i,1,1,1); yy = dof%y(i,1,1,1); zz = dof%z(i,1,1,1)
       sx = sin(xx); sy = sin(yy); sz = sin(zz)
       cx = cos(xx); cy = cos(yy); cz = cos(zz)
       m  = mu0 * (1.0_rp + a * sx * sy * sz)
       qv = cx + cy + cz
       dmx = mu0 * a * cx * sy * sz
       dmy = mu0 * a * sx * cy * sz
       dmz = mu0 * a * sx * sy * cz

       u%x(i,1,1,1) = sx; v%x(i,1,1,1) = sy; w%x(i,1,1,1) = sz
       mu%x(i,1,1,1) = m
       rho%x(i,1,1,1) = 1.0_rp
       qt%x(i,1,1,1) = qv
       p%x(i,1,1,1) = 0.0_rp
       fx%x(i,1,1,1) = 0.0_rp; fy%x(i,1,1,1) = 0.0_rp; fz%x(i,1,1,1) = 0.0_rp

       au%x(i,1,1,1) = -c43 * m * sx + dmx * (2.0_rp * cx - c23 * qv)
       av%x(i,1,1,1) = -c43 * m * sy + dmy * (2.0_rp * cy - c23 * qv)
       aw%x(i,1,1,1) = -c43 * m * sz + dmz * (2.0_rp * cz - c23 * qv)
    end do

    ! ---- Test A: momentum residual returns weak div(tau) ----
    call vres%compute(ax_vel, u, v, w, ur, vr, wr, p, fx, fy, fz, &
         coef, msh, Xh, mu, rho, qt, 0.0_rp, 1.0_rp, n)
    call gs_h%op(ur, GS_OP_ADD)
    call gs_h%op(vr, GS_OP_ADD)
    call gs_h%op(wr, GS_OP_ADD)
    lv = 0.0_rp
    do i = 1, n
       ur%x(i,1,1,1) = ur%x(i,1,1,1) * coef%Binv(i,1,1,1) - au%x(i,1,1,1)
       vr%x(i,1,1,1) = vr%x(i,1,1,1) * coef%Binv(i,1,1,1) - av%x(i,1,1,1)
       wr%x(i,1,1,1) = wr%x(i,1,1,1) * coef%Binv(i,1,1,1) - aw%x(i,1,1,1)
       lv = max(lv, abs(ur%x(i,1,1,1)), abs(vr%x(i,1,1,1)), abs(wr%x(i,1,1,1)))
    end do
    l2e = glsc3(ur%x, ur%x, coef%B, n) + glsc3(vr%x, vr%x, coef%B, n) &
        + glsc3(wr%x, wr%x, coef%B, n)
    l2a = glsc3(au%x, au%x, coef%B, n) + glsc3(av%x, av%x, coef%B, n) &
        + glsc3(aw%x, aw%x, coef%B, n)
    l2e = sqrt(l2e / l2a)

    ! ---- Test B: pressure residual returns weak div(div tau) ----
    call bc_prs%init_from_components(coef); call bc_prs%finalize()
    call bc_sym%init_from_components(coef); call bc_sym%finalize()

    ! Shipped pressure residual: p=0, f=0, bd=0, rho=1 -> p_res = weak div(div tau).
    call pres%compute(p, pr, u, v, w, u, v, w, fx, fy, fz, coef, gs_h, &
         bc_prs, bc_sym, ax_prs, 0.0_rp, 1.0_rp, mu, rho, qt, event)

    ! Discrete divergence of the analytic stress vector, same post-processing.
    call cdtp(t1%x, au%x, coef%drdx, coef%dsdx, coef%dtdx, coef)
    call cdtp(t2%x, av%x, coef%drdy, coef%dsdy, coef%dtdy, coef)
    call cdtp(t3%x, aw%x, coef%drdz, coef%dsdz, coef%dtdz, coef)
    do i = 1, n
       la%x(i,1,1,1) = t1%x(i,1,1,1) + t2%x(i,1,1,1) + t3%x(i,1,1,1)
    end do
    call gs_h%op(pr, GS_OP_ADD)
    call gs_h%op(la, GS_OP_ADD)
    lp = 0.0_rp
    do i = 1, n
       pr%x(i,1,1,1) = (pr%x(i,1,1,1) - la%x(i,1,1,1)) * coef%Binv(i,1,1,1)
       la%x(i,1,1,1) = la%x(i,1,1,1) * coef%Binv(i,1,1,1)
       lp = max(lp, abs(pr%x(i,1,1,1)))
    end do
    l2a = glsc3(la%x, la%x, coef%B, n)
    ! relative L2 for pressure
    block
      real(kind=rp) :: pe
      pe = glsc3(pr%x, pr%x, coef%B, n)
      write(*,'(I5, 4ES19.9)') lx, lv, l2e, lp, sqrt(pe / l2a)
    end block

    deallocate(ax_vel)
    deallocate(ax_prs)
  end subroutine run_resolution

end program mms_lowmach
