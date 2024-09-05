!>
!> USE:
!>
!>  type(tamg_hierarchy_t) :: amg
!>  type(tamg_solver_t) :: amg_solver
!>
!>  call amg%init(ax, Xh, coef, msh, gs_h, 3)
!>  call amg_solver%init(amg, niter)
!>
!>  call amg_solver%solve(x%x, f, n)
!>
module tree_amg_multigrid
  use math
  use utils
  use tree_amg
  use tree_amg_matvec
  use tree_amg_smoother
  implicit none
  private

  type, public :: tamg_solver_t
    type(tamg_hierarchy_t), pointer :: amg
    type(amg_cheby_t), allocatable :: smoo(:)
    integer :: nlvls
    integer :: max_iter
  contains
    procedure, pass(this) :: init => tamg_mg_init
    procedure, pass(this) :: solve => tamg_mg_solve
  end type tamg_solver_t

contains

  subroutine tamg_mg_init(this, amg, max_iter)
    class(tamg_solver_t), intent(inout), target :: this
    class(tamg_hierarchy_t), intent(inout), target :: amg
    integer, intent(in) :: max_iter
    integer :: lvl, nlvls, n

    this%amg => amg

    this%max_iter = max_iter

    nlvls = amg%nlvls
    this%nlvls = nlvls
    if (this%nlvls .gt. amg%nlvls) then
      call neko_error("Requested number multigrid levels is greater than the initialized AMG levels")
    end if
    !allocate(this%smoo(0:(nlvls-1)))
    allocate(this%smoo(0:(nlvls)))
    do lvl = 0, nlvls-1
      n = amg%lvl(lvl+1)%fine_lvl_dofs
      call this%smoo(lvl)%init(n ,lvl, 10)
    end do

  end subroutine tamg_mg_init

  subroutine tamg_mg_solve(this, z, r, n)
    integer, intent(in) :: n
    class(tamg_solver_t), intent(inout) :: this
    real(kind=rp), dimension(n), intent(inout) :: z
    real(kind=rp), dimension(n), intent(inout) :: r
    integer :: iter, max_iter

    max_iter = this%max_iter

    z = 0d0

    do iter = 1, max_iter
      !print *, "MG iter:", iter
      call tamg_mg_cycle(z, r, n, 0, this%amg, this)
    end do
  end subroutine tamg_mg_solve

  recursive subroutine tamg_mg_cycle(x, b, n, lvl, amg, mgstuff)
    integer, intent(in) :: n
    real(kind=rp), intent(inout) :: x(n)
    real(kind=rp), intent(inout) :: b(n)
    type(tamg_hierarchy_t), intent(inout) :: amg
    type(tamg_solver_t), intent(inout) :: mgstuff
    integer, intent(in) :: lvl
    real(kind=rp) :: r(n)
    real(kind=rp) :: rc(n)
    real(kind=rp) :: tmp(n)
    integer :: iter, num_iter
    integer :: sit, max_lvl
    integer :: i, cyt

    r = 0d0
    rc = 0d0
    tmp = 0d0

    sit = 10
    max_lvl = mgstuff%nlvls-1

    !>----------<!
    !> SMOOTH   <!
    !>----------<!
    call mgstuff%smoo(lvl)%solve(x,b, n, amg, sit)
    if (lvl .eq. max_lvl) then !> Is coarsest grid.
      return
    end if

    !>----------<!
    !> Residual <!
    !>----------<!
    call calc_resid(r,x,b,amg,lvl,n)

    !>----------<!
    !> Restrict <!
    !>----------<!
    if (lvl .eq. 0) then
      call average_duplicates(r,amg,lvl,n)
    end if
    call amg%interp_f2c(rc, r, lvl+1)

    !>-------------------<!
    !> Call Coarse solve <!
    !>-------------------<!
    tmp = 0d0
    call tamg_mg_cycle(tmp, rc, amg%lvl(lvl+1)%nnodes, lvl+1, amg, mgstuff)

    !>----------<!
    !> Project  <!
    !>----------<!
    call amg%interp_c2f(r, tmp, lvl+1)
    if (lvl .eq. 0) then
      call average_duplicates(r,amg,lvl,n)
    end if

    !>----------<!
    !> Correct  <!
    !>----------<!
    call add2(x, r, n)

    !>----------<!
    !> SMOOTH   <!
    !>----------<!
    call mgstuff%smoo(lvl)%solve(x,b, n, amg, sit)

    !>----------<!
    !> Residual <!
    !>----------<!
    !call calc_resid(r,x,b,amg,lvl,n)!> TODO: for debug
    !print *, "LVL:",lvl, "POST RESID:", sqrt(glsc2(r, r, n))
  end subroutine tamg_mg_cycle

  subroutine my_matvec(vec_out, vec_in, n, lvl, amg)
    integer, intent(in) :: n
    real(kind=rp), intent(inout) :: vec_out(n)
    real(kind=rp), intent(inout) :: vec_in(n)
    integer, intent(in) :: lvl
    type(tamg_hierarchy_t), intent(inout) :: amg
    integer :: lvl_start

    lvl_start = amg%nlvls

    call amg%matvec(vec_out, vec_in, lvl)
  end subroutine my_matvec

  subroutine calc_resid(r, x, b, amg, lvl, n)
    integer, intent(in) :: n
    real(kind=rp), intent(inout) :: r(n)
    real(kind=rp), intent(inout) :: x(n)
    real(kind=rp), intent(inout) :: b(n)
    type(tamg_hierarchy_t), intent(inout) :: amg
    integer, intent(in) :: lvl
    integer :: i
    r = 0d0
    call my_matvec(r, x, n, lvl, amg)
    do  i = 1, n
      r(i) = b(i) - r(i)
    end do
  end subroutine calc_resid

  subroutine average_duplicates(U, amg, lvl, n)
    integer, intent(in) :: n
    real(kind=rp), intent(inout) :: U(n)
    type(tamg_hierarchy_t), intent(inout) :: amg
    integer, intent(in) :: lvl
    integer :: i
    call amg%gs_h%op(U, n, GS_OP_ADD)
    do  i = 1, n
      U(i) = U(i) * amg%coef%mult(i,1,1,1)
    end do
  end subroutine average_duplicates
end module tree_amg_multigrid
