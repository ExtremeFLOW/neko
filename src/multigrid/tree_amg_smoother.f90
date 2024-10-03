module tree_amg_smoother
  use tree_amg
  use num_types
  use utils
  use math
  use krylov, only : ksp_monitor_t
  use bc, only: bc_list_t, bc_list_apply
  use gather_scatter, only : gs_t, GS_OP_ADD
  implicit none
  private

  type, public :: amg_cheby_t
     real(kind=rp), allocatable :: d(:)
     real(kind=rp), allocatable :: w(:)
     real(kind=rp), allocatable :: r(:)
     real(kind=rp) :: tha, dlt
     integer :: lvl
     integer :: n
     integer :: power_its = 250
     integer :: max_iter = 10
     logical :: recompute_eigs = .true.
   contains
     procedure, pass(this) :: init => amg_cheby_init
     procedure, pass(this) :: solve => amg_cheby_solve
     procedure, pass(this) :: comp_eig => amg_cheby_power
  end type amg_cheby_t

contains

  !> Initialise a standard solver
  subroutine amg_cheby_init(this, n, lvl, max_iter)
    class(amg_cheby_t), intent(inout), target :: this
    integer, intent(in) :: n
    integer, intent(in) :: lvl
    integer, intent(in) :: max_iter

    allocate(this%d(n))
    allocate(this%w(n))
    allocate(this%r(n))
    this%n = n
    this%lvl = lvl
    this%max_iter = max_iter
    this%recompute_eigs = .true.
    print *, "INIT SMOO ON LVL", lvl

  end subroutine amg_cheby_init

  subroutine amg_cheby_power(this, amg, n)
    class(amg_cheby_t), intent(inout) :: this
    type(tamg_hierarchy_t), intent(inout) :: amg
    integer, intent(in) :: n
    real(kind=rp) :: lam, b, a, rn
    !real(kind=rp) :: boost = 1.2_rp
    real(kind=rp) :: boost = 1.1_rp
    real(kind=rp) :: lam_factor = 30.0_rp
    real(kind=rp) :: wtw, dtw, dtd
    integer :: i
    associate(w => this%w, d => this%d, coef => amg%coef, gs_h => amg%gs_h, msh=>amg%msh, Xh=>amg%Xh, blst=>amg%blst)

      print *, "COMP EIGS on lvl", this%lvl, "n", n
      do i = 1, n
        !TODO: replace with a better way to initialize power method
        !call random_number(rn)
        !d(i) = rn + 10.0_rp
        d(i) = sin(real(i))
      end do
      if (this%lvl .eq. 0) then
        call gs_h%op(d, n, GS_OP_ADD)!TODO
        call bc_list_apply(blst, d, n)
      end if

      !Power method to get lamba max
      do i = 1, this%power_its
        w = 0d0
        call amg%matvec(w, d, this%lvl)

        if (this%lvl .eq. 0) then
          !call bc_list_apply(blst, w, n)
          wtw = glsc3(w, coef%mult, w, n)
        else
          wtw = glsc2(w, w, n)
        end if

        call cmult2(d, w, 1.0_rp/sqrt(wtw), n)

        if (this%lvl .eq. 0) then
          !call bc_list_apply(blst, d, n)
        end if
      end do

      w = 0d0
      call amg%matvec(w, d, this%lvl)
      if (this%lvl .eq. 0) then
        !call bc_list_apply(blst, w, n)
      end if

      if (this%lvl .eq. 0) then
        dtw = glsc3(d, coef%mult, w, n)
        dtd = glsc3(d, coef%mult, d, n)
      else
        dtw = glsc2(d, w, n)
        dtd = glsc2(d, d, n)
      end if
      lam = dtw / dtd
      !print *, "LAM:", lam
      b = lam * boost
      a = lam / lam_factor
      this%tha = (b+a)/2.0_rp
      this%dlt = (b-a)/2.0_rp

      this%recompute_eigs = .false.
    end associate
  end subroutine amg_cheby_power

  !> A chebyshev preconditioner
  !> From Saad's iterative methods textbook
  subroutine amg_cheby_solve(this, x, f, n, amg, niter)
    class(amg_cheby_t), intent(inout) :: this
    integer, intent(in) :: n
    real(kind=rp), dimension(n), intent(inout) :: x
    real(kind=rp), dimension(n), intent(inout) :: f
    class(tamg_hierarchy_t), intent(inout) :: amg
    type(ksp_monitor_t) :: ksp_results
    integer, optional, intent(in) :: niter
    integer :: iter, max_iter
    real(kind=rp) :: rtr, rnorm
    real(kind=rp) :: rhok, rhokp1, s1, thet, delt, tmp1, tmp2

    if (this%recompute_eigs) then
       call this%comp_eig(amg, n)
    end if

    if (present(niter)) then
       max_iter = niter
    else
       max_iter = this%max_iter
    end if

    associate( w => this%w, r => this%r, d => this%d, blst=>amg%blst)
      call copy(r, f, n)
      w = 0d0
      call amg%matvec(w, x, this%lvl)
      call sub2(r, w, n)

      thet = this%tha
      delt = this%dlt
      s1 = thet / delt
      rhok = 1.0_rp / s1

      call copy(d, r, n)
      call cmult(d, 1.0_rp/thet, n)

      do iter = 1, max_iter
        call add2(x,d,n)

        w = 0d0
        call amg%matvec(w, d, this%lvl)
        call sub2(r, w, n)

        rhokp1 = 1.0_rp / (2.0_rp * s1 - rhok)
        tmp1 = rhokp1 * rhok
        tmp2 = 2.0_rp * rhokp1 / delt
        rhok = rhokp1

        call cmult(d, tmp1, n)
        call add2s2(d, r, tmp2, n)
      end do

    end associate
  end subroutine amg_cheby_solve

end module tree_amg_smoother
