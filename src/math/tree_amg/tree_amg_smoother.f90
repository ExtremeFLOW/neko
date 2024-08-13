module tree_amg_smoother
  use tree_amg
  use num_types
  use utils
  use math
  use krylov, only : ksp_monitor_t
  implicit none
  private

  type, public :: amg_cheby_t
     real(kind=rp), allocatable :: d(:)
     real(kind=rp), allocatable :: w(:)
     real(kind=rp), allocatable :: r(:)
     real(kind=rp) :: tha, dlt
     integer :: lvl
     integer :: n
     integer :: power_its = 15
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
    real(kind=rp) :: boost = 1.2_rp
    real(kind=rp) :: lam_factor = 30.0_rp
    real(kind=rp) :: wtw, dtw, dtd
    integer :: i
    associate(w => this%w, d => this%d, coef => amg%coef, gs_h => amg%gs_h, msh=>amg%msh, Xh=>amg%Xh)

      print *, "COMP EIGS on lvl", this%lvl, "n", n
      do i = 1, n
        !TODO: replace with a better way to initialize power method
        call random_number(rn)
        d(i) = rn + 10.0_rp
      end do
      !!if (this%lvl .eq. 1) then
      !!  call gs_h%op(d, n, GS_OP_ADD)TODO
      !!  !--call bc_list_apply(blst, d, n)
      !!end if

      !Power method to get lamba max
      do i = 1, this%power_its
        w = 0d0
        call amg%matvec(w, d, this%lvl)

        if (this%lvl .eq. 1) then
          wtw = glsc3(w, coef%mult, w, n)
        else
          wtw = glsc2(w, w, n)
        end if
        call cmult2(d, w, 1.0_rp/sqrt(wtw), n)
      end do

      w = 0d0
      call amg%matvec(w, d, this%lvl)

      if (this%lvl .eq. 1) then
        dtw = glsc3(d, coef%mult, w, n)
        dtd = glsc3(d, coef%mult, d, n)
      else
        dtw = glsc2(d, w, n)
        dtd = glsc2(d, d, n)
      end if
      lam = dtw / dtd
      b = lam * boost
      a = lam / lam_factor
      this%tha = (b+a)/2.0_rp
      this%dlt = (b-a)/2.0_rp

      this%recompute_eigs = .false.
    end associate
  end subroutine amg_cheby_power

  !> A chebyshev preconditioner
  subroutine amg_cheby_solve(this, x, f, n, amg, niter)
    class(amg_cheby_t), intent(inout) :: this
    integer, intent(in) :: n
    real(kind=rp), dimension(n), intent(inout) :: x
    real(kind=rp), dimension(n), intent(inout) :: f
    class(tamg_hierarchy_t), intent(inout) :: amg
    type(ksp_monitor_t) :: ksp_results
    integer, optional, intent(in) :: niter
    integer :: iter, max_iter
    real(kind=rp) :: a, b, rtr, rnorm

    if (this%recompute_eigs) then
       call this%comp_eig(amg, n)
    end if

    if (present(niter)) then
       max_iter = niter
    else
       max_iter = this%max_iter
    end if

    associate( w => this%w, r => this%r, d => this%d)
      ! calculate residual
      call copy(r, f, n)
      w = 0d0
      call amg%matvec(w, x, this%lvl)
      call sub2(r, w, n)

      rtr = glsc2(r, r, n)
      rnorm = sqrt(rtr)
      ksp_results%res_start = rnorm
      ksp_results%res_final = rnorm
      ksp_results%iter = 0

      ! First iteration
      call copy(w, r, n) ! PRECOND
      call copy(d, w, n)
      a = 2.0_rp / this%tha
      call add2s2(x, d, a, n)! x = x + a*d

      ! Rest of the iterations
      do iter = 2, max_iter
        ! calculate residual
        call copy(r, f, n)
        w = 0d0
        call amg%matvec(w, x, this%lvl)
        call sub2(r, w, n)

        call copy(w, r, n)! PRECOND

        b = (this%dlt * a / 2.0_rp)**2
        a = 1.0_rp / (this%tha - b)
        call add2s1(d, w, b, n)! d = w + b*d

        call add2s2(x, d, a, n)! x = x + a*d
      end do

      ! calculate residual
      call copy(r, f, n)
      w = 0d0
      call amg%matvec(w, x, this%lvl)
      call sub2(r, w, n)
      rtr = glsc2(r, r, n)
      rnorm = sqrt(rtr)
      ksp_results%res_final = rnorm
      ksp_results%iter = iter
    end associate
  end subroutine amg_cheby_solve

end module tree_amg_smoother
