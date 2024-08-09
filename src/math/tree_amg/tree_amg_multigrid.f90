module tree_amg_multigrid
  use math
  use tree_amg_matvec
  implicit none

contains

  recursive subroutine tamg_mg_cycle(x, b, n, lvl, amg)
    integer, intent(in) :: n
    real(kind=rp), intent(inout) :: x(n)
    real(kind=rp), intent(inout) :: b(n)
    type(tamg_hierarchy_t), intent(inout) :: amg
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

    sit = 4
    max_lvl = 1

    !>----------<!
    !> Residual <! NOT NEEDED
    !>----------<!
    call calc_resid(r,x,b,amg,lvl,n)
    if (lvl .eq. 1) then
      call average_duplicates(r,amg,lvl,n)
    end if
    print *, "LVL:",lvl, "INIT RESID:", sqrt(glsc2(r, r, n))

    !>----------<!
    !> SMOOTH   <!
    !>----------<!
    call cheby(x,b,amg,lvl,n,sit)

    !>----------<!
    !> Residual <!
    !>----------<!
    call calc_resid(r,x,b,amg,lvl,n)
    print *, "LVL:",lvl, "PRE RESID:", sqrt(glsc2(r, r, n))

    if (lvl .eq. max_lvl) then!TODO: move when done debugging
      return
    end if

    !>----------<!
    !> Restrict <!
    !>----------<!
    if (lvl .eq. 1) then
      call average_duplicates(r,amg,lvl,n)
    end if
    call restriction_operator(rc, r, lvl+1, amg)

    !>-------------------<!
    !> Call Coarse solve <!
    !>-------------------<!
    tmp = 0d0
    call tamg_mg_cycle(tmp, rc, amg, lvl+1, n)

    !>----------<!
    !> Project  <!
    !>----------<!
    call prolongation_operator(r, tmp, lvl+1, amg)
    if (lvl .eq. 1) then
      call average_duplicates(r,amg,lvl,n)
    end if

    !>----------<!
    !> Correct  <!
    !>----------<!
    call add2(x, r, n)

    !>----------<!
    !> SMOOTH   <!
    !>----------<!
    call cheby(x,b,amg,lvl,n,sit)

    !>----------<!
    !> Residual <!
    !>----------<!
    call calc_resid(r,x,b,amg,lvl,n)
    print *, "LVL:",lvl, "POST RESID:", sqrt(glsc2(r, r, n))
  end subroutine tamg_mg_cycle

  subroutine my_matvec(vec_out, vec_in, n, lvl, amg)
    integer, intent(in) :: n
    real(kind=rp), intent(inout) :: vec_out(n)
    real(kind=rp), intent(inout) :: vec_in(n)
    integer, intent(in) :: lvl
    type(tamg_hierarchy_t), intent(inout) :: amg

    call matvec_operator(vec_out, vec_in, 1, lvl, amg)
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
    call my_matvec(r, x, lvl, amg, n)
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
    call 
    do  i = 1, n
      U(i) = U(i) * amg%coef%mult(i,1,1,1)
    end do
  end subroutine average_duplicates
end module tree_amg_multigrid
