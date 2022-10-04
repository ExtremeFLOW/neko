module rhs_maker_sx
  use rhs_maker
  implicit none
  private
  
  type, public, extends(rhs_maker_sumab_t) :: rhs_maker_sumab_sx_t
   contains
     procedure, nopass :: compute_fluid => rhs_maker_sumab_sx
  end type rhs_maker_sumab_sx_t

  type, public, extends(rhs_maker_ext_t) ::  rhs_maker_ext_sx_t
   contains
     procedure, nopass :: compute_fluid => rhs_maker_ext_sx
     procedure, nopass :: compute_scalar => scalar_rhs_maker_ext_sx
  end type rhs_maker_ext_sx_t

  type, public, extends(rhs_maker_bdf_t) :: rhs_maker_bdf_sx_t
   contains
     procedure, nopass :: compute_fluid => rhs_maker_bdf_sx
  end type rhs_maker_bdf_sx_t
  
contains

  subroutine rhs_maker_sumab_sx(u, v, w, uu, vv, ww, uulag, vvlag, wwlag, ab, nab)
    type(field_t), intent(inout) :: u,v, w
    type(field_t), intent(inout) :: uu, vv, ww
    type(field_series_t), intent(inout) :: uulag, vvlag, wwlag
    real(kind=rp), dimension(3), intent(in) :: ab
    integer, intent(in) :: nab
    integer :: i, n

    n = uu%dof%size()

    do i = 1, n
       u%x(i,1,1,1) = ab(1) * uu%x(i,1,1,1) + ab(2) * uulag%lf(1)%x(i,1,1,1)
       v%x(i,1,1,1) = ab(1) * vv%x(i,1,1,1) + ab(2) * vvlag%lf(1)%x(i,1,1,1)
       w%x(i,1,1,1) = ab(1) * ww%x(i,1,1,1) + ab(2) * wwlag%lf(1)%x(i,1,1,1)
    end do

    if (nab .eq. 3) then
       do i = 1, n
          u%x(i,1,1,1) = u%x(i,1,1,1) + ab(3) * uulag%lf(2)%x(i,1,1,1)
          v%x(i,1,1,1) = v%x(i,1,1,1) + ab(3) * vvlag%lf(2)%x(i,1,1,1)
          w%x(i,1,1,1) = w%x(i,1,1,1) + ab(3) * wwlag%lf(2)%x(i,1,1,1)
       end do
    end if
    
  end subroutine rhs_maker_sumab_sx

  subroutine rhs_maker_ext_sx(temp1, temp2, temp3, fx_lag, fy_lag, fz_lag, &
                             fx_laglag, fy_laglag, fz_laglag, fx, fy, fz, &
                             rho, ext_coeffs, n)
    type(field_t), intent(inout) :: temp1, temp2, temp3
    type(field_t), intent(inout) :: fx_lag, fy_lag, fz_lag
    type(field_t), intent(inout) :: fx_laglag, fy_laglag, fz_laglag
    real(kind=rp), intent(inout) :: rho, ext_coeffs(10)
    integer, intent(in) :: n
    real(kind=rp), intent(inout) :: fx(n), fy(n), fz(n)
    integer :: i

    do i = 1, n
       temp1%x(i,1,1,1) = ext_coeffs(2) * fx_lag%x(i,1,1,1) + &
                          ext_coeffs(3) * fx_laglag%x(i,1,1,1)
       temp2%x(i,1,1,1) = ext_coeffs(2) * fy_lag%x(i,1,1,1) + &
                          ext_coeffs(3) * fy_laglag%x(i,1,1,1)
       temp3%x(i,1,1,1) = ext_coeffs(2) * fz_lag%x(i,1,1,1) + &
                          ext_coeffs(3) * fz_laglag%x(i,1,1,1)
    end do

    do i = 1, n
       fx_laglag%x(i,1,1,1) = fx_lag%x(i,1,1,1)
       fy_laglag%x(i,1,1,1) = fy_lag%x(i,1,1,1)
       fz_laglag%x(i,1,1,1) = fz_lag%x(i,1,1,1)
       fx_lag%x(i,1,1,1) = fx(i)
       fy_lag%x(i,1,1,1) = fy(i)
       fz_lag%x(i,1,1,1) = fz(i)
    end do

    do i = 1, n
       fx(i) = (ext_coeffs(1) * fx(i) + temp1%x(i,1,1,1)) * rho
       fy(i) = (ext_coeffs(1) * fy(i) + temp2%x(i,1,1,1)) * rho
       fz(i) = (ext_coeffs(1) * fz(i) + temp3%x(i,1,1,1)) * rho
    end do
    
  end subroutine rhs_maker_ext_sx

  subroutine scalar_rhs_maker_ext_sx(temp, fs_lag, fs_laglag, fs, rho, ext_coeffs, &
                                n)
    type(field_t), intent(inout) :: temp
    type(field_t), intent(inout) :: fs_lag
    type(field_t), intent(inout) :: fs_laglag
    real(kind=rp), intent(inout) :: rho, ext_coeffs(10)
    integer, intent(in) :: n
    real(kind=rp), intent(inout) :: fs(n)
    integer :: i

    do i = 1, n
       temp%x(i,1,1,1) = ext_coeffs(2) * fs_lag%x(i,1,1,1) + &
                          ext_coeffs(3) * fs_laglag%x(i,1,1,1)
    end do

    do i = 1, n
       fs_laglag%x(i,1,1,1) = fs_lag%x(i,1,1,1)
       fs_lag%x(i,1,1,1) = fs(i)
    end do

    do i = 1, n
       fs(i) = (ext_coeffs(1) * fs(i) + temp%x(i,1,1,1)) * rho
    end do
    
  end subroutine scalar_rhs_maker_ext_sx

  subroutine rhs_maker_bdf_sx(ta1, ta2, ta3, tb1, tb2, tb3, &
                               ulag, vlag, wlag, bfx, bfy, bfz, &
                               u, v, w, B, rho, dt, bd, nbd, n)    
    integer, intent(in) :: n, nbd
    type(field_t), intent(inout) :: ta1, ta2, ta3
    type(field_t), intent(in) :: u, v, w
    type(field_t), intent(inout) :: tb1, tb2, tb3
    type(field_series_t), intent(in) :: ulag, vlag, wlag        
    real(kind=rp), intent(inout) :: bfx(n), bfy(n), bfz(n)
    real(kind=rp), intent(in) :: B(n)
    real(kind=rp), intent(in) :: dt, rho, bd(10)
    integer :: i, ilag

    do i = 1, n
       tb1%x(i,1,1,1) = u%x(i,1,1,1) * B(i) * bd(2)
       tb2%x(i,1,1,1) = v%x(i,1,1,1) * B(i) * bd(2)
       tb3%x(i,1,1,1) = w%x(i,1,1,1) * B(i) * bd(2)
    end do

    do ilag = 2, nbd
       do i = 1, n
          ta1%x(i,1,1,1) = ulag%lf(ilag-1)%x(i,1,1,1) * B(i) * bd(ilag+1)
          ta2%x(i,1,1,1) = vlag%lf(ilag-1)%x(i,1,1,1) * B(i) * bd(ilag+1)
          ta3%x(i,1,1,1) = wlag%lf(ilag-1)%x(i,1,1,1) * B(i) * bd(ilag+1)
       end do

       do i = 1, n
          tb1%x(i,1,1,1) = tb1%x(i,1,1,1) + ta1%x(i,1,1,1)
          tb2%x(i,1,1,1) = tb2%x(i,1,1,1) + ta2%x(i,1,1,1)
          tb3%x(i,1,1,1) = tb3%x(i,1,1,1) + ta3%x(i,1,1,1)
       end do
    end do

    do i = 1, n
       bfx(i) = bfx(i) + tb1%x(i,1,1,1) * (rho / dt)
       bfy(i) = bfy(i) + tb2%x(i,1,1,1) * (rho / dt)
       bfz(i) = bfz(i) + tb3%x(i,1,1,1) * (rho / dt)
    end do

  end subroutine rhs_maker_bdf_sx

  subroutine scalar_rhs_maker_bdf_sx(ta1, tb1, slag, bfs, s, B, rho, dt, bd, nbd, n)    
    integer, intent(in) :: n, nbd
    type(field_t), intent(inout) :: ta1
    type(field_t), intent(in) :: s
    type(field_t), intent(inout) :: tb1
    type(field_series_t), intent(in) :: slag
    real(kind=rp), intent(inout) :: bfs(n)
    real(kind=rp), intent(in) :: B(n)
    real(kind=rp), intent(in) :: dt, rho, bd(10)
    integer :: i, ilag

    do i = 1, n
       tb1%x(i,1,1,1) = s%x(i,1,1,1) * B(i) * bd(2)
    end do

    do ilag = 2, nbd
       do i = 1, n
          ta1%x(i,1,1,1) = slag%lf(ilag-1)%x(i,1,1,1) * B(i) * bd(ilag+1)
       end do

       do i = 1, n
          tb1%x(i,1,1,1) = tb1%x(i,1,1,1) + ta1%x(i,1,1,1)
       end do
    end do

    do i = 1, n
       bfs(i) = bfs(i) + tb1%x(i,1,1,1) * (rho / dt)
    end do

  end subroutine scalar_rhs_maker_bdf_sx

end module rhs_maker_sx

