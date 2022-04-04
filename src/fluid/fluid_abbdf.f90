module fluid_abbdf
  use num_types
  use field_series, only : field_series_t
  use field, only : field_t
  implicit none 

  !> Abstract type to sum up AB/BDF contributions
  type, abstract :: fluid_sumab_t
   contains
     procedure(fluid_sumab), nopass, deferred :: compute
  end type fluid_sumab_t

  !> Abstract type to sum up contributions to kth order extrapolation scheme
  type, abstract :: fluid_makeabf_t
   contains
     procedure(fluid_makeabf), nopass, deferred :: compute
  end type fluid_makeabf_t

  !> Abstract type to add contributions to F from lagged BD terms
  type, abstract :: fluid_makebdf_t
   contains
     procedure(fluid_makebdf), nopass, deferred :: compute
  end type fluid_makebdf_t

  abstract interface
     subroutine fluid_sumab(u, v, w, uu, vv, ww, uulag, vvlag, wwlag, ab, nab)
       import field_t
       import field_series_t
       import rp
       type(field_t), intent(inout) :: u,v, w
       type(field_t), intent(inout) :: uu, vv, ww
       type(field_series_t), intent(inout) :: uulag, vvlag, wwlag
       real(kind=rp), dimension(3), intent(in) :: ab
       integer, intent(in) :: nab
     end subroutine fluid_sumab
  end interface

  abstract interface
     subroutine fluid_makeabf(ta1, ta2, ta3, abx1, aby1, abz1, &
                              abx2, aby2, abz2, bfx, bfy, bfz, rho, ab, n)
       import field_t
       import rp
       type(field_t), intent(inout) :: ta1, ta2, ta3
       type(field_t), intent(inout) :: abx1, aby1, abz1
       type(field_t), intent(inout) :: abx2, aby2, abz2
       real(kind=rp), intent(inout) :: rho, ab(10)
       integer, intent(in) :: n
       real(kind=rp), intent(inout) :: bfx(n), bfy(n), bfz(n)
     end subroutine fluid_makeabf
  end interface

  abstract interface
     subroutine fluid_makebdf(ta1, ta2, ta3, tb1, tb2, tb3, &
                              ulag, vlag, wlag, bfx, bfy, bfz, &
                              u, v, w, B, rho, dt, bd, nbd, n)
       import field_series_t
       import field_t
       import rp
       integer, intent(in) :: n, nbd
       type(field_t), intent(inout) :: ta1, ta2, ta3
       type(field_t), intent(in) :: u, v, w
       type(field_t), intent(inout) :: tb1, tb2, tb3
       type(field_series_t), intent(in) :: ulag, vlag, wlag        
       real(kind=rp), intent(inout) :: bfx(n), bfy(n), bfz(n)
       real(kind=rp), intent(in) :: B(n)
       real(kind=rp), intent(in) :: dt, rho, bd(10)
     end subroutine fluid_makebdf
  end interface

end module fluid_abbdf
