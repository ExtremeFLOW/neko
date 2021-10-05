!> Fast Diagonalization 
module fdm_cpu
  use num_types
  use tensor_cpu
  implicit none

contains

  subroutine fdm_do_fast_cpu(e, r, s, d, nl, ldim, nelv)
    integer, intent(in) :: nl, nelv, ldim
    real(kind=rp), intent(inout) :: e(nl**ldim, nelv)
    real(kind=rp), intent(inout) :: r(nl**ldim, nelv)
    real(kind=rp), intent(inout) :: s(nl*nl,2,ldim, nelv)
    real(kind=rp), intent(inout) :: d(nl**ldim, nelv)    
    integer ::  ie, nn, i

    nn = nl**ldim
    if(.not. ldim .eq. 3) then
       do ie = 1, nelv
          call tnsr2d_el_cpu(e(1,ie), nl, r(1,ie), nl, s(1,2,1,ie), s(1,1,2,ie))
          do i = 1, nn
             r(i,ie) = d(i,ie) * e(i,ie)
          end do
          call tnsr2d_el_cpu(e(1,ie), nl, r(1,ie), nl, s(1,1,1,ie), s(1,2,2,ie))
       end do
    else
       do ie = 1, nelv
          call tnsr3d_el_cpu(e(1,ie), nl, r(1,ie), nl, &
               s(1,2,1,ie), s(1,1,2,ie), s(1,1,3,ie))
          do i = 1, nn
             r(i,ie) = d(i,ie) * e(i,ie)
          end do
          call tnsr3d_el_cpu(e(1,ie), nl, r(1,ie), nl, &
               s(1,1,1,ie), s(1,2,2,ie), s(1,2,3,ie))
       end do
    end if
  end subroutine fdm_do_fast_cpu
  
end module fdm_cpu
