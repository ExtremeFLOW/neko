module user
  use neko
  use ieee_arithmetic, only: ieee_is_nan

  implicit none

contains

  ! Register user defined functions (see user_intf.f90)
  subroutine user_setup(usr)
    type(user_t), intent(inout) :: usr

    usr%fluid_user_ic => user_ic
    
  end subroutine user_setup

  ! User defined initial condition
  subroutine user_ic(u, v, w, p, params)
    type(field_t), intent(inout) :: u
    type(field_t), intent(inout) :: v
    type(field_t), intent(inout) :: w
    type(field_t), intent(inout) :: p
    type(json_file), intent(inout) :: params
    integer :: iel, ix,iy,iz
    real(kind=rp) :: fcoeff(3), xl(2)



    do iel = 1, u%msh%nelv
    do iz = 1, u%Xh%lz
    do iy = 1, u%Xh%ly
    do ix = 1, u%Xh%lx
    	   xl(1) = u%dof%x(ix,iy,iz,iel)
    	   xl(2) = u%dof%y(ix,iy,iz,iel)
         fcoeff(1)=  3.0e4
         fcoeff(2)= -1.5e3
         fcoeff(3)=  0.5e5
       	u%x(ix,iy,iz,iel) = math_ran_dst(ix,iy,iz,iel,xl,fcoeff)*1.0e-08
         fcoeff(1)=  2.3e4
         fcoeff(2)=  2.3e3
         fcoeff(3)= -2.0e5
       	v%x(ix,iy,iz,iel) = math_ran_dst(ix,iy,iz,iel,xl,fcoeff)*1.0e-08
       	w%x(ix,iy,iz,iel) = 0.0
    end do
    end do
    end do
    end do

  end subroutine user_ic



  !=======================================================================
!> @brief Give random distribution depending on position
!! @ingroup math
!! @details The original Nek5000 random number generator is implementted
!!  in @ref ran1. This totally ad-hoc random number generator below
!!  could be preferable to the original one for the simple reason that it
!!  gives the same initial cindition independent of the number of
!!  processors, which is important for code verification.
!! @param[in] ix,iy,iz     GLL point index
!! @param[in] ieg          global element number
!! @param[in] xl           physical point coordinates
!! @param[in] fcoeff       function coefficients
!! @return  random distribution
      real function math_ran_dst(ix,iy,iz,ieg,xl,fcoeff)
      implicit none


      ! argument list
      integer ix,iy,iz,ieg
      real(kind=rp) :: fcoeff(3), xl(2)
!-----------------------------------------------------------------------
      math_ran_dst = fcoeff(1)*(ieg+xl(1)*sin(xl(2))) + &
          fcoeff(2)*ix*iy + fcoeff(3)*ix
      math_ran_dst = 1.e3*sin(math_ran_dst)
      math_ran_dst = 1.e3*sin(math_ran_dst)
      math_ran_dst = cos(math_ran_dst)

      return
      end function math_ran_dst

end module user
