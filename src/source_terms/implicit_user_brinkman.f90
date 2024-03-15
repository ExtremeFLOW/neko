module implicit_user_brinkman

	use field, only: field_t 
	use field_registry, only: neko_field_registry
	use coefs, only: coef_t
   use user_intf, only: user_t
   use num_types, only: rp
	private

	public :: implicit_brinkman

	type, public :: implicit_user_brinkman_t
		type(field_t), pointer :: chi
	contains
		procedure, pass(this) :: init => implicit_user_brinkman_init
		procedure, pass(this) :: compute => implicit_user_brinkman_compute
		procedure(implicit_brinkman), nopass, pointer :: compute_user_brinkman => null()		
   end type implicit_user_brinkman_t

	 ! Harry ---------------------------------------
  !> Abstract interface for user implicit Brinkman forcing
  abstract interface
     subroutine implicit_brinkman(chi, t)
       import field_t
       import rp
       type(field_t), intent(inout) :: chi
    	 real(kind=rp), intent(in) :: t
     end subroutine implicit_brinkman
  end interface
	 ! ---------------------------------------------

	contains

	subroutine implicit_user_brinkman_init(this, coef, user_implicit_brinkman)
		class(implicit_user_brinkman_t), intent(inout) :: this
		type(coef_t), intent(in) :: coef
		procedure(implicit_brinkman) :: user_implicit_brinkman

		call neko_field_registry%add_field(coef%dof, "chi")
		this%chi => neko_field_registry%get_field("chi")

		this%compute_user_brinkman => user_implicit_brinkman

	end subroutine implicit_user_brinkman_init
	
	subroutine implicit_user_brinkman_compute(this, t)
	implicit none
		class(implicit_user_brinkman_t), intent(inout) :: this
		!real(kind=rp), intent(in) :: t
		real(kind=rp) :: t
		integer :: i

		!call user%compute(this%chi, t)
		! VICTOR WILL FIGURE THIS OUT

		call this%compute_user_brinkman(this%chi, t)

!		do i=1, this%chi%dof%size()
!		  if(this%chi%dof%y(i,1,1,1).lt.0.5) then
!		  	 this%chi%x(i,1,1,1) = 0.
!		  else
!		  	 this%chi%x(i,1,1,1) = 100000.
!		  endif
!		enddo  

	end subroutine implicit_user_brinkman_compute

end module implicit_user_brinkman
