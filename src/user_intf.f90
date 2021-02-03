!> Interfaces for user interaction with NEKO
module user_intf
  use field
  use source
  use parameters
  implicit none

  !> Abstract interface for user defined initial conditions
  abstract interface
     subroutine useric(u, v, w, p, params)
       import field_t
       import param_t
       type(field_t), intent(inout) :: u
       type(field_t), intent(inout) :: v
       type(field_t), intent(inout) :: w
       type(field_t), intent(inout) :: p
       type(param_t), intent(inout) :: params
     end subroutine useric
  end interface

  type :: user_t
     procedure(useric), nopass, pointer :: fluid_usr_ic => null()
     procedure(source_term_pw), nopass, pointer :: fluid_usr_f => null()
   contains
     procedure, pass(u) :: init => user_intf_init
  end type user_t
  
contains

  !> User interface initialization
  subroutine user_intf_init(u)
    class(user_t), intent(inout) :: u

    if (.not. associated(u%fluid_usr_ic)) then
       u%fluid_usr_ic => dummy_user_ic
    end if

    if (.not. associated(u%fluid_usr_f)) then
       u%fluid_usr_f => dummy_user_f
    end if
    
  end subroutine user_intf_init

  
  !
  ! Below is the dummy user interface
  ! when running in pure turboNEKO mode
  !

  !> Dummy user initial condition
  subroutine dummy_user_ic(u, v, w, p, params)
    type(field_t), intent(inout) :: u
    type(field_t), intent(inout) :: v
    type(field_t), intent(inout) :: w
    type(field_t), intent(inout) :: p
    type(param_t), intent(inout) :: params
    call neko_error('Dummy user defined initial condition set')    
  end subroutine dummy_user_ic

  !> Dummy user forcing
  subroutine dummy_user_f(u, v, w, j, k, l, e)
    real(kind=dp), intent(inout) :: u
    real(kind=dp), intent(inout) :: v
    real(kind=dp), intent(inout) :: w
    integer, intent(inout) :: j
    integer, intent(inout) :: k
    integer, intent(inout) :: l
    integer, intent(inout) :: e
    call neko_error('Dummy user defined forcing set')    
  end subroutine dummy_user_f
  
end module user_intf
