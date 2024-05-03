!> Defines structs that are used... Dont know if we should keep it though.
module structs
  use num_types, only : rp, dp
  implicit none
  private

  type, public :: struct_curve_t
     real(kind=dp) :: curve_data(5,12)
     integer :: curve_type(12)
     integer :: el_idx
  end type struct_curve_t

  !> Pointer to array
  type, public :: array_ptr_t
     real(kind=rp), pointer :: ptr(:)
  end type array_ptr_t
end module structs

