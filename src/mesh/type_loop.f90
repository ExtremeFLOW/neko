module base
  implicit none
  private

  public :: base_t

  type, abstract :: base_t
     integer :: id
  end type base_t
end module base


module type_loop
  use base, only : base_t
  implicit none
  
  type, abstract :: loop_1_t
     class(base_t), pointer :: ptr => null()
  end type loop_1_t

  type, extends(loop_1_t), abstract :: desc_t
  end type desc_t

  type :: alloc_t
     class(desc_t), allocatable :: obj
  end type alloc_t

  type, extends(base_t), abstract :: loop_2_t
     type(alloc_t) , dimension(:), allocatable :: cmp
  end type loop_2_t

end module type_loop
