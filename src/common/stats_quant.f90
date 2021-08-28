!> Defines a statistical quantity
module stats_quant
  use num_types
  implicit none

  !> Abstract type defining a statistical quantity
  type, abstract :: stats_quant_t
   contains
     procedure(stats_quant_update), pass(this), deferred :: update
  end type stats_quant_t

  !> Abstract interface for updating/adding data to a quantitiy
  abstract interface
     subroutine stats_quant_update(this, k)
       import :: stats_quant_t
       import dp
       class(stats_quant_t), intent(inout) :: this
       real(kind=dp), intent(in) :: k
     end subroutine stats_quant_update
  end interface
  
end module stats_quant
