!> Defines a statistical quantity
module stats_quant
  use num_types
  implicit none

  !> Abstract type defining a statistical quantity
  type, abstract :: stat_quant_t
   contains
     procedure(stat_quant_update), pass(this), deferred :: update
  end type stat_quant_t

  !> Abstract interface for updating/adding data to a quantitiy
  abstract interface
     subroutine stat_quant_update(this, k)
       import :: stat_quant_t
       import dp
       class(stat_quant_t), intent(inout) :: this
       real(kind=dp), intent(in) :: k
     end subroutine stat_quant_update
  end interface
  
end module stats_quant
