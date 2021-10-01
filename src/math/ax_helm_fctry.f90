module ax_helm_fctry
  use neko_config
  use ax_product
  use ax_helm_xsmm
  use ax_helm_sx
  use ax_helm
  implicit none

contains

  subroutine ax_helm_factory(Ax)
    class(ax_t), allocatable, intent(inout) :: Ax
    
    if (allocated(Ax)) then
       deallocate(Ax)
    end if

    if (NEKO_BCKND_SX .eq. 1) then
       allocate(ax_helm_sx_t::Ax)
    else if (NEKO_BCKND_XSMM .eq. 1) then
       allocate(ax_helm_xsmm_t::Ax)
    else
       allocate(ax_helm_t::Ax)
    end if

  end subroutine ax_helm_factory


end module ax_helm_fctry
