!> Routines to interpolate between different spaces
module interpolation
  use speclib
  use device
  use utils
  use math
  use fast3d
  use tensor
  use space
  use device
  use, intrinsic :: iso_c_binding
  implicit none
  private
  
  type, public :: interpolator_t
     type(space_t), pointer :: Xh
     type(space_t), pointer :: Yh
     real(kind=rp), allocatable :: Xh_to_Yh(:,:), Xh_to_YhT(:,:)
     real(kind=rp), allocatable :: Yh_to_Xh(:,:), Yh_to_XhT(:,:)
     type(c_ptr) :: Xh_Yh_d = C_NULL_PTR
     type(c_ptr) :: Xh_YhT_d = C_NULL_PTR
     type(c_ptr) :: Yh_Xh_d = C_NULL_PTR
     type(c_ptr) :: Yh_XhT_d = C_NULL_PTR

   contains
     procedure, pass(this) :: init => interp_init
     procedure, pass(this) :: free => interp_free
     procedure, pass(this) :: map => interpolate
  end type interpolator_t
  
contains
  
  subroutine interp_init(this, Xh, Yh)
    class(interpolator_t), intent(inout) :: this
    type(space_t), intent(inout), target :: Xh
    type(space_t), intent(inout), target :: Yh
    integer :: deg_derivate

    call this%free()

    allocate(this%Xh_to_Yh(Yh%lx,Xh%lx))
    allocate(this%Xh_to_YhT(Xh%lx,Yh%lx))
    allocate(this%Yh_to_Xh(Xh%lx,Yh%lx))
    allocate(this%Yh_to_XhT(Yh%lx,Xh%lx))
    if (Xh%t .eq. GLL .and. Yh%t .eq. GLL) then
    else if ((Xh%t .eq. GL .and. Yh%t .eq. GLL) .or. &
         (Yh%t .eq. GL .and. Xh%t .eq. GLL)) then
    else
       call neko_error('Unsupported interpolation')
    end if
    deg_derivate = 0
    call setup_intp(this%Xh_to_Yh, this%Xh_to_YhT, &
         Yh%zg, Xh%zg, Yh%lx, Xh%lx, deg_derivate)
    call setup_intp(this%Yh_to_Xh, this%Yh_to_XhT, &
         Xh%zg, Yh%zg, Xh%lx, Yh%lx, deg_derivate)

    this%Xh => Xh
    this%Yh => Yh
    if ((NEKO_BCKND_HIP .eq. 1) .or. (NEKO_BCKND_CUDA .eq. 1) .or. &
         (NEKO_BCKND_OPENCL .eq. 1)) then
       call device_map(this%Xh_to_Yh, this%Xh_Yh_d, Yh%lx*Xh%lx)
       call device_map(this%Xh_to_YhT, this%Xh_YhT_d, Yh%lx*Xh%lx)
       call device_map(this%Yh_to_Xh, this%Yh_Xh_d, Yh%lx*Xh%lx)
       call device_map(this%Yh_to_XhT, this%Yh_XhT_d, Yh%lx*Xh%lx)
       call device_memcpy(this%Xh_to_Yh, this%Xh_Yh_d, Yh%lx*Xh%lx, HOST_TO_DEVICE)
       call device_memcpy(this%Xh_to_YhT, this%Xh_YhT_d, Yh%lx*Xh%lx, HOST_TO_DEVICE)
       call device_memcpy(this%Yh_to_Xh, this%Yh_Xh_d, Yh%lx*Xh%lx, HOST_TO_DEVICE)
       call device_memcpy(this%Yh_to_XhT, this%Yh_XhT_d, Yh%lx*Xh%lx, HOST_TO_DEVICE)
    end if

  end subroutine interp_init

  subroutine interp_free(this)
    class(interpolator_t), intent(inout) :: this

    if (allocated(this%Xh_to_Yh)) then
       deallocate(this%Xh_to_Yh)
    end if
    if (allocated(this%Xh_to_YhT)) then
       deallocate(this%Xh_to_YhT)
    end if
    if (allocated(this%Yh_to_Xh)) then
       deallocate(this%Yh_to_Xh)
    end if
    if (allocated(this%Yh_to_XhT)) then
       deallocate(this%Yh_to_XhT)
    end if
    if (c_associated(this%Yh_Xh_d)) then
       call device_free(this%Yh_Xh_d)
    end if
    if (c_associated(this%Yh_XhT_d)) then
       call device_free(this%Yh_XhT_d)
    end if
    if (c_associated(this%Xh_Yh_d)) then
       call device_free(this%Xh_Yh_d)
    end if
    if (c_associated(this%Xh_YhT_d)) then
       call device_free(this%Xh_YhT_d)
    end if

  end subroutine interp_free

  !> Interpolates array x -> y in to_space
  subroutine interpolate(this, y, x, nel,to_space)
    class(interpolator_t), intent(inout) :: this
    integer :: nel
    type(space_t) :: to_space
    real(kind=rp), intent(inout) :: x(1,nel)
    real(kind=rp), intent(inout) :: y(1,nel)
    if (to_space .eq. this%Yh) then
       call tnsr3d(y, this%Yh%lx, x, &
                   this%Xh%lx,this%Yh_to_XhT, &
                   this%Yh_to_Xh, this%Yh_to_Xh, nel)
    else if (to_space .eq. this%Xh) then
       call tnsr3d(y, this%Xh%lx, x, &
                   this%Yh%lx,this%Yh_to_Xh, &
                   this%Yh_to_XhT, this%Yh_to_XhT, nel)
    else
       call neko_error('Invalid interpolation')
    end if
  end subroutine interpolate

end module interpolation
