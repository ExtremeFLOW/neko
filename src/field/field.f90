!> Defines a field
!
module field
  use neko_config
  use num_types
  use device
  use math
  use mesh
  use space
  use dofmap
  use, intrinsic :: iso_c_binding
  implicit none
  
  type field_t
     real(kind=rp), allocatable :: x(:,:,:,:) !< Field data
     
     type(space_t), pointer :: Xh   !< Function space \f$ X_h \f$
     type(mesh_t), pointer :: msh   !< Mesh
     type(dofmap_t), pointer :: dof !< Dofmap

     logical :: internal_dofmap = .false. !< Does the field have an own dofmap
     character(len=80) :: name            !< Name of the field
     type(c_ptr) :: x_d
  end type field_t

  interface field_init
     module procedure field_init_external_dof, field_init_internal_dof
  end interface field_init
  
  interface assignment(=)
     module procedure field_assign_field, field_assign_scalar
  end interface assignment(=)

  interface field_add
     module procedure field_add_field, field_add_scalar
  end interface field_add

contains

  !> Initialize a field @a f on the mesh @a msh using an internal dofmap
  subroutine field_init_internal_dof(f, msh, space, fld_name)
    type(field_t), intent(inout) :: f       !< Field to be initialized
    type(mesh_t), target, intent(in) :: msh !< underlying mesh of the field
    type(space_t), target, intent(in) :: space !< Function space for the field
    character(len=*), optional :: fld_name     !< Name of the field

    call field_free(f)

    f%Xh => space
    f%msh => msh

    allocate(f%dof)
    f%dof = dofmap_t(f%msh, f%Xh)
    f%internal_dofmap = .true.
    
    if (present(fld_name)) then
       call field_init_common(f, fld_name)
    else
       call field_init_common(f)
    end if
    
  end subroutine field_init_internal_dof

  !> Initialize a field @a f on the mesh @a msh using an internal dofmap
  subroutine field_init_external_dof(f, dof, fld_name)
    type(field_t), intent(inout) :: f       !< Field to be initialized
    type(dofmap_t), target, intent(in) :: dof  !< External dofmap for the field
    character(len=*), optional :: fld_name     !< Name of the field

    call field_free(f)

    f%Xh => dof%Xh
    f%msh => dof%msh
    f%dof => dof

    if (present(fld_name)) then
       call field_init_common(f, fld_name)
    else
       call field_init_common(f)
    end if
    
  end subroutine field_init_external_dof

  !> Initialize a field @a f 
  subroutine field_init_common(f, fld_name)
    type(field_t), intent(inout) :: f       !< Field to be initialized
    character(len=*), optional :: fld_name  !< Name of the field
    integer :: ierr
    integer :: lx, ly, lz, nelv, n

    lx = f%Xh%lx
    ly = f%Xh%ly
    lz = f%Xh%lz
    nelv = f%msh%nelv
    
    if (.not. allocated(f%x)) then
       allocate(f%x(lx, ly, lz, nelv), stat = ierr)        
       f%x = 0d0
    end if

    if (present(fld_name)) then
       f%name = fld_name
    else
       f%name = "Field"
    end if

    if ((NEKO_BCKND_HIP .eq. 1) .or. (NEKO_BCKND_CUDA .eq. 1)) then
       n = lx * ly * lz * nelv           
       call device_map(f%x, f%x_d, n)
    end if
    
  end subroutine field_init_common

  !> Deallocate a field @a f
  subroutine field_free(f)
    type(field_t), intent(inout) :: f
    
    if (allocated(f%x)) then
       deallocate(f%x)
    end if

    if (f%internal_dofmap) then
       deallocate(f%dof)
       f%internal_dofmap = .false.
    end if
    
    nullify(f%msh)
    nullify(f%Xh)
    nullify(f%dof)

    if (c_associated(f%x_d)) then
       call device_free(f%x_d)
    end if

  end subroutine field_free

  !> Assignment \f$ F = G \f$
  !! @note @a F will be initialized if it has a different size than
  !! @a G or it's not allocated
  subroutine field_assign_field(f, g)
    type(field_t), intent(inout) :: f
    type(field_t), intent(in) :: g
    integer :: n

    if (allocated(f%x) .and. (f%Xh .ne. g%Xh)) then
       call field_free(f)
    end if
    
    f%Xh =>g%Xh
    f%msh => g%msh
    f%dof => g%dof
    
    
    f%Xh%lx = g%Xh%lx
    f%Xh%ly = g%Xh%ly
    f%Xh%lz = g%Xh%lz
    
    if (.not. allocated(f%x)) then
       allocate(f%x(f%Xh%lx, f%Xh%ly, f%Xh%lz, f%msh%nelv))
    end if
    
    n = f%msh%nelv * f%Xh%lx * f%Xh%ly * f%Xh%lz
    call copy(f%x, g%x, n)
    
  end subroutine field_assign_field

  !> Assignment \f$ F = a \f$
  subroutine field_assign_scalar(f, a)
    type(field_t), intent(inout) :: f
    real(kind=rp), intent(in) :: a
    integer :: n, i, j, k, l

    n = f%msh%nelv * f%Xh%lx * f%Xh%ly * f%Xh%lz
    do i = 1, f%msh%nelv
       do l = 1, f%Xh%lz
          do k = 1, f%Xh%ly
             do j = 1, f%Xh%lx
                f%x(j, k, l, i) = a
             end do
          end do
       end do
    end do

  end subroutine field_assign_scalar
  
  !> Add \f$ F(u_1, u_2, ... , u_n) =
  !! F(u_1, u_2, ... , u_n) + G(u_1, u_2, ... , u_n) \f$
  !! @note Component wise
  subroutine field_add_field(f, g)
    type(field_t), intent(inout) :: f
    type(field_t), intent(inout) :: g
    integer :: n

    n = f%msh%nelv * f%Xh%lx * f%Xh%ly * f%Xh%lz
    call add2(f%x, g%x, n)

  end subroutine field_add_field


  !> Add \f$ F(u_1, u_2, ... , u_n) =
  !! F(u_1, u_2, ... , u_n) + a \f$
  subroutine field_add_scalar(f, a)
    type(field_t), intent(inout) :: f
    real(kind=rp), intent(inout) :: a
    integer :: n

    n = f%msh%nelv * f%Xh%lx * f%Xh%ly * f%Xh%lz
    call cadd(f%x, a, n)

  end subroutine field_add_scalar

end module field

