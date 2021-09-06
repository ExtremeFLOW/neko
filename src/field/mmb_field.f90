module mmb_field
  use mamba
  use field
  use num_types
  implicit none


  type, public, extends(field_t) ::  mmb_field_t
     type(mmbArray) :: mba_x
     type(mmbLayout) :: layout
     type(mmbTileIterator) :: mba_it
   contains
     procedure, pass(this) :: free => mmb_field_free
     procedure, pass(this) :: init => mmb_field_init_internal_dof, &
          mmb_field_init_external_dof
  end type mmb_field_t

  interface mmb_field_init
     module procedure mmb_field_init_external_dof, mmb_field_init_internal_dof
  end interface mmb_field_init

  
contains

  subroutine mmb_field_init_internal_dof(this, msh, space, fld_name)
    class(mmb_field_t), intent(inout) :: this
    type(mesh_t), target, intent(in) :: msh 
    type(space_t), target, intent(in) :: space 
    character(len=*), optional :: fld_name

    call mmb_field_free(this)

    this%Xh => space
    this%msh => msh

    allocate(this%dof)
    this%dof = dofmap_t(this%msh, this%Xh)
    this%internal_dofmap = .true.
    
    if (present(fld_name)) then
       call mmb_field_init_common(this, fld_name)
    else
       call mmb_field_init_common(this)
    end if

  end subroutine mmb_field_init_internal_dof

  subroutine mmb_field_init_external_dof(this, dof, fld_name)
    class(mmb_field_t), intent(inout) :: this
    type(dofmap_t), target, intent(in) :: dof
    character(len=*), optional :: fld_name

    call mmb_field_free(this)
    
    this%Xh => dof%Xh
    this%msh => dof%msh
    this%dof => dof

    if (present(fld_name)) then
       call mmb_field_init_common(this, fld_name)
    else
       call mmb_field_init_common(this)
    end if
    
  end subroutine mmb_field_init_external_dof

  subroutine mmb_field_init_common(f, fld_name)
    type(mmb_field_t), intent(inout) :: f
    character(len=*), optional :: fld_name  
    integer :: ierr
    integer :: lx, ly, lz, nelv
    integer(mmbErrorKind) :: err
    integer(mmbIndexKind), dimension(4) :: dims
  
    
    lx = f%Xh%lx
    ly = f%Xh%ly
    lz = f%Xh%lz
    nelv = f%msh%nelv
    
    if (.not. allocated(f%x)) then
       allocate(f%x(lx, ly, lz, nelv), stat = ierr)        
       f%x = 0d0
    end if

   

    dims = [lx, ly, lz, nelv]
    call mmb_layout_create_regular_nd(int(storage_size(1.0)/8, mmbSizeKind), &
         4_mmbSizeKind, MMB_COLMAJOR, mmb_layout_padding_create_zero(),&
         f%layout, err)
    if (err .ne. MMB_OK) call neko_error('Mamba layout fail')

    call mmb_array_create_wrapped(f%x, dims, f%layout, &
         dram_interface, MMB_READ_WRITE, f%mba_x, err)
    if (err .ne. MMB_OK) call neko_error('Mamba array wrap. fail')

    call mmb_array_tile(f%mba_x, dims, err)

    call mmb_tile_iterator_create(f%mba_x, f%mba_it, err)

    if (present(fld_name)) then
       f%name = fld_name
    else
       f%name = "MAMBA Field"
    end if

  end subroutine mmb_field_init_common


  subroutine mmb_field_free(this)
    class(mmb_field_t), intent(inout) :: this

    if (this%internal_dofmap) then
       deallocate(this%dof)
       this%internal_dofmap = .false.
    end if
    
    nullify(this%msh)
    nullify(this%Xh)
    nullify(this%dof)

    call mmb_array_destroy(this%mba_x)
   
  end subroutine mmb_field_free
  
end module mmb_field
