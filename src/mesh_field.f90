!> Defines a mesh field
!! @details A mesh field is a scalar integer cell based field (\f$ dQ_0 \f$)
module mesh_field
  use num_types
  use mesh
  implicit none
  
  !> @todo Add support for different data types
  type mesh_fld_t
     integer, allocatable :: data(:) !< Data
     type(mesh_t), pointer :: msh !< Mesh
  end type mesh_fld_t

contains
  
  subroutine mesh_field_init(fld, msh)
    type(mesh_fld_t), intent(inout) :: fld
    type(mesh_t), target, intent(in) :: msh

    call mesh_field_free(fld)

    fld%msh => msh
    if (.not. allocated(fld%data)) then
       allocate(fld%data(msh%nelv))
    end if

    fld%data = 0
  end subroutine mesh_field_init

  subroutine mesh_field_free(fld)
    type(mesh_fld_t), intent(inout) :: fld

    if (allocated(fld%data)) then
       deallocate(fld%data)
    end if

    nullify(fld%msh)
  end subroutine mesh_field_free

end module mesh_field
