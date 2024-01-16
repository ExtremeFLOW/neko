program evaluate_mesh 
  use neko
  use data_processor

  implicit none

  !> --------------------
  !> Declaration phase 
  !> --------------------
  
  !> Declare objects
  type(data_processor_t) :: pd

  !> Declare fields
  type(field_t) :: work_field
  type(field_t) :: work_field2
  type(field_t), pointer :: u,v,w
  type(field_t), pointer :: tau_x, tau_y
  
  !> declare variables
  integer :: i,j,k,n, reader, probe
  real(kind=rp) :: Ra
  real(kind=rp) :: Pr
  real(kind=rp) :: mu

  !> --------------------
  !> Initialization phase 
  !> --------------------

  !> Initialize neko 
  call neko_init 

  !> Initialize the parameters
  call init_params(pd%params)
 
  !> Initialize the file interpolator object and read first file
  call init_readers(pd)
  
  !> Initialize work fields and other misc items
  call work_field%init(pd%fld_ioc(1)%dof)
  call work_field2%init(pd%fld_ioc(1)%dof)
  
  
  !> --------------------
  !> Compute phase phase
  !> --------------------
  call check_mesh_holes(pd%fld_ioc(1)%coef, pd%fld_ioc(1)%dof, &
                        pd%fld_ioc(1)%msh, pd%fld_ioc(1)%Xh, &
                        pd%fld_ioc(1)%gs_h)



  !> --------------------
  !> Finalization phase
  !> --------------------
  
  call work_field%free()
  call work_field2%free()


  !> Finalize neko
  call neko_finalize

end program evaluate_mesh
