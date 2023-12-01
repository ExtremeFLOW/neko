program post_process
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
 
  !> Read from case file
  call json_get(pd%params, 'case.non_dimensional_quantities.Ra', Ra)
  call json_get(pd%params, 'case.non_dimensional_quantities.Pr', Pr)
  mu = sqrt(pr/ra)

  !> Initialize the file interpolator object and read first file
  call init_readers(pd)
  
  !> Initialize work fields and other misc items
  call work_field%init(pd%fld_ioc(1)%dof)
  call work_field2%init(pd%fld_ioc(1)%dof)
  call neko_field_registry%add_field(pd%fld_ioc(1)%dof, "tau_x")
  call neko_field_registry%add_field(pd%fld_ioc(1)%dof, "tau_y")
  
  !> Initialize the standalone simulation components 
  !! Do this after all variables you want in registry have been added
  call init_probes(pd)

  !> --------------------
  !> Compute phase
  !> --------------------
  u => neko_field_registry%get_field("u")
  v => neko_field_registry%get_field("v")
  tau_x => neko_field_registry%get_field("tau_x")
  tau_y => neko_field_registry%get_field("tau_y")


  do reader = 1, size(pd%fld_ioc)

     do i=0, pd%fld_ioc(reader)%number_of_files -1

        !!> --------------------
        !!> Read data
        !!> --------------------
        !!  Read the data in the next file but do not do it for i=0 since the file data is already in memory from init 
        if (i.ne.0) call pd%fld_ioc(reader)%step()
  
        !!> --------------------
        !!> Calculate shear stress in relevant directions
        !!> --------------------
        call calculate_shear_stress(tau_x, u, pd%fld_ioc(1)%coef, "z", mu)  
        call calculate_shear_stress(tau_y, v, pd%fld_ioc(1)%coef, "z", mu)  
       
        !!> --------------------
        !!> Interpolate selected fields in registry
        !!> --------------------
        do probe = 1, size(pd%pb)
           ! here assume one time for all the probes interpolated
           call pd%pb(probe)%compute_(pd%fld_ioc(1)%t, pd%fld_ioc(1)%file_counter)
        end do
 
     end do
  end do
 
  !> --------------------
  !> Finalization phase
  !> --------------------
  
  call work_field%free()
  call work_field2%free()

  !> Finalize the probes
  do probe = 1, size(pd%pb)
     call pd%pb(probe)%free()
  end do

  !> Finalize neko
  call neko_finalize

end program post_process
