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
  type(field_t), pointer :: uzt, eps_t, eps_k
  type(field_t), pointer :: s, ss, s_rms
  
  !> declare variables
  integer :: i,j,k,n, reader, probe
  real(kind=rp) :: Ra
  real(kind=rp) :: Pr

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

  !> Initialize the file interpolator object and read first file
  call init_readers(pd)
  
  !> Initialize work fields and other misc items
  call work_field%init(pd%fld_ioc(1)%dof)
  call work_field2%init(pd%fld_ioc(1)%dof)
  call neko_field_registry%add_field(pd%fld_ioc(1)%dof, "t_rms")
  
  !> Initialize the standalone simulation components 
  !! Do this after all variables you want in registry have been added
  call init_probes(pd)

  !> --------------------
  !> Compute phase
  !> --------------------
  

  do reader = 1, size(pd%fld_ioc)

     do i=0, pd%fld_ioc(reader)%number_of_files -1

        !!> --------------------
        !!> Read data
        !!> --------------------
        !!  Read the data in the next file but do not do it for i=0 since the file data is already in memory from init 
        if (i.ne.0) call pd%fld_ioc(reader)%step()
  
        !!> --------------------
        !!> Average the quantities of all files in all readers
        !!> --------------------
        call pd%fld_ioc(reader)%average_registry(pd%fld_ioc(reader)%t, sync=.true.)
        
     end do
  end do

  !!> --------------------
  !!> Put all averages in the registry for ease of computing
  !!> --------------------
  !> For all the readers, mode the data from the averaging object to their corresponding
  !! fields in the registry
  do reader = 1, size(pd%fld_ioc)
     !> Move the data from mean_fields to the registry so it can be interpolated
     call pd%fld_ioc(reader)%put_averages_in_registry()
  end do

  !!> --------------------
  !!> Calculate anything you want to interpolate from the averages
  !!> --------------------
  s   => neko_field_registry%get_field("t")
  ss => neko_field_registry%get_field("tt")
  s_rms => neko_field_registry%get_field("t_rms")
  call calculate_rms(s_rms, s, ss, work_field) 

  !!> --------------------
  !!> Interpolate the averaged quantities
  !!> --------------------
  do probe = 1, size(pd%pb)
     ! here assume one time for all the probes interpolated
     call pd%pb(probe)%compute_(pd%fld_ioc(1)%t, pd%fld_ioc(1)%file_counter)
  end do
  
  !!> --------------------
  !!> Calculate nusselt
  !!> --------------------
  uzt   => neko_field_registry%get_field("uzt")
  eps_t => neko_field_registry%get_field("eps_t")
  eps_k => neko_field_registry%get_field("eps_k")
  call calculate_nusselt_from_fields(uzt, eps_t, eps_k, work_field, work_field2, &
                                     pd%fld_ioc(1)%coef, Ra, Pr) 

  !!> --------------------
  !!> Find the ratio between the mesh spacing and kolmogorov scale
  !!> --------------------
  call compare_mesh_kolmogorov(eps_t, eps_k, work_field, work_field2, &
                               pd%fld_ioc(1)%coef, pd%fld_ioc(1)%dof, &
                               pd%fld_ioc(1)%msh, pd%fld_ioc(1)%Xh, Ra, Pr)

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
