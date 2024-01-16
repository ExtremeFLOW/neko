module post_processing_library
  use neko
  use json_utils
  use fld_io_controller, only : fld_io_controller_t

  implicit none
    
contains
  
  subroutine calculate_rms(s_rms, s, ss, work_field) 
     type(field_t), intent(inout) :: s, ss, s_rms
     type(field_t), intent(inout) :: work_field

     !> Compute t_rms
     !Do it also in the gpu to verufy that all is good
     if (NEKO_BCKND_DEVICE .eq. 1) then 
        !  First get mean(t)^2
        call device_col3(work_field%x_d, s%x_d,  s%x_d, s%dof%size())
        ! then substract mean(tt) - mean(t)^2
        call device_sub3(s_rms%x_d, ss%x_d,  work_field%x_d, s%dof%size())
     else
        !  First get mean(t)^2
        call col3(work_field%x, s%x,  s%x, s%dof%size())
        ! then substract mean(tt) - mean(t)^2
        call sub3(s_rms%x, ss%x,  work_field%x, s%dof%size())
     end if

  end subroutine calculate_rms


  subroutine calculate_nusselt_from_fields(uzt, eps_t, eps_k, work_field, work_field2, coef, ra, pr) 
     type(coef_t), intent(inout) :: coef
     type(field_t), intent(inout) :: uzt, eps_t, eps_k
     type(field_t), intent(inout) :: work_field
     type(field_t), intent(inout) :: work_field2
     real(kind=rp), intent(inout) :: ra
     real(kind=rp), intent(inout) :: pr

     real(kind=rp) :: lambda, bar_uzt, nu_uzt
     real(kind=rp) :: mu, bar_eps_t, bar_eps_k, nu_eps_t, nu_eps_k
     !> declare variables
     integer :: i,j,k,n
 
     mu = sqrt(pr/ra)
     lambda = 1_rp/sqrt(ra*pr)
  
     !> Integrate thermal dissipation
     bar_eps_t = 0_rp
     call average_from_weights(bar_eps_t, eps_t, &
                     work_field, coef%B, coef%volume)
     !! Calculate nusselt from thermal dissipaton
     nu_eps_t = 1_rp/lambda*bar_eps_t
  
     !> Integrate kinetic
     bar_eps_k = 0_rp
     call average_from_weights(bar_eps_k, eps_k, &
                     work_field, coef%B, coef%volume)
  
     !! Calculate nusselt from kinetic dissipaton
     nu_eps_k = 1_rp + bar_eps_k/((mu**3_rp)*ra/(pr**2_rp))
     
     !> Integrate convective current
     bar_uzt = 0_rp
     call average_from_weights(bar_uzt, uzt, &
                     work_field, coef%B, coef%volume)
  
     !! Calculate nusselt from convective current
     nu_uzt = 1_rp + sqrt(ra)*bar_uzt

     write(*,*) " field in registry nu_eps_t,  nu_eps_k and nu_uzt are:", nu_eps_t, nu_eps_k, nu_uzt

  end subroutine calculate_nusselt_from_fields

  subroutine average_from_weights(avrg, field, work_field, weights, sum_weights)
     real(kind=rp), intent(inout) :: avrg
     type(field_t), intent(inout) :: field
     type(field_t), intent(inout) :: work_field
     real(kind=rp), dimension(field%Xh%lxyz,field%msh%nelv), intent(inout) :: weights
     real(kind=rp), intent(inout) :: sum_weights
     integer :: n
     type(c_ptr) :: weights_d
    
     n = field%dof%size()

     if (NEKO_BCKND_DEVICE .eq. 1) then 
        weights_d = device_get_ptr(weights)
        call device_col3(work_field%x_d,field%x_d,weights_d,n)      
        avrg = device_glsum(work_field%x_d,n)                  
        avrg = avrg / abs(sum_weights)            
     else
        call col3(work_field%x,field%x,weights,n)      
        avrg = glsum(work_field%x,n)                  
        avrg = avrg / abs(sum_weights)            
     end if

  end subroutine average_from_weights

  function assert_gpu_cpu_synch(field, work_field, tol) result(report)
     type(field_t), intent(inout) :: field
     type(field_t), intent(inout) :: work_field
     real(kind=rp), intent(inout) :: tol
     logical :: report
     integer :: i, n

     n = field%dof%size()

     !> Make sure the work_fields are zeroed
     call rzero(work_field%x,n)
     call device_rzero(work_field%x_d,n)

     !> Copy the contents of the gpu to the work field
     call device_copy(work_field%x_d , field%x_d, n)

     !> Copy the information from the gpu to the cpu
     call device_memcpy(work_field%x, work_field%x_d, field%dof%size(), &
                        DEVICE_TO_HOST, sync=.true.)

     !> Now assert
     report = .true.
     do i=1, n
        if (abs(field%x(i,1,1,1)-work_field%x(i,1,1,1)) .gt. tol ) then
           report = .false.
        end if
     end do
  end function assert_gpu_cpu_synch

  function assert_cpu_count_entries_less_than_tol(field, work_field, tol) result(report)
     type(field_t), intent(inout) :: field
     type(field_t), intent(inout) :: work_field
     real(kind=rp), intent(inout) :: tol
     integer :: report
     integer :: i, n

     n = field%dof%size()
    
     !> Make sure the work_fields are zeroed
     call rzero(work_field%x,n)
     call device_rzero(work_field%x_d,n)

     !> Copy the contents of the cpu to the work field
     call copy(work_field%x , field%x, n)

     !> Now assert
     report = 0
     do i=1, n
        if ( work_field%x(i,1,1,1) .lt. tol ) then
           report = report + 1
        end if
     end do

  end function assert_cpu_count_entries_less_than_tol

  function assert_gpu_count_entries_less_than_tol(field, work_field, tol) result(report)
     type(field_t), intent(inout) :: field
     type(field_t), intent(inout) :: work_field
     real(kind=rp), intent(inout) :: tol
     integer :: report
     integer :: i, n

     n = field%dof%size()
    
     !> Make sure the work_fields are zeroed
     call rzero(work_field%x,n)
     call device_rzero(work_field%x_d,n)

     !> Copy the contents of the gpu to the work field
     call device_copy(work_field%x_d , field%x_d, n)
    
     !> Copy the information from the gpu to the cpu
     call device_memcpy(work_field%x, work_field%x_d, n, &
                        DEVICE_TO_HOST, sync=.true.)

     !> Now assert
     report = 0
     do i=1, n
        if ( work_field%x(i,1,1,1) .lt. tol ) then
           report = report + 1
        end if
     end do
  end function assert_gpu_count_entries_less_than_tol

  subroutine report_status_of_field(field, work_field)
     type(field_t), intent(inout) :: field
     type(field_t), intent(inout) :: work_field
     logical :: sync
     integer :: count_cpu, count_gpu
     real(kind=rp) tol

     !logical :: assert_gpu_cpu_synch
     !integer :: assert_cpu_count_entries_less_than_tol
     !integer :: assert_gpu_count_entries_less_than_tol

     tol = 1e-8
     sync = assert_gpu_cpu_synch(field, work_field, tol)
     tol = -1e-7
     count_cpu = assert_cpu_count_entries_less_than_tol(field, work_field, tol)
     count_gpu = assert_cpu_count_entries_less_than_tol(field, work_field, tol)

     if (sync) then
        write(*,*) "The field is syncrhonized in cpu/gpu"
     else
        write(*,*) "The field is NOT syncrhonized in cpu/gpu"
     end if
       
     write(*,*) "The number of entries with tolerance less than ", tol, "is:"
     write(*,*) "cpu: ", count_cpu
     write(*,*) "gpu: ", count_gpu
     write(*,*) "The number of entries in the field is ", field%dof%size()
     !write(*,*) "---------------------------------------------------------------"

  end subroutine report_status_of_field

  subroutine sync_field(field)
    type(field_t), intent(inout) :: field
    integer :: n

    n = field%dof%size()
    
    if (NEKO_BCKND_DEVICE .eq. 1) then 
       call device_memcpy(field%x, &
                          field%x_d, &
                          n,DEVICE_TO_HOST,sync=.true.)
    end if

  end subroutine sync_field
  
  subroutine sync_field_cpu_to_gpu(field)
    type(field_t), intent(inout) :: field
    integer :: n

    n = field%dof%size()
    
    if (NEKO_BCKND_DEVICE .eq. 1) then 
       call device_memcpy(field%x, &
                          field%x_d, &
                          n,HOST_TO_DEVICE,sync=.true.)
    end if

  end subroutine sync_field_cpu_to_gpu
  
  subroutine lcsum(a, b, lx, nelv)
     integer, intent(in) :: lx
     integer, intent(in) :: nelv
     real(kind=rp), intent(inout) :: a(nelv)
     real(kind=rp), intent(inout) :: b(lx,lx,lx,nelv)
     real(kind=rp) :: suma
     integer :: i, e
        
     !> Get local volumes
     do e = 1, nelv
        suma = 0_rp
        do i = 1, lx**3
           suma = suma + b(i,1,1,e) 
        end do
        a(e) = suma
     end do
  end subroutine lcsum


  subroutine check_mesh_holes(coef, dof, msh, Xh, gs) 
     type(coef_t), intent(inout) :: coef
     type(dofmap_t), intent(inout) :: dof
     type(mesh_t), intent(inout) :: msh
     type(space_t), intent(inout) :: Xh
     type(gs_t), intent(inout) :: gs

     type(field_t), target :: x, y, z
     type(field_t), target :: x_cpy, y_cpy, z_cpy
     type(file_t) :: file_obj
     type(field_list_t) :: field_list
     real(kind=rp) :: time, avrg_x, avrg_y, avrg_z
     
     integer :: i,j,k,n,e
     
     !> Initialize fields that will hold the data 

     call x%init(dof, "x")
     call y%init(dof, "y")
     call z%init(dof, "z")
     call x_cpy%init(dof, "x_cpy")
     call y_cpy%init(dof, "y_cpy")
     call z_cpy%init(dof, "z_cpy")

     ! Copy the data to the fields
     do e = 1, msh%nelv
        do i=1, Xh%lx
           do j=1, Xh%ly
              do k=1, Xh%lz
                 x%x(i,j,k,e) = dof%x(i,j,k,e)
                 y%x(i,j,k,e) = dof%y(i,j,k,e)
                 z%x(i,j,k,e) = dof%z(i,j,k,e)
                 x_cpy%x(i,j,k,e) = dof%x(i,j,k,e)
                 y_cpy%x(i,j,k,e) = dof%y(i,j,k,e)
                 z_cpy%x(i,j,k,e) = dof%z(i,j,k,e)
              end do
           end do
        end do
     end do
     
     call sync_field_cpu_to_gpu(x)
     call sync_field_cpu_to_gpu(y)
     call sync_field_cpu_to_gpu(z)
     call sync_field_cpu_to_gpu(x_cpy)
     call sync_field_cpu_to_gpu(y_cpy)
     call sync_field_cpu_to_gpu(z_cpy)
    
     !> Perform dssum
     call gs%op(x%x, x%dof%size(), GS_OP_ADD)
     call gs%op(y%x, y%dof%size(), GS_OP_ADD)
     call gs%op(z%x, z%dof%size(), GS_OP_ADD)
     if (NEKO_BCKND_DEVICE .eq. 1) then
        call device_col2(x%x_d, coef%mult_d, x%dof%size())
        call device_col2(y%x_d, coef%mult_d, y%dof%size())
        call device_col2(z%x_d, coef%mult_d, z%dof%size())
     else
        call col2(x%x, coef%mult, x%dof%size())
        call col2(y%x, coef%mult, y%dof%size())
        call col2(z%x, coef%mult, z%dof%size())
     end if 
    
     !> Substract arrays to see if they are the same
     if (NEKO_BCKND_DEVICE .eq. 1) then 
        call device_sub2(x%x_d, x_cpy%x_d, x%dof%size())
        call device_sub2(y%x_d, y_cpy%x_d, y%dof%size())
        call device_sub2(z%x_d, z_cpy%x_d, z%dof%size())
     else
        call sub2(x%x, x_cpy%x, x%dof%size())
        call sub2(y%x, y_cpy%x, y%dof%size())
        call sub2(z%x, z_cpy%x, z%dof%size())
     end if

     !> Syncrhonize the fields with the cpu
     call sync_field(x)
     call sync_field(y)
     call sync_field(z)

     !> Get the sums. These are not averages
     avrg_x = 0_rp
     avrg_y = 0_rp
     avrg_z = 0_rp
     if (NEKO_BCKND_DEVICE .eq. 1) then 
        avrg_x = device_glsum(x%x_d,x%dof%size())                  
        avrg_y = device_glsum(y%x_d,y%dof%size())               
        avrg_z = device_glsum(z%x_d,z%dof%size())                 
     else
        avrg_x = glsum(x%x,x%dof%size())                  
        avrg_y = glsum(y%x,y%dof%size())               
        avrg_z = glsum(z%x,z%dof%size())                 
     end if
     
     if (pe_rank.eq.0) then
             write(*,*) "The averages are: ", avrg_x, avrg_y, avrg_z
     end if

     !> ----------------------------------------------------------------
     !> Write the fields
     file_obj = file_t("mesh_holes_file.fld")
     allocate(field_list%fields(3))
     field_list%fields(1)%f => x
     field_list%fields(2)%f => y
     field_list%fields(3)%f => z 
     time = 0_rp
     call file_obj%write(field_list,time)
     deallocate(field_list%fields)
     call file_free(file_obj)

  end subroutine check_mesh_holes




  subroutine compare_mesh_kolmogorov(eps_t, eps_k, work_field, work_field2, coef, dof, msh, Xh, ra, pr) 
     type(coef_t), intent(inout) :: coef
     type(dofmap_t), intent(inout) :: dof
     type(mesh_t), intent(inout) :: msh
     type(space_t), intent(inout) :: Xh
     type(field_t), intent(inout) :: eps_t, eps_k
     type(field_t), intent(inout) :: work_field
     type(field_t), intent(inout) :: work_field2
     real(kind=rp), intent(inout) :: ra
     real(kind=rp), intent(inout) :: pr


     type(field_t), pointer :: dx, dy, dz, dr, ds, dt
     type(field_t), pointer :: dx_e, dy_e, dz_e
     type(file_t) :: file_obj
     type(field_list_t) :: field_list
     !
     real(kind=rp) :: suma
     real(kind=rp) :: lambda
     real(kind=rp) :: mu, bar_eps_t, bar_eps_k, eta_t, eta_k
     real(kind=rp) :: time
     integer :: i,j,k,n,e

     ! local averages
     real(kind=rp) :: dx_e_nelv(msh%nelv) 
     real(kind=rp) :: dy_e_nelv(msh%nelv) 
     real(kind=rp) :: dz_e_nelv(msh%nelv) 
     real(kind=rp) :: vol_e(msh%nelv) 
     type(c_ptr) :: dx_e_nelv_d = C_NULL_PTR
     type(c_ptr) :: dy_e_nelv_d = C_NULL_PTR
     type(c_ptr) :: dz_e_nelv_d = C_NULL_PTR
     type(c_ptr) :: vol_e_d = C_NULL_PTR
     
     !> Get more non dimensional parameters
     mu = sqrt(pr/ra)
     lambda = 1_rp/sqrt(ra*pr)
     
     !> ----------------------------------------------------------------
     !> get kolmogorov scale
     
     !> Integrate thermal dissipation
     bar_eps_t = 0_rp
     call average_from_weights(bar_eps_t, eps_t, &
                     work_field, coef%B, coef%volume)
  
     !> Integrate kinetic dissipation
     bar_eps_k = 0_rp
     call average_from_weights(bar_eps_k, eps_k, &
                     work_field, coef%B, coef%volume)
   
     !> Calculate global kolmogorov scale
     eta_k = 0_rp
     eta_k = ((mu**3_rp)/(bar_eps_k))**(0.25_rp)
    
     !> Calculate global batchelor scale
     eta_t = 0_rp
     eta_t = ((lambda**3_rp)/(bar_eps_k))**(0.25_rp)

     !> ----------------------------------------------------------------
     !> Calculate the grid spacing in the xyz

     call neko_field_registry%add_field(dof, "dx")
     call neko_field_registry%add_field(dof, "dy")
     call neko_field_registry%add_field(dof, "dz")
     call neko_field_registry%add_field(dof, "dx_e")
     call neko_field_registry%add_field(dof, "dy_e")
     call neko_field_registry%add_field(dof, "dz_e")
     call neko_field_registry%add_field(dof, "dr")
     call neko_field_registry%add_field(dof, "ds")
     call neko_field_registry%add_field(dof, "dt")
     dx => neko_field_registry%get_field("dx")
     dy => neko_field_registry%get_field("dy")
     dz => neko_field_registry%get_field("dz") 
     dx_e => neko_field_registry%get_field("dx_e")
     dy_e => neko_field_registry%get_field("dy_e")
     dz_e => neko_field_registry%get_field("dz_e") 
     dr => neko_field_registry%get_field("dr")
     ds => neko_field_registry%get_field("ds")
     dt => neko_field_registry%get_field("dt")

     ! get drst
     do e = 1, msh%nelv
        do i=1, Xh%lx
           do j=1, Xh%ly
              do k=1, Xh%lz
                 dr%x(i,j,k,e) = 1/Xh%dr_inv(i)
                 ds%x(i,j,k,e) = 1/Xh%ds_inv(j)
                 dt%x(i,j,k,e) = 1/Xh%dt_inv(k)
              end do
           end do
        end do
     end do
     
     call sync_field_cpu_to_gpu(dr)
     call sync_field_cpu_to_gpu(ds)
     call sync_field_cpu_to_gpu(dt)
     
     !!> get dx**2 = (dx/dr * dr)**2 + (dx/ds * ds)**2
     if (NEKO_BCKND_DEVICE .eq. 1) then
        call device_rzero(dx%x_d, dx%dof%size())
        call device_col3(work_field%x_d, coef%dxdr_d, dr%x_d, dx%dof%size())
        call device_addcol3(dx%x_d, work_field%x_d, work_field%x_d, dx%dof%size())
        call device_col3(work_field%x_d, coef%dxds_d, ds%x_d, dx%dof%size())
        call device_addcol3(dx%x_d, work_field%x_d, work_field%x_d, dx%dof%size())

        call device_rzero(dy%x_d, dx%dof%size())
        call device_col3(work_field%x_d, coef%dydr_d, dr%x_d, dx%dof%size())
        call device_addcol3(dy%x_d, work_field%x_d, work_field%x_d, dx%dof%size())
        call device_col3(work_field%x_d, coef%dyds_d, ds%x_d, dx%dof%size())
        call device_addcol3(dy%x_d, work_field%x_d, work_field%x_d, dx%dof%size())
        
        call device_rzero(dz%x_d, dx%dof%size())
        call device_col3(work_field%x_d, coef%dzdt_d, dt%x_d, dx%dof%size())
        call device_addcol3(dz%x_d, work_field%x_d, work_field%x_d, dx%dof%size())
     else
        call rzero(dx%x, dx%dof%size())
        call col3(work_field%x, coef%dxdr, dr%x, dx%dof%size())
        call addcol3(dx%x, work_field%x, work_field%x, dx%dof%size())
        call col3(work_field%x, coef%dxds, ds%x, dx%dof%size())
        call addcol3(dx%x, work_field%x, work_field%x, dx%dof%size())

        call rzero(dy%x, dx%dof%size())
        call col3(work_field%x, coef%dydr, dr%x, dx%dof%size())
        call addcol3(dy%x, work_field%x, work_field%x, dx%dof%size())
        call col3(work_field%x, coef%dyds, ds%x, dx%dof%size())
        call addcol3(dy%x, work_field%x, work_field%x, dx%dof%size())
        
        call rzero(dz%x, dx%dof%size())
        call col3(work_field%x, coef%dzdt, dt%x, dx%dof%size())
        call addcol3(dz%x, work_field%x, work_field%x, dx%dof%size())
     end if
    
   
     !> Syncrhonize the fields with the cpu
     call sync_field(dx)
     call sync_field(dy)
     call sync_field(dz)

     !> Take the square root and scale in the cpu
     do i = 1, dx%dof%size()
        dx%x(i,1,1,1) = sqrt(dx%x(i,1,1,1))/eta_k
        dy%x(i,1,1,1) = sqrt(dy%x(i,1,1,1))/eta_k
        dz%x(i,1,1,1) = sqrt(dz%x(i,1,1,1))/eta_k
     end do

     !> ----------------------------------------------------------------
     !> Write the fields
     file_obj = file_t("mesh_eta.fld")
     allocate(field_list%fields(3))
     field_list%fields(1)%f => dx
     field_list%fields(2)%f => dy
     field_list%fields(3)%f => dz 
     time = 0_rp
     call file_obj%write(field_list,time)
     deallocate(field_list%fields)
     call file_free(file_obj)


     !> ----------------------------------------------------------------
     !> Get the averages per element to compare to spectral error indicator
     
     !> Syncrhonize the fields with the cpu
     call sync_field_cpu_to_gpu(dx)
     call sync_field_cpu_to_gpu(dy)
     call sync_field_cpu_to_gpu(dz)

     if (NEKO_BCKND_DEVICE .eq. 1) then
        ! Map the pointers
        call device_map(dx_e_nelv,  dx_e_nelv_d,  msh%nelv)
        call device_map(dy_e_nelv,  dy_e_nelv_d,  msh%nelv)
        call device_map(dz_e_nelv,  dz_e_nelv_d,  msh%nelv)
        call device_map(vol_e,  vol_e_d,  msh%nelv)


        !> Get local volumes
        call device_lcsum(vol_e_d,coef%B_d, &
                            Xh%lx,msh%nelv)
        !! put it in the cpu
        call device_memcpy(vol_e, vol_e_d, msh%nelv, &
                           DEVICE_TO_HOST, sync=.true.)

        !> Get the weighted spacings
        call device_col3(dx_e%x_d, dx%x_d, coef%B_d, dx%dof%size())
        call device_col3(dy_e%x_d, dy%x_d, coef%B_d, dx%dof%size())
        call device_col3(dz_e%x_d, dz%x_d, coef%B_d, dx%dof%size())

        !> Now get the local integrals
        call device_lcsum(dx_e_nelv_d, dx_e%x_d, Xh%lx, msh%nelv)
        call device_lcsum(dy_e_nelv_d, dy_e%x_d, Xh%lx, msh%nelv)
        call device_lcsum(dz_e_nelv_d, dz_e%x_d, Xh%lx, msh%nelv)
        !! put it in the cpu
        call device_memcpy(dx_e_nelv, dx_e_nelv_d, msh%nelv, &
                           DEVICE_TO_HOST, sync=.true.)
        call device_memcpy(dy_e_nelv, dy_e_nelv_d, msh%nelv, &
                           DEVICE_TO_HOST, sync=.true.)
        call device_memcpy(dz_e_nelv, dz_e_nelv_d, msh%nelv, &
                           DEVICE_TO_HOST, sync=.true.)                   

     else
        
        !> Get local volumes
        call lcsum(vol_e,coef%B, &
                            Xh%lx,msh%nelv)
        
        !> Get the weighted spacings
        call col3(dx_e%x, dx%x, coef%B, dx%dof%size())
        call col3(dy_e%x, dy%x, coef%B, dx%dof%size())
        call col3(dz_e%x, dz%x, coef%B, dx%dof%size())
        
        !> Get the local integrals
        call lcsum(dx_e_nelv, dx_e%x, Xh%lx, msh%nelv)
        call lcsum(dy_e_nelv, dy_e%x, Xh%lx, msh%nelv)
        call lcsum(dz_e_nelv, dz_e%x, Xh%lx, msh%nelv)

     end if

     !> Being on the cpu, perform the local average operations
     do e = 1, msh%nelv
        dx_e_nelv(e) = dx_e_nelv(e)/vol_e(e)
        dy_e_nelv(e) = dy_e_nelv(e)/vol_e(e)
        dz_e_nelv(e) = dz_e_nelv(e)/vol_e(e)
     end do

     !> Put the averages in the field format to write it
     do e = 1, msh%nelv
        do i=1, Xh%lx
           do j=1, Xh%ly
              do k=1, Xh%lz
                 dx_e%x(i,j,k,e) = dx_e_nelv(e)
                 dy_e%x(i,j,k,e) = dy_e_nelv(e)
                 dz_e%x(i,j,k,e) = dz_e_nelv(e)
              end do
           end do
        end do
     end do

     !> Write the averages
     file_obj = file_t("mesh_eta_avrg.fld")
     allocate(field_list%fields(3))
     field_list%fields(1)%f => dx_e
     field_list%fields(2)%f => dy_e
     field_list%fields(3)%f => dz_e
     time = 0_rp
     call file_obj%write(field_list,time)
     deallocate(field_list%fields)
     call file_free(file_obj)

  end subroutine compare_mesh_kolmogorov
 
 !> @params
 !! dir is the wall normal direction 
  subroutine calculate_shear_stress(tau, u, coef, dir, mu)
    type(field_t), intent(inout) :: tau
    type(field_t), intent(inout) :: u
    type(coef_t), intent(inout) :: coef
    character(len=1), intent(in)  :: dir
    real(kind=rp), intent(inout)  :: mu
    integer :: n
    
    n = u%dof%size()
   
    !> Get the derivative in the relevant direction 
    !> dudxyz detects if there is a GPU pointer
    if (trim(dir) .eq. "x") then
       call dudxyz (tau%x, u%x, coef%drdx, coef%dsdx, coef%dtdx, coef)
    else if (trim(dir) .eq. "y") then
       call dudxyz (tau%x, u%x, coef%drdy, coef%dsdy, coef%dtdy, coef)
    else if (trim(dir) .eq. "z") then
       call dudxyz (tau%x, u%x, coef%drdz, coef%dsdz, coef%dtdz, coef)
    else
       write(*,*) "Direction not currently supported"
    end if
   
    !> Multiply by mu 
    if (NEKO_BCKND_DEVICE .eq. 1) then 
       call device_cmult(tau%x_d, mu, n)
    else
       call cmult(tau%x, mu, n)
    end if
     
  end subroutine calculate_shear_stress

end module post_processing_library
