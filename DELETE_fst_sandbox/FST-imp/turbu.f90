module turbu
    use neko
    use spec
    use global_params
    implicit none


contains

  !---------------------------------------------------------------------- 
    
  subroutine make_turbu(u)

    implicit none

    type(field_t), intent(in) :: u

    integer :: k,i,j
    integer :: shellno
    integer :: seed
    real(kind=rp) :: ue,ve,we
    real(kind=rp) :: uamp,vamp,wamp
    real(kind=rp) :: amp, bb1(fst_modes, 3), dlx, dly, dlz
    real(kind=rp) :: u_hat(fst_modes, 3), u_hat_p(fst_modes, 3)
    character(len=LOG_SIZE) :: log_buf

    ! call neko_log%message("calculating dlx")
    ! write(log_buf, *) "u%dof%x is of shape", shape(u%dof%x)
    ! call neko_log%message(log_buf)
    ! call print_int("u%xh%lx * u%xh%ly * u%xh%lz * u%msh%nelv",u%xh%lx * u%xh%ly * u%xh%lz * u%msh%nelv)
    dlx = glmax(u%dof%x,u%xh%lx * u%xh%ly * u%xh%lz * u%msh%nelv) - &
          glmin(u%dof%x,u%xh%lx * u%xh%ly * u%xh%lz * u%msh%nelv)
    ! call neko_log%message("done calculating dlx")

    ! call neko_log%message("calculating dly")
    dly = glmax(u%dof%y,u%xh%lx * u%xh%ly * u%xh%lz * u%msh%nelv) - &
          glmin(u%dof%y,u%xh%lx * u%xh%ly * u%xh%lz * u%msh%nelv)
    ! call neko_log%message("done calculating dly")

    ! call neko_log%message("calculating dlz")
    dlz = glmax(u%dof%z,u%xh%lx * u%xh%ly * u%xh%lz * u%msh%nelv) - &
          glmin(u%dof%z,u%xh%lx * u%xh%ly * u%xh%lz * u%msh%nelv)
    ! call neko_log%message("done calculating dlz")

          
    if ( pe_rank.eq.0 ) then
            
      seed = -143        

      call print_param("==== RAN", ran2(seed))
      call print_param("==== RAN", ran2(seed))
      call print_param("==== RAN", ran2(seed))
      call print_param("==== RAN", ran2(seed))
      call print_param("==== RAN", ran2(seed))
      ! call exit
      ! call print_param("==== PI" , pi)

      if (write_files) open(unit=137,form='formatted',file='bb.txt')
      
      call neko_log%message("calling spec_s")
      call spec_s(u, dlx, dly, dlz) ! get isotropically distributed wavenumbers in spheres 
      call neko_log%message("done spec_s")
      
      seed = -143

      do k=1,u%msh%gdim
        do i=1,fst_modes

        !write (*,*) "SEED=_", seed
                bb(i,k) = ran2(seed)*2.0*pi        ! random phase shift
        !write (*,*) "SEED=_", seed
                bb1(i,k) = 2.0*ran2(seed)-1.0      ! random amplitude
        !write (*,*) "SEED=_", seed

          if (write_files) write(137,*) bb(i,1), bb1(i,1)
          ! write(6,*) 'BB', bb(i,1)
        enddo
      enddo

      if (write_files) close(137)
      ! write(6,*) 'FST - Random amplitude generated'
      call neko_log%message("FST - Random amplitude generated")
      
      !     make sure that continuity is enforced by ensuring u_vec.k_vec=(0 0 0)
      do i=1,k_length
         do j= 1,u%msh%gdim
            u_hat(i,j)=bb1(i,j)
         enddo

         do j=1,u%msh%gdim
            u_hat_p(i,j) = u_hat(i,j) &
                 - (u_hat(i,1)*k_num_all(i,1) &
                 + u_hat(i,2)*k_num_all(i,2) &
                 + u_hat(i,3)*k_num_all(i,3)) &
                 * k_num_all(i,j) &
                 / (k_num_all(i,1)**2 &
                 +  k_num_all(i,2)**2 &
                 +  k_num_all(i,3)**2)
         enddo

         do j=1,u%msh%gdim
            u_hat_pn(i,j) = u_hat_p(i,j) &
                 / sqrt(u_hat_p(i,1)**2 &
                 + u_hat_p(i,2)**2 &
                 + u_hat_p(i,3)**2)
         enddo
      enddo

      call neko_log%message('FST - Amplitudes projection done')

      !           Check energy in individual modes
      ue=0.
      ve=0.
      we=0.
      !           Also write the modes
      if (write_files) open(file='fst_spectrum.csv',unit=13)
      if (write_files) write(13,'(9(A, ","),A)') 'ShellNo','kx','ky','kz', &
            'u_amp','v_amp','w_amp','u_hat_pn1','u_hat_pn2', 'u_hat_pn3'
      do i=1,k_length
        shellno = shell(i)
        amp = shell_amp(shellno)

        !write (*,*) "AMP: ", amp

        uamp = u_hat_pn(i,1)*amp
        vamp = u_hat_pn(i,2)*amp
        wamp = u_hat_pn(i,3)*amp
        
        if (write_files) write(13,'(9(g0, ","), g0)') shellno,k_num_all(i,1),k_num_all(i,2), &
        k_num_all(i,3),uamp,vamp,wamp, u_hat_pn(i,1), u_hat_pn(i,2), u_hat_pn(i,3)

        ue =  ue + ((uamp)**2)/2.
        ve =  ve + ((vamp)**2)/2.
        we =  we + ((wamp)**2)/2.
      enddo
      
      if (write_files) close(13)
      
      write(log_buf,'(A18,10x,E12.5E2)') 'FST - Energy in u',ue 
      call neko_log%message(log_buf)      
      write(log_buf,'(A18,10x,E12.5E2)') 'FST - Energy in v',ve 
      call neko_log%message(log_buf)      
      write(log_buf,'(A18,10x,E12.5E2)') 'FST - Energy in w',we 
      call neko_log%message(log_buf)      
      write(log_buf,'(A20,8x,E12.5E2)') 'FST - Estimated tke', &
                                    (ue+ve+we)/2.
      call neko_log%message(log_buf)                                    
      write(log_buf,'(A19,9x,E12.5E2)') 'FST - Estimated ti', &
      sqrt((ue+ve+we)/3.)
      call neko_log%message(log_buf)                                  
    
    end if
                                 
    return
  end subroutine make_turbu
  !----------------------------------------------------------------------       
  
end module turbu
