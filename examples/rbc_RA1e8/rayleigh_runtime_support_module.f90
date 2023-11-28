module rbc
  use neko     
  use time_based_controller, only : time_based_controller_t
  use mean_field, only : mean_field_t
  use, intrinsic :: iso_c_binding
  implicit none
  private

  type, public :: rbc_t 

     !> Store the simulation time ste
     integer :: tstep
     type(coef_t), pointer :: coef => null()

     !> Timers
     real(kind=dp) :: start_time, end_time
     
     !> Flow quantities
     type(field_t) :: work_field ! Field to perform operations
     type(field_t) :: work_field2 ! Field to perform operations
     type(field_t) :: tt ! T ^2
     type(field_t) :: uxt ! u_z * T
     type(field_t) :: uyt ! u_z * T
     type(field_t) :: uzt ! u_z * T
     type(field_t) :: dtdx ! Derivative of scalar wrt x
     type(field_t) :: dtdy ! Detivative of scalar wrt y
     type(field_t) :: dtdz ! Derivative of scalar wrt z
     type(field_t) :: dtdn ! Derivative of scalar wrt normal
     type(field_t) :: div_dtdX ! divergence of heat flux
     type(field_t) :: dudx ! Derivative wrt x
     type(field_t) :: dudy ! Detivative wrt y
     type(field_t) :: dudz ! Derivative wrt z
     type(field_t) :: dvdx ! Derivative wrt x
     type(field_t) :: dvdy ! Detivative wrt y
     type(field_t) :: dvdz ! Derivative wrt z
     type(field_t) :: dwdx ! Derivative wrt x
     type(field_t) :: dwdy ! Detivative wrt y
     type(field_t) :: dwdz ! Derivative wrt z
     type(field_t) :: eps_k ! kinetic dissipation
     type(field_t) :: eps_t ! thermal dissipation
     
     !> Mean fields     
     type(mean_field_t) :: mean_t     ! Derivative wrt z
     type(mean_field_t) :: mean_tt     ! Derivative wrt z
     type(mean_field_t) :: mean_uxt   ! Derivative wrt z
     type(mean_field_t) :: mean_uyt   ! Derivative wrt z
     type(mean_field_t) :: mean_uzt   ! Derivative wrt z
     type(mean_field_t) :: mean_dtdx     ! Derivative wrt z
     type(mean_field_t) :: mean_dtdy     ! Derivative wrt z
     type(mean_field_t) :: mean_dtdz     ! Derivative wrt z
     type(mean_field_t) :: mean_dtdn  ! Derivative wrt z
     type(mean_field_t) :: mean_eps_k ! Detivative wrt y
     type(mean_field_t) :: mean_eps_t ! Derivative wrt z
     type(mean_field_t) :: mean_speri_u  ! Derivative wrt z
     type(mean_field_t) :: mean_speri_v  ! Derivative wrt z
     type(mean_field_t) :: mean_speri_w  ! Derivative wrt z
     real(kind = rp):: t_last_sample = 0_rp
 
     !> Support fields for boundary/facet 
     type(field_t) :: mass_area_top ! mass matrix for area on top wall
     type(field_t) :: mass_area_bot ! mass matrix for area on bottom wall
     type(field_t) :: mass_area_side ! mass matrix for area on top wall
     type(field_t) :: bdry_mask ! mass matrix for area on top wall
     type(field_t) :: normal_x ! normal vector in x (non zero at the boundaries)
     type(field_t) :: normal_y ! normal vector in y (non zero at the boundaries)
     type(field_t) :: normal_z ! normal vector in z (non zero at the boundaries)
     type(field_t) :: bm1 ! mass matrix for area on bottom wall
  
     !> Real values from intigral quantities
     real(kind=rp) :: bar_uzt ! Volume integral
     real(kind=rp) :: bar_div_dtdX ! Volume integral
     real(kind=rp) :: top_wall_bar_dtdz ! area integral
     real(kind=rp) :: bot_wall_bar_dtdz ! area integral
     real(kind=rp) :: top_wall_area ! area of top and bottom wall
     real(kind=rp) :: bot_wall_area ! area of top and bottom wall
     real(kind=rp) :: side_wall_area ! area of top and bottom wall
     real(kind=rp) :: heat_bot_wall
     real(kind=rp) :: heat_top_wall
     real(kind=rp) :: heat_side_wall
     real(kind=rp) :: heat_balance
     real(kind=rp) :: nu_eps_k
     real(kind=rp) :: nu_eps_t
     real(kind=rp) :: bar_eps_k
     real(kind=rp) :: bar_eps_t
     real(kind=rp) :: eta_k !Globally averaged kolmogorov scale
     real(kind=rp) :: eta_t !Globally averaged Batchelor scale
     real(kind=rp) :: tke ! turbulent kinetic energy in space
     real(kind=rp) :: uuprime ! turbulent kinetic energy in x
     real(kind=rp) :: vvprime ! "" in y
     real(kind=rp) :: wwprime ! "" in z
  
     !> field list for outputs
     type(field_list_t) :: area_l
     type(field_list_t) :: dtdX_l
     type(field_list_t) :: eps_l
     type(field_list_t) :: mean_fields_l
     logical :: data_in_stats = .false.

     !> Variables to write extra files
     character(len=NEKO_FNAME_LEN) :: fname
  
     !> Files
     type(file_t) :: mf_dtdn
     type(file_t) :: mf_dtdX
     type(file_t) :: mf_eps
     type(file_t) :: mf_mean_fields
     type(file_t) :: mf_area
     type(file_t) :: mf_bm1

     !> List of elements and facets in uper and lower boundary
     type(stack_i4t4_t) :: top_wall_facet
     type(stack_i4t4_t) :: bot_wall_facet
     type(stack_i4t4_t) :: side_wall_facet
     type(stack_i4t4_t) :: wall_facet

     !> Spectral error indicator
     type(spectral_error_indicator_t) :: spectral_error_indicator
    
     !> Calculation controllers
     type(time_based_controller_t) :: sample_control
     type(time_based_controller_t) :: field_write_control
     type(time_based_controller_t) :: data_stream_control
        
     !> for mean
     real(kind = rp):: eps_t_mean = 0_rp

     !> for data stream
     type(data_streamer_t) :: dstream
     logical :: stream_data = .false.

   contains
     !> Constructor
     procedure, pass(this) :: init => rbc_init
     !> Destructor
     procedure, pass(this) :: free => rbc_free
     !> calculate quantities
     procedure, pass(this) :: calculate => rbc_calculate
     !> Sync the values in the cpu to those of GPU
     procedure, pass(this) :: sync => rbc_sync_fromGPU
     !> Sync the values in the cpu to those of GPU
     procedure, pass(this) :: write_fields => rbc_write_fields
     !> Sync the values in the cpu to those of GPU
     procedure, pass(this) :: write_stats => rbc_write_stats
     !> Sync the values in the cpu to those of GPU
     procedure, pass(this) :: update_stats => rbc_update_stats
     !> Integrate relevant quantities
     procedure, pass(this) :: get_integral_quantities => rbc_get_integral_quantities


  end type rbc_t

contains

  subroutine rbc_init(this, t, u, v, w, p, coef, params)
    class(rbc_t), intent(inout), target :: this
    real(kind=rp) :: t
    type(field_t), intent(inout) :: u
    type(field_t), intent(inout) :: v
    type(field_t), intent(inout) :: w
    type(field_t), intent(inout) :: p
    type(coef_t), intent(inout), target :: coef
    type(json_file), intent(inout) :: params
    
    ! Local variables 
    character(len=NEKO_FNAME_LEN) :: fname
    type(field_t), pointer :: s

    !> Parameters to populate the list of elements and facets
    real(kind=rp) :: stack_size(1) !Do it like this so it is cast
    real(kind=rp) :: stack_size_global
    integer :: e, facet, typ, e_counter, facet_counter, i,j,k,n
    real(kind=rp) :: normal(3)
    real(kind=rp) :: area_mass, area_rank(1),area_global(1)
    type(tuple4_i4_t), pointer :: wall_facet_list(:)
    type(tuple4_i4_t), pointer :: sidewall_facet_list(:)
    real(kind=rp) :: voll(1), voll_temp(1)
    integer :: lx, ly, lz, nones
    real(kind=rp) :: ones(u%Xh%lx,u%Xh%ly,6,u%msh%nelv)
    type(tuple4_i4_t)  :: facet_el_type_0

    !> Controller data
    real(kind=rp)     :: write_par
    character(len=:), allocatable  :: write_control
    real(kind=rp)     :: T_end
    integer :: how_many_samples
    logical :: found

    !> Get te pointer to coef
    this%coef => coef

    s => neko_field_registry%get_field('s')
    
    !> Initialize the spectral error indicator
    call this%spectral_error_indicator%init(u,v,w,coef)
    
    !> Initialize variables related to nusselt calculation
    call this%work_field%init( u%dof, 'work_field')
    call this%work_field2%init( u%dof, 'work_field2')
    call this%tt%init( u%dof, 'tt')
    call this%uxt%init( u%dof, 'uxt')
    call this%uyt%init( u%dof, 'uyt')
    call this%uzt%init( u%dof, 'uzt')
    call this%dtdx%init( u%dof, 'dtdx')
    call this%dtdy%init( u%dof, 'dtdy')
    call this%dtdz%init( u%dof, 'dtdz')
    call this%dtdn%init( u%dof, 'dtdn')
    
    call this%dudx%init( u%dof, 'dudx')
    call this%dudy%init( u%dof, 'dudy')
    call this%dudz%init( u%dof, 'dudz')
    call this%dvdx%init( u%dof, 'dvdx')
    call this%dvdy%init( u%dof, 'dvdy')
    call this%dvdz%init( u%dof, 'dvdz')
    call this%dwdx%init( u%dof, 'dwdx')
    call this%dwdy%init( u%dof, 'dwdy')
    call this%dwdz%init( u%dof, 'dwdz')


    call this%eps_k%init( u%dof, 'eps_k')
    call this%eps_t%init( u%dof, 'eps_t')

    call this%div_dtdX%init( u%dof, 'div_dtdX')
    call this%mass_area_top%init( u%dof, 'mat')
    call this%mass_area_bot%init( u%dof, 'mab')
    call this%mass_area_side%init( u%dof, 'masd')
    call this%bdry_mask%init( u%dof, 'bdry_mask')
    call this%bm1%init( u%dof, 'mass_mat')
    call this%normal_x%init( u%dof, 'normal_x')
    call this%normal_y%init( u%dof, 'normal_y')
    call this%normal_z%init( u%dof, 'normal_z')
    
    !> Mean fields
    call this%mean_t%init( s, 'mean_t')
    call this%mean_tt%init( this%tt, 'mean_tt')
    call this%mean_uxt%init( this%uxt, 'mean_uxt')
    call this%mean_uyt%init( this%uyt, 'mean_uyt')
    call this%mean_uzt%init( this%uzt, 'mean_uzt')
    call this%mean_dtdx%init( this%dtdx, 'mean_dtdx')
    call this%mean_dtdy%init( this%dtdy, 'mean_dtdy')
    call this%mean_dtdz%init( this%dtdz, 'mean_dtdz')
    call this%mean_dtdn%init( this%dtdn, 'mean_dtdn')
    call this%mean_eps_k%init( this%eps_k, 'mean_eps_k')
    call this%mean_eps_t%init( this%eps_t, 'mean_eps_t')
    call this%mean_speri_u%init( this%spectral_error_indicator%u_hat, 'mean_speri_u')
    call this%mean_speri_v%init( this%spectral_error_indicator%v_hat, 'mean_speri_v')
    call this%mean_speri_w%init( this%spectral_error_indicator%w_hat, 'mean_speri_w')
    
    !> Initialize the controllers
    call json_get(params, 'case.end_time', T_end)
    !! Sample controller
    call json_get(params, 'case.rbc_sampler.output_control', &
                                     write_control)
    call json_get(params, 'case.rbc_sampler.calc_frequency', &
                                     write_par)    
    call this%sample_control%init(T_end, write_control, &
                                      write_par)
    if (this%sample_control%nsteps .eq. 0) then
       this%sample_control%nexecutions = &
               int(t / this%sample_control%time_interval) + 1
    end if
    !! Field writer controller
    call json_get(params, 'case.rbc_field_writer.output_control', &
                                     write_control)
    call json_get(params, 'case.rbc_field_writer.calc_frequency', &
                                     write_par)    
    call this%field_write_control%init(T_end, write_control, &
                                      write_par)
    if (this%field_write_control%nsteps .eq. 0) then
       this%field_write_control%nexecutions = &
               int(t / this%field_write_control%time_interval) + 1
    end if 
    !! Field writer controller
    call json_get(params, 'case.rbc_field_writer.output_control', &
                                     write_control)
    call json_get(params, 'case.rbc_field_writer.calc_frequency', &
                                     write_par)    
    call this%field_write_control%init(T_end, write_control, &
                                      write_par)
    if (this%field_write_control%nsteps .eq. 0) then
       this%field_write_control%nexecutions = &
               int(t / this%field_write_control%time_interval) + 1
    end if
    !! Data stream controller
    call json_get_or_default(params, 'case.rbc_data_streamer.stream_data',&
                                   found, .false.)
    if (found .eqv. .true.) then
       call json_get(params, 'case.rbc_data_streamer.output_control', &
                                     write_control)
       call json_get(params, 'case.rbc_data_streamer.calc_frequency', &
                                     write_par)    
       call this%data_stream_control%init(T_end, write_control, &
                                      write_par)
       if (this%data_stream_control%nsteps .eq. 0) then
          this%data_stream_control%nexecutions = &
                  int(t / this%data_stream_control%time_interval) + 1
       end if
    end if



    !> Initialize the sampling interval for stats
    if (t .gt. 0) then
       if (pe_rank .eq. 0) then 
          write(*,*) 'restarting the value for last stat sample '
       end if
       !if (this%sample_control%nsteps .eq. 0) then 
       !   !> find how many samples have been taken
       !   how_many_samples = floor(t/this%sample_control%time_interval)
       !   !> Update the time since the last sample based on the information
       !   this%t_last_sample = how_many_samples * &
       !                     this%sample_control%time_interval
       !else if (this%sample_control%nsteps .gt. 0) then
       !   this%t_last_sample = t
       !end if
       
       if (pe_rank .eq. 0) then 
          write(*,*) 'by construction last user stat sample was last time step of previous run'
       end if
       this%t_last_sample = t

       if (pe_rank .eq. 0) then 
          write(*,*) 'last sample was at t= ', this%t_last_sample
       end if

    end if

    !> Initialize the files
    fname = 'dtdX.fld'
    this%mf_dtdX =  file_t(fname)
    call this%mf_dtdX%set_counter(this%field_write_control%nexecutions)
    call this%mf_dtdX%set_start_counter(this%field_write_control%nexecutions)

    fname = 'eps.fld'
    this%mf_eps =  file_t(fname)
    call this%mf_eps%set_counter(this%field_write_control%nexecutions)
    call this%mf_eps%set_start_counter(this%field_write_control%nexecutions)
    
    fname = 'mean_fields_usr.fld'
    this%mf_mean_fields =  file_t(fname)
    call this%mf_mean_fields%set_counter(this%field_write_control%nexecutions)
    call this%mf_mean_fields%set_start_counter(this%field_write_control%nexecutions)
    
    fname = 'dtdn.fld'
    this%mf_dtdn =  file_t(fname)
    call this%mf_dtdn%set_counter(this%field_write_control%nexecutions)
    call this%mf_dtdn%set_start_counter(this%field_write_control%nexecutions)

    fname = 'area.fld'
    this%mf_area =  file_t(fname)
    fname = 'bm1.fld'
    this%mf_bm1 =  file_t(fname)


    !> Initialize the lists
    allocate(this%area_l%fields(3))
    this%area_l%fields(1)%f => this%mass_area_top 
    this%area_l%fields(2)%f => this%mass_area_bot
    this%area_l%fields(3)%f => this%mass_area_side 

    allocate(this%dtdX_l%fields(3))
    this%dtdX_l%fields(1)%f => this%dtdx 
    this%dtdX_l%fields(2)%f => this%dtdy 
    this%dtdX_l%fields(3)%f => this%dtdz

    allocate(this%eps_l%fields(2))
    this%eps_l%fields(1)%f => this%eps_t
    this%eps_l%fields(2)%f => this%eps_k

    allocate(this%mean_fields_l%fields(14))
    this%mean_fields_l%fields(1)%f => this%mean_t%mf
    this%mean_fields_l%fields(2)%f => this%mean_tt%mf
    this%mean_fields_l%fields(3)%f => this%mean_uxt%mf
    this%mean_fields_l%fields(4)%f => this%mean_uyt%mf
    this%mean_fields_l%fields(5)%f => this%mean_uzt%mf
    this%mean_fields_l%fields(6)%f => this%mean_dtdx%mf
    this%mean_fields_l%fields(7)%f => this%mean_dtdy%mf
    this%mean_fields_l%fields(8)%f => this%mean_dtdz%mf
    this%mean_fields_l%fields(9)%f => this%mean_dtdn%mf
    this%mean_fields_l%fields(10)%f => this%mean_eps_t%mf
    this%mean_fields_l%fields(11)%f => this%mean_eps_k%mf
    this%mean_fields_l%fields(12)%f => this%mean_speri_u%mf
    this%mean_fields_l%fields(13)%f => this%mean_speri_v%mf
    this%mean_fields_l%fields(14)%f => this%mean_speri_w%mf


    !> Initialize list to identify relevant facets in boundaries
    call this%wall_facet%init()
    call this%top_wall_facet%init()
    call this%bot_wall_facet%init()
    call this%side_wall_facet%init()

    !> Populate the list with upper and lower and wall facets 
    do e = 1, u%msh%nelv !Go over all elements
      do facet = 1, 6 ! Go over all facets of hex element
        normal = coef%get_normal(1,1,1,e,facet) ! Get facet normal
        typ =  u%msh%facet_type(facet,e)             ! Get facet type
        if (typ.ne.0) then !then it is a boundary facet
          if ((normal(3)).ge.0.999) then !then it is on the top plate
            facet_el_type_0%x = (/facet, e, typ, 0/)
            call this%wall_facet%push(facet_el_type_0)
            call this%top_wall_facet%push(facet_el_type_0)
          else if ((normal(3)).le.(-0.999)) then !then it is on the bottom
            facet_el_type_0%x = (/facet, e, typ, 0/)
            call this%wall_facet%push(facet_el_type_0)
            call this%bot_wall_facet%push(facet_el_type_0)
          else !then it is on the sidewall
            facet_el_type_0%x = (/facet, e, typ, 0/)
            call this%wall_facet%push(facet_el_type_0)
            call this%side_wall_facet%push(facet_el_type_0)
          end if
        end if
      end do
    end do
    !! Determine the number of facets in the walls
    stack_size = real(this%wall_facet%top_,kind=rp)
    stack_size_global = glsum(stack_size,1)
    if (pe_rank .eq. 0) then
       write(*,*) 'Facets at the wall in this rank: ', stack_size
       write(*,*) 'Facets at the wall global: ', stack_size_global
    end if

    ! Map from facet data to field data to have an easy way to perform operations with col3
    !! Area mass matrices
    call map_from_facets_to_field(this%mass_area_top, coef%area, this%top_wall_facet)
    call map_from_facets_to_field(this%mass_area_bot, coef%area, this%bot_wall_facet)
    call map_from_facets_to_field(this%mass_area_side, coef%area, this%side_wall_facet)
    !! Normal vectors in the required facets
    call map_from_facets_to_field(this%normal_x, coef%nx, this%wall_facet)
    call map_from_facets_to_field(this%normal_y, coef%ny, this%wall_facet)
    call map_from_facets_to_field(this%normal_z, coef%nz, this%wall_facet)
    
    !! Create a mask with the boundary facet data
    n = size(coef%B)
    nones = u%Xh%lx*u%Xh%ly*6*u%msh%nelv
    call rone(ones, nones)
    call map_from_facets_to_field(this%bdry_mask, ones, this%wall_facet)

    this%top_wall_area = glsum(this%mass_area_top%x,n)
    this%bot_wall_area = glsum(this%mass_area_bot%x,n)
    this%side_wall_area = glsum(this%mass_area_side%x,n)

    !! Verify validity of averaging routines
    this%bar_uzt = 0_rp
    this%top_wall_bar_dtdz = 0_rp
    this%bot_wall_bar_dtdz = 0_rp
    call average_from_weights(this%bar_uzt, this%bdry_mask, &
                        this%work_field, this%mass_area_side%x, this%side_wall_area)
    call average_from_weights(this%top_wall_bar_dtdz, this%bdry_mask, &
                        this%work_field, this%mass_area_top%x, this%top_wall_area)
    call average_from_weights(this%bot_wall_bar_dtdz, this%bdry_mask, &
                        this%work_field, this%mass_area_bot%x, this%bot_wall_area)
    !! Write it out
    if (pe_rank .eq. 0) then
       write(*,*) 'Area of top wall * normal= ', this%top_wall_area
       write(*,*) 'Area of bottom wall * normal= ', this%bot_wall_area
       write(*,*) 'Area of side wall= ', this%side_wall_area
       write(*,*) 'Validity check of averaging #1 (should be 1): ', this%bar_uzt 
       write(*,*) 'Validity check of averaging #2 (should be 1): ', this%top_wall_bar_dtdz 
       write(*,*) 'Validity check of averaging #3 (should be 1): ', this%bot_wall_bar_dtdz 
    end if 

    !> Store volume mass matrix for post processing
    n = size(coef%B) 
    call rzero(this%bm1%x,u%dof%size())
    call copy(this%bm1%x,coef%B,n)
    voll = glsum(this%bm1%x,n)
    if (pe_rank .eq. 0) then
       write(*,*) 'volume=', voll
    end if 

    !> Perform IO
    call this%mf_area%write(this%area_l,t)
    call this%mf_bm1%write(this%bm1,t)      



    !> Initialize data streamer
    call json_get_or_default(params, 'case.rbc_data_streamer.stream_data',&
                                   found, .false.)
    if (found .eqv. .true.) then
       !> Initialize the streamer object
       call this%dstream%init(coef, 1)
       this%stream_data = .true.
    end if



  end subroutine rbc_init


  subroutine rbc_free(this)
    class(rbc_t), intent(inout), target :: this

    ! Finalize variables related to nusselt calculation
    call this%work_field%free()
    call this%work_field2%free()
    call this%tt%free()
    call this%uxt%free()
    call this%uyt%free()
    call this%uzt%free()
    call this%dtdx%free()
    call this%dtdy%free()
    call this%dtdz%free()
    call this%dtdn%free()
    call this%div_dtdX%free()

    call this%dudx%free()
    call this%dudy%free()
    call this%dudz%free()
    call this%dvdx%free()
    call this%dvdy%free()
    call this%dvdz%free()
    call this%dwdx%free()
    call this%dwdy%free()
    call this%dwdz%free()

    call this%eps_k%free()
    call this%eps_t%free()
    
    !> Mean fields
    call this%mean_t%free()
    call this%mean_tt%free()
    call this%mean_uxt%free()
    call this%mean_uyt%free()
    call this%mean_uzt%free()
    call this%mean_dtdx%free()
    call this%mean_dtdy%free()
    call this%mean_dtdz%free()
    call this%mean_dtdn%free()
    call this%mean_eps_k%free()
    call this%mean_eps_t%free()
    call this%mean_speri_u%free()
    call this%mean_speri_v%free()
    call this%mean_speri_w%free()
    
    call this%mass_area_top%free()
    call this%mass_area_bot%free()
    call this%mass_area_side%free()
    call this%bdry_mask%free()
    call this%bm1%free()
    call this%normal_x%free()
    call this%normal_y%free()
    call this%normal_z%free()
    
    ! Finilize list that contains uper and lower wall facets
    call this%wall_facet%free()
    call this%top_wall_facet%free()
    call this%bot_wall_facet%free()
    call this%side_wall_facet%free()
 
    ! Finilize the field lists
    if (allocated(this%area_l%fields)) then
       deallocate(this%area_l%fields)       
    end if
    if (allocated(this%dtdX_l%fields)) then
       deallocate(this%dtdX_l%fields)       
    end if
    if (allocated(this%eps_l%fields)) then
       deallocate(this%eps_l%fields)       
    end if
    if (allocated(this%mean_fields_l%fields)) then
       deallocate(this%mean_fields_l%fields)       
    end if
    
    !> Finalize the spectral error indicator
    call this%spectral_error_indicator%free()

    if (this%stream_data .eqv. .true.) then
       !> finalize the streamer object
       call this%dstream%free()
    end if


  end subroutine rbc_free
  
  subroutine rbc_write_fields(this, t)
    class(rbc_t), intent(inout), target :: this
    real(kind=rp) :: t

    !> Write instantaneous quantitiess
    call this%mf_eps%write(this%eps_l,t)
    call this%mf_dtdn%write(this%dtdn,t)
    call this%mf_dtdX%write(this%dtdX_l,t)
    call this%spectral_error_indicator%write_as_field(t)

  end subroutine rbc_write_fields

  subroutine rbc_write_stats(this, t)
    class(rbc_t), intent(inout), target :: this
    real(kind=rp) :: t
    
    !> Write and reset mean fields
    call this%mf_mean_fields%write(this%mean_fields_l,t) 
    call this%mean_t%reset() 
    call this%mean_tt%reset() 
    call this%mean_uxt%reset() 
    call this%mean_uyt%reset() 
    call this%mean_uzt%reset() 
    call this%mean_dtdx%reset() 
    call this%mean_dtdy%reset() 
    call this%mean_dtdz%reset() 
    call this%mean_dtdn%reset() 
    call this%mean_eps_k%reset()   
    call this%mean_eps_t%reset() 
    call this%mean_speri_u%reset() 
    call this%mean_speri_v%reset() 
    call this%mean_speri_w%reset()
    this%data_in_stats = .false. 

  end subroutine rbc_write_stats

 
  subroutine rbc_calculate(this, t, tstep, coef, params, Ra, Pr, get_spec_err_ind)
    class(rbc_t), intent(inout), target :: this
    real(kind=rp), intent(in) :: t
    integer, intent(in) :: tstep
    type(coef_t), intent(inout) :: coef
    type(json_file), intent(inout) :: params
    real(kind=rp), intent(inout) :: Ra
    real(kind=rp), intent(inout) :: Pr
    logical, intent(inout) :: get_spec_err_ind
    
    !> Input fields
    type(field_t), pointer :: u
    type(field_t), pointer :: v
    type(field_t), pointer :: w
    type(field_t), pointer :: p
    type(field_t), pointer :: s
    
    !> Misc parameters
    integer :: n,ntot, i,j,k, facet, e, lx,ly,lz
    integer :: index(4)
    type (tuple_i4_t) :: facet_el
    real(kind=rp) :: normal(3)
    real(kind=rp) :: voll(1), voll_temp(1)
    logical :: verify_bc = .false. ! write boundary conditions
    real(kind=rp) :: dt


    !> Get the needed fields
    u => neko_field_registry%get_field('u')
    v => neko_field_registry%get_field('v')
    w => neko_field_registry%get_field('w')
    p => neko_field_registry%get_field('p')
    s => neko_field_registry%get_field('s')

    !> Misc parameters
    n = size(coef%B)
    ntot = coef%dof%size()
    lx = u%Xh%lx
    ly = u%Xh%ly
    lz = u%Xh%lz

    if (pe_rank == 0) write(*,*) 'Performing user calculations of fields'
    
    ! ==================================================================
    ! ===================== "GPU computations" =========================
    ! ==================================================================
    
    this%start_time = MPI_WTIME()
    
    ! ======================== "Primitives" ============================

    !> Derivatives of temperature
    call dudxyz (this%dtdx%x, s%x, coef%drdx, coef%dsdx, coef%dtdx, coef)
    call dudxyz (this%dtdy%x, s%x, coef%drdy, coef%dsdy, coef%dtdy, coef)
    call dudxyz (this%dtdz%x, s%x, coef%drdz, coef%dsdz, coef%dtdz, coef)
    
    !> Derivatives of velocity
    call dudxyz (this%dudx%x, u%x, coef%drdx, coef%dsdx, coef%dtdx, coef)
    call dudxyz (this%dudy%x, u%x, coef%drdy, coef%dsdy, coef%dtdy, coef)
    call dudxyz (this%dudz%x, u%x, coef%drdz, coef%dsdz, coef%dtdz, coef) 
    call dudxyz (this%dvdx%x, v%x, coef%drdx, coef%dsdx, coef%dtdx, coef)
    call dudxyz (this%dvdy%x, v%x, coef%drdy, coef%dsdy, coef%dtdy, coef)
    call dudxyz (this%dvdz%x, v%x, coef%drdz, coef%dsdz, coef%dtdz, coef)
    call dudxyz (this%dwdx%x, w%x, coef%drdx, coef%dsdx, coef%dtdx, coef)
    call dudxyz (this%dwdy%x, w%x, coef%drdy, coef%dsdy, coef%dtdy, coef)
    call dudxyz (this%dwdz%x, w%x, coef%drdz, coef%dsdz, coef%dtdz, coef)


    !> Temperature squared
    if (NEKO_BCKND_DEVICE .eq. 1) then 
       call device_col3(this%tt%x_d,s%x_d,s%x_d,n)               
    else
       call col3(this%tt%x,s%x,s%x,n)                          
    end if

    
    ! ======================== "Derived" ============================
    
    !> Convective heat flux
    if (NEKO_BCKND_DEVICE .eq. 1) then 
       call device_col3(this%uxt%x_d,u%x_d,s%x_d,n)               
       call device_col3(this%uyt%x_d,v%x_d,s%x_d,n)               
       call device_col3(this%uzt%x_d,w%x_d,s%x_d,n)               
    else
       call col3(this%uxt%x,u%x,s%x,n)                          
       call col3(this%uyt%x,v%x,s%x,n)                          
       call col3(this%uzt%x,w%x,s%x,n)                          
    end if
    
    !> Derivatives of temperature normal to the boundaries
    call project_to_vector(this%dtdn, &
                           this%dtdx, this%dtdy, this%dtdz, &
                           this%normal_x%x, this%normal_y%x, this%normal_z%x)

    !> Thermal dissipation
    call calculate_thermal_dissipation(this%eps_t, this%dtdx,this%dtdy, &
                                       this%dtdz, this%work_field,Ra,Pr)
    
    !> Kinetic dissipation
    call calculate_kinetic_dissipation(this%eps_k, &
                                       this%dudx,this%dudy,this%dudz, &
                                       this%dvdx,this%dvdy,this%dvdz, &
                                       this%dwdx,this%dwdy,this%dwdz, &
                                       this%work_field,Ra,Pr)

    !> Divergence of the heat flux
    call divergence_of_field(this%div_dtdX, &
                             this%dtdx, this%dtdy, this%dtdz, &
                             this%work_field, coef)
   
    this%end_time = MPI_WTIME()
    if (pe_rank == 0) write(*,*) 'Elapsed time for calculate in GPU(s) = ', &
                                  this%end_time - this%start_time

    ! ==================================================================
    ! ===================== "CPU computations" =========================
    ! ==================================================================
    
    if (get_spec_err_ind) then

       this%start_time = MPI_WTIME()

       !> Get the spectral error indicators for the mesh
       call this%spectral_error_indicator%get_indicators(coef)
    
       this%end_time = MPI_WTIME()
       if (pe_rank == 0) write(*,*) 'Elapsed time for calculate in CPU(s) = ', &
                                  this%end_time - this%start_time
    end if

  end subroutine rbc_calculate


  subroutine rbc_update_stats(this, t)
    class(rbc_t), intent(inout), target :: this
    real(kind=rp), intent(in) :: t

    this%data_in_stats = .true. !> Indicate that there is data in buffer
    call this%mean_t%update(t-this%t_last_sample)
    call this%mean_tt%update(t-this%t_last_sample)
    call this%mean_uxt%update(t-this%t_last_sample)
    call this%mean_uyt%update(t-this%t_last_sample)
    call this%mean_uzt%update(t-this%t_last_sample)
    call this%mean_dtdx%update(t-this%t_last_sample)
    call this%mean_dtdy%update(t-this%t_last_sample)
    call this%mean_dtdz%update(t-this%t_last_sample)
    call this%mean_dtdn%update(t-this%t_last_sample)
    call this%mean_speri_u%update(t-this%t_last_sample)
    call this%mean_speri_v%update(t-this%t_last_sample)
    call this%mean_speri_w%update(t-this%t_last_sample)
    call this%mean_eps_k%update(t-this%t_last_sample)
    call this%mean_eps_t%update(t-this%t_last_sample)
    this%t_last_sample = t !< Update the last time that we updated

  end subroutine rbc_update_stats

  
  subroutine rbc_get_integral_quantities(this, t, tstep, coef, Ra, Pr)
    class(rbc_t), intent(inout), target :: this
    real(kind=rp), intent(in) :: t
    integer, intent(in) :: tstep
    type(coef_t), intent(inout) :: coef 
    real(kind=rp), intent(in) :: Ra
    real(kind=rp), intent(in) :: Pr
    real(kind=rp) :: lambda
    real(kind=rp) :: mu
    
    !> Input fields
    type(field_t), pointer :: u
    type(field_t), pointer :: v
    type(field_t), pointer :: w
    type(field_t), pointer :: p
    type(field_t), pointer :: s
    
    !> Get the needed fields
    u => neko_field_registry%get_field('u')
    v => neko_field_registry%get_field('v')
    w => neko_field_registry%get_field('w')
    p => neko_field_registry%get_field('p')
    s => neko_field_registry%get_field('s')

   
    mu = sqrt(Pr/Ra)
    lambda = 1_rp/sqrt(Ra*Pr)

    
    ! ==================================================================
    ! ===================== " Calculations " ===========================
    ! ==================================================================

    !> Integrate convective current
    !! from vertical temperature derivative
    this%bar_uzt = 0_rp
    call average_from_weights(this%bar_uzt, this%uzt, &
                        this%work_field, coef%B, coef%volume)
    
    !> Integrate thermal dissipation
    this%bar_eps_t = 0_rp
    call average_from_weights(this%bar_eps_t, this%eps_t, &
                        this%work_field, coef%B, coef%volume)

    !> Integrate kinetic dissipation
    this%bar_eps_k = 0_rp
    call average_from_weights(this%bar_eps_k, this%eps_k, &
                        this%work_field, coef%B, coef%volume)

    !> Calculate nusselt from thermal dissipaton
    this%nu_eps_t = 1_rp/lambda*this%bar_eps_t
            
    !> Calculate nusselt from kinetic (the following two formulations are the same)
    !this%nu_eps_k = 1 + 1_rp/lambda*this%bar_eps_k 
    this%nu_eps_k = 1_rp + this%bar_eps_k/((mu**3_rp)*Ra/(Pr**2_rp))
    
    !> Calculate global kolmogorov scale
    this%eta_k = 0_rp
    this%eta_k = ((mu**3_rp)/(this%bar_eps_k))**(0.25_rp)
    
    !> Calculate global batchelor scale
    this%eta_t = 0_rp
    this%eta_t = ((lambda**3_rp)/(this%bar_eps_k))**(0.25_rp)


    !> Integrate vertical temperature derivatives at plates
    this%top_wall_bar_dtdz = 0_rp
    this%bot_wall_bar_dtdz = 0_rp
    call average_from_weights(this%top_wall_bar_dtdz, this%dtdz, &
                        this%work_field, this%mass_area_top%x, &
                        this%top_wall_area)
    call average_from_weights(this%bot_wall_bar_dtdz, this%dtdz, &
                        this%work_field, this%mass_area_bot%x, &
                        this%bot_wall_area)

       
    !> Calculate the mean heatflux at the walls 
    call average_from_weights(this%heat_bot_wall, this%dtdn, &
                              this%work_field, this%mass_area_bot%x, &
                              this%bot_wall_area)
    call average_from_weights(this%heat_top_wall, this%dtdn, &
                              this%work_field, this%mass_area_top%x, &
                              this%top_wall_area)
    call average_from_weights(this%heat_side_wall, this%dtdn, &
                              this%work_field, this%mass_area_side%x, &
                              this%side_wall_area)

    !> Calculate the balance
    this%heat_balance = this%heat_bot_wall + this%heat_top_wall &
                        + this%heat_side_wall
   
    !> Inegrate the divergence of heatflux
    call average_from_weights(this%bar_div_dtdX, this%div_dtdX, &
                        this%work_field, coef%B, coef%volume)


    !> Get turbulent kinetic energy
    this%uuprime = rms_in_space( u, this%work_field, this%work_field2, coef)
    this%uuprime = this%uuprime**2

    this%vvprime = rms_in_space( v, this%work_field, this%work_field2, coef)
    this%vvprime = this%vvprime**2

    this%wwprime = rms_in_space( w, this%work_field, this%work_field2, coef)
    this%wwprime = this%wwprime**2

    this%tke = 0.5_rp*(this%uuprime + this%vvprime + this%wwprime)
                            
    ! ==================================================================
    ! ===================== " IO operations " ==========================
    ! ==================================================================
           
    ! Write them, ideally put this in a buffer            
    if (pe_rank .eq. 0) then
       open(10,file="nusselt.txt",position="append")
       write(10,*) t,'', this%bar_uzt, '', this%top_wall_bar_dtdz, '', &
                   this%bot_wall_bar_dtdz
       close(10)

       open(20,file="bc_heat_balance.txt",position="append")
       write(20,*) t,'', this%heat_top_wall, '', this%heat_bot_wall, '', &
                   this%heat_side_wall, '', this%heat_balance, '', &
                   this%bar_div_dtdX
       close(20)
       
       open(30,file="nu_eps.txt",position="append")
       write(30,*) t,'', this%nu_eps_t, '', this%nu_eps_k
       close(30)
       
       open(40,file="eta.txt",position="append")
       write(40,*) t,'', this%eta_t, '', this%eta_k, &
                     '', this%bar_eps_t, '', this%bar_eps_k
       close(40)
       
       open(50,file="tke.txt",position="append")
       write(50,*) t,'', this%tke, '', this%uuprime, &
                     '', this%vvprime, '', this%wwprime
       close(50)
    end if


  end subroutine rbc_get_integral_quantities


  subroutine rbc_sync_fromGPU(this)
    class(rbc_t), intent(inout), target :: this
    integer :: n
  
    n = this%dudx%dof%size()
    
    if (pe_rank == 0) write(*,*) 'Peforming user field sync with GPU-CPU'
    
    if (NEKO_BCKND_DEVICE .eq. 1) then 

       call device_memcpy(this%work_field%x, &
                          this%work_field%x_d, &
                          n,DEVICE_TO_HOST,sync=.true.)

       call device_memcpy(this%uzt%x, &
                          this%uzt%x_d, &
                          n,DEVICE_TO_HOST,sync=.true.)

       call device_memcpy(this%dtdx%x, &
                          this%dtdx%x_d, &
                          n,DEVICE_TO_HOST,sync=.true.)

       call device_memcpy(this%dtdy%x, &
                          this%dtdy%x_d, &
                          n,DEVICE_TO_HOST,sync=.true.)

       call device_memcpy(this%dtdz%x, &
                          this%dtdz%x_d, &
                          n,DEVICE_TO_HOST,sync=.true.)
       
       call device_memcpy(this%div_dtdX%x, &
                          this%div_dtdX%x_d, &
                          n,DEVICE_TO_HOST,sync=.true.)
       
       call device_memcpy(this%dtdn%x, &
                          this%dtdn%x_d, &
                          n,DEVICE_TO_HOST,sync=.true.)
       
       call device_memcpy(this%dudx%x, &
                          this%dudx%x_d, &
                          n,DEVICE_TO_HOST,sync=.true.)
       
       call device_memcpy(this%dudy%x, &
                          this%dudy%x_d, &
                          n,DEVICE_TO_HOST,sync=.true.)
       
       call device_memcpy(this%dudz%x, &
                          this%dudz%x_d, &
                          n,DEVICE_TO_HOST,sync=.true.)
       
       call device_memcpy(this%dvdx%x, &
                          this%dvdx%x_d, &
                          n,DEVICE_TO_HOST,sync=.true.)
       
       call device_memcpy(this%dvdy%x, &
                          this%dvdy%x_d, &
                          n,DEVICE_TO_HOST,sync=.true.)
       
       call device_memcpy(this%dvdz%x, &
                          this%dvdz%x_d, &
                          n,DEVICE_TO_HOST,sync=.true.)
       
       call device_memcpy(this%dwdx%x, &
                          this%dwdx%x_d, &
                          n,DEVICE_TO_HOST,sync=.true.)
       
       call device_memcpy(this%dwdy%x, &
                          this%dwdy%x_d, &
                          n,DEVICE_TO_HOST,sync=.true.)
       
       call device_memcpy(this%dwdz%x, &
                          this%dwdz%x_d, &
                          n,DEVICE_TO_HOST,sync=.true.)
       
       call device_memcpy(this%eps_k%x, &
                          this%eps_k%x_d, &
                          n,DEVICE_TO_HOST,sync=.true.)
       
       call device_memcpy(this%eps_t%x, &
                          this%eps_t%x_d, &
                          n,DEVICE_TO_HOST,sync=.true.)
       
       !> Mean quantities
       
       call device_memcpy(this%mean_t%mf%x, &
                          this%mean_t%mf%x_d, &
                          n,DEVICE_TO_HOST,sync=.true.)
       
       call device_memcpy(this%mean_tt%mf%x, &
                          this%mean_tt%mf%x_d, &
                          n,DEVICE_TO_HOST,sync=.true.)
       
       call device_memcpy(this%mean_uxt%mf%x, &
                          this%mean_uxt%mf%x_d, &
                          n,DEVICE_TO_HOST,sync=.true.)
       
       call device_memcpy(this%mean_uyt%mf%x, &
                          this%mean_uyt%mf%x_d, &
                          n,DEVICE_TO_HOST,sync=.true.)
       
       call device_memcpy(this%mean_uzt%mf%x, &
                          this%mean_uzt%mf%x_d, &
                          n,DEVICE_TO_HOST,sync=.true.)
       
       call device_memcpy(this%mean_dtdx%mf%x, &
                          this%mean_dtdx%mf%x_d, &
                          n,DEVICE_TO_HOST,sync=.true.)
       
       call device_memcpy(this%mean_dtdy%mf%x, &
                          this%mean_dtdy%mf%x_d, &
                          n,DEVICE_TO_HOST,sync=.true.)
       
       call device_memcpy(this%mean_dtdz%mf%x, &
                          this%mean_dtdz%mf%x_d, &
                          n,DEVICE_TO_HOST,sync=.true.)
       
       call device_memcpy(this%mean_dtdn%mf%x, &
                          this%mean_dtdn%mf%x_d, &
                          n,DEVICE_TO_HOST,sync=.true.)
       
       call device_memcpy(this%mean_eps_k%mf%x, &
                          this%mean_eps_k%mf%x_d, &
                          n,DEVICE_TO_HOST,sync=.true.)
       
       call device_memcpy(this%mean_eps_t%mf%x, &
                          this%mean_eps_t%mf%x_d, &
                          n,DEVICE_TO_HOST,sync=.true.)

       call device_memcpy(this%mean_speri_u%mf%x, &
                          this%mean_speri_u%mf%x_d, &
                          n,DEVICE_TO_HOST,sync=.true.)

       call device_memcpy(this%mean_speri_v%mf%x, &
                          this%mean_speri_v%mf%x_d, &
                          n,DEVICE_TO_HOST,sync=.true.)

       call device_memcpy(this%mean_speri_w%mf%x, &
                          this%mean_speri_w%mf%x_d, &
                          n,DEVICE_TO_HOST,sync=.true.)
       
    end if

  end subroutine rbc_sync_fromGPU



  !==================================================
  !> Supporting subroutines
  !==================================================

  subroutine map_from_facets_to_field(u, area, wall_facet)
    !> this is needed because the facet data is arrayed differently than filed data
    !> u(lx,ly,lz,e) is the field where you want to put the facet data
    !> area(lx,ly,facet,e) is the array in facet form
    !> wall_facet is the list of facets that should be mapped
    type(field_t), intent(inout) :: u
    real(kind=rp) :: area(u%Xh%lx,u%Xh%ly,6,u%msh%nelv)
    type(stack_i4t4_t) :: wall_facet

    !> Parameters to populate the list of elements and facets
    integer :: e, facet, e_counter, facet_counter, i,j,k,n
    type(tuple4_i4_t), pointer :: wall_facet_list(:)
    integer :: lx, ly, lz

    wall_facet_list => wall_facet%array()
    !! Initialize as 0. This serves as a mask for nodes that are not integrated
    call rzero(u%x,u%dof%size())
    !! Fill up the top and bottom mass matrices
    lx = u%Xh%lx
    ly = u%Xh%ly
    lz = u%Xh%lz
    do e_counter = 1, wall_facet%top_
      e = wall_facet_list(e_counter)%x(2)
      facet = wall_facet_list(e_counter)%x(1)
      select case(facet)
         case(1)
            do j = 1 , ly
               do k = 1 , lz
                  u%x(1,j,k,e) = area(j,k,facet,e) 
               end do
            end do
         case(2)
            do j = 1 , ly
               do k = 1 , lz
                  u%x(lx,j,k,e) = area(j,k,facet,e) 
               end do
            end do
         case(3)
            do i = 1 , lx
               do k = 1 , lz
                  u%x(i,1,k,e) = area(i,k,facet,e)    
               end do
            end do
         case(4)
            do i = 1 , lx
               do k = 1 , lz
                  u%x(i,ly,k,e) = area(i,k,facet,e)       
               end do
            end do
         case(5)
            do i = 1 , lx
               do j = 1 , ly
                  u%x(i,j,1,e) = area(i,j,facet,e)
               end do
            end do
         case(6)
            do i = 1 , lx
               do j = 1 , ly
                  u%x(i,j,lz,e) = area(i,j,facet,e)
               end do
            end do
      end select
    end do


    ! Put it also on device
    if (NEKO_BCKND_DEVICE .eq. 1) then 
       call device_memcpy(u%x,u%x_d, &
                          u%dof%size(),HOST_TO_DEVICE)
    end if

  end subroutine map_from_facets_to_field


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


  subroutine project_to_vector(proj, uu, vv, ww, nx, ny, nz)
    type(field_t), intent(inout) :: proj
    type(field_t), intent(inout) :: uu
    type(field_t), intent(inout) :: vv
    type(field_t), intent(inout) :: ww
    real(kind=rp), dimension(proj%Xh%lxyz,proj%msh%nelv), intent(inout) :: nx
    real(kind=rp), dimension(proj%Xh%lxyz,proj%msh%nelv), intent(inout) :: ny
    real(kind=rp), dimension(proj%Xh%lxyz,proj%msh%nelv), intent(inout) :: nz
    integer :: n
    type(c_ptr) :: nx_d
    type(c_ptr) :: ny_d
    type(c_ptr) :: nz_d
    
    n = proj%dof%size()

    ! Find the normal temperature derivative at the boundaries
    if (NEKO_BCKND_DEVICE .eq. 1) then 
       nx_d = device_get_ptr(nx)
       ny_d = device_get_ptr(ny)
       nz_d = device_get_ptr(nz)
       call device_rzero(proj%x_d,n)
       call device_vdot3(proj%x_d, uu%x_d,vv%x_d,ww%x_d, &
                         nx_d,ny_d,nz_d, n) 
    else
       call rzero(proj%x,n)
       call vdot3(proj%x, uu%x,vv%x,ww%x, &
                         nx,ny,nz, n) 
    end if

  end subroutine project_to_vector


  subroutine divergence_of_field(div_u, u, v, w, work_field, coef)
    type(field_t), intent(inout) :: div_u
    type(field_t), intent(inout) :: u
    type(field_t), intent(inout) :: v
    type(field_t), intent(inout) :: w
    type(field_t), intent(inout) :: work_field
    type(coef_t), intent(inout) :: coef
    integer :: n
    
    n = div_u%dof%size()
    
    !dudx in div_u
    call dudxyz (div_u%x, u%x, coef%drdx, coef%dsdx, coef%dtdx, coef)
    !dvdy in work_field
    call dudxyz (work_field%x, v%x, coef%drdy, coef%dsdy, coef%dtdy, coef)
    ! Add dudx + dvdy    
    if (NEKO_BCKND_DEVICE .eq. 1) then 
          call device_add2(div_u%x_d,work_field%x_d,n)
    else
          call add2(div_u%x,work_field%x,n)
    end if
    ! dwdz in work_dielf
    call dudxyz (work_field%x, w%x, coef%drdz, coef%dsdz, coef%dtdz, coef)
    ! Add (dudx + dvdy) + dwdz 
    if (NEKO_BCKND_DEVICE .eq. 1) then 
          call device_add2(div_u%x_d,work_field%x_d,n)
    else
          call add2(div_u%x,work_field%x,n)
    end if
     
  end subroutine divergence_of_field


  subroutine calculate_thermal_dissipation(eps_t, dtdx,dtdy,dtdz, work_field,Ra,Pr)
    type(field_t), intent(inout) :: eps_t
    type(field_t), intent(inout) :: dtdx
    type(field_t), intent(inout) :: dtdy
    type(field_t), intent(inout) :: dtdz
    type(field_t), intent(inout) :: work_field
    real(kind=rp), intent(in) :: Ra
    real(kind=rp), intent(in) :: Pr
     
    real(kind=rp) :: lambda
    integer :: n

    n = eps_t%dof%size()
    lambda = 1_rp/sqrt(Ra*Pr)

    !Thermal dissipation epsT=(grad T)**2/sqrt(RaPr)
    !jorg did it as epsT=((dtdx)**2+(dtdy)**2+(dtdz)**2)/sqrt(RaPr)
    if (NEKO_BCKND_DEVICE .eq. 1) then 
        call device_rzero(eps_t%x_d,n)
        call device_addsqr2s2(eps_t%x_d, dtdx%x_d, lambda, n)
        call device_addsqr2s2(eps_t%x_d, dtdy%x_d, lambda, n)
        call device_addsqr2s2(eps_t%x_d, dtdz%x_d, lambda, n)
    else
        call rzero(eps_t%x,n)
        call addsqr2s2(eps_t%x, dtdx%x, lambda, n)
        call addsqr2s2(eps_t%x, dtdy%x, lambda, n)
        call addsqr2s2(eps_t%x, dtdz%x, lambda, n)
  
    end if

  end subroutine calculate_thermal_dissipation
  
  
  
  subroutine calculate_kinetic_dissipation(eps_k, dudx,dudy,dudz, &
                                                  dvdx,dvdy,dvdz, &
                                                  dwdx,dwdy,dwdz, &
                                                  work_field,Ra,Pr)
    type(field_t), intent(inout) :: eps_k
    type(field_t), intent(inout) :: dudx
    type(field_t), intent(inout) :: dudy
    type(field_t), intent(inout) :: dudz
    type(field_t), intent(inout) :: dvdx
    type(field_t), intent(inout) :: dvdy
    type(field_t), intent(inout) :: dvdz
    type(field_t), intent(inout) :: dwdx
    type(field_t), intent(inout) :: dwdy
    type(field_t), intent(inout) :: dwdz
    type(field_t), intent(inout) :: work_field
    real(kind=rp), intent(in) :: Ra
    real(kind=rp), intent(in) :: Pr
     
    real(kind=rp) :: mu_half, c1 = 1_rp, c2 = 1_rp
    integer :: n

    n = eps_k%dof%size()
    mu_half = sqrt(Pr/Ra)*0.5_rp ! 0.5*mu

    !Energy dissipation epsv=0.5*(du_i/dx_j+du_j/dx_i)**2*sqrt(Pr/Ra)
    ! Some explanation https://www.cfd-online.com/Wiki/Introduction_to_turbulence/Turbulence_kinetic_energy
    if (NEKO_BCKND_DEVICE .eq. 1) then 
        call device_rzero(eps_k%x_d,n)
        
        call device_add3s2(work_field%x_d,dudx%x_d,dudx%x_d,c1,c1,n)
        call device_addsqr2s2(eps_k%x_d, work_field%x_d, mu_half, n)
        call device_add3s2(work_field%x_d,dudy%x_d,dvdx%x_d,c1,c1,n)
        call device_addsqr2s2(eps_k%x_d, work_field%x_d, mu_half, n)
        call device_add3s2(work_field%x_d,dudz%x_d,dwdx%x_d,c1,c1,n)
        call device_addsqr2s2(eps_k%x_d, work_field%x_d, mu_half, n)

        call device_add3s2(work_field%x_d,dvdx%x_d,dudy%x_d,c1,c1,n)
        call device_addsqr2s2(eps_k%x_d, work_field%x_d, mu_half, n)
        call device_add3s2(work_field%x_d,dvdy%x_d,dvdy%x_d,c1,c1,n)
        call device_addsqr2s2(eps_k%x_d, work_field%x_d, mu_half, n)
        call device_add3s2(work_field%x_d,dvdz%x_d,dwdy%x_d,c1,c1,n)
        call device_addsqr2s2(eps_k%x_d, work_field%x_d, mu_half, n)

        call device_add3s2(work_field%x_d,dwdx%x_d,dudz%x_d,c1,c1,n)
        call device_addsqr2s2(eps_k%x_d, work_field%x_d, mu_half, n)
        call device_add3s2(work_field%x_d,dwdy%x_d,dvdz%x_d,c1,c1,n)
        call device_addsqr2s2(eps_k%x_d, work_field%x_d, mu_half, n)
        call device_add3s2(work_field%x_d,dwdz%x_d,dwdz%x_d,c1,c1,n)
        call device_addsqr2s2(eps_k%x_d, work_field%x_d, mu_half, n)
    else
        call rzero(eps_k%x,n)
        
        call add3s2(work_field%x,dudx%x,dudx%x,c1,c1,n)
        call addsqr2s2(eps_k%x, work_field%x, mu_half, n)
        call add3s2(work_field%x,dudy%x,dvdx%x,c1,c1,n)
        call addsqr2s2(eps_k%x, work_field%x, mu_half, n)
        call add3s2(work_field%x,dudz%x,dwdx%x,c1,c1,n)
        call addsqr2s2(eps_k%x, work_field%x, mu_half, n)

        call add3s2(work_field%x,dvdx%x,dudy%x,c1,c1,n)
        call addsqr2s2(eps_k%x, work_field%x, mu_half, n)
        call add3s2(work_field%x,dvdy%x,dvdy%x,c1,c1,n)
        call addsqr2s2(eps_k%x, work_field%x, mu_half, n)
        call add3s2(work_field%x,dvdz%x,dwdy%x,c1,c1,n)
        call addsqr2s2(eps_k%x, work_field%x, mu_half, n)

        call add3s2(work_field%x,dwdx%x,dudz%x,c1,c1,n)
        call addsqr2s2(eps_k%x, work_field%x, mu_half, n)
        call add3s2(work_field%x,dwdy%x,dvdz%x,c1,c1,n)
        call addsqr2s2(eps_k%x, work_field%x, mu_half, n)
        call add3s2(work_field%x,dwdz%x,dwdz%x,c1,c1,n)
        call addsqr2s2(eps_k%x, work_field%x, mu_half, n)

    end if

  end subroutine calculate_kinetic_dissipation


  function rms_in_space( u, work_field, work_field2, coef) result(urms)
    type(field_t), intent(inout) :: u
    type(field_t), intent(inout) :: work_field
    type(field_t), intent(inout) :: work_field2
    type(coef_t), intent(inout) :: coef
    real(kind=rp) :: urms
    real(kind=rp) :: u_mean
    integer :: n
    real(kind=rp) :: c1 = 1_rp, c2 = 1_rp

    n = u%dof%size()
    
    !get the mean of velocity in space
    u_mean = 0_rp    
    call average_from_weights(u_mean, u, &
                        work_field2, coef%B, coef%volume)


    if (NEKO_BCKND_DEVICE .eq. 1) then 
        
       call device_rzero(work_field%x_d,n)
       call device_rzero(work_field2%x_d,n)

       ! Substract the mean from the field
       call device_copy(work_field2%x_d, u%x_d, n)
       call device_cadd(work_field2%x_d, -1_rp * u_mean, n)

       ! Square the value
       call device_addsqr2s2(work_field%x_d, work_field2%x_d, c1, n)

    else
       
       call rzero(work_field%x,n)
       call rzero(work_field2%x,n)
       
       ! Substract the mean from the field
       call copy(work_field2%x, u%x, n)
       call cadd(work_field2%x, -1_rp * u_mean, n)

       ! Square the value
       call addsqr2s2(work_field%x, work_field2%x, c1 , n)
       
    end if

    !Average in space again
    urms = 0_rp    
    call average_from_weights(urms, work_field, &
                        work_field2, coef%B, coef%volume)

    !get rms
    urms = sqrt(urms)

  end function rms_in_space

end module rbc

