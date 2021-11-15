!> NEKO parameters
module parameters
  use num_types
  implicit none  

  type param_t
     integer :: nsamples        !< Number of samples
     logical :: output_bdry     !< Output boundary markings
     logical :: output_part     !< Output partitions
     logical :: output_chkp     !< Output checkpoints
     real(kind=rp) :: dt        !< time-step size               
     real(kind=rp) :: T_end     !< Final time
     real(kind=rp) :: rho       !< Density \f$ \rho \f$
     real(kind=rp) :: mu        !< Dynamic viscosity \f$ \mu \f$
     real(kind=rp) :: Re        !< Reynolds number
     real(kind=rp), dimension(3) :: uinf !< Free-stream velocity \f$ u_\infty \f$
     real(kind=rp) :: abstol_vel  !< Tolerance for velocity solver
     real(kind=rp) :: abstol_prs  !< Tolerance for pressure solver
     character(len=20) :: ksp_vel !< Krylov solver for velocity 
     character(len=20) :: ksp_prs !< Krylov solver for pressure
     character(len=20) :: pc_vel  !< Precon for velocity solver
     character(len=20) :: pc_prs  !< Precon for pressure solver
     character(len=20) :: fluid_inflow !< Fluid inflow condition
     integer :: vol_flow_dir !< Direction of forced volume flow x=1, y=2, z=3
     logical :: avflow       !< If we should use the averaged flow for vol_flow
     logical :: loadb        !< Load-balancing
     real(kind=rp) :: flow_rate !< Volume flow speed
     integer :: proj_dim     !< Projection space for pressure solution
     integer :: time_order   !< Order of the time stepping
     character(len=8) :: jlimit !< Job limit in HH:MM:SS
     character(len=80) :: restart_file !< Checkpoint filename
     real(kind=rp) :: stats_begin      !< Start time for statistics
     logical :: stats_mean_flow        !< Mean flow statistics
     logical :: output_mean_flow       !< Output mean flow field
     logical :: stats_mean_sqr_flow    !< Mean squared flow statistics
     logical :: output_mean_sqr_flow   !< Output mean squared flow field
     character(len=1024) :: output_dir !< Output directory
     real(kind=rp) :: delta !< Boundary layer thickness \f$ \delta \f$
     character(len=10) :: blasius_approx !< Type of approximate Blasius profile
  end type param_t

  type param_io_t
     type(param_t) p
   contains
     procedure  :: param_read
     generic :: read(formatted) => param_read
  end type param_io_t

  interface write(formatted)
     module procedure param_write
  end interface write(formatted)
  
contains

  subroutine param_read(param, unit, iotype, v_list, iostat, iomsg)
    class(param_io_t), intent(inout) ::  param
    integer(kind=4), intent(in) :: unit
    character(len=*), intent(in) :: iotype
    integer, intent(in) :: v_list(:)
    integer(kind=4), intent(out) :: iostat
    character(len=*), intent(inout) :: iomsg

    integer :: nsamples = 0
    logical :: output_bdry = .false.
    logical :: output_part = .false.
    logical :: output_chkp = .false.
    real(kind=rp) :: dt = 0d0
    real(kind=rp) :: T_end = 0d0
    real(kind=rp) :: rho = 1d0
    real(kind=rp) :: mu = 1d0
    real(kind=rp) :: Re = 1d0
    real(kind=rp), dimension(3) :: uinf = (/ 0d0, 0d0, 0d0 /)
    real(kind=rp) :: abstol_vel = 1d-9
    real(kind=rp) :: abstol_prs = 1d-9
    real(kind=rp) :: delta = 1d0
    character(len=20) :: ksp_vel = 'cg'
    character(len=20) :: ksp_prs = 'gmres'
    character(len=20) :: pc_vel = 'jacobi'
    character(len=20) :: pc_prs = 'hsmg'
    character(len=20) :: fluid_inflow = 'default'
    integer :: vol_flow_dir = 0
    logical :: avflow = .true.
    logical :: loadb = .false.
    real(kind=rp) :: flow_rate = 0d0
    integer :: proj_dim = 20
    integer :: time_order = 3
    character(len=8) :: jlimit = '00:00:00'
    character(len=80) :: restart_file = ''
    real(kind=rp) :: stats_begin = 0d0
    logical :: stats_mean_flow = .false.
    logical :: output_mean_flow = .false.
    logical :: stats_mean_sqr_flow = .false.
    logical :: output_mean_sqr_flow = .false.
    character(len=1024) :: output_dir = ''
    character(len=10) :: blasius_approx = 'sin'
    
    namelist /NEKO_PARAMETERS/ nsamples, output_bdry, output_part, output_chkp, &
         dt, T_end, rho, mu, Re, uinf, abstol_vel, abstol_prs, ksp_vel, ksp_prs, &
         pc_vel, pc_prs, fluid_inflow, vol_flow_dir, loadb, avflow, flow_rate, &
         proj_dim, time_order, jlimit, restart_file, stats_begin, &
         stats_mean_flow, output_mean_flow, stats_mean_sqr_flow, &
         output_mean_sqr_flow, output_dir, delta, blasius_approx

    read(unit, nml=NEKO_PARAMETERS, iostat=iostat, iomsg=iomsg)

    param%p%output_bdry = output_bdry
    param%p%output_part = output_part
    param%p%output_chkp = output_chkp
    param%p%nsamples = nsamples
    param%p%dt = dt
    param%p%T_end = T_end
    param%p%rho = rho
    param%p%mu = mu
    param%p%Re = Re
    param%p%uinf = uinf
    param%p%abstol_vel = abstol_vel
    param%p%abstol_prs = abstol_prs
    param%p%ksp_vel = ksp_vel
    param%p%ksp_prs = ksp_prs
    param%p%pc_vel = pc_vel
    param%p%pc_prs = pc_prs
    param%p%fluid_inflow = fluid_inflow
    param%p%vol_flow_dir = vol_flow_dir
    param%p%avflow = avflow
    param%p%loadb = loadb
    param%p%flow_rate = flow_rate
    param%p%proj_dim = proj_dim
    param%p%time_order = time_order
    param%p%jlimit = jlimit
    param%p%restart_file = restart_file
    param%p%stats_begin = stats_begin
    param%p%stats_mean_flow = stats_mean_flow
    param%p%output_mean_flow = output_mean_flow
    param%p%stats_mean_sqr_flow = stats_mean_sqr_flow
    param%p%output_mean_sqr_flow = output_mean_sqr_flow
    param%p%output_dir = output_dir
    param%p%delta = delta
    param%p%blasius_approx = blasius_approx

  end subroutine param_read

  subroutine param_write(param, unit, iotype, v_list, iostat, iomsg)
    class(param_io_t), intent(in) ::  param
    integer(kind=4), intent(in) :: unit
    character(len=*), intent(in) :: iotype
    integer, intent(in) :: v_list(:)
    integer(kind=4), intent(out) :: iostat
    character(len=*), intent(inout) :: iomsg

    real(kind=rp) :: dt, T_End, rho, mu, Re, abstol_vel, abstol_prs, flow_rate
    real(kind=rp) :: stats_begin, delta
    character(len=20) :: ksp_vel, ksp_prs, pc_vel, pc_prs, fluid_inflow
    real(kind=rp), dimension(3) :: uinf
    logical :: output_part, output_bdry, output_chkp
    logical :: avflow, loadb, stats_mean_flow, output_mean_flow
    logical :: stats_mean_sqr_flow, output_mean_sqr_flow
    integer :: nsamples, vol_flow_dir, proj_dim, time_order
    character(len=8) :: jlimit    
    character(len=80) :: restart_file
    character(len=1024) :: output_dir
    character(len=10) :: blasius_approx

    namelist /NEKO_PARAMETERS/ nsamples, output_bdry, output_part, output_chkp, &
         dt, T_end, rho, mu, Re, uinf, abstol_vel, abstol_prs, ksp_vel, ksp_prs, &
         pc_vel, pc_prs, fluid_inflow, vol_flow_dir, avflow, loadb, flow_rate, &
         proj_dim, time_order, jlimit, restart_file, stats_begin, &
         stats_mean_flow, output_mean_flow, stats_mean_sqr_flow, &
         output_mean_sqr_flow, output_dir, delta, blasius_approx

    nsamples = param%p%nsamples
    output_bdry = param%p%output_bdry
    output_part = param%p%output_part
    output_chkp = param%p%output_chkp
    dt = param%p%dt
    T_end = param%p%T_end
    rho = param%p%rho
    mu = param%p%mu
    Re = param%p%Re
    uinf = param%p%uinf
    abstol_vel = param%p%abstol_vel
    abstol_prs = param%p%abstol_prs
    ksp_vel = param%p%ksp_vel
    ksp_prs = param%p%ksp_prs
    pc_vel = param%p%pc_vel
    pc_prs = param%p%pc_prs
    fluid_inflow = param%p%fluid_inflow
    vol_flow_dir = param%p%vol_flow_dir
    avflow = param%p%avflow
    loadb = param%p%loadb
    flow_rate = param%p%flow_rate
    proj_dim = param%p%proj_dim
    time_order = param%p%time_order
    jlimit = param%p%jlimit
    restart_file = param%p%restart_file
    stats_begin = param%p%stats_begin
    stats_mean_flow = param%p%stats_mean_flow
    output_mean_flow = param%p%output_mean_flow
    stats_mean_sqr_flow = param%p%stats_mean_sqr_flow
    output_mean_sqr_flow = param%p%output_mean_sqr_flow
    output_dir = param%p%output_dir
    delta = param%p%delta
    blasius_approx = param%p%blasius_approx
    
    write(unit, nml=NEKO_PARAMETERS, iostat=iostat, iomsg=iomsg)

        
  end subroutine param_write

  
end module parameters

