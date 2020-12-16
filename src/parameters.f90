!> NEKO parameters
module parameters
  use num_types
  implicit none  

  type param_t
     real(kind=dp) :: dt        !< time-step size     
     integer :: nsteps          !< Number of time-stpes     
     real(kind=dp) :: rho       !< Density \f$ \rho \f$
     real(kind=dp) :: mu        !< Dynamic viscosity \f$ \mu \f$
     real(kind=dp) :: Re        !< Reynolds number
     real(kind=dp), dimension(3) :: uinf !< Free-stream velocity \f$ u_\infty \f$
     logical :: output_bdry              !< Output boundary markings
  end type param_t

  type param_io_t
     type(param_t) p
   contains
     procedure  :: param_read
     generic :: read(formatted) => param_read
  end type param_io_t

  interface write(formatted)
     module procedure :: param_write
  end interface write(formatted)
  
contains

  subroutine param_read(param, unit, iotype, v_list, iostat, iomsg)
    class(param_io_t), intent(inout) ::  param
    integer(kind=4), intent(in) :: unit
    character(len=*), intent(in) :: iotype
    integer, intent(in) :: v_list(:)
    integer(kind=4), intent(out) :: iostat
    character(len=*), intent(inout) :: iomsg

    real(kind=dp) :: dt = 0d0
    integer :: nsteps = 0
    real(kind=dp) :: rho = 1d0
    real(kind=dp) :: mu = 1d0
    real(kind=dp) :: Re = 1d0
    real(kind=dp), dimension(3) :: uinf = (/ 0d0, 0d0, 0d0 /)
    logical :: output_bdry = .false.
    namelist /NEKO_PARAMETERS/ dt, nsteps, rho, mu, Re, uinf, output_bdry

    read(unit, nml=NEKO_PARAMETERS, iostat=iostat, iomsg=iomsg)
    param%p%dt = dt
    param%p%nsteps = nsteps 
    param%p%rho = rho 
    param%p%mu = mu
    param%p%Re = Re
    param%p%uinf = uinf
    param%p%output_bdry = output_bdry

  end subroutine param_read

  subroutine param_write(param, unit, iotype, v_list, iostat, iomsg)
    class(param_io_t), intent(in) ::  param
    integer(kind=4), intent(in) :: unit
    character(len=*), intent(in) :: iotype
    integer, intent(in) :: v_list(:)
    integer(kind=4), intent(out) :: iostat
    character(len=*), intent(inout) :: iomsg

    real(kind=dp) :: dt, rho, mu, Re
    real(kind=dp), dimension(3) :: uinf
    logical :: output_bdry
    integer :: nsteps
    namelist /NEKO_PARAMETERS/ dt, nsteps, rho, mu, Re, uinf, output_bdry

    dt = param%p%dt
    nsteps = param%p%nsteps
    rho = param%p%rho  
    mu = param%p%mu
    Re = param%p%Re
    uinf = param%p%uinf
    output_bdry = param%p%output_bdry
    
    write(unit, nml=NEKO_PARAMETERS, iostat=iostat, iomsg=iomsg)

        
  end subroutine param_write

  
end module parameters

