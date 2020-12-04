!> NEKO parameters
module parameters
  use num_types
  implicit none  

  type param_t
     real(kind=dp) :: dt        !< time-step size     
     integer :: nsteps          !< Number of time-stpes     
     real(kind=dp) :: rho       !< Density     
     real(kind=dp) :: mu        !< Dynamic viscosity   
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
    namelist /NEKO_PARAMETERS/ dt, nsteps, rho, mu

    read(unit, nml=NEKO_PARAMETERS, iostat=iostat)
    param%p%dt = dt
    param%p%nsteps = nsteps 
    param%p%rho = rho 
    param%p%mu = mu 
        
  end subroutine param_read

  subroutine param_write(param, unit, iotype, v_list, iostat, iomsg)
    class(param_io_t), intent(in) ::  param
    integer(kind=4), intent(in) :: unit
    character(len=*), intent(in) :: iotype
    integer, intent(in) :: v_list(:)
    integer(kind=4), intent(out) :: iostat
    character(len=*), intent(inout) :: iomsg

    real(kind=dp) :: dt, rho, mu
    integer :: nsteps
    namelist /NEKO_PARAMETERS/ dt, nsteps, rho, mu

    dt = param%p%dt
    nsteps = param%p%nsteps
    rho = param%p%rho  
    mu = param%p%mu 

    write(unit, nml=NEKO_PARAMETERS)
        
  end subroutine param_write

  
end module parameters

