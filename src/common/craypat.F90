!> Interface to CrayPat F77 API
module craypat
  use iso_c_binding
  use stack
  use utils
  implicit none

  type(stack_i4_t), private :: region_depth
  logical, private :: craypat_on = .false.
#ifdef CRAYPAT
  include 'pat_apif.h'

contains

  !> Turn on CrayPat recording
  subroutine craypat_record_start
    integer :: ierr
    call region_depth%init()
    call PAT_record(PAT_STATE_ON, ierr)
    craypat_on = .true.
  end subroutine craypat_record_start

  !> Turn off CrayPat recording
  subroutine craypat_record_stop
    integer :: ierr
    call PAT_record(PAT_STATE_OFF, ierr)
    craypat_on = .false.
  end subroutine craypat_record_stop

  !> Start a CrayPat region
  subroutine craypat_region_begin(name)
    character(kind=c_char,len=*) :: name
    integer :: ierr, region_id
    
    if (craypat_on) then
       !> @todo Don't hardcode region names...
       if (name .eq. 'Time-Step') then
          region_id = 1
       else if(name .eq. 'Pressure') then
          region_id = 2
       else if (name .eq. 'Velocity') then
          region_id = 3
       else if (name .eq. 'gather-scatter') then
          region_id = 4
       else ! User defined region
          region_id = 99 
       end if
       
       call region_depth%push(region_id)
       call PAT_region_begin(region_id, name, ierr)
    end if
    
  end subroutine craypat_region_begin

  !> End a CrayPat region
  subroutine craypat_region_end
    integer :: ierr

    if (craypat_on) then
       if (region_depth%size() .le. 0) then
          call neko_error('Invalid CrayPat region id')
       end if
       call PAT_region_end(region_depth%pop(), ierr)
    end if
    
  end subroutine craypat_region_end

#endif
 
end module craypat
