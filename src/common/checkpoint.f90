!> Defines a checkpoint
module checkpoint
  use num_types
  use field
  use utils
  implicit none
  private

  type, public :: chkp_t
     type(field_t), pointer :: u => null()
     type(field_t), pointer :: v => null()
     type(field_t), pointer :: w => null()
     type(field_t), pointer :: p => null()

     !
     ! Optional payload
     !
     real(kind=rp), pointer :: ulag(:,:,:,:,:) => null()
     real(kind=rp), pointer :: vlag(:,:,:,:,:) => null()
     real(kind=rp), pointer :: wlag(:,:,:,:,:) => null()

     real(kind=dp) :: t         !< Restart time (valid after load)
   contains
     procedure, pass(this) :: init => chkp_init
     procedure, pass(this) :: add_lag => chkp_add_lag
     procedure, pass(this) :: restart_time => chkp_restart_time
     final :: chkp_free
  end type chkp_t

contains

  !> Initialize checkpoint structure with mandatory data
  subroutine chkp_init(this, u, v, w, p)
    class(chkp_t), intent(inout) :: this
    type(field_t), intent(in), target :: u
    type(field_t), intent(in), target :: v
    type(field_t), intent(in), target :: w
    type(field_t), intent(in), target :: p

    ! Check that all velocity components are defined on the same
    ! function space
    if ( u%Xh .ne. v%Xh .or. &
         u%Xh .ne. w%Xh ) then
       call neko_error('Different function spaces for each velocity component')
    end if

    ! Check that both velocity and pressure is defined on the same mesh
    if ( u%msh%nelv .ne. p%msh%nelv ) then
       call neko_error('Velocity and pressure defined on different meshes')
    end if
    
    this%u => u
    this%v => v
    this%w => w
    this%p => p

    this%t = 0d0
    
  end subroutine chkp_init

  !> Reset checkpoint
  subroutine chkp_free(this)
    type(chkp_t), intent(inout) :: this

    nullify(this%u)
    nullify(this%v)
    nullify(this%w)
    nullify(this%p)

    nullify(this%ulag)
    nullify(this%vlag)
    nullify(this%wlag)
    
  end subroutine chkp_free

  !> Add lagged velocity terms
  subroutine chkp_add_lag(this, ulag, vlag, wlag)
    class(chkp_t), intent(inout) :: this
    real(kind=rp), target :: ulag(:,:,:,:,:)
    real(kind=rp), target :: vlag(:,:,:,:,:)
    real(kind=rp), target :: wlag(:,:,:,:,:)

    this%ulag => ulag
    this%vlag => vlag
    this%wlag => wlag
    
  end subroutine chkp_add_lag

  !> Return restart time from a loaded checkpoint
  pure function chkp_restart_time(this) result(rtime)
    class(chkp_t), intent(in) :: this
    real(kind=dp) :: rtime

    rtime = this%t
  end function chkp_restart_time
  
end module checkpoint
