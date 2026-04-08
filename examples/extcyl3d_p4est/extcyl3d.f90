! Two-dimensional flow past circular cylinder cavity
!
! Note that the domain is actually 3D with width one element. In order
! to prevent any instability in the z direction, the w velocity is
! set to zero at every step. This is needed for higher Reynolds numbers.
!
module user
  use neko
  use amr_reconstruct, only : amr_reconstruct_t, amr_flg_none, amr_flg_h_ref, &
       amr_flg_h_crs
  use fld_file, only : fld_file_t
  use fld_file_output, only : fld_file_output_t
  implicit none

  ! Global user variables
  type(fld_file_output_t), save :: tst_fld
  real(rp) :: t

contains

  ! Register user-defined functions (see user_intf.f90)
  subroutine user_setup(user)
    type(user_t), intent(inout) :: user
    user%initialize => user_initialize
    user%compute => compute
    user%amr_refine_flag => amr_refine_flag
    user%amr_reconstruct => amr_reconstructing
  end subroutine user_setup

  ! User-defined initialization called just before time loop starts
  ! I use it to set up AMR refinement output
  subroutine user_initialize(time)
    type(time_state_t), intent(in) :: time
    real(kind=rp) :: t
    type(field_t), pointer :: u, v, w, p

    u => neko_registry%get_field('u')
    v => neko_registry%get_field('v')
    w => neko_registry%get_field('w')
    p => neko_registry%get_field('p')

    ! initialise file output for debugging
    call tst_fld%init(rp, "testing_refine", 4)
    call tst_fld%fields%assign_to_ptr(1, p)
    call tst_fld%fields%assign_to_ptr(2, u)
    call tst_fld%fields%assign_to_ptr(3, v)
    call tst_fld%fields%assign_to_ptr(4, w)
    select type(file => tst_fld%file_%file_type)
    type is (fld_file_t)
       file%write_mesh = .true.
    end select

  end subroutine user_initialize

  ! User-defined routine called at the end of every time step
  subroutine compute(time)
    type(time_state_t), intent(in) :: time

    type(field_t), pointer :: w

    ! set the w component to zero to avoid any 3D instability
    ! in this quasi-2D flow
    w => neko_registry%get_field("w")
    call field_rzero(w)

  end subroutine compute

  subroutine amr_refine_flag(time, ref_mark, ifrefine)
    type(time_state_t), intent(in) :: time
    integer, dimension(:), intent(inout) :: ref_mark
    logical, intent(inout) :: ifrefine

    if (time%tstep .eq. 4) then
!       if (pe_rank .eq. 0) then
!          ref_mark(:) = amr_flg_h_ref
!       else
!          ref_mark(:) = amr_flg_none
!       end if
       ref_mark(:) = amr_flg_h_ref
       ifrefine = .true.
    else if (time%tstep .eq. 5) then
       ref_mark(:) = amr_flg_h_crs
!       ifrefine = .true.
       ifrefine = .false.
    else
       ref_mark(:) = amr_flg_none
       ifrefine = .false.
    end if

    ! for debugging
    if (ifrefine) then
       t = time%t
       call tst_fld%sample(t)
    end if

  end subroutine amr_refine_flag

  subroutine amr_reconstructing(reconstruct, counter, tstep)
    type(amr_reconstruct_t), intent(inout) :: reconstruct
    integer, intent(in) :: counter, tstep

    ! for debugging
    call tst_fld%sample(t)

  end subroutine amr_reconstructing

end module user
