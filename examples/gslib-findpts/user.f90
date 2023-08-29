! For more details:
!
! HPTS as implemented in Nek5000:
! --> https://github.com/Nek5000/Nek5000/blob/2c81113187a54e3ab0401858b68a7e0dafb1df72/core/postpro.f#L1510
! GSLIB Documentation and variables:
! --> https://github.com/Nek5000/gslib/blob/39d1baae8f4bfebe3ebca6a234dcc8ba1ee5edc7/src/findpts.c#L86

module user
  use neko
  implicit none

  character(len=LOG_SIZE) :: log_buf ! For logging status

  type(probes_t) :: pb

  ! output variables
  type(file_t) :: fout
  type(matrix_t) :: mat_out

contains

  ! Register user defined functions (see user_intf.f90)
  subroutine user_setup(u)
    type(user_t), intent(inout) :: u
    u%user_init_modules     => initialize ! Initialize parameters
    u%user_check            => check
    u%user_finalize_modules => finalize ! Finalize
  end subroutine user_setup

  subroutine initialize(t, u, v, w, p, coef, params)
    implicit none

    real(kind=rp) :: t
    type(field_t), intent(inout) :: u
    type(field_t), intent(inout) :: v
    type(field_t), intent(inout) :: w
    type(field_t), intent(inout) :: p
    type(coef_t), intent(inout) :: coef
    type(json_file), intent(inout) :: params

    integer :: i
    type(file_t) :: input_file
    type(matrix_t) :: mat_coords
    type(vector_t) :: vec_header

    input_file = file_t('input.csv')

    call input_file%read(pb)
    call pb%setup(coef) ! basically findpts_setup
    call pb%map(coef) ! map x,y,z -> r,s,t

    call pb%show()
!!$    call pb%debug()

    !
    ! Initialize the output
    !
    fout = file_t("output.csv")
    call mat_out%init(pb%n_probes, pb%n_fields)

    !
    ! Write coordinates in output file (imitate nek5000)
    !
    call mat_coords%init(pb%n_probes,3)
    call transpose(mat_coords%x, pb%n_probes, pb%xyz, 3)
    call fout%write(mat_coords)
    call mat_coords%free

  end subroutine initialize
! usrcheck, this is called at the end of every time step

  subroutine check(t, tstep,u, v, w, p, coef, param)
    real(kind=rp), intent(in) :: t
    integer, intent(in) :: tstep
    type(coef_t), intent(inout) :: coef
    type(field_t), intent(inout) :: u
    type(field_t), intent(inout) :: v
    type(field_t), intent(inout) :: w
    type(field_t), intent(inout) :: p
    type(json_file), intent(inout) :: param

    integer :: ierr

    if (mod(tstep, 5) .ne. 0) return

    !
    ! Interpolation: x, y, z and u
    !
    call pb%interpolate(u, v, w, p)

    !
    ! Output to file
    !
    call transpose(mat_out%x, pb%n_probes, pb%out_fields, pb%n_fields)
    call fout%write(mat_out, t)

  end subroutine check

  ! Free relevant objects
  subroutine finalize(t, params)
    real(kind=rp) :: t
    type(json_file), intent(inout) :: params

    call pb%free
    call mat_out%free
    call file_free(fout)

  end subroutine finalize


end module user
