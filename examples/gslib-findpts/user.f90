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

  ! output variables
  type(file_t) :: fout
  type(matrix_t) :: mat_out

contains

  ! Register user defined functions (see user_intf.f90)
  subroutine user_setup(u)
    type(user_t), intent(inout) :: u
    u%user_init_modules     => initialize ! Initialize parameters
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

    type(matrix_t) :: mat_coords

    real(kind=rp) :: delta = 0.04_rp

    ! ******** # of POINTS *******
    integer, parameter :: N = 5
    ! ****************************
    integer, parameter :: nfields = 4
    integer :: handle
    real(kind=rp) :: xyz_raw(3,N)
    integer :: proc(N), elid(N), dist(N), rcode(N)
    real(kind=rp) :: rst_raw(3*N)
    real(kind=rp) :: interpolated_field(nfields,N)

    integer :: i, lx, ly, lz, ierr, nelv, j
    real(kind=rp) :: tol_dist ! Copied from Nek5000 hpts, tolerance for border points
    tol_dist = 5.0d-6 ! Set from Nek5000

    lx = coef%xh%lx
    ly = coef%xh%ly
    lz = coef%xh%lz
    nelv = coef%msh%nelv

    !
    ! Input probes
    !
    xyz_raw(:,1) = (/-0.9,-0.5,-0.5/)
    xyz_raw(:,2) = (/-0.8,-0.5,-0.5/)
    xyz_raw(:,3) = (/-0.7,-0.5,-0.5/)
    xyz_raw(:,4) = (/-0.6,-0.5,-0.5/)
    xyz_raw(:,5) = (/-0.5,-0.5,-0.5/)
!!$    xyz_raw(:,6) = (/-0.3,-0.5,-0.5/)
!!$    xyz_raw(:,7) = (/-0.2,-0.5,-0.5/)
!!$    xyz_raw(:,8) = (/-0.6,-0.2,-0.5/)

    do j = 1,N
       do i = 1,3
          xyz_raw(i,j) = xyz_raw(i,j) + delta
       end do
    end do

    !
    ! Setup, see findpts.c for description of parameters
    !
    call fgslib_findpts_setup(handle, &
         NEKO_COMM, pe_size, &
         coef%msh%gdim, &
         coef%dof%x, coef%dof%y, coef%dof%z, &
         lx, ly, lz, &
         nelv, & ! in Nek: uses "nelt" = # of elts per processor (on the temperature mesh)
         2*lx, 2*ly, 2*lz, & ! Mesh size for bounding box computation
         0.01, lx*ly*lz*nelv, lx*ly*lz*nelv, 128, 5e-13)

    rcode = 0
    elid = 0
    proc = 0

    !
    ! MAPPING from xyz to rst
    !
    call fgslib_findpts(handle, &
         rcode, 1, &
         proc, 1, &
         elid, 1, &
         rst_raw, coef%msh%gdim, &
         dist, 1, &
         xyz_raw(1,1), coef%msh%gdim, &
         xyz_raw(2,1), coef%msh%gdim, &
         xyz_raw(3,1), coef%msh%gdim, N)

    !
    ! Quick check to see who owns what?
    ! Will show, for each point: process, process owner, element owner, error code
    !
    write (*, *) pe_rank, "/", proc, "/" , elid, "/", rcode
!!$ 100 format(((I2," "),A1,5(I2," "),A1,5(I10," "),A1,5(I1," ")))

    !
    ! Final check to see if there are any problems
    !
    do i=1,N
       if (rcode(i) .eq. 1) then
          if (dist(i) .gt. tol_dist) then
             call neko_warning("Point on boundary or outside the mesh!")
          end if
       end if

       if (rcode(i) .eq. 2) call neko_warning("Point not within mesh!")
    end do

    !
    ! INTERPOLATION
    !
    call fgslib_findpts_eval(handle, interpolated_field(1,1), nfields, &
         rcode, 1, &
         proc, 1, &
         elid, 1, &
         rst_raw, coef%msh%gdim, &
         N, u%x(1,1,1,1))

    call fgslib_findpts_eval(handle, interpolated_field(2,1), nfields, &
         rcode, 1, &
         proc, 1, &
         elid, 1, &
         rst_raw, coef%msh%gdim, &
         N, v%x(1,1,1,1))

    call fgslib_findpts_eval(handle, interpolated_field(3,1), nfields, &
         rcode, 1, &
         proc, 1, &
         elid, 1, &
         rst_raw, coef%msh%gdim, &
         N, w%x(1,1,1,1))

    call fgslib_findpts_eval(handle, interpolated_field(4,1), nfields, &
         rcode, 1, &
         proc, 1, &
         elid, 1, &
         rst_raw, coef%msh%gdim, &
         N, p%x(1,1,1,1))

    !
    ! Initialize the output
    !
    fout = file_t("probes_output.csv")
    call mat_out%init(N, nfields)

    !
    ! Write coordinates in output file (imitate nek5000)
    !
    call mat_coords%init(N,3)
    call transpose(mat_coords%x, N, xyz_raw, 3)
    call fout%write(mat_coords)
    call mat_coords%free ! We don't free fout as it will be used in usercheck

    !
    ! Output to file
    !
    call transpose(mat_out%x, N, interpolated_field, nfields)
    call fout%write(mat_out, t)

  end subroutine initialize

  ! Free relevant objects
  subroutine finalize(t, params)
    real(kind=rp) :: t
    type(json_file), intent(inout) :: params

    call fgslib_findpts_free(handle)
    call mat_out%free
    call file_free(fout)

  end subroutine finalize


end module user
