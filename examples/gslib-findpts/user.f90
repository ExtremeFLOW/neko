! For more details:
!
! HPTS as implemented in Nek5000:
! --> https://github.com/Nek5000/Nek5000/blob/2c81113187a54e3ab0401858b68a7e0dafb1df72/core/postpro.f#L1510
! GSLIB Documentation and variables:
! --> https://github.com/Nek5000/gslib/blob/39d1baae8f4bfebe3ebca6a234dcc8ba1ee5edc7/src/findpts.c#L86

module user
  use neko
  implicit none
  
  !> Log variable
  character(len=LOG_SIZE) :: log_buf ! For logging status

  !> Probe type
  type(probes_t) :: pb

  !> Output variables
  type(file_t) :: fout
  type(matrix_t) :: mat_out

  !> Case IO parameters  
  integer            :: n_fields
  character(len=:), allocatable  :: pts_in_filename
  character(len=:), allocatable  :: pts_out_filename

contains

  !> Register user defined functions (see user_intf.f90)
  subroutine user_setup(u)
    type(user_t), intent(inout) :: u
    u%user_init_modules     => initialize
    u%user_check            => check
    u%user_finalize_modules => finalize
  end subroutine user_setup

  !> Initialize relevant fields
  subroutine initialize(t, u, v, w, p, coef, params)
    implicit none
    real(kind=rp) :: t
    type(field_t), intent(inout) :: u
    type(field_t), intent(inout) :: v
    type(field_t), intent(inout) :: w
    type(field_t), intent(inout) :: p
    type(coef_t), intent(inout) :: coef
    type(json_file), intent(inout) :: params
    
    !> Support variables 
    integer :: i
    type(file_t) :: input_file
    type(matrix_t) :: mat_coords
    type(vector_t) :: vec_header

    !> Read from json how many fields to interpolate
    call json_get(params, 'case.probes.n_fields', n_fields)
    call json_get(params, 'case.probes.pts_in_filename', pts_in_filename) 
    call json_get(params, 'case.probes.pts_out_filename', pts_out_filename) 

    !> Define the input file object
    input_file = file_t(trim(pts_in_filename))

    !> Read file and initialize the probe object
    !! Pre initialize the number of fields according to case file
    pb%n_fields = n_fields
    !! Read the data and initialize/allocate-memory in probe object
    call input_file%read(pb)
    !! Initialize data that depends on the case file in probe object
    call pb%usr_init(t, params)
    !! Perform the set up of gslib_findpts
    call pb%setup(coef)
    !! Find the probes in the mesh. Map from xyz -> rst
    call pb%map(coef)

    !> Write a summary of the probe info
    call pb%show()

    !> Initialize the output
    fout = file_t(trim(pts_out_filename))
    call mat_out%init(pb%n_probes, pb%n_fields)

    !> Write coordinates in output file (imitate nek5000)
    !! Initialize the arrays
    call mat_coords%init(pb%n_probes,3)
    !! Array them as rows
    call transpose(mat_coords%x, pb%n_probes, pb%xyz, 3)
    !! Write the data to the file
    call fout%write(mat_coords)
    !! Free the memory 
    call mat_coords%free

  end subroutine initialize

  subroutine check(t, tstep,u, v, w, p, coef, param)
    real(kind=rp), intent(in) :: t
    integer, intent(in) :: tstep
    type(coef_t), intent(inout) :: coef
    type(field_t), intent(inout) :: u
    type(field_t), intent(inout) :: v
    type(field_t), intent(inout) :: w
    type(field_t), intent(inout) :: p
    type(json_file), intent(inout) :: param
    logical :: write_output = .false.

    integer :: ierr
    integer :: i, il

    !> Interpolate the desired fields
    call pb%interpolate(t,tstep, write_output)
    !! Write if the interpolate function returs write_output=.true.
    if (write_output) then
       call transpose(mat_out%x, pb%n_probes, pb%out_fields, pb%n_fields)
       call fout%write(mat_out, t)
       write_output = .false.
    end if

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
