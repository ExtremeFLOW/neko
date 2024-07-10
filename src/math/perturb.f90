module perturb
  use pcs_f
  use logger
  use num_types
  use comm
  use math
  use device_math, only: device_copy
  use device
  use iso_c_binding
  implicit none

  contains

  subroutine perturb_init_opts_char(opts, inputchar)
    character(len=LOG_SIZE) :: inputchar
    type(pcs_struct), target :: opts

    opts%fpopts_ptr = c_loc(opts%fpopts)

    opts%fpopts%round = CPFLOAT_RND_NE 
    opts%fpopts%flip=CPFLOAT_SOFTERR_NO !> bitflips not
    opts%fpopts%subnormal=CPFLOAT_SUBN_USE !> use subnormal
    !> Uestion whether this is correct, I think this eliminates the impoartance of teh emax and emin...
    opts%fpopts%explim = CPFLOAT_EXPRANGE_TARG !> use exponent from target format
    opts%fpopts%saturation=CPFLOAT_SAT_NO !> if use staturation arithmetic
    opts%fpopts%saturation=CPFLOAT_SAT_USE !> use staturation arithmetic
    opts%fpopts%format = c_null_char
    opts%fpopts%bitseed = c_null_ptr
    opts%fpopts%randseedf = c_null_ptr
    opts%fpopts%randseed= c_null_ptr

    if (trim(inputchar) .eq. 'fp64') then
       opts%oper = PCS_CPFLOAT
       opts%fpopts%precision = 52
    else if (trim(inputchar) .eq. 'fp32') then
       opts%oper = PCS_CPFLOAT
       opts%fpopts%precision = 24 !Bits in the significand + 1.
       opts%fpopts%emax = 127
       opts%fpopts%emin = -126
       opts%fpopts%explim = CPFLOAT_EXPRANGE_TARG !> use emax and emin
    else if (trim(inputchar) .eq. 'fp16') then
       opts%oper = PCS_CPFLOAT
       opts%fpopts%precision = 11 !Bits in the significand + 1.
       opts%fpopts%emax = 15
       opts%fpopts%emin = -14
       opts%fpopts%explim = CPFLOAT_EXPRANGE_TARG !> use emax and emin
    else if (trim(inputchar) .eq. 'e4m3') then
       opts%oper = PCS_CPFLOAT
       opts%fpopts%precision = 4 !Bits in the significand + 1.
       opts%fpopts%emax = 7
       opts%fpopts%emin = -6
       opts%fpopts%explim = CPFLOAT_EXPRANGE_TARG !> use emax and emin
    else if (trim(inputchar) .eq. 'e5m2') then
       opts%oper = PCS_CPFLOAT
       opts%fpopts%precision = 3 !Bits in the significand + 1.
       opts%fpopts%emax = 15
       opts%fpopts%emin = -14
       opts%fpopts%explim = CPFLOAT_EXPRANGE_TARG !> use emax and emin
    else if (trim(inputchar) .eq. 'fp64_stor') then
       opts%oper = PCS_CPFLOAT
       opts%fpopts%precision = 52
       opts%fpopts%explim = CPFLOAT_EXPRANGE_STOR !> use exponent from storage format
    else if (trim(inputchar) .eq. 'fp32_stor') then
       opts%oper = PCS_CPFLOAT
       opts%fpopts%precision = 24 !Bits in the significand + 1.
       opts%fpopts%emax = 127
       opts%fpopts%emin = -126
       opts%fpopts%explim = CPFLOAT_EXPRANGE_STOR !> use exponent from storage format
    else if (trim(inputchar) .eq. 'fp16_stor') then
       opts%oper = PCS_CPFLOAT
       opts%fpopts%precision = 11 !Bits in the significand + 1.
       opts%fpopts%emax = 15
       opts%fpopts%emin = -14
       opts%fpopts%explim = CPFLOAT_EXPRANGE_STOR !> use exponent from storage format
    else if (trim(inputchar) .eq. 'e4m3_stor') then
       opts%oper = PCS_CPFLOAT
       opts%fpopts%precision = 4 !Bits in the significand + 1.
       opts%fpopts%emax = 7
       opts%fpopts%emin = -6
       opts%fpopts%explim = CPFLOAT_EXPRANGE_STOR !> use exponent from storage format
    else if (trim(inputchar) .eq. 'e5m2_stor') then
       opts%oper = PCS_CPFLOAT
       opts%fpopts%precision = 3 !Bits in the significand + 1.
       opts%fpopts%emax = 15
       opts%fpopts%emin = -14
       opts%fpopts%explim = CPFLOAT_EXPRANGE_STOR !> use exponent from storage format
    end if


    if (pe_rank .eq. 0) then
       write(*,*) 'PCS validate', validate_pcs_struct(opts)
       write(*,*) 'PCSSettings operations', opts%oper

       write(*,*) 'PCS rounding',opts%fpopts%round  
       write(*,*)'PCS precision', opts%fpopts%precision !Bits in the significand + 1.
       write(*,*)'PCS emax', opts%fpopts%emax
       write(*,*)'PCS emin', opts%fpopts%emin
       write(*,*) 'PCS flip',opts%fpopts%flip !> bitflips not
       write(*,*) 'PCS subnormal',opts%fpopts%subnormal !> use subnormal
       write(*,*) 'PCS explim ',opts%fpopts%explim  !> use exponent from storage format
       write(*,*) 'PCS saturation',opts%fpopts%saturation !> if use s
    end if
    

  end subroutine perturb_init_opts_char

  subroutine perturb_vector(x, y, n, opts )
    integer, intent(in) :: n
    real(kind=rp), dimension(n), intent(in) :: y  
    real(kind=rp), dimension(n), intent(out) :: x
    type(pcs_struct), target,  intent(inout) :: opts
    integer :: ierr

    if ((opts%fpopts%precision .ne. 52)) then 
        ierr = pcs(x,y,int(n,8),opts)
    else
        call copy(x,y,n) 
    end if

  end subroutine perturb_vector

  subroutine perturb_vector_device(x_d, y_d,temp, n, opts )
    integer, intent(in) :: n
    type(c_ptr), intent(inout) :: x_d, y_d
    type(pcs_struct), intent(in) :: opts
    real(kind=rp), dimension(n), intent(inout) :: temp
    integer :: ierr
    if ((opts%fpopts%precision .ne. 52)) then 
       call device_memcpy(temp, y_d, n, DEVICE_TO_HOST, sync=.true.)
       ierr = pcs(temp, temp,int(n,8),opts)
       call device_memcpy(temp, x_d, n, HOST_TO_DEVICE, sync=.true.)
   else
       call device_copy(x_d,y_d,n)
   end if

  end subroutine perturb_vector_device

end module perturb
