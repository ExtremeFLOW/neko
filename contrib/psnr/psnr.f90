
module calculate_psnr
  use neko
  implicit none

contains

  subroutine calculate_psnr_sp(compressed, original, size)
    real(kind=sp), allocatable, intent(in) :: compressed(:)
    real(kind=sp), allocatable, intent(in) :: original(:)
    integer, intent(in) :: size
    real(kind=sp) :: minVal, maxVal, peakVal, mse, psnr, sum
    integer :: i

    peakVal = -huge(peakVal)

    minVal = huge(minVal)
    maxVal = -huge(maxVal)

    mse = 0.0_sp
    psnr = 0.0_sp

    do i = 1, size
       minVal = min(minVal, real(compressed(i), sp))
       maxVal = max(maxVal, real(compressed(i), sp))
       sum = real((compressed(i) - original(i)), sp)
       mse = mse + sum * sum
    end do
    mse = mse / real(size, sp)

    peakVal = max(peakVal, maxVal - minVal)
    psnr = 20.0_sp * log10(peakVal / (2.0_sp * sqrt(mse)))

    write(*,*) "MSE, PSNR : ", mse, psnr
  end subroutine calculate_psnr_sp

  subroutine calculate_psnr_dp(original, compressed, size)
    real(kind=dp), allocatable, intent(in) :: original(:)
    real(kind=dp), allocatable, intent(in) :: compressed(:)
    integer, intent(in) :: size
    real(kind=dp) :: minVal, maxVal, peakVal, mse, psnr, sum
    integer :: i

    peakVal = -huge(peakVal)

    minVal = huge(minVal)
    maxVal = -huge(maxVal)

    mse = 0.0_dp
    psnr = 0.0_dp

    do i = 1, size
       minVal = min(minVal, real(original(i), dp))
       maxVal = max(maxVal, real(original(i), dp))
       sum = real((original(i) - compressed(i)), dp)
       mse = mse + sum * sum
    end do
    mse = mse / real(size, dp)

    peakVal = max(peakVal, maxVal - minVal)
    psnr = 20.0_dp * log10(peakVal / (2.0_dp * sqrt(mse)))

    write(*,*) "MSE, PSNR : ", mse, psnr

  end subroutine calculate_psnr_dp

end module calculate_psnr

!> Program to convert calculate MSE and PSNR of compressed to original data from adios2 file format.
program psnr
  use neko
  use calculate_psnr
  implicit none

  character(len=NEKO_FNAME_LEN) :: inputchar, original_fname, compressed_fname
  type(file_t) :: original_file, compressed_file
  type(fld_file_data_t) :: original_data
  type(fld_file_data_t) :: compressed_data
  integer :: argc, file_precision, i, j, n_scalars
  logical :: dp_precision

  argc = command_argument_count()

  if ((argc .lt. 3) .or. (argc .gt. 3)) then
     if (pe_rank .eq. 0) then
        write(*,*) 'Usage: ./psnr compresser_field.bp original_field.bp precision'
        write(*,*) 'Example command: ./psnr compresser_field.bp original_field.bp .false.'
        write(*,*) 'Calculate MSE and PSNR of compressed to original data'
     end if
     stop
  end if

  call neko_init

  call get_command_argument(1, inputchar)
  read(inputchar, fmt='(A)') original_fname
  call get_command_argument(2, inputchar)
  read(inputchar, fmt='(A)') compressed_fname
  call get_command_argument(3, inputchar)
  read(inputchar, *) dp_precision

  if (dp_precision) then
     file_precision = dp
  else
     file_precision = sp
  end if

  if (pe_rank .eq. 0) write(*,*) 'Reading file:', trim(original_fname)
  original_file = file_t(trim(original_fname),precision=file_precision)
  call original_data%init()
  call original_file%read(original_data)

  if (pe_rank .eq. 0) write(*,*) 'Reading file:', trim(compressed_fname)
  compressed_file = file_t(trim(compressed_fname),precision=file_precision)
  call compressed_data%init()
  call compressed_file%read(compressed_data)

  do i = 1, compressed_data%meta_nsamples-1
     if (pe_rank .eq. 0) write(*,*) 'Reading file:', i+1, trim(original_fname)
     call original_file%read(original_data)

     if (pe_rank .eq. 0) write(*,*) 'Reading file:', i+1, trim(compressed_fname)
     call compressed_file%read(compressed_data)

     n_scalars = compressed_data%size() - 5
     if (compressed_data%u%n .eq. original_data%u%n) then
        write(*,*) "Error in velocity components"
        call calculate_psnr_dp(original_data%u%x, compressed_data%u%x, &
             compressed_data%u%n)
        call calculate_psnr_dp(original_data%v%x, compressed_data%v%x, &
             compressed_data%v%n)
        call calculate_psnr_dp(original_data%w%x, compressed_data%w%x, &
             compressed_data%w%n)
     end if
     if (compressed_data%u%n .eq. original_data%p%n) then
        write(*,*) "Error in pressure"
        call calculate_psnr_dp(original_data%p%x, compressed_data%p%x, &
             compressed_data%p%n)
     end if
     if (compressed_data%u%n .eq. original_data%t%n) then
        write(*,*) "Error in temperature"
        call calculate_psnr_dp(original_data%t%x, compressed_data%t%x, &
             compressed_data%t%n)
     end if
     do j = 1, n_scalars
        if (compressed_data%s(j)%n .eq. original_data%s(j)%n) then
           write(*,*) "Error in scalar ", j
           call calculate_psnr_dp(original_data%s(j)%x, compressed_data%s(j)%x, &
                compressed_data%s(j)%n)
        end if
     end do

  end do

  if (pe_rank .eq. 0) write(*,*) 'Done'

  call neko_finalize

end program psnr

