module tensor_cpu
  use num_types, only : rp
  use mxm_wrapper, only : mxm
  implicit none
  private

  public :: tnsr2d_el_cpu, tnsr3d_el_cpu, tnsr3d_cpu, tnsr1_3d_cpu

contains

  !> Two-dimensional tensor product on a single element.
  !! @details Applies \f$ v = (A \otimes B^T) u \f$: the one-dimensional
  !! operator A along the first reference direction and B along the second,
  !! taking the element from order @a nu to order @a nv. Routed through `mxm`,
  !! so it uses whichever matrix-matrix backend the build selected.
  !! @param v Result, nv*nv values.
  !! @param nv Output polynomial order.
  !! @param u Input, nu*nu values.
  !! @param nu Input polynomial order.
  !! @param A Operator along the first reference direction, (nv, nu).
  !! @param Bt Transposed operator along the second direction, (nu, nv).
  subroutine tnsr2d_el_cpu(v, nv, u, nu, A, Bt)
    integer, intent(in) :: nv, nu
    real(kind=rp), intent(inout) :: v(nv*nv), u(nu*nu)
    real(kind=rp), intent(inout) :: A(nv, nu), Bt(nu, nv)
    real(kind=rp) :: work(0:nu**2*nv)

    call mxm(A, nv, u, nu, work, nu)
    call mxm(work, nv, Bt, nu, v, nv)

  end subroutine tnsr2d_el_cpu

  !> Three-dimensional tensor product on a single element.
  !! @details Applies \f$ v = (A \otimes B^T \otimes C^T) u \f$, one
  !! one-dimensional operator per reference direction, taking the element from
  !! order @a nu to order @a nv. Dispatches on the shape: @a nv = @a nu selects
  !! the square kernels (filtering and SVV, where the order is unchanged),
  !! @a nv = 1 the single-point kernels (evaluating an element at one location,
  !! as local_interpolation does), and anything else the generic kernel.
  !! @note Contains no OpenMP and is marked `!OCL SERIAL`, because callers
  !! invoke it per element from inside their own parallel region.
  !! @param v Result, nv*nv*nv values.
  !! @param nv Output polynomial order.
  !! @param u Input, nu*nu*nu values.
  !! @param nu Input polynomial order.
  !! @param A Operator along the first reference direction, (nv, nu).
  !! @param Bt Transposed operator along the second direction, (nu, nv).
  !! @param Ct Transposed operator along the third direction, (nu, nv).
  subroutine tnsr3d_el_cpu(v, nv, u, nu, A, Bt, Ct)
    integer, intent(in) :: nv, nu
    real(kind=rp), intent(inout) :: v(nv*nv*nv), u(nu*nu*nu)
    real(kind=rp), intent(inout) :: A(nv, nu), Bt(nu, nv), Ct(nu, nv)

    if (nv .eq. nu) then
       select case (nv)
       case (14)
          call tnsr3d_el_n14_cpu(v, u, A, Bt, Ct)
       case (13)
          call tnsr3d_el_n13_cpu(v, u, A, Bt, Ct)
       case (12)
          call tnsr3d_el_n12_cpu(v, u, A, Bt, Ct)
       case (11)
          call tnsr3d_el_n11_cpu(v, u, A, Bt, Ct)
       case (10)
          call tnsr3d_el_n10_cpu(v, u, A, Bt, Ct)
       case (9)
          call tnsr3d_el_n9_cpu(v, u, A, Bt, Ct)
       case (8)
          call tnsr3d_el_n8_cpu(v, u, A, Bt, Ct)
       case (7)
          call tnsr3d_el_n7_cpu(v, u, A, Bt, Ct)
       case (6)
          call tnsr3d_el_n6_cpu(v, u, A, Bt, Ct)
       case (5)
          call tnsr3d_el_n5_cpu(v, u, A, Bt, Ct)
       case (4)
          call tnsr3d_el_n4_cpu(v, u, A, Bt, Ct)
       case (3)
          call tnsr3d_el_n3_cpu(v, u, A, Bt, Ct)
       case (2)
          call tnsr3d_el_n2_cpu(v, u, A, Bt, Ct)
       case default
          call tnsr3d_el_n_cpu(v, u, A, Bt, Ct, nv)
       end select
    else if (nv .eq. 1) then
       select case (nu)
       case (4)
          call tnsr3d_el_1_4_cpu(v, u, A, Bt, Ct)
       case (6)
          call tnsr3d_el_1_6_cpu(v, u, A, Bt, Ct)
       case (8)
          call tnsr3d_el_1_8_cpu(v, u, A, Bt, Ct)
       case (10)
          call tnsr3d_el_1_10_cpu(v, u, A, Bt, Ct)
       case (12)
          call tnsr3d_el_1_12_cpu(v, u, A, Bt, Ct)
       case default
          call tnsr3d_el_1_nu_cpu(v, u, nu, A, Bt, Ct)
       end select
    else
       call tnsr3d_el_nvnu_cpu(v, nv, u, nu, A, Bt, Ct)
    end if

  end subroutine tnsr3d_el_cpu

  !> Single-element tensor product for arbitrary @a nu and @a nv.
  !! @details The fallback of tnsr3d_el_cpu, used when neither the square nor
  !! the single-point shape applies. Contracts one direction at a time through
  !! two work buffers.
  !! @param v Result, nv*nv*nv values.
  !! @param nv Output polynomial order.
  !! @param u Input, nu*nu*nu values.
  !! @param nu Input polynomial order.
  !! @param v Result, nv*nv*nv values.
  !! @param u Input, nu*nu*nu values.
  !! @param A Operator along the first reference direction, (nv, nu).
  !! @param Bt Transposed operator along the second direction, (nu, nv).
  !! @param Ct Transposed operator along the third direction, (nu, nv).
  !OCL SERIAL
  subroutine tnsr3d_el_nvnu_cpu(v, nv, u, nu, A, Bt, Ct)
    integer, intent(in) :: nv, nu
    real(kind=rp), intent(inout) :: v(nv*nv*nv), u(nu*nu*nu)
    real(kind=rp), intent(inout) :: A(nv, nu), Bt(nu, nv), Ct(nu, nv)
    real(kind=rp) :: work(nu**2*nv), work2(nu*nv**2)
    real(kind=rp) :: tmp
    integer :: i, j, k, l, nunu, nvnu, nvnv
    integer :: ii, jj
    nvnu = nv * nu
    nunu = nu * nu
    nvnv = nv * nv

    do j = 1, nunu
       do i = 1, nv
          ii = i + nv * (j - 1)
          tmp = 0.0_rp
          do k = 1, nu
             tmp = tmp + A(i,k) * u(k + nu * (j - 1))
          end do
          work(ii) = tmp
       end do
    end do

    do i = 1, nu
       do j = 1, nv
          do l = 1, nv
             ii = l + nv * (j - 1) + nvnv * (i - 1)
             tmp = 0.0_rp
             do k = 1, nu
                jj = l + nv * (k - 1) + nvnu * (i - 1)
                tmp = tmp + work(jj) * Bt(k,j)
             end do
             work2(ii) = tmp
          end do
       end do
    end do

    do j = 1, nv
       do i = 1, nvnv
          jj = i + nvnv * (j - 1)
          tmp = 0.0_rp
          do k = 1, nu
             ii = i + nvnv * (k - 1)
             tmp = tmp + work2(ii) * Ct(k, j)
          end do
          v(jj) = tmp
       end do
    end do

  end subroutine tnsr3d_el_nvnu_cpu

  !> Single-element evaluation at one point, for arbitrary @a nu.
  !! @details The @a nv = 1 case of tnsr3d_el_cpu: each operator is a single
  !! row, so the three contractions collapse the element to one value. Used to
  !! interpolate a field to a scattered point inside an element.
  !! @param v Result, a single value.
  !! @param u Input, nu*nu*nu values.
  !! @param nu Input polynomial order.
  !! @param A Operator along the first reference direction, (1, nu).
  !! @param Bt Transposed operator along the second direction, (nu, 1).
  !! @param Ct Transposed operator along the third direction, (nu, 1).
  !OCL SERIAL
  subroutine tnsr3d_el_1_nu_cpu(v, u, nu, A, Bt, Ct)
    integer, intent(in) :: nu
    real(kind=rp), intent(inout) :: v(1)
    real(kind=rp), intent(in) :: u(nu*nu*nu)
    real(kind=rp), intent(in) :: A(1, nu), Bt(nu, 1), Ct(nu, 1)
    real(kind=rp) :: work(nu**2), work2(nu)
    real(kind=rp) :: tmp
    integer :: i, j, k, nunu
    integer :: jj
    nunu = nu * nu

    do j = 1, nunu
       tmp = 0.0_rp
       do k = 1, nu
          tmp = tmp + A(1,k) * u(k + nu * (j - 1))
       end do
       work(j) = tmp
    end do

    do i = 1, nu
       tmp = 0.0_rp
       do k = 1, nu
          jj = k + nu * (i - 1)
          tmp = tmp + work(jj) * Bt(k,1)
       end do
       work2(i) = tmp
    end do

    tmp = 0.0_rp
    do k = 1, nu
       tmp = tmp + work2(k) * Ct(k, 1)
    end do
    v(1) = tmp

  end subroutine tnsr3d_el_1_nu_cpu

  !> Single-element evaluation at one point, specialised for @a nu = 4.
  !! @details Unrolled form of tnsr3d_el_1_nu_cpu; see tnsr3d_el_cpu for the
  !! operation and the argument meanings.
  !OCL SERIAL
  subroutine tnsr3d_el_1_4_cpu(v, u, A, Bt, Ct)
    integer, parameter :: n = 4
    integer, parameter :: nn = n**2
    real(kind=rp), intent(inout) :: v(1)
    real(kind=rp), intent(in) :: u(n*n*n)
    real(kind=rp), intent(in) :: A(1,n), Bt(n,1), Ct(n,1)
    real(kind=rp) :: work(n**2), work2(n)
    integer :: i, j

    do j = 1, nn
       work(j) = A(1,1) * u(1 + n * (j - 1)) &
            + A(1,2) * u(2 + n * (j - 1)) &
            + A(1,3) * u(3 + n * (j - 1)) &
            + A(1,4) * u(4 + n * (j - 1))
    end do

    do i = 1, n
       work2(i) = work(1 + n * (i - 1)) * Bt(1,1) &
            + work(2 + n * (i - 1)) * Bt(2,1) &
            + work(3 + n * (i - 1)) * Bt(3,1) &
            + work(4 + n * (i - 1)) * Bt(4,1)
    end do

    v(1) = work2(1) * Ct(1, 1) &
         + work2(2) * Ct(2, 1) &
         + work2(3) * Ct(3, 1) &
         + work2(4) * Ct(4, 1)

  end subroutine tnsr3d_el_1_4_cpu

  !> Single-element evaluation at one point, specialised for @a nu = 6.
  !! @details Unrolled form of tnsr3d_el_1_nu_cpu; see tnsr3d_el_cpu for the
  !! operation and the argument meanings.
  !OCL SERIAL
  subroutine tnsr3d_el_1_6_cpu(v, u, A, Bt, Ct)
    integer, parameter :: n = 6
    integer, parameter :: nn = n**2
    real(kind=rp), intent(inout) :: v(1)
    real(kind=rp), intent(in) :: u(n*n*n)
    real(kind=rp), intent(in) :: A(1,n), Bt(n,1), Ct(n,1)
    real(kind=rp) :: work(n**2), work2(n)
    integer :: i, j

    do j = 1, nn
       work(j) = A(1,1) * u(1 + n * (j - 1)) &
            + A(1,2) * u(2 + n * (j - 1)) &
            + A(1,3) * u(3 + n * (j - 1)) &
            + A(1,4) * u(4 + n * (j - 1)) &
            + A(1,5) * u(5 + n * (j - 1)) &
            + A(1,6) * u(6 + n * (j - 1))
    end do

    do i = 1, n
       work2(i) = work(1 + n * (i - 1)) * Bt(1,1) &
            + work(2 + n * (i - 1)) * Bt(2,1) &
            + work(3 + n * (i - 1)) * Bt(3,1) &
            + work(4 + n * (i - 1)) * Bt(4,1) &
            + work(5 + n * (i - 1)) * Bt(5,1) &
            + work(6 + n * (i - 1)) * Bt(6,1)
    end do

    v(1) = work2(1) * Ct(1, 1) &
         + work2(2) * Ct(2, 1) &
         + work2(3) * Ct(3, 1) &
         + work2(4) * Ct(4, 1) &
         + work2(5) * Ct(5, 1) &
         + work2(6) * Ct(6, 1)

  end subroutine tnsr3d_el_1_6_cpu

  !> Single-element evaluation at one point, specialised for @a nu = 8.
  !! @details Unrolled form of tnsr3d_el_1_nu_cpu; see tnsr3d_el_cpu for the
  !! operation and the argument meanings.
  !OCL SERIAL
  subroutine tnsr3d_el_1_8_cpu(v, u, A, Bt, Ct)
    integer, parameter :: n = 8
    integer, parameter :: nn = n**2
    real(kind=rp), intent(inout) :: v(1)
    real(kind=rp), intent(in) :: u(n*n*n)
    real(kind=rp), intent(in) :: A(1,n), Bt(n,1), Ct(n,1)
    real(kind=rp) :: work(n**2), work2(n)
    integer :: i, j

    do j = 1, nn
       work(j) = A(1,1) * u(1 + n * (j - 1)) &
            + A(1,2) * u(2 + n * (j - 1)) &
            + A(1,3) * u(3 + n * (j - 1)) &
            + A(1,4) * u(4 + n * (j - 1)) &
            + A(1,5) * u(5 + n * (j - 1)) &
            + A(1,6) * u(6 + n * (j - 1)) &
            + A(1,7) * u(7 + n * (j - 1)) &
            + A(1,8) * u(8 + n * (j - 1))
    end do

    do i = 1, n
       work2(i) = work(1 + n * (i - 1)) * Bt(1,1) &
            + work(2 + n * (i - 1)) * Bt(2,1) &
            + work(3 + n * (i - 1)) * Bt(3,1) &
            + work(4 + n * (i - 1)) * Bt(4,1) &
            + work(5 + n * (i - 1)) * Bt(5,1) &
            + work(6 + n * (i - 1)) * Bt(6,1) &
            + work(7 + n * (i - 1)) * Bt(7,1) &
            + work(8 + n * (i - 1)) * Bt(8,1)
    end do

    v(1) = work2(1) * Ct(1, 1) &
         + work2(2) * Ct(2, 1) &
         + work2(3) * Ct(3, 1) &
         + work2(4) * Ct(4, 1) &
         + work2(5) * Ct(5, 1) &
         + work2(6) * Ct(6, 1) &
         + work2(7) * Ct(7, 1) &
         + work2(8) * Ct(8, 1)


  end subroutine tnsr3d_el_1_8_cpu

  !> Single-element evaluation at one point, specialised for @a nu = 10.
  !! @details Unrolled form of tnsr3d_el_1_nu_cpu; see tnsr3d_el_cpu for the
  !! operation and the argument meanings.
  !OCL SERIAL
  subroutine tnsr3d_el_1_10_cpu(v, u, A, Bt, Ct)
    integer, parameter :: n = 10
    integer, parameter :: nn = n**2
    real(kind=rp), intent(inout) :: v(1)
    real(kind=rp), intent(in) :: u(n*n*n)
    real(kind=rp), intent(in) :: A(1,n), Bt(n,1), Ct(n,1)
    real(kind=rp) :: work(n**2), work2(n)
    integer :: i, j

    do j = 1, nn
       work(j) = A(1,1) * u(1 + n * (j - 1)) &
            + A(1,2) * u(2 + n * (j - 1)) &
            + A(1,3) * u(3 + n * (j - 1)) &
            + A(1,4) * u(4 + n * (j - 1)) &
            + A(1,5) * u(5 + n * (j - 1)) &
            + A(1,6) * u(6 + n * (j - 1)) &
            + A(1,7) * u(7 + n * (j - 1)) &
            + A(1,8) * u(8 + n * (j - 1)) &
            + A(1,9) * u(9 + n * (j - 1)) &
            + A(1,10) * u(10 + n * (j - 1))
    end do

    do i = 1, n
       work2(i) = work(1 + n * (i - 1)) * Bt(1,1) &
            + work(2 + n * (i - 1)) * Bt(2,1) &
            + work(3 + n * (i - 1)) * Bt(3,1) &
            + work(4 + n * (i - 1)) * Bt(4,1) &
            + work(5 + n * (i - 1)) * Bt(5,1) &
            + work(6 + n * (i - 1)) * Bt(6,1) &
            + work(7 + n * (i - 1)) * Bt(7,1) &
            + work(8 + n * (i - 1)) * Bt(8,1) &
            + work(9 + n * (i - 1)) * Bt(9,1) &
            + work(10 + n * (i - 1)) * Bt(10,1)
    end do

    v(1) = work2(1) * Ct(1, 1) &
         + work2(2) * Ct(2, 1) &
         + work2(3) * Ct(3, 1) &
         + work2(4) * Ct(4, 1) &
         + work2(5) * Ct(5, 1) &
         + work2(6) * Ct(6, 1) &
         + work2(7) * Ct(7, 1) &
         + work2(8) * Ct(8, 1) &
         + work2(9) * Ct(9, 1) &
         + work2(10) * Ct(10, 1)

  end subroutine tnsr3d_el_1_10_cpu

  !> Single-element evaluation at one point, specialised for @a nu = 12.
  !! @details Unrolled form of tnsr3d_el_1_nu_cpu; see tnsr3d_el_cpu for the
  !! operation and the argument meanings.
  !OCL SERIAL
  subroutine tnsr3d_el_1_12_cpu(v, u, A, Bt, Ct)
    integer, parameter :: n = 12
    integer, parameter :: nn = n**2
    real(kind=rp), intent(inout) :: v(1)
    real(kind=rp), intent(in) :: u(n*n*n)
    real(kind=rp), intent(in) :: A(1,n), Bt(n,1), Ct(n,1)
    real(kind=rp) :: work(n**2), work2(n)
    integer :: i, j

    do j = 1, nn
       work(j) = A(1,1) * u(1 + n * (j - 1)) &
            + A(1,2) * u(2 + n * (j - 1)) &
            + A(1,3) * u(3 + n * (j - 1)) &
            + A(1,4) * u(4 + n * (j - 1)) &
            + A(1,5) * u(5 + n * (j - 1)) &
            + A(1,6) * u(6 + n * (j - 1)) &
            + A(1,7) * u(7 + n * (j - 1)) &
            + A(1,8) * u(8 + n * (j - 1)) &
            + A(1,9) * u(9 + n * (j - 1)) &
            + A(1,10) * u(10 + n * (j - 1)) &
            + A(1,11) * u(11 + n * (j - 1)) &
            + A(1,12) * u(12 + n * (j - 1))
    end do

    do i = 1, n
       work2(i) = work(1 + n * (i - 1)) * Bt(1,1) &
            + work(2 + n * (i - 1)) * Bt(2,1) &
            + work(3 + n * (i - 1)) * Bt(3,1) &
            + work(4 + n * (i - 1)) * Bt(4,1) &
            + work(5 + n * (i - 1)) * Bt(5,1) &
            + work(6 + n * (i - 1)) * Bt(6,1) &
            + work(7 + n * (i - 1)) * Bt(7,1) &
            + work(8 + n * (i - 1)) * Bt(8,1) &
            + work(9 + n * (i - 1)) * Bt(9,1) &
            + work(10 + n * (i - 1)) * Bt(10,1) &
            + work(11 + n * (i - 1)) * Bt(11,1) &
            + work(12 + n * (i - 1)) * Bt(12,1)
    end do

    v(1) = work2(1) * Ct(1, 1) &
         + work2(2) * Ct(2, 1) &
         + work2(3) * Ct(3, 1) &
         + work2(4) * Ct(4, 1) &
         + work2(5) * Ct(5, 1) &
         + work2(6) * Ct(6, 1) &
         + work2(7) * Ct(7, 1) &
         + work2(8) * Ct(8, 1) &
         + work2(9) * Ct(9, 1) &
         + work2(10) * Ct(10, 1) &
         + work2(11) * Ct(11, 1) &
         + work2(12) * Ct(12, 1)

  end subroutine tnsr3d_el_1_12_cpu

  !> Single-element tensor product with @a nv = @a nu, for arbitrary order.
  !! @details The order-preserving case of tnsr3d_el_cpu, which is what
  !! filtering and the one-sided SVV operators need. The fallback for orders
  !! without an unrolled kernel below.
  !! @param v Result, n*n*n values.
  !! @param u Input, n*n*n values.
  !! @param A Operator along the first reference direction, (n, n).
  !! @param Bt Transposed operator along the second direction, (n, n).
  !! @param Ct Transposed operator along the third direction, (n, n).
  !! @param n Polynomial order, the same on input and output.
  !OCL SERIAL
  subroutine tnsr3d_el_n_cpu(v, u, A, Bt, Ct, n)
    integer, intent(in) :: n
    real(kind=rp), intent(inout) :: v(n*n*n), u(n*n*n)
    real(kind=rp), intent(inout) :: A(n,n), Bt(n,n), Ct(n,n)
    real(kind=rp) :: work(n**3), work2(n**3), tmp
    integer :: i, j, l, k
    integer :: ii, jj, nn

    nn = n**2

    do j = 1, nn
       do i = 1, n
          ii = i + n * (j - 1)
          tmp = 0.0_rp
          do k = 1, n
             tmp = tmp + A(i,k) * u(k + n * (j - 1))
          end do
          work(ii) = tmp
       end do
    end do

    do i = 1, n
       do j = 1, n
          do l = 1, n
             ii = l + n * (j - 1) + nn * (i - 1)
             tmp = 0.0_rp
             do k = 1, n
                tmp = tmp + work(l + n * (k - 1) + nn * (i - 1)) * Bt(k,j)
             end do
             work2(ii) = tmp
          end do
       end do
    end do

    do j = 1, n
       do i = 1, nn
          jj = i + nn * (j - 1)
          tmp = 0.0_rp
          do k = 1, n
             tmp = tmp + work2(i + nn * (k - 1)) * Ct(k, j)
          end do
          v(jj) = tmp
       end do
    end do

  end subroutine tnsr3d_el_n_cpu

  !> Single-element tensor product specialised for @a nv = @a nu = 14.
  !! @details Unrolled form of tnsr3d_el_n_cpu; see tnsr3d_el_cpu for the
  !! operation and the argument meanings.
  !OCL SERIAL
  subroutine tnsr3d_el_n14_cpu(v, u, A, Bt, Ct)
    integer, parameter :: n = 14
    integer, parameter :: nn = n**2
    real(kind=rp), intent(inout) :: v(n*n*n), u(n*n*n)
    real(kind=rp), intent(inout) :: A(n,n), Bt(n,n), Ct(n,n)
    real(kind=rp) :: work(n**3), work2(n**3)
    integer :: i, j, l
    integer :: ii, jj

    do j = 1, nn
       do i = 1, n
          ii = i + n * (j - 1)
          work(ii) = A(i,1) * u(1 + n * (j - 1)) &
               + A(i,2) * u(2 + n * (j - 1)) &
               + A(i,3) * u(3 + n * (j - 1)) &
               + A(i,4) * u(4 + n * (j - 1)) &
               + A(i,5) * u(5 + n * (j - 1)) &
               + A(i,6) * u(6 + n * (j - 1)) &
               + A(i,7) * u(7 + n * (j - 1)) &
               + A(i,8) * u(8 + n * (j - 1)) &
               + A(i,9) * u(9 + n * (j - 1)) &
               + A(i,10) * u(10 + n * (j - 1)) &
               + A(i,11) * u(11 + n * (j - 1)) &
               + A(i,12) * u(12 + n * (j - 1)) &
               + A(i,13) * u(13 + n * (j - 1)) &
               + A(i,14) * u(14 + n * (j - 1))
       end do
    end do

    do i = 1, n
       do j = 1, n
          do l = 1, n
             ii = l + n * (j - 1) + nn * (i - 1)
             work2(ii) = work(l + n * (1 - 1) + nn * (i - 1)) * Bt(1,j) &
                  + work(l + n * (2 - 1) + nn * (i - 1)) * Bt(2,j) &
                  + work(l + n * (3 - 1) + nn * (i - 1)) * Bt(3,j) &
                  + work(l + n * (4 - 1) + nn * (i - 1)) * Bt(4,j) &
                  + work(l + n * (5 - 1) + nn * (i - 1)) * Bt(5,j) &
                  + work(l + n * (6 - 1) + nn * (i - 1)) * Bt(6,j) &
                  + work(l + n * (7 - 1) + nn * (i - 1)) * Bt(7,j) &
                  + work(l + n * (8 - 1) + nn * (i - 1)) * Bt(8,j) &
                  + work(l + n * (9 - 1) + nn * (i - 1)) * Bt(9,j) &
                  + work(l + n * (10 - 1) + nn * (i - 1)) * Bt(10,j) &
                  + work(l + n * (11 - 1) + nn * (i - 1)) * Bt(11,j) &
                  + work(l + n * (12 - 1) + nn * (i - 1)) * Bt(12,j) &
                  + work(l + n * (13 - 1) + nn * (i - 1)) * Bt(13,j) &
                  + work(l + n * (14 - 1) + nn * (i - 1)) * Bt(14,j)
          end do
       end do
    end do

    do j = 1, n
       do i = 1, nn
          jj = i + nn * (j - 1)
          v(jj) = work2(i + nn * (1 - 1)) * Ct(1, j) &
               + work2(i + nn * (2 - 1)) * Ct(2, j) &
               + work2(i + nn * (3 - 1)) * Ct(3, j) &
               + work2(i + nn * (4 - 1)) * Ct(4, j) &
               + work2(i + nn * (5 - 1)) * Ct(5, j) &
               + work2(i + nn * (6 - 1)) * Ct(6, j) &
               + work2(i + nn * (7 - 1)) * Ct(7, j) &
               + work2(i + nn * (8 - 1)) * Ct(8, j) &
               + work2(i + nn * (9 - 1)) * Ct(9, j) &
               + work2(i + nn * (10 - 1)) * Ct(10, j) &
               + work2(i + nn * (11 - 1)) * Ct(11, j) &
               + work2(i + nn * (12 - 1)) * Ct(12, j) &
               + work2(i + nn * (13 - 1)) * Ct(13, j) &
               + work2(i + nn * (14 - 1)) * Ct(14, j)
       end do
    end do

  end subroutine tnsr3d_el_n14_cpu

  !> Single-element tensor product specialised for @a nv = @a nu = 13.
  !! @details Unrolled form of tnsr3d_el_n_cpu; see tnsr3d_el_cpu for the
  !! operation and the argument meanings.
  !OCL SERIAL
  subroutine tnsr3d_el_n13_cpu(v, u, A, Bt, Ct)
    integer, parameter :: n = 13
    integer, parameter :: nn = n**2
    real(kind=rp), intent(inout) :: v(n*n*n), u(n*n*n)
    real(kind=rp), intent(inout) :: A(n,n), Bt(n,n), Ct(n,n)
    real(kind=rp) :: work(n**3), work2(n**3)
    integer :: i, j, l
    integer :: ii, jj

    do j = 1, nn
       do i = 1, n
          ii = i + n * (j - 1)
          work(ii) = A(i,1) * u(1 + n * (j - 1)) &
               + A(i,2) * u(2 + n * (j - 1)) &
               + A(i,3) * u(3 + n * (j - 1)) &
               + A(i,4) * u(4 + n * (j - 1)) &
               + A(i,5) * u(5 + n * (j - 1)) &
               + A(i,6) * u(6 + n * (j - 1)) &
               + A(i,7) * u(7 + n * (j - 1)) &
               + A(i,8) * u(8 + n * (j - 1)) &
               + A(i,9) * u(9 + n * (j - 1)) &
               + A(i,10) * u(10 + n * (j - 1)) &
               + A(i,11) * u(11 + n * (j - 1)) &
               + A(i,12) * u(12 + n * (j - 1)) &
               + A(i,13) * u(13 + n * (j - 1))
       end do
    end do

    do i = 1, n
       do j = 1, n
          do l = 1, n
             ii = l + n * (j - 1) + nn * (i - 1)
             work2(ii) = work(l + n * (1 - 1) + nn * (i - 1)) * Bt(1,j) &
                  + work(l + n * (2 - 1) + nn * (i - 1)) * Bt(2,j) &
                  + work(l + n * (3 - 1) + nn * (i - 1)) * Bt(3,j) &
                  + work(l + n * (4 - 1) + nn * (i - 1)) * Bt(4,j) &
                  + work(l + n * (5 - 1) + nn * (i - 1)) * Bt(5,j) &
                  + work(l + n * (6 - 1) + nn * (i - 1)) * Bt(6,j) &
                  + work(l + n * (7 - 1) + nn * (i - 1)) * Bt(7,j) &
                  + work(l + n * (8 - 1) + nn * (i - 1)) * Bt(8,j) &
                  + work(l + n * (9 - 1) + nn * (i - 1)) * Bt(9,j) &
                  + work(l + n * (10 - 1) + nn * (i - 1)) * Bt(10,j) &
                  + work(l + n * (11 - 1) + nn * (i - 1)) * Bt(11,j) &
                  + work(l + n * (12 - 1) + nn * (i - 1)) * Bt(12,j) &
                  + work(l + n * (13 - 1) + nn * (i - 1)) * Bt(13,j)
          end do
       end do
    end do

    do j = 1, n
       do i = 1, nn
          jj = i + nn * (j - 1)
          v(jj) = work2(i + nn * (1 - 1)) * Ct(1, j) &
               + work2(i + nn * (2 - 1)) * Ct(2, j) &
               + work2(i + nn * (3 - 1)) * Ct(3, j) &
               + work2(i + nn * (4 - 1)) * Ct(4, j) &
               + work2(i + nn * (5 - 1)) * Ct(5, j) &
               + work2(i + nn * (6 - 1)) * Ct(6, j) &
               + work2(i + nn * (7 - 1)) * Ct(7, j) &
               + work2(i + nn * (8 - 1)) * Ct(8, j) &
               + work2(i + nn * (9 - 1)) * Ct(9, j) &
               + work2(i + nn * (10 - 1)) * Ct(10, j) &
               + work2(i + nn * (11 - 1)) * Ct(11, j) &
               + work2(i + nn * (12 - 1)) * Ct(12, j) &
               + work2(i + nn * (13 - 1)) * Ct(13, j)
       end do
    end do

  end subroutine tnsr3d_el_n13_cpu

  !> Single-element tensor product specialised for @a nv = @a nu = 12.
  !! @details Unrolled form of tnsr3d_el_n_cpu; see tnsr3d_el_cpu for the
  !! operation and the argument meanings.
  !OCL SERIAL
  subroutine tnsr3d_el_n12_cpu(v, u, A, Bt, Ct)
    integer, parameter :: n = 12
    integer, parameter :: nn = n**2
    real(kind=rp), intent(inout) :: v(n*n*n), u(n*n*n)
    real(kind=rp), intent(inout) :: A(n,n), Bt(n,n), Ct(n,n)
    real(kind=rp) :: work(n**3), work2(n**3)
    integer :: i, j, l
    integer :: ii, jj

    do j = 1, nn
       do i = 1, n
          ii = i + n * (j - 1)
          work(ii) = A(i,1) * u(1 + n * (j - 1)) &
               + A(i,2) * u(2 + n * (j - 1)) &
               + A(i,3) * u(3 + n * (j - 1)) &
               + A(i,4) * u(4 + n * (j - 1)) &
               + A(i,5) * u(5 + n * (j - 1)) &
               + A(i,6) * u(6 + n * (j - 1)) &
               + A(i,7) * u(7 + n * (j - 1)) &
               + A(i,8) * u(8 + n * (j - 1)) &
               + A(i,9) * u(9 + n * (j - 1)) &
               + A(i,10) * u(10 + n * (j - 1)) &
               + A(i,11) * u(11 + n * (j - 1)) &
               + A(i,12) * u(12 + n * (j - 1))
       end do
    end do

    do i = 1, n
       do j = 1, n
          do l = 1, n
             ii = l + n * (j - 1) + nn * (i - 1)
             work2(ii) = work(l + n * (1 - 1) + nn * (i - 1)) * Bt(1,j) &
                  + work(l + n * (2 - 1) + nn * (i - 1)) * Bt(2,j) &
                  + work(l + n * (3 - 1) + nn * (i - 1)) * Bt(3,j) &
                  + work(l + n * (4 - 1) + nn * (i - 1)) * Bt(4,j) &
                  + work(l + n * (5 - 1) + nn * (i - 1)) * Bt(5,j) &
                  + work(l + n * (6 - 1) + nn * (i - 1)) * Bt(6,j) &
                  + work(l + n * (7 - 1) + nn * (i - 1)) * Bt(7,j) &
                  + work(l + n * (8 - 1) + nn * (i - 1)) * Bt(8,j) &
                  + work(l + n * (9 - 1) + nn * (i - 1)) * Bt(9,j) &
                  + work(l + n * (10 - 1) + nn * (i - 1)) * Bt(10,j) &
                  + work(l + n * (11 - 1) + nn * (i - 1)) * Bt(11,j) &
                  + work(l + n * (12 - 1) + nn * (i - 1)) * Bt(12,j)
          end do
       end do
    end do

    do j = 1, n
       do i = 1, nn
          jj = i + nn * (j - 1)
          v(jj) = work2(i + nn * (1 - 1)) * Ct(1, j) &
               + work2(i + nn * (2 - 1)) * Ct(2, j) &
               + work2(i + nn * (3 - 1)) * Ct(3, j) &
               + work2(i + nn * (4 - 1)) * Ct(4, j) &
               + work2(i + nn * (5 - 1)) * Ct(5, j) &
               + work2(i + nn * (6 - 1)) * Ct(6, j) &
               + work2(i + nn * (7 - 1)) * Ct(7, j) &
               + work2(i + nn * (8 - 1)) * Ct(8, j) &
               + work2(i + nn * (9 - 1)) * Ct(9, j) &
               + work2(i + nn * (10 - 1)) * Ct(10, j) &
               + work2(i + nn * (11 - 1)) * Ct(11, j) &
               + work2(i + nn * (12 - 1)) * Ct(12, j)
       end do
    end do

  end subroutine tnsr3d_el_n12_cpu

  !> Single-element tensor product specialised for @a nv = @a nu = 11.
  !! @details Unrolled form of tnsr3d_el_n_cpu; see tnsr3d_el_cpu for the
  !! operation and the argument meanings.
  !OCL SERIAL
  subroutine tnsr3d_el_n11_cpu(v, u, A, Bt, Ct)
    integer, parameter :: n = 11
    integer, parameter :: nn = n**2
    real(kind=rp), intent(inout) :: v(n*n*n), u(n*n*n)
    real(kind=rp), intent(inout) :: A(n,n), Bt(n,n), Ct(n,n)
    real(kind=rp) :: work(n**3), work2(n**3)
    integer :: i, j, l
    integer :: ii, jj

    do j = 1, nn
       do i = 1, n
          ii = i + n * (j - 1)
          work(ii) = A(i,1) * u(1 + n * (j - 1)) &
               + A(i,2) * u(2 + n * (j - 1)) &
               + A(i,3) * u(3 + n * (j - 1)) &
               + A(i,4) * u(4 + n * (j - 1)) &
               + A(i,5) * u(5 + n * (j - 1)) &
               + A(i,6) * u(6 + n * (j - 1)) &
               + A(i,7) * u(7 + n * (j - 1)) &
               + A(i,8) * u(8 + n * (j - 1)) &
               + A(i,9) * u(9 + n * (j - 1)) &
               + A(i,10) * u(10 + n * (j - 1)) &
               + A(i,11) * u(11 + n * (j - 1))
       end do
    end do

    do i = 1, n
       do j = 1, n
          do l = 1, n
             ii = l + n * (j - 1) + nn * (i - 1)
             work2(ii) = work(l + n * (1 - 1) + nn * (i - 1)) * Bt(1,j) &
                  + work(l + n * (2 - 1) + nn * (i - 1)) * Bt(2,j) &
                  + work(l + n * (3 - 1) + nn * (i - 1)) * Bt(3,j) &
                  + work(l + n * (4 - 1) + nn * (i - 1)) * Bt(4,j) &
                  + work(l + n * (5 - 1) + nn * (i - 1)) * Bt(5,j) &
                  + work(l + n * (6 - 1) + nn * (i - 1)) * Bt(6,j) &
                  + work(l + n * (7 - 1) + nn * (i - 1)) * Bt(7,j) &
                  + work(l + n * (8 - 1) + nn * (i - 1)) * Bt(8,j) &
                  + work(l + n * (9 - 1) + nn * (i - 1)) * Bt(9,j) &
                  + work(l + n * (10 - 1) + nn * (i - 1)) * Bt(10,j) &
                  + work(l + n * (11 - 1) + nn * (i - 1)) * Bt(11,j)
          end do
       end do
    end do

    do j = 1, n
       do i = 1, nn
          jj = i + nn * (j - 1)
          v(jj) = work2(i + nn * (1 - 1)) * Ct(1, j) &
               + work2(i + nn * (2 - 1)) * Ct(2, j) &
               + work2(i + nn * (3 - 1)) * Ct(3, j) &
               + work2(i + nn * (4 - 1)) * Ct(4, j) &
               + work2(i + nn * (5 - 1)) * Ct(5, j) &
               + work2(i + nn * (6 - 1)) * Ct(6, j) &
               + work2(i + nn * (7 - 1)) * Ct(7, j) &
               + work2(i + nn * (8 - 1)) * Ct(8, j) &
               + work2(i + nn * (9 - 1)) * Ct(9, j) &
               + work2(i + nn * (10 - 1)) * Ct(10, j) &
               + work2(i + nn * (11 - 1)) * Ct(11, j)
       end do
    end do

  end subroutine tnsr3d_el_n11_cpu

  !> Single-element tensor product specialised for @a nv = @a nu = 10.
  !! @details Unrolled form of tnsr3d_el_n_cpu; see tnsr3d_el_cpu for the
  !! operation and the argument meanings.
  !OCL SERIAL
  subroutine tnsr3d_el_n10_cpu(v, u, A, Bt, Ct)
    integer, parameter :: n = 10
    integer, parameter :: nn = n**2
    real(kind=rp), intent(inout) :: v(n*n*n), u(n*n*n)
    real(kind=rp), intent(inout) :: A(n,n), Bt(n,n), Ct(n,n)
    real(kind=rp) :: work(n**3), work2(n**3)
    integer :: i, j, l
    integer :: ii, jj

    do j = 1, nn
       do i = 1, n
          ii = i + n * (j - 1)
          work(ii) = A(i,1) * u(1 + n * (j - 1)) &
               + A(i,2) * u(2 + n * (j - 1)) &
               + A(i,3) * u(3 + n * (j - 1)) &
               + A(i,4) * u(4 + n * (j - 1)) &
               + A(i,5) * u(5 + n * (j - 1)) &
               + A(i,6) * u(6 + n * (j - 1)) &
               + A(i,7) * u(7 + n * (j - 1)) &
               + A(i,8) * u(8 + n * (j - 1)) &
               + A(i,9) * u(9 + n * (j - 1)) &
               + A(i,10) * u(10 + n * (j - 1))
       end do
    end do

    do i = 1, n
       do j = 1, n
          do l = 1, n
             ii = l + n * (j - 1) + nn * (i - 1)
             work2(ii) = work(l + n * (1 - 1) + nn * (i - 1)) * Bt(1,j) &
                  + work(l + n * (2 - 1) + nn * (i - 1)) * Bt(2,j) &
                  + work(l + n * (3 - 1) + nn * (i - 1)) * Bt(3,j) &
                  + work(l + n * (4 - 1) + nn * (i - 1)) * Bt(4,j) &
                  + work(l + n * (5 - 1) + nn * (i - 1)) * Bt(5,j) &
                  + work(l + n * (6 - 1) + nn * (i - 1)) * Bt(6,j) &
                  + work(l + n * (7 - 1) + nn * (i - 1)) * Bt(7,j) &
                  + work(l + n * (8 - 1) + nn * (i - 1)) * Bt(8,j) &
                  + work(l + n * (9 - 1) + nn * (i - 1)) * Bt(9,j) &
                  + work(l + n * (10 - 1) + nn * (i - 1)) * Bt(10,j)
          end do
       end do
    end do

    do j = 1, n
       do i = 1, nn
          jj = i + nn * (j - 1)
          v(jj) = work2(i + nn * (1 - 1)) * Ct(1, j) &
               + work2(i + nn * (2 - 1)) * Ct(2, j) &
               + work2(i + nn * (3 - 1)) * Ct(3, j) &
               + work2(i + nn * (4 - 1)) * Ct(4, j) &
               + work2(i + nn * (5 - 1)) * Ct(5, j) &
               + work2(i + nn * (6 - 1)) * Ct(6, j) &
               + work2(i + nn * (7 - 1)) * Ct(7, j) &
               + work2(i + nn * (8 - 1)) * Ct(8, j) &
               + work2(i + nn * (9 - 1)) * Ct(9, j) &
               + work2(i + nn * (10 - 1)) * Ct(10, j)
       end do
    end do

  end subroutine tnsr3d_el_n10_cpu

  !> Single-element tensor product specialised for @a nv = @a nu = 9.
  !! @details Unrolled form of tnsr3d_el_n_cpu; see tnsr3d_el_cpu for the
  !! operation and the argument meanings.
  !OCL SERIAL
  subroutine tnsr3d_el_n9_cpu(v, u, A, Bt, Ct)
    integer, parameter :: n = 9
    integer, parameter :: nn = n**2
    real(kind=rp), intent(inout) :: v(n*n*n), u(n*n*n)
    real(kind=rp), intent(inout) :: A(n,n), Bt(n,n), Ct(n,n)
    real(kind=rp) :: work(n**3), work2(n**3)
    integer :: i, j, l
    integer :: ii, jj

    do j = 1, nn
       do i = 1, n
          ii = i + n * (j - 1)
          work(ii) = A(i,1) * u(1 + n * (j - 1)) &
               + A(i,2) * u(2 + n * (j - 1)) &
               + A(i,3) * u(3 + n * (j - 1)) &
               + A(i,4) * u(4 + n * (j - 1)) &
               + A(i,5) * u(5 + n * (j - 1)) &
               + A(i,6) * u(6 + n * (j - 1)) &
               + A(i,7) * u(7 + n * (j - 1)) &
               + A(i,8) * u(8 + n * (j - 1)) &
               + A(i,9) * u(9 + n * (j - 1))
       end do
    end do

    do i = 1, n
       do j = 1, n
          do l = 1, n
             ii = l + n * (j - 1) + nn * (i - 1)
             work2(ii) = work(l + n * (1 - 1) + nn * (i - 1)) * Bt(1,j) &
                  + work(l + n * (2 - 1) + nn * (i - 1)) * Bt(2,j) &
                  + work(l + n * (3 - 1) + nn * (i - 1)) * Bt(3,j) &
                  + work(l + n * (4 - 1) + nn * (i - 1)) * Bt(4,j) &
                  + work(l + n * (5 - 1) + nn * (i - 1)) * Bt(5,j) &
                  + work(l + n * (6 - 1) + nn * (i - 1)) * Bt(6,j) &
                  + work(l + n * (7 - 1) + nn * (i - 1)) * Bt(7,j) &
                  + work(l + n * (8 - 1) + nn * (i - 1)) * Bt(8,j) &
                  + work(l + n * (9 - 1) + nn * (i - 1)) * Bt(9,j)
          end do
       end do
    end do

    do j = 1, n
       do i = 1, nn
          jj = i + nn * (j - 1)
          v(jj) = work2(i + nn * (1 - 1)) * Ct(1, j) &
               + work2(i + nn * (2 - 1)) * Ct(2, j) &
               + work2(i + nn * (3 - 1)) * Ct(3, j) &
               + work2(i + nn * (4 - 1)) * Ct(4, j) &
               + work2(i + nn * (5 - 1)) * Ct(5, j) &
               + work2(i + nn * (6 - 1)) * Ct(6, j) &
               + work2(i + nn * (7 - 1)) * Ct(7, j) &
               + work2(i + nn * (8 - 1)) * Ct(8, j) &
               + work2(i + nn * (9 - 1)) * Ct(9, j)
       end do
    end do

  end subroutine tnsr3d_el_n9_cpu

  !> Single-element tensor product specialised for @a nv = @a nu = 8.
  !! @details Unrolled form of tnsr3d_el_n_cpu; see tnsr3d_el_cpu for the
  !! operation and the argument meanings.
  !OCL SERIAL
  subroutine tnsr3d_el_n8_cpu(v, u, A, Bt, Ct)
    integer, parameter :: n = 8
    integer, parameter :: nn = n**2
    real(kind=rp), intent(inout) :: v(n*n*n), u(n*n*n)
    real(kind=rp), intent(inout) :: A(n,n), Bt(n,n), Ct(n,n)
    real(kind=rp) :: work(n**3), work2(n**3)
    integer :: i, j, l
    integer :: ii, jj

    do j = 1, nn
       do i = 1, n
          ii = i + n * (j - 1)
          work(ii) = A(i,1) * u(1 + n * (j - 1)) &
               + A(i,2) * u(2 + n * (j - 1)) &
               + A(i,3) * u(3 + n * (j - 1)) &
               + A(i,4) * u(4 + n * (j - 1)) &
               + A(i,5) * u(5 + n * (j - 1)) &
               + A(i,6) * u(6 + n * (j - 1)) &
               + A(i,7) * u(7 + n * (j - 1)) &
               + A(i,8) * u(8 + n * (j - 1))
       end do
    end do

    do i = 1, n
       do j = 1, n
          do l = 1, n
             ii = l + n * (j - 1) + nn * (i - 1)
             work2(ii) = work(l + n * (1 - 1) + nn * (i - 1)) * Bt(1,j) &
                  + work(l + n * (2 - 1) + nn * (i - 1)) * Bt(2,j) &
                  + work(l + n * (3 - 1) + nn * (i - 1)) * Bt(3,j) &
                  + work(l + n * (4 - 1) + nn * (i - 1)) * Bt(4,j) &
                  + work(l + n * (5 - 1) + nn * (i - 1)) * Bt(5,j) &
                  + work(l + n * (6 - 1) + nn * (i - 1)) * Bt(6,j) &
                  + work(l + n * (7 - 1) + nn * (i - 1)) * Bt(7,j) &
                  + work(l + n * (8 - 1) + nn * (i - 1)) * Bt(8,j)
          end do
       end do
    end do

    do j = 1, n
       do i = 1, nn
          jj = i + nn * (j - 1)
          v(jj) = work2(i + nn * (1 - 1)) * Ct(1, j) &
               + work2(i + nn * (2 - 1)) * Ct(2, j) &
               + work2(i + nn * (3 - 1)) * Ct(3, j) &
               + work2(i + nn * (4 - 1)) * Ct(4, j) &
               + work2(i + nn * (5 - 1)) * Ct(5, j) &
               + work2(i + nn * (6 - 1)) * Ct(6, j) &
               + work2(i + nn * (7 - 1)) * Ct(7, j) &
               + work2(i + nn * (8 - 1)) * Ct(8, j)
       end do
    end do

  end subroutine tnsr3d_el_n8_cpu

  !> Single-element tensor product specialised for @a nv = @a nu = 7.
  !! @details Unrolled form of tnsr3d_el_n_cpu; see tnsr3d_el_cpu for the
  !! operation and the argument meanings.
  !OCL SERIAL
  subroutine tnsr3d_el_n7_cpu(v, u, A, Bt, Ct)
    integer, parameter :: n = 7
    integer, parameter :: nn = n**2
    real(kind=rp), intent(inout) :: v(n*n*n), u(n*n*n)
    real(kind=rp), intent(inout) :: A(n,n), Bt(n,n), Ct(n,n)
    real(kind=rp) :: work(n**3), work2(n**3)
    integer :: i, j, l
    integer :: ii, jj

    do j = 1, nn
       do i = 1, n
          ii = i + n * (j - 1)
          work(ii) = A(i,1) * u(1 + n * (j - 1)) &
               + A(i,2) * u(2 + n * (j - 1)) &
               + A(i,3) * u(3 + n * (j - 1)) &
               + A(i,4) * u(4 + n * (j - 1)) &
               + A(i,5) * u(5 + n * (j - 1)) &
               + A(i,6) * u(6 + n * (j - 1)) &
               + A(i,7) * u(7 + n * (j - 1))
       end do
    end do

    do i = 1, n
       do j = 1, n
          do l = 1, n
             ii = l + n * (j - 1) + nn * (i - 1)
             work2(ii) = work(l + n * (1 - 1) + nn * (i - 1)) * Bt(1,j) &
                  + work(l + n * (2 - 1) + nn * (i - 1)) * Bt(2,j) &
                  + work(l + n * (3 - 1) + nn * (i - 1)) * Bt(3,j) &
                  + work(l + n * (4 - 1) + nn * (i - 1)) * Bt(4,j) &
                  + work(l + n * (5 - 1) + nn * (i - 1)) * Bt(5,j) &
                  + work(l + n * (6 - 1) + nn * (i - 1)) * Bt(6,j) &
                  + work(l + n * (7 - 1) + nn * (i - 1)) * Bt(7,j)
          end do
       end do
    end do

    do j = 1, n
       do i = 1, nn
          jj = i + nn * (j - 1)
          v(jj) = work2(i + nn * (1 - 1)) * Ct(1, j) &
               + work2(i + nn * (2 - 1)) * Ct(2, j) &
               + work2(i + nn * (3 - 1)) * Ct(3, j) &
               + work2(i + nn * (4 - 1)) * Ct(4, j) &
               + work2(i + nn * (5 - 1)) * Ct(5, j) &
               + work2(i + nn * (6 - 1)) * Ct(6, j) &
               + work2(i + nn * (7 - 1)) * Ct(7, j)
       end do
    end do

  end subroutine tnsr3d_el_n7_cpu

  !> Single-element tensor product specialised for @a nv = @a nu = 6.
  !! @details Unrolled form of tnsr3d_el_n_cpu; see tnsr3d_el_cpu for the
  !! operation and the argument meanings.
  !OCL SERIAL
  subroutine tnsr3d_el_n6_cpu(v, u, A, Bt, Ct)
    integer, parameter :: n = 6
    integer, parameter :: nn = n**2
    real(kind=rp), intent(inout) :: v(n*n*n), u(n*n*n)
    real(kind=rp), intent(inout) :: A(n,n), Bt(n,n), Ct(n,n)
    real(kind=rp) :: work(n**3), work2(n**3)
    integer :: i, j, l
    integer :: ii, jj

    do j = 1, nn
       do i = 1, n
          ii = i + n * (j - 1)
          work(ii) = A(i,1) * u(1 + n * (j - 1)) &
               + A(i,2) * u(2 + n * (j - 1)) &
               + A(i,3) * u(3 + n * (j - 1)) &
               + A(i,4) * u(4 + n * (j - 1)) &
               + A(i,5) * u(5 + n * (j - 1)) &
               + A(i,6) * u(6 + n * (j - 1))
       end do
    end do

    do i = 1, n
       do j = 1, n
          do l = 1, n
             ii = l + n * (j - 1) + nn * (i - 1)
             work2(ii) = work(l + n * (1 - 1) + nn * (i - 1)) * Bt(1,j) &
                  + work(l + n * (2 - 1) + nn * (i - 1)) * Bt(2,j) &
                  + work(l + n * (3 - 1) + nn * (i - 1)) * Bt(3,j) &
                  + work(l + n * (4 - 1) + nn * (i - 1)) * Bt(4,j) &
                  + work(l + n * (5 - 1) + nn * (i - 1)) * Bt(5,j) &
                  + work(l + n * (6 - 1) + nn * (i - 1)) * Bt(6,j)
          end do
       end do
    end do

    do j = 1, n
       do i = 1, nn
          jj = i + nn * (j - 1)
          v(jj) = work2(i + nn * (1 - 1)) * Ct(1, j) &
               + work2(i + nn * (2 - 1)) * Ct(2, j) &
               + work2(i + nn * (3 - 1)) * Ct(3, j) &
               + work2(i + nn * (4 - 1)) * Ct(4, j) &
               + work2(i + nn * (5 - 1)) * Ct(5, j) &
               + work2(i + nn * (6 - 1)) * Ct(6, j)
       end do
    end do

  end subroutine tnsr3d_el_n6_cpu

  !> Single-element tensor product specialised for @a nv = @a nu = 5.
  !! @details Unrolled form of tnsr3d_el_n_cpu; see tnsr3d_el_cpu for the
  !! operation and the argument meanings.
  !OCL SERIAL
  subroutine tnsr3d_el_n5_cpu(v, u, A, Bt, Ct)
    integer, parameter :: n = 5
    integer, parameter :: nn = n**2
    real(kind=rp), intent(inout) :: v(n*n*n), u(n*n*n)
    real(kind=rp), intent(inout) :: A(n,n), Bt(n,n), Ct(n,n)
    real(kind=rp) :: work(n**3), work2(n**3)
    integer :: i, j, l
    integer :: ii, jj

    do j = 1, nn
       do i = 1, n
          ii = i + n * (j - 1)
          work(ii) = A(i,1) * u(1 + n * (j - 1)) &
               + A(i,2) * u(2 + n * (j - 1)) &
               + A(i,3) * u(3 + n * (j - 1)) &
               + A(i,4) * u(4 + n * (j - 1)) &
               + A(i,5) * u(5 + n * (j - 1))
       end do
    end do

    do i = 1, n
       do j = 1, n
          do l = 1, n
             ii = l + n * (j - 1) + nn * (i - 1)
             work2(ii) = work(l + n * (1 - 1) + nn * (i - 1)) * Bt(1,j) &
                  + work(l + n * (2 - 1) + nn * (i - 1)) * Bt(2,j) &
                  + work(l + n * (3 - 1) + nn * (i - 1)) * Bt(3,j) &
                  + work(l + n * (4 - 1) + nn * (i - 1)) * Bt(4,j) &
                  + work(l + n * (5 - 1) + nn * (i - 1)) * Bt(5,j)
          end do
       end do
    end do

    do j = 1, n
       do i = 1, nn
          jj = i + nn * (j - 1)
          v(jj) = work2(i + nn * (1 - 1)) * Ct(1, j) &
               + work2(i + nn * (2 - 1)) * Ct(2, j) &
               + work2(i + nn * (3 - 1)) * Ct(3, j) &
               + work2(i + nn * (4 - 1)) * Ct(4, j) &
               + work2(i + nn * (5 - 1)) * Ct(5, j)
       end do
    end do

  end subroutine tnsr3d_el_n5_cpu

  !> Single-element tensor product specialised for @a nv = @a nu = 4.
  !! @details Unrolled form of tnsr3d_el_n_cpu; see tnsr3d_el_cpu for the
  !! operation and the argument meanings.
  !OCL SERIAL
  subroutine tnsr3d_el_n4_cpu(v, u, A, Bt, Ct)
    integer, parameter :: n = 4
    integer, parameter :: nn = n**2
    real(kind=rp), intent(inout) :: v(n*n*n), u(n*n*n)
    real(kind=rp), intent(inout) :: A(n,n), Bt(n,n), Ct(n,n)
    real(kind=rp) :: work(n**3), work2(n**3)
    integer :: i, j, l
    integer :: ii, jj

    do j = 1, nn
       do i = 1, n
          ii = i + n * (j - 1)
          work(ii) = A(i,1) * u(1 + n * (j - 1)) &
               + A(i,2) * u(2 + n * (j - 1)) &
               + A(i,3) * u(3 + n * (j - 1)) &
               + A(i,4) * u(4 + n * (j - 1))
       end do
    end do

    do i = 1, n
       do j = 1, n
          do l = 1, n
             ii = l + n * (j - 1) + nn * (i - 1)
             work2(ii) = work(l + n * (1 - 1) + nn * (i - 1)) * Bt(1,j) &
                  + work(l + n * (2 - 1) + nn * (i - 1)) * Bt(2,j) &
                  + work(l + n * (3 - 1) + nn * (i - 1)) * Bt(3,j) &
                  + work(l + n * (4 - 1) + nn * (i - 1)) * Bt(4,j)
          end do
       end do
    end do

    do j = 1, n
       do i = 1, nn
          jj = i + nn * (j - 1)
          v(jj) = work2(i + nn * (1 - 1)) * Ct(1, j) &
               + work2(i + nn * (2 - 1)) * Ct(2, j) &
               + work2(i + nn * (3 - 1)) * Ct(3, j) &
               + work2(i + nn * (4 - 1)) * Ct(4, j)
       end do
    end do

  end subroutine tnsr3d_el_n4_cpu

  !> Single-element tensor product specialised for @a nv = @a nu = 3.
  !! @details Unrolled form of tnsr3d_el_n_cpu; see tnsr3d_el_cpu for the
  !! operation and the argument meanings.
  !OCL SERIAL
  subroutine tnsr3d_el_n3_cpu(v, u, A, Bt, Ct)
    integer, parameter :: n = 3
    integer, parameter :: nn = n**2
    real(kind=rp), intent(inout) :: v(n*n*n), u(n*n*n)
    real(kind=rp), intent(inout) :: A(n,n), Bt(n,n), Ct(n,n)
    real(kind=rp) :: work(n**3), work2(n**3)
    integer :: i, j, l
    integer :: ii, jj

    do j = 1, nn
       do i = 1, n
          ii = i + n * (j - 1)
          work(ii) = A(i,1) * u(1 + n * (j - 1)) &
               + A(i,2) * u(2 + n * (j - 1)) &
               + A(i,3) * u(3 + n * (j - 1))
       end do
    end do

    do i = 1, n
       do j = 1, n
          do l = 1, n
             ii = l + n * (j - 1) + nn * (i - 1)
             work2(ii) = work(l + n * (1 - 1) + nn * (i - 1)) * Bt(1,j) &
                  + work(l + n * (2 - 1) + nn * (i - 1)) * Bt(2,j) &
                  + work(l + n * (3 - 1) + nn * (i - 1)) * Bt(3,j)
          end do
       end do
    end do

    do j = 1, n
       do i = 1, nn
          jj = i + nn * (j - 1)
          v(jj) = work2(i + nn * (1 - 1)) * Ct(1, j) &
               + work2(i + nn * (2 - 1)) * Ct(2, j) &
               + work2(i + nn * (3 - 1)) * Ct(3, j)
       end do
    end do

  end subroutine tnsr3d_el_n3_cpu

  !> Single-element tensor product specialised for @a nv = @a nu = 2.
  !! @details Unrolled form of tnsr3d_el_n_cpu; see tnsr3d_el_cpu for the
  !! operation and the argument meanings.
  !OCL SERIAL
  subroutine tnsr3d_el_n2_cpu(v, u, A, Bt, Ct)
    integer, parameter :: n = 2
    integer, parameter :: nn = n**2
    real(kind=rp), intent(inout) :: v(n*n*n), u(n*n*n)
    real(kind=rp), intent(inout) :: A(n,n), Bt(n,n), Ct(n,n)
    real(kind=rp) :: work(n**3), work2(n**3)
    integer :: i, j, l
    integer :: ii, jj

    do j = 1, nn
       do i = 1, n
          ii = i + n * (j - 1)
          work(ii) = A(i,1) * u(1 + n * (j - 1)) &
               + A(i,2) * u(2 + n * (j - 1))
       end do
    end do

    do i = 1, n
       do j = 1, n
          do l = 1, n
             ii = l + n * (j - 1) + nn * (i - 1)
             work2(ii) = work(l + n * (1 - 1) + nn * (i - 1)) * Bt(1,j) &
                  + work(l + n * (2 - 1) + nn * (i - 1)) * Bt(2,j)
          end do
       end do
    end do

    do j = 1, n
       do i = 1, nn
          jj = i + nn * (j - 1)
          v(jj) = work2(i + nn * (1 - 1)) * Ct(1, j) &
               + work2(i + nn * (2 - 1)) * Ct(2, j)
       end do
    end do

  end subroutine tnsr3d_el_n2_cpu

  !> Three-dimensional tensor product over a list of elements.
  !! @details The batched counterpart of tnsr3d_el_cpu: applies
  !! \f$ v = (A \otimes B^T \otimes C^T) u \f$ to every element of @a u,
  !! taking each from order @a nu to order @a nv. This is the entry point used
  !! for dealiasing (@a nu and @a nv related by the 3/2 rule) and for the phmg
  !! and hsmg grid transfers. Dispatches on @a nu to an unrolled kernel where
  !! one exists, otherwise to the generic tnsr3d_nvnu_cpu.
  !! @note Unlike the `_el_` routines this opens its own parallel region, over
  !! the element loop.
  !! @param v Result, nv*nv*nv per element.
  !! @param nv Output polynomial order.
  !! @param u Input, nu*nu*nu per element.
  !! @param nu Input polynomial order.
  !! @param A Operator along the first reference direction, (nv, nu).
  !! @param Bt Transposed operator along the second direction, (nu, nv).
  !! @param Ct Transposed operator along the third direction, (nu, nv).
  !! @param nelv Number of elements.
  subroutine tnsr3d_cpu(v, nv, u, nu, A, Bt, Ct, nelv)
    integer, intent(in) :: nv, nu, nelv
    real(kind=rp), intent(inout) :: v(nv*nv*nv, nelv)
    real(kind=rp), intent(in) :: u(nu*nu*nu, nelv)
    real(kind=rp), intent(in) :: A(nv, nu), Bt(nu, nv), Ct(nu, nv)

    if (nu .eq. 2 .and. nv .eq. 4) then
       call tnsr3d_nu2nv4_cpu(v, u, A, Bt, Ct, nelv)
    else if (nu .eq. 4) then
       call tnsr3d_nu4_cpu(v, nv, u, A, Bt, Ct, nelv)
    else if (nu .eq. 8) then
       call tnsr3d_nu8_cpu(v, nv, u, A, Bt, Ct, nelv)
    else if (nu .eq. 12) then
       call tnsr3d_nu12_cpu(v, nv, u, A, Bt, Ct, nelv)
    else
       call tnsr3d_nvnu_cpu(v, nv, u, nu, A, Bt, Ct, nelv)
    end if

  end subroutine tnsr3d_cpu

  !> Generic tensor-product evaluation for arbitrary @a nu, @a nv.
  !! @details Each of the three contractions sums over @a k, the index that
  !! strides through memory. Written with @a k innermost the reduction is a
  !! dependent scalar chain over strided loads, which neither SVE nor any other
  !! unit-stride SIMD unit can do anything with. The @a k loop is therefore
  !! hoisted out and peeled at @a k = 1, leaving an innermost loop that runs
  !! over the unit-stride index of both operand and result with a
  !! loop-invariant scalar multiplier, i.e. a plain axpy. The summation order
  !! over @a k is unchanged, so results match the reduction form. This is the
  !! same shape the unrolled tnsr3d_nu4_cpu already has, expressed without
  !! knowing @a nu at compile time.
  !! See tnsr3d_cpu for the argument meanings.
  subroutine tnsr3d_nvnu_cpu(v, nv, u, nu, A, Bt, Ct, nelv)
    integer, intent(in) :: nv, nu, nelv
    real(kind=rp), intent(inout) :: v(nv*nv*nv, nelv)
    real(kind=rp), intent(in) :: u(nu*nu*nu, nelv)
    real(kind=rp), intent(in) :: A(nv, nu), Bt(nu, nv), Ct(nu, nv)
    real(kind=rp) :: work(nu**2*nv), work2(nu*nv**2), c
    integer :: ie, i, j, k, l, ii, jj, kk
    integer :: nunu, nvnu, nvnv

    nvnu = nv * nu
    nunu = nu * nu
    nvnv = nv * nv

    !$omp parallel do private(ie, i, j, k, l, ii, jj, kk, c, work, work2)
    do ie = 1, nelv
       do j = 1, nunu
          ii = nv * (j - 1)
          kk = nu * (j - 1)
          c = u(1 + kk, ie)
          do i = 1, nv
             work(i + ii) = A(i,1) * c
          end do
          do k = 2, nu
             c = u(k + kk, ie)
             !OCL NORECURRENCE, NOVREC, NOALIAS
             !DIR$ IVDEP
             !GCC$ ivdep
             do i = 1, nv
                work(i + ii) = work(i + ii) + A(i,k) * c
             end do
          end do
       end do

       do i = 1, nu
          do j = 1, nv
             ii = nv * (j - 1) + nvnv * (i - 1)
             jj = nvnu * (i - 1)
             c = Bt(1,j)
             do l = 1, nv
                work2(l + ii) = work(l + jj) * c
             end do
             do k = 2, nu
                c = Bt(k,j)
                kk = nv * (k - 1) + nvnu * (i - 1)
                !OCL NORECURRENCE, NOVREC, NOALIAS
                !DIR$ IVDEP
                !GCC$ ivdep
                do l = 1, nv
                   work2(l + ii) = work2(l + ii) + work(l + kk) * c
                end do
             end do
          end do
       end do

       do j = 1, nv
          jj = nvnv * (j - 1)
          c = Ct(1,j)
          do i = 1, nvnv
             v(i + jj, ie) = work2(i) * c
          end do
          do k = 2, nu
             c = Ct(k,j)
             ii = nvnv * (k - 1)
             !OCL NORECURRENCE, NOVREC, NOALIAS
             !DIR$ IVDEP
             !GCC$ ivdep
             do i = 1, nvnv
                v(i + jj, ie) = v(i + jj, ie) + work2(i + ii) * c
             end do
          end do
       end do
    end do
    !$omp end parallel do

  end subroutine tnsr3d_nvnu_cpu

  !> Batched tensor product specialised for @a nu = 2, @a nv = 4.
  !! @details Prolongation from the coarsest multigrid level. See tnsr3d_cpu
  !! for the operation and the argument meanings.
  subroutine tnsr3d_nu2nv4_cpu(v, u, A, Bt, Ct, nelv)
    integer, parameter :: nu = 2
    integer, parameter :: nv = 4
    integer, parameter :: nunu = 4
    integer, parameter :: nvnu = 8
    integer, parameter :: nvnv = 16
    integer, intent(in) :: nelv
    real(kind=rp), intent(inout) :: v(nv*nv*nv, nelv)
    real(kind=rp), intent(in) :: u(nu*nu*nu, nelv)
    real(kind=rp), intent(in) :: A(nv, nu), Bt(nu, nv), Ct(nu, nv)
    real(kind=rp) :: work(nu**2*nv), work2(nu*nv**2), tmp
    integer :: ie, i, j, k, l, ii, jj

    !$omp parallel do private(ie, i, j, k, l, ii, jj, tmp, work, work2)
    do ie = 1, nelv
       do j = 1, nunu
          do i = 1, nv
             ii = i + nv * (j - 1)
             work(ii) = A(i,1) * u(1 + nu * (j - 1), ie) &
                  + A(i,2) * u(2 + nu * (j - 1), ie)
          end do
       end do

       do i = 1, nu
          do j = 1, nv
             do l = 1, nv
                ii = l + nv * (j - 1) + nvnv * (i - 1)
                tmp = 0.0_rp
                do k = 1, nu
                   jj = l + nv * (k - 1) + nvnu * (i - 1)
                   tmp = tmp + work(jj) * Bt(k,j)
                end do
                work2(ii) = tmp
             end do
          end do
       end do

       do j = 1, nv
          do i = 1, nvnv
             jj = i + nvnv * (j - 1)
             v(jj, ie) = work2(i + nvnv * (1 - 1)) * Ct(1, j) &
                  + work2(i + nvnv * (2 - 1)) * Ct(2, j)
          end do
       end do
    end do
    !$omp end parallel do

  end subroutine tnsr3d_nu2nv4_cpu

  !> Batched tensor product specialised for @a nu = 4, generic in @a nv.
  !! @details Covers the phmg transfers to and from an lx = 4 level. See
  !! tnsr3d_cpu for the operation and the argument meanings.
  subroutine tnsr3d_nu4_cpu(v, nv, u, A, Bt, Ct, nelv)
    integer, parameter :: nu = 4
    integer, parameter :: nunu = 16
    integer, intent(in) :: nv, nelv
    real(kind=rp), intent(inout) :: v(nv*nv*nv, nelv)
    real(kind=rp), intent(in) :: u(nu*nu*nu, nelv)
    real(kind=rp), intent(in) :: A(nv, nu), Bt(nu, nv), Ct(nu, nv)
    real(kind=rp) :: work(nu**2*nv), work2(nu*nv**2), tmp
    integer :: ie, i, j, k, l, ii, jj
    integer :: nvnu, nvnv

    nvnu = nv * nu
    nvnv = nv * nv

    !$omp parallel do private(ie, i, j, k, l, ii, jj, tmp, work, work2)
    do ie = 1, nelv
       do j = 1, nunu
          do i = 1, nv
             ii = i + nv * (j - 1)
             work(ii) = A(i,1) * u(1 + nu * (j - 1), ie) &
                  + A(i,2) * u(2 + nu * (j - 1), ie) &
                  + A(i,3) * u(3 + nu * (j - 1), ie) &
                  + A(i,4) * u(4 + nu * (j - 1), ie)
          end do
       end do

       do i = 1, nu
          do j = 1, nv
             do l = 1, nv
                ii = l + nv * (j - 1) + nvnv * (i - 1)
                tmp = 0.0_rp
                do k = 1, nu
                   jj = l + nv * (k - 1) + nvnu * (i - 1)
                   tmp = tmp + work(jj) * Bt(k,j)
                end do
                work2(ii) = tmp
             end do
          end do
       end do

       do j = 1, nv
          do i = 1, nvnv
             jj = i + nvnv * (j - 1)
             v(jj, ie) = work2(i + nvnv * (1 - 1)) * Ct(1, j) &
                  + work2(i + nvnv * (2 - 1)) * Ct(2, j) &
                  + work2(i + nvnv * (3 - 1)) * Ct(3, j) &
                  + work2(i + nvnv * (4 - 1)) * Ct(4, j)
          end do
       end do
    end do
    !$omp end parallel do

  end subroutine tnsr3d_nu4_cpu

  !> Tensor-product evaluation specialised for @a nu = 8, generic in @a nv.
  !! @details Covers the 3/2-rule dealiasing pair for lx = 8
  !! (@a nu = 8 to @a nv = 12) and the phmg transfers that start from an
  !! lx = 8 level. All three contractions are unrolled over @a nu, so the
  !! innermost loop runs over the unit-stride index with no reduction chain and
  !! the result is written once rather than accumulated over @a nu passes.
  !! See tnsr3d_cpu for the operation and the argument meanings.
  subroutine tnsr3d_nu8_cpu(v, nv, u, A, Bt, Ct, nelv)
    integer, parameter :: nu = 8
    integer, parameter :: nunu = 64
    integer, intent(in) :: nv, nelv
    real(kind=rp), intent(inout) :: v(nv*nv*nv, nelv)
    real(kind=rp), intent(in) :: u(nu*nu*nu, nelv)
    real(kind=rp), intent(in) :: A(nv, nu), Bt(nu, nv), Ct(nu, nv)
    real(kind=rp) :: work(nu**2*nv), work2(nu*nv**2)
    integer :: ie, i, j, l, ii, jj
    integer :: nvnu, nvnv

    nvnu = nv * nu
    nvnv = nv * nv

    !$omp parallel do private(ie, i, j, l, ii, jj, work, work2)
    do ie = 1, nelv
       do j = 1, nunu
          do i = 1, nv
             ii = i + nv * (j - 1)
             work(ii) = A(i,1) * u(1 + nu * (j - 1), ie) &
                  + A(i,2) * u(2 + nu * (j - 1), ie) &
                  + A(i,3) * u(3 + nu * (j - 1), ie) &
                  + A(i,4) * u(4 + nu * (j - 1), ie) &
                  + A(i,5) * u(5 + nu * (j - 1), ie) &
                  + A(i,6) * u(6 + nu * (j - 1), ie) &
                  + A(i,7) * u(7 + nu * (j - 1), ie) &
                  + A(i,8) * u(8 + nu * (j - 1), ie)
          end do
       end do

       do i = 1, nu
          do j = 1, nv
             ii = nv * (j - 1) + nvnv * (i - 1)
             jj = nvnu * (i - 1)
             do l = 1, nv
                work2(l + ii) = work(l + jj) * Bt(1,j) &
                     + work(l + nv + jj) * Bt(2,j) &
                     + work(l + 2 * nv + jj) * Bt(3,j) &
                     + work(l + 3 * nv + jj) * Bt(4,j) &
                     + work(l + 4 * nv + jj) * Bt(5,j) &
                     + work(l + 5 * nv + jj) * Bt(6,j) &
                     + work(l + 6 * nv + jj) * Bt(7,j) &
                     + work(l + 7 * nv + jj) * Bt(8,j)
             end do
          end do
       end do

       do j = 1, nv
          do i = 1, nvnv
             jj = i + nvnv * (j - 1)
             v(jj, ie) = work2(i + nvnv * (1 - 1)) * Ct(1, j) &
                  + work2(i + nvnv * (2 - 1)) * Ct(2, j) &
                  + work2(i + nvnv * (3 - 1)) * Ct(3, j) &
                  + work2(i + nvnv * (4 - 1)) * Ct(4, j) &
                  + work2(i + nvnv * (5 - 1)) * Ct(5, j) &
                  + work2(i + nvnv * (6 - 1)) * Ct(6, j) &
                  + work2(i + nvnv * (7 - 1)) * Ct(7, j) &
                  + work2(i + nvnv * (8 - 1)) * Ct(8, j)
          end do
       end do
    end do
    !$omp end parallel do

  end subroutine tnsr3d_nu8_cpu

  !> Tensor-product evaluation specialised for @a nu = 12, generic in @a nv.
  !! @details Covers the 3/2-rule dealiasing pair for lx = 12
  !! (@a nu = 12 to @a nv = 18) and the reverse direction of the lx = 8 pair
  !! (@a nu = 12 to @a nv = 8). Unrolled over @a nu on the same grounds as
  !! tnsr3d_nu8_cpu.
  !! See tnsr3d_cpu for the operation and the argument meanings.
  subroutine tnsr3d_nu12_cpu(v, nv, u, A, Bt, Ct, nelv)
    integer, parameter :: nu = 12
    integer, parameter :: nunu = 144
    integer, intent(in) :: nv, nelv
    real(kind=rp), intent(inout) :: v(nv*nv*nv, nelv)
    real(kind=rp), intent(in) :: u(nu*nu*nu, nelv)
    real(kind=rp), intent(in) :: A(nv, nu), Bt(nu, nv), Ct(nu, nv)
    real(kind=rp) :: work(nu**2*nv), work2(nu*nv**2)
    integer :: ie, i, j, l, ii, jj
    integer :: nvnu, nvnv

    nvnu = nv * nu
    nvnv = nv * nv

    !$omp parallel do private(ie, i, j, l, ii, jj, work, work2)
    do ie = 1, nelv
       do j = 1, nunu
          do i = 1, nv
             ii = i + nv * (j - 1)
             work(ii) = A(i,1) * u(1 + nu * (j - 1), ie) &
                  + A(i,2) * u(2 + nu * (j - 1), ie) &
                  + A(i,3) * u(3 + nu * (j - 1), ie) &
                  + A(i,4) * u(4 + nu * (j - 1), ie) &
                  + A(i,5) * u(5 + nu * (j - 1), ie) &
                  + A(i,6) * u(6 + nu * (j - 1), ie) &
                  + A(i,7) * u(7 + nu * (j - 1), ie) &
                  + A(i,8) * u(8 + nu * (j - 1), ie) &
                  + A(i,9) * u(9 + nu * (j - 1), ie) &
                  + A(i,10) * u(10 + nu * (j - 1), ie) &
                  + A(i,11) * u(11 + nu * (j - 1), ie) &
                  + A(i,12) * u(12 + nu * (j - 1), ie)
          end do
       end do

       do i = 1, nu
          do j = 1, nv
             ii = nv * (j - 1) + nvnv * (i - 1)
             jj = nvnu * (i - 1)
             do l = 1, nv
                work2(l + ii) = work(l + jj) * Bt(1,j) &
                     + work(l + nv + jj) * Bt(2,j) &
                     + work(l + 2 * nv + jj) * Bt(3,j) &
                     + work(l + 3 * nv + jj) * Bt(4,j) &
                     + work(l + 4 * nv + jj) * Bt(5,j) &
                     + work(l + 5 * nv + jj) * Bt(6,j) &
                     + work(l + 6 * nv + jj) * Bt(7,j) &
                     + work(l + 7 * nv + jj) * Bt(8,j) &
                     + work(l + 8 * nv + jj) * Bt(9,j) &
                     + work(l + 9 * nv + jj) * Bt(10,j) &
                     + work(l + 10 * nv + jj) * Bt(11,j) &
                     + work(l + 11 * nv + jj) * Bt(12,j)
             end do
          end do
       end do

       do j = 1, nv
          do i = 1, nvnv
             jj = i + nvnv * (j - 1)
             v(jj, ie) = work2(i + nvnv * (1 - 1)) * Ct(1, j) &
                  + work2(i + nvnv * (2 - 1)) * Ct(2, j) &
                  + work2(i + nvnv * (3 - 1)) * Ct(3, j) &
                  + work2(i + nvnv * (4 - 1)) * Ct(4, j) &
                  + work2(i + nvnv * (5 - 1)) * Ct(5, j) &
                  + work2(i + nvnv * (6 - 1)) * Ct(6, j) &
                  + work2(i + nvnv * (7 - 1)) * Ct(7, j) &
                  + work2(i + nvnv * (8 - 1)) * Ct(8, j) &
                  + work2(i + nvnv * (9 - 1)) * Ct(9, j) &
                  + work2(i + nvnv * (10 - 1)) * Ct(10, j) &
                  + work2(i + nvnv * (11 - 1)) * Ct(11, j) &
                  + work2(i + nvnv * (12 - 1)) * Ct(12, j)
          end do
       end do
    end do
    !$omp end parallel do

  end subroutine tnsr3d_nu12_cpu

  !> Three-dimensional tensor product applied in place.
  !! @details As tnsr3d_cpu, but @a v is both input and output: each element is
  !! read at order @a nu and written back at order @a nv. Elements are walked
  !! backwards when @a nv > @a nu, so that a growing element never overwrites
  !! input it has not read yet.
  !! @param v Input and result, nu*nu*nu per element on entry, nv*nv*nv on exit.
  !! @param nv Output polynomial order.
  !! @param nu Input polynomial order.
  !! @param A Operator along the first reference direction, (nv, nu).
  !! @param Bt Transposed operator along the second direction, (nu, nv).
  !! @param Ct Transposed operator along the third direction, (nu, nv).
  !! @param nelv Number of elements.
  subroutine tnsr1_3d_cpu(v, nv, nu, A, Bt, Ct, nelv)
    integer, intent(in) :: nv, nu, nelv
    real(kind=rp), intent(inout) :: v(nv*nv*nv*nelv)
    real(kind=rp), intent(inout) :: A(nv, nu), Bt(nu, nv), Ct(nu, nv)

    if (nu .eq. 4 .and. nv .eq. 2) then
       call tnsr1_3d_nu4nv2_cpu(v, A, Bt, Ct, nelv)
    else
       call tnsr1_3d_nvnu_cpu(v, nv, nu, A, Bt, Ct, nelv)
    end if

  end subroutine tnsr1_3d_cpu

  !> In-place tensor product for arbitrary @a nu and @a nv.
  !! @details The fallback of tnsr1_3d_cpu, including the reversed element
  !! traversal when the order grows. See tnsr1_3d_cpu for the argument
  !! meanings.
  subroutine tnsr1_3d_nvnu_cpu(v, nv, nu, A, Bt, Ct, nelv)
    integer, intent(in) :: nv, nu, nelv
    real(kind=rp), intent(inout) :: v(nv*nv*nv*nelv)
    real(kind=rp), intent(inout) :: A(nv, nu), Bt(nu, nv), Ct(nu, nv)
    real(kind=rp) :: work(nu**2*nv), work2(nu*nv**2)
    integer :: e, e0, ee, es, iu, iv, nu3, nv3
    integer :: i, j, k, l, ii, jj, kk
    integer :: nunu, nvnu, nvnv
    real(kind=rp) :: tmp

    nvnu = nv * nu
    nunu = nu * nu
    nvnv = nv * nv

    e0 = 1
    es = 1
    ee = nelv

    if (nv .gt. nu) then
       e0 = nelv
       es = -1
       ee = 1
    end if

    nu3 = nu**3
    nv3 = nv**3

    !$omp parallel do private(e,iu,iv,i,j,k,l,ii,jj,kk,work,work2,tmp)
    do e = e0, ee, es
       iu = (e-1)*nu3
       iv = (e-1)*nv3

       do j = 1, nunu
          do i = 1, nv
             ii = i + nv * (j - 1)
             tmp = 0.0_rp
             do k = 1, nu
                kk = k + nu * (j - 1) + iu
                tmp = tmp + A(i,k) * v(kk)
             end do
             work(ii) = tmp
          end do
       end do

       do i = 1, nu
          do j = 1, nv
             do l = 1, nv
                ii = l + nv * (j - 1) + nvnv * (i - 1)
                tmp = 0.0_rp
                do k = 1, nu
                   jj = l + nv * (k - 1) + nvnu * (i - 1)
                   tmp = tmp + work(jj) * Bt(k,j)
                end do
                work2(ii) = tmp
             end do
          end do
       end do

       do j = 1, nv
          do i = 1, nvnv
             jj = i + nvnv * (j - 1) + iv
             tmp = 0.0_rp
             do k = 1, nu
                ii = i + nvnv * (k - 1)
                tmp = tmp + work2(ii) * Ct(k, j)
             end do
             v(jj) = tmp
          end do
       end do
    end do
    !$omp end parallel do
  end subroutine tnsr1_3d_nvnu_cpu

  !> In-place tensor product specialised for @a nu = 4, @a nv = 2.
  !! @details See tnsr1_3d_cpu for the operation and the argument meanings.
  subroutine tnsr1_3d_nu4nv2_cpu(v, A, Bt, Ct, nelv)
    integer, parameter :: nu = 4
    integer, parameter :: nv = 2
    integer, parameter :: nunu = 16
    integer, parameter :: nvnu = 8
    integer, parameter :: nvnv = 4
    integer, parameter :: nununu = 64
    integer, parameter :: nvnvnv = 8
    integer, intent(in) :: nelv
    real(kind=rp), intent(inout) :: v(nv*nv*nv*nelv)
    real(kind=rp), intent(inout) :: A(nv, nu), Bt(nu, nv), Ct(nu, nv)
    real(kind=rp) :: work(nu**2*nv), work2(nu*nv**2)
    integer :: e, iu, iv
    integer :: i, j, k, l, ii, jj
    real(kind=rp) :: tmp

    !$omp parallel do private(e,iu,iv,i,j,k,l,ii,jj, work, work2, tmp)
    do e = 1, nelv
       iu = (e-1)*nununu
       iv = (e-1)*nvnvnv

       do j = 1, nunu
          do i = 1, nv
             ii = i + nv * (j - 1)
             work(ii) = A(i,1) * v(1 + nu * (j - 1) + iu) &
                  + A(i,2) * v(2 + nu * (j - 1) + iu) &
                  + A(i,3) * v(3 + nu * (j - 1) + iu) &
                  + A(i,4) * v(4 + nu * (j - 1) + iu)
          end do
       end do

       do i = 1, nu
          do j = 1, nv
             do l = 1, nv
                ii = l + nv * (j - 1) + nvnv * (i - 1)
                tmp = 0.0_rp
                do k = 1, nu
                   jj = l + nv * (k - 1) + nvnu * (i - 1)
                   tmp = tmp + work(jj) * Bt(k,j)
                end do
                work2(ii) = tmp
             end do
          end do
       end do

       do j = 1, nv
          do i = 1, nvnv
             jj = i + nvnv * (j - 1) + iv
             v(jj) = work2(i + nvnv * (1 - 1)) * Ct(1, j) &
                  + work2(i + nvnv * (2 - 1)) * Ct(2, j) &
                  + work2(i + nvnv * (3 - 1)) * Ct(3, j) &
                  + work2(i + nvnv * (4 - 1)) * Ct(4, j)

          end do
       end do
    end do
    !$omp end parallel do
  end subroutine tnsr1_3d_nu4nv2_cpu

end module tensor_cpu
