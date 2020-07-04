subroutine cg(x, f, g, c, r, w, p, z, n, msk, niter)
  use num_types
  use field
  implicit none

  type(field_t), intent(inout) :: x
  type(field_t), intent(inout) :: msk
  real(kind=dp), intent(inout) :: g(6, x%Xh%lx, x%Xh%ly, x%Xh%lz, x%msh%nelv)
  real(kind=dp), intent(inout), dimension(n) :: f
  real(kind=dp), intent(inout), dimension(n) :: c
  real(kind=dp), intent(inout), dimension(n) :: r
  real(kind=dp), intent(inout), dimension(n) :: w
  real(kind=dp), intent(inout), dimension(n) :: p
  real(kind=dp), intent(inout), dimension(n) :: z
  integer, intent(inout) :: n
  integer, intent(inout) :: niter

  real(kind=dp) :: rnorm, rtr, rtz2, rtz1, beta, pap, alpha, alphm
  integer :: iter

  call rzero(x%x, n)
  call copy(r, f, n)
  call col2(r, msk%x, n)

  pap = 0d0
  rtz1 = 1d0
  
  rnorm = sqrt(glsc3(r, c, r, n))

  do iter = 1, niter
     !sovleM

     rtz2 = rtz1
     rtz1 = glsc3(r, c, z, n)

     beta = rtz1 / rtz2
     if (iter .eq. 1) beta = 0d0
     call add2s1(p, z, beta, n)
     
     !call ax
     pap = glsc3(w, c, p, n)

     alpha = rtz1/pap
     alphm = -alpha
     call add2s2(x%x, p, alpha, n)
     call add2s2(r, w, alphm, n)

     rtr = glsc3(r, c, r, n)
     !!!
     rnorm = sqrt(rtr)
  end do
  
end subroutine cg

subroutine solveM(z, r, n)
  use num_types
  use math
  implicit none
  real(kind=dp), intent(inout), dimension(n) :: z
  real(kind=dp), intent(inout), dimension(n) :: r
  integer, intent(inout) :: n

  call copy(z, r, n)
  
end subroutine solveM


