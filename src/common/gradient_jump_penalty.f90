! Copyright (c) 2024, The Neko Authors
! All rights reserved.
!
! Redistribution and use in source and binary forms, with or without
! modification, are permitted provided that the following conditions
! are met:
!
!   * Redistributions of source code must retain the above copyright
!     notice, this list of conditions and the following disclaimer.
!
!   * Redistributions in binary form must reproduce the above
!     copyright notice, this list of conditions and the following
!     disclaimer in the documentation and/or other materials provided
!     with the distribution.
!
!   * Neither the name of the authors nor the names of its
!     contributors may be used to endorse or promote products derived
!     from this software without specific prior written permission.
!
! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
! "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
! LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
! FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
! COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
! INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
! BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
! LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
! CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
! LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
! ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
! POSSIBILITY OF SUCH DAMAGE.
!
!
!> Implements `gradient_jump_penalty_t`.
module gradient_jump_penalty
  use num_types, only : rp
  use speclib, only : pnleg, pndleg
  use utils, only : neko_error
  use math
  use point, only : point_t
  use field, only : field_t
  use dofmap , only : dofmap_t
  use neko_config, only : NEKO_BCKND_DEVICE
  use coefs, only : coef_t
  use scratch_registry, only : neko_scratch_registry
  use element, only : element_t
  use hex
  use quad
  use operators, only : dudxyz
  use gs_ops, only : GS_OP_ADD
  use tuple, only : tuple_i4_t

  implicit none
  private

  !> Implements the gradient jump penalty.
  !! @note Reference DOI: 10.1016/j.cma.2021.114200
  type, public :: gradient_jump_penalty_t
     !> Coefficient of the penalty.
     real(kind=rp) :: tau
     !> Polynomial order
     integer :: p
     !> Collocation point per dimension (p + 1)
     integer :: lx
     !> Penalty terms
     type(field_t), pointer :: penalty
     !> Work array to store penalty terms on each element
     real(kind=rp), allocatable :: penalty_el(:, :, :)
     !> SEM coefficients.
     type(coef_t), pointer :: coef => null()
     !> Gradient of the quatity of interest
     type(field_t), pointer :: grad_1, grad_2, grad_3 
     !> Gradient jump on the elementary interface (zero inside each element)
     real(kind=rp), allocatable, dimension(:, :, :, :) :: G
     !> 3 parts of the flux of the quantity
     real(kind=rp), allocatable, dimension(:, :, :, :) :: flux1, flux2, flux3
     !> 3 parts of the flux of the volumetric flow
     real(kind=rp), allocatable, dimension(:, :, :, :) :: volflux1, volflux2, volflux3
     !> The absolute flux of the volumetric flow (volflux1 + volflux2 + volflux3)
     real(kind=rp), allocatable, dimension(:, :, :, :) :: absvolflux
     !> Expanded array of facet normal (zero inside each element)
     real(kind=rp), pointer, dimension(:, :, :, :) :: n1, n2, n3
     !> Number of facet in elements and its maximum
     integer, allocatable :: n_facet(:)
     integer :: n_facet_max
     !> Length scale for element regarding a facet
     real(kind=rp), allocatable :: h(:, :)
     !> Polynomial evaluated at collocation points
     real(kind=rp), allocatable :: phi(:, :)
     !> The first derivative of polynomial at two ends of the interval
     real(kind=rp), allocatable :: dphidxi(:, :)
     !> could be build at initialisation as a tuple_i4_t. (facet, nelv) -> (neigh_facet, neigh_el)
     type(tuple_i4_t), allocatable, dimension(:, :) :: facet_neigh_index

  contains
     !> Constructor.
     procedure, pass(this) :: init => gradient_jump_penalty_init
     !> Destructor
     procedure, pass(this) :: free => gradient_jump_penalty_free
     !> Compute gradient jump penalty term.
     procedure, pass(this) :: compute => gradient_jump_penalty_compute
     !> Perform gradient jump penalty term.
     procedure, pass(this) :: perform => gradient_jump_penalty_perform

  end type gradient_jump_penalty_t

contains
  !> Constructor.
  !! @param dofmap SEM map of degrees of freedom.
  !! @param coef SEM coefficients.
  subroutine gradient_jump_penalty_init(this, dofmap, coef)
    implicit none
    class(gradient_jump_penalty_t), intent(inout) :: this
    type(dofmap_t), intent(in) :: dofmap
    type(coef_t), target, intent(in) :: coef
    
    class(element_t), pointer :: ep
    integer :: temp_indices(12), i, j
    real(kind=rp), allocatable :: zg(:) ! Quadrature points
    
    call this%free()

    this%p = dofmap%xh%lx - 1
    this%lx = dofmap%xh%lx
    if (this%p .gt. 1) then
       this%tau = -0.8_rp * (this%p + 1) ** (-4.0_rp)
    else
       this%tau = -0.02_rp
    end if

    this%coef => coef

    allocate(this%penalty_el(this%p + 1, this%p + 1 , this%p + 1))
    
    allocate(this%n_facet(this%coef%msh%nelv))
    do i = 1, this%coef%msh%nelv
       ep => this%coef%msh%elements(i)%e
       select type(ep)
       type is (hex_t)
          this%n_facet(i) = 6
       type is (quad_t)
          call neko_error("Only Hexahedral element is supported now for gradient jump penalty")
       end select
    end do
    this%n_facet_max = maxval(this%n_facet)

    allocate(this%h(this%n_facet_max, this%coef%msh%nelv))
    do i = 1, this%coef%msh%nelv
       ep => this%coef%msh%elements(i)%e
       select type(ep)
       type is (hex_t)
          call eval_h_hex(this%h(:, i), ep)
       type is (quad_t)
          call neko_error("Gradient jump penalty error: mesh size evaluation is not supported for quad_t")
       end select
    end do
    
    allocate(zg(this%lx))
    allocate(this%phi(this%p + 1, this%p + 1))
    allocate(this%dphidxi(this%p + 1, this%p + 1))

    zg = dofmap%xh%zg(:,1)
    do i = 1, dofmap%xh%lx
       do j = 1, dofmap%xh%lx
          this%dphidxi(j,i) = this%coef%Xh%dx(j,i)
       end do
    end do

    call neko_scratch_registry%request_field(this%penalty, temp_indices(1))
    call neko_scratch_registry%request_field(this%grad_1, temp_indices(2))
    call neko_scratch_registry%request_field(this%grad_2, temp_indices(3))
    call neko_scratch_registry%request_field(this%grad_3, temp_indices(4))

    allocate(this%G(this%lx, this%lx, this%n_facet_max, this%coef%msh%nelv))
    allocate(this%flux1(this%lx, this%lx, this%n_facet_max, this%coef%msh%nelv))
    allocate(this%flux2(this%lx, this%lx, this%n_facet_max, this%coef%msh%nelv))
    allocate(this%flux3(this%lx, this%lx, this%n_facet_max, this%coef%msh%nelv))
    allocate(this%volflux1(this%lx, this%lx, this%n_facet_max, this%coef%msh%nelv))
    allocate(this%volflux2(this%lx, this%lx, this%n_facet_max, this%coef%msh%nelv))
    allocate(this%volflux3(this%lx, this%lx, this%n_facet_max, this%coef%msh%nelv))
    allocate(this%absvolflux(this%lx, this%lx, this%n_facet_max, this%coef%msh%nelv))

    this%n1 => this%coef%nx
    this%n2 => this%coef%ny
    this%n3 => this%coef%nz

    ! build an 2d array to find the neighbouring facet of a specified facet
    allocate(this%facet_neigh_index(this%n_facet_max, this%coef%msh%nelv))
    do i = 1, this%coef%msh%nelv
       ep => this%coef%msh%elements(i)%e
       select type(ep)
       type is (hex_t)
          call facet_matching_hex(this, i)
       type is (quad_t)
          call neko_error("Gradient jump penalty error: neighbouring facet matching is not supported for quad_t")
       end select 
    end do

  end subroutine gradient_jump_penalty_init
  
  !> Matching neighbouring facet for hexahedral mesh
  !!!!! Now only works for 1-core run, MPI needs more discussion
  !! @param i The local index of element
  subroutine facet_matching_hex(this, i)
     class(gradient_jump_penalty_t), intent(inout) :: this
     integer, intent(in) :: i

     integer :: j, k, glb_neigh_el_id, local_neigh_el_id

     do j = 1, 6
        glb_neigh_el_id = this%coef%msh%facet_neigh(j, i)
        if (glb_neigh_el_id .gt. this%coef%msh%nelv) then
           call neko_error("facet matching is not implemented for MPI")
        end if
        local_neigh_el_id = glb_neigh_el_id - this%coef%msh%offset_el
        do k = 1, 6
           if (this%coef%msh%facet_neigh(k, local_neigh_el_id) &
              .eq. i + this%coef%msh%offset_el) then
              this%facet_neigh_index(j, i)%x = (/  k, local_neigh_el_id/)
           end if
        end do
     end do

  end subroutine facet_matching_hex
  
  !> Evaluate h for each element for hexahedral mesh
  !! @param h_el The length scale of an element
  !! @param ep The pointer to the element
  subroutine eval_h_hex(h_el, ep)
    real(kind=rp), intent(inout) :: h_el(6)
    type(hex_t), pointer, intent(in) :: ep
    
    integer :: i
    type(point_t), pointer :: p1, p2, p3, p4, p5, p6, p7, p8

   !  !! strategy 1: use the diameter of the hexahedral
   !  do i = 1, 6
   !     h_el(i) = ep%diameter()
   !  end do

    !! strategy 2: hard code it, only works for cuboid mesh
    !! should be refined for distorted mesh as well
    p1 => ep%p(1)
    p2 => ep%p(2)
    p3 => ep%p(3)
    p4 => ep%p(4)
    p5 => ep%p(5)
    p6 => ep%p(6)
    p7 => ep%p(7)
    p8 => ep%p(8)
    h_el(:) = 0.0_rp
    do i = 1, NEKO_HEX_GDIM
       h_el(1) = h_el(1) + (p1%x(i) - p2%x(i))**2
       h_el(2) = h_el(2) + (p1%x(i) - p2%x(i))**2
       h_el(3) = h_el(3) + (p1%x(i) - p3%x(i))**2
       h_el(4) = h_el(4) + (p1%x(i) - p3%x(i))**2
       h_el(5) = h_el(5) + (p1%x(i) - p5%x(i))**2
       h_el(6) = h_el(6) + (p1%x(i) - p5%x(i))**2
    end do
    do i = 1, 6
       h_el(i) = sqrt(h_el(i))
    end do

  end subroutine eval_h_hex

  !> Destructor for the gradient_jump_penalty_t class.
  subroutine gradient_jump_penalty_free(this)
    implicit none
    class(gradient_jump_penalty_t), intent(inout) :: this

    nullify(this%penalty)
    if (allocated(this%penalty_el)) then
       deallocate(this%penalty_el)
    end if
    if (allocated(this%h)) then
       deallocate(this%h)
    end if
    if (allocated(this%n_facet)) then
       deallocate(this%n_facet)
    end if
    if (allocated(this%phi)) then
       deallocate(this%phi)
    end if
    if (allocated(this%dphidxi)) then
       deallocate(this%dphidxi)
    end if
    if (allocated(this%G)) then
       deallocate(this%G)
    end if
    if (allocated(this%flux1)) then
       deallocate(this%flux1)
    end if
    if (allocated(this%flux2)) then
       deallocate(this%flux2)
    end if
    if (allocated(this%flux3)) then
       deallocate(this%flux3)
    end if
     if (allocated(this%volflux1)) then
       deallocate(this%volflux1)
    end if
    if (allocated(this%volflux2)) then
       deallocate(this%volflux2)
    end if
    if (allocated(this%volflux3)) then
       deallocate(this%volflux3)
    end if
    if (allocated(this%absvolflux)) then
       deallocate(this%absvolflux)
    end if
    nullify(this%coef)
    nullify(this%grad_1)
    nullify(this%grad_2)
    nullify(this%grad_3)
    nullify(this%n1)
    nullify(this%n2)
    nullify(this%n3)

  end subroutine gradient_jump_penalty_free

  !> Compute the gradient jump penalty term.
  !! @param u x-velocity
  !! @param v y-velocity
  !! @param w z-velocity
  !! @param s quantity of interest
  subroutine gradient_jump_penalty_compute(this, u, v, w, s)
    class(gradient_jump_penalty_t), intent(inout) :: this
    type(field_t), intent(in) :: u, v, w, s

    class(element_t), pointer :: ep
    integer :: i

    call G_compute(this, s)
    call absvolflux_compute(this, u, v, w)
    do i = 1, this%coef%msh%nelv
       ep => this%coef%msh%elements(i)%e
       select type(ep)
       type is (hex_t)
          call gradient_jump_penalty_compute_hex_el(this, u, v, w, s, i)
          this%penalty%x(:, :, :, i) = this%penalty_el
       type is (quad_t)
          call neko_error("Only Hexahedral element is supported now for gradient jump penalty")
       end select
    end do

  end subroutine gradient_jump_penalty_compute

  !> Assign the gradient jump penalty term.
  !! @param f A field object to store RHS terms in the weak form equation.
  subroutine gradient_jump_penalty_perform(this, f)
    class(gradient_jump_penalty_t), intent(inout) :: this
    type(field_t), intent(inout) :: f

    call add2(f%x, this%penalty%x, this%coef%dof%size())

  end subroutine gradient_jump_penalty_perform

  !> Compute the gradient jump penalty term for a single hexatedral element.
  !> <tau * h^2 * abs(u .dot. n) * G * phij * phik * dphi_idxi * dxidn>
  !! @param u x-velocity
  !! @param v y-velocity
  !! @param w z-velocity
  !! @param s The quantity of interest
  !! @param i_el The index of the element
  subroutine gradient_jump_penalty_compute_hex_el(this, u, v, w, s, i_el)
    class(gradient_jump_penalty_t), intent(inout) :: this
    type(field_t), intent(in) :: u, v, w, s
    integer, intent(in) :: i_el

    real(kind=rp) :: integrant_facet(this%p + 1, this%p + 1)
    integer :: i, j, k, l

    do i = 1, this%p + 1
       do j = 1, this%p + 1
          do k = 1, this%p + 1
             this%penalty_el(i, j, k) = 0.0_rp
             do l = 1, 6
                this%penalty_el(i, j, k) = this%penalty_el(i, j, k) + &
                                  weak_integrate_over_facet(this, i, j, k, l, i_el)
             end do
          end do
       end do  
    end do

  end subroutine gradient_jump_penalty_compute_hex_el

  !> Compute the average of the flux over facets
  !! @param s The quantity of interest
  subroutine G_compute(this, s)
    class(gradient_jump_penalty_t), intent(inout) :: this
    type(field_t), intent(in) :: s

    integer :: i, j
    real(kind=rp) :: tmp(this%lx, this%lx, this%n_facet_max, this%coef%msh%nelv)
    type(tuple_i4_t) :: i_neigh

    call dudxyz(this%grad_1%x, s%x, this%coef%drdx, this%coef%dsdx, this%coef%dtdx, this%coef)
    call dudxyz(this%grad_2%x, s%x, this%coef%drdy, this%coef%dsdy, this%coef%dtdy, this%coef)
    call dudxyz(this%grad_3%x, s%x, this%coef%drdz, this%coef%dsdz, this%coef%dtdz, this%coef)

    call pick_facet_value_hex(this%flux1, this%grad_1%x, this%lx, this%coef%msh%nelv)
    call pick_facet_value_hex(this%flux2, this%grad_2%x, this%lx, this%coef%msh%nelv)
    call pick_facet_value_hex(this%flux3, this%grad_3%x, this%lx, this%coef%msh%nelv)

    call col2(this%flux1, this%n1, size(this%n1))
    call col2(this%flux2, this%n2, size(this%n1))
    call col2(this%flux3, this%n3, size(this%n1))

    call add3(this%G, this%flux1, this%flux2, size(this%n1))
    call add2(this%G, this%flux3, size(this%n1))
    
    do i = 1, this%coef%msh%nelv
       do j = 1, this%n_facet_max
          i_neigh = this%facet_neigh_index(j, i)
          tmp(:, :, j, i) = this%G(:, :, j, i) + this%G(:, :, i_neigh%x(1), i_neigh%x(2))
       end do
    end do
    this%G = tmp

  end subroutine G_compute

  !> Compute the average of the flux over facets
  !! @param s The quantity of interest
  subroutine absvolflux_compute(this, u, v, w)
    class(gradient_jump_penalty_t), intent(inout) :: this
    type(field_t), intent(in) :: u, v, w

    integer :: i
    
    call pick_facet_value_hex(this%volflux1, u%x, this%lx, this%coef%msh%nelv)
    call pick_facet_value_hex(this%volflux2, v%x, this%lx, this%coef%msh%nelv)
    call pick_facet_value_hex(this%volflux3, w%x, this%lx, this%coef%msh%nelv)

    call col3(this%volflux1, u%x, this%n1, size(this%n1))
    call col3(this%volflux2, v%x, this%n2, size(this%n1))
    call col3(this%volflux3, w%x, this%n3, size(this%n1))

    call add3(this%absvolflux, this%volflux1, this%volflux2, size(this%n1))
    call add2(this%absvolflux, this%volflux3, size(this%n1))

    do i = 1, size(this%n1)
       this%absvolflux(i, 1, 1, 1) = abs(this%absvolflux(i, 1, 1, 1))
    end do

  end subroutine absvolflux_compute

  !> Pick facet values of a field
  !! @param f_facet The data on facets
  !! @param f_field A particular field
  !! @param lx      The number of collocation points on one dimension
  !! @param nelv    The number of element
  subroutine pick_facet_value_hex(f_facet, f_field, lx, nelv)
    integer, intent(in) :: lx, nelv
    real(kind=rp), intent(in) :: f_field(lx, lx, lx, nelv)
    real(kind=rp), intent(inout) :: f_facet(lx, lx, 6, nelv)

    f_facet(:, :, 1, :) = f_field(1, :, :, :)
    f_facet(:, :, 2, :) = f_field(lx, :, :, :)
    f_facet(:, :, 3, :) = f_field(:, 1, :, :)
    f_facet(:, :, 4, :) = f_field(:, lx, :, :)
    f_facet(:, :, 5, :) = f_field(:, :, 1, :)
    f_facet(:, :, 6, :) = f_field(:, :, lx, :)

  end subroutine pick_facet_value_hex

  !> Integrate over a facet in weak from
  !! @param f Integrant
  pure function weak_integrate_over_facet(this, i, j, k, facet_index, i_el) result(f_int)
    class(gradient_jump_penalty_t), intent(in) :: this
   !  real(kind=rp), intent(in) :: f(this%p + 1, this%p + 1)
    integer, intent(in) :: i, j, k, facet_index, i_el
    real(kind=rp) :: f_int
    
    real(kind=rp) :: f(this%p + 1, this%p + 1), dphidxi, dxidn
    integer :: i_pt, j_pt

   select case (facet_index) ! Identify the facet indexing
    case(1)
        i_pt = j
        j_pt = k
        dphidxi = this%dphidxi(1, i)
        dxidn = this%coef%drdx(1, j, k, i_el) + &
                this%coef%drdy(1, j, k, i_el) + &
                this%coef%drdz(1, j, k, i_el)
    case(2)
        i_pt = j
        j_pt = k
        dphidxi = this%dphidxi(this%p + 1, i)
        dxidn = this%coef%drdx(this%p + 1, j, k, i_el) + &
                this%coef%drdy(this%p + 1, j, k, i_el) + &
                this%coef%drdz(this%p + 1, j, k, i_el)
    case(3)
        i_pt = i
        j_pt = k
        dphidxi = this%dphidxi(1, j)
        dxidn = this%coef%dsdx(i, 1, k, i_el) + &
                this%coef%dsdy(i, 1, k, i_el) + &
                this%coef%dsdz(i, 1, k, i_el)
    case(4)
        i_pt = i
        j_pt = k
        dphidxi = this%dphidxi(this%p + 1, j)
        dxidn = this%coef%dsdx(i, this%p + 1, k, i_el) + &
                this%coef%dsdy(i, this%p + 1, k, i_el) + &
                this%coef%dsdz(i, this%p + 1, k, i_el)
    case(5)
        i_pt = i
        j_pt = j
        dphidxi = this%dphidxi(1, k)
        dxidn = this%coef%dtdx(i, j, 1, i_el) + &
                this%coef%dtdy(i, j, 1, i_el) + &
                this%coef%dtdz(i, j, 1, i_el)
    case(6)
        i_pt = i
        j_pt = j
        dphidxi = this%dphidxi(this%p + 1, k)
        dxidn = this%coef%dtdx(i, j, this%p + 1, i_el) + &
                this%coef%dtdy(i, j, this%p + 1, i_el) + &
                this%coef%dtdz(i, j, this%p + 1, i_el)
    end select
    f_int = this%absvolflux(i_pt, j_pt, facet_index, i_el) * &
            this%G(i_pt, j_pt, facet_index, i_el) * &
            dphidxi * dxidn 
    f_int = f_int * this%tau * this%h(facet_index, i_el) ** 2
    f_int = f_int * this%coef%area(i_pt, j_pt, facet_index, i_el)

  end function weak_integrate_over_facet

end module gradient_jump_penalty