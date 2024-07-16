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
  use utils, only : neko_error
  use math
  use point, only : point_t
  use field, only : field_t
  use dofmap , only : dofmap_t
  use neko_config, only : NEKO_BCKND_DEVICE
  use coefs, only : coef_t
  use element, only : element_t
  use hex
  use quad
  use operators, only : dudxyz
  use gs_ops, only : GS_OP_ADD
  use space
  use gather_scatter, only : gs_t
  use device
  use device_math

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
     real(kind=rp), allocatable, dimension(:, :, :, :) :: penalty
     type(c_ptr) :: penalty_d = C_NULL_PTR
     !> work array storing integrant of the penalty
     real(kind=rp), allocatable, dimension(:, :, :, :) :: penalty_facet
     type(c_ptr) :: penalty_facet_d = C_NULL_PTR
     !> SEM coefficients.
     type(coef_t), pointer :: coef => null()
     !> Gradient of the quatity of interest
     real(kind=rp), allocatable, dimension(:, :, :, :) :: grad1, &
                                                          grad2, grad3
     type(c_ptr) :: grad1_d = C_NULL_PTR
     type(c_ptr) :: grad2_d = C_NULL_PTR
     type(c_ptr) :: grad3_d = C_NULL_PTR
     !> Gradient jump on the elementary interface (zero inside each element)
     real(kind=rp), allocatable, dimension(:, :, :, :) :: G
     type(c_ptr) :: G_d = C_NULL_PTR
     !> 3 parts of the flux of the quantity
     real(kind=rp), allocatable, dimension(:, :, :, :) :: flux1, flux2, flux3
     type(c_ptr) :: flux1_d = C_NULL_PTR
     type(c_ptr) :: flux2_d = C_NULL_PTR
     type(c_ptr) :: flux3_d = C_NULL_PTR
     !> 3 parts of the flux of the volumetric flow
     real(kind=rp), allocatable, dimension(:, :, :, :) :: volflux1, &
                                                          volflux2, volflux3
     type(c_ptr) :: volflux1_d = C_NULL_PTR
     type(c_ptr) :: volflux2_d = C_NULL_PTR
     type(c_ptr) :: volflux3_d = C_NULL_PTR
     !> The absolute flux of the volumetric flow
     real(kind=rp), allocatable, dimension(:, :, :, :) :: absvolflux
     type(c_ptr) :: absvolflux_d = C_NULL_PTR
     !> Expanded array of facet normal (zero inside each element)
     real(kind=rp), allocatable, dimension(:, :, :, :) :: n1, n2, n3
     type(c_ptr) :: n1_d = C_NULL_PTR
     type(c_ptr) :: n2_d = C_NULL_PTR
     type(c_ptr) :: n3_d = C_NULL_PTR
     !> Facet factor: quantities in the integrant fix with time stepping
     real(kind=rp), allocatable, dimension(:, :, :, :) :: facet_factor
     type(c_ptr) :: facet_factor_d = C_NULL_PTR
     !> Number of facet in elements and its maximum
     integer, allocatable :: n_facet(:)
     integer :: n_facet_max
     !> Length scale for element regarding a facet
     real(kind=rp), allocatable, dimension(:, :, :, :) :: h2
     !> The first derivative of polynomial at two ends of the interval
     real(kind=rp), allocatable :: dphidxi(:, :)
     !> gather-scattering operation related variables
     type(space_t) :: Xh_GJP !< needed to init gs
     type(dofmap_t) :: dm_GJP !< needed to init gs
     type(gs_t) :: gs_GJP
     !> Size of fields of size lx ** 3 * nelv
     integer :: n
     !> Size of fields of size lx ** 3 * nelv
     integer :: n_large

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
    integer :: i, j
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
    this%n = this%lx ** 3 * this%coef%msh%nelv
    this%n_large = (this%lx + 2) ** 3 * this%coef%msh%nelv

    allocate(this%n_facet(this%coef%msh%nelv))
    do i = 1, this%coef%msh%nelv
       ep => this%coef%msh%elements(i)%e
       select type (ep)
       type is (hex_t)
          this%n_facet(i) = 6
       type is (quad_t)
          call neko_error("Only Hexahedral element is &
                           &supported now for gradient jump penalty")
       end select
    end do
    this%n_facet_max = maxval(this%n_facet)

    allocate(this%h2(this%lx + 2, this%lx + 2, &
                        this%lx + 2, this%coef%msh%nelv))

    do i = 1, this%coef%msh%nelv
       ep => this%coef%msh%elements(i)%e
       select type (ep)
       type is (hex_t)
          call eval_h2_hex(this%h2(:, :, :, i), this%lx, ep)
       type is (quad_t)
          call neko_error("Gradient jump penalty error: mesh size &
                            &evaluation is not supported for quad_t")
       end select
    end do

    allocate(zg(this%lx))
    allocate(this%dphidxi(this%lx, this%lx))

    zg = dofmap%xh%zg(:,1)
    do i = 1, dofmap%xh%lx
       do j = 1, dofmap%xh%lx
          this%dphidxi(j,i) = this%coef%Xh%dx(j,i)
       end do
    end do

    allocate(this%penalty(this%lx, this%lx, this%lx, this%coef%msh%nelv))
    allocate(this%grad1(this%lx, this%lx, this%lx, this%coef%msh%nelv))
    allocate(this%grad2(this%lx, this%lx, this%lx, this%coef%msh%nelv))
    allocate(this%grad3(this%lx, this%lx, this%lx, this%coef%msh%nelv))

    allocate(this%penalty_facet(this%lx + 2, this%lx + 2, &
                    this%lx + 2, this%coef%msh%nelv))
    allocate(this%G(this%lx + 2, this%lx + 2, &
                    this%lx + 2, this%coef%msh%nelv))
    allocate(this%flux1(this%lx + 2, this%lx + 2, &
                        this%lx + 2, this%coef%msh%nelv))
    allocate(this%flux2(this%lx + 2, this%lx + 2, &
                        this%lx + 2, this%coef%msh%nelv))
    allocate(this%flux3(this%lx + 2, this%lx + 2, &
                        this%lx + 2, this%coef%msh%nelv))
    allocate(this%volflux1(this%lx + 2, this%lx + 2, &
                        this%lx + 2, this%coef%msh%nelv))
    allocate(this%volflux2(this%lx + 2, this%lx + 2, &
                        this%lx + 2, this%coef%msh%nelv))
    allocate(this%volflux3(this%lx + 2, this%lx + 2, &
                        this%lx + 2, this%coef%msh%nelv))
    allocate(this%absvolflux(this%lx + 2, this%lx + 2, &
                        this%lx + 2, this%coef%msh%nelv))
    allocate(this%n1(this%lx + 2, this%lx + 2, &
                        this%lx + 2, this%coef%msh%nelv))
    allocate(this%n2(this%lx + 2, this%lx + 2, &
                        this%lx + 2, this%coef%msh%nelv))
    allocate(this%n3(this%lx + 2, this%lx + 2, &
                        this%lx + 2, this%coef%msh%nelv))

    ! Extract facets' normals
    this%n1(1, 2: this%lx + 1, 2: this%lx + 1, :) = this%coef%nx(:, :, 1, :)
    this%n2(1, 2: this%lx + 1, 2: this%lx + 1, :) = this%coef%ny(:, :, 1, :)
    this%n3(1, 2: this%lx + 1, 2: this%lx + 1, :) = this%coef%nz(:, :, 1, :)

    this%n1(this%lx + 2, 2: this%lx + 1, 2: this%lx + 1, :) = &
                                                    this%coef%nx(:, :, 2, :)
    this%n2(this%lx + 2, 2: this%lx + 1, 2: this%lx + 1, :) = &
                                                    this%coef%ny(:, :, 2, :)
    this%n3(this%lx + 2, 2: this%lx + 1, 2: this%lx + 1, :) = &
                                                    this%coef%nz(:, :, 2, :)

    this%n1(2: this%lx + 1, 1, 2: this%lx + 1, :) = this%coef%nx(:, :, 3, :)
    this%n2(2: this%lx + 1, 1, 2: this%lx + 1, :) = this%coef%ny(:, :, 3, :)
    this%n3(2: this%lx + 1, 1, 2: this%lx + 1, :) = this%coef%nz(:, :, 3, :)

    this%n1(2: this%lx + 1, this%lx + 2, 2: this%lx + 1, :) = &  
                                                    this%coef%nx(:, :, 4, :)
    this%n2(2: this%lx + 1, this%lx + 2, 2: this%lx + 1, :) = &
                                                    this%coef%ny(:, :, 4, :)
    this%n3(2: this%lx + 1, this%lx + 2, 2: this%lx + 1, :) = &
                                                    this%coef%nz(:, :, 4, :)

    this%n1(2: this%lx + 1, 2: this%lx + 1, 1, :) = this%coef%nx(:, :, 5, :)
    this%n2(2: this%lx + 1, 2: this%lx + 1, 1, :) = this%coef%ny(:, :, 5, :)
    this%n3(2: this%lx + 1, 2: this%lx + 1, 1, :) = this%coef%nz(:, :, 5, :)

    this%n1(2: this%lx + 1, 2: this%lx + 1, this%lx + 2, :) = &
                                                    this%coef%nx(:, :, 6, :)
    this%n2(2: this%lx + 1, 2: this%lx + 1, this%lx + 2, :) = &
                                                    this%coef%ny(:, :, 6, :)
    this%n3(2: this%lx + 1, 2: this%lx + 1, this%lx + 2, :) = &
                                                    this%coef%nz(:, :, 6, :)

    ! Assemble facet factor
    call facet_factor_init(this)

    ! Initialize Gather-Scatter
    call this%Xh_GJP%init(GLL, this%lx+2, this%lx+2, this%lx+2)
    this%dm_GJP = dofmap_t(this%coef%msh, this%Xh_GJP)
    call this%gs_GJP%init(this%dm_GJP)

    ! Initialize pointers for GPU
    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_map(this%penalty, this%penalty_d, this%n)
       call device_map(this%grad1, this%grad1_d, this%n)
       call device_map(this%grad2, this%grad2_d, this%n)
       call device_map(this%grad3, this%grad3_d, this%n)
       
       call device_map(this%penalty_facet, this%penalty_facet_d, this%n_large)
       call device_map(this%G, this%G_d, this%n_large)
       call device_map(this%flux1, this%flux1_d, this%n_large)
       call device_map(this%flux2, this%flux2_d, this%n_large)
       call device_map(this%flux3, this%flux3_d, this%n_large)

       call device_map(this%volflux1, this%volflux1_d, this%n_large)
       call device_map(this%volflux2, this%volflux2_d, this%n_large)
       call device_map(this%volflux3, this%volflux3_d, this%n_large)
       call device_map(this%absvolflux, this%absvolflux_d, this%n_large)

       call device_map(this%n1, this%n1_d, this%n_large)
       call device_map(this%n2, this%n2_d, this%n_large)
       call device_map(this%n3, this%n3_d, this%n_large)
       call device_map(this%facet_factor, this%facet_factor_d, this%n_large)

       call device_memcpy(this%n1, this%n1_d, this%n_large, &
                          HOST_TO_DEVICE, sync=.false.)
       call device_memcpy(this%n2, this%n2_d, this%n_large, &
                          HOST_TO_DEVICE, sync=.false.)
       call device_memcpy(this%n3, this%n3_d, this%n_large, &
                          HOST_TO_DEVICE, sync=.false.)
       call device_memcpy(this%facet_factor, this%facet_factor_d, this%n_large,&
                          HOST_TO_DEVICE, sync = .false.)

    end if

  end subroutine gradient_jump_penalty_init
  
  !> Evaluate h^2 for each element for hexahedral mesh
  !! @param h2_el The sqaure of the length scale of an element
  !! @param ep The pointer to the element
  subroutine eval_h2_hex(h2_el, n, ep)
    integer, intent(in) :: n
    real(kind=rp), intent(inout) :: h2_el(n + 2, n + 2, n + 2)
    type(hex_t), pointer, intent(in) :: ep
    
    integer :: i
    type(point_t), pointer :: p1, p2, p3, p4, p5, p6, p7, p8

    p1 => ep%p(1)
    p2 => ep%p(2)
    p3 => ep%p(3)
    p4 => ep%p(4)
    p5 => ep%p(5)
    p6 => ep%p(6)
    p7 => ep%p(7)
    p8 => ep%p(8)
    h2_el = 0.0_rp
    
    h2_el(1, 2 : n + 1, 2 : n + 1) = dist2_facets_hex(p1, p5, p7, p3, &
                                                     p2, p6, p8, p4)
    h2_el(n + 2, 2 : n + 1, 2 : n + 1) = h2_el(1, 2, 2)

    h2_el(2 : n + 1, 1, 2 : n + 1) = dist2_facets_hex(p1, p2, p6, p5, &
                                                     p3, p4, p8, p7)
    h2_el(2 : n + 1, n + 2, 2 : n + 1) = h2_el(2, 1, 2)

    h2_el(2 : n + 1, 2 : n + 1, 1) = dist2_facets_hex(p1, p2, p4, p3, &
                                                     p5, p6, p8, p7)
    h2_el(2 : n + 1, 2 : n + 1, n + 2) = h2_el(2, 2, 1)

  end subroutine eval_h2_hex

  !> "Distrance" of two facets of a hexahedral element
  !! @param p11, p12, p13, p14 Vertices defining facet 1
  !! @param p21, p22, p23, p24 Vertices defining facet 2
  function dist2_facets_hex(p11, p12, p13, p14, &
                           p21, p22, p23, p24) result(dist2)
  type(point_t), intent(in) :: p11, p12, p13, p14
  type(point_t), intent(in) :: p21, p22, p23, p24
  real(kind=rp) :: dist2

  dist2 = 0.0_rp

  dist2 = dist2 + dist_vertex_to_plane_3d(p11, p21, p22, p23)
  dist2 = dist2 + dist_vertex_to_plane_3d(p11, p21, p22, p24)
  dist2 = dist2 + dist_vertex_to_plane_3d(p11, p21, p23, p24)
  dist2 = dist2 + dist_vertex_to_plane_3d(p11, p22, p23, p24)

  dist2 = dist2 + dist_vertex_to_plane_3d(p12, p21, p22, p23)
  dist2 = dist2 + dist_vertex_to_plane_3d(p12, p21, p22, p24)
  dist2 = dist2 + dist_vertex_to_plane_3d(p12, p21, p23, p24)
  dist2 = dist2 + dist_vertex_to_plane_3d(p12, p22, p23, p24)

  dist2 = dist2 + dist_vertex_to_plane_3d(p13, p21, p22, p23)
  dist2 = dist2 + dist_vertex_to_plane_3d(p13, p21, p22, p24)
  dist2 = dist2 + dist_vertex_to_plane_3d(p13, p21, p23, p24)
  dist2 = dist2 + dist_vertex_to_plane_3d(p13, p22, p23, p24)

  dist2 = dist2 + dist_vertex_to_plane_3d(p14, p21, p22, p23)
  dist2 = dist2 + dist_vertex_to_plane_3d(p14, p21, p22, p24)
  dist2 = dist2 + dist_vertex_to_plane_3d(p14, p21, p23, p24)
  dist2 = dist2 + dist_vertex_to_plane_3d(p14, p22, p23, p24)

  dist2 = dist2 / 16.0_rp
  dist2 = dist2 * dist2

  end function dist2_facets_hex

  !> Distance from a vertex to a plane in a 3d space
  !! @param pv The vertex
  !! @param p1, p2, p3 Points defining a plane
  function dist_vertex_to_plane_3d(pv, p1, p2, p3) result(dist)
  type(point_t), intent(in) :: pv, p1, p2, p3
  !> Vectors connecting p1&p2, p1&p3 and pv&p1
  real(kind=rp), dimension(3) :: u12, u13, uv1
  real(kind=rp) :: norm_u12, norm_u13
  !> The vector normal to the plane
  real(kind=rp) :: un(3)
  !> The result of this function
  real(kind=rp) :: dist
  integer :: i
  
  ! Set up u12 and u13
  norm_u12 = 0.0_rp
  norm_u13 = 0.0_rp
  do i = 1, 3
     u12(i) = p2%x(i) - p1%x(i)
     norm_u12 = norm_u12 + u12(i) * u12(i)
     u13(i) = p3%x(i) - p1%x(i)
     norm_u13 = norm_u13 + u13(i) * u13(i)
     uv1(i) = pv%x(i) - p1%x(i)
  end do
  norm_u12 = sqrt(norm_u12)
  norm_u13 = sqrt(norm_u13)

  ! Normalized cross product of u12 and u13 to get un
  un(1) = u12(2) * u13(3) - u12(3) * u13(2)
  un(2) = - u12(1) * u13(3) + u12(3) * u13(1)
  un(3) = u12(1) * u13(2) - u12(2) * u13(1)
  do i = 1, 3
     un(i) = un(i) / norm_u12 / norm_u13
  end do

  ! Project of uv1 onto un, noting un is a unit vector
  dist = 0.0_rp
  do i = 1, 3
     dist = dist + uv1(i) * un(i)
  end do

  end function dist_vertex_to_plane_3d

  !> Initialize the facet factor array
  subroutine facet_factor_init(this)
    class(gradient_jump_penalty_t), intent(inout) :: this
    !> work array
    real(kind=rp) :: wa(this%lx, this%lx, this%lx, this%coef%msh%nelv)

    allocate(this%facet_factor(this%lx + 2, this%lx + 2, &
                               this%lx + 2, this%coef%msh%nelv))

    associate(facet_factor => this%facet_factor, &
              lx => this%lx, area => this%coef%area, &
              nelv => this%coef%msh%nelv, &
              jacinv => this%coef%jacinv, n => this%n, &
              n_large => this%n_large, h2 => this%h2, &
              tau => this%tau)
    
    ! Assemble facet_factor for facet 1 and 2
    call add4(wa, this%coef%drdx, this%coef%drdy, this%coef%drdz, n)
    call col2(wa, jacinv, n)
    facet_factor(1, 2: lx + 1, 2: lx + 1, :) = -1.0_rp * wa(1, :, :, :)
    facet_factor(lx + 2, 2: lx + 1, 2: lx + 1, :) = wa(lx, :, :, :)

    ! Assemble facet_factor for facet 3 and 4
    call add4(wa, this%coef%dsdx, this%coef%dsdy, this%coef%dsdz, n)
    call col2(wa, jacinv, n)
    facet_factor(2: lx + 1, 1, 2: lx + 1, :) = -1.0_rp * wa(:, 1, :, :)
    facet_factor(2: lx + 1, lx + 2, 2: lx + 1, :) = wa(:, lx, :, :)

    ! Assemble facet_factor for facet 5 and 6
    call add4(wa, this%coef%dsdx, this%coef%dsdy, this%coef%dsdz, n)
    call col2(wa, jacinv, n)
    facet_factor(2: lx + 1, 2: lx + 1, 1, :) = -1.0_rp * wa(:, :, 1, :)
    facet_factor(2: lx + 1, 2: lx + 1, lx + 2, :) = wa(:, :, lx, :)

    ! Multiplied by h^2 * tau
    call col2(facet_factor, h2, n_large)
    call cmult(facet_factor, tau, this%n_large)

    ! Multiplied by the quadrant weight
    call col2(facet_factor(1, 2: lx + 1, 2: lx + 1, :), &
              area(:, :, 1, :), lx * lx * nelv)
    call col2(facet_factor(lx + 2, 2: lx + 1, 2: lx + 1, :), &
              area(:, :, 2, :), lx * lx * nelv)
    call col2(facet_factor(2: lx + 1, 1, 2: lx + 1, :), &
              area(:, :, 3, :), lx * lx * nelv)
    call col2(facet_factor(2: lx + 1, lx + 2, 2: lx + 1, :), &
              area(:, :, 4, :), lx * lx * nelv)
    call col2(facet_factor(2: lx + 1, 2: lx + 1, 1, :), &
              area(:, :, 5, :), lx * lx * nelv)
    call col2(facet_factor(2: lx + 1, 2: lx + 1, lx + 2, :), &
              area(:, :, 6, :), lx * lx * nelv)

    end associate
  end subroutine facet_factor_init

  !> Destructor for the gradient_jump_penalty_t class.
  subroutine gradient_jump_penalty_free(this)
    implicit none
    class(gradient_jump_penalty_t), intent(inout) :: this

    if (c_associated(this%penalty_d)) then
       call device_free(this%penalty_d)
    end if
    if (c_associated(this%grad1_d)) then
       call device_free(this%grad1_d)
    end if
    if (c_associated(this%grad2_d)) then
       call device_free(this%grad2_d)
    end if
    if (c_associated(this%grad3_d)) then
       call device_free(this%grad3_d)
    end if
    if (c_associated(this%penalty_facet_d)) then
       call device_free(this%penalty_facet_d)
    end if
    if (c_associated(this%G_d)) then
       call device_free(this%G_d)
    end if
    if (c_associated(this%flux1_d)) then
       call device_free(this%flux1_d)
    end if
    if (c_associated(this%flux2_d)) then
       call device_free(this%flux2_d)
    end if
    if (c_associated(this%flux3_d)) then
       call device_free(this%flux3_d)
    end if
    if (c_associated(this%volflux1_d)) then
       call device_free(this%volflux1_d)
    end if
    if (c_associated(this%volflux2_d)) then
       call device_free(this%volflux2_d)
    end if
    if (c_associated(this%volflux3_d)) then
       call device_free(this%volflux3_d)
    end if
    if (c_associated(this%absvolflux_d)) then
       call device_free(this%absvolflux_d)
    end if
    if (c_associated(this%n1_d)) then
       call device_free(this%n1_d)
    end if
    if (c_associated(this%n2_d)) then
       call device_free(this%n2_d)
    end if
    if (c_associated(this%n3_d)) then
       call device_free(this%n3_d)
    end if
    if (c_associated(this%facet_factor_d)) then
       call device_free(this%facet_factor_d)
    end if

    if (allocated(this%penalty)) then
       deallocate(this%penalty)
    end if
    if (allocated(this%grad1)) then
       deallocate(this%grad1)
    end if
    if (allocated(this%grad2)) then
       deallocate(this%grad2)
    end if
    if (allocated(this%grad3)) then
       deallocate(this%grad3)
    end if
    if (allocated(this%h2)) then
       deallocate(this%h2)
    end if
    if (allocated(this%n_facet)) then
       deallocate(this%n_facet)
    end if
    if (allocated(this%dphidxi)) then
       deallocate(this%dphidxi)
    end if
    if (allocated(this%penalty_facet)) then
       deallocate(this%penalty_facet)
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
    if (allocated(this%n1)) then
       deallocate(this%n1)
    end if
    if (allocated(this%n2)) then
       deallocate(this%n2)
    end if
    if (allocated(this%n3)) then
       deallocate(this%n3)
    end if
    if (allocated(this%facet_factor)) then
       deallocate(this%facet_factor)
    end if

    nullify(this%coef)
    
    call this%Xh_GJP%free()
    call this%gs_GJP%free()

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

    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_col3(this%penalty_facet_d, this%absvolflux_d, this%G_d, this%n_large)
       call device_col2(this%penalty_facet_d, this%facet_factor_d, this%n_large)
    else
       call col3(this%penalty_facet, this%absvolflux, this%G, this%n_large)
       call col2(this%penalty_facet, this%facet_factor, this%n_large)
    end if
    
    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_memcpy(this%penalty_facet, this%penalty_facet_d, this%n_large, &
                          DEVICE_TO_HOST, sync=.true.)
    end if

    do i = 1, this%coef%msh%nelv
       ep => this%coef%msh%elements(i)%e
       select type (ep)
       type is (hex_t)
          call gradient_jump_penalty_compute_hex_el(this, this%penalty_facet, i)
       type is (quad_t)
          call neko_error("Only Hexahedral element is supported &
                                       now for gradient jump penalty")
       end select
    end do

    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_memcpy(this%penalty, this%penalty_d, this%n, &
                          HOST_TO_DEVICE, sync = .true.)
    end if

  end subroutine gradient_jump_penalty_compute

  !> Assign the gradient jump penalty term.
  !! @param f A field object to store RHS terms in the weak form equation.
  subroutine gradient_jump_penalty_perform(this, f)
    class(gradient_jump_penalty_t), intent(inout) :: this
    type(field_t), intent(inout) :: f

    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_add2(f%x_d, this%penalty_d, this%coef%dof%size())
    else
       call add2(f%x, this%penalty, this%coef%dof%size())
    end if

  end subroutine gradient_jump_penalty_perform

  !> Compute the gradient jump penalty term for a single hexatedral element.
  !> <tau * h^2 * absvolflux * G * phij * phik * dphi_idxi * dxidn>
  !! @param wa work array containing tau * absvolflux * G * dxidn
  !! @param i_el The index of the element
  subroutine gradient_jump_penalty_compute_hex_el(this, wa, i_el)
    class(gradient_jump_penalty_t), intent(inout) :: this
    integer, intent(in) :: i_el
    real(kind=rp), intent(in) :: wa(this%lx + 2, this%lx + 2, &
                        this%lx + 2, this%coef%msh%nelv)

    integer :: i, j, k

    associate(lx => this%lx, &
              penalty => this%penalty, dphidxi => this%dphidxi)

    do i = 1, lx
       do j = 1, lx
          do k = 1, lx
             penalty(i, j, k, i_el) = & 
               wa(1, j + 1, k + 1, i_el) * &
                  dphidxi(1, i) + &
               wa(lx + 2, j + 1, k + 1, i_el) * &
                  dphidxi(lx, i) + &
               wa(i + 1, 1, k + 1, i_el) * &
                  dphidxi(1, j) + &
               wa(i + 1, lx + 2, k + 1, i_el) * &
                  dphidxi(lx, j) + &
               wa(i + 1, j + 1, 1, i_el) * &
                  dphidxi(1, k) + &
               wa(i + 1, j + 1, lx + 2, i_el) * &
                  dphidxi(lx, k)
          end do
       end do
    end do

    end associate
  end subroutine gradient_jump_penalty_compute_hex_el

  !> Compute the average of the flux over facets
  !! @param s The quantity of interest
  subroutine G_compute(this, s)
    class(gradient_jump_penalty_t), intent(inout) :: this
    type(field_t), intent(in) :: s

    call dudxyz(this%grad1, s%x, this%coef%drdx, &
                this%coef%dsdx, this%coef%dtdx, this%coef)
    call dudxyz(this%grad2, s%x, this%coef%drdy, &
                this%coef%dsdy, this%coef%dtdy, this%coef)
    call dudxyz(this%grad3, s%x, this%coef%drdz, &
                this%coef%dsdz, this%coef%dtdz, this%coef)

    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_memcpy(this%grad1, this%grad1_d, this%n, &
                          DEVICE_TO_HOST, sync = .false.)
       call device_memcpy(this%grad2, this%grad2_d, this%n, &
                          DEVICE_TO_HOST, sync = .false.)
       call device_memcpy(this%grad3, this%grad3_d, this%n, &
                          DEVICE_TO_HOST, sync = .true.)
    end if

    call pick_facet_value_hex(this%flux1, this%grad1, &
                              this%lx, this%coef%msh%nelv)
    call pick_facet_value_hex(this%flux2, this%grad2, &
                              this%lx, this%coef%msh%nelv)
    call pick_facet_value_hex(this%flux3, this%grad3, &
                              this%lx, this%coef%msh%nelv)
    
    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_memcpy(this%flux1, this%flux1_d, this%n_large, &
                          HOST_TO_DEVICE, sync = .false.)
       call device_memcpy(this%flux2, this%flux2_d, this%n_large, &
                          HOST_TO_DEVICE, sync = .false.)
       call device_memcpy(this%flux3, this%flux3_d, this%n_large, &
                          HOST_TO_DEVICE, sync = .true.)
    end if
     
    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_col2(this%flux1_d, this%n1_d, this%n_large)
       call device_col2(this%flux2_d, this%n2_d, this%n_large)
       call device_col2(this%flux3_d, this%n3_d, this%n_large)
       call device_add3s2(this%G_d, this%flux1_d, this%flux2_d, &
                          1.0_rp, 1.0_rp, this%n_large)
       call device_add2(this%G_d, this%flux3_d, this%n_large)
    else
       call col2(this%flux1, this%n1, this%n_large)
       call col2(this%flux2, this%n2, this%n_large)
       call col2(this%flux3, this%n3, this%n_large)
       call add3(this%G, this%flux1, this%flux2, this%n_large)
       call add2(this%G, this%flux3, this%n_large)
    end if
    
    call this%gs_GJP%op(this%G, this%n_large, GS_OP_ADD)

  end subroutine G_compute

  !> Compute the average of the volumetric flux over facets
  !! @param u x-velocity
  !! @param v y-velocity
  !! @param w z-velocity
  subroutine absvolflux_compute(this, u, v, w)
    class(gradient_jump_penalty_t), intent(inout) :: this
    type(field_t), intent(in) :: u, v, w

    integer :: i
    
    call pick_facet_value_hex(this%volflux1, u%x, this%lx, this%coef%msh%nelv)
    call pick_facet_value_hex(this%volflux2, v%x, this%lx, this%coef%msh%nelv)
    call pick_facet_value_hex(this%volflux3, w%x, this%lx, this%coef%msh%nelv)

    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_memcpy(this%volflux1, this%volflux1_d, this%n_large, &
                          HOST_TO_DEVICE, sync = .false.)
       call device_memcpy(this%volflux2, this%volflux2_d, this%n_large, &
                          HOST_TO_DEVICE, sync = .false.)
       call device_memcpy(this%volflux3, this%volflux3_d, this%n_large, &
                          HOST_TO_DEVICE, sync = .true.)
    end if
    
    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_col2(this%volflux1_d, this%n1_d, this%n_large)
       call device_col2(this%volflux2_d, this%n2_d, this%n_large)
       call device_col2(this%volflux3_d, this%n3_d, this%n_large)
       call device_add3s2(this%absvolflux_d, this%volflux1_d, this%volflux2_d, &
                          1.0_rp, 1.0_rp, this%n_large)
       call device_add2(this%absvolflux_d, this%volflux3_d, this%n_large)
    else
       call col2(this%volflux1, this%n1, this%n_large)
       call col2(this%volflux2, this%n2, this%n_large)
       call col2(this%volflux3, this%n3, this%n_large)
       call add3(this%absvolflux, this%volflux1, this%volflux2, this%n_large)
       call add2(this%absvolflux, this%volflux3, this%n_large)
    end if

    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_memcpy(this%absvolflux, this%absvolflux_d, this%n_large, &
                          DEVICE_TO_HOST, sync = .true.)
    end if
    do i = 1, this%n_large
       this%absvolflux(i, 1, 1, 1) = abs(this%absvolflux(i, 1, 1, 1))
    end do
    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_memcpy(this%absvolflux, this%absvolflux_d, this%n_large, &
                          HOST_TO_DEVICE, sync = .true.)
    end if

  end subroutine absvolflux_compute

  !> Pick facet values of a field
  !! @param f_facet The data on facets
  !! @param f_field A particular field
  !! @param lx      The number of collocation points on one dimension
  !! @param nelv    The number of element
  subroutine pick_facet_value_hex(f_facet, f_field, lx, nelv)
    integer, intent(in) :: lx, nelv
    real(kind=rp), intent(in) :: f_field(lx, lx, lx, nelv)
    real(kind=rp), intent(inout) :: f_facet(lx + 2, lx + 2, lx + 2, nelv)

    call copy(f_facet(1, 2: lx + 1, 2: lx + 1, :), &
              f_field(1, :, :, :), lx * lx * nelv)
    call copy(f_facet(lx + 2, 2: lx + 1, 2: lx + 1, :), &
              f_field(lx, :, :, :), lx * lx * nelv)
    call copy(f_facet(2: lx + 1, 1, 2: lx + 1, :), &
              f_field(:, 1, :, :), lx * lx * nelv)
    call copy(f_facet(2: lx + 1, lx + 2, 2: lx + 1, :), &
              f_field(:, lx, :, :), lx * lx * nelv)
    call copy(f_facet(2: lx + 1, 2: lx + 1, 1, :), &
              f_field(:, :, 1, :), lx * lx * nelv)
    call copy(f_facet(2: lx + 1, 2: lx + 1, lx + 2, :), &
              f_field(:, :, lx, :), lx * lx * nelv)

  end subroutine pick_facet_value_hex

end module gradient_jump_penalty