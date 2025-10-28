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
  use json_utils, only : json_get_or_default
  use json_module, only : json_file
  use math, only: add2, col2, col3, invcol2, add3, copy, absval
  use point, only : point_t
  use field, only : field_t
  use dofmap , only : dofmap_t
  use neko_config, only : NEKO_BCKND_DEVICE
  use coefs, only : coef_t
  use element, only : element_t
  use hex, only : hex_t
  use quad, only : quad_t
  use gs_ops, only : GS_OP_ADD
  use space, only : space_t, GLL
  use gather_scatter, only : gs_t
  use device, only : device_map, HOST_TO_DEVICE, device_memcpy, device_free
  use device_math, only : device_col2, device_add2, device_add2s2, &
       device_col3, device_add3s2, device_invcol2, device_absval
  use device_gradient_jump_penalty, only : device_pick_facet_value_hex, &
       device_gradient_jump_penalty_finalize
  use source_term, only : source_term_t
  use field_list, only : field_list_t
  use field_registry, only : neko_field_registry
  use time_state, only : time_state_t
  use operators, only : dudxyz
  use, intrinsic :: iso_c_binding, only : c_ptr, C_NULL_PTR, c_associated

  implicit none
  private

  !> Implements the gradient jump penalty.
  !! @note Reference DOI: 10.1016/j.cma.2021.114200
  type, public, extends(source_term_t):: gradient_jump_penalty_t
     !> The solution fields for which to apply the penalty.
     type(field_list_t) :: s_fields
     !> The x component of velocity
     type(field_t), pointer :: u
     !> The y component of velocity
     type(field_t), pointer :: v
     !> The z component of velocity
     type(field_t), pointer :: w

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
     type(c_ptr) :: dphidxi_d = C_NULL_PTR
     !> gather-scattering operation related variables
     type(space_t) :: Xh_GJP !< needed to init gs
     type(dofmap_t) :: dm_GJP !< needed to init gs
     type(gs_t) :: gs_GJP
     !> Size of fields of size lx ** 3 * nelv
     integer :: n
     !> Size of fields of size lx ** 3 * nelv
     integer :: n_large

   contains
     !> Constructor from JSON.
     procedure, pass(this) :: init => gradient_jump_penalty_init
     !> Constructor from components.
     procedure, pass(this) :: init_from_components => &
          gradient_jump_penalty_init_from_components
     !> Destructor.
     procedure, pass(this) :: free => gradient_jump_penalty_free
     !> Compute gradient jump penalty term.
     procedure, pass(this) :: compute_single => &
          gradient_jump_penalty_compute_single
     !> Perform gradient jump penalty term.
     procedure, pass(this) :: compute_ => gradient_jump_penalty_compute

  end type gradient_jump_penalty_t

contains
  !> Constructor.
  !! @param json The JSON object for the source.
  !! @param fields A list of fields for adding the source values.
  !! @param coef The SEM coeffs.
  !! @param a, b Coefficients to determine tau
  !! @param variable_name The name of the variable for this source term.
  subroutine gradient_jump_penalty_init(this, json, fields, coef, variable_name)
    implicit none
    class(gradient_jump_penalty_t), intent(inout) :: this
    type(json_file), intent(inout) :: json
    type(coef_t), target, intent(in) :: coef
    character(len=*), intent(in) :: variable_name
    type(field_list_t), intent(in), target :: fields

    real(kind=rp) :: start_time, end_time
    real(kind=rp) :: a, b

    call json_get_or_default(json, "start_time", start_time, 0.0_rp)
    call json_get_or_default(json, "end_time", end_time, huge(0.0_rp))

    if ((coef%Xh%lx - 1) .eq. 1) then
       call json_get_or_default(json, 'tau', a, 0.02_rp)
       b = 0.0_rp
    else
       call json_get_or_default(json, 'scaling_factor', a, 0.8_rp)
       call json_get_or_default(json, 'scaling_exponent', b, 4.0_rp)
    end if

    call this%init_from_components(coef, fields, start_time, end_time, a, b, &
         variable_name)

  end subroutine gradient_jump_penalty_init

  !> Constructor from components.
  !! @param fields A list of fields for adding the source values.
  !! @param coef The SEM coeffs.
  !! @param a, b Coefficients to determine tau
  !! @param start_time When to start adding the source term.
  !! @param end_time When to stop adding the source term.
  !! @param variable_name The name of the variable for this source term.
  subroutine gradient_jump_penalty_init_from_components(this, coef, fields, &
       start_time, end_time, a, b, variable_name)
    class(gradient_jump_penalty_t), intent(inout) :: this
    type(coef_t), target, intent(in) :: coef
    type(field_list_t), intent(in), target :: fields
    real(kind=rp), intent(in) :: start_time
    real(kind=rp), intent(in) :: end_time
    real(kind=rp), intent(in) :: a, b
    character(len=*), intent(in) :: variable_name

    integer :: i, j, k, l
    real(kind=rp), allocatable :: zg(:) ! Quadrature points
    real(kind=rp) :: normal(3)

    call this%free()

    call this%init_base(fields, coef, start_time, end_time)

    this%u => neko_field_registry%get_field("u")
    this%v => neko_field_registry%get_field("v")
    this%w => neko_field_registry%get_field("w")

    if (fields%size() .eq. 1) then
       call this%s_fields%init(1)
       call this%s_fields%assign(1, &
            neko_field_registry%get_field(variable_name))
    else if (fields%size() .eq. 3) then
       call this%s_fields%init(3)
       call this%s_fields%assign(1, this%u)
       call this%s_fields%assign(2, this%v)
       call this%s_fields%assign(3, this%w)
    else
       call neko_error("The GJP source assumes either 3 or 1 RHS fields.")
    end if




    this%p = coef%dof%Xh%lx - 1
    this%lx = coef%dof%Xh%lx

    if (this%p .gt. 1) then
       this%tau = -a * (this%p + 1) ** (-b)
    else
       this%tau = -a
    end if

    this%n = this%lx ** 3 * this%coef%msh%nelv
    this%n_large = (this%lx + 2) ** 3 * this%coef%msh%nelv

    allocate(this%n_facet(this%coef%msh%nelv))
    do i = 1, this%coef%msh%nelv
       select type (ep => this%coef%msh%elements(i)%e)
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
       select type (ep => this%coef%msh%elements(i)%e)
       type is (hex_t)
          call eval_h2_hex(this%h2(:, :, :, i), this%lx, i, this%coef)
       type is (quad_t)
          call neko_error("Gradient jump penalty error: mesh size &
          &evaluation is not supported for quad_t")
       end select
    end do

    allocate(zg(this%lx))
    allocate(this%dphidxi(this%lx, this%lx))

    zg = coef%Xh%zg(:,1)
    do i = 1, coef%Xh%lx
       do j = 1, coef%Xh%lx
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
    do i = 1, this%coef%msh%nelv
       do j = 1, 6 ! for hexahedral elements
          do k = 1, this%lx
             do l = 1, this%lx
                select case (j)
                case (1)
                   normal = this%coef%get_normal(1, l, k, i, j)
                   this%n1(1, l + 1, k + 1, i) = normal(1)
                   this%n2(1, l + 1, k + 1, i) = normal(2)
                   this%n3(1, l + 1, k + 1, i) = normal(3)
                case (2)
                   normal = this%coef%get_normal(1, l, k, i, j)
                   this%n1(this%lx + 2, l + 1, k + 1, i) = normal(1)
                   this%n2(this%lx + 2, l + 1, k + 1, i) = normal(2)
                   this%n3(this%lx + 2, l + 1, k + 1, i) = normal(3)
                case (3)
                   normal = this%coef%get_normal(l, 1, k, i, j)
                   this%n1(l + 1, 1, k + 1, i) = normal(1)
                   this%n2(l + 1, 1, k + 1, i) = normal(2)
                   this%n3(l + 1, 1, k + 1, i) = normal(3)
                case (4)
                   normal = this%coef%get_normal(l, 1, k, i, j)
                   this%n1(l + 1, this%lx + 2, k + 1, i) = normal(1)
                   this%n2(l + 1, this%lx + 2, k + 1, i) = normal(2)
                   this%n3(l + 1, this%lx + 2, k + 1, i) = normal(3)
                case (5)
                   normal = this%coef%get_normal(l, k, 1, i, j)
                   this%n1(l + 1, k + 1, 1, i) = normal(1)
                   this%n2(l + 1, k + 1, 1, i) = normal(2)
                   this%n3(l + 1, k + 1, 1, i) = normal(3)
                case (6)
                   normal = this%coef%get_normal(l, k, 1, i, j)
                   this%n1(l + 1, k + 1, this%lx + 2, i) = normal(1)
                   this%n2(l + 1, k + 1, this%lx + 2, i) = normal(2)
                   this%n3(l + 1, k + 1, this%lx + 2, i) = normal(3)
                case default
                   call neko_error("The face index is not correct")
                end select
             end do
          end do
       end do
    end do

    ! Assemble facet factor
    call facet_factor_init(this)

    ! Initialize Gather-Scatter
    call this%Xh_GJP%init(GLL, this%lx+2, this%lx+2, this%lx+2)
    call this%dm_GJP%init(this%coef%msh, this%Xh_GJP)
    call this%gs_GJP%init(this%dm_GJP)

    ! Initialize pointers for device
    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_map(this%dphidxi, this%dphidxi_d, &
            this%lx * this%lx)
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

       call device_memcpy(this%dphidxi, this%dphidxi_d, &
            this%lx * this%lx, &
            HOST_TO_DEVICE, sync = .false.)
       call device_memcpy(this%n1, this%n1_d, this%n_large, &
            HOST_TO_DEVICE, sync = .false.)
       call device_memcpy(this%n2, this%n2_d, this%n_large, &
            HOST_TO_DEVICE, sync = .false.)
       call device_memcpy(this%n3, this%n3_d, this%n_large, &
            HOST_TO_DEVICE, sync = .false.)
       call device_memcpy(this%facet_factor, this%facet_factor_d, this%n_large,&
            HOST_TO_DEVICE, sync = .false.)

    end if

  end subroutine gradient_jump_penalty_init_from_components

  !> Evaluate h^2 for each element for hexahedral mesh
  !! @param h2_el The square of the length scale of an element
  !! @param ep The pointer to the element
  subroutine eval_h2_hex(h2_el, n, i, coef)
    integer, intent(in) :: n, i
    type(coef_t), pointer, intent(in) :: coef
    real(kind=rp), intent(inout) :: h2_el(n + 2, n + 2, n + 2)

    type(dofmap_t), pointer :: dm
    integer :: j, k, l

    dm => coef%dof

    h2_el = 0.0_rp

    do j = 1, 6
       do k = 1, n
          do l = 1, n
             select case (j)
             case (1)
                h2_el(1, l + 1, k + 1) = &
                     dist2_quadrature_hex(l, k, j, i, n, dm, coef)
             case (2)
                h2_el(n + 2, l + 1, k + 1) = &
                     dist2_quadrature_hex(l, k, j, i, n, dm, coef)
             case (3)
                h2_el(l + 1, 1, k + 1) = &
                     dist2_quadrature_hex(l, k, j, i, n, dm, coef)
             case (4)
                h2_el(l + 1, n + 2, k + 1) = &
                     dist2_quadrature_hex(l, k, j, i, n, dm, coef)
             case (5)
                h2_el(l + 1, k + 1, 1) = &
                     dist2_quadrature_hex(l, k, j, i, n, dm, coef)
             case (6)
                h2_el(l + 1, k + 1, n + 2) = &
                     dist2_quadrature_hex(l, k, j, i, n, dm, coef)
             case default
                call neko_error("The face index is not correct")
             end select
          end do
       end do
    end do

  end subroutine eval_h2_hex

  function dist2_quadrature_hex(l, k, j, i, n, dm, coef) result(dist2)
    integer, intent(in) :: l, k, j, i, n
    type(dofmap_t), pointer, intent(in) :: dm
    type(coef_t), pointer, intent(in) :: coef
    real(kind=rp) :: dist2, dist_1, dist_2

    real(kind=rp) :: x1, y1, z1, x2, y2, z2
    real(kind=rp) :: normal1(3), normal2(3)
    real(kind=rp) :: n11, n12, n13, n21, n22, n23
    real(kind=rp) :: v1, v2, v3

    dist2 = 0.0_rp
    select case (j)
    case (1)
       normal1 = coef%get_normal(1, l, k, i, 1)
       n11 = normal1(1)
       n12 = normal1(2)
       n13 = normal1(3)
       x1 = dm%x(1, l, k, i)
       y1 = dm%y(1, l, k, i)
       z1 = dm%z(1, l, k, i)
       normal2 = coef%get_normal(1, l, k, i, 2)
       n21 = normal2(1)
       n22 = normal2(2)
       n23 = normal2(3)
       x2 = dm%x(n, l, k, i)
       y2 = dm%y(n, l, k, i)
       z2 = dm%z(n, l, k, i)
    case (2)
       ! now the facet pair share the same value for h
       ! but just let it be here for furture possible changes
       normal1 = coef%get_normal(1, l, k, i, 2)
       n11 = normal1(1)
       n12 = normal1(2)
       n13 = normal1(3)
       x1 = dm%x(n, l, k, i)
       y1 = dm%y(n, l, k, i)
       z1 = dm%z(n, l, k, i)
       normal2 = coef%get_normal(1, l, k, i, 1)
       n21 = normal2(1)
       n22 = normal2(2)
       n23 = normal2(3)
       x2 = dm%x(1, l, k, i)
       y2 = dm%y(1, l, k, i)
       z2 = dm%z(1, l, k, i)
    case (3)
       normal1 = coef%get_normal(1, l, k, i, 3)
       n11 = normal1(1)
       n12 = normal1(2)
       n13 = normal1(3)
       x1 = dm%x(l, 1, k, i)
       y1 = dm%y(l, 1, k, i)
       z1 = dm%z(l, 1, k, i)
       normal2 = coef%get_normal(1, l, k, i, 4)
       n21 = normal2(1)
       n22 = normal2(2)
       n23 = normal2(3)
       x2 = dm%x(l, n, k, i)
       y2 = dm%y(l, n, k, i)
       z2 = dm%z(l, n, k, i)
    case (4)
       normal1 = coef%get_normal(1, l, k, i, 4)
       n11 = normal1(1)
       n12 = normal1(2)
       n13 = normal1(3)
       x1 = dm%x(l, n, k, i)
       y1 = dm%y(l, n, k, i)
       z1 = dm%z(l, n, k, i)
       normal2 = coef%get_normal(1, l, k, i, 3)
       n21 = normal2(1)
       n22 = normal2(2)
       n23 = normal2(3)
       x2 = dm%x(l, 1, k, i)
       y2 = dm%y(l, 1, k, i)
       z2 = dm%z(l, 1, k, i)
    case (5)
       normal1 = coef%get_normal(1, l, k, i, 5)
       n11 = normal1(1)
       n12 = normal1(2)
       n13 = normal1(3)
       x1 = dm%x(l, k, 1, i)
       y1 = dm%y(l, k, 1, i)
       z1 = dm%z(l, k, 1, i)
       normal2 = coef%get_normal(1, l, k, i, 6)
       n21 = normal2(1)
       n22 = normal2(2)
       n23 = normal2(3)
       x2 = dm%x(l, k, n, i)
       y2 = dm%y(l, k, n, i)
       z2 = dm%z(l, k, n, i)
    case (6)
       normal1 = coef%get_normal(1, l, k, i, 6)
       n11 = normal1(1)
       n12 = normal1(2)
       n13 = normal1(3)
       x1 = dm%x(l, k, n, i)
       y1 = dm%y(l, k, n, i)
       z1 = dm%z(l, k, n, i)
       normal2 = coef%get_normal(1, l, k, i, 5)
       n21 = normal2(1)
       n22 = normal2(2)
       n23 = normal2(3)
       x2 = dm%x(l, k, 1, i)
       y2 = dm%y(l, k, 1, i)
       z2 = dm%z(l, k, 1, i)
    case default
       call neko_error("The face index is not correct")
    end select

    ! get the vector from the quadrature point to the one on the other side
    v1 = x2 - x1
    v2 = y2 - y1
    v3 = z2 - z1
    ! Project onto tabsvolflhe facet-normal direction of the point
    dist_1 = v1*n11 + v2*n12 + v3*n13
    dist_2 = - (v1*n21 + v2*n22 + v3*n23)

    dist2 = ((dist_1 + dist_2)/2.0_rp)*((dist_1 + dist_2)/2.0_rp)

  end function dist2_quadrature_hex

  !> Initialize the facet factor array
  subroutine facet_factor_init(this)
    class(gradient_jump_penalty_t), intent(inout) :: this
    integer :: i, j, k, l
    real(kind=rp) :: area_tmp

    allocate(this%facet_factor(this%lx + 2, this%lx + 2, &
         this%lx + 2, this%coef%msh%nelv))

    associate(facet_factor => this%facet_factor, &
         coef => this%coef, &
         lx => this%lx, &
         nelv => this%coef%msh%nelv, &
         jacinv => this%coef%jacinv, h2 => this%h2, &
         tau => this%tau, n1 => this%n1, &
         n2 => this%n2, n3 => this%n3)

      do i = 1, nelv
         do j = 1, 6 ! for hexahedral elementsh2
            do k = 1, lx
               do l = 1, lx
                  select case (j)
                  case (1)
                     area_tmp = coef%get_area(1, l, k, i, j)
                     facet_factor(1, l + 1, k + 1, i) = area_tmp * tau * &
                          h2(1, l + 1, k + 1, i) * &
                          (n1(1, l + 1, k + 1, i) * coef%drdx(1, l, k, i) + &
                          n2(1, l + 1, k + 1, i) * coef%drdy(1, l, k, i) + &
                          n3(1, l + 1, k + 1, i) * coef%drdz(1, l, k, i) ) &
                          * jacinv(1, l, k, i)
                  case (2)
                     area_tmp = coef%get_area(1, l, k, i, j)
                     facet_factor(lx + 2, l + 1, k + 1, i) = area_tmp * tau * &
                          h2(lx + 2, l + 1, k + 1, i) * &
                          (n1(lx + 2, l + 1, k + 1, i) * &
                          coef%drdx(lx, l, k, i) + &
                          n2(lx + 2, l + 1, k + 1, i) * &
                          coef%drdy(lx, l, k, i) + &
                          n3(lx + 2, l + 1, k + 1, i) * &
                          coef%drdz(lx, l, k, i) ) &
                          * jacinv(lx, l, k, i)
                  case (3)
                     area_tmp = coef%get_area(l, 1, k, i, j)
                     facet_factor(l + 1, 1, k + 1, i) = area_tmp * tau * &
                          h2(l + 1, 1, k + 1, i) * &
                          (n1(l + 1, 1, k + 1, i) * coef%dsdx(l, 1, k, i) + &
                          n2(l + 1, 1, k + 1, i) * coef%dsdy(l, 1, k, i) + &
                          n3(l + 1, 1, k + 1, i) * coef%dsdz(l, 1, k, i) ) &
                          * jacinv(l, 1, k, i)
                  case (4)
                     area_tmp = coef%get_area(l, 1, k, i, j)
                     facet_factor(l + 1, lx + 2, k + 1, i) = area_tmp * tau * &
                          h2(l + 1, lx + 2, k + 1, i) * &
                          (n1(l + 1, lx + 2, k + 1, i) * &
                          coef%dsdx(l, lx, k, i) + &
                          n2(l + 1, lx + 2, k + 1, i) * &
                          coef%dsdy(l, lx, k, i) + &
                          n3(l + 1, lx + 2, k + 1, i) * &
                          coef%dsdz(l, lx, k, i) ) &
                          * jacinv(l, lx, k, i)
                  case (5)
                     area_tmp = coef%get_area(l, k, 1, i, j)
                     facet_factor(l + 1, k + 1, 1, i) = area_tmp * tau * &
                          h2(l + 1, k + 1, 1, i) * &
                          (n1(l + 1, k + 1, 1, i) * coef%dtdx(l, k, 1, i) + &
                          n2(l + 1, k + 1, 1, i) * coef%dtdy(l, k, 1, i) + &
                          n3(l + 1, k + 1, 1, i) * coef%dtdz(l, k, 1, i) ) &
                          * jacinv(l, k, 1, i)
                  case (6)
                     area_tmp = coef%get_area(l, k, 1, i, j)
                     facet_factor(l + 1, k + 1, lx + 2, i) = area_tmp * tau * &
                          h2(l + 1, k + 1, lx + 2, i) * &
                          ( &
                          n1(l + 1, k + 1, lx + 2, i) * &
                          coef%dtdx(l, k, lx, i) + &
                          n2(l + 1, k + 1, lx + 2, i) * &
                          coef%dtdy(l, k, lx, i) + &
                          n3(l + 1, k + 1, lx + 2, i) * &
                          coef%dtdz(l, k, lx, i) &
                          ) &
                          * jacinv(l, k, lx, i)
                  case default
                     call neko_error("The face index is not correct")
                  end select
               end do
            end do
         end do
      end do

    end associate
  end subroutine facet_factor_init

  !> Destructor for the gradient_jump_penalty_t class.
  subroutine gradient_jump_penalty_free(this)
    implicit none
    class(gradient_jump_penalty_t), intent(inout) :: this

    call this%free_base

    if (c_associated(this%dphidxi_d)) then
       call device_free(this%dphidxi_d)
    end if
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

    nullify(this%u)
    nullify(this%v)
    nullify(this%w)

    call this%s_fields%free()

    call this%Xh_GJP%free()
    call this%gs_GJP%free()
    call this%dm_GJP%free()

  end subroutine gradient_jump_penalty_free

  !> Compute the gradient jump penalty term for a single field.
  !! @param s Quantity of interest.
  subroutine gradient_jump_penalty_compute_single(this, s)
    class(gradient_jump_penalty_t), intent(inout) :: this
    type(field_t), intent(in) :: s

    class(element_t), pointer :: ep
    integer :: i

    call G_compute(this, s)
    call absvolflux_compute(this, this%u, this%v, this%w)

    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_col3(this%penalty_facet_d, this%absvolflux_d, this%G_d, &
            this%n_large)
       call device_col2(this%penalty_facet_d, this%facet_factor_d, this%n_large)
       call device_gradient_jump_penalty_finalize(this%penalty_d, &
            this%penalty_facet_d, &
            this%dphidxi_d, &
            this%lx, this%coef%msh%nelv)
    else
       call col3(this%penalty_facet, this%absvolflux, this%G, this%n_large)
       call col2(this%penalty_facet, this%facet_factor, this%n_large)
       call gradient_jump_penalty_finalize(this%penalty, this%penalty_facet, &
            this%dphidxi, &
            this%lx, this%coef%msh%nelv)
    end if

  end subroutine gradient_jump_penalty_compute_single

  !> Assign the gradient jump penalty term.
  !! @param time The time state.
  subroutine gradient_jump_penalty_compute(this, time)
    class(gradient_jump_penalty_t), intent(inout) :: this
    type(time_state_t), intent(in) :: time
    integer :: i, n_fields, n

    n_fields = this%fields%size()
    n = this%coef%dof%size()



    do i = 1, n_fields

       call this%compute_single(this%s_fields%items(i)%ptr)

       if (NEKO_BCKND_DEVICE .eq. 1) then
          call device_invcol2(this%penalty_d, this%coef%B_d, n)
          call device_add2(this%fields%x_d(i), this%penalty_d, n)
       else
          call invcol2(this%penalty, this%coef%B, n)
          call add2(this%fields%items(i)%ptr%x, this%penalty, n)
       end if
    end do

  end subroutine gradient_jump_penalty_compute

  !> Interface of finalizing the gradient jump penalty term.
  !> <tau * h^2 * absvolflux * G * phij * phik * dphi_idxi * dxidn>
  !! @param penalty Gradient Jump Penalty array
  !! @param wa Work array containing tau * absvolflux * G * dxidn
  !! @param dphidxi The first derivative of polynomial
  !! @param lx Order of polynomial plus one
  !! @param nelv Number of elements
  subroutine gradient_jump_penalty_finalize(penalty, wa, dphidxi, lx, nelv)
    integer, intent(in) :: lx, nelv
    real(kind=rp), intent(inout) :: penalty(lx, lx, lx, nelv)
    real(kind=rp), intent(in) :: wa(lx + 2, lx + 2, lx + 2, nelv)
    real(kind=rp), intent(in) :: dphidxi(lx, lx)

    call gradient_jump_penalty_finalize_hex(penalty, wa, dphidxi, lx, nelv)

  end subroutine gradient_jump_penalty_finalize

  !> Finalizinge the gradient jump penalty term for hexahedral elements.
  !> <tau * h^2 * absvolflux * G * phij * phik * dphi_idxi * dxidn>
  !! @param penalty Gradient Jump Penalty array
  !! @param wa Work array containing tau * absvolflux * G * dxidn
  !! @param dphidxi The first derivative of polynomial
  !! @param lx Order of polynomial plus one
  !! @param nelv Number of elements
  subroutine gradient_jump_penalty_finalize_hex(penalty, wa, dphidxi, lx, nelv)
    integer, intent(in) :: lx, nelv
    real(kind=rp), intent(inout) :: penalty(lx, lx, lx, nelv)
    real(kind=rp), intent(in) :: wa(lx + 2, lx + 2, lx + 2, nelv)
    real(kind=rp), intent(in) :: dphidxi(lx, lx)

    integer :: i, j, k

    do i = 1, lx
       do j = 1, lx
          do k = 1, lx
             penalty(i, j, k, :) = &
                  wa(1, j + 1, k + 1, :) * &
                  dphidxi(1, i) + &
                  wa(lx + 2, j + 1, k + 1, :) * &
                  dphidxi(lx, i) + &
                  wa(i + 1, 1, k + 1, :) * &
                  dphidxi(1, j) + &
                  wa(i + 1, lx + 2, k + 1, :) * &
                  dphidxi(lx, j) + &
                  wa(i + 1, j + 1, 1, :) * &
                  dphidxi(1, k) + &
                  wa(i + 1, j + 1, lx + 2, :) * &
                  dphidxi(lx, k)
          end do
       end do
    end do

  end subroutine gradient_jump_penalty_finalize_hex

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
       call device_pick_facet_value_hex(this%flux1_d, this%grad1_d, &
            this%lx, this%coef%msh%nelv)
       call device_pick_facet_value_hex(this%flux2_d, this%grad2_d, &
            this%lx, this%coef%msh%nelv)
       call device_pick_facet_value_hex(this%flux3_d, this%grad3_d, &
            this%lx, this%coef%msh%nelv)
       call device_col2(this%flux1_d, this%n1_d, this%n_large)
       call device_col2(this%flux2_d, this%n2_d, this%n_large)
       call device_col2(this%flux3_d, this%n3_d, this%n_large)
       call device_add3s2(this%G_d, this%flux1_d, this%flux2_d, &
            1.0_rp, 1.0_rp, this%n_large)
       call device_add2(this%G_d, this%flux3_d, this%n_large)
    else
       call pick_facet_value_hex(this%flux1, this%grad1, &
            this%lx, this%coef%msh%nelv)
       call pick_facet_value_hex(this%flux2, this%grad2, &
            this%lx, this%coef%msh%nelv)
       call pick_facet_value_hex(this%flux3, this%grad3, &
            this%lx, this%coef%msh%nelv)
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

    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_pick_facet_value_hex(this%volflux1_d, u%x_d, this%lx, &
            this%coef%msh%nelv)
       call device_pick_facet_value_hex(this%volflux2_d, v%x_d, this%lx, &
            this%coef%msh%nelv)
       call device_pick_facet_value_hex(this%volflux3_d, w%x_d, this%lx, &
            this%coef%msh%nelv)
       call device_col2(this%volflux1_d, this%n1_d, this%n_large)
       call device_col2(this%volflux2_d, this%n2_d, this%n_large)
       call device_col2(this%volflux3_d, this%n3_d, this%n_large)
       call device_add3s2(this%absvolflux_d, this%volflux1_d, &
            this%volflux2_d, 1.0_rp, 1.0_rp, this%n_large)
       call device_add2(this%absvolflux_d, this%volflux3_d, this%n_large)
       call device_absval(this%absvolflux_d, this%n_large)
    else
       call pick_facet_value_hex(this%volflux1, u%x, &
            this%lx, this%coef%msh%nelv)
       call pick_facet_value_hex(this%volflux2, v%x, &
            this%lx, this%coef%msh%nelv)
       call pick_facet_value_hex(this%volflux3, w%x, &
            this%lx, this%coef%msh%nelv)
       call col2(this%volflux1, this%n1, this%n_large)
       call col2(this%volflux2, this%n2, this%n_large)
       call col2(this%volflux3, this%n3, this%n_large)
       call add3(this%absvolflux, this%volflux1, this%volflux2, this%n_large)
       call add2(this%absvolflux, this%volflux3, this%n_large)
       call absval(this%absvolflux, this%n_large)
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
