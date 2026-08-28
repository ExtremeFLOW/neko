! Copyright (c) 2026, The Neko Authors
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
!     copyright notice, this list of conditions and the following disclaimer
!     in the documentation and/or other materials provided with the
!     distribution.
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
!> Sparse cumulative-SBP graph on the collocated Euler GLL nodes.
module euler_gll_graph
  use mpi_f08, only : MPI_Allreduce, MPI_LAND, MPI_LOGICAL, MPI_MAX
  use num_types, only : rp
  use field, only : field_t
  use coefs, only : coef_t
  use gather_scatter, only : gs_t
  use gs_ops, only : GS_OP_ADD
  use comm, only : NEKO_COMM, MPI_REAL_PRECISION
  use utils, only : neko_error
  use logger, only : LOG_SIZE
  implicit none
  private

  !> Element-local coordinate-neighbour multigraph.
  type, public :: euler_gll_graph_t
     logical :: initialized = .false.
     integer :: lx = 0
     integer :: ly = 0
     integer :: lz = 0
     integer :: nelv = 0
     integer :: n_directions = 0
     integer :: edges_per_element = 0
     integer :: n_edges = 0
     integer, allocatable :: left(:,:)
     integer, allocatable :: right(:,:)
     integer, allocatable :: direction(:)
     integer, allocatable :: element_edge_start(:)
     real(kind=rp), allocatable :: coefficient(:,:)
     real(kind=rp), allocatable :: diagonal_coefficient(:,:,:,:,:)
     !> One-dimensional GLL quadrature weights used by reconstruction.
     real(kind=rp), allocatable :: weight_x(:)
     real(kind=rp), allocatable :: weight_y(:)
     real(kind=rp), allocatable :: weight_z(:)
     !> Element-local work arrays used by correction reconstruction.
     real(kind=rp), allocatable :: reconstruction_work_1(:,:,:,:)
     real(kind=rp), allocatable :: reconstruction_work_2(:,:,:,:)
     type(field_t) :: mass
     type(field_t) :: directional_degree(3)
   contains
     procedure, pass(this) :: init => euler_gll_graph_init
     procedure, pass(this) :: free => euler_gll_graph_free
     procedure, pass(this) :: local_flux_divergence => &
          euler_gll_graph_local_flux_divergence
     procedure, pass(this) :: high_order_flux_divergence => &
          euler_gll_graph_high_order_flux_divergence
     procedure, pass(this) :: incidence => euler_gll_graph_incidence
     procedure, pass(this) :: directional_correction_flux => &
          euler_gll_graph_directional_correction_flux
     procedure, pass(this) :: periodic_facets_covered => &
          euler_gll_graph_periodic_facets_covered
  end type euler_gll_graph_t

contains

  !> Build the sparse tensor-product cumulative-SBP graph.
  !! @param this Euler GLL graph to initialise.
  !! @param coef SEM coefficients defining the geometry and quadrature.
  !! @param gs Gather-scatter operator for assembled nodal quantities.
  subroutine euler_gll_graph_init(this, coef, gs)
    class(euler_gll_graph_t), intent(inout) :: this
    type(coef_t), target, intent(in) :: coef
    type(gs_t), intent(inout) :: gs
    real(kind=rp) :: metric_left(3), metric_right(3)
    real(kind=rp) :: transverse_weight
    real(kind=rp) :: local_error, global_error, local_scale, global_scale
    real(kind=rp), allocatable :: coefficient_correction(:)
    character(len=LOG_SIZE) :: message
    character(len=48) :: name
    integer :: direction, edge, e, i, ierr, j, k

    call this%free()

    this%lx = coef%Xh%lx
    this%ly = coef%Xh%ly
    this%lz = coef%Xh%lz
    this%nelv = coef%msh%nelv
    this%n_directions = 2
    if (this%lz .gt. 1) this%n_directions = 3
    this%edges_per_element = (this%lx - 1) * this%ly * this%lz + &
         this%lx * (this%ly - 1) * this%lz
    if (this%lz .gt. 1) then
       this%edges_per_element = this%edges_per_element + &
            this%lx * this%ly * (this%lz - 1)
    end if
    this%n_edges = this%edges_per_element * this%nelv

    allocate(this%left(4, this%n_edges))
    allocate(this%right(4, this%n_edges))
    allocate(this%direction(this%n_edges))
    allocate(this%coefficient(3, this%n_edges))
    allocate(this%element_edge_start(this%nelv + 1))
    this%element_edge_start = 1
    allocate(this%diagonal_coefficient(3, this%lx, this%ly, this%lz, &
         this%nelv))
    this%diagonal_coefficient = 0.0_rp
    allocate(this%weight_x(this%lx))
    allocate(this%weight_y(this%ly))
    allocate(this%weight_z(this%lz))
    this%weight_x = coef%Xh%wx
    this%weight_y = coef%Xh%wy
    this%weight_z = coef%Xh%wz
    allocate(this%reconstruction_work_1(this%lx, this%ly, this%lz, &
         this%nelv))
    this%reconstruction_work_1 = 0.0_rp
    allocate(this%reconstruction_work_2(this%lx, this%ly, this%lz, &
         this%nelv))
    this%reconstruction_work_2 = 0.0_rp

    edge = 0
    do e = 1, this%nelv
       this%element_edge_start(e) = edge + 1

       do k = 1, this%lz
          do j = 1, this%ly
             transverse_weight = 0.5_rp * coef%Xh%wy(j) * coef%Xh%wz(k)
             do i = 1, this%lx - 1
                edge = edge + 1
                metric_left = [coef%drdx(i,j,k,e), &
                     coef%drdy(i,j,k,e), coef%drdz(i,j,k,e)]
                metric_right = [coef%drdx(i + 1,j,k,e), &
                     coef%drdy(i + 1,j,k,e), &
                     coef%drdz(i + 1,j,k,e)]
                call euler_gll_graph_set_edge(this, edge, i, j, k, &
                     i + 1, j, k, e, 1, &
                     0.5_rp * transverse_weight * &
                     (metric_left + metric_right))
             end do
             metric_left = [coef%drdx(1,j,k,e), &
                  coef%drdy(1,j,k,e), coef%drdz(1,j,k,e)]
             metric_right = [coef%drdx(this%lx,j,k,e), &
                  coef%drdy(this%lx,j,k,e), &
                  coef%drdz(this%lx,j,k,e)]
             this%diagonal_coefficient(:,1,j,k,e) = &
                  this%diagonal_coefficient(:,1,j,k,e) - &
                  transverse_weight * metric_left
             this%diagonal_coefficient(:,this%lx,j,k,e) = &
                  this%diagonal_coefficient(:,this%lx,j,k,e) + &
                  transverse_weight * metric_right
          end do
       end do

       do k = 1, this%lz
          do i = 1, this%lx
             transverse_weight = 0.5_rp * coef%Xh%wx(i) * coef%Xh%wz(k)
             do j = 1, this%ly - 1
                edge = edge + 1
                metric_left = [coef%dsdx(i,j,k,e), &
                     coef%dsdy(i,j,k,e), coef%dsdz(i,j,k,e)]
                metric_right = [coef%dsdx(i,j + 1,k,e), &
                     coef%dsdy(i,j + 1,k,e), &
                     coef%dsdz(i,j + 1,k,e)]
                call euler_gll_graph_set_edge(this, edge, i, j, k, &
                     i, j + 1, k, e, 2, &
                     0.5_rp * transverse_weight * &
                     (metric_left + metric_right))
             end do
             metric_left = [coef%dsdx(i,1,k,e), &
                  coef%dsdy(i,1,k,e), coef%dsdz(i,1,k,e)]
             metric_right = [coef%dsdx(i,this%ly,k,e), &
                  coef%dsdy(i,this%ly,k,e), &
                  coef%dsdz(i,this%ly,k,e)]
             this%diagonal_coefficient(:,i,1,k,e) = &
                  this%diagonal_coefficient(:,i,1,k,e) - &
                  transverse_weight * metric_left
             this%diagonal_coefficient(:,i,this%ly,k,e) = &
                  this%diagonal_coefficient(:,i,this%ly,k,e) + &
                  transverse_weight * metric_right
          end do
       end do

       if (this%lz .gt. 1) then
          do j = 1, this%ly
             do i = 1, this%lx
                transverse_weight = 0.5_rp * coef%Xh%wx(i) * &
                     coef%Xh%wy(j)
                do k = 1, this%lz - 1
                   edge = edge + 1
                   metric_left = [coef%dtdx(i,j,k,e), &
                        coef%dtdy(i,j,k,e), coef%dtdz(i,j,k,e)]
                   metric_right = [coef%dtdx(i,j,k + 1,e), &
                        coef%dtdy(i,j,k + 1,e), &
                        coef%dtdz(i,j,k + 1,e)]
                   call euler_gll_graph_set_edge(this, edge, i, j, k, &
                        i, j, k + 1, e, 3, &
                        0.5_rp * transverse_weight * &
                        (metric_left + metric_right))
                end do
                metric_left = [coef%dtdx(i,j,1,e), &
                     coef%dtdy(i,j,1,e), coef%dtdz(i,j,1,e)]
                metric_right = [coef%dtdx(i,j,this%lz,e), &
                     coef%dtdy(i,j,this%lz,e), &
                     coef%dtdz(i,j,this%lz,e)]
                this%diagonal_coefficient(:,i,j,1,e) = &
                     this%diagonal_coefficient(:,i,j,1,e) - &
                     transverse_weight * metric_left
                this%diagonal_coefficient(:,i,j,this%lz,e) = &
                     this%diagonal_coefficient(:,i,j,this%lz,e) + &
                     transverse_weight * metric_right
             end do
          end do
       end if
       this%element_edge_start(e + 1) = edge + 1
    end do

    if (edge .ne. this%n_edges) then
       call neko_error('Incorrect Euler GLL graph edge count')
    end if

    allocate(coefficient_correction(this%n_edges))
    do direction = 1, 3
       this%reconstruction_work_1 = &
            this%diagonal_coefficient(direction,:,:,:,:)
       do edge = 1, this%n_edges
          associate(a => this%left(:,edge), b => this%right(:,edge))
            this%reconstruction_work_1(a(1),a(2),a(3),a(4)) = &
                 this%reconstruction_work_1(a(1),a(2),a(3),a(4)) + &
                 this%coefficient(direction,edge)
            this%reconstruction_work_1(b(1),b(2),b(3),b(4)) = &
                 this%reconstruction_work_1(b(1),b(2),b(3),b(4)) - &
                 this%coefficient(direction,edge)
          end associate
       end do
       this%reconstruction_work_1 = -this%reconstruction_work_1
       call euler_gll_graph_reconstruct(this, coefficient_correction, &
            this%reconstruction_work_1, local_error, local_scale)
       this%coefficient(direction,:) = &
            this%coefficient(direction,:) + coefficient_correction
    end do
    deallocate(coefficient_correction)

    local_error = 0.0_rp
    local_scale = 0.0_rp
    do direction = 1, 3
       this%reconstruction_work_1 = &
            this%diagonal_coefficient(direction,:,:,:,:)
       do edge = 1, this%n_edges
          associate(a => this%left(:,edge), b => this%right(:,edge))
            this%reconstruction_work_1(a(1),a(2),a(3),a(4)) = &
                 this%reconstruction_work_1(a(1),a(2),a(3),a(4)) + &
                 this%coefficient(direction,edge)
            this%reconstruction_work_1(b(1),b(2),b(3),b(4)) = &
                 this%reconstruction_work_1(b(1),b(2),b(3),b(4)) - &
                 this%coefficient(direction,edge)
          end associate
       end do
       local_error = max(local_error, &
            maxval(abs(this%reconstruction_work_1)))
       local_scale = max(local_scale, &
            maxval(abs(this%diagonal_coefficient(direction,:,:,:,:))), &
            maxval(abs(this%coefficient(direction,:))))
    end do
    call MPI_Allreduce(local_error, global_error, 1, MPI_REAL_PRECISION, &
         MPI_MAX, NEKO_COMM, ierr)
    call MPI_Allreduce(local_scale, global_scale, 1, MPI_REAL_PRECISION, &
         MPI_MAX, NEKO_COMM, ierr)
    if (global_error .gt. 1024.0_rp * epsilon(1.0_rp) * &
         max(1.0_rp, global_scale)) then
       call neko_error('Euler GLL graph metric projection failed')
    end if

    call this%mass%init(coef%dof, 'euler_gll_graph_mass')
    do direction = 1, 3
       write(name, '(A,I0)') 'euler_gll_graph_directional_degree_', &
            direction
       call this%directional_degree(direction)%init(coef%dof, trim(name))
    end do
    this%mass%x = coef%B
    call gs%op(this%mass, GS_OP_ADD)
    if (any(this%mass%x .le. 0.0_rp)) then
       write(message, '(A,ES13.6)') &
            'Euler GLL graph requires positive assembled masses, minimum: ', &
            minval(this%mass%x)
       call neko_error(trim(message))
    end if

    do direction = 1, 3
       this%directional_degree(direction)%x = 0.0_rp
    end do
    do edge = 1, this%n_edges
       associate(a => this%left(:,edge), b => this%right(:,edge))
         direction = this%direction(edge)
         this%directional_degree(direction)%x(a(1),a(2),a(3),a(4)) = &
              this%directional_degree(direction)%x(a(1),a(2),a(3),a(4)) + &
              1.0_rp
         this%directional_degree(direction)%x(b(1),b(2),b(3),b(4)) = &
              this%directional_degree(direction)%x(b(1),b(2),b(3),b(4)) + &
              1.0_rp
       end associate
    end do
    do direction = 1, 3
       call gs%op(this%directional_degree(direction), GS_OP_ADD)
    end do
    this%initialized = .true.
  end subroutine euler_gll_graph_init

  !> Free graph storage.
  subroutine euler_gll_graph_free(this)
    class(euler_gll_graph_t), intent(inout) :: this
    integer :: direction

    if (allocated(this%left)) deallocate(this%left)
    if (allocated(this%right)) deallocate(this%right)
    if (allocated(this%direction)) deallocate(this%direction)
    if (allocated(this%coefficient)) deallocate(this%coefficient)
    if (allocated(this%element_edge_start)) then
       deallocate(this%element_edge_start)
    end if
    if (allocated(this%diagonal_coefficient)) then
       deallocate(this%diagonal_coefficient)
    end if
    if (allocated(this%weight_x)) deallocate(this%weight_x)
    if (allocated(this%weight_y)) deallocate(this%weight_y)
    if (allocated(this%weight_z)) deallocate(this%weight_z)
    if (allocated(this%reconstruction_work_1)) then
       deallocate(this%reconstruction_work_1)
    end if
    if (allocated(this%reconstruction_work_2)) then
       deallocate(this%reconstruction_work_2)
    end if
    call this%mass%free()
    do direction = 1, 3
       call this%directional_degree(direction)%free()
    end do
    this%initialized = .false.
    this%lx = 0
    this%ly = 0
    this%lz = 0
    this%nelv = 0
    this%n_directions = 0
    this%edges_per_element = 0
    this%n_edges = 0
  end subroutine euler_gll_graph_free

  !> Apply the element-local sparse derivative to a Cartesian flux.
  subroutine euler_gll_graph_local_flux_divergence(this, output, flux_x, &
       flux_y, flux_z)
    class(euler_gll_graph_t), intent(in) :: this
    real(kind=rp), intent(out) :: output(:,:,:,:)
    real(kind=rp), intent(in) :: flux_x(:,:,:,:), flux_y(:,:,:,:)
    real(kind=rp), intent(in) :: flux_z(:,:,:,:)
    real(kind=rp) :: flux_left, flux_right
    integer :: edge, e, i, j, k

    call euler_gll_graph_assert_initialized(this)
    output = 0.0_rp
    do e = 1, this%nelv
       do k = 1, this%lz
          do j = 1, this%ly
             do i = 1, this%lx
                output(i,j,k,e) = &
                     this%diagonal_coefficient(1,i,j,k,e) * flux_x(i,j,k,e) + &
                     this%diagonal_coefficient(2,i,j,k,e) * flux_y(i,j,k,e) + &
                     this%diagonal_coefficient(3,i,j,k,e) * flux_z(i,j,k,e)
             end do
          end do
       end do
    end do

    do edge = 1, this%n_edges
       associate(a => this%left(:,edge), b => this%right(:,edge), &
            c => this%coefficient(:,edge))
         flux_left = c(1) * flux_x(a(1),a(2),a(3),a(4)) + &
              c(2) * flux_y(a(1),a(2),a(3),a(4)) + &
              c(3) * flux_z(a(1),a(2),a(3),a(4))
         flux_right = c(1) * flux_x(b(1),b(2),b(3),b(4)) + &
              c(2) * flux_y(b(1),b(2),b(3),b(4)) + &
              c(3) * flux_z(b(1),b(2),b(3),b(4))
         output(a(1),a(2),a(3),a(4)) = &
              output(a(1),a(2),a(3),a(4)) + flux_right
         output(b(1),b(2),b(3),b(4)) = &
              output(b(1),b(2),b(3),b(4)) - flux_left
       end associate
    end do
  end subroutine euler_gll_graph_local_flux_divergence

  !> Apply the mass-weighted conservative SEM divergence.
  !! @param output Mass-weighted divergence.
  !! @param flux_x Cartesian x-component of the flux.
  !! @param flux_y Cartesian y-component of the flux.
  !! @param flux_z Cartesian z-component of the flux.
  !! @param coef SEM coefficients.
  subroutine euler_gll_graph_high_order_flux_divergence(this, output, flux_x, &
       flux_y, flux_z, coef)
    class(euler_gll_graph_t), intent(in) :: this
    real(kind=rp), intent(out) :: output(:,:,:,:)
    real(kind=rp), intent(in) :: flux_x(:,:,:,:), flux_y(:,:,:,:)
    real(kind=rp), intent(in) :: flux_z(:,:,:,:)
    type(coef_t), intent(in) :: coef
    real(kind=rp) :: contravariant_flux, quadrature_weight
    integer :: e, i, j, k, l

    call euler_gll_graph_assert_initialized(this)
    output = 0.0_rp
    do e = 1, this%nelv
       do k = 1, this%lz
          do j = 1, this%ly
             do i = 1, this%lx
                quadrature_weight = coef%Xh%wx(i) * coef%Xh%wy(j) * &
                     coef%Xh%wz(k)
                do l = 1, this%lx
                   contravariant_flux = &
                        coef%drdx(l,j,k,e) * flux_x(l,j,k,e) + &
                        coef%drdy(l,j,k,e) * flux_y(l,j,k,e) + &
                        coef%drdz(l,j,k,e) * flux_z(l,j,k,e)
                   output(i,j,k,e) = output(i,j,k,e) + &
                        quadrature_weight * coef%Xh%dx(i,l) * &
                        contravariant_flux
                end do
                do l = 1, this%ly
                   contravariant_flux = &
                        coef%dsdx(i,l,k,e) * flux_x(i,l,k,e) + &
                        coef%dsdy(i,l,k,e) * flux_y(i,l,k,e) + &
                        coef%dsdz(i,l,k,e) * flux_z(i,l,k,e)
                   output(i,j,k,e) = output(i,j,k,e) + &
                        quadrature_weight * coef%Xh%dy(j,l) * &
                        contravariant_flux
                end do
                if (this%lz .gt. 1) then
                   do l = 1, this%lz
                      contravariant_flux = &
                           coef%dtdx(i,j,l,e) * flux_x(i,j,l,e) + &
                           coef%dtdy(i,j,l,e) * flux_y(i,j,l,e) + &
                           coef%dtdz(i,j,l,e) * flux_z(i,j,l,e)
                      output(i,j,k,e) = output(i,j,k,e) + &
                           quadrature_weight * coef%Xh%dz(k,l) * &
                           contravariant_flux
                   end do
                end if
             end do
          end do
       end do
    end do
  end subroutine euler_gll_graph_high_order_flux_divergence

  !> Apply the element-local incidence matrix to edge values.
  subroutine euler_gll_graph_incidence(this, output, edge_values)
    class(euler_gll_graph_t), intent(in) :: this
    real(kind=rp), intent(out) :: output(:,:,:,:)
    real(kind=rp), intent(in) :: edge_values(:)
    integer :: edge

    call euler_gll_graph_assert_initialized(this)
    if (size(edge_values) .ne. this%n_edges) then
       call neko_error('Incorrect Euler GLL graph edge-vector size')
    end if
    output = 0.0_rp
    do edge = 1, this%n_edges
       associate(a => this%left(:,edge), b => this%right(:,edge))
         output(a(1),a(2),a(3),a(4)) = &
              output(a(1),a(2),a(3),a(4)) + edge_values(edge)
         output(b(1),b(2),b(3),b(4)) = &
              output(b(1),b(2),b(3),b(4)) - edge_values(edge)
       end associate
    end do
  end subroutine euler_gll_graph_incidence

  !> Reconstruct an element-local correction with weighted tensor sweeps.
  !! The input must have zero unweighted sum on every element. Any roundoff
  !! defect is left as the tensor product of the normalised GLL weights and is
  !! therefore visible in the caller's divergence-residual diagnostic.
  !! @param edge_values Reconstructed oriented edge values.
  !! @param residual Element-local nodal residual to reconstruct.
  !! @param compatibility_error Maximum relative element-sum defect.
  !! @param residual_scale Maximum element residual scale.
  subroutine euler_gll_graph_reconstruct(this, edge_values, residual, &
       compatibility_error, residual_scale)
    class(euler_gll_graph_t), intent(in) :: this
    real(kind=rp), intent(out) :: edge_values(:)
    real(kind=rp), intent(in) :: residual(:,:,:,:)
    real(kind=rp), intent(out) :: compatibility_error, residual_scale
    real(kind=rp) :: cumulative, element_scale, element_sum
    real(kind=rp) :: line_sum, plane_sum, value
    real(kind=rp) :: weight_x, weight_y, weight_z
    real(kind=rp) :: weight_sum_x, weight_sum_y, weight_sum_z
    integer :: edge, e, i, j, k

    if (size(edge_values) .ne. this%n_edges) then
       call neko_error('Incorrect Euler GLL reconstruction edge size')
    end if

    edge_values = 0.0_rp
    compatibility_error = 0.0_rp
    residual_scale = 0.0_rp
    weight_sum_x = sum(this%weight_x)
    weight_sum_y = sum(this%weight_y)
    weight_sum_z = sum(this%weight_z)

    do e = 1, this%nelv
       element_sum = sum(residual(:,:,:,e))
       element_scale = sum(abs(residual(:,:,:,e)))
       residual_scale = max(residual_scale, element_scale)
       compatibility_error = max(compatibility_error, &
            abs(element_sum) / max(1.0_rp, element_scale))

       edge = this%element_edge_start(e) - 1
       do k = 1, this%lz
          do j = 1, this%ly
             line_sum = sum(residual(:,j,k,e))
             cumulative = 0.0_rp
             do i = 1, this%lx - 1
                weight_x = this%weight_x(i) / weight_sum_x
                value = residual(i,j,k,e) - weight_x * line_sum
                cumulative = cumulative + value
                edge = edge + 1
                edge_values(edge) = cumulative
             end do
          end do
       end do

       do k = 1, this%lz
          plane_sum = sum(residual(:,:,k,e))
          do i = 1, this%lx
             weight_x = this%weight_x(i) / weight_sum_x
             cumulative = 0.0_rp
             do j = 1, this%ly - 1
                weight_y = this%weight_y(j) / weight_sum_y
                line_sum = sum(residual(:,j,k,e))
                value = weight_x * (line_sum - weight_y * plane_sum)
                cumulative = cumulative + value
                edge = edge + 1
                edge_values(edge) = cumulative
             end do
          end do
       end do

       if (this%lz .gt. 1) then
          do j = 1, this%ly
             weight_y = this%weight_y(j) / weight_sum_y
             do i = 1, this%lx
                weight_x = this%weight_x(i) / weight_sum_x
                cumulative = 0.0_rp
                do k = 1, this%lz - 1
                   weight_z = this%weight_z(k) / weight_sum_z
                   plane_sum = sum(residual(:,:,k,e))
                   value = weight_x * weight_y * &
                        (plane_sum - weight_z * element_sum)
                   cumulative = cumulative + value
                   edge = edge + 1
                   edge_values(edge) = cumulative
                end do
             end do
          end do
       end if
    end do
  end subroutine euler_gll_graph_reconstruct

  !> Reconstruct the conservative inviscid high-to-low correction.
  !! @param edge_values Reconstructed oriented edge correction.
  !! @param flux_x Cartesian x-component of the flux.
  !! @param flux_y Cartesian y-component of the flux.
  !! @param flux_z Cartesian z-component of the flux.
  !! @param coef SEM coefficients.
  !! @param compatibility_error Maximum relative element-sum defect.
  subroutine euler_gll_graph_directional_correction_flux(this, edge_values, &
       flux_x, flux_y, flux_z, coef, compatibility_error)
    class(euler_gll_graph_t), intent(inout) :: this
    real(kind=rp), intent(out) :: edge_values(:)
    real(kind=rp), intent(in) :: flux_x(:,:,:,:), flux_y(:,:,:,:)
    real(kind=rp), intent(in) :: flux_z(:,:,:,:)
    type(coef_t), intent(in) :: coef
    real(kind=rp), intent(out) :: compatibility_error
    real(kind=rp) :: residual_scale

    call euler_gll_graph_assert_initialized(this)
    if (size(edge_values) .ne. this%n_edges) then
       call neko_error('Incorrect Euler directional correction edge size')
    end if

    call this%local_flux_divergence(this%reconstruction_work_1, flux_x, &
         flux_y, flux_z)
    call this%high_order_flux_divergence(this%reconstruction_work_2, flux_x, &
         flux_y, flux_z, coef)
    this%reconstruction_work_1 = this%reconstruction_work_1 - &
         this%reconstruction_work_2
    call euler_gll_graph_reconstruct(this, edge_values, &
         this%reconstruction_work_1, compatibility_error, residual_scale)
  end subroutine euler_gll_graph_directional_correction_flux

  !> Check that every physical exterior facet has periodic metadata.
  logical function euler_gll_graph_periodic_facets_covered(this, coef) &
       result(covered)
    class(euler_gll_graph_t), intent(in) :: this
    type(coef_t), intent(in) :: coef
    logical :: facet_covered, local_covered
    integer :: e, f, i, ierr, n_facets

    call euler_gll_graph_assert_initialized(this)
    n_facets = merge(4, 6, coef%msh%gdim .eq. 2)
    local_covered = .true.
    do e = 1, coef%msh%nelv
       do f = 1, n_facets
          if (coef%msh%facet_neigh(f,e) .eq. 0) then
             facet_covered = .false.
             do i = 1, coef%msh%periodic%size
                if (coef%msh%periodic%facet_el(i)%x(1) .eq. f .and. &
                     coef%msh%periodic%facet_el(i)%x(2) .eq. e) then
                   facet_covered = .true.
                   exit
                end if
             end do
             if (.not. facet_covered) local_covered = .false.
          end if
       end do
    end do
    call MPI_Allreduce(local_covered, covered, 1, MPI_LOGICAL, MPI_LAND, &
         NEKO_COMM, ierr)
  end function euler_gll_graph_periodic_facets_covered

  !> Count graph components by distributed label propagation.
  !! This verification operation is not used in the timestep path.

  !> Store one oriented edge occurrence.
  subroutine euler_gll_graph_set_edge(this, edge, i1, j1, k1, i2, j2, k2, &
       element, direction, coefficient)
    class(euler_gll_graph_t), intent(inout) :: this
    integer, intent(in) :: edge, i1, j1, k1, i2, j2, k2, element, direction
    real(kind=rp), intent(in) :: coefficient(3)

    this%left(:,edge) = [i1, j1, k1, element]
    this%right(:,edge) = [i2, j2, k2, element]
    this%direction(edge) = direction
    this%coefficient(:,edge) = coefficient
  end subroutine euler_gll_graph_set_edge

  !> Require a constructed graph before an operator application.
  subroutine euler_gll_graph_assert_initialized(this)
    class(euler_gll_graph_t), intent(in) :: this

    if (.not. this%initialized) then
       call neko_error('Euler GLL graph is not initialised')
    end if
  end subroutine euler_gll_graph_assert_initialized

end module euler_gll_graph
