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
  use mpi_f08, only : MPI_Allreduce, MPI_LAND, MPI_LOGICAL
  use num_types, only : rp
  use field, only : field_t
  use coefs, only : coef_t
  use gather_scatter, only : gs_t
  use gs_ops, only : GS_OP_ADD
  use comm, only : NEKO_COMM
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
     real(kind=rp), allocatable :: reconstruction_work_1(:,:,:,:)
     type(field_t) :: mass
     type(field_t) :: directional_degree(3)
   contains
     procedure, pass(this) :: init => euler_gll_graph_init
     procedure, pass(this) :: free => euler_gll_graph_free
     procedure, pass(this) :: local_flux_divergence => &
          euler_gll_graph_local_flux_divergence
     procedure, pass(this) :: incidence => euler_gll_graph_incidence
     procedure, pass(this) :: directional_correction_flux => &
          euler_gll_graph_directional_correction_flux
     procedure, pass(this) :: periodic_facets_covered => &
          euler_gll_graph_periodic_facets_covered
  end type euler_gll_graph_t

contains

  !> Build the sparse tensor-product cumulative-SBP graph.
  subroutine euler_gll_graph_init(this, coef, gs)
    class(euler_gll_graph_t), intent(inout) :: this
    type(coef_t), target, intent(in) :: coef
    type(gs_t), intent(inout) :: gs
    real(kind=rp) :: metric(3), transverse_weight
    character(len=LOG_SIZE) :: message
    character(len=48) :: name
    integer :: direction, edge, e, i, j, k

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
    allocate(this%reconstruction_work_1(this%lx, this%ly, this%lz, &
         this%nelv))
    this%reconstruction_work_1 = 0.0_rp

    edge = 0
    do e = 1, this%nelv
       this%element_edge_start(e) = edge + 1

       metric = [coef%drdx(1,1,1,e), &
            coef%drdy(1,1,1,e), coef%drdz(1,1,1,e)]
       do k = 1, this%lz
          do j = 1, this%ly
             transverse_weight = 0.5_rp * coef%Xh%wy(j) * coef%Xh%wz(k)
             do i = 1, this%lx - 1
                edge = edge + 1
                call euler_gll_graph_set_edge(this, edge, i, j, k, &
                     i + 1, j, k, e, 1, transverse_weight * metric)
             end do
             this%diagonal_coefficient(:,1,j,k,e) = &
                  this%diagonal_coefficient(:,1,j,k,e) - &
                  transverse_weight * metric
             this%diagonal_coefficient(:,this%lx,j,k,e) = &
                  this%diagonal_coefficient(:,this%lx,j,k,e) + &
                  transverse_weight * metric
          end do
       end do

       metric = [coef%dsdx(1,1,1,e), &
            coef%dsdy(1,1,1,e), coef%dsdz(1,1,1,e)]
       do k = 1, this%lz
          do i = 1, this%lx
             transverse_weight = 0.5_rp * coef%Xh%wx(i) * coef%Xh%wz(k)
             do j = 1, this%ly - 1
                edge = edge + 1
                call euler_gll_graph_set_edge(this, edge, i, j, k, &
                     i, j + 1, k, e, 2, transverse_weight * metric)
             end do
             this%diagonal_coefficient(:,i,1,k,e) = &
                  this%diagonal_coefficient(:,i,1,k,e) - &
                  transverse_weight * metric
             this%diagonal_coefficient(:,i,this%ly,k,e) = &
                  this%diagonal_coefficient(:,i,this%ly,k,e) + &
                  transverse_weight * metric
          end do
       end do

       if (this%lz .gt. 1) then
          metric = [coef%dtdx(1,1,1,e), &
               coef%dtdy(1,1,1,e), coef%dtdz(1,1,1,e)]
          do j = 1, this%ly
             do i = 1, this%lx
                transverse_weight = 0.5_rp * coef%Xh%wx(i) * &
                     coef%Xh%wy(j)
                do k = 1, this%lz - 1
                   edge = edge + 1
                   call euler_gll_graph_set_edge(this, edge, i, j, k, &
                        i, j, k + 1, e, 3, transverse_weight * metric)
                end do
                this%diagonal_coefficient(:,i,j,1,e) = &
                     this%diagonal_coefficient(:,i,j,1,e) - &
                     transverse_weight * metric
                this%diagonal_coefficient(:,i,j,this%lz,e) = &
                     this%diagonal_coefficient(:,i,j,this%lz,e) + &
                     transverse_weight * metric
             end do
          end do
       end if
       this%element_edge_start(e + 1) = edge + 1
    end do

    if (edge .ne. this%n_edges) then
       call neko_error('Incorrect Euler GLL graph edge count')
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

    if (.not. this%initialized) return

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
    if (allocated(this%reconstruction_work_1)) then
       deallocate(this%reconstruction_work_1)
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

  !> Apply and gather-scatter the sparse derivative.

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

  !> Apply the transpose of the element-local incidence matrix.

  !> Reconstruct an element-local correction with weighted tensor sweeps.
  !! The input must have zero unweighted sum on every element. Any roundoff
  !! defect is left as the tensor product of the normalised GLL weights and is
  !! therefore visible in the caller's divergence-residual diagnostic.

  !> Reconstruct the inviscid high-to-low correction without mixing tensor
  !! directions. Each reference-direction residual is inverted only on the
  !! coordinate edges that produced it. On affine elements every coordinate
  !! line is independently compatible, so no global graph solve is required.
  subroutine euler_gll_graph_directional_correction_flux(this, edge_values, &
       flux_x, flux_y, flux_z, coef, compatibility_error)
    class(euler_gll_graph_t), intent(inout) :: this
    real(kind=rp), intent(out) :: edge_values(:)
    real(kind=rp), intent(in) :: flux_x(:,:,:,:), flux_y(:,:,:,:)
    real(kind=rp), intent(in) :: flux_z(:,:,:,:)
    type(coef_t), intent(in) :: coef
    real(kind=rp), intent(out) :: compatibility_error
    real(kind=rp) :: cumulative, derivative_x, derivative_y, derivative_z
    real(kind=rp) :: high_residual, line_scale
    integer :: edge, edge_first, edge_last, e, i, j, k, l

    call euler_gll_graph_assert_initialized(this)
    if (size(edge_values) .ne. this%n_edges) then
       call neko_error('Incorrect Euler directional correction edge size')
    end if

    edge_values = 0.0_rp
    compatibility_error = 0.0_rp

    ! Reference-x contribution.
    this%reconstruction_work_1 = 0.0_rp
    do edge = 1, this%n_edges
       if (this%direction(edge) .ne. 1) cycle
       associate(a => this%left(:,edge), b => this%right(:,edge), &
            c => this%coefficient(:,edge))
         this%reconstruction_work_1(a(1),a(2),a(3),a(4)) = &
              this%reconstruction_work_1(a(1),a(2),a(3),a(4)) + &
              c(1) * flux_x(b(1),b(2),b(3),b(4)) + &
              c(2) * flux_y(b(1),b(2),b(3),b(4)) + &
              c(3) * flux_z(b(1),b(2),b(3),b(4))
         this%reconstruction_work_1(b(1),b(2),b(3),b(4)) = &
              this%reconstruction_work_1(b(1),b(2),b(3),b(4)) - &
              c(1) * flux_x(a(1),a(2),a(3),a(4)) - &
              c(2) * flux_y(a(1),a(2),a(3),a(4)) - &
              c(3) * flux_z(a(1),a(2),a(3),a(4))
       end associate
    end do
    do e = 1, this%nelv
       edge = this%element_edge_start(e) - 1
       do k = 1, this%lz
          do j = 1, this%ly
             edge_first = edge + 1
             edge_last = edge + this%lx - 1
             associate(c_first => this%coefficient(:,edge_first), &
                  c_last => this%coefficient(:,edge_last))
               this%reconstruction_work_1(1,j,k,e) = &
                    this%reconstruction_work_1(1,j,k,e) - &
                    c_first(1) * flux_x(1,j,k,e) - &
                    c_first(2) * flux_y(1,j,k,e) - &
                    c_first(3) * flux_z(1,j,k,e)
               this%reconstruction_work_1(this%lx,j,k,e) = &
                    this%reconstruction_work_1(this%lx,j,k,e) + &
                    c_last(1) * flux_x(this%lx,j,k,e) + &
                    c_last(2) * flux_y(this%lx,j,k,e) + &
                    c_last(3) * flux_z(this%lx,j,k,e)
             end associate
             line_scale = 0.0_rp
             do i = 1, this%lx
                derivative_x = 0.0_rp
                derivative_y = 0.0_rp
                derivative_z = 0.0_rp
                do l = 1, this%lx
                   derivative_x = derivative_x + coef%Xh%dx(i,l) * &
                        flux_x(l,j,k,e)
                   derivative_y = derivative_y + coef%Xh%dx(i,l) * &
                        flux_y(l,j,k,e)
                   derivative_z = derivative_z + coef%Xh%dx(i,l) * &
                        flux_z(l,j,k,e)
                end do
                high_residual = coef%B(i,j,k,e) * &
                     coef%jacinv(i,j,k,e) * ( &
                     coef%drdx(i,j,k,e) * derivative_x + &
                     coef%drdy(i,j,k,e) * derivative_y + &
                     coef%drdz(i,j,k,e) * derivative_z)
                line_scale = line_scale + &
                     abs(this%reconstruction_work_1(i,j,k,e)) + &
                     abs(high_residual)
                this%reconstruction_work_1(i,j,k,e) = &
                     this%reconstruction_work_1(i,j,k,e) - high_residual
             end do
             line_scale = max(1.0_rp, line_scale)
             compatibility_error = max(compatibility_error, &
                  abs(sum(this%reconstruction_work_1(:,j,k,e))) / line_scale)
             cumulative = 0.0_rp
             do i = 1, this%lx - 1
                edge = edge + 1
                cumulative = cumulative + &
                     this%reconstruction_work_1(i,j,k,e)
                edge_values(edge) = cumulative
             end do
          end do
       end do
    end do

    ! Reference-y contribution.
    this%reconstruction_work_1 = 0.0_rp
    do edge = 1, this%n_edges
       if (this%direction(edge) .ne. 2) cycle
       associate(a => this%left(:,edge), b => this%right(:,edge), &
            c => this%coefficient(:,edge))
         this%reconstruction_work_1(a(1),a(2),a(3),a(4)) = &
              this%reconstruction_work_1(a(1),a(2),a(3),a(4)) + &
              c(1) * flux_x(b(1),b(2),b(3),b(4)) + &
              c(2) * flux_y(b(1),b(2),b(3),b(4)) + &
              c(3) * flux_z(b(1),b(2),b(3),b(4))
         this%reconstruction_work_1(b(1),b(2),b(3),b(4)) = &
              this%reconstruction_work_1(b(1),b(2),b(3),b(4)) - &
              c(1) * flux_x(a(1),a(2),a(3),a(4)) - &
              c(2) * flux_y(a(1),a(2),a(3),a(4)) - &
              c(3) * flux_z(a(1),a(2),a(3),a(4))
       end associate
    end do
    do e = 1, this%nelv
       edge = this%element_edge_start(e) - 1 + &
            (this%lx - 1) * this%ly * this%lz
       do k = 1, this%lz
          do i = 1, this%lx
             edge_first = edge + 1
             edge_last = edge + this%ly - 1
             associate(c_first => this%coefficient(:,edge_first), &
                  c_last => this%coefficient(:,edge_last))
               this%reconstruction_work_1(i,1,k,e) = &
                    this%reconstruction_work_1(i,1,k,e) - &
                    c_first(1) * flux_x(i,1,k,e) - &
                    c_first(2) * flux_y(i,1,k,e) - &
                    c_first(3) * flux_z(i,1,k,e)
               this%reconstruction_work_1(i,this%ly,k,e) = &
                    this%reconstruction_work_1(i,this%ly,k,e) + &
                    c_last(1) * flux_x(i,this%ly,k,e) + &
                    c_last(2) * flux_y(i,this%ly,k,e) + &
                    c_last(3) * flux_z(i,this%ly,k,e)
             end associate
             line_scale = 0.0_rp
             do j = 1, this%ly
                derivative_x = 0.0_rp
                derivative_y = 0.0_rp
                derivative_z = 0.0_rp
                do l = 1, this%ly
                   derivative_x = derivative_x + coef%Xh%dy(j,l) * &
                        flux_x(i,l,k,e)
                   derivative_y = derivative_y + coef%Xh%dy(j,l) * &
                        flux_y(i,l,k,e)
                   derivative_z = derivative_z + coef%Xh%dy(j,l) * &
                        flux_z(i,l,k,e)
                end do
                high_residual = coef%B(i,j,k,e) * &
                     coef%jacinv(i,j,k,e) * ( &
                     coef%dsdx(i,j,k,e) * derivative_x + &
                     coef%dsdy(i,j,k,e) * derivative_y + &
                     coef%dsdz(i,j,k,e) * derivative_z)
                line_scale = line_scale + &
                     abs(this%reconstruction_work_1(i,j,k,e)) + &
                     abs(high_residual)
                this%reconstruction_work_1(i,j,k,e) = &
                     this%reconstruction_work_1(i,j,k,e) - high_residual
             end do
             line_scale = max(1.0_rp, line_scale)
             compatibility_error = max(compatibility_error, &
                  abs(sum(this%reconstruction_work_1(i,:,k,e))) / line_scale)
             cumulative = 0.0_rp
             do j = 1, this%ly - 1
                edge = edge + 1
                cumulative = cumulative + &
                     this%reconstruction_work_1(i,j,k,e)
                edge_values(edge) = cumulative
             end do
          end do
       end do
    end do

    if (this%lz .gt. 1) then
       ! Reference-z contribution.
       this%reconstruction_work_1 = 0.0_rp
       do edge = 1, this%n_edges
          if (this%direction(edge) .ne. 3) cycle
          associate(a => this%left(:,edge), b => this%right(:,edge), &
               c => this%coefficient(:,edge))
            this%reconstruction_work_1(a(1),a(2),a(3),a(4)) = &
                 this%reconstruction_work_1(a(1),a(2),a(3),a(4)) + &
                 c(1) * flux_x(b(1),b(2),b(3),b(4)) + &
                 c(2) * flux_y(b(1),b(2),b(3),b(4)) + &
                 c(3) * flux_z(b(1),b(2),b(3),b(4))
            this%reconstruction_work_1(b(1),b(2),b(3),b(4)) = &
                 this%reconstruction_work_1(b(1),b(2),b(3),b(4)) - &
                 c(1) * flux_x(a(1),a(2),a(3),a(4)) - &
                 c(2) * flux_y(a(1),a(2),a(3),a(4)) - &
                 c(3) * flux_z(a(1),a(2),a(3),a(4))
          end associate
       end do
       do e = 1, this%nelv
          edge = this%element_edge_start(e) - 1 + &
               (this%lx - 1) * this%ly * this%lz + &
               this%lx * (this%ly - 1) * this%lz
          do j = 1, this%ly
             do i = 1, this%lx
                edge_first = edge + 1
                edge_last = edge + this%lz - 1
                associate(c_first => this%coefficient(:,edge_first), &
                     c_last => this%coefficient(:,edge_last))
                  this%reconstruction_work_1(i,j,1,e) = &
                       this%reconstruction_work_1(i,j,1,e) - &
                       c_first(1) * flux_x(i,j,1,e) - &
                       c_first(2) * flux_y(i,j,1,e) - &
                       c_first(3) * flux_z(i,j,1,e)
                  this%reconstruction_work_1(i,j,this%lz,e) = &
                       this%reconstruction_work_1(i,j,this%lz,e) + &
                       c_last(1) * flux_x(i,j,this%lz,e) + &
                       c_last(2) * flux_y(i,j,this%lz,e) + &
                       c_last(3) * flux_z(i,j,this%lz,e)
                end associate
                line_scale = 0.0_rp
                do k = 1, this%lz
                   derivative_x = 0.0_rp
                   derivative_y = 0.0_rp
                   derivative_z = 0.0_rp
                   do l = 1, this%lz
                      derivative_x = derivative_x + coef%Xh%dz(k,l) * &
                           flux_x(i,j,l,e)
                      derivative_y = derivative_y + coef%Xh%dz(k,l) * &
                           flux_y(i,j,l,e)
                      derivative_z = derivative_z + coef%Xh%dz(k,l) * &
                           flux_z(i,j,l,e)
                   end do
                   high_residual = coef%B(i,j,k,e) * &
                        coef%jacinv(i,j,k,e) * ( &
                        coef%dtdx(i,j,k,e) * derivative_x + &
                        coef%dtdy(i,j,k,e) * derivative_y + &
                        coef%dtdz(i,j,k,e) * derivative_z)
                   line_scale = line_scale + &
                        abs(this%reconstruction_work_1(i,j,k,e)) + &
                        abs(high_residual)
                   this%reconstruction_work_1(i,j,k,e) = &
                        this%reconstruction_work_1(i,j,k,e) - high_residual
                end do
                line_scale = max(1.0_rp, line_scale)
                compatibility_error = max(compatibility_error, &
                     abs(sum(this%reconstruction_work_1(i,j,:,e))) / &
                     line_scale)
                cumulative = 0.0_rp
                do k = 1, this%lz - 1
                   edge = edge + 1
                   cumulative = cumulative + &
                        this%reconstruction_work_1(i,j,k,e)
                   edge_values(edge) = cumulative
                end do
             end do
          end do
       end do
    end if
  end subroutine euler_gll_graph_directional_correction_flux

  !> Apply one element-local coordinate-neighbour Laplacian.

  !> Apply and assemble the coordinate-neighbour graph Laplacian.

  !> Return assembled row and column sums of the graph coefficient matrix.

  !> Sum one replicated nodal array with each assembled node counted once.

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
