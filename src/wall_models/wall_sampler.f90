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
!> Defines the abstract interface for wall-model field samplers.
module wall_sampler
  use field, only : field_t
  use coefs, only : coef_t
  use vector, only : vector_t
  use json_module, only : json_file
  use user_intf, only : user_t
  use utils, only : neko_error
  implicit none
  private

  !> Base type for sampling solution fields at points associated with wall
  !! nodes. Samples belonging to one wall node are stored contiguously.
  !! @details The design supports multiple sampling points per wall node.
  type, abstract, public :: wall_sampler_t
     !> Number of wall nodes local to this rank.
     integer :: n_nodes = 0
     !> Number of sampling points associated with each wall node.
     integer :: n_samples = 0
     !> True when sampling values are supplied by a user callback.
     logical :: user_values = .false.
     !> Wall-normal distance of every sampling point. The layout is
     !! `(node - 1) * n_samples + sample`.
     type(vector_t) :: h
   contains
     !> Initialise state common to all fully constructed wall samplers.
     procedure, pass(this) :: init_base => wall_sampler_init_base
     !> Parse sampler-specific configuration from JSON.
     procedure(wall_sampler_init), pass(this), deferred :: init
     !> Complete sampler setup after geometric and wall-node data are known.
     procedure(wall_sampler_finalize), pass(this), deferred :: finalize
     !> Sample the solution field at the sampling points.
     procedure(wall_sampler_sample), pass(this), deferred :: sample
     !> Release sampler resources.
     procedure(wall_sampler_free), pass(this), deferred :: free
  end type wall_sampler_t

  abstract interface
     !> Parse sampler-specific configuration from JSON.
     !! @param json Sampler configuration data.
     subroutine wall_sampler_init(this, json)
       import wall_sampler_t, json_file
       class(wall_sampler_t), intent(inout) :: this
       type(json_file), intent(inout) :: json
     end subroutine wall_sampler_init

     !> Complete sampler setup after geometric and wall-node data are known.
     !! @param coef SEM coefficients.
     !! @param msk Mask selecting local wall nodes.
     !! @param facet Facet index for each selected wall node.
     !! @param n_x X-component of wall-normal vectors at wall nodes.
     !! @param n_y Y-component of wall-normal vectors at wall nodes.
     !! @param n_z Z-component of wall-normal vectors at wall nodes.
     !! @param bc_name Name of the owning boundary condition when user values
     !! are requested.
     !! @param user User interface providing sampling callbacks when user
     !! values are requested.
     subroutine wall_sampler_finalize(this, coef, msk, facet, n_x, n_y, n_z, &
          bc_name, user)
       import wall_sampler_t, coef_t, vector_t, user_t
       class(wall_sampler_t), intent(inout) :: this
       type(coef_t), intent(in) :: coef
       integer, intent(in) :: msk(0:)
       integer, intent(in) :: facet(0:)
       type(vector_t), intent(in) :: n_x, n_y, n_z
       character(len=*), optional, intent(in) :: bc_name
       type(user_t), target, optional, intent(in) :: user
     end subroutine wall_sampler_finalize

     !> Evaluate and store sampled field values at all sampling points.
     !! @param field Field data to be sampled.
     !! @param values Output sampled values in sampler ordering.
     subroutine wall_sampler_sample(this, field, values)
       import wall_sampler_t, field_t, vector_t
       class(wall_sampler_t), intent(inout) :: this
       type(field_t), intent(inout) :: field
       type(vector_t), intent(inout) :: values
     end subroutine wall_sampler_sample

     !> Release sampler resources.
     subroutine wall_sampler_free(this)
       import wall_sampler_t
       class(wall_sampler_t), intent(inout) :: this
     end subroutine wall_sampler_free
  end interface

contains

  !> Initialise state common to all fully constructed wall samplers.
  !! @param n_nodes Number of local wall nodes.
  !! @param n_samples Number of samples at each wall node.
  !! @param h Wall-normal distances in sampler ordering.
  subroutine wall_sampler_init_base(this, n_nodes, n_samples, h)
    class(wall_sampler_t), intent(inout) :: this
    integer, intent(in) :: n_nodes
    integer, intent(in) :: n_samples
    type(vector_t), intent(in) :: h

    if (n_nodes < 1 .or. n_samples < 1) then
       call neko_error('Wall sampler dimensions must be positive')
    end if
    if (h%size() /= n_nodes * n_samples) then
       call neko_error('Wall sampler distances have an invalid size')
    end if

    call this%h%free()
    this%n_nodes = n_nodes
    this%n_samples = n_samples
    this%user_values = .false.
    this%h = h
  end subroutine wall_sampler_init_base

end module wall_sampler
