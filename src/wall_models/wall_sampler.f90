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
     !> Wall-normal distance of every sampling point. The layout is
     !! `(node - 1) * n_samples + sample`.
     type(vector_t) :: h
   contains
     !> Parse sampler-specific configuration.
     procedure(wall_sampler_init), pass(this), deferred :: init
     !> Final initialization step, mirroring the wall_model_t interface.
     procedure(wall_sampler_finalize), pass(this), deferred :: finalize
     !> Sample the solution field at the sampling points.
     procedure(wall_sampler_sample), pass(this), deferred :: sample
     !> Destructor.
     procedure(wall_sampler_free), pass(this), deferred :: free
  end type wall_sampler_t

  abstract interface
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
     subroutine wall_sampler_finalize(this, coef, msk, facet, n_x, n_y, n_z)
       import wall_sampler_t, coef_t, vector_t
       class(wall_sampler_t), intent(inout) :: this
       type(coef_t), intent(in) :: coef
       integer, intent(in) :: msk(0:)
       integer, intent(in) :: facet(0:)
       type(vector_t), intent(in) :: n_x, n_y, n_z
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

     !> Destructor
     subroutine wall_sampler_free(this)
       import wall_sampler_t
       class(wall_sampler_t), intent(inout) :: this
     end subroutine wall_sampler_free
  end interface

end module wall_sampler
