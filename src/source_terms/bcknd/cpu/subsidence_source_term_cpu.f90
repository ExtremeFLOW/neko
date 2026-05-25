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
!> Implements the cpu kernel for the  `subsidence_source_term_t` type.
! Copyright (c) 2026, The Neko Authors
! All rights reserved.
!
!> Implements the CPU kernel for the `subsidence_source_term_t` type.
module subsidence_source_term_cpu
  use num_types, only : rp
  use field_list, only : field_list_t
  use field, only : field_t
  use scratch_registry, only : neko_scratch_registry
  use coefs, only : coef_t
  use operators, only : grad
  use gs_ops, only : GS_OP_ADD
  use math, only : col2
  implicit none
  private

  public :: subsidence_source_term_compute_cpu

contains

  !> Computes and aggregates subsidence force terms onto the velocity RHS fields.
  !! The routine computes directional derivatives of `u` and `v` along the
  !! provided `vertical_dir` and adds w_sub * d(phi)/d(vertical_dir) to the
  !! corresponding RHS fields in `fields`.
  subroutine subsidence_source_term_compute_cpu(u, v, fields, coef, w_sub, vertical_dir)
    type(field_t), intent(in) :: u, v, w_sub
    type(field_list_t), intent(inout) :: fields
    type(coef_t), intent(in) :: coef
    real(kind=rp), intent(in) :: vertical_dir(3)

    integer :: n, i
    integer :: temp_indices(3)
    type(field_t), pointer :: fu, fv
    type(field_t), pointer :: dx, dy, dz
    real(kind=rp) :: deriv

    n = fields%item_size(1)

    fu => fields%get_by_index(1)
    fv => fields%get_by_index(2)

    ! Request three scratch fields for gradients
    call neko_scratch_registry%request_field(dx, temp_indices(1), .false.)
    call neko_scratch_registry%request_field(dy, temp_indices(2), .false.)
    call neko_scratch_registry%request_field(dz, temp_indices(3), .false.)

    ! Compute gradients of u
    call grad(dx%x, dy%x, dz%x, u%x, coef)

    call coef%gs_h%op(dx, GS_OP_ADD)
    call coef%gs_h%op(dy, GS_OP_ADD)
    call coef%gs_h%op(dz, GS_OP_ADD)

    call col2(dx%x, coef%mult, u%dof%size())
    call col2(dy%x, coef%mult, u%dof%size())
    call col2(dz%x, coef%mult, u%dof%size())

    ! Add subsidence contribution to U RHS: fu += w_sub * (d/dz_dir u)
    do concurrent (i = 1:u%dof%size())
       deriv = vertical_dir(1) * dx%x(i,1,1,1) + &
               vertical_dir(2) * dy%x(i,1,1,1) + &
               vertical_dir(3) * dz%x(i,1,1,1)
       fu%x(i,1,1,1) = fu%x(i,1,1,1) + w_sub%x(i,1,1,1) * deriv
    end do

    ! Compute gradients of v (reuse scratch fields)
    call grad(dx%x, dy%x, dz%x, v%x, coef)

    call coef%gs_h%op(dx, GS_OP_ADD)
    call coef%gs_h%op(dy, GS_OP_ADD)
    call coef%gs_h%op(dz, GS_OP_ADD)

    call col2(dx%x, coef%mult, v%dof%size())
    call col2(dy%x, coef%mult, v%dof%size())
    call col2(dz%x, coef%mult, v%dof%size())

    do concurrent (i = 1:v%dof%size())
       deriv = vertical_dir(1) * dx%x(i,1,1,1) + &
               vertical_dir(2) * dy%x(i,1,1,1) + &
               vertical_dir(3) * dz%x(i,1,1,1)
       fv%x(i,1,1,1) = fv%x(i,1,1,1) + w_sub%x(i,1,1,1) * deriv
    end do

    call neko_scratch_registry%relinquish_field(temp_indices)

  end subroutine subsidence_source_term_compute_cpu

end module subsidence_source_term_cpu
