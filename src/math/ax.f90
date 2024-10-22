! Copyright (c) 2020-2021, The Neko Authors
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
!> Defines a Matrix-vector product
module ax_product
  use num_types, only : rp
  use coefs, only : coef_t
  use space, only : space_t
  use mesh, only : mesh_t
  implicit none
  private

  !> Base type for a matrix-vector product providing \f$ Ax \f$
  type, public, abstract :: ax_t
   contains
     procedure(ax_compute), nopass, deferred :: compute
     procedure(ax_compute_vector), pass(this), deferred :: compute_vector
  end type ax_t

  interface
     !> Factory routine for the a Helmholtz problem matrix-vector product.
     !! The selection is based on the compute backend.
     !! @param object The matrix-vector product type to be allocated.
     !! @param full_formulation Whether to use the formulation with the full
     !! viscous stress tensor, not assuming constant material properties.
     module subroutine ax_helm_factory(object, full_formulation)
       class(ax_t), allocatable, intent(inout) :: object
       logical, intent(in) :: full_formulation
     end subroutine ax_helm_factory
  end interface

  public :: ax_helm_factory
  
  !> Abstract interface for computing\f$ Ax \f$ inside a Krylov method
  !!
  !! @param w Vector of size @a (lx,ly,lz,nelv).
  !! @param u Vector of size @a (lx,ly,lz,nelv).
  !! @param coef Coefficients.
  !! @param msh Mesh.
  !! @param Xh Function space \f$ X_h \f$.
  abstract interface
     subroutine ax_compute(w, u, coef, msh, Xh)
       import space_t
       import mesh_t
       import coef_t
       import ax_t
       import rp
       implicit none
       type(space_t), intent(inout) :: Xh
       type(mesh_t), intent(inout) :: msh
       type(coef_t), intent(inout) :: coef
       real(kind=rp), intent(inout) :: w(Xh%lx, Xh%ly, Xh%lz, msh%nelv)
       real(kind=rp), intent(inout) :: u(Xh%lx, Xh%ly, Xh%lz, msh%nelv)
     end subroutine ax_compute
  end interface

  !> Abstract interface for computing\f$ Ax \f$ inside a Krylov method,
  !! taking 3 components of a vector field in a coupled manner.
  !! @param au Result for the first component of the vector.
  !! @param av Result for the first component of the vector.
  !! @param aw Result for the first component of the vector.
  !! @param u The first component of the vector.
  !! @param v The second component of the vector.
  !! @param w The third component of the vector.
  !! @param coef Coefficients.
  !! @param msh Mesh.
  !! @param Xh Function space \f$ X_h \f$.
  abstract interface
     subroutine ax_compute_vector(this, au, av, aw, u, v, w, coef, msh, Xh)
       import space_t
       import mesh_t
       import coef_t
       import ax_t
       import rp
       implicit none
       class(ax_t), intent(in) :: this
       type(space_t), intent(inout) :: Xh
       type(mesh_t), intent(inout) :: msh
       type(coef_t), intent(inout) :: coef
       real(kind=rp), intent(inout) :: au(Xh%lx, Xh%ly, Xh%lz, msh%nelv)
       real(kind=rp), intent(inout) :: av(Xh%lx, Xh%ly, Xh%lz, msh%nelv)
       real(kind=rp), intent(inout) :: aw(Xh%lx, Xh%ly, Xh%lz, msh%nelv)
       real(kind=rp), intent(inout) :: u(Xh%lx, Xh%ly, Xh%lz, msh%nelv)
       real(kind=rp), intent(inout) :: v(Xh%lx, Xh%ly, Xh%lz, msh%nelv)
       real(kind=rp), intent(inout) :: w(Xh%lx, Xh%ly, Xh%lz, msh%nelv)
     end subroutine ax_compute_vector
  end interface

end module ax_product
