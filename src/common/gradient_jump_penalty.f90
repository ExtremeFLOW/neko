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
  use field, only : field_t
  use dofmap , only : dofmap_t
  use neko_config, only : NEKO_BCKND_DEVICE
  use coefs, only : coef_t
  use scratch_registry, only : neko_scratch_registry
  use element, only : element_t
  use hex
  use quad

  implicit none
  private

  !> Implements the gradient jump penalty.
  !! @note Reference DOI: 10.1016/j.cma.2021.114200
  type, public :: gradient_jump_penalty_t
     !> Coefficient of the penalty.
     real :: tau
     !> Polynomial order
     integer :: p
     !> Penalty terms
     real(kind=rp), allocatable :: penalty(:, :, :, :)
     !> SEM coefficients.
     type(coef_t), pointer :: coef => null()
     !> Gradient jump on the elementary interface (zero inside each element)
     type(field_t), pointer :: G
     !> Gradient of the quantity
     type(field_t), pointer :: grad1, grad2, grad3
     !> Number of facet in elements and its maximum
     integer, allocatable :: n_facet(:)
     integer :: n_facet_max
     !> Length scale for element regarding a facet
     real(kind=rp), allocatable :: h(:)
     !> Polynomial Quadrature
     real(kind=rp), allocatable :: w(:), w2(:,:)
     !> Polynomial evaluated at collocation points
     real(kind=rp), allocatable :: phi(:,:)
     !> The first derivative of polynomial at two ends of the interval
     real(kind=rp), allocatable :: dphidxi(:,:)

  contains
     !> Constructor.
     procedure, pass(this) :: init => gradient_jump_penalty_init
     !> Destructor
     procedure, pass(this) :: free => gradient_jump_penalty_free
    !  !> Compute gradient jump penalty term.
    !  procedure, pass(this) :: compute => gradient_jump_penalty_compute

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
    integer :: temp_indices(4), i, j
    real(kind=rp), allocatable :: zg(:) ! Quadrature points
    
    call this%free()

    this%p = dofmap%xh%lx - 1
    if (this%p .gt. 1) then
       this%tau = 0.8_rp * (this%p + 1) ** (-4)
    else
       this%tau = 0.02_rp
    end if
    this%coef => coef

    allocate(this%penalty(this%p + 1, this%p + 1 , this%p + 1 , this%coef%msh%nelv))
    
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

    allocate(this%h(this%coef%msh%nelv))
    do i = 1, this%coef%msh%nelv
       ep => this%coef%msh%elements(i)%e
       select type(ep)
       type is (hex_t)
          call eval_h_hex(this%h(i), ep)
       type is (quad_t)
          call neko_error("Gradient jump penalty error: mesh size evaluation is not supported for quad_t")
       end select
    end do
    
    allocate(zg(dofmap%xh%lx))
    allocate(this%w(dofmap%xh%lx))
    allocate(this%w2(dofmap%xh%lx, dofmap%xh%lx))
    allocate(this%phi(this%p + 1, this%p + 1))
    allocate(this%dphidxi(this%p + 1, this%p + 1))

    zg = dofmap%xh%zg(:,1)
    this%w = dofmap%xh%wx
    do i = 1, dofmap%xh%lx
       do j = 1, dofmap%xh%lx
          this%w2(i, j) = this%w(i) * this%w(j)
       end do
    end do

    do i = 1, dofmap%xh%lx
       do j = 1, dofmap%xh%lx
          this%phi(j,i) = pnleg(zg(j),i-1)
          this%dphidxi(j,i) = pndleg(zg(j),i-1)
       end do
    end do
    write(*,*) this%dphidxi(:,3)


    call neko_scratch_registry%request_field(this%G, temp_indices(1))
    call neko_scratch_registry%request_field(this%grad1, temp_indices(2))
    call neko_scratch_registry%request_field(this%grad2, temp_indices(3))
    call neko_scratch_registry%request_field(this%grad3, temp_indices(4))

  end subroutine gradient_jump_penalty_init
  
  !> Evaluate h for each element
  subroutine eval_h_hex(h_el, ep)
  real(kind=rp), intent(inout) :: h_el
  type(hex_t), pointer, intent(in) :: ep
  
  !! todo: estimation of the length scale of the mesh could be more elegant
  h_el = ep%diameter()

  end subroutine eval_h_hex

  !> Destructor for the gradient_jump_penalty_t class.
  subroutine gradient_jump_penalty_free(this)
    implicit none
    class(gradient_jump_penalty_t), intent(inout) :: this
    
    if (allocated(this%penalty)) then
       deallocate(this%penalty)
    end if
    if (allocated(this%h)) then
       deallocate(this%h)
    end if
    if (allocated(this%n_facet)) then
       deallocate(this%n_facet)
    end if
    if (allocated(this%w)) then
       deallocate(this%w)
    end if
    if (allocated(this%w2)) then
       deallocate(this%w2)
    end if
    if (allocated(this%phi)) then
       deallocate(this%phi)
    end if
    if (allocated(this%dphidxi)) then
       deallocate(this%dphidxi)
    end if
    nullify(this%coef)
    nullify(this%G)
    nullify(this%grad1)
    nullify(this%grad2)
    nullify(this%grad3)
    
  end subroutine gradient_jump_penalty_free

!   !> Compute eddy viscosity.
!   !! @param t The time value.
!   !! @param tstep The current time-step.
!   subroutine gradient_jump_penalty_compute(this, t, tstep)
!     class(grad_jump_penalty_t), intent(inout) :: this
!     real(kind=rp), intent(in) :: t
!     integer, intent(in) :: tstep

!     call set_ds_filt(this%test_filter)

!     if (NEKO_BCKND_DEVICE .eq. 1) then
!         call neko_error("Dynamic Smagorinsky model not implemented on accelarators.")
!     else
!         call dynamic_smagorinsky_compute_cpu(t, tstep, this%coef, this%nut, &
!                                 this%delta, this%c_dyn, this%test_filter, &
!                                 this%mij, this%lij, this%num, this%den)
!     end if

!   end subroutine gradient_jump_penalty_compute

  ! !> Integrate over a facet
  ! !! @param f Integrant
  ! subroutine integrate_over_facet(f)

  ! end subroutine integrate_over_facet

end module gradient_jump_penalty