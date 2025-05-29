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
!> Implements utilities for the TreeAMG hierarchy structure.
module tree_amg_utils
  use num_types, only : rp
  use tree_amg, only : tamg_hierarchy_t
  use gather_scatter, only : GS_OP_ADD
  implicit none
  private

  public :: tamg_sample_matrix_val

contains

  !> Sample the values in a matix (expensive, use with caution)
  !! @param val Return the matrix value in A(i,j)
  !! @param amg The TreeAMG object (the matrix A)
  !! @param lvl Current level of the cycle
  !! @param i The row index to sample A(i,j)
  !! @param j The column index to sample A(i,j)
  subroutine tamg_sample_matrix_val(val, amg, lvl, i, j)
    real(kind=rp), intent(inout) :: val
    type(tamg_hierarchy_t), intent(inout) :: amg
    integer, intent(in) :: lvl, i, j
    real(kind=rp), allocatable :: u(:), v(:)
    integer :: n
    n = amg%lvl(lvl+1)%fine_lvl_dofs
    allocate( u(n) )
    allocate( v(n) )
    u = 0.0_rp
    v = 0.0_rp

    u(i) = 1.0_rp

    !> Get the ghosts to agree on the value of 1.0
    if (lvl.eq.0) call amg%gs_h%op(u, n, GS_OP_ADD)

    call amg%matvec(v, u, lvl)

    val = v(j)
    deallocate( u )
    deallocate( v )
  end subroutine tamg_sample_matrix_val

  !> never call this function
  !> also this does not account for duplicates/ghosts in dofs
  subroutine tamg_print_matrix(amg, lvl)
    type(tamg_hierarchy_t), intent(inout) :: amg
    integer, intent(in) :: lvl
    integer :: n, i, j
    real(kind=rp) :: val
    real(kind=rp), allocatable :: u(:), v(:)

    n = amg%lvl(lvl+1)%fine_lvl_dofs
    allocate( u(n) )
    allocate( v(n) )

    do i = 1, n
       do j = 1, n
          u = 0.0_rp
          v = 0.0_rp

          u(i) = 1.0_rp
          !> Get the ghosts to agree on the value of 1.0
          if (lvl.eq.0) call amg%gs_h%op(u, n, GS_OP_ADD)

          call amg%matvec(v, u, lvl)

          val = v(j)

          print *, j, i, val
       end do
    end do

    deallocate( u )
    deallocate( v )
  end subroutine tamg_print_matrix

  !> never call this function
  !> also this does not account for duplicates/ghosts in dofs
  subroutine tamg_print_restriction_matrix(amg, lvl)
    type(tamg_hierarchy_t), intent(inout) :: amg
    integer, intent(in) :: lvl
    integer :: n, m, i, j
    real(kind=rp) :: val
    real(kind=rp), allocatable :: u(:), v(:)

    n = amg%lvl(lvl+1)%fine_lvl_dofs
    m = amg%lvl(lvl+2)%fine_lvl_dofs
    allocate( u(n) )
    allocate( v(m) )

    do i = 1, n
       do j = 1, m
          u = 0.0_rp
          v = 0.0_rp

          u(i) = 1.0_rp
          !> Get the ghosts to agree on the value of 1.0
          if (lvl.eq.0) call amg%gs_h%op(u, n, GS_OP_ADD)

          call amg%interp_f2c(v, u, lvl+1)

          val = v(j)

          print *, j, i, val
       end do
    end do

    deallocate( u )
    deallocate( v )
  end subroutine tamg_print_restriction_matrix

  !> never call this function
  !> also this does not account for duplicates/ghosts in dofs
  subroutine tamg_print_prolongation_matrix(amg, lvl)
    type(tamg_hierarchy_t), intent(inout) :: amg
    integer, intent(in) :: lvl
    integer :: n, m, i, j
    real(kind=rp) :: val
    real(kind=rp), allocatable :: u(:), v(:)

    n = amg%lvl(lvl+1)%fine_lvl_dofs
    m = amg%lvl(lvl+2)%fine_lvl_dofs
    allocate( u(m) )
    allocate( v(n) )

    do i = 1, m
       do j = 1, n
          u = 0.0_rp
          v = 0.0_rp

          u(i) = 1.0_rp
          !---!> Get the ghosts to agree on the value of 1.0
          !---if (lvl.eq.0) call amg%gs_h%op(u, n, GS_OP_ADD)

          call amg%interp_c2f(v, u, lvl+1)

          val = v(j)

          print *, j, i, val
       end do
    end do

    deallocate( u )
    deallocate( v )
  end subroutine tamg_print_prolongation_matrix


end module tree_amg_utils
