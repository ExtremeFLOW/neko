! Copyright (c) 2018-2023, The Neko Authors
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
!> Edge alignment
module alignment_edge
  use num_types, only : i4, i8, dp
  use utils, only : neko_error
  use alignment, only : alignment_t
  implicit none
  private

  public :: alignment_edge_t

  !> number of operations different from identity
  integer, parameter :: NEKO_EDGE_NOPERATION = 1

  !> Base type for polytope alignment
  type, extends(alignment_t) :: alignment_edge_t
   contains
     !> Set edge specific operation number
     procedure, pass(this) :: init => alignment_edge_init
     !> array transformation
     procedure, pass(this) :: trans_i4 => transfrorm_edge_i4
     procedure, pass(this) :: trans_i8 => transfrorm_edge_i8
     procedure, pass(this) :: trans_dp => transfrorm_edge_dp
     !> general transformation
     generic :: trans => trans_i4, trans_i8, trans_dp
     !> inverse array transformation
     procedure, pass(this) :: trans_inv_i4 => transfrorm_edge_i4
     procedure, pass(this) :: trans_inv_i8 => transfrorm_edge_i8
     procedure, pass(this) :: trans_inv_dp => transfrorm_edge_dp
     !> general transformation
     generic :: trans_inv => trans_inv_i4, trans_inv_i8, trans_inv_dp
  end type alignment_edge_t

contains
  !> @brief Set edge specific operation number
  subroutine  alignment_edge_init(this)
    class(alignment_edge_t), intent(inout) :: this
    call this%set_nop(NEKO_EDGE_NOPERATION)
    return
  end subroutine alignment_edge_init

  !> @brief Transform single integer array rank 1
  !! @parameter[in]     ifbnd    do we include boundary points
  !! @parameter[inout]  edg      edge data
  subroutine transfrorm_edge_i4(this, ifbnd, edg)
    class(alignment_edge_t), intent(in) :: this
    logical, intent(in) :: ifbnd
    integer(i4), dimension(:), intent(inout) :: edg
    ! local variables
    integer(i4) :: algn, sz, istart, il, itmp1, itmp2
    integer(i4) :: iedg

    ! check alignment type; zero means identity; nothing to do
    algn = this%algn()
    if (algn /= 0) then
       ! edge size
       sz = size(edg, 1, i4)
       ! do we work on boundary points?
       if (ifbnd) then
          istart = 1
       else
          istart = 2
       end if
       ! apply transformations
       ! P - permutation
       select case(algn)
       case(1) ! P
          itmp1 = sz + 1
          do il = istart, sz/2
             itmp2 = itmp1 - il
             iedg = edg(il)
             edg(il) = edg(itmp2)
             edg(itmp2) = iedg
          end do
       case default
          call neko_error('Edge alignment not initialised properly')
       end select
    end if

    return
  end subroutine transfrorm_edge_i4

  !> @brief Transform double integer array rank 1
  !! @parameter[in]     ifbnd    do we include boundary points
  !! @parameter[inout]  edg      edge data
  subroutine transfrorm_edge_i8(this, ifbnd, edg)
    class(alignment_edge_t), intent(in) :: this
    logical, intent(in) :: ifbnd
    integer(i8), dimension(:), intent(inout) :: edg
    ! local variables
    integer(i4) :: algn, sz, istart, il, itmp1, itmp2
    integer(i8) :: iedg

    ! check alignment type; zero means identity; nothing to do
    algn = this%algn()
    if (algn /= 0) then
       ! edge size
       sz = size(edg, 1, i4)
       ! do we work on boundary points?
       if (ifbnd) then
          istart = 1
       else
          istart = 2
       end if
       ! apply transformations
       ! P - permutation
       select case(algn)
       case(1) ! P
          itmp1 = sz + 1
          do il = istart, sz/2
             itmp2 = itmp1 - il
             iedg = edg(il)
             edg(il) = edg(itmp2)
             edg(itmp2) = iedg
          end do
       case default
          call neko_error('Edge alignment not initialised properly')
       end select
    end if

    return
  end subroutine transfrorm_edge_i8

  !> @brief Transform double real array rank 1
  !! @parameter[in]     ifbnd    do we include boundary points
  !! @parameter[inout]  edg      edge data
  subroutine transfrorm_edge_dp(this, ifbnd, edg)
    class(alignment_edge_t), intent(in) :: this
    logical, intent(in) :: ifbnd
    real(dp), dimension(:), intent(inout) :: edg
    ! local variables
    integer(i4) :: algn, sz, istart, il, itmp1, itmp2
    real(dp) :: redg

    ! check alignment type; zero means identity; nothing to do
    algn = this%algn()
    if (algn /= 0) then
       ! edge size
       sz = size(edg, 1, i4)
       ! do we work on boundary points?
       if (ifbnd) then
          istart = 1
       else
          istart = 2
       end if
       ! apply transformations
       ! P - permutation
       select case(algn)
       case(1) ! P
          itmp1 = sz + 1
          do il = istart, sz/2
             itmp2 = itmp1 - il
             redg = edg(il)
             edg(il) = edg(itmp2)
             edg(itmp2) = redg
          end do
       case default
          call neko_error('Edge alignment not initialised properly')
       end select
    end if

    return
  end subroutine transfrorm_edge_dp

end module alignment_edge
