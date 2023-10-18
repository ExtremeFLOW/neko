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
module alignment_quad
  use num_types, only : i4, i8, dp
  use utils, only : neko_error
  use alignment, only : alignment_t
  implicit none
  private

  public :: alignment_quad_t

  !> number of operations different from identity
  integer, parameter :: NEKO_QUAD_NOPERATION = 7

  !> Base type for polytope alignment
  type, extends(alignment_t) :: alignment_quad_t
   contains
     !> Set quad specific operation number
     procedure, pass(this) :: init => alignment_quad_init
     !> array transformation
     procedure, pass(this) :: trans_i4 => transfrorm_quad_i4
     procedure, pass(this) :: trans_i8 => transfrorm_quad_i8
     procedure, pass(this) :: trans_dp => transfrorm_quad_dp
     !> general transformation
     generic :: trans => trans_i4, trans_i8, trans_dp
     !> inverse array transformation
     procedure, pass(this) :: trans_inv_i4 => transfrorm_inv_quad_i4
     procedure, pass(this) :: trans_inv_i8 => transfrorm_inv_quad_i8
     procedure, pass(this) :: trans_inv_dp => transfrorm_inv_quad_dp
     !> general transformation
     generic :: trans_inv => trans_inv_i4, trans_inv_i8, trans_inv_dp
  end type alignment_quad_t

contains
  !> @brief Set quad specific operation number
  subroutine  alignment_quad_init(this)
    class(alignment_quad_t), intent(inout) :: this
    call this%set_nop(NEKO_QUAD_NOPERATION)
    return
  end subroutine alignment_quad_init

  !> @brief Transform single integer array rank 2
  !! @parameter[in]     ifbnd    do we include boundary point
  !! @parameter[inout]  fcs      face data
  !! @parameter[inout]  work     work space
  subroutine transfrorm_quad_i4(this, ifbnd, fcs, work)
    class(alignment_quad_t), intent(in) :: this
    logical, intent(in) :: ifbnd
    integer(i4), dimension(:, :), intent(inout) :: fcs
    integer(i4), dimension(:), intent(inout) :: work
    ! local variables
    integer(i4) :: algn, sz, istart, iend, il, jl
    integer(i4) :: iface

    ! check alignment type; zero means identity; nothing to do
    algn = this%algn()
    if (algn /= 0) then
       ! face size
       ! I assume face is a square matrix and work has the same size
       ! possible place for check
       sz = size(fcs, 1, i4)

       ! do we work on boundary points?
       if (ifbnd) then
          istart = 1
          iend = sz
       else
          istart = 2
          iend = sz - 1
       end if
       ! apply transformations
       ! T - transpose
       ! PX - column permutation
       ! PY - row permutation
       select case(algn)
       case(1) ! T
          do jl = istart, iend ! T
             do il = istart, jl -1
                iface = fcs(il, jl)
                fcs(il, jl) = fcs(jl, il)
                fcs(jl, il) = iface
             end do
          end do
       case(2) ! PX
          do jl = istart, iend ! PX
             do il = istart, sz/2
                iface = fcs(il, jl)
                fcs(il, jl) = fcs(sz + 1 - il, jl)
                fcs(sz + 1 - il, jl) = iface
             end do
          end do
       case(3) ! PXT = TPY
          do jl = istart, iend ! PX
             do il = istart, sz/2
                iface = fcs(il, jl)
                fcs(il, jl) = fcs(sz + 1 - il, jl)
                fcs(sz + 1 - il, jl) = iface
             end do
          end do
          do jl = istart, iend ! T
             do il = istart, jl -1
                iface = fcs(il, jl)
                fcs(il, jl) = fcs(jl, il)
                fcs(jl, il) = iface
             end do
          end do
       case(4) ! PYT = TPX
          do jl = istart, sz/2 ! PY
             work(istart:iend) = fcs(istart:iend, jl)
             fcs(istart:iend, jl) = fcs(istart:iend, sz + 1 - jl)
             fcs(istart:iend, sz + 1 - jl) = work(istart:iend)
          end do
          ! or
          !do jl = istart, sz/2 ! PY
          !   do il = istart, iend
          !      iface = fcs(il, jl)
          !      fcs(il, jl) = fcs(il, sz + 1 - jl)
          !      fcs(il, sz + 1 - jl) = iface
          !   end do
          !end do
          do jl = istart, iend ! T
             do il = istart, jl -1
                iface = fcs(il, jl)
                fcs(il, jl) = fcs(jl, il)
                fcs(jl, il) = iface
             end do
          end do
       case(5) ! PY
          do jl = istart, sz/2 ! PY
             work(istart:iend) = fcs(istart:iend, jl)
             fcs(istart:iend, jl) = fcs(istart:iend, sz + 1 - jl)
             fcs(istart:iend, sz + 1 - jl) = work(istart:iend)
          end do
          ! or
          !do jl = istart, sz/2 ! PY
          !   do il = istart, iend
          !      iface = fcs(il, jl)
          !      fcs(il, jl) = fcs(il, sz + 1 - jl)
          !      fcs(il, sz + 1 - jl) = iface
          !   end do
          !end do
       case(6) ! PXPYT = PYPXT = TPXPY = TPYPX
          do jl = istart, iend ! PX
             do il = istart, sz/2
                iface = fcs(il, jl)
                fcs(il, jl) = fcs(sz + 1 - il, jl)
                fcs(sz + 1 - il, jl) = iface
             end do
          end do
          do jl = istart, sz/2 ! PY
             work(istart:iend) = fcs(istart:iend, jl)
             fcs(istart:iend, jl) = fcs(istart:iend, sz + 1 - jl)
             fcs(istart:iend, sz + 1 - jl) = work(istart:iend)
          end do
          ! or
          !do jl = istart, sz/2 ! PY
          !   do il = istart, iend
          !      iface = fcs(il, jl)
          !      fcs(il, jl) = fcs(il, sz + 1 - jl)
          !      fcs(il, sz + 1 - jl) = iface
          !   end do
          !end do
          do jl = istart, iend ! T
             do il = istart, jl -1
                iface = fcs(il, jl)
                fcs(il, jl) = fcs(jl, il)
                fcs(jl, il) = iface
             end do
          end do
       case(7) ! PXPY = PYPX
          do jl = istart, iend ! PX
             do il = istart, sz/2
                iface = fcs(il, jl)
                fcs(il, jl) = fcs(sz + 1 - il, jl)
                fcs(sz + 1 - il, jl) = iface
             end do
          end do
          do jl = istart, sz/2 ! PY
             work(istart:iend) = fcs(istart:iend, jl)
             fcs(istart:iend, jl) = fcs(istart:iend, sz + 1 - jl)
             fcs(istart:iend, sz + 1 - jl) = work(istart:iend)
          end do
       case default
          call neko_error('Quad alignment not initialised properly')
       end select
    end if

    return
  end subroutine transfrorm_quad_i4

  !> @brief Transform double integer array rank 2
  !! @parameter[in]     ifbnd    do we include boundary point
  !! @parameter[inout]  fcs      face data
  !! @parameter[inout]  work     work space
  subroutine transfrorm_quad_i8(this, ifbnd, fcs, work)
    class(alignment_quad_t), intent(in) :: this
    logical, intent(in) :: ifbnd
    integer(i8), dimension(:, :), intent(inout) :: fcs
    integer(i8), dimension(:), intent(inout) :: work
    ! local variables
    integer(i4) :: algn, sz, istart, iend, il, jl
    integer(i8) :: iface

    ! check alignment type; zero means identity; nothing to do
    algn = this%algn()
    if (algn /= 0) then
       ! face size
       ! I assume face is a square matrix and work has the same size
       ! possible place for check
       sz = size(fcs, 1, i4)

       ! do we work on boundary points?
       if (ifbnd) then
          istart = 1
          iend = sz
       else
          istart = 2
          iend = sz - 1
       end if
       ! apply transformations
       ! T - transpose
       ! PX - column permutation
       ! PY - row permutation
       select case(algn)
       case(1) ! T
          do jl = istart, iend ! T
             do il = istart, jl -1
                iface = fcs(il, jl)
                fcs(il, jl) = fcs(jl, il)
                fcs(jl, il) = iface
             end do
          end do
       case(2) ! PX
          do jl = istart, iend ! PX
             do il = istart, sz/2
                iface = fcs(il, jl)
                fcs(il, jl) = fcs(sz + 1 - il, jl)
                fcs(sz + 1 - il, jl) = iface
             end do
          end do
       case(3) ! PXT = TPY
          do jl = istart, iend ! PX
             do il = istart, sz/2
                iface = fcs(il, jl)
                fcs(il, jl) = fcs(sz + 1 - il, jl)
                fcs(sz + 1 - il, jl) = iface
             end do
          end do
          do jl = istart, iend ! T
             do il = istart, jl -1
                iface = fcs(il, jl)
                fcs(il, jl) = fcs(jl, il)
                fcs(jl, il) = iface
             end do
          end do
       case(4) ! PYT = TPX
          do jl = istart, sz/2 ! PY
             work(istart:iend) = fcs(istart:iend, jl)
             fcs(istart:iend, jl) = fcs(istart:iend, sz + 1 - jl)
             fcs(istart:iend, sz + 1 - jl) = work(istart:iend)
          end do
          ! or
          !do jl = istart, sz/2 ! PY
          !   do il = istart, iend
          !      iface = fcs(il, jl)
          !      fcs(il, jl) = fcs(il, sz + 1 - jl)
          !      fcs(il, sz + 1 - jl) = iface
          !   end do
          !end do
          do jl = istart, iend ! T
             do il = istart, jl -1
                iface = fcs(il, jl)
                fcs(il, jl) = fcs(jl, il)
                fcs(jl, il) = iface
             end do
          end do
       case(5) ! PY
          do jl = istart, sz/2 ! PY
             work(istart:iend) = fcs(istart:iend, jl)
             fcs(istart:iend, jl) = fcs(istart:iend, sz + 1 - jl)
             fcs(istart:iend, sz + 1 - jl) = work(istart:iend)
          end do
          ! or
          !do jl = istart, sz/2 ! PY
          !   do il = istart, iend
          !      iface = fcs(il, jl)
          !      fcs(il, jl) = fcs(il, sz + 1 - jl)
          !      fcs(il, sz + 1 - jl) = iface
          !   end do
          !end do
       case(6) ! PXPYT = PYPXT = TPXPY = TPYPX
          do jl = istart, iend ! PX
             do il = istart, sz/2
                iface = fcs(il, jl)
                fcs(il, jl) = fcs(sz + 1 - il, jl)
                fcs(sz + 1 - il, jl) = iface
             end do
          end do
          do jl = istart, sz/2 ! PY
             work(istart:iend) = fcs(istart:iend, jl)
             fcs(istart:iend, jl) = fcs(istart:iend, sz + 1 - jl)
             fcs(istart:iend, sz + 1 - jl) = work(istart:iend)
          end do
          ! or
          !do jl = istart, sz/2 ! PY
          !   do il = istart, iend
          !      iface = fcs(il, jl)
          !      fcs(il, jl) = fcs(il, sz + 1 - jl)
          !      fcs(il, sz + 1 - jl) = iface
          !   end do
          !end do
          do jl = istart, iend ! T
             do il = istart, jl -1
                iface = fcs(il, jl)
                fcs(il, jl) = fcs(jl, il)
                fcs(jl, il) = iface
             end do
          end do
       case(7) ! PXPY = PYPX
          do jl = istart, iend ! PX
             do il = istart, sz/2
                iface = fcs(il, jl)
                fcs(il, jl) = fcs(sz + 1 - il, jl)
                fcs(sz + 1 - il, jl) = iface
             end do
          end do
          do jl = istart, sz/2 ! PY
             work(istart:iend) = fcs(istart:iend, jl)
             fcs(istart:iend, jl) = fcs(istart:iend, sz + 1 - jl)
             fcs(istart:iend, sz + 1 - jl) = work(istart:iend)
          end do
       case default
          call neko_error('Quad alignment not initialised properly')
       end select
    end if

    return
  end subroutine transfrorm_quad_i8

  !> @brief Transform double real array rank 2
  !! @parameter[in]     ifbnd    do we include boundary point
  !! @parameter[inout]  fcs      face data
  !! @parameter[inout]  work     work space
  subroutine transfrorm_quad_dp(this, ifbnd, fcs, work)
    class(alignment_quad_t), intent(in) :: this
    logical, intent(in) :: ifbnd
    real(dp), dimension(:, :), intent(inout) :: fcs
    real(dp), dimension(:), intent(inout) :: work
    ! local variables
    integer(i4) :: algn, sz, istart, iend, il, jl
    real(dp) :: rface

    ! check alignment type; zero means identity; nothing to do
    algn = this%algn()
    if (algn /= 0) then
       ! face size
       ! I assume face is a square matrix and work has the same size
       ! possible place for check
       sz = size(fcs, 1, i4)

       ! do we work on boundary points?
       if (ifbnd) then
          istart = 1
          iend = sz
       else
          istart = 2
          iend = sz - 1
       end if
       ! apply transformations
       ! T - transpose
       ! PX - column permutation
       ! PY - row permutation
       select case(algn)
       case(1) ! T
          do jl = istart, iend ! T
             do il = istart, jl -1
                rface = fcs(il, jl)
                fcs(il, jl) = fcs(jl, il)
                fcs(jl, il) = rface
             end do
          end do
       case(2) ! PX
          do jl = istart, iend ! PX
             do il = istart, sz/2
                rface = fcs(il, jl)
                fcs(il, jl) = fcs(sz + 1 - il, jl)
                fcs(sz + 1 - il, jl) = rface
             end do
          end do
       case(3) ! PXT = TPY
          do jl = istart, iend ! PX
             do il = istart, sz/2
                rface = fcs(il, jl)
                fcs(il, jl) = fcs(sz + 1 - il, jl)
                fcs(sz + 1 - il, jl) = rface
             end do
          end do
          do jl = istart, iend ! T
             do il = istart, jl -1
                rface = fcs(il, jl)
                fcs(il, jl) = fcs(jl, il)
                fcs(jl, il) = rface
             end do
          end do
       case(4) ! PYT = TPX
          do jl = istart, sz/2 ! PY
             work(istart:iend) = fcs(istart:iend, jl)
             fcs(istart:iend, jl) = fcs(istart:iend, sz + 1 - jl)
             fcs(istart:iend, sz + 1 - jl) = work(istart:iend)
          end do
          ! or
          !do jl = istart, sz/2 ! PY
          !   do il = istart, iend
          !      rface = fcs(il, jl)
          !      fcs(il, jl) = fcs(il, sz + 1 - jl)
          !      fcs(il, sz + 1 - jl) = rface
          !   end do
          !end do
          do jl = istart, iend ! T
             do il = istart, jl -1
                rface = fcs(il, jl)
                fcs(il, jl) = fcs(jl, il)
                fcs(jl, il) = rface
             end do
          end do
       case(5) ! PY
          do jl = istart, sz/2 ! PY
             work(istart:iend) = fcs(istart:iend, jl)
             fcs(istart:iend, jl) = fcs(istart:iend, sz + 1 - jl)
             fcs(istart:iend, sz + 1 - jl) = work(istart:iend)
          end do
          ! or
          !do jl = istart, sz/2 ! PY
          !   do il = istart, iend
          !      rface = fcs(il, jl)
          !      fcs(il, jl) = fcs(il, sz + 1 - jl)
          !      fcs(il, sz + 1 - jl) = rface
          !   end do
          !end do
       case(6) ! PXPYT = PYPXT = TPXPY = TPYPX
          do jl = istart, iend ! PX
             do il = istart, sz/2
                rface = fcs(il, jl)
                fcs(il, jl) = fcs(sz + 1 - il, jl)
                fcs(sz + 1 - il, jl) = rface
             end do
          end do
          do jl = istart, sz/2 ! PY
             work(istart:iend) = fcs(istart:iend, jl)
             fcs(istart:iend, jl) = fcs(istart:iend, sz + 1 - jl)
             fcs(istart:iend, sz + 1 - jl) = work(istart:iend)
          end do
          ! or
          !do jl = istart, sz/2 ! PY
          !   do il = istart, iend
          !      rface = fcs(il, jl)
          !      fcs(il, jl) = fcs(il, sz + 1 - jl)
          !      fcs(il, sz + 1 - jl) = rface
          !   end do
          !end do
          do jl = istart, iend ! T
             do il = istart, jl -1
                rface = fcs(il, jl)
                fcs(il, jl) = fcs(jl, il)
                fcs(jl, il) = rface
             end do
          end do
       case(7) ! PXPY = PYPX
          do jl = istart, iend ! PX
             do il = istart, sz/2
                rface = fcs(il, jl)
                fcs(il, jl) = fcs(sz + 1 - il, jl)
                fcs(sz + 1 - il, jl) = rface
             end do
          end do
          do jl = istart, sz/2 ! PY
             work(istart:iend) = fcs(istart:iend, jl)
             fcs(istart:iend, jl) = fcs(istart:iend, sz + 1 - jl)
             fcs(istart:iend, sz + 1 - jl) = work(istart:iend)
          end do
       case default
          call neko_error('Quad alignment not initialised properly')
       end select
    end if

    return
  end subroutine transfrorm_quad_dp

  !> @brief Inverse transform single integer array rank 2
  !! @parameter[in]     ifbnd    do we include boundary point
  !! @parameter[inout]  fcs      face data
  !! @parameter[inout]  work     work space
  subroutine transfrorm_inv_quad_i4(this, ifbnd, fcs, work)
    class(alignment_quad_t), intent(in) :: this
    logical, intent(in) :: ifbnd
    integer(i4), dimension(:, :), intent(inout) :: fcs
    integer(i4), dimension(:), intent(inout) :: work
    ! local variables
    integer(i4) :: algn, sz, istart, iend, il, jl
    integer(i4) :: iface

    ! check alignment type; zero means identity; nothing to do
    algn = this%algn()
    if (algn /= 0) then
       ! face size
       ! I assume face is a square matrix and work has the same size
       ! possible place for check
       sz = size(fcs, 1, i4)

       ! do we work on boundary points?
       if (ifbnd) then
          istart = 1
          iend = sz
       else
          istart = 2
          iend = sz - 1
       end if
       ! apply transformations
       ! T - transpose
       ! PX - column permutation
       ! PY - row permutation
       select case(algn)
       case(1) ! T
          do jl = istart, iend ! T
             do il = istart, jl -1
                iface = fcs(il, jl)
                fcs(il, jl) = fcs(jl, il)
                fcs(jl, il) = iface
             end do
          end do
       case(2) ! PX
          do jl = istart, iend ! PX
             do il = istart, sz/2
                iface = fcs(il, jl)
                fcs(il, jl) = fcs(sz + 1 - il, jl)
                fcs(sz + 1 - il, jl) = iface
             end do
          end do
       case(3) ! PYT = TPX
          do jl = istart, sz/2 ! PY
             work(istart:iend) = fcs(istart:iend, jl)
             fcs(istart:iend, jl) = fcs(istart:iend, sz + 1 - jl)
             fcs(istart:iend, sz + 1 - jl) = work(istart:iend)
          end do
          ! or
          !do jl = istart, sz/2 ! PY
          !   do il = istart, iend
          !      iface = fcs(il, jl)
          !      fcs(il, jl) = fcs(il, sz + 1 - jl)
          !      fcs(il, sz + 1 - jl) = iface
          !   end do
          !end do
          do jl = istart, iend ! T
             do il = istart, jl -1
                iface = fcs(il, jl)
                fcs(il, jl) = fcs(jl, il)
                fcs(jl, il) = iface
             end do
          end do
       case(4) ! PXT = TPY
          do jl = istart, iend ! PX
             do il = istart, sz/2
                iface = fcs(il, jl)
                fcs(il, jl) = fcs(sz + 1 - il, jl)
                fcs(sz + 1 - il, jl) = iface
             end do
          end do
          do jl = istart, iend ! T
             do il = istart, jl -1
                iface = fcs(il, jl)
                fcs(il, jl) = fcs(jl, il)
                fcs(jl, il) = iface
             end do
          end do
       case(5) ! PY
          do jl = istart, sz/2 ! PY
             work(istart:iend) = fcs(istart:iend, jl)
             fcs(istart:iend, jl) = fcs(istart:iend, sz + 1 - jl)
             fcs(istart:iend, sz + 1 - jl) = work(istart:iend)
          end do
          ! or
          !do jl = istart, sz/2 ! PY
          !   do il = istart, iend
          !      iface = fcs(il, jl)
          !      fcs(il, jl) = fcs(il, sz + 1 - jl)
          !      fcs(il, sz + 1 - jl) = iface
          !   end do
          !end do
       case(6) ! PXPYT = PYPXT = TPXPY = TPYPX
          do jl = istart, iend ! PX
             do il = istart, sz/2
                iface = fcs(il, jl)
                fcs(il, jl) = fcs(sz + 1 - il, jl)
                fcs(sz + 1 - il, jl) = iface
             end do
          end do
          do jl = istart, sz/2 ! PY
             work(istart:iend) = fcs(istart:iend, jl)
             fcs(istart:iend, jl) = fcs(istart:iend, sz + 1 - jl)
             fcs(istart:iend, sz + 1 - jl) = work(istart:iend)
          end do
          ! or
          !do jl = istart, sz/2 ! PY
          !   do il = istart, iend
          !      iface = fcs(il, jl)
          !      fcs(il, jl) = fcs(il, sz + 1 - jl)
          !      fcs(il, sz + 1 - jl) = iface
          !   end do
          !end do
          do jl = istart, iend ! T
             do il = istart, jl -1
                iface = fcs(il, jl)
                fcs(il, jl) = fcs(jl, il)
                fcs(jl, il) = iface
             end do
          end do
       case(7) ! PXPY = PYPX
          do jl = istart, iend ! PX
             do il = istart, sz/2
                iface = fcs(il, jl)
                fcs(il, jl) = fcs(sz + 1 - il, jl)
                fcs(sz + 1 - il, jl) = iface
             end do
          end do
          do jl = istart, sz/2 ! PY
             work(istart:iend) = fcs(istart:iend, jl)
             fcs(istart:iend, jl) = fcs(istart:iend, sz + 1 - jl)
             fcs(istart:iend, sz + 1 - jl) = work(istart:iend)
          end do
       case default
          call neko_error('Quad alignment not initialised properly')
       end select
    end if

    return
  end subroutine transfrorm_inv_quad_i4

  !> @brief Inverse transform double integer array rank 2
  !! @parameter[in]     ifbnd    do we include boundary point
  !! @parameter[inout]  fcs      face data
  !! @parameter[inout]  work     work space
  subroutine transfrorm_inv_quad_i8(this, ifbnd, fcs, work)
    class(alignment_quad_t), intent(in) :: this
    logical, intent(in) :: ifbnd
    integer(i8), dimension(:, :), intent(inout) :: fcs
    integer(i8), dimension(:), intent(inout) :: work
    ! local variables
    integer(i4) :: algn, sz, istart, iend, il, jl
    integer(i8) :: iface

    ! check alignment type; zero means identity; nothing to do
    algn = this%algn()
    if (algn /= 0) then
       ! face size
       ! I assume face is a square matrix and work has the same size
       ! possible place for check
       sz = size(fcs, 1, i4)

       ! do we work on boundary points?
       if (ifbnd) then
          istart = 1
          iend = sz
       else
          istart = 2
          iend = sz - 1
       end if
       ! apply transformations
       ! T - transpose
       ! PX - column permutation
       ! PY - row permutation
       select case(algn)
       case(1) ! T
          do jl = istart, iend ! T
             do il = istart, jl -1
                iface = fcs(il, jl)
                fcs(il, jl) = fcs(jl, il)
                fcs(jl, il) = iface
             end do
          end do
       case(2) ! PX
          do jl = istart, iend ! PX
             do il = istart, sz/2
                iface = fcs(il, jl)
                fcs(il, jl) = fcs(sz + 1 - il, jl)
                fcs(sz + 1 - il, jl) = iface
             end do
          end do
       case(3) ! PYT = TPX
          do jl = istart, sz/2 ! PY
             work(istart:iend) = fcs(istart:iend, jl)
             fcs(istart:iend, jl) = fcs(istart:iend, sz + 1 - jl)
             fcs(istart:iend, sz + 1 - jl) = work(istart:iend)
          end do
          ! or
          !do jl = istart, sz/2 ! PY
          !   do il = istart, iend
          !      iface = fcs(il, jl)
          !      fcs(il, jl) = fcs(il, sz + 1 - jl)
          !      fcs(il, sz + 1 - jl) = iface
          !   end do
          !end do
          do jl = istart, iend ! T
             do il = istart, jl -1
                iface = fcs(il, jl)
                fcs(il, jl) = fcs(jl, il)
                fcs(jl, il) = iface
             end do
          end do
       case(4) ! PXT = TPY
          do jl = istart, iend ! PX
             do il = istart, sz/2
                iface = fcs(il, jl)
                fcs(il, jl) = fcs(sz + 1 - il, jl)
                fcs(sz + 1 - il, jl) = iface
             end do
          end do
          do jl = istart, iend ! T
             do il = istart, jl -1
                iface = fcs(il, jl)
                fcs(il, jl) = fcs(jl, il)
                fcs(jl, il) = iface
             end do
          end do
       case(5) ! PY
          do jl = istart, sz/2 ! PY
             work(istart:iend) = fcs(istart:iend, jl)
             fcs(istart:iend, jl) = fcs(istart:iend, sz + 1 - jl)
             fcs(istart:iend, sz + 1 - jl) = work(istart:iend)
          end do
          ! or
          !do jl = istart, sz/2 ! PY
          !   do il = istart, iend
          !      iface = fcs(il, jl)
          !      fcs(il, jl) = fcs(il, sz + 1 - jl)
          !      fcs(il, sz + 1 - jl) = iface
          !   end do
          !end do
       case(6) ! PXPYT = PYPXT = TPXPY = TPYPX
          do jl = istart, iend ! PX
             do il = istart, sz/2
                iface = fcs(il, jl)
                fcs(il, jl) = fcs(sz + 1 - il, jl)
                fcs(sz + 1 - il, jl) = iface
             end do
          end do
          do jl = istart, sz/2 ! PY
             work(istart:iend) = fcs(istart:iend, jl)
             fcs(istart:iend, jl) = fcs(istart:iend, sz + 1 - jl)
             fcs(istart:iend, sz + 1 - jl) = work(istart:iend)
          end do
          ! or
          !do jl = istart, sz/2 ! PY
          !   do il = istart, iend
          !      iface = fcs(il, jl)
          !      fcs(il, jl) = fcs(il, sz + 1 - jl)
          !      fcs(il, sz + 1 - jl) = iface
          !   end do
          !end do
          do jl = istart, iend ! T
             do il = istart, jl -1
                iface = fcs(il, jl)
                fcs(il, jl) = fcs(jl, il)
                fcs(jl, il) = iface
             end do
          end do
       case(7) ! PXPY = PYPX
          do jl = istart, iend ! PX
             do il = istart, sz/2
                iface = fcs(il, jl)
                fcs(il, jl) = fcs(sz + 1 - il, jl)
                fcs(sz + 1 - il, jl) = iface
             end do
          end do
          do jl = istart, sz/2 ! PY
             work(istart:iend) = fcs(istart:iend, jl)
             fcs(istart:iend, jl) = fcs(istart:iend, sz + 1 - jl)
             fcs(istart:iend, sz + 1 - jl) = work(istart:iend)
          end do
       case default
          call neko_error('Quad alignment not initialised properly')
       end select
    end if

    return
  end subroutine transfrorm_inv_quad_i8

  !> @brief Inverse transform double real array rank 2
  !! @parameter[in]     ifbnd    do we include boundary point
  !! @parameter[inout]  fcs      face data
  !! @parameter[inout]  work     work space
  subroutine transfrorm_inv_quad_dp(this, ifbnd, fcs, work)
    class(alignment_quad_t), intent(in) :: this
    logical, intent(in) :: ifbnd
    real(dp), dimension(:, :), intent(inout) :: fcs
    real(dp), dimension(:), intent(inout) :: work
    ! local variables
    integer(i4) :: algn, sz, istart, iend, il, jl
    real(dp) :: rface

    ! check alignment type; zero means identity; nothing to do
    algn = this%algn()
    if (algn /= 0) then
       ! face size
       ! I assume face is a square matrix and work has the same size
       ! possible place for check
       sz = size(fcs, 1, i4)

       ! do we work on boundary points?
       if (ifbnd) then
          istart = 1
          iend = sz
       else
          istart = 2
          iend = sz - 1
       end if
       ! apply transformations
       ! T - transpose
       ! PX - column permutation
       ! PY - row permutation
       select case(algn)
       case(1) ! T
          do jl = istart, iend ! T
             do il = istart, jl -1
                rface = fcs(il, jl)
                fcs(il, jl) = fcs(jl, il)
                fcs(jl, il) = rface
             end do
          end do
       case(2) ! PX
          do jl = istart, iend ! PX
             do il = istart, sz/2
                rface = fcs(il, jl)
                fcs(il, jl) = fcs(sz + 1 - il, jl)
                fcs(sz + 1 - il, jl) = rface
             end do
          end do
       case(3) ! PYT = TPX
          do jl = istart, sz/2 ! PY
             work(istart:iend) = fcs(istart:iend, jl)
             fcs(istart:iend, jl) = fcs(istart:iend, sz + 1 - jl)
             fcs(istart:iend, sz + 1 - jl) = work(istart:iend)
          end do
          ! or
          !do jl = istart, sz/2 ! PY
          !   do il = istart, iend
          !      rface = fcs(il, jl)
          !      fcs(il, jl) = fcs(il, sz + 1 - jl)
          !      fcs(il, sz + 1 - jl) = rface
          !   end do
          !end do
          do jl = istart, iend ! T
             do il = istart, jl -1
                rface = fcs(il, jl)
                fcs(il, jl) = fcs(jl, il)
                fcs(jl, il) = rface
             end do
          end do
       case(4) ! PXT = TPY
          do jl = istart, iend ! PX
             do il = istart, sz/2
                rface = fcs(il, jl)
                fcs(il, jl) = fcs(sz + 1 - il, jl)
                fcs(sz + 1 - il, jl) = rface
             end do
          end do
          do jl = istart, iend ! T
             do il = istart, jl -1
                rface = fcs(il, jl)
                fcs(il, jl) = fcs(jl, il)
                fcs(jl, il) = rface
             end do
          end do
       case(5) ! PY
          do jl = istart, sz/2 ! PY
             work(istart:iend) = fcs(istart:iend, jl)
             fcs(istart:iend, jl) = fcs(istart:iend, sz + 1 - jl)
             fcs(istart:iend, sz + 1 - jl) = work(istart:iend)
          end do
          ! or
          !do jl = istart, sz/2 ! PY
          !   do il = istart, iend
          !      rface = fcs(il, jl)
          !      fcs(il, jl) = fcs(il, sz + 1 - jl)
          !      fcs(il, sz + 1 - jl) = rface
          !   end do
          !end do
       case(6) ! PXPYT = PYPXT = TPXPY = TPYPX
          do jl = istart, iend ! PX
             do il = istart, sz/2
                rface = fcs(il, jl)
                fcs(il, jl) = fcs(sz + 1 - il, jl)
                fcs(sz + 1 - il, jl) = rface
             end do
          end do
          do jl = istart, sz/2 ! PY
             work(istart:iend) = fcs(istart:iend, jl)
             fcs(istart:iend, jl) = fcs(istart:iend, sz + 1 - jl)
             fcs(istart:iend, sz + 1 - jl) = work(istart:iend)
          end do
          ! or
          !do jl = istart, sz/2 ! PY
          !   do il = istart, iend
          !      rface = fcs(il, jl)
          !      fcs(il, jl) = fcs(il, sz + 1 - jl)
          !      fcs(il, sz + 1 - jl) = rface
          !   end do
          !end do
          do jl = istart, iend ! T
             do il = istart, jl -1
                rface = fcs(il, jl)
                fcs(il, jl) = fcs(jl, il)
                fcs(jl, il) = rface
             end do
          end do
       case(7) ! PXPY = PYPX
          do jl = istart, iend ! PX
             do il = istart, sz/2
                rface = fcs(il, jl)
                fcs(il, jl) = fcs(sz + 1 - il, jl)
                fcs(sz + 1 - il, jl) = rface
             end do
          end do
          do jl = istart, sz/2 ! PY
             work(istart:iend) = fcs(istart:iend, jl)
             fcs(istart:iend, jl) = fcs(istart:iend, sz + 1 - jl)
             fcs(istart:iend, sz + 1 - jl) = work(istart:iend)
          end do
       case default
          call neko_error('Quad alignment not initialised properly')
       end select
    end if

    return
  end subroutine transfrorm_inv_quad_dp

end module alignment_quad
