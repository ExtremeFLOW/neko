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
!> Generic Gather-scatter backend for NEC Vector Engines
module gs_sx
  use num_types
  use gs_bcknd
  use gs_ops
  implicit none
  private

  !> Gather-scatter backend for NEC SX-Aurora
  type, public, extends(gs_bcknd_t) :: gs_sx_t
     real(kind=rp), allocatable :: local_wrk(:)
     real(kind=rp), allocatable :: shared_wrk(:)
     integer :: nlocal
     integer :: nshared
   contains
     procedure, pass(this) :: init => gs_sx_init
     procedure, pass(this) :: free => gs_sx_free
     procedure, pass(this) :: gather => gs_gather_sx
     procedure, pass(this) :: scatter => gs_scatter_sx
  end type gs_sx_t
  
contains
  
  !> SX backend initialisation
  subroutine gs_sx_init(this, nlocal, nshared, nlcl_blks, nshrd_blks)
    class(gs_sx_t), intent(inout) :: this
    integer, intent(in) :: nlocal
    integer, intent(in) :: nshared
    integer, intent(in) :: nlcl_blks
    integer, intent(in) :: nshrd_blks
       
    call this%free()

    this%nlocal = nlocal
    this%nshared = nshared

    allocate(this%local_wrk(nlocal))
    allocate(this%shared_wrk(nshared))
    
  end subroutine gs_sx_init

  !> SX backend deallocation
  subroutine gs_sx_free(this)
    class(gs_sx_t), intent(inout) :: this

    if (allocated(this%local_wrk)) then
       deallocate(this%local_wrk)
    end if

    if (allocated(this%shared_wrk)) then
       deallocate(this%shared_wrk)
    end if

    this%nlocal = 0
    this%nshared = 0
    
  end subroutine gs_sx_free

  !> Gather kernel
  subroutine gs_gather_sx(this, v, m, o, dg, u, n, gd, nb, b, op, shrd)
    integer, intent(inout) :: m
    integer, intent(inout) :: n
    integer, intent(inout) :: nb
    class(gs_sx_t), intent(inout) :: this
    real(kind=rp), dimension(m), intent(inout) :: v
    integer, dimension(m), intent(inout) :: dg
    real(kind=rp), dimension(n), intent(inout) :: u
    integer, dimension(m), intent(inout) :: gd
    integer, dimension(nb), intent(inout) :: b
    integer, intent(inout) :: o
    integer, intent(inout) :: op
    logical, intent(in) :: shrd

    if (.not. shrd) then
       associate(w=>this%local_wrk)
         select case(op)
         case (GS_OP_ADD)
            call gs_gather_kernel_add(v, m, o, dg, u, n, gd, nb, b, w)
         case (GS_OP_MUL)
            call gs_gather_kernel_mul(v, m, o, dg, u, n, gd, nb, b, w)
         case (GS_OP_MIN)
            call gs_gather_kernel_min(v, m, o, dg, u, n, gd, nb, b, w)
         case (GS_OP_MAX)
            call gs_gather_kernel_max(v, m, o, dg, u, n, gd, nb, b, w)
         end select
       end associate
    else if (shrd) then
       associate(w=>this%shared_wrk)
         select case(op)
         case (GS_OP_ADD)
            call gs_gather_kernel_add(v, m, o, dg, u, n, gd, nb, b, w)
         case (GS_OP_MUL)
            call gs_gather_kernel_mul(v, m, o, dg, u, n, gd, nb, b, w)
         case (GS_OP_MIN)
            call gs_gather_kernel_min(v, m, o, dg, u, n, gd, nb, b, w)
         case (GS_OP_MAX)
            call gs_gather_kernel_max(v, m, o, dg, u, n, gd, nb, b, w)
         end select
       end associate
    end if
    
  end subroutine gs_gather_sx
 
  !> Gather kernel for addition of data
  !! \f$ v(dg(i)) = v(dg(i)) + u(gd(i)) \f$
  subroutine gs_gather_kernel_add(v, m, o, dg, u, n, gd, nb, b, w)
    integer, intent(in) :: m
    integer, intent(in) :: n
    integer, intent(in) :: nb
    real(kind=rp), dimension(m), intent(inout) :: v
    real(kind=rp), dimension(m), intent(inout) :: w
    integer, dimension(m), intent(inout) :: dg
    real(kind=rp), dimension(n), intent(inout) :: u
    integer, dimension(m), intent(inout) :: gd
    integer, dimension(nb), intent(inout) :: b
    integer, intent(in) :: o
    integer :: i
    real(kind=rp) :: tmp

    v = 0d0 
    do i = 1, abs(o) - 1
       w(i) = u(gd(i))
    end do

    do i = 1, abs(o) - 1
       v(dg(i)) = v(dg(i)) + w(i)
    end do
    
    if (o .lt. 0) then
       do i = abs(o), m
          v(dg(i)) = u(gd(i))
       end do
    else
       do i = o, m, 2
          tmp  = u(gd(i)) + u(gd(i+1))
          v(dg(i)) = tmp
       end do
    end if
    
  end subroutine gs_gather_kernel_add

  !> Gather kernel for multiplication of data
  !! \f$ v(dg(i)) = v(dg(i)) \cdot u(gd(i)) \f$
  subroutine gs_gather_kernel_mul(v, m, o, dg, u, n, gd, nb, b, w)
    integer, intent(in) :: m
    integer, intent(in) :: n
    integer, intent(in) :: nb
    real(kind=rp), dimension(m), intent(inout) :: v
    real(kind=rp), dimension(m), intent(inout) :: w
    integer, dimension(m), intent(inout) :: dg
    real(kind=rp), dimension(n), intent(inout) :: u
    integer, dimension(m), intent(inout) :: gd
    integer, dimension(nb), intent(inout) :: b
    integer, intent(in) :: o
    integer :: i
    real(kind=rp) :: tmp
    
    do i = 1, abs(o) - 1
       w(i) = u(gd(i))
    end do

    do i = 1, abs(o) - 1
       v(dg(i)) = v(dg(i)) * w(i)
    end do
    
    if (o .lt. 0) then
       do i = abs(o), m
          v(dg(i)) = u(gd(i))
       end do
    else
       do i = o, m, 2
          tmp  = u(gd(i)) * u(gd(i+1))
          v(dg(i)) = tmp
       end do
    end if
    
  end subroutine gs_gather_kernel_mul
  
  !> Gather kernel for minimum of data
  !! \f$ v(dg(i)) = \min(v(dg(i)), u(gd(i))) \f$
  subroutine gs_gather_kernel_min(v, m, o, dg, u, n, gd, nb, b, w)
    integer, intent(in) :: m
    integer, intent(in) :: n
    integer, intent(in) :: nb
    real(kind=rp), dimension(m), intent(inout) :: v
    real(kind=rp), dimension(m), intent(inout) :: w
    integer, dimension(m), intent(inout) :: dg
    real(kind=rp), dimension(n), intent(inout) :: u
    integer, dimension(m), intent(inout) :: gd
    integer, dimension(nb), intent(inout) :: b
    integer, intent(in) :: o
    integer :: i
    real(kind=rp) :: tmp

    do i = 1, abs(o) - 1
       w(i) = u(gd(i))
    end do

    do i = 1, abs(o) - 1
       v(dg(i)) = min(v(dg(i)), w(i))
    end do
    
    if (o .lt. 0) then
       do i = abs(o), m
          v(dg(i)) = u(gd(i))
       end do
    else
       do i = o, m, 2
          tmp  = min(u(gd(i)), u(gd(i+1)))
          v(dg(i)) = tmp
       end do
    end if
    
  end subroutine gs_gather_kernel_min

  !> Gather kernel for maximum of data
  !! \f$ v(dg(i)) = \max(v(dg(i)), u(gd(i))) \f$
  subroutine gs_gather_kernel_max(v, m, o, dg, u, n, gd, nb, b, w)
    integer, intent(in) :: m
    integer, intent(in) :: n
    integer, intent(in) :: nb
    real(kind=rp), dimension(m), intent(inout) :: v
    real(kind=rp), dimension(m), intent(inout) :: w
    integer, dimension(m), intent(inout) :: dg
    real(kind=rp), dimension(n), intent(inout) :: u
    integer, dimension(m), intent(inout) :: gd
    integer, dimension(nb), intent(inout) :: b
    integer, intent(in) :: o
    integer :: i
    real(kind=rp) :: tmp

    do i = 1, abs(o) - 1
       w(i) = u(gd(i))
    end do

    do i = 1, abs(o) - 1
       v(dg(i)) = max(v(dg(i)), w(i))
    end do
    
    if (o .lt. 0) then
       do i = abs(o), m
          v(dg(i)) = u(gd(i))
       end do
    else
       do i = o, m, 2
          tmp  = max(u(gd(i)), u(gd(i+1)))
          v(dg(i)) = tmp
       end do
    end if
    
  end subroutine gs_gather_kernel_max

  !> Scatter kernel  @todo Make the kernel abstract
  subroutine gs_scatter_sx(this, v, m, dg, u, n, gd, nb, b, shrd)
    integer, intent(in) :: m
    integer, intent(in) :: n
    integer, intent(in) :: nb
    class(gs_sx_t), intent(inout) :: this
    real(kind=rp), dimension(m), intent(inout) :: v
    integer, dimension(m), intent(inout) :: dg
    real(kind=rp), dimension(n), intent(inout) :: u
    integer, dimension(m), intent(inout) :: gd
    integer, dimension(nb), intent(inout) :: b
    logical, intent(in) :: shrd
        
    if (.not. shrd) then
       call gs_scatter_kernel(v, m, dg, u, n, gd, nb, b, this%local_wrk)
    else if (shrd) then
       call gs_scatter_kernel(v, m, dg, u, n, gd, nb, b, this%shared_wrk)
    end if

  end subroutine gs_scatter_sx

  !> Scatter kernel \f$ u(gd(i) = v(dg(i)) \f$
  subroutine gs_scatter_kernel(v, m, dg, u, n, gd, nb, b, w)
    integer, intent(in) :: m
    integer, intent(in) :: n
    integer, intent(in) :: nb
    real(kind=rp), dimension(m), intent(inout) :: v
    real(kind=rp), dimension(m), intent(inout) :: w
    integer, dimension(m), intent(inout) :: dg
    real(kind=rp), dimension(n), intent(inout) :: u
    integer, dimension(m), intent(inout) :: gd
    integer, dimension(nb), intent(inout) :: b
    integer :: i
    
    !NEC$ IVDEP
    do i = 1, m
       w(i) = v(dg(i))
    end do

    !NEC$ IVDEP
    do i = 1, m
       u(gd(i)) = w(i)
    end do
    
  end subroutine gs_scatter_kernel

end module gs_sx
