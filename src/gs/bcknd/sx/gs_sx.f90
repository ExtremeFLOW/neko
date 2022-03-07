! Copyright (c) 2020-2022, The Neko Authors
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
     real(kind=rp), allocatable :: local_wrk1(:)
     real(kind=rp), allocatable :: local_wrk2(:)
     real(kind=rp), allocatable :: local_wrk3(:)
     real(kind=rp), allocatable :: shared_wrk1(:)
     real(kind=rp), allocatable :: shared_wrk2(:)
     real(kind=rp), allocatable :: shared_wrk3(:)
     integer :: nlocal
     integer :: nshared
   contains
     procedure, pass(this) :: init => gs_sx_init
     procedure, pass(this) :: free => gs_sx_free
     procedure, pass(this) :: gather => gs_gather_sx
     procedure, pass(this) :: scatter => gs_scatter_sx
     procedure, pass(this) :: gather_many => gs_gather_many_sx
     procedure, pass(this) :: scatter_many => gs_scatter_many_sx
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

    allocate(this%local_wrk1(nlocal))
    allocate(this%local_wrk2(nlocal))
    allocate(this%local_wrk3(nlocal))
    
    allocate(this%shared_wrk1(nshared))
    allocate(this%shared_wrk2(nshared))    
    allocate(this%shared_wrk3(nshared))
    
  end subroutine gs_sx_init

  !> SX backend deallocation
  subroutine gs_sx_free(this)
    class(gs_sx_t), intent(inout) :: this

    if (allocated(this%local_wrk1)) then
       deallocate(this%local_wrk1)
    end if

    if (allocated(this%local_wrk2)) then
       deallocate(this%local_wrk2)
    end if

    if (allocated(this%local_wrk3)) then
       deallocate(this%local_wrk3)
    end if

    if (allocated(this%shared_wrk1)) then
       deallocate(this%shared_wrk1)
    end if

    if (allocated(this%shared_wrk2)) then
       deallocate(this%shared_wrk2)
    end if

    if (allocated(this%shared_wrk3)) then
       deallocate(this%shared_wrk3)
    end if

    this%nlocal = 0
    this%nshared = 0
    
  end subroutine gs_sx_free

  !> Gather kernel
  subroutine gs_gather_sx(this, v, m, o, dg, u, n, gd, nb, b, op)
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
    integer :: op

    if (this%nlocal .eq. m) then
       associate(w=>this%local_wrk1)
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
    else if (this%nshared .eq. m) then
       associate(w=>this%shared_wrk1)
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
    integer :: i, j, k, blk_len
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
    integer :: i, j, k, blk_len
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
    integer :: i, j, k, blk_len
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
    integer :: i, j, k, blk_len
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
  subroutine gs_scatter_sx(this, v, m, dg, u, n, gd, nb, b)
    integer, intent(in) :: m
    integer, intent(in) :: n
    integer, intent(in) :: nb
    class(gs_sx_t), intent(inout) :: this
    real(kind=rp), dimension(m), intent(inout) :: v
    integer, dimension(m), intent(inout) :: dg
    real(kind=rp), dimension(n), intent(inout) :: u
    integer, dimension(m), intent(inout) :: gd
    integer, dimension(nb), intent(inout) :: b
        
    if (this%nlocal .eq. m) then
       call gs_scatter_kernel(v, m, dg, u, n, gd, nb, b, this%local_wrk1)
    else if (this%nshared .eq. m) then
       call gs_scatter_kernel(v, m, dg, u, n, gd, nb, b, this%shared_wrk1)
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
    integer :: i, j, k, blk_len
    real(kind=rp) :: tmp
    
    !NEC$ IVDEP
    do i = 1, m
       w(i) = v(dg(i))
    end do

    !NEC$ IVDEP
    do i = 1, m
       u(gd(i)) = w(i)
    end do
    
  end subroutine gs_scatter_kernel

  !> Gather kernel many
  subroutine gs_gather_many_sx(this, v1, v2, v3, m, o, dg, &
                                     u1, u2, u3, n, gd, nb, b, op)
    integer, intent(inout) :: m
    integer, intent(inout) :: n
    integer, intent(inout) :: nb
    class(gs_sx_t), intent(inout) :: this
    real(kind=rp), dimension(m), intent(inout) :: v1
    real(kind=rp), dimension(m), intent(inout) :: v2
    real(kind=rp), dimension(m), intent(inout) :: v3
    integer, dimension(m), intent(inout) :: dg
    real(kind=rp), dimension(n), intent(inout) :: u1
    real(kind=rp), dimension(n), intent(inout) :: u2
    real(kind=rp), dimension(n), intent(inout) :: u3
    integer, dimension(m), intent(inout) :: gd
    integer, dimension(nb), intent(inout) :: b
    integer, intent(inout) :: o
    integer :: op

    if (this%nlocal .eq. m) then
       associate(w1=>this%local_wrk1, w2=>this%local_wrk2, w3=>this%local_wrk3)
         select case(op)
         case (GS_OP_ADD)
            call gs_gather_many_kernel_add(v1, v2, v3, m, o, dg, &
                                           u1, u2, u3, n, gd, nb, b, w1, w2, w3)
         case (GS_OP_MUL)
            call gs_gather_many_kernel_mul(v1, v2, v3, m, o, dg, &
                                           u1, u2, u3, n, gd, nb, b, w1, w2, w3)
         case (GS_OP_MIN)
            call gs_gather_many_kernel_min(v1, v2, v3, m, o, dg, &
                                           u1, u2, u3, n, gd, nb, b, w1, w2, w3)
         case (GS_OP_MAX)
            call gs_gather_many_kernel_max(v1, v2, v3, m, o, dg, &
                                           u1, u2, u3, n, gd, nb, b, w1, w2, w3)
         end select
       end associate
    else if (this%nshared .eq. m) then
       associate(w1=>this%shared_wrk1, w2=>this%shared_wrk2,&
                 w3=>this%shared_wrk3)
         select case(op)
         case (GS_OP_ADD)
            call gs_gather_many_kernel_add(v1, v2, v3,  m, o, dg, &
                                           u1, u2, u3, n, gd, nb, b, w1, w2, w3)
         case (GS_OP_MUL)
            call gs_gather_many_kernel_mul(v1, v2, v3, m, o, dg, &
                                           u1, u2, u3, n, gd, nb, b, w1, w2, w3)
         case (GS_OP_MIN)
            call gs_gather_many_kernel_min(v1, v2, v3, m, o, dg, &
                                           u1, u2, u3, n, gd, nb, b, w1, w2, w3)
         case (GS_OP_MAX)
            call gs_gather_many_kernel_max(v1, v2, v3, m, o, dg, &
                                           u1, u2, u3, n, gd, nb, b, w1, w2, w3)
         end select
       end associate
    end if
    
  end subroutine gs_gather_many_sx
 
  !> Gather kernel for addition of data (many version)
  !! \f$ v_n(dg(i)) = v_n(dg(i)) + u_n(gd(i)) \f$
  subroutine gs_gather_many_kernel_add(v1, v2, v3, m, o, dg, &
                                       u1, u2, u3, n, gd, nb, b, w1, w2, w3)
    integer, intent(in) :: m
    integer, intent(in) :: n
    integer, intent(in) :: nb
    real(kind=rp), dimension(m), intent(inout) :: v1
    real(kind=rp), dimension(m), intent(inout) :: v2
    real(kind=rp), dimension(m), intent(inout) :: v3
    real(kind=rp), dimension(m), intent(inout) :: w1
    real(kind=rp), dimension(m), intent(inout) :: w2
    real(kind=rp), dimension(m), intent(inout) :: w3
    integer, dimension(m), intent(inout) :: dg
    real(kind=rp), dimension(n), intent(inout) :: u1
    real(kind=rp), dimension(n), intent(inout) :: u2
    real(kind=rp), dimension(n), intent(inout) :: u3
    integer, dimension(m), intent(inout) :: gd
    integer, dimension(nb), intent(inout) :: b
    integer, intent(in) :: o
    integer :: i, j, k, blk_len
    real(kind=rp) :: tmp1, tmp2, tmp3

    v1 = 0d0
    v2 = 0d0
    v3 = 0d0 
    do i = 1, abs(o) - 1
       w1(i) = u1(gd(i))
       w2(i) = u2(gd(i))
       w3(i) = u3(gd(i))
    end do

    do i = 1, abs(o) - 1
       v1(dg(i)) = v1(dg(i)) + w1(i)
       v2(dg(i)) = v2(dg(i)) + w2(i)
       v3(dg(i)) = v3(dg(i)) + w3(i)
    end do
    
    if (o .lt. 0) then
       do i = abs(o), m
          v1(dg(i)) = u1(gd(i))
          v2(dg(i)) = u2(gd(i))
          v3(dg(i)) = u3(gd(i))
       end do
    else
       do i = o, m, 2
          tmp1  = u1(gd(i)) + u1(gd(i+1))
          tmp2  = u2(gd(i)) + u2(gd(i+1))
          tmp3  = u3(gd(i)) + u3(gd(i+1))
          v1(dg(i)) = tmp1
          v2(dg(i)) = tmp2
          v3(dg(i)) = tmp3
       end do
    end if
    
  end subroutine gs_gather_many_kernel_add

  !> Gather kernel for multiplication of data (many version)
  !! \f$ v_n(dg(i)) = v_n(dg(i)) \cdot u_n(gd(i)) \f$
  subroutine gs_gather_many_kernel_mul(v1, v2, v3, m, o, dg, &
                                       u1, u2, u3, n, gd, nb, b, w1, w2, w3)
    integer, intent(in) :: m
    integer, intent(in) :: n
    integer, intent(in) :: nb
    real(kind=rp), dimension(m), intent(inout) :: v1
    real(kind=rp), dimension(m), intent(inout) :: v2
    real(kind=rp), dimension(m), intent(inout) :: v3
    real(kind=rp), dimension(m), intent(inout) :: w1
    real(kind=rp), dimension(m), intent(inout) :: w2
    real(kind=rp), dimension(m), intent(inout) :: w3
    integer, dimension(m), intent(inout) :: dg
    real(kind=rp), dimension(n), intent(inout) :: u1
    real(kind=rp), dimension(n), intent(inout) :: u2
    real(kind=rp), dimension(n), intent(inout) :: u3
    integer, dimension(m), intent(inout) :: gd
    integer, dimension(nb), intent(inout) :: b
    integer, intent(in) :: o
    integer :: i, j, k, blk_len
    real(kind=rp) :: tmp1, tmp2, tmp3
    
    do i = 1, abs(o) - 1
       w1(i) = u1(gd(i))
       w2(i) = u2(gd(i))
       w3(i) = u3(gd(i))
    end do

    do i = 1, abs(o) - 1
       v1(dg(i)) = v1(dg(i)) * w1(i)
       v2(dg(i)) = v2(dg(i)) * w2(i)
       v3(dg(i)) = v3(dg(i)) * w3(i)
    end do
    
    if (o .lt. 0) then
       do i = abs(o), m
          v1(dg(i)) = u1(gd(i))
          v2(dg(i)) = u2(gd(i))
          v3(dg(i)) = u3(gd(i))
       end do
    else
       do i = o, m, 2
          tmp1  = u1(gd(i)) * u1(gd(i+1))
          tmp2  = u2(gd(i)) * u2(gd(i+1))
          tmp3  = u3(gd(i)) * u3(gd(i+1))
          v1(dg(i)) = tmp1
          v2(dg(i)) = tmp2
          v3(dg(i)) = tmp3
       end do
    end if
    
  end subroutine gs_gather_many_kernel_mul
  
  !> Gather kernel for minimum of data (many version)
  !! \f$ v_n(dg(i)) = \min(v_n(dg(i)), u_n(gd(i))) \f$
  subroutine gs_gather_many_kernel_min(v1, v2, v3, m, o, dg, &
                                       u1, u2, u3, n, gd, nb, b, w1, w2, w3)
    integer, intent(in) :: m
    integer, intent(in) :: n
    integer, intent(in) :: nb
    real(kind=rp), dimension(m), intent(inout) :: v1
    real(kind=rp), dimension(m), intent(inout) :: v2
    real(kind=rp), dimension(m), intent(inout) :: v3
    real(kind=rp), dimension(m), intent(inout) :: w1
    real(kind=rp), dimension(m), intent(inout) :: w2
    real(kind=rp), dimension(m), intent(inout) :: w3
    integer, dimension(m), intent(inout) :: dg
    real(kind=rp), dimension(n), intent(inout) :: u1
    real(kind=rp), dimension(n), intent(inout) :: u2
    real(kind=rp), dimension(n), intent(inout) :: u3
    integer, dimension(m), intent(inout) :: gd
    integer, dimension(nb), intent(inout) :: b
    integer, intent(in) :: o
    integer :: i, j, k, blk_len
    real(kind=rp) :: tmp1, tmp2, tmp3

    do i = 1, abs(o) - 1
       w1(i) = u1(gd(i))
       w2(i) = u2(gd(i))
       w3(i) = u3(gd(i))
    end do

    do i = 1, abs(o) - 1
       v1(dg(i)) = min(v1(dg(i)), w1(i))
       v2(dg(i)) = min(v2(dg(i)), w2(i))
       v3(dg(i)) = min(v3(dg(i)), w3(i))
    end do
    
    if (o .lt. 0) then
       do i = abs(o), m
          v1(dg(i)) = u1(gd(i))
          v2(dg(i)) = u2(gd(i))
          v3(dg(i)) = u3(gd(i))
       end do
    else
       do i = o, m, 2
          tmp1  = min(u1(gd(i)), u1(gd(i+1)))
          tmp2  = min(u2(gd(i)), u2(gd(i+1)))
          tmp3  = min(u3(gd(i)), u3(gd(i+1)))
          v1(dg(i)) = tmp1
          v2(dg(i)) = tmp2
          v3(dg(i)) = tmp3
       end do
    end if
    
  end subroutine gs_gather_many_kernel_min

  !> Gather kernel for maximum of data (many version)
  !! \f$ v_n(dg(i)) = \max(v_n(dg(i)), u_n(gd(i))) \f$
  subroutine gs_gather_many_kernel_max(v1, v2, v3, m, o, dg, &
                                       u1, u2, u3, n, gd, nb, b, w1, w2, w3)
    integer, intent(in) :: m
    integer, intent(in) :: n
    integer, intent(in) :: nb
    real(kind=rp), dimension(m), intent(inout) :: v1
    real(kind=rp), dimension(m), intent(inout) :: v2
    real(kind=rp), dimension(m), intent(inout) :: v3
    real(kind=rp), dimension(m), intent(inout) :: w1
    real(kind=rp), dimension(m), intent(inout) :: w2
    real(kind=rp), dimension(m), intent(inout) :: w3
    integer, dimension(m), intent(inout) :: dg
    real(kind=rp), dimension(n), intent(inout) :: u1
    real(kind=rp), dimension(n), intent(inout) :: u2
    real(kind=rp), dimension(n), intent(inout) :: u3
    integer, dimension(m), intent(inout) :: gd
    integer, dimension(nb), intent(inout) :: b
    integer, intent(in) :: o
    integer :: i, j, k, blk_len
    real(kind=rp) :: tmp1, tmp2, tmp3

    do i = 1, abs(o) - 1
       w1(i) = u1(gd(i))
       w2(i) = u2(gd(i))
       w3(i) = u3(gd(i))
    end do

    do i = 1, abs(o) - 1
       v1(dg(i)) = max(v1(dg(i)), w1(i))
       v2(dg(i)) = max(v2(dg(i)), w2(i))
       v3(dg(i)) = max(v3(dg(i)), w3(i))
    end do
    
    if (o .lt. 0) then
       do i = abs(o), m
          v1(dg(i)) = u1(gd(i))
          v2(dg(i)) = u2(gd(i))
          v3(dg(i)) = u3(gd(i))
       end do
    else
       do i = o, m, 2
          tmp1 = max(u1(gd(i)), u1(gd(i+1)))
          tmp2 = max(u2(gd(i)), u2(gd(i+1)))
          tmp3 = max(u3(gd(i)), u3(gd(i+1)))
          v1(dg(i)) = tmp1
          v2(dg(i)) = tmp2
          v3(dg(i)) = tmp3
       end do
    end if
    
  end subroutine gs_gather_many_kernel_max

  !> Scatter kernel many @todo Make the kernel abstract 
  subroutine gs_scatter_many_sx(this, v1, v2, v3, m, dg, &
                                      u1, u2, u3, n, gd, nb, b)
    integer, intent(in) :: m
    integer, intent(in) :: n
    integer, intent(in) :: nb
    class(gs_sx_t), intent(inout) :: this
    real(kind=rp), dimension(m), intent(inout) :: v1
    real(kind=rp), dimension(m), intent(inout) :: v2
    real(kind=rp), dimension(m), intent(inout) :: v3
    integer, dimension(m), intent(inout) :: dg
    real(kind=rp), dimension(n), intent(inout) :: u1
    real(kind=rp), dimension(n), intent(inout) :: u2
    real(kind=rp), dimension(n), intent(inout) :: u3
    integer, dimension(m), intent(inout) :: gd
    integer, dimension(nb), intent(inout) :: b
        
    if (this%nlocal .eq. m) then
       associate(w1=>this%local_wrk1, w2=>this%local_wrk2, w3=>this%local_wrk3)
         call gs_scatter_many_kernel(v1, v2, v3, m, dg, &
                                     u1, u2, u3, n, gd, nb, b, w1, w2, w3)
       end associate
    else if (this%nshared .eq. m) then
       associate(w1=>this%shared_wrk1, w2=>this%shared_wrk2, &
                 w3=>this%shared_wrk3)
         call gs_scatter_many_kernel(v1, v2, v3, m, dg, &
                                     u1, u2, u3, n, gd, nb, b, w1, w2, w3)
       end associate
    end if

  end subroutine gs_scatter_many_sx

  !> Scatter many kernel \f$ u_n(gd(i) = v_n(dg(i)) \f$ 
  subroutine gs_scatter_many_kernel(v1, v2, v3, m, dg, &
                                    u1, u2, u3, n, gd, nb, b, w1, w2, w3)
    integer, intent(in) :: m
    integer, intent(in) :: n
    integer, intent(in) :: nb
    real(kind=rp), dimension(m), intent(inout) :: v1
    real(kind=rp), dimension(m), intent(inout) :: v2
    real(kind=rp), dimension(m), intent(inout) :: v3
    real(kind=rp), dimension(m), intent(inout) :: w1
    real(kind=rp), dimension(m), intent(inout) :: w2
    real(kind=rp), dimension(m), intent(inout) :: w3
    integer, dimension(m), intent(inout) :: dg
    real(kind=rp), dimension(n), intent(inout) :: u1
    real(kind=rp), dimension(n), intent(inout) :: u2
    real(kind=rp), dimension(n), intent(inout) :: u3
    integer, dimension(m), intent(inout) :: gd
    integer, dimension(nb), intent(inout) :: b
    integer :: i, j, k, blk_len
    real(kind=rp) :: tmp
    
    !NEC$ IVDEP
    do i = 1, m
       w1(i) = v1(dg(i))
       w2(i) = v2(dg(i))
       w3(i) = v3(dg(i))
    end do

    !NEC$ IVDEP
    do i = 1, m
       u1(gd(i)) = w1(i)
       u2(gd(i)) = w2(i)
       u3(gd(i)) = w3(i)
    end do
    
  end subroutine gs_scatter_many_kernel
  
end module gs_sx
