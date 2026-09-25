! Copyright (c) 2025-2026, The Neko Authors
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
!> Fortran bindings to SHMEM's C API
module shmem
  use, intrinsic :: iso_c_binding
  implicit none

#ifdef HAVE_OPENSHMEM

  enum, bind(c)
    enumerator :: SHMEM_MAX_NAME_LEN = 256
  end enum

  enum, bind(c)
    enumerator :: SHMEM_CMP_EQ = 0
    enumerator :: SHMEM_CMP_NE = 1
    enumerator :: SHMEM_CMP_GT = 2
    enumerator :: SHMEM_CMP_LE = 3
    enumerator :: SHMEM_CMP_LT = 4
    enumerator :: SHMEM_CMP_GE = 5
  end enum

  enum, bind(c)
    enumerator :: SHMEM_THREAD_SINGLE = 0
    enumerator :: SHMEM_THREAD_FUNNELED = 1
    enumerator :: SHMEM_THREAD_SERIALIZED = 2
    enumerator :: SHMEM_THREAD_MULTIPLE= 3
  end enum

  enum, bind(c)
    enumerator :: SHMEM_CTX_PRIVATE = 1
    enumerator :: SHMEM_CTX_SERIALIZED = 2
    enumerator :: SHMEM_CTX_NOSTORE = 4
  end enum

  enum, bind(c)
    enumerator :: SHMEM_CTX_LOW_LATENCY = int(Z'100')
    enumerator :: SHMEM_CTX_DEDICATED = int(Z'200')
    enumerator :: SHMEM_CTX_BEST_EFFORT = int(Z'400')
    enumerator :: SHMEM_CTX_BULK_DATA = int(Z'800')
  end enum

  enum, bind(c)
    enumerator :: SHMEM_MALLOC_ATOMICS_REMOTE = 1
    enumerator :: SHMEM_MALLOC_SIGNAL_REMOTE = 2
  end enum

  enum, bind(c)
    enumerator :: SHMEM_SIGNAL_SET = 1
    enumerator :: SHMEM_SIGNAL_ADD = 2
  end enum

  !
  ! Library Setup
  !

  interface
     subroutine shmem_init() bind (c, name='shmem_init')
     end subroutine shmem_init
  end interface

  interface
     subroutine shmem_finalize() bind (c, name='shmem_finalize')
     end subroutine shmem_finalize
  end interface

  interface
     integer(c_int) function shmem_me_pe() bind(c, name='shmem_my_pe')
       use, intrinsic :: iso_c_binding
     end function shmem_me_pe
  end interface

  interface
     integer(c_int) function shmem_n_pes() bind(c, name='shmem_n_pes')
       use, intrinsic :: iso_c_binding
     end function shmem_n_pes
  end interface

  !
  ! Library Query
  !

  interface
     subroutine shmem_info_get_name(name) &
          bind(c, name = 'shmem_info_get_name')
       use, intrinsic :: iso_c_binding
       import SHMEM_MAX_NAME_LEN
       implicit none
       character(kind=c_char, len=SHMEM_MAX_NAME_LEN) :: name
     end subroutine shmem_info_get_name
  end interface

  interface
     subroutine shmem_info_get_version(major, minor) &
          bind(c, name = 'shmem_info_get_version')
       use, intrinsic :: iso_c_binding
       integer(c_int) :: major
       integer(c_int) :: minor
     end subroutine shmem_info_get_version
  end interface

  interface
     integer(c_int) function shmem_pe_accessible(pe) &
          bind(c, name = 'shmem_pe_accessible')
       use, intrinsic :: iso_c_binding
       integer(c_int), value :: pe
     end function shmem_pe_accessible
  end interface

  interface
     integer(c_int) function shmem_addr_accessible(addr, pe) &
          bind(c, name = 'shmem_addr_accessible')
       use, intrinsic :: iso_c_binding
       type(c_ptr) :: addr
       integer(c_int), value :: pe
     end function shmem_addr_accessible
  end interface

  !
  ! Memory Mangement
  !

  interface
     type(c_ptr) function shmem_malloc(size) &
          bind(c, name = 'shmem_malloc')
       use, intrinsic :: iso_c_binding
       integer(c_size_t), value :: size
     end function shmem_malloc
  end interface

  interface
     subroutine shmem_free(ptr) &
          bind(c, name = 'shmem_free')
       use, intrinsic :: iso_c_binding
       type(c_ptr), value :: ptr
     end subroutine shmem_free
  end interface

  interface
     type(c_ptr) function shmem_realloc(ptr, size) &
          bind(c, name = 'shmem_realloc')
       use, intrinsic :: iso_c_binding
       type(c_ptr), value :: ptr
       integer(c_size_t), value :: size
     end function shmem_realloc
  end interface

  interface
     type(c_ptr) function shmem_align(alignment, size) &
          bind(c, name = 'shmem_align')
       use, intrinsic :: iso_c_binding
       integer(c_size_t), value :: alignment
       integer(c_size_t), value :: size
     end function shmem_align
  end interface

  interface
     type(c_ptr) function shmem_calloc(count, size) &
          bind(c, name = 'shmem_calloc')
       use, intrinsic :: iso_c_binding
       integer(c_size_t), value :: count
       integer(c_size_t), value :: size
     end function shmem_calloc
  end interface

  interface
     type(c_ptr) function shmem_malloc_with_hints(size, hints) &
          bind(c, name = 'shmem_malloc_with_hints')
       use, intrinsic :: iso_c_binding
       integer(c_size_t), value :: size
       integer(c_long), value :: hints
     end function shmem_malloc_with_hints
  end interface

  interface
     type(c_ptr) function shmem_ptr(dest, pe) &
          bind(c, name = 'shmem_ptr')
       use, intrinsic :: iso_c_binding
       type(c_ptr), value :: dest
       integer(c_int), value :: pe
     end function shmem_ptr
  end interface

  !
  ! Thread Support
  !

  interface
     integer(c_int) function shmem_init_thread(requested, provided) &
          bind(c, name = 'shmem_init_thread')
       use, intrinsic :: iso_c_binding
       integer(c_int), value :: requested
       integer(c_int) :: provided
     end function shmem_init_thread
  end interface

  interface
     subroutine shmem_query_thread(provided) &
          bind(c, name = 'shmem_query_thread')
       use, intrinsic :: iso_c_binding
       integer(c_int) :: provided
     end subroutine shmem_query_thread
  end interface

  !
  ! Communication Management
  !

  interface
     integer(c_int) function shmem_ctx_create(options, ctx) &
          bind(c, name = 'shmem_ctx_create')
       use, intrinsic :: iso_c_binding
       integer(c_long), value :: options
       type(c_ptr), value :: ctx
     end function shmem_ctx_create
  end interface

  interface
     subroutine shmem_ctx_destroy(ctx) &
          bind(c, name = 'shmem_ctx_destroy')
       use, intrinsic :: iso_c_binding
       type(c_ptr), value :: ctx
     end subroutine shmem_ctx_destroy
  end interface

  !
  ! Remote Memory Access
  !

  interface
     subroutine shmem_putmem(dest, source, nelems, pe) &
          bind(c, name = 'shmem_putmem')
       use, intrinsic :: iso_c_binding
       type(c_ptr), value :: dest
       type(c_ptr), value :: source
       integer(c_size_t), value :: nelems
       integer(c_int), value :: pe
     end subroutine shmem_putmem
  end interface

  interface
     subroutine shmem_ctx_putmem(ctx, dest, source, nelems, pe) &
          bind(c, name = 'shmem_ctx_putmem')
       use, intrinsic :: iso_c_binding
       type(c_ptr), value :: ctx
       type(c_ptr), value :: dest
       type(c_ptr), value :: source
       integer(c_size_t), value :: nelems
       integer(c_int), value :: pe
     end subroutine shmem_ctx_putmem
  end interface

  interface
     subroutine shmem_getmem(dest, source, nelems, pe) &
          bind(c, name = 'shmem_getmem')
       use, intrinsic :: iso_c_binding
       type(c_ptr), value :: dest
       type(c_ptr), value :: source
       integer(c_size_t), value :: nelems
       integer(c_int), value :: pe
     end subroutine shmem_getmem
  end interface

  interface
     subroutine shmem_ctx_getmem(ctx, dest, source, nelems, pe) &
          bind(c, name = 'shmem_ctx_getmem')
       use, intrinsic :: iso_c_binding
       type(c_ptr), value :: ctx
       type(c_ptr), value :: dest
       type(c_ptr), value :: source
       integer(c_size_t), value :: nelems
       integer(c_int), value :: pe
     end subroutine shmem_ctx_getmem
  end interface

  interface
     subroutine shmem_putmem_nbi(dest, source, nelems, pe) &
          bind(c, name = 'shmem_putmem_nbi')
       use, intrinsic :: iso_c_binding
       type(c_ptr), value :: dest
       type(c_ptr), value :: source
       integer(c_size_t), value :: nelems
       integer(c_int), value :: pe
     end subroutine shmem_putmem_nbi
  end interface

  interface
     subroutine shmem_ctx_putmem_nbi(ctx, dest, source, nelems, pe) &
          bind(c, name = 'shmem_ctx_putmem_nbi')
       use, intrinsic :: iso_c_binding
       type(c_ptr), value :: ctx
       type(c_ptr), value :: dest
       type(c_ptr), value :: source
       integer(c_size_t), value :: nelems
       integer(c_int), value :: pe
     end subroutine shmem_ctx_putmem_nbi
  end interface

  interface
     subroutine shmem_getmem_nbi(dest, source, nelems, pe) &
          bind(c, name = 'shmem_getmem_nbi')
       use, intrinsic :: iso_c_binding
       type(c_ptr), value :: dest
       type(c_ptr), value :: source
       integer(c_size_t), value :: nelems
       integer(c_int), value :: pe
     end subroutine shmem_getmem_nbi
  end interface

  interface
     subroutine shmem_ctx_getmem_nbi(ctx, dest, source, nelems, pe) &
          bind(c, name = 'shmem_ctx_getmem_nbi')
       use, intrinsic :: iso_c_binding
       type(c_ptr), value :: ctx
       type(c_ptr), value :: dest
       type(c_ptr), value :: source
       integer(c_size_t), value :: nelems
       integer(c_int), value :: pe
     end subroutine shmem_ctx_getmem_nbi
  end interface

  !
  ! Signaling Operations (OpenSHMEM 1.5)
  !

  interface
     subroutine shmem_putmem_signal_nbi(dest, source, nelems, sig_addr, &
          signal, sig_op, pe) bind(c, name = 'shmem_putmem_signal_nbi')
       use, intrinsic :: iso_c_binding
       type(c_ptr), value :: dest
       type(c_ptr), value :: source
       integer(c_size_t), value :: nelems
       type(c_ptr), value :: sig_addr
       integer(c_int64_t), value :: signal
       integer(c_int), value :: sig_op
       integer(c_int), value :: pe
     end subroutine shmem_putmem_signal_nbi
  end interface

  interface
     integer(c_int64_t) function shmem_signal_wait_until(sig_addr, cmp, &
          cmp_value) bind(c, name = 'shmem_signal_wait_until')
       use, intrinsic :: iso_c_binding
       type(c_ptr), value :: sig_addr
       integer(c_int), value :: cmp
       integer(c_int64_t), value :: cmp_value
     end function shmem_signal_wait_until
  end interface

  interface
     subroutine shmem_uint64_atomic_set(dest, val, pe) &
          bind(c, name = 'shmem_uint64_atomic_set')
       use, intrinsic :: iso_c_binding
       type(c_ptr), value :: dest
       integer(c_int64_t), value :: val
       integer(c_int), value :: pe
     end subroutine shmem_uint64_atomic_set
  end interface

  interface
     subroutine shmem_uint64_wait_until(ivar, cmp, cmp_value) &
          bind(c, name = 'shmem_uint64_wait_until')
       use, intrinsic :: iso_c_binding
       type(c_ptr), value :: ivar
       integer(c_int), value :: cmp
       integer(c_int64_t), value :: cmp_value
     end subroutine shmem_uint64_wait_until
  end interface

  !
  ! Collective Operations
  !

  interface
     integer(c_int) function shmem_alltoallmem(team, dest, source, nelems) &
          bind(c, name = 'shmem_alltoall_mem')
       use, intrinsic :: iso_c_binding
       type(c_ptr), value :: team
       type(c_ptr), value :: dest
       type(c_ptr), value :: source
       integer(c_size_t), value :: nelems
     end function shmem_alltoallmem
  end interface

  interface
     subroutine shmem_barrier(PE_start, logPE_stride, PE_size, pSync) &
          bind(c, name = 'shmem_barrier')
       use, intrinsic :: iso_c_binding
       integer(c_int), value :: PE_start
       integer(c_int), value :: logPE_stride
       integer(c_int), value :: PE_size
       integer(c_long) :: pSync
     end subroutine shmem_barrier
  end interface

  interface
     subroutine shmem_barrier_all() &
          bind(c, name = 'shmem_barrier_all')
     end subroutine shmem_barrier_all
  end interface

  interface
     integer(c_int) function shmem_broadcastmem(team, dest, source, &
          nelems, PE_root) bind(c, name = 'shmem_broadcastmem')
       use, intrinsic :: iso_c_binding
       type(c_ptr), value :: team
       type(c_ptr), value :: dest
       type(c_ptr), value :: source
       integer(c_size_t), value :: nelems
       integer(c_int), value :: PE_root
     end function shmem_broadcastmem
  end interface

  interface
     integer(c_int) function shmem_collectmem(team, dest, source, nelems) &
          bind(c, name = 'shmem_collectmem')
       use, intrinsic :: iso_c_binding
       type(c_ptr), value :: team
       type(c_ptr), value :: dest
       type(c_ptr), value :: source
       integer(c_size_t), value :: nelems
     end function shmem_collectmem
  end interface

  interface
     integer(c_int) function shmem_fcollectmem(team, dest, source, nelems) &
          bind(c, name = 'shmem_fcollectmem')
       use, intrinsic :: iso_c_binding
       type(c_ptr), value :: team
       type(c_ptr), value :: dest
       type(c_ptr), value :: source
       integer(c_size_t), value :: nelems
     end function shmem_fcollectmem
  end interface

  interface
     subroutine shmem_sync(team) &
          bind(c, name = 'shmem_sync')
       use, intrinsic :: iso_c_binding
       type(c_ptr), value :: team
     end subroutine shmem_sync
  end interface

  interface
     subroutine shmem_sync_all() &
          bind(c, name = 'shmem_sync_all')
     end subroutine shmem_sync_all
  end interface

  interface
     subroutine shmem_fence() &
          bind(c, name = 'shmem_fence')
     end subroutine shmem_fence
  end interface

  interface
     subroutine shmem_ctx_fence(ctx) &
          bind(c, name = 'shmem_ctx_fence')
       use, intrinsic :: iso_c_binding
       type(c_ptr), value :: ctx
     end subroutine shmem_ctx_fence
  end interface

  !
  ! Memory ordering
  !

  interface
     subroutine shmem_quiet() &
          bind(c, name = 'shmem_quiet')
     end subroutine shmem_quiet
  end interface

  interface
     subroutine shmem_ctx_quiet(ctx) &
          bind(c, name = 'shmem_ctx_quiet')
       use, intrinsic :: iso_c_binding
       type(c_ptr), value :: ctx
     end subroutine shmem_ctx_quiet
  end interface

  !
  ! Distributed Locks
  !

  interface
     subroutine shmem_clear_lock(lock) &
          bind(c, name = 'shmem_clear_lock')
       use, intrinsic :: iso_c_binding
       integer(c_long) :: lock
     end subroutine shmem_clear_lock
  end interface

  interface
     subroutine shmem_set_lock(lock) &
          bind(c, name = 'shmem_set_lock')
       use, intrinsic :: iso_c_binding
       integer(c_long) :: lock
     end subroutine shmem_set_lock
  end interface

  interface
     integer(c_int) function shmem_test_lock(lock) &
          bind(c, name = 'shmem_test_lock')
       use, intrinsic :: iso_c_binding
       integer(c_long) :: lock
     end function shmem_test_lock
  end interface

#endif

end module shmem
