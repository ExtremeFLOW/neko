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
!> Implements multigrid using the TreeAMG hierarchy structure.
!> USE:
!>
!>  type(tamg_hierarchy_t) :: amg
!>  type(tamg_solver_t) :: amg_solver
!>
!>  call amg%init(ax, Xh, coef, msh, gs_h, 3)
!>  call amg_solver%init(amg, niter)
!>
!>  call amg_solver%solve(x%x, f, n)
!>
module tree_amg_multigrid
  use num_types, only: rp
  use utils, only : neko_error, neko_warning
  use math, only : add2, rzero, glsc2, sub3, col2
  use device_math, only : device_rzero, device_col2, device_add2, device_sub3, &
       device_glsc2
  use comm
  use mpi_f08, only: MPI_Allreduce, MPI_MIN, MPI_IN_PLACE, MPI_INTEGER
  use coefs, only : coef_t
  use mesh, only : mesh_t
  use space, only : space_t
  use ax_product, only: ax_t
  use bc_list, only : bc_list_t
  use gather_scatter, only : gs_t, GS_OP_ADD
  use tree_amg, only : tamg_hierarchy_t, tamg_lvl_init, tamg_node_init
  use tree_amg_aggregate, only : aggregate_finest_level, aggregate_greedy, &
       aggregate_end, aggregate_pairs
  use tree_amg_smoother, only : amg_cheby_t
  use logger, only : neko_log, LOG_SIZE
  use device, only: device_map, device_free, device_memcpy, HOST_TO_DEVICE, &
       device_get_ptr
  use neko_config, only: NEKO_BCKND_DEVICE
  use, intrinsic :: iso_c_binding
  implicit none
  private

  type :: tamg_wrk_t
     real(kind=rp), allocatable :: r(:)
     real(kind=rp), allocatable :: rc(:)
     real(kind=rp), allocatable :: tmp(:)
     type(c_ptr) :: r_d = C_NULL_PTR
     type(c_ptr) :: rc_d = C_NULL_PTR
     type(c_ptr) :: tmp_d = C_NULL_PTR
  end type tamg_wrk_t

  !> Type for the TreeAMG solver
  type, public :: tamg_solver_t
     type(tamg_hierarchy_t), allocatable :: amg
     type(amg_cheby_t), allocatable :: smoo(:)
     type(tamg_wrk_t), allocatable :: wrk(:)
     !type(amg_jacobi_t), allocatable :: jsmoo(:)
     integer :: nlvls
     integer :: max_iter
   contains
     procedure, pass(this) :: init => tamg_mg_init
     procedure, pass(this) :: solve => tamg_mg_solve
     procedure, pass(this) :: free =>tamg_mg_free
  end type tamg_solver_t

contains

  !> Initialization of the TreeAMG multigrid solver
  !! @param ax Finest level matvec operator
  !! @param Xh Finest level field
  !! @param coef Finest level coeff thing
  !! @param msh Finest level mesh information
  !! @param gs_h Finest level gather scatter operator
  !! @param nlvls Number of levels for the TreeAMG hierarchy
  !! @param blst Finest level BC list
  !! @param max_iter Number of AMG iterations
  subroutine tamg_mg_init(this, ax, Xh, coef, msh, gs_h, nlvls, blst, &
       max_iter, cheby_degree)
    class(tamg_solver_t), intent(inout), target :: this
    class(ax_t), target, intent(in) :: ax
    type(space_t), target, intent(in) :: Xh
    type(coef_t), target, intent(in) :: coef
    type(mesh_t), target, intent(in) :: msh
    type(gs_t), target, intent(in) :: gs_h
    type(bc_list_t), target, intent(in) :: blst
    integer, intent(in) :: nlvls
    integer, intent(in) :: max_iter
    integer, intent(in) :: cheby_degree
    integer :: lvl, n, mlvl, target_num_aggs
    integer, allocatable :: agg_nhbr(:,:), nhbr_tmp(:,:)
    character(len=LOG_SIZE) :: log_buf
    integer :: glb_min_target_aggs
    logical :: use_greedy_agg

    call neko_log%section('AMG')

    write(log_buf, '(A28,I2,A8)') 'Creating AMG hierarchy with', &
         nlvls, 'levels.'
    call neko_log%message(log_buf)

    allocate( this%amg )
    call this%amg%init(ax, Xh, coef, msh, gs_h, nlvls, blst)

    ! Aggregation
    use_greedy_agg = .true.
    ! Create level 1 (neko elements are level 0)
    call aggregate_finest_level(this%amg, Xh%lx, Xh%ly, Xh%lz, msh%nelv)

    ! Create the remaining levels
    allocate( agg_nhbr, SOURCE = msh%facet_neigh )
    do mlvl = 2, nlvls-1
       ! estimate number of aggregates
       if (use_greedy_agg) then
          target_num_aggs = this%amg%lvl(mlvl-1)%nnodes / 8
       else
          target_num_aggs = this%amg%lvl(mlvl-1)%nnodes / 2
       end if

       glb_min_target_aggs = target_num_aggs
       call MPI_Allreduce(MPI_IN_PLACE, glb_min_target_aggs, 1, &
            MPI_INTEGER, MPI_MIN, NEKO_COMM)
       if (glb_min_target_aggs .lt. 4 ) then
          call neko_warning( &
               "TAMG: Too many levels. Not enough DOFs for coarsest grid.")
          this%amg%nlvls = mlvl
          exit
       end if

       if (use_greedy_agg) then
          call print_preagg_info( mlvl, glb_min_target_aggs, 1)
          call aggregate_greedy( this%amg, mlvl, target_num_aggs, agg_nhbr, nhbr_tmp)
       else
          call print_preagg_info( mlvl, glb_min_target_aggs, 2)
          call aggregate_pairs( this%amg, mlvl, target_num_aggs, agg_nhbr, nhbr_tmp)
       end if

       agg_nhbr = nhbr_tmp
       deallocate( nhbr_tmp )
    end do
    deallocate( agg_nhbr )

    ! Create the end point
    call aggregate_end(this%amg, this%amg%nlvls)

    this%max_iter = max_iter

    this%nlvls = this%amg%nlvls
    if (this%nlvls .gt. this%amg%nlvls) then
       call neko_error( &
            "Requested number multigrid levels &
       & is greater than the initialized AMG levels")
    end if

    ! Initialize relaxation methods
    allocate(this%smoo(0:(this%amg%nlvls)))
    do lvl = 0, this%amg%nlvls-1
       n = this%amg%lvl(lvl+1)%fine_lvl_dofs
       call this%smoo(lvl)%init(n, lvl, cheby_degree)
    end do

    ! Allocate work space on each level
    allocate(this%wrk(this%amg%nlvls))
    do lvl = 1, this%amg%nlvls
       n = this%amg%lvl(lvl)%fine_lvl_dofs
       allocate( this%wrk(lvl)%r(n) )
       allocate( this%wrk(lvl)%rc(n) )
       allocate( this%wrk(lvl)%tmp(n) )
       if (NEKO_BCKND_DEVICE .eq. 1) then
          call device_map( this%wrk(lvl)%r, this%wrk(lvl)%r_d, n)
          call device_map( this%wrk(lvl)%rc, this%wrk(lvl)%rc_d, n)
          call device_map( this%wrk(lvl)%tmp, this%wrk(lvl)%tmp_d, n)
       end if
    end do

    !allocate(this%jsmoo(0:(this%amg%nlvls)))
    !do lvl = 0, this%amg%nlvls-1
    !  n = this%amg%lvl(lvl+1)%fine_lvl_dofs
    !  call this%jsmoo(lvl)%init(n ,lvl, cheby_degree)
    !end do

    ! Create index mapping between levels
    call fill_lvl_map(this%amg)

    call neko_log%end_section()

  end subroutine tamg_mg_init

  !> free tree amg solver object
  subroutine tamg_mg_free(this)
    class(tamg_solver_t), intent(inout), target :: this
    integer :: i
    if (allocated(this%amg)) then
       call this%amg%free()
       deallocate(this%amg)
    end if
    if (allocated(this%smoo)) then
       do i = 1, size(this%smoo)
          call this%smoo(i)%free()
       end do
       deallocate(this%smoo)
    end if
    if (allocated(this%wrk)) then
       do i = 1, size(this%wrk)
          if (NEKO_BCKND_DEVICE .eq. 1) then
             call device_free(this%wrk(i)%r_d)
             call device_free(this%wrk(i)%rc_d)
             call device_free(this%wrk(i)%tmp_d)
          end if
          if (allocated(this%wrk(i)%r)) deallocate(this%wrk(i)%r)
          if (allocated(this%wrk(i)%rc)) deallocate(this%wrk(i)%rc)
          if (allocated(this%wrk(i)%tmp)) deallocate(this%wrk(i)%tmp)
       end do
    end if
  end subroutine tamg_mg_free


  !> Solver function for the TreeAMG solver object
  !! @param z The solution to be returned
  !! @param r The right-hand side
  !! @param n Number of dofs
  subroutine tamg_mg_solve(this, z, r, n)
    integer, intent(in) :: n
    class(tamg_solver_t), intent(inout) :: this
    real(kind=rp), dimension(n), intent(inout) :: z
    real(kind=rp), dimension(n), intent(inout) :: r
    type(c_ptr) :: z_d
    type(c_ptr) :: r_d
    integer :: iter, max_iter
    logical :: zero_initial_guess

    max_iter = this%max_iter

    if (NEKO_BCKND_DEVICE .eq. 1) then
       z_d = device_get_ptr(z)
       r_d = device_get_ptr(r)
       ! Zero out the initial guess becuase we do not handle null spaces very well...
       call device_rzero(z_d, n)
       zero_initial_guess = .true.
       ! Call the amg cycle
       do iter = 1, max_iter
          call tamg_mg_cycle_d(z, r, z_d, r_d, n, 0, this%amg, this, &
               zero_initial_guess)
          zero_initial_guess = .false.
       end do
    else
       ! Zero out the initial guess becuase we do not handle null spaces very well...
       call rzero(z, n)
       zero_initial_guess = .true.
       ! Call the amg cycle
       do iter = 1, max_iter
          call tamg_mg_cycle(z, r, n, 0, this%amg, this, &
               zero_initial_guess)
          zero_initial_guess = .false.
       end do
    end if
  end subroutine tamg_mg_solve


  !> Recrsive multigrid cycle for the TreeAMG solver object
  !! @param x The solution to be returned
  !! @param b The right-hand side
  !! @param n Number of dofs
  !! @param lvl Current level of the cycle
  !! @param amg The TreeAMG object
  !! @param mgstuff The Solver object. TODO: rename this
  recursive subroutine tamg_mg_cycle(x, b, n, lvl, amg, mgstuff, &
       zero_initial_guess)
    integer, intent(in) :: n
    real(kind=rp), intent(inout) :: x(n)
    real(kind=rp), intent(inout) :: b(n)
    type(tamg_hierarchy_t), intent(inout) :: amg
    type(tamg_solver_t), intent(inout) :: mgstuff
    logical, intent(in) :: zero_initial_guess
    integer, intent(in) :: lvl
    real(kind=rp) :: r(n)
    real(kind=rp) :: rc(n)
    real(kind=rp) :: tmp(n)
    integer :: iter, num_iter
    integer :: max_lvl
    integer :: i, cyt
    max_lvl = mgstuff%nlvls-1
    !!----------!!
    !! SMOOTH   !!
    !!----------!!
    call mgstuff%smoo(lvl)%solve(x, b, n, amg, &
         zero_initial_guess)
    if (lvl .eq. max_lvl) then !> Is coarsest grid.
       return
    end if
    !!----------!!
    !! Residual !!
    !!----------!!
    call calc_resid(r, x, b, amg, lvl, n)
    !!----------!!
    !! Restrict !!
    !!----------!!
    call amg%interp_f2c(rc, r, lvl+1)
    !!-------------------!!
    !! Call Coarse solve !!
    !!-------------------!!
    call rzero(tmp, n)
    call tamg_mg_cycle(tmp, rc, amg%lvl(lvl+1)%nnodes, lvl+1, amg, mgstuff, &
         .true.)
    !!----------!!
    !! Project  !!
    !!----------!!
    call amg%interp_c2f(r, tmp, lvl+1)
    !!----------!!
    !! Correct  !!
    !!----------!!
    call add2(x, r, n)
    !!----------!!
    !! SMOOTH   !!
    !!----------!!
    call mgstuff%smoo(lvl)%solve(x,b, n, amg)
  end subroutine tamg_mg_cycle

  !> Recrsive multigrid cycle for the TreeAMG solver object on device
  !! @param x The solution to be returned
  !! @param b The right-hand side
  !! @param n Number of dofs
  !! @param lvl Current level of the cycle
  !! @param amg The TreeAMG object
  !! @param mgstuff The Solver object. TODO: rename this
  recursive subroutine tamg_mg_cycle_d(x, b, x_d, b_d, n, lvl, amg, mgstuff, &
       zero_initial_guess)
    integer, intent(in) :: n
    real(kind=rp), intent(inout) :: x(n)
    real(kind=rp), intent(inout) :: b(n)
    type(c_ptr) :: x_d
    type(c_ptr) :: b_d
    type(tamg_hierarchy_t), intent(inout) :: amg
    type(tamg_solver_t), intent(inout) :: mgstuff
    logical, intent(in) :: zero_initial_guess
    integer, intent(in) :: lvl
    integer :: iter, num_iter
    integer :: max_lvl
    integer :: i, cyt
    max_lvl = mgstuff%nlvls-1
    !!----------!!
    !! SMOOTH   !!
    !!----------!!
    call mgstuff%smoo(lvl)%device_solve(x, b, x_d, b_d, n, amg, &
         zero_initial_guess)
    if (lvl .eq. max_lvl) then !> Is coarsest grid.
       return
    end if
    associate( r => mgstuff%wrk(lvl+1)%r, r_d => mgstuff%wrk(lvl+1)%r_d, &
         rc => mgstuff%wrk(lvl+1)%rc, rc_d => mgstuff%wrk(lvl+1)%rc_d, &
         tmp => mgstuff%wrk(lvl+1)%tmp, tmp_d => mgstuff%wrk(lvl+1)%tmp_d )
      !!----------!!
      !! Residual !!
      !!----------!!
      call amg%device_matvec(r, x, r_d, x_d, lvl)
      call device_sub3(r_d, b_d, r_d, n)
      !!----------!!
      !! Restrict !!
      !!----------!!
      call amg%interp_f2c_d(rc_d, r_d, lvl+1)
      !!-------------------!!
      !! Call Coarse solve !!
      !!-------------------!!
      call device_rzero(tmp_d, n)
      call tamg_mg_cycle_d(tmp, rc, tmp_d, rc_d, &
           amg%lvl(lvl+1)%nnodes, lvl+1, amg, mgstuff, .true.)
      !!----------!!
      !! Project  !!
      !!----------!!
      call amg%interp_c2f_d(r_d, tmp_d, lvl+1, r)
      !!----------!!
      !! Correct  !!
      !!----------!!
      call device_add2(x_d, r_d, n)
      !!----------!!
      !! SMOOTH   !!
      !!----------!!
      call mgstuff%smoo(lvl)%device_solve(x, b, x_d, b_d, n, amg)
    end associate
  end subroutine tamg_mg_cycle_d


  !> Wrapper function to calculate residyal
  !! @param r The residual to be returned
  !! @param x The current solution
  !! @param b The right-hand side
  !! @param amg The TreeAMG object
  !! @param lvl Current level of the cycle
  !! @param n Number of dofs
  subroutine calc_resid(r, x, b, amg, lvl, n)
    integer, intent(in) :: n
    real(kind=rp), intent(inout) :: r(n)
    real(kind=rp), intent(inout) :: x(n)
    real(kind=rp), intent(inout) :: b(n)
    type(tamg_hierarchy_t), intent(inout) :: amg
    integer, intent(in) :: lvl
    integer :: i
    call amg%matvec(r, x, lvl)
    call sub3(r, b, r, n)
  end subroutine calc_resid


  subroutine print_preagg_info(lvl, nagg, agg_type)
    integer, intent(in) :: lvl, nagg, agg_type
    character(len=LOG_SIZE) :: log_buf
    !TODO: calculate min and max agg size
    if (agg_type .eq. 1) then
       write(log_buf, '(A8,I2,A31)') '-- level', lvl, &
            '-- Calling Greedy Aggregation'
    else if (agg_type .eq. 2) then
       write(log_buf, '(A8,I2,A33)') '-- level', lvl, &
            '-- Calling Pairwise Aggregation'
    else
       write(log_buf, '(A8,I2,A31)') '-- level', lvl, &
            '-- UNKNOWN Aggregation'
    end if
    call neko_log%message(log_buf)
    write(log_buf, '(A33,I6)') 'Target Aggregates:', nagg
    call neko_log%message(log_buf)
  end subroutine print_preagg_info

  subroutine print_resid_info(r, x, b, r_d, x_d, b_d, amg, lvl, n)
    integer, intent(in) :: lvl, n
    real(kind=rp), intent(inout) :: r(n)
    real(kind=rp), intent(inout) :: x(n)
    real(kind=rp), intent(inout) :: b(n)
    type(c_ptr) :: r_d
    type(c_ptr) :: x_d
    type(c_ptr) :: b_d
    type(tamg_hierarchy_t), intent(inout) :: amg
    real(kind=rp) :: val
    character(len=LOG_SIZE) :: log_buf

    call amg%device_matvec(r, x, r_d, x_d, lvl)
    call device_sub3(r_d, b_d, r_d, n)
    val = device_glsc2(r_d, r_d, n)

    write(log_buf, '(A33,I6,F12.6)') 'tAMG resid:', lvl, val
    call neko_log%message(log_buf)
  end subroutine print_resid_info

  !> Create index mapping between levels and directly to finest level
  !! @param amg The tamg hierarchy
  subroutine fill_lvl_map(amg)
    type(tamg_hierarchy_t), intent(inout) :: amg
    integer :: i, j, k, l, nid, n
    do j = 1, amg%lvl(1)%nnodes
       do k = 1, amg%lvl(1)%nodes(j)%ndofs
          nid = amg%lvl(1)%nodes(j)%dofs(k)
          amg%lvl(1)%map_finest2lvl(nid) = amg%lvl(1)%nodes(j)%gid
       end do
    end do
    n = size(amg%lvl(1)%map_finest2lvl)
    do l = 2, amg%nlvls
       do i = 1, n
          nid = amg%lvl(l-1)%map_finest2lvl(i)
          do j = 1, amg%lvl(l)%nnodes
             do k = 1, amg%lvl(l)%nodes(j)%ndofs
                if (nid .eq. amg%lvl(l)%nodes(j)%dofs(k)) then
                   amg%lvl(l)%map_finest2lvl(i) = amg%lvl(l)%nodes(j)%gid
                end if
             end do
          end do
       end do
    end do
    if (NEKO_BCKND_DEVICE .eq. 1) then
       do l = 1, amg%nlvls
          amg%lvl(l)%map_finest2lvl(0) = n
          call device_memcpy( amg%lvl(l)%map_finest2lvl, &
               amg%lvl(l)%map_finest2lvl_d, n, &
               HOST_TO_DEVICE, .true.)
          call device_memcpy( amg%lvl(l)%map_f2c, &
               amg%lvl(l)%map_f2c_d, amg%lvl(l)%fine_lvl_dofs+1, &
               HOST_TO_DEVICE, .true.)
       end do
    end if
  end subroutine fill_lvl_map
end module tree_amg_multigrid
