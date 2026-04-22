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
  use num_types, only : rp
  use utils, only : neko_error, neko_warning
  use math, only : add2, rzero, glsc2, col2, copy, add2s1
  use device_math, only : device_rzero, device_col2, device_add2, device_sub3, &
       device_glsc2, device_copy
  use comm
  use mpi_f08, only : MPI_Allreduce, MPI_MIN, MPI_IN_PLACE, MPI_INTEGER
  use coefs, only : coef_t
  use mesh, only : mesh_t
  use space, only : space_t
  use ax_product, only : ax_t
  use scalar_bc_resolver, only : scalar_bc_resolver_t
  use gather_scatter, only : gs_t, GS_OP_ADD
  use tree_amg, only : tamg_hierarchy_t, tamg_lvl_init, tamg_node_init
  use tree_amg_aggregate, only : aggregate_finest_level, aggregate_greedy, &
       aggregate_end, aggregate_pairs
  use tree_amg_smoother, only : amg_cheby_t
  use profiler, only : profiler_start_region, profiler_end_region
  use logger, only : neko_log, LOG_SIZE
  use device, only : device_map, device_free, device_memcpy, HOST_TO_DEVICE, &
       device_get_ptr
  use neko_config, only : NEKO_BCKND_DEVICE
  use, intrinsic :: iso_c_binding
  implicit none
  private

  type :: tamg_wrk_t
     integer :: n = -1
     real(kind=rp), allocatable :: r(:)
     real(kind=rp), allocatable :: b(:)
     real(kind=rp), allocatable :: x(:)
     type(c_ptr) :: r_d = C_NULL_PTR
     type(c_ptr) :: b_d = C_NULL_PTR
     type(c_ptr) :: x_d = C_NULL_PTR
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
     procedure, pass(this) :: free => tamg_mg_free
     procedure, private, pass(this) :: mg_cycle => tamg_mg_cycle
     procedure, private, pass(this) :: mg_cycle_d => tamg_mg_cycle_d
  end type tamg_solver_t

contains

  !> Initialization of the TreeAMG multigrid solver
  !! @param ax Finest level matvec operator
  !! @param Xh Finest level field
  !! @param coef Finest level coeff thing
  !! @param msh Finest level mesh information
  !! @param gs_h Finest level gather scatter operator
  !! @param nlvls Number of levels for the TreeAMG hierarchy
  !! @param bc_resolver Finest level BC resolver
  !! @param max_iter Number of AMG iterations
  subroutine tamg_mg_init(this, ax, Xh, coef, msh, gs_h, nlvls, bc_resolver, &
       max_iter, cheby_degree)
    class(tamg_solver_t), intent(inout), target :: this
    class(ax_t), target, intent(in) :: ax
    type(space_t), target, intent(in) :: Xh
    type(coef_t), target, intent(in) :: coef
    type(mesh_t), target, intent(in) :: msh
    type(gs_t), target, intent(in) :: gs_h
    type(scalar_bc_resolver_t), target, intent(in) :: bc_resolver
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
    call this%amg%init(ax, Xh, coef, msh, gs_h, nlvls, bc_resolver)

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
          call aggregate_greedy(this%amg, mlvl, target_num_aggs, agg_nhbr, &
               nhbr_tmp)
       else
          call print_preagg_info( mlvl, glb_min_target_aggs, 2)
          call aggregate_pairs(this%amg, mlvl, target_num_aggs, agg_nhbr, &
               nhbr_tmp)
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
    allocate(this%wrk(0:(this%amg%nlvls)))
    do lvl = 0, this%amg%nlvls-1
       n = this%amg%lvl(lvl+1)%fine_lvl_dofs
       this%wrk(lvl)%n = n
       allocate( this%wrk(lvl)%r(n) )
       allocate( this%wrk(lvl)%b(n) )
       allocate( this%wrk(lvl)%x(n) )
       if (NEKO_BCKND_DEVICE .eq. 1) then
          call device_map( this%wrk(lvl)%r, this%wrk(lvl)%r_d, n)
          call device_map( this%wrk(lvl)%b, this%wrk(lvl)%b_d, n)
          call device_map( this%wrk(lvl)%x, this%wrk(lvl)%x_d, n)
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
       do i = 0, (size(this%smoo)-1)
          call this%smoo(i)%free()
       end do
       deallocate(this%smoo)
    end if
    if (allocated(this%wrk)) then
       do i = 0, (size(this%wrk)-1)
          if (NEKO_BCKND_DEVICE .eq. 1) then
             call device_free(this%wrk(i)%r_d)
             call device_free(this%wrk(i)%b_d)
             call device_free(this%wrk(i)%x_d)
          end if
          if (allocated(this%wrk(i)%r)) deallocate(this%wrk(i)%r)
          if (allocated(this%wrk(i)%b)) deallocate(this%wrk(i)%b)
          if (allocated(this%wrk(i)%x)) deallocate(this%wrk(i)%x)
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
       ! Zero out the initial guess because we do not handle null spaces very
       ! well.
       call device_rzero(this%wrk(0)%x_d, n)
       call device_copy(this%wrk(0)%b_d, r_d, n)
       zero_initial_guess = .true.
       ! Call the amg cycle
       do iter = 1, max_iter
          call this%mg_cycle_d(zero_initial_guess)
          zero_initial_guess = .false.
       end do
       call device_copy(z_d, this%wrk(0)%x_d, n)
    else
       ! Zero out the initial guess because we do not handle null spaces very
       ! well.
       call rzero(this%wrk(0)%x, n)
       call copy(this%wrk(0)%b, r, n)
       zero_initial_guess = .true.
       ! Call the amg cycle
       do iter = 1, max_iter
          call this%mg_cycle(zero_initial_guess)
          zero_initial_guess = .false.
       end do
       call copy(z, this%wrk(0)%x, n)
    end if
  end subroutine tamg_mg_solve


  !> multigrid cycle for the TreeAMG solver object
  !! @param this The Solver object.
  !! @param zero_initial_guess Flag for when the initial guess is zero
  subroutine tamg_mg_cycle(this, zero_initial_guess)
    class(tamg_solver_t), intent(inout), target :: this
    logical, intent(inout) :: zero_initial_guess
    character(len=2) :: lvl_name
    integer :: max_lvl, lvl

    max_lvl = this%nlvls-1
    ! Loop down hierarchy. Fine to coarse
    do lvl = 0, max_lvl-1
       write(lvl_name, '(I0)') lvl
       call profiler_start_region( "AMG_level_" // trim(lvl_name))
       associate(x => this%wrk(lvl)%x, b => this%wrk(lvl)%b, &
            r => this%wrk(lvl)%r, n => this%wrk(lvl)%n)
         !!----------!!
         !! SMOOTH   !!
         !!----------!!
         call this%smoo(lvl)%solve(x, b, n, this%amg, &
              zero_initial_guess)
         !!----------!!
         !! Residual !!
         !!----------!!
         call calc_resid(r, x, b, this%amg, lvl, n)
         !!----------!!
         !! Restrict !!
         !!----------!!
         call this%amg%interp_f2c(this%wrk(lvl+1)%b, r, lvl+1)

         call rzero(this%wrk(lvl+1)%x, this%wrk(lvl+1)%n)
         zero_initial_guess = .true.
       end associate
       call profiler_end_region( "AMG_level_" // trim(lvl_name))
    end do
    write(lvl_name, '(I0)') max_lvl
    call profiler_start_region( "AMG_level_" // trim(lvl_name))
    !!-------------------!!
    !! Call Coarse solve !!
    !!-------------------!!
    call this%smoo(max_lvl)%solve(this%wrk(max_lvl)%x, &
         this%wrk(max_lvl)%b, this%amg%lvl(max_lvl)%nnodes, this%amg, &
         zero_initial_guess)
    call profiler_end_region( "AMG_level_" // trim(lvl_name))

    zero_initial_guess = .false.
    ! Loop up hierarchy. Coarse to fine
    do lvl = max_lvl-1, 0, -1
       write(lvl_name, '(I0)') lvl
       call profiler_start_region( "AMG_level_" // trim(lvl_name))
       associate(x => this%wrk(lvl)%x, b => this%wrk(lvl)%b, &
            r => this%wrk(lvl)%r, n => this%wrk(lvl)%n)
         !!----------!!
         !! Project  !!
         !!----------!!
         call this%amg%interp_c2f(r, this%wrk(lvl+1)%x, lvl+1)
         !!----------!!
         !! Correct  !!
         !!----------!!
         call add2(x, r, n)
         !!----------!!
         !! SMOOTH   !!
         !!----------!!
         call this%smoo(lvl)%solve(x, b, n, this%amg)
       end associate
       call profiler_end_region( "AMG_level_" // trim(lvl_name))
    end do
  end subroutine tamg_mg_cycle

  !> multigrid cycle for the TreeAMG solver object on device
  !! @param this The Solver object.
  !! @param zero_initial_guess Flag for when the initial guess is zero
  subroutine tamg_mg_cycle_d(this, zero_initial_guess)
    class(tamg_solver_t), intent(inout), target :: this
    logical, intent(inout) :: zero_initial_guess
    character(len=2) :: lvl_name
    integer :: max_lvl, lvl

    max_lvl = this%nlvls-1
    ! Loop down hierarchy. Fine to coarse
    do lvl = 0, max_lvl-1
       write(lvl_name, '(I0)') lvl
       call profiler_start_region( "AMG_level_" // trim(lvl_name))
       associate(x => this%wrk(lvl)%x, x_d => this%wrk(lvl)%x_d, &
            b => this%wrk(lvl)%b, b_d => this%wrk(lvl)%b_d, &
            r => this%wrk(lvl)%r, r_d => this%wrk(lvl)%r_d, &
            n => this%wrk(lvl)%n)
         !!----------!!
         !! SMOOTH   !!
         !!----------!!
         call this%smoo(lvl)%device_solve(x, b, x_d, b_d, n, this%amg, &
              zero_initial_guess)
         !!----------!!
         !! Residual !!
         !!----------!!
         call this%amg%device_matvec(r, x, r_d, x_d, lvl)
         call device_sub3(r_d, b_d, r_d, n)
         !!----------!!
         !! Restrict !!
         !!----------!!
         call this%amg%interp_f2c_d(this%wrk(lvl+1)%b_d, r_d, lvl+1)

         call device_rzero(this%wrk(lvl+1)%x_d, this%wrk(lvl+1)%n)
         zero_initial_guess = .true.
       end associate
       call profiler_end_region( "AMG_level_" // trim(lvl_name))
    end do
    write(lvl_name, '(I0)') max_lvl
    call profiler_start_region( "AMG_level_" // trim(lvl_name))
    !!-------------------!!
    !! Call Coarse solve !!
    !!-------------------!!
    call this%smoo(max_lvl)%device_solve( &
         this%wrk(max_lvl)%x, this%wrk(max_lvl)%b, &
         this%wrk(max_lvl)%x_d, this%wrk(max_lvl)%b_d, &
         this%amg%lvl(max_lvl)%nnodes, this%amg, &
         zero_initial_guess)
    call profiler_end_region( "AMG_level_" // trim(lvl_name))

    zero_initial_guess = .false.
    ! Loop up hierarchy. Coarse to fine
    do lvl = max_lvl-1, 0, -1
       write(lvl_name, '(I0)') lvl
       call profiler_start_region( "AMG_level_" // trim(lvl_name))
       associate(x => this%wrk(lvl)%x, x_d => this%wrk(lvl)%x_d, &
            b => this%wrk(lvl)%b, b_d => this%wrk(lvl)%b_d, &
            r => this%wrk(lvl)%r, r_d => this%wrk(lvl)%r_d, &
            n => this%wrk(lvl)%n)
         !!----------!!
         !! Project  !!
         !!----------!!
         call this%amg%interp_c2f_d(r_d, this%wrk(lvl+1)%x_d, lvl+1, r)
         !!----------!!
         !! Correct  !!
         !!----------!!
         call device_add2(x_d, r_d, n)
         !!----------!!
         !! SMOOTH   !!
         !!----------!!
         call this%smoo(lvl)%device_solve(x, b, x_d, b_d, n, this%amg)
       end associate
       call profiler_end_region( "AMG_level_" // trim(lvl_name))
    end do
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
    call add2s1(r, b, -1.0_rp, n)
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
