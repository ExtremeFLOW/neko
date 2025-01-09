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
  use num_types
  use utils
  use math
  use device_math
  use comm
  use coefs, only : coef_t
  use mesh, only : mesh_t
  use space, only : space_t
  use ax_product, only: ax_t
  use bc_list, only : bc_list_t
  use gather_scatter, only : gs_t, GS_OP_ADD
  use tree_amg, only : tamg_hierarchy_t, tamg_lvl_init, tamg_node_init
  use tree_amg_aggregate
  use tree_amg_smoother
  use logger, only : neko_log, LOG_SIZE
  use device, only: device_map, device_free, c_ptr, C_NULL_PTR, device_memcpy, HOST_TO_DEVICE
  use neko_config, only: NEKO_BCKND_DEVICE
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
    procedure, pass(this) :: device_solve => tamg_mg_solve_device
  end type tamg_solver_t

contains

  !> Initialization of the TreeAMG multigrid solver
  !! @param ax Finest level matvec operator
  !! @param Xh Finest level field
  !! @param coef Finest level coeff thing
  !! @param msh Finest level mesh information
  !! @param gs_h Finest level gather scatter operator
  !! @param nlvls_in Number of levels for the TreeAMG hierarchy
  !! @param blst Finest level BC list
  !! @param max_iter Number of AMG iterations
  subroutine tamg_mg_init(this, ax, Xh, coef, msh, gs_h, nlvls_in, blst, max_iter)
    class(tamg_solver_t), intent(inout), target :: this
    class(ax_t), target, intent(in) :: ax
    type(space_t),target, intent(in) :: Xh
    type(coef_t), target, intent(in) :: coef
    type(mesh_t), target, intent(in) :: msh
    type(gs_t), target, intent(in) :: gs_h
    type(bc_list_t), target, intent(in) :: blst
    integer, intent(in) :: nlvls_in
    integer, intent(in) :: max_iter
    integer :: nlvls, lvl, n, cheby_degree, env_len, mlvl, target_num_aggs
    integer, allocatable :: agg_nhbr(:,:), asdf(:,:)
    character(len=255) :: env_cheby_degree, env_mlvl
    character(len=LOG_SIZE) :: log_buf

    call neko_log%section('AMG')
    
    call get_environment_variable("NEKO_TAMG_MAX_LVL", &
         env_mlvl, env_len)
    if (env_len .eq. 0) then
       !yeah...
       nlvls = nlvls_in
    else
       read(env_mlvl(1:env_len), *) mlvl
       nlvls = mlvl
    end if

    write(log_buf, '(A28,I2,A8)') 'Creating AMG hierarchy with', nlvls, 'levels.'
    call neko_log%message(log_buf)

    allocate( this%amg )
    call this%amg%init(ax, Xh, coef, msh, gs_h, nlvls, blst)

    !> Create level 1 (neko elements are level 0)
    call aggregate_finest_level(this%amg, Xh%lx, Xh%ly, Xh%lz, msh%nelv)

    !> Create the remaining levels
    allocate( agg_nhbr, SOURCE=msh%facet_neigh )
    do mlvl = 2, nlvls-1
      target_num_aggs = this%amg%lvl(mlvl-1)%nnodes / 8
      call print_preagg_info( mlvl, target_num_aggs)
      if ( target_num_aggs .lt. 4 ) then
        call neko_error("TAMG: Too many levels. Not enough DOFs for coarsest grid.")
      end if
      call aggregate_greedy( this%amg, mlvl, target_num_aggs, agg_nhbr, asdf)
      agg_nhbr = asdf
      deallocate( asdf )
    end do
    deallocate( agg_nhbr )

    !> Create the end point
    call aggregate_end(this%amg, nlvls)

    this%max_iter = max_iter

    this%nlvls = this%amg%nlvls!TODO: read from parameter
    if (this%nlvls .gt. this%amg%nlvls) then
      call neko_error("Requested number multigrid levels is greater than the initialized AMG levels")
    end if

    call get_environment_variable("NEKO_TAMG_CHEBY_DEGREE", &
         env_cheby_degree, env_len)
    if (env_len .eq. 0) then
       cheby_degree = 10
    else
       read(env_cheby_degree(1:env_len), *) cheby_degree
    end if
    
    allocate(this%smoo(0:(nlvls)))
    do lvl = 0, nlvls-1
      n = this%amg%lvl(lvl+1)%fine_lvl_dofs
      call this%smoo(lvl)%init(n ,lvl, cheby_degree)
    end do

    !> Allocate work space on each level
    allocate(this%wrk(nlvls))
    do lvl = 1, nlvls
      n = this%amg%lvl(lvl)%fine_lvl_dofs
      allocate( this%wrk(lvl)%r(n) )
      allocate( this%wrk(lvl)%rc(n) )
      allocate( this%wrk(lvl)%tmp(n) )
      if (NEKO_BCKND_DEVICE .eq. 1) then
        call device_map(  this%wrk(lvl)%r, this%wrk(lvl)%r_d, n)
        call device_map(  this%wrk(lvl)%rc, this%wrk(lvl)%rc_d, n)
        call device_map(  this%wrk(lvl)%tmp, this%wrk(lvl)%tmp_d, n)
      end if
    end do

    !allocate(this%jsmoo(0:(nlvls)))
    !do lvl = 0, nlvls-1
    !  n = this%amg%lvl(lvl+1)%fine_lvl_dofs
    !  call this%jsmoo(lvl)%init(n ,lvl, cheby_degree)
    !end do

    call fill_lvl_map(this%amg)

    call neko_log%end_section()
    
  end subroutine tamg_mg_init
 

  !> Solver function for the TreeAMG solver object
  !! @param z The solution to be returned
  !! @param r The right-hand side
  !! @param n Number of dofs
  subroutine tamg_mg_solve(this, z, r, n)
    integer, intent(in) :: n
    class(tamg_solver_t), intent(inout) :: this
    real(kind=rp), dimension(n), intent(inout) :: z
    real(kind=rp), dimension(n), intent(inout) :: r
    integer :: iter, max_iter

    max_iter = this%max_iter

    ! Zero out the initial guess becuase we do not handle null spaces very well...
    z = 0d0

    ! Call the amg cycle
    do iter = 1, max_iter
      call tamg_mg_cycle(z, r, n, 0, this%amg, this)
    end do
  end subroutine tamg_mg_solve

  !> Solver function for the TreeAMG solver object
  !! @param z The solution to be returned
  !! @param r The right-hand side
  !! @param n Number of dofs
  subroutine tamg_mg_solve_device(this, z, r, z_d, r_d, n)
    integer, intent(in) :: n
    class(tamg_solver_t), intent(inout) :: this
    real(kind=rp), dimension(n), intent(inout) :: z
    real(kind=rp), dimension(n), intent(inout) :: r
    type(c_ptr) :: z_d
    type(c_ptr) :: r_d
    integer :: iter, max_iter

    max_iter = this%max_iter

    ! Zero out the initial guess becuase we do not handle null spaces very well...
    call device_rzero(z_d, n)

    ! Call the amg cycle
    do iter = 1, max_iter
      call tamg_mg_cycle_d(z, r, z_d, r_d, n, 0, this%amg, this)
    end do
  end subroutine tamg_mg_solve_device

  !> Recrsive multigrid cycle for the TreeAMG solver object
  !! @param x The solution to be returned
  !! @param b The right-hand side
  !! @param n Number of dofs
  !! @param lvl Current level of the cycle
  !! @param amg The TreeAMG object
  !! @param mgstuff The Solver object. TODO: rename this
  recursive subroutine tamg_mg_cycle(x, b, n, lvl, amg, mgstuff)
    integer, intent(in) :: n
    real(kind=rp), intent(inout) :: x(n)
    real(kind=rp), intent(inout) :: b(n)
    type(tamg_hierarchy_t), intent(inout) :: amg
    type(tamg_solver_t), intent(inout) :: mgstuff
    integer, intent(in) :: lvl
    real(kind=rp) :: r(n)
    real(kind=rp) :: rc(n)
    real(kind=rp) :: tmp(n)
    integer :: iter, num_iter
    integer :: max_lvl
    integer :: i, cyt

    r = 0d0
    rc = 0d0
    tmp = 0d0

    max_lvl = mgstuff%nlvls-1

    !call calc_resid(r,x,b,amg,lvl,n)!> TODO: for debug
    !print *, "LVL:",lvl, "PRE RESID:", sqrt(glsc2(r, r, n))
    !>----------<!
    !> SMOOTH   <!
    !>----------<!
    call mgstuff%smoo(lvl)%solve(x,b, n, amg)
    !call mgstuff%jsmoo(lvl)%solve(x,b, n, amg)
    if (lvl .eq. max_lvl) then !> Is coarsest grid.
      return
    end if

    !>----------<!
    !> Residual <!
    !>----------<!
    call calc_resid(r,x,b,amg,lvl,n)

    !>----------<!
    !> Restrict <!
    !>----------<!
    if (lvl .eq. 0) then
      call average_duplicates(r,amg,lvl,n)
    end if
    call amg%interp_f2c(rc, r, lvl+1)

    !>-------------------<!
    !> Call Coarse solve <!
    !>-------------------<!
    tmp = 0d0
    call tamg_mg_cycle(tmp, rc, amg%lvl(lvl+1)%nnodes, lvl+1, amg, mgstuff)

    !>----------<!
    !> Project  <!
    !>----------<!
    call amg%interp_c2f(r, tmp, lvl+1)
    if (lvl .eq. 0) then
      call average_duplicates(r,amg,lvl,n)
    end if

    !>----------<!
    !> Correct  <!
    !>----------<!
    call add2(x, r, n)

    !>----------<!
    !> SMOOTH   <!
    !>----------<!
    call mgstuff%smoo(lvl)%solve(x,b, n, amg)
    !call mgstuff%jsmoo(lvl)%solve(x,b, n, amg)

    !>----------<!
    !> Residual <!
    !>----------<!
    !call calc_resid(r,x,b,amg,lvl,n)!> TODO: for debug
    !print *, "LVL:",lvl, "POST RESID:", sqrt(glsc2(r, r, n))
  end subroutine tamg_mg_cycle

  !> Recrsive multigrid cycle for the TreeAMG solver object on device
  !! @param x The solution to be returned
  !! @param b The right-hand side
  !! @param n Number of dofs
  !! @param lvl Current level of the cycle
  !! @param amg The TreeAMG object
  !! @param mgstuff The Solver object. TODO: rename this
  recursive subroutine tamg_mg_cycle_d(x, b, x_d, b_d, n, lvl, amg, mgstuff)
    integer, intent(in) :: n
    real(kind=rp), intent(inout) :: x(n)
    real(kind=rp), intent(inout) :: b(n)
    type(c_ptr) :: x_d
    type(c_ptr) :: b_d
    type(tamg_hierarchy_t), intent(inout) :: amg
    type(tamg_solver_t), intent(inout) :: mgstuff
    integer, intent(in) :: lvl
    integer :: iter, num_iter
    integer :: max_lvl
    integer :: i, cyt
    max_lvl = mgstuff%nlvls-1
    !>----------<!
    !> SMOOTH   <!
    !>----------<!
    call mgstuff%smoo(lvl)%device_solve(x, b, x_d, b_d, n, amg)
    if (lvl .eq. max_lvl) then !> Is coarsest grid.
      return
    end if
  associate( r => mgstuff%wrk(lvl+1)%r, r_d => mgstuff%wrk(lvl+1)%r_d, &
             rc => mgstuff%wrk(lvl+1)%rc, rc_d => mgstuff%wrk(lvl+1)%rc_d, &
             tmp => mgstuff%wrk(lvl+1)%tmp, tmp_d => mgstuff%wrk(lvl+1)%tmp_d )
    !>----------<!
    !> Residual <!
    !>----------<!
    !--call device_rzero(r_d, n)
    call amg%device_matvec(r, x, r_d, x_d, lvl)
    call device_sub3(r_d, b_d, r_d, n)
    !>----------<!
    !> Restrict <!
    !>----------<!
    if (lvl .eq. 0) then
      call amg%gs_h%op(r, n, GS_OP_ADD)
      call device_col2(r_d, amg%coef%mult_d, n)
    end if
    call amg%interp_f2c_d(rc_d, r_d, lvl+1)
    !>-------------------<!
    !> Call Coarse solve <!
    !>-------------------<!
    call device_rzero(tmp_d, n)
    call tamg_mg_cycle_d(tmp, rc, tmp_d, rc_d, amg%lvl(lvl+1)%nnodes, lvl+1, amg, mgstuff)
    !>----------<!
    !> Project  <!
    !>----------<!
    call amg%interp_c2f_d(r_d, tmp_d, lvl+1)
    if (lvl .eq. 0) then
      call amg%gs_h%op(r, n, GS_OP_ADD)
      call device_col2(r_d, amg%coef%mult_d, n)
    end if
    !>----------<!
    !> Correct  <!
    !>----------<!
    call device_add2(x_d, r_d, n)
  end associate
    !>----------<!
    !> SMOOTH   <!
    !>----------<!
    call mgstuff%smoo(lvl)%device_solve(x, b, x_d, b_d, n, amg)
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
    r = 0d0
    call amg%matvec(r, x, lvl)
    do  i = 1, n
      r(i) = b(i) - r(i)
    end do
  end subroutine calc_resid

  !> Wrapper function to gather scatter and average the duplicates
  !! @param U The target array
  !! @param amg The TreeAMG object
  !! @param lvl Current level of the cycle
  !! @param n Number of dofs
  subroutine average_duplicates(U, amg, lvl, n)
    integer, intent(in) :: n
    real(kind=rp), intent(inout) :: U(n)
    type(tamg_hierarchy_t), intent(inout) :: amg
    integer, intent(in) :: lvl
    integer :: i
    call amg%gs_h%op(U, n, GS_OP_ADD)
    do  i = 1, n
      U(i) = U(i) * amg%coef%mult(i,1,1,1)
    end do
  end subroutine average_duplicates

  subroutine print_preagg_info(lvl,nagg)
    integer, intent(in) :: lvl,nagg
    character(len=LOG_SIZE) :: log_buf
    !TODO: calculate min and max agg size
    write(log_buf, '(A8,I2,A31)') '-- level',lvl,'-- Calling Greedy Aggregation'
    call neko_log%message(log_buf)
    write(log_buf, '(A33,I6)') 'Target Aggregates:',nagg
    call neko_log%message(log_buf)
  end subroutine print_preagg_info

  subroutine fill_lvl_map(amg)
    type(tamg_hierarchy_t), intent(inout) :: amg
    integer :: i, j, k, l, nid, n
    do j = 1, amg%lvl(1)%nnodes
      do k = 1, amg%lvl(1)%nodes(j)%ndofs
        nid = amg%lvl(1)%nodes(j)%dofs(k)
        amg%lvl(1)%map_f2c_dof(nid) = amg%lvl(1)%nodes(j)%gid
      end do
    end do
    n = size(amg%lvl(1)%map_f2c_dof)
    do l = 2, amg%nlvls
      do i = 1, n
        nid = amg%lvl(l-1)%map_f2c_dof(i)
        do j = 1, amg%lvl(l)%nnodes
          do k = 1, amg%lvl(l)%nodes(j)%ndofs
            if (nid .eq. amg%lvl(l)%nodes(j)%dofs(k)) then
              amg%lvl(l)%map_f2c_dof(i) = amg%lvl(l)%nodes(j)%gid
            end if
          end do
        end do
      end do
    end do
    if (NEKO_BCKND_DEVICE .eq. 1) then
      do l = 1, amg%nlvls
        amg%lvl(l)%map_f2c_dof(0) = n
        call device_memcpy( amg%lvl(l)%map_f2c_dof, amg%lvl(l)%f2c_d, n, HOST_TO_DEVICE, .true.)
        call device_memcpy( amg%lvl(l)%nodes_gids, amg%lvl(l)%nodes_gids_d, amg%lvl(l)%fine_lvl_dofs+1, HOST_TO_DEVICE, .true.)
      end do
    end if
  end subroutine
end module tree_amg_multigrid
