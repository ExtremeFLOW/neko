! Copyright (c) 2026, The Neko Authors
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
!> Checks of a Krylov solver's operator and preconditioner, and of the
!! residual it reports.
!!
!! `krylov_check_setup` measures, with deterministic test vectors, whether
!! the (masked, assembled) operator is symmetric and positive, and whether
!! the preconditioner is symmetric, positive, fixed and linear. It then
!! compares these properties with what the chosen Krylov method requires:
!! the CG family needs a symmetric positive definite operator and a
!! symmetric, positive, fixed and linear preconditioner; BiCGStab needs a
!! fixed and linear preconditioner; flexible GMRES needs nothing. A mismatch
!! is reported as a warning, and the simulation continues.
!!
!! The work fields of both checks are taken from the scratch registry, which
!! the time loop fills with more fields than the checks need, so the checks
!! add no memory. Fields are allocated locally only when the registry is not
!! set up on the solver's dofmap.
!!
!! `krylov_check_residual` recomputes the true residual \f$ \|f - Ax\| \f$
!! after a solve and compares it with the residual the solver reported, in
!! the solver's own norm. The two differ when a solver's recursively updated
!! residual drifts (PipeCG, or GMRES after many iterations).
module krylov_check
  use num_types, only : rp, i8
  use krylov, only : ksp_t, ksp_monitor_t
  use precon, only : pc_t
  use ax_product, only : ax_t
  use coefs, only : coef_t
  use gather_scatter, only : gs_t
  use gs_ops, only : GS_OP_ADD
  use scalar_bc_projector, only : scalar_bc_projector_t
  use vector_bc_projector, only : vector_bc_projector_t, &
       segregated_vector_bc_projector_t
  use field, only : field_t, field_ptr_t
  use scratch_registry, only : neko_scratch_registry
  use field_math, only : field_copy, field_add2, field_sub2
  use math, only : glsc3, col2, glmin, glmax, NEKO_EPS
  use operators, only : ortho
  use device, only : device_memcpy, HOST_TO_DEVICE
  use device_math, only : device_glsc3, device_col2
  use opr_device, only : device_ortho
  use neko_config, only : NEKO_BCKND_DEVICE
  use logger, only : neko_log, LOG_SIZE
  use utils, only : neko_warning
  implicit none
  private

  public :: krylov_check_setup, krylov_check_residual, krylov_check_due, &
       krylov_requirements, krylov_check_compatible

  !> Check a solver's operator and preconditioner, for a scalar or a
  !! (component-wise) vector system.
  interface krylov_check_setup
     module procedure krylov_check_setup_scalar, krylov_check_setup_vector
  end interface krylov_check_setup

  !> Compare the true with the reported residual, for a scalar or a vector
  !! system.
  interface krylov_check_residual
     module procedure krylov_check_residual_scalar, &
          krylov_check_residual_vector
  end interface krylov_check_residual

  !> Number of test-vector pairs
  integer, parameter :: NPAIRS = 3
  !> Relative deviation of the true from the reported residual that is
  !! reported as a warning
  real(kind=rp), parameter :: RESIDUAL_WARN = 0.2_rp

contains

  !> The properties a Krylov method requires of its operator and
  !! preconditioner.
  !! @param ksp_type Solver type name, as given to the factory.
  !! @param need_spd Whether the operator must be symmetric positive definite.
  !! @param need_sym Whether the preconditioner must be symmetric and positive.
  !! @param need_fixed Whether the preconditioner must be fixed and linear.
  !! @param known Whether the solver type is known to this table.
  pure subroutine krylov_requirements(ksp_type, need_spd, need_sym, &
       need_fixed, known)
    character(len=*), intent(in) :: ksp_type
    logical, intent(out) :: need_spd, need_sym, need_fixed, known

    need_spd = .false.
    need_sym = .false.
    need_fixed = .false.
    known = .true.
    select case (trim(ksp_type))
    case ('cg', 'pipecg', 'fused_cg', 'cacg', 'coupled_cg', &
         'fused_coupled_cg')
       need_spd = .true.
       need_sym = .true.
       need_fixed = .true.
    case ('bicgstab', 'coupled_bicgstab')
       need_fixed = .true.
    case ('gmres')
       ! flexible GMRES: the preconditioner may vary
    case default
       known = .false.
    end select
  end subroutine krylov_requirements

  !> Whether measured properties satisfy a Krylov method's requirements.
  !! @param ksp_type Solver type name.
  !! @param op_sym Operator is symmetric.
  !! @param op_pos Operator is positive.
  !! @param pc_sym Preconditioner is symmetric.
  !! @param pc_pos Preconditioner is positive.
  !! @param pc_fixed Preconditioner is fixed (same output for the same input).
  !! @param pc_linear Preconditioner is linear.
  !! @param reason Why the pairing is not compatible, empty if it is.
  pure subroutine krylov_check_compatible(ksp_type, op_sym, op_pos, pc_sym, &
       pc_pos, pc_fixed, pc_linear, compatible, reason)
    character(len=*), intent(in) :: ksp_type
    logical, intent(in) :: op_sym, op_pos, pc_sym, pc_pos, pc_fixed, pc_linear
    logical, intent(out) :: compatible
    character(len=*), intent(out) :: reason
    logical :: need_spd, need_sym, need_fixed, known

    call krylov_requirements(ksp_type, need_spd, need_sym, need_fixed, known)
    compatible = .true.
    reason = ''
    if (.not. known) return

    if (need_spd .and. .not. (op_sym .and. op_pos)) then
       compatible = .false.
       reason = trim(ksp_type) // ' requires a symmetric positive ' // &
            'definite operator'
    else if (need_sym .and. .not. (pc_sym .and. pc_pos)) then
       compatible = .false.
       reason = trim(ksp_type) // ' requires a symmetric positive ' // &
            'preconditioner'
    else if (need_fixed .and. .not. (pc_fixed .and. pc_linear)) then
       compatible = .false.
       reason = trim(ksp_type) // ' requires a fixed, linear preconditioner'
    end if
  end subroutine krylov_check_compatible

  !> Whether the true residual is to be checked at this time step: always at
  !! the first step, then every @a interval steps (0 disables the periodic
  !! check).
  pure function krylov_check_due(tstep, interval) result(due)
    integer, intent(in) :: tstep, interval
    logical :: due

    due = (tstep .eq. 1)
    if (interval .gt. 0) due = due .or. (mod(tstep, interval) .eq. 0)
  end function krylov_check_due

  !> Measure the properties of a solver's operator and preconditioner and
  !! compare them with what the solver requires.
  !! @param name Name of the equation, for the log.
  !! @param ksp The Krylov solver, with its preconditioner set.
  !! @param Ax The operator.
  !! @param coef Coefficients (mesh, function space, multiplicity).
  !! @param gs Gather-scatter of the solver's function space.
  !! @param bclst Projector onto the strong boundary conditions (mask).
  !! @param demean Whether the system has a constant null space that is
  !! removed by de-meaning (pure Neumann problems).
  !! @param pc_type Preconditioner type name, for the log (optional).
  !! @note Uses four work fields and about a dozen preconditioner
  !! applications.
  subroutine krylov_check_setup_scalar(name, ksp, Ax, coef, gs, bclst, &
       demean, pc_type)
    character(len=*), intent(in) :: name
    class(ksp_t), intent(inout) :: ksp
    class(ax_t), intent(inout) :: Ax
    type(coef_t), intent(inout) :: coef
    type(gs_t), intent(inout) :: gs
    type(scalar_bc_projector_t), intent(inout) :: bclst
    logical, intent(in) :: demean
    character(len=*), intent(in), optional :: pc_type
    type(field_ptr_t) :: w(4)
    type(field_t), pointer :: x, y, t1, t2
    integer :: w_idx(4)
    logical :: w_scratch
    integer(kind=i8) :: glb_n
    integer :: n, k
    real(kind=rp) :: tol, a1, a2, nx, ny, nax, nmx, nmy
    real(kind=rp) :: asym_a, cos_a, asym_m, cos_m, fixed_m, lin_m
    logical :: op_sym, op_pos, pc_sym, pc_pos, pc_fixed, pc_linear
    logical :: compatible
    character(len=LOG_SIZE) :: log_buf, reason, pc_name

    if (.not. associated(ksp%M)) return

    n = coef%dof%size()
    glb_n = int(coef%msh%glb_nelv, i8) * int(coef%Xh%lxyz, i8)
    ! Round-off in the relative measures grows with the number of points;
    ! the cap keeps the threshold meaningful in single precision
    tol = 1.0e3_rp * NEKO_EPS * sqrt(real(glb_n, rp))
    tol = min(max(tol, 1.0e-12_rp), 1.0e-3_rp)

    pc_name = 'preconditioner'
    if (present(pc_type)) pc_name = trim(pc_type)

    call work_fields_get(w, w_idx, coef, w_scratch)
    x => w(1)%ptr
    y => w(2)%ptr
    t1 => w(3)%ptr
    t2 => w(4)%ptr

    asym_a = 0.0_rp
    cos_a = huge(1.0_rp)
    asym_m = 0.0_rp
    cos_m = huge(1.0_rp)
    fixed_m = 0.0_rp
    lin_m = 0.0_rp
    do k = 1, NPAIRS
       call test_vector(x, 2 * k - 1)
       call test_vector(y, 2 * k)
       nx = sqrt(dot(x, x))
       ny = sqrt(dot(y, y))

       ! Operator: symmetry <y, Ax> = <x, Ay> and positivity <x, Ax> > 0
       call apply_a(t1, x)
       nax = sqrt(dot(t1, t1))
       a1 = dot(y, t1)
       cos_a = min(cos_a, dot(x, t1) / (nax * nx))
       call apply_a(t1, y)
       a2 = dot(x, t1)
       asym_a = max(asym_a, abs(a1 - a2) / (nax * ny))

       ! Preconditioner: symmetry and positivity, t1 = Mx
       call ksp%M%solve(t1%x, x%x, n)
       nmx = sqrt(dot(t1, t1))
       a1 = dot(y, t1)
       cos_m = min(cos_m, dot(x, t1) / (nmx * nx))

       ! Fixedness: the same input twice
       call ksp%M%solve(t2%x, x%x, n)
       call field_sub2(t2, t1, n)
       fixed_m = max(fixed_m, sqrt(dot(t2, t2)) / nmx)

       ! t2 = My
       call ksp%M%solve(t2%x, y%x, n)
       nmy = sqrt(dot(t2, t2))
       a2 = dot(x, t2)
       asym_m = max(asym_m, abs(a1 - a2) / (nmx * ny))

       ! Linearity: M(x + y) against Mx + My; x is not needed afterwards
       call field_add2(t1, t2, n)
       call field_add2(x, y, n)
       call ksp%M%solve(t2%x, x%x, n)
       call field_sub2(t2, t1, n)
       lin_m = max(lin_m, sqrt(dot(t2, t2)) / sqrt(nmx**2 + nmy**2))
    end do

    op_sym = asym_a .lt. tol
    op_pos = cos_a .gt. 0.0_rp
    pc_sym = asym_m .lt. tol
    pc_pos = cos_m .gt. 0.0_rp
    pc_fixed = fixed_m .lt. tol
    pc_linear = lin_m .lt. tol

    call neko_log%section('Solver check: ' // trim(name))
    write(log_buf, '(A,A,ES9.2,A,A,F6.3,A)') 'Operator    : ', &
         trim(merge('symmetric    ', 'not symmetric', op_sym)) // ' (', &
         asym_a, '), ', trim(merge('positive    ', 'not positive', op_pos)) &
         // ' (', cos_a, ')'
    call neko_log%message(log_buf)
    write(log_buf, '(A12,A,A,ES9.2,A,A,ES9.2,A)') pc_name, ': ', &
         trim(merge('symmetric    ', 'not symmetric', pc_sym)) // ' (', &
         asym_m, '), ', trim(merge('positive    ', 'not positive', pc_pos)) &
         // ' (', cos_m, ')'
    call neko_log%message(log_buf)
    write(log_buf, '(A14,A,ES9.2,A,A,ES9.2,A)') ' ', &
         trim(merge('fixed    ', 'not fixed', pc_fixed)) // ' (', &
         fixed_m, '), ', trim(merge('linear    ', 'not linear', pc_linear)) &
         // ' (', lin_m, ')'
    call neko_log%message(log_buf)
    if (pc_pos .and. cos_m .lt. 1.0e-2_rp) then
       call neko_log%message('Note: the preconditioner is close to ' // &
            'indefinite; expect slow convergence')
    end if

    call krylov_check_compatible(trim(ksp%type_name), op_sym, op_pos, pc_sym, &
         pc_pos, pc_fixed, pc_linear, compatible, reason)
    if (compatible) then
       call neko_log%message('Pairing     : ' // trim(ksp%type_name) // &
            ' + ' // trim(pc_name) // ' is compatible')
    else
       call neko_warning(trim(name) // ' solver: ' // trim(reason) // &
            ', which ' // trim(pc_name) // ' is not here. ' // &
            'The solve may stall or converge to a wrong solution; ' // &
            'gmres accepts any preconditioner.')
    end if
    call neko_log%end_section()

    call work_fields_put(w, w_idx, w_scratch)

  contains

    !> A deterministic, continuous, masked test vector, built from smooth
    !! functions of the coordinates so that it does not depend on the
    !! partitioning.
    subroutine test_vector(v, seed)
      type(field_t), intent(inout) :: v
      integer, intent(in) :: seed
      real(kind=rp) :: xmin, xmax, ymin, ymax, zmin, zmax
      real(kind=rp) :: kx, ky, kz, px, py, pz
      integer :: i

      ! Wave numbers from the bounding box, a few waves per direction
      call bounding_box(xmin, xmax, ymin, ymax, zmin, zmax)
      kx = 2.0_rp * acos(-1.0_rp) * real(2 + seed, rp) / max(xmax - xmin, &
           epsilon(1.0_rp))
      ky = 2.0_rp * acos(-1.0_rp) * real(3 + seed, rp) / max(ymax - ymin, &
           epsilon(1.0_rp))
      kz = 2.0_rp * acos(-1.0_rp) * real(1 + seed, rp) / max(zmax - zmin, &
           epsilon(1.0_rp))
      do i = 1, n
         px = coef%dof%x%x(i,1,1,1) - xmin
         py = coef%dof%y%x(i,1,1,1) - ymin
         pz = coef%dof%z%x(i,1,1,1) - zmin
         v%x(i,1,1,1) = sin(kx * px + 0.3_rp * seed) * cos(ky * py) &
              + 0.5_rp * cos(kz * pz + 0.7_rp * seed) * sin(ky * py + 1.0_rp) &
              + 0.25_rp * sin(0.37_rp * kx * px * ky * py / &
              max(kx * (xmax - xmin), 1.0_rp) + seed)
      end do
      if (NEKO_BCKND_DEVICE .eq. 1) then
         call device_memcpy(v%x, v%x_d, n, HOST_TO_DEVICE, sync = .true.)
      end if
      ! Continuous (gather-scatter sum, multiplicity weights), masked at the
      ! strong boundaries, and de-meaned for a pure Neumann problem
      call gs%op(v, GS_OP_ADD)
      if (NEKO_BCKND_DEVICE .eq. 1) then
         call device_col2(v%x_d, coef%mult_d, n)
      else
         call col2(v%x, coef%mult, n)
      end if
      call bclst%apply(v%x, n)
      if (demean) then
         if (NEKO_BCKND_DEVICE .eq. 1) then
            call device_ortho(v%x_d, glb_n, n)
         else
            call ortho(v%x, glb_n, n)
         end if
      end if
    end subroutine test_vector

    !> Global bounding box of the GLL points.
    subroutine bounding_box(xmin, xmax, ymin, ymax, zmin, zmax)
      real(kind=rp), intent(out) :: xmin, xmax, ymin, ymax, zmin, zmax
      xmin = glmin(coef%dof%x%x, n)
      xmax = glmax(coef%dof%x%x, n)
      ymin = glmin(coef%dof%y%x, n)
      ymax = glmax(coef%dof%y%x, n)
      zmin = glmin(coef%dof%z%x, n)
      zmax = glmax(coef%dof%z%x, n)
    end subroutine bounding_box

    !> The assembled, masked operator.
    subroutine apply_a(w, v)
      type(field_t), intent(inout) :: w, v
      call Ax%compute(w%x, v%x, coef, coef%msh, coef%Xh)
      call gs%op(w, GS_OP_ADD)
      call bclst%apply(w%x, n)
    end subroutine apply_a

    !> The solver's inner product (multiplicity weighted).
    function dot(a, b) result(s)
      type(field_t), intent(in) :: a, b
      real(kind=rp) :: s
      if (NEKO_BCKND_DEVICE .eq. 1) then
         s = device_glsc3(a%x_d, coef%mult_d, b%x_d, n)
      else
         s = glsc3(a%x, coef%mult, b%x, n)
      end if
    end function dot

  end subroutine krylov_check_setup_scalar

  !> Check the solver of a vector system whose components are solved one
  !! at a time (segregated boundary conditions), using the x component.
  !! Coupled vector solvers are not checked.
  subroutine krylov_check_setup_vector(name, ksp, Ax, coef, gs, bclst, &
       pc_type)
    character(len=*), intent(in) :: name
    class(ksp_t), intent(inout) :: ksp
    class(ax_t), intent(inout) :: Ax
    type(coef_t), intent(inout) :: coef
    type(gs_t), intent(inout) :: gs
    class(vector_bc_projector_t), intent(inout) :: bclst
    character(len=*), intent(in), optional :: pc_type

    select type (bclst)
    type is (segregated_vector_bc_projector_t)
       call krylov_check_setup_scalar(name, ksp, Ax, coef, gs, bclst%x, &
            .false., pc_type)
    end select
  end subroutine krylov_check_setup_vector

  !> Recompute the true residual of a solve and compare it with the
  !! residual the solver reported.
  !! @param name Name of the equation, for the log.
  !! @param Ax The operator.
  !! @param x The solution.
  !! @param f The right-hand side the solver was given (assembled, masked).
  !! @param coef Coefficients.
  !! @param gs Gather-scatter.
  !! @param bclst Projector onto the strong boundary conditions (mask).
  !! @param results The solver's monitor, with the reported final residual.
  !! @param tol The solver's absolute tolerance; differences well below it
  !! are not reported.
  subroutine krylov_check_residual_scalar(name, Ax, x, f, coef, gs, bclst, &
       results, tol)
    character(len=*), intent(in) :: name
    class(ax_t), intent(inout) :: Ax
    type(field_t), intent(inout) :: x, f
    type(coef_t), intent(inout) :: coef
    type(gs_t), intent(inout) :: gs
    type(scalar_bc_projector_t), intent(inout) :: bclst
    type(ksp_monitor_t), intent(in) :: results
    real(kind=rp), intent(in) :: tol
    type(field_ptr_t) :: w(1)
    type(field_t), pointer :: r
    integer :: w_idx(1), n
    logical :: w_scratch

    n = coef%dof%size()
    call work_fields_get(w, w_idx, coef, w_scratch)
    r => w(1)%ptr
    call Ax%compute(r%x, x%x, coef, coef%msh, coef%Xh)
    call gs%op(r, GS_OP_ADD)
    call bclst%apply(r%x, n)
    call field_sub2(r, f, n)
    call report_residual(name, residual_norm(r, coef), results, tol)
    call work_fields_put(w, w_idx, w_scratch)
  end subroutine krylov_check_residual_scalar

  !> Recompute the true residuals of a vector solve, with the vector
  !! boundary-condition projector applied to the three components together.
  !! @param names Names of the three components, for the log.
  !! @param results The three solver monitors.
  !! @param tol The solver's absolute tolerance.
  subroutine krylov_check_residual_vector(names, Ax, x, y, z, fx, fy, fz, &
       coef, gs, bclst, results, tol)
    character(len=*), intent(in) :: names(3)
    class(ax_t), intent(inout) :: Ax
    type(field_t), intent(inout) :: x, y, z, fx, fy, fz
    type(coef_t), intent(inout) :: coef
    type(gs_t), intent(inout) :: gs
    class(vector_bc_projector_t), intent(inout) :: bclst
    type(ksp_monitor_t), intent(in) :: results(3)
    real(kind=rp), intent(in) :: tol
    type(field_ptr_t) :: w(3)
    type(field_t), pointer :: rx, ry, rz
    integer :: w_idx(3), n
    logical :: w_scratch

    n = coef%dof%size()
    call work_fields_get(w, w_idx, coef, w_scratch)
    rx => w(1)%ptr
    ry => w(2)%ptr
    rz => w(3)%ptr
    call Ax%compute(rx%x, x%x, coef, coef%msh, coef%Xh)
    call Ax%compute(ry%x, y%x, coef, coef%msh, coef%Xh)
    call Ax%compute(rz%x, z%x, coef, coef%msh, coef%Xh)
    call gs%op(rx, GS_OP_ADD)
    call gs%op(ry, GS_OP_ADD)
    call gs%op(rz, GS_OP_ADD)
    call bclst%apply(rx%x, ry%x, rz%x, n)
    call field_sub2(rx, fx, n)
    call field_sub2(ry, fy, n)
    call field_sub2(rz, fz, n)
    call report_residual(names(1), residual_norm(rx, coef), results(1), tol)
    call report_residual(names(2), residual_norm(ry, coef), results(2), tol)
    call report_residual(names(3), residual_norm(rz, coef), results(3), tol)
    call work_fields_put(w, w_idx, w_scratch)
  end subroutine krylov_check_residual_vector

  !> Take `size(w)` work fields: from the scratch registry when it is set up
  !! on the dofmap of @a coef (the time loop allocates and reuses these, so
  !! nothing is added), otherwise allocated here.
  subroutine work_fields_get(w, idx, coef, scratch)
    type(field_ptr_t), intent(inout) :: w(:)
    integer, intent(out) :: idx(:)
    type(coef_t), intent(in) :: coef
    logical, intent(out) :: scratch
    integer :: i

    scratch = associated(neko_scratch_registry%dof, coef%dof)
    do i = 1, size(w)
       if (scratch) then
          call neko_scratch_registry%request_field(w(i)%ptr, idx(i), .false.)
       else
          idx(i) = 0
          allocate(w(i)%ptr)
          call w(i)%ptr%init(coef%dof)
       end if
    end do
  end subroutine work_fields_get

  !> Return the work fields taken by work_fields_get.
  subroutine work_fields_put(w, idx, scratch)
    type(field_ptr_t), intent(inout) :: w(:)
    integer, intent(inout) :: idx(:)
    logical, intent(in) :: scratch
    integer :: i

    do i = 1, size(w)
       if (scratch) then
          call neko_scratch_registry%relinquish_field(idx(i))
       else
          call w(i)%ptr%free()
          deallocate(w(i)%ptr)
       end if
       nullify(w(i)%ptr)
    end do
  end subroutine work_fields_put

  !> The solver's residual norm, \f$ (\sum_i m_i r_i^2 / V)^{1/2} \f$.
  function residual_norm(r, coef) result(res)
    type(field_t), intent(in) :: r
    type(coef_t), intent(in) :: coef
    real(kind=rp) :: res
    integer :: n

    n = coef%dof%size()
    if (NEKO_BCKND_DEVICE .eq. 1) then
       res = sqrt(device_glsc3(r%x_d, coef%mult_d, r%x_d, n) / coef%volume)
    else
       res = sqrt(glsc3(r%x, coef%mult, r%x, n) / coef%volume)
    end if
  end function residual_norm

  !> Log the true and the reported residual, and warn when they differ.
  subroutine report_residual(name, true_res, results, tol)
    character(len=*), intent(in) :: name
    real(kind=rp), intent(in) :: true_res
    type(ksp_monitor_t), intent(in) :: results
    real(kind=rp), intent(in) :: tol
    real(kind=rp) :: reported, ref
    character(len=LOG_SIZE) :: log_buf
    character(len=256) :: warn_buf

    reported = results%res_final
    write(log_buf, '(A,A,ES12.4,A,ES12.4)') trim(name), &
         ' true residual: ', true_res, ', reported: ', reported
    call neko_log%message(log_buf)

    ! A difference matters relative to the larger of the two, and not when
    ! both are well below the tolerance or at round-off of the right-hand
    ! side (a zero right-hand side gives a reported residual of zero)
    ref = max(reported, true_res, 0.1_rp * tol, &
         1.0e3_rp * NEKO_EPS * results%res_start)
    if (abs(true_res - reported) .gt. RESIDUAL_WARN * ref) then
       write(warn_buf, '(A,A,ES10.3,A,ES10.3,A)') trim(name), &
            ' solver: the true residual ', true_res, &
            ' differs from the reported one ', reported, &
            '; the stopping test does not measure the true residual'
       call neko_warning(trim(warn_buf))
    end if
  end subroutine report_residual

end module krylov_check
