!> Jacobi preconditioner accelerator backend
module device_jacobi
  use precon
  use coefs
  use num_types
  use device_math
  use gather_scatter
  use, intrinsic :: iso_c_binding
  implicit none
  private

  !> Defines a jacobi preconditioner
  type, public, extends(pc_t) :: device_jacobi_t
     real(kind=rp), allocatable :: d(:,:,:,:)
     type(c_ptr) :: d_d
     type(gs_t), pointer :: gs_h
     type(dofmap_t), pointer :: dof
     type(coef_t), pointer :: coef
   contains
     procedure, pass(this) :: init => device_jacobi_init
     procedure, pass(this) :: free => device_jacobi_free
     procedure, pass(this) :: solve => device_jacobi_solve
     procedure, pass(this) :: update => device_jacobi_update
  end type device_jacobi_t

  interface
     subroutine hip_jacobi_update(d_d, dxt_d, dyt_d, dzt_d, &
          G11_d, G22_D, G33_d, G12_d, G13_d, G23_d, nelv, lx) &
          bind(c, name='hip_jacobi_update')
       use, intrinsic :: iso_c_binding
       type(c_ptr), value :: d_d, dxt_d, dyt_d, dzt_d
       type(c_ptr), value :: G11_d, G22_D, G33_d, G12_d, G13_d, G23_d
       integer(c_int) :: nelv, lx
     end subroutine hip_jacobi_update
  end interface
  
contains
  
  subroutine device_jacobi_init(this, coef, dof, gs_h)
    class(device_jacobi_t), intent(inout) :: this
    type(coef_t), intent(inout), target :: coef
    type(dofmap_t), intent(inout), target :: dof
    type(gs_t), intent(inout), target :: gs_h

    call this%free()

    this%gs_h => gs_h
    this%dof => dof
    this%coef => coef

    allocate(this%d(dof%Xh%lx,dof%Xh%ly,dof%Xh%lz, dof%msh%nelv))

    call device_map(this%d, this%d_d, size(this%d))

    call device_jacobi_update(this)

  end subroutine device_jacobi_init

  subroutine device_jacobi_free(this)
    class(device_jacobi_t), intent(inout) :: this

    if (c_associated(this%d_d)) then
       call device_free(this%d_d)
       this%d_d = C_NULL_PTR
    end if
    
    if (allocated(this%d)) then
       deallocate(this%d)
    end if

    nullify(this%dof)
    nullify(this%gs_h)
    nullify(this%coef)
  end subroutine device_jacobi_free

  !> The jacobi preconditioner \f$ J z = r \f$
  !! \f$ z = J^{-1}r\f$ where \f$ J^{-1} ~= 1/diag(A) \f$
  subroutine device_jacobi_solve(this, z, r, n)
    integer, intent(inout) :: n
    class(device_jacobi_t), intent(inout) :: this
    real(kind=rp), dimension(n), intent(inout) :: z
    real(kind=rp), dimension(n), intent(inout) :: r
    type(c_ptr) :: z_d, r_d
    
    z_d = device_get_ptr(z, n)
    r_d = device_get_ptr(r, n)
    
    call device_col3(z_d, r_d, this%d_d, n)
    
  end subroutine device_jacobi_solve

  subroutine device_jacobi_update(this)
    class(device_jacobi_t), intent(inout) :: this
    integer :: i, j, k, l, e, lz, ly, lx
    associate(dof => this%dof, coef => this%coef, Xh => this%dof%Xh, &
         gs_h => this%gs_h, nelv => this%dof%msh%nelv)

      lx = Xh%lx
      ly = Xh%ly
      lz = Xh%lz

      call hip_jacobi_update(this%d_d, Xh%dxt_d, Xh%dyt_d, Xh%dzt_d, &
                             coef%G11_d, coef%G22_d, coef%G33_d, &
                             coef%G12_d, coef%G13_d, coef%G23_d, &
                             nelv, lx)

      call device_col2(this%d_d, coef%h1_d, coef%dof%n_dofs)

      if (coef%ifh2) then
         call device_addcol3(this%d_d, coef%h2_d, coef%B_d, coef%dof%n_dofs)
      end if
      
      call gs_op_vector(gs_h, this%d, dof%n_dofs, GS_OP_ADD)

      if (.not. coef%ifh2) then
         call device_col2(this%d_d, coef%mult_d, coef%dof%n_dofs)
      end if
      
      call invcol1(this%d, dof%n_dofs)
    end associate
  end subroutine device_jacobi_update

end module device_jacobi
