!> Krylov preconditioner
module jacobi
  use math
  use utils
  use precon
  use ax_product
  use gather_scatter
  implicit none
  
  !> Defines a canonical Krylov preconditioner
  type, public, extends(pc_t) :: jacobi_t
     real(kind=rp), allocatable :: d(:,:,:,:)
     type(gs_t), pointer :: gs_h
     type(dofmap_t), pointer :: dof
     type(coef_t), pointer :: coef
  contains
     procedure, pass(this) :: init => jacobi_init
     procedure, pass(this) :: free => jacobi_free
     procedure, pass(this) :: solve => jacobi_solve
     procedure, pass(this) :: update => jacobi_update
  end type jacobi_t

  !> Abstract interface for solving \f$ M z = r \f$
  !!
  !! @param z vector of length @a n
  !! @param r vector of length @a n
contains
  subroutine jacobi_init(this, coef, dof, gs_h)
    class(jacobi_t), intent(inout) :: this
    type(coef_t), intent(inout), target :: coef
    type(dofmap_t), intent(inout), target :: dof
    type(gs_t), intent(inout), target :: gs_h
    
    call this%free()
    this%gs_h => gs_h
    this%dof => dof
    this%coef => coef
    allocate(this%d(dof%Xh%lx,dof%Xh%ly,dof%Xh%lz, dof%msh%nelv))
    call jacobi_update(this)

  end subroutine jacobi_init

  subroutine jacobi_free(this)
    class(jacobi_t), intent(inout) :: this
    if (allocated(this%d)) then
      deallocate(this%d)
    end if
    nullify(this%dof)
    nullify(this%gs_h)
    nullify(this%coef)
  end subroutine jacobi_free

  !> The jacobi preconditioner \f$ J z = r \f$
  !! \f$ z = J^{-1}r\f$ where \f$ J^{-1} ~= 1/diag(A) \f$
  subroutine jacobi_solve(this, z, r, n)
    integer, intent(inout) :: n
    class(jacobi_t), intent(inout) :: this
    real(kind=rp), dimension(n), intent(inout) :: z
    real(kind=rp), dimension(n), intent(inout) :: r
    call col3(z,r,this%d,n)
  end subroutine jacobi_solve



  subroutine jacobi_update(this)
    class(jacobi_t), intent(inout) :: this
    integer :: i, j, k, l, e, lz, ly, lx
    associate(dof => this%dof, coef => this%coef, gs_h => this%gs_h)

    lx = dof%Xh%lx
    ly = dof%Xh%ly
    lz = dof%Xh%lz

    this%d = 0d0
    
    do e=1,dof%msh%nelv
      do l=1,lx
         do k=1,lz
            do j=1,ly
               do i=1,lx
                  this%d(i,j,k,e) = this%d(i,j,k,e) + &
                                    coef%G1(l,j,k,e) * dof%Xh%dxt(i,l)**2
               end do
            end do
         end do
      end do
      do l=1,ly
         do k=1,lz
            do j=1,ly
               do i=1,lx
                  this%d(i,j,k,e) = this%d(i,j,k,e) + &
                                    coef%G2(i,l,k,e) * dof%Xh%dyt(j,l)**2
               end do
            end do
         end do
      end do
      do l=1,lz
         do k=1,lz
            do j=1,ly
               do i=1,lx
                  this%d(i,j,k,e) = this%d(i,j,k,e) + &
                                    coef%G3(i,j,l,e) * dof%Xh%dzt(k,l)**2
               end do
            end do
         end do
     end do

     if (dof%msh%dfrmd_el(e)) then
        do j=1,ly,ly-1
           do k=1,lz,lz-1
              this%d(1,j,k,e) = this%d(1,j,k,e) &
                              + coef%G4(1,j,k,e) * dof%Xh%dxt(1,1)*dof%Xh%dyt(j,j) &
                              + coef%G5(1,j,k,e) * dof%Xh%dxt(1,1)*dof%Xh%dzt(k,k)
              this%d(lx,j,k,e) = this%d(lx,j,k,e) &
                               + coef%G4(lx,j,k,e) * dof%Xh%dxt(lx,lx)*dof%Xh%dyt(j,j) &
                               + coef%G5(lx,j,k,e) * dof%Xh%dxt(lx,lx)*dof%Xh%dzt(k,k)
           end do
        end do

        do i=1,lx,lx-1
           do k=1,lz,lz-1
              this%d(i,1,k,e) = this%d(i,1,k,e) &
                              + coef%G4(i,1,k,e) * dof%Xh%dyt(1,1)*dof%Xh%dxt(i,i) &
                              + coef%G6(i,1,k,e) * dof%Xh%dyt(1,1)*dof%Xh%dzt(k,k)
              this%d(i,ly,k,e) = this%d(i,ly,k,e) &
                               + coef%G4(i,ly,k,e) * dof%Xh%dyt(ly,ly)*dof%Xh%dxt(i,i) &
                               + coef%G6(i,ly,k,e) * dof%Xh%dyt(ly,ly)*dof%Xh%dzt(k,k)
           end do
        end do
        do i=1,lx,lx-1
           do j=1,ly,ly-1
              this%d(i,j,1,e) = this%d(i,j,1,e) &
                              + coef%G5(i,j,1,e) * dof%Xh%dzt(1,1)*dof%Xh%dxt(i,i) &
                              + coef%G6(i,j,1,e) * dof%Xh%dzt(1,1)*dof%Xh%dyt(j,j)
              this%d(i,j,lz,e) = this%d(i,j,lz,e) &
                               + coef%G5(i,j,lz,e) * dof%Xh%dzt(lz,lz)*dof%Xh%dxt(i,i) &
                               + coef%G6(i,j,lz,e) * dof%Xh%dzt(lz,lz)*dof%Xh%dyt(j,j)
           end do
        end do
     end if
   end do 
   call col2(this%d,coef%h1,coef%dof%n_dofs)
   if (coef%ifh2) call addcol3(this%d,coef%h2,coef%B,coef%dof%n_dofs)
   call gs_op_vector(gs_h, this%d, dof%n_dofs, GS_OP_ADD)
   if (.not. coef%ifh2) call col2(this%d, coef%mult, coef%dof%n_dofs)
   call invcol1(this%d,dof%n_dofs)
   end associate
  end subroutine jacobi_update
  
 end module jacobi
