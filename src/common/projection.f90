! Copyright (c) 2008-2020, UCHICAGO ARGONNE, LLC. 
!
! The UChicago Argonne, LLC as Operator of Argonne National
! Laboratory holds copyright in the Software. The copyright holder
! reserves all rights except those expressly granted to licensees,
! and U.S. Government license rights.
!
! Redistribution and use in source and binary forms, with or without
! modification, are permitted provided that the following conditions
! are met:
!
! 1. Redistributions of source code must retain the above copyright
! notice, this list of conditions and the disclaimer below.
!
! 2. Redistributions in binary form must reproduce the above copyright
! notice, this list of conditions and the disclaimer (as noted below)
! in the documentation and/or other materials provided with the
! distribution.
!
! 3. Neither the name of ANL nor the names of its contributors
! may be used to endorse or promote products derived from this software
! without specific prior written permission.
!
! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS 
! "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT 
! LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS 
! FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL 
! UCHICAGO ARGONNE, LLC, THE U.S. DEPARTMENT OF 
! ENERGY OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, 
! SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED 
! TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, 
! DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY 
! THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT 
! (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE 
! OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
!
! Additional BSD Notice
! ---------------------
! 1. This notice is required to be provided under our contract with
! the U.S. Department of Energy (DOE). This work was produced at
! Argonne National Laboratory under Contract 
! No. DE-AC02-06CH11357 with the DOE.
!
! 2. Neither the United States Government nor UCHICAGO ARGONNE, 
! LLC nor any of their employees, makes any warranty, 
! express or implied, or assumes any liability or responsibility for the
! accuracy, completeness, or usefulness of any information, apparatus,
! product, or process disclosed, or represents that its use would not
! infringe privately-owned rights.
!
! 3. Also, reference herein to any specific commercial products, process, 
! or services by trade name, trademark, manufacturer or otherwise does 
! not necessarily constitute or imply its endorsement, recommendation, 
! or favoring by the United States Government or UCHICAGO ARGONNE LLC. 
! The views and opinions of authors expressed 
! herein do not necessarily state or reflect those of the United States 
! Government or UCHICAGO ARGONNE, LLC, and shall 
! not be used for advertising or product endorsement purposes.
!
!> Project x onto X, the space of old solutions and back again
!! @note In this code we assume that the matrix project for the
!! pressure Ax does not vary in time.
module projection
  use num_types
  use math
  use coefs
  use ax_product
  use bc
  use comm
  use gather_scatter
  use neko_config
  use device
  use device_math
  use, intrinsic :: iso_c_binding
  implicit none
  private

  type, public ::  projection_t
     real(kind=rp), allocatable :: xx(:,:)
     real(kind=rp), allocatable :: bb(:,:)
     real(kind=rp), allocatable :: xbar(:)
     type(c_ptr), allocatable :: xx_d(:)
     type(c_ptr), allocatable :: bb_d(:)
     type(c_ptr) :: xbar_d = C_NULL_PTR
     type(c_ptr) :: alpha_d = C_NULL_PTR
     type(c_ptr) :: xx_d_d = C_NULL_PTR
     type(c_ptr) :: bb_d_d = C_NULL_PTR
     integer :: m, L
     real(kind=rp) :: tol = 1d-7
   contains
     procedure, pass(this) :: project_on => bcknd_project_on
     procedure, pass(this) :: project_back => bcknd_project_back
     procedure, pass(this) :: init => projection_init
     procedure, pass(this) :: free => projection_free
  end type projection_t

contains

  subroutine projection_init(this, n, L)
    class(projection_t), target, intent(inout) :: this
    integer, intent(in) :: n
    integer, optional, intent(in) :: L
    integer :: i
    integer(c_size_t) :: ptr_size
    type(c_ptr) :: ptr
    
    call this%free()
    
    if (present(L)) then
       this%L = L
    else
       this%L = 20
    end if

    this%m = 0 

    allocate(this%xx(n,this%L))
    allocate(this%bb(n,this%L))
    allocate(this%xbar(n))
    allocate(this%xx_d(this%L))
    allocate(this%bb_d(this%L))
    call rzero(this%xbar,n)
    do i = 1, this%L
       call rzero(this%xx(1,i),n)
       call rzero(this%bb(1,i),n)
    end do
    if ((NEKO_BCKND_HIP .eq. 1) .or. (NEKO_BCKND_CUDA .eq. 1) .or. &
         (NEKO_BCKND_OPENCL .eq. 1)) then
       
       call device_map(this%xbar, this%xbar_d,n)
       call device_alloc(this%alpha_d, int(rp*this%L,c_size_t))

       do i = 1, this%L
          this%xx_d(i) = C_NULL_PTR
          call device_map_r1(this%xx(:,i), this%xx_d(i), n)
          this%bb_d(i) = C_NULL_PTR
          call device_map_r1(this%bb(:,i), this%bb_d(i), n)
       end do

       ptr_size = 8*this%L
       call device_alloc(this%xx_d_d, ptr_size)
       ptr = c_loc(this%xx_d)
       call device_memcpy(ptr,this%xx_d_d, ptr_size, HOST_TO_DEVICE)
       call device_alloc(this%bb_d_d, ptr_size)
       ptr = c_loc(this%bb_d)
       call device_memcpy(ptr,this%bb_d_d, ptr_size, HOST_TO_DEVICE)
    end if


  end subroutine projection_init

  subroutine projection_free(this)
    class(projection_t), intent(inout) :: this
    integer :: i
    if (allocated(this%xx)) then
       deallocate(this%xx)
    end if
    if (allocated(this%bb)) then
       deallocate(this%bb)
    end if
    if (allocated(this%xbar)) then
       deallocate(this%xbar)
    end if
    if (allocated(this%xx_d)) then
       do i = 1, this%L
          if (c_associated(this%xx_d(i))) then
             call device_free(this%xx_d(i))
          end if
       end do
    end if
    if (c_associated(this%xx_d_d)) then
       call device_free(this%xx_d_d)
    end if
    if (c_associated(this%xbar_d)) then
       call device_free(this%xbar_d)
    end if
    if (c_associated(this%alpha_d)) then
       call device_free(this%alpha_d)
    end if
    if (allocated(this%bb_d)) then
       do i = 1, this%L
          if (c_associated(this%bb_d(i))) then
             call device_free(this%bb_d(i))
          end if
       end do
    end if
    if (c_associated(this%bb_d_d)) then
       call device_free(this%bb_d_d)
    end if

  end subroutine projection_free

  subroutine bcknd_project_on(this, b, coef, n)
    class(projection_t), intent(inout) :: this
    integer, intent(inout) :: n
    class(coef_t), intent(inout) :: coef   
    real(kind=rp), intent(inout), dimension(n) :: b 
    if ((NEKO_BCKND_HIP .eq. 1) .or. (NEKO_BCKND_CUDA .eq. 1) .or. &
         (NEKO_BCKND_OPENCL .eq. 1)) then
       call device_project_on(this, b, coef, n)
    else
       call cpu_project_on(this, b, coef, n)
    end if
  end subroutine bcknd_project_on
  subroutine bcknd_project_back(this,x,Ax,coef, bclst, gs_h, n)
    class(projection_t) :: this
    integer, intent(inout) :: n
    class(Ax_t), intent(inout) :: Ax    
    class(coef_t), intent(inout) :: coef   
    class(bc_list_t), intent(inout) :: bclst
    type(gs_t), intent(inout) :: gs_h
    real(kind=rp), intent(inout), dimension(n) :: x 
    type(c_ptr) :: x_d

    this%m = min(this%m+1,this%L)
    
    if ((NEKO_BCKND_HIP .eq. 1) .or. (NEKO_BCKND_CUDA .eq. 1) .or. &
         (NEKO_BCKND_OPENCL .eq. 1)) then
       x_d = device_get_ptr(x,n)
       if (this%m .gt. 0) call device_add2(x_d,this%xbar_d,n)      ! Restore desired solution
       if (this%m .eq. this%L) this%m = 1
       call device_copy(this%xx_d(this%m),x_d,n)   ! Update (X,B)

    else
       if (this%m.gt.0) call add2(x,this%xbar,n)      ! Restore desired solution
       call copy        (this%xx(1,this%m),x,n)   ! Update (X,B)
    end if

    call Ax%compute(this%bb(1,this%m), x, coef, coef%msh, coef%Xh)
    call gs_op_vector(gs_h, this%bb(1,this%m), n, GS_OP_ADD)
    call bc_list_apply_scalar(bclst, this%bb(1,this%m), n)

    if ((NEKO_BCKND_HIP .eq. 1) .or. (NEKO_BCKND_CUDA .eq. 1) .or. &
         (NEKO_BCKND_OPENCL .eq. 1)) then
       call device_proj_ortho(this, this%xx_d, this%bb_d, coef%mult_d, n)

    else
       call cpu_proj_ortho  (this,this%xx,this%bb,coef%mult,n) 
    end if
  end subroutine bcknd_project_back
 


  subroutine cpu_project_on(this, b, coef, n)
    class(projection_t), intent(inout) :: this
    integer, intent(inout) :: n
    class(coef_t), intent(inout) :: coef   
    real(kind=rp), intent(inout), dimension(n) :: b 
    integer :: i, j, k, ierr
    real(kind=rp) :: work(this%L), alpha(this%L)

    associate(xbar => this%xbar, xx => this%xx, &
              bb => this%bb)
    
      if (this%m .le. 0) return

      !First round of CGS
      call rzero(alpha, this%m)
      do i = 1, n, NEKO_BLK_SIZE
         j = min(NEKO_BLK_SIZE, n-i+1)
         do k = 1, this%m 
            alpha(k) = alpha(k) + vlsc3(xx(i,k), coef%mult(i,1,1,1), b(i), j)
         end do
      end do
      
      !First one outside loop to avoid zeroing xbar and bbar
      call MPI_Allreduce(MPI_IN_PLACE, alpha, this%m, &
           MPI_REAL_PRECISION, MPI_SUM, NEKO_COMM, ierr)
      
      call rzero(work, this%m)
      
      do i = 1, n, NEKO_BLK_SIZE
         j = min(NEKO_BLK_SIZE, n-i+1)
         call cmult2(xbar(i), xx(i,1), alpha(1), j)
         call add2s2(b(i), bb(i,1), -alpha(1), j)
         do k = 2,this%m
            call add2s2(xbar(i), xx(i,k), alpha(k), j)
            call add2s2(b(i), bb(i,k), -alpha(k), j)
         end do
         !Second round of CGS
         do k = 1, this%m
            work(k) = work(k) + vlsc3(xx(i,k), coef%mult(i,1,1,1), b(i), j)
         end do
      end do
      
      call MPI_Allreduce(work, alpha, this%m, &
           MPI_REAL_PRECISION, MPI_SUM, NEKO_COMM, ierr)
      
      do i = 1, n, NEKO_BLK_SIZE
         j = min(NEKO_BLK_SIZE, n-i+1)
         do k = 1,this%m
            call add2s2(xbar(i), xx(i,k), alpha(k), j)
            call add2s2(b(i), bb(i,k), -alpha(k), j)
         end do
      end do
    end associate
  end subroutine cpu_project_on

  subroutine device_project_on(this, b, coef, n)
    class(projection_t), intent(inout) :: this
    integer, intent(inout) :: n
    class(coef_t), intent(inout) :: coef   
    real(kind=rp), intent(inout), dimension(n) :: b 
    integer :: i, j, k, ierr
    real(kind=rp) :: work(this%L), alpha(this%L)
    type(c_ptr) :: b_d
    b_d = device_get_ptr(b, n)

    associate(xbar_d => this%xbar_d, xx_d => this%xx_d, xx_d_d => this%xx_d_d, &
              bb_d => this%bb_d, bb_d_d => this%bb_d_d, alpha_d => this%alpha_d)
    
      if (this%m .le. 0) return


      call device_glsc3_many(alpha,b_d,xx_d_d,coef%mult_d,this%m,n)
      call device_memcpy(alpha, alpha_d, this%m, HOST_TO_DEVICE) 
      call device_rzero(xbar_d, n)
      call device_add2s2_many(xbar_d, xx_d_d, alpha_d, this%m, n)
      call device_cmult(alpha_d, -1.0_rp, this%m)
      call device_add2s2_many(b_d, bb_d_d, alpha_d, this%m, n)
      call device_glsc3_many(alpha,b_d,xx_d_d,coef%mult_d,this%m,n)
 
      call device_memcpy(alpha, alpha_d, this%m, HOST_TO_DEVICE) 
      call device_add2s2_many(xbar_d, xx_d_d, alpha_d, this%m, n)
      call device_cmult(alpha_d, -1.0_rp, this%m)
      call device_add2s2_many(b_d, bb_d_d, alpha_d, this%m, n)
      
    end associate
  end subroutine device_project_on

  !This is a lot more primitive than on the CPU
  subroutine device_proj_ortho(this, xx_d, bb_d, w_d, n)
    type(projection_t)  :: this
    integer, intent(inout) :: n
    type(c_ptr), dimension(this%L) :: xx_d, bb_d
    type(c_ptr) :: w_d
    real(kind=rp) :: nrm, scl
    real(kind=rp) :: alpha(this%L)

    associate(m => this%m,  xx_d_d => this%xx_d_d, &
              bb_d_d => this%bb_d_d, alpha_d => this%alpha_d)
      
      if(m .le. 0) return

      call device_glsc3_many(alpha,bb_d(m),xx_d_d,w_d,m,n)
      nrm = sqrt(alpha(m))
      call cmult(alpha, -1.0_rp,m)
      call device_memcpy(alpha, alpha_d, this%m, HOST_TO_DEVICE) 
      call device_add2s2_many(xx_d(m),xx_d_d,alpha_d,m-1,n)
      call device_add2s2_many(bb_d(m),bb_d_d,alpha_d,m-1,n)
    
      call device_glsc3_many(alpha,bb_d(m),xx_d_d,w_d,m,n)
      call cmult(alpha, -1.0_rp,m)
      call device_memcpy(alpha, alpha_d, m, HOST_TO_DEVICE) 
      call device_add2s2_many(xx_d(m),xx_d_d,alpha_d,m-1,n)
      call device_add2s2_many(bb_d(m),bb_d_d,alpha_d,m-1,n)
      call device_glsc3_many(alpha,bb_d(m),xx_d_d,w_d,m,n)
      alpha(m) = device_glsc3(xx_d(m), w_d, bb_d(m), n)
      alpha(m) = sqrt(alpha(m))

      if(alpha(m) .gt. this%tol*nrm) then !New vector is linearly independent    
         scl = 1.0_rp / alpha(m) 
         call device_cmult(xx_d(m), scl, n)   
         call device_cmult(bb_d(m), scl, n)   

        
      else !New vector is not linearly independent, forget about it
         if(pe_rank .eq. 0) then
            call neko_warning('New vector not linearly indepependent!')
         end if
         m = m - 1 !Remove column
      endif

    end associate
    
  end subroutine device_proj_ortho
  
  
  subroutine cpu_proj_ortho(this, xx, bb, w, n)
    type(projection_t)  :: this
    integer, intent(inout) :: n
    real(kind=rp), dimension(n, this%L), intent(inout) :: xx, bb
    real(kind=rp), dimension(n), intent(inout) :: w
    real(kind=rp) :: nrm, scl1, scl2, c, s
    real(kind=rp) :: work(this%L), alpha(this%L), beta(this%L)
    integer :: i, j, k, h, ierr

    associate(m => this%m)
      
      if(m .le. 0) return !No vectors to ortho-normalize 

      ! AX = B
      ! Calculate dx, db: dx = x-XX^Tb, db=b-BX^Tb     
      call rzero(alpha, m)
      do i = 1, n, NEKO_BLK_SIZE
         j = min(NEKO_BLK_SIZE, n-i+1)
         do k = 1, m !First round CGS
            alpha(k) = alpha(k) + 0.5_rp * (vlsc3(xx(i,k), w(i), bb(i,m), j) &
                 + vlsc3(bb(i,k), w(i), xx(i,m), j))
         end do
      end do
      
      call MPI_Allreduce(MPI_IN_PLACE, alpha, this%m, &
           MPI_REAL_PRECISION, MPI_SUM, NEKO_COMM, ierr)
      
      nrm = sqrt(alpha(m)) !Calculate A-norm of new vector
      
      call rzero(beta,m)
    
      do i = 1, n, NEKO_BLK_SIZE
         j = min(NEKO_BLK_SIZE, n-i+1)
         do k = 1,m-1
            call add2s2(xx(i,m), xx(i,k), -alpha(k), j)
            call add2s2(bb(i,m), bb(i,k), -alpha(k), j)
            beta(k) = beta(k) + 0.5_rp * (vlsc3(xx(i,k), w(i), bb(i,m), j) &
                 + vlsc3(bb(i,k), w(i), xx(i,m), j))
         end do
      end do
    
      call MPI_Allreduce(MPI_IN_PLACE, beta, this%m-1, &
           MPI_REAL_PRECISION, MPI_SUM, NEKO_COMM, ierr)

      alpha(m) = 0.0_rp
      
      do i = 1, n, NEKO_BLK_SIZE
         j = min(NEKO_BLK_SIZE,n-i+1)
         do k = 1, m-1
            call add2s2(xx(i,m), xx(i,k), -beta(k), j)
            call add2s2(bb(i,m), bb(i,k), -beta(k), j)
         end do
         alpha(m) = alpha(m) + vlsc3(xx(i,m), w(i), bb(i,m), j)
      end do
      do k = 1, m-1
         alpha(k) = alpha(k) + beta(k)
      end do
      
      !alpha(m) = glsc3(xx(1,m), w, bb(1,m), n) 
      call MPI_Allreduce(MPI_IN_PLACE, alpha(m), 1, &
           MPI_REAL_PRECISION, MPI_SUM, NEKO_COMM, ierr)
      alpha(m) = sqrt(alpha(m))
      !dx and db now stored in last column of xx and bb

      if(alpha(m) .gt. this%tol*nrm) then !New vector is linearly independent    
         !Normalize dx and db
         scl1 = 1.0_rp / alpha(m) 
         call cmult(xx(1,m), scl1, n)   
         call cmult(bb(1,m), scl1, n)   

         !We want to throw away the oldest information
         !The below propagates newest information to first vector.
         !This will make the first vector a scalar 
         !multiple of x.
         do k = m, 2, -1
            h = k - 1   
            call givens_rotation(alpha(h), alpha(k), c, s, nrm)
            alpha(h) = nrm     
            do i = 1, n !Apply rotation to xx and bb
               scl1 = c*xx(i,h) + s*xx(i,k)
               xx(i,k) = -s*xx(i,h) + c*xx(i,k)
               xx(i,h) = scl1       
               scl2 = c*bb(i,h) + s*bb(i,k)
               bb(i,k) = -s*bb(i,h) + c*bb(i,k)    
               bb(i,h) = scl2        
            end do
         end do
         
      else !New vector is not linearly independent, forget about it
         k = m !location of rank deficient column
         if(pe_rank .eq. 0) then
            call neko_warning('New vector not linearly indepependent!')
         end if
         m = m - 1 !Remove column
      endif

    end associate
    
  end subroutine cpu_proj_ortho
  
  subroutine givens_rotation(a, b, c, s, r)
    real(kind=rp), intent(inout) :: a, b, c, s, r
    real(kind=rp) ::  h, d

    if(b .ne. 0.0_rp) then
       h = hypot(a, b) 
       d = 1.0_rp / h
       c = abs(a) * d
       s = sign(d, a) * b
       r = sign(1.0_rp, a) * h
    else
       c = 1.0_rp
       s = 0.0_rp
       r = a
    endif
      
    return
  end subroutine givens_rotation

end module projection
