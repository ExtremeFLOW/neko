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
!> Defines inflow dirichlet conditions
module usr_inflow
  use num_types
  use coefs
  use inflow
  use device
  use device_inhom_dirichlet
  use utils
  implicit none
  private
  
  !> User defined dirichlet condition for inlet (vector valued)
  type, public, extends(inflow_t) :: usr_inflow_t
     type(coef_t), pointer :: c => null()
     procedure(usr_inflow_eval), nopass, pointer :: eval => null()
     procedure(usr_scalar_bc_eval), nopass, pointer :: eval_scalar_bc => null()
     type(c_ptr), private :: usr_x_d = C_NULL_PTR
     type(c_ptr), private :: usr_y_d = C_NULL_PTR
     type(c_ptr), private :: usr_z_d = C_NULL_PTR
   contains
     procedure, pass(this) :: apply_scalar => usr_inflow_apply_scalar
     procedure, pass(this) :: apply_vector => usr_inflow_apply_vector
     procedure, pass(this) :: validate => usr_inflow_validate
     procedure, pass(this) :: set_coef => usr_inflow_set_coef
     procedure, pass(this) :: set_eval => usr_inflow_set_eval
     procedure, pass(this) :: set_scalar_bc => usr_set_scalar_bc
     procedure, pass(this) :: apply_vector_dev => usr_inflow_apply_vector_dev
     procedure, pass(this) :: apply_scalar_dev => usr_inflow_apply_scalar_dev
     final :: usr_inflow_free
  end type usr_inflow_t

  !> Abstract interface defining a user defined inflow condition (pointwise)
  abstract interface
     subroutine usr_inflow_eval(u, v, w, x, y, z, nx, ny, nz, ix, iy, iz, ie)
       import rp
       real(kind=rp), intent(inout) :: u
       real(kind=rp), intent(inout) :: v
       real(kind=rp), intent(inout) :: w
       real(kind=rp), intent(in) :: x
       real(kind=rp), intent(in) :: y
       real(kind=rp), intent(in) :: z
       real(kind=rp), intent(in) :: nx
       real(kind=rp), intent(in) :: ny
       real(kind=rp), intent(in) :: nz
       integer, intent(in) :: ix
       integer, intent(in) :: iy
       integer, intent(in) :: iz
       integer, intent(in) :: ie
     end subroutine usr_inflow_eval
  end interface
  !> Abstract interface defining a user defined scalar boundary condition (pointwise)
  !! JUst imitating inflow for now, but we should update this
  !! Probably passing the whole field, params, coef, etcetc would be good
  abstract interface
     subroutine usr_scalar_bc_eval(s, x, y, z, nx, ny, nz, ix, iy, iz, ie)
       import rp
       real(kind=rp), intent(inout) :: s
       real(kind=rp), intent(in) :: x
       real(kind=rp), intent(in) :: y
       real(kind=rp), intent(in) :: z
       real(kind=rp), intent(in) :: nx
       real(kind=rp), intent(in) :: ny
       real(kind=rp), intent(in) :: nz
       integer, intent(in) :: ix
       integer, intent(in) :: iy
       integer, intent(in) :: iz
       integer, intent(in) :: ie
     end subroutine usr_scalar_bc_eval
  end interface


  public :: usr_inflow_eval, usr_scalar_bc_eval

contains
     
  subroutine usr_inflow_free(this)
    type(usr_inflow_t), intent(inout) :: this

    if (c_associated(this%usr_x_d)) then
       call device_free(this%usr_x_d)
    end if
    
    if (c_associated(this%usr_y_d)) then
       call device_free(this%usr_y_d)
    end if

    if (c_associated(this%usr_z_d)) then
       call device_free(this%usr_z_d)
    end if
    
  end subroutine usr_inflow_free
  
  !> Scalar apply
  !! Just imitating inflow for now, but we should look this over
  subroutine usr_inflow_apply_scalar(this, x, n)
    class(usr_inflow_t), intent(inout) :: this
    integer, intent(in) :: n
    real(kind=rp), intent(inout),  dimension(n) :: x
    integer :: i, m, k, idx(4), facet

    associate(xc => this%c%dof%x, yc => this%c%dof%y, zc => this%c%dof%z, &
         nx => this%c%nx, ny => this%c%ny, nz => this%c%nz, &
         lx => this%c%Xh%lx)
      m = this%msk(0)
      print *, m
      do i = 1, m
         k = this%msk(i)
         facet = this%facet(i)
         idx = nonlinear_index(k, lx, lx, lx)
         select case(facet)
         case(1,2)          
            call this%eval_scalar_bc(x(k), &
                 xc(idx(1), idx(2), idx(3), idx(4)), &
                 yc(idx(1), idx(2), idx(3), idx(4)), &
                 zc(idx(1), idx(2), idx(3), idx(4)), &
                 nx(idx(2), idx(3), facet, idx(4)), &
                 ny(idx(2), idx(3), facet, idx(4)), &
                 nz(idx(2), idx(3), facet, idx(4)), &
                 idx(1), idx(2), idx(3), idx(4))
         case(3,4)
            call this%eval_scalar_bc(x(k), &
                 xc(idx(1), idx(2), idx(3), idx(4)), &
                 yc(idx(1), idx(2), idx(3), idx(4)), &
                 zc(idx(1), idx(2), idx(3), idx(4)), &       
                 nx(idx(1), idx(3), facet, idx(4)), &
                 ny(idx(1), idx(3), facet, idx(4)), &
                 nz(idx(1), idx(3), facet, idx(4)), &
                 idx(1), idx(2), idx(3), idx(4))
         case(5,6)
            call this%eval_scalar_bc(x(k), &
                 xc(idx(1), idx(2), idx(3), idx(4)), &
                 yc(idx(1), idx(2), idx(3), idx(4)), &
                 zc(idx(1), idx(2), idx(3), idx(4)), &                     
                 nx(idx(1), idx(2), facet, idx(4)), &
                 ny(idx(1), idx(2), facet, idx(4)), &
                 nz(idx(1), idx(2), facet, idx(4)), &
                 idx(1), idx(2), idx(3), idx(4))
         end select
      end do
    end associate
  end subroutine usr_inflow_apply_scalar

  !> Apply user defined inflow conditions (vector valued)
  subroutine usr_inflow_apply_vector(this, x, y, z, n)
    class(usr_inflow_t), intent(inout) :: this
    integer, intent(in) :: n
    real(kind=rp), intent(inout),  dimension(n) :: x
    real(kind=rp), intent(inout),  dimension(n) :: y
    real(kind=rp), intent(inout),  dimension(n) :: z
    integer :: i, m, k, idx(4), facet

    associate(xc => this%c%dof%x, yc => this%c%dof%y, zc => this%c%dof%z, &
         nx => this%c%nx, ny => this%c%ny, nz => this%c%nz, &
         lx => this%c%Xh%lx)
      m = this%msk(0)
      do i = 1, m
         k = this%msk(i)
         facet = this%facet(i)
         idx = nonlinear_index(k, lx, lx, lx)
         select case(facet)
         case(1,2)          
            call this%eval(x(k), y(k), z(k), &
                 xc(idx(1), idx(2), idx(3), idx(4)), &
                 yc(idx(1), idx(2), idx(3), idx(4)), &
                 zc(idx(1), idx(2), idx(3), idx(4)), &
                 nx(idx(2), idx(3), facet, idx(4)), &
                 ny(idx(2), idx(3), facet, idx(4)), &
                 nz(idx(2), idx(3), facet, idx(4)), &
                 idx(1), idx(2), idx(3), idx(4))
         case(3,4)
            call this%eval(x(k), y(k), z(k), &
                 xc(idx(1), idx(2), idx(3), idx(4)), &
                 yc(idx(1), idx(2), idx(3), idx(4)), &
                 zc(idx(1), idx(2), idx(3), idx(4)), &       
                 nx(idx(1), idx(3), facet, idx(4)), &
                 ny(idx(1), idx(3), facet, idx(4)), &
                 nz(idx(1), idx(3), facet, idx(4)), &
                 idx(1), idx(2), idx(3), idx(4))
         case(5,6)
            call this%eval(x(k), y(k), z(k), &
                 xc(idx(1), idx(2), idx(3), idx(4)), &
                 yc(idx(1), idx(2), idx(3), idx(4)), &
                 zc(idx(1), idx(2), idx(3), idx(4)), &                     
                 nx(idx(1), idx(2), facet, idx(4)), &
                 ny(idx(1), idx(2), facet, idx(4)), &
                 nz(idx(1), idx(2), facet, idx(4)), &
                 idx(1), idx(2), idx(3), idx(4))
         end select
      end do
    end associate
    
  end subroutine usr_inflow_apply_vector

  subroutine usr_inflow_apply_vector_dev(this, x_d, y_d, z_d)
    class(usr_inflow_t), intent(inout), target :: this
    type(c_ptr) :: x_d
    type(c_ptr) :: y_d
    type(c_ptr) :: z_d
    integer :: i, m, k, idx(4), facet
    integer(c_size_t) :: s
    real(kind=rp), allocatable :: x(:)
    real(kind=rp), allocatable :: y(:)
    real(kind=rp), allocatable :: z(:)

    associate(xc => this%c%dof%x, yc => this%c%dof%y, zc => this%c%dof%z, &
         nx => this%c%nx, ny => this%c%ny, nz => this%c%nz, &
         lx => this%c%Xh%lx, usr_x_d => this%usr_x_d, usr_y_d => this%usr_y_d, &
         usr_z_d => this%usr_z_d)

      m = this%msk(0)


      ! Pretabulate values during first call to apply
      if (.not. c_associated(usr_x_d)) then
         allocate(x(m), y(m), z(m)) ! Temp arrays
         
         s = m*rp

         call device_alloc(usr_x_d, s)
         call device_alloc(usr_y_d, s)
         call device_alloc(usr_z_d, s)

         associate(xc => this%c%dof%x, yc => this%c%dof%y, zc => this%c%dof%z, &
                   nx => this%c%nx, ny => this%c%ny, nz => this%c%nz, &
                   lx => this%c%Xh%lx)
           do i = 1, m
              k = this%msk(i)
              facet = this%facet(i)
              idx = nonlinear_index(k, lx, lx, lx)
              select case(facet)
              case(1,2)          
                 call this%eval(x(i), y(i), z(i), &
                      xc(idx(1), idx(2), idx(3), idx(4)), &
                      yc(idx(1), idx(2), idx(3), idx(4)), &
                      zc(idx(1), idx(2), idx(3), idx(4)), &
                      nx(idx(2), idx(3), facet, idx(4)), &
                      ny(idx(2), idx(3), facet, idx(4)), &
                      nz(idx(2), idx(3), facet, idx(4)), &
                      idx(1), idx(2), idx(3), idx(4))
              case(3,4)
                 call this%eval(x(i), y(i), z(i), &
                      xc(idx(1), idx(2), idx(3), idx(4)), &
                      yc(idx(1), idx(2), idx(3), idx(4)), &
                      zc(idx(1), idx(2), idx(3), idx(4)), &       
                      nx(idx(1), idx(3), facet, idx(4)), &
                      ny(idx(1), idx(3), facet, idx(4)), &
                      nz(idx(1), idx(3), facet, idx(4)), &
                      idx(1), idx(2), idx(3), idx(4))
              case(5,6)
                 call this%eval(x(i), y(i), z(i), &
                      xc(idx(1), idx(2), idx(3), idx(4)), &
                      yc(idx(1), idx(2), idx(3), idx(4)), &
                      zc(idx(1), idx(2), idx(3), idx(4)), &                     
                      nx(idx(1), idx(2), facet, idx(4)), &
                      ny(idx(1), idx(2), facet, idx(4)), &
                      nz(idx(1), idx(2), facet, idx(4)), &
                      idx(1), idx(2), idx(3), idx(4))
              end select
           end do
         end associate
 
        
         call device_memcpy(x, usr_x_d, m, HOST_TO_DEVICE)
         call device_memcpy(y, usr_y_d, m, HOST_TO_DEVICE)
         call device_memcpy(z, usr_z_d, m, HOST_TO_DEVICE)

         deallocate(x, y, z)
      end if

      call device_inhom_dirichlet_apply_vector(this%msk_d, x_d, y_d, z_d, &
           usr_x_d, usr_y_d, usr_z_d, m)
      
    end associate

  end subroutine usr_inflow_apply_vector_dev

  !> No-op scalar apply (device version)
  subroutine usr_inflow_apply_scalar_dev(this, x_d)
    class(usr_inflow_t), intent(inout), target :: this
    type(c_ptr) :: x_d
    integer :: i, m, k, idx(4), facet
    integer(c_size_t) :: s
    real(kind=rp), allocatable :: x(:)
      m = this%msk(0)
      print *, m, m*rp

    associate(xc => this%c%dof%x, yc => this%c%dof%y, zc => this%c%dof%z, &
         nx => this%c%nx, ny => this%c%ny, nz => this%c%nz, &
         lx => this%c%Xh%lx, usr_x_d => this%usr_x_d)


      ! Pretabulate values during first call to apply
      if (.not. c_associated(usr_x_d)) then
         allocate(x(m)) ! Temp arrays
         s = m*rp

         call device_alloc(this%usr_x_d, s)

         do i = 1, m
            k = this%msk(i)
            facet = this%facet(i)
            idx = nonlinear_index(k, lx, lx, lx)
            select case(facet)
            case(1,2)          
               call this%eval_scalar_bc(x(i), &
                    xc(idx(1), idx(2), idx(3), idx(4)), &
                    yc(idx(1), idx(2), idx(3), idx(4)), &
                    zc(idx(1), idx(2), idx(3), idx(4)), &
                    nx(idx(2), idx(3), facet, idx(4)), &
                    ny(idx(2), idx(3), facet, idx(4)), &
                    nz(idx(2), idx(3), facet, idx(4)), &
                    idx(1), idx(2), idx(3), idx(4))
            case(3,4)
               call this%eval_scalar_bc(x(i), &
                    xc(idx(1), idx(2), idx(3), idx(4)), &
                    yc(idx(1), idx(2), idx(3), idx(4)), &
                    zc(idx(1), idx(2), idx(3), idx(4)), &       
                    nx(idx(1), idx(3), facet, idx(4)), &
                    ny(idx(1), idx(3), facet, idx(4)), &
                    nz(idx(1), idx(3), facet, idx(4)), &
                    idx(1), idx(2), idx(3), idx(4))
            case(5,6)
               call this%eval_scalar_bc(x(i), &
                    xc(idx(1), idx(2), idx(3), idx(4)), &
                    yc(idx(1), idx(2), idx(3), idx(4)), &
                    zc(idx(1), idx(2), idx(3), idx(4)), &                     
                    nx(idx(1), idx(2), facet, idx(4)), &
                    ny(idx(1), idx(2), facet, idx(4)), &
                    nz(idx(1), idx(2), facet, idx(4)), &
                    idx(1), idx(2), idx(3), idx(4))
            end select
         end do
 
        
         call device_memcpy(x, this%usr_x_d, m, HOST_TO_DEVICE)

         deallocate(x)
      end if

      call device_inhom_dirichlet_apply_scalar(this%msk_d, x_d, &
                                               this%usr_x_d, m)
    end associate


  end subroutine usr_inflow_apply_scalar_dev
  


  !> Assign coefficients (facet normals etc)
  subroutine usr_inflow_set_coef(this, c)
    class(usr_inflow_t), intent(inout) :: this
    type(coef_t), target, intent(inout) :: c
    this%c => c
  end subroutine usr_inflow_set_coef

  !> Assign user provided eval function
  subroutine usr_inflow_set_eval(this, usr_eval)
    class(usr_inflow_t), intent(inout) :: this
    procedure(usr_inflow_eval) :: usr_eval
    this%eval => usr_eval
  end subroutine usr_inflow_set_eval

  !> Assign user provided eval function
  subroutine usr_set_scalar_bc(this, user_scalar_bc)
    class(usr_inflow_t), intent(inout) :: this
    procedure(usr_scalar_bc_eval) :: user_scalar_bc
    this%eval_scalar_bc => user_scalar_bc
  end subroutine usr_set_scalar_bc


  !> Validate user inflow condition
  subroutine usr_inflow_validate(this)
    class(usr_inflow_t), intent(inout) :: this
    logical :: valid

    valid = .true. ! Assert it's going to be ok...    
    if (.not. associated(this%c)) then
       call neko_warning('Missing coefficients')
       valid = .false.       
    end if

    if (.not. associated(this%eval)) then
       call neko_warning('Missing eval function')
       valid = .false.
    end if

    if (.not. valid) then
       call neko_error('Invalid user defined inflow condition')
    end if
    
  end subroutine usr_inflow_validate
  
end module usr_inflow
