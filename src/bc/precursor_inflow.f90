! Copyright (c) 2021, The Neko Authors
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
!> Defines a dirichlet condition for inlet from a precursor simulation (vector valued)
module precursor_inflow
  use num_types
  use coefs
  use utils
  use inflow
  use device
  use device_inhom_dirichlet
!  use flow_profile
  use, intrinsic :: iso_c_binding
  implicit none
  private

  !> Blasius profile for inlet (vector valued)
  type, public, extends(inflow_t) :: precursor_inflow_t
     type(coef_t), pointer :: c => null()
!     real(kind=rp) :: delta
!     procedure(blasius_profile), nopass, pointer :: bla => null()
     type(c_ptr), private :: blax_d = C_NULL_PTR
     type(c_ptr), private :: blay_d = C_NULL_PTR
     type(c_ptr), private :: blaz_d = C_NULL_PTR
   contains
     procedure, pass(this) :: apply_scalar => precursor_inflow_apply_scalar
     procedure, pass(this) :: apply_vector => precursor_inflow_apply_vector
     procedure, pass(this) :: apply_scalar_dev => precursor_inflow_apply_scalar_dev
     procedure, pass(this) :: apply_vector_dev => precursor_inflow_apply_vector_dev
!     procedure, pass(this) :: set_params => blasius_set_params
     procedure, pass(this) :: set_coef => precursor_inflow_set_coef
  end type precursor_inflow_t

contains

  subroutine precursor_inflow_free(this)
    type(precursor_inflow_t), intent(inout) :: this

!    nullify(this%bla)

    if (c_associated(this%blax_d)) then
       call device_free(this%blax_d)
    end if

    if (c_associated(this%blay_d)) then
       call device_free(this%blay_d)
    end if

    if (c_associated(this%blaz_d)) then
       call device_free(this%blaz_d)
    end if
    
  end subroutine precursor_inflow_free
  
  !> No-op scalar apply
  subroutine precursor_inflow_apply_scalar(this, x, n, t, tstep)
    class(precursor_inflow_t), intent(inout) :: this
    integer, intent(in) :: n
    real(kind=rp), intent(inout),  dimension(n) :: x
    real(kind=rp), intent(in), optional :: t
    integer, intent(in), optional :: tstep
  end subroutine precursor_inflow_apply_scalar

  !> No-op scalar apply (device version)
  subroutine precursor_inflow_apply_scalar_dev(this, x_d, t, tstep)
    class(precursor_inflow_t), intent(inout), target :: this
    type(c_ptr) :: x_d
    real(kind=rp), intent(in), optional :: t
    integer, intent(in), optional :: tstep
  end subroutine precursor_inflow_apply_scalar_dev
  
  !> Apply preursor inflow condition (vector valued)
  subroutine precursor_inflow_apply_vector(this, x, y, z, n, t, tstep)
    class(precursor_inflow_t), intent(inout) :: this
    integer, intent(in) :: n
    real(kind=rp), intent(inout),  dimension(n) :: x
    real(kind=rp), intent(inout),  dimension(n) :: y
    real(kind=rp), intent(inout),  dimension(n) :: z
    real(kind=rp), intent(in), optional :: t
    integer, intent(in), optional :: tstep
    integer :: i, m, k, idx(4), facet
    real(kind=rp), allocatable :: yz_map(:)
    integer :: nely, nelz
    integer :: first_gll, last_gll, first_idx(4), last_idx(4)
    integer :: first_facet, last_facet
    real(kind=rp), allocatable :: map(:)  
    real(kind=rp), allocatable :: y_simson(:,:,:), z_simson(:,:,:)
    real(kind=rp) :: min_y,max_y, min_z, max_z  
    integer :: ny1_simson, nz1_simson, nely_simson, nelz_simson
    integer :: l

    associate(xc => this%c%dof%x, yc => this%c%dof%y, zc => this%c%dof%z, &
         nx => this%c%nx, ny => this%c%ny, nz => this%c%nz, &
         lx => this%c%Xh%lx)

!!   Later a better way to automatically get this should be implemented
!    nely = 10 ! Total no. elements in y
!    nelz = 10 ! Total no. elements in z
!     
!    !if (tstep==1)  write(*,*) 'Mask:', this%msk(0), mask_i, size(this%msk)
!    ! Create the map on all processors in the 1st time step. We have tcounter1 as this routine is called per gll point and we only
!    ! want tis part to be called once. 
!    !if (tstep==1.and.tcounter1==0.and.pe_rank==0)  then 
!!    if (tstep==1.and.tcounter1==0)  then 
!    if (tstep==1)  then 
!     allocate(yz_ele_plane(2,nely*nelz)) ! y,z gll coordinates of lower left corner of all element
!     allocate(yz_map(1:nely*nelz))
!     open(unit=10,file='vel_plane_ele',status='old',form='unformatted',access="stream")
!     read(10) yz_ele_plane ! size(2,nely*nelz). (1, nely*nelz) is y coord and (2,nely*nelz) is z coord
!     close(10)
!     ! allocate(yz_map(1:nely*nelz))
!     yz_map(1:nely*nelz) = 0 ! Shows which glbal element no. corr. to each element in inflow vel. plane
!     ! Create map to find global element nos. that corr. to elements on vel. plane
!     ! Check on all cpus, and bcast to all, and add face element global number in an array ordered correctly as the vel. plane.
!     ! Later we can use findloc to get the index on vel slice corr. to a given global ele. no.
!     do yz_ele=1,nely*nelz ! Loop over elements in vel. plane
!      do l = 1,msh%labeled_zones(1)%size ! Loop over elements in Neko inflow plane 
!       el_and_facet = msh%labeled_zones(1)%facet_el(l)
!       facet = el_and_facet%x(1)
!       e = el_and_facet%x(2)
!       if (index_is_on_facet(1,1,1,lx,ly,lz, facet)) then
!!        if ((y(1,1,1,e).eq.yz_ele_plane(1,yz_ele)).and.(z(1,1,1,e).eq.yz_ele_plane(2,yz_ele))) then
!        if ((abs(y(1,1,1,e)-yz_ele_plane(1,yz_ele))<1e-8).and.(abs(z(1,1,1,e)-yz_ele_plane(2,yz_ele))<1e-8)) then
!            write(*,001) 'Plane:',yz_ele,yz_ele_plane(1,yz_ele),yz_ele_plane(2,yz_ele),&
!                & ' Neko:',e+msh%offset_el,y(1,1,1,e),z(1,1,1,e)
!            yz_map(yz_ele) = e+msh%offset_el ! Assigning the global element number in Neko inlet corresponding to this element in the velocity plane
!        end if      
!       end if
!      end do  
!     end do
!
!
!
!        write(*,*) 'Mask:', 'pe_rank=',pe_rank, yz_map(:)!, yz_ele_plane(1,1)
!        
!        do i = 1, 1!this%msk(0)
!           k = this%msk(i)
!           idx = nonlinear_index(k, lx, lx, lx)
!           loc = findloc(yz_map,idx(4)+msh%offset_el,dim=1) ! Finding index (loc) of the Neko element on the SIMSON inflow plane 
!           write(*,*) 'idx(4)=', idx(4), ', offset=', msh%offset_el, &
!&                ', idx(4)+msh%offset_el=',idx(4)+msh%offset_el, ', yz_map(loc)=', yz_map(loc), ', loc=',loc
!        end do
!        
!!        tcounter1 = tcounter1+1
!    end if

      if (tstep==1) then
!      write(*,*) 'Printing mask: gl_ele, i,x,y,z'    
      m = this%msk(0)
      allocate(map(1:m))
      min_y = 10000000.0_rp; min_z = 10000000.0_rp
      max_y = 0.0_rp; max_z = 0.0_rp
      ny1_simson = 6; nz1_simson = 6; ! GLL points within each element on the SIMSON mesh
      nely_simson = 10; nelz_simson = 10 ! No. of elements in y and z on SIMSON mesh
      ! Arrays to read in SIMSON mesh with the same structure as velocity planes
      allocate(y_simson(1:ny1_simson,1:nz1_simson,1:nely_simson*nelz_simson))
      allocate(z_simson(1:ny1_simson,1:nz1_simson,1:nely_simson*nelz_simson))
      y_simson(1:ny1_simson,1:nz1_simson,1:nely_simson*nelz_simson) = 0.0_rp
      z_simson(1:ny1_simson,1:nz1_simson,1:nely_simson*nelz_simson) = 0.0_rp
      open(unit=11,file='/scratch/ronith/Neko_1/neko/examples/test_cyl_boundary_layer/y_neko',&
              status='old',form='unformatted',access="stream")
      read(11) y_simson ! size(1:ny1_simson, 1:nz1_simson, 1:nely_simson*nelz_simson)
      close(11)
      open(unit=12,file='/scratch/ronith/Neko_1/neko/examples/test_cyl_boundary_layer/z_neko',& 
              status='old',form='unformatted',access="stream")
      read(12) z_simson ! size(1:ny1_simson, 1:nz1_simson, 1:nely_simson*nelz_simson)
      close(12)



      do i = 1, m
         k = this%msk(i)
         facet = this%facet(i)
         idx = nonlinear_index(k, lx, lx, lx)
!         write(*,*) idx(4)+this%msh%offset_el, i, xc(idx(1), idx(2), idx(3), idx(4)), &
!             yc(idx(1), idx(2), idx(3), idx(4)), &
!             zc(idx(1), idx(2), idx(3), idx(4))

!          if (i==1) write(*,*) 'First:', idx(4)+this%msh%offset_el, i, xc(idx(1), idx(2), idx(3), idx(4)), &
!             yc(idx(1), idx(2), idx(3), idx(4)), &
!             zc(idx(1), idx(2), idx(3), idx(4))
!          if (i==m) write(*,*) 'Last:', idx(4)+this%msh%offset_el, i, xc(idx(1), idx(2), idx(3), idx(4)), &
!             yc(idx(1), idx(2), idx(3), idx(4)), &
!             zc(idx(1), idx(2), idx(3), idx(4))

         select case(facet)
         case(1,2)
!            x(k) = 2.0_rp !this%bla(zc(idx(1), idx(2), idx(3), idx(4)), &
!!                 this%delta, this%x(1))
!            y(k) = 0.0_rp 
!            z(k) = 0.0_rp
!          loc_z(i)=zc(idx(1), idx(2), idx(3), idx(4))
          if (yc(idx(1), idx(2), idx(3), idx(4)).le.min_y) min_y=yc(idx(1), idx(2), idx(3), idx(4))
          if (yc(idx(1), idx(2), idx(3), idx(4)).ge.max_y) max_y=yc(idx(1), idx(2), idx(3), idx(4))
          if (zc(idx(1), idx(2), idx(3), idx(4)).le.min_z) min_z=zc(idx(1), idx(2), idx(3), idx(4))
          if (zc(idx(1), idx(2), idx(3), idx(4)).ge.max_z) max_z=zc(idx(1), idx(2), idx(3), idx(4))
!         case(3,4)
!            x(k) = 0.0_rp 
!            y(k) = 2.0_rp !this%bla(xc(idx(1), idx(2), idx(3), idx(4)), &
!!                 this%delta, this%x(2))
!            z(k) = 0.0_rp
!         case(5,6)
!            x(k) = 0.0_rp
!            y(k) = 0.0_rp
!            z(k) = 2.0_rp !this%bla(yc(idx(1), idx(2), idx(3), idx(4)), &
!!                 this%delta, this%x(3))
         end select            
      end do

!      write(*,*) 'Minval(loc_z)=', minval(loc_z),'Maxval(loc_z)=', maxval(loc_z), ', Maxloc(loc_z)=',maxloc(loc_z), ', m=',m
      write(*,*) 'min y=',min_y,', max y=',max_y,', min z=',min_z,', max z=',max_z

      do i = 1,1! m
         k = this%msk(i)
         facet = this%facet(i)
         idx = nonlinear_index(k, lx, lx, lx)

         do l = 1, 1!nely*nelz
            write(*,*) 'y_simson(:,:,l)',y_simson(:,:,l)
            write(*,*) 'y_simson(:,l,l)',y_simson(:,l,l)
            write(*,*) 'y_simson(l,:,l)',y_simson(l,:,l)
            write(*,*) 'y_simson(:,:,:)',y_simson(:,:,:)
            write(*,*) 'z_simson(:,:,:)',z_simson(:,:,:)
         end do ! do l = 1, nely*nelz 


       end do ! do i = 1, m  
    

!      first_gll=this%msk(1)
!      first_idx=nonlinear_index(first_gll, lx, lx, lx)
!      first_facet=this%facet(1)
!      last_gll=this%msk(this%msk(0))
!      last_idx=nonlinear_index(last_gll, lx, lx, lx)
!      last_facet=this%facet(this%msk(0))
!      select case(first_facet)
!          case(1,2)
!      write(*,*) 'Min y coord:',yc(first_idx(1), first_idx(2), first_idx(3), first_idx(4)), &
!                 'and min z coord:',zc(first_idx(1), first_idx(2), first_idx(3), first_idx(4))
!      end select
!      select case(last_facet)
!          case(1,2)
!      write(*,*) 'Max y coord:', yc(last_idx(1), last_idx(2), last_idx(3), last_idx(4)), &
!                  'and max z coord:', zc(last_idx(1), last_idx(2), last_idx(3), last_idx(4))
!      end select


!      write(*,*) 'Min and max y coord:',yc(first_idx(1), first_idx(2), first_idx(3), first_idx(4)), &
!                         yc(last_idx(1), last_idx(2), last_idx(3), last_idx(4))
!      write(*,*) 'Min and max z coord:',zc(first_idx(1), first_idx(2), first_idx(3), first_idx(4)), &
!                         zc(last_idx(1), last_idx(2), last_idx(3), last_idx(4))
     end if !if (tstep==1) then



      m = this%msk(0)
      do i = 1, m
         k = this%msk(i)
         facet = this%facet(i)
         idx = nonlinear_index(k, lx, lx, lx)
         select case(facet)
         case(1,2)
            x(k) = 2.0_rp !this%bla(zc(idx(1), idx(2), idx(3), idx(4)), &
!                 this%delta, this%x(1))
            y(k) = 0.0_rp 
            z(k) = 0.0_rp
         case(3,4)
            x(k) = 0.0_rp 
            y(k) = 2.0_rp !this%bla(xc(idx(1), idx(2), idx(3), idx(4)), &
!                 this%delta, this%x(2))
            z(k) = 0.0_rp
         case(5,6)
            x(k) = 0.0_rp
            y(k) = 0.0_rp
            z(k) = 2.0_rp !this%bla(yc(idx(1), idx(2), idx(3), idx(4)), &
!                 this%delta, this%x(3))
         end select            
      end do
    end associate
  end subroutine precursor_inflow_apply_vector

  !> Apply precursor inflow condition (vector valued) (device version)
  subroutine precursor_inflow_apply_vector_dev(this, x_d, y_d, z_d, t, tstep)
    class(precursor_inflow_t), intent(inout), target :: this
    type(c_ptr) :: x_d
    type(c_ptr) :: y_d
    type(c_ptr) :: z_d
    real(kind=rp), intent(in), optional :: t
    integer, intent(in), optional :: tstep
    integer :: i, m, k, idx(4), facet
    integer(c_size_t) :: s
    real(kind=rp), allocatable :: bla_x(:), bla_y(:), bla_z(:)

    associate(xc => this%c%dof%x, yc => this%c%dof%y, zc => this%c%dof%z, &
         nx => this%c%nx, ny => this%c%ny, nz => this%c%nz, &
         lx => this%c%Xh%lx, blax_d => this%blax_d, blay_d => this%blay_d, &
         blaz_d => this%blaz_d)

      m = this%msk(0)


      ! Pretabulate values during first call to apply
      if (.not. c_associated(blax_d)) then
         allocate(bla_x(m), bla_y(m), bla_z(m)) ! Temp arrays

         if (rp .eq. REAL32) then
            s = m * 4
         else if (rp .eq. REAL64) then
            s = m * 8
         end if

         call device_alloc(blax_d, s)
         call device_alloc(blay_d, s)
         call device_alloc(blaz_d, s)
         
         do i = 1, m
            k = this%msk(i)
            facet = this%facet(i)
            idx = nonlinear_index(k, lx, lx, lx)
            select case(facet)
            case(1,2)
               bla_x(i) = 2.0_rp !this%bla(zc(idx(1), idx(2), idx(3), idx(4)), &
!                    this%delta, this%x(1))
               bla_y(i) = 0.0_rp 
               bla_z(i) = 0.0_rp
            case(3,4)
               bla_x(i) = 0.0_rp 
               bla_y(i) = 2.0-rp !this%bla(xc(idx(1), idx(2), idx(3), idx(4)), &
!                    this%delta, this%x(2))
               bla_z(i) = 0.0_rp
            case(5,6)
               bla_x(i) = 0.0_rp
               bla_y(i) = 0.0_rp
               bla_z(i) = 2.0_rp !this%bla(yc(idx(1), idx(2), idx(3), idx(4)), &
!                    this%delta, this%x(3))
            end select
         end do

         call device_memcpy(bla_x, blax_d, m, HOST_TO_DEVICE)
         call device_memcpy(bla_y, blay_d, m, HOST_TO_DEVICE)
         call device_memcpy(bla_z, blaz_d, m, HOST_TO_DEVICE)

         deallocate(bla_x, bla_y, bla_z)
      end if

      call device_inhom_dirichlet_apply_vector(this%msk_d, x_d, y_d, z_d, &
           blax_d, blay_d, blaz_d, m)
      
    end associate

  end subroutine precursor_inflow_apply_vector_dev

!  !> Set Blasius parameters
!  subroutine blasius_set_params(this, delta, type)
!    class(blasius_t), intent(inout) :: this
!    real(kind=rp) :: delta
!    character(len=*) :: type
!    this%delta = delta
!    
!    select case(trim(type))
!    case('linear')
!       this%bla => blasius_linear
!    case('quadratic')
!       this%bla => blasius_quadratic
!    case('cubic')
!       this%bla => blasius_cubic
!    case('quartic')
!       this%bla => blasius_quartic
!    case('sin')
!       this%bla => blasius_sin
!    case default
!       call neko_error('Invalid Blasius approximation')
!    end select
!  end subroutine blasius_set_params

  !> Assign coefficients (facet normals etc)
  subroutine precursor_inflow_set_coef(this, c)
    class(precursor_inflow_t), intent(inout) :: this
    type(coef_t), target, intent(inout) :: c
    this%c => c
  end subroutine precursor_inflow_set_coef
     
end module precursor_inflow
