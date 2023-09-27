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
  use comm, only: pe_rank
!  use flow_profile
  use, intrinsic :: iso_c_binding
  implicit none
  private

  !> Blasius profile for inlet (vector valued)
  type, public, extends(inflow_t) :: precursor_inflow_t
     type(coef_t), pointer :: c => null()
     !> Map between GLL points on SIMSON inflow plane to Neko inflow plane
     integer, allocatable :: yz_map(:) 
!     real(kind=rp) :: delta
!     procedure(blasius_profile), nopass, pointer :: bla => null()
     type(c_ptr), private :: blax_d = C_NULL_PTR
     type(c_ptr), private :: blay_d = C_NULL_PTR
     type(c_ptr), private :: blaz_d = C_NULL_PTR
   contains
     !> Constructor  
     procedure, pass(this) :: init => precursor_inflow_init  
     !> Destructor
    !  procedure, pass(this) :: free => precursor_inflow_free  
     procedure, pass(this) :: apply_scalar => precursor_inflow_apply_scalar
     procedure, pass(this) :: apply_vector => precursor_inflow_apply_vector
     procedure, pass(this) :: apply_scalar_dev => precursor_inflow_apply_scalar_dev
     procedure, pass(this) :: apply_vector_dev => precursor_inflow_apply_vector_dev
!     procedure, pass(this) :: set_params => blasius_set_params
     procedure, pass(this) :: set_coef => precursor_inflow_set_coef
  end type precursor_inflow_t

contains

  subroutine precursor_inflow_init(this, x, y, z, n)
!  subroutine precursor_inflow_init(this)
   type(precursor_inflow_t), intent(inout) :: this
    integer, intent(in) :: n
    real(kind=rp), intent(inout),  dimension(n) :: x
    real(kind=rp), intent(inout),  dimension(n) :: y
    real(kind=rp), intent(inout),  dimension(n) :: z
    !real(kind=rp), intent(in), optional :: t
    !integer, intent(in), optional :: tstep
    integer :: i, m, k, idx(4), facet
    !integer, allocatable :: yz_map(:)
    !real(kind=rp), allocatable :: u_in_glob1(:), u_in_glob2(:)
    !real(kind=rp), allocatable :: u_in_loc(:,:)
    integer :: nely, nelz
    integer :: first_gll, last_gll, first_idx(4), last_idx(4)
    integer :: first_facet, last_facet
    real(kind=rp), allocatable :: y_simson(:,:,:), z_simson(:,:,:)
    real(kind=rp) :: y_min, y_max, z_min, z_max  
    real(kind=rp) :: y_coord, z_coord
    integer :: ny1_simson, nz1_simson, nely_simson, nelz_simson, nel_yz_simson
    integer :: el_simson
    real(kind=rp) :: y_min_el_simson, y_max_el_simson, z_min_el_simson, z_max_el_simson
    real(kind=rp), allocatable :: flat_y(:), flat_z(:)
    integer :: kk, j, map

    associate(xc => this%c%dof%x, yc => this%c%dof%y, zc => this%c%dof%z, &
         nx => this%c%nx, ny => this%c%ny, nz => this%c%nz, &
         lx => this%c%Xh%lx)

! Setting up the map at tstep==1
   !   if (tstep==1) then
!      write(*,*) 'Printing mask: gl_ele, i,x,y,z'    
      m = this%msk(0)
      allocate(this%yz_map(1:m))
      y_min = 10000000.0_rp; z_min = 10000000.0_rp
      y_max = 0.0_rp; z_max = 0.0_rp
      y_min_el_simson = 0.0_rp
      y_max_el_simson = 0.0_rp
      z_min_el_simson = 0.0_rp
      z_max_el_simson = 0.0_rp
! Later I can take these input constants through the json file      
      ny1_simson = 8; nz1_simson = 8; ! GLL points within each element on the SIMSON mesh
      nely_simson = 10; nelz_simson = 10 ! No. of elements in y and z on SIMSON mesh
      nel_yz_simson = nely_simson*nelz_simson

      ! Arrays to read in SIMSON mesh with the same structure as velocity planes
      allocate(y_simson(1:ny1_simson,1:nz1_simson,1:nely_simson*nelz_simson))
      allocate(z_simson(1:ny1_simson,1:nz1_simson,1:nely_simson*nelz_simson))
      y_simson(1:ny1_simson,1:nz1_simson,1:nely_simson*nelz_simson) = 0.0_rp
      z_simson(1:ny1_simson,1:nz1_simson,1:nely_simson*nelz_simson) = 0.0_rp
      ! Global/Whole inflow velocity plane
      !allocate(u_in_glob1(1:ny1_simson*nz1_simson*nel_yz_simson))
      !allocate(u_in_glob2(1:ny1_simson*nz1_simson*nel_yz_simson))
      ! Local/rank-specific inflow velocity plane
      !allocate(u_in_loc(1:2,1:m))
      !u_in_glob1(:) = 0.0_rp
      !u_in_glob2(:) = 0.0_rp
      !u_in_loc(1:2,1:m) = 0.0_rp

      open(unit=11,file='/scratch/ronith/Neko_1/neko/examples/test_cyl_boundary_layer/y_neko',&
              status='old',form='unformatted',access="stream")
      read(11) y_simson ! size(1:ny1_simson, 1:nz1_simson, 1:nely_simson*nelz_simson)
      close(11)
      open(unit=12,file='/scratch/ronith/Neko_1/neko/examples/test_cyl_boundary_layer/z_neko',&
              status='old',form='unformatted',access="stream")
      read(12) z_simson ! size(1:ny1_simson, 1:nz1_simson, 1:nely_simson*nelz_simson)
      close(12)
     

      ! Finding y,z min/max coords. of gll points lying on the inlet plane within each rank 
      ! that contains part of inlet plane
      do i = 1, m
         k = this%msk(i)
         facet = this%facet(i)
         idx = nonlinear_index(k, lx, lx, lx)
         y_coord = yc(idx(1), idx(2), idx(3), idx(4))
         z_coord = zc(idx(1), idx(2), idx(3), idx(4))
!         select case(facet)
!         case(1,2)
          if (y_coord.le.y_min) y_min = y_coord
          if (y_coord.ge.y_max) y_max = y_coord
          if (z_coord.le.z_min) z_min = z_coord
          if (z_coord.ge.z_max) z_max = z_coord
!         end select            
      end do

!      write(*,*) 'Minval(loc_z)=', minval(loc_z),'Maxval(loc_z)=', maxval(loc_z), ', Maxloc(loc_z)=',maxloc(loc_z), ', m=',m
      write(*,*) 'Neko: pe_rank=',pe_rank,', min y=',y_min,', max y=',y_max,', min z=',z_min,', max z=',z_max
      map=0
      ! Map each gll point on the local inlet plane to a corresponding point on the SIMSON plane
      do i = 1, m
         k = this%msk(i)
         !facet = this%facet(i)
         idx = nonlinear_index(k, lx, lx, lx)
         y_coord = yc(idx(1), idx(2), idx(3), idx(4))
         z_coord = zc(idx(1), idx(2), idx(3), idx(4))

         do el_simson = 1, nel_yz_simson
            !write(*,*) 'Max y-coord in ele=100: maxval(y_simson(:,8,100))', maxval(y_simson(:,nz1_simson,100)) ! Use this
            !write(*,*) 'Min y-coord in ele=100: minval(y_simson(:,8,100))', minval(y_simson(:,nz1_simson,100)) ! Use this
            !write(*,*) 'Max z-coord in ele=100: maxval(z_simson(8,:,100))', maxval(z_simson(ny1_simson,:,100)) ! Use this
            !write(*,*) 'Min z-coord in ele=100: minval(z_simson(8,:,100))', minval(z_simson(ny1_simson,:,100)) ! Use this
            
            ! Finding min/max within each element on SIMSON plane to determine which elements can be mapped to the current rank
            y_min_el_simson = minval(y_simson(:,nz1_simson,el_simson))
            y_max_el_simson = maxval(y_simson(:,nz1_simson,el_simson))
            z_min_el_simson = minval(z_simson(ny1_simson,:,el_simson))
            z_max_el_simson = maxval(z_simson(ny1_simson,:,el_simson))
  
            ! If this SIMSON element lies within the current rank
            if (((y_min_el_simson.ge.y_min).and.(y_max_el_simson.le.y_max)) &
                                           .and. & 
                ((z_min_el_simson.ge.z_min).and.(z_max_el_simson.le.z_max))) then
            do kk = 1, nz1_simson
               do j = 1, ny1_simson
                  if ( (abs(y_simson(j,kk,el_simson)-y_coord).le.1e-07_rp).and. &
                      (abs(z_simson(j,kk,el_simson)-z_coord).le.1e-07_rp)  ) then
                  !if ( ((y_simson(j,k,el_simson)==y_coord)).and.(z_simson(j,k,el_simson)==z_coord)  ) then
                     ! Save the flattened 1D index of the 3D array y_simson(1:ny1_simson, 1:nz1_simson,1_nel_yz_simson) as map
                     ! The velocity files are saved in binary as a vector of the same form from the matlab file
                     this%yz_map(i) = el_simson + nel_yz_simson * (kk-1) + (nel_yz_simson*nz1_simson) * (j-1)
!                     write(*,*) 'Test inside:', y_coord, y_simson(j,k,el_simson), yz_map(i)
                     map = map+1
                  end if  
               end do ! j =1, ny1_simson 
            end do ! kk = 1, nz1_simson
            end if
          end do ! do el_simson = 1, nel_yz_simson
       end do ! do i = 1, m 
 
       write(*,*) 'Total mapped points=',map

    end associate
  end subroutine precursor_inflow_init


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
    integer :: i, m, k!, idx(4), facet
!    integer, allocatable :: yz_map(:)
    real(kind=rp), allocatable :: u_in_glob1(:), u_in_glob2(:)
    real(kind=rp), allocatable :: u_in_loc(:,:)
    integer :: nely, nelz
!    integer :: first_gll, last_gll, first_idx(4), last_idx(4)
!    integer :: first_facet, last_facet
!    real(kind=rp), allocatable :: y_simson(:,:,:), z_simson(:,:,:)
!    real(kind=rp) :: y_min, y_max, z_min, z_max  
!    real(kind=rp) :: y_coord, z_coord
!    integer :: ny1_simson, nz1_simson, nely_simson, nelz_simson, nel_yz_simson
!    integer :: el_simson
!    real(kind=rp) :: y_min_el_simson, y_max_el_simson, z_min_el_simson, z_max_el_simson
!    real(kind=rp), allocatable :: flat_y(:), flat_z(:)
    integer :: kk, j, map

    associate(xc => this%c%dof%x, yc => this%c%dof%y, zc => this%c%dof%z, &
         nx => this%c%nx, ny => this%c%ny, nz => this%c%nz, &
         lx => this%c%Xh%lx)

! Setting up the map at tstep==1
   !   if (tstep==1) then
!      write(*,*) 'Printing mask: gl_ele, i,x,y,z'    
      m = this%msk(0)
!      allocate(yz_map(1:m))
!      y_min = 10000000.0_rp; z_min = 10000000.0_rp
!      y_max = 0.0_rp; z_max = 0.0_rp
!      y_min_el_simson = 0.0_rp
!      y_max_el_simson = 0.0_rp
!      z_min_el_simson = 0.0_rp
!      z_max_el_simson = 0.0_rp
!! Later I can take these input constants through the json file      
!      ny1_simson = 8; nz1_simson = 8; ! GLL points within each element on the SIMSON mesh
!      nely_simson = 10; nelz_simson = 10 ! No. of elements in y and z on SIMSON mesh
!      nel_yz_simson = nely_simson*nelz_simson
!
!      ! Arrays to read in SIMSON mesh with the same structure as velocity planes
!      allocate(y_simson(1:ny1_simson,1:nz1_simson,1:nely_simson*nelz_simson))
!      allocate(z_simson(1:ny1_simson,1:nz1_simson,1:nely_simson*nelz_simson))
!      y_simson(1:ny1_simson,1:nz1_simson,1:nely_simson*nelz_simson) = 0.0_rp
!      z_simson(1:ny1_simson,1:nz1_simson,1:nely_simson*nelz_simson) = 0.0_rp
      ! Global/Whole inflow velocity plane
      allocate(u_in_glob1(1:ny1_simson*nz1_simson*nel_yz_simson))
      allocate(u_in_glob2(1:ny1_simson*nz1_simson*nel_yz_simson))
!      ! Local/rank-specific inflow velocity plane
!      allocate(u_in_loc(1:2,1:m))
      u_in_glob1(:) = 0.0_rp
      u_in_glob2(:) = 0.0_rp
!      u_in_loc(1:2,1:m) = 0.0_rp
!
!      open(unit=11,file='/scratch/ronith/Neko_1/neko/examples/test_cyl_boundary_layer/y_neko',&
!              status='old',form='unformatted',access="stream")
!      read(11) y_simson ! size(1:ny1_simson, 1:nz1_simson, 1:nely_simson*nelz_simson)
!      close(11)
!      open(unit=12,file='/scratch/ronith/Neko_1/neko/examples/test_cyl_boundary_layer/z_neko',&
!              status='old',form='unformatted',access="stream")
!      read(12) z_simson ! size(1:ny1_simson, 1:nz1_simson, 1:nely_simson*nelz_simson)
!      close(12)
!     
!
!      ! Finding y,z min/max coords. of gll points lying on the inlet plane within each rank 
!      ! that contains part of inlet plane
!      do i = 1, m
!         k = this%msk(i)
!         facet = this%facet(i)
!         idx = nonlinear_index(k, lx, lx, lx)
!         y_coord = yc(idx(1), idx(2), idx(3), idx(4))
!         z_coord = zc(idx(1), idx(2), idx(3), idx(4))
!!         select case(facet)
!!         case(1,2)
!          if (y_coord.le.y_min) y_min = y_coord
!          if (y_coord.ge.y_max) y_max = y_coord
!          if (z_coord.le.z_min) z_min = z_coord
!          if (z_coord.ge.z_max) z_max = z_coord
!!         end select            
!      end do
!
!!      write(*,*) 'Minval(loc_z)=', minval(loc_z),'Maxval(loc_z)=', maxval(loc_z), ', Maxloc(loc_z)=',maxloc(loc_z), ', m=',m
!      write(*,*) 'Neko: pe_rank=',pe_rank,', min y=',y_min,', max y=',y_max,', min z=',z_min,', max z=',z_max
!      map=0
!      ! Map each gll point on the local inlet plane to a corresponding point on the SIMSON plane
!      do i = 1, m
!         k = this%msk(i)
!         !facet = this%facet(i)
!         idx = nonlinear_index(k, lx, lx, lx)
!         y_coord = yc(idx(1), idx(2), idx(3), idx(4))
!         z_coord = zc(idx(1), idx(2), idx(3), idx(4))
!
!         do el_simson = 1, nel_yz_simson
!            !write(*,*) 'Max y-coord in ele=100: maxval(y_simson(:,8,100))', maxval(y_simson(:,nz1_simson,100)) ! Use this
!            !write(*,*) 'Min y-coord in ele=100: minval(y_simson(:,8,100))', minval(y_simson(:,nz1_simson,100)) ! Use this
!            !write(*,*) 'Max z-coord in ele=100: maxval(z_simson(8,:,100))', maxval(z_simson(ny1_simson,:,100)) ! Use this
!            !write(*,*) 'Min z-coord in ele=100: minval(z_simson(8,:,100))', minval(z_simson(ny1_simson,:,100)) ! Use this
!            
!            ! Finding min/max within each element on SIMSON plane to determine which elements can be mapped to the current rank
!            y_min_el_simson = minval(y_simson(:,nz1_simson,el_simson))
!            y_max_el_simson = maxval(y_simson(:,nz1_simson,el_simson))
!            z_min_el_simson = minval(z_simson(ny1_simson,:,el_simson))
!            z_max_el_simson = maxval(z_simson(ny1_simson,:,el_simson))
!  
!            ! If this SIMSON element lies within the current rank
!            if (((y_min_el_simson.ge.y_min).and.(y_max_el_simson.le.y_max)) &
!                                           .and. & 
!                ((z_min_el_simson.ge.z_min).and.(z_max_el_simson.le.z_max))) then
!            do kk = 1, nz1_simson
!               do j = 1, ny1_simson
!                  if ( (abs(y_simson(j,kk,el_simson)-y_coord).le.1e-07_rp).and. &
!                      (abs(z_simson(j,kk,el_simson)-z_coord).le.1e-07_rp)  ) then
!                  !if ( ((y_simson(j,k,el_simson)==y_coord)).and.(z_simson(j,k,el_simson)==z_coord)  ) then
!                     ! Save the flattened 1D index of the 3D array y_simson(1:ny1_simson, 1:nz1_simson,1_nel_yz_simson) as map
!                     ! The velocity files are saved in binary as a vector of the same form from the matlab file
!                     yz_map(i) = el_simson + nel_yz_simson * (kk-1) + (nel_yz_simson*nz1_simson) * (j-1)
!!                     write(*,*) 'Test inside:', y_coord, y_simson(j,k,el_simson), yz_map(i)
!                     map = map+1
!                  end if  
!               end do ! j =1, ny1_simson 
!            end do ! kk = 1, nz1_simson
!            end if
!          end do ! do el_simson = 1, nel_yz_simson
!       end do ! do i = 1, m 
 
       write(*,*) 'Total mapped points=',map

       open(unit=13,file='/scratch/ronith/Neko_1/neko/examples/test_cyl_boundary_layer/buin1D_0',&
              status='old',form='unformatted',access="stream")
       read(13) u_in_glob1(:) ! 3D array flattened into 1D vec. of size(1:ny1_simson*nz1_simson*nely_simson*nelz_simson)
       close(13)

       !write(*,*) u_in_glob1(:)

       open(unit=14,file='/scratch/ronith/Neko_1/neko/examples/test_cyl_boundary_layer/buin1D_1',&
              status='old',form='unformatted',access="stream")
       read(14) u_in_glob2(:) ! 3D array flattened into 1D vec. of size(1:ny1_simson*nz1_simson*nely_simson*nelz_simson)
       close(14)

       !write(*,*) u_in_glob2(:)
       !write(*,*) 'u_in_loc(54,1000)' ,u_in_loc(54,1000)
       !write(*,*) 'yz_map(54)',yz_map(54)
       !write(*,*) 'u_in_global(yz_map(54))',u_in_glob1(yz_map(54))
       !write(*,*) 'u_in_loc(1,54)',u_in_loc(1,54)
       !write(*,*) 'shape(u_in_glob1)', shape(u_in_glob1)
       !write(*,*) 'shpae(u_in_loc)',shape(u_in_loc)
       !write(*,*) 'shape(yz_map)',shape(yz_map)

       ! Extract mapped local velocities into a new vector which will be used for temp. interpolation as well
       do i = 1, m
         !do el_simson = 1, nel_yz_simson
         !   do k = 1, nz1_simson
         !      do j = 1, ny1_simson
                  ! Use the flattened 1D index of the 3D array y_simson(1:ny1_simson, 1:nz1_simson,1_nel_yz_simson) as map
        !write(*,*) 'i=',i, '/',m,', ' ,yz_map(i),'/',ny1_simson*nz1_simson*nel_yz_simson
         u_in_loc(1,i) = u_in_glob1(yz_map(i))
         u_in_loc(2,i) = u_in_glob2(yz_map(i))
         !      end do ! j =1, ny1_simson
         !   end do ! k =1, nz1_simson   
         !end do ! el_simson = 1, nel_yz_simson
        end do ! i =1, m


     !allocate(flat_y(1:m))
     !allocate(flat_z(1:m))
     !!flat_y = pack(y_simson,.true.) 

     !do el_simson = 1, nel_yz_simson
     !   do k = 1, nz1_simson
     !      do j = 1, ny1_simson
     !         flat_y( el_simson + nel_yz_simson * (k-1) + (nel_yz_simson*nz1_simson) * (j-1)) = y_simson(j,k,el_simson)
     !         flat_z( el_simson + nel_yz_simson * (k-1) + (nel_yz_simson*nz1_simson) * (j-1)) = z_simson(j,k,el_simson)
     !      end do
     !   end do
     !end do

     !do i=1, 100!m
     ! k = this%msk(i)
     ! idx = nonlinear_index(k, lx, lx, lx)
     ! y_coord = yc(idx(1), idx(2), idx(3), idx(4))
     ! z_coord = zc(idx(1), idx(2), idx(3), idx(4))
     ! write(*,*) 'Test y:',y_coord, flat_y(yz_map(i)), yz_map(i) 
     ! write(*,*) 'Test z:',z_coord, flat_z(yz_map(i)), yz_map(i) 
     !end do


  !   end if !if (tstep==1) then

      write(*,*) 'size of yz_map:',size(yz_map)
      write(*,*) 'shape of yz_map:',shape(yz_map)
      write(*,*) 'size of u_in_loc:',size(u_in_loc)
      write(*,*) 'shape of u_in_loc:',shape(u_in_loc)

      m = this%msk(0)
      do i = 1, m
         k = this%msk(i)
         facet = this%facet(i)
         idx = nonlinear_index(k, lx, lx, lx)
         select case(facet)
         case(1,2)
            write(*,*) 'tstep=',tstep,', i=',i,', k=',k
            x(k) = u_in_loc(1,i) !2.0_rp !this%bla(zc(idx(1), idx(2), idx(3), idx(4)), &
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
               bla_y(i) = 2.0_rp !this%bla(xc(idx(1), idx(2), idx(3), idx(4)), &
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
