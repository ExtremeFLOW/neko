!> Maps a 3D dofmap to a 2D spectral element grid.

module map_2d
  use num_types, only: rp
  use space, only: space_t
  use dofmap, only: dofmap_t
  use map_1d
  use gather_scatter
  use mesh, only: mesh_t
  use device
  use utils
  use comm
  use field_list
  use coefs, only: coef_t
  use field_list, only: field_list_t
  use matrix, only: matrix_t
  use vector, only: vector_ptr_t
  use logger, only: neko_log, LOG_SIZE
  use utils, only: neko_error, neko_warning
  use math, only: glmax, glmin, glimax, glimin, relcmp, cmult, &
                  add2s1, col2, copy, rzero
  use neko_mpi_types
  use fld_file_data
  use field
  use, intrinsic :: iso_c_binding
  implicit none
  private
  
  type, public :: map_2d_t
     integer :: nelv_2d = 0 !< Number of elements in 2D mesh on this rank
     integer :: glb_nelv_2d = 0 !< global number of elements in 2d
     integer :: offset_el_2d = 0 !< element offset for this rank
     integer :: lxy = 0 !< number of gll points per 2D element
     integer :: n_2d = 0 !< total number of gll points (nelv_2d*lxy)
     integer, allocatable :: idx_2d(:) !< Mapping of GLL point from 3D to 2D
     integer, allocatable :: el_idx_2d(:) !< Mapping of element in 3D to 2D
     type(map_1d_t) :: map_1d !< 1D map in normal direction to 2D plane
     type(mesh_t), pointer :: msh !< 3D mesh
     type(dofmap_t), pointer :: dof => null() !< 3D dofmap
     type(coef_t), pointer :: coef => null() !< 3D SEM coefs
     type(field_t) :: u !< Work array 1, naming based on when they are used
     type(field_t) :: old_u !< Work array 2
     type(field_t) :: avg_u !< Work array 3
     type(field_t) :: el_heights !< Weight elements by their size in integral
     integer :: dir !< direction normal to 2D plane
     real(kind=rp) :: domain_height !< total height of 3D domain
   contains
     procedure, pass(this) :: init_int => map_2d_init
     procedure, pass(this) :: init_char => map_2d_init_char
     generic :: init => init_int, init_char
     procedure, pass(this) :: average_file => map_2d_average
     procedure, pass(this) :: average_list => map_2d_average_field_list
     generic :: average => average_list, average_file
  end type map_2d_t

contains

  subroutine map_2d_init(this,coef, dir, tol)
    class(map_2d_t), intent(inout) :: this
    type(coef_t), intent(inout), target :: coef
    integer, intent(in) :: dir
    real(kind=rp), intent(in) :: tol
    real(kind=rp) :: el_dim(3,3), el_h
    integer :: i, e, j, ierr, k, lx, lxy, n
    call this%map_1d%init(coef,dir,tol)
    this%msh => coef%msh
    this%coef => coef
    this%dof => coef%dof
    
    call this%u%init(this%dof)
    call this%old_u%init(this%dof)
    call this%avg_u%init(this%dof)
    call this%el_heights%init(this%dof)
    this%dir = dir

    n = this%dof%size()
    lx = this%dof%Xh%lx
    lxy = this%dof%Xh%lxy
    this%nelv_2d = 0
    do i = 1, this%msh%nelv
       if (this%map_1d%el_lvl(i) .eq. 1) then
          this%nelv_2d = this%nelv_2d + 1
       end if
    end do
    this%glb_nelv_2d = 0
    call MPI_Allreduce(this%nelv_2d, this%glb_nelv_2d, 1, &
         MPI_INTEGER, MPI_SUM, NEKO_COMM, ierr)

    this%offset_el_2d = 0
    call MPI_Exscan(this%nelv_2d, this%offset_el_2d, 1, &
         MPI_INTEGER, MPI_SUM, NEKO_COMM, ierr)
    allocate(this%el_idx_2d(this%nelv_2d)) 
    do i = 1, this%nelv_2d
       this%el_idx_2d(i) = this%offset_el_2d + i
    end do
    this%n_2d = this%nelv_2d*lxy
    allocate(this%idx_2d(this%n_2d))

    j = 0
    do e = 1, this%msh%nelv
       if (this%map_1d%el_lvl(e) .eq. 1) then
          if (this%map_1d%dir_el(e) .eq. 1) then
             do i = 1, lxy
                this%idx_2d(j*lxy+i) = linear_index(1,i,1,e,lx,lx,lx)
             end do
          end if
          if (this%map_1d%dir_el(e) .eq. 2) then
             do i = 1, lx
                do k = 1, lx
                   this%idx_2d(j*lxy+k+lx*(i-1)) = linear_index(k,1,i,e,lx,lx,lx)
                end do
             end do
          end if
          if (this%map_1d%dir_el(e) .eq. 3) then
             do i = 1, lxy
                this%idx_2d(j*lxy+i) = linear_index(i,1,1,e,lx,lx,lx)
             end do
          end if
          j = j +1
       end if
    end do
    do i = 1, this%msh%nelv
       !find height in hom-dir
       !direction in local coords (r,s,t) that is hom is stored in map_1d%dir_el
       !set element to height
       !we assume elements are stacked on eachother...
       el_dim(1,:) = abs(this%msh%elements(i)%e%pts(1)%p%x-&
                     this%msh%elements(i)%e%pts(2)%p%x)
       el_dim(2,:) = abs(this%msh%elements(i)%e%pts(1)%p%x-&
                     this%msh%elements(i)%e%pts(3)%p%x)
       el_dim(3,:) = abs(this%msh%elements(i)%e%pts(1)%p%x-&
                     this%msh%elements(i)%e%pts(5)%p%x)
       ! 1 corresponds to x, 2 to y, 3 to z
       el_h = el_dim(this%map_1d%dir_el(i),dir)
       this%el_heights%x(:,:,:,i) = el_h
    end do
    !Need to compute mapping for 3d to 2d
    !Does order matter? Think its ok as long as values written in same order

    call copy(this%u%x,this%el_heights%x,n)
    call copy(this%old_u%x,this%el_heights%x,n)
    call copy(this%avg_u%x,this%el_heights%x,n)
    call perform_global_summation(this%u, this%avg_u, this%old_u, &
         this%map_1d%n_el_lvls, &
         this%map_1d%dir_el,this%coef%gs_h, this%coef%mult, this%msh%nelv, lx)
    this%domain_height = this%u%x(1,1,1,1)

  end subroutine map_2d_init
   
  subroutine map_2d_init_char(this, coef,  dir, tol)
    class(map_2d_t) :: this
    type(coef_t), intent(inout), target :: coef
    character(len=*), intent(in) :: dir
    real(kind=rp), intent(in) :: tol
    integer :: idir

    if (trim(dir) .eq. 'x') then
       idir = 1
    else if (trim(dir) .eq. 'y') then
       idir = 2
    else if (trim(dir) .eq. 'z') then
       idir = 3
    else
       call neko_error('Direction not supported, map_2d')
    end if

    call this%init(coef,idir,tol)

  end subroutine map_2d_init_char


  !> Computes average if field list in one direction and outputs 2D field
  !! with averaged values. 
  !! @param fld_data2D output 2D averages
  !! @param field_list list of fields to be averaged
  subroutine map_2d_average_field_list(this,fld_data2D,fld_data3D)
    class(map_2d_t), intent(inout) :: this
    type(fld_file_data_t), intent(inout) :: fld_data2D
    type(field_list_t), intent(inout) :: fld_data3D
    real(kind=rp), pointer, dimension(:,:,:,:) :: x_ptr, y_ptr
    
    type(vector_ptr_t), allocatable :: fields2d(:)
    integer :: n_2d, n
    integer :: i, j, lx, lxy, e

    call fld_data2D%init(this%nelv_2d,this%offset_el_2d)
    fld_data2D%gdim = 2
    fld_data2D%lx = this%dof%Xh%lx
    fld_data2D%ly = this%dof%Xh%ly
    fld_data2D%lz = 1
    lx = this%dof%Xh%lx
    fld_data2D%glb_nelv = this%glb_nelv_2d
    lxy = fld_data2D%lx*fld_data2D%ly
    n_2d = lxy*this%nelv_2d
    n = this%dof%size()
    call fld_data2D%x%init(n_2d)
    call fld_data2D%y%init(n_2d)
    allocate(fld_data2D%idx(this%nelv_2d))

    if (this%dir .eq. 1) then
       x_ptr => this%dof%z
       y_ptr => this%dof%y
    end if
    if (this%dir .eq. 2) then
        x_ptr => this%dof%x
        y_ptr => this%dof%z
    end if
    if (this%dir .eq. 3) then
        x_ptr => this%dof%x
        y_ptr => this%dof%y
    end if
    do j = 1, this%nelv_2d
       fld_data2d%idx(j) = this%el_idx_2d(j)
    end do
    do j = 1, n_2d
       fld_data2d%x%x(j) = x_ptr(this%idx_2d(j),1,1,1)
       fld_data2d%y%x(j) = y_ptr(this%idx_2d(j),1,1,1)
    end do
    allocate(fields2D(fld_data3D%size()))

    call fld_data2D%init_n_fields(fld_data3D%size(),n_2d)
    this%u = 0.0_rp
    this%old_u = 0.0_rp
    this%avg_u = 0.0_rp
    do i = 1, fld_data3D%size()
       call copy(this%old_u%x,fld_data3D%items(i)%ptr%x,n)
       call perform_local_summation(this%u,this%old_u, &
            this%el_heights, this%domain_height, &
            this%map_1d%dir_el, this%coef, this%msh%nelv, lx)
       call copy(this%old_u%x,this%u%x,n)
       call copy(this%avg_u%x,this%u%x,n)
       call perform_global_summation(this%u, this%avg_u, & 
            this%old_u, this%map_1d%n_el_lvls, &
            this%map_1d%dir_el,this%coef%gs_h, this%coef%mult, &
            this%msh%nelv, lx)
       call copy(fld_data3D%items(i)%ptr%x,this%u%x,n)
    end do
    call fld_data2D%get_list(fields2d,fld_data2D%size())
    do i = 1, fld_data3D%size()
       do j = 1, n_2d
          fields2d(i)%ptr%x(j) = fld_data3D%items(i)%ptr%x(this%idx_2d(j),1,1,1)
       end do
    end do

  end subroutine map_2d_average_field_list
 
  
  !> Computes average if field list in one direction and outputs 2D field
  !! with averaged values. 
  !! @param fld_data2D output 2D averages
  !! @param fld_data3D fld_file_data of fields to be averaged
  subroutine map_2d_average(this,fld_data2D,fld_data3D)
    class(map_2d_t), intent(inout) :: this
    type(fld_file_data_t), intent(inout) :: fld_data2D
    type(fld_file_data_t), intent(inout) :: fld_data3D
    real(kind=rp), pointer, dimension(:,:,:,:) :: x_ptr, y_ptr
    
    type(vector_ptr_t), allocatable :: fields3d(:), fields2d(:)
    integer :: n_2d, n
    integer :: i, j, lx, lxy, e

    call fld_data2D%init(this%nelv_2d,this%offset_el_2d)
    fld_data2D%gdim = 2
    fld_data2D%lx = this%dof%Xh%lx
    fld_data2D%ly = this%dof%Xh%ly
    fld_data2D%lz = 1
    lx = this%dof%Xh%lx
    fld_data2D%glb_nelv = this%glb_nelv_2d
    lxy = fld_data2D%lx*fld_data2D%ly
    n_2d = lxy*this%nelv_2d
    n = this%dof%size()
    call fld_data2D%x%init(n_2d)
    call fld_data2D%y%init(n_2d)
    allocate(fld_data2D%idx(n_2d))

    if (this%dir .eq. 1) then
       x_ptr => this%dof%z
       y_ptr => this%dof%y
    end if
    if (this%dir .eq. 2) then
        x_ptr => this%dof%x
        y_ptr => this%dof%z
    end if
    if (this%dir .eq. 3) then
        x_ptr => this%dof%x
        y_ptr => this%dof%y
    end if
    do j = 1, n_2d
       fld_data2d%idx(j) = this%idx_2d(j)
       fld_data2d%x%x(j) = x_ptr(this%idx_2d(j),1,1,1)
       fld_data2d%y%x(j) = y_ptr(this%idx_2d(j),1,1,1)
    end do
    allocate(fields3D(fld_data3D%size()))
    allocate(fields2D(fld_data3D%size()))

    call fld_data2D%init_n_fields(fld_data3D%size(),n_2d)
    call fld_data3D%get_list(fields3D,fld_data3D%size())
    
    this%u = 0.0_rp
    this%old_u = 0.0_rp
    this%avg_u = 0.0_rp
    do i = 1, fld_data3D%size()
       call copy(this%old_u%x,fields3D(i)%ptr%x,n)
       call perform_local_summation(this%u,this%old_u,&
            this%el_heights, this%domain_height, &
            this%map_1d%dir_el, this%coef, this%msh%nelv, lx)
       call copy(this%old_u%x,this%u%x,n)
       call copy(this%avg_u%x,this%u%x,n)
       call perform_global_summation(this%u, this%avg_u, &
            this%old_u, this%map_1d%n_el_lvls, &
            this%map_1d%dir_el,this%coef%gs_h,&
            this%coef%mult, this%msh%nelv, lx)
       call copy(fields3D(i)%ptr%x,this%u%x,n)
    end do
    call fld_data2D%get_list(fields2d,fld_data2D%size())
    do i = 1, fld_data3D%size()
       do j = 1, n_2d
          fields2d(i)%ptr%x(j) = fields3D(i)%ptr%x(this%idx_2d(j))
       end do
    end do
  end subroutine map_2d_average
 
  subroutine perform_global_summation(u, avg_u, old_u, n_levels, &
                                     hom_dir_el, gs_h, mult, nelv, lx)
    type(field_t), intent(inout) :: u, avg_u, old_u
    type(gs_t), intent(inout) :: gs_h
    integer, intent(in) :: n_levels, nelv, lx
    integer, intent(in) :: hom_dir_el(nelv)
    real(kind=rp), intent(in) :: mult(nelv*lx**3)
    real(kind=rp) :: temp_el(lx,lx,lx)
    integer :: n, i, j, e, k
    type(c_ptr) :: ptr
  
    n = u%dof%size()
  
    do i = 1, n_levels-1
       !compute average
       if (NEKO_BCKND_DEVICE .eq. 1) &
            call device_memcpy(u%x, u%x_d, n, &
                               HOST_TO_DEVICE, sync=.true.)
       call gs_h%op(u,GS_OP_ADD)
       if (NEKO_BCKND_DEVICE .eq. 1) &
            call device_memcpy(u%x, u%x_d, n, DEVICE_TO_HOST, sync=.true.)
       call col2(u%x,mult,n)
       !Assumes sugarcube I think
       do e = 1, nelv
          temp_el = 2.0*u%x(:,:,:,e)-old_u%x(:,:,:,e)
          if (hom_dir_el(e) .eq. 1) then
             u%x(1,:,:,e) = temp_el(lx,:,:)
             avg_u%x(1,:,:,e) = avg_u%x(1,:,:,e)+temp_el(1,:,:)
             u%x(lx,:,:,e) = temp_el(1,:,:)
          else if (hom_dir_el(e) .eq. 2) then
             u%x(:,1,:,e) = temp_el(:,lx,:)
             avg_u%x(:,1,:,e) = avg_u%x(:,1,:,e)+temp_el(:,1,:)
             u%x(:,lx,:,e) = temp_el(:,1,:)
          else if (hom_dir_el(e) .eq. 3) then
             u%x(:,:,1,e) = temp_el(:,:,lx)
             avg_u%x(:,:,1,e) = avg_u%x(:,:,1,e)+temp_el(:,:,1)
             u%x(:,:,lx,e) = temp_el(:,:,1)
          end if
          old_u%x(:,:,:,e) = u%x(:,:,:,e)
       end do
    end do
    do e = 1, nelv
       do i = 1,lx
          do j = 1, lx
             if (hom_dir_el(e) .eq. 1) then
                u%x(:,i,j,e) = avg_u%x(1,i,j,e)
             else if (hom_dir_el(e) .eq. 2) then
                u%x(i,:,j,e) = avg_u%x(i,1,j,e)
             else if (hom_dir_el(e) .eq. 3) then
                u%x(i,j,:,e) = avg_u%x(i,j,1,e)
             end if
          end do
       end do
    end do
  end subroutine perform_global_summation
  
  subroutine perform_local_summation(u_out, u, el_heights, domain_height,&
                                     hom_dir_el, coef, nelv, lx)
    type(field_t), intent(inout) :: u, u_out, el_heights
    type(coef_t), intent(inout) :: coef
    integer, intent(in) :: nelv, lx
    integer, intent(in) :: hom_dir_el(nelv)
    real(kind=rp), intent(in) :: domain_height
    real(kind=rp) :: wt
    integer :: n, i, j, e, k
  
    n = nelv*lx**3
  
    call col2(u%x,el_heights%x,n)
    call cmult(u%x, 1.0_rp/(2.0*domain_height),n)
    call rzero(u_out%x,n)
  
  
    do e = 1, nelv
       do i = 1,lx
          do j = 1, lx
             do k = 1, lx
                wt = coef%Xh%wx(k)
                if (hom_dir_el(e) .eq. 1) then
                   u_out%x(1,i,j,e) = u_out%x(1,i,j,e)+wt*u%x(k,i,j,e)
                else if (hom_dir_el(e) .eq. 2) then
                   u_out%x(i,1,j,e) = u_out%x(i,1,j,e)+wt*u%x(i,k,j,e)
                else if (hom_dir_el(e) .eq. 3) then
                   u_out%x(i,j,1,e) = u_out%x(i,j,1,e)+wt*u%x(i,j,k,e)
                end if
             end do
          end do
       end do
    end do
  
    do e = 1, nelv
       do i = 1,lx
          do j = 1, lx
             if (hom_dir_el(e) .eq. 1) then
                u_out%x(:,i,j,e) = u_out%x(1,i,j,e)
             else if (hom_dir_el(e) .eq. 2) then
                u_out%x(i,:,j,e) = u_out%x(i,1,j,e)
             else if (hom_dir_el(e) .eq. 3) then
                u_out%x(i,j,:,e) = u_out%x(i,j,1,e)
             end if
          end do
       end do
    end do
  end subroutine perform_local_summation

end module map_2d
