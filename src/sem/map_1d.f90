
!> Creates a 1d GLL point map along a specified direction based on the connectivity in the mesh.
module map_1d
  use num_types, only: rp
  use space, only: space_t
  use dofmap, only: dofmap_t
  use gather_scatter
  use mesh, only: mesh_t
  use device
  use comm
  use coefs, only: coef_t
  use field_list, only: field_list_t
  use matrix, only: matrix_t
  use vector, only: vector_ptr_t
  use logger, only: neko_log, LOG_SIZE
  use utils, only: neko_error, neko_warning
  use math, only: glmax, glmin, glimax, glimin, relcmp, cmult, add2s1, col2
  use neko_mpi_types
  use, intrinsic :: iso_c_binding
  implicit none
  private
  !> Type that encapsulates a mapping from each gll point in the mesh
  !! to its corresponding (global) GLL point index in one direction.
  !! @remark Could also be rather easily extended to say polar coordinates
  !! as well ( I think). Martin Karp
  type, public :: map_1d_t
     !> Checks whether the specified direction is in the r,s, or t
     !! direction for each element.
     integer, allocatable :: dir_el(:)
     !> Checks which level an element belongs to.
     integer, allocatable :: el_lvl(:)
     !> Checks which level or id in the 1D GLL mapping each point in the dofmap is.
     integer, allocatable :: pt_lvl(:,:,:,:)
     !> Number of elements stacked on top of eachother in the specified direction
     integer :: n_el_lvls
     !> Number of total gll levels
     integer :: n_gll_lvls
     !> Dofmap
     !> Dofmap
     type(dofmap_t), pointer :: dof => null()
     !> SEM coefs
     type(coef_t), pointer :: coef => null()
     !> Mesh
     type(mesh_t), pointer :: msh => null()
     !> The specified direction in which we create the 1D mapping
     integer :: dir
     !> Tolerance for the mesh
     real(kind=rp) :: tol = 1e-7
     !> Volume per level in the 1d grid
     real(kind=rp), allocatable :: volume_per_gll_lvl(:)
   contains
     !> Constructor
     procedure, pass(this) :: init_int => map_1d_init

     procedure, pass(this) :: init_char => map_1d_init_char
     generic :: init => init_int, init_char
     !> Destructor
     procedure, pass(this) :: free => map_1d_free
     !> Average field list along planes
     procedure, pass(this) :: average_planes_fld_lst =>  map_1d_average_field_list

     procedure, pass(this) :: average_planes_vec_ptr =>  map_1d_average_vector_ptr
     generic :: average_planes => average_planes_fld_lst, average_planes_vec_ptr
  end type map_1d_t


contains

  subroutine map_1d_init(this, coef,  dir, tol)
    class(map_1d_t) :: this
    type(coef_t), intent(inout), target :: coef
    integer, intent(in) :: dir
    real(kind=rp), intent(in) :: tol
    integer :: nelv, lx, n, i, e, lvl, ierr
    real(kind=rp), contiguous, pointer :: line(:,:,:,:)
    real(kind=rp), allocatable :: min_vals(:,:,:,:)
    real(kind=rp), allocatable :: min_temp(:,:,:,:)
    type(c_ptr) :: min_vals_d = c_null_ptr
    real(kind=rp) :: el_dim(3,3), glb_min, glb_max, el_min

    call this%free()

    if (NEKO_BCKND_DEVICE .eq. 1) then
       if (pe_rank .eq. 0) then
          call neko_warning('map_1d does not copy indices to device,'// &
                            ' but ok if used on cpu and for io')
       end if
    end if

    this%dir = dir
    this%dof => coef%dof
    this%coef => coef
    this%msh => coef%msh
    nelv = this%msh%nelv
    lx = this%dof%Xh%lx
    n = this%dof%size()

    if (dir .eq. 1) then
       line => this%dof%x
    else if (dir .eq. 2) then
       line => this%dof%y
    else if(dir .eq. 3) then
       line => this%dof%z
    else
       call neko_error('Invalid dir for geopmetric comm')
    end if
    allocate(this%dir_el(nelv))
    allocate(this%el_lvl(nelv))
    allocate(this%pt_lvl(lx, lx, lx, nelv))
    allocate(min_vals(lx, lx, lx, nelv))
    allocate(min_temp(lx, lx, lx, nelv))
    call MPI_BARRIER(NEKO_COMM)
    if (NEKO_BCKND_DEVICE .eq. 1) then

       call device_map(min_vals,min_vals_d,n)
    end if
    call MPI_BARRIER(NEKO_COMM)

    do i = 1, nelv
       !store which direction r,s,t corresponds to speciifed direction, x,y,z
       !we assume elements are stacked on each other...
       ! Check which one of the normalized vectors are closest to dir
       ! If we want to incorporate other directions, we should look here
       el_dim(1,:) = abs(this%msh%elements(i)%e%pts(1)%p%x - &
                     this%msh%elements(i)%e%pts(2)%p%x)
       el_dim(1,:) = el_dim(1,:)/norm2(el_dim(1,:))
       el_dim(2,:) = abs(this%msh%elements(i)%e%pts(1)%p%x - &
                     this%msh%elements(i)%e%pts(3)%p%x)
       el_dim(2,:) = el_dim(2,:)/norm2(el_dim(2,:))
       el_dim(3,:) = abs(this%msh%elements(i)%e%pts(1)%p%x - &
                     this%msh%elements(i)%e%pts(5)%p%x)
       el_dim(3,:) = el_dim(3,:)/norm2(el_dim(3,:))
       ! Checks which directions in rst the xyz corresponds to
       ! 1 corresponds to r, 2 to s, 3 to t and are stored in dir_el
       this%dir_el(i) = maxloc(el_dim(:,this%dir),dim=1)
    end do
    glb_min =  glmin(line,n)
    glb_max =  glmax(line,n)

    i = 1
    this%el_lvl = -1
    ! Check what the mimum value in each element and put in min_vals
    do e = 1, nelv
       el_min = minval(line(:,:,:,e))
       min_vals(:,:,:,e) = el_min
       ! Check if this element is on the bottom, in this case assign el_lvl = i = 1
       if (relcmp(el_min, glb_min, this%tol)) then
          if(this%el_lvl(e) .eq. -1) this%el_lvl(e) = i
       end if
    end do
    ! While loop where at each iteation the global maximum value 
    ! propagates down one level.
    ! When the minumum value has propagated to the highest level this stops.
    ! Only works when the bottom plate of the domain is flat.
    do while (.not. relcmp(glmax(min_vals,n), glb_min, this%tol))
       i = i + 1
       do e = 1, nelv
          !Sets the value at the bottom of each element to glb_max
          if (this%dir_el(e) .eq. 1) then
             if (line(1,1,1,e) .gt. line(lx,1,1,e)) then
                min_vals(lx,:,:,e) = glb_max
             else
                min_vals(1,:,:,e) = glb_max
             end if
          end if
          if (this%dir_el(e) .eq. 2) then
             if (line(1,1,1,e) .gt. line(1,lx,1,e)) then
                min_vals(:,lx,:,e) = glb_max
             else
                min_vals(:,1,:,e) = glb_max
             end if
          end if
          if (this%dir_el(e) .eq. 3) then
             if (line(1,1,1,e) .gt. line(1,1,lx,e)) then
                min_vals(:,:,lx,e) = glb_max
             else
                min_vals(:,:,1,e) = glb_max
             end if
          end if
       end do
       !Make sketchy min as GS_OP_MIN is not supported with device mpi
       min_temp = min_vals
       if (NEKO_BCKND_DEVICE .eq. 1) &
         call device_memcpy(min_vals, min_vals_d, n,&
                            HOST_TO_DEVICE, sync=.false.)
       !Propagates the minumum value along the element boundary.
       call coef%gs_h%op(min_vals, n, GS_OP_ADD)
       if (NEKO_BCKND_DEVICE .eq. 1) &
          call device_memcpy(min_vals, min_vals_d, n,&
                             DEVICE_TO_HOST, sync=.true.)
       !Obtain average along boundary

       call col2(min_vals, coef%mult, n)
       call cmult(min_temp, -1.0_rp, n)
       call add2s1(min_vals, min_temp, 2.0_rp, n)


       !Checks the new minimum value on each element
       !Assign this value to all points in this element in min_val
       !If the element has not already been assinged a level,
       !and it has obtained the minval, set el_lvl = i
       do e = 1, nelv
          el_min = minval(min_vals(:,:,:,e))
          min_vals(:,:,:,e) = el_min
          if (relcmp(el_min, glb_min, this%tol)) then
             if (this%el_lvl(e) .eq. -1) this%el_lvl(e) = i
          end if
       end do
    end do
    this%n_el_lvls = glimax(this%el_lvl,nelv)
    this%n_gll_lvls = this%n_el_lvls*lx
    !Numbers the points in each element based on the element level
    !and its orientation
    do e = 1, nelv
       do i = 1, lx
          lvl = lx * (this%el_lvl(e) - 1) + i
          if (this%dir_el(e) .eq. 1) then
             if (line(1,1,1,e) .gt. line(lx,1,1,e)) then
                this%pt_lvl(lx-i+1,:,:,e) = lvl
             else
                this%pt_lvl(i,:,:,e) = lvl
             end if
          end if
          if (this%dir_el(e) .eq. 2) then
             if (line(1,1,1,e) .gt. line(1,lx,1,e)) then
                this%pt_lvl(:,lx-i+1,:,e) = lvl
             else
                this%pt_lvl(:,i,:,e) = lvl
             end if
          end if
          if (this%dir_el(e) .eq. 3) then
             if (line(1,1,1,e) .gt. line(1,1,lx,e)) then
                this%pt_lvl(:,:,lx-i+1,e) = lvl
             else
                this%pt_lvl(:,:,i,e) = lvl
             end if
          end if
       end do
    end do
    if(allocated(min_vals)) deallocate(min_vals)
    if(c_associated(min_vals_d)) call device_free(min_vals_d)
    if(allocated(min_temp)) deallocate(min_temp)
    allocate(this%volume_per_gll_lvl(this%n_gll_lvls))
    this%volume_per_gll_lvl =0.0_rp
    do i = 1, n
       this%volume_per_gll_lvl(this%pt_lvl(i,1,1,1)) = &
       this%volume_per_gll_lvl(this%pt_lvl(i,1,1,1)) + coef%B(i,1,1,1)
    end do
    call MPI_Allreduce(MPI_IN_PLACE,this%volume_per_gll_lvl, this%n_gll_lvls, &
         MPI_REAL_PRECISION, MPI_SUM, NEKO_COMM, ierr)

  end subroutine map_1d_init

  subroutine map_1d_init_char(this, coef,  dir, tol)
    class(map_1d_t) :: this
    type(coef_t), intent(inout), target :: coef
    character(len=*), intent(in) :: dir
    real(kind=rp), intent(in) :: tol
    integer :: idir

    if (trim(dir) .eq. 'yz' .or. trim(dir) .eq. 'zy') then
       idir = 1
    else if (trim(dir) .eq. 'xz' .or. trim(dir) .eq. 'zx') then
       idir = 2
    else if (trim(dir) .eq. 'xy' .or. trim(dir) .eq. 'yx') then
       idir = 3
    else
       call neko_error('homogenous direction not supported')
    end if

    call  this%init(coef,idir,tol)

  end subroutine map_1d_init_char

  subroutine map_1d_free(this)
    class(map_1d_t) :: this

    if(allocated(this%dir_el)) deallocate(this%dir_el)
    if(allocated(this%el_lvl)) deallocate(this%el_lvl)
    if(allocated(this%pt_lvl)) deallocate(this%pt_lvl)
    if(associated(this%dof)) nullify(this%dof)
    if(associated(this%msh)) nullify(this%msh)
    if(associated(this%coef)) nullify(this%coef)
    if(allocated(this%volume_per_gll_lvl)) deallocate(this%volume_per_gll_lvl)
    this%dir = 0
    this%n_el_lvls = 0
    this%n_gll_lvls = 0

  end subroutine map_1d_free


  !> Computes average if field list in two directions and outputs matrix
  !! with averaged values
  !! avg_planes contains coordinates in first row, avg. of fields in the rest
  !! @param avg_planes output averages
  !! @param field_list list of fields to be averaged
  subroutine map_1d_average_field_list(this, avg_planes, field_list)
    class(map_1d_t), intent(inout) :: this
    type(field_list_t), intent(inout) :: field_list
    type(matrix_t), intent(inout) :: avg_planes
    integer :: n, ierr, j, i
    real(kind=rp) :: coord
    call avg_planes%free()
    call avg_planes%init(this%n_gll_lvls, field_list%size()+1)
    avg_planes = 0.0_rp
    !ugly way of getting coordinates, computes average
    n = this%dof%size()
    do i = 1, n
       if (this%dir .eq. 1) coord = this%dof%x(i,1,1,1)
       if (this%dir .eq. 2) coord = this%dof%y(i,1,1,1)
       if (this%dir .eq. 3) coord = this%dof%z(i,1,1,1)
       avg_planes%x(this%pt_lvl(i,1,1,1),1) = &
       avg_planes%x(this%pt_lvl(i,1,1,1),1) + coord*this%coef%B(i,1,1,1) &
       /this%volume_per_gll_lvl(this%pt_lvl(i,1,1,1))
    end do
    do j = 2, field_list%size() + 1
       do i = 1, n
          avg_planes%x(this%pt_lvl(i,1,1,1),j) = &
          avg_planes%x(this%pt_lvl(i,1,1,1),j) + &
          field_list%items(j-1)%ptr%x(i,1,1,1)*this%coef%B(i,1,1,1) &
          /this%volume_per_gll_lvl(this%pt_lvl(i,1,1,1))
       end do
    end do
    if (pe_size .gt. 1) then
       call MPI_Allreduce(MPI_IN_PLACE, avg_planes%x, (field_list%size()+1)*this%n_gll_lvls, &
            MPI_REAL_PRECISION, MPI_SUM, NEKO_COMM, ierr)
    end if


  end subroutine map_1d_average_field_list
  
  !> Computes average if vector_pt in two directions and outputs matrix
  !! with averaged values
  !! avg_planes contains coordinates in first row, avg. of fields in the rest
  !! @param avg_planes output averages
  !! @param vector_pts to vectors to be averaged
  subroutine map_1d_average_vector_ptr(this, avg_planes, vector_ptr)
    class(map_1d_t), intent(inout) :: this
    !Observe is an array...
    type(vector_ptr_t), intent(inout) :: vector_ptr(:)
    type(matrix_t), intent(inout) :: avg_planes
    integer :: n, ierr, j, i
    real(kind=rp) :: coord

    call avg_planes%free()
    call avg_planes%init(this%n_gll_lvls,size(vector_ptr)+1)
    !ugly way of getting coordinates, computes average
    avg_planes = 0.0_rp

    n = this%dof%size()
    do i = 1, n
       if (this%dir .eq. 1) coord = this%dof%x(i,1,1,1)
       if (this%dir .eq. 2) coord = this%dof%y(i,1,1,1)
       if (this%dir .eq. 3) coord = this%dof%z(i,1,1,1)
       avg_planes%x(this%pt_lvl(i,1,1,1),1) = &
       avg_planes%x(this%pt_lvl(i,1,1,1),1) + coord*this%coef%B(i,1,1,1) &
       /this%volume_per_gll_lvl(this%pt_lvl(i,1,1,1))
    end do
    do j = 2, size(vector_ptr)+1
       do i = 1, n
          avg_planes%x(this%pt_lvl(i,1,1,1),j) = &
          avg_planes%x(this%pt_lvl(i,1,1,1),j) + &
          vector_ptr(j-1)%ptr%x(i)*this%coef%B(i,1,1,1) &
          /this%volume_per_gll_lvl(this%pt_lvl(i,1,1,1))
       end do
    end do
    call MPI_Allreduce(MPI_IN_PLACE,avg_planes%x, (size(vector_ptr)+1)*this%n_gll_lvls, &
       MPI_REAL_PRECISION, MPI_SUM, NEKO_COMM, ierr)


  end subroutine map_1d_average_vector_ptr

end module map_1d
