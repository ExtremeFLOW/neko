
!> Creates a 1d GLL point map along aspecified direction based on the connectivity in the mesh.
module map_1d
  use num_types, only: rp
  use space, only: space_t
  use dofmap, only: dofmap_t
  use gather_scatter
  use mesh, only: mesh_t
  use comm
  use logger, only: neko_log, LOG_SIZE
  use utils, only: neko_error
  use math, only: glmax, glmin, glimax, glimin, relcmp
  use neko_mpi_types
  use, intrinsic :: iso_c_binding
  implicit none
  private
  !> Type that encapsulates a mapping from each gll point in the mesh 
  !! to its corresponding (global) GLL point index in one direction.
  !! I envision this could also be rather easily extended to say polar coordinates as well.
  type, public :: map_1d_t
     !> Checks whether the specified direction is in the r,s, or t direction for each element.
     integer, allocatable :: dir_el(:)
     !> Checks which level an element belongs to.
     integer, allocatable :: el_lvl(:)
     !> Checks which level or id in the 1D GLL mapping each point in the dofmap is.
     integer, allocatable :: pt_lvl(:,:,:,:) 
     !> Number of elements stacked on top of eachother in the specified direction
     integer :: n_el_lvls
     !> Dofmap
     type(dofmap_t), pointer :: dof => null()
     !> Mesh
     type(mesh_t), pointer :: msh => null()
     !> The specified direction in which we create the 1D mapping
     integer :: dir
     !> Tolerance for the mesh
     real(kind=rp) :: tol = 1e-7
   contains
     procedure, pass(this) :: init => map_1d_init
     procedure, pass(this) :: free => map_1d_free
  end type map_1d_t


contains
  
  subroutine map_1d_init(this, dof, gs, dir, tol)
    class(map_1d_t) :: this
    type(dofmap_t), target, intent(in) :: dof
    type(gs_t), intent(inout) :: gs
    integer, intent(in) :: dir
    real(kind=rp), intent(in) :: tol
    integer :: nelv, lx, n, i, e, lvl
    real(kind=rp), pointer :: line(:,:,:,:)
    real(kind=rp), allocatable :: min_vals(:,:,:,:)
    real(kind=rp) :: el_dim(3,3), glb_min, glb_max, el_min
    call this%free()

    if (NEKO_BCKND_DEVICE .eq. 1) then
       call neko_error('map_1d not yet supported on device')
    end if
    
    this%dir = dir
    this%dof => dof
    this%msh => dof%msh
    nelv = this%msh%nelv
    lx = this%dof%Xh%lx
    n = dof%size()
  
    if (dir .eq. 1) then
        line => dof%x
    else if (dir .eq. 2) then
        line => dof%y
    else if(dir .eq. 3) then
        line => dof%z
    else
        call neko_error('Invalid dir for geopmetric comm')
    end if
    
    allocate(this%dir_el(nelv))
    allocate(this%el_lvl(nelv))
    allocate(min_vals(lx, lx, lx, nelv))
    allocate(this%pt_lvl(lx, lx, lx, nelv))

    do i = 1, nelv
       !store which direction r,s,t corresponds to speciifed direction, x,y,z
       !we assume elements are stacked on eachother...
       ! Check which one of the normalized vectors are closest to dir
       ! If we want to incorporate other directions, we should look here
       el_dim(1,:) = abs(this%msh%elements(i)%e%pts(1)%p%x-this%msh%elements(i)%e%pts(2)%p%x)
       el_dim(1,:) = el_dim(1,:)/norm2(el_dim(1,:))
       el_dim(2,:) = abs(this%msh%elements(i)%e%pts(1)%p%x-this%msh%elements(i)%e%pts(3)%p%x)
       el_dim(2,:) = el_dim(2,:)/norm2(el_dim(2,:))
       el_dim(3,:) = abs(this%msh%elements(i)%e%pts(1)%p%x-this%msh%elements(i)%e%pts(5)%p%x)
       el_dim(3,:) = el_dim(3,:)/norm2(el_dim(3,:))
       ! Checks which directions in rst the xyz corresponds to
       ! 1 corresponds to r, 2 to s, 3 to t and are stored in dir_el
       this%dir_el(i) = maxloc(el_dim(:,this%dir),dim=1)
    end do
    glb_min =  glmin(line,n)
    glb_max =  glmax(line,n)

    i = 1
    this%el_lvl = -1
    do e = 1, nelv
       el_min = minval(line(:,:,:,e))
       min_vals(:,:,:,e) = el_min
       if (relcmp(el_min,glb_min,this%tol)) then
          if(this%el_lvl(e) .eq. -1) this%el_lvl(e) = i
       end if
    end do
    do while (.not. relcmp(glmax(min_vals,n), glb_min, this%tol))
      i = i + 1
      do e = 1, nelv
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
      call gs%op(min_vals,n,GS_OP_MIN)
      do e = 1, nelv
         el_min = minval(min_vals(:,:,:,e))
         min_vals(:,:,:,e) = el_min
         if (relcmp(el_min, glb_min, this%tol)) then
            if (this%el_lvl(e) .eq. -1) this%el_lvl(e) = i
         end if
      end do
    end do
    this%n_el_lvls = glimax(this%el_lvl,nelv)
    write(*,*) 'number of element levels', this%n_el_lvls
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

  end subroutine map_1d_init

  subroutine map_1d_free(this)
    class(map_1d_t) :: this

    if(allocated(this%dir_el)) deallocate(this%dir_el)
    if(allocated(this%el_lvl)) deallocate(this%el_lvl)
    if(allocated(this%pt_lvl)) deallocate(this%pt_lvl)
    if(associated(this%dof)) nullify(this%dof)
    if(associated(this%msh)) nullify(this%msh)
    this%dir = 0
    this%n_el_lvls = 0

  end subroutine map_1d_free

end module map_1d
