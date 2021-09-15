!> Defines a function space
module space
  use num_types
  use speclib
  use utils
  use math
  use mamba
  implicit none

  integer, parameter :: GL = 0, GLL = 1, GJ = 2

  type space_t
     integer :: lx              !< Polynomial dimension in x-direction
     integer :: ly              !< Polynomial dimension in y-direction
     integer :: lz              !< Polynomial dimension in z-direction
     integer :: lxy             !< Number of points in xy-plane
     integer :: lyz             !< Number of points in yz-plane
     integer :: lxz             !< Number of points in xz-plane
     integer :: lxyz            !< Number of points in xyz-block
     
     real(kind=rp), allocatable :: zg(:,:) !< Quadrature points
     
     real(kind=rp), allocatable :: dr_inv(:) !< 1/dist quadrature points
     real(kind=rp), allocatable :: ds_inv(:) !< 1/dist quadrature points
     real(kind=rp), allocatable :: dt_inv(:) !< 1/dist quadrature points

     real(kind=rp), allocatable :: wx(:)   !< Quadrature weights
     real(kind=rp), allocatable :: wy(:)   !< Quadrature weights
     real(kind=rp), allocatable :: wz(:)   !< Quadrature weights

     real(kind=rp), allocatable :: w3(:,:,:)

     real(kind=rp), allocatable :: dx(:,:) !< Derivative operator
     real(kind=rp), allocatable :: dy(:,:) !< Derivative operator
     real(kind=rp), allocatable :: dz(:,:) !< Derivative operator

     real(kind=rp), allocatable :: dxt(:,:) !< Derivative operator
     real(kind=rp), allocatable :: dyt(:,:) !< Derivative operator
     real(kind=rp), allocatable :: dzt(:,:) !< Derivative operator

     type(mmbArray) :: mba_dx
     type(mmbArray) :: mba_dxt
     type(mmbLayout) :: layout
     type(mmbTileIterator) :: mba_dx_it
     type(mmbTileIterator) :: mba_dxt_it
  end type space_t

  interface operator(.eq.)
     module procedure space_eq
  end interface operator(.eq.)

  interface operator(.ne.)
     module procedure space_ne
  end interface operator(.ne.)
  
contains

  !> Initialize a function space @a s with given polynomial dimensions
  subroutine space_init(s, t, lx, ly, lz)
    type(space_t), intent(inout) :: s
    integer, intent(in) :: t            !< Quadrature type
    integer, intent(in) :: lx           !< Polynomial dimension in x-direction
    integer, intent(in) :: ly           !< Polynomial dimension in y-direction
    integer, optional, intent(in) :: lz !< Polynomial dimension in z-direction
    integer :: ix, iy, iz
    integer(mmbErrorKind) :: err
    integer(mmbIndexKind), dimension(2) :: dims

    call space_free(s)

    s%lx = lx
    s%ly = ly
    if (present(lz)) then
       s%lz = lz
       if (lx .ne. ly .or. lx .ne. lz) then
          call neko_error("Unsupported polynomial dimension")
       end if
    else
       if (lx .ne. ly) then
          call neko_error("Unsupported polynomial dimension")
       end if
       s%lz = 1
    end if
    s%lxy = s%ly*s%lx
    s%lyz = s%ly*s%lz
    s%lxz = s%lx*s%lz
    s%lxyz = s%lx*s%ly*s%lz

    allocate(s%zg(lx, 3))

    allocate(s%wx(s%lx))
    allocate(s%wy(s%ly))
    allocate(s%wz(s%lz))
    
    allocate(s%dr_inv(s%lx))
    allocate(s%ds_inv(s%ly))
    allocate(s%dt_inv(s%lz))

    allocate(s%w3(s%lx, s%ly, s%lz))

    allocate(s%dx(s%lx, s%lx))
    allocate(s%dy(s%ly, s%ly))
    allocate(s%dz(s%lz, s%lz))

    allocate(s%dxt(s%lx, s%lx))
    allocate(s%dyt(s%ly, s%ly))
    allocate(s%dzt(s%lz, s%lz))
    
    if (t .eq. GLL) then
       call zwgll(s%zg(1,1), s%wx, s%lx)
       call zwgll(s%zg(1,2), s%wy, s%ly)
       if (s%lz .gt. 1) then
          call zwgll(s%zg(1,3), s%wz, s%lz)
       else
          s%zg(:,3) = 0d0
          s%wz = 1d0
       end if
    else if (t .eq. GL) then
       call zwgl(s%zg(1,1), s%wx, s%lx)
       call zwgl(s%zg(1,2), s%wy, s%ly)
       if (s%lz .gt. 1) then
          call zwgl(s%zg(1,3), s%wz, s%lz)
       else
          s%zg(:,3) = 0d0
          s%wz = 1d0
       end if
    else
       call neko_error("Invalid quadrature rule")
    end if

    do iz = 1, s%lz
       do iy = 1, s%ly
          do ix = 1, s%lx
             s%w3(ix, iy, iz) = s%wx(ix) * s%wy(iy) * s%wz(iz)
          end do
       end do
    end do

    call dgll(s%dx, s%dxt, s%zg(1,1), s%lx, s%lx)
    call dgll(s%dy, s%dyt, s%zg(1,2), s%ly, s%ly)
    if (s%lz .gt. 1) then
       call dgll(s%dz, s%dzt, s%zg(1,3), s%lz, s%lz)
    else
       s%dz = 0d0
       s%dzt = 0d0
    end if
    
    call space_compute_dist(s%dr_inv, s%zg(1,1), s%lx)
    call space_compute_dist(s%ds_inv, s%zg(1,2), s%ly)
    if (s%lz .gt. 1) then
       call space_compute_dist(s%dt_inv, s%zg(1,3), s%lz)
    else
       s%dt_inv = 0d0
    end if
    
    dims = [s%lx, s%lx]
    call mmb_layout_create_regular_nd(int(storage_size(1.0_rp)/8, mmbSizeKind), &
         2_mmbSizeKind, MMB_COLMAJOR, mmb_layout_padding_create_zero(),&
         s%layout, err)

    call mmb_array_create_wrapped(s%dx, dims, s%layout, &
         dram_interface, MMB_READ_WRITE, s%mba_dx, err)
    
    call mmb_array_create_wrapped(s%dxt, dims, s%layout, &
         dram_interface, MMB_READ_WRITE, s%mba_dxt, err)

    call mmb_array_tile(s%mba_dx, dims, err)
    call mmb_array_tile(s%mba_dxt, dims, err)

    call mmb_tile_iterator_create(s%mba_dx, s%mba_dx_it, err)
    call mmb_tile_iterator_create(s%mba_dxt, s%mba_dxt_it, err)

  end subroutine space_init
   
  !> Deallocate a space @a s
  subroutine space_free(s)
    type(space_t), intent(inout) :: s
    
    if (allocated(s%zg)) then
       deallocate(s%zg)
    end if

    if (allocated(s%wx)) then
       deallocate(s%wx)
    end if

    if (allocated(s%wy)) then
       deallocate(s%wy)
    end if

    if (allocated(s%wz)) then
       deallocate(s%wz)
    end if

    if (allocated(s%w3)) then
       deallocate(s%w3)
    end if

    if (allocated(s%dx)) then
       deallocate(s%dx)
    end if

    if (allocated(s%dy)) then
       deallocate(s%dy)
    end if

    if (allocated(s%dz)) then
       deallocate(s%dz)
    end if

    if (allocated(s%dxt)) then
       deallocate(s%dxt)
    end if

    if (allocated(s%dyt)) then
       deallocate(s%dyt)
    end if

    if (allocated(s%dzt)) then
       deallocate(s%dzt)
    end if
    
    if (allocated(s%dr_inv)) then
       deallocate(s%dr_inv)
    end if
    
    if (allocated(s%ds_inv)) then
       deallocate(s%ds_inv)
    end if
    
    if (allocated(s%dt_inv)) then
       deallocate(s%dt_inv)
    end if

  end subroutine space_free

  !> Check if \f$ X_h = Y_H \f$
  !! @note this only checks the polynomial dimensions
  pure function space_eq(Xh, Yh) result(res)
    type(space_t), intent(in) :: Xh
    type(space_t), intent(in) :: Yh
    logical :: res

    if ( (Xh%lx .eq. Yh%lx) .and. &
         (Xh%ly .eq. Yh%ly) .and. &
         (Xh%lz .eq. Yh%lz) ) then
       res = .true.
    else
       res = .false.
    end if
    
  end function space_eq

  !> Check if \f$ X_h \ne Y_H \f$
  !! @note this only checks the polynomial dimensions
  pure function space_ne(Xh, Yh) result(res)
    type(space_t), intent(in) :: Xh
    type(space_t), intent(in) :: Yh
    logical :: res

    if ( (Xh%lx .eq. Yh%lx) .and. &
         (Xh%ly .eq. Yh%ly) .and. &
         (Xh%lz .eq. Yh%lz) ) then
       res = .false.
    else
       res = .true.
    end if
    
  end function space_ne
  
  subroutine space_compute_dist(dx, x, lx)
    integer, intent(in) :: lx
    real(kind=rp), intent(inout) :: dx(lx), x(lx)
    integer :: i
    dx(1) = x(2) - x(1)
    do i = 2, lx - 1
       dx(i) = 0.5*(x(i+1) - x(i-1))
    enddo
    dx(lx) = x(lx) - x(lx-1)
    call invcol1(dx, lx)
  end subroutine space_compute_dist

end module space
