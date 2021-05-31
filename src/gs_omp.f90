!> Generic Gather-scatter backend for accelerators using OpenMP target offload
module gs_omp
  use num_types
  use gs_bcknd
  use gs_ops
  implicit none
  private

  !> Gather-scatter backend for accelerators 
  type, public, extends(gs_bcknd_t) :: gs_omp_t
     real(kind=rp), allocatable :: local_wrk(:)
     real(kind=rp), allocatable :: shared_wrk(:)
     integer :: nlocal
     integer :: nshared
   contains
     procedure, pass(this) :: init => gs_omp_init
     procedure, pass(this) :: free => gs_omp_free
     procedure, pass(this) :: gather => gs_gather_omp
     procedure, pass(this) :: scatter => gs_scatter_omp
  end type gs_omp_t
  
contains
  
  !> Dummy backend initialisation
  subroutine gs_omp_init(this, nlocal, nshared)
    class(gs_omp_t), intent(inout) :: this
    integer, intent(in) :: nlocal
    integer, intent(in) :: nshared

    call this%free()

    this%nlocal = nlocal
    this%nshared = nshared

    allocate(this%local_wrk(nlocal))
    allocate(this%shared_wrk(nshared))

    !$omp target enter data map(alloc: this%local_wrk, this%shared_wrk)
    
  end subroutine gs_omp_init

  !> Dummy backend deallocation
  subroutine gs_omp_free(this)
    class(gs_omp_t), intent(inout) :: this

    if (allocated(this%local_wrk)) then
       !$omp target exit data map(delete: this%local_wrk)
       deallocate(this%local_wrk)
    end if

    if (allocated(this%shared_wrk)) then
       !$omp target exit data map(delete: this%shared_wrk)
       deallocate(this%shared_wrk)
    end if

    this%nlocal = 0
    this%nshared = 0
    
  end subroutine gs_omp_free

  !> Gather kernel
  subroutine gs_gather_omp(this, v, m, o, dg, u, n, gd, nb, b, op)
    integer, intent(inout) :: m
    integer, intent(inout) :: n
    integer, intent(inout) :: nb
    class(gs_omp_t), intent(inout) :: this
    real(kind=rp), dimension(m), intent(inout) :: v
    integer, dimension(m), intent(inout) :: dg
    real(kind=rp), dimension(n), intent(inout) :: u
    integer, dimension(m), intent(inout) :: gd
    integer, dimension(nb), intent(inout) :: b
    integer, intent(inout) :: o
    integer :: op

    if (this%nlocal .eq. m) then
       associate(w=>this%local_wrk)
         select case(op)
         case (GS_OP_ADD)
            call gs_gather_kernel_add(v, m, o, dg, u, n, gd, nb, b, w)
         case (GS_OP_MUL)
            call gs_gather_kernel_mul(v, m, o, dg, u, n, gd, nb, b, w)
         case (GS_OP_MIN)
            call gs_gather_kernel_min(v, m, o, dg, u, n, gd, nb, b, w)
         case (GS_OP_MAX)
            call gs_gather_kernel_max(v, m, o, dg, u, n, gd, nb, b, w)
         end select
       end associate
    else if (this%nshared .eq. m) then
       associate(w=>this%shared_wrk)
         select case(op)
         case (GS_OP_ADD)
            call gs_gather_kernel_add(v, m, o, dg, u, n, gd, nb, b, w)
         case (GS_OP_MUL)
            call gs_gather_kernel_mul(v, m, o, dg, u, n, gd, nb, b, w)
         case (GS_OP_MIN)
            call gs_gather_kernel_min(v, m, o, dg, u, n, gd, nb, b, w)
         case (GS_OP_MAX)
            call gs_gather_kernel_max(v, m, o, dg, u, n, gd, nb, b, w)
         end select
       end associate
    end if
    
  end subroutine gs_gather_omp
 
  !> Gather kernel for addition of data
  !! \f$ v(dg(i)) = v(dg(i)) + u(gd(i)) \f$
  subroutine gs_gather_kernel_add(v, m, o, dg, u, n, gd, nb, b, w)
    integer, intent(in) :: m
    integer, intent(in) :: n
    integer, intent(in) :: nb
    real(kind=rp), dimension(m), intent(inout) :: v
    real(kind=rp), dimension(m), intent(inout) :: w
    integer, dimension(m), intent(inout) :: dg
    real(kind=rp), dimension(n), intent(inout) :: u
    integer, dimension(m), intent(inout) :: gd
    integer, dimension(nb), intent(inout) :: b
    integer, intent(in) :: o
    integer :: i, j, k, blk_len
    real(kind=rp) :: tmp

    !$omp target teams distribute parallel do
    do i = 1, abs(o) - 1
       v(i) = 0d0
       w(i) = u(gd(i))
    end do
    !$omp end target teams distribute parallel do

    !$omp target teams distribute parallel do
    do i = 1, abs(o) - 1
       !$omp atomic
       v(dg(i)) = v(dg(i)) + w(i)
    end do
    !$omp end target teams distribute parallel do
    
    if (o .lt. 0) then
       !$omp target teams distribute parallel do
       do i = abs(o), m
          v(dg(i)) = u(gd(i))
       end do
       !$omp end target teams distribute parallel do
    else
       !$omp target teams distribute parallel do private(tmp)
       do i = o, m, 2
          tmp  = u(gd(i)) + u(gd(i+1))
          v(dg(i)) = tmp
       end do
       !$omp end target teams distribute parallel do
    end if
    
  end subroutine gs_gather_kernel_add

  !> Gather kernel for multiplication of data
  !! \f$ v(dg(i)) = v(dg(i)) \cdot u(gd(i)) \f$
  subroutine gs_gather_kernel_mul(v, m, o, dg, u, n, gd, nb, b, w)
    integer, intent(in) :: m
    integer, intent(in) :: n
    integer, intent(in) :: nb
    real(kind=rp), dimension(m), intent(inout) :: v
    real(kind=rp), dimension(m), intent(inout) :: w
    integer, dimension(m), intent(inout) :: dg
    real(kind=rp), dimension(n), intent(inout) :: u
    integer, dimension(m), intent(inout) :: gd
    integer, dimension(nb), intent(inout) :: b
    integer, intent(in) :: o
    integer :: i, j, k, blk_len
    real(kind=rp) :: tmp
    
    !$omp target teams distribute parallel do
    do i = 1, abs(o) - 1
       w(i) = u(gd(i))
       v(i) = 1d0
    end do
    !$omp end target teams distribute parallel do

    !$omp target teams distribute parallel do
    do i = 1, abs(o) - 1
       !$omp atomic
       v(dg(i)) = v(dg(i)) * w(i)
    end do
    !$omp end target teams distribute parallel do
    
    if (o .lt. 0) then
       !$omp target teams distribute parallel do
       do i = abs(o), m          
          v(dg(i)) = u(gd(i))
       end do
       !$omp end target teams distribute parallel do
    else
       !$omp target teams distribute parallel do private(tmp)
       do i = o, m, 2
          tmp  = u(gd(i)) * u(gd(i+1))
          v(dg(i)) = tmp
       end do
       !$omp end target teams distribute parallel do
    end if
    
  end subroutine gs_gather_kernel_mul
  
  !> Gather kernel for minimum of data
  !! \f$ v(dg(i)) = \min(v(dg(i)), u(gd(i))) \f$
  subroutine gs_gather_kernel_min(v, m, o, dg, u, n, gd, nb, b, w)
    integer, intent(in) :: m
    integer, intent(in) :: n
    integer, intent(in) :: nb
    real(kind=rp), dimension(m), intent(inout) :: v
    real(kind=rp), dimension(m), intent(inout) :: w
    integer, dimension(m), intent(inout) :: dg
    real(kind=rp), dimension(n), intent(inout) :: u
    integer, dimension(m), intent(inout) :: gd
    integer, dimension(nb), intent(inout) :: b
    integer, intent(in) :: o
    integer :: i, j, k, blk_len
    real(kind=rp) :: tmp

    !$omp target teams distribute parallel do
    do i = 1, abs(o) - 1
       w(i) = u(gd(i))
       v(i) = 1d32
    end do
    !$omp end target teams distribute parallel do

    !$omp target teams distribute parallel do
    do i = 1, abs(o) - 1
       !$omp atomic
       v(dg(i)) = min(v(dg(i)), w(i))
    end do
    !$omp end target teams distribute parallel do
    
    if (o .lt. 0) then
       !$omp target teams distribute parallel do
       do i = abs(o), m
          v(dg(i)) = u(gd(i))
       end do
       !$omp end target teams distribute parallel do
    else
       !$omp target teams distribute parallel do private(tmp)
       do i = o, m, 2          
          tmp  = min(u(gd(i)), u(gd(i+1)))
          v(dg(i)) = tmp
       end do
       !$omp end target teams distribute parallel do
    end if
    
  end subroutine gs_gather_kernel_min

  !> Gather kernel for maximum of data
  !! \f$ v(dg(i)) = \max(v(dg(i)), u(gd(i))) \f$
  subroutine gs_gather_kernel_max(v, m, o, dg, u, n, gd, nb, b, w)
    integer, intent(in) :: m
    integer, intent(in) :: n
    integer, intent(in) :: nb
    real(kind=rp), dimension(m), intent(inout) :: v
    real(kind=rp), dimension(m), intent(inout) :: w
    integer, dimension(m), intent(inout) :: dg
    real(kind=rp), dimension(n), intent(inout) :: u
    integer, dimension(m), intent(inout) :: gd
    integer, dimension(nb), intent(inout) :: b
    integer, intent(in) :: o
    integer :: i, j, k, blk_len
    real(kind=rp) :: tmp

    !$omp target teams distribute parallel do
    do i = 1, abs(o) - 1
       w(i) = u(gd(i))
       v(i) = 1d-32
    end do
    !$omp end target teams distribute parallel do

    !$omp target teams distribute parallel do
    do i = 1, abs(o) - 1
       !$omp atomic
       v(dg(i)) = max(v(dg(i)), w(i))
    end do
    !$omp end target teams distribute parallel do
    
    if (o .lt. 0) then
       !$omp target teams distribute parallel do
       do i = abs(o), m
          v(dg(i)) = u(gd(i))
       end do
       !$omp end target teams distribute parallel do
    else
       !$omp target teams distribute parallel do private(tmp)
       do i = o, m, 2
          tmp  = max(u(gd(i)), u(gd(i+1)))
          v(dg(i)) = tmp
       end do
       !$omp end target teams distribute parallel do
    end if
    
  end subroutine gs_gather_kernel_max

  !> Scatter kernel  @todo Make the kernel abstract
  subroutine gs_scatter_omp(this, v, m, dg, u, n, gd, nb, b)
    integer, intent(in) :: m
    integer, intent(in) :: n
    integer, intent(in) :: nb
    class(gs_omp_t), intent(inout) :: this
    real(kind=rp), dimension(m), intent(inout) :: v
    integer, dimension(m), intent(inout) :: dg
    real(kind=rp), dimension(n), intent(inout) :: u
    integer, dimension(m), intent(inout) :: gd
    integer, dimension(nb), intent(inout) :: b
        
    if (this%nlocal .eq. m) then
       call gs_scatter_kernel(v, m, dg, u, n, gd, nb, b, this%local_wrk)
    else if (this%nshared .eq. m) then
       call gs_scatter_kernel(v, m, dg, u, n, gd, nb, b, this%shared_wrk)
    end if

  end subroutine gs_scatter_omp

  !> Scatter kernel \f$ u(gd(i) = v(dg(i)) \f$
  subroutine gs_scatter_kernel(v, m, dg, u, n, gd, nb, b, w)
    integer, intent(in) :: m
    integer, intent(in) :: n
    integer, intent(in) :: nb
    real(kind=rp), dimension(m), intent(inout) :: v
    real(kind=rp), dimension(m), intent(inout) :: w
    integer, dimension(m), intent(inout) :: dg
    real(kind=rp), dimension(n), intent(inout) :: u
    integer, dimension(m), intent(inout) :: gd
    integer, dimension(nb), intent(inout) :: b
    integer :: i, j, k, blk_len
    real(kind=rp) :: tmp

    !$omp target teams distribute parallel do
    do i = 1, m
       w(i) = v(dg(i))
    end do
    !$omp end target teams distribute parallel do

    !$omp target teams distribute parallel do
    do i = 1, m
       u(gd(i)) = w(i)
    end do
    !$omp end target teams distribute parallel do
    
  end subroutine gs_scatter_kernel

end module gs_omp
