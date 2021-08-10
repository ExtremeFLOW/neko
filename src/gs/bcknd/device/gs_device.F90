!> Generic Gather-scatter backend for accelerators
module gs_device
  use num_types
  use gs_bcknd
  use device    
  use gs_ops
  use, intrinsic :: iso_c_binding
  implicit none
  private

  !> Gather-scatter backend for offloading devices
  type, public, extends(gs_bcknd_t) :: gs_device_t
     real(kind=rp), allocatable :: local_wrk(:)  !< Scratch array for local
     real(kind=rp), allocatable :: shared_wrk(:) !< Scratch array for shared
     type(c_ptr) :: local_gs_d      !< Device ptr for local gs-ops
     type(c_ptr) :: local_dof_gs_d  !< Device ptr for local dof to gs mapping
     type(c_ptr) :: local_gs_dof_d  !< Device ptr for local gs to dof mapping
     type(c_ptr) :: shared_gs_d     !< Device ptr for shared gs-ops
     type(c_ptr) :: shared_dof_gs_d !< Device ptr for shared dof to gs mapping
     type(c_ptr) :: shared_gs_dof_d !< Device ptr for shared gs to dof mapping
     type(c_ptr) :: local_wrk_d     !< Device ptr for local scratch array
     type(c_ptr) :: shared_wrk_d    !< Device ptr for shared scratch array
     integer :: nlocal              
     integer :: nshared
   contains
     procedure, pass(this) :: init => gs_device_init
     procedure, pass(this) :: free => gs_device_free
     procedure, pass(this) :: gather => gs_gather_device
     procedure, pass(this) :: scatter => gs_scatter_device
  end type gs_device_t

#ifdef HAVE_HIP
  interface
     subroutine hip_gather_kernel(v, m, o, dg, u, n, gd, w, op) &
          bind(c, name='hip_gather_kernel')
       use, intrinsic :: iso_c_binding
       implicit none
       integer(c_int) :: m, n, o, op
       type(c_ptr), value :: v, w, u, dg, gd
     end subroutine hip_gather_kernel
  end interface

  interface
     subroutine hip_scatter_kernel(v, m, dg, u, n, gd, w) &
          bind(c, name='hip_scatter_kernel')
       use, intrinsic :: iso_c_binding
       implicit none
       integer(c_int) :: m, n
       type(c_ptr), value :: v, w, u, dg, gd
     end subroutine hip_scatter_kernel
  end interface
#endif

contains
  
  !> Accelerator backend initialisation
  subroutine gs_device_init(this, nlocal, nshared)
    class(gs_device_t), intent(inout) :: this
    integer, intent(in) :: nlocal
    integer, intent(in) :: nshared

    call this%free()

    this%nlocal = nlocal
    this%nshared = nshared

    allocate(this%local_wrk(nlocal))
    allocate(this%shared_wrk(nshared))

    call device_map(this%local_wrk, this%local_wrk_d, nlocal)
    call device_map(this%shared_wrk, this%shared_wrk_d, nshared)

    this%local_gs_d = C_NULL_PTR
    this%local_dof_gs_d = C_NULL_PTR
    this%local_gs_dof_d = C_NULL_PTR
    this%shared_gs_d = C_NULL_PTR
    this%shared_dof_gs_d = C_NULL_PTR
    this%shared_gs_dof_d = C_NULL_PTR    
      
  end subroutine gs_device_init

  !> Dummy backend deallocation
  subroutine gs_device_free(this)
    class(gs_device_t), intent(inout) :: this

    if (allocated(this%local_wrk)) then
       deallocate(this%local_wrk)
    end if

    if (allocated(this%shared_wrk)) then
       deallocate(this%shared_wrk)
    end if

    if (c_associated(this%local_gs_d)) then
       call device_free(this%local_gs_d)
    end if

    if (c_associated(this%local_dof_gs_d)) then
       call device_free(this%local_dof_gs_d)
    end if

    if (c_associated(this%local_gs_dof_d)) then
       call device_free(this%local_gs_dof_d)
    end if

    if (c_associated(this%local_wrk_d)) then
       call device_free(this%local_wrk_d)
    end if

    if (c_associated(this%shared_wrk_d)) then
       call device_free(this%shared_wrk_d)
    end if   

    this%nlocal = 0
    this%nshared = 0
    
  end subroutine gs_device_free

  !> Gather kernel
  subroutine gs_gather_device(this, v, m, o, dg, u, n, gd, nb, b, op)
    integer, intent(inout) :: m
    integer, intent(inout) :: n
    integer, intent(inout) :: nb
    class(gs_device_t), intent(inout) :: this
    real(kind=rp), dimension(m), intent(inout) :: v
    integer, dimension(m), intent(inout) :: dg
    real(kind=rp), dimension(n), intent(inout) :: u
    integer, dimension(m), intent(inout) :: gd
    integer, dimension(nb), intent(inout) :: b
    integer, intent(inout) :: o
    integer :: op
    type(c_ptr) :: u_d

    u_d = device_get_ptr(u, n)
        
    if (this%nlocal .eq. m) then       
       associate(v_d=>this%local_gs_d, dg_d=>this%local_dof_gs_d, &
            gd_d=>this%local_gs_dof_d, w_d=>this%local_wrk_d)

         if (.not. c_associated(v_d)) then
            call device_map(v, v_d, m)
         end if

         if (.not. c_associated(dg_d)) then
            call device_map(dg, dg_d, m)
            call device_memcpy(dg, dg_d, m, HOST_TO_DEVICE)
         end if

         if (.not. c_associated(gd_d)) then
            call device_map(gd, gd_d, m)
            call device_memcpy(gd, gd_d, m, HOST_TO_DEVICE)
         end if
         
#ifdef HAVE_HIP
         call hip_gather_kernel(v_d, m, o, dg_d, u_d, n, gd_d, w_d, op)
#else
         call neko_error('No device backend configured')
#endif
         
       end associate
    else if (this%nshared .eq. m) then
       associate(v_d=>this%shared_gs_d, dg_d=>this%shared_dof_gs_d, &
            gd_d=>this%shared_gs_dof_d, w_d=>this%shared_wrk_d)

         if (.not. c_associated(v_d)) then
            call device_map(v, v_d, m)
         end if

         if (.not. c_associated(dg_d)) then
            call device_map(dg, dg_d, m)
            call device_memcpy(dg, dg_d, m, HOST_TO_DEVICE)
         end if

         if (.not. c_associated(gd_d)) then
            call device_map(gd, gd_d, m)
            call device_memcpy(gd, gd_d, m, HOST_TO_DEVICE)
         end if
#ifdef HAVE_HIP   
         call hip_gather_kernel(v_d, m, o, dg_d, u_d, n, gd_d, w_d, op)
#else
         call neko_error('No device backend configured')
#endif
         
       end associate
    end if

  end subroutine gs_gather_device
 
  !> Scatter kernel
  subroutine gs_scatter_device(this, v, m, dg, u, n, gd, nb, b)
    integer, intent(in) :: m
    integer, intent(in) :: n
    integer, intent(in) :: nb
    class(gs_device_t), intent(inout) :: this
    real(kind=rp), dimension(m), intent(inout) :: v
    integer, dimension(m), intent(inout) :: dg
    real(kind=rp), dimension(n), intent(inout) :: u
    integer, dimension(m), intent(inout) :: gd
    integer, dimension(nb), intent(inout) :: b
    type(c_ptr) :: u_d

    u_d = device_get_ptr(u, n)

    if (this%nlocal .eq. m) then
       associate(v_d=>this%local_gs_d, dg_d=>this%local_dof_gs_d, &
            gd_d=>this%local_gs_dof_d, w_d=>this%local_wrk_d)
#ifdef HAVE_HIP
         call hip_scatter_kernel(v_d, m, dg_d, u_d, n, gd_d, w_d)
#else
         call neko_error('No device backend configured')
#endif
       end associate
    else if (this%nshared .eq. m) then
       associate(v_d=>this%shared_gs_d, dg_d=>this%shared_dof_gs_d, &
            gd_d=>this%shared_gs_dof_d, w_d=>this%shared_wrk_d)
#ifdef HAVE_HIP
         call hip_scatter_kernel(v_d, m, dg_d, u_d, n, gd_d, w_d)
#else
         call neko_error('No device backend configured')
#endif
       end associate
    end if

  end subroutine gs_scatter_device

end module gs_device
