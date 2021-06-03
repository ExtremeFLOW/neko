!> Coefficients 
module coefs
  use gather_scatter
  use neko_config
  use num_types
  use space  
  use math
  use mesh
  use device
  use mxm_wrapper
  use, intrinsic :: iso_c_binding
  implicit none
  private
  
  !> Coefficients defined on a given (mesh, \f$ X_h \f$) tuple
  type, public :: coef_t     
     real(kind=rp), allocatable :: G1(:,:,:,:) !< Geometric data
     real(kind=rp), allocatable :: G2(:,:,:,:) !< Geometric data
     real(kind=rp), allocatable :: G3(:,:,:,:) !< Geometric data
     real(kind=rp), allocatable :: G4(:,:,:,:) !< Geometric data
     real(kind=rp), allocatable :: G5(:,:,:,:) !< Geometric data
     real(kind=rp), allocatable :: G6(:,:,:,:) !< Geometric data

     real(kind=rp), allocatable :: mult(:,:,:,:) !< Multiplicity
     ! generate mapping data between element and reference element 
     !! \f$ dx/dr, dy/dr, dz/dr \f$
     !! \f$ dx/ds, dy/ds, dz/ds \f$
     !! \f$ dx/dt, dy/dt, dz/dt \f$
     real(kind=rp), allocatable :: dxdr(:,:,:,:), dydr(:,:,:,:), dzdr(:,:,:,:) 
     real(kind=rp), allocatable :: dxds(:,:,:,:), dyds(:,:,:,:), dzds(:,:,:,:)
     real(kind=rp), allocatable :: dxdt(:,:,:,:), dydt(:,:,:,:), dzdt(:,:,:,:) 
     !< \f$ dr/dx, dr/dy, dr/dz \f$
     !! \f$ ds/dx, ds/dy, ds/dz \f$
     !! \f$ dt/dx, dt/dy, dt/dz \f$
     real(kind=rp), allocatable :: drdx(:,:,:,:), drdy(:,:,:,:), drdz(:,:,:,:) 
     real(kind=rp), allocatable :: dsdx(:,:,:,:), dsdy(:,:,:,:), dsdz(:,:,:,:)
     real(kind=rp), allocatable :: dtdx(:,:,:,:), dtdy(:,:,:,:), dtdz(:,:,:,:) 
     
     real(kind=rp), allocatable :: h1(:,:,:,:) 
     real(kind=rp), allocatable :: h2(:,:,:,:)
     logical :: ifh2
     
     real(kind=rp), allocatable :: jac(:,:,:,:) !< Jacobian
     real(kind=rp), allocatable :: jacinv(:,:,:,:) !< Inverted Jacobian
     real(kind=rp), allocatable :: B(:,:,:,:) !< Mass matrix/volume matrix
     real(kind=rp), allocatable :: Binv(:,:,:,:) !< Inverted Mass matrix/volume matrix

     real(kind=rp), allocatable :: area(:,:,:,:) !< Facet area
     real(kind=rp), allocatable :: nx(:,:,:,:)   !< x-direction of facet normal
     real(kind=rp), allocatable :: ny(:,:,:,:)   !< y-direction of facet normal
     real(kind=rp), allocatable :: nz(:,:,:,:)   !< z-direction of facet normal
     
     real(kind=rp) :: volume
     
     type(space_t), pointer :: Xh => null()
     type(mesh_t), pointer :: msh => null()
     type(dofmap_t), pointer :: dof => null()
     type(gs_t), pointer :: gs_h=> null()

     !
     ! Device pointers (if present)
     ! 
     
     type(c_ptr) :: G1_d, G2_d, G3_d, G4_d, G5_d, G6_d
     type(c_ptr) :: dxdr_d, dydr_d, dzdr_d
     type(c_ptr) :: dxds_d, dyds_d, dzds_d
     type(c_ptr) :: dxdt_d, dydt_d, dzdt_d
     type(c_ptr) :: drdx_d, drdy_d, drdz_d
     type(c_ptr) :: dsdx_d, dsdy_d, dsdz_d
     type(c_ptr) :: dtdx_d, dtdy_d, dtdz_d
     type(c_ptr) :: mult_d, h1_d, h2_d
     type(c_ptr) :: jac_d, jacinv_d, B_d, Binv_d
     type(c_ptr) :: area_d, nx_d, ny_d, nz_d

  end type coef_t

  public :: coef_init, coef_free
  
contains

  !> Initialize coefficients
  subroutine coef_init(coef, gs_h)
    type(coef_t), intent(inout) :: coef
    type(gs_t), intent(inout), target :: gs_h
    integer :: n, m
    integer(c_size_t) :: len
    call coef_free(coef)
    
    coef%msh => gs_h%dofmap%msh
    coef%Xh => gs_h%dofmap%Xh
    coef%dof => gs_h%dofmap
    coef%gs_h => gs_h

    !
    ! Allocate arrays for geometric data
    !
    !>@todo Be clever and try to avoid allocating zeroed geom. factors
    allocate(coef%G1(coef%Xh%lx, coef%Xh%ly, coef%Xh%lz, coef%msh%nelv))
    allocate(coef%G2(coef%Xh%lx, coef%Xh%ly, coef%Xh%lz, coef%msh%nelv))
    allocate(coef%G3(coef%Xh%lx, coef%Xh%ly, coef%Xh%lz, coef%msh%nelv))
    allocate(coef%G4(coef%Xh%lx, coef%Xh%ly, coef%Xh%lz, coef%msh%nelv))
    allocate(coef%G5(coef%Xh%lx, coef%Xh%ly, coef%Xh%lz, coef%msh%nelv))
    allocate(coef%G6(coef%Xh%lx, coef%Xh%ly, coef%Xh%lz, coef%msh%nelv))
    
    allocate(coef%dxdr(coef%Xh%lx, coef%Xh%ly, coef%Xh%lz, coef%msh%nelv))
    allocate(coef%dxds(coef%Xh%lx, coef%Xh%ly, coef%Xh%lz, coef%msh%nelv))
    allocate(coef%dxdt(coef%Xh%lx, coef%Xh%ly, coef%Xh%lz, coef%msh%nelv))
    
    allocate(coef%dydr(coef%Xh%lx, coef%Xh%ly, coef%Xh%lz, coef%msh%nelv))
    allocate(coef%dyds(coef%Xh%lx, coef%Xh%ly, coef%Xh%lz, coef%msh%nelv))
    allocate(coef%dydt(coef%Xh%lx, coef%Xh%ly, coef%Xh%lz, coef%msh%nelv))
    
    allocate(coef%dzdr(coef%Xh%lx, coef%Xh%ly, coef%Xh%lz, coef%msh%nelv))
    allocate(coef%dzds(coef%Xh%lx, coef%Xh%ly, coef%Xh%lz, coef%msh%nelv))
    allocate(coef%dzdt(coef%Xh%lx, coef%Xh%ly, coef%Xh%lz, coef%msh%nelv))
    
    allocate(coef%drdx(coef%Xh%lx, coef%Xh%ly, coef%Xh%lz, coef%msh%nelv))
    allocate(coef%dsdx(coef%Xh%lx, coef%Xh%ly, coef%Xh%lz, coef%msh%nelv))
    allocate(coef%dtdx(coef%Xh%lx, coef%Xh%ly, coef%Xh%lz, coef%msh%nelv))
    
    allocate(coef%drdy(coef%Xh%lx, coef%Xh%ly, coef%Xh%lz, coef%msh%nelv))
    allocate(coef%dsdy(coef%Xh%lx, coef%Xh%ly, coef%Xh%lz, coef%msh%nelv))
    allocate(coef%dtdy(coef%Xh%lx, coef%Xh%ly, coef%Xh%lz, coef%msh%nelv))
    
    allocate(coef%drdz(coef%Xh%lx, coef%Xh%ly, coef%Xh%lz, coef%msh%nelv))
    allocate(coef%dsdz(coef%Xh%lx, coef%Xh%ly, coef%Xh%lz, coef%msh%nelv))
    allocate(coef%dtdz(coef%Xh%lx, coef%Xh%ly, coef%Xh%lz, coef%msh%nelv))
    
    allocate(coef%jac(coef%Xh%lx, coef%Xh%ly, coef%Xh%lz, coef%msh%nelv))
    allocate(coef%jacinv(coef%Xh%lx, coef%Xh%ly, coef%Xh%lz, coef%msh%nelv))
    
    allocate(coef%area(coef%Xh%lx, coef%Xh%ly, 6, coef%msh%nelv))
    allocate(coef%nx(coef%Xh%lx, coef%Xh%ly, 6, coef%msh%nelv))
    allocate(coef%ny(coef%Xh%lx, coef%Xh%ly, 6, coef%msh%nelv))
    allocate(coef%nz(coef%Xh%lx, coef%Xh%ly, 6, coef%msh%nelv))
    
    allocate(coef%B(coef%Xh%lx, coef%Xh%ly, coef%Xh%lz, coef%msh%nelv))
    allocate(coef%Binv(coef%Xh%lx, coef%Xh%ly, coef%Xh%lz, coef%msh%nelv))

    allocate(coef%h1(coef%Xh%lx, coef%Xh%ly, coef%Xh%lz, coef%msh%nelv))
    allocate(coef%h2(coef%Xh%lx, coef%Xh%ly, coef%Xh%lz, coef%msh%nelv))

    allocate(coef%mult(coef%Xh%lx, coef%Xh%ly, coef%Xh%lz, coef%msh%nelv))

    !
    ! Setup device memory (if present)
    !
    
    n = coef%Xh%lx * coef%Xh%ly * coef%Xh%lz * coef%msh%nelv
    if ((NEKO_BCKND_HIP .eq. 1) .or. (NEKO_BCKND_CUDA .eq. 1)) then
       len = n * 8              !<@todo Respect sizeof(real(rp))
       call device_alloc(coef%G1_d, len)
       call device_alloc(coef%G2_d, len)
       call device_alloc(coef%G3_d, len)
       call device_alloc(coef%G4_d, len)
       call device_alloc(coef%G5_d, len)
       call device_alloc(coef%G6_d, len)
       call device_associate(coef%G1, coef%G1_d, n)
       call device_associate(coef%G2, coef%G2_d, n)
       call device_associate(coef%G3, coef%G3_d, n)
       call device_associate(coef%G4, coef%G4_d, n)
       call device_associate(coef%G5, coef%G5_d, n)
       call device_associate(coef%G6, coef%G6_d, n)
       
       call device_alloc(coef%dxdr_d, len)
       call device_alloc(coef%dydr_d, len)
       call device_alloc(coef%dzdr_d, len)
       call device_associate(coef%dxdr, coef%dxdr_d, n)
       call device_associate(coef%dydr, coef%dydr_d, n)
       call device_associate(coef%dzdr, coef%dzdr_d, n)

       call device_alloc(coef%dxds_d, len)
       call device_alloc(coef%dyds_d, len)
       call device_alloc(coef%dzds_d, len)
       call device_associate(coef%dxds, coef%dxds_d, n)
       call device_associate(coef%dyds, coef%dyds_d, n)
       call device_associate(coef%dzds, coef%dzds_d, n)

       call device_alloc(coef%dxdt_d, len)
       call device_alloc(coef%dydt_d, len)
       call device_alloc(coef%dzdt_d, len)
       call device_associate(coef%dxdt, coef%dxdt_d, n)
       call device_associate(coef%dydt, coef%dydt_d, n)
       call device_associate(coef%dzdt, coef%dzdt_d, n)

       call device_alloc(coef%drdx_d, len)
       call device_alloc(coef%drdy_d, len)
       call device_alloc(coef%drdz_d, len)
       call device_associate(coef%drdx, coef%drdx_d, n)
       call device_associate(coef%drdy, coef%drdy_d, n)
       call device_associate(coef%drdz, coef%drdz_d, n)

       call device_alloc(coef%dsdx_d, len)
       call device_alloc(coef%dsdy_d, len)
       call device_alloc(coef%dsdz_d, len)
       call device_associate(coef%dsdx, coef%dsdx_d, n)
       call device_associate(coef%dsdy, coef%dsdy_d, n)
       call device_associate(coef%dsdz, coef%dsdz_d, n)

       call device_alloc(coef%dtdx_d, len)
       call device_alloc(coef%dtdy_d, len)
       call device_alloc(coef%dtdz_d, len)
       call device_associate(coef%dtdx, coef%dtdx_d, n)
       call device_associate(coef%dtdy, coef%dtdy_d, n)
       call device_associate(coef%dtdz, coef%dtdz_d, n)

       call device_alloc(coef%mult_d, len)
       call device_associate(coef%mult, coef%mult_d, n)
       
       call device_alloc(coef%h1_d, len)       
       call device_associate(coef%h1, coef%h1_d, n)

       call device_alloc(coef%h2_d, len)
       call device_associate(coef%h2, coef%h2_d, n)

       call device_alloc(coef%jac_d, len)
       call device_alloc(coef%jacinv_d, len)
       call device_alloc(coef%B_d, len)
       call device_alloc(coef%Binv_d, len)
       call device_associate(coef%jac, coef%jac_d, n)
       call device_associate(coef%jacinv, coef%jacinv_d, n)
       call device_associate(coef%B, coef%B_d, n)
       call device_associate(coef%Binv, coef%Binv_d, n)

       m = coef%Xh%lx * coef%Xh%ly * 6 * coef%msh%nelv
       len = m * 8
       call device_alloc(coef%area_d, len)
       call device_alloc(coef%nx_d, len)
       call device_alloc(coef%ny_d, len)
       call device_alloc(coef%nz_d, len)
       call device_associate(coef%area, coef%area_d, m)
       call device_associate(coef%nx, coef%nx_d, m)
       call device_associate(coef%ny, coef%ny_d, m)
       call device_associate(coef%nz, coef%nz_d, m)
    else
       coef%G1_d = C_NULL_PTR
       coef%G2_d = C_NULL_PTR
       coef%G3_d = C_NULL_PTR
       coef%G4_d = C_NULL_PTR
       coef%G5_d = C_NULL_PTR
       coef%G6_d = C_NULL_PTR
       coef%dxdr_d = C_NULL_PTR
       coef%dydr_d = C_NULL_PTR
       coef%dzdr_d = C_NULL_PTR
       coef%dxds_d = C_NULL_PTR
       coef%dyds_d = C_NULL_PTR
       coef%dzds_d = C_NULL_PTR
       coef%dxdt_d = C_NULL_PTR
       coef%dydt_d = C_NULL_PTR
       coef%dzdt_d = C_NULL_PTR
       coef%drdx_d = C_NULL_PTR
       coef%drdy_d = C_NULL_PTR
       coef%drdz_d = C_NULL_PTR
       coef%dsdx_d = C_NULL_PTR
       coef%dsdy_d = C_NULL_PTR
       coef%dsdz_d = C_NULL_PTR
       coef%dtdx_d = C_NULL_PTR
       coef%dtdy_d = C_NULL_PTR
       coef%dtdz_d = C_NULL_PTR
       coef%mult_d = C_NULL_PTR
       coef%h1_d = C_NULL_PTR
       coef%h2_d = C_NULL_PTR
       coef%jac_d = C_NULL_PTR
       coef%jacinv_d = C_NULL_PTR
       coef%B_d = C_NULL_PTR
       coef%Binv_d = C_NULL_PTR
       coef%area_d = C_NULL_PTR
       coef%nx_d = C_NULL_PTR
       coef%ny_d = C_NULL_PTR
       coef%nz_d = C_NULL_PTR       
  end if

    call coef_generate_dxyzdrst(coef)
    
    call coef_generate_geo(coef)

    call coef_generate_area_and_normal(coef)

    call coef_generate_mass(coef)
    
    ! This is a placeholder, just for now
    ! We can probably find a prettier solution
    call rone(coef%h1,n)
    call rone(coef%h2,n)
    coef%ifh2 = .false.

    if ((NEKO_BCKND_HIP .eq. 1) .or. (NEKO_BCKND_CUDA .eq. 1)) then
       call device_memcpy(coef%h1, coef%h2_d, n, HOST_TO_DEVICE)
       call device_memcpy(coef%h2, coef%h2_d, n, HOST_TO_DEVICE)
    end if
    
    !
    ! Set up multiplicity
    !
    call rone(coef%mult, n)

    !>  @todo cleanup once we have device math in place
    if ((NEKO_BCKND_HIP .eq. 1) .or. (NEKO_BCKND_CUDA .eq. 1)) then
       call device_memcpy(coef%mult, coef%mult_d, n, HOST_TO_DEVICE)
    end if
       
    call gs_op_vector(gs_h, coef%mult, n, GS_OP_ADD)

    if ((NEKO_BCKND_HIP .eq. 1) .or. (NEKO_BCKND_CUDA .eq. 1)) then
       call device_memcpy(coef%mult, coef%mult_d, n, DEVICE_TO_HOST)
    end if
    
    call invcol1(coef%mult, n)

    if ((NEKO_BCKND_HIP .eq. 1) .or. (NEKO_BCKND_CUDA .eq. 1)) then
       call device_memcpy(coef%mult, coef%mult_d, n, HOST_TO_DEVICE)
    end if
    
  end subroutine coef_init

  !> Deallocate coefficients
  subroutine coef_free(coef)
    type(coef_t), intent(inout) :: coef

    if (allocated(coef%G1)) then
       deallocate(coef%G1)
    end if

    if (allocated(coef%G2)) then
       deallocate(coef%G2)
    end if

    if (allocated(coef%G3)) then
       deallocate(coef%G3)
    end if

    if (allocated(coef%G4)) then
       deallocate(coef%G4)
    end if

    if (allocated(coef%G5)) then
       deallocate(coef%G5)
    end if
    
    if (allocated(coef%G6)) then
       deallocate(coef%G6)
    end if

    if (allocated(coef%mult)) then
       deallocate(coef%mult)
    end if
    
    if (allocated(coef%B)) then
       deallocate(coef%B)
    end if
    
    if (allocated(coef%Binv)) then
       deallocate(coef%Binv)
    end if
    
    if(allocated(coef%dxdr)) then
       deallocate(coef%dxdr)
    end if
    
    if(allocated(coef%dxds)) then
       deallocate(coef%dxds)
    end if
    
    if(allocated(coef%dxdt)) then
       deallocate(coef%dxdt)
    end if
    
    if(allocated(coef%dydr)) then
       deallocate(coef%dydr)
    end if
    
    if(allocated(coef%dyds)) then
       deallocate(coef%dyds)
    end if
    
    if(allocated(coef%dydt)) then
       deallocate(coef%dydt)
    end if
    
    if(allocated(coef%dzdr)) then
       deallocate(coef%dzdr)
    end if
    
    if(allocated(coef%dzds)) then
       deallocate(coef%dzds)
    end if
    
    if(allocated(coef%dzdt)) then
       deallocate(coef%dzdt)
    end if
    
    if(allocated(coef%drdx)) then
       deallocate(coef%drdx)
    end if
    
    if(allocated(coef%dsdx)) then
       deallocate(coef%dsdx)
    end if
    
    if(allocated(coef%dtdx)) then
       deallocate(coef%dtdx)
    end if
    
    if(allocated(coef%drdy)) then
       deallocate(coef%drdy)
    end if
    
    if(allocated(coef%dsdy)) then
       deallocate(coef%dsdy)
    end if
    
    if(allocated(coef%dtdy)) then
       deallocate(coef%dtdy)
    end if
    
    if(allocated(coef%drdz)) then
       deallocate(coef%drdz)
    end if
    
    if(allocated(coef%dsdz)) then
       deallocate(coef%dsdz)
    end if
    
    if(allocated(coef%dtdz)) then
       deallocate(coef%dtdz)
    end if
    
    if(allocated(coef%jac)) then
       deallocate(coef%jac)
    end if
    
    if(allocated(coef%jacinv)) then
       deallocate(coef%jacinv)
    end if
    
    if(allocated(coef%h1)) then
       deallocate(coef%h1)
    end if
    
    if(allocated(coef%h2)) then
       deallocate(coef%h2)
    end if

    if (allocated(coef%area)) then
       deallocate(coef%area)
    end if

    if (allocated(coef%nx)) then
       deallocate(coef%nx)
    end if

    if (allocated(coef%ny)) then
       deallocate(coef%ny)
    end if

    if (allocated(coef%nz)) then
       deallocate(coef%nz)
    end if
    
    nullify(coef%msh)
    nullify(coef%Xh)
    nullify(coef%dof)

    !
    ! Cleanup the device (if present)
    !
    
    if (c_associated(coef%G1_d)) then
       call device_free(coef%G1_d)
    end if

    if (c_associated(coef%G2_d)) then
       call device_free(coef%G2_d)
    end if

    if (c_associated(coef%G3_d)) then
       call device_free(coef%G3_d)
    end if

    if (c_associated(coef%G4_d)) then
       call device_free(coef%G4_d)
    end if

    if (c_associated(coef%G5_d)) then
       call device_free(coef%G5_d)
    end if

    if (c_associated(coef%G6_d)) then
       call device_free(coef%G6_d)
    end if

    if (c_associated(coef%dxdr_d)) then
       call device_Free(coef%dxdr_d)
    end if

    if (c_associated(coef%dydr_d)) then
       call device_Free(coef%dydr_d)
    end if

    if (c_associated(coef%dzdr_d)) then
       call device_Free(coef%dzdr_d)
    end if

    if (c_associated(coef%dxds_d)) then
       call device_Free(coef%dxds_d)
    end if

    if (c_associated(coef%dyds_d)) then
       call device_free(coef%dyds_d)
    end if

    if (c_associated(coef%dzds_d)) then
       call device_free(coef%dzds_d)
    end if
    
    if (c_associated(coef%dxdt_d)) then
       call device_free(coef%dxdt_d)
    end if

    if (c_associated(coef%dydt_d)) then
       call device_free(coef%dydt_d)
    end if

    if (c_associated(coef%dzdt_d)) then
       call device_free(coef%dzdt_d)
    end if

    if (c_associated(coef%drdx_d)) then
       call device_Free(coef%drdx_d)
    end if

    if (c_associated(coef%drdy_d)) then
       call device_Free(coef%drdy_d)
    end if

    if (c_associated(coef%drdz_d)) then
       call device_Free(coef%drdz_d)
    end if

    if (c_associated(coef%dsdx_d)) then
       call device_Free(coef%dsdx_d)
    end if

    if (c_associated(coef%dsdy_d)) then
       call device_free(coef%dsdy_d)
    end if

    if (c_associated(coef%dsdz_d)) then
       call device_free(coef%dsdz_d)
    end if
    
    if (c_associated(coef%dtdx_d)) then
       call device_free(coef%dtdx_d)
    end if

    if (c_associated(coef%dtdy_d)) then
       call device_free(coef%dtdy_d)
    end if

    if (c_associated(coef%dtdz_d)) then
       call device_free(coef%dtdz_d)
    end if
    
    if (c_associated(coef%mult_d)) then
       call device_free(coef%mult_d)
    end if

    if (c_associated(coef%h1_d)) then
       call device_free(coef%h1_d)
    end if
    
    if (c_associated(coef%h2_d)) then
       call device_free(coef%h2_d)
    end if

    if (c_associated(coef%jac_d)) then
       call device_free(coef%jac_d)
    end if

    if (c_associated(coef%jacinv_d)) then
       call device_free(coef%jacinv_d)
    end if
    
    if (c_associated(coef%B_d)) then
       call device_free(coef%B_d)
    end if
    
    if (c_associated(coef%Binv_d)) then
       call device_free(coef%Binv_d)
    end if

    if (c_associated(coef%area_d)) then
       call device_free(coef%area_d)
    end if
    
    if (c_associated(coef%nx_d)) then
       call device_free(coef%nx_d)
    end if

    if (c_associated(coef%ny_d)) then
       call device_free(coef%ny_d)
    end if

    if (c_associated(coef%nz_d)) then
       call device_Free(coef%nz_d)
    end if

  end subroutine coef_free

  subroutine coef_generate_dxyzdrst(c)
    type(coef_t), intent(inout) :: c
    integer :: e,i,lxy,lyz
    
    lxy=c%Xh%lx*c%Xh%ly
    lyz=c%Xh%ly*c%Xh%lz
    do e =1,c%msh%nelv
       call mxm(c%Xh%dx,c%Xh%lx,c%dof%x(1,1,1,e),c%Xh%lx,c%dxdr(1,1,1,e),lyz)
       call mxm(c%Xh%dx,c%Xh%lx,c%dof%y(1,1,1,e),c%Xh%lx,c%dydr(1,1,1,e),lyz)
       call mxm(c%Xh%dx,c%Xh%lx,c%dof%z(1,1,1,e),c%Xh%lx,c%dzdr(1,1,1,e),lyz)

       DO i=1,c%Xh%lz
          call mxm(c%dof%x(1,1,i,e),c%Xh%lx,c%Xh%dyt,c%Xh%ly,c%dxds(1,1,i,e),c%Xh%ly)
          call mxm(c%dof%y(1,1,i,e),c%Xh%lx,c%Xh%dyt,c%Xh%ly,c%dyds(1,1,i,e),c%Xh%ly)
          call mxm(c%dof%z(1,1,i,e),c%Xh%lx,c%Xh%dyt,c%Xh%ly,c%dzds(1,1,i,e),c%Xh%ly)
       end do
       !> We actually take 2d into account, wow, need to do that for the rest.
       if(c%msh%gdim .eq. 3) then
          call mxm(c%dof%x(1,1,1,e),lxy,c%Xh%dzt,c%Xh%lz,c%dxdt(1,1,1,e),c%Xh%lz)
          call mxm(c%dof%y(1,1,1,e),lxy,c%Xh%dzt,c%Xh%lz,c%dydt(1,1,1,e),c%Xh%lz)
          call mxm(c%dof%z(1,1,1,e),lxy,c%Xh%dzt,c%Xh%lz,c%dzdt(1,1,1,e),c%Xh%lz)
       else
          call rzero(c%dxdt(1,1,1,e),lxy)
          call rzero(c%dydt(1,1,1,e),lxy)
          call rone(c%dzdt(1,1,1,e),lxy)
       endif
    end do

    if (c%msh%gdim .eq. 2) then
       call rzero   (c%jac,c%dof%n_dofs)
       call addcol3 (c%jac,c%dxdr,c%dyds,c%dof%n_dofs)
       call subcol3 (c%jac,c%dxds,c%dydr,c%dof%n_dofs)
       call copy    (c%drdx,c%dyds,c%dof%n_dofs)
       call copy    (c%drdy,c%dxds,c%dof%n_dofs)
       call chsign  (c%drdy,c%dof%n_dofs)
       call copy    (c%dsdx,c%dydr,c%dof%n_dofs)
       call chsign  (c%dsdx,c%dof%n_dofs)
       call copy    (c%dsdy,c%dxdr,c%dof%n_dofs)
       call rzero   (c%drdz,c%dof%n_dofs)
       call rzero   (c%dsdz,c%dof%n_dofs)
       call RONE    (c%dtdz,c%dof%n_dofs)
   else
       call rzero   (c%jac,c%dof%n_dofs)
       call addcol4 (c%jac,c%dxdr,c%dyds,c%dzdt,c%dof%n_dofs)
       call addcol4 (c%jac,c%dxdt,c%dydr,c%dzds,c%dof%n_dofs)
       call addcol4 (c%jac,c%dxds,c%dydt,c%dzdr,c%dof%n_dofs)
       call subcol4 (c%jac,c%dxdr,c%dydt,c%dzds,c%dof%n_dofs)
       call subcol4 (c%jac,c%dxds,c%dydr,c%dzdt,c%dof%n_dofs)
       call subcol4 (c%jac,c%dxdt,c%dyds,c%dzdr,c%dof%n_dofs)
       call ascol5  (c%drdx,c%dyds,c%dzdt,c%dydt,c%dzds,c%dof%n_dofs)
       call ascol5  (c%drdy,c%dxdt,c%dzds,c%dxds,c%dzdt,c%dof%n_dofs)
       call ascol5  (c%drdz,c%dxds,c%dydt,c%dxdt,c%dyds,c%dof%n_dofs)
       call ascol5  (c%dsdx,c%dydt,c%dzdr,c%dydr,c%dzdt,c%dof%n_dofs)
       call ascol5  (c%dsdy,c%dxdr,c%dzdt,c%dxdt,c%dzdr,c%dof%n_dofs)
       call ascol5  (c%dsdz,c%dxdt,c%dydr,c%dxdr,c%dydt,c%dof%n_dofs)
       call ascol5  (c%dtdx,c%dydr,c%dzds,c%dyds,c%dzdr,c%dof%n_dofs)
       call ascol5  (c%dtdy,c%dxds,c%dzdr,c%dxdr,c%dzds,c%dof%n_dofs)
       call ascol5  (c%dtdz,c%dxdr,c%dyds,c%dxds,c%dydr,c%dof%n_dofs)
    end if

    call invers2(c%jacinv,c%jac,c%dof%n_dofs)

    !>  @todo cleanup once we have device math in place
    if ((NEKO_BCKND_HIP .eq. 1) .or. (NEKO_BCKND_CUDA .eq. 1)) then
       call device_memcpy(c%dxds, c%dxds_d, c%dof%n_dofs, HOST_TO_DEVICE)
       call device_memcpy(c%dydr, c%dyds_d, c%dof%n_dofs, HOST_TO_DEVICE)
       call device_memcpy(c%dzdr, c%dzds_d, c%dof%n_dofs, HOST_TO_DEVICE)
       call device_memcpy(c%dxds, c%dxds_d, c%dof%n_dofs, HOST_TO_DEVICE)
       call device_memcpy(c%dyds, c%dyds_d, c%dof%n_dofs, HOST_TO_DEVICE)
       call device_memcpy(c%dzds, c%dzds_d, c%dof%n_dofs, HOST_TO_DEVICE)
       call device_memcpy(c%dxdt, c%dxdt_d, c%dof%n_dofs, HOST_TO_DEVICE)
       call device_memcpy(c%dydt, c%dydt_d, c%dof%n_dofs, HOST_TO_DEVICE)
       call device_memcpy(c%dzdt, c%dzdt_d, c%dof%n_dofs, HOST_TO_DEVICE)       
       call device_memcpy(c%drdx, c%drdx_d, c%dof%n_dofs, HOST_TO_DEVICE)
       call device_memcpy(c%drdy, c%drdy_d, c%dof%n_dofs, HOST_TO_DEVICE)
       call device_memcpy(c%drdz, c%drdz_d, c%dof%n_dofs, HOST_TO_DEVICE)
       call device_memcpy(c%dsdx, c%dsdx_d, c%dof%n_dofs, HOST_TO_DEVICE)
       call device_memcpy(c%dsdy, c%dsdy_d, c%dof%n_dofs, HOST_TO_DEVICE)
       call device_memcpy(c%dsdz, c%dsdz_d, c%dof%n_dofs, HOST_TO_DEVICE)
       call device_memcpy(c%dtdx, c%dtdx_d, c%dof%n_dofs, HOST_TO_DEVICE)
       call device_memcpy(c%dtdy, c%dtdy_d, c%dof%n_dofs, HOST_TO_DEVICE)
       call device_memcpy(c%dtdz, c%dtdz_d, c%dof%n_dofs, HOST_TO_DEVICE)       
       call device_memcpy(c%jac, c%jac_d, c%dof%n_dofs, HOST_TO_DEVICE)
       call device_memcpy(c%jacinv, c%jacinv_d, c%dof%n_dofs, HOST_TO_DEVICE)
    end if

  end subroutine coef_generate_dxyzdrst

  !> Generate geometric data for the given mesh
  !! @note Current implementation assumes regular shaped hex elements
  subroutine coef_generate_geo(c)
    type(coef_t), intent(inout) :: c
    integer :: e, lxyz

    lxyz = c%Xh%lx * c%Xh%ly * c%Xh%lz
    
    if(c%msh%gdim .eq. 2) then
      call vdot2(c%G1,c%drdx,c%drdy,c%drdx,c%drdy,c%dof%n_dofs)
      call vdot2(c%G2,c%dsdx,c%dsdy,c%dsdx,c%dsdy,c%dof%n_dofs)
      call vdot2(c%G4,c%drdx,c%drdy,c%dsdx,c%dsdy,c%dof%n_dofs)
      call  col2(c%G1,c%jacinv,c%dof%n_dofs)
      call  col2(c%G2,c%jacinv,c%dof%n_dofs)
      call  col2(c%G4,c%jacinv,c%dof%n_dofs)
      call rzero(c%G3,c%dof%n_dofs)
      call rzero(c%G5,c%dof%n_dofs)
      call rzero(c%G6,c%dof%n_dofs)
    else
      call vdot3(c%G1,c%drdx,c%drdy,c%drdz,c%drdx,c%drdy,c%drdz,c%dof%n_dofs)
      call vdot3(c%G2,c%dsdx,c%dsdy,c%dsdz,c%dsdx,c%dsdy,c%dsdz,c%dof%n_dofs)
      call vdot3(c%G3,c%dtdx,c%dtdy,c%dtdz,c%dtdx,c%dtdy,c%dtdz,c%dof%n_dofs)
      call vdot3(c%G4,c%drdx,c%drdy,c%drdz,c%dsdx,c%dsdy,c%dsdz,c%dof%n_dofs)
      call vdot3(c%G5,c%drdx,c%drdy,c%drdz,c%dtdx,c%dtdy,c%dtdz,c%dof%n_dofs)
      call vdot3(c%G6,c%dsdx,c%dsdy,c%dsdz,c%dtdx,c%dtdy,c%dtdz,c%dof%n_dofs)
      
      call col2(c%G1,c%jacinv,c%dof%n_dofs)
      call col2(c%G2,c%jacinv,c%dof%n_dofs)
      call col2(c%G3,c%jacinv,c%dof%n_dofs)
      call col2(c%G4,c%jacinv,c%dof%n_dofs)
      call col2(c%G5,c%jacinv,c%dof%n_dofs)
      call col2(c%G6,c%jacinv,c%dof%n_dofs)
    end if
    do e=1,c%msh%nelv
       call col2(c%G1(1,1,1,e),c%Xh%w3,lxyz)
       call col2(c%G2(1,1,1,e),c%Xh%w3,lxyz)
       call col2(c%G4(1,1,1,e),c%Xh%w3,lxyz)
       if (c%msh%gdim .eq. 3) then
         call col2(c%G3(1,1,1,e),c%Xh%w3,lxyz)
         call col2(c%G5(1,1,1,e),c%Xh%w3,lxyz)
         call col2(c%G6(1,1,1,e),c%Xh%w3,lxyz)
       end if
    end do

    !>  @todo cleanup once we have device math in place
    if ((NEKO_BCKND_HIP .eq. 1) .or. (NEKO_BCKND_CUDA .eq. 1)) then
       call device_memcpy(c%G1, c%G1_d, c%dof%n_dofs, HOST_TO_DEVICE)
       call device_memcpy(c%G2, c%G2_d, c%dof%n_dofs, HOST_TO_DEVICE)
       call device_memcpy(c%G3, c%G3_d, c%dof%n_dofs, HOST_TO_DEVICE)
       call device_memcpy(c%G4, c%G4_d, c%dof%n_dofs, HOST_TO_DEVICE)
       call device_memcpy(c%G5, c%G5_d, c%dof%n_dofs, HOST_TO_DEVICE)
       call device_memcpy(c%G6, c%G6_d, c%dof%n_dofs, HOST_TO_DEVICE)
    end if
  end subroutine coef_generate_geo
 
  !> Generate mass matrix B for the given mesh and space
  !! @note This is also a stapleholder, we need to go through the coef class properly.
  subroutine coef_generate_mass(c)
    type(coef_t), intent(inout) :: c
    integer :: e, j, k, l, lxyz
    
    lxyz = c%Xh%lx * c%Xh%ly * c%Xh%lz
    
    call rone(c%B,c%dof%n_dofs)
    do e=1,c%msh%nelv
       ! Here we need to handle things differently for axis symmetric elements
       call col3(c%B(1,1,1,e),c%jac(1,1,1,e),c%Xh%w3,lxyz)
    end do
    
    call copy(c%Binv,c%B,c%dof%n_dofs)

    !>  @todo cleanup once we have device math in place
    if ((NEKO_BCKND_HIP .eq. 1) .or. (NEKO_BCKND_CUDA .eq. 1)) then
       call device_memcpy(c%Binv, c%Binv_d, c%dof%n_dofs, HOST_TO_DEVICE)
    end if
    
    call gs_op_vector(c%gs_h,c%Binv, c%dof%n_dofs,GS_OP_ADD)

    !>  @todo cleanup once we have device math in place
    if ((NEKO_BCKND_HIP .eq. 1) .or. (NEKO_BCKND_CUDA .eq. 1)) then
       call device_memcpy(c%Binv, c%Binv_d, c%dof%n_dofs, DEVICE_TO_HOST)
    end if

    call invcol1(c%Binv,c%dof%n_dofs)

    !>  @todo cleanup once we have device math in place
    if ((NEKO_BCKND_HIP .eq. 1) .or. (NEKO_BCKND_CUDA .eq. 1)) then
       call device_memcpy(c%B, c%B_d, c%dof%n_dofs, HOST_TO_DEVICE)
       call device_memcpy(c%Binv, c%Binv_d, c%dof%n_dofs, HOST_TO_DEVICE)
    end if

    c%volume = glsum(c%B,c%dof%n_dofs)
  end subroutine coef_generate_mass

  !> Generate facet area and surface normals
  subroutine coef_generate_area_and_normal(coef)
    type(coef_t), intent(inout) :: coef
    real(kind=rp), allocatable :: a(:,:,:,:)
    real(kind=rp), allocatable :: b(:,:,:,:)
    real(kind=rp), allocatable :: c(:,:,:,:)
    real(kind=rp), allocatable :: dot(:,:,:,:)
    integer :: n, e, j, k, l, lx
    real(kind=rp) :: weight, len
    n = coef%dof%n_dofs
    lx = coef%Xh%lx
    
    allocate(a(coef%Xh%lx, coef%Xh%lx, coef%Xh%lx, coef%msh%nelv))
    allocate(b(coef%Xh%lx, coef%Xh%lx, coef%Xh%lx, coef%msh%nelv))
    allocate(c(coef%Xh%lx, coef%Xh%lx, coef%Xh%lx, coef%msh%nelv))
    allocate(dot(coef%Xh%lx, coef%Xh%lx, coef%Xh%lx, coef%msh%nelv))

    call vcross(a,b,c, coef%dxds, coef%dyds, coef%dzds, &
         coef%dxdt, coef%dydt, coef%dzdt, n)
    call vdot3(dot, a, b, c, a, b, c, n)

    do e = 1, coef%msh%nelv
       do k = 1, coef%Xh%lx
          do j = 1, coef%Xh%lx
             weight = coef%Xh%wy(j) * coef%Xh%wz(k)
             coef%area(j, k, 2, e) = sqrt(dot(lx, j, k, e)) * weight
             coef%area(j, k, 1, e) = sqrt(dot(1, j, k, e)) * weight
             coef%nx(j,k, 1, e) = -A(1, j, k, e)
             coef%nx(j,k, 2, e) =  A(lx, j, k, e)
             coef%ny(j,k, 1, e) = -B(1, j, k, e)
             coef%ny(j,k, 2, e) =  B(lx, j, k, e)
             coef%nz(j,k, 1, e) = -C(1, j, k, e)
             coef%nz(j,k, 2, e) =  C(lx, j, k, e)
          end do
       end do
    end do

    call vcross(a,b,c, coef%dxdr, coef%dydr, coef%dzdr, &
         coef%dxdt, coef%dydt, coef%dzdt, n)
    call vdot3(dot, a, b, c, a, b, c, n)
    
    do e = 1, coef%msh%nelv
       do k = 1, coef%Xh%lx
          do j = 1, coef%Xh%lx
             weight = coef%Xh%wx(j) * coef%Xh%wz(k)
             coef%area(j, k, 3, e) = sqrt(dot(j, 1, k, e)) * weight
             coef%area(j, k, 4, e) = sqrt(dot(j, lx, k, e)) * weight
             coef%nx(j,k, 3, e) =  A(j, 1, k, e)
             coef%nx(j,k, 4, e) = -A(j, lx, k, e)
             coef%ny(j,k, 3, e) =  B(j, 1, k, e)
             coef%ny(j,k, 4, e) = -B(j, lx, k, e)
             coef%nz(j,k, 3, e) =  C(j, 1, k, e)
             coef%nz(j,k, 4, e) = -C(j, lx, k, e)             
          end do
       end do
    end do


    call vcross(a,b,c, coef%dxdr, coef%dydr, coef%dzdr, &
         coef%dxds, coef%dyds, coef%dzds, n)
    call vdot3(dot, a, b, c, a, b, c, n)
    
    do e = 1, coef%msh%nelv
       do k = 1, coef%Xh%lx
          do j = 1, coef%Xh%lx
             weight = coef%Xh%wx(j) * coef%Xh%wy(k)
             coef%area(j, k, 5, e) = sqrt(dot(j, k, 1, e)) * weight
             coef%area(j, k, 6, e) = sqrt(dot(j, j, lx, e)) * weight
             coef%nx(j,k, 5, e) = -A(j, k, 1, e)
             coef%nx(j,k, 6, e) =  A(j, k, lx, e)
             coef%ny(j,k, 5, e) = -B(j, k, 1, e)
             coef%ny(j,k, 6, e) =  B(j, k, lx, e)
             coef%nz(j,k, 5, e) = -C(j, k, 1, e)
             coef%nz(j,k, 6, e) =  C(j, k, lx, e)             
          end do
       end do
    end do

    ! Normalize
    n = size(coef%nz)
    do j = 1, n
       len = sqrt(coef%nx(j,1,1,1)**2 + &
            coef%ny(j,1,1,1)**2 + coef%nz(j,1,1,1)**2)
       if (len .gt. NEKO_EPS) then
          coef%nx(j,1,1,1) = coef%nx(j,1,1,1) / len
          coef%ny(j,1,1,1) = coef%ny(j,1,1,1) / len
          coef%nz(j,1,1,1) = coef%nz(j,1,1,1) / len
       end if
    end do

    deallocate(dot)
    deallocate(c)
    deallocate(b)
    deallocate(a)
    !>  @todo cleanup once we have device math in place
    if ((NEKO_BCKND_HIP .eq. 1) .or. (NEKO_BCKND_CUDA .eq. 1)) then
       n = size(coef%area)
       call device_memcpy(coef%area, coef%area_d, n, HOST_TO_DEVICE)
       call device_memcpy(coef%nx, coef%nx_d, n, HOST_TO_DEVICE)
       call device_memcpy(coef%ny, coef%ny_d, n, HOST_TO_DEVICE)
       call device_memcpy(coef%nz, coef%nz_d, n, HOST_TO_DEVICE)
    end if
    
  end subroutine coef_generate_area_and_normal
  
end module coefs
