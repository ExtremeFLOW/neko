! Copyright (c) 2020-2023, The Neko Authors
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
!> Coefficients 
module coefs
  use gather_scatter
  use neko_config
  use num_types
  use space  
  use math
  use mesh
  use device
  use device_coef
  use mxm_wrapper
  use, intrinsic :: iso_c_binding
  implicit none
  private
  
  !> Coefficients defined on a given (mesh, \f$ X_h \f$) tuple.
  !! Arrays use indices (i,j,k,e): element e, local coordinate (i,j,k).
  type, public :: coef_t
     !> Geometric factors \f$ G_{11} \f$
     real(kind=rp), allocatable :: G11(:,:,:,:)
     !> Geometric factors \f$ G_{22} \f$
     real(kind=rp), allocatable :: G22(:,:,:,:)
     !> Geometric factors \f$ G_{33} \f$
     real(kind=rp), allocatable :: G33(:,:,:,:)
     !> Geometric factors \f$ G_{12} \f$
     real(kind=rp), allocatable :: G12(:,:,:,:)
     !> Geometric factors \f$ G_{13} \f$
     real(kind=rp), allocatable :: G13(:,:,:,:)
     !> Geometric factors \f$ G_{23} \f$
     real(kind=rp), allocatable :: G23(:,:,:,:)

     real(kind=rp), allocatable :: mult(:,:,:,:) !< Multiplicity
     !> generate mapping data between element and reference element 
     !! \f$ dx/dr, dy/dr, dz/dr \f$
     !! \f$ dx/ds, dy/ds, dz/ds \f$
     !! \f$ dx/dt, dy/dt, dz/dt \f$
     real(kind=rp), allocatable :: dxdr(:,:,:,:), dydr(:,:,:,:), dzdr(:,:,:,:) 
     real(kind=rp), allocatable :: dxds(:,:,:,:), dyds(:,:,:,:), dzds(:,:,:,:)
     real(kind=rp), allocatable :: dxdt(:,:,:,:), dydt(:,:,:,:), dzdt(:,:,:,:) 
     !> \f$ dr/dx, dr/dy, dr/dz \f$
     !! \f$ ds/dx, ds/dy, ds/dz \f$
     !! \f$ dt/dx, dt/dy, dt/dz \f$
     real(kind=rp), allocatable :: drdx(:,:,:,:), drdy(:,:,:,:), drdz(:,:,:,:) 
     real(kind=rp), allocatable :: dsdx(:,:,:,:), dsdy(:,:,:,:), dsdz(:,:,:,:)
     real(kind=rp), allocatable :: dtdx(:,:,:,:), dtdy(:,:,:,:), dtdz(:,:,:,:) 
     
     real(kind=rp), allocatable :: h1(:,:,:,:) !< Stiffness scaling
     real(kind=rp), allocatable :: h2(:,:,:,:) !< Mass scaling
     logical :: ifh2 !< True if h2 .ne. 0
     
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
     
     type(c_ptr) :: G11_d = C_NULL_PTR
     type(c_ptr) :: G22_d = C_NULL_PTR
     type(c_ptr) :: G33_d = C_NULL_PTR
     type(c_ptr) :: G12_d = C_NULL_PTR
     type(c_ptr) :: G13_d = C_NULL_PTR
     type(c_ptr) :: G23_d = C_NULL_PTR
     type(c_ptr) :: dxdr_d = C_NULL_PTR
     type(c_ptr) :: dydr_d = C_NULL_PTR
     type(c_ptr) :: dzdr_d = C_NULL_PTR
     type(c_ptr) :: dxds_d = C_NULL_PTR
     type(c_ptr) :: dyds_d = C_NULL_PTR
     type(c_ptr) :: dzds_d = C_NULL_PTR
     type(c_ptr) :: dxdt_d = C_NULL_PTR
     type(c_ptr) :: dydt_d = C_NULL_PTR
     type(c_ptr) :: dzdt_d = C_NULL_PTR
     type(c_ptr) :: drdx_d = C_NULL_PTR
     type(c_ptr) :: drdy_d = C_NULL_PTR
     type(c_ptr) :: drdz_d = C_NULL_PTR
     type(c_ptr) :: dsdx_d = C_NULL_PTR
     type(c_ptr) :: dsdy_d = C_NULL_PTR
     type(c_ptr) :: dsdz_d = C_NULL_PTR
     type(c_ptr) :: dtdx_d = C_NULL_PTR
     type(c_ptr) :: dtdy_d = C_NULL_PTR
     type(c_ptr) :: dtdz_d = C_NULL_PTR
     type(c_ptr) :: mult_d = C_NULL_PTR
     type(c_ptr) :: h1_d = C_NULL_PTR
     type(c_ptr) :: h2_d = C_NULL_PTR
     type(c_ptr) :: jac_d = C_NULL_PTR
     type(c_ptr) :: jacinv_d = C_NULL_PTR
     type(c_ptr) :: B_d = C_NULL_PTR
     type(c_ptr) :: Binv_d = C_NULL_PTR
     type(c_ptr) :: area_d = C_NULL_PTR
     type(c_ptr) :: nx_d = C_NULL_PTR
     type(c_ptr) :: ny_d = C_NULL_PTR
     type(c_ptr) :: nz_d = C_NULL_PTR

  end type coef_t

  interface coef_init
     module procedure coef_init_empty, coef_init_all
  end interface coef_init
  
  public :: coef_init, coef_free, coef_get_normal
  
contains
  
  !> Initialize empty coefs for a space and a mesh 
  subroutine coef_init_empty(coef, Xh, msh)
    type(coef_t), intent(inout) :: coef
    type(space_t), intent(inout), target :: Xh
    type(mesh_t), intent(inout), target :: msh
    integer :: n    
    call coef_free(coef)
    coef%msh => msh
    coef%Xh => Xh
    
    allocate(coef%drdx(coef%Xh%lx, coef%Xh%ly, coef%Xh%lz, coef%msh%nelv))
    allocate(coef%dsdx(coef%Xh%lx, coef%Xh%ly, coef%Xh%lz, coef%msh%nelv))
    allocate(coef%dtdx(coef%Xh%lx, coef%Xh%ly, coef%Xh%lz, coef%msh%nelv))
    
    allocate(coef%drdy(coef%Xh%lx, coef%Xh%ly, coef%Xh%lz, coef%msh%nelv))
    allocate(coef%dsdy(coef%Xh%lx, coef%Xh%ly, coef%Xh%lz, coef%msh%nelv))
    allocate(coef%dtdy(coef%Xh%lx, coef%Xh%ly, coef%Xh%lz, coef%msh%nelv))
    
    allocate(coef%drdz(coef%Xh%lx, coef%Xh%ly, coef%Xh%lz, coef%msh%nelv))
    allocate(coef%dsdz(coef%Xh%lx, coef%Xh%ly, coef%Xh%lz, coef%msh%nelv))
    allocate(coef%dtdz(coef%Xh%lx, coef%Xh%ly, coef%Xh%lz, coef%msh%nelv))
    
    !
    ! Setup device memory (if present)
    !
    
    n = coef%Xh%lx * coef%Xh%ly * coef%Xh%lz * coef%msh%nelv
    if ((NEKO_BCKND_HIP .eq. 1) .or. (NEKO_BCKND_CUDA .eq. 1) .or. &
        (NEKO_BCKND_OPENCL .eq. 1)) then

       call device_map(coef%drdx, coef%drdx_d, n)
       call device_map(coef%drdy, coef%drdy_d, n)
       call device_map(coef%drdz, coef%drdz_d, n)

       call device_map(coef%dsdx, coef%dsdx_d, n)
       call device_map(coef%dsdy, coef%dsdy_d, n)
       call device_map(coef%dsdz, coef%dsdz_d, n)

       call device_map(coef%dtdx, coef%dtdx_d, n)
       call device_map(coef%dtdy, coef%dtdy_d, n)
       call device_map(coef%dtdz, coef%dtdz_d, n)
       
    end if
    
  end subroutine coef_init_empty

  !> Initialize coefficients
  subroutine coef_init_all(coef, gs_h)
    type(coef_t), intent(inout) :: coef
    type(gs_t), intent(inout), target :: gs_h
    integer :: n, m
    call coef_free(coef)
    
    coef%msh => gs_h%dofmap%msh
    coef%Xh => gs_h%dofmap%Xh
    coef%dof => gs_h%dofmap
    coef%gs_h => gs_h

    !
    ! Allocate arrays for geometric data
    !
    !>@todo Be clever and try to avoid allocating zeroed geom. factors
    allocate(coef%G11(coef%Xh%lx, coef%Xh%ly, coef%Xh%lz, coef%msh%nelv))
    allocate(coef%G22(coef%Xh%lx, coef%Xh%ly, coef%Xh%lz, coef%msh%nelv))
    allocate(coef%G33(coef%Xh%lx, coef%Xh%ly, coef%Xh%lz, coef%msh%nelv))
    allocate(coef%G12(coef%Xh%lx, coef%Xh%ly, coef%Xh%lz, coef%msh%nelv))
    allocate(coef%G13(coef%Xh%lx, coef%Xh%ly, coef%Xh%lz, coef%msh%nelv))
    allocate(coef%G23(coef%Xh%lx, coef%Xh%ly, coef%Xh%lz, coef%msh%nelv))
    
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
    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_map(coef%G11, coef%G11_d, n)
       call device_map(coef%G22, coef%G22_d, n)
       call device_map(coef%G33, coef%G33_d, n)
       call device_map(coef%G12, coef%G12_d, n)
       call device_map(coef%G13, coef%G13_d, n)
       call device_map(coef%G23, coef%G23_d, n)
       
       call device_map(coef%dxdr, coef%dxdr_d, n)
       call device_map(coef%dydr, coef%dydr_d, n)
       call device_map(coef%dzdr, coef%dzdr_d, n)

       call device_map(coef%dxds, coef%dxds_d, n)
       call device_map(coef%dyds, coef%dyds_d, n)
       call device_map(coef%dzds, coef%dzds_d, n)
       
       call device_map(coef%dxdt, coef%dxdt_d, n)
       call device_map(coef%dydt, coef%dydt_d, n)
       call device_map(coef%dzdt, coef%dzdt_d, n)

       call device_map(coef%drdx, coef%drdx_d, n)
       call device_map(coef%drdy, coef%drdy_d, n)
       call device_map(coef%drdz, coef%drdz_d, n)

       call device_map(coef%dsdx, coef%dsdx_d, n)
       call device_map(coef%dsdy, coef%dsdy_d, n)
       call device_map(coef%dsdz, coef%dsdz_d, n)

       call device_map(coef%dtdx, coef%dtdx_d, n)
       call device_map(coef%dtdy, coef%dtdy_d, n)
       call device_map(coef%dtdz, coef%dtdz_d, n)

       call device_map(coef%mult, coef%mult_d, n)
       call device_map(coef%h1, coef%h1_d, n)
       call device_map(coef%h2, coef%h2_d, n)

       call device_map(coef%jac, coef%jac_d, n)
       call device_map(coef%jacinv, coef%jacinv_d, n)
       call device_map(coef%B, coef%B_d, n)
       call device_map(coef%Binv, coef%Binv_d, n)

       m = coef%Xh%lx * coef%Xh%ly * 6 * coef%msh%nelv
       
       call device_map(coef%area, coef%area_d, m)
       call device_map(coef%nx, coef%nx_d, m)
       call device_map(coef%ny, coef%ny_d, m)
       call device_map(coef%nz, coef%nz_d, m)
    end if

    call coef_generate_dxyzdrst(coef)
    
    call coef_generate_geo(coef)

    call coef_generate_area_and_normal(coef)

    call coef_generate_mass(coef)
    
    ! This is a placeholder, just for now
    ! We can probably find a prettier solution
    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_rone(coef%h1_d, n)
       call device_rone(coef%h2_d, n)
       call device_memcpy(coef%h1, coef%h1_d, n, DEVICE_TO_HOST)
       call device_memcpy(coef%h2, coef%h2_d, n, DEVICE_TO_HOST)
    else
       call rone(coef%h1,n)
       call rone(coef%h2,n)
    end if

    coef%ifh2 = .false.
    
    !
    ! Set up multiplicity
    !
    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_rone(coef%mult_d, n)
    else
       call rone(coef%mult, n)
    end if
       
    call gs_op(gs_h, coef%mult, n, GS_OP_ADD)

    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_invcol1(coef%mult_d, n)
       call device_memcpy(coef%mult, coef%mult_d, n, DEVICE_TO_HOST)
    else    
       call invcol1(coef%mult, n)
    end if

  end subroutine coef_init_all

  !> Deallocate coefficients
  subroutine coef_free(coef)
    type(coef_t), intent(inout) :: coef

    if (allocated(coef%G11)) then
       deallocate(coef%G11)
    end if

    if (allocated(coef%G22)) then
       deallocate(coef%G22)
    end if

    if (allocated(coef%G33)) then
       deallocate(coef%G33)
    end if

    if (allocated(coef%G12)) then
       deallocate(coef%G12)
    end if

    if (allocated(coef%G13)) then
       deallocate(coef%G13)
    end if
    
    if (allocated(coef%G23)) then
       deallocate(coef%G23)
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
    
    if (c_associated(coef%G11_d)) then
       call device_free(coef%G11_d)
    end if

    if (c_associated(coef%G22_d)) then
       call device_free(coef%G22_d)
    end if

    if (c_associated(coef%G33_d)) then
       call device_free(coef%G33_d)
    end if

    if (c_associated(coef%G12_d)) then
       call device_free(coef%G12_d)
    end if

    if (c_associated(coef%G13_d)) then
       call device_free(coef%G13_d)
    end if

    if (c_associated(coef%G23_d)) then
       call device_free(coef%G23_d)
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
    integer :: e,i,lxy,lyz,ntot
    
    lxy=c%Xh%lx*c%Xh%ly
    lyz=c%Xh%ly*c%Xh%lz
    ntot = c%dof%size()
       
    associate(drdx => c%drdx, drdy => c%drdy, drdz => c%drdz, &
         dsdx => c%dsdx, dsdy => c%dsdy, dsdz => c%dsdz, &
         dtdx => c%dtdx, dtdy => c%dtdy, dtdz => c%dtdz, &
         dxdr => c%dxdr, dydr => c%dydr, dzdr => c%dzdr, &
         dxds => c%dxds, dyds => c%dyds, dzds => c%dzds, &
         dxdt => c%dxdt, dydt => c%dydt, dzdt => c%dzdt, &
         dx => c%Xh%dx, dy => c%Xh%dy, dz => c%Xh%dz, &
         x => c%dof%x, y => c%dof%y, z => c%dof%z, &
         lx => c%Xh%lx, ly => c%Xh%ly, lz => c%Xh%lz, &
         dyt => c%Xh%dyt, dzt => c%Xh%dzt, &
         jacinv => c%jacinv, jac => c%jac)

        if (NEKO_BCKND_DEVICE .eq. 1) then

         call device_coef_generate_dxydrst(c%drdx_d, c%drdy_d, c%drdz_d, &
              c%dsdx_d, c%dsdy_d, c%dsdz_d, c%dtdx_d, c%dtdy_d, c%dtdz_d, &
              c%dxdr_d, c%dydr_d, c%dzdr_d, c%dxds_d, c%dyds_d, c%dzds_d, &
              c%dxdt_d, c%dydt_d, c%dzdt_d, c%Xh%dx_d, c%Xh%dy_d, c%Xh%dz_d, &
              c%dof%x_d, c%dof%y_d, c%dof%z_d, c%jacinv_d, c%jac_d, &
              c%Xh%lx, c%msh%nelv)

         call device_memcpy(dxdr, c%dxdr_d, ntot, DEVICE_TO_HOST)
         call device_memcpy(dydr, c%dydr_d, ntot, DEVICE_TO_HOST)
         call device_memcpy(dzdr, c%dzdr_d, ntot, DEVICE_TO_HOST)
         call device_memcpy(dxds, c%dxds_d, ntot, DEVICE_TO_HOST)
         call device_memcpy(dyds, c%dyds_d, ntot, DEVICE_TO_HOST)
         call device_memcpy(dzds, c%dzds_d, ntot, DEVICE_TO_HOST)
         call device_memcpy(dxdt, c%dxdt_d, ntot, DEVICE_TO_HOST)
         call device_memcpy(dydt, c%dydt_d, ntot, DEVICE_TO_HOST)
         call device_memcpy(dzdt, c%dzdt_d, ntot, DEVICE_TO_HOST)
         call device_memcpy(drdx, c%drdx_d, ntot, DEVICE_TO_HOST)
         call device_memcpy(drdy, c%drdy_d, ntot, DEVICE_TO_HOST)
         call device_memcpy(drdz, c%drdz_d, ntot, DEVICE_TO_HOST)
         call device_memcpy(dsdx, c%dsdx_d, ntot, DEVICE_TO_HOST)
         call device_memcpy(dsdy, c%dsdy_d, ntot, DEVICE_TO_HOST)
         call device_memcpy(dsdz, c%dsdz_d, ntot, DEVICE_TO_HOST)
         call device_memcpy(dtdx, c%dtdx_d, ntot, DEVICE_TO_HOST)
         call device_memcpy(dtdy, c%dtdy_d, ntot, DEVICE_TO_HOST)
         call device_memcpy(dtdz, c%dtdz_d, ntot, DEVICE_TO_HOST)
         call device_memcpy(jac, c%jac_d, ntot, DEVICE_TO_HOST)
         call device_memcpy(jacinv, c%jacinv_d, ntot, DEVICE_TO_HOST)

      else
         do e = 1, c%msh%nelv
            call mxm(dx, lx, x(1,1,1,e), lx, dxdr(1,1,1,e), lyz)
            call mxm(dx, lx, y(1,1,1,e), lx, dydr(1,1,1,e), lyz)
            call mxm(dx, lx, z(1,1,1,e), lx, dzdr(1,1,1,e), lyz)
            
            do i = 1, lz
               call mxm(x(1,1,i,e), lx, dyt, ly, dxds(1,1,i,e), ly)
               call mxm(y(1,1,i,e), lx, dyt, ly, dyds(1,1,i,e), ly)
               call mxm(z(1,1,i,e), lx, dyt, ly, dzds(1,1,i,e), ly)
            end do
       
            ! We actually take 2d into account, wow, need to do that for the rest.
            if(c%msh%gdim .eq. 3) then
               call mxm(x(1,1,1,e), lxy, dzt, lz, dxdt(1,1,1,e), lz)
               call mxm(y(1,1,1,e), lxy, dzt, lz, dydt(1,1,1,e), lz)
               call mxm(z(1,1,1,e), lxy, dzt, lz, dzdt(1,1,1,e), lz)
            else
               call rzero(dxdt(1,1,1,e), lxy)
               call rzero(dydt(1,1,1,e), lxy)
               call rone(dzdt(1,1,1,e), lxy)
            end if
         end do

         if (c%msh%gdim .eq. 2) then
            call rzero   (jac, ntot)
            call addcol3 (jac, dxdr, dyds, ntot)
            call subcol3 (jac, dxds, dydr, ntot)
            call copy    (drdx, dyds, ntot)
            call copy    (drdy, dxds, ntot)
            call chsign  (drdy, ntot)
            call copy    (dsdx, dydr, ntot)
            call chsign  (dsdx, ntot)
            call copy    (dsdy, dxdr, ntot)
            call rzero   (drdz, ntot)
            call rzero   (dsdz, ntot)
            call rone    (dtdz, ntot)
         else

            do i = 1, ntot               
               c%jac(i, 1, 1, 1) = 0.0_rp
            end do

            do i = 1, ntot               
               c%jac(i, 1, 1, 1) = c%jac(i, 1, 1, 1) + ( c%dxdr(i, 1, 1, 1)  &
                                 * c%dyds(i, 1, 1, 1) * c%dzdt(i, 1, 1, 1) )

               c%jac(i, 1, 1, 1) = c%jac(i, 1, 1, 1) + ( c%dxdt(i, 1, 1, 1)  &
                                 * c%dydr(i, 1, 1, 1) * c%dzds(i, 1, 1, 1) )

               c%jac(i, 1, 1, 1) = c%jac(i, 1, 1, 1) + ( c%dxds(i, 1, 1, 1)  &
                                 * c%dydt(i, 1, 1, 1) * c%dzdr(i, 1, 1, 1) )
            end do

            do i = 1, ntot
               c%jac(i, 1, 1, 1) = c%jac(i, 1, 1, 1) - ( c%dxdr(i, 1, 1, 1)  &
                                 * c%dydt(i, 1, 1, 1) * c%dzds(i, 1, 1, 1) )

               c%jac(i, 1, 1, 1) = c%jac(i, 1, 1, 1) - ( c%dxds(i, 1, 1, 1)  &
                                 * c%dydr(i, 1, 1, 1) * c%dzdt(i, 1, 1, 1) )

               c%jac(i, 1, 1, 1) = c%jac(i, 1, 1, 1) - ( c%dxdt(i, 1, 1, 1)  &
                                 * c%dyds(i, 1, 1, 1) * c%dzdr(i, 1, 1, 1) )
            end do

            do i = 1, ntot
               c%drdx(i, 1, 1, 1) = c%dyds(i, 1, 1, 1) * c%dzdt(i, 1, 1, 1) &
                                  - c%dydt(i, 1, 1, 1) * c%dzds(i, 1, 1, 1)

               c%drdy(i, 1, 1, 1) = c%dxdt(i, 1, 1, 1) * c%dzds(i, 1, 1, 1) &
                                  - c%dxds(i, 1, 1, 1) * c%dzdt(i, 1, 1, 1)
               
               c%drdz(i, 1, 1, 1) = c%dxds(i, 1, 1, 1) * c%dydt(i, 1, 1, 1) &
                                  - c%dxdt(i, 1, 1, 1) * c%dyds(i, 1, 1, 1)
            end do

            do i = 1, ntot
               c%dsdx(i, 1, 1, 1) = c%dydt(i, 1, 1, 1) * c%dzdr(i, 1, 1, 1) &
                                  - c%dydr(i, 1, 1, 1) * c%dzdt(i, 1, 1, 1)

               c%dsdy(i, 1, 1, 1) = c%dxdr(i, 1, 1, 1) * c%dzdt(i, 1, 1, 1) &
                                  - c%dxdt(i, 1, 1, 1) * c%dzdr(i, 1, 1, 1)
               
               c%dsdz(i, 1, 1, 1) = c%dxdt(i, 1, 1, 1) * c%dydr(i, 1, 1, 1) &
                                  - c%dxdr(i, 1, 1, 1) * c%dydt(i, 1, 1, 1)
            end do

            do i = 1, ntot 
               c%dtdx(i, 1, 1, 1) = c%dydr(i, 1, 1, 1) * c%dzds(i, 1, 1, 1) &
                                  - c%dyds(i, 1, 1, 1) * c%dzdr(i, 1, 1, 1)

               c%dtdy(i, 1, 1, 1) = c%dxds(i, 1, 1, 1) * c%dzdr(i, 1, 1, 1) &
                                  - c%dxdr(i, 1, 1, 1) * c%dzds(i, 1, 1, 1)
               
               c%dtdz(i, 1, 1, 1) = c%dxdr(i, 1, 1, 1) * c%dyds(i, 1, 1, 1) &
                                  - c%dxds(i, 1, 1, 1) * c%dydr(i, 1, 1, 1)
            end do

         end if
         call invers2(jacinv, jac, ntot)
      end if
    end associate
    
  end subroutine coef_generate_dxyzdrst
  
  !> Generate geometric data for the given mesh
  !! @note Current implementation assumes regular shaped hex elements
  subroutine coef_generate_geo(c)
    type(coef_t), intent(inout) :: c
    integer :: e, i, lxyz, ntot

    lxyz = c%Xh%lx * c%Xh%ly * c%Xh%lz
    ntot = c%dof%size()

    if (NEKO_BCKND_DEVICE .eq. 1) then

       call device_coef_generate_geo(c%G11_d, c%G12_d, c%G13_d, &
                                     c%G22_d, c%G23_d, c%G33_d, &
                                     c%drdx_d, c%drdy_d, c%drdz_d, &
                                     c%dsdx_d, c%dsdy_d, c%dsdz_d, &
                                     c%dtdx_d, c%dtdy_d, c%dtdz_d, &
                                     c%jacinv_d, c%Xh%w3_d, c%msh%nelv, &
                                     c%Xh%lx, c%msh%gdim)

       call device_memcpy(c%G11, c%G11_d, ntot, DEVICE_TO_HOST)
       call device_memcpy(c%G22, c%G22_d, ntot, DEVICE_TO_HOST)
       call device_memcpy(c%G33, c%G33_d, ntot, DEVICE_TO_HOST)
       call device_memcpy(c%G12, c%G12_d, ntot, DEVICE_TO_HOST)
       call device_memcpy(c%G13, c%G13_d, ntot, DEVICE_TO_HOST)
       call device_memcpy(c%G23, c%G23_d, ntot, DEVICE_TO_HOST)
       
    else
       if(c%msh%gdim .eq. 2) then

          do i = 1, ntot
             c%G11(i, 1, 1, 1) = c%drdx(i, 1, 1, 1) * c%drdx(i, 1, 1, 1) &
                               + c%drdy(i, 1, 1, 1) * c%drdy(i, 1, 1, 1) 

             c%G22(i, 1, 1, 1) = c%dsdx(i, 1, 1, 1) * c%dsdx(i, 1, 1, 1) &
                               + c%dsdy(i, 1, 1, 1) * c%dsdy(i, 1, 1, 1)

             c%G12(i, 1, 1, 1) = c%drdx(i, 1, 1, 1) * c%dsdx(i, 1, 1, 1) &
                               + c%drdy(i, 1, 1, 1) * c%dsdy(i, 1, 1, 1) 
          end do

          do i = 1, ntot
             c%G11(i, 1, 1, 1) = c%G11(i, 1, 1, 1) * c%jacinv(i, 1, 1, 1)
             c%G22(i, 1, 1, 1) = c%G22(i, 1, 1, 1) * c%jacinv(i, 1, 1, 1)
             c%G12(i, 1, 1, 1) = c%G12(i, 1, 1, 1) * c%jacinv(i, 1, 1, 1)
             c%G33(i, 1, 1, 1) = 0.0_rp
             c%G13(i, 1, 1, 1) = 0.0_rp
             c%G23(i, 1, 1, 1) = 0.0_rp
          end do

          do e = 1, c%msh%nelv
             do i = 1, lxyz
               c%G11(i,1,1,e) = c%G11(i,1,1,e) * c%Xh%w3(i,1,1)
               c%G22(i,1,1,e) = c%G22(i,1,1,e) * c%Xh%w3(i,1,1)
               c%G12(i,1,1,e) = c%G12(i,1,1,e) * c%Xh%w3(i,1,1)
            end do
         end do
            
       else

          do i = 1, ntot
             c%G11(i, 1, 1, 1) = c%drdx(i, 1, 1, 1) * c%drdx(i, 1, 1, 1) &
                               + c%drdy(i, 1, 1, 1) * c%drdy(i, 1, 1, 1) &
                               + c%drdz(i, 1, 1, 1) * c%drdz(i, 1, 1, 1)

             c%G22(i, 1, 1, 1) = c%dsdx(i, 1, 1, 1) * c%dsdx(i, 1, 1, 1) &
                               + c%dsdy(i, 1, 1, 1) * c%dsdy(i, 1, 1, 1) &
                               + c%dsdz(i, 1, 1, 1) * c%dsdz(i, 1, 1, 1)

             c%G33(i, 1, 1, 1) = c%dtdx(i, 1, 1, 1) * c%dtdx(i, 1, 1, 1) &
                               + c%dtdy(i, 1, 1, 1) * c%dtdy(i, 1, 1, 1) &
                               + c%dtdz(i, 1, 1, 1) * c%dtdz(i, 1, 1, 1)
          end do

          do i = 1, ntot
             c%G11(i, 1, 1, 1) = c%G11(i, 1, 1, 1) * c%jacinv(i, 1, 1, 1)
             c%G22(i, 1, 1, 1) = c%G22(i, 1, 1, 1) * c%jacinv(i, 1, 1, 1)
             c%G33(i, 1, 1, 1) = c%G33(i, 1, 1, 1) * c%jacinv(i, 1, 1, 1)
          end do

          do i = 1, ntot
             c%G12(i, 1, 1, 1) = c%drdx(i, 1, 1, 1) * c%dsdx(i, 1, 1, 1) &
                               + c%drdy(i, 1, 1, 1) * c%dsdy(i, 1, 1, 1) &
                               + c%drdz(i, 1, 1, 1) * c%dsdz(i, 1, 1, 1)

             c%G13(i, 1, 1, 1) = c%drdx(i, 1, 1, 1) * c%dtdx(i, 1, 1, 1) &
                               + c%drdy(i, 1, 1, 1) * c%dtdy(i, 1, 1, 1) &
                               + c%drdz(i, 1, 1, 1) * c%dtdz(i, 1, 1, 1)

             c%G23(i, 1, 1, 1) = c%dsdx(i, 1, 1, 1) * c%dtdx(i, 1, 1, 1) &
                               + c%dsdy(i, 1, 1, 1) * c%dtdy(i, 1, 1, 1) &
                               + c%dsdz(i, 1, 1, 1) * c%dtdz(i, 1, 1, 1)
          end do

          do i = 1, ntot
             c%G12(i, 1, 1, 1) = c%G12(i, 1, 1, 1) * c%jacinv(i, 1, 1, 1)
             c%G13(i, 1, 1, 1) = c%G13(i, 1, 1, 1) * c%jacinv(i, 1, 1, 1)
             c%G23(i, 1, 1, 1) = c%G23(i, 1, 1, 1) * c%jacinv(i, 1, 1, 1)
          end do

          do e = 1, c%msh%nelv
             do i = 1, lxyz
                c%G11(i,1,1,e) = c%G11(i,1,1,e) * c%Xh%w3(i,1,1)
                c%G22(i,1,1,e) = c%G22(i,1,1,e) * c%Xh%w3(i,1,1)
                c%G12(i,1,1,e) = c%G12(i,1,1,e) * c%Xh%w3(i,1,1)
                
                c%G33(i,1,1,e) = c%G33(i,1,1,e) * c%Xh%w3(i,1,1)
                c%G13(i,1,1,e) = c%G13(i,1,1,e) * c%Xh%w3(i,1,1)
                c%G23(i,1,1,e) = c%G23(i,1,1,e) * c%Xh%w3(i,1,1)
             end do
          end do

       end if
    end if
    
  end subroutine coef_generate_geo
 
  !> Generate mass matrix B for the given mesh and space
  !! @note This is also a stapleholder, we need to go through the coef class properly.
  subroutine coef_generate_mass(c)
    type(coef_t), intent(inout) :: c
    integer :: e, i, lxyz, ntot
    
    lxyz = c%Xh%lx * c%Xh%ly * c%Xh%lz
    ntot = c%dof%size()
    
    !> @todo rewrite this nest into a device kernel
    do e = 1, c%msh%nelv
       ! Here we need to handle things differently for axis symmetric elements
       do i = 1, lxyz
          c%B(i,1,1,e) = c%jac(i,1,1,e) * c%Xh%w3(i,1,1)
          c%Binv(i,1,1,e) = c%B(i,1,1,e)
       end do
    end do

    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_memcpy(c%B, c%B_d, ntot, HOST_TO_DEVICE)
       call device_memcpy(c%Binv, c%Binv_d, ntot, HOST_TO_DEVICE)
    end if
    
    call gs_op(c%gs_h, c%Binv, ntot, GS_OP_ADD)

    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_invcol1(c%Binv_d, ntot)
       call device_memcpy(c%Binv, c%Binv_d, ntot, DEVICE_TO_HOST)
    else
       call invcol1(c%Binv, ntot)
    end if

    !>  @todo cleanup once we have device math in place
    if (NEKO_BCKND_DEVICE .eq. 1) then
       c%volume = device_glsum(c%B_d, ntot)
    else
       c%volume = glsum(c%B, ntot)
    end if

  end subroutine coef_generate_mass

  pure function coef_get_normal(coef, i, j, k, e, facet) result(normal)
    type(coef_t), intent(in) :: coef
    integer, intent(in) :: i, j, k, e, facet
    real(kind=rp) :: normal(3)
      
    select case (facet)               
      case(1,2)
        normal(1) = coef%nx(j, k, facet, e)
        normal(2) = coef%ny(j, k, facet, e)
        normal(3) = coef%nz(j, k, facet, e)
      case(3,4)
        normal(1) = coef%nx(i, k, facet, e)
        normal(2) = coef%ny(i, k, facet, e)
        normal(3) = coef%nz(i, k, facet, e)
      case(5,6)
        normal(1) = coef%nx(i, j, facet, e)
        normal(2) = coef%ny(i, j, facet, e)
        normal(3) = coef%nz(i, j, facet, e)
      end select
  end function coef_get_normal

  !> Generate facet area and surface normals
  subroutine coef_generate_area_and_normal(coef)
    type(coef_t), intent(inout) :: coef
    real(kind=rp), allocatable :: a(:,:,:,:)
    real(kind=rp), allocatable :: b(:,:,:,:)
    real(kind=rp), allocatable :: c(:,:,:,:)
    real(kind=rp), allocatable :: dot(:,:,:,:)
    integer :: n, e, i, j, k, lx
    real(kind=rp) :: weight, len
    n = coef%dof%size()
    lx = coef%Xh%lx
    
    allocate(a(coef%Xh%lx, coef%Xh%lx, coef%Xh%lx, coef%msh%nelv))
    allocate(b(coef%Xh%lx, coef%Xh%lx, coef%Xh%lx, coef%msh%nelv))
    allocate(c(coef%Xh%lx, coef%Xh%lx, coef%Xh%lx, coef%msh%nelv))
    allocate(dot(coef%Xh%lx, coef%Xh%lx, coef%Xh%lx, coef%msh%nelv))

    ! ds x dt
    do i = 1, n
       a(i, 1, 1, 1) = coef%dyds(i, 1, 1, 1) * coef%dzdt(i, 1, 1, 1) &
                     - coef%dzds(i, 1, 1, 1) * coef%dydt(i, 1, 1, 1)
 
       b(i, 1, 1, 1) = coef%dzds(i, 1, 1, 1) * coef%dxdt(i, 1, 1, 1) &
                     - coef%dxds(i, 1, 1, 1) * coef%dzdt(i, 1, 1, 1)

       c(i, 1, 1, 1) = coef%dxds(i, 1, 1, 1) * coef%dydt(i, 1, 1, 1) &
                     - coef%dyds(i, 1, 1, 1) * coef%dxdt(i, 1, 1, 1)
    end do

    do i = 1, n
       dot(i, 1, 1, 1) = a(i, 1, 1, 1) * a(i, 1, 1, 1) &
                       + b(i, 1, 1, 1) * b(i, 1, 1, 1) &
                       + c(i, 1, 1, 1) * c(i, 1, 1, 1)
    end do

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

    ! dr x dt
    do i = 1, n
       a(i, 1, 1, 1) = coef%dydr(i, 1, 1, 1) * coef%dzdt(i, 1, 1, 1) &
                     - coef%dzdr(i, 1, 1, 1) * coef%dydt(i, 1, 1, 1)
 
       b(i, 1, 1, 1) = coef%dzdr(i, 1, 1, 1) * coef%dxdt(i, 1, 1, 1) &
                     - coef%dxdr(i, 1, 1, 1) * coef%dzdt(i, 1, 1, 1)

       c(i, 1, 1, 1) = coef%dxdr(i, 1, 1, 1) * coef%dydt(i, 1, 1, 1) &
                     - coef%dydr(i, 1, 1, 1) * coef%dxdt(i, 1, 1, 1)
    end do

    do i = 1, n
       dot(i, 1, 1, 1) = a(i, 1, 1, 1) * a(i, 1, 1, 1) &
                       + b(i, 1, 1, 1) * b(i, 1, 1, 1) &
                       + c(i, 1, 1, 1) * c(i, 1, 1, 1)
    end do

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

    ! dr x ds
    do i = 1, n
       a(i, 1, 1, 1) = coef%dydr(i, 1, 1, 1) * coef%dzds(i, 1, 1, 1) &
                     - coef%dzdr(i, 1, 1, 1) * coef%dyds(i, 1, 1, 1)
 
       b(i, 1, 1, 1) = coef%dzdr(i, 1, 1, 1) * coef%dxds(i, 1, 1, 1) &
                     - coef%dxdr(i, 1, 1, 1) * coef%dzds(i, 1, 1, 1)

       c(i, 1, 1, 1) = coef%dxdr(i, 1, 1, 1) * coef%dyds(i, 1, 1, 1) &
                     - coef%dydr(i, 1, 1, 1) * coef%dxds(i, 1, 1, 1)
    end do

    do i = 1, n
       dot(i, 1, 1, 1) = a(i, 1, 1, 1) * a(i, 1, 1, 1) &
                       + b(i, 1, 1, 1) * b(i, 1, 1, 1) &
                       + c(i, 1, 1, 1) * c(i, 1, 1, 1)
    end do

    do e = 1, coef%msh%nelv
       do k = 1, coef%Xh%lx
          do j = 1, coef%Xh%lx
             weight = coef%Xh%wx(j) * coef%Xh%wy(k)
             coef%area(j, k, 5, e) = sqrt(dot(j, k, 1, e)) * weight
             coef%area(j, k, 6, e) = sqrt(dot(j, k, lx, e)) * weight
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
    do j = 1, size(coef%nz)
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
    if (NEKO_BCKND_DEVICE .eq. 1) then
       n = size(coef%area)
       call device_memcpy(coef%area, coef%area_d, n, HOST_TO_DEVICE)
       call device_memcpy(coef%nx, coef%nx_d, n, HOST_TO_DEVICE)
       call device_memcpy(coef%ny, coef%ny_d, n, HOST_TO_DEVICE)
       call device_memcpy(coef%nz, coef%nz_d, n, HOST_TO_DEVICE)
    end if
    
  end subroutine coef_generate_area_and_normal
  
end module coefs
