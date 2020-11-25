!> Coefficients 
module coefs
  use gather_scatter
  use num_types
  use space  
  use math
  use mesh
  implicit none
  private
  
  !> Coefficients defined on a given (mesh, \f$ X_h \f$) tuple
  type, public :: coef_t
     
     real(kind=dp), allocatable :: G1(:,:,:,:) !< Geometric data
     real(kind=dp), allocatable :: G2(:,:,:,:) !< Geometric data
     real(kind=dp), allocatable :: G3(:,:,:,:) !< Geometric data
     real(kind=dp), allocatable :: G4(:,:,:,:) !< Geometric data
     real(kind=dp), allocatable :: G5(:,:,:,:) !< Geometric data
     real(kind=dp), allocatable :: G6(:,:,:,:) !< Geometric data

     real(kind=dp), allocatable :: mult(:,:,:,:) !< Multiplicity
     
     type(space_t), pointer :: Xh => null()
     type(mesh_t), pointer :: msh => null()
  end type coef_t

  public :: coef_init, coef_free
  
contains

  !> Initialize coefficients
  subroutine coef_init(coef, gs_h)
    type(coef_t), intent(inout) :: coef
    type(gs_t), intent(inout), target :: gs_h
    integer :: n

    call coef_free(coef)
    
    coef%msh => gs_h%dofmap%msh
    coef%Xh => gs_h%dofmap%Xh

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

    call coef_generate_geom(coef)

    !
    ! Set up multiplicity
    !
    n = coef%Xh%lx * coef%Xh%ly * coef%Xh%lz * coef%msh%nelv
    allocate(coef%mult(coef%Xh%lx, coef%Xh%ly, coef%Xh%lz, coef%msh%nelv))
    call rone(coef%mult, n)   
    call gs_op_vector(gs_h, coef%mult, n, GS_OP_ADD)
    call invcol1(coef%mult, n)
    
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

    coef%msh => null()
    coef%Xh => null()
    
  end subroutine coef_free

  !> Generate geometric data for the given mesh
  !! @note Current implementation assumes regular shaped hex elements
  subroutine coef_generate_geom(c)
    type(coef_t), intent(inout) :: c
    integer :: e, j, k, l
    
    do e = 1, c%msh%nelv
       do l = 1, c%Xh%lz
          do k = 1, c%Xh%ly
             do j = 1, c%Xh%lx
                c%G1(j, k, l, e) = c%Xh%wx(j) * c%Xh%wx(k) * c%Xh%wx(l)
                c%G2(j, k, l, e) = 0d0
                c%G3(j, k, l, e) = 0d0
                c%G4(j, k, l, e) = c%Xh%wx(j) * c%Xh%wx(k) * c%Xh%wx(l)
                c%G5(j, k, l, e) = 0d0
                c%G6(j, k, l, e) = c%Xh%wx(j) * c%Xh%wx(k) * c%Xh%wx(l)
             end do
          end do
       end do
    end do
    
  end subroutine coef_generate_geom
  
  
end module coefs
