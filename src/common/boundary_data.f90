! Copyright (c) 2026, The Neko Authors
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
!> Implements the `boundary_data_t` type.
module boundary_data
  use num_types, only : rp
  use coefs, only : coef_t
  use dirichlet, only : dirichlet_t
  use field, only : field_t
  use vector, only : vector_t
  use registry, only : neko_registry
  use mesh, only : NEKO_MSH_MAX_ZLBLS
  use field_math, only : field_rzero
  use vector_math, only : vector_masked_gather_copy_0, &
       vector_masked_scatter_copy_0, vector_face_masked_gather_copy_0, &
       vector_glsc2, vector_glsc3, vector_glsum, &
       vector_vdot3, vector_subcol3, vector_col3, vector_cmult, vector_copy
  use neko_config, only : NEKO_BCKND_DEVICE
  use utils, only : neko_error
  use comm, only : NEKO_COMM
  use mpi_f08, only : MPI_Allreduce, MPI_INTEGER, MPI_SUM
  implicit none
  private

  !> Collects data on the boundary points of one or more labelled zones
  !! and perform some bounary operations.
  type, public :: boundary_data_t
     !> SEM coefficients.
     type(coef_t), pointer :: coef => null()
     !> Boundary mask built from the chosen zones.
     type(dirichlet_t) :: bc
     !> Labelled zones included.
     integer, allocatable :: zone_indices(:)
     !> Number of boundary points on this rank.
     integer :: n_local = 0
     !> Number of boundary points over all ranks.
     integer :: n_global = 0
     !> Coordinates at the boundary points.
     type(vector_t) :: x, y, z
     !> Unit normal components at the boundary points.
     type(vector_t) :: n_x, n_y, n_z
     !> Surface quadrature weights at the boundary points.
     type(vector_t) :: area
     !> Work vector.
     type(vector_t), private :: work
     !> Whether the stored normals point out of the wall into the fluid. When
     !! false they keep the `coef` convention, pointing out of the fluid.
     logical :: outward_normals = .true.
   contains
     !> Constructor.
     procedure, pass(this) :: init => boundary_data_init
     !> Destructor.
     procedure, pass(this) :: free => boundary_data_free
     !> Re-gather the geometry factors at the boundary points.
     procedure, pass(this) :: update_geometry => &
          boundary_data_update_geometry
     !> Sample a quantity into a boundary vector.
     generic :: get => get_vector_by_name, get_vector_by_field
     procedure, private, pass(this) :: get_vector_by_name => &
          boundary_data_get_vector_by_name
     procedure, private, pass(this) :: get_vector_by_field => &
          boundary_data_get_vector_by_field
     !> Scatter boundary values back into a full field.
     generic :: scatter => scatter_to_field_by_vector, &
          scatter_to_field_by_name
     procedure, private, pass(this) :: scatter_to_field_by_vector => &
          boundary_data_scatter_to_field_by_vector
     procedure, private, pass(this) :: scatter_to_field_by_name => &
          boundary_data_scatter_to_field_by_name
     !> Whether a name refers to a geometry attribute.
     procedure, nopass :: is_geometry => boundary_data_is_geometry
     !> Surface integral of a quantity over the zones.
     generic :: integrate => integrate_by_name, integrate_by_field
     procedure, private, pass(this) :: integrate_by_name => &
          boundary_data_integrate_by_name
     procedure, private, pass(this) :: integrate_by_field => &
          boundary_data_integrate_by_field
     !> Area weighted average of a quantity over the zones.
     generic :: average => average_by_name, average_by_field
     procedure, private, pass(this) :: average_by_name => &
          boundary_data_average_by_name
     procedure, private, pass(this) :: average_by_field => &
          boundary_data_average_by_field
     !> Total surface area of the zones.
     procedure, pass(this) :: surface_area => boundary_data_surface_area
     !> Area weighted geometric centre of the zones.
     procedure, pass(this) :: centroid => boundary_data_centroid
     !> Unweighted mean of the boundary point coordinates.
     procedure, pass(this) :: point_average => boundary_data_point_average
     !> Surface integral of a vector quantity, component by component.
     generic :: integrate_vector => integrate_vector_by_name, &
          integrate_vector_by_field
     procedure, private, pass(this) :: integrate_vector_by_name => &
          boundary_data_integrate_vector_by_name
     procedure, private, pass(this) :: integrate_vector_by_field => &
          boundary_data_integrate_vector_by_field
     !> Surface integral of a scalar times the normal.
     generic :: integrate_normal => integrate_normal_by_name, &
          integrate_normal_by_field
     procedure, private, pass(this) :: integrate_normal_by_name => &
          boundary_data_integrate_normal_by_name
     procedure, private, pass(this) :: integrate_normal_by_field => &
          boundary_data_integrate_normal_by_field
     !> Flux of a vector quantity through the zones.
     generic :: flux => flux_by_name, flux_by_field
     procedure, private, pass(this) :: flux_by_name => &
          boundary_data_flux_by_name
     procedure, private, pass(this) :: flux_by_field => &
          boundary_data_flux_by_field
     !> Wall tangential part of a vector.
     generic :: tangential => tangential_inplace, tangential_split
     procedure, private, pass(this) :: tangential_inplace => &
          boundary_data_tangential_inplace
     procedure, private, pass(this) :: tangential_split => &
          boundary_data_tangential_split
     !> Wall normal part of a vector.
     generic :: normal => normal_inplace, normal_split
     procedure, private, pass(this) :: normal_inplace => &
          boundary_data_normal_inplace
     procedure, private, pass(this) :: normal_split => &
          boundary_data_normal_split
  end type boundary_data_t

contains

  !> Build the boundary point mask and gather the geometry.
  !! @param coef The SEM coefficients.
  !! @param zone_indices The labelled zones to include.
  !! @param outward_normals Whether the normals should point out of the wall
  !! into the fluid, which is the default. Pass false to keep the `coef`
  !! convention, where they point out of the fluid domain.
  subroutine boundary_data_init(this, coef, zone_indices, outward_normals)
    class(boundary_data_t), intent(inout) :: this
    type(coef_t), intent(inout), target :: coef
    integer, intent(in) :: zone_indices(:)
    logical, intent(in), optional :: outward_normals
    integer :: i, ierr

    call this%free()

    this%coef => coef
    this%outward_normals = .true.
    if (present(outward_normals)) this%outward_normals = outward_normals

    if (size(zone_indices) .eq. 0) then
       call neko_error("boundary_data: at least one zone index is required")
    end if
    do i = 1, size(zone_indices)
       if (zone_indices(i) .lt. 1 .or. &
            zone_indices(i) .gt. NEKO_MSH_MAX_ZLBLS) then
          call neko_error("boundary_data: zone index out of range")
       end if
    end do

    allocate(this%zone_indices(size(zone_indices)))
    this%zone_indices = zone_indices

    ! The face-node mask (facet_node_msk) is used.
    call this%bc%init_base(this%coef)
    this%bc%zone_indices = this%zone_indices
    do i = 1, size(this%zone_indices)
       call this%bc%mark_zone( &
            this%coef%dof%msh%labeled_zones(this%zone_indices(i)))
    end do
    call this%bc%finalize()

    this%n_local = this%bc%facet_node_msk(0)

    call MPI_Allreduce(this%n_local, this%n_global, 1, MPI_INTEGER, &
         MPI_SUM, NEKO_COMM, ierr)

    if (this%n_global .eq. 0) then
       call neko_error("boundary_data: the requested zones contain no " // &
            "boundary points")
    end if

    call this%x%init(this%n_local)
    call this%y%init(this%n_local)
    call this%z%init(this%n_local)
    call this%n_x%init(this%n_local)
    call this%n_y%init(this%n_local)
    call this%n_z%init(this%n_local)
    call this%area%init(this%n_local)
    call this%work%init(this%n_local)

    ! Populate the geometry once at construction.
    call this%update_geometry()

  end subroutine boundary_data_init

  !> Destructor.
  subroutine boundary_data_free(this)
    class(boundary_data_t), intent(inout) :: this

    call this%bc%free()

    call this%x%free()
    call this%y%free()
    call this%z%free()
    call this%n_x%free()
    call this%n_y%free()
    call this%n_z%free()
    call this%area%free()
    call this%work%free()

    if (allocated(this%zone_indices)) deallocate(this%zone_indices)

    this%n_local = 0
    this%n_global = 0
    nullify(this%coef)

  end subroutine boundary_data_free

  !> Re-gather the coordinates, normals and surface weights.
  subroutine boundary_data_update_geometry(this)
    class(boundary_data_t), intent(inout) :: this
    integer :: n

    if (this%n_local .le. 0) return

    n = this%coef%dof%size()

    call vector_masked_gather_copy_0(this%x, this%coef%dof%x, &
         this%bc%facet_node_msk, n, this%n_local)
    call vector_masked_gather_copy_0(this%y, this%coef%dof%y, &
         this%bc%facet_node_msk, n, this%n_local)
    call vector_masked_gather_copy_0(this%z, this%coef%dof%z, &
         this%bc%facet_node_msk, n, this%n_local)

    call vector_face_masked_gather_copy_0(this%n_x, this%coef%nx, &
         this%bc%facet_node_msk, this%bc%facet, this%coef%Xh%lx, &
         this%coef%Xh%ly, this%coef%Xh%lz, this%n_local)
    call vector_face_masked_gather_copy_0(this%n_y, this%coef%ny, &
         this%bc%facet_node_msk, this%bc%facet, this%coef%Xh%lx, &
         this%coef%Xh%ly, this%coef%Xh%lz, this%n_local)
    call vector_face_masked_gather_copy_0(this%n_z, this%coef%nz, &
         this%bc%facet_node_msk, this%bc%facet, this%coef%Xh%lx, &
         this%coef%Xh%ly, this%coef%Xh%lz, this%n_local)

    call vector_face_masked_gather_copy_0(this%area, this%coef%area, &
         this%bc%facet_node_msk, this%bc%facet, this%coef%Xh%lx, &
         this%coef%Xh%ly, this%coef%Xh%lz, this%n_local)

    if (this%outward_normals) then
       call vector_cmult(this%n_x, -1.0_rp)
       call vector_cmult(this%n_y, -1.0_rp)
       call vector_cmult(this%n_z, -1.0_rp)
    end if

  end subroutine boundary_data_update_geometry


  !> Whether `name` refers to one of the geometry attributes.
  !! @param name The requested name.
  pure function boundary_data_is_geometry(name) result(is_geom)
    character(len=*), intent(in) :: name
    logical :: is_geom

    select case (trim(name))
    case ("x", "y", "z", "n_x", "n_y", "n_z", "area")
       is_geom = .true.
    case default
       is_geom = .false.
    end select

  end function boundary_data_is_geometry

  !> Sample a named quantity at the boundary points into a vector.
  !! @details The name is either one of the geometry attributes `x`, `y`, `z`,
  !! `n_x`, `n_y`, `n_z`, `area`, or the name of a field in the registry.
  !! @param name The requested quantity.
  !! @param v The vector to fill.
  subroutine boundary_data_get_vector_by_name(this, name, v)
    class(boundary_data_t), intent(inout) :: this
    character(len=*), intent(in) :: name
    type(vector_t), intent(inout) :: v
    type(field_t), pointer :: f

    if (v%size() .ne. this%n_local) then
       call v%free()
       call v%init(this%n_local)
    end if
    if (this%n_local .le. 0) return

    select case (trim(name))
    case ("x")
       call vector_copy(v, this%x)
    case ("y")
       call vector_copy(v, this%y)
    case ("z")
       call vector_copy(v, this%z)
    case ("n_x")
       call vector_copy(v, this%n_x)
    case ("n_y")
       call vector_copy(v, this%n_y)
    case ("n_z")
       call vector_copy(v, this%n_z)
    case ("area")
       call vector_copy(v, this%area)
    case default
       if (.not. neko_registry%field_exists(trim(name))) then
          call neko_error("boundary_data: '" // trim(name) // &
               "' is neither a geometry attribute (x, y, z, n_x, n_y, " // &
               "n_z, area) nor a field in the registry")
       end if
       f => neko_registry%get_field_by_name(trim(name))
       call this%get_vector_by_field(f, v)
    end select

  end subroutine boundary_data_get_vector_by_name

  !> Sample a field at the boundary points into a vector.
  !! @param f The source field, which must live on the same dofmap.
  !! @param v The vector to fill.
  subroutine boundary_data_get_vector_by_field(this, f, v)
    class(boundary_data_t), intent(inout) :: this
    type(field_t), intent(in) :: f
    type(vector_t), intent(inout) :: v

    if (.not. associated(f%dof, this%coef%dof)) then
       call neko_error("boundary_data: the field '" // trim(f%name) // &
            "' is on a different dofmap than the boundary mask, so the " // &
            "masked indices do not apply to it")
    end if

    if (v%size() .ne. this%n_local) then
       call v%free()
       call v%init(this%n_local)
    end if
    if (this%n_local .le. 0) return

    call vector_masked_gather_copy_0(v, f%x, this%bc%facet_node_msk, &
         this%coef%dof%size(), this%n_local)

  end subroutine boundary_data_get_vector_by_field

  !> Scatter a boundary vector back into a full field.
  !! @param v The boundary point values.
  !! @param f The destination field.
  !! @param clear If true, zero the whole field before writing the
  !! boundary points; defaults to false.
  subroutine boundary_data_scatter_to_field_by_vector(this, v, f, clear)
    class(boundary_data_t), intent(inout) :: this
    type(vector_t), intent(in) :: v
    type(field_t), intent(inout) :: f
    logical, intent(in), optional :: clear
    logical :: clear_
    integer :: n

    if (.not. associated(f%dof, this%coef%dof)) then
       call neko_error("boundary_data: the destination field is on a " // &
            "different dofmap than the boundary mask")
    end if

    clear_ = .false.
    if (present(clear)) clear_ = clear

    n = this%coef%dof%size()

    if (clear_) call field_rzero(f)

    call vector_masked_scatter_copy_0(f%x, v, this%bc%facet_node_msk, n, &
         this%n_local)

  end subroutine boundary_data_scatter_to_field_by_vector

  !> Scatter a named boundary quantity into a full field.
  !! @param name The requested quantity.
  !! @param f The destination field.
  !! @param clear If true, zero the whole field before writing the
  !! boundary points; defaults to false.
  subroutine boundary_data_scatter_to_field_by_name(this, name, f, clear)
    class(boundary_data_t), intent(inout) :: this
    character(len=*), intent(in) :: name
    type(field_t), intent(inout) :: f
    logical, intent(in), optional :: clear
    logical :: clear_

    clear_ = .false.
    if (present(clear)) clear_ = clear

    call this%get(trim(name), this%work)
    call this%scatter_to_field_by_vector(this%work, f, clear_)

  end subroutine boundary_data_scatter_to_field_by_name

  !> Keep only the wall tangential part of a vector at the boundary points,
  !! in place. Each component is overwritten with `v - (v . n) n`.
  !! @note The result does not depend on the sign of the normal, since the
  !! normal enters twice, so it is the same whether the stored normals are the
  !! outward or not.
  !! @param vx The x component, overwritten with its tangential part.
  !! @param vy The y component, overwritten with its tangential part.
  !! @param vz The z component, overwritten with its tangential part.
  subroutine boundary_data_tangential_inplace(this, vx, vy, vz)
    class(boundary_data_t), intent(inout) :: this
    type(vector_t), intent(inout) :: vx, vy, vz

    if (vx%size() .lt. this%n_local .or. vy%size() .lt. this%n_local .or. &
         vz%size() .lt. this%n_local) then
       call neko_error("boundary_data: the vectors passed to tangential " // &
            "are shorter than the number of boundary points")
    end if

    call vector_vdot3(this%work, vx, vy, vz, this%n_x, this%n_y, this%n_z)
    call vector_subcol3(vx, this%work, this%n_x)
    call vector_subcol3(vy, this%work, this%n_y)
    call vector_subcol3(vz, this%work, this%n_z)

  end subroutine boundary_data_tangential_inplace

  !> Wall tangential part of a vector at the boundary points, out of place.
  !! The input `(vx, vy, vz)` is left unchanged and `(tx, ty, tz)` receive
  !! `v - (v . n) n`.
  !! @note The result does not depend on the sign of the normal, since the
  !! normal enters twice, so it is the same whether the stored normals are the
  !! outward or not.
  !! @param vx The input x component.
  !! @param vy The input y component.
  !! @param vz The input z component.
  !! @param tx The output x component of the tangential part.
  !! @param ty The output y component of the tangential part.
  !! @param tz The output z component of the tangential part.
  subroutine boundary_data_tangential_split(this, vx, vy, vz, tx, ty, tz)
    class(boundary_data_t), intent(inout) :: this
    type(vector_t), intent(in) :: vx, vy, vz
    type(vector_t), intent(inout) :: tx, ty, tz

    if (vx%size() .lt. this%n_local .or. vy%size() .lt. this%n_local .or. &
         vz%size() .lt. this%n_local .or. tx%size() .lt. this%n_local .or. &
         ty%size() .lt. this%n_local .or. tz%size() .lt. this%n_local) then
       call neko_error("boundary_data: the vectors passed to tangential " // &
            "are shorter than the number of boundary points")
    end if

    call vector_vdot3(this%work, vx, vy, vz, this%n_x, this%n_y, this%n_z)
    call vector_copy(tx, vx)
    call vector_copy(ty, vy)
    call vector_copy(tz, vz)
    call vector_subcol3(tx, this%work, this%n_x)
    call vector_subcol3(ty, this%work, this%n_y)
    call vector_subcol3(tz, this%work, this%n_z)

  end subroutine boundary_data_tangential_split

  !> Keep only the wall normal part of a vector at the boundary points, in
  !! place. Each component is overwritten with `(v . n) n`.
  !! @note The result does not depend on the sign of the normal, since the
  !! normal enters twice, so it is the same whether the stored normals are the
  !! outward or not.
  !! @param vx The x component, overwritten with its normal part.
  !! @param vy The y component, overwritten with its normal part.
  !! @param vz The z component, overwritten with its normal part.
  subroutine boundary_data_normal_inplace(this, vx, vy, vz)
    class(boundary_data_t), intent(inout) :: this
    type(vector_t), intent(inout) :: vx, vy, vz

    if (vx%size() .lt. this%n_local .or. vy%size() .lt. this%n_local .or. &
         vz%size() .lt. this%n_local) then
       call neko_error("boundary_data: the vectors passed to normal " // &
            "are shorter than the number of boundary points")
    end if

    call vector_vdot3(this%work, vx, vy, vz, this%n_x, this%n_y, this%n_z)
    call vector_col3(vx, this%work, this%n_x)
    call vector_col3(vy, this%work, this%n_y)
    call vector_col3(vz, this%work, this%n_z)

  end subroutine boundary_data_normal_inplace

  !> Wall normal part of a vector at the boundary points, out of place. The
  !! input `(vx, vy, vz)` is left unchanged and `(nx, ny, nz)` receive
  !! `(v . n) n`.
  !! @note The result does not depend on the sign of the normal, since the
  !! normal enters twice, so it is the same whether the stored normals are the
  !! outward or not.
  !! @param vx The input x component.
  !! @param vy The input y component.
  !! @param vz The input z component.
  !! @param nx The output x component of the normal part.
  !! @param ny The output y component of the normal part.
  !! @param nz The output z component of the normal part.
  subroutine boundary_data_normal_split(this, vx, vy, vz, nx, ny, nz)
    class(boundary_data_t), intent(inout) :: this
    type(vector_t), intent(in) :: vx, vy, vz
    type(vector_t), intent(inout) :: nx, ny, nz

    if (vx%size() .lt. this%n_local .or. vy%size() .lt. this%n_local .or. &
         vz%size() .lt. this%n_local .or. nx%size() .lt. this%n_local .or. &
         ny%size() .lt. this%n_local .or. nz%size() .lt. this%n_local) then
       call neko_error("boundary_data: the vectors passed to normal " // &
            "are shorter than the number of boundary points")
    end if

    call vector_vdot3(this%work, vx, vy, vz, this%n_x, this%n_y, this%n_z)
    call vector_col3(nx, this%work, this%n_x)
    call vector_col3(ny, this%work, this%n_y)
    call vector_col3(nz, this%work, this%n_z)

  end subroutine boundary_data_normal_split

  !> Surface integral of a named quantity over the zones.
  !! @param name The quantity to integrate.
  function boundary_data_integrate_by_name(this, name) result(val)
    class(boundary_data_t), intent(inout) :: this
    character(len=*), intent(in) :: name
    real(kind=rp) :: val

    call this%get(trim(name), this%work)
    val = vector_glsc2(this%work, this%area)

  end function boundary_data_integrate_by_name

  !> Surface integral of a field over the zones.
  !! @param f The field to integrate, on the same dofmap.
  function boundary_data_integrate_by_field(this, f) result(val)
    class(boundary_data_t), intent(inout) :: this
    type(field_t), intent(in) :: f
    real(kind=rp) :: val

    call this%get(f, this%work)
    val = vector_glsc2(this%work, this%area)

  end function boundary_data_integrate_by_field

  !> Total surface area of the zones.
  function boundary_data_surface_area(this) result(val)
    class(boundary_data_t), intent(inout) :: this
    real(kind=rp) :: val

    val = vector_glsum(this%area)

  end function boundary_data_surface_area

  !> Area weighted average of a named quantity over the zones.
  !! @param name The quantity to average.
  function boundary_data_average_by_name(this, name) result(val)
    class(boundary_data_t), intent(inout) :: this
    character(len=*), intent(in) :: name
    real(kind=rp) :: val
    real(kind=rp) :: total

    total = this%surface_area()
    if (total .le. 0.0_rp) then
       call neko_error("boundary_data: the total surface area is not " // &
            "positive, cannot form an area weighted average")
    end if
    val = this%integrate(trim(name)) / total

  end function boundary_data_average_by_name

  !> Area weighted average of a field over the zones.
  !! @param f The field to average.
  function boundary_data_average_by_field(this, f) result(val)
    class(boundary_data_t), intent(inout) :: this
    type(field_t), intent(in) :: f
    real(kind=rp) :: val
    real(kind=rp) :: total

    total = this%surface_area()
    if (total .le. 0.0_rp) then
       call neko_error("boundary_data: the total surface area is not " // &
            "positive, cannot form an area weighted average")
    end if
    val = this%integrate(f) / total

  end function boundary_data_average_by_field

  !> Area weighted geometric centre of the zones.
  function boundary_data_centroid(this) result(c)
    class(boundary_data_t), intent(inout) :: this
    real(kind=rp) :: c(3)
    real(kind=rp) :: total

    total = vector_glsum(this%area)
    if (total .le. 0.0_rp) then
       call neko_error("boundary_data: the total surface area is not " // &
            "positive, cannot form a centroid")
    end if

    c(1) = vector_glsc2(this%x, this%area) / total
    c(2) = vector_glsc2(this%y, this%area) / total
    c(3) = vector_glsc2(this%z, this%area) / total

  end function boundary_data_centroid

  !> Unweighted mean of the boundary point coordinates.
  function boundary_data_point_average(this) result(c)
    class(boundary_data_t), intent(inout) :: this
    real(kind=rp) :: c(3)

    c(1) = vector_glsum(this%x) / this%n_global
    c(2) = vector_glsum(this%y) / this%n_global
    c(3) = vector_glsum(this%z) / this%n_global

  end function boundary_data_point_average

  !> Surface integral of a vector quantity, component by component.
  !! @param fx The x component.
  !! @param fy The y component.
  !! @param fz The z component.
  function boundary_data_integrate_vector_by_field(this, fx, fy, fz) &
       result(f)
    class(boundary_data_t), intent(inout) :: this
    type(field_t), intent(in) :: fx, fy, fz
    real(kind=rp) :: f(3)

    f(1) = this%integrate(fx)
    f(2) = this%integrate(fy)
    f(3) = this%integrate(fz)

  end function boundary_data_integrate_vector_by_field

  !> Surface integral of a named vector quantity, component by component.
  !! @param fx The name of the x component.
  !! @param fy The name of the y component.
  !! @param fz The name of the z component.
  function boundary_data_integrate_vector_by_name(this, fx, fy, fz) result(f)
    class(boundary_data_t), intent(inout) :: this
    character(len=*), intent(in) :: fx, fy, fz
    real(kind=rp) :: f(3)

    f(1) = this%integrate(trim(fx))
    f(2) = this%integrate(trim(fy))
    f(3) = this%integrate(trim(fz))

  end function boundary_data_integrate_vector_by_name

  !> Flux of a vector quantity through the zones.
  !! @note The sign follows the stored normals. With the default outward
  !! convention a positive value is a flux out of the wall into the fluid.
  !! @param u The x component.
  !! @param v The y component.
  !! @param w The z component.
  function boundary_data_flux_by_field(this, u, v, w) result(q)
    class(boundary_data_t), intent(inout) :: this
    type(field_t), intent(in) :: u, v, w
    real(kind=rp) :: q

    call this%get(u, this%work)
    q = vector_glsc3(this%work, this%n_x, this%area)
    call this%get(v, this%work)
    q = q + vector_glsc3(this%work, this%n_y, this%area)
    call this%get(w, this%work)
    q = q + vector_glsc3(this%work, this%n_z, this%area)

  end function boundary_data_flux_by_field

  !> Flux of a named vector quantity through the zones.
  !! @param u The name of the x component.
  !! @param v The name of the y component.
  !! @param w The name of the z component.
  function boundary_data_flux_by_name(this, u, v, w) result(q)
    class(boundary_data_t), intent(inout) :: this
    character(len=*), intent(in) :: u, v, w
    real(kind=rp) :: q

    call this%get(trim(u), this%work)
    q = vector_glsc3(this%work, this%n_x, this%area)
    call this%get(trim(v), this%work)
    q = q + vector_glsc3(this%work, this%n_y, this%area)
    call this%get(trim(w), this%work)
    q = q + vector_glsc3(this%work, this%n_z, this%area)

  end function boundary_data_flux_by_name

  !> Surface integral of a scalar quantity times the normal.
  !! @note The sign follows the stored normals, so it changes with
  !! `outward_normals`.
  !! @param f The scalar field.
  function boundary_data_integrate_normal_by_field(this, f) result(fn)
    class(boundary_data_t), intent(inout) :: this
    type(field_t), intent(in) :: f
    real(kind=rp) :: fn(3)

    call this%get(f, this%work)
    fn(1) = vector_glsc3(this%work, this%n_x, this%area)
    fn(2) = vector_glsc3(this%work, this%n_y, this%area)
    fn(3) = vector_glsc3(this%work, this%n_z, this%area)

  end function boundary_data_integrate_normal_by_field

  !> Surface integral of a named scalar quantity times the normal.
  !! @param name The scalar quantity.
  function boundary_data_integrate_normal_by_name(this, name) result(fn)
    class(boundary_data_t), intent(inout) :: this
    character(len=*), intent(in) :: name
    real(kind=rp) :: fn(3)

    call this%get(trim(name), this%work)
    fn(1) = vector_glsc3(this%work, this%n_x, this%area)
    fn(2) = vector_glsc3(this%work, this%n_y, this%area)
    fn(3) = vector_glsc3(this%work, this%n_z, this%area)

  end function boundary_data_integrate_normal_by_name

end module boundary_data
