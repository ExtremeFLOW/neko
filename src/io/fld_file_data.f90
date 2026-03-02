!> Simple module to handle fld file series.
!! Provides an interface to the different fields sotred in a fld file
!! Also provides simple functions to scale and add different fld files.
!! An example of using this module is shown in contrib/average_fields.f90
!! The fld_file_data_t should dynamically update each time one reads a new fld
!! file.
!! Martin Karp 1/2-2023
module fld_file_data
  use num_types, only : rp
  use math, only : cmult, add2
  use vector, only : vector_t, vector_ptr_t
  use interpolation, only : interpolator_t
  use field, only : field_t
  use field_list, only : field_list_t
  use logger, only : neko_log, LOG_SIZE
  use device, only : HOST_TO_DEVICE
  use dofmap, only : dofmap_t
  use space, only : space_t, GLL
  use global_interpolation, only : global_interpolation_t
  use utils, only : neko_error, NEKO_FNAME_LEN, extract_fld_file_index
  use mesh, only : mesh_t
  implicit none
  private

  type, public :: fld_file_data_t
     type(vector_t) :: x !< x-coordinates.
     type(vector_t) :: y !< y-coordinates.
     type(vector_t) :: z !< z-coordinates.
     type(vector_t) :: u !< x-velocity field.
     type(vector_t) :: v !< y-velocity field.
     type(vector_t) :: w !< z-velocity field.
     type(vector_t) :: p !< Pressure field.
     type(vector_t) :: t !< Temperature field.
     integer, allocatable :: idx(:) !< element indices.
     type(vector_t), allocatable :: s(:) !< Numbered scalar fields.
     integer :: gdim !< Spatial dimensions.
     integer :: n_scalars = 0 !< Number of numbered scalar fields.
     real(kind=rp) :: time = 0.0 !< Time of latest sample/read.
     integer :: glb_nelv = 0 !< Global number of elements.
     integer :: nelv = 0 !< Number of elements on this rank.
     integer :: offset_el = 0 !< Element offset for this rank.
     integer :: lx = 0 !< Number of GLL points in x.
     integer :: ly = 0 !< Number of GLL points in y.
     integer :: lz = 0 !< Number of GLL points in z.
     integer :: t_counter = 0 !< counter of samples.
     ! meta file information (if any)
     !> Number of samples specified in .nek5000 file.
     integer :: meta_nsamples = 0
     integer :: meta_start_counter = 0 !< Index of the first field file.
     !> name of fld series as specified in .nek5000 (meta) file.
     character(len=1024) :: fld_series_fname
   contains
     procedure, pass(this) :: init => fld_file_data_init
     procedure, pass(this) :: free => fld_file_data_free
     procedure, pass(this) :: scale => fld_file_data_scale
     procedure, pass(this) :: add => fld_file_data_add
     procedure, pass(this) :: size => fld_file_data_size
     procedure, pass(this) :: get_list => fld_file_data_get_list
     procedure, pass(this) :: init_same => fld_file_data_init_same
     procedure, pass(this) :: init_n_fields => fld_file_data_init_n_fields
     !> Generates a global_interpolation object to interpolate the fld data.
     procedure, pass(this) :: generate_interpolator => &
          fld_file_data_generate_interpolator
     !> Reads an fld file and import fields, with/without interpolation.
     procedure, pass(this) :: import_fields => fld_file_data_import_fields
  end type fld_file_data_t

contains

  !> Imports fields from an fld_file_data object, potentially with
  !! interpolation.
  !! @param this fld_data object. Must already be initialized, no
  !! checks are done.
  !! @param u The field on which to import the u component of the fld data.
  !! @param v The field on which to import the v component of the fld data.
  !! @param w The field on which to import the w component of the fld data.
  !! @param p The field on which to import the pressure field of the fld data.
  !! @param t The field on which to import the temperature field of the fld
  !! data.
  !! @param s_target_list Field list containing the fields on which to import
  !! the scalar fields of the fld data. Unless a list of target indices is
  !! provided in `s_index_list`, assigns field at position `i` in the list
  !! to scalar `i` in the fld file.
  !! @param s_index_list The list of scalars indices from which to load the
  !! fields provided in `s_target_list`. Must have the same size as
  !! `s_target_list`. For example, s_index_list = (/2,3/) will load scalar #2
  !! in `s_target_list%items(1)` and scalar #3 in `s_target_list%items(2)`.
  !! Index  0 corresponds to temperature by default. Therefore using
  !! `s_index_list = (/0/)` is equivalent to using the argument `t=...`.
  !! @param interpolate Whether or not to interpolate the fld data.
  !! @param tolerance If interpolation is enabled, the tolerance to use for the
  !! point finding.
  !! @note If interpolation is disabled, space-to-space interpolation is still
  !! performed within each element to allow for seamless change of polynomial
  !! order for the same given mesh.
  !! @attention No data movement between CPU and GPU is done in this
  !! subroutine. The required data must be copied manually beforehand (see
  !! import_field_utils.f90).
  subroutine fld_file_data_import_fields(this, u, v, w, p, t, &
       s_target_list, s_index_list, interpolate, tolerance)
    class(fld_file_data_t), intent(inout) :: this
    type(field_t), pointer, intent(inout), optional :: u,v,w,p,t
    type(field_list_t), intent(inout), optional :: s_target_list
    integer, intent(in), optional :: s_index_list(:)
    logical, intent(in), optional :: interpolate
    real(kind=rp), intent(in), optional :: tolerance

    integer :: i

    ! ---- For the mesh to mesh interpolation
    logical :: mesh_mismatch
    logical :: interpolate_
    type(global_interpolation_t) :: global_interp
    type(dofmap_t), pointer :: dof
    type(mesh_t) , pointer :: msh
    ! -----

    ! ---- For space to space interpolation
    type(space_t) :: prev_Xh
    type(space_t) , pointer :: Xh
    type(interpolator_t) :: space_interp
    ! ----

    character(len=LOG_SIZE) :: log_buf

    ! ---- Default values
    interpolate_ = .false.
    if (present(interpolate)) interpolate_ = interpolate
    ! ----

    !
    ! Handle the passing of arguments and pointers
    !
    dof => null()
    msh => null()
    Xh => null()

    if (present(u)) then
       dof => u%dof; msh => u%msh; Xh => u%Xh
    else if (present(v)) then
       dof => v%dof; msh => v%msh; Xh => v%Xh
    else if (present(w)) then
       dof => w%dof; msh => w%msh; Xh => w%Xh
    else if (present(p)) then
       dof => p%dof; msh => p%msh; Xh => p%Xh
    else if (present(t)) then
       dof => t%dof; msh => t%msh; Xh => t%Xh
    else if (present(s_target_list)) then
       if (s_target_list%size() .eq. 0) then
          call neko_error("Scalar target list is empty")
       else
          dof => s_target_list%items(1)%ptr%dof
          msh => s_target_list%items(1)%ptr%msh
          Xh => s_target_list%items(1)%ptr%Xh
       end if
    else
       call neko_error("At least one field must be passed")
    end if

    !
    ! Check that the data in the fld file matches the current case.
    ! Note that this is a safeguard and there are corner cases where
    ! two different meshes have the same dimension and same # of elements
    ! but this should be enough to cover obvious cases.
    !
    mesh_mismatch = (this%glb_nelv .ne. msh%glb_nelv .or. &
         this%gdim .ne. msh%gdim)

    if (mesh_mismatch .and. .not. interpolate_) then
       call neko_error("The fld file must match the current mesh! &
       &Use 'interpolate': 'true' to enable interpolation.")
    else if (.not. mesh_mismatch .and. interpolate_) then
       call neko_log%warning("You have activated interpolation but you might &
       &still be using the same mesh.")
    end if

    !
    ! Handling of interpolation and I/O
    !
    if (interpolate_) then

       ! Throw error if dof or msh are not specified
       ! This should never be thrown, but just in case.
       if (.not. associated(dof) .or. .not. associated(msh)) then
          call neko_error("both dof and msh must be associated")
       end if

       ! Generates an interpolator object and performs the point search
       call this%generate_interpolator(global_interp, dof, msh, &
            tolerance = tolerance)

       ! Evaluate all the fields
       if (present(u)) call global_interp%evaluate(u%x(:,1,1,1), this%u%x, &
            on_host = .false.)
       if (present(v)) call global_interp%evaluate(v%x(:,1,1,1), this%v%x, &
            on_host = .false.)
       if (present(w)) call global_interp%evaluate(w%x(:,1,1,1), this%w%x, &
            on_host = .false.)
       if (present(p)) call global_interp%evaluate(p%x(:,1,1,1), this%p%x, &
            on_host = .false.)
       if (present(t)) call global_interp%evaluate(t%x(:,1,1,1), this%t%x, &
            on_host = .false.)
       if (present(s_target_list)) then

          ! If the index list exists, use it as a "mask"
          if (present(s_index_list)) then
             do i = 1, size(s_index_list)
                ! Take care that if we set i=0 we want temperature
                if (s_index_list(i) .eq. 0) then
                   call global_interp%evaluate(s_target_list%x(i), &
                        this%t%x, on_host = .false.)
                else
                   ! For scalar fields, require indices in 1:this%n_scalars
                   if (s_index_list(i) < 1 .or. &
                           s_index_list(i) > this%n_scalars) then
                      call neko_error("s_index_list entry out of bounds")
                   end if
                   call global_interp%evaluate(s_target_list%x(i), &
                        this%s(s_index_list(i))%x, on_host = .false.)
                end if
             end do

             ! otherwise, just copy element-to-element
          else
             do i = 1, s_target_list%size()
                call global_interp%evaluate(s_target_list%x(i), this%s(i)%x, &
                     on_host = .false.)
             end do
          end if ! present s_index_list
       end if ! present s_tgt

       call global_interp%free()

    else ! No interpolation, but potentially just from different spaces

       ! throw an error is the space is not passed
       if (.not. associated(Xh)) call neko_error("Xh is not associated")

       ! Build a space_t object from the data in the fld file
       call prev_Xh%init(GLL, this%lx, this%ly, this%lz)
       call space_interp%init(Xh, prev_Xh)

       ! Do the space-to-space interpolation
       if (present(u)) call space_interp%map(u%x, this%u%x, this%nelv, Xh)
       if (present(v)) call space_interp%map(v%x, this%v%x, this%nelv, Xh)
       if (present(w)) call space_interp%map(w%x, this%w%x, this%nelv, Xh)
       if (present(p)) call space_interp%map(p%x, this%p%x, this%nelv, Xh)
       if (present(t)) call space_interp%map(t%x, this%t%x, this%nelv, Xh)
       if (present(s_target_list)) then

          ! If the index list exists, use it as a "mask"
          if (present(s_index_list)) then
             do i = 1, size(s_index_list)

                ! 0 means we want temperature
                if (s_index_list(i) .eq. 0) then
                   call space_interp%map(s_target_list%x(i), &
                        this%t%x, this%nelv, Xh)
                else
                   call space_interp%map(s_target_list%x(i), &
                        this%s(s_index_list(i))%x, this%nelv, Xh)
                end if
             end do

             ! otherwise, just copy element-to-element
          else
             do i = 1, s_target_list%size()
                call space_interp%map(s_target_list%x(i), this%s(i)%x, &
                     this%nelv, Xh)
             end do
          end if ! present s_index_list
       end if ! present s_tgt

       call space_interp%free
       call prev_Xh%free

    end if ! if interpolate

    nullify(dof)
    nullify(Xh)
    nullify(msh)

  end subroutine fld_file_data_import_fields

  !> Initializes a fld_file_data object.
  !! @param nelv Number of elements (on this rank).
  !! @param offset_el Element offset for this rank.
  subroutine fld_file_data_init(this, nelv, offset_el)
    class(fld_file_data_t), intent(inout) :: this
    integer, intent(in), optional :: nelv, offset_el

    call this%free()
    if (present(nelv)) this%nelv = nelv
    if (present(offset_el)) this%offset_el = offset_el

  end subroutine fld_file_data_init

  !> Get the number of initialized fields in this fld file.
  function fld_file_data_size(this) result(i)
    class(fld_file_data_t) :: this
    integer :: i
    i = 0
    if (this%u%size() .gt. 0) i = i + 1
    if (this%v%size() .gt. 0) i = i + 1
    if (this%w%size() .gt. 0) i = i + 1
    if (this%p%size() .gt. 0) i = i + 1
    if (this%t%size() .gt. 0) i = i + 1
    i = i + this%n_scalars

  end function fld_file_data_size

  !> Genereate same fields as in another fld_file
  subroutine fld_file_data_init_same(this, fld_file, n)
    class(fld_file_data_t), target, intent(inout) :: this
    class(fld_file_data_t), target, intent(in) :: fld_file
    integer, intent(in) :: n
    integer :: i, j

    if (fld_file%u%size() .gt. 0) then
       call this%u%init(n)
    end if
    if (fld_file%v%size() .gt. 0) then
       call this%v%init(n)
    end if
    if (fld_file%w%size() .gt. 0) then
       call this%w%init(n)
    end if
    if (fld_file%p%size() .gt. 0) then
       call this%p%init(n)
    end if
    if (fld_file%t%size() .gt. 0) then
       call this%t%init(n)
    end if
    this%n_scalars = fld_file%n_scalars
    allocate(this%s(fld_file%n_scalars))
    do j = 1, fld_file%n_scalars
       call this%s(j)%init(n)
    end do

  end subroutine fld_file_data_init_same

  !> Generate same fields as in another fld_file
  subroutine fld_file_data_init_n_fields(this, n_fields, n)
    class(fld_file_data_t), target, intent(inout) :: this
    integer, intent(in) :: n, n_fields
    integer :: i, j


    if (n_fields .gt. 0) then
       call this%u%init(n)
    end if
    if (n_fields .gt. 1) then
       call this%v%init(n)
    end if
    if (n_fields .gt. 2) then
       call this%w%init(n)
    end if
    if (n_fields .gt. 3) then
       call this%p%init(n)
    end if
    if (n_fields .gt. 4) then
       call this%t%init(n)
    end if
    if (n_fields .gt. 5) then
       this%n_scalars = n_fields-5
       allocate(this%s(this%n_scalars))
       do j = 1, this%n_scalars
          call this%s(j)%init(n)
       end do
    end if

  end subroutine fld_file_data_init_n_fields

  !> Get a list with pointers to the fields in the fld file
  subroutine fld_file_data_get_list(this, ptr_list, n)
    class(fld_file_data_t), target, intent(in) :: this
    integer, intent(in) :: n
    integer :: i, j
    type(vector_ptr_t), intent(inout) :: ptr_list(n)
    i = 1
    if (this%u%size() .gt. 0) then
       ptr_list(i)%ptr => this%u
       i = i + 1
    end if
    if (this%v%size() .gt. 0) then
       ptr_list(i)%ptr => this%v
       i = i + 1
    end if
    if (this%w%size() .gt. 0) then
       ptr_list(i)%ptr => this%w
       i = i + 1
    end if
    if (this%p%size() .gt. 0) then
       ptr_list(i)%ptr => this%p
       i = i + 1
    end if
    if (this%t%size() .gt. 0) then
       ptr_list(i)%ptr => this%t
       i = i + 1
    end if
    do j = 1, this%n_scalars
       ptr_list(i)%ptr => this%s(j)
       i = i +1
    end do

  end subroutine fld_file_data_get_list



  !> Scale the values stored in this fld_file_data
  subroutine fld_file_data_scale(this, c)
    class(fld_file_data_t), intent(inout) :: this
    real(kind=rp), intent(in) :: c
    integer :: i

    if (this%u%size() .gt. 0) call cmult(this%u%x, c, this%u%size())
    if (this%v%size() .gt. 0) call cmult(this%v%x, c, this%v%size())
    if (this%w%size() .gt. 0) call cmult(this%w%x, c, this%w%size())
    if (this%p%size() .gt. 0) call cmult(this%p%x, c, this%p%size())
    if (this%t%size() .gt. 0) call cmult(this%t%x, c, this%t%size())

    do i = 1, this%n_scalars
       if (this%s(i)%size() .gt. 0) call cmult(this%s(i)%x, c, this%s(i)%size())
    end do

  end subroutine fld_file_data_scale

  !> Add the values in another fld file to this
  subroutine fld_file_data_add(this, other)
    class(fld_file_data_t), intent(inout) :: this
    class(fld_file_data_t), intent(in) :: other
    integer :: i, n

    if (this%u%size() .gt. 0) call add2(this%u%x, other%u%x, this%u%size())
    if (this%v%size() .gt. 0) call add2(this%v%x, other%v%x, this%v%size())
    if (this%w%size() .gt. 0) call add2(this%w%x, other%w%x, this%w%size())
    if (this%p%size() .gt. 0) call add2(this%p%x, other%p%x, this%p%size())
    if (this%t%size() .gt. 0) call add2(this%t%x, other%t%x, this%t%size())

    do i = 1, this%n_scalars
       if (this%s(i)%size() .gt. 0) call add2(this%s(i)%x, other%s(i)%x, &
            this%s(i)%size())
    end do
  end subroutine fld_file_data_add

  !> Deallocate fld file data type
  subroutine fld_file_data_free(this)
    class(fld_file_data_t), intent(inout) :: this
    integer :: i
    call this%x%free()
    call this%y%free()
    call this%z%free()
    call this%u%free()
    call this%v%free()
    call this%w%free()
    call this%p%free()
    call this%t%free()
    if (allocated(this%s)) then
       do i = 1, this%n_scalars
          call this%s(i)%free()
       end do
       deallocate(this%s)
    end if
    this%n_scalars = 0
    this%time = 0.0
    this%glb_nelv = 0
    this%nelv = 0
    this%offset_el = 0
    this%lx = 0
    this%ly = 0
    this%lz = 0
    this%t_counter = 0
    this%meta_nsamples = 0
    this%meta_start_counter = 0
    if (allocated(this%idx)) deallocate(this%idx)
  end subroutine fld_file_data_free

  !> Generates a global_interpolation object to interpolate the fld data.
  !! @param global_interp Global interpolation object which will be initialized
  !! and ready to be used for interpolation.
  !! @param to_dof Dofmap on which to interpolate.
  !! @param to_msh Mesh on which to interpolate.
  !! @param tolerance Tolerance for the newton iterations.
  subroutine fld_file_data_generate_interpolator(this, global_interp, to_dof, &
       to_msh, tolerance)
    class(fld_file_data_t), intent(in) :: this
    type(global_interpolation_t), intent(inout) :: global_interp
    type(dofmap_t), intent(in), target :: to_dof
    type(mesh_t), intent(in), target :: to_msh
    real(kind=rp), intent(in) :: tolerance

    ! --- variables for interpolation
    type(space_t) :: fld_Xh
    real(kind=rp), allocatable :: x_coords(:,:,:,:), y_coords(:,:,:,:), &
         z_coords(:,:,:,:)
    real(kind=rp) :: center_x, center_y, center_z
    integer :: e, i
    ! ---

    type(space_t), pointer :: to_Xh
    to_Xh => to_dof%Xh

    ! Safeguard in case we didn't read mesh information
    if (.not. allocated(this%x%x) .or. &
         .not. allocated(this%y%x) .or. &
         .not. allocated(this%z%x)) call neko_error("Unable to retrieve &
    &mesh information from fld data.")

    ! Create a space based on the fld data
    call fld_Xh%init(GLL, this%lx, this%ly, this%lz)

    ! These are the coordinates of our current dofmap
    ! that we use for the interpolation
    allocate(x_coords(to_Xh%lx, to_Xh%ly, to_Xh%lz, to_msh%nelv))
    allocate(y_coords(to_Xh%lx, to_Xh%ly, to_Xh%lz, to_msh%nelv))
    allocate(z_coords(to_Xh%lx, to_Xh%ly, to_Xh%lz, to_msh%nelv))

    !> To ensure that each point is within an element
    !! Remedies issue with points on the boundary
    !! Technically gives each point a slightly different value
    !! but still within the specified tolerance
    do e = 1, to_msh%nelv
       center_x = 0d0
       center_y = 0d0
       center_z = 0d0
       do i = 1, to_Xh%lxyz
          center_x = center_x + to_dof%x(i, 1, 1, e)
          center_y = center_y + to_dof%y(i, 1, 1, e)
          center_z = center_z + to_dof%z(i, 1, 1, e)
       end do
       center_x = center_x / to_Xh%lxyz
       center_y = center_y / to_Xh%lxyz
       center_z = center_z / to_Xh%lxyz
       do i = 1, to_Xh%lxyz
          x_coords(i, 1, 1, e) = to_dof%x(i, 1, 1, e) - &
               tolerance * (to_dof%x(i, 1, 1, e) - center_x)
          y_coords(i, 1, 1, e) = to_dof%y(i, 1, 1, e) - &
               tolerance * (to_dof%y(i, 1, 1, e) - center_y)
          z_coords(i, 1, 1, e) = to_dof%z(i, 1, 1, e) - &
               tolerance * (to_dof%z(i, 1, 1, e) - center_z)
       end do
    end do

    ! The initialization is done based on the variables created from
    ! fld data
    call global_interp%init(this%x%x, this%y%x, this%z%x, this%gdim, &
         this%nelv, fld_Xh, tol = tolerance)
    call global_interp%find_points(x_coords, y_coords, z_coords, &
         to_dof%size())

    deallocate(x_coords)
    deallocate(y_coords)
    deallocate(z_coords)
    call fld_Xh%free()

  end subroutine fld_file_data_generate_interpolator

end module fld_file_data
