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
 !> Abstract type for mesh element class
module polytope_mesh
  use num_types, only : i4, dp
  use utils, only : neko_error
  use polytope, only : polytope_t
  use polytope_actualisation, only : polytope_actualisation_t
  use point, only : point_ptr, point_t
  implicit none
  private

  public :: mesh_object_t, polytope_mesh_t, mesh_element_t

  !> Single element object allocatable space
  type :: mesh_object_t
     class(polytope_actualisation_t), allocatable :: obj
  end type mesh_object_t

  !> Base type for a mesh polytope.
  !! @details This is an abstract type building on actualisation and point data.
  !! It contains both topology and geometry information and provides the highest
  !! dimension objects in the mesh called elements in a spectral element method.
  type, extends(polytope_t), abstract :: polytope_mesh_t
     !> Geometrical dimension
     integer(i4), private :: gdim_ = -1
     !> Polytope facets
     type(mesh_object_t), dimension(:), allocatable :: facet
     !> Polytope ridges
     type(mesh_object_t), dimension(:), allocatable :: ridge
     !> Polytope peaks
     type(mesh_object_t), dimension(:), allocatable ::  peak
     !> Number of points
     integer(i4), private :: npts_
     !> Element vertices position
     type(point_ptr), dimension(:), allocatable :: pts
   contains
     !> Free aligned polytope and interpolation data
     procedure, pass(this) :: free => polytope_free
     !> Initialise general data
     procedure, pass(this) :: init_dat => polytope_init_data
     !> Return geometrical dimension
     procedure, pass(this) :: gdim => polytope_gdim_get
     !> Return number of points
     procedure, pass(this) :: npts => polytope_npts_get
     !> Return a pointer to the polytope facets (faces or edges)
     procedure, pass(this) :: fct => polytope_fct_ptr
     !> Return a pointer to the polytope ridges (vertices or edges)
     procedure, pass(this) :: rdg => polytope_rdg_ptr
     !> Return a pointer to the polytope peaks (vertices)
     procedure, pass(this) :: pek => polytope_pek_ptr
     !> Return a pointer to the polytope vertex points
     procedure, pass(this) :: pnt => polytope_pnt_ptr
     !> Is polytope self-periodic?
     procedure, pass(this) :: selfp => polytope_self_periodic
     !> Return facets shared by polytopes
     procedure, pass(this) :: fct_share => polytope_facet_share
     !> Return ridges shared by polytopes
     procedure, pass(this) :: rdg_share => polytope_ridge_share
     !> Return peaks shared by polytopes
     procedure, pass(this) :: pek_share => polytope_peak_share
     !> Return boundary flag for a given facet
     procedure, pass(this) :: fct_bnd => polytope_fct_bnd
     !> Return communication id for a given facet
     procedure, pass(this) :: fct_gsid => polytope_fct_gsid
     !> Return communication id for a given ridge
     procedure, pass(this) :: rdg_gsid => polytope_rdg_gsid
     !> Return communication id for a given peak
     procedure, pass(this) :: pek_gsid => polytope_pek_gsid
     !> Return facet alignment
     procedure, pass(this) :: falgn => polytope_fct_algn
     !> Just required
     procedure, pass(this) :: bnd => polytope_int_dummy
     !> Just required
     procedure, pass(this) :: gsid => polytope_int_dummy
     !> Initialise a topology polytope
     procedure(polytope_mesh_init), pass(this), deferred :: init
     !> Return element diameter
     procedure(polytope_mesh_diameter), pass(this), deferred :: diameter
     !> Return element centroid
     procedure(polytope_mesh_centroid), pass(this), deferred :: centroid
     !> Return facet @a r and @s local directions with respect to the element
     procedure(polytope_fct_dir), pass(this), deferred :: fct_dir
     !> Return ridge @a r local direction with respect to the element
     procedure(polytope_rdg_dir), pass(this), deferred :: rdg_dir
  end type polytope_mesh_t

  !> Single mesh element allocatable space
  type :: mesh_element_t
     class(polytope_mesh_t), allocatable :: el
  end type mesh_element_t

  !> Abstract interface to initialise a polytope with geometry information
  !! @parameter[in]   id       polytope id
  !! @parameter[in]   nfct     number of facets
  !! @parameter[in]   fct      polytope facets
  !! @parameter[in]   npts     number of points
  !! @parameter[in]   pts      points
  !! @parameter[in]   gdim     geometrical dimension
  !! @parameter[in]   nrdg     number of hanging ridges
  !! @parameter[in]   rdg_hng  ridge hanging flag
  abstract interface
     subroutine polytope_mesh_init(this, id, nfct, fct, npts, pts, gdim, nrdg, &
          & rdg_hng)
       import i4
       import polytope_mesh_t
       import mesh_object_t
       import point_ptr
       class(polytope_mesh_t), intent(inout) :: this
       integer(i4), intent(in) :: id, nfct, npts, gdim, nrdg
       type(mesh_object_t), dimension(nfct), intent(inout) :: fct
       type(point_ptr), dimension(npts), intent(in) :: pts
       integer(i4), dimension(2, 3), intent(in) :: rdg_hng
     end subroutine polytope_mesh_init
  end interface

  !> Abstract interface to get element diameter
  !! @return res
  abstract interface
     function polytope_mesh_diameter(this) result(res)
       import dp
       import polytope_mesh_t
       class(polytope_mesh_t), intent(in) :: this
       real(dp) :: res
     end function polytope_mesh_diameter
  end interface

  !> Abstract interface to get element centroid
  !! @return res
  abstract interface
     function polytope_mesh_centroid(this) result(res)
       import polytope_mesh_t
       import point_t
       class(polytope_mesh_t), intent(in) :: this
       type(point_t) :: res
     end function polytope_mesh_centroid
  end interface

  !> Abstract interface to get @a r and @a s facet local directions
  !! @parameter[in]   pos          facet position
  !! @parameter[out]  dirr, dirs   local directions
  abstract interface
     subroutine polytope_fct_dir(this, pos, dirr, dirs)
       import i4
       import polytope_mesh_t
       class(polytope_mesh_t), intent(in) :: this
       integer(i4), intent(in) :: pos
       integer(i4), intent(out) :: dirr, dirs
     end subroutine polytope_fct_dir
  end interface

  !> Abstract interface to get @a r ridge local direction
  !! @parameter[in]   pos          ridge position
  !! @parameter[out]  dirr         local direction
  abstract interface
     subroutine polytope_rdg_dir(this, pos, dirr)
       import i4
       import polytope_mesh_t
       class(polytope_mesh_t), intent(in) :: this
       integer(i4), intent(in) :: pos
       integer(i4), intent(out) :: dirr
     end subroutine polytope_rdg_dir
  end interface

contains

  !> Free polytope data
  subroutine polytope_free(this)
    class(polytope_mesh_t), intent(inout) :: this
    integer(i4) :: il
    this%gdim_ = -1
    if (allocated(this%facet)) then
       do il = 1, this%nfacet
          call this%facet(il)%obj%free()
          deallocate(this%facet(il)%obj)
       end do
       deallocate(this%facet)
    end if
    if (allocated(this%ridge)) then
       do il = 1, this%nridge
          call this%ridge(il)%obj%free()
          deallocate(this%ridge(il)%obj)
       end do
       deallocate(this%ridge)
    end if
    if (allocated(this%peak)) then
       do il = 1, this%npeak
          call this%peak(il)%obj%free()
          deallocate(this%peak(il)%obj)
       end do
       deallocate(this%peak)
    end if
    call this%freep()
    if (allocated(this%pts)) then
       do il = 1, this%npts_
          this%pts(il)%p => null()
       end do
       deallocate(this%pts)
    end if
    this%npts_ = -1
  end subroutine polytope_free

  !> Initialise general data
  !! @parameter[in]   pltp   polytope
  !! @parameter[in]   ifalgn if non identity alignment
  !! @parameter[in]   ifint  interpolation flag
  !! @parameter[in]   hng    hanging information
  !! @parameter[in]   pos    position in the higher order element
  subroutine polytope_init_data(this, gdim, npts)
    class(polytope_mesh_t), intent(inout) :: this
    integer(i4), intent(in) :: gdim, npts

    this%gdim_ = gdim
    this%npts_ = npts

  end subroutine polytope_init_data

  !> @brief Get geometrical dimension
  !! @return   gdim
  pure function polytope_gdim_get(this) result(gdim)
    class(polytope_mesh_t), intent(in) :: this
    integer(i4) :: gdim
    gdim = this%gdim_
  end function polytope_gdim_get

  !> @brief Get number of points
  !! @return   npts
  pure function polytope_npts_get(this) result(npts)
    class(polytope_mesh_t), intent(in) :: this
    integer(i4) :: npts
    npts = this%npts_
  end function polytope_npts_get

  !> @brief Return pointer to the polytope facet
  !! @parameter[in]   pos   polytope element position
  !! @return ptr
  function polytope_fct_ptr(this, pos) result(ptr)
    class(polytope_mesh_t), intent(in) :: this
    integer(i4), intent(in) :: pos
    class(polytope_t), pointer :: ptr
    if ((pos > 0) .and. (pos <= this%nfacet)) then
       ptr => this%facet(pos)%obj%polytope
    else
       call neko_error('Wrong facet number for mesh objects.')
    end if
  end function polytope_fct_ptr

  !> @brief Return pointer to the polytope ridge
  !! @parameter[in]   pos   polytope element position
  !! @return ptr
  function polytope_rdg_ptr(this, pos) result(ptr)
    class(polytope_mesh_t), intent(in) :: this
    integer(i4), intent(in) :: pos
    class(polytope_t), pointer :: ptr
    if ((pos > 0) .and. (pos <= this%nridge)) then
       ptr => this%ridge(pos)%obj%polytope
    else
       call neko_error('Wrong ridge number for mesh objects.')
    end if
  end function polytope_rdg_ptr

  !> @brief Return pointer to the polytope peak
  !! @parameter[in]   pos   polytope element position
  !! @return ptr
  function polytope_pek_ptr(this, pos) result(ptr)
    class(polytope_mesh_t), intent(in) :: this
    integer(i4), intent(in) :: pos
    class(polytope_t), pointer :: ptr
    if ((pos > 0) .and. (pos <= this%npeak)) then
       ptr => this%peak(pos)%obj%polytope
    else
       call neko_error('Wrong peak number for mesh objects.')
    end if
  end function polytope_pek_ptr

  !> @brief Return pointer to the polytope vertex point
  !! @parameter[in]   pos   polytope point position
  !! @return ptr
  function polytope_pnt_ptr(this, pos) result(ptr)
    class(polytope_mesh_t), intent(in) :: this
    integer(i4), intent(in) :: pos
    type(point_t), pointer :: ptr
    if ((pos > 0) .and. (pos <= this%npts_)) then
       ptr => this%pts(pos)%p
    else
       call neko_error('Wrong point number for mesh objects.')
    end if
  end function polytope_pnt_ptr

  !> @brief Check if polytope is self-periodic
  !! @return   selfp
  function polytope_self_periodic(this) result(selfp)
    class(polytope_mesh_t), intent(in) :: this
    logical :: selfp
    integer(i4) :: il, jl, itmp, algn

    ! count self periodic facets (faces or edges)
    itmp = 0
    do il = 1, this%nfacet - 1
       do jl = il + 1, this%nfacet
          selfp = this%facet(il)%obj%equal(this%facet(jl)%obj%polytope)
          if (selfp) itmp = itmp + 1
       end do
    end do
    ! count self periodic ridges (edges or vertices)
    itmp = 0
    do il = 1, this%nridge - 1
       do jl = il + 1, this%nridge
          selfp = this%ridge(il)%obj%equal(this%ridge(jl)%obj%polytope)
          if (selfp) itmp = itmp + 1
       end do
    end do
    ! count self periodic peaks (vertices)
    do il = 1, this%npeak - 1
       do jl = il + 1, this%npeak
          selfp = (this%peak(il)%obj%polytope%id() == &
               & this%peak(jl)%obj%polytope%id())
          if (selfp) itmp = itmp + 1
       end do
    end do
    if (itmp == 0) then
       selfp = .false.
    else
       selfp = .true.
    end if
  end function polytope_self_periodic

  !> @brief Return positions of facets shared by polytopes
  !! @note Polytopes can be self-periodic
  !! @parameter[in]   other   second polytope
  !! @parameter[out]  ishare  number of shared facets
  !! @parameter[out]  facetp  integer position of shared facets
  subroutine polytope_facet_share(this, other, ishare, facetp)
    class(polytope_mesh_t), intent(in) :: this, other
    integer(i4), intent(out) :: ishare
    integer(i4), dimension(:, :), allocatable, intent(out) :: facetp
    integer(i4) :: il, jl

    allocate(facetp(2, this%nfacet * other%nfacet))
    ishare = 0
    facetp(:, :) = 0
    do il = 1, this%nfacet
       do jl = 1, other%nfacet
          if (this%facet(il)%obj%equal(other%facet(jl)%obj%polytope)) then
             ishare = ishare + 1
             facetp(1, ishare) = il
             facetp(2, ishare) = jl
          end if
       end do
    end do
  end subroutine polytope_facet_share

  !> @brief Return positions of ridges shared by polytopes
  !! @note Polytopes can be self-periodic
  !! @parameter[in]   other   second polytope
  !! @parameter[out]  ishare  number of shared facets
  !! @parameter[out]  facetp  integer position of shared facets
  subroutine polytope_ridge_share(this, other, ishare, facetp)
    class(polytope_mesh_t), intent(in) :: this, other
    integer(i4), intent(out) :: ishare
    integer(i4), dimension(:, :), allocatable, intent(out) :: facetp
    integer(i4) :: il, jl

    allocate(facetp(2, this%nridge * other%nridge))
    ishare = 0
    facetp(:, :) = 0
    do il = 1, this%nridge
       do jl = 1, other%nridge
          if (this%ridge(il)%obj%equal(other%ridge(jl)%obj%polytope)) then
             ishare = ishare + 1
             facetp(1, ishare) = il
             facetp(2, ishare) = jl
          end if
       end do
    end do
  end subroutine polytope_ridge_share

  !> @brief Return positions of peaks (vertices) shared by polytopes
  !! @note Plytopes can be self-periodic
  !! @parameter[in]   other   second polytope
  !! @parameter[out]  ishare  number of shared vertices
  !! @parameter[out]  ridgep  integer position of shared vertices
  pure subroutine polytope_peak_share(this, other, ishare, ridgep)
    class(polytope_mesh_t), intent(in) :: this, other
    integer(i4), intent(out) :: ishare
    integer(i4), dimension(:, :), allocatable, intent(out) :: ridgep
    integer(i4) :: il, jl

    allocate(ridgep(2, this%npeak * other%npeak))
    ishare = 0
    ridgep(:, :) = 0
    do il = 1, this%npeak
       do jl = 1, other%npeak
          if (this%peak(il)%obj%polytope%id() == &
               & other%peak(jl)%obj%polytope%id()) then
             ishare = ishare + 1
             ridgep(1, ishare) = il
             ridgep(2, ishare) = jl
          end if
       end do
    end do
  end subroutine polytope_peak_share

  !> @brief Return boundary flag of the polytope facet
  !! @parameter[in]   pos   polytope facet position
  !! @return bnd
  function polytope_fct_bnd(this, pos) result(bnd)
    class(polytope_mesh_t), intent(in) :: this
    integer(i4), intent(in) :: pos
    integer(i4) :: bnd
    if ((pos > 0) .and. (pos <= this%nfacet)) then
       bnd = this%facet(pos)%obj%polytope%bnd()
    else
       call neko_error('Wrong facet number for mesh objects boundary.')
    end if
  end function polytope_fct_bnd

  !> @brief Return communication id of the polytope facet
  !! @parameter[in]   pos   polytope facet position
  !! @return gsid
  function polytope_fct_gsid(this, pos) result(gsid)
    class(polytope_mesh_t), intent(in) :: this
    integer(i4), intent(in) :: pos
    integer(i4) :: gsid
    if ((pos > 0) .and. (pos <= this%nfacet)) then
       gsid = this%facet(pos)%obj%polytope%gsid()
    else
       call neko_error('Wrong facet number for mesh objects gsid.')
    end if
  end function polytope_fct_gsid

  !> @brief Return communication id of the polytope ridge
  !! @parameter[in]   pos   polytope ridge position
  !! @return gsid
  function polytope_rdg_gsid(this, pos) result(gsid)
    class(polytope_mesh_t), intent(in) :: this
    integer(i4), intent(in) :: pos
    integer(i4) :: gsid
    if ((pos > 0) .and. (pos <= this%nridge)) then
       gsid = this%ridge(pos)%obj%polytope%gsid()
    else
       call neko_error('Wrong ridge number for mesh objects gsid.')
    end if
  end function polytope_rdg_gsid

  !> @brief Return communication id of the polytope peak
  !! @parameter[in]   pos   polytope peak position
  !! @return gsid
  function polytope_pek_gsid(this, pos) result(gsid)
    class(polytope_mesh_t), intent(in) :: this
    integer(i4), intent(in) :: pos
    integer(i4) :: gsid
    if ((pos > 0) .and. (pos <= this%npeak)) then
       gsid = this%peak(pos)%obj%polytope%gsid()
    else
       call neko_error('Wrong peak number for mesh objects gsid.')
    end if
  end function polytope_pek_gsid

  !> Return facet alignment
  !! @return algn
  function polytope_fct_algn(this, pos) result(algn)
    class(polytope_mesh_t), intent(in) :: this
    integer(i4), intent(in) :: pos
    integer(i4) :: algn
    algn = this%facet(pos)%obj%algn_op%algn()
  end function polytope_fct_algn

  !> Dummy routine pretending to extract integer property
  !! @return intp
  function polytope_int_dummy(this) result(intp)
    class(polytope_mesh_t), intent(in) :: this
    integer(i4) :: intp
    intp = - 1
    call neko_error('This property should not be used for mesh objects.')
  end function polytope_int_dummy

end module polytope_mesh
