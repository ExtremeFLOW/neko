! Copyright (c) 2020-2021, The Neko Authors
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
!> Dirichlet condition on axis aligned plane in the non normal direction
module non_normal
  use json_module, only : json_file
  use symmetry, only : symmetry_t
  use num_types, only : rp
  use tuple, only : tuple_i4_t
  use coefs, only : coef_t
  use utils, only : neko_error
  use neko_config, only : NEKO_BCKND_DEVICE
  use logger, only : neko_log, LOG_SIZE, NEKO_LOG_VERBOSE
  use amr_reconstruct, only : amr_reconstruct_t
  use, intrinsic :: iso_c_binding
  implicit none
  private

  !> Dirichlet condition in non normal direction of a plane.
  !! @warning Only works for axis-aligned plane boundaries.
  type, public, extends(symmetry_t) :: non_normal_t
   contains
     !> Constructor.
     procedure, pass(this) :: init => non_normal_init
     !> Constructor from components
     procedure, pass(this) :: init_from_components => &
          non_normal_init_from_components
     !> Destructor.
     procedure, pass(this) :: free => non_normal_free
     !> Finalize.
     procedure, pass(this) :: finalize => non_normal_finalize
     !> AMR restart
     procedure, pass(this) :: amr_restart => non_normal_amr_restart
  end type non_normal_t

contains

  !> Constructor
  !! @param[in] coef The SEM coefficients.
  !! @param[inout] json The JSON object configuring the boundary condition.
  subroutine non_normal_init(this, coef, json)
    class(non_normal_t), target, intent(inout) :: this
    type(coef_t), target, intent(in) :: coef
    type(json_file), intent(inout) ::json

    call this%init_from_components(coef)
  end subroutine non_normal_init

  !> Constructor from components.
  !! @param[in] coef The SEM coefficients.
  subroutine non_normal_init_from_components(this, coef)
    class(non_normal_t), target, intent(inout) :: this
    type(coef_t), target, intent(in) :: coef

    call this%free()
    call this%symmetry_t%init_from_components(coef)

  end subroutine non_normal_init_from_components

  !> Finalize
  subroutine non_normal_finalize(this, only_facets)
    class(non_normal_t), target, intent(inout) :: this
    logical, optional, intent(in) :: only_facets
    logical :: only_facets_ = .false.
    integer :: i, j, k, l
    type(tuple_i4_t), pointer :: bfp(:)
    real(kind=rp) :: sx, sy, sz
    real(kind=rp), parameter :: TOL = 1d-3
    type(tuple_i4_t) :: bc_facet
    integer :: facet, el

    if (present(only_facets)) then
       only_facets_ = only_facets
    else
       only_facets_ = this%only_facets
    end if

    associate(c => this%coef, nx => this%coef%nx, ny => this%coef%ny, &
         nz => this%coef%nz)
      bfp => this%marked_facet%array()
      do i = 1, this%marked_facet%size()
         bc_facet = bfp(i)
         facet = bc_facet%x(1)
         el = bc_facet%x(2)
         call this%get_normal_axis(sx, sy, sz, facet, el)

         if (sx .lt. TOL) then
            call this%bc_y%mark_facet(facet, el)
            call this%bc_z%mark_facet(facet, el)
         end if

         if (sy .lt. TOL) then
            call this%bc_x%mark_facet(facet, el)
            call this%bc_z%mark_facet(facet, el)
         end if

         if (sz .lt. TOL) then
            call this%bc_y%mark_facet(facet, el)
            call this%bc_x%mark_facet(facet, el)
         end if
      end do
    end associate
    call this%bc_x%finalize(only_facets_)
    call this%bc_y%finalize(only_facets_)
    call this%bc_z%finalize(only_facets_)

    call this%finalize_base(only_facets_)
  end subroutine non_normal_finalize

  !> Destructor
  subroutine non_normal_free(this)
    class(non_normal_t), target, intent(inout) :: this

    call this%symmetry_t%free()
  end subroutine non_normal_free

  !> AMR restart
  !! @param[inout]  reconstruct   data reconstruction type
  !! @param[in]     counter       restart counter
  !! @param[in]     tstep         time step
  subroutine non_normal_amr_restart(this, reconstruct, counter, tstep)
    class(non_normal_t), intent(inout) :: this
    type(amr_reconstruct_t), intent(inout) :: reconstruct
    integer, intent(in) :: counter, tstep
    character(len=LOG_SIZE) :: log_buf
    integer :: il

    ! Was this component already restarted?
    if (this%counter .eq. counter) return

    this%counter = counter

    ! For defined zone indices perform full reconstruction including
    ! finalisation. If zones are missing just prepare for collecting
    ! facets. Do not forget to finalise those bc later.
    if (allocated(this%zone_indices)) then
       if (allocated(this%type)) then
          log_buf = 'Reconstruct Non Normal: '//trim(this%type)
       else
          log_buf = 'Reconstruct Non Normal'
       end if
       call neko_log%section(log_buf, NEKO_LOG_VERBOSE)

       ! reconstruct dofmap; No problem, as AMR restart prevents recursive
       ! reconstructions
       if (associated(this%dof)) call this%dof%amr_restart(reconstruct, &
            counter, tstep)
       ! reconstruct coef; No problem, as AMR restart prevents recursive
       ! reconstructions
       if (associated(this%coef)) call this%coef%amr_restart(reconstruct, &
            counter, tstep)

       if (NEKO_BCKND_DEVICE .eq. 1) then
          ! added utils module; could be removed
          call neko_error('Non normal:: Nothing done for device.')
       end if

       ! free space
       if (allocated(this%msk)) deallocate(this%msk)
       if (allocated(this%facet)) deallocate(this%facet)
!       call this%marked_facet%free()
!       call this%marked_facet%init()
       call this%marked_facet%clear()

       ! Clean all bc components
       ! there should be no allocated zones
       if (allocated(this%bc_x%zone_indices) .or. &
            allocated(this%bc_y%zone_indices) .or. &
            allocated(this%bc_z%zone_indices)) &
            call neko_error('Non normal:: component zone_indices allocated')
       call this%bc_x%amr_restart(reconstruct, counter, tstep)
       call this%bc_y%amr_restart(reconstruct, counter, tstep)
       call this%bc_z%amr_restart(reconstruct, counter, tstep)

       this%iffinalised = .false.

       ! get zones
       do il = 1, size(this%zone_indices)
          call this%mark_zone(this%coef%msh%labeled_zones(&
               this%zone_indices(il)))
       end do
       call this%finalize()

       call neko_log%end_section(lvl = NEKO_LOG_VERBOSE)
    else
       if (allocated(this%type)) then
          log_buf = 'Clean Non Normal: '//trim(this%type)
       else
          log_buf = 'Clean Non Normal'
       end if
       call neko_log%section(log_buf, NEKO_LOG_VERBOSE)

       ! reconstruct dofmap; No problem, as AMR restart prevents recursive
       ! reconstructions
       if (associated(this%dof)) call this%dof%amr_restart(reconstruct, &
            counter, tstep)
       ! reconstruct coef; No problem, as AMR restart prevents recursive
       ! reconstructions
       if (associated(this%coef)) call this%coef%amr_restart(reconstruct, &
            counter, tstep)

       if (NEKO_BCKND_DEVICE .eq. 1) then
          ! added utils module; could be removed
          call neko_error('Non normal:: Nothing done for device.')
       end if

       ! free space
       if (allocated(this%msk)) deallocate(this%msk)
       if (allocated(this%facet)) deallocate(this%facet)
!       call this%marked_facet%free()
!       call this%marked_facet%init()
       call this%marked_facet%clear()

       ! Clean all bc components
       ! there should be no allocated zones
       if (allocated(this%bc_x%zone_indices) .or. &
            allocated(this%bc_y%zone_indices) .or. &
            allocated(this%bc_z%zone_indices)) &
            call neko_error('Non normal:: component zone_indices allocated')
       call this%bc_x%amr_restart(reconstruct, counter, tstep)
       call this%bc_y%amr_restart(reconstruct, counter, tstep)
       call this%bc_z%amr_restart(reconstruct, counter, tstep)

       this%iffinalised = .false.

       call neko_log%end_section(lvl = NEKO_LOG_VERBOSE)
    end if

  end subroutine non_normal_amr_restart

end module non_normal
