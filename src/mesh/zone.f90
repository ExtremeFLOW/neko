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
!> Defines a zone as a subset of facets in a mesh
module zone
  use tuple, only : tuple_i4_t, tuple4_i4_t
  use stack, only : stack_i4t2_t, stack_i4t4_t
  use utils
  implicit none
  private
  
  type, public :: zone_t
     type(tuple_i4_t), allocatable :: facet_el(:)
     integer :: size = 0
     logical, private :: finalized = .false.
     type(stack_i4t2_t), private :: scratch
   contains
     procedure, pass(z) :: init => zone_init
     procedure, pass(z) :: free => zone_free
     procedure, pass(z) :: finalize => zone_finalize
     procedure, pass(z) :: add_facet => zone_add_facet
  end type zone_t

  type, public, extends(zone_t) :: zone_periodic_t
     type(tuple_i4_t), allocatable :: p_facet_el(:)
     type(stack_i4t2_t), private :: p_scratch
     type(tuple4_i4_t), allocatable :: p_ids(:)  !< Periodic ids, same for each periodic point
     type(stack_i4t4_t), private :: p_id_scratch
     type(tuple4_i4_t), allocatable :: org_ids(:)  !< Original ids point ids
     type(stack_i4t4_t), private :: org_id_scratch
   contains
     procedure, pass(z) :: init => zone_periodic_init
     procedure, pass(z) :: free => zone_periodic_free
     procedure, pass(z) :: finalize => zone_periodic_finalize
     procedure, pass(z) :: add_periodic_facet => zone_periodic_add_facet
  end type zone_periodic_t
 
contains

  !> Initialize a zone
  subroutine zone_init(z, size)
    class(zone_t), intent(inout) :: z
    integer, optional :: size

    call zone_free(z)

    if (present(size)) then
       call z%scratch%init(size)
    else
       call z%scratch%init()
    end if
    
  end subroutine zone_init

  !> Deallocate a zone
  subroutine zone_free(z)
    class(zone_t), intent(inout) :: z
    if (allocated(z%facet_el)) then
       deallocate(z%facet_el)
    end if

    z%finalized = .false.
    z%size = 0

    call z%scratch%free()
    
  end subroutine zone_free

  !> Finalize a zone list
  !! @details Create a static list of (facet,el) tuples
  subroutine zone_finalize(z)
    class(zone_t), intent(inout) :: z
    type(tuple_i4_t), pointer :: tp(:)
    integer :: i
    
    if (.not. z%finalized) then

       allocate(z%facet_el(z%scratch%size()))
       
       tp => z%scratch%array()
       do i = 1, z%scratch%size()
          z%facet_el(i) = tp(i)
       end do

       z%size = z%scratch%size()

       call z%scratch%clear()

       z%finalized = .true.
       
    end if
    
  end subroutine zone_finalize

  !> Add a (facet, el) tuple to an unfinalized zone
  subroutine zone_add_facet(z, facet, el)
    class(zone_t), intent(inout) :: z
    integer, intent(in) :: facet   !< Facet in the zone
    integer, intent(in) :: el      !< Element  in the zone
    type(tuple_i4_t) :: t

    if (z%finalized) then
       call neko_error('Zone already finalized')
    end if

    t%x = (/ facet, el /)
    call z%scratch%push(t)
    
  end subroutine zone_add_facet

    !> Initialize a periodic zone
  subroutine zone_periodic_init(z, size)
    class(zone_periodic_t), intent(inout) :: z
    integer, optional :: size

    call z%free()

    if (present(size)) then
       call zone_init(z, size)
       call z%p_scratch%init(size)
       call z%p_id_scratch%init(size)
       call z%org_id_scratch%init(size)
    else
       call zone_init(z)
       call z%p_scratch%init()
       call z%p_id_scratch%init()
       call z%org_id_scratch%init()
    end if
    
  end subroutine zone_periodic_init

  !> Deallocate a zone
  subroutine zone_periodic_free(z)
    class(zone_periodic_t), intent(inout) :: z

    call zone_free(z)

    if (allocated(z%p_facet_el)) then
       deallocate(z%p_facet_el)
    end if

    if (allocated(z%p_ids)) then
       deallocate(z%p_ids)
    end if
    if (allocated(z%org_ids)) then
       deallocate(z%org_ids)
    end if

    call z%p_scratch%free()
    call z%p_id_scratch%free()
    call z%org_id_scratch%free()
    
  end subroutine zone_periodic_free

  !> Finalize a periodic zone list
  !! @details Create a static list of (facet,el) tuples
  subroutine zone_periodic_finalize(z)
    class(zone_periodic_t), intent(inout) :: z
    type(tuple_i4_t), pointer :: tp(:)
    type(tuple4_i4_t), pointer :: tp2(:)
    type(tuple4_i4_t), pointer :: tp3(:)
    integer :: i
    
    if (.not. z%finalized) then

       call zone_finalize(z)

       if (z%size .ne. z%p_scratch%size()) then
          call neko_error('Zone size mismatch')
       end if

       allocate(z%p_facet_el(z%size))
       allocate(z%p_ids(z%size))
       allocate(z%org_ids(z%size))
       
       tp => z%p_scratch%array()
       do i = 1, z%size
          z%p_facet_el(i) = tp(i)
       end do
       tp2 => z%p_id_scratch%array()
       do i = 1, z%size
          z%p_ids(i) = tp2(i)
       end do
       tp3 => z%org_id_scratch%array()
       do i = 1, z%size
          z%org_ids(i) = tp3(i)
       end do

       call z%p_scratch%clear()
       call z%p_id_scratch%clear()
       call z%org_id_scratch%clear()

    end if
    
  end subroutine zone_periodic_finalize

  !> Add a (facet, el) tuple to an unfinalized zone
  subroutine zone_periodic_add_facet(z, facet, el, p_facet, p_el, pids, org_ids)
    class(zone_periodic_t), intent(inout) :: z
    integer, intent(in) :: facet   !< Facet in the zone
    integer, intent(in) :: el      !< Element  in the zone
    integer, intent(in) :: p_facet !< Facet at periodic length
    integer, intent(in) :: p_el    !< Element at periodic length
    integer, intent(in) :: pids(4) !< Periodic id of points
    integer, intent(in) :: org_ids(4) !< Original id of points
    type(tuple_i4_t) :: t
    type(tuple4_i4_t) :: t2
    type(tuple4_i4_t) :: t3

    if (z%finalized) then
       call neko_error('Zone already finalized')
    end if

    call z%add_facet(facet, el)

    t%x = (/ p_facet, p_el /)
    call z%p_scratch%push(t)
    t2%x = pids
    call z%p_id_scratch%push(t2)
    t3%x = org_ids
    call z%org_id_scratch%push(t3)
    
  end subroutine zone_periodic_add_facet

end module zone
