! Copyright (c) 2021, The Neko Authors
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
!> Redistribution routines
module redist
  use mesh_field, only : mesh_fld_t, mesh_field_init, mesh_field_free
  use neko_mpi_types
  use mpi_f08
  use htable, only : htable_i4_t
  use point, only : point_t
  use stack, only : stack_i4_t, stack_nh_t, stack_nc_t, stack_nz_t
  use curve, only : curve_t
  use comm
  use mesh, only : mesh_t, NEKO_MSH_MAX_ZLBLS
  use nmsh, only : nmsh_hex_t, nmsh_zone_t, nmsh_curve_el_t
  use facet_zone, only : facet_zone_t, facet_zone_periodic_t
  use element, only : element_t
  implicit none
  private

  public :: redist_mesh

contains

  !> Redistribute a mesh @a msh according to new partitions
  subroutine redist_mesh(msh, parts)
    type(mesh_t), intent(inout), target :: msh    !< Mesh
    type(mesh_fld_t), intent(in) :: parts         !< Partitions
    type(stack_nh_t), allocatable :: new_mesh_dist(:)
    type(stack_nz_t), allocatable :: new_zone_dist(:)
    type(stack_nc_t), allocatable :: new_curve_dist(:)
    type(nmsh_hex_t) :: el
    class(nmsh_hex_t), pointer :: np(:)
    type(nmsh_hex_t), allocatable :: recv_buf_msh(:)
    class(nmsh_zone_t), pointer :: zp(:)
    type(nmsh_zone_t), allocatable :: recv_buf_zone(:)
    type(nmsh_curve_el_t), pointer :: cp(:)
    type(nmsh_curve_el_t), allocatable :: recv_buf_curve(:)
    class(element_t), pointer :: ep
    integer, allocatable :: recv_buf_idx(:), send_buf_idx(:)
    type(MPI_Status) :: status
    integer :: i, j, k, ierr, max_recv_idx, label
    integer :: src, dst, recv_size, gdim, tmp, new_el_idx, new_pel_idx
    integer :: max_recv(3)
    type(point_t) :: p(8)
    type(htable_i4_t) :: el_map, glb_map
    type(stack_i4_t) :: pe_lst


    !
    ! Reset possible periodic ids
    !
    call msh%reset_periodic_ids()

    !
    ! Extract new zone distributions
    !

    allocate(new_zone_dist(0:pe_size - 1))
    do i = 0, pe_size - 1
       call new_zone_dist(i)%init()
    end do

    call redist_zone(msh, msh%wall, 1, parts, new_zone_dist)
    call redist_zone(msh, msh%inlet, 2, parts, new_zone_dist)
    call redist_zone(msh, msh%outlet, 3, parts, new_zone_dist)
    call redist_zone(msh, msh%sympln, 4, parts, new_zone_dist)
    call redist_zone(msh, msh%periodic, 5, parts, new_zone_dist)

    do j = 1, NEKO_MSH_MAX_ZLBLS
       label = j
       call redist_zone(msh, msh%labeled_zones(j), 7, parts, &
            new_zone_dist, label)
    end do

    !
    ! Extract new curve info. distributions
    !
    allocate(new_curve_dist(0:pe_size - 1))
    do i = 0, pe_size - 1
       call new_curve_dist(i)%init()
    end do

    call redist_curve(msh, msh%curve, parts, new_curve_dist)

    !
    ! Extract new mesh distribution
    !

    allocate(new_mesh_dist(0:(pe_size - 1)))
    do i = 0, pe_size - 1
       call new_mesh_dist(i)%init()
    end do

    do i = 1, msh%nelv
       ep => msh%elements(i)%e
       el%el_idx = ep%id()
       do j = 1, 8
          el%v(j)%v_idx = ep%pts(j)%p%id()
          el%v(j)%v_xyz = ep%pts(j)%p%x
       end do
       call new_mesh_dist(parts%data(i))%push(el)
    end do


    gdim = msh%gdim
    call msh%free()

    max_recv = 0
    do i = 0, pe_size - 1
       max_recv(1) = max(max_recv(1), new_mesh_dist(i)%size())
       max_recv(2) = max(max_recv(2), new_zone_dist(i)%size())
       max_recv(3) = max(max_recv(3), new_curve_dist(i)%size())
    end do

    call MPI_Allreduce(MPI_IN_PLACE, max_recv, 3, MPI_INTEGER, &
         MPI_MAX, NEKO_COMM, ierr)
    allocate(recv_buf_msh(max_recv(1)))
    allocate(recv_buf_zone(max_recv(2)))
    allocate(recv_buf_curve(max_recv(3)))

    do i = 1, pe_size - 1
       src = modulo(pe_rank - i + pe_size, pe_size)
       dst = modulo(pe_rank + i, pe_size)

       ! We should use the %array() procedure, which works great for
       ! GNU, Intel and NEC, but it breaks horribly on Cray when using
       ! certain data types
       select type (nmd_array => new_mesh_dist(dst)%data)
       type is (nmsh_hex_t)
          call MPI_Sendrecv(nmd_array, &
               new_mesh_dist(dst)%size(), MPI_NMSH_HEX, dst, 0, recv_buf_msh, &
               max_recv(1), MPI_NMSH_HEX, src, 0, NEKO_COMM, status, ierr)
       end select
       call MPI_Get_count(status, MPI_NMSH_HEX, recv_size, ierr)

       do j = 1, recv_size
          call new_mesh_dist(pe_rank)%push(recv_buf_msh(j))
       end do

       ! We should use the %array() procedure, which works great for
       ! GNU, Intel and NEC, but it breaks horribly on Cray when using
       ! certain data types
       select type (nzd_array => new_zone_dist(dst)%data)
       type is (nmsh_zone_t)
          call MPI_Sendrecv(nzd_array, &
            new_zone_dist(dst)%size(), MPI_NMSH_ZONE, dst, 1, recv_buf_zone,&
            max_recv(2), MPI_NMSH_ZONE, src, 1, NEKO_COMM, status, ierr)
       end select
       call MPI_Get_count(status, MPI_NMSH_ZONE, recv_size, ierr)

       do j = 1, recv_size
          call new_zone_dist(pe_rank)%push(recv_buf_zone(j))
       end do

       call MPI_Sendrecv(new_curve_dist(dst)%array(), &
            new_curve_dist(dst)%size(), MPI_NMSH_CURVE, dst, 2, recv_buf_curve,&
            max_recv(3), MPI_NMSH_CURVE, src, 2, NEKO_COMM, status, ierr)
       call MPI_Get_count(status, MPI_NMSH_CURVE, recv_size, ierr)

       do j = 1, recv_size
          call new_curve_dist(pe_rank)%push(recv_buf_curve(j))
       end do
    end do

    deallocate(recv_buf_msh)
    deallocate(recv_buf_zone)
    deallocate(recv_buf_curve)

    do i = 0, pe_size - 1
       if (i .ne. pe_rank) then
          call new_mesh_dist(i)%free
          call new_zone_dist(i)%free
          call new_curve_dist(i)%free
       end if
    end do

    !
    ! Create a new mesh based on the given distribution
    !
    call msh%init(gdim, new_mesh_dist(pe_rank)%size())

    call el_map%init(new_mesh_dist(pe_rank)%size())
    call glb_map%init(new_mesh_dist(pe_rank)%size())

    ! We should use the %array() procedure, which works great for
    ! GNU, Intel and NEC, but it breaks horribly on Cray when using
    ! certain data types
    select type (np => new_mesh_dist(pe_rank)%data)
    type is (nmsh_hex_t)
       do i = 1, new_mesh_dist(pe_rank)%size()
          do j = 1, 8
             p(j) = point_t(np(i)%v(j)%v_xyz, np(i)%v(j)%v_idx)
          end do
          call msh%add_element(i, &
               p(1), p(2), p(3), p(4), p(5), p(6), p(7), p(8))

          if (el_map%get(np(i)%el_idx, tmp) .gt. 0) then
             ! Old glb to new local
             tmp = i
             call el_map%set(np(i)%el_idx, tmp)

             ! Old glb to new glb
             tmp = msh%elements(i)%e%id()
             call glb_map%set(np(i)%el_idx,  tmp)
          else
             call neko_error('Global element id already defined')
          end if
       end do
    end select
    call new_mesh_dist(pe_rank)%free()


    !
    ! Figure out new mesh distribution (necessary for periodic zones)
    !
    call pe_lst%init()

    ! We should use the %array() procedure, which works great for
    ! GNU, Intel and NEC, but it breaks horribly on Cray when using
    ! certain data types
    select type(zp => new_zone_dist(pe_rank)%data)
    type is (nmsh_zone_t)
       do i = 1, new_zone_dist(pe_rank)%size()
          if (zp(i)%type .eq. 5) then
             call pe_lst%push(zp(i)%p_e)
          end if
       end do
    end select

    max_recv_idx = 2 * pe_lst%size()
    call MPI_Allreduce(MPI_IN_PLACE, max_recv_idx, 1, MPI_INTEGER, &
         MPI_MAX, NEKO_COMM, ierr)
    allocate(recv_buf_idx(max_recv_idx))
    allocate(send_buf_idx(max_recv_idx))

    do i = 1, pe_size - 1
       src = modulo(pe_rank - i + pe_size, pe_size)
       dst = modulo(pe_rank + i, pe_size)

       ! We should use the %array() procedure, which works great for
       ! GNU, Intel and NEC, but it breaks horribly on Cray when using
       ! certain data types
       select type (pe_lst_array => pe_lst%data)
       type is (integer)
          call MPI_Sendrecv(pe_lst_array, &
               pe_lst%size(), MPI_INTEGER, dst, 0, recv_buf_idx, &
               max_recv_idx, MPI_INTEGER, src, 0, NEKO_COMM, status, ierr)
       end select
       call MPI_Get_count(status, MPI_INTEGER, recv_size, ierr)

       k = 0
       do j = 1, recv_size
          if (glb_map%get(recv_buf_idx(j), tmp) .eq. 0) then
             k = k + 1
             send_buf_idx(k) = recv_buf_idx(j)
             k = k + 1
             send_buf_idx(k) = tmp
          end if
       end do

       call MPI_Sendrecv(send_buf_idx, k, MPI_INTEGER, src, 1, &
            recv_buf_idx, max_recv_idx, MPI_INTEGER, dst, 1, &
            NEKO_COMM, status, ierr)
       call MPI_Get_count(status, MPI_INTEGER, recv_size, ierr)

       do j = 1, recv_size, 2
          call glb_map%set(recv_buf_idx(j),  recv_buf_idx(j+1))
       end do
    end do
    deallocate(recv_buf_idx)
    deallocate(send_buf_idx)
    call pe_lst%free()

    !
    ! Add zone data for new mesh distribution
    !
    zp => new_zone_dist(pe_rank)%array()
    do i = 1, new_zone_dist(pe_rank)%size()
       if (el_map%get(zp(i)%e, new_el_idx) .gt. 0) then
          call neko_error('Missing element after redistribution')
       end if
       select case(zp(i)%type)
       case(1)
          call msh%mark_wall_facet(zp(i)%f, new_el_idx)
       case(2)
          call msh%mark_inlet_facet(zp(i)%f, new_el_idx)
       case(3)
          call msh%mark_outlet_facet(zp(i)%f, new_el_idx)
       case(4)
          call msh%mark_sympln_facet(zp(i)%f, new_el_idx)
       case(5)
          if (glb_map%get(zp(i)%p_e, new_pel_idx) .gt. 0) then
             call neko_error('Missing periodic element after redistribution')
          end if

          call msh%mark_periodic_facet(zp(i)%f, new_el_idx, &
               zp(i)%p_f, new_pel_idx, zp(i)%glb_pt_ids)
       case(7)
          call msh%mark_labeled_facet(zp(i)%f, new_el_idx, zp(i)%p_f)
       end select
    end do
    do i = 1, new_zone_dist(pe_rank)%size()
       if (el_map%get(zp(i)%e, new_el_idx) .gt. 0) then
          call neko_error('Missing element after redistribution')
       end if
       select case(zp(i)%type)
       case(5)
          if (glb_map%get(zp(i)%p_e, new_pel_idx) .gt. 0) then
             call neko_error('Missing periodic element after redistribution')
          end if

          call msh%apply_periodic_facet(zp(i)%f, new_el_idx, &
               zp(i)%p_f, new_pel_idx, zp(i)%glb_pt_ids)
       end select
    end do

    call new_zone_dist(pe_rank)%free()

    !
    ! Add curve element information for new mesh distribution
    !
    cp => new_curve_dist(pe_rank)%array()
    do i = 1, new_curve_dist(pe_rank)%size()
       if (el_map%get(cp(i)%e, new_el_idx) .gt. 0) then
          call neko_error('Missing element after redistribution')
       end if
       call msh%mark_curve_element(new_el_idx, cp(i)%curve_data, cp(i)%type)
    end do
    call new_curve_dist(pe_rank)%free()

    call msh%finalize()

  end subroutine redist_mesh

  !> Fill redistribution list for zone data
  subroutine redist_zone(msh, z, type, parts, new_dist, label)
    type(mesh_t), intent(inout) :: msh
    class(facet_zone_t), intent(in) :: z
    integer, intent(in) :: type
    type(mesh_fld_t), intent(in) :: parts
    type(stack_nz_t), intent(inout), allocatable :: new_dist(:)
    integer, intent(in), optional :: label
    type(nmsh_zone_t) :: nmsh_zone
    integer :: i, j, zone_el, lbl

    if (present(label)) then
       lbl = label
    else
       lbl = 0
    end if

    select type(zp => z)
    type is (facet_zone_periodic_t)
       do i = 1, zp%size
          zone_el =  zp%facet_el(i)%x(2)
          nmsh_zone%e = zp%facet_el(i)%x(2) + msh%offset_el
          nmsh_zone%f = zp%facet_el(i)%x(1)
          nmsh_zone%p_e = zp%p_facet_el(i)%x(2)
          nmsh_zone%p_f = zp%p_facet_el(i)%x(1)
          nmsh_zone%glb_pt_ids = zp%p_ids(i)%x
          nmsh_zone%type = type
          call new_dist(parts%data(zone_el))%push(nmsh_zone)
       end do
    type is (facet_zone_t)
       do i = 1, zp%size
          zone_el =  zp%facet_el(i)%x(2)
          nmsh_zone%e = zp%facet_el(i)%x(2) + msh%offset_el
          nmsh_zone%f = zp%facet_el(i)%x(1)
          nmsh_zone%p_f = lbl ! Labels are encoded in the periodic facet...
          nmsh_zone%type = type
          call new_dist(parts%data(zone_el))%push(nmsh_zone)
       end do
    end select
  end subroutine redist_zone

  !> Fill redistribution list for zone data
  subroutine redist_curve(msh, c, parts, new_dist)
    type(mesh_t), intent(inout) :: msh
    type(curve_t), intent(in) :: c
    type(mesh_fld_t), intent(in) :: parts
    type(stack_nc_t), intent(inout), allocatable :: new_dist(:)
    type(nmsh_curve_el_t) :: nmsh_curve
    integer :: i, j, curve_el

    do i = 1, c%size
       curve_el = c%curve_el(i)%el_idx
       nmsh_curve%e = curve_el + msh%offset_el
       nmsh_curve%curve_data = c%curve_el(i)%curve_data
       nmsh_curve%type = c%curve_el(i)%curve_type
       call new_dist(parts%data(curve_el))%push(nmsh_curve)
    end do

  end subroutine redist_curve

end module redist
