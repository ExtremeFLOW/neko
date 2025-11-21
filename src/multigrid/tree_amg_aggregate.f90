! Copyright (c) 2024, The Neko Authors
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
!> Implements an aggregation for TreeAMG hierarchy structure.
module tree_amg_aggregate
  use tree_amg, only : tamg_lvl_init, tamg_node_init, tamg_hierarchy_t
  use utils, only : neko_error, linear_index
  use num_types, only : rp, dp, i8
  use comm, only : NEKO_COMM, MPI_REAL_PRECISION, pe_size
  use mpi_f08, only : MPI_Allreduce, MPI_INTEGER, MPI_INTEGER8, &
       MPI_MIN, MPI_MAX, MPI_SUM
  use mesh, only : mesh_t
  use logger, only : neko_log, LOG_SIZE
  implicit none

  type, public :: tamg_agg_monitor_t
     ! Summary info
     integer :: level
     integer :: num_dofs
     integer :: num_aggs
     ! aggregation progress info
     integer :: phase1_naggs
     integer :: phase1_ndof_aggregated
     integer :: phase2_naggs
     integer :: phase2_ndof_aggregated
  end type tamg_agg_monitor_t

contains

  !> Aggregaiton on finest level
  !! Aggregates all dofs in an element into a single aggregate
  !! @param tamg TreeAMG hierarchy data structure being aggregated
  !! @param lx Number of dofs in x direction per element
  !! @param ly Number of dofs in y direction per element
  !! @param lz Number of dofs in z direction per element
  !! @param ne Number of elements
  subroutine aggregate_finest_level(tamg, lx, ly, lz, ne)
    type(tamg_hierarchy_t), intent(inout) :: tamg
    integer, intent(in) :: lx, ly, lz, ne
    integer :: i, j, k, l, nl, nt
    integer :: lid, gid_ptr
    integer :: lvl_id
    lvl_id = 1
    ! Count things
    nl = lx*ly*lz
    nt = nl*ne
    ! Allocate. For finest level, each aggregate is a node.
    call tamg_lvl_init(tamg%lvl(lvl_id), lvl_id, ne, nt)
    gid_ptr = 1
    do l = 1, ne
       call tamg_node_init( tamg%lvl(lvl_id)%nodes(l), l, nl)
       ! Fill the nodes
       lid = 0
       do k = 1, lz
          do j = 1, ly
             do i = 1, lx
                lid = lid + 1
                tamg%lvl(lvl_id)%nodes(l)%dofs(lid) = linear_index(i, j, k, l, &
                     lx, ly, lz)

                tamg%lvl(lvl_id)%map_f2c(linear_index(i,j,k,l,lx,ly,lz)) = l
                gid_ptr = gid_ptr + 1
             end do
          end do
       end do
    end do

    call aggregation_monitor_finest(lvl_id,nt,ne)

  end subroutine aggregate_finest_level

  !> First pass of a greedy aggregation
  !! Loop through all dofs and aggregate on dof that has all unaggregated neighbors
  !! @param naggs The number of aggregates that have been created
  !! @param max_aggs The maximum number of aggregates to allow to be created
  !! @param facet_neigh Dof adjacency information
  !! @param offset_el Offset for facet_neigh
  !! @param n_facet Max number of adjecnt dofs (ex. 6 if element is a cube)
  !! @param is_aggregated Array containing aggregate info. Maps from dof id to agg id
  !! @param aggregate_size Array containing the size of each aggregate
  subroutine agg_greedy_first_pass(naggs, max_aggs, n_elements, &
       facet_neigh, offset_el, n_facet, is_aggregated, aggregate_size)
    integer, intent(inout):: naggs
    integer, intent(in) :: max_aggs, n_elements
    integer, intent(in) :: facet_neigh(:,:)
    integer, intent(in) :: offset_el, n_facet
    integer, intent(inout) :: is_aggregated(:)
    integer, allocatable, intent(inout) :: aggregate_size(:)
    integer, allocatable :: as_tmp(:)
    integer, allocatable :: rand_order(:)
    real(kind=dp) :: random_value, r
    integer :: i, side, nhbr, j, tmp
    logical :: no_nhbr_agg

    ! Initialize a random permutation
    allocate( rand_order( n_elements ) )
    do i = 1, n_elements
       rand_order(i) = i
    end do
    ! Shuffle rand_order using Fisher-Yates algorithm
    do i = n_elements, 2, -1
       call random_number(r)
       j = int(r * real(i, kind=rp)) + 1
       tmp = rand_order(i)
       rand_order(i) = rand_order(j)
       rand_order(j) = tmp
    end do

!-----!    do while (naggs .le. max_aggs)
!-----!      call random_number(random_value)
!-----!      i = floor(random_value * n_elements + 1)
    !do i = 1, n_elements! THE NON-RANDOM VERSION...
    do tmp = 1, n_elements
       i = rand_order(tmp)
       if (is_aggregated(i) .eq. -1) then
          ! Check to see if any of the points neighbors are aggregated
          no_nhbr_agg = .true.
          do side = 1, n_facet
             nhbr = facet_neigh(side, i) - offset_el
             if ((nhbr .gt. 0).and.(nhbr .le. n_elements)) then! if nhbr exists
                if (is_aggregated(nhbr) .ne. -1) then
                   no_nhbr_agg = .false.
                end if
             end if
          end do! side

          ! if no neighbors are aggregated, create new aggregate
          if (no_nhbr_agg) then
             naggs = naggs + 1
             is_aggregated(i) = naggs
             if(size(aggregate_size).lt.naggs) then
                allocate(as_tmp(naggs + 20))
                as_tmp(1:size(aggregate_size)) = aggregate_size
                call move_alloc(as_tmp, aggregate_size)
             end if
             aggregate_size(naggs) = 1
             do side = 1, n_facet
                nhbr = facet_neigh(side, i) - offset_el
                if ((nhbr .gt. 0).and.(nhbr .le. n_elements)) then! if nhbr exists
                   if (is_aggregated(nhbr) .eq. -1) then
                      is_aggregated(nhbr) = naggs
                      aggregate_size(naggs) = aggregate_size(naggs) + 1
                   end if
                end if
             end do! side
          end if! no_nhbr_agg
       end if! is_aggregated(i)
    end do
  end subroutine agg_greedy_first_pass

  !> Second pass of a greedy aggregation
  !! Loop through all unaggregated dofs and add them to a neighboring aggregate.
  !! If no neighboring aggregates, create a new aggregate.
  !! @param naggs The number of aggregates that have ben created
  !! @param max_aggs The maximum number of aggregates to allow to be created
  !! @param facet_neigh Dof adjacency information
  !! @param offset_el Offset for facet_neigh
  !! @param n_facet Max number of adjecnt dofs (ex. 6 if element is a cube)
  !! @param is_aggregated Array containing aggregate info. Maps from dof id to agg id
  !! @param aggregate_size Array containing the size of each aggregate
  subroutine agg_greedy_second_pass(naggs, max_aggs, n_elements, &
       facet_neigh, offset_el, n_facet, is_aggregated, aggregate_size)
    integer, intent(inout):: naggs
    integer, intent(in) :: max_aggs, n_elements
    integer, intent(in) :: facet_neigh(:,:)
    integer, intent(in) :: offset_el, n_facet
    integer, intent(inout) :: is_aggregated(:)
    integer, intent(inout) :: aggregate_size(:)
    integer :: i, side, nhbr
    integer :: tnt_agg, tst_agg, tnt_size, tst_size

    ! Add remaining unaggregated nodes to aggregates
    do i = 1, n_elements
       if (is_aggregated(i) .eq. -1) then
          ! dof i is unaggregated. Check neighbors, add to smallest neighbor
          tnt_agg = -1
          tnt_size = 999!TODO: replace with large number
          tst_agg = -1
          tst_size = 999!TODO: replace with large number
          do side = 1, n_facet
             nhbr = facet_neigh(side, i) - offset_el
             if ((nhbr .gt. 0).and.(nhbr .le. n_elements)) then! if nhbr exists
                if (is_aggregated(nhbr) .ne. -1) then
                   tst_agg = is_aggregated(nhbr)
                   tst_size = aggregate_size(tst_agg)
                   if (tst_size .lt. tnt_size) then
                      tnt_size = tst_size
                      tnt_agg = tst_agg
                   end if
                end if
             end if
          end do

          if (tnt_agg .ne. -1) then
             ! if neighbor aggregate found add to that aggregate
             is_aggregated(i) = tnt_agg
             aggregate_size(tnt_agg) = aggregate_size(tnt_agg) + 1
          else
             ! if none of the neignbors are aggregated. might as well make a new aggregate
             naggs = naggs + 1
             if (naggs .gt. size(aggregate_size)) then!TODO: another movealoc here? Error? the max_aggs needs to change though...
                call neko_error("Creating too many aggregates... something might be wrong... try increasing max_aggs")
             end if
             is_aggregated(i) = naggs
             aggregate_size(naggs) = 1
             ! Add neighbors to aggregate if unaggregated
             do side = 1, n_facet
                nhbr = facet_neigh(side, i) - offset_el
                if ((nhbr .gt. 0).and.(nhbr .le. n_elements)) then! if nhbr exists
                   if (is_aggregated(nhbr) .eq. -1) then
                      aggregate_size(naggs) = aggregate_size(naggs) + 1
                      is_aggregated(nhbr) = naggs
                   end if
                end if
             end do
          end if

       end if
    end do
  end subroutine agg_greedy_second_pass

  !> Create information on which aggregates are "adjacent" to eachother
  !! @param agg_nhbr Aggregate adjacency information (same structure as facet_neigh)
  !! @param n_agg_nhbr The max number of aggregate neighbors over all aggregates
  !! @param n_elements The number of elements (dofs)
  !! @param facet_neigh Dof adjacency information
  !! @param offset_el Offset for facet_neigh
  !! @param n_facet Max number of adjecnt dofs (ex. 6 if element is a cube)
  !! @param is_aggregated Array containing aggregate info. Maps from dof id to agg id
  subroutine agg_fill_nhbr_info(agg_nhbr, n_agg_nhbr, n_elements, &
       facet_neigh, offset_el, n_facet, is_aggregated)
    integer, allocatable, intent(inout) :: agg_nhbr(:,:)
    integer, intent(inout) :: n_agg_nhbr
    integer, intent(in) :: n_elements
    integer, intent(in) :: facet_neigh(:,:)
    integer, intent(in) :: offset_el, n_facet
    integer, intent(inout) :: is_aggregated(:)
    integer :: i, j, side, nhbr, tst_agg, tnt_agg, n_nhbr_loc
    logical :: agg_added
    integer, allocatable :: agg_nhbr_counter(:)

    n_agg_nhbr = 0
    n_nhbr_loc = 0

    allocate(agg_nhbr_counter( maxval(is_aggregated)), source=0)

    ! Count max possible neighbors to an aggregate
    do i = 1, n_elements!TODO: this is the lazy expensive way...
       tnt_agg = is_aggregated(i)
       n_nhbr_loc = 0
       do side = 1, n_facet
          nhbr = facet_neigh(side,i) - offset_el
          if ((nhbr .gt. 0).and.(nhbr .le. n_elements)) then! if nhbr exists
             tst_agg = is_aggregated(nhbr)
             if (tst_agg .le. 0) call neko_error("Unaggregated element detected. We do not want to handle that here...")
             if (tst_agg .ne. tnt_agg) then
                agg_nhbr_counter(tnt_agg) = agg_nhbr_counter(tnt_agg) + 1
             end if
          end if
       end do
       n_agg_nhbr = max(n_agg_nhbr, agg_nhbr_counter(tnt_agg))
    end do

    ! Allocate for neighbor info
    allocate(agg_nhbr(n_agg_nhbr, maxval(is_aggregated)), source=-1)

    ! fill neighbor info
    do i = 1, n_elements!TODO: this is the lazy expensive way...
       tnt_agg = is_aggregated(i)
       do side = 1, n_facet
          nhbr = facet_neigh(side,i) - offset_el
          if ((nhbr .gt. 0).and.(nhbr .le. n_elements)) then! if nhbr exists
             tst_agg = is_aggregated(nhbr)
             if (tst_agg .le. 0) call neko_error("Unaggregated element detected. We do not want to handle that here...")
             if (tst_agg .ne. tnt_agg) then
                agg_added = .false.
                do j = 1, n_agg_nhbr
                   if ((agg_nhbr(j,tnt_agg) .eq. tst_agg)) then
                      agg_added = .true.
                   else if ((agg_nhbr(j,tnt_agg).eq.-1).and.(.not.agg_added)) then
                      agg_nhbr(j,tnt_agg) = tst_agg
                      agg_added = .true.
                      n_agg_nhbr = max(n_agg_nhbr, j)
                   end if
                end do! j
                if (.not.agg_added) call neko_error("Aggregates have too many neighbors... probably. Or some other error.")
             end if
          end if
       end do! side
    end do! i
  end subroutine agg_fill_nhbr_info

  !> Aggregates dofs based on adjacent dofs
  !! @param tamg TreeAMG hierarchy data structure being aggregated
  !! @param lvl_id The level id for which aggregates are being created
  !! @param max_aggs Target number of aggregates to create on level
  !! @param facet_neigh Input array that contains adjacency of dofs on level
  !! @param agg_nhbr Output array that contains adjacency of created aggregates
  subroutine aggregate_greedy(tamg, lvl_id, max_aggs, facet_neigh, agg_nhbr)
    type(tamg_hierarchy_t), intent(inout) :: tamg
    integer, intent(in) :: lvl_id
    integer, intent(in) :: max_aggs
    integer, intent(in) :: facet_neigh(:,:)
    integer, intent(inout), allocatable :: agg_nhbr(:,:)
    integer, allocatable :: is_aggregated(:)
    integer, allocatable :: aggregate_size(:)
    integer :: n_elements, naggs, n_facet, offset_el
    integer :: i, j, l, ntot, n_agg_facet, gid_ptr

    if (lvl_id .lt. 2) then
       call neko_error("For now, can only use greedy agg after elms have been aggregated to points (level 1)")
    else if (lvl_id .eq. 2) then
       n_facet = 6 ! NEKO elements are hexes, thus have 6 face neighbors
       offset_el = tamg%msh%offset_el
    else
       n_facet = size(facet_neigh, 1)
       offset_el = 0
    end if

    naggs = 0
    n_elements = tamg%lvl(lvl_id-1)%nnodes
    allocate( is_aggregated( n_elements ) )
    allocate( aggregate_size( max_aggs ) )

    ! fill with false
    is_aggregated = -1

    ! Fill with large number
    aggregate_size = huge(i)!999999

    ! First pass of greedy aggregation.
    call agg_greedy_first_pass(naggs, max_aggs, n_elements, &
         facet_neigh, offset_el, n_facet, &
         is_aggregated, aggregate_size)

    call aggregation_monitor_phase1(lvl_id, n_elements, naggs, is_aggregated)

    ! Second pass of greedy aggregation, adding unaggregated dofs to neighboring aggregates.
    call agg_greedy_second_pass(naggs, max_aggs, n_elements, &
         facet_neigh, offset_el, n_facet, &
         is_aggregated, aggregate_size)

    call aggregation_monitor_phase2(lvl_id, n_elements, naggs, is_aggregated)

    if (.true.) then! if needed on next level...
       call agg_fill_nhbr_info(agg_nhbr, n_agg_facet, n_elements, &
            facet_neigh, offset_el, n_facet, is_aggregated)
    end if

    ! count things
    ntot = 0
    do l = 1, naggs
       ntot = ntot + aggregate_size(l)
    end do
    ! Allocate and fill lvl and nodes
    call tamg_lvl_init( tamg%lvl(lvl_id), lvl_id, naggs, ntot)
    ntot = 0
    gid_ptr = 1
    do l = 1, naggs
       call tamg_node_init( tamg%lvl(lvl_id)%nodes(l), l, aggregate_size(l))
       j = 0
       do i = 1, n_elements!TODO: this is the lazy expensive way...
          if (is_aggregated(i) .eq. l) then
             j = j+1
             tamg%lvl(lvl_id)%nodes(l)%dofs(j) = i

             tamg%lvl(lvl_id)%map_f2c(i) = l
             gid_ptr = gid_ptr + 1
          end if
       end do
       if (j .ne. tamg%lvl(lvl_id)%nodes(l)%ndofs) then
          call neko_error("Aggregation problem. Not enough dofs in node.")
       end if
       ntot = ntot + aggregate_size(l)
    end do

    call aggregation_monitor_final(lvl_id,ntot,naggs)

    deallocate( is_aggregated )
    deallocate( aggregate_size )
  end subroutine aggregate_greedy

  !> Aggregate all dofs to a single point to form a tree-like structure.
  !! @param tamg TreeAMG hierarchy data structure being aggregated
  !! @param lvl_id The level id for which aggregates are being created
  subroutine aggregate_end( tamg, lvl_id)
    type(tamg_hierarchy_t), intent(inout) :: tamg
    integer, intent(in) :: lvl_id
    integer :: nt, i
    ! link all branches together at a point
    nt = tamg%lvl(lvl_id-1)%nnodes
    ! Allocate lvl
    call tamg_lvl_init( tamg%lvl(lvl_id), lvl_id, 1, nt)

    ! Allocate node
    call tamg_node_init( tamg%lvl(lvl_id)%nodes(1), 1, nt)
    ! Fill node
    do i = 1, tamg%lvl(lvl_id-1)%nnodes
       tamg%lvl(lvl_id)%nodes(1)%dofs(i) = tamg%lvl(lvl_id-1)%nodes(i)%gid

       tamg%lvl(lvl_id)%map_f2c(i) = 1
    end do

  end subroutine aggregate_end

  subroutine aggregation_monitor_finest(lvl,ndof,nagg)
    integer, intent(in) :: lvl,ndof,nagg
    integer :: na_max, na_min, na_avg, na_sum
    character(len=LOG_SIZE) :: log_buf

    write(log_buf, '(A8,I2,A37)') '-- level',lvl,'-- Aggregation: Element-as-Aggregate'
    !write(log_buf, '(A44)') 'Aggregation: Element-as-Aggregate'
    call neko_log%message(log_buf)

    call MPI_ALLREDUCE(nagg, na_max, 1, MPI_INTEGER, MPI_MAX, NEKO_COMM)
    call MPI_ALLREDUCE(nagg, na_min, 1, MPI_INTEGER, MPI_MIN, NEKO_COMM)
    call MPI_ALLREDUCE(nagg, na_sum, 1, MPI_INTEGER, MPI_SUM, NEKO_COMM)
    na_avg = na_sum / pe_size
    write(log_buf, '(A35,I6,A1,I6,A1,I6,A1)') 'Number of Aggregates: (',na_min,',',na_avg,',',na_max,')'
    call neko_log%message(log_buf)

  end subroutine aggregation_monitor_finest

  subroutine aggregation_monitor_phase1(lvl,ndof,nagg,is_aggregated)
    integer, intent(in) :: lvl,ndof,nagg
    integer, intent(in) :: is_aggregated(:)
    integer :: num_aggregated, i, na_max, na_min, na_avg, na_sum
    real(kind=rp) :: agg_frac
    integer(kind=i8) :: glb_dof, glb_aggd, loc_dof, loc_aggd
    character(len=LOG_SIZE) :: log_buf
    num_aggregated = 0
    do i = 1, ndof
       if (is_aggregated(i) .gt. -1) then
          num_aggregated = num_aggregated + 1
       end if
    end do

    !write(log_buf, '(A8,I2,A24)') '-- level',lvl,'-- Aggregation: phase1'
    write(log_buf, '(A27)') 'Aggregation: phase1'
    call neko_log%message(log_buf)

    loc_aggd = int(num_aggregated, i8)
    loc_dof = int(ndof, i8)
    call MPI_ALLREDUCE(loc_dof, glb_dof, 1, MPI_INTEGER8, MPI_SUM, NEKO_COMM)
    call MPI_ALLREDUCE(loc_aggd, glb_aggd, 1, MPI_INTEGER8, MPI_SUM, NEKO_COMM)
    agg_frac = real(glb_aggd,rp) / real(glb_dof,rp)

    call MPI_ALLREDUCE(nagg, na_max, 1, MPI_INTEGER, MPI_MAX, NEKO_COMM)
    call MPI_ALLREDUCE(nagg, na_min, 1, MPI_INTEGER, MPI_MIN, NEKO_COMM)
    call MPI_ALLREDUCE(nagg, na_sum, 1, MPI_INTEGER, MPI_SUM, NEKO_COMM)
    na_avg = na_sum / pe_size
    write(log_buf, '(A35,I6,A1,I6,A1,I6,A1)') 'Number of Aggregates: (',na_min,',',na_avg,',',na_max,')'
    call neko_log%message(log_buf)

    call MPI_ALLREDUCE(num_aggregated, na_max, 1, MPI_INTEGER, MPI_MAX, NEKO_COMM)
    call MPI_ALLREDUCE(num_aggregated, na_min, 1, MPI_INTEGER, MPI_MIN, NEKO_COMM)
    call MPI_ALLREDUCE(num_aggregated, na_sum, 1, MPI_INTEGER, MPI_SUM, NEKO_COMM)
    na_avg = na_sum / pe_size
    write(log_buf, '(A35,I6,A1,I6,A1,I6,A1,F6.2)') 'Aggregated: (',na_min,',',na_avg,',',na_max,')',(agg_frac*100)
    call neko_log%message(log_buf)
  end subroutine aggregation_monitor_phase1

  subroutine aggregation_monitor_phase2(lvl,ndof,nagg,is_aggregated)
    integer, intent(in) :: lvl,ndof,nagg
    integer, intent(in) :: is_aggregated(:)
    integer :: num_aggregated, i, na_max, na_min, na_avg, na_sum
    real(kind=rp) :: agg_frac
    integer(kind=i8) :: glb_dof, glb_aggd, loc_dof, loc_aggd
    character(len=LOG_SIZE) :: log_buf
    num_aggregated = 0
    do i = 1, ndof
       if (is_aggregated(i) .gt. -1) then
          num_aggregated = num_aggregated + 1
       end if
    end do
    !write(log_buf, '(A8,I2,A24)') '-- level',lvl,'-- Aggregation: phase2'
    write(log_buf, '(A27)') 'Aggregation: phase2'
    call neko_log%message(log_buf)

    loc_aggd = int(num_aggregated, i8)
    loc_dof = int(ndof, i8)
    call MPI_ALLREDUCE(loc_dof, glb_dof, 1, MPI_INTEGER8, MPI_SUM, NEKO_COMM)
    call MPI_ALLREDUCE(loc_aggd, glb_aggd, 1, MPI_INTEGER8, MPI_SUM, NEKO_COMM)
    agg_frac = real(glb_aggd,rp) / real(glb_dof,rp)

    call MPI_ALLREDUCE(nagg, na_max, 1, MPI_INTEGER, MPI_MAX, NEKO_COMM)
    call MPI_ALLREDUCE(nagg, na_min, 1, MPI_INTEGER, MPI_MIN, NEKO_COMM)
    call MPI_ALLREDUCE(nagg, na_sum, 1, MPI_INTEGER, MPI_SUM, NEKO_COMM)
    na_avg = na_sum / pe_size
    write(log_buf, '(A35,I6,A1,I6,A1,I6,A1)') 'Number of Aggregates: (',na_min,',',na_avg,',',na_max,')'
    call neko_log%message(log_buf)

    call MPI_ALLREDUCE(num_aggregated, na_max, 1, MPI_INTEGER, MPI_MAX, NEKO_COMM)
    call MPI_ALLREDUCE(num_aggregated, na_min, 1, MPI_INTEGER, MPI_MIN, NEKO_COMM)
    call MPI_ALLREDUCE(num_aggregated, na_sum, 1, MPI_INTEGER, MPI_SUM, NEKO_COMM)
    na_avg = na_sum / pe_size
    write(log_buf, '(A35,I6,A1,I6,A1,I6,A1,F6.2)') 'Aggregated: (',na_min,',',na_avg,',',na_max,')',(agg_frac*100)
    call neko_log%message(log_buf)
  end subroutine aggregation_monitor_phase2

  subroutine aggregation_monitor_final(lvl,ndof,nagg)
    integer, intent(in) :: lvl,ndof,nagg
    character(len=LOG_SIZE) :: log_buf
    !TODO: calculate min and max agg size
    !write(log_buf, '(A8,I2,A23,I6)') '-- level',lvl,'-- Aggregation: Done.', nagg
    write(log_buf, '(A26,I6)') 'Aggregation: Done.', nagg
    call neko_log%message(log_buf)
  end subroutine aggregation_monitor_final

  !> Aggregates pairs of dofs based on adjacent dofs
  !! @param tamg TreeAMG hierarchy data structure being aggregated
  !! @param lvl_id The level id for which aggregates are being created
  !! @param max_aggs Target number of aggregates to create on level
  !! @param facet_neigh Input array that contains adjacency of dofs on level
  !! @param agg_nhbr Output array that contains adjacency of created aggregates
  subroutine aggregate_pairs(tamg, lvl_id, max_aggs, facet_neigh, agg_nhbr)
    type(tamg_hierarchy_t), intent(inout) :: tamg
    integer, intent(in) :: lvl_id
    integer, intent(in) :: max_aggs
    integer, intent(in) :: facet_neigh(:,:)
    integer, intent(inout), allocatable :: agg_nhbr(:,:)
    integer, allocatable :: is_aggregated(:)
    integer, allocatable :: aggregate_size(:)
    integer, allocatable :: rand_order(:)
    integer :: n_elements, naggs, n_facet, offset_el
    integer :: i, j, l, ntot, n_agg_facet, gid_ptr, n_agg_nhbr
    integer :: n_pairs, tmp
    integer :: side, nhbr, nhbr_id
    real(kind=rp) :: nhbr_msr, nhbr_tst, r

    if (lvl_id .lt. 2) then
       call neko_error("For now, can only use greedy agg after elms have been aggregated to points (level 1)")
    else if (lvl_id .eq. 2) then
       n_facet = 6 !> NEKO elements are hexes, thus have 6 face neighbors
       offset_el = tamg%msh%offset_el
    else
       n_facet = size(facet_neigh, 1)
       offset_el = 0
    end if

    n_elements = tamg%lvl(lvl_id-1)%nnodes
    n_pairs = n_elements / 2
    allocate( is_aggregated( n_elements ) )
    allocate( aggregate_size( n_elements ) )
    ! fill with false
    is_aggregated = -1
    ! fill with large number
    aggregate_size = huge(i)

    ! Initialize a random permutation
    allocate( rand_order( n_elements ) )
    do i = 1, n_elements
       rand_order(i) = i
    end do
    ! Shuffle rand_order using Fisher-Yates algorithm
    do i = n_elements, 2, -1
       call random_number(r)
       j = int(r * real(i,kind=rp)) + 1
       tmp = rand_order(i)
       rand_order(i) = rand_order(j)
       rand_order(j) = tmp
    end do

    naggs = 0
    ! first pass of pair agg
    do tmp = 1, n_elements
       i = rand_order(tmp)
       if (is_aggregated(i) .eq. -1) then
          nhbr_id = -1
          nhbr_msr = -1.0_rp
          nhbr_tst = -1.0_rp
          do side = 1, n_facet
             nhbr = facet_neigh(side, i) - offset_el
             if ((nhbr .gt. 0).and.(nhbr .le. n_elements)) then ! nhbr exists
                if (is_aggregated(nhbr) .eq. -1) then
                   nhbr_tst = 1.0_rp! insert desired metric here
                   if (nhbr_tst .gt. nhbr_msr) then! if nhbr has goodest metric
                      nhbr_id = nhbr
                      nhbr_msr = nhbr_tst
                   end if
                end if! is_aggregated(nhbr)
             end if! nhbr exists
          end do! side

          if (nhbr_id .ne. -1) then
             naggs = naggs + 1
             is_aggregated(i) = naggs
             is_aggregated(nhbr_id) = naggs
             aggregate_size(naggs) = 2
          else! singleton, in theory we want to avoid
             naggs = naggs + 1
             is_aggregated(i) = naggs
             aggregate_size(naggs) = 1
          endif
       end if! is_aggregated(i)
    end do
    call agg_fill_nhbr_info( agg_nhbr, n_agg_nhbr, n_elements, &
         facet_neigh, offset_el, n_facet, is_aggregated)

    ! count things
    ntot = n_elements
    ! Allocate and fill lvl and nodes
    call tamg_lvl_init( tamg%lvl(lvl_id), lvl_id, naggs, ntot)
    ntot = 0
    gid_ptr = 1
    do l = 1, naggs
       call tamg_node_init( tamg%lvl(lvl_id)%nodes(l), l, aggregate_size(l))
       j = 0
       do i = 1, n_elements!TODO: this is the lazy expensive way...
          if (is_aggregated(i) .eq. l) then
             j = j+1
             tamg%lvl(lvl_id)%nodes(l)%dofs(j) = i

             tamg%lvl(lvl_id)%map_f2c(i) = l
             gid_ptr = gid_ptr + 1
          end if
       end do
       if (j .ne. tamg%lvl(lvl_id)%nodes(l)%ndofs) then
          call neko_error("Aggregation problem. Not enough dofs in node.")
       end if
       ntot = ntot + aggregate_size(l)
    end do
    call aggregation_monitor_final(lvl_id,ntot,naggs)
    deallocate( is_aggregated )
    deallocate( aggregate_size )
  end subroutine aggregate_pairs

end module tree_amg_aggregate
