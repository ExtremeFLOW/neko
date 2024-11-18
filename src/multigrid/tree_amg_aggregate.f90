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
!> Implements an aggregation for TreeAMG hierarchy structure.
module tree_amg_aggregate
  use tree_amg
  use utils
  use num_types
  use mesh, only : mesh_t
  implicit none

contains

  !> Aggregaiton on finest level
  !> Aggregates all dofs in an element into a single aggregate
  !! @param tamg TreeAMG hierarchy data structure being aggregated
  !! @param lx Number of dofs in x direction per element
  !! @param ly Number of dofs in y direction per element
  !! @param lz Number of dofs in z direction per element
  !! @param ne Number of elements
  subroutine aggregate_finest_level(tamg, lx, ly, lz, ne)
    type(tamg_hierarchy_t), intent(inout) :: tamg
    integer, intent(in) :: lx, ly, lz, ne
    integer :: i, j, k, l, nl, nt
    integer :: lid
    integer :: lvl_id
    lvl_id = 1
    !> Allocate. For finest level, each aggregate is a node.
    call tamg_lvl_init(tamg%lvl(lvl_id), lvl_id, ne)
    !> Count things
    nl = lx*ly*lz
    nt = nl*ne
    do l = 1, ne
      call tamg_node_init( tamg%lvl(lvl_id)%nodes(l), l, nl)
      !> Fill the nodes
      lid = 0
      do k = 1, lz
        do j = 1, ly
          do i = 1, lx
            lid = lid + 1
            tamg%lvl(lvl_id)%nodes(l)%dofs(lid) = linear_index(i,j,k,l,lx,ly,lz)
          end do
        end do
      end do
    end do
    tamg%lvl(lvl_id)%fine_lvl_dofs = nt
    allocate(tamg%lvl(lvl_id)%wrk_in( nt ))
    allocate(tamg%lvl(lvl_id)%wrk_out( nt ))
    print *, "work allocated on lvl", lvl_id, "dofs", nt
  end subroutine aggregate_finest_level

  !> First pass of a greedy aggregation
  !> Loop through all dofs and aggregate on dof that has all unaggregated neighbors
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
    integer, intent(inout) :: aggregate_size(:)
    real(kind=dp) :: random_value
    integer :: i, side, nhbr
    logical :: no_nhbr_agg

!-----!    do while (naggs .le. max_aggs)
!-----!      call random_number(random_value)
!-----!      i = floor(random_value * n_elements + 1)
    do i = 1, n_elements!> THE NON-RANDOM VERSION...
      if (is_aggregated(i) .eq. -1) then
        !> Check to see if any of the points neighbors are aggregated
        no_nhbr_agg = .true.
        do side = 1, n_facet
          nhbr = facet_neigh(side, i) - offset_el
          if ((nhbr .gt. 0).and.(nhbr .le. n_elements)) then!> if nhbr exists
            if (is_aggregated(nhbr) .ne. -1) then
              no_nhbr_agg = .false.
            end if
          end if
        end do! side

        !> if no neighbors are aggregated, create new aggregate
        if (no_nhbr_agg) then
          naggs = naggs + 1
          is_aggregated(i) = naggs
          aggregate_size(naggs) = 1
          do side = 1, n_facet
            nhbr = facet_neigh(side, i) - offset_el
            if ((nhbr .gt. 0).and.(nhbr .le. n_elements)) then!> if nhbr exists
              if (is_aggregated(nhbr) .eq. -1) then
                is_aggregated(nhbr) = naggs
                aggregate_size(naggs) = aggregate_size(naggs) + 1
              end if
            end if
          end do! side
        end if! no_nhbr_agg
      end if! is_aggregated(i)
    end do
    print *, "done with first pass of aggregation"
  end subroutine agg_greedy_first_pass

  !> Second pass of a greedy aggregation
  !> Loop through all unaggregated dofs and add them to a neighboring aggregate.
  !> If no neighboring aggregates, create a new aggregate.
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

    !> Add remaining unaggregated nodes to aggregates
    do i = 1, n_elements
      if (is_aggregated(i) .eq. -1) then
        !> dof i is unaggregated. Check neighbors, add to smallest neighbor
        tnt_agg = -1
        tnt_size = 999!TODO: replace with large number
        tst_agg = -1
        tst_size = 999!TODO: replace with large number
        do side = 1, n_facet
          nhbr = facet_neigh(side, i) - offset_el
          if ((nhbr .gt. 0).and.(nhbr .le. n_elements)) then!> if nhbr exists
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
          !> if neighbor aggregate found add to that aggregate
          is_aggregated(i) = tnt_agg
          aggregate_size(tnt_agg) = aggregate_size(tnt_agg) + 1
        else
          !> if none of the neignbors are aggregated. might as well make a new aggregate
          naggs = naggs + 1
          if (naggs .gt. max_aggs*2) then
            print *, "AGGS:", naggs
            call neko_error("Creating too many aggregates... something might be wrong... try increasing max_aggs")
          end if
          is_aggregated(i) = naggs
          aggregate_size(naggs) = 1
          !> Add neighbors to aggregate if unaggregated
          do side = 1, n_facet
            nhbr = facet_neigh(side, i) - offset_el
            if ((nhbr .gt. 0).and.(nhbr .le. n_elements)) then!> if nhbr exists
              if (is_aggregated(nhbr) .eq. -1) then
                aggregate_size(naggs) = aggregate_size(naggs) + 1
                is_aggregated(nhbr) = naggs
              end if
            end if
          end do
        end if

      end if
    end do
    print *, "done with second pass of aggregation: number of aggregates", naggs
  end subroutine agg_greedy_second_pass

  !> Create information on which aggregates are "adjacent" to eachother
  !! @param agg_nhbr Aggregate adjacency information (same structure as facet_neigh)
  !! @param n_agg_nhbr The max number of aggregate neighbors over all aggregates
  !! @param n_elements The number of elements (dofs)
  !! @param facet_neigh Dof adjacency information
  !! @param offset_el Offset for facet_neigh
  !! @param n_facet Max number of adjecnt dofs (ex. 6 if element is a cube)
  !! @param is_aggregated Array containing aggregate info. Maps from dof id to agg id
  !! @param aggregate_size Array containing the size of each aggregate
  subroutine agg_fill_nhbr_info(agg_nhbr, n_agg_nhbr, n_elements, &
      facet_neigh, offset_el, n_facet, is_aggregated, aggregate_size)
    integer, intent(inout) :: agg_nhbr(:,:)
    integer, intent(inout) :: n_agg_nhbr
    integer, intent(in) :: n_elements
    integer, intent(in) :: facet_neigh(:,:)
    integer, intent(in) :: offset_el, n_facet
    integer, intent(inout) :: is_aggregated(:)
    integer, intent(inout) :: aggregate_size(:)
    integer :: i, j, side, nhbr, tst_agg, tnt_agg
    logical :: agg_added

    n_agg_nhbr = 0

    do i = 1, n_elements!TODO: this is the lazy expensive way...
      tnt_agg = is_aggregated(i)
      do side = 1, n_facet
        nhbr = facet_neigh(side,i) - offset_el
        if ((nhbr .gt. 0).and.(nhbr .le. n_elements)) then!> if nhbr exists
          tst_agg = is_aggregated(nhbr)
          if (tst_agg .le. 0) call neko_error("Unaggregated element detected. We do not want to handle that here...")
          if (tst_agg .ne. tnt_agg) then
            agg_added = .false.
            do j = 1, 20!TODO: this hard-coded value
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
    integer :: i, j, l, ntot, n_agg_facet

    if (lvl_id .lt. 2) then
      call neko_error("For now, can only use greedy agg after elms have been aggregated to points (level 1)")
    else if (lvl_id .eq. 2) then
      n_facet = 6 !> NEKO elements are hexes, thus have 6 face neighbors
      offset_el = tamg%msh%offset_el
    else
      n_facet = 20!TODO: this hard-coded value. how many neighbors can there be?
      offset_el = 0
    end if

    naggs = 0
    n_elements = tamg%lvl(lvl_id-1)%nnodes
    allocate( is_aggregated( n_elements ) )
    allocate( aggregate_size( max_aggs*2 ) )

    !> fill with false
    is_aggregated = -1

    !> Fill with large number
    aggregate_size = huge(i)!999999

    !> First pass of greedy aggregation.
    call agg_greedy_first_pass(naggs, max_aggs, n_elements, &
      facet_neigh, offset_el, n_facet, &
      is_aggregated, aggregate_size)

    !> Second pass of greedy aggregation, adding unaggregated dofs to neighboring aggregates.
    call agg_greedy_second_pass(naggs, max_aggs, n_elements, &
      facet_neigh, offset_el, n_facet, &
      is_aggregated, aggregate_size)

    if (.true.) then!> if needed on next level...
      allocate( agg_nhbr(20, max_aggs*2) )!TODO: this hard-coded n_facet value (20)...
      agg_nhbr = -1
      call agg_fill_nhbr_info(agg_nhbr, n_agg_facet, n_elements, &
        facet_neigh, offset_el, n_facet, &
        is_aggregated, aggregate_size)
    end if

    !> Allocate and fill lvl and nodes
    ntot = 0
    call tamg_lvl_init( tamg%lvl(lvl_id), lvl_id, naggs)
    do l = 1, naggs
      call tamg_node_init( tamg%lvl(lvl_id)%nodes(l), l, aggregate_size(l))
      j = 0
      do i = 1, n_elements!TODO: this is the lazy expensive way...
        if (is_aggregated(i) .eq. l) then
          j = j+1
          tamg%lvl(lvl_id)%nodes(l)%dofs(j) = i
        end if
      end do
      if (j .ne. tamg%lvl(lvl_id)%nodes(l)%ndofs) then
        print *, j, tamg%lvl(lvl_id)%nodes(l)%ndofs
        call neko_error("Aggregation problem. Not enough dofs in node.")
      end if
      ntot = ntot + aggregate_size(l)
    end do
    tamg%lvl(lvl_id)%fine_lvl_dofs = ntot
    allocate( tamg%lvl(lvl_id)%wrk_in( ntot ) )
    allocate( tamg%lvl(lvl_id)%wrk_out( ntot ) )
    print *, "work allocated on lvl", lvl_id, "dofs", ntot

    deallocate( is_aggregated )
    deallocate( aggregate_size )
  end subroutine aggregate_greedy

  !> Aggregates elements based on face-adjacent elements
  !! @param tamg TreeAMG hierarchy data structure being aggregated
  !! @param max_aggs Target number of aggregates to create on level
  !! @param agg_nhbr Output that tracks adjacency of aggregates for use on next level
  subroutine aggregate_elm(tamg, max_aggs, agg_nhbr)
    type(tamg_hierarchy_t), intent(inout) :: tamg
    integer, intent(in) :: max_aggs
    integer, intent(inout), allocatable :: agg_nhbr(:,:)
    integer :: naggs
    real(kind=dp) :: random_value
    integer :: i,j,k,l
    integer :: side, nhbr, n_elements
    integer :: lvl_id, ntot
    integer, allocatable :: is_aggregated(:)
    integer, allocatable :: aggregate_size(:)
    integer :: tnt_agg, tst_agg
    integer :: tnt_size, tst_size
    integer :: aa
    logical :: agg_added

    lvl_id = 2

    associate( msh => tamg%msh )

    naggs = 0
    n_elements = msh%nelv
    allocate( is_aggregated( n_elements ) )
    is_aggregated = -1!> fill with false
    allocate( aggregate_size( max_aggs*2 ) )
    aggregate_size = 999999!> Fill with large number
    !allocate( agg_nhbr(10, max_aggs*2) )
    allocate( agg_nhbr(20, max_aggs*2) )
    agg_nhbr = -1

    print *, "n_elements", msh%nelv, msh%offset_el
    do while (naggs .le. max_aggs)
      call random_number(random_value)
      i = floor(random_value * n_elements + 1)
      if (is_aggregated(i) .eq. -1) then!> if not aggregated, create new aggregate
        naggs = naggs + 1
        is_aggregated(i) = naggs
        aggregate_size(naggs) = 1
        !> Add neighbors to aggregate if unaggregated
        do side = 1, 6!> loop through neighbors
          nhbr = msh%facet_neigh(side, i) - msh%offset_el
          if ((nhbr .gt. 0).and.(nhbr .le. msh%nelv)) then!> if nhbr exists
            if (is_aggregated(nhbr) .eq. -1) then!> if nhbr unaggregated
              aggregate_size(naggs) = aggregate_size(naggs) + 1
              is_aggregated(nhbr) = naggs
            end if
          end if
        end do
      endif
    end do
    print *, "done with first pass of aggregation"

    !> Add remaining unaggregated nodes to aggregates
    do i = 1, n_elements
      if (is_aggregated(i) .eq. -1) then
        !> dof i is unaggregated. Check neighbors, add to smallest neighbor
        tnt_agg = -1
        tnt_size = 999!TODO: replace with large number
        tst_agg = -1
        tst_size = 999!TODO: replace with large number
        do side = 1, 6
          nhbr = msh%facet_neigh(side, i) - msh%offset_el
          if ((nhbr .gt. 0).and.(nhbr .le. msh%nelv)) then
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
          !> if neighbor aggregate found add to that aggregate
          is_aggregated(i) = tnt_agg
          aggregate_size(tnt_agg) = aggregate_size(tnt_agg) + 1
        else
          !> if none of the neignbors are aggregated. might as well make a new aggregate
          naggs = naggs + 1
          if (naggs .gt. max_aggs*2) then
            print *, "AGGS:", naggs
            call neko_error("I did not think of a way to handle creating too many aggregates... increase max_aggs")
          end if
          is_aggregated(i) = naggs
          aggregate_size(naggs) = 1
          !> Add neighbors to aggregate if unaggregated
          do side = 1, 6
            nhbr = msh%facet_neigh(side, i) - msh%offset_el
            if ((nhbr .gt. 0).and.(nhbr .le. msh%nelv)) then
              if (is_aggregated(nhbr) .eq. -1) then
                aggregate_size(naggs) = aggregate_size(naggs) + 1
                is_aggregated(nhbr) = naggs
              end if
            end if
          end do
        end if

      end if
    end do
    print *, "done with second pass of aggregation: number of aggregates", naggs
    do i = 1, n_elements!TODO: this is the lazy expensive way...
      tnt_agg = is_aggregated(i)
      do side = 1, 6
        nhbr = msh%facet_neigh(side,i) - msh%offset_el
        if ((nhbr .gt. 0).and.(nhbr .le. msh%nelv)) then
          tst_agg = is_aggregated(nhbr)
          if (tst_agg .le. 0) then
            call neko_error("Unaggregated element detected. We do not want to handle that here...")
          end if
          if (tst_agg .ne. tnt_agg) then
            agg_added = .false.
            do j = 1, 20
              if ((agg_nhbr(j,tnt_agg) .eq. tst_agg)) then
                agg_added = .true.
              else if ((agg_nhbr(j,tnt_agg).eq.-1).and.(.not.agg_added)) then
                agg_nhbr(j,tnt_agg) = tst_agg
                agg_added = .true.
              end if
            end do
          end if
        end if
      end do
    end do

    !> Allocate and fill lvl and nodes
    ntot = 0
    call tamg_lvl_init( tamg%lvl(lvl_id), lvl_id, naggs)
    do l = 1, naggs
      call tamg_node_init( tamg%lvl(lvl_id)%nodes(l), l, aggregate_size(l))
      j = 0
      do i = 1, n_elements!TODO: this is the lazy expensive way...
        if (is_aggregated(i) .eq. l) then
          j = j+1
          tamg%lvl(lvl_id)%nodes(l)%dofs(j) = i
        end if
      end do
      if (j .ne. tamg%lvl(lvl_id)%nodes(l)%ndofs) then
        print *, j, tamg%lvl(lvl_id)%nodes(l)%ndofs
        call neko_error("Aggregation problem. Not enough dofs in node.")
      end if
      ntot = ntot + aggregate_size(l)
    end do
    tamg%lvl(lvl_id)%fine_lvl_dofs = ntot
    allocate( tamg%lvl(lvl_id)%wrk_in( ntot ) )
    allocate( tamg%lvl(lvl_id)%wrk_out( ntot ) )
    print *, "work allocated on lvl", lvl_id, "dofs", ntot

    end associate
  end subroutine aggregate_elm

  !> Aggregates dofs based on adjacent dofs
  !! @param tamg TreeAMG hierarchy data structure being aggregated
  !! @param max_aggs Target number of aggregates to create on level
  !! @param lvl_id The level id for which aggregates are being created
  !! @param agg_nhbr Input array that contains adjacency of aggregates
  subroutine aggregate_general(tamg, max_aggs, lvl_id, agg_nhbr)
    type(tamg_hierarchy_t), intent(inout) :: tamg
    integer, intent(in) :: max_aggs
    integer, intent(in) :: lvl_id
    integer, intent(inout), allocatable :: agg_nhbr(:,:)
    integer :: naggs
    real(kind=dp) :: random_value
    integer :: i,j,k,l
    integer :: side, nhbr, n_elements
    integer :: ntot
    integer, allocatable :: is_aggregated(:)
    integer, allocatable :: aggregate_size(:)
    integer :: tnt_agg, tst_agg
    integer :: tnt_size, tst_size
    integer :: aa
    logical :: agg_added

    naggs = 0
    n_elements = tamg%lvl(lvl_id-1)%nnodes
    allocate( is_aggregated( n_elements ) )
    is_aggregated = -1!> fill with false
    allocate( aggregate_size( max_aggs*2 ) )
    aggregate_size = 999999!> Fill with large number
    !!allocate( agg_nhbr(10, max_aggs*2) )
    !!agg_nhbr = -1

    do while (naggs .le. max_aggs)
      call random_number(random_value)
      i = floor(random_value * n_elements + 1)
      if (is_aggregated(i) .eq. -1) then!> if not aggregated, create new aggregate
        naggs = naggs + 1
        is_aggregated(i) = naggs
        aggregate_size(naggs) = 1
        !> Add neighbors to aggregate if unaggregated
        do side = 1, 20!> loop through neighbors
          nhbr = agg_nhbr(side, i)
          if (nhbr .gt. 0) then!> if nhbr exists
            if (is_aggregated(nhbr) .eq. -1) then!> if nhbr unaggregated
              aggregate_size(naggs) = aggregate_size(naggs) + 1
              is_aggregated(nhbr) = naggs
            end if
          end if
        end do
      endif
    end do
    print *, "done with first pass of aggregation"

    !> Add remaining unaggregated nodes to aggregates
    do i = 1, n_elements
      if (is_aggregated(i) .eq. -1) then
        !> dof i is unaggregated. Check neighbors, add to smallest neighbor
        tnt_agg = -1
        tnt_size = 999!TODO: replace with large number
        tst_agg = -1
        tst_size = 999!TODO: replace with large number
        do side = 1, 20
          nhbr = agg_nhbr(side, i)
          if (nhbr .gt. 0) then
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
          !> if neighbor aggregate found add to that aggregate
          is_aggregated(i) = tnt_agg
          aggregate_size(tnt_agg) = aggregate_size(tnt_agg) + 1
        else
          !> if none of the neignbors are aggregated. might as well make a new aggregate
          naggs = naggs + 1
          if (naggs .gt. max_aggs*2) then
            print *, "AGGS:", naggs
            call neko_error("I did not think of a way to handle creating too many aggregates... increase max_aggs")
          end if
          is_aggregated(i) = naggs
          aggregate_size(naggs) = 1
          !> Add neighbors to aggregate if unaggregated
          do side = 1, 20
            nhbr = agg_nhbr(side, i)
            if (nhbr .gt. 0) then
              if (is_aggregated(nhbr) .eq. -1) then
                aggregate_size(naggs) = aggregate_size(naggs) + 1
                is_aggregated(nhbr) = naggs
              end if
            end if
          end do
        end if

      end if
    end do
    print *, "done with second pass of aggregation: number of aggregates", naggs
    !do i = 1, n_elements!TODO: this is the lazy expensive way...
    !  tnt_agg = is_aggregated(i)
    !  do side = 1, 10
    !    nhbr = msh_blah(side,i)
    !    tst_agg = is_aggregated(nhbr)
    !    if (tst_agg .ne. tnt_agg) then
    !      agg_added = .false.
    !      do j = 1, 10
    !        if ((agg_nhbr(tnt_agg,j) .eq. tst_agg)) then
    !          agg_added = .true.
    !        else if ((agg_nhbr(tnt_agg,j) .ne. tst_agg).and.(agg_nhbr(tnt_agg,j).eq.-1).and.(.not.agg_added)) then
    !          agg_nhbr(tnt_agg,j) = tst_agg
    !          agg_added = .true.
    !        end if
    !      end do
    !    end if
    !  end do
    !end do

    !> Allocate and fill lvl and nodes
    ntot = 0
    call tamg_lvl_init( tamg%lvl(lvl_id), lvl_id, naggs)
    do l = 1, naggs
      call tamg_node_init( tamg%lvl(lvl_id)%nodes(l), l, aggregate_size(l))
      j = 0
      do i = 1, n_elements!TODO: this is the lazy expensive way...
        if (is_aggregated(i) .eq. l) then
          j = j+1
          tamg%lvl(lvl_id)%nodes(l)%dofs(j) = tamg%lvl(lvl_id-1)%nodes(i)%gid
        end if
      end do
      if (j .ne. tamg%lvl(lvl_id)%nodes(l)%ndofs) then
        call neko_error("Aggregation problem. Not enough dofs in node.")
      end if
      ntot = ntot + aggregate_size(l)
    end do
    tamg%lvl(lvl_id)%fine_lvl_dofs = ntot
    allocate( tamg%lvl(lvl_id)%wrk_in( ntot ) )
    allocate( tamg%lvl(lvl_id)%wrk_out( ntot ) )
    print *, "work allocated on lvl", lvl_id, "dofs", ntot

  end subroutine aggregate_general

  !> Aggregate all dofs to a single point to form a tree-like structure.
  !! @param tamg TreeAMG hierarchy data structure being aggregated
  !! @param lvl_id The level id for which aggregates are being created
  subroutine aggregate_end( tamg, lvl_id)
    type(tamg_hierarchy_t), intent(inout) :: tamg
    integer, intent(in) :: lvl_id
    integer :: nt, i
    !> link all branches together at a point
    call tamg_lvl_init( tamg%lvl(lvl_id), lvl_id, 1)
    !> Allocate lvl
    nt = tamg%lvl(lvl_id-1)%nnodes
    tamg%lvl(lvl_id)%fine_lvl_dofs = nt
    allocate( tamg%lvl(lvl_id)%wrk_in( nt ) )
    allocate( tamg%lvl(lvl_id)%wrk_out( nt ) )
    print *, "work allocated on lvl", lvl_id, "dofs", nt

    !> Allocate node
    call tamg_node_init( tamg%lvl(lvl_id)%nodes(1), 1, nt)
    !> Fill node
    do i = 1, tamg%lvl(lvl_id-1)%nnodes
      tamg%lvl(lvl_id)%nodes(1)%dofs(i) = tamg%lvl(lvl_id-1)%nodes(i)%gid
    end do
  end subroutine aggregate_end
end module tree_amg_aggregate
