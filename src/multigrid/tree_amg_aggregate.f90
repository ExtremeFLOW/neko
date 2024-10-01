module tree_amg_aggregate
  use tree_amg
  use utils
  use num_types
  use mesh, only : mesh_t
  implicit none

contains

  ! Fill finest level
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
  end subroutine aggregate_finest_level

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
    allocate( agg_nhbr(10, max_aggs*2) )
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
            do j = 1, 10
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

    end associate
  end subroutine aggregate_elm

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
        do side = 1, 10!> loop through neighbors
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
        do side = 1, 10
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
          do side = 1, 10
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

  end subroutine aggregate_general

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

    !> Allocate node
    call tamg_node_init( tamg%lvl(lvl_id)%nodes(1), 1, nt)
    !> Fill node
    do i = 1, tamg%lvl(lvl_id-1)%nnodes
      tamg%lvl(lvl_id)%nodes(1)%dofs(i) = tamg%lvl(lvl_id-1)%nodes(i)%gid
    end do
  end subroutine aggregate_end
end module tree_amg_aggregate
