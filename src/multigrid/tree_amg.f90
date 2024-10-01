module tree_amg
  use num_types
  use utils
  use math
  use coefs, only : coef_t
  use mesh, only : mesh_t
  use space, only : space_t
  use ax_product, only: ax_t
  use bc, only: bc_list_t, bc_list_apply
  use gather_scatter, only : gs_t, GS_OP_ADD
  implicit none
  private

  type, private :: tamg_node_t
    logical :: isleaf = .true.
    integer :: gid = -1
    integer :: lvl = -1
    integer :: ndofs = 0
    integer, allocatable :: dofs(:)
    real(kind=rp) :: xyz(3)
    real(kind=rp), allocatable :: interp_r(:)
    real(kind=rp), allocatable :: interp_p(:)
  end type tamg_node_t

  type, private :: tamg_lvl_t
    integer :: lvl = -1
    integer :: nnodes = 0
    type(tamg_node_t), allocatable :: nodes(:)
    integer :: fine_lvl_dofs = 0
    real(kind=rp), allocatable :: wrk_in(:)
    real(kind=rp), allocatable :: wrk_out(:)
  end type tamg_lvl_t

  type, public :: tamg_hierarchy_t
    integer :: nlvls
    type(tamg_lvl_t), allocatable :: lvl(:)

    !Things needed to do finest level matvec
    class(ax_t), pointer :: ax
    type(mesh_t), pointer :: msh
    type(space_t), pointer :: Xh
    type(coef_t), pointer :: coef
    type(gs_t), pointer :: gs_h
    type(bc_list_t), pointer :: blst

  contains
    procedure, pass(this) :: init => tamg_init
    procedure, pass(this) :: matvec => tamg_matvec
    procedure, pass(this) :: matvec_impl => tamg_matvec_impl
    procedure, pass(this) :: interp_f2c => tamg_restriction_operator
    procedure, pass(this) :: interp_c2f => tamg_prolongation_operator
  end type tamg_hierarchy_t

contains

  subroutine tamg_init(this, ax, Xh, coef, msh, gs_h, nlvls, blst)
    class(tamg_hierarchy_t), target, intent(inout) :: this
    class(ax_t), target, intent(in) :: ax
    type(space_t),target, intent(in) :: Xh
    type(coef_t), target, intent(in) :: coef
    type(mesh_t), target, intent(in) :: msh
    type(gs_t), target, intent(in) :: gs_h
    integer, intent(in) :: nlvls
    type(bc_list_t), target, intent(in) :: blst
    integer, allocatable :: agg_nhbr(:,:)

    this%ax => ax
    this%msh => msh
    this%Xh => Xh
    this%coef => coef
    this%gs_h => gs_h
    this%blst => blst

    this%nlvls = nlvls
    allocate( this%lvl(this%nlvls) )

    call aggregate_finest_level(this, Xh%lx, Xh%ly, Xh%lz, msh%nelv)
    print *, "Calling lazy aggregation"
    print *, "-- target aggregates:", (msh%nelv/8)
    call aggregate_elm(this, (msh%nelv/8), agg_nhbr)
    print *, "-- Aggregation done. Aggregates:", this%lvl(2)%nnodes

    if (nlvls .gt. 3) then
      print *, "Calling lazy aggregation"
      print *, "-- target aggregates:", (this%lvl(2)%nnodes/8)
      call aggregate_general(this, (this%lvl(2)%nnodes/8), 3, agg_nhbr)
      print *, "-- Aggregation done. Aggregates:", this%lvl(3)%nnodes
    end if

    call aggregate_end(this, nlvls)

  end subroutine tamg_init

  subroutine tamg_lvl_init(tamg_lvl, lvl, nnodes)
    type(tamg_lvl_t), intent(inout) :: tamg_lvl
    integer, intent(in) :: lvl
    integer, intent(in) :: nnodes

    tamg_lvl%lvl = lvl
    tamg_lvl%nnodes = nnodes
    allocate( tamg_lvl%nodes(tamg_lvl%nnodes) )
  end subroutine tamg_lvl_init

  subroutine tamg_node_init(node, gid, ndofs)
    type(tamg_node_t), intent(inout) :: node
    integer, intent(in) :: gid
    integer, intent(in) :: ndofs

    node%gid = gid
    node%ndofs = ndofs
    allocate( node%dofs( node%ndofs) )
    node%dofs = -1
    allocate( node%interp_r( node%ndofs) )
    allocate( node%interp_p( node%ndofs) )
    node%interp_r = 1.0
    node%interp_p = 1.0
  end subroutine tamg_node_init

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

  recursive subroutine tamg_matvec(this, vec_out, vec_in, lvl_out)
    class(tamg_hierarchy_t), intent(inout) :: this
    real(kind=rp), intent(inout) :: vec_out(:)
    real(kind=rp), intent(inout) :: vec_in(:)
    integer, intent(in) :: lvl_out
    integer :: i, n, e
    call this%matvec_impl(vec_out, vec_in, this%nlvls, lvl_out)
    !call tamg_matvec_flat_impl(this, vec_out, vec_in, this%nlvls, lvl_out)
  end subroutine tamg_matvec

  recursive subroutine tamg_matvec_impl(this, vec_out, vec_in, lvl, lvl_out)
    class(tamg_hierarchy_t), intent(inout) :: this
    real(kind=rp), intent(inout) :: vec_out(:)
    real(kind=rp), intent(inout) :: vec_in(:)
    integer, intent(in) :: lvl
    integer, intent(in) :: lvl_out
    integer :: i, n, e

    vec_out = 0d0

    if (lvl .eq. 0) then !> isleaf true
      n = size(vec_in)
      !> Call local finite element assembly
      call this%gs_h%op(vec_in, n, GS_OP_ADD)
      do i = 1, n
        vec_in(i) = vec_in(i) * this%coef%mult(i,1,1,1)
      end do
      !>
      call this%ax%compute(vec_out, vec_in, this%coef, this%msh, this%Xh)
      !>
      call this%gs_h%op(vec_out, n, GS_OP_ADD)
      call bc_list_apply(this%blst, vec_out, n)
      !>
    else !> pass down through hierarchy
      if (lvl_out .ge. lvl) then
        !> lvl is finer than desired output
        !> project input vector to finer grid
        associate( wrk_in => this%lvl(lvl)%wrk_in, wrk_out => this%lvl(lvl)%wrk_out)
        wrk_in = 0d0
        wrk_out = 0d0
        do n = 1, this%lvl(lvl)%nnodes
          associate (node => this%lvl(lvl)%nodes(n))
          do i = 1, node%ndofs
            wrk_in( node%dofs(i) ) = wrk_in( node%dofs(i) ) + vec_in( node%gid ) * node%interp_p( i )
          end do
          end associate
        end do

        call this%matvec_impl(wrk_out, wrk_in, lvl-1, lvl_out)

        !> restrict to coarser grid
        do n = 1, this%lvl(lvl)%nnodes
          associate (node => this%lvl(lvl)%nodes(n))
          do i = 1, node%ndofs
            vec_out( node%gid ) = vec_out(node%gid ) + wrk_out( node%dofs(i) ) * node%interp_r( i )
          end do
          end associate
        end do
        end associate
      else if (lvl_out .lt. lvl) then
        !> lvl is coarser then desired output. Continue down tree
        call this%matvec_impl(vec_out, vec_in, lvl-1, lvl_out)
      else
        print *, "ERROR"
      end if
    end if
  end subroutine tamg_matvec_impl


  recursive subroutine tamg_matvec_flat_impl(this, vec_out, vec_in, lvl_blah, lvl_out)
    class(tamg_hierarchy_t), intent(inout) :: this
    real(kind=rp), intent(inout) :: vec_out(:)
    real(kind=rp), intent(inout) :: vec_in(:)
    integer, intent(in) :: lvl_blah
    integer, intent(in) :: lvl_out
    integer :: i, n, e, nn, lvl

    vec_out = 0d0
    do lvl = 1, lvl_out+1
      associate( wrk_in => this%lvl(lvl)%wrk_in, wrk_out => this%lvl(lvl)%wrk_out)
        wrk_in = 0d0
        wrk_out = 0d0
      end associate
    end do

    !> copy to make things easier to think about
    do n = 1, this%lvl(lvl_out+1)%nnodes
      do i = 1, this%lvl(lvl_out+1)%nodes(n)%ndofs
        associate( wrk_in => this%lvl(lvl_out+1)%wrk_in, node => this%lvl(lvl_out+1)%nodes(n))
          wrk_in( node%dofs(i) ) = vec_in( node%dofs(i) )
        end associate
      end do!i
    end do!n

    if(lvl_out .gt. 0) then
      do n = 1, this%lvl(lvl_out)%nnodes!> this loop is independent
        do lvl = lvl_out, 1, -1

          associate( wrk_in => this%lvl(lvl)%wrk_in, wrk_out => this%lvl(lvl)%wrk_out)
          do nn = 1, this%lvl(lvl)%nnodes
            associate (node => this%lvl(lvl)%nodes(nn))
            do i = 1, node%ndofs
              wrk_in( node%dofs(i) ) = wrk_in( node%dofs(i) ) + this%lvl(lvl+1)%wrk_in( node%gid ) * node%interp_p( i )
            end do!i
            end associate
          end do!nn
          end associate

        end do!lvl
      end do!n
    end if

    associate( wrk_in => this%lvl(1)%wrk_in, wrk_out => this%lvl(1)%wrk_out)
    !> Do finest level matvec
    n = size(wrk_in)
    !> Call local finite element assembly
    call this%gs_h%op(wrk_in, n, GS_OP_ADD)
    do i = 1, n
      wrk_in(i) = wrk_in(i) * this%coef%mult(i,1,1,1)
    end do
    !>
    call this%ax%compute(wrk_out, wrk_in, this%coef, this%msh, this%Xh)
    !>
    call this%gs_h%op(wrk_out, n, GS_OP_ADD)
    call bc_list_apply(this%blst, wrk_out, n)
    !>
    end associate


    if(lvl_out .gt. 0) then
      do n = 1, this%lvl(lvl_out)%nnodes!> this loop is independent
        do lvl = 2, lvl_out+1

          associate( wrk_in => this%lvl(lvl)%wrk_in, wrk_out => this%lvl(lvl)%wrk_out)
          do nn = 1, this%lvl(lvl)%nnodes
            associate (node => this%lvl(lvl)%nodes(nn))
            do i = 1, node%ndofs
              wrk_out( node%gid ) = wrk_out( node%gid ) + this%lvl(lvl-1)%wrk_out( node%dofs(i) ) * node%interp_r( i )
            end do!i
            end associate
          end do!nn
          end associate

        end do!lvl
      end do!n
    end if

    !> copy here to make things easier to think about
    do n = 1, this%lvl(lvl_out+1)%nnodes
      do i = 1, this%lvl(lvl_out+1)%nodes(n)%ndofs
        associate( wrk_out => this%lvl(lvl_out+1)%wrk_out, node => this%lvl(lvl_out+1)%nodes(n))
          vec_out( node%dofs(i) ) = wrk_out( node%dofs(i) )
        end associate
      end do!i
    end do!n
  end subroutine tamg_matvec_flat_impl


  subroutine tamg_restriction_operator(this, vec_out, vec_in, lvl)
    class(tamg_hierarchy_t), intent(inout) :: this
    real(kind=rp), intent(inout) :: vec_out(:)
    real(kind=rp), intent(inout) :: vec_in(:)
    integer, intent(in) :: lvl
    integer :: i, n

    vec_out = 0d0
    do n = 1, this%lvl(lvl)%nnodes
      associate (node => this%lvl(lvl)%nodes(n))
      do i = 1, node%ndofs
        vec_out( node%gid ) = vec_out( node%gid ) + vec_in( node%dofs(i) ) * node%interp_r( i )
      end do
      end associate
    end do
  end subroutine tamg_restriction_operator

  subroutine tamg_prolongation_operator(this, vec_out, vec_in, lvl)
    class(tamg_hierarchy_t), intent(inout) :: this
    real(kind=rp), intent(inout) :: vec_out(:)
    real(kind=rp), intent(inout) :: vec_in(:)
    integer, intent(in) :: lvl
    integer :: i, n

    vec_out = 0d0
    do n = 1, this%lvl(lvl)%nnodes
      associate (node => this%lvl(lvl)%nodes(n))
      do i = 1, node%ndofs
        vec_out( node%dofs(i) ) = vec_out( node%dofs(i) ) + vec_in( node%gid ) * node%interp_p( i )
      end do
      end associate
    end do
  end subroutine tamg_prolongation_operator

end module tree_amg
