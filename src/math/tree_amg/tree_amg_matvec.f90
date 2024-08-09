module tree_amg_matvec
  use tree_amg
  use num_types
  use gather_scatter, only : gs_t, GS_OP_ADD
  implicit none

contains

  subroutine restriction_operator(vec_out, vec_in, lvl, tamg)
    real(kind=rp), intent(inout) :: vec_out(:)
    real(kind=rp), intent(inout) :: vec_in(:)
    integer, intent(in) :: lvl
    type(tamg_hierarchy_t), intent(inout) :: tamg
    integer :: i, n

    vec_out = 0d0
    do n = 1, tamg%lvl(lvl)%nnodes
      associate (node => tamg%lvl(lvl)%nodes(n))
      do i = 1, node%ndofs
        vec_out( node%gid ) = vec_out( node%gid ) + vec_in( node%dofs(i) ) * node%interp_r( i )
      end do
      end associate
    end do
  end subroutine restriction_operator

  subroutine prolongation_operator(vec_out, vec_in, lvl, tamg)
    real(kind=rp), intent(inout) :: vec_out(:)
    real(kind=rp), intent(inout) :: vec_in(:)
    integer, intent(in) :: lvl
    type(tamg_hierarchy_t), intent(inout) :: tamg
    integer :: i, n

    vec_out = 0d0
    do n = 1, tamg%lvl(lvl)%nnodes
      associate (node => tamg%lvl(lvl)%nodes(n))
      do i = 1, node%ndofs
        vec_out( node%dofs(i) ) = vec_out( node%dofs(i) ) + vec_in( node%gid ) * node%interp_p( i )
      end do
      end associate
    end do
  end subroutine prolongation_operator


  recursive subroutine matvec_operator(vec_out, vec_in, lvl, lvl_out, tamg)
    real(kind=rp), intent(inout) :: vec_out(:)
    real(kind=rp), intent(inout) :: vec_in(:)
    integer, intent(in) :: lvl
    integer, intent(in) :: lvl_out
    type(tamg_hierarchy_t), intent(inout) :: tamg
    integer :: i, n, e

    vec_out = 0d0

    if (lvl .eq. 0) then !> isleaf true
      !> Call local finite element assembly
      call tamg%gs_h%op(vec_in, size(vec_in), GS_OP_ADD)
      do i = 1, size(vec_in)
        vec_in(i) = vec_in(i) * tamg%coef%mult(i,1,1,1)
      end do
      !>
      call tamg%ax%compute(vec_out, vec_in, tamg%coef, tamg%msh, tamg%Xh)
      !>
      call tamg%gs_h%op(vec_out, size(vec_out), GS_OP_ADD)
      !>
    else !> pass down through hierarchy
      if (lvl_out .ge. lvl) then
        !> lvl is finer than desired output
        !> project input vector to finer grid
        associate( wrk_in => tamg%lvl(lvl)%wrk_in, wrk_out => tamg%lvl(lvl)%wrk_out)
        wrk_in = 0d0
        wrk_out = 0d0
        do n = 1, tamg%lvl(lvl)%nnodes
          associate (node => tamg%lvl(lvl)%nodes(n))
          do i = 1, node%ndofs
            wrk_in( node%dofs(i) ) = wrk_in( node%dofs(i) ) + vec_in( node%gid ) * node%interp_p( i )
          end do
          end associate
        end do

        call matvec_operator(wrk_out, wrk_in, lvl-1, lvl_out, tamg)

        !> restrict to coarser grid
        do n = 1, tamg%lvl(lvl)%nnodes
          associate (node => tamg%lvl(lvl)%nodes(n))
          do i = 1, node%ndofs
            vec_out( node%gid ) = vec_out(node%gid ) + wrk_out( node%dofs(i) ) * node%interp_r( i )
          end do
          end associate
        end do
        end associate
      else if (lvl_out .lt. lvl) then
        !> lvl is coarser then desired output. Continue down tree
        call matvec_operator(vec_out, vec_in, lvl-1, lvl_out, tamg)
      else
        print *, "ERROR"
      end if
    end if
  end subroutine matvec_operator
end module tree_amg_matvec
