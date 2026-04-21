! Copyright (c) 2019-2026, The Neko Authors
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
!> Implementation for the nonconforming faces/edges during communication
module gs_interp_cpu
  use num_types, only : i4, i8, rp
  use mxm_wrapper, only: mxm
  use tensor, only : transpose
  use mesh_conn, only : mesh_conn_t
  use amr_interpolate, only : amr_interpolate_t, amr_nchildren
  use field, only : field_t
  use gs_interp, only : gs_interp_t
  use amr_reconstruct, only : amr_reconstruct_t

  implicit none
  private

  !> Type for face/edge interpolation
  type, public, extends(gs_interp_t) :: gs_interp_cpu_t
     !> Face interpolation operator
     real(rp), allocatable, dimension(:, :, :, :) :: jm_fcs
     !> Face inverse interpolation operator
     real(rp), allocatable, dimension(:, :, :, :) :: jm_fcsi
     !> Edge interpolation operator
     real(rp), allocatable, dimension(:, :, :) :: jm_edg
     !> Edge inverse interpolation operator
     real(rp), allocatable, dimension(:, :, :) :: jm_edgi
     !> Mask to zero children's faces/edges
     real(rp), allocatable, dimension(:, :, :, :) :: zero_msk
     !> Face global multiplicity for J^T operator (after action)
     real(rp), allocatable, dimension(:, :, :) :: mult_fcs_jt
     !> Inverse of face global multiplicity for J^T operator (before action)
     real(rp), allocatable, dimension(:, :, :) :: mult_fcs_jt_inv
     !> Face global multiplicity for J^-1 operator (after action)
     real(rp), allocatable, dimension(:, :, :) :: mult_fcs_ji
     !> Inverse of face global multiplicity for J^-1 operator (before action)
     real(rp), allocatable, dimension(:, :, :) :: mult_fcs_ji_inv
     !> Face global multiplicity for H1 operator (after action)
     real(rp), allocatable, dimension(:, :, :) :: mult_fcs_h1
     !> Inverse of face global multiplicity for H1 operator (before action)
     real(rp), allocatable, dimension(:, :, :) :: mult_fcs_h1_inv
     !> Edge global multiplicity for J^T operator (after action)
     real(rp), allocatable, dimension(:, :) :: mult_edg_jt
     !> Inverse of edge global multiplicity for J^T operator (before action)
     real(rp), allocatable, dimension(:, :) :: mult_edg_jt_inv
     !> Edge global multiplicity for J^-1 operator (after action)
     real(rp), allocatable, dimension(:, :) :: mult_edg_ji
     !> Inverse of edge global multiplicity for J^-1 operator (before action)
     real(rp), allocatable, dimension(:, :) :: mult_edg_ji_inv
     !> Edge global multiplicity for H1 operator (after action)
     real(rp), allocatable, dimension(:, :) :: mult_edg_h1
     !> Inverse of edg global multiplicity for H1 operator (before action)
     real(rp), allocatable, dimension(:, :) :: mult_edg_h1_inv
     !> Face work array
     real(rp), allocatable, dimension(:, :) :: face_tmp
     !> Edge work array
     real(rp), allocatable, dimension(:) :: edge_tmp
     !> Face storage for interpolation
     real(rp), allocatable, dimension(:, :, :) :: facein, faceout, facetmp
     !> Edge storage for interpolation
     real(rp), allocatable, dimension(:, :) :: edgein, edgeout
   contains
     !> Initialise type
     procedure, pass(this) :: init => gs_interp_cpu_init
     !> Initialise multiplicity arrays
     procedure, pass(this) :: init_mult => gs_interp_cpu_init_mult
     !> Free type
     procedure, pass(this) :: free => gs_interp_cpu_free
     !> Perform face/edge interpolation
     procedure, pass(this) :: apply_j => gs_interp_cpu_apply_j
     !> Perform inverse face/edge interpolation
     procedure, pass(this) :: apply_ji => gs_interp_cpu_apply_ji
     !> Perform transposed face/edge interpolation
     procedure, pass(this) :: apply_jt => gs_interp_cpu_apply_jt
     !> Zero children faces/edges
     procedure, pass(this) :: zero_children => gs_interp_cpu_zero_children
     !> Remove multiplicity for J^T
     procedure, pass(this) :: remove_mult_jt => gs_interp_cpu_remove_mult_jt
     !> Remove multiplicity for J^-1
     procedure, pass(this) :: remove_mult_ji => gs_interp_cpu_remove_mult_ji
     !> Remove multiplicity for H1
     procedure, pass(this) :: remove_mult_h1 => gs_interp_cpu_remove_mult_h1
     !> Add multiplicity for J^T
     procedure, pass(this) :: add_mult_jt => gs_interp_cpu_add_mult_jt
     !> Add multiplicity for J^-1
     procedure, pass(this) :: add_mult_ji => gs_interp_cpu_add_mult_ji
     !> Add multiplicity for H1
     procedure, pass(this) :: add_mult_h1 => gs_interp_cpu_add_mult_h1
     !> AMR restart
     procedure, pass(this) :: amr_restart => gs_interp_cpu_amr_restart
  end type gs_interp_cpu_t

contains
  !> Initialise gs interpolation type
  !! @param[in]  lx    polynomial order + 1
  !! @param[in]  conn  mesh connectivity
  subroutine gs_interp_cpu_init(this, lx, conn)
    class(gs_interp_cpu_t), intent(inout) :: this
    integer, intent(in) :: lx
    type(mesh_conn_t), target, intent(in) :: conn

    call this%free()

    call this%init_base(lx, conn)

    ! Work space
    allocate(this%face_tmp(this%lx, this%lx), this%edge_tmp(this%lx))
    if (this%nhang_fcs .gt. 0) &
         allocate(this%facein(this%lx, this%lx, this%nhang_fcs), &
         this%faceout(this%lx, this%lx, this%nhang_fcs), &
         this%facetmp(this%lx, this%lx, this%nhang_fcs))
    if (this%nhang_edg .gt. 0) &
         allocate(this%edgein(this%lx, this%nhang_edg), &
         this%edgeout(this%lx, this%nhang_edg))

    ! interpolation operators
    call gs_interp_cpu_init_operator(this)

    ! multiplicity is set from gs_t, as it requires complete gs setup

  end subroutine gs_interp_cpu_init

  !> Initialise face/edge interpolation operators
  subroutine gs_interp_cpu_init_operator(this)
    type(gs_interp_cpu_t), intent(inout) :: this
    integer :: il, jl, lposx, lposy, itmp

    if (this%ifhang) then
       ! face operators
       if (this%nhang_fcs .gt. 0) then
          allocate(this%jm_fcs(this%lx, this%lx, 2, this%nhang_fcs), &
               this%jm_fcsi(this%lx, this%lx, 2, this%nhang_fcs))

          do il =1, this%nhang_fcs
             lposx = mod(this%hang_fcs_pos(il), 2) + 1
             lposy = this%hang_fcs_pos(il) / 2 + 1
             select case(this%hang_fcs(il))
             case(1, 3)
                this%jm_fcs(:, :, 1, il) = this%interpolate%x_cr2fn(:, :, lposx)
                this%jm_fcs(:, :, 2, il) = this%interpolate%z_cr2fnT(:, :, &
                     lposy)
                this%jm_fcsi(:, :, 1, il) = this%interpolate%x_fn2cr(:, :, &
                     lposx)
                this%jm_fcsi(:, :, 2, il) = this%interpolate%z_fn2crT(:, :, &
                     lposy)
             case(2, 4)
                this%jm_fcs(:, :, 1, il) = this%interpolate%y_cr2fn(:, :, lposx)
                this%jm_fcs(:, :, 2, il) = this%interpolate%z_cr2fnT(:, :, &
                     lposy)
                this%jm_fcsi(:, :, 1, il) = this%interpolate%y_fn2cr(:, :, &
                     lposx)
                this%jm_fcsi(:, :, 2, il) = this%interpolate%z_fn2crT(:, :, &
                     lposy)
             case(5, 6)
                this%jm_fcs(:, :, 1, il) = this%interpolate%x_cr2fn(:, :, lposx)
                this%jm_fcs(:, :, 2, il) = this%interpolate%y_cr2fnT(:, :, &
                     lposy)
                this%jm_fcsi(:, :, 1, il) = this%interpolate%x_fn2cr(:, :, &
                     lposx)
                this%jm_fcsi(:, :, 2, il) = this%interpolate%y_fn2crT(:, :, &
                     lposy)
             end select
          end do
       end if

       ! edge operators
       if (this%nhang_edg .gt. 0) then
          allocate(this%jm_edg(this%lx, this%lx, this%nhang_edg), &
               this%jm_edgi(this%lx, this%lx, this%nhang_edg))

          do il =1, this%nhang_edg
             lposx = this%hang_edg_pos(il) + 1
             ! face unrelated edges only
             if (lposx .eq. 1 .or. lposx .eq. 2) then
                select case(this%hang_edg(il))
                case(1: 4)
                   this%jm_edg(:, :, il) = this%interpolate%x_cr2fn(:, :, lposx)
                   this%jm_edgi(:, :, il) = this%interpolate%x_fn2cr(:, :, &
                        lposx)
                case(5: 8)
                   this%jm_edg(:, :, il) = this%interpolate%y_cr2fn(:, :, lposx)
                   this%jm_edgi(:, :, il) = this%interpolate%y_fn2cr(:, :, &
                        lposx)
                case(9: 12)
                   this%jm_edg(:, :, il) = this%interpolate%z_cr2fn(:, :, lposx)
                   this%jm_edgi(:, :, il) = this%interpolate%z_fn2cr(:, :, &
                        lposx)
                end select
             end if
          end do
       end if

       ! Set mask to remove children's faces/edges
       allocate(this%zero_msk(this%lx, this%lx, this%lx, this%nhang_el))
       this%zero_msk(:, :, :, :) = 1.0_rp
       this%face_tmp(:, :) = 0.0_rp
       this%edge_tmp(:) = 0.0_rp
       do il = 1, this%nhang_el
          ! faces
          itmp = this%hang_fcs_off(il + 1) - this%hang_fcs_off(il)
          if (itmp .gt. 0) then
             do jl = this%hang_fcs_off(il), this%hang_fcs_off(il + 1) - 1
                call vector_to_face(this%zero_msk(:, :, :, il), this%face_tmp, &
                     this%hang_fcs(jl), this%lx)
             end do
          end if

          ! edges
          itmp = this%hang_edg_off(il + 1) - this%hang_edg_off(il)
          if (itmp .gt. 0) then
             do jl = this%hang_edg_off(il), this%hang_edg_off(il + 1) - 1
                call vector_to_edge(this%zero_msk(:, :, :, il), this%edge_tmp, &
                     this%hang_edg(jl), this%lx)
             end do
          end if
       end do

    end if

  end subroutine gs_interp_cpu_init_operator

  !> Initialise multiplicity arrays
  !! @param[in]  mult_jt  multiplicity array for J^T for the whole mesh
  !! @param[in]  mult_ji  multiplicity array for J^-1 for the whole mesh
  !! @param[in]  mult_h1  multiplicity array for H1 for the whole mesh
  subroutine gs_interp_cpu_init_mult(this, mult_jt, mult_ji, mult_h1)
    class(gs_interp_cpu_t), intent(inout) :: this
    real(rp), dimension(:, :, :, :) , intent(in) :: mult_jt, mult_ji, mult_h1
    integer :: il, jl, kl, lposx, lposy, itmp
    real(rp), parameter :: one = 1.0_rp

    associate(lx => this%lx, nhang_fcs => this%nhang_fcs, &
         nhang_edg => this%nhang_edg)
      if (this%ifhang) then
         ! face operators
         if (nhang_fcs .gt. 0) then
            ! Allocate face multiplicity arrays
            allocate(this%mult_fcs_jt(lx, lx, nhang_fcs), &
                 this%mult_fcs_jt_inv(lx, lx, nhang_fcs), &
                 this%mult_fcs_ji(lx, lx, nhang_fcs), &
                 this%mult_fcs_ji_inv(lx, lx, nhang_fcs), &
                 this%mult_fcs_h1(lx, lx, nhang_fcs), &
                 this%mult_fcs_h1_inv(lx, lx, nhang_fcs))

            ! extract face multiplicity
            do il = 1, this%nhang_el
               itmp = this%hang_fcs_off(il + 1) - this%hang_fcs_off(il)
               if (itmp .gt. 0) then
                  do jl = this%hang_fcs_off(il), this%hang_fcs_off(il + 1) - 1
                     call face_to_vector(mult_jt(:, :, :, this%hang_el(il)), &
                          this%mult_fcs_jt(:, :, jl), this%hang_fcs(jl), lx)
                     call face_to_vector(mult_ji(:, :, :, this%hang_el(il)), &
                          this%mult_fcs_ji(:, :, jl), this%hang_fcs(jl), lx)
                     call face_to_vector(mult_h1(:, :, :, this%hang_el(il)), &
                          this%mult_fcs_h1(:, :, jl), this%hang_fcs(jl), lx)
                  end do
               end if
            end do
            ! get inverse multiplicity before operator action
            do concurrent(il = 1 : lx, jl = 1 : lx, kl = 1: nhang_fcs)
               this%mult_fcs_jt_inv(il, jl, kl) = one / &
                    this%mult_fcs_jt(il, jl, kl)
               this%mult_fcs_ji_inv(il, jl, kl) = one / &
                    this%mult_fcs_ji(il, jl, kl)
               this%mult_fcs_h1_inv(il, jl, kl) = one / &
                    this%mult_fcs_h1(il, jl, kl)
            end do
            ! get multiplicity after operator action
            ! this does not work for lx = 2
            do il = 1, nhang_fcs
               lposx = mod(this%hang_fcs_pos(il), 2) + 1
               lposy = this%hang_fcs_pos(il) / 2 + 1
               call mult_fill_fcs(lx, lposx, lposy, this%mult_fcs_jt(:, :, il))
               call mult_fill_fcs(lx, lposx, lposy, this%mult_fcs_ji(:, :, il))
               call mult_fill_fcs(lx, lposx, lposy, this%mult_fcs_h1(:, :, il))
            end do
         end if

         ! edge operators
         if (nhang_edg .gt. 0) then
            ! Allocate face multiplicity arrays
            allocate(this%mult_edg_jt(lx, nhang_edg), &
                 this%mult_edg_jt_inv(lx, nhang_edg), &
                 this%mult_edg_ji(lx, nhang_edg), &
                 this%mult_edg_ji_inv(lx, nhang_edg), &
                 this%mult_edg_h1(lx, nhang_edg), &
                 this%mult_edg_h1_inv(lx, nhang_edg))

            ! extract edge multiplicity
            do il = 1, this%nhang_el
               itmp = this%hang_edg_off(il + 1) - this%hang_edg_off(il)
               if (itmp .gt. 0) then
                  do jl = this%hang_edg_off(il), this%hang_edg_off(il + 1) - 1
                     call edge_to_vector(mult_jt(:, :, :, this%hang_el(il)), &
                          this%mult_edg_jt(:, jl), this%hang_edg(jl), lx)
                     call edge_to_vector(mult_ji(:, :, :, this%hang_el(il)), &
                          this%mult_edg_ji(:, jl), this%hang_edg(jl), lx)
                     call edge_to_vector(mult_h1(:, :, :, this%hang_el(il)), &
                          this%mult_edg_h1(:, jl), this%hang_edg(jl), lx)
                  end do
               end if
            end do
            ! get inverse multiplicity before operator action
            do concurrent(il = 1 : lx, kl = 1: nhang_edg)
               this%mult_edg_jt_inv(il, kl) = one / this%mult_edg_jt(il, kl)
               this%mult_edg_ji_inv(il, kl) = one / this%mult_edg_ji(il, kl)
               this%mult_edg_h1_inv(il, kl) = one / this%mult_edg_h1(il, kl)
            end do
            ! get multiplicity after operator action
            ! this does not work for lx = 2
            do il = 1, nhang_edg
               lposx = this%hang_edg_pos(il) + 1
               call mult_fill_edg(lx, lposx, this%mult_edg_jt(:, il))
               call mult_fill_edg(lx, lposx, this%mult_edg_ji(:, il))
               call mult_fill_edg(lx, lposx, this%mult_edg_h1(:, il))
            end do
         end if
      end if

    end associate

  end subroutine gs_interp_cpu_init_mult

  !> Fill face multiplicity array
  subroutine mult_fill_fcs(lx, lposx, lposy, mult)
    integer, intent(in) :: lx, lposx, lposy
    real(rp), dimension(lx, lx), intent(inout) :: mult
    real(rp) :: rtmp

    ! Vertex doesn't change; correct edges and face
    select case (lposx)
    case (1)
       select case (lposy)
       case (1)
          rtmp = mult(1, 2) ! edge 1
          mult(1, 2 : lx)  = rtmp
          rtmp = mult(2, 1) ! edge 2
          mult(2 : lx, 1)  = rtmp
          rtmp = mult(2, 2) ! face
          mult(2 : lx, 2 : lx)  = rtmp
       case (2)
          rtmp = mult(1, lx -1) ! edge 1
          mult(1, 1 : lx - 1)  = rtmp
          rtmp = mult(2, lx) ! edge 2
          mult(2 : lx, lx)  = rtmp
          rtmp = mult(2, lx - 1) ! face
          mult(2 : lx, 1 : lx - 1)  = rtmp
       end select
    case (2)
       select case (lposy)
       case (1)
          rtmp = mult(lx, 2) ! edge 1
          mult(lx, 2 : lx)  = rtmp
          rtmp = mult(lx - 1, 1) ! edge 2
          mult(1 : lx - 1, 1)  = rtmp
          rtmp = mult(lx - 1, 2) ! face
          mult(1 : lx - 1, 2 : lx)  = rtmp
       case (2)
          rtmp = mult(lx, lx - 1) ! edge 1
          mult(lx, 1 : lx - 1)  = rtmp
          rtmp = mult(lx - 1, lx) ! edge 2
          mult(1 : lx - 1, lx)  = rtmp
          rtmp = mult(lx - 1, lx - 1) ! face
          mult(1 : lx - 1, 1 : lx - 1)  = rtmp
       end select
    end select
  end subroutine mult_fill_fcs

  !> Fill edge multiplicity array
  subroutine mult_fill_edg(lx, lposx, mult)
    integer, intent(in) :: lx, lposx
    real(rp), dimension(lx), intent(inout) :: mult
    real(rp) :: rtmp

    ! Vertex doesn't change; correct edge
    select case (lposx)
    case (1)
       rtmp = mult(2) ! edge
       mult(2 : lx)  = rtmp
    case (2)
       rtmp = mult(lx -1) ! edge
       mult(1 : lx - 1)  = rtmp
    end select
  end subroutine mult_fill_edg

  !> Free gs interpolation type
  subroutine gs_interp_cpu_free(this)
    class(gs_interp_cpu_t), intent(inout) :: this

    call this%free_base()

    if (allocated(this%jm_fcs)) deallocate(this%jm_fcs)
    if (allocated(this%jm_fcsi)) deallocate(this%jm_fcsi)
    if (allocated(this%jm_edg)) deallocate(this%jm_edg)
    if (allocated(this%jm_edgi)) deallocate(this%jm_edgi)
    if (allocated(this%zero_msk)) deallocate(this%zero_msk)
    if (allocated(this%mult_fcs_jt)) deallocate(this%mult_fcs_jt)
    if (allocated(this%mult_fcs_jt_inv)) deallocate(this%mult_fcs_jt_inv)
    if (allocated(this%mult_fcs_ji)) deallocate(this%mult_fcs_ji)
    if (allocated(this%mult_fcs_ji_inv)) deallocate(this%mult_fcs_ji_inv)
    if (allocated(this%mult_fcs_h1)) deallocate(this%mult_fcs_h1)
    if (allocated(this%mult_fcs_h1_inv)) deallocate(this%mult_fcs_h1_inv)
    if (allocated(this%mult_edg_jt)) deallocate(this%mult_edg_jt)
    if (allocated(this%mult_edg_jt_inv)) deallocate(this%mult_edg_jt_inv)
    if (allocated(this%mult_edg_ji)) deallocate(this%mult_edg_ji)
    if (allocated(this%mult_edg_ji_inv)) deallocate(this%mult_edg_ji_inv)
    if (allocated(this%mult_edg_h1)) deallocate(this%mult_edg_h1)
    if (allocated(this%mult_edg_h1_inv)) deallocate(this%mult_edg_h1_inv)
    if (allocated(this%face_tmp)) deallocate(this%face_tmp)
    if (allocated(this%edge_tmp)) deallocate(this%edge_tmp)
    if (allocated(this%facein)) deallocate(this%facein)
    if (allocated(this%faceout)) deallocate(this%faceout)
    if (allocated(this%facetmp)) deallocate(this%facetmp)
    if (allocated(this%edgein)) deallocate(this%edgein)
    if (allocated(this%edgeout)) deallocate(this%edgeout)

  end subroutine gs_interp_cpu_free

  !> Perform face/edge interpolation
  !! @param[inout]  field    field for face interpolation
  subroutine gs_interp_cpu_apply_j(this, field)
    class(gs_interp_cpu_t), intent(inout) :: this
    type(field_t), intent(inout) :: field
    integer :: il, itmp

    if (this%ifhang) then
       do il = 1, this%nhang_el
          ! face interpolation
          itmp = this%hang_fcs_off(il + 1) - this%hang_fcs_off(il)
          if (itmp .gt. 0) call gs_interp_cpu_apply_j_elem_face(this%lx, &
                  field%x(:, :, :, this%hang_el(il)), itmp, &
                  this%hang_fcs(this%hang_fcs_off(il) : &
                  this%hang_fcs_off(il + 1) - 1), &
                  this%jm_fcs(:, :, :, this%hang_fcs_off(il) : &
                  this%hang_fcs_off(il + 1) - 1), this%facein(:, :, 1 : itmp), &
                  this%faceout(:, :, 1 : itmp), this%facetmp(:, :, 1 : itmp))

          ! edge interpolation
          itmp = this%hang_edg_off(il + 1) - this%hang_edg_off(il)
          if (itmp .gt. 0) call gs_interp_cpu_apply_j_elem_edge(this%lx, &
                  field%x(:, :, :, this%hang_el(il)), itmp, &
                  this%hang_edg(this%hang_edg_off(il) : &
                  this%hang_edg_off(il + 1) - 1), &
                  this%jm_edg(:, :, this%hang_edg_off(il) : &
                  this%hang_edg_off(il + 1) - 1), this%edgein(:, 1 : itmp), &
                  this%edgeout(:, 1 : itmp))
       end do
    end if

  end subroutine gs_interp_cpu_apply_j

  !> Perform face interpolation in a single element for j
  !! @param[in]     lx        number of points in 1D
  !! @param[inout]  elem      field element
  !! @param[in]     nface     face number
  !! @param[in]     facelist  list of faces
  !! @param[in]     operator  interpolation operators
  !! @param[inout]  facein, faceout, facetmp   work arrays
  subroutine gs_interp_cpu_apply_j_elem_face(lx, elem, nface, facelist, &
       operator, facein, faceout, facetmp)
    integer, intent(in) :: lx, nface
    real(rp), dimension(lx, lx, lx), intent(inout) :: elem
    integer, dimension(nface), intent(in) :: facelist
    real(rp), dimension(lx, lx, 2, nface), intent(in) :: operator
    real(rp), dimension(lx, lx, nface), intent(inout) :: facein, faceout, &
         facetmp
    integer :: il

    ! extract faces
    do il = 1, nface
       call face_to_vector(elem, facein(:, :, il), facelist(il), lx)
    end do

    ! interpolate faces
    do il = 1, nface
       call mxm(operator(:, :, 1, il), lx, facein(:, :, il), lx, &
            facetmp(:, :, il), lx)
       call mxm(facetmp(:, :, il) , lx, operator(:, :, 2, il), lx, &
            faceout(:, :, il), lx)
    end do

    ! put faces back
    do il = 1, nface
       call vector_to_face(elem, faceout(:, :, il), facelist(il), lx)
    end do

  end subroutine gs_interp_cpu_apply_j_elem_face

  !> Perform edge interpolation in a single element for j
  !! @param[in]     lx        number of points in 1D
  !! @param[inout]  elem      field element
  !! @param[in]     nedge     edge number
  !! @param[in]     edgelist  list of edges
  !! @param[in]     operator  interpolation operators
  !! @param[inout]  edgein, edgeout    work arrays
  subroutine gs_interp_cpu_apply_j_elem_edge(lx, elem, nedge, edgelist, &
       operator, edgein, edgeout)
    integer, intent(in) :: lx, nedge
    real(rp), dimension(lx, lx, lx), intent(inout) :: elem
    integer, dimension(nedge), intent(in) :: edgelist
    real(rp), dimension(lx, lx, nedge), intent(in) :: operator
    real(rp), dimension(lx, 1, nedge), intent(inout) :: edgein, edgeout
    integer :: il

    ! extract edges
    do il = 1, nedge
       call edge_to_vector(elem, edgein(:, 1, il), edgelist(il), lx)
    end do

    ! interpolate edges
    do il = 1, nedge
       call mxm(operator(:, :, il), lx, edgein(:, :, il), lx, &
            edgeout(:, :, il), 1)
    end do

    ! put edges back
    do il = 1, nedge
       call vector_to_edge(elem, edgeout(:, 1, il), edgelist(il), lx)
    end do

  end subroutine gs_interp_cpu_apply_j_elem_edge

  !> Perform inverse face/edge interpolation
  !! @param[inout]  field    field for face interpolation
  subroutine gs_interp_cpu_apply_ji(this, field)
    class(gs_interp_cpu_t), intent(inout) :: this
    type(field_t), intent(inout) :: field
    integer :: il, itmp

    if (this%ifhang) then
       do il = 1, this%nhang_el
          ! face interpolation
          itmp = this%hang_fcs_off(il + 1) - this%hang_fcs_off(il)
          if (itmp .gt. 0) call gs_interp_cpu_apply_ji_elem_face(this%lx, &
                  field%x(:, :, :, this%hang_el(il)), itmp, &
                  this%hang_fcs(this%hang_fcs_off(il) : &
                  this%hang_fcs_off(il + 1) - 1), &
                  this%jm_fcsi(:, :, :, this%hang_fcs_off(il) : &
                  this%hang_fcs_off(il + 1) - 1), this%interpolate%fc_mult, &
                  this%facein(:, :, 1 : itmp), this%faceout(:, :, 1 : itmp), &
                  this%facetmp(:, :, 1 : itmp))

          ! edge interpolation
          itmp = this%hang_edg_off(il + 1) - this%hang_edg_off(il)
          if (itmp .gt. 0) call gs_interp_cpu_apply_ji_elem_edge(this%lx, &
                  field%x(:, :, :, this%hang_el(il)), itmp, &
                  this%hang_edg(this%hang_edg_off(il) : &
                  this%hang_edg_off(il + 1) - 1), &
                  this%jm_edgi(:, :, this%hang_edg_off(il) : &
                  this%hang_edg_off(il + 1) - 1), this%interpolate%ed_mult, &
                  this%edgein(:, 1 : itmp), this%edgeout(:, 1 : itmp))
       end do
    end if

  end subroutine gs_interp_cpu_apply_ji

  !> Perform inverse face interpolation in a single element for j
  !! @param[in]     lx        number of points in 1D
  !! @param[inout]  elem      field element
  !! @param[in]     nface     face number
  !! @param[in]     facelist  list of faces
  !! @param[in]     operator  interpolation operators
  !! @param[in]     mult      multiplicity
  !! @param[inout]  facein, faceout, facetmp   work arrays
  subroutine gs_interp_cpu_apply_ji_elem_face(lx, elem, nface, facelist, &
       operator, mult, facein, faceout, facetmp)
    integer, intent(in) :: lx, nface
    real(rp), dimension(lx, lx, lx), intent(inout) :: elem
    integer, dimension(nface), intent(in) :: facelist
    real(rp), dimension(lx, lx, 2, nface), intent(in) :: operator
    real(rp), dimension(lx, lx), intent(in) :: mult
    real(rp), dimension(lx, lx, nface), intent(inout) :: facein, faceout, &
         facetmp
    integer :: il

    ! extract faces
    do il = 1, nface
       call face_to_vector(elem, facein(:, :, il), facelist(il), lx)
    end do

    ! interpolate faces
    do il = 1, nface
       call mxm(operator(:, :, 1, il), lx, facein(:, :, il), lx, &
            facetmp(:, :, il), lx)
       call mxm(facetmp(:, :, il) , lx, operator(:, :, 2, il), lx, &
            faceout(:, :, il), lx)
       faceout(:, :, il) = faceout(:, :, il) * mult(:, :)
    end do

    ! put faces back
    do il = 1, nface
       call vector_to_face(elem, faceout(:, :, il), facelist(il), lx)
    end do

  end subroutine gs_interp_cpu_apply_ji_elem_face

  !> Perform inverse edge interpolation in a single element for j
  !! @param[in]     lx        number of points in 1D
  !! @param[inout]  elem      field element
  !! @param[in]     nedge     edge number
  !! @param[in]     edgelist  list of edges
  !! @param[in]     operator  interpolation operators
  !! @param[in]     mult      multiplicity
  !! @param[inout]  edgein, edgeout    work arrays
  subroutine gs_interp_cpu_apply_ji_elem_edge(lx, elem, nedge, edgelist, &
       operator, mult, edgein, edgeout)
    integer, intent(in) :: lx, nedge
    real(rp), dimension(lx, lx, lx), intent(inout) :: elem
    integer, dimension(nedge), intent(in) :: edgelist
    real(rp), dimension(lx, lx, nedge), intent(in) :: operator
    real(rp), dimension(lx), intent(in) :: mult
    real(rp), dimension(lx, 1, nedge), intent(inout) :: edgein, edgeout
    integer :: il

    ! extract edges
    do il = 1, nedge
       call edge_to_vector(elem, edgein(:, 1, il), edgelist(il), lx)
    end do

    ! interpolate edges
    do il = 1, nedge
       call mxm(operator(:, :, il), lx, edgein(:, :, il), lx, &
            edgeout(:, :, il), 1)
       edgeout(:, 1, il) = edgeout(:, 1, il) * mult(:)
    end do

    ! put edges back
    do il = 1, nedge
       call vector_to_edge(elem, edgeout(:, 1, il), edgelist(il), lx)
    end do

  end subroutine gs_interp_cpu_apply_ji_elem_edge

  !> Perform transposed face/edge interpolation
  !! @param[inout]  field    field for face interpolation
  subroutine gs_interp_cpu_apply_jt(this, field)
    class(gs_interp_cpu_t), intent(inout) :: this
    type(field_t), intent(inout) :: field
    integer :: il, itmp

    if (this%ifhang) then
       do il = 1, this%nhang_el
          ! face interpolation
          itmp = this%hang_fcs_off(il + 1) - this%hang_fcs_off(il)
          if (itmp .gt. 0) call gs_interp_cpu_apply_jt_elem_face(this%lx, &
                  field%x(:, :, :, this%hang_el(il)), itmp, &
                  this%hang_fcs(this%hang_fcs_off(il) : &
                  this%hang_fcs_off(il + 1) - 1), &
                  this%jm_fcs(:, :, :, this%hang_fcs_off(il) : &
                  this%hang_fcs_off(il + 1) - 1), this%facein(:, :, 1 : itmp), &
                  this%faceout(:, :, 1 : itmp), this%facetmp(:, :, 1 : itmp))

          ! edge interpolation
          itmp = this%hang_edg_off(il + 1) - this%hang_edg_off(il)
          if (itmp .gt. 0) call gs_interp_cpu_apply_jt_elem_edge(this%lx, &
                  field%x(:, :, :, this%hang_el(il)), itmp, &
                  this%hang_edg(this%hang_edg_off(il) : &
                  this%hang_edg_off(il + 1) - 1), &
                  this%jm_edg(:, :, this%hang_edg_off(il) : &
                  this%hang_edg_off(il + 1) - 1), this%edgein(:, 1 : itmp), &
                  this%edgeout(:, 1 : itmp))
       end do
    end if

  end subroutine gs_interp_cpu_apply_jt

  !> Perform face interpolation in a single element for j transposed
  !! @param[in]     lx        number of points in 1D
  !! @param[inout]  elem      field element
  !! @param[in]     nface     face number
  !! @param[in]     facelist  list of faces
  !! @param[in]     operator  interpolation operators
  !! @param[inout]  facein, faceout, facetmp   work arrays
  subroutine gs_interp_cpu_apply_jt_elem_face(lx, elem, nface, facelist, &
       operator, facein, faceout, facetmp)
    integer, intent(in) :: lx, nface
    real(rp), dimension(lx, lx, lx), intent(inout) :: elem
    integer, dimension(nface), intent(in) :: facelist
    real(rp), dimension(lx, lx, 2, nface), intent(in) :: operator
    real(rp), dimension(lx, lx, nface), intent(inout) :: facein, faceout, &
         facetmp
    integer :: il

    ! extract faces
    ! NOTICE ORIGINAL ftovec_0l
    do il = 1, nface
       call face_to_vector(elem, facein(:, :, il), facelist(il), lx)
    end do

    ! interpolate faces
    do il = 1, nface
       call transpose(facein(:, :, il), lx)
       call mxm(operator(:, :, 2, il), lx, facein(:, :, il), lx, &
            facetmp(:, :, il), lx)
       call mxm(facetmp(:, :, il) , lx, operator(:, :, 1, il), lx, &
            faceout(:, :, il), lx)
       call transpose(faceout(:, :, il), lx)
    end do

    ! put faces back
    ! NOTICE ORIGINAL vectof_addl
    do il = 1, nface
       call vector_to_face(elem, faceout(:, :, il), facelist(il), lx)
    end do

  end subroutine gs_interp_cpu_apply_jt_elem_face

  !> Perform edge interpolation in a single element for j transposed
  !! @param[in]     lx        number of points in 1D
  !! @param[inout]  elem      field element
  !! @param[in]     nedge     edge number
  !! @param[in]     edgelist  list of edges
  !! @param[in]     operator  interpolation operators
  !! @param[inout]  edgein, edgeout    work arrays
  subroutine gs_interp_cpu_apply_jt_elem_edge(lx, elem, nedge, edgelist, &
       operator, edgein, edgeout)
    integer, intent(in) :: lx, nedge
    real(rp), dimension(lx, lx, lx), intent(inout) :: elem
    integer, dimension(nedge), intent(in) :: edgelist
    real(rp), dimension(lx, lx, nedge), intent(in) :: operator
    real(rp), dimension(1, lx, nedge), intent(inout) :: edgein, edgeout
    integer :: il

    ! extract edges
    do il = 1, nedge
       call edge_to_vector(elem, edgein(1, :, il), edgelist(il), lx)
    end do

    ! interpolate edges
    do il = 1, nedge
       call mxm(edgein(:, :, il), 1, operator(:, :, il), lx, &
            edgeout(:, :, il), lx)
    end do

    ! put edges back
    do il = 1, nedge
       call vector_to_edge(elem, edgeout(1, :, il), edgelist(il), lx)
    end do

  end subroutine gs_interp_cpu_apply_jt_elem_edge

  !> Zero children's faces/edges
  !! @param[inout]  field    field for face interpolation
  subroutine gs_interp_cpu_zero_children(this, field)
    class(gs_interp_cpu_t), intent(inout) :: this
    type(field_t), intent(inout) :: field
    integer :: il

    if (this%ifhang) then
       do il = 1, this%nhang_el
          field%x(:, :, :, this%hang_el(il)) = &
               field%x(:, :, :, this%hang_el(il)) * this%zero_msk(:, :, :, il)
       end do
    end if

  end subroutine gs_interp_cpu_zero_children

  !> Remove multiplicity for J^T
  !! @param[inout]  field    field for face interpolation
  subroutine gs_interp_cpu_remove_mult_jt(this, field)
    class(gs_interp_cpu_t), intent(inout) :: this
    type(field_t), intent(inout) :: field
    integer :: il, jl, itmp

    if (this%ifhang) then
       do il = 1, this%nhang_el
          ! face multiplicity
          itmp = this%hang_fcs_off(il + 1) - this%hang_fcs_off(il)
          if (itmp .gt. 0) then
             do jl = this%hang_fcs_off(il), this%hang_fcs_off(il + 1) - 1
                call face_to_vector(field%x(:, :, :, this%hang_el(il)), &
                     this%face_tmp, this%hang_fcs(jl), this%lx)
                this%face_tmp(:, :) = this%face_tmp(:, :) * &
                     this%mult_fcs_jt_inv(:, :, jl)
                call vector_to_face(field%x(:, :, :, this%hang_el(il)), &
                     this%face_tmp, this%hang_fcs(jl), this%lx)
             end do
          end if

          ! edge multiplicity
          itmp = this%hang_edg_off(il + 1) - this%hang_edg_off(il)
          if (itmp .gt. 0) then
             do jl = this%hang_edg_off(il), this%hang_edg_off(il + 1) - 1
                call edge_to_vector(field%x(:, :, :, this%hang_el(il)), &
                     this%edge_tmp, this%hang_edg(jl), this%lx)
                this%edge_tmp(:) = this%edge_tmp(:) * &
                     this%mult_edg_jt_inv(:, jl)
                call vector_to_edge(field%x(:, :, :, this%hang_el(il)), &
                     this%edge_tmp, this%hang_edg(jl), this%lx)
             end do
          end if
       end do
    end if

  end subroutine gs_interp_cpu_remove_mult_jt

  !> Remove multiplicity for J^-1
  !! @param[inout]  field    field for face interpolation
  subroutine gs_interp_cpu_remove_mult_ji(this, field)
    class(gs_interp_cpu_t), intent(inout) :: this
    type(field_t), intent(inout) :: field
    integer :: il, jl, itmp

    if (this%ifhang) then
       do il = 1, this%nhang_el
          ! face multiplicity
          itmp = this%hang_fcs_off(il + 1) - this%hang_fcs_off(il)
          if (itmp .gt. 0) then
             do jl = this%hang_fcs_off(il), this%hang_fcs_off(il + 1) - 1
                call face_to_vector(field%x(:, :, :, this%hang_el(il)), &
                     this%face_tmp, this%hang_fcs(jl), this%lx)
                this%face_tmp(:, :) = this%face_tmp(:, :) * &
                     this%mult_fcs_ji_inv(:, :, jl)
                call vector_to_face(field%x(:, :, :, this%hang_el(il)), &
                     this%face_tmp, this%hang_fcs(jl), this%lx)
             end do
          end if

          ! edge multiplicity
          itmp = this%hang_edg_off(il + 1) - this%hang_edg_off(il)
          if (itmp .gt. 0) then
             do jl = this%hang_edg_off(il), this%hang_edg_off(il + 1) - 1
                call edge_to_vector(field%x(:, :, :, this%hang_el(il)), &
                     this%edge_tmp, this%hang_edg(jl), this%lx)
                this%edge_tmp(:) = this%edge_tmp(:) * &
                     this%mult_edg_ji_inv(:, jl)
                call vector_to_edge(field%x(:, :, :, this%hang_el(il)), &
                     this%edge_tmp, this%hang_edg(jl), this%lx)
             end do
          end if
       end do
    end if

  end subroutine gs_interp_cpu_remove_mult_ji

  !> Remove multiplicity for H1
  !! @param[inout]  field    field for face interpolation
  subroutine gs_interp_cpu_remove_mult_h1(this, field)
    class(gs_interp_cpu_t), intent(inout) :: this
    type(field_t), intent(inout) :: field
    integer :: il, jl, itmp

    if (this%ifhang) then
       do il = 1, this%nhang_el
          ! face multiplicity
          itmp = this%hang_fcs_off(il + 1) - this%hang_fcs_off(il)
          if (itmp .gt. 0) then
             do jl = this%hang_fcs_off(il), this%hang_fcs_off(il + 1) - 1
                call face_to_vector(field%x(:, :, :, this%hang_el(il)), &
                     this%face_tmp, this%hang_fcs(jl), this%lx)
                this%face_tmp(:, :) = this%face_tmp(:, :) * &
                     this%mult_fcs_h1_inv(:, :, jl)
                call vector_to_face(field%x(:, :, :, this%hang_el(il)), &
                     this%face_tmp, this%hang_fcs(jl), this%lx)
             end do
          end if

          ! edge multiplicity
          itmp = this%hang_edg_off(il + 1) - this%hang_edg_off(il)
          if (itmp .gt. 0) then
             do jl = this%hang_edg_off(il), this%hang_edg_off(il + 1) - 1
                call edge_to_vector(field%x(:, :, :, this%hang_el(il)), &
                     this%edge_tmp, this%hang_edg(jl), this%lx)
                this%edge_tmp(:) = this%edge_tmp(:) * &
                     this%mult_edg_h1_inv(:, jl)
                call vector_to_edge(field%x(:, :, :, this%hang_el(il)), &
                     this%edge_tmp, this%hang_edg(jl), this%lx)
             end do
          end if
       end do
    end if

  end subroutine gs_interp_cpu_remove_mult_h1

  !> Add multiplicity for J^T
  !! @param[inout]  field    field for face interpolation
  subroutine gs_interp_cpu_add_mult_jt(this, field)
    class(gs_interp_cpu_t), intent(inout) :: this
    type(field_t), intent(inout) :: field
    integer :: il, jl, itmp

    if (this%ifhang) then
       do il = 1, this%nhang_el
          ! face multiplicity
          itmp = this%hang_fcs_off(il + 1) - this%hang_fcs_off(il)
          if (itmp .gt. 0) then
             do jl = this%hang_fcs_off(il), this%hang_fcs_off(il + 1) - 1
                call face_to_vector(field%x(:, :, :, this%hang_el(il)), &
                     this%face_tmp, this%hang_fcs(jl), this%lx)
                this%face_tmp(:, :) = this%face_tmp(:, :) * &
                     this%mult_fcs_jt(:, :, jl)
                call vector_to_face(field%x(:, :, :, this%hang_el(il)), &
                     this%face_tmp, this%hang_fcs(jl), this%lx)
             end do
          end if

          ! edge multiplicity
          itmp = this%hang_edg_off(il + 1) - this%hang_edg_off(il)
          if (itmp .gt. 0) then
             do jl = this%hang_edg_off(il), this%hang_edg_off(il + 1) - 1
                call edge_to_vector(field%x(:, :, :, this%hang_el(il)), &
                     this%edge_tmp, this%hang_edg(jl), this%lx)
                this%edge_tmp(:) = this%edge_tmp(:) * &
                     this%mult_edg_jt(:, jl)
                call vector_to_edge(field%x(:, :, :, this%hang_el(il)), &
                     this%edge_tmp, this%hang_edg(jl), this%lx)
             end do
          end if
       end do
    end if

  end subroutine gs_interp_cpu_add_mult_jt

  !> Add multiplicity for J^-1
  !! @param[inout]  field    field for face interpolation
  subroutine gs_interp_cpu_add_mult_ji(this, field)
    class(gs_interp_cpu_t), intent(inout) :: this
    type(field_t), intent(inout) :: field
    integer :: il, jl, itmp

    if (this%ifhang) then
       do il = 1, this%nhang_el
          ! face multiplicity
          itmp = this%hang_fcs_off(il + 1) - this%hang_fcs_off(il)
          if (itmp .gt. 0) then
             do jl = this%hang_fcs_off(il), this%hang_fcs_off(il + 1) - 1
                call face_to_vector(field%x(:, :, :, this%hang_el(il)), &
                     this%face_tmp, this%hang_fcs(jl), this%lx)
                this%face_tmp(:, :) = this%face_tmp(:, :) * &
                     this%mult_fcs_ji(:, :, jl)
                call vector_to_face(field%x(:, :, :, this%hang_el(il)), &
                     this%face_tmp, this%hang_fcs(jl), this%lx)
             end do
          end if

          ! edge multiplicity
          itmp = this%hang_edg_off(il + 1) - this%hang_edg_off(il)
          if (itmp .gt. 0) then
             do jl = this%hang_edg_off(il), this%hang_edg_off(il + 1) - 1
                call edge_to_vector(field%x(:, :, :, this%hang_el(il)), &
                     this%edge_tmp, this%hang_edg(jl), this%lx)
                this%edge_tmp(:) = this%edge_tmp(:) * &
                     this%mult_edg_ji(:, jl)
                call vector_to_edge(field%x(:, :, :, this%hang_el(il)), &
                     this%edge_tmp, this%hang_edg(jl), this%lx)
             end do
          end if
       end do
    end if

  end subroutine gs_interp_cpu_add_mult_ji

  !> Add multiplicity for H1
  !! @param[inout]  field    field for face interpolation
  subroutine gs_interp_cpu_add_mult_h1(this, field)
    class(gs_interp_cpu_t), intent(inout) :: this
    type(field_t), intent(inout) :: field
    integer :: il, jl, itmp

    if (this%ifhang) then
       do il = 1, this%nhang_el
          ! face multiplicity
          itmp = this%hang_fcs_off(il + 1) - this%hang_fcs_off(il)
          if (itmp .gt. 0) then
             do jl = this%hang_fcs_off(il), this%hang_fcs_off(il + 1) - 1
                call face_to_vector(field%x(:, :, :, this%hang_el(il)), &
                     this%face_tmp, this%hang_fcs(jl), this%lx)
                this%face_tmp(:, :) = this%face_tmp(:, :) * &
                     this%mult_fcs_h1(:, :, jl)
                call vector_to_face(field%x(:, :, :, this%hang_el(il)), &
                     this%face_tmp, this%hang_fcs(jl), this%lx)
             end do
          end if

          ! edge multiplicity
          itmp = this%hang_edg_off(il + 1) - this%hang_edg_off(il)
          if (itmp .gt. 0) then
             do jl = this%hang_edg_off(il), this%hang_edg_off(il + 1) - 1
                call edge_to_vector(field%x(:, :, :, this%hang_el(il)), &
                     this%edge_tmp, this%hang_edg(jl), this%lx)
                this%edge_tmp(:) = this%edge_tmp(:) * &
                     this%mult_edg_h1(:, jl)
                call vector_to_edge(field%x(:, :, :, this%hang_el(il)), &
                     this%edge_tmp, this%hang_edg(jl), this%lx)
             end do
          end if
       end do
    end if

  end subroutine gs_interp_cpu_add_mult_h1

  !> AMR restart
  !! @param[inout]  reconstruct   data reconstruction type
  !! @param[in]     counter       restart counter
  !! @param[in]     tstep         time step
  subroutine gs_interp_cpu_amr_restart(this, reconstruct, counter, tstep)
    class(gs_interp_cpu_t), intent(inout) :: this
    type(amr_reconstruct_t), intent(inout) :: reconstruct
    integer, intent(in) :: counter, tstep

    ! Was this component already restarted?
    if (this%counter .eq. counter) return

    this%counter = counter

    if (allocated(this%jm_fcs)) deallocate(this%jm_fcs)
    if (allocated(this%jm_fcsi)) deallocate(this%jm_fcsi)
    if (allocated(this%jm_edg)) deallocate(this%jm_edg)
    if (allocated(this%jm_edgi)) deallocate(this%jm_edgi)
    if (allocated(this%zero_msk)) deallocate(this%zero_msk)
    if (allocated(this%mult_edg_jt)) deallocate(this%mult_edg_jt)
    if (allocated(this%mult_edg_jt_inv)) deallocate(this%mult_edg_jt_inv)
    if (allocated(this%mult_edg_ji)) deallocate(this%mult_edg_ji)
    if (allocated(this%mult_edg_ji_inv)) deallocate(this%mult_edg_ji_inv)
    if (allocated(this%mult_edg_h1)) deallocate(this%mult_edg_h1)
    if (allocated(this%mult_edg_h1_inv)) deallocate(this%mult_edg_h1_inv)
    if (allocated(this%facein)) deallocate(this%facein)
    if (allocated(this%faceout)) deallocate(this%faceout)
    if (allocated(this%facetmp)) deallocate(this%facetmp)
    if (allocated(this%edgein)) deallocate(this%edgein)
    if (allocated(this%edgeout)) deallocate(this%edgeout)

    call this%amr_restart_base()

    ! Work space
    if (this%nhang_fcs .gt. 0) &
         allocate(this%facein(this%lx, this%lx, this%nhang_fcs), &
         this%faceout(this%lx, this%lx, this%nhang_fcs), &
         this%facetmp(this%lx, this%lx, this%nhang_fcs))
    if (this%nhang_edg .gt. 0) &
         allocate(this%edgein(this%lx, this%nhang_edg), &
         this%edgeout(this%lx, this%nhang_edg))

    ! interpolation operators
    call gs_interp_cpu_init_operator(this)

    ! multiplicity will be reconstructed from gs_t, as communicator has to be
    ! finalised at this point

  end subroutine gs_interp_cpu_amr_restart

  ! What follows is hex specific with assumption all three dimensions are
  ! equal. Not sure this should not be located here.
  ! Just a temporary solution

  !> Extract face from the hex element; symmetric notation
  !> @param[in]     elem    element vector
  !> @param[out]    face    element face
  !> @param[in]     iface   face number
  !> @param[in]     nx      1D size
  pure subroutine face_to_vector(elem, face, iface, nx)
    integer, intent(in) :: iface, nx
    real(rp), dimension(nx, nx, nx), intent(in) :: elem
    real(rp), dimension(nx, nx), intent(out) :: face

    select case(iface)
    case(1)
       face(:, :) = elem(1, :, :)
    case(2)
       face(:, :) = elem(nx, :, :)
    case(3)
       face(:, :) = elem(:, 1, :)
    case(4)
       face(:, :) = elem(:, nx, :)
    case(5)
       face(:, :) = elem(:, :, 1)
    case(6)
       face(:, :) = elem(:, :, nx)
    end select

  end subroutine face_to_vector

  !> Fill face of the hex element with data; symmetric notation
  !> @param[out]    elem    element vector
  !> @param[in]     face    element face
  !> @param[in]     iface   face number
  !> @param[in]     nx      1D size
  pure subroutine vector_to_face(elem, face, iface, nx)
    integer, intent(in) :: iface, nx
    real(rp), dimension(nx, nx, nx), intent(out) :: elem
    real(rp), dimension(nx, nx), intent(in) :: face

    select case(iface)
    case(1)
       elem(1, :, :) = face(:, :)
    case(2)
       elem(nx, :, :) = face(:, :)
    case(3)
       elem(:, 1, :) = face(:, :)
    case(4)
       elem(:, nx, :) = face(:, :)
    case(5)
       elem(:, :, 1) = face(:, :)
    case(6)
       elem(:, :, nx) = face(:, :)
    end select

  end subroutine vector_to_face

  !> Extract edge from the hex element; symmetric notation
  !> @param[in]     elem    element vector
  !> @param[out]    edge    element edge
  !> @param[in]     iedge   edge number
  !> @param[in]     nx      1D size
  pure subroutine edge_to_vector(elem, edge, iedge, nx)
    integer, intent(in) :: iedge, nx
    real(rp), dimension(nx, nx, nx), intent(in) :: elem
    real(rp), dimension(nx), intent(out) :: edge

    select case(iedge)
    case(1)
       edge(:) = elem(:, 1, 1)
    case(2)
       edge(:) = elem(:, nx, 1)
    case(3)
       edge(:) = elem(:, 1, nx)
    case(4)
       edge(:) = elem(:, nx, nx)
    case(5)
       edge(:) = elem(1, :, 1)
    case(6)
       edge(:) = elem(nx, :, 1)
    case(7)
       edge(:) = elem(1, :, nx)
    case(8)
       edge(:) = elem(nx, :, nx)
    case(9)
       edge(:) = elem(1, 1, :)
    case(10)
       edge(:) = elem(nx, 1, :)
    case(11)
       edge(:) = elem(1, nx, :)
    case(12)
       edge(:) = elem(nx, nx, :)
    end select

  end subroutine edge_to_vector

  !> Fill edge of the hex element with data; symmetric notation
  !> @param[out]    elem    element vector
  !> @param[in]     edge    element edge
  !> @param[in]     iedge   edge number
  !> @param[in]     nx      1D size
  pure subroutine vector_to_edge(elem, edge, iedge, nx)
    integer, intent(in) :: iedge, nx
    real(rp), dimension(nx, nx, nx), intent(out) :: elem
    real(rp), dimension(nx), intent(in) :: edge

    select case(iedge)
    case(1)
       elem(:, 1, 1) = edge(:)
    case(2)
       elem(:, nx, 1) = edge(:)
    case(3)
       elem(:, 1, nx) = edge(:)
    case(4)
       elem(:, nx, nx) = edge(:)
    case(5)
       elem(1, :, 1) = edge(:)
    case(6)
       elem(nx, :, 1) = edge(:)
    case(7)
       elem(1, :, nx) = edge(:)
    case(8)
       elem(nx, :, nx) = edge(:)
    case(9)
       elem(1, 1, :) = edge(:)
    case(10)
       elem(nx, 1, :) = edge(:)
    case(11)
       elem(1, nx, :) = edge(:)
    case(12)
       elem(nx, nx, :) = edge(:)
    end select

  end subroutine vector_to_edge

end module gs_interp_cpu
