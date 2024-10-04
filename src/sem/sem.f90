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
!> Implements the `sem_t` type.
module sem
  use neko_config
  use mesh, only : mesh_t
  use space, only : space_t
  use gather_scatter, only : gs_t
  use coefs, only : coef_t
  use dofmap, only : dofmap_t
  use num_types, only : rp
  use utils, only : neko_error
  use file, only : file_t
  use mesh, only : mesh_t
  implicit none
  private

  !> Wrapper type containing the entire SEM backbone.
  type, public :: sem_t
     !> Approximation space.
     type(space_t) :: Xh 
     !> Computational mesh.
     type(mesh_t) :: msh
     !> Direct stiffness summation kernels (gather-scatter)
     !! Often denoted QQ^T in matrix notation in the literature.
     type(gs_t) :: gs
     !> Map of degrees of freedom.
     type(dofmap_t) :: dofmap
     !> Coeffients, including transormation metrics.
     type(coef_t) :: coef
   contains
     !> Constructor
     procedure, pass(this) :: init => sem_init
     !> Destructor
     procedure, pass(this) :: free => sem_free


  end type sem_t

contains
  !> Constructor.
  subroutine sem_init(this, mesh_file, quadrature, lx, ly, lz)
    class(sem_t), intent(inout) :: this
    character(len=*), intent(in) :: mesh_file
    integer, intent(in) :: quadrature   !< Quadrature type
    integer, intent(in) :: lx           !< Polynomial dimension in x-direction
    integer, intent(in) :: ly           !< Polynomial dimension in y-direction
    integer, intent(in) :: lz !< Polynomial dimension in z-direction

    type(file_t) :: msh_file

    call this%free()

    call this%Xh%init(quadrature, lx, ly, lz)
    msh_file = file_t(mesh_file)

    call msh_file%read(this%msh)

    call this%dofmap%init(this%msh, this%Xh)
    call this%gs%init(this%dofmap)
    call this%coef%init(this%gs)

  end subroutine sem_init

  !> Destructor.
  subroutine sem_free(this)
    class(sem_t), intent(inout) :: this
    
    call this%Xh%free()
    call this%msh%free()
    call this%dofmap%free()
    call this%coef%free()
    call this%gs%free()

  end subroutine sem_free

end module sem
