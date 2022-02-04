! Copyright (c) 2019-2021, The Neko Authors
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
!> Defines a mesh field
!! @details A mesh field is a scalar integer cell based field (\f$ dQ_0 \f$)
module mesh_field
  use num_types
  use mesh
  implicit none
  
  !> @todo Add support for different data types
  type mesh_fld_t
     integer, allocatable :: data(:) !< Data
     type(mesh_t), pointer :: msh    !< Mesh
     character(len=80) :: name
  end type mesh_fld_t

contains
  
  subroutine mesh_field_init(fld, msh, fld_name)
    type(mesh_fld_t), intent(inout) :: fld
    type(mesh_t), target, intent(in) :: msh
    character(len=*), optional :: fld_name 

    call mesh_field_free(fld)

    fld%msh => msh
    if (.not. allocated(fld%data)) then
       allocate(fld%data(msh%nelv))
    end if

    if (present(fld_name)) then
       fld%name = fld_name
    else
       fld%name = 'MeshField'
    end if

    fld%data = 0
  end subroutine mesh_field_init

  subroutine mesh_field_free(fld)
    type(mesh_fld_t), intent(inout) :: fld

    if (allocated(fld%data)) then
       deallocate(fld%data)
    end if

    nullify(fld%msh)
  end subroutine mesh_field_free

end module mesh_field
