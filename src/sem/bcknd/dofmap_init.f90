! Copyright (c) 2020-2026, The Neko Authors
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
!> Submodule containing the mapping based initialization of dofmaps
submodule (dofmap) dofmap_init_m
  use interpolation, only : interpolator_t

contains

  !> Constructor.
  !! @param dof The existing dofmap to initialize from.
  !! @param Xh The SEM function space.
  module subroutine dofmap_init_and_map(this, dof, Xh)
    class(dofmap_t), intent(inout) :: this
    type(dofmap_t), intent(in) :: dof
    type(space_t), intent(in) :: Xh
    type(interpolator_t) :: interpolator

    ! Initialize as usual
    call this%init_from_mesh(dof%msh, Xh)

    ! Interpolate if needed
    if (dof%Xh%lxyz .ne. this%Xh%lxyz) then
       call interpolator%init(this%Xh, dof%Xh)

       if (NEKO_BCKND_DEVICE .eq. 1) then
          call interpolator%map(this%x_d, dof%x_d, this%msh%nelv, this%Xh)
          call interpolator%map(this%y_d, dof%y_d, this%msh%nelv, this%Xh)
          call interpolator%map(this%z_d, dof%z_d, this%msh%nelv, this%Xh)
       else
          call interpolator%map_host(this%x, dof%x, this%msh%nelv, this%Xh)
          call interpolator%map_host(this%y, dof%y, this%msh%nelv, this%Xh)
          call interpolator%map_host(this%z, dof%z, this%msh%nelv, this%Xh)
       end if

       call interpolator%free()

    else
       if (NEKO_BCKND_DEVICE .eq. 1) then
          call device_copy(this%x_d, dof%x_d, this%ntot)
          call device_copy(this%y_d, dof%y_d, this%ntot)
          call device_copy(this%z_d, dof%z_d, this%ntot)
       else
          call copy(this%x, dof%x, this%ntot)
          call copy(this%y, dof%y, this%ntot)
          call copy(this%z, dof%z, this%ntot)
       end if

    end if

  end subroutine dofmap_init_and_map

end submodule dofmap_init_m
