! Copyright (c) 2019-2023, The Neko Authors
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
!> Implements a dynamic stack ADT
!! @details a stack storing values @a data of an arbitrary type
module stack
  use num_types
  use nmsh
  use utils, only : neko_error, neko_warning
  use point, only : point_t
  use structs, only : struct_curve_t
  use math, only : NEKO_M_LN2
  use tuple, only : tuple_t, tuple_i4_t, tuple4_i4_t, tuple_i4r8_t, tuple_2i4r8_t
  implicit none
  private
  
  integer, parameter :: NEKO_STACK_SIZE_T = 32

  !> Base type for a stack
  type, abstract, private :: stack_t
     class(*),  allocatable :: data(:)
     integer :: top_
     integer :: size_
   contains
     procedure, non_overridable, pass(this) :: init => stack_init
     procedure, non_overridable, pass(this) :: free => stack_free
     procedure, non_overridable, pass(this) :: clear => stack_clear
     procedure, non_overridable, pass(this) :: size => stack_size
     procedure, non_overridable, pass(this) :: push => stack_push
  end type stack_t

  !> Integer based stack
  type, public, extends(stack_t) :: stack_i4_t
   contains
     procedure, public, pass(this) :: pop => stack_i4_pop
     procedure, public, pass(this) :: array => stack_i4_data
  end type stack_i4_t

  !> Integer*8 based stack
  type, public, extends(stack_t) :: stack_i8_t
   contains
     procedure, public, pass(this) :: pop => stack_i8_pop
     procedure, public, pass(this) :: array => stack_i8_data
  end type stack_i8_t

  !> Double precision based stack
  type, public, extends(stack_t) :: stack_r8_t
   contains
     procedure, public, pass(this) :: pop => stack_r8_pop
     procedure, public, pass(this) :: array => stack_r8_data
  end type stack_r8_t

  !> Integer 2-tuple based stack
  type, public, extends(stack_t) :: stack_i4t2_t
   contains
     procedure, public, pass(this) :: pop => stack_i4t2_pop
     procedure, public, pass(this) :: array => stack_i4t2_data
  end type stack_i4t2_t

  !> Integer 4-tuple based stack
  type, public, extends(stack_t) :: stack_i4t4_t
   contains
     procedure, public, pass(this) :: pop => stack_i4t4_pop
     procedure, public, pass(this) :: array => stack_i4t4_data
  end type stack_i4t4_t

  !> Mixed integer-double precision 2-tuple based stack
  type, public, extends(stack_t) :: stack_i4r8t2_t
   contains
     procedure, public, pass(this) :: pop => stack_i4r8t2_pop
     procedure, public, pass(this) :: array => stack_i4r8t2_data
  end type stack_i4r8t2_t

  !> Mixed integer-double precision 3-tuple based stack
  type, public, extends(stack_t) :: stack_2i4r8t3_t
   contains
     procedure, public, pass(this) :: pop => stack_2i4r8t3_pop
     procedure, public, pass(this) :: array => stack_2i4r8t3_data
  end type stack_2i4r8t3_t
  
  !> Curved element stack
  type, public, extends(stack_t) :: stack_curve_t
   contains
     procedure, public, pass(this) :: pop => stack_curve_element_pop
     procedure, public, pass(this) :: array => stack_curve_element_data
  end type stack_curve_t

  !> Neko quad element based stack
  type, public, extends(stack_t) :: stack_nq_t
   contains
     procedure, public, pass(this) :: pop => stack_nq_pop
     procedure, public, pass(this) :: array => stack_nq_data
  end type stack_nq_t

  !> Neko hex element based stack
  type, public, extends(stack_t) :: stack_nh_t
   contains
     procedure, public, pass(this) :: pop => stack_nh_pop
     procedure, public, pass(this) :: array => stack_nh_data
  end type stack_nh_t

  !> Neko zone based stack
  type, public, extends(stack_t) :: stack_nz_t
   contains
     procedure, public, pass(this) :: pop => stack_nz_pop
     procedure, public, pass(this) :: array => stack_nz_data
  end type stack_nz_t

  !> Neko curve info based stack
  type, public, extends(stack_t) :: stack_nc_t
   contains
     procedure, public, pass(this) :: pop => stack_nc_pop
     procedure, public, pass(this) :: array => stack_nc_data
  end type stack_nc_t

  !> Point based stack
  type, public, extends(stack_t) :: stack_pt_t
   contains
     procedure, public, pass(this) :: pop => stack_pt_pop
     procedure, public, pass(this) :: array => stack_pt_data
  end type stack_pt_t

contains

  !> Initialize a stack of arbitrary type 
  subroutine stack_init(this, size)
    class(stack_t), intent(inout) :: this 
    integer, optional :: size !< Initial size of the stack
    integer :: size_t

    if (present(size)) then
       if (size .gt. 0) then
          size_t = size
       else
          call neko_warning('Invalid stack size, using default')
          size_t = NEKO_STACK_SIZE_T
       end if
    else
       size_t = NEKO_STACK_SIZE_T
    end if

    this%size_ = ishft(1, ceiling(log(real(size_t, rp)) / NEKO_M_LN2))
    this%top_ = 0
    select type(this)
    type is(stack_i4_t)
       allocate(integer::this%data(this%size_))
    type is(stack_i8_t)
       allocate(integer(i8)::this%data(this%size_))
    type is (stack_r8_t)
       allocate(double precision::this%data(this%size_))
    type is (stack_i4t2_t)
       allocate(tuple_i4_t::this%data(this%size_))
    type is (stack_i4t4_t)
       allocate(tuple4_i4_t::this%data(this%size_))
    type is (stack_i4r8t2_t)
       allocate(tuple_i4r8_t::this%data(this%size_))
    type is (stack_2i4r8t3_t)
       allocate(tuple_2i4r8_t::this%data(this%size_))
    type is (stack_curve_t)
       allocate(struct_curve_t::this%data(this%size_))
    type is (stack_nq_t)
       allocate(nmsh_quad_t::this%data(this%size_))
    type is (stack_nh_t)
       allocate(nmsh_hex_t::this%data(this%size_))
    type is (stack_nz_t)
       allocate(nmsh_zone_t::this%data(this%size_))
    type is (stack_nc_t)
       allocate(nmsh_curve_el_t::this%data(this%size_))
    type is (stack_pt_t)
       allocate(point_t::this%data(this%size_))
    class default
       call neko_error('Invalid data type')
    end select

  end subroutine stack_init
  
  !> Destroy a stack
  subroutine stack_free(this)
    class(stack_t), intent(inout) :: this
    
    if (allocated(this%data)) then
       deallocate(this%data)
       this%size_ = 0 
       this%top_ = 0
    end if    

  end subroutine stack_free

  !> Clear all entries of a stack
  subroutine stack_clear(this)
    class(stack_t), intent(inout) :: this
    this%top_ = 0
  end subroutine stack_clear

  !> Return number of entries in the stack
  pure function stack_size(this) result(size)
    class(stack_t), intent(in) :: this
    integer :: size
    size = this%top_
  end function stack_size

  !> Push data onto the stack
  subroutine stack_push(this, data)
    class(stack_t), target, intent(inout) :: this
    class(*), intent(inout) :: data !< Arbitrary typed data (same type as stack)
    class(*), allocatable :: tmp(:)
    integer :: i

    if (this%top_ .eq. this%size_) then
       this%size_ = ishft(this%size_, 1)
       select type(data)
       type is(integer)
          allocate(integer::tmp(this%size_))
       type is(integer(i8))
          allocate(integer(i8)::tmp(this%size_))
       type is(double precision)          
          allocate(double precision::tmp(this%size_))
       type is(tuple_i4_t)
          allocate(tuple_i4_t::tmp(this%size_))
       type is(tuple4_i4_t)
          allocate(tuple4_i4_t::tmp(this%size_))
       type is(tuple_i4r8_t)
          allocate(tuple_i4r8_t::tmp(this%size_))
       type is(tuple_2i4r8_t)
          allocate(tuple_2i4r8_t::tmp(this%size_))
       type is(struct_curve_t)
          allocate(struct_curve_t::tmp(this%size_))
       type is (nmsh_quad_t)
          allocate(nmsh_quad_t::tmp(this%size_))
       type is (nmsh_hex_t)
          allocate(nmsh_hex_t::tmp(this%size_))
       type is (nmsh_zone_t)
          allocate(nmsh_zone_t::tmp(this%size_))
       type is (nmsh_curve_el_t)
          allocate(nmsh_curve_el_t::tmp(this%size_))
       type is (point_t)
          allocate(point_t::tmp(this%size_))
       class default
          call neko_error('Invalid data type (stack_push)')
       end select
       
       select type(tmp)
       type is (integer)
          select type(sdp=>this%data)
          type is (integer)
             tmp(1:this%top_) = sdp
          end select
       type is (integer(i8))
          select type(sdp=>this%data)
          type is (integer(i8))
             tmp(1:this%top_) = sdp
          end select
       type is (double precision)
          select type(sdp=>this%data)
          type is (double precision)
             tmp(1:this%top_) = sdp
          end select
       type is (tuple_i4_t)
          select type(sdp=>this%data)
          type is (tuple_i4_t)
             do i = 1, this%top_
                tmp(i) = sdp(i)
             end do
          end select
       type is (tuple4_i4_t)
          select type(sdp=>this%data)
          type is (tuple4_i4_t)
             do i = 1, this%top_
                tmp(i) = sdp(i)
             end do
          end select
       type is (tuple_i4r8_t)
          select type(sdp=>this%data)
          type is (tuple_i4r8_t)
             do i = 1, this%top_
                tmp(i) = sdp(i)
             end do
          end select
       type is (tuple_2i4r8_t)
          select type(sdp=>this%data)
          type is (tuple_2i4r8_t)
             do i = 1, this%top_
                tmp(i) = sdp(i)
             end do
          end select
       type is (struct_curve_t)
          select type(sdp=>this%data)
          type is (struct_curve_t)
             tmp(1:this%top_) = sdp
          end select
       type is (nmsh_quad_t)
          select type(sdp=>this%data)
          type is(nmsh_quad_t)
             tmp(1:this%top_) = sdp
          end select
       type is (nmsh_hex_t)
          select type(sdp=>this%data)
          type is(nmsh_hex_t)
             tmp(1:this%top_) = sdp
          end select
       type is (nmsh_zone_t)
          select type(sdp=>this%data)
          type is(nmsh_zone_t)
             tmp(1:this%top_) = sdp
          end select
       type is (nmsh_curve_el_t)
          select type(sdp=>this%data)
          type is(nmsh_curve_el_t)
             tmp(1:this%top_) = sdp
          end select
       type is (point_t)
          select type(sdp=>this%data)
          type is(point_t)
             tmp(1:this%top_) = sdp
          end select
       class default
          call neko_error('Invalid data type (stack_push tmp)')
       end select
       call move_alloc(tmp, this%data)
    end if
    
    this%top_ = this%top_ + 1

    select type(sdp=>this%data)
    type is (integer)
       select type(data)
       type is (integer)
          sdp(this%top_) = data
       end select
    type is (integer(i8))
       select type(data)
       type is (integer(i8))
          sdp(this%top_) = data
       end select
    type is (double precision)
       select type(data)
       type is (double precision)
          sdp(this%top_) = data
       end select
    type is (tuple_i4_t)
       select type(data)
       type is (tuple_i4_t)
          sdp(this%top_) = data
       end select
    type is (tuple4_i4_t)
       select type(data)
       type is (tuple4_i4_t)
          sdp(this%top_) = data
       end select
    type is (tuple_i4r8_t)
       select type(data)
       type is (tuple_i4r8_t)
          sdp(this%top_) = data
       end select
    type is (tuple_2i4r8_t)
       select type(data)
       type is (tuple_2i4r8_t)
          sdp(this%top_) = data
       end select
    type is (struct_curve_t)
       select type(data)
       type is (struct_curve_t)
          sdp(this%top_) = data
       end select
    type is (nmsh_quad_t)
       select type(data)
       type is (nmsh_quad_t)
          sdp(this%top_) = data
       end select
    type is (nmsh_hex_t)
       select type(data)
       type is (nmsh_hex_t)
          sdp(this%top_) = data
       end select
    type is (nmsh_zone_t)
       select type(data)
       type is (nmsh_zone_t)
          sdp(this%top_) = data
       end select
    type is (nmsh_curve_el_t)
       select type(data)
       type is (nmsh_curve_el_t)
          sdp(this%top_) = data
       end select
    type is (point_t)
       select type(data)
       type is (point_t)
          sdp(this%top_) = data
       end select
    class default
       call neko_error('Invalid data type in stack (stack_push)')
    end select
  end subroutine stack_push

  !> Pop an integer of the stack
  function stack_i4_pop(this) result(data)
    class(stack_i4_t), target, intent(inout) :: this
    integer :: data

    select type (sdp=>this%data)
    type is (integer)       
       data = sdp(this%top_)
    class default
       call neko_error('Invalid data type (i4 pop)')
    end select
    this%top_ = this%top_ - 1
  end function stack_i4_pop

  !> Return a pointer to the internal integer array
  function stack_i4_data(this) result(data)
    class(stack_i4_t), target, intent(inout) :: this
    integer, pointer :: data(:)

    select type (sdp=>this%data)
    type is (integer)       
       data => sdp
    class default
       call neko_error('Invalid data type (i4 array)')
    end select
  end function stack_i4_data

  !> Pop an integer*8 of the stack
  function stack_i8_pop(this) result(data)
    class(stack_i8_t), target, intent(inout) :: this
    integer(kind=i8) :: data

    select type (sdp=>this%data)
    type is (integer(i8))       
       data = sdp(this%top_)
    class default
       call neko_error('Invalid data type (i8 pop)')
    end select
    this%top_ = this%top_ - 1
  end function stack_i8_pop

  !> Return a pointer to the internal integer*8 array
  function stack_i8_data(this) result(data)
    class(stack_i8_t), target, intent(inout) :: this
    integer(kind=i8), pointer :: data(:)

    select type (sdp=>this%data)
    type is (integer(i8))       
       data => sdp
    class default
       call neko_error('Invalid data type (i8 array)')
    end select
  end function stack_i8_data

  !> Pop a double precision value of the stack
  function stack_r8_pop(this) result(data)
    class(stack_r8_t), target, intent(inout) :: this
    real(kind=dp) :: data
    
    select type (sdp=>this%data)
    type is (double precision)       
       data = sdp(this%top_)
    class default
       call neko_error('Invalid data type (r8 pop)')
    end select
    this%top_ = this%top_ -1
  end function stack_r8_pop

  !> Return a pointer to the internal double precision array 
  function stack_r8_data(this) result(data)
    class(stack_r8_t), target, intent(inout) :: this
    real(kind=dp), pointer :: data(:)

    select type (sdp=>this%data)
    type is (double precision)       
       data => sdp
    class default
       call neko_error('Invalid data type (r8 array)')
    end select
  end function stack_r8_data

  !> Pop an integer 2-tuple of the stack
  function stack_i4t2_pop(this) result(data)
    class(stack_i4t2_t), target, intent(inout) :: this
    type(tuple_i4_t) :: data
    
    select type (sdp=>this%data)
    type is (tuple_i4_t)       
       data = sdp(this%top_)
    class default
       call neko_error('Invalid data type (i4t2 pop)')
    end select
    this%top_ = this%top_ -1
  end function stack_i4t2_pop

  !> Return a pointer to the interal 2-tuple array
  function stack_i4t2_data(this) result(data)
    class(stack_i4t2_t), target, intent(inout) :: this
    type(tuple_i4_t), pointer :: data(:)

    select type (sdp=>this%data)
    type is (tuple_i4_t)       
       data => sdp
    class default
       call neko_error('Invalid data type (i4t2 array)')
    end select
  end function stack_i4t2_data

  !> Pop an integer 4-tuple of the stack
  function stack_i4t4_pop(this) result(data)
    class(stack_i4t4_t), target, intent(inout) :: this
    type(tuple4_i4_t) :: data
    
    select type (sdp=>this%data)
    type is (tuple4_i4_t)       
       data = sdp(this%top_)
    class default
       call neko_error('Invalid data type (i4t4 pop)')
    end select
    this%top_ = this%top_ -1
  end function stack_i4t4_pop

  !> Return a pointer to the internal 4-tuple array
  function stack_i4t4_data(this) result(data)
    class(stack_i4t4_t), target, intent(inout) :: this
    type(tuple4_i4_t), pointer :: data(:)

    select type (sdp=>this%data)
    type is (tuple4_i4_t)       
       data => sdp
    class default
       call neko_error('Invalid data type (i4t4 array)')
    end select
  end function stack_i4t4_data

  !> Pop a mixed integer-double precision  2-tuple of the stack
  function stack_i4r8t2_pop(this) result(data)
    class(stack_i4r8t2_t), target, intent(inout) :: this
    type(tuple_i4r8_t) :: data
    
    select type (sdp=>this%data)
    type is (tuple_i4r8_t)       
       data = sdp(this%top_)
    class default
       call neko_error('Invalid data type (i4r8t2 pop)')
    end select
    this%top_ = this%top_ -1
  end function stack_i4r8t2_pop

  !> Return a pointer to the internal 2-tuple array
  function stack_i4r8t2_data(this) result(data)
    class(stack_i4r8t2_t), target, intent(inout) :: this
    type(tuple_i4r8_t), pointer :: data(:)

    select type (sdp=>this%data)
    type is (tuple_i4r8_t)       
       data => sdp
    class default
       call neko_error('Invalid data type (i4r8t2 array)')
    end select
  end function stack_i4r8t2_data

  !> Pop a mixed integer-double precision  3-tuple of the stack
  function stack_2i4r8t3_pop(this) result(data)
    class(stack_2i4r8t3_t), target, intent(inout) :: this
    type(tuple_2i4r8_t) :: data
    
    select type (sdp=>this%data)
    type is (tuple_2i4r8_t)       
       data = sdp(this%top_)
    class default
       call neko_error('Invalid data type (i4r8t2 pop)')
    end select
    this%top_ = this%top_ -1
  end function stack_2i4r8t3_pop

  !> Return a pointer to the internal 2-tuple array
  function stack_2i4r8t3_data(this) result(data)
    class(stack_2i4r8t3_t), target, intent(inout) :: this
    type(tuple_2i4r8_t), pointer :: data(:)

    select type (sdp=>this%data)
    type is (tuple_2i4r8_t)       
       data => sdp
    class default
       call neko_error('Invalid data type (i4r8t2 array)')
    end select
  end function stack_2i4r8t3_data
 
  !> Pop a curve element of the stack
  function stack_curve_element_pop(this) result(data)
    class(stack_curve_t), target, intent(inout) :: this
    type(struct_curve_t) :: data
    
    select type (sdp=>this%data)
    type is (struct_curve_t)       
       data = sdp(this%top_)
    class default
       call neko_error('Invalid data type (curve pop)')
    end select
    this%top_ = this%top_ -1
  end function stack_curve_element_pop

  !> Return a pointer to the internal curve element array
  function stack_curve_element_data(this) result(data)
    class(stack_curve_t), target, intent(inout) :: this
    type(struct_curve_t), pointer :: data(:)

    select type (sdp=>this%data)
    type is (struct_curve_t)       
       data => sdp
    class default
       call neko_error('Invalid data type (curve array)')
    end select
  end function stack_curve_element_data

  !> Pop a Neko quad element of the stack
  function stack_nq_pop(this) result(data)
    class(stack_nq_t), target, intent(inout) :: this
    type(nmsh_quad_t) :: data

    select type (sdp=>this%data)
    type is (nmsh_quad_t)       
       data = sdp(this%top_)
    class default
       call neko_error('Invalid data type (nq pop)')
    end select
    this%top_ = this%top_ -1
  end function stack_nq_pop

  !> Return a pointer to the internal Neko quad array
  function stack_nq_data(this) result(data)
    class(stack_nq_t), target, intent(inout) :: this
    type(nmsh_quad_t), pointer :: data(:)

    select type (sdp=>this%data)
    type is (nmsh_quad_t)       
       data => sdp
    class default
       call neko_error('Invalid data type (nq array)')
    end select
  end function stack_nq_data

  !> Pop a Neko hex element of the stack
  function stack_nh_pop(this) result(data)
    class(stack_nh_t), target, intent(inout) :: this
    type(nmsh_hex_t) :: data

    select type (sdp=>this%data)
    type is (nmsh_hex_t)       
       data = sdp(this%top_)
    class default
       call neko_error('Invalid data type (nh pop)')
    end select
    this%top_ = this%top_ -1
  end function stack_nh_pop

  !> Return a pointer to the internal Neko quad array
  function stack_nh_data(this) result(data)
    class(stack_nh_t), target, intent(inout) :: this
    type(nmsh_hex_t), pointer :: data(:)

    select type (sdp => this%data)
    type is (nmsh_hex_t)       
       data => sdp
    class default
       call neko_error('Invalid data type (nh array)')
    end select
  end function stack_nh_data

  !> Pop a Neko zone of the stack
  function stack_nz_pop(this) result(data)
    class(stack_nz_t), target, intent(inout) :: this
    type(nmsh_zone_t) :: data

    select type (sdp=>this%data)
    type is (nmsh_zone_t)       
       data = sdp(this%top_)
    class default
       call neko_error('Invalid data type (nz pop)')
    end select
    this%top_ = this%top_ -1
  end function stack_nz_pop

  !> Return a pointer to the internal Neko zone array
  function stack_nz_data(this) result(data)
    class(stack_nz_t), target, intent(inout) :: this
    type(nmsh_zone_t), pointer :: data(:)

    select type (sdp=>this%data)
    type is (nmsh_zone_t)       
       data => sdp
    class default
       call neko_error('Invalid data type (nz array)')
    end select
  end function stack_nz_data

  !> Pop a Neko curve info of the stack
  function stack_nc_pop(this) result(data)
    class(stack_nc_t), target, intent(inout) :: this
    type(nmsh_curve_el_t) :: data

    select type (sdp=>this%data)
    type is (nmsh_curve_el_t)       
       data = sdp(this%top_)
    class default
       call neko_error('Invalid data type (nc pop)')
    end select
    this%top_ = this%top_ -1
  end function stack_nc_pop

  !> Return a pointer to the internal Neko curve info array
  function stack_nc_data(this) result(data)
    class(stack_nc_t), target, intent(inout) :: this
    type(nmsh_curve_el_t), pointer :: data(:)

    select type (sdp=>this%data)
    type is (nmsh_curve_el_t)       
       data => sdp
    class default
       call neko_error('Invalid data type (nc array)')
    end select
  end function stack_nc_data

  !> Pop a point of the stack
  function stack_pt_pop(this) result(data)
    class(stack_pt_t), target, intent(inout) :: this
    type(point_t) :: data

    select type (sdp=>this%data)
    type is (point_t)       
       data = sdp(this%top_)
    class default
       call neko_error('Invalid data type (point pop)')
    end select
    this%top_ = this%top_ -1
  end function stack_pt_pop

  !> Return a pointer to the internal point array
  function stack_pt_data(this) result(data)
    class(stack_pt_t), target, intent(inout) :: this
    type(point_t), pointer :: data(:)

    select type (sdp=>this%data)
    type is (point_t)       
       data => sdp
    class default
       call neko_error('Invalid data type (point array)')
    end select
  end function stack_pt_data
  
end module stack
