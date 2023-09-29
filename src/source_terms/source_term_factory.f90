! Copyright (c) 2023, The Neko Authors
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
!
!> Defines a factory subroutine for source terms.
module source_term_fctry
  use source_term, only : source_term_t
  use const_source_term, only : const_source_term_t
  use json_module, only : json_file
  use json_utils, only : json_get
  use field_list, only : field_list_t
  use utils, only : neko_error
  use coefs, only : coef_t
  implicit none
  private
  
  public :: source_term_factory
  
  contains

  !> Source term factory. Both constructs and initializes the object.
  !! @param json JSON object initializing the source term.
  subroutine source_term_factory(source_term, json, fields, coef)
       class(source_term_t), allocatable, intent(inout) :: source_term
       type(json_file), intent(inout) :: json
       type(field_list_t), intent(inout) :: fields
       type(coef_t), intent(inout) :: coef
       character(len=:), allocatable :: source_type
              
       call json_get(json, "type", source_type)

       if (trim(source_type) .eq. "constant") then 
          allocate(const_source_term_t::source_term)
       else
           call neko_error('Unknown source term '//trim(source_type))
       end if
       
       ! Initialize
       call source_term%init(json, fields, coef)

  end subroutine source_term_factory

end module source_term_fctry
