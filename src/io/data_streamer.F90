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
!> Implements type data_streamer_t.
module data_streamer
  use num_types, only: rp, c_rp
  use field, only: field_t
  use coefs, only: coef_t
  use utils, only: neko_warning
  use device
  use comm
  use neko_mpi_types
  use neko_config
  use, intrinsic :: iso_c_binding
  implicit none
  private

  !> Provides access to data streaming by interfacing with c++
  !! ADIOS2 subroutines.
  !! @details
  !! Adios2 is an API that allows for easy coupling of codes
  !! through data streaming and gives the posibility to perform
  !! other IO operations such as data compression, etc.
  !! This type wraps and interfaces the needed calls to allow
  !! the use of the c++ routines that ultimately expose the data
  !! from neko to any executable that counts with a proper reader.
  type, public :: data_streamer_t
     !> Define if the execution is asyncrhonous
     integer :: if_asynch
     !> global element numbers
     integer, allocatable :: lglel(:)
   contains
     !> Constructor
     procedure, pass(this) :: init => data_streamer_init
     !> Destructor
     procedure, pass(this) :: free => data_streamer_free
     !> Stream data
     procedure, pass(this) :: stream => data_streamer_stream
     !> Stream back the data
     procedure, pass(this) :: recieve => data_streamer_recieve

  end type data_streamer_t

contains

  !> Constructor
  !! Wraps the adios2 set-up.
  !! @param coef Type that contains geometrical information
  !! on the case.
  !! @param if_asynch Controls whether the asyncrhonous executions
  !! is to be enabled.
  subroutine data_streamer_init(this, coef)
    class(data_streamer_t), intent(inout) :: this
    type(coef_t), intent(inout) :: coef
    integer :: nelb, nelv, nelgv, npts, gdim

    !Assign the set up parameters
    nelv  = coef%msh%nelv
    npts  = coef%Xh%lx*coef%Xh%ly*coef%Xh%lz
    nelgv = coef%msh%glb_nelv
    nelb  = coef%msh%offset_el
    gdim = coef%msh%gdim

#ifdef HAVE_ADIOS2
    call fortran_adios2_initialize(npts, nelv, nelb, nelgv, gdim, NEKO_COMM)
#else
    call neko_warning('Is not being built with ADIOS2 support.')
    call neko_warning('Not able to use stream/compression functionality')
#endif


  end subroutine data_streamer_init

  !> Destructor
  !! wraps the adios2 finalize routine. Closes insitu writer
  subroutine data_streamer_free(this)
    class(data_streamer_t), intent(inout) :: this

#ifdef HAVE_ADIOS2
    call fortran_adios2_finalize()
#else
    call neko_warning('Is not being built with ADIOS2 support.')
    call neko_warning('Not able to use stream/compression functionality')
#endif

  end subroutine data_streamer_free
   
  !> streamer
  !! @param fld array of shape field%x
  subroutine data_streamer_stream(this, fld)
    class(data_streamer_t), intent(inout) :: this
    real(kind=rp), dimension(:,:,:,:), intent(inout) :: fld

#ifdef HAVE_ADIOS2
    call fortran_adios2_stream(fld)
#else
    call neko_warning('Is not being built with ADIOS2 support.')
    call neko_warning('Not able to use stream/compression functionality')
#endif

  end subroutine data_streamer_stream
  
  !> reciever
  !! @param fld array of shape field%x
  subroutine data_streamer_recieve(this, fld)
    class(data_streamer_t), intent(inout) :: this
    real(kind=rp), dimension(:,:,:,:), intent(inout) :: fld

#ifdef HAVE_ADIOS2
    call fortran_adios2_recieve(fld)
#else
    call neko_warning('Is not being built with ADIOS2 support.')
    call neko_warning('Not able to use stream/compression functionality')
#endif

  end subroutine data_streamer_recieve


#ifdef HAVE_ADIOS2

  !> Interface to adios2_initialize in c++.
  !! @details This routine interfaces with c++ routine that set up adios2
  !! if streaming, the global array to pair writer and reader is opened.
  !! @param npts number of points per element
  !! @param nelv number of elements in this rank
  !! @param nelb number of elements in ranks before this one
  !! @param nelgv total number of elements in velocity mesh
  !! @param gdim dimension (2d or 3d)
  !! @param comm simulation communicator
  subroutine fortran_adios2_initialize(npts, nelv, nelb, nelgv, gdim, comm)
    use, intrinsic :: ISO_C_BINDING
    implicit none
    integer, intent(in) :: npts, nelv, nelb, nelgv, gdim
    type(MPI_COMM) :: comm

    interface
       !> C-definition is: void adios2_initialize_(const int *nval,
       !! const int *nelvin,const int *nelb, const int *nelgv,
       !! const int *nelgt, const double *xml,const double *yml,
       !! const double *zml, const int *if_asynchronous,
       !! const int *comm_int)
       subroutine c_adios2_initialize(npts, nelv, nelb, nelgv, gdim, &
                                      comm) bind(C,name="adios2_initialize_")
         use, intrinsic :: ISO_C_BINDING
         import c_rp
         implicit none
         integer(kind=C_INT) :: npts
         integer(kind=C_INT) :: nelv
         integer(kind=C_INT) :: nelb
         integer(kind=C_INT) :: nelgv
         integer(kind=C_INT) :: gdim
         type(*) :: comm
       end subroutine c_adios2_initialize
    end interface

    call c_adios2_initialize(npts, nelv, nelb, nelgv, gdim, comm)
  end subroutine fortran_adios2_initialize

  !> Interface to adios2_finalize in c++.
  !! closes any writer openned at initialization time
  subroutine fortran_adios2_finalize()
    use, intrinsic :: ISO_C_BINDING
    implicit none

    interface
       !> C-definition is: void adios2_finalize_()
       subroutine c_adios2_finalize() bind(C,name="adios2_finalize_")
         use, intrinsic :: ISO_C_BINDING
         implicit none
       end subroutine c_adios2_finalize
    end interface

    call c_adios2_finalize()
  end subroutine fortran_adios2_finalize
  
  !> Interface to adios2_stream in c++.
  !! @details This routine communicates the data to a global array that
  !! is accessed by a data processor. The operations do not write to disk.
  !! data is communicated with mpi.
  !! @param fld array of shape field%x
  subroutine fortran_adios2_stream(fld)
    use, intrinsic :: ISO_C_BINDING
    implicit none
    real(kind=rp), dimension(:,:,:,:), intent(inout) :: fld

    interface
       !> C-definition is: void adios2_stream_(const double *fld)
       subroutine c_adios2_stream(fld) &
                                  bind(C,name="adios2_stream_")
         use, intrinsic :: ISO_C_BINDING
         import c_rp
         implicit none
         real(kind=c_rp), intent(INOUT) :: fld(*)
       end subroutine c_adios2_stream
    end interface

    call c_adios2_stream(fld)
  end subroutine fortran_adios2_stream
  
  !> Interface to adios2_recieve in ci++.
  !! @details This routine communicates the data to a global array that
  !! is accessed by a data processor. The operations do not write to disk.
  !! data is communicated with mpi.
  !! @param fld array of shape field%x
  subroutine fortran_adios2_recieve(fld)
    use, intrinsic :: ISO_C_BINDING
    implicit none
    real(kind=rp), dimension(:,:,:,:), intent(inout) :: fld

    interface
       !> C-definition is: void adios2_stream_(const double *fld)
       subroutine c_adios2_recieve(fld) &
                                  bind(C,name="adios2_recieve_")
         use, intrinsic :: ISO_C_BINDING
         import c_rp
         implicit none
         real(kind=c_rp), intent(INOUT) :: fld(*)
       end subroutine c_adios2_recieve
    end interface

    call c_adios2_recieve(fld)
  end subroutine fortran_adios2_recieve

#endif

end module data_streamer
