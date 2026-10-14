! Copyright (c) 2026, The Neko Authors
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
! Version, passed in by configure; "unknown" when compiled on its own
#ifndef PACKAGE_VERSION
#define PACKAGE_VERSION "unknown"
#endif

!> Kinds, messages and error handling for gmsh2nmsh
!! @details gmsh2nmsh does not depend on the Neko library, so that it can be
!! compiled on its own with any Fortran 2008 compiler.
module gmsh2nmsh_util
  use, intrinsic :: iso_fortran_env, only : int32, int64, real64, &
       error_unit, output_unit
  implicit none
  private

  integer, public, parameter :: i4 = int32
  integer, public, parameter :: i8 = int64
  integer, public, parameter :: dp = real64

  public :: fatal, warn, progress, int_str, real_str

contains

  !> Print an error message and stop
  !! @param msg Message to print before stopping
  subroutine fatal(msg)
    character(len=*), intent(in) :: msg

    flush(output_unit)
    write(error_unit, '(A)') ' *** ERROR: ' // trim(msg) // ' ***'
    flush(error_unit)
    stop 1

  end subroutine fatal

  !> Print a warning
  !! @param msg Message to print
  subroutine warn(msg)
    character(len=*), intent(in) :: msg

    write(output_unit, '(A)') ' *** WARNING: ' // trim(msg) // ' ***'

  end subroutine warn

  !> Print a progress line: the label, a row of dots and the value
  !! @details The output is flushed, so that progress is visible at once
  !! also when it is redirected.
  !! @param label Text at the start of the line
  !! @param value Text after the dots
  subroutine progress(label, value)
    character(len=*), intent(in) :: label, value
    integer, parameter :: WIDTH = 33
    integer :: ndots

    ndots = max(3, WIDTH - len_trim(label) - 1)
    write(output_unit, '(A)') ' ' // trim(label) // ' ' // &
         repeat('.', ndots) // ' ' // trim(value)
    flush(output_unit)

  end subroutine progress

  !> An integer as text
  !! @param i The integer
  function int_str(i) result(str)
    integer, intent(in) :: i
    character(len=:), allocatable :: str
    character(len=24) :: buf

    write(buf, '(I0)') i
    str = trim(buf)

  end function int_str

  !> A real as short text, without trailing zeros
  !! @param x The real
  function real_str(x) result(str)
    real(kind=dp), intent(in) :: x
    character(len=:), allocatable :: str
    character(len=32) :: buf
    integer :: i

    if (abs(x) .lt. tiny(1.0_dp)) then
       str = '0'
       return
    end if

    if (abs(x) .ge. 1.0e-3_dp .and. abs(x) .lt. 1.0e6_dp) then
       write(buf, '(F24.6)') x
       buf = adjustl(buf)
       i = len_trim(buf)
       do while (buf(i:i) .eq. '0')
          buf(i:i) = ' '
          i = i - 1
       end do
       if (buf(i:i) .eq. '.') buf(i:i) = ' '
    else
       write(buf, '(ES10.3)') x
    end if
    str = trim(adjustl(buf))

  end function real_str

end module gmsh2nmsh_util

!> Reader for meshes in the Gmsh MSH format
!! @details Supports the MSH 2.2 and 4.1 formats, in ASCII and binary form.
!! Only the data needed to build a Neko mesh is kept: the node coordinates,
!! first and second order hexahedra (3D) or quadrilaterals (2D), and the
!! boundary facets (quadrilaterals or lines) of physical groups.
module gmsh2nmsh_reader
  use gmsh2nmsh_util, only : i4, i8, dp, fatal, progress, int_str
  implicit none
  private

  !> Maximum length of a physical group name
  integer, public, parameter :: GMSH_NAME_LEN = 128

  !> Slot of each of the first 20 nodes of a Gmsh hexahedron
  !! @details Slots 1-8 hold the vertices in Neko's symmetric order, slots
  !! 9-20 the midpoints of the 12 edges in the curve edge order used by
  !! `dofmap_xyzquad` (the Nek5000 convention). The face and volume centre
  !! nodes of 27-node hexahedra are not used.
  integer, public, parameter :: HEX_G2S(20) = [1, 2, 4, 3, 5, 6, 8, 7, &
       9, 12, 17, 10, 18, 11, 19, 20, 13, 16, 14, 15]

  !> Slot of each of the first 8 nodes of a Gmsh quadrilateral
  !! @details Slots 1-4 hold the vertices in symmetric order, slots 5-8 the
  !! edge midpoints in curve edge order (bottom, right, top, left).
  integer, public, parameter :: QUAD_G2S(8) = [1, 2, 4, 3, 5, 6, 7, 8]

  ! Element categories
  integer, parameter :: KIND_IGNORE = 0
  integer, parameter :: KIND_LINE = 1
  integer, parameter :: KIND_QUAD = 2
  integer, parameter :: KIND_HEX = 3

  !> Number of nodes of the Gmsh element types 1 to 19
  integer, parameter :: MSH_NNODES(19) = [2, 3, 4, 4, 8, 6, 5, 3, 6, 9, &
       10, 27, 18, 14, 1, 8, 20, 15, 13]

  !> Size of the buffer used when parsing ASCII data
  integer, parameter :: STREAM_BUF_LEN = 4194304

  !> Number of nodes or elements decoded per binary read
  integer, parameter :: BIN_CHUNK = 262144

  !> Buffered reader for files mixing ASCII and binary data
  type :: msh_stream_t
     integer :: unit = -1 !< Fortran unit of the file, -1 when closed
     integer(kind=i8) :: fsize = 0 !< File size in bytes
     integer(kind=i8) :: pos = 1 !< Next byte to consume
     integer(kind=i8) :: buf_start = 1 !< File position of buf(1:1)
     integer :: buf_len = 0 !< Number of valid bytes in buf
     character(len=:), allocatable :: buf !< Read buffer
     logical :: binary = .false. !< Read binary data
     logical :: swap = .false. !< Swap the byte order of data
     integer :: size_t_len = 8 !< Byte length of size_t values
   contains
     procedure, pass(this) :: free => msh_stream_free
  end type msh_stream_t

  !> First physical tag of the entities of one dimension (MSH 4)
  type :: entity_list_t
     integer :: n = 0 !< Number of entities
     integer, allocatable :: tag(:) !< Entity tags
     !> First physical tag of each entity, 0 if it is in no group
     integer, allocatable :: phys(:)
   contains
     procedure, pass(this) :: free => entity_list_free
  end type entity_list_t

  !> Map from Gmsh node tags to node indices
  type :: node_map_t
     integer(kind=i8) :: tag_min = 1 !< Smallest node tag
     !> Index of the node with tag tag_min + i - 1, 0 if there is none
     integer, allocatable :: idx(:)
   contains
     procedure, pass(this) :: free => node_map_free
  end type node_map_t

  !> Mesh data extracted from a Gmsh file
  type, public :: gmsh_mesh_t
     real(kind=dp) :: version = 0.0_dp !< MSH format version
     logical :: binary = .false. !< The file is binary
     integer :: nnodes = 0 !< Number of nodes
     !> Node coordinates (3, nnodes)
     real(kind=dp), allocatable :: xyz(:,:)
     integer :: nphys = 0 !< Number of physical names
     integer, allocatable :: phys_dim(:) !< Dimension of each physical name
     integer, allocatable :: phys_tag(:) !< Tag of each physical name
     !> Physical names, without quotes
     character(len=GMSH_NAME_LEN), allocatable :: phys_name(:)
     integer :: nhex = 0 !< Number of hexahedra
     !> Number of hexahedra with edge nodes
     integer :: nhex_ho = 0
     !> Node indices of the hexahedra (20, nhex), see HEX_G2S
     integer, allocatable :: hex(:,:)
     integer :: nquad = 0 !< Number of quadrilaterals
     !> Number of quadrilaterals with edge nodes
     integer :: nquad_ho = 0
     !> Node indices of the quadrilaterals (8, nquad), see QUAD_G2S
     integer, allocatable :: quad(:,:)
     !> Physical tag of each quadrilateral, 0 if none
     integer, allocatable :: quad_phys(:)
     !> Number of lines in a physical group
     integer :: nline = 0
     !> Node indices of the lines (3, nline): end points and midpoint
     integer, allocatable :: line(:,:)
     integer, allocatable :: line_phys(:) !< Physical tag of each line
     !> Number of curves and surfaces in more than one physical group
     integer :: n_multi_phys = 0
     !> Number of volumes in the geometry (MSH 4 only, else 0)
     integer :: n_volumes = 0
   contains
     procedure, pass(this) :: free => gmsh_mesh_free
  end type gmsh_mesh_t

  public :: gmsh_read

contains

  !> Close the file of a stream and release its buffer
  subroutine msh_stream_free(this)
    class(msh_stream_t), intent(inout) :: this
    logical :: is_open

    if (this%unit .ne. -1) then
       inquire(unit = this%unit, opened = is_open)
       if (is_open) close(this%unit)
    end if
    this%unit = -1
    if (allocated(this%buf)) deallocate(this%buf)
    this%fsize = 0
    this%pos = 1
    this%buf_start = 1
    this%buf_len = 0
    this%binary = .false.
    this%swap = .false.
    this%size_t_len = 8

  end subroutine msh_stream_free

  !> Deallocate an entity list
  subroutine entity_list_free(this)
    class(entity_list_t), intent(inout) :: this

    if (allocated(this%tag)) deallocate(this%tag)
    if (allocated(this%phys)) deallocate(this%phys)
    this%n = 0

  end subroutine entity_list_free

  !> Deallocate a node map
  subroutine node_map_free(this)
    class(node_map_t), intent(inout) :: this

    if (allocated(this%idx)) deallocate(this%idx)
    this%tag_min = 1

  end subroutine node_map_free

  !> Deallocate all data of a Gmsh mesh
  subroutine gmsh_mesh_free(this)
    class(gmsh_mesh_t), intent(inout) :: this

    if (allocated(this%xyz)) deallocate(this%xyz)
    if (allocated(this%phys_dim)) deallocate(this%phys_dim)
    if (allocated(this%phys_tag)) deallocate(this%phys_tag)
    if (allocated(this%phys_name)) deallocate(this%phys_name)
    if (allocated(this%hex)) deallocate(this%hex)
    if (allocated(this%quad)) deallocate(this%quad)
    if (allocated(this%quad_phys)) deallocate(this%quad_phys)
    if (allocated(this%line)) deallocate(this%line)
    if (allocated(this%line_phys)) deallocate(this%line_phys)
    this%version = 0.0_dp
    this%binary = .false.
    this%nnodes = 0
    this%nphys = 0
    this%nhex = 0
    this%nhex_ho = 0
    this%nquad = 0
    this%nquad_ho = 0
    this%nline = 0
    this%n_multi_phys = 0
    this%n_volumes = 0

  end subroutine gmsh_mesh_free

  !> Read a Gmsh mesh file
  !! @param fname Name of the .msh file
  !! @param gm Extracted mesh data
  subroutine gmsh_read(fname, gm)
    character(len=*), intent(in) :: fname
    type(gmsh_mesh_t), intent(inout) :: gm
    type(msh_stream_t) :: s
    type(entity_list_t) :: ent(0:3)
    type(node_map_t) :: nmap
    character(len=256) :: name
    character(len=16) :: str
    logical :: ok, have_format, have_nodes
    integer :: d

    call gm%free()
    call stream_open(s, fname)

    have_format = .false.
    have_nodes = .false.
    do
       call stream_section(s, name, ok)
       if (.not. ok) exit

       if (.not. have_format .and. trim(name) .ne. 'MeshFormat') then
          call fatal('Not a Gmsh mesh file: ' // trim(fname))
       end if

       select case (trim(name))
       case ('MeshFormat')
          call read_format(s, gm)
          have_format = .true.
          write(str, '(F4.1)') gm%version
          if (gm%binary) then
             call progress('Reading ' // fname, 'MSH ' // &
                  trim(adjustl(str)) // ', binary')
          else
             call progress('Reading ' // fname, 'MSH ' // &
                  trim(adjustl(str)) // ', ASCII')
          end if
       case ('PhysicalNames')
          call read_physical_names(s, gm)
       case ('Entities')
          call read_entities(s, gm, ent)
       case ('PartitionedEntities')
          call fatal('Partitioned Gmsh meshes are not supported, ' // &
               'save the mesh without partitions')
       case ('Nodes')
          if (gm%version .lt. 3.0_dp) then
             call read_nodes_v2(s, gm, nmap)
          else
             call read_nodes_v4(s, gm, nmap)
          end if
          have_nodes = .true.
          call progress('  Nodes', int_str(gm%nnodes))
       case ('Elements')
          if (.not. have_nodes) then
             call fatal('Gmsh file has $Elements before $Nodes')
          end if
          if (gm%version .lt. 3.0_dp) then
             call read_elements_v2(s, gm, nmap)
          else
             call read_elements_v4(s, gm, nmap, ent)
          end if
       case default
          call stream_skip_section(s, name)
          cycle
       end select

       call stream_end_section(s, name)
    end do

    call s%free()
    do d = 0, 3
       call ent(d)%free()
    end do
    call nmap%free()

    if (.not. have_nodes) then
       call fatal('No $Nodes section found in ' // trim(fname))
    end if

  end subroutine gmsh_read

  !> Read the $MeshFormat section
  !! @param s Input stream
  !! @param gm Mesh receiving the version and file type
  subroutine read_format(s, gm)
    type(msh_stream_t), intent(inout) :: s
    type(gmsh_mesh_t), intent(inout) :: gm
    integer(kind=i4) :: one(1)
    integer :: file_type, data_size
    character(len=32) :: str

    gm%version = stream_real(s)
    file_type = int(stream_int(s))
    data_size = int(stream_int(s))

    if (gm%version .lt. 2.0_dp .or. gm%version .ge. 5.0_dp .or. &
         (gm%version .ge. 3.0_dp .and. gm%version .lt. 4.05_dp)) then
       write(str, '(F0.1)') gm%version
       call fatal('Unsupported Gmsh format version ' // &
            trim(str) // ', save the mesh in format 4.1 or 2.2')
    end if

    gm%binary = file_type .eq. 1
    s%binary = gm%binary
    if (gm%version .ge. 4.0_dp) then
       if (data_size .ne. 4 .and. data_size .ne. 8) then
          call fatal('Unsupported size_t length in Gmsh file')
       end if
       s%size_t_len = data_size
    end if

    if (gm%binary) then
       call stream_skip_line(s)
       call stream_read_i4(s, one)
       if (one(1) .ne. 1) then
          if (swap_i4(one(1)) .eq. 1) then
             s%swap = .true.
          else
             call fatal('Cannot determine the byte order of the ' // &
                  'binary Gmsh file')
          end if
       end if
    end if

  end subroutine read_format

  !> Read the $PhysicalNames section (always ASCII)
  !! @param s Input stream
  !! @param gm Mesh receiving the physical names
  subroutine read_physical_names(s, gm)
    type(msh_stream_t), intent(inout) :: s
    type(gmsh_mesh_t), intent(inout) :: gm
    character(len=1024) :: str
    integer :: i, q
    logical :: ok

    if (allocated(gm%phys_dim)) then
       deallocate(gm%phys_dim, gm%phys_tag, gm%phys_name)
    end if
    gm%nphys = int(stream_int(s))
    allocate(gm%phys_dim(gm%nphys), gm%phys_tag(gm%nphys))
    allocate(gm%phys_name(gm%nphys))

    do i = 1, gm%nphys
       gm%phys_dim(i) = int(stream_int(s))
       gm%phys_tag(i) = int(stream_int(s))
       call stream_line(s, str, ok)
       str = adjustl(str)
       if (str(1:1) .eq. '"') then
          str = str(2:)
          q = index(str, '"', back = .true.)
          if (q .gt. 0) str = str(1:q-1)
       end if
       gm%phys_name(i) = trim(str)
    end do

  end subroutine read_physical_names

  !> Read the $Entities section (MSH 4), keeping the first physical tag of
  !! each entity
  !! @param s Input stream
  !! @param gm Mesh receiving the volume and multi-group counts
  !! @param ent Entity lists of each dimension, filled on return
  subroutine read_entities(s, gm, ent)
    type(msh_stream_t), intent(inout) :: s
    type(gmsh_mesh_t), intent(inout) :: gm
    type(entity_list_t), intent(inout) :: ent(0:3)
    integer(kind=i8) :: cnt(0:3), np, nb, j
    real(kind=dp) :: bbox(6)
    integer :: d, i, nr, ph

    do d = 0, 3
       cnt(d) = read_size_t(s)
    end do

    do d = 0, 3
       call ent(d)%free()
       ent(d)%n = int(cnt(d))
       if (d .eq. 3) gm%n_volumes = ent(d)%n
       allocate(ent(d)%tag(ent(d)%n), ent(d)%phys(ent(d)%n))
       nr = 6
       if (d .eq. 0) nr = 3
       do i = 1, ent(d)%n
          ent(d)%tag(i) = read_int(s)
          call read_reals(s, bbox(1:nr))
          np = read_size_t(s)
          ent(d)%phys(i) = 0
          do j = 1, np
             ph = read_int(s)
             if (j .eq. 1) ent(d)%phys(i) = ph
          end do
          if (np .gt. 1 .and. (d .eq. 1 .or. d .eq. 2)) then
             gm%n_multi_phys = gm%n_multi_phys + 1
          end if
          if (d .gt. 0) then
             nb = read_size_t(s)
             do j = 1, nb
                ph = read_int(s)
             end do
          end if
       end do
    end do

  end subroutine read_entities

  !> Read the $Nodes section of a MSH 4 file
  !! @param s Input stream
  !! @param gm Mesh receiving the node coordinates
  !! @param nmap Node tag map, built on return
  subroutine read_nodes_v4(s, gm, nmap)
    type(msh_stream_t), intent(inout) :: s
    type(gmsh_mesh_t), intent(inout) :: gm
    type(node_map_t), intent(inout) :: nmap
    integer(kind=i8) :: nblocks, nnodes, b, nn, k, done, m, i, tag_bound
    integer(kind=i8), allocatable :: tags(:)
    real(kind=dp), allocatable :: cbuf(:)
    integer :: edim, etag, param, npar, nc, j

    nblocks = read_size_t(s)
    nnodes = read_size_t(s)
    ! The tag range is recomputed from the tags themselves
    tag_bound = read_size_t(s)
    tag_bound = read_size_t(s)

    if (nnodes .gt. int(huge(1), i8)) then
       call fatal('Too many nodes in Gmsh file')
    end if
    gm%nnodes = int(nnodes)
    allocate(gm%xyz(3, gm%nnodes), tags(gm%nnodes))
    if (s%binary) allocate(cbuf(6 * BIN_CHUNK))

    k = 0
    do b = 1, nblocks
       edim = read_int(s)
       etag = read_int(s)
       param = read_int(s)
       nn = read_size_t(s)
       npar = 0
       if (param .ne. 0) npar = edim
       nc = 3 + npar
       if (k + nn .gt. nnodes) then
          call fatal('Inconsistent $Nodes section in Gmsh file')
       end if

       if (s%binary) then
          if (nn .gt. 0) call stream_read_size_t(s, tags(k+1:k+nn))
          done = 0
          do while (done .lt. nn)
             m = min(int(BIN_CHUNK, i8), nn - done)
             call stream_read_r8(s, cbuf(1:nc*m))
             do i = 1, m
                gm%xyz(:, k + done + i) = cbuf(nc*(i - 1) + 1:nc*(i - 1) + 3)
             end do
             done = done + m
          end do
       else
          do i = 1, nn
             tags(k + i) = stream_int(s)
          end do
          do i = 1, nn
             do j = 1, 3
                gm%xyz(j, k + i) = stream_real(s)
             end do
             do j = 1, npar
                call stream_skip_token(s)
             end do
          end do
       end if
       k = k + nn
    end do

    if (k .ne. nnodes) then
       call fatal('Inconsistent $Nodes section in Gmsh file')
    end if

    call build_node_map(tags, nmap)
    deallocate(tags)
    if (allocated(cbuf)) deallocate(cbuf)

  end subroutine read_nodes_v4

  !> Read the $Nodes section of a MSH 2 file
  !! @param s Input stream
  !! @param gm Mesh receiving the node coordinates
  !! @param nmap Node tag map, built on return
  subroutine read_nodes_v2(s, gm, nmap)
    type(msh_stream_t), intent(inout) :: s
    type(gmsh_mesh_t), intent(inout) :: gm
    type(node_map_t), intent(inout) :: nmap
    integer(kind=i8), parameter :: REC_LEN = 28
    integer(kind=i8) :: n, done, m, i, p
    integer(kind=i8), allocatable :: tags(:)
    character(len=:), allocatable :: rec
    integer(kind=i4) :: t4
    real(kind=dp) :: x3(3)

    n = stream_int(s)
    if (n .gt. int(huge(1), i8)) then
       call fatal('Too many nodes in Gmsh file')
    end if
    gm%nnodes = int(n)
    allocate(gm%xyz(3, gm%nnodes), tags(gm%nnodes))

    if (s%binary) then
       ! Each node is an int tag followed by three doubles
       call stream_skip_line(s)
       allocate(character(len=REC_LEN*BIN_CHUNK) :: rec)
       done = 0
       do while (done .lt. n)
          m = min(int(BIN_CHUNK, i8), n - done)
          call stream_read_chars(s, rec(1:REC_LEN*m))
          do i = 1, m
             p = REC_LEN * (i - 1) + 1
             t4 = transfer(rec(p:p+3), t4)
             x3 = transfer(rec(p+4:p+27), x3)
             if (s%swap) then
                t4 = swap_i4(t4)
                x3 = swap_r8(x3)
             end if
             tags(done + i) = int(t4, i8)
             gm%xyz(:, done + i) = x3
          end do
          done = done + m
       end do
    else
       do i = 1, n
          tags(i) = stream_int(s)
          gm%xyz(1, i) = stream_real(s)
          gm%xyz(2, i) = stream_real(s)
          gm%xyz(3, i) = stream_real(s)
       end do
    end if

    call build_node_map(tags, nmap)
    deallocate(tags)
    if (allocated(rec)) deallocate(rec)

  end subroutine read_nodes_v2

  !> Read the $Elements section of a MSH 4 file
  !! @param s Input stream
  !! @param gm Mesh receiving the elements
  !! @param nmap Node tag map
  !! @param ent Entity lists giving the physical tags
  subroutine read_elements_v4(s, gm, nmap, ent)
    type(msh_stream_t), intent(inout) :: s
    type(gmsh_mesh_t), intent(inout) :: gm
    type(node_map_t), intent(in) :: nmap
    type(entity_list_t), intent(in) :: ent(0:3)
    integer(kind=i8) :: nblocks, ne, b, done, m, i, off, tmp
    integer(kind=i8), allocatable :: buf(:)
    integer(kind=i8) :: ntags(27)
    integer :: edim, etag, etype, ekind, nn, phys, j

    nblocks = read_size_t(s)
    tmp = read_size_t(s)
    tmp = read_size_t(s)
    tmp = read_size_t(s)

    do b = 1, nblocks
       edim = read_int(s)
       etag = read_int(s)
       etype = read_int(s)
       ne = read_size_t(s)
       ekind = element_kind(etype)
       nn = MSH_NNODES(etype)
       phys = entity_phys(ent, edim, etag)
       call reserve(gm, ekind, phys, ne)

       if (s%binary) then
          if (.not. allocated(buf)) allocate(buf(28 * BIN_CHUNK))
          done = 0
          do while (done .lt. ne)
             m = min(int(BIN_CHUNK, i8), ne - done)
             call stream_read_size_t(s, buf(1:(nn + 1)*m))
             do i = 1, m
                off = (nn + 1) * (i - 1) + 1
                call add_element(gm, ekind, buf(off+1:off+nn), phys, nmap)
             end do
             done = done + m
          end do
       else
          do i = 1, ne
             tmp = stream_int(s)
             do j = 1, nn
                ntags(j) = stream_int(s)
             end do
             call add_element(gm, ekind, ntags(1:nn), phys, nmap)
          end do
       end if
    end do

    if (allocated(buf)) deallocate(buf)

  end subroutine read_elements_v4

  !> Read the $Elements section of a MSH 2 file
  !! @param s Input stream
  !! @param gm Mesh receiving the elements
  !! @param nmap Node tag map
  subroutine read_elements_v2(s, gm, nmap)
    type(msh_stream_t), intent(inout) :: s
    type(gmsh_mesh_t), intent(inout) :: gm
    type(node_map_t), intent(in) :: nmap
    integer(kind=i8) :: n, k, tmp
    integer(kind=i8) :: ntags(27)
    integer(kind=i4) :: hdr(3)
    integer(kind=i4), allocatable :: buf(:)
    integer :: etype, nfollow, ntag, ekind, nn, phys, reclen
    integer :: done, m, i, j, off

    n = stream_int(s)

    if (s%binary) then
       ! Blocks of elements of one type, each preceded by a header with the
       ! type, the number of elements and the number of tags
       call stream_skip_line(s)
       k = 0
       do while (k .lt. n)
          call stream_read_i4(s, hdr)
          etype = hdr(1)
          nfollow = hdr(2)
          ntag = hdr(3)
          ekind = element_kind(etype)
          nn = MSH_NNODES(etype)
          reclen = 1 + ntag + nn
          call reserve(gm, ekind, 1, int(nfollow, i8))
          done = 0
          do while (done .lt. nfollow)
             m = min(BIN_CHUNK, nfollow - done)
             if (allocated(buf)) then
                if (size(buf) .lt. reclen * m) deallocate(buf)
             end if
             if (.not. allocated(buf)) allocate(buf(reclen * m))
             call stream_read_i4(s, buf(1:reclen*m))
             do i = 1, m
                off = reclen * (i - 1)
                phys = 0
                if (ntag .ge. 1) phys = buf(off + 2)
                do j = 1, nn
                   ntags(j) = int(buf(off + 1 + ntag + j), i8)
                end do
                call add_element(gm, ekind, ntags(1:nn), phys, nmap)
             end do
             done = done + m
          end do
          k = k + nfollow
       end do
    else
       do k = 1, n
          tmp = stream_int(s)
          etype = int(stream_int(s))
          ntag = int(stream_int(s))
          phys = 0
          do j = 1, ntag
             tmp = stream_int(s)
             if (j .eq. 1) phys = int(tmp)
          end do
          ekind = element_kind(etype)
          nn = MSH_NNODES(etype)
          do j = 1, nn
             ntags(j) = stream_int(s)
          end do
          call add_element(gm, ekind, ntags(1:nn), phys, nmap)
       end do
    end if

    if (allocated(buf)) deallocate(buf)

  end subroutine read_elements_v2

  !> Classify a Gmsh element type, stopping on types Neko cannot use
  !! @param etype Gmsh element type number
  function element_kind(etype) result(ekind)
    integer, intent(in) :: etype
    integer :: ekind
    character(len=16) :: str

    ekind = KIND_IGNORE
    select case (etype)
    case (5, 12, 17)
       ekind = KIND_HEX
    case (3, 10, 16)
       ekind = KIND_QUAD
    case (1, 8)
       ekind = KIND_LINE
    case (15)
       ekind = KIND_IGNORE
    case (2, 9, 20, 21, 22, 23, 24, 25)
       call fatal('The Gmsh file contains triangles. Neko needs ' // &
            'a mesh of quadrilaterals (2D) or hexahedra (3D)')
    case (4, 11, 29, 30, 31)
       call fatal('The Gmsh file contains tetrahedra. Neko needs ' // &
            'an all-hexahedral mesh')
    case (6, 13, 18)
       call fatal('The Gmsh file contains prisms. Neko needs ' // &
            'an all-hexahedral mesh')
    case (7, 14, 19)
       call fatal('The Gmsh file contains pyramids. Neko needs ' // &
            'an all-hexahedral mesh')
    case (26, 27, 28, 36, 37, 38, 39, 40, 41, 92, 93)
       call fatal('The Gmsh file contains elements of order 3 or ' // &
            'higher. Only first and second order meshes are supported')
    case default
       write(str, '(I0)') etype
       call fatal('Unsupported Gmsh element type ' // trim(str))
    end select

  end function element_kind

  !> First physical tag of an entity, 0 if it is in no physical group
  !! @param ent Entity lists of each dimension
  !! @param dim Entity dimension
  !! @param tag Entity tag
  pure function entity_phys(ent, dim, tag) result(phys)
    type(entity_list_t), intent(in) :: ent(0:3)
    integer, intent(in) :: dim, tag
    integer :: phys
    integer :: i

    phys = 0
    if (dim .lt. 0 .or. dim .gt. 3) return
    if (.not. allocated(ent(dim)%tag)) return
    do i = 1, ent(dim)%n
       if (ent(dim)%tag(i) .eq. tag) then
          phys = ent(dim)%phys(i)
          return
       end if
    end do

  end function entity_phys

  !> Make room for @a n more elements of kind @a ekind
  !! @param gm Mesh whose element arrays are grown
  !! @param ekind Element kind
  !! @param phys Physical tag, lines without one are not stored
  !! @param n Number of elements to make room for
  subroutine reserve(gm, ekind, phys, n)
    type(gmsh_mesh_t), intent(inout) :: gm
    integer, intent(in) :: ekind, phys
    integer(kind=i8), intent(in) :: n
    integer(kind=i8) :: need

    select case (ekind)
    case (KIND_HEX)
       need = int(gm%nhex, i8) + n
    case (KIND_QUAD)
       need = int(gm%nquad, i8) + n
    case (KIND_LINE)
       if (phys .le. 0) return
       need = int(gm%nline, i8) + n
    case default
       return
    end select

    if (need .gt. int(huge(1), i8)) then
       call fatal('Too many elements in Gmsh file')
    end if

    select case (ekind)
    case (KIND_HEX)
       call grow2(gm%hex, 20, int(need))
    case (KIND_QUAD)
       call grow2(gm%quad, 8, int(need))
       call grow1(gm%quad_phys, int(need))
    case (KIND_LINE)
       call grow2(gm%line, 3, int(need))
       call grow1(gm%line_phys, int(need))
    end select

  end subroutine reserve

  !> Store one element given by its Gmsh node tags
  !! @param gm Mesh receiving the element
  !! @param ekind Element kind
  !! @param tags Gmsh node tags of the element
  !! @param phys Physical tag of the element, 0 if none
  !! @param nmap Node tag map
  subroutine add_element(gm, ekind, tags, phys, nmap)
    type(gmsh_mesh_t), intent(inout) :: gm
    integer, intent(in) :: ekind
    integer(kind=i8), intent(in) :: tags(:)
    integer, intent(in) :: phys
    type(node_map_t), intent(in) :: nmap
    integer :: i, n

    select case (ekind)
    case (KIND_HEX)
       n = gm%nhex + 1
       call grow2(gm%hex, 20, n)
       gm%hex(:, n) = 0
       do i = 1, min(size(tags), 20)
          gm%hex(HEX_G2S(i), n) = node_index(nmap, tags(i))
       end do
       if (size(tags) .gt. 8) gm%nhex_ho = gm%nhex_ho + 1
       gm%nhex = n
    case (KIND_QUAD)
       n = gm%nquad + 1
       call grow2(gm%quad, 8, n)
       call grow1(gm%quad_phys, n)
       gm%quad(:, n) = 0
       do i = 1, min(size(tags), 8)
          gm%quad(QUAD_G2S(i), n) = node_index(nmap, tags(i))
       end do
       gm%quad_phys(n) = phys
       if (size(tags) .gt. 4) gm%nquad_ho = gm%nquad_ho + 1
       gm%nquad = n
    case (KIND_LINE)
       ! Lines only matter as boundary facets of 2D meshes
       if (phys .le. 0) return
       n = gm%nline + 1
       call grow2(gm%line, 3, n)
       call grow1(gm%line_phys, n)
       gm%line(:, n) = 0
       do i = 1, min(size(tags), 3)
          gm%line(i, n) = node_index(nmap, tags(i))
       end do
       gm%line_phys(n) = phys
       gm%nline = n
    end select

  end subroutine add_element

  !> Build the map from node tags to node indices
  !! @param tags Node tags in the order of the nodes
  !! @param nmap Map, built on return
  subroutine build_node_map(tags, nmap)
    integer(kind=i8), intent(in) :: tags(:)
    type(node_map_t), intent(inout) :: nmap
    integer(kind=i8) :: tmin, tmax, range, j
    integer :: i

    if (size(tags) .eq. 0) then
       call fatal('The Gmsh file contains no nodes')
    end if

    tmin = minval(tags)
    tmax = maxval(tags)
    if (tmin .lt. 1) then
       call fatal('Invalid node tag in Gmsh file')
    end if
    range = tmax - tmin + 1
    if (range .gt. int(huge(1), i8) .or. &
         range .gt. 16_i8 * size(tags, kind=i8) + 10000000_i8) then
       call fatal('The node tags of the Gmsh file are too sparse, ' // &
            'renumber the mesh in Gmsh before saving it')
    end if

    call nmap%free()
    allocate(nmap%idx(range))
    nmap%idx = 0
    nmap%tag_min = tmin
    do i = 1, size(tags)
       j = tags(i) - tmin + 1
       if (nmap%idx(j) .ne. 0) then
          call fatal('Duplicate node tag in Gmsh file')
       end if
       nmap%idx(j) = i
    end do

  end subroutine build_node_map

  !> Node index of a node tag
  !! @param nmap Node tag map
  !! @param tag Gmsh node tag
  function node_index(nmap, tag) result(idx)
    type(node_map_t), intent(in) :: nmap
    integer(kind=i8), intent(in) :: tag
    integer :: idx
    integer(kind=i8) :: j
    character(len=32) :: str

    j = tag - nmap%tag_min + 1
    idx = 0
    if (j .ge. 1 .and. j .le. size(nmap%idx, kind=i8)) idx = nmap%idx(j)
    if (idx .eq. 0) then
       write(str, '(I0)') tag
       call fatal('Element refers to unknown node ' // trim(str))
    end if

  end function node_index

  !> Make sure a 2D array has at least @a need columns
  !! @param a Array to grow
  !! @param nrow Number of rows
  !! @param need Number of columns needed
  subroutine grow2(a, nrow, need)
    integer, allocatable, intent(inout) :: a(:,:)
    integer, intent(in) :: nrow, need
    integer, allocatable :: tmp(:,:)
    integer :: cap

    if (.not. allocated(a)) then
       allocate(a(nrow, max(need, 1024)))
    else if (size(a, 2) .lt. need) then
       cap = int(min(int(huge(1), i8), &
            max(int(need, i8), 2_i8 * size(a, 2, kind=i8))))
       allocate(tmp(nrow, cap))
       tmp(:, 1:size(a, 2)) = a
       call move_alloc(tmp, a)
    end if

  end subroutine grow2

  !> Make sure a 1D array has at least @a need entries
  !! @param a Array to grow
  !! @param need Number of entries needed
  subroutine grow1(a, need)
    integer, allocatable, intent(inout) :: a(:)
    integer, intent(in) :: need
    integer, allocatable :: tmp(:)
    integer :: cap

    if (.not. allocated(a)) then
       allocate(a(max(need, 1024)))
    else if (size(a) .lt. need) then
       cap = int(min(int(huge(1), i8), &
            max(int(need, i8), 2_i8 * size(a, kind=i8))))
       allocate(tmp(cap))
       tmp(1:size(a)) = a
       call move_alloc(tmp, a)
    end if

  end subroutine grow1

  ! ---------------------------------------------------------------------
  ! Value readers that follow the ASCII/binary mode of MSH 4 sections
  ! ---------------------------------------------------------------------

  !> Read an int
  !! @param s Input stream
  function read_int(s) result(v)
    type(msh_stream_t), intent(inout) :: s
    integer :: v
    integer(kind=i4) :: b(1)

    if (s%binary) then
       call stream_read_i4(s, b)
       v = int(b(1))
    else
       v = int(stream_int(s))
    end if

  end function read_int

  !> Read a size_t
  !! @param s Input stream
  function read_size_t(s) result(v)
    type(msh_stream_t), intent(inout) :: s
    integer(kind=i8) :: v
    integer(kind=i8) :: b(1)

    if (s%binary) then
       call stream_read_size_t(s, b)
       v = b(1)
    else
       v = stream_int(s)
    end if

  end function read_size_t

  !> Read an array of doubles
  !! @param s Input stream
  !! @param v Values read
  subroutine read_reals(s, v)
    type(msh_stream_t), intent(inout) :: s
    real(kind=dp), intent(out) :: v(:)
    integer :: i

    if (s%binary) then
       call stream_read_r8(s, v)
    else
       do i = 1, size(v)
          v(i) = stream_real(s)
       end do
    end if

  end subroutine read_reals

  ! ---------------------------------------------------------------------
  ! Buffered stream
  ! ---------------------------------------------------------------------

  !> Open a file as a byte stream
  !! @param s Stream to open
  !! @param fname Name of the file
  subroutine stream_open(s, fname)
    type(msh_stream_t), intent(inout) :: s
    character(len=*), intent(in) :: fname
    integer :: ierr

    open(newunit = s%unit, file = trim(fname), access = 'stream', &
         form = 'unformatted', status = 'old', action = 'read', &
         iostat = ierr)
    if (ierr .ne. 0) then
       call fatal('Cannot open ' // trim(fname))
    end if
    inquire(unit = s%unit, size = s%fsize)
    allocate(character(len=STREAM_BUF_LEN) :: s%buf)
    s%pos = 1
    s%buf_start = 1
    s%buf_len = 0

  end subroutine stream_open

  !> Make sure the byte at the current position is buffered
  !! @return .false. at the end of the file
  !! @param s Stream
  function stream_fill(s) result(ok)
    type(msh_stream_t), intent(inout) :: s
    logical :: ok
    integer :: n, ierr

    ok = .false.
    if (s%pos .gt. s%fsize) return
    if (s%pos .lt. s%buf_start .or. &
         s%pos .ge. s%buf_start + int(s%buf_len, i8)) then
       n = int(min(int(STREAM_BUF_LEN, i8), s%fsize - s%pos + 1_i8))
       read(s%unit, pos = s%pos, iostat = ierr) s%buf(1:n)
       if (ierr .ne. 0) call fatal('Read error in Gmsh file')
       s%buf_start = s%pos
       s%buf_len = n
    end if
    ok = .true.

  end function stream_fill

  !> Whitespace test
  !! @param c Character to test
  pure function is_space(c) result(res)
    character, intent(in) :: c
    logical :: res

    res = c .eq. ' ' .or. c .eq. achar(9) .or. c .eq. achar(10) .or. &
         c .eq. achar(13)

  end function is_space

  !> Read the next whitespace separated token
  !! @param s Stream
  !! @param tok The token, truncated to len(tok)
  !! @param n Length of the token, 0 at the end of the file
  subroutine stream_token(s, tok, n)
    type(msh_stream_t), intent(inout) :: s
    character(len=*), intent(out) :: tok
    integer, intent(out) :: n
    integer :: i

    tok = ''
    n = 0

    ! Skip leading whitespace
    do
       if (.not. stream_fill(s)) return
       i = int(s%pos - s%buf_start) + 1
       do while (i .le. s%buf_len)
          if (.not. is_space(s%buf(i:i))) exit
          i = i + 1
       end do
       s%pos = s%buf_start + int(i - 1, i8)
       if (i .le. s%buf_len) exit
    end do

    ! Collect the token
    do
       if (.not. stream_fill(s)) return
       i = int(s%pos - s%buf_start) + 1
       do while (i .le. s%buf_len)
          if (is_space(s%buf(i:i))) exit
          n = n + 1
          if (n .le. len(tok)) tok(n:n) = s%buf(i:i)
          i = i + 1
       end do
       s%pos = s%buf_start + int(i - 1, i8)
       if (i .le. s%buf_len) return
    end do

  end subroutine stream_token

  !> Skip the next token
  !! @param s Stream
  subroutine stream_skip_token(s)
    type(msh_stream_t), intent(inout) :: s
    character(len=1) :: tok
    integer :: n

    call stream_token(s, tok, n)
    if (n .eq. 0) call fatal('Unexpected end of Gmsh file')

  end subroutine stream_skip_token

  !> Read the next token as an integer
  !! @param s Stream
  function stream_int(s) result(v)
    type(msh_stream_t), intent(inout) :: s
    integer(kind=i8) :: v
    character(len=64) :: tok
    integer :: n, i, i0, d
    logical :: neg

    call stream_token(s, tok, n)
    if (n .eq. 0) call fatal('Unexpected end of Gmsh file')
    if (n .gt. len(tok)) call fatal('Invalid integer in Gmsh file')

    v = 0
    neg = tok(1:1) .eq. '-'
    i0 = 1
    if (tok(1:1) .eq. '-' .or. tok(1:1) .eq. '+') i0 = 2
    if (i0 .gt. n) call fatal('Invalid integer in Gmsh file')
    do i = i0, n
       d = iachar(tok(i:i)) - iachar('0')
       if (d .lt. 0 .or. d .gt. 9) then
          call fatal('Invalid integer in Gmsh file: ' // tok(1:n))
       end if
       v = 10_i8 * v + int(d, i8)
    end do
    if (neg) v = -v

  end function stream_int

  !> Read the next token as a real
  !! @param s Stream
  function stream_real(s) result(v)
    type(msh_stream_t), intent(inout) :: s
    real(kind=dp) :: v
    character(len=64) :: tok
    integer :: n, ierr

    call stream_token(s, tok, n)
    if (n .eq. 0) call fatal('Unexpected end of Gmsh file')
    read(tok(1:min(n, len(tok))), *, iostat = ierr) v
    if (ierr .ne. 0 .or. n .gt. len(tok)) then
       call fatal('Invalid real in Gmsh file: ' // trim(tok))
    end if

  end function stream_real

  !> Read the rest of the current line, without the line terminator
  !! @param s Stream
  !! @param str The line, truncated to len(str)
  !! @param ok .false. at the end of the file
  subroutine stream_line(s, str, ok)
    type(msh_stream_t), intent(inout) :: s
    character(len=*), intent(out) :: str
    logical, intent(out) :: ok
    integer :: i, n
    logical :: found

    str = ''
    n = 0
    ok = stream_fill(s)
    if (.not. ok) return

    found = .false.
    do while (.not. found)
       if (.not. stream_fill(s)) exit
       i = int(s%pos - s%buf_start) + 1
       do while (i .le. s%buf_len)
          if (s%buf(i:i) .eq. achar(10)) then
             found = .true.
             exit
          end if
          if (s%buf(i:i) .ne. achar(13)) then
             n = n + 1
             if (n .le. len(str)) str(n:n) = s%buf(i:i)
          end if
          i = i + 1
       end do
       if (found) then
          s%pos = s%buf_start + int(i, i8)
       else
          s%pos = s%buf_start + int(i - 1, i8)
       end if
    end do

  end subroutine stream_line

  !> Skip the rest of the current line
  !! @param s Stream
  subroutine stream_skip_line(s)
    type(msh_stream_t), intent(inout) :: s
    character(len=1) :: str
    logical :: ok

    call stream_line(s, str, ok)

  end subroutine stream_skip_line

  !> Find the next section header
  !! @param s Stream
  !! @param name Section name without the leading $
  !! @param ok .false. at the end of the file
  subroutine stream_section(s, name, ok)
    type(msh_stream_t), intent(inout) :: s
    character(len=*), intent(out) :: name
    logical, intent(out) :: ok
    character(len=256) :: str

    name = ''
    do
       call stream_line(s, str, ok)
       if (.not. ok) return
       str = adjustl(str)
       if (len_trim(str) .eq. 0) cycle
       if (str(1:1) .ne. '$') then
          call fatal('Malformed Gmsh file, expected a section ' // &
               'header but found: ' // trim(str))
       end if
       name = str(2:)
       return
    end do

  end subroutine stream_section

  !> Consume the end marker of section @a name
  !! @param s Stream
  !! @param name Section name without the leading $
  subroutine stream_end_section(s, name)
    type(msh_stream_t), intent(inout) :: s
    character(len=*), intent(in) :: name
    character(len=256) :: str
    logical :: ok

    do
       call stream_line(s, str, ok)
       if (.not. ok) then
          call fatal('Missing $End' // trim(name) // ' in Gmsh file')
       end if
       str = adjustl(str)
       if (len_trim(str) .eq. 0) cycle
       if (trim(str) .ne. '$End' // trim(name)) then
          call fatal('Malformed $' // trim(name) // &
               ' section in Gmsh file')
       end if
       return
    end do

  end subroutine stream_end_section

  !> Skip a section that is not needed
  !! @param s Stream
  !! @param name Section name without the leading $
  subroutine stream_skip_section(s, name)
    type(msh_stream_t), intent(inout) :: s
    character(len=*), intent(in) :: name
    character(len=256) :: str
    logical :: ok

    do
       call stream_line(s, str, ok)
       if (.not. ok) then
          call fatal('Missing $End' // trim(name) // ' in Gmsh file')
       end if
       if (trim(adjustl(str)) .eq. '$End' // trim(name)) return
    end do

  end subroutine stream_skip_section

  !> Read raw bytes
  !! @param s Stream
  !! @param str Bytes read, len(str) of them
  subroutine stream_read_chars(s, str)
    type(msh_stream_t), intent(inout) :: s
    character(len=*), intent(out) :: str
    integer :: ierr

    if (len(str) .eq. 0) return
    read(s%unit, pos = s%pos, iostat = ierr) str
    if (ierr .ne. 0) then
       call fatal('Unexpected end of binary data in Gmsh file')
    end if
    s%pos = s%pos + int(len(str), i8)

  end subroutine stream_read_chars

  !> Read binary 4-byte integers
  !! @param s Stream
  !! @param a Values read
  subroutine stream_read_i4(s, a)
    type(msh_stream_t), intent(inout) :: s
    integer(kind=i4), intent(out) :: a(:)
    integer :: ierr

    if (size(a) .eq. 0) return
    read(s%unit, pos = s%pos, iostat = ierr) a
    if (ierr .ne. 0) then
       call fatal('Unexpected end of binary data in Gmsh file')
    end if
    s%pos = s%pos + 4_i8 * size(a, kind=i8)
    if (s%swap) a = swap_i4(a)

  end subroutine stream_read_i4

  !> Read binary doubles
  !! @param s Stream
  !! @param a Values read
  subroutine stream_read_r8(s, a)
    type(msh_stream_t), intent(inout) :: s
    real(kind=dp), intent(out) :: a(:)
    integer :: ierr

    if (size(a) .eq. 0) return
    read(s%unit, pos = s%pos, iostat = ierr) a
    if (ierr .ne. 0) then
       call fatal('Unexpected end of binary data in Gmsh file')
    end if
    s%pos = s%pos + 8_i8 * size(a, kind=i8)
    if (s%swap) a = swap_r8(a)

  end subroutine stream_read_r8

  !> Read binary size_t values
  !! @param s Stream
  !! @param a Values read
  subroutine stream_read_size_t(s, a)
    type(msh_stream_t), intent(inout) :: s
    integer(kind=i8), intent(out) :: a(:)
    integer(kind=i4), allocatable :: tmp(:)
    integer :: ierr

    if (size(a) .eq. 0) return
    if (s%size_t_len .eq. 8) then
       read(s%unit, pos = s%pos, iostat = ierr) a
       if (ierr .ne. 0) then
          call fatal('Unexpected end of binary data in Gmsh file')
       end if
       s%pos = s%pos + 8_i8 * size(a, kind=i8)
       if (s%swap) a = swap_i8(a)
    else
       allocate(tmp(size(a)))
       call stream_read_i4(s, tmp)
       a = int(tmp, i8)
       deallocate(tmp)
    end if

  end subroutine stream_read_size_t

  !> Reverse the byte order of a 4-byte integer
  !! @param x Value to swap
  elemental function swap_i4(x) result(y)
    integer(kind=i4), intent(in) :: x
    integer(kind=i4) :: y
    character(len=4) :: c, r
    integer :: k

    c = transfer(x, c)
    do k = 1, 4
       r(k:k) = c(5-k:5-k)
    end do
    y = transfer(r, y)

  end function swap_i4

  !> Reverse the byte order of an 8-byte integer
  !! @param x Value to swap
  elemental function swap_i8(x) result(y)
    integer(kind=i8), intent(in) :: x
    integer(kind=i8) :: y
    character(len=8) :: c, r
    integer :: k

    c = transfer(x, c)
    do k = 1, 8
       r(k:k) = c(9-k:9-k)
    end do
    y = transfer(r, y)

  end function swap_i8

  !> Reverse the byte order of a double
  !! @param x Value to swap
  elemental function swap_r8(x) result(y)
    real(kind=dp), intent(in) :: x
    real(kind=dp) :: y
    character(len=8) :: c, r
    integer :: k

    c = transfer(x, c)
    do k = 1, 8
       r(k:k) = c(9-k:9-k)
    end do
    y = transfer(r, y)

  end function swap_r8

end module gmsh2nmsh_reader

!> Convert a Gmsh mesh to Neko's nmsh format
!! @details Reads a mesh of first or second order hexahedra (3D) or
!! quadrilaterals (2D) in Gmsh's MSH 2.2 or 4.1 format, ASCII or binary,
!! and writes it as a Neko mesh:
!! - boundary facets in physical groups become labeled zones,
!! - pairs of physical groups can be made periodic (pure translations),
!! - curved edges of second order elements are kept as midpoint curves,
!! - left-handed elements are mirrored to right-handed ones.
!!
!! The nmsh file is written directly, with the layout of Neko's
!! nmsh_file_write, so the tool needs neither MPI nor the Neko library.
program gmsh2nmsh
  use gmsh2nmsh_util, only : i4, i8, dp, fatal, warn, progress, &
       int_str, real_str
  use gmsh2nmsh_reader, only : gmsh_mesh_t, gmsh_read, GMSH_NAME_LEN
  implicit none

  !> Maximum length of a file name
  integer, parameter :: FNAME_LEN = 1024
  !> Number of labeled zones Neko supports (NEKO_MSH_MAX_ZLBLS)
  integer, parameter :: MAX_ZONES = 20
  !> nmsh zone types
  integer, parameter :: PERIODIC_ZONE = 5
  integer, parameter :: LABELED_ZONE = 7
  !> Symmetric number of each vertex of an nmsh element record, which
  !! stores the vertices in cyclic order
  integer, parameter :: CYC_TO_SYM(8) = [1, 2, 4, 3, 5, 6, 8, 7]

  !> Vertices of the facets of a hexahedron (symmetric numbering)
  integer, parameter :: HEX_FACE(4, 6) = reshape([1, 5, 7, 3, 2, 6, 8, 4, &
       1, 2, 6, 5, 3, 4, 8, 7, 1, 2, 4, 3, 5, 6, 8, 7], [4, 6])
  !> Vertices of the facets of a quadrilateral
  integer, parameter :: QUAD_FACE(2, 4) = reshape([1, 3, 2, 4, 1, 2, 3, 4], &
       [2, 4])
  !> End vertices of the curve edges of a hexahedron
  integer, parameter :: HEX_EDGE(2, 12) = reshape([1, 2, 2, 4, 4, 3, 3, 1, &
       5, 6, 6, 8, 8, 7, 7, 5, 1, 5, 2, 6, 4, 8, 3, 7], [2, 12])
  !> End vertices of the curve edges of a quadrilateral
  integer, parameter :: QUAD_EDGE(2, 4) = reshape([1, 2, 2, 4, 4, 3, 3, 1], &
       [2, 4])
  !> Slot permutation that mirrors a hexahedron in t
  integer, parameter :: HEX_FLIP(20) = [5, 6, 7, 8, 1, 2, 3, 4, &
       13, 14, 15, 16, 9, 10, 11, 12, 17, 18, 19, 20]
  !> Slot permutation that mirrors a quadrilateral in s
  integer, parameter :: QUAD_FLIP(8) = [3, 4, 1, 2, 7, 6, 5, 8]
  !> Default relative midpoint offset above which an edge is curved
  real(kind=dp), parameter :: DEFAULT_CURVE_TOL = 1.0e-4_dp

  ! Command line
  character(len=FNAME_LEN) :: input_fname, output_fname
  character(len=4096) :: periodic_spec
  real(kind=dp) :: user_tol, curve_tol
  logical :: user_tol_set, curve_tol_set, linear

  type(gmsh_mesh_t) :: gm

  ! Elements: corner and edge midpoint node indices per element, see the
  ! slot layout in gmsh2nmsh_reader
  integer :: gdim, nel, ncrn, nedge, nfacet, nfv
  integer, allocatable :: conn(:,:)
  integer, allocatable :: facet_vtx(:,:), edge_vtx(:,:), flip(:)
  integer :: nflip, nbad

  ! Compact vertex ids of the element corners, per Gmsh node
  integer :: nvert
  integer, allocatable :: vid(:)

  ! Boundary facets from physical groups
  integer :: nbf, n_bf_unknown, n_bf_dup, n_bf_conflict
  integer, allocatable :: bf_key(:,:), bf_phys(:), bf_nmatch(:), bf_grp(:)
  integer, allocatable :: bf_el(:,:), bf_side(:,:)
  logical, allocatable :: bf_keep(:)

  ! Boundary physical groups, sorted by tag
  integer :: ngrp
  integer, allocatable :: grp_tag(:), grp_label(:), grp_partner(:)
  integer, allocatable :: grp_nfacets(:), grp_ninternal(:)
  character(len=GMSH_NAME_LEN), allocatable :: grp_name(:)
  logical :: renumbered

  ! Zone records, see mark_zones
  integer :: nlab, nper
  integer, allocatable :: lab_e(:), lab_f(:), lab_label(:)
  integer, allocatable :: per_e(:), per_f(:), per_pe(:), per_pf(:)

  ! Union-find over vertex ids, used to merge periodic vertices
  integer, allocatable :: parent(:)

  ! Curved elements
  integer :: ncurved_el, ncurved_edges

  ! Clock at the start of the conversion
  integer(kind=i8) :: clock_start

  call print_banner()
  call parse_arguments()
  call print_settings()

  call system_clock(clock_start)
  call gmsh_read(trim(input_fname), gm)
  call setup_elements()
  call number_vertices()
  call print_elements()
  call match_boundary_facets()
  call setup_groups()
  call print_groups()
  call mark_zones()
  call count_curves()
  call print_curves()
  call write_nmsh(trim(output_fname))

  call print_table()
  call print_warnings()

  call gm%free()
  call free_mesh_data()

contains

  !> Deallocate the arrays of the converted mesh
  subroutine free_mesh_data()

    if (allocated(conn)) deallocate(conn)
    if (allocated(facet_vtx)) deallocate(facet_vtx)
    if (allocated(edge_vtx)) deallocate(edge_vtx)
    if (allocated(flip)) deallocate(flip)
    if (allocated(vid)) deallocate(vid)
    if (allocated(bf_key)) deallocate(bf_key)
    if (allocated(bf_phys)) deallocate(bf_phys)
    if (allocated(bf_nmatch)) deallocate(bf_nmatch)
    if (allocated(bf_grp)) deallocate(bf_grp)
    if (allocated(bf_el)) deallocate(bf_el)
    if (allocated(bf_side)) deallocate(bf_side)
    if (allocated(bf_keep)) deallocate(bf_keep)
    if (allocated(grp_tag)) deallocate(grp_tag)
    if (allocated(grp_label)) deallocate(grp_label)
    if (allocated(grp_partner)) deallocate(grp_partner)
    if (allocated(grp_nfacets)) deallocate(grp_nfacets)
    if (allocated(grp_ninternal)) deallocate(grp_ninternal)
    if (allocated(grp_name)) deallocate(grp_name)
    if (allocated(lab_e)) deallocate(lab_e)
    if (allocated(lab_f)) deallocate(lab_f)
    if (allocated(lab_label)) deallocate(lab_label)
    if (allocated(per_e)) deallocate(per_e)
    if (allocated(per_f)) deallocate(per_f)
    if (allocated(per_pe)) deallocate(per_pe)
    if (allocated(per_pf)) deallocate(per_pf)
    if (allocated(parent)) deallocate(parent)

  end subroutine free_mesh_data

  !> Print the banner with the version and the compiler
  subroutine print_banner()
    use, intrinsic :: iso_fortran_env, only : compiler_version

    write(*, '(A)') ''
    write(*, '(A)') ' N E K O  Gmsh to nmsh converter, version ' // &
         PACKAGE_VERSION
    write(*, '(A)') ' Built with ' // compiler_version()
    write(*, '(A)') ''

  end subroutine print_banner

  !> Print the input, the output and the options in use
  subroutine print_settings()

    write(*, '(A)') ' Input:      ' // trim(input_fname)
    write(*, '(A)') ' Output:     ' // trim(output_fname)
    if (len_trim(periodic_spec) .gt. 0) then
       write(*, '(A)') ' Periodic:   ' // trim(adjustl(periodic_spec))
    end if
    if (user_tol_set) then
       write(*, '(A)') ' Tolerance:  ' // real_str(user_tol)
    end if
    if (curve_tol_set) then
       write(*, '(A)') ' Curve tol:  ' // real_str(curve_tol)
    end if
    if (linear) then
       write(*, '(A)') ' Linear:     yes'
    end if
    write(*, '(A)') ''

  end subroutine print_settings

  !> Print the command line usage
  subroutine print_usage()
    write(*, '(A)') 'Usage: gmsh2nmsh <mesh.msh> [<mesh.nmsh>] [options]'
    write(*, '(A)') ''
    write(*, '(A)') 'Convert a Gmsh mesh (MSH 4.1 or 2.2, ASCII or ' // &
         'binary) of first or second'
    write(*, '(A)') 'order hexahedra (3D) or quadrilaterals (2D) to ' // &
         'the Neko nmsh format.'
    write(*, '(A)') 'Boundary facets in physical groups become labeled ' // &
         'zones.'
    write(*, '(A)') ''
    write(*, '(A)') 'Options:'
    write(*, '(A)') '  --periodic=PAIRS   Make pairs of physical groups ' // &
         'periodic, given by tag'
    write(*, '(A)') '                     or name, e.g. "1:2" or ' // &
         '"inlet:outlet,front:back"'
    write(*, '(A)') '  --tol=VALUE        Absolute periodic matching ' // &
         'tolerance (default:'
    write(*, '(A)') '                     1e-6 times the shortest facet ' // &
         'edge of the pair)'
    write(*, '(A)') '  --curve-tol=VALUE  Midpoint offset, relative to ' // &
         'the edge length,'
    write(*, '(A)') '                     above which an edge is stored ' // &
         'as curved (default: 1e-4)'
    write(*, '(A)') '  --linear           Ignore high-order nodes, ' // &
         'write straight elements'
    write(*, '(A)') '  -h, --help         Show this help'
  end subroutine print_usage

  !> Parse the command line
  subroutine parse_arguments()
    character(len=FNAME_LEN) :: arg, val
    character(len=80) :: suffix
    integer :: argc, i, npos, ios

    periodic_spec = ''
    user_tol = 0.0_dp
    user_tol_set = .false.
    curve_tol = DEFAULT_CURVE_TOL
    curve_tol_set = .false.
    linear = .false.
    npos = 0

    argc = command_argument_count()
    if (argc .lt. 1) then
       call print_usage()
       stop
    end if

    i = 1
    do while (i .le. argc)
       call get_command_argument(i, arg)
       if (trim(arg) .eq. '-h' .or. trim(arg) .eq. '--help') then
          call print_usage()
          stop
       else if (option_value(arg, '--periodic', i, argc, val)) then
          periodic_spec = trim(periodic_spec) // ' ' // trim(val)
       else if (option_value(arg, '--tol', i, argc, val)) then
          read(val, *, iostat = ios) user_tol
          if (ios .ne. 0 .or. user_tol .le. 0.0_dp) then
             call fatal('Invalid value for --tol: ' // trim(val))
          end if
          user_tol_set = .true.
       else if (option_value(arg, '--curve-tol', i, argc, val)) then
          read(val, *, iostat = ios) curve_tol
          if (ios .ne. 0 .or. curve_tol .lt. 0.0_dp) then
             call fatal('Invalid value for --curve-tol: ' // trim(val))
          end if
          curve_tol_set = .true.
       else if (trim(arg) .eq. '--linear') then
          linear = .true.
       else if (arg(1:2) .eq. '--') then
          call fatal('Unknown option ' // trim(arg))
       else
          npos = npos + 1
          if (npos .eq. 1) then
             input_fname = arg
          else if (npos .eq. 2) then
             output_fname = arg
          else
             call fatal('Too many arguments, see gmsh2nmsh --help')
          end if
       end if
       i = i + 1
    end do

    if (npos .lt. 1) then
       call fatal('No input file given, see gmsh2nmsh --help')
    end if

    suffix = file_suffix(input_fname)
    if (trim(suffix) .ne. 'msh') then
       call fatal('The input file must be a Gmsh .msh file')
    end if

    if (npos .ge. 2) then
       suffix = file_suffix(output_fname)
       if (trim(suffix) .ne. 'nmsh') then
          call fatal('The output file must end in .nmsh')
       end if
    else
       output_fname = input_fname(1:scan(trim(input_fname), '.', &
            back = .true.)) // 'nmsh'
    end if

  end subroutine parse_arguments

  !> Text after the last dot of a file name
  !! @param fname File name
  function file_suffix(fname) result(suffix)
    character(len=*), intent(in) :: fname
    character(len=80) :: suffix
    integer :: i

    i = scan(trim(fname), '.', back = .true.)
    suffix = ''
    if (i .gt. 0) suffix = fname(i+1:len_trim(fname))

  end function file_suffix

  !> Match an option given as "name=value" or "name value"
  !! @param arg Current command line argument
  !! @param name Option name including the leading dashes
  !! @param i Index of the current argument, advanced if the value follows it
  !! @param argc Number of command line arguments
  !! @param val Option value
  function option_value(arg, name, i, argc, val) result(match)
    character(len=*), intent(in) :: arg, name
    integer, intent(inout) :: i
    integer, intent(in) :: argc
    character(len=*), intent(out) :: val
    logical :: match
    integer :: n

    n = len_trim(name)
    match = .false.
    val = ''
    if (trim(arg) .eq. name) then
       if (i .ge. argc) call fatal('Missing value for ' // name)
       i = i + 1
       call get_command_argument(i, val)
       match = .true.
    else if (len_trim(arg) .gt. n + 1) then
       if (arg(1:n+1) .eq. name // '=') then
          val = arg(n+2:)
          match = .true.
       end if
    end if

  end function option_value

  !> Pick the elements, put them in Neko order and make them right-handed
  subroutine setup_elements()
    real(kind=dp) :: c(3, 8)
    integer :: e

    if (gm%nhex .gt. 0) then
       gdim = 3
       nel = gm%nhex
       ncrn = 8
       nedge = 12
       nfacet = 6
       nfv = 4
       conn = gm%hex(:, 1:nel)
       deallocate(gm%hex)
       facet_vtx = HEX_FACE
       edge_vtx = HEX_EDGE
       flip = HEX_FLIP
    else if (gm%nquad .gt. 0 .and. gm%n_volumes .gt. 0) then
       call fatal('The geometry has volumes but the file contains ' // &
            'no hexahedra. When physical groups are defined, Gmsh only ' // &
            'saves elements in them: add the volume to a physical group, ' // &
            'or set Mesh.SaveAll = 1')
    else if (gm%nquad .gt. 0) then
       gdim = 2
       nel = gm%nquad
       ncrn = 4
       nedge = 4
       nfacet = 4
       nfv = 2
       conn = gm%quad(:, 1:nel)
       deallocate(gm%quad, gm%quad_phys)
       gm%nquad = 0
       facet_vtx = QUAD_FACE
       edge_vtx = QUAD_EDGE
       flip = QUAD_FLIP
    else
       call fatal('No hexahedra or quadrilaterals found. When ' // &
            'physical groups are defined, Gmsh only saves elements in ' // &
            'them: add the volume (3D) or surface (2D) to a physical ' // &
            'group, or set Mesh.SaveAll = 1')
    end if

    nflip = 0
    nbad = 0
    c = 0.0_dp
    do e = 1, nel
       call check_degenerate(e)
       call corner_coords(e, c)
       if (center_jacobian(c) .lt. 0.0_dp) then
          conn(:, e) = conn(flip, e)
          nflip = nflip + 1
          call corner_coords(e, c)
       end if
       if (min_corner_jacobian(c) .le. 0.0_dp) nbad = nbad + 1
    end do

  end subroutine setup_elements

  !> Stop on elements with repeated vertices
  !! @param e Element index
  subroutine check_degenerate(e)
    integer, intent(in) :: e
    integer :: a, b
    character(len=16) :: str

    do a = 1, ncrn - 1
       do b = a + 1, ncrn
          if (conn(a, e) .eq. conn(b, e)) then
             write(str, '(I0)') e
             call fatal('Element ' // trim(str) // ' has a ' // &
                  'repeated vertex, collapsed elements are not supported')
          end if
       end do
    end do

  end subroutine check_degenerate

  !> Corner coordinates of element @a e in symmetric order
  !! @param e Element index
  !! @param c Corner coordinates, column k for symmetric vertex k
  subroutine corner_coords(e, c)
    integer, intent(in) :: e
    real(kind=dp), intent(inout) :: c(3, 8)
    integer :: k

    do k = 1, ncrn
       c(:, k) = gm%xyz(:, conn(k, e))
    end do

  end subroutine corner_coords

  !> Triple product a . (b x c)
  !! @param a First vector
  !! @param b Second vector
  !! @param c Third vector
  pure function triple(a, b, c) result(v)
    real(kind=dp), intent(in) :: a(3), b(3), c(3)
    real(kind=dp) :: v

    v = a(1) * (b(2) * c(3) - b(3) * c(2)) + &
         a(2) * (b(3) * c(1) - b(1) * c(3)) + &
         a(3) * (b(1) * c(2) - b(2) * c(1))

  end function triple

  !> Sign carrying Jacobian of the (bi/tri)linear map at the element centre
  !! @param c Corner coordinates in symmetric order
  function center_jacobian(c) result(det)
    real(kind=dp), intent(in) :: c(3, 8)
    real(kind=dp) :: det
    real(kind=dp) :: dr(3), ds(3), dt(3)

    if (gdim .eq. 3) then
       dr = (c(:,2) - c(:,1)) + (c(:,4) - c(:,3)) + &
            (c(:,6) - c(:,5)) + (c(:,8) - c(:,7))
       ds = (c(:,3) - c(:,1)) + (c(:,4) - c(:,2)) + &
            (c(:,7) - c(:,5)) + (c(:,8) - c(:,6))
       dt = (c(:,5) - c(:,1)) + (c(:,6) - c(:,2)) + &
            (c(:,7) - c(:,3)) + (c(:,8) - c(:,4))
       det = triple(dr, ds, dt)
    else
       dr = (c(:,2) - c(:,1)) + (c(:,4) - c(:,3))
       ds = (c(:,3) - c(:,1)) + (c(:,4) - c(:,2))
       det = dr(1) * ds(2) - dr(2) * ds(1)
    end if

  end function center_jacobian

  !> Smallest Jacobian of the (bi/tri)linear map over the element corners
  !! @param c Corner coordinates in symmetric order
  function min_corner_jacobian(c) result(det)
    real(kind=dp), intent(in) :: c(3, 8)
    real(kind=dp) :: det
    real(kind=dp) :: er(3), es(3), et(3)
    integer :: i, j, k

    det = huge(1.0_dp)
    if (gdim .eq. 3) then
       do k = 0, 1
          do j = 0, 1
             do i = 0, 1
                er = c(:, cidx(1, j, k)) - c(:, cidx(0, j, k))
                es = c(:, cidx(i, 1, k)) - c(:, cidx(i, 0, k))
                et = c(:, cidx(i, j, 1)) - c(:, cidx(i, j, 0))
                det = min(det, triple(er, es, et))
             end do
          end do
       end do
    else
       do j = 0, 1
          do i = 0, 1
             er = c(:, cidx(1, j, 0)) - c(:, cidx(0, j, 0))
             es = c(:, cidx(i, 1, 0)) - c(:, cidx(i, 0, 0))
             det = min(det, er(1) * es(2) - er(2) * es(1))
          end do
       end do
    end if

  end function min_corner_jacobian

  !> Symmetric vertex number of corner (i, j, k)
  !! @param i Position along r, 0 or 1
  !! @param j Position along s, 0 or 1
  !! @param k Position along t, 0 or 1
  pure function cidx(i, j, k) result(idx)
    integer, intent(in) :: i, j, k
    integer :: idx

    idx = 1 + i + 2 * j + 4 * k

  end function cidx

  !> Give the element corners compact vertex ids
  subroutine number_vertices()
    real(kind=dp) :: zmin, zmax, ext
    integer :: e, k, n

    allocate(vid(gm%nnodes))
    vid = 0
    nvert = 0
    do e = 1, nel
       do k = 1, ncrn
          n = conn(k, e)
          if (vid(n) .eq. 0) then
             nvert = nvert + 1
             vid(n) = nvert
          end if
       end do
    end do

    ! A 2D mesh is extruded in z by Neko, so it has to be planar in z
    if (gdim .eq. 2) then
       zmin = huge(1.0_dp)
       zmax = -huge(1.0_dp)
       ext = 0.0_dp
       do n = 1, gm%nnodes
          if (vid(n) .eq. 0) cycle
          zmin = min(zmin, gm%xyz(3, n))
          zmax = max(zmax, gm%xyz(3, n))
          ext = max(ext, maxval(abs(gm%xyz(1:2, n))))
       end do
       if (zmax - zmin .gt. 1.0e-6_dp * max(ext, tiny(1.0_dp))) then
          call fatal('No hexahedra found, and the quadrilaterals ' // &
               'do not lie in a plane z = constant as a 2D mesh must. ' // &
               'If this is a 3D mesh, add the volume to a physical ' // &
               'group or set Mesh.SaveAll = 1')
       end if
    end if

  end subroutine number_vertices

  !> Sort the first @a n entries of a facet key
  !! @param key Facet key
  !! @param n Number of entries to sort
  pure subroutine sort_key(key, n)
    integer, intent(inout) :: key(4)
    integer, intent(in) :: n
    integer :: i, j, t

    do i = 2, n
       t = key(i)
       j = i - 1
       do while (j .ge. 1)
          if (key(j) .le. t) exit
          key(j + 1) = key(j)
          j = j - 1
       end do
       key(j + 1) = t
    end do

  end subroutine sort_key

  !> Bucket of a key of four integers
  !! @param key Facet key
  !! @param nb Number of buckets
  pure function hash_key(key, nb) result(h)
    integer, intent(in) :: key(4), nb
    integer :: h
    integer(kind=i8) :: s

    s = modulo(int(key(1), i8), 1000003_i8) * 73856093_i8 + &
         modulo(int(key(2), i8), 1000003_i8) * 19349663_i8 + &
         modulo(int(key(3), i8), 1000003_i8) * 83492791_i8 + &
         modulo(int(key(4), i8), 1000003_i8) * 50331653_i8
    h = int(modulo(s, int(nb, i8))) + 1

  end function hash_key

  !> Bucket of a cell of a uniform grid
  !! @param cell Cell indices
  !! @param nb Number of buckets
  pure function hash_cell(cell, nb) result(h)
    integer(kind=i8), intent(in) :: cell(3)
    integer, intent(in) :: nb
    integer :: h
    integer(kind=i8) :: s

    s = modulo(cell(1), 1000003_i8) * 73856093_i8 + &
         modulo(cell(2), 1000003_i8) * 19349663_i8 + &
         modulo(cell(3), 1000003_i8) * 83492791_i8
    h = int(modulo(s, int(nb, i8))) + 1

  end function hash_cell

  !> Boundary facet with a given key, 0 if none
  !! @param key Facet key
  !! @param h Bucket of the key
  !! @param head First facet of each bucket
  !! @param nxt Next facet in the same bucket
  pure function find_key(key, h, head, nxt) result(j)
    integer, intent(in) :: key(4), h
    integer, intent(in) :: head(:), nxt(:)
    integer :: j

    j = head(h)
    do while (j .ne. 0)
       if (all(bf_key(:, j) .eq. key)) exit
       j = nxt(j)
    end do

  end function find_key

  !> Find the element facets of the boundary facets in physical groups
  subroutine match_boundary_facets()
    integer, allocatable :: head(:), nxt(:)
    integer :: key(4), corner(4)
    integer :: nraw, i, j, k, e, f, h, nb, phys
    logical :: unknown

    if (gdim .eq. 3) then
       nraw = gm%nquad
    else
       nraw = gm%nline
    end if

    allocate(bf_key(4, max(nraw, 1)), bf_phys(max(nraw, 1)))
    nbf = 0
    n_bf_unknown = 0
    corner = 0
    do i = 1, nraw
       if (gdim .eq. 3) then
          phys = gm%quad_phys(i)
          corner(1:4) = gm%quad(1:4, i)
       else
          phys = gm%line_phys(i)
          corner(1:2) = gm%line(1:2, i)
       end if
       if (phys .le. 0) cycle

       key = 0
       unknown = .false.
       do k = 1, nfv
          key(k) = vid(corner(k))
          if (key(k) .eq. 0) unknown = .true.
       end do
       if (unknown) then
          n_bf_unknown = n_bf_unknown + 1
          cycle
       end if
       call sort_key(key, nfv)
       nbf = nbf + 1
       bf_key(:, nbf) = key
       bf_phys(nbf) = phys
    end do

    ! Hash the facets, keeping the first of any duplicates
    nb = 2 * nbf + 1
    allocate(head(nb), nxt(max(nbf, 1)), bf_keep(nbf))
    head = 0
    nxt = 0
    n_bf_dup = 0
    n_bf_conflict = 0
    do i = 1, nbf
       h = hash_key(bf_key(:, i), nb)
       j = find_key(bf_key(:, i), h, head, nxt)
       if (j .ne. 0) then
          bf_keep(i) = .false.
          n_bf_dup = n_bf_dup + 1
          if (bf_phys(j) .ne. bf_phys(i)) n_bf_conflict = n_bf_conflict + 1
       else
          bf_keep(i) = .true.
          nxt(i) = head(h)
          head(h) = i
       end if
    end do

    ! Look up every element facet
    allocate(bf_nmatch(nbf), bf_el(2, nbf), bf_side(2, nbf))
    bf_nmatch = 0
    bf_el = 0
    bf_side = 0
    do e = 1, nel
       do f = 1, nfacet
          key = 0
          do k = 1, nfv
             key(k) = vid(conn(facet_vtx(k, f), e))
          end do
          call sort_key(key, nfv)
          j = find_key(key, hash_key(key, nb), head, nxt)
          if (j .eq. 0) cycle
          bf_nmatch(j) = bf_nmatch(j) + 1
          if (bf_nmatch(j) .gt. 2) then
             call fatal('A facet is shared by more than two ' // &
                  'elements, the mesh is not conforming')
          end if
          bf_el(bf_nmatch(j), j) = e
          bf_side(bf_nmatch(j), j) = f
       end do
    end do

    deallocate(head, nxt)

  end subroutine match_boundary_facets

  !> Index of a physical name of the boundary dimension, 0 if none
  !! @param tag Physical tag
  function find_phys_name(tag) result(idx)
    integer, intent(in) :: tag
    integer :: idx
    integer :: i

    idx = 0
    do i = 1, gm%nphys
       if (gm%phys_dim(i) .eq. gdim - 1 .and. gm%phys_tag(i) .eq. tag) then
          idx = i
          return
       end if
    end do

  end function find_phys_name

  !> Collect the boundary groups, resolve periodic pairs and assign labels
  subroutine setup_groups()
    character(len=GMSH_NAME_LEN), allocatable :: tok(:)
    integer, allocatable :: tmp(:)
    integer :: i, j, g, ga, gb, nlabeled, ntok, t
    logical :: fits

    ! Distinct physical tags of the matched facets
    allocate(tmp(max(nbf, 1)))
    ngrp = 0
    do i = 1, nbf
       if (.not. bf_keep(i) .or. bf_nmatch(i) .eq. 0) cycle
       if (any(tmp(1:ngrp) .eq. bf_phys(i))) cycle
       ngrp = ngrp + 1
       tmp(ngrp) = bf_phys(i)
    end do
    ! Insertion sort, the number of groups is small
    do i = 2, ngrp
       t = tmp(i)
       j = i - 1
       do while (j .ge. 1)
          if (tmp(j) .le. t) exit
          tmp(j + 1) = tmp(j)
          j = j - 1
       end do
       tmp(j + 1) = t
    end do

    allocate(grp_tag(ngrp), grp_label(ngrp), grp_partner(ngrp))
    allocate(grp_nfacets(ngrp), grp_ninternal(ngrp), grp_name(ngrp))
    grp_tag = tmp(1:ngrp)
    deallocate(tmp)
    grp_label = 0
    grp_partner = 0
    grp_nfacets = 0
    grp_ninternal = 0
    do g = 1, ngrp
       j = find_phys_name(grp_tag(g))
       grp_name(g) = ''
       if (j .gt. 0) grp_name(g) = gm%phys_name(j)
    end do

    allocate(bf_grp(nbf))
    bf_grp = 0
    do i = 1, nbf
       if (.not. bf_keep(i) .or. bf_nmatch(i) .eq. 0) cycle
       do g = 1, ngrp
          if (grp_tag(g) .eq. bf_phys(i)) exit
       end do
       bf_grp(i) = g
       grp_nfacets(g) = grp_nfacets(g) + 1
       if (bf_nmatch(i) .eq. 2) grp_ninternal(g) = grp_ninternal(g) + 1
    end do

    ! Periodic pairs
    call split_tokens(periodic_spec, tok, ntok)
    if (mod(ntok, 2) .ne. 0) then
       call fatal('--periodic expects pairs of physical groups')
    end if
    do i = 1, ntok, 2
       ga = resolve_group(tok(i))
       gb = resolve_group(tok(i + 1))
       if (ga .eq. gb) then
          call fatal('A physical group cannot be periodic with ' // &
               'itself: ' // trim(tok(i)))
       end if
       if (grp_partner(ga) .ne. 0 .or. grp_partner(gb) .ne. 0) then
          call fatal('A physical group is in more than one ' // &
               'periodic pair')
       end if
       grp_partner(ga) = gb
       grp_partner(gb) = ga
    end do
    deallocate(tok)

    ! Labels: the physical tags if they are valid zone indices, otherwise
    ! the groups are numbered 1, 2, ... in order of their tags
    nlabeled = count(grp_partner .eq. 0)
    if (nlabeled .gt. MAX_ZONES) then
       call fatal('More boundary physical groups than Neko ' // &
            'supports as labeled zones')
    end if
    fits = .true.
    do g = 1, ngrp
       if (grp_partner(g) .ne. 0) cycle
       if (grp_tag(g) .lt. 1 .or. grp_tag(g) .gt. MAX_ZONES) then
          fits = .false.
       end if
    end do
    renumbered = .not. fits
    j = 0
    do g = 1, ngrp
       if (grp_partner(g) .ne. 0) cycle
       j = j + 1
       if (fits) then
          grp_label(g) = grp_tag(g)
       else
          grp_label(g) = j
       end if
    end do

  end subroutine setup_groups

  !> Split a periodic specification into tokens
  !! @param spec Specification as given on the command line
  !! @param tok The tokens
  !! @param ntok Number of tokens
  subroutine split_tokens(spec, tok, ntok)
    character(len=*), intent(in) :: spec
    character(len=GMSH_NAME_LEN), allocatable, intent(out) :: tok(:)
    integer, intent(out) :: ntok
    character(len=len(spec)) :: str
    integer :: i, i0, pass

    str = spec
    do i = 1, len(str)
       if (index(',;:()' // achar(9), str(i:i)) .gt. 0) str(i:i) = ' '
    end do

    ! First pass counts, second pass stores
    do pass = 1, 2
       ntok = 0
       i = 1
       do while (i .le. len_trim(str))
          if (str(i:i) .eq. ' ') then
             i = i + 1
             cycle
          end if
          i0 = i
          do while (i .le. len(str))
             if (str(i:i) .eq. ' ') exit
             i = i + 1
          end do
          ntok = ntok + 1
          if (pass .eq. 2) tok(ntok) = str(i0:i-1)
       end do
       if (pass .eq. 1) allocate(tok(ntok))
    end do

  end subroutine split_tokens

  !> Group index of a physical group given by tag or name
  !! @param token Physical tag or name of the group
  function resolve_group(token) result(g)
    character(len=*), intent(in) :: token
    integer :: g
    integer :: i, tag, ios

    tag = -1
    if (verify(trim(token), '0123456789') .eq. 0) then
       read(token, *, iostat = ios) tag
    else
       do i = 1, gm%nphys
          if (gm%phys_dim(i) .eq. gdim - 1 .and. &
               trim(gm%phys_name(i)) .eq. trim(token)) then
             tag = gm%phys_tag(i)
             exit
          end if
       end do
    end if

    do g = 1, ngrp
       if (grp_tag(g) .eq. tag) return
    end do
    call fatal('Periodic group ' // trim(token) // ' is not a ' // &
         'boundary physical group of the mesh')

  end function resolve_group

  !> Collect the labeled zone records and make the periodic pairs
  subroutine mark_zones()
    integer :: i, m, g, n

    n = 2 * max(nbf, 1)
    allocate(lab_e(n), lab_f(n), lab_label(n))
    nlab = 0
    do i = 1, nbf
       if (bf_grp(i) .eq. 0) cycle
       g = bf_grp(i)
       if (grp_partner(g) .ne. 0) cycle
       do m = 1, bf_nmatch(i)
          nlab = nlab + 1
          lab_e(nlab) = bf_el(m, i)
          lab_f(nlab) = bf_side(m, i)
          lab_label(nlab) = grp_label(g)
       end do
    end do

    allocate(per_e(n), per_f(n), per_pe(n), per_pf(n))
    nper = 0
    allocate(parent(nvert))
    do i = 1, nvert
       parent(i) = i
    end do

    do g = 1, ngrp
       if (grp_partner(g) .gt. g) call make_periodic(g, grp_partner(g))
    end do

  end subroutine mark_zones

  !> Add a periodic facet record
  !! @param f Facet
  !! @param e Element
  !! @param pf Periodic facet
  !! @param pe Periodic element
  subroutine add_periodic(f, e, pf, pe)
    integer, intent(in) :: f, e, pf, pe

    nper = nper + 1
    per_f(nper) = f
    per_e(nper) = e
    per_pf(nper) = pf
    per_pe(nper) = pe

  end subroutine add_periodic
  !> Boundary facets of group @a g
  !! @param g Group index
  !! @param list Indices of the boundary facets of the group
  subroutine group_facets(g, list)
    integer, intent(in) :: g
    integer, allocatable, intent(out) :: list(:)
    integer :: i, n

    allocate(list(grp_nfacets(g)))
    n = 0
    do i = 1, nbf
       if (bf_grp(i) .ne. g) cycle
       n = n + 1
       list(n) = i
    end do

  end subroutine group_facets

  !> Corners, centroid and vertex ids of boundary facet @a b
  !! @param b Boundary facet index
  !! @param x Corner coordinates
  !! @param c Centroid
  !! @param v Vertex ids of the corners
  !! @param emin Shortest facet edge so far, updated
  !! @param scale Largest coordinate magnitude so far, updated
  subroutine facet_geometry(b, x, c, v, emin, scale)
    integer, intent(in) :: b
    real(kind=dp), intent(out) :: x(3, nfv), c(3)
    integer, intent(out) :: v(nfv)
    real(kind=dp), intent(inout) :: emin, scale
    integer :: k, node

    do k = 1, nfv
       node = conn(facet_vtx(k, bf_side(1, b)), bf_el(1, b))
       x(:, k) = gm%xyz(:, node)
       if (gdim .eq. 2) x(3, k) = 0.0_dp
       v(k) = vid(node)
       scale = max(scale, maxval(abs(x(:, k))))
    end do
    c = sum(x, dim = 2) / real(nfv, dp)

    if (nfv .eq. 2) then
       emin = min(emin, norm2(x(:, 2) - x(:, 1)))
    else
       do k = 1, nfv
          emin = min(emin, norm2(x(:, mod(k, nfv) + 1) - x(:, k)))
       end do
    end if

  end subroutine facet_geometry

  !> Match the corners of two facets under a translation
  !! @param xa Corners of the first facet
  !! @param xb Corners of the second facet
  !! @param offset Translation from the first facet to the second
  !! @param tol Matching tolerance
  !! @param pair Corner of @a xb matching each corner of @a xa
  function match_corners(xa, xb, offset, tol, pair) result(match)
    real(kind=dp), intent(in) :: xa(3, nfv), xb(3, nfv), offset(3), tol
    integer, intent(out) :: pair(4)
    logical :: match
    integer :: k, m, n

    match = .false.
    pair = 0
    do k = 1, nfv
       n = 0
       do m = 1, nfv
          if (norm2(xb(:, m) - xa(:, k) - offset) .le. tol) then
             n = n + 1
             pair(k) = m
          end if
       end do
       if (n .ne. 1) return
    end do
    match = .true.

  end function match_corners

  !> Make the groups @a ga and @a gb periodic
  !! @details The two groups must be related by a translation, which is
  !! taken as the difference of their mean facet centroids. Each facet of
  !! @a ga is matched to one facet of @a gb with a hash of the centroids,
  !! then both directions are marked as periodic facets and the matched
  !! vertices are merged.
  !! @param ga Group index of the first side
  !! @param gb Group index of the second side
  subroutine make_periodic(ga, gb)
    integer, intent(in) :: ga, gb
    integer, allocatable :: fa(:), fb(:), va(:,:), vb(:,:)
    integer, allocatable :: head(:), nxt(:)
    integer(kind=i8), allocatable :: cell(:,:)
    real(kind=dp), allocatable :: xa(:,:,:), xb(:,:,:), ca(:,:), cb(:,:)
    logical, allocatable :: used(:)
    real(kind=dp) :: offset(3), t(3), emin, scale, tol, h
    integer(kind=i8) :: c0(3), cc(3)
    integer :: pair(4), pair_m(4)
    integer :: na, nb, i, j, k, jm, nfound, dx, dy, dz
    integer :: ea, fa_side, eb, fb_side
    character(len=512) :: msg

    if (grp_ninternal(ga) .gt. 0 .or. grp_ninternal(gb) .gt. 0) then
       call fatal('Periodic group ' // trim(group_label_str(ga)) // &
            ' or ' // trim(group_label_str(gb)) // ' contains interior ' // &
            'facets')
    end if
    call group_facets(ga, fa)
    call group_facets(gb, fb)
    na = size(fa)
    if (size(fb) .ne. na) then
       write(msg, '(A,I0,A,I0,A)') 'Periodic groups ' // &
            trim(group_label_str(ga)) // ' and ' // &
            trim(group_label_str(gb)) // ' have ', na, ' and ', size(fb), &
            ' facets, they must match one to one'
       call fatal(trim(msg))
    end if

    allocate(xa(3, nfv, na), xb(3, nfv, na), ca(3, na), cb(3, na))
    allocate(va(nfv, na), vb(nfv, na))
    emin = huge(1.0_dp)
    scale = 0.0_dp
    do i = 1, na
       call facet_geometry(fa(i), xa(:,:,i), ca(:,i), va(:,i), emin, scale)
       call facet_geometry(fb(i), xb(:,:,i), cb(:,i), vb(:,i), emin, scale)
    end do

    if (emin .le. 0.0_dp) then
       call fatal('Periodic group ' // trim(group_label_str(ga)) // &
            ' or ' // trim(group_label_str(gb)) // ' has a facet with ' // &
            'coincident vertices')
    end if

    offset = (sum(cb, dim = 2) - sum(ca, dim = 2)) / real(na, dp)
    if (user_tol_set) then
       tol = user_tol
    else
       tol = max(1.0e-6_dp * emin, 1.0e-12_dp * scale)
    end if
    h = max(0.5_dp * emin, 2.0_dp * tol)

    ! Bucket the facets of gb by the grid cell of their centroid
    nb = 2 * na + 1
    allocate(head(nb), nxt(na), cell(3, na), used(na))
    head = 0
    used = .false.
    do j = 1, na
       cell(:, j) = floor(cb(:, j) / h, kind=i8)
       k = hash_cell(cell(:, j), nb)
       nxt(j) = head(k)
       head(k) = j
    end do

    do i = 1, na
       t = ca(:, i) + offset
       c0 = floor(t / h, kind=i8)
       nfound = 0
       jm = 0
       pair_m = 0
       do dz = -1, 1
          do dy = -1, 1
             do dx = -1, 1
                cc = c0 + int([dx, dy, dz], i8)
                j = head(hash_cell(cc, nb))
                do while (j .ne. 0)
                   if (all(cell(:, j) .eq. cc)) then
                      if (norm2(cb(:, j) - t) .le. tol) then
                         if (match_corners(xa(:,:,i), xb(:,:,j), offset, &
                              tol, pair)) then
                            nfound = nfound + 1
                            jm = j
                            pair_m = pair
                         end if
                      end if
                   end if
                   j = nxt(j)
                end do
             end do
          end do
       end do

       if (nfound .ne. 1 .or. used(max(jm, 1))) then
          write(msg, '(A,3(1X,ES12.5),A,ES9.2,A)') 'No unique periodic ' // &
               'partner in group ' // trim(group_label_str(gb)) // &
               ' for the facet of group ' // trim(group_label_str(ga)) // &
               ' at', ca(:, i), ' (tolerance', tol, '). The groups ' // &
               'must be related by a translation, see --tol'
          call fatal(trim(msg))
       end if
       used(jm) = .true.

       ea = bf_el(1, fa(i))
       fa_side = bf_side(1, fa(i))
       eb = bf_el(1, fb(jm))
       fb_side = bf_side(1, fb(jm))
       call add_periodic(fa_side, ea, fb_side, eb)
       call add_periodic(fb_side, eb, fa_side, ea)

       do k = 1, nfv
          call uf_union(va(k, i), vb(pair_m(k), jm))
       end do
    end do

    call progress('Periodic ' // trim(group_label_str(ga)) // ' -> ' // &
         trim(group_label_str(gb)), 'offset ' // real_str(offset(1)) // &
         ' ' // real_str(offset(2)) // ' ' // real_str(offset(3)) // &
         ', tolerance ' // real_str(tol))

    deallocate(fa, fb, va, vb, xa, xb, ca, cb, head, nxt, cell, used)

  end subroutine make_periodic

  !> Root of the set containing @a x0, with path compression
  !! @param x0 Vertex id
  function uf_find(x0) result(r)
    integer, intent(in) :: x0
    integer :: r
    integer :: x, nx

    r = x0
    do while (parent(r) .ne. r)
       r = parent(r)
    end do
    x = x0
    do while (parent(x) .ne. r)
       nx = parent(x)
       parent(x) = r
       x = nx
    end do

  end function uf_find

  !> Merge the sets containing @a a and @a b, keeping the smaller root
  !! @param a Vertex id in the first set
  !! @param b Vertex id in the second set
  subroutine uf_union(a, b)
    integer, intent(in) :: a, b
    integer :: ra, rb

    ra = uf_find(a)
    rb = uf_find(b)
    if (ra .lt. rb) then
       parent(rb) = ra
    else if (rb .lt. ra) then
       parent(ra) = rb
    end if

  end subroutine uf_union

  !> Midpoint curves of the curved edges of element @a e
  !! @return .true. if the element has a curved edge
  !! @param e Element index
  !! @param curve_data Curve data of each edge, the midpoint in entries 1 to 3
  !! @param curve_type Curve type of each edge, 4 for a midpoint curve, else 0
  function element_curve(e, curve_data, curve_type) result(curved)
    integer, intent(in) :: e
    real(kind=dp), intent(out) :: curve_data(5, 12)
    integer, intent(out) :: curve_type(12)
    logical :: curved
    real(kind=dp) :: a(3), b(3), xm(3)
    integer :: k, m

    curve_data = 0.0_dp
    curve_type = 0
    curved = .false.
    if (linear) return

    do k = 1, nedge
       m = conn(ncrn + k, e)
       if (m .eq. 0) cycle
       a = gm%xyz(:, conn(edge_vtx(1, k), e))
       b = gm%xyz(:, conn(edge_vtx(2, k), e))
       xm = gm%xyz(:, m)
       if (gdim .eq. 2) then
          a(3) = 0.0_dp
          b(3) = 0.0_dp
          xm(3) = 0.0_dp
       end if
       if (norm2(xm - 0.5_dp * (a + b)) .gt. curve_tol * norm2(b - a)) then
          ! Neko curve type 4: quadratic edge through its midpoint
          curve_type(k) = 4
          curve_data(1:3, k) = xm
          curved = .true.
       end if
    end do

  end function element_curve

  !> Count the curved elements and edges
  subroutine count_curves()
    real(kind=dp) :: curve_data(5, 12)
    integer :: curve_type(12)
    integer :: e

    ncurved_el = 0
    ncurved_edges = 0
    do e = 1, nel
       if (element_curve(e, curve_data, curve_type)) then
          ncurved_el = ncurved_el + 1
          ncurved_edges = ncurved_edges + count(curve_type .ne. 0)
       end if
    end do

  end subroutine count_curves

  !> Write the mesh in the nmsh format
  !! @details The layout is that of Neko's nmsh_file_write: the number of
  !! elements and the dimension, the elements with their vertices in cyclic
  !! order, the zones (periodic facets first, then the labeled zones in
  !! order of their index) and the curved elements. Integers are 4 bytes,
  !! reals 8 bytes, without padding. Elements keep the original vertex ids,
  !! periodic facets store the merged ids of their vertices.
  !! @param fname Name of the nmsh file
  subroutine write_nmsh(fname)
    character(len=*), intent(in) :: fname
    real(kind=dp) :: x(3, 8), curve_data(5, 12)
    integer :: ids(8), pids(4), curve_type(12)
    integer :: unit, ierr, nv, e, i, j, k, node
    integer(kind=i8) :: clock_end, clock_rate
    character(len=16) :: time_str

    open(newunit = unit, file = fname, access = 'stream', &
         form = 'unformatted', status = 'replace', action = 'write', &
         iostat = ierr)
    if (ierr .ne. 0) call fatal('Cannot open ' // fname // ' for writing')

    nv = 2**gdim
    write(unit) int(nel, i4), int(gdim, i4)
    do e = 1, nel
       do j = 1, nv
          node = conn(CYC_TO_SYM(j), e)
          ids(j) = vid(node)
          x(:, j) = gm%xyz(:, node)
          if (gdim .eq. 2) x(3, j) = 0.0_dp
       end do
       write(unit) int(e, i4), (int(ids(j), i4), x(:, j), j = 1, nv)
    end do

    write(unit) int(nper + nlab, i4)
    do i = 1, nper
       pids = 0
       do k = 1, nfv
          pids(k) = uf_find(vid(conn(facet_vtx(k, per_f(i)), per_e(i))))
       end do
       write(unit) int([per_e(i), per_f(i), per_pe(i), per_pf(i), pids, &
            PERIODIC_ZONE], i4)
    end do
    do k = 1, MAX_ZONES
       do i = 1, nlab
          if (lab_label(i) .ne. k) cycle
          write(unit) int([lab_e(i), lab_f(i), 0, k, 0, 0, 0, 0, &
               LABELED_ZONE], i4)
       end do
    end do

    write(unit) int(ncurved_el, i4)
    do e = 1, nel
       if (element_curve(e, curve_data, curve_type)) then
          write(unit) int(e, i4), curve_data, int(curve_type, i4)
       end if
    end do

    close(unit)

    call system_clock(clock_end, clock_rate)
    write(time_str, '(F12.1)') real(clock_end - clock_start, dp) / &
         real(clock_rate, dp)
    call progress('Writing ' // fname, 'done, ' // &
         trim(adjustl(time_str)) // ' s in total')

  end subroutine write_nmsh
  !> Name of a group for messages
  !! @param g Group index
  function group_label_str(g) result(str)
    integer, intent(in) :: g
    character(len=GMSH_NAME_LEN + 16) :: str

    if (len_trim(grp_name(g)) .gt. 0) then
       write(str, '(A,A,A)') "'", trim(grp_name(g)), "'"
    else
       write(str, '(I0)') grp_tag(g)
    end if

  end function group_label_str

  !> Print the elements and vertices
  subroutine print_elements()
    character(len=:), allocatable :: str
    integer :: nho

    if (gdim .eq. 3) then
       nho = gm%nhex_ho
       str = int_str(nel) // ' hexahedra'
    else
       nho = gm%nquad_ho
       str = int_str(nel) // ' quadrilaterals, 2D'
    end if
    if (nho .eq. nel) then
       str = str // ', all second order'
    else if (nho .gt. 0) then
       str = str // ', ' // int_str(nho) // ' second order'
    end if
    call progress('  Elements', str)
    call progress('  Vertices', int_str(nvert))
    if (nflip .gt. 0) then
       call progress('  Reoriented', int_str(nflip) // &
            ' left-handed elements')
    end if
    deallocate(str)

  end subroutine print_elements

  !> Print the number of matched boundary facets
  subroutine print_groups()

    if (ngrp .eq. 0) then
       call progress('Matching boundary facets', 'no physical groups')
    else if (ngrp .eq. 1) then
       call progress('Matching boundary facets', &
            int_str(sum(grp_nfacets)) // ' facets in 1 physical group')
    else
       call progress('Matching boundary facets', &
            int_str(sum(grp_nfacets)) // ' facets in ' // int_str(ngrp) // &
            ' physical groups')
    end if

  end subroutine print_groups

  !> Print the number of curved elements
  subroutine print_curves()

    if (linear) then
       call progress('Curved elements', 'none, --linear')
    else if (ncurved_el .eq. 0) then
       call progress('Curved elements', 'none')
    else
       call progress('Curved elements', int_str(ncurved_el) // &
            ' elements, ' // int_str(ncurved_edges) // ' edges')
    end if

  end subroutine print_curves

  !> Print the boundary groups with their zones
  subroutine print_table()
    character(len=32) :: zone
    integer :: g

    if (ngrp .eq. 0) return

    write(*, '(A)') ''
    write(*, '(A)') '   Tag     Facets  Zone        Physical group'
    do g = 1, ngrp
       if (grp_partner(g) .ne. 0) then
          write(zone, '(A,I0)') 'periodic ', grp_tag(grp_partner(g))
       else
          write(zone, '(I0)') grp_label(g)
       end if
       write(*, '(I6,1X,I10,2X,A,A)') grp_tag(g), grp_nfacets(g), &
            zone(1:12), trim(grp_name(g))
    end do
    write(*, '(A)') ''

  end subroutine print_table

  !> Print warnings about the conversion
  subroutine print_warnings()
    integer :: nunmatched, ninternal

    if (ngrp .eq. 0) then
       write(*, '(A)') ''
       call warn('No boundary physical groups found, the mesh ' // &
            'has no labeled zones')
    end if
    if (renumbered) then
       call warn('Physical tags outside [1, 20], the labeled ' // &
            'zones are numbered in order of the tags (see the table)')
    end if
    if (nflip .gt. 0) then
       call warn('Left-handed elements were mirrored')
    end if
    if (nbad .gt. 0) then
       call warn('Some elements have a non-positive Jacobian ' // &
            'at a vertex, check the mesh quality')
    end if
    if (gm%n_multi_phys .gt. 0 .or. n_bf_conflict .gt. 0) then
       call warn('Some boundary facets are in more than one ' // &
            'physical group, the first group is used')
    end if
    nunmatched = n_bf_unknown
    ninternal = 0
    if (nbf .gt. 0) then
       nunmatched = nunmatched + count(bf_keep .and. bf_nmatch .eq. 0)
       ninternal = count(bf_keep .and. bf_nmatch .eq. 2)
    end if
    if (nunmatched .gt. 0) then
       call warn('Some facets in physical groups are not ' // &
            'facets of any element and were ignored')
    end if
    if (ninternal .gt. 0) then
       call warn('Some facets in physical groups are interior ' // &
            'facets, they are labeled on both sides')
    end if

  end subroutine print_warnings
end program gmsh2nmsh
