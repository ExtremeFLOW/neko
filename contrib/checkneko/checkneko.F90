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

!> Kinds, text helpers, JSON access and messages for checkneko
!! @details checkneko does not depend on the Neko library, only on
!! json-fortran, so that it can be built with any Fortran 2008 compiler, also
!! on its own for a machine Neko is not built for.
module checkneko_util
  use, intrinsic :: iso_fortran_env, only : int32, int64, real64, &
       error_unit, output_unit
  use json_module, only : json_core, json_value, json_IK, json_RK, json_LK, &
       json_null, json_object, json_array, json_logical, json_integer, &
       json_real, json_string
  implicit none
  private

  integer, public, parameter :: i4 = int32
  integer, public, parameter :: i8 = int64
  integer, public, parameter :: dp = real64

  !> Severity of a message
  integer, public, parameter :: SEV_ERROR = 1, SEV_WARNING = 2

  !> What a message from the schema validation is about. The kinds are used
  !! to pick the alternative of a oneOf or anyOf the user most likely meant.
  integer, public, parameter :: MK_OTHER = 0, MK_TYPE = 1, MK_CONST = 2, &
       MK_ENUM = 3, MK_REQUIRED = 4, MK_UNKNOWN_KEY = 5, MK_NOT = 6

  !> A string of any length, for arrays of strings
  type, public :: str_t
     character(len=:), allocatable :: s
  end type str_t

  !> A message about the case file
  type, public :: msg_t
     integer :: sev = SEV_ERROR
     integer :: kind = MK_OTHER
     !> Where in the case file, e.g. case.fluid.boundary_conditions[2]
     character(len=:), allocatable :: path
     character(len=:), allocatable :: text
     !> The value the message is about, which gives the line number
     type(json_value), pointer :: node => null()
     !> The value as text
     character(len=:), allocatable :: got
     !> The accepted values (MK_CONST, MK_ENUM) or types (MK_TYPE)
     type(str_t), allocatable :: allowed(:)
  end type msg_t

  !> A list of messages
  type, public :: msg_list_t
     type(msg_t), allocatable :: m(:)
     integer :: n = 0
   contains
     procedure, pass(this) :: add => msg_list_add
     procedure, pass(this) :: truncate => msg_list_truncate
  end type msg_list_t

  !> The JSON core used to access the case file and the schemas
  type(json_core), public :: jcore

  !> The messages about the case file being checked
  type(msg_list_t), public :: msgs

  public :: fatal, exit_with, warn, progress, report, int_str, real_text, &
       quoted, &
       describe, render, type_phrase, value_list, edit_distance, &
       typo_limit, is_whole, str_eq, file_exists, join_path, file_suffix, &
       jtype, jnchildren, jname, jchild, jelement, jfirst, jnext, jstr, &
       jint, jreal, jbool, is_number, node_at, path_key, path_index, &
       parent_path, add_str

contains

  !> Print an error about the tool itself, not about the case, and stop
  !! @param msg Message to print before stopping
  subroutine fatal(msg)
    character(len=*), intent(in) :: msg

    flush(output_unit)
    write(error_unit, '(A)') 'checkneko: ' // trim(msg)
    call exit_with(1)

  end subroutine fatal

  !> Stop with an exit status
  !! @details Unlike STOP, this does not print the status, which some
  !! compilers do. The status is how a script sees the result.
  !! @param status The exit status
  subroutine exit_with(status)
    use, intrinsic :: iso_c_binding, only : c_int
    integer, intent(in) :: status
    interface
       subroutine c_exit(code) bind(c, name = 'exit')
         import :: c_int
         integer(kind=c_int), value :: code
       end subroutine c_exit
    end interface

    flush(output_unit)
    flush(error_unit)
    call c_exit(int(status, c_int))

  end subroutine exit_with

  !> Print a warning about the tool itself, not about the case
  !! @param msg Message to print
  subroutine warn(msg)
    character(len=*), intent(in) :: msg

    flush(output_unit)
    write(error_unit, '(A)') 'checkneko: warning: ' // trim(msg)

  end subroutine warn

  !> Print a progress line: the label, a row of dots and the value
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

  !> Add a message about the case file
  !! @param sev Severity, SEV_ERROR or SEV_WARNING
  !! @param path Where in the case file
  !! @param text The message
  !! @param node The value the message is about, may be null
  !! @param kind What the message is about, see MK_OTHER
  !! @param got The value as text
  !! @param allowed The accepted values or types
  subroutine report(sev, path, text, node, kind, got, allowed)
    integer, intent(in) :: sev
    character(len=*), intent(in) :: path, text
    type(json_value), pointer, intent(in) :: node
    integer, intent(in), optional :: kind
    character(len=*), intent(in), optional :: got
    type(str_t), intent(in), optional :: allowed(:)
    type(msg_t) :: msg

    msg%sev = sev
    msg%path = path
    msg%text = text
    msg%node => node
    if (present(kind)) msg%kind = kind
    if (present(got)) msg%got = got
    if (present(allowed)) msg%allowed = allowed
    call msgs%add(msg)

  end subroutine report

  !> Append a message to a list
  !! @param msg The message
  subroutine msg_list_add(this, msg)
    class(msg_list_t), intent(inout) :: this
    type(msg_t), intent(in) :: msg
    type(msg_t), allocatable :: tmp(:)

    if (.not. allocated(this%m)) allocate(this%m(64))
    if (this%n .eq. size(this%m)) then
       allocate(tmp(2 * size(this%m)))
       tmp(1:this%n) = this%m(1:this%n)
       call move_alloc(tmp, this%m)
    end if
    this%n = this%n + 1
    this%m(this%n) = msg

  end subroutine msg_list_add

  !> Drop the messages after the first n
  !! @param n Number of messages to keep
  subroutine msg_list_truncate(this, n)
    class(msg_list_t), intent(inout) :: this
    integer, intent(in) :: n

    this%n = min(this%n, n)

  end subroutine msg_list_truncate

  !> Append a string to an array of strings, unless it is there already
  !! @param list The array
  !! @param s The string
  subroutine add_str(list, s)
    type(str_t), allocatable, intent(inout) :: list(:)
    character(len=*), intent(in) :: s
    type(str_t), allocatable :: tmp(:)
    integer :: i, n

    n = 0
    if (allocated(list)) then
       n = size(list)
       do i = 1, n
          if (str_eq(list(i)%s, s)) return
       end do
    end if
    allocate(tmp(n + 1))
    if (n .gt. 0) tmp(1:n) = list
    tmp(n + 1)%s = s
    call move_alloc(tmp, list)

  end subroutine add_str

  !> An integer as text
  !! @param i The integer
  function int_str(i) result(str)
    integer, intent(in) :: i
    character(len=:), allocatable :: str
    character(len=24) :: buf

    write(buf, '(I0)') i
    str = trim(buf)

  end function int_str

  !> A real as text that reads as a real, e.g. 7.0 rather than 7
  !! @param x The real
  function real_text(x) result(str)
    real(kind=dp), intent(in) :: x
    character(len=:), allocatable :: str
    character(len=40) :: buf
    integer :: i, e

    if (abs(x) .lt. tiny(x)) then
       str = '0.0'
       return
    end if

    if (abs(x) .ge. 1.0e-3_dp .and. abs(x) .lt. 1.0e12_dp) then
       write(buf, '(F32.10)') x
       buf = adjustl(buf)
       i = len_trim(buf)
       do while (buf(i:i) .eq. '0' .and. buf(i-1:i-1) .ne. '.')
          buf(i:i) = ' '
          i = i - 1
       end do
       str = trim(buf)
    else
       write(buf, '(ES24.15E3)') x
       buf = adjustl(buf)
       e = index(buf, 'E')
       i = e - 1
       do while (buf(i:i) .eq. '0' .and. buf(i-1:i-1) .ne. '.')
          i = i - 1
       end do
       str = buf(1:i) // trim(buf(e:))
    end if

  end function real_text

  !> A string in double quotes, shortened if it is long
  !! @param s The string
  function quoted(s) result(str)
    character(len=*), intent(in) :: s
    character(len=:), allocatable :: str
    integer, parameter :: MAX_LEN = 60

    if (len(s) .gt. MAX_LEN) then
       str = '"' // s(1:MAX_LEN - 3) // '..."'
    else
       str = '"' // s // '"'
    end if

  end function quoted

  !> A JSON value as it would be written in the case file
  !! @param p The value
  function render(p) result(str)
    type(json_value), pointer, intent(in) :: p
    character(len=:), allocatable :: str

    select case (jtype(p))
    case (json_string)
       str = quoted(jstr(p))
    case (json_integer)
       str = int_str(jint(p))
    case (json_real)
       str = real_text(jreal(p))
    case (json_logical)
       if (jbool(p)) then
          str = 'true'
       else
          str = 'false'
       end if
    case (json_null)
       str = 'null'
    case (json_object)
       str = '{...}'
    case (json_array)
       str = '[...]'
    case default
       str = '?'
    end select

  end function render

  !> A JSON value described for a message, e.g. "the integer 3"
  !! @param p The value
  function describe(p) result(str)
    type(json_value), pointer, intent(in) :: p
    character(len=:), allocatable :: str

    select case (jtype(p))
    case (json_string)
       str = 'the string ' // render(p)
    case (json_integer)
       str = 'the integer ' // render(p)
    case (json_real)
       str = 'the real number ' // render(p)
    case (json_logical)
       str = render(p)
    case (json_null)
       str = 'null'
    case (json_object)
       str = 'an object'
    case (json_array)
       str = 'an array'
    case default
       str = 'a value'
    end select

  end function describe

  !> Schema type names as a phrase, e.g. "an integer or a string"
  !! @param names The type names, "constant" for the name of a constant
  function type_phrase(names) result(str)
    type(str_t), intent(in) :: names(:)
    character(len=:), allocatable :: str, word
    integer :: i

    str = ''
    do i = 1, size(names)
       select case (names(i)%s)
       case ('object')
          word = 'an object'
       case ('array')
          word = 'an array'
       case ('string')
          word = 'a string'
       case ('integer')
          word = 'an integer'
       case ('number')
          word = 'a number'
       case ('boolean')
          word = 'true or false'
       case ('null')
          word = 'null'
       case ('constant')
          word = 'the name of a constant'
       case default
          word = names(i)%s
       end select
       if (i .eq. 1) then
          str = word
       else if (i .eq. size(names)) then
          str = str // ' or ' // word
       else
          str = str // ', ' // word
       end if
    end do

  end function type_phrase

  !> Values as a comma separated list
  !! @param values The values, already rendered
  function value_list(values) result(str)
    type(str_t), intent(in) :: values(:)
    character(len=:), allocatable :: str
    integer :: i

    str = ''
    do i = 1, size(values)
       if (i .gt. 1) str = str // ', '
       str = str // values(i)%s
    end do

  end function value_list

  !> Edit distance between two strings, ignoring case
  !! @details The optimal string alignment distance: insertions, deletions,
  !! substitutions and transpositions of neighbours all cost one.
  !! @param a First string
  !! @param b Second string
  function edit_distance(a, b) result(dist)
    character(len=*), intent(in) :: a, b
    integer :: dist
    integer, allocatable :: d(:,:)
    character(len=len(a)) :: la
    character(len=len(b)) :: lb
    integer :: i, j, m, n, cost

    la = lower(a)
    lb = lower(b)
    m = len(a)
    n = len(b)
    allocate(d(0:m, 0:n))
    do i = 0, m
       d(i, 0) = i
    end do
    do j = 0, n
       d(0, j) = j
    end do
    do j = 1, n
       do i = 1, m
          cost = 1
          if (la(i:i) .eq. lb(j:j)) cost = 0
          d(i, j) = min(d(i - 1, j) + 1, d(i, j - 1) + 1, &
               d(i - 1, j - 1) + cost)
          if (i .gt. 1 .and. j .gt. 1) then
             if (la(i:i) .eq. lb(j - 1:j - 1) .and. &
                  la(i - 1:i - 1) .eq. lb(j:j)) then
                d(i, j) = min(d(i, j), d(i - 2, j - 2) + 1)
             end if
          end if
       end do
    end do
    dist = d(m, n)

  end function edit_distance

  !> Whether a real is a whole number
  !! @param x The real
  elemental function is_whole(x) result(whole)
    real(kind=dp), intent(in) :: x
    logical :: whole

    whole = abs(x - aint(x)) .le. 0.0_dp

  end function is_whole

  !> Largest edit distance at which a name counts as a misspelling
  !! @param n Length of the name
  pure function typo_limit(n) result(limit)
    integer, intent(in) :: n
    integer :: limit

    if (n .le. 4) then
       limit = 1
    else if (n .le. 9) then
       limit = 2
    else
       limit = 3
    end if

  end function typo_limit

  !> A string in lower case
  !! @param s The string
  pure function lower(s) result(l)
    character(len=*), intent(in) :: s
    character(len=len(s)) :: l
    integer :: i, c

    l = s
    do i = 1, len(s)
       c = iachar(s(i:i))
       if (c .ge. iachar('A') .and. c .le. iachar('Z')) then
          l(i:i) = achar(c + 32)
       end if
    end do

  end function lower

  !> Exact string equality, unlike == which ignores trailing blanks
  !! @param a First string
  !! @param b Second string
  pure function str_eq(a, b) result(eq)
    character(len=*), intent(in) :: a, b
    logical :: eq

    eq = len(a) .eq. len(b)
    if (eq) eq = a .eq. b

  end function str_eq

  !> Whether a file exists
  !! @param fname File name
  function file_exists(fname) result(exists)
    character(len=*), intent(in) :: fname
    logical :: exists

    inquire(file = fname, exist = exists)

  end function file_exists

  !> A file name relative to a directory
  !! @param dir The directory, "." for the current directory
  !! @param fname The file name, returned as it is if absolute
  function join_path(dir, fname) result(path)
    character(len=*), intent(in) :: dir, fname
    character(len=:), allocatable :: path

    if (len(fname) .gt. 0) then
       if (fname(1:1) .eq. '/') then
          path = fname
          return
       end if
    end if
    if (str_eq(dir, '.') .or. len(dir) .eq. 0) then
       path = fname
    else if (dir(len(dir):len(dir)) .eq. '/') then
       path = dir // fname
    else
       path = dir // '/' // fname
    end if

  end function join_path

  !> Text after the last dot of a file name, empty if there is none
  !! @param fname File name
  function file_suffix(fname) result(suffix)
    character(len=*), intent(in) :: fname
    character(len=:), allocatable :: suffix
    integer :: i

    i = scan(fname, '.', back = .true.)
    if (i .eq. 0 .or. i .lt. scan(fname, '/', back = .true.)) then
       suffix = ''
    else
       suffix = fname(i + 1:)
    end if

  end function file_suffix

  !> Type of a JSON value, one of json_object, json_integer, ...
  !! @param p The value
  function jtype(p) result(t)
    type(json_value), pointer, intent(in) :: p
    integer :: t
    integer(kind=json_IK) :: vt

    t = 0
    if (.not. associated(p)) return
    call jcore%info(p, var_type = vt)
    t = int(vt)

  end function jtype

  !> Number of entries of an array or keys of an object
  !! @param p The value
  function jnchildren(p) result(n)
    type(json_value), pointer, intent(in) :: p
    integer :: n
    integer(kind=json_IK) :: nc

    n = 0
    if (.not. associated(p)) return
    call jcore%info(p, n_children = nc)
    n = int(nc)

  end function jnchildren

  !> Key of a value in an object
  !! @param p The value
  function jname(p) result(name)
    type(json_value), pointer, intent(in) :: p
    character(len=:), allocatable :: name

    call jcore%info(p, name = name)
    if (.not. allocated(name)) name = ''

  end function jname

  !> Value of a key in an object, null if the key is not there
  !! @details Gives the first one if the key appears more than once, as Neko
  !! does.
  !! @param p The object
  !! @param name The key
  function jchild(p, name) result(c)
    type(json_value), pointer, intent(in) :: p
    character(len=*), intent(in) :: name
    type(json_value), pointer :: c
    logical :: found

    c => null()
    if (jtype(p) .ne. json_object) return
    call jcore%get_child(p, name, c, found)
    if (.not. found) c => null()

  end function jchild

  !> Entry of an array, null if there is no such entry
  !! @param p The array
  !! @param i Index of the entry, starting at 1
  function jelement(p, i) result(c)
    type(json_value), pointer, intent(in) :: p
    integer, intent(in) :: i
    type(json_value), pointer :: c
    logical :: found

    c => null()
    if (jtype(p) .ne. json_array .or. i .lt. 1) return
    call jcore%get_child(p, int(i, json_IK), c, found)
    if (.not. found) c => null()

  end function jelement

  !> First entry of an array or object, null if it is empty
  !! @param p The array or object
  function jfirst(p) result(c)
    type(json_value), pointer, intent(in) :: p
    type(json_value), pointer :: c

    c => null()
    if (jtype(p) .ne. json_object .and. jtype(p) .ne. json_array) return
    call jcore%get_child(p, c)

  end function jfirst

  !> Next entry of the array or object that holds a value
  !! @param p The value
  function jnext(p) result(c)
    type(json_value), pointer, intent(in) :: p
    type(json_value), pointer :: c

    c => null()
    if (.not. associated(p)) return
    call jcore%get_next(p, c)

  end function jnext

  !> Value of a string, empty for other types
  !! @param p The value
  function jstr(p) result(s)
    type(json_value), pointer, intent(in) :: p
    character(len=:), allocatable :: s

    if (jtype(p) .eq. json_string) then
       call jcore%get(p, s)
    end if
    if (.not. allocated(s)) s = ''

  end function jstr

  !> Value of an integer, zero for other types
  !! @param p The value
  function jint(p) result(i)
    type(json_value), pointer, intent(in) :: p
    integer :: i
    integer(kind=json_IK) :: v

    i = 0
    if (jtype(p) .eq. json_integer) then
       call jcore%get(p, v)
       i = int(v)
    end if

  end function jint

  !> Value of an integer or real as a real, zero for other types
  !! @param p The value
  function jreal(p) result(x)
    type(json_value), pointer, intent(in) :: p
    real(kind=dp) :: x
    real(kind=json_RK) :: v

    x = 0.0_dp
    if (jtype(p) .eq. json_real) then
       call jcore%get(p, v)
       x = real(v, dp)
    else if (jtype(p) .eq. json_integer) then
       x = real(jint(p), dp)
    end if

  end function jreal

  !> Value of a boolean, false for other types
  !! @param p The value
  function jbool(p) result(b)
    type(json_value), pointer, intent(in) :: p
    logical :: b
    logical(kind=json_LK) :: v

    b = .false.
    if (jtype(p) .eq. json_logical) then
       call jcore%get(p, v)
       b = v
    end if

  end function jbool

  !> Whether a value is an integer or a real
  !! @param p The value
  function is_number(p) result(num)
    type(json_value), pointer, intent(in) :: p
    logical :: num

    num = jtype(p) .eq. json_integer .or. jtype(p) .eq. json_real

  end function is_number

  !> Value at a path of object keys, such as case.fluid.scheme
  !! @param root The value the path starts from
  !! @param path Keys separated by dots
  function node_at(root, path) result(p)
    type(json_value), pointer, intent(in) :: root
    character(len=*), intent(in) :: path
    type(json_value), pointer :: p
    integer :: i, j

    p => root
    i = 1
    do while (associated(p) .and. i .le. len(path))
       j = index(path(i:), '.')
       if (j .eq. 0) then
          p => jchild(p, path(i:))
          exit
       end if
       p => jchild(p, path(i:i + j - 2))
       i = i + j
    end do

  end function node_at

  !> Path of a key of the object at a path
  !! @param path Path of the object
  !! @param key The key
  function path_key(path, key) result(p)
    character(len=*), intent(in) :: path, key
    character(len=:), allocatable :: p

    if (len(path) .eq. 0) then
       p = key
    else
       p = path // '.' // key
    end if

  end function path_key

  !> Path of an entry of the array at a path
  !! @param path Path of the array
  !! @param i Index of the entry, starting at 1
  function path_index(path, i) result(p)
    character(len=*), intent(in) :: path
    integer, intent(in) :: i
    character(len=:), allocatable :: p

    p = path // '[' // int_str(i) // ']'

  end function path_index

  !> Path of the object or array that holds the value at a path
  !! @param path The path
  function parent_path(path) result(p)
    character(len=*), intent(in) :: path
    character(len=:), allocatable :: p
    integer :: i

    i = scan(path, '.[', back = .true.)
    if (i .eq. 0) then
       p = ''
    else
       p = path(1:i - 1)
    end if

  end function parent_path

end module checkneko_util

!> A small regular expression matcher for the "pattern" schema keyword
!! @details Supports the part of ECMA-262 regular expressions the schemas
!! need: literals and escapes, ".", character classes with ranges, groups,
!! alternation, the quantifiers *, +, ?, {n}, {n,} and {n,m}, and the anchors
!! ^ and $. The pattern is compiled into a Thompson NFA, which is simulated
!! without backtracking. As for "pattern", a match may start anywhere.
module checkneko_regex
  implicit none
  private

  integer, parameter :: OP_CHAR = 1, OP_SPLIT = 2, OP_JUMP = 3, &
       OP_BOL = 4, OP_EOL = 5, OP_MATCH = 6

  character, parameter :: BSLASH = achar(92)

  !> A state of the NFA
  type :: rx_state_t
     integer :: op = OP_JUMP
     integer :: out1 = 0
     integer :: out2 = 0
     !> Characters accepted by an OP_CHAR state
     logical :: set(0:255) = .false.
  end type rx_state_t

  !> A part of the NFA under construction: its first state and its exits,
  !! which are yet to be connected. An exit is encoded as 2*s for out1 and
  !! 2*s + 1 for out2 of state s.
  type :: frag_t
     integer :: start = 0
     integer, allocatable :: outs(:)
  end type frag_t

  !> A compiled regular expression
  type, public :: regex_t
     type(rx_state_t), allocatable :: st(:)
     integer :: n = 0
     integer :: start = 0
     !> Why the pattern could not be compiled, empty if it could
     character(len=:), allocatable :: err
   contains
     procedure, pass(this) :: compile => regex_compile
     procedure, pass(this) :: search => regex_search
  end type regex_t

  ! The pattern being compiled and the position in it
  character(len=:), allocatable :: pat
  integer :: pos

contains

  !> Compile a pattern, see err for whether it worked
  !! @param pattern The pattern
  subroutine regex_compile(this, pattern)
    class(regex_t), intent(inout) :: this
    character(len=*), intent(in) :: pattern
    type(frag_t) :: f
    integer :: m

    if (allocated(this%st)) deallocate(this%st)
    allocate(this%st(64))
    this%n = 0
    this%err = ''
    pat = pattern
    pos = 1

    f = parse_alt(this)
    if (len(this%err) .eq. 0 .and. pos .le. len(pat)) then
       this%err = 'unmatched )'
    end if
    if (len(this%err) .gt. 0) return

    m = new_state(this, OP_MATCH)
    call patch(this, f%outs, m)
    this%start = f%start

  end subroutine regex_compile

  !> Whether the pattern matches anywhere in a string
  !! @param text The string
  function regex_search(this, text) result(found)
    class(regex_t), intent(in) :: this
    character(len=*), intent(in) :: text
    logical :: found
    integer, allocatable :: clist(:), nlist(:), mark(:)
    integer :: nc, nn, gen, i, k, s, c

    found = .false.
    if (len(this%err) .gt. 0 .or. this%n .eq. 0) return

    allocate(clist(this%n), nlist(this%n), mark(this%n))
    mark = 0
    gen = 1
    nc = 0
    call add(clist, nc, this%start, 1)
    do i = 1, len(text)
       if (has_match(clist, nc)) then
          found = .true.
          return
       end if
       gen = gen + 1
       nn = 0
       c = iand(ichar(text(i:i)), 255)
       do k = 1, nc
          s = clist(k)
          if (this%st(s)%op .eq. OP_CHAR) then
             if (this%st(s)%set(c)) call add(nlist, nn, this%st(s)%out1, i + 1)
          end if
       end do
       ! A match may start at any position
       call add(nlist, nn, this%start, i + 1)
       clist(1:nn) = nlist(1:nn)
       nc = nn
    end do
    found = has_match(clist, nc)

  contains

    !> Add a state and the states reachable from it without input
    !! @param list The list of states
    !! @param nl Number of states in the list
    !! @param s The state
    !! @param p Position in the text, for the anchors
    recursive subroutine add(list, nl, s, p)
      integer, intent(inout) :: list(:)
      integer, intent(inout) :: nl
      integer, intent(in) :: s, p

      if (s .eq. 0) return
      if (mark(s) .eq. gen) return
      mark(s) = gen
      select case (this%st(s)%op)
      case (OP_SPLIT)
         call add(list, nl, this%st(s)%out1, p)
         call add(list, nl, this%st(s)%out2, p)
      case (OP_JUMP)
         call add(list, nl, this%st(s)%out1, p)
      case (OP_BOL)
         if (p .eq. 1) call add(list, nl, this%st(s)%out1, p)
      case (OP_EOL)
         if (p .eq. len(text) + 1) call add(list, nl, this%st(s)%out1, p)
      case default
         nl = nl + 1
         list(nl) = s
      end select
    end subroutine add

    !> Whether a list of states holds the final state
    !! @param list The list of states
    !! @param nl Number of states in the list
    function has_match(list, nl) result(hm)
      integer, intent(in) :: list(:)
      integer, intent(in) :: nl
      logical :: hm
      integer :: j

      hm = .false.
      do j = 1, nl
         if (this%st(list(j))%op .eq. OP_MATCH) hm = .true.
      end do
    end function has_match

  end function regex_search

  !> Parse alternatives separated by |
  recursive function parse_alt(this) result(f)
    class(regex_t), intent(inout) :: this
    type(frag_t) :: f, g
    integer :: s

    allocate(f%outs(0))
    f = parse_concat(this)
    do while (pos .le. len(pat) .and. len(this%err) .eq. 0)
       if (pat(pos:pos) .ne. '|') exit
       pos = pos + 1
       g = parse_concat(this)
       if (len(this%err) .gt. 0) return
       s = new_state(this, OP_SPLIT)
       this%st(s)%out1 = f%start
       this%st(s)%out2 = g%start
       f%start = s
       f%outs = [f%outs, g%outs]
    end do

  end function parse_alt

  !> Parse a sequence of quantified atoms
  recursive function parse_concat(this) result(f)
    class(regex_t), intent(inout) :: this
    type(frag_t) :: f, g
    logical :: have

    allocate(f%outs(0))
    have = .false.
    do while (pos .le. len(pat) .and. len(this%err) .eq. 0)
       if (pat(pos:pos) .eq. '|' .or. pat(pos:pos) .eq. ')') exit
       g = parse_repeat(this)
       if (len(this%err) .gt. 0) return
       if (have) then
          call patch(this, f%outs, g%start)
          f%outs = g%outs
       else
          f = g
          have = .true.
       end if
    end do
    if (.not. have) f = single(this, OP_JUMP)

  end function parse_concat

  !> Parse an atom and its quantifier, if any
  recursive function parse_repeat(this) result(f)
    class(regex_t), intent(inout) :: this
    type(frag_t) :: f, g, c
    integer :: p0, nmin, nmax, k
    character :: q

    p0 = pos
    f = parse_atom(this)
    if (len(this%err) .gt. 0 .or. pos .gt. len(pat)) return

    q = pat(pos:pos)
    if (q .eq. '*') then
       pos = pos + 1
       f = star(this, f)
    else if (q .eq. '+') then
       pos = pos + 1
       f = plus(this, f)
    else if (q .eq. '?') then
       pos = pos + 1
       f = quest(this, f)
    else if (q .eq. '{') then
       ! Not a valid count: the { is a literal, parsed as the next atom
       if (.not. parse_counts(nmin, nmax)) return
       if (nmin .eq. 0) then
          g = single(this, OP_JUMP)
       else
          g = f
          do k = 2, nmin
             c = copy_atom(this, p0)
             call patch(this, g%outs, c%start)
             g%outs = c%outs
          end do
       end if
       if (nmax .lt. 0) then
          c = copy_atom(this, p0)
          c = star(this, c)
          call patch(this, g%outs, c%start)
          g%outs = c%outs
       else
          do k = 1, nmax - nmin
             c = copy_atom(this, p0)
             c = quest(this, c)
             call patch(this, g%outs, c%start)
             g%outs = c%outs
          end do
       end if
       f = g
    else
       return
    end if

    ! A lazy quantifier matches the same strings
    if (pos .le. len(pat)) then
       if (pat(pos:pos) .eq. '?') pos = pos + 1
    end if

  end function parse_repeat

  !> Parse the count of a {n}, {n,} or {n,m} quantifier
  !! @details On success pos is moved past the closing brace.
  !! @param nmin Minimum count
  !! @param nmax Maximum count, -1 if there is none
  function parse_counts(nmin, nmax) result(ok)
    integer, intent(out) :: nmin, nmax
    logical :: ok
    integer :: p

    ok = .false.
    nmin = 0
    nmax = -1
    p = pos + 1
    if (.not. read_int(p, nmin)) return
    if (p .gt. len(pat)) return
    if (pat(p:p) .eq. ',') then
       p = p + 1
       if (p .gt. len(pat)) return
       if (pat(p:p) .ne. '}') then
          if (.not. read_int(p, nmax)) return
          if (nmax .lt. nmin) return
       end if
    else
       nmax = nmin
    end if
    if (p .gt. len(pat)) return
    if (pat(p:p) .ne. '}') return
    pos = p + 1
    ok = .true.

  end function parse_counts

  !> Read a decimal number from the pattern
  !! @param p Position, moved past the digits
  !! @param val The number
  function read_int(p, val) result(ok)
    integer, intent(inout) :: p
    integer, intent(out) :: val
    logical :: ok

    ok = .false.
    val = 0
    do while (p .le. len(pat))
       if (pat(p:p) .lt. '0' .or. pat(p:p) .gt. '9') exit
       val = 10 * val + (iachar(pat(p:p)) - iachar('0'))
       p = p + 1
       ok = .true.
    end do

  end function read_int

  !> Parse the atom at a position again, giving a separate copy of it
  !! @param p0 Position of the atom
  recursive function copy_atom(this, p0) result(f)
    class(regex_t), intent(inout) :: this
    integer, intent(in) :: p0
    type(frag_t) :: f
    integer :: saved

    saved = pos
    pos = p0
    f = parse_atom(this)
    pos = saved

  end function copy_atom

  !> Parse a group, a character class, an escape or a single character
  recursive function parse_atom(this) result(f)
    class(regex_t), intent(inout) :: this
    type(frag_t) :: f
    logical :: set(0:255)
    character :: c

    allocate(f%outs(0))
    c = pat(pos:pos)
    select case (c)
    case ('(')
       pos = pos + 1
       if (pos .le. len(pat)) then
          if (pat(pos:pos) .eq. '?') then
             if (pos + 1 .le. len(pat)) then
                if (pat(pos + 1:pos + 1) .eq. ':') then
                   pos = pos + 2
                else
                   this%err = 'lookaround and named groups are not supported'
                   return
                end if
             else
                this%err = 'missing )'
                return
             end if
          end if
       end if
       f = parse_alt(this)
       if (len(this%err) .gt. 0) return
       if (pos .gt. len(pat)) then
          this%err = 'missing )'
          return
       end if
       pos = pos + 1
    case ('[')
       call parse_class(this, set)
       if (len(this%err) .gt. 0) return
       f = single(this, OP_CHAR, set)
    case ('.')
       pos = pos + 1
       set = .true.
       set(10) = .false.
       set(13) = .false.
       f = single(this, OP_CHAR, set)
    case ('^')
       pos = pos + 1
       f = single(this, OP_BOL)
    case ('$')
       pos = pos + 1
       f = single(this, OP_EOL)
    case ('*', '+', '?')
       this%err = 'nothing to repeat'
    case default
       if (c .eq. BSLASH) then
          call parse_escape(this, set, .false.)
          if (len(this%err) .gt. 0) return
       else
          set = .false.
          set(code(c)) = .true.
          pos = pos + 1
       end if
       f = single(this, OP_CHAR, set)
    end select

  end function parse_atom

  !> Parse a character class such as [0-5] or [^a-z_]
  !! @param set The characters in the class
  subroutine parse_class(this, set)
    class(regex_t), intent(inout) :: this
    logical, intent(out) :: set(0:255)
    logical :: item(0:255), negate
    integer :: lo, hi

    set = .false.
    pos = pos + 1
    negate = .false.
    if (pos .le. len(pat)) then
       if (pat(pos:pos) .eq. '^') then
          negate = .true.
          pos = pos + 1
       end if
    end if

    do
       if (pos .gt. len(pat)) then
          this%err = 'missing ]'
          return
       end if
       if (pat(pos:pos) .eq. ']') exit

       if (pat(pos:pos) .eq. BSLASH) then
          call parse_escape(this, item, .true.)
          if (len(this%err) .gt. 0) return
          lo = single_char(item)
       else
          item = .false.
          lo = code(pat(pos:pos))
          item(lo) = .true.
          pos = pos + 1
       end if

       ! A range, unless the - is the last character of the class
       if (lo .ge. 0 .and. pos + 1 .le. len(pat)) then
          if (pat(pos:pos) .eq. '-' .and. pat(pos + 1:pos + 1) .ne. ']') then
             pos = pos + 1
             if (pat(pos:pos) .eq. BSLASH) then
                call parse_escape(this, item, .true.)
                if (len(this%err) .gt. 0) return
                hi = single_char(item)
             else
                hi = code(pat(pos:pos))
                pos = pos + 1
             end if
             if (hi .lt. lo) then
                this%err = 'invalid range in a character class'
                return
             end if
             item = .false.
             item(lo:hi) = .true.
          end if
       end if
       set = set .or. item
    end do
    pos = pos + 1

    if (negate) set = .not. set

  end subroutine parse_class

  !> Parse an escape sequence, starting at the backslash
  !! @param set The characters the escape stands for
  !! @param in_class Whether the escape is inside a character class
  subroutine parse_escape(this, set, in_class)
    class(regex_t), intent(inout) :: this
    logical, intent(out) :: set(0:255)
    logical, intent(in) :: in_class
    character :: c

    set = .false.
    pos = pos + 1
    if (pos .gt. len(pat)) then
       this%err = 'trailing backslash'
       return
    end if
    c = pat(pos:pos)
    pos = pos + 1

    select case (c)
    case ('d')
       set(48:57) = .true.
    case ('D')
       set = .true.
       set(48:57) = .false.
    case ('w')
       call word_set(set)
    case ('W')
       call word_set(set)
       set = .not. set
    case ('s')
       call space_set(set)
    case ('S')
       call space_set(set)
       set = .not. set
    case ('n')
       set(10) = .true.
    case ('t')
       set(9) = .true.
    case ('r')
       set(13) = .true.
    case ('f')
       set(12) = .true.
    case ('v')
       set(11) = .true.
    case ('0')
       set(0) = .true.
    case ('b')
       if (in_class) then
          set(8) = .true.
       else
          this%err = 'word boundaries are not supported'
       end if
    case default
       if (is_alnum(c)) then
          this%err = 'unsupported escape ' // BSLASH // c
       else
          set(code(c)) = .true.
       end if
    end select

  end subroutine parse_escape

  !> Characters matched by \w
  !! @param set The characters, which are added to
  subroutine word_set(set)
    logical, intent(inout) :: set(0:255)

    set(48:57) = .true.
    set(65:90) = .true.
    set(97:122) = .true.
    set(95) = .true.

  end subroutine word_set

  !> Characters matched by \s
  !! @param set The characters, which are added to
  subroutine space_set(set)
    logical, intent(inout) :: set(0:255)

    set(9:13) = .true.
    set(32) = .true.

  end subroutine space_set

  !> Code of the only character in a set, -1 if it has more or none
  !! @param set The characters
  function single_char(set) result(c)
    logical, intent(in) :: set(0:255)
    integer :: c, i

    c = -1
    if (count(set) .ne. 1) return
    do i = 0, 255
       if (set(i)) c = i
    end do

  end function single_char

  !> Whether a character is a letter or a digit
  !! @param c The character
  function is_alnum(c) result(alnum)
    character, intent(in) :: c
    logical :: alnum

    alnum = (c .ge. 'a' .and. c .le. 'z') .or. (c .ge. 'A' .and. c .le. 'Z') &
         .or. (c .ge. '0' .and. c .le. '9')

  end function is_alnum

  !> Code of a character, 0 to 255
  !! @param c The character
  function code(c) result(i)
    character, intent(in) :: c
    integer :: i

    i = iand(ichar(c), 255)

  end function code

  !> Add a state to the NFA
  !! @param op What the state does
  function new_state(this, op) result(s)
    class(regex_t), intent(inout) :: this
    integer, intent(in) :: op
    integer :: s
    type(rx_state_t), allocatable :: tmp(:)

    if (this%n .eq. size(this%st)) then
       allocate(tmp(2 * size(this%st)))
       tmp(1:this%n) = this%st(1:this%n)
       call move_alloc(tmp, this%st)
    end if
    this%n = this%n + 1
    s = this%n
    this%st(s)%op = op
    this%st(s)%out1 = 0
    this%st(s)%out2 = 0
    this%st(s)%set = .false.

  end function new_state

  !> Connect exits to a state
  !! @param outs The exits
  !! @param target The state
  subroutine patch(this, outs, target)
    class(regex_t), intent(inout) :: this
    integer, intent(in) :: outs(:)
    integer, intent(in) :: target
    integer :: k, s

    do k = 1, size(outs)
       s = outs(k) / 2
       if (mod(outs(k), 2) .eq. 0) then
          this%st(s)%out1 = target
       else
          this%st(s)%out2 = target
       end if
    end do

  end subroutine patch

  !> A fragment of a single state with one exit
  !! @param op What the state does
  !! @param set Characters accepted by an OP_CHAR state
  function single(this, op, set) result(f)
    class(regex_t), intent(inout) :: this
    integer, intent(in) :: op
    logical, intent(in), optional :: set(0:255)
    type(frag_t) :: f
    integer :: s

    s = new_state(this, op)
    if (present(set)) this%st(s)%set = set
    f%start = s
    f%outs = [2 * s]

  end function single

  !> A fragment repeated zero or more times
  !! @param a The fragment
  function star(this, a) result(f)
    class(regex_t), intent(inout) :: this
    type(frag_t), intent(in) :: a
    type(frag_t) :: f
    integer :: s

    s = new_state(this, OP_SPLIT)
    this%st(s)%out1 = a%start
    call patch(this, a%outs, s)
    f%start = s
    f%outs = [2 * s + 1]

  end function star

  !> A fragment repeated one or more times
  !! @param a The fragment
  function plus(this, a) result(f)
    class(regex_t), intent(inout) :: this
    type(frag_t), intent(in) :: a
    type(frag_t) :: f
    integer :: s

    s = new_state(this, OP_SPLIT)
    this%st(s)%out1 = a%start
    call patch(this, a%outs, s)
    f%start = a%start
    f%outs = [2 * s + 1]

  end function plus

  !> An optional fragment
  !! @param a The fragment
  function quest(this, a) result(f)
    class(regex_t), intent(inout) :: this
    type(frag_t), intent(in) :: a
    type(frag_t) :: f
    integer :: s

    s = new_state(this, OP_SPLIT)
    this%st(s)%out1 = a%start
    f%start = s
    f%outs = [a%outs, 2 * s + 1]

  end function quest

end module checkneko_regex

!> Validation of a case file against the JSON schemas in doc/schemas
!! @details Implements the part of JSON Schema 2020-12 the schemas use, with
!! Neko's distinction between JSON integers and reals: "integer" accepts only
!! integers such as 10, not 10.0, and "x-neko-real" rejects integers where
!! Neko reads a real. Keywords the validator does not implement are reported
!! when the schemas are loaded, so that they cannot be ignored silently.
!!
!! Besides the messages, the validation records annotations: the keys each
!! object is declared to have, which the checks use to spot misspelled keys
!! in objects that allow any key, and the strings that name a constant from
!! case.constants, with the kind of value Neko looks the constant up as.
module checkneko_schema
  use json_module, only : json_value, json_null, json_object, json_array, &
       json_logical, json_integer, json_real, json_string
  use checkneko_util
  use checkneko_regex, only : regex_t
  implicit none
  private

  !> What an annotation records
  integer, public, parameter :: AK_DECLARED = 1, AK_CONSTREF = 2

  !> How a constant is looked up: as an integer, a real, either, or as an
  !! array of integers or reals
  integer, public, parameter :: CK_ANY = 0, CK_INTEGER = 1, CK_REAL = 2, &
       CK_NUMBER = 3, CK_INTEGER_ARRAY = 4, CK_REAL_ARRAY = 5

  !> Information collected while validating
  type, public :: ann_t
     integer :: kind = AK_DECLARED
     !> The object (AK_DECLARED) or the string naming a constant (AK_CONSTREF)
     type(json_value), pointer :: node => null()
     character(len=:), allocatable :: path
     !> The declared key, for AK_DECLARED
     character(len=:), allocatable :: name
     !> How the constant is looked up, for AK_CONSTREF
     integer :: ck = CK_ANY
     integer :: min_items = -1
     integer :: max_items = -1
  end type ann_t

  !> A list of annotations
  type, public :: ann_list_t
     type(ann_t), allocatable :: a(:)
     integer :: n = 0
   contains
     procedure, pass(this) :: add => ann_list_add
     procedure, pass(this) :: truncate => ann_list_truncate
  end type ann_list_t

  !> The annotations of the case file being checked
  type(ann_list_t), public :: anns

  !> A loaded schema
  type :: doc_t
     character(len=:), allocatable :: id
     type(json_value), pointer :: root => null()
  end type doc_t

  !> A compiled "pattern"
  type :: pattern_t
     character(len=:), allocatable :: src
     type(regex_t) :: rx
  end type pattern_t

  !> Id of the schema of the whole case file
  character(len=*), parameter :: ROOT_ID = 'urn:neko:schema:case-file'
  !> Prefix of all schema ids, the rest gives the file name
  character(len=*), parameter :: ID_PREFIX = 'urn:neko:schema:'

  type(doc_t), allocatable :: docs(:)
  integer :: ndocs = 0
  type(pattern_t), allocatable :: patterns(:)
  integer :: npatterns = 0
  !> Directory the schemas are read from
  character(len=:), allocatable :: schema_dir

  public :: schema_init, schema_validate, schema_free, schema_count

contains

  !> Load the schemas, starting from the one of the whole case file
  !! @details Every schema reachable through $ref is loaded and checked for
  !! keywords the validator does not support.
  !! @param dir Directory to read the schemas from, e.g. doc/schemas
  subroutine schema_init(dir)
    character(len=*), intent(in) :: dir
    integer :: k

    schema_dir = dir
    allocate(docs(64))
    allocate(patterns(8))
    ndocs = 0
    npatterns = 0
    k = doc_index(ROOT_ID)

  end subroutine schema_init

  !> Number of loaded schemas
  function schema_count() result(n)
    integer :: n

    n = ndocs

  end function schema_count

  !> Free the schemas
  subroutine schema_free()
    integer :: k

    do k = 1, ndocs
       call jcore%destroy(docs(k)%root)
    end do
    if (allocated(docs)) deallocate(docs)
    if (allocated(patterns)) deallocate(patterns)
    ndocs = 0
    npatterns = 0

  end subroutine schema_free

  !> Validate a case file, adding messages to msgs and annotations to anns
  !! @param inst The case file
  subroutine schema_validate(inst)
    type(json_value), pointer, intent(in) :: inst
    integer :: k

    k = doc_index(ROOT_ID)
    call validate(inst, '', docs(k)%root, k)

  end subroutine schema_validate

  !> Index of a schema, which is loaded the first time it is asked for
  !! @param id Id of the schema
  recursive function doc_index(id) result(k)
    character(len=*), intent(in) :: id
    integer :: k
    type(json_value), pointer :: root
    type(doc_t), allocatable :: tmp(:)
    character(len=:), allocatable :: fname, msg
    logical :: ok

    do k = 1, ndocs
       if (str_eq(docs(k)%id, id)) return
    end do

    root => null()
    fname = join_path(schema_dir, id_file(id))
    if (.not. file_exists(fname)) then
       call fatal('cannot find the schema ' // id // ' (' // fname // ')')
    end if
    call jcore%load(file = fname, p = root)
    if (jcore%failed()) then
       call jcore%check_for_errors(ok, msg)
       call fatal('cannot parse the schema ' // id // ': ' // msg)
    end if

    if (ndocs .eq. size(docs)) then
       allocate(tmp(2 * size(docs)))
       tmp(1:ndocs) = docs(1:ndocs)
       call move_alloc(tmp, docs)
    end if
    ndocs = ndocs + 1
    k = ndocs
    docs(k)%id = id
    docs(k)%root => root

    call check_keywords(root, k)

  end function doc_index

  !> File of a schema, relative to the schema directory
  !! @details urn:neko:schema:simcomps:probes is simcomps/probes.schema.json
  !! @param id Id of the schema
  function id_file(id) result(fname)
    character(len=*), intent(in) :: id
    character(len=:), allocatable :: fname
    integer :: i

    fname = id
    if (len(id) .gt. len(ID_PREFIX)) then
       if (id(1:len(ID_PREFIX)) .eq. ID_PREFIX) fname = id(len(ID_PREFIX) + 1:)
    end if
    do i = 1, len(fname)
       if (fname(i:i) .eq. ':') fname(i:i) = '/'
    end do
    fname = fname // '.schema.json'

  end function id_file

  !> Check that a schema uses only keywords the validator implements
  !! @details Also loads the schemas its $refs point to, checks that the
  !! $refs resolve and compiles its patterns.
  !! @param s The schema
  !! @param doc Index of the schema document s is part of
  recursive subroutine check_keywords(s, doc)
    type(json_value), pointer, intent(in) :: s
    integer, intent(in) :: doc
    type(json_value), pointer :: c, e, t
    character(len=:), allocatable :: key
    integer :: tdoc, k

    if (jtype(s) .ne. json_object) return

    c => jfirst(s)
    do while (associated(c))
       key = jname(c)
       select case (key)
       case ('properties', '$defs')
          e => jfirst(c)
          do while (associated(e))
             call check_keywords(e, doc)
             e => jnext(e)
          end do
       case ('allOf', 'anyOf', 'oneOf')
          e => jfirst(c)
          do while (associated(e))
             call check_keywords(e, doc)
             e => jnext(e)
          end do
       case ('not', 'if', 'then', 'else', 'items', 'additionalProperties')
          call check_keywords(c, doc)
       case ('$ref')
          call resolve(jstr(c), doc, t, tdoc)
       case ('pattern')
          k = pattern_index(jstr(c))
          if (len(patterns(k)%rx%err) .gt. 0) then
             call warn('cannot compile the pattern ' // jstr(c) // &
                  ' in the schema ' // docs(doc)%id // ': ' // &
                  patterns(k)%rx%err // ', it is not checked')
          end if
       case ('type', 'enum', 'const', 'required', 'dependentRequired', &
            'minItems', 'maxItems', 'minimum', 'maximum', &
            'exclusiveMinimum', 'exclusiveMaximum', 'minLength', &
            'maxLength', 'x-neko-real')
          ! Validated
       case ('$schema', '$id', '$comment', 'title', 'description', &
            'default', 'examples', 'format', 'deprecated')
          ! Annotations only
       case default
          call warn('the schema ' // docs(doc)%id // ' uses the keyword "' &
               // key // '", which checkneko does not support, so ' // &
               'it is ignored')
       end select
       c => jnext(c)
    end do

  end subroutine check_keywords

  !> Resolve a $ref, such as urn:neko:schema:common#/$defs/constRef
  !! @param ref The reference
  !! @param doc Index of the schema document the reference is in
  !! @param t The schema referred to
  !! @param tdoc Index of the schema document t is part of
  recursive subroutine resolve(ref, doc, t, tdoc)
    character(len=*), intent(in) :: ref
    integer, intent(in) :: doc
    type(json_value), pointer, intent(out) :: t
    integer, intent(out) :: tdoc
    character(len=:), allocatable :: frag, seg
    integer :: i, p, q, idx, ios

    i = index(ref, '#')
    if (i .eq. 0) then
       tdoc = doc_index(ref)
       frag = ''
    else if (i .eq. 1) then
       tdoc = doc
       frag = ref(2:)
    else
       tdoc = doc_index(ref(1:i - 1))
       frag = ref(i + 1:)
    end if
    t => docs(tdoc)%root

    if (len(frag) .eq. 0) return
    if (frag(1:1) .ne. '/') then
       call fatal('cannot resolve $ref ' // ref // ' in ' // docs(doc)%id)
    end if

    ! A JSON pointer, with ~1 for / and ~0 for ~ in its segments
    p = 2
    do while (p .le. len(frag) + 1)
       q = index(frag(p:), '/')
       if (q .eq. 0) then
          seg = frag(p:)
          p = len(frag) + 2
       else
          seg = frag(p:p + q - 2)
          p = p + q
       end if
       seg = unescape_pointer(seg)
       if (jtype(t) .eq. json_object) then
          t => jchild(t, seg)
       else if (jtype(t) .eq. json_array) then
          read(seg, *, iostat = ios) idx
          if (ios .ne. 0) then
             t => null()
          else
             t => jelement(t, idx + 1)
          end if
       else
          t => null()
       end if
       if (.not. associated(t)) then
          call fatal('cannot resolve $ref ' // ref // ' in ' // docs(doc)%id)
       end if
    end do

  end subroutine resolve

  !> A segment of a JSON pointer with its escapes replaced
  !! @param seg The segment
  function unescape_pointer(seg) result(s)
    character(len=*), intent(in) :: seg
    character(len=:), allocatable :: s
    integer :: i

    s = ''
    i = 1
    do while (i .le. len(seg))
       if (seg(i:i) .eq. '~' .and. i .lt. len(seg)) then
          if (seg(i + 1:i + 1) .eq. '1') then
             s = s // '/'
          else
             s = s // '~'
          end if
          i = i + 2
       else
          s = s // seg(i:i)
          i = i + 1
       end if
    end do

  end function unescape_pointer

  !> Index of a compiled pattern, compiled the first time it is asked for
  !! @param src The pattern
  function pattern_index(src) result(k)
    character(len=*), intent(in) :: src
    integer :: k
    type(pattern_t), allocatable :: tmp(:)

    do k = 1, npatterns
       if (str_eq(patterns(k)%src, src)) return
    end do
    if (npatterns .eq. size(patterns)) then
       allocate(tmp(2 * size(patterns)))
       tmp(1:npatterns) = patterns(1:npatterns)
       call move_alloc(tmp, patterns)
    end if
    npatterns = npatterns + 1
    k = npatterns
    patterns(k)%src = src
    call patterns(k)%rx%compile(src)

  end function pattern_index

  !> Validate a value against a schema
  !! @param inst The value
  !! @param path Where the value is in the case file
  !! @param sch The schema
  !! @param doc Index of the schema document sch is part of
  recursive subroutine validate(inst, path, sch, doc)
    type(json_value), pointer, intent(in) :: inst
    character(len=*), intent(in) :: path
    type(json_value), pointer, intent(in) :: sch
    integer, intent(in) :: doc

    select case (jtype(sch))
    case (json_logical)
       if (.not. jbool(sch)) call not_allowed(inst, path)
       return
    case (json_object)
    case default
       return
    end select

    call kw_ref(inst, path, sch, doc)
    call kw_type(inst, path, sch)
    call kw_const_enum(inst, path, sch)
    call kw_number(inst, path, sch)
    call kw_string(inst, path, sch)
    call kw_object(inst, path, sch, doc)
    call kw_array(inst, path, sch, doc)
    call kw_combinators(inst, path, sch, doc)
    call kw_conditional(inst, path, sch, doc)

  end subroutine validate

  !> Report a value the schema false rejects
  !! @param inst The value
  !! @param path Where the value is
  subroutine not_allowed(inst, path)
    type(json_value), pointer, intent(in) :: inst
    character(len=*), intent(in) :: path
    integer :: i

    i = scan(path, '.[', back = .true.)
    if (i .gt. 0) then
       if (path(i:i) .eq. '.') then
          call report(SEV_ERROR, path, 'the key ' // quoted(path(i + 1:)) // &
               ' is not supported here', inst, MK_UNKNOWN_KEY)
          return
       end if
    else if (len(path) .gt. 0) then
       call report(SEV_ERROR, path, 'the key ' // quoted(path) // &
            ' is not supported here', inst, MK_UNKNOWN_KEY)
       return
    end if
    call report(SEV_ERROR, path, 'this value is not allowed here', inst)

  end subroutine not_allowed

  !> $ref: validate against the schema referred to
  !! @param inst The value
  !! @param path Where the value is in the case file
  !! @param sch The schema
  !! @param doc Index of the schema document sch is part of
  recursive subroutine kw_ref(inst, path, sch, doc)
    type(json_value), pointer, intent(in) :: inst
    character(len=*), intent(in) :: path
    type(json_value), pointer, intent(in) :: sch
    integer, intent(in) :: doc
    type(json_value), pointer :: p, t
    character(len=:), allocatable :: ref
    integer :: tdoc, n0, j, k
    logical :: const_ref
    type(ann_t) :: ann

    p => jchild(sch, '$ref')
    if (.not. associated(p)) return
    ref = jstr(p)
    call resolve(ref, doc, t, tdoc)

    const_ref = .false.
    if (len(ref) .ge. 15) const_ref = ref(len(ref) - 14:) .eq. '/$defs/constRef'

    n0 = msgs%n
    call validate(inst, path, t, tdoc)

    if (.not. const_ref) return
    if (msgs%n .eq. n0) then
       ann%kind = AK_CONSTREF
       ann%node => inst
       ann%path = path
       call anns%add(ann)
    else
       ! Say "the name of a constant" rather than "a string"
       do j = n0 + 1, msgs%n
          if (msgs%m(j)%kind .ne. MK_TYPE) cycle
          if (.not. allocated(msgs%m(j)%allowed)) cycle
          do k = 1, size(msgs%m(j)%allowed)
             if (msgs%m(j)%allowed(k)%s .eq. 'string') then
                msgs%m(j)%allowed(k)%s = 'constant'
             end if
          end do
          msgs%m(j)%text = 'expected ' // &
               type_phrase(msgs%m(j)%allowed) // ', got ' // describe(inst)
       end do
    end if

  end subroutine kw_ref

  !> type
  !! @param inst The value
  !! @param path Where the value is in the case file
  !! @param sch The schema
  subroutine kw_type(inst, path, sch)
    type(json_value), pointer, intent(in) :: inst
    character(len=*), intent(in) :: path
    type(json_value), pointer, intent(in) :: sch
    type(json_value), pointer :: p, e
    type(str_t), allocatable :: names(:)
    logical :: ok

    p => jchild(sch, 'type')
    if (.not. associated(p)) return

    ok = .false.
    if (jtype(p) .eq. json_string) then
       call add_str(names, jstr(p))
       ok = type_matches(inst, jstr(p))
    else
       e => jfirst(p)
       do while (associated(e))
          call add_str(names, jstr(e))
          ok = ok .or. type_matches(inst, jstr(e))
          e => jnext(e)
       end do
    end if
    if (.not. ok) call type_error(inst, path, names)

  end subroutine kw_type

  !> Whether a value is of a schema type
  !! @param inst The value
  !! @param name The type
  function type_matches(inst, name) result(match)
    type(json_value), pointer, intent(in) :: inst
    character(len=*), intent(in) :: name
    logical :: match

    select case (name)
    case ('object')
       match = jtype(inst) .eq. json_object
    case ('array')
       match = jtype(inst) .eq. json_array
    case ('string')
       match = jtype(inst) .eq. json_string
    case ('boolean')
       match = jtype(inst) .eq. json_logical
    case ('null')
       match = jtype(inst) .eq. json_null
    case ('integer')
       match = jtype(inst) .eq. json_integer
    case ('number')
       match = is_number(inst)
    case default
       match = .true.
    end select

  end function type_matches

  !> Report a value of the wrong type
  !! @param inst The value
  !! @param path Where the value is
  !! @param names The accepted types
  subroutine type_error(inst, path, names)
    type(json_value), pointer, intent(in) :: inst
    character(len=*), intent(in) :: path
    type(str_t), intent(in) :: names(:)

    call report(SEV_ERROR, path, 'expected ' // type_phrase(names) // &
         ', got ' // describe(inst) // integer_hint(inst, names), inst, &
         MK_TYPE, describe(inst), names)

  end subroutine type_error

  !> A hint for a real such as 7.0 where Neko reads an integer
  !! @param inst The value
  !! @param names The accepted types
  function integer_hint(inst, names) result(hint)
    type(json_value), pointer, intent(in) :: inst
    type(str_t), intent(in) :: names(:)
    character(len=:), allocatable :: hint
    logical :: int_only
    real(kind=dp) :: x
    integer :: i

    hint = ''
    if (jtype(inst) .ne. json_real) return
    int_only = .false.
    do i = 1, size(names)
       if (names(i)%s .eq. 'integer') int_only = .true.
       if (names(i)%s .eq. 'number') return
    end do
    x = jreal(inst)
    if (int_only .and. is_whole(x) .and. abs(x) .lt. 2.0e9_dp) then
       hint = ' (write ' // int_str(nint(x)) // ', without a decimal point)'
    end if

  end function integer_hint

  !> const and enum
  !! @param inst The value
  !! @param path Where the value is in the case file
  !! @param sch The schema
  subroutine kw_const_enum(inst, path, sch)
    type(json_value), pointer, intent(in) :: inst
    character(len=*), intent(in) :: path
    type(json_value), pointer, intent(in) :: sch
    type(json_value), pointer :: p, e
    type(str_t), allocatable :: allowed(:)
    logical :: ok

    p => jchild(sch, 'const')
    if (associated(p)) then
       if (.not. json_equal(inst, p)) then
          call add_str(allowed, render(p))
          call report(SEV_ERROR, path, 'must be ' // render(p) // ', got ' &
               // render(inst), inst, MK_CONST, render(inst), allowed)
          deallocate(allowed)
       end if
    end if

    p => jchild(sch, 'enum')
    if (associated(p)) then
       ok = .false.
       e => jfirst(p)
       do while (associated(e))
          call add_str(allowed, render(e))
          if (json_equal(inst, e)) ok = .true.
          e => jnext(e)
       end do
       if (.not. ok .and. allocated(allowed)) then
          call report(SEV_ERROR, path, choice_text(render(inst), allowed), &
               inst, MK_ENUM, render(inst), allowed)
       end if
    end if

  end subroutine kw_const_enum

  !> Text for a value that is not one of the accepted ones
  !! @param got The value, rendered
  !! @param allowed The accepted values, rendered
  function choice_text(got, allowed) result(text)
    character(len=*), intent(in) :: got
    type(str_t), intent(in) :: allowed(:)
    character(len=:), allocatable :: text
    integer :: i, best, d, dbest

    if (size(allowed) .eq. 1) then
       text = 'must be ' // allowed(1)%s // ', got ' // got
    else
       text = 'must be one of ' // value_list(allowed) // ', got ' // got
    end if

    ! Suggest the closest string, as for a misspelled keyword
    if (len(got) .lt. 2 .or. got(1:1) .ne. '"') return
    best = 0
    dbest = huge(1)
    do i = 1, size(allowed)
       if (allowed(i)%s(1:1) .ne. '"') cycle
       d = edit_distance(got, allowed(i)%s)
       if (d .lt. dbest) then
          dbest = d
          best = i
       end if
    end do
    if (best .gt. 0 .and. dbest .le. typo_limit(len(got) - 2) .and. &
         dbest .lt. len(got) - 2) then
       text = text // ' (did you mean ' // allowed(best)%s // '?)'
    end if

  end function choice_text

  !> Whether two JSON values are equal, with 1 equal to 1.0 as in JSON Schema
  !! @param a First value
  !! @param b Second value
  recursive function json_equal(a, b) result(eq)
    type(json_value), pointer, intent(in) :: a, b
    logical :: eq
    type(json_value), pointer :: ea, eb

    eq = .false.
    if (is_number(a) .and. is_number(b)) then
       eq = .not. (jreal(a) .lt. jreal(b) .or. jreal(a) .gt. jreal(b))
       return
    end if
    if (jtype(a) .ne. jtype(b)) return

    select case (jtype(a))
    case (json_string)
       eq = str_eq(jstr(a), jstr(b))
    case (json_logical)
       eq = jbool(a) .eqv. jbool(b)
    case (json_null)
       eq = .true.
    case (json_array)
       if (jnchildren(a) .ne. jnchildren(b)) return
       ea => jfirst(a)
       eb => jfirst(b)
       do while (associated(ea))
          if (.not. json_equal(ea, eb)) return
          ea => jnext(ea)
          eb => jnext(eb)
       end do
       eq = .true.
    case (json_object)
       if (jnchildren(a) .ne. jnchildren(b)) return
       ea => jfirst(a)
       do while (associated(ea))
          eb => jchild(b, jname(ea))
          if (.not. associated(eb)) return
          if (.not. json_equal(ea, eb)) return
          ea => jnext(ea)
       end do
       eq = .true.
    end select

  end function json_equal

  !> minimum, maximum, exclusiveMinimum, exclusiveMaximum and x-neko-real
  !! @param inst The value
  !! @param path Where the value is in the case file
  !! @param sch The schema
  subroutine kw_number(inst, path, sch)
    type(json_value), pointer, intent(in) :: inst
    character(len=*), intent(in) :: path
    type(json_value), pointer, intent(in) :: sch
    type(json_value), pointer :: p
    real(kind=dp) :: x

    if (.not. is_number(inst)) return
    x = jreal(inst)

    p => jchild(sch, 'minimum')
    if (associated(p)) then
       if (x .lt. jreal(p)) call report(SEV_ERROR, path, 'must be at least ' &
            // render(p) // ', got ' // render(inst), inst)
    end if
    p => jchild(sch, 'exclusiveMinimum')
    if (associated(p)) then
       if (x .le. jreal(p)) call report(SEV_ERROR, path, &
            'must be greater than ' // render(p) // ', got ' // &
            render(inst), inst)
    end if
    p => jchild(sch, 'maximum')
    if (associated(p)) then
       if (x .gt. jreal(p)) call report(SEV_ERROR, path, 'must be at most ' &
            // render(p) // ', got ' // render(inst), inst)
    end if
    p => jchild(sch, 'exclusiveMaximum')
    if (associated(p)) then
       if (x .ge. jreal(p)) call report(SEV_ERROR, path, &
            'must be less than ' // render(p) // ', got ' // &
            render(inst), inst)
    end if

    ! Neko reads the value as a real, which it does not accept an integer for
    p => jchild(sch, 'x-neko-real')
    if (associated(p)) then
       if (jbool(p) .and. jtype(inst) .eq. json_integer) then
          call report(SEV_ERROR, path, 'expected a real number, got ' // &
               describe(inst) // ' (write ' // render(inst) // '.0)', inst)
       end if
    end if

  end subroutine kw_number

  !> minLength, maxLength and pattern
  !! @param inst The value
  !! @param path Where the value is in the case file
  !! @param sch The schema
  subroutine kw_string(inst, path, sch)
    type(json_value), pointer, intent(in) :: inst
    character(len=*), intent(in) :: path
    type(json_value), pointer, intent(in) :: sch
    type(json_value), pointer :: p
    character(len=:), allocatable :: s
    integer :: n, i, k

    if (jtype(inst) .ne. json_string) return
    s = jstr(inst)

    ! Length in characters, not counting UTF-8 continuation bytes
    n = 0
    do i = 1, len(s)
       if (iand(ichar(s(i:i)), 192) .ne. 128) n = n + 1
    end do

    p => jchild(sch, 'minLength')
    if (associated(p)) then
       if (n .lt. jint(p)) then
          if (jint(p) .eq. 1) then
             call report(SEV_ERROR, path, 'must not be empty', inst)
          else
             call report(SEV_ERROR, path, 'must be at least ' // &
                  int_str(jint(p)) // ' characters long', inst)
          end if
       end if
    end if
    p => jchild(sch, 'maxLength')
    if (associated(p)) then
       if (n .gt. jint(p)) call report(SEV_ERROR, path, &
            'must be at most ' // int_str(jint(p)) // ' characters long', inst)
    end if

    p => jchild(sch, 'pattern')
    if (associated(p)) then
       k = pattern_index(jstr(p))
       if (len(patterns(k)%rx%err) .eq. 0) then
          if (.not. patterns(k)%rx%search(s)) then
             call report(SEV_ERROR, path, render(inst) // &
                  ' does not have the expected form (the pattern ' // &
                  jstr(p) // ')', inst)
          end if
       end if
    end if

  end subroutine kw_string

  !> properties, additionalProperties, required and dependentRequired
  !! @param inst The value
  !! @param path Where the value is in the case file
  !! @param sch The schema
  !! @param doc Index of the schema document sch is part of
  recursive subroutine kw_object(inst, path, sch, doc)
    type(json_value), pointer, intent(in) :: inst
    character(len=*), intent(in) :: path
    type(json_value), pointer, intent(in) :: sch
    integer, intent(in) :: doc
    type(json_value), pointer :: props, addl, req, dep, c, sub, e, d
    character(len=:), allocatable :: key, name
    type(ann_t) :: ann

    if (jtype(inst) .ne. json_object) return

    props => jchild(sch, 'properties')
    addl => jchild(sch, 'additionalProperties')
    req => jchild(sch, 'required')
    dep => jchild(sch, 'dependentRequired')

    if (associated(props)) then
       e => jfirst(props)
       do while (associated(e))
          ann%kind = AK_DECLARED
          ann%node => inst
          ann%path = path
          ann%name = jname(e)
          call anns%add(ann)
          e => jnext(e)
       end do
    end if

    ! Every key, also one that appears twice
    c => jfirst(inst)
    do while (associated(c))
       key = jname(c)
       sub => jchild(props, key)
       if (associated(sub)) then
          call validate(c, path_key(path, key), sub, doc)
       else if (associated(addl)) then
          if (jtype(addl) .eq. json_logical .and. .not. jbool(addl)) then
             call report(SEV_ERROR, path_key(path, key), 'unknown key ' // &
                  quoted(key) // suggestion(key, props), c, MK_UNKNOWN_KEY)
          else
             call validate(c, path_key(path, key), addl, doc)
          end if
       end if
       c => jnext(c)
    end do

    if (associated(req)) then
       e => jfirst(req)
       do while (associated(e))
          name = jstr(e)
          if (.not. associated(jchild(inst, name))) then
             call report(SEV_ERROR, path, 'missing the required key ' // &
                  quoted(name) // misspelled_as(name, inst, props), inst, &
                  MK_REQUIRED, allowed = [str_t(quoted(name))])
          end if
          e => jnext(e)
       end do
    end if

    if (associated(dep)) then
       d => jfirst(dep)
       do while (associated(d))
          if (associated(jchild(inst, jname(d)))) then
             e => jfirst(d)
             do while (associated(e))
                name = jstr(e)
                if (.not. associated(jchild(inst, name))) then
                   call report(SEV_ERROR, path, 'missing the key ' // &
                        quoted(name) // ', which is required when ' // &
                        quoted(jname(d)) // ' is given', inst, MK_REQUIRED)
                end if
                e => jnext(e)
             end do
          end if
          d => jnext(d)
       end do
    end if

  end subroutine kw_object

  !> A suggestion for a misspelled key, from the declared keys
  !! @param key The key
  !! @param props The declared keys, as the properties of a schema
  function suggestion(key, props) result(text)
    character(len=*), intent(in) :: key
    type(json_value), pointer, intent(in) :: props
    character(len=:), allocatable :: text, best
    type(json_value), pointer :: e
    integer :: d, dbest

    text = ''
    if (.not. associated(props)) return
    dbest = huge(1)
    e => jfirst(props)
    do while (associated(e))
       d = edit_distance(key, jname(e))
       if (d .lt. dbest) then
          dbest = d
          best = jname(e)
       end if
       e => jnext(e)
    end do
    if (dbest .le. typo_limit(len(key)) .and. dbest .lt. len(key)) then
       text = ' (did you mean ' // quoted(best) // '?)'
    end if

  end function suggestion

  !> A hint that a missing key is there with a misspelled name
  !! @param name The missing key
  !! @param inst The object
  !! @param props The declared keys, as the properties of a schema
  function misspelled_as(name, inst, props) result(text)
    character(len=*), intent(in) :: name
    type(json_value), pointer, intent(in) :: inst, props
    character(len=:), allocatable :: text, key
    type(json_value), pointer :: c
    integer :: d, dbest

    text = ''
    dbest = huge(1)
    c => jfirst(inst)
    do while (associated(c))
       key = jname(c)
       if (.not. associated(jchild(props, key))) then
          d = edit_distance(name, key)
          if (d .lt. dbest .and. d .le. typo_limit(len(name)) .and. &
               d .lt. len(name)) then
             dbest = d
             text = ' (is ' // quoted(key) // ' a misspelling of it?)'
          end if
       end if
       c => jnext(c)
    end do

  end function misspelled_as

  !> items, minItems and maxItems
  !! @param inst The value
  !! @param path Where the value is in the case file
  !! @param sch The schema
  !! @param doc Index of the schema document sch is part of
  recursive subroutine kw_array(inst, path, sch, doc)
    type(json_value), pointer, intent(in) :: inst
    character(len=*), intent(in) :: path
    type(json_value), pointer, intent(in) :: sch
    integer, intent(in) :: doc
    type(json_value), pointer :: items, pmin, pmax, e
    integer :: n, i

    if (jtype(inst) .ne. json_array) return
    n = jnchildren(inst)

    items => jchild(sch, 'items')
    if (associated(items)) then
       e => jfirst(inst)
       i = 0
       do while (associated(e))
          i = i + 1
          call validate(e, path_index(path, i), items, doc)
          e => jnext(e)
       end do
    end if

    pmin => jchild(sch, 'minItems')
    pmax => jchild(sch, 'maxItems')
    if (associated(pmin) .and. associated(pmax)) then
       if (jint(pmin) .eq. jint(pmax) .and. n .ne. jint(pmin)) then
          call report(SEV_ERROR, path, 'must have exactly ' // &
               entries(jint(pmin)) // ', got ' // int_str(n), inst)
          return
       end if
    end if
    if (associated(pmin)) then
       if (n .lt. jint(pmin)) call report(SEV_ERROR, path, &
            'must have at least ' // entries(jint(pmin)) // ', got ' // &
            int_str(n), inst)
    end if
    if (associated(pmax)) then
       if (n .gt. jint(pmax)) call report(SEV_ERROR, path, &
            'must have at most ' // entries(jint(pmax)) // ', got ' // &
            int_str(n), inst)
    end if

  end subroutine kw_array

  !> A number of array entries as text
  !! @param n The number
  function entries(n) result(text)
    integer, intent(in) :: n
    character(len=:), allocatable :: text

    if (n .eq. 1) then
       text = '1 entry'
    else
       text = int_str(n) // ' entries'
    end if

  end function entries

  !> allOf, anyOf, oneOf and not
  !! @param inst The value
  !! @param path Where the value is in the case file
  !! @param sch The schema
  !! @param doc Index of the schema document sch is part of
  recursive subroutine kw_combinators(inst, path, sch, doc)
    type(json_value), pointer, intent(in) :: inst
    character(len=*), intent(in) :: path
    type(json_value), pointer, intent(in) :: sch
    integer, intent(in) :: doc
    type(json_value), pointer :: p, e
    integer :: n0, a0
    logical :: passed

    p => jchild(sch, 'allOf')
    if (associated(p)) then
       e => jfirst(p)
       do while (associated(e))
          call validate(inst, path, e, doc)
          e => jnext(e)
       end do
    end if

    p => jchild(sch, 'anyOf')
    if (associated(p)) call alternatives(inst, path, p, doc, .false.)

    p => jchild(sch, 'oneOf')
    if (associated(p)) call alternatives(inst, path, p, doc, .true.)

    p => jchild(sch, 'not')
    if (associated(p)) then
       n0 = msgs%n
       a0 = anns%n
       call validate(inst, path, p, doc)
       passed = msgs%n .eq. n0
       call msgs%truncate(n0)
       call anns%truncate(a0)
       if (passed) call report(SEV_ERROR, path, not_text(p), inst, MK_NOT)
    end if

  end subroutine kw_combinators

  !> Text for a value that matches a schema it must not match
  !! @param p The schema of the not
  function not_text(p) result(text)
    type(json_value), pointer, intent(in) :: p
    character(len=:), allocatable :: text
    type(json_value), pointer :: req

    req => jchild(p, 'required')
    if (jnchildren(p) .eq. 1 .and. jnchildren(req) .eq. 1) then
       text = 'the key ' // quoted(jstr(jfirst(req))) // &
            ' is not allowed here'
    else
       text = 'this combination of settings is not allowed'
    end if

  end function not_text

  !> anyOf and oneOf: at least or exactly one of the schemas must match
  !! @details If none matches, the messages of the alternative the user most
  !! likely meant are kept, see summarise.
  !! @param inst The value
  !! @param path Where the value is
  !! @param arr The schemas
  !! @param doc Index of the schema document arr is part of
  !! @param exactly_one Whether exactly one must match (oneOf)
  recursive subroutine alternatives(inst, path, arr, doc, exactly_one)
    type(json_value), pointer, intent(in) :: inst
    character(len=*), intent(in) :: path
    type(json_value), pointer, intent(in) :: arr
    integer, intent(in) :: doc
    logical, intent(in) :: exactly_one
    type(msg_list_t), allocatable :: br(:)
    integer, allocatable :: score(:)
    type(json_value), pointer :: e
    integer :: nb, k, j, n0, a0, ab, npass, kpass

    nb = jnchildren(arr)
    if (nb .eq. 0) return
    allocate(br(nb), score(nb))
    score = 0
    n0 = msgs%n
    a0 = anns%n
    npass = 0
    kpass = 0

    e => jfirst(arr)
    k = 0
    do while (associated(e))
       k = k + 1
       ab = anns%n
       call validate(inst, path, e, doc)
       if (msgs%n .eq. n0) then
          npass = npass + 1
          if (kpass .eq. 0) kpass = k
       else
          score(k) = declared_keys(inst, ab)
          do j = n0 + 1, msgs%n
             call br(k)%add(msgs%m(j))
          end do
          call msgs%truncate(n0)
          call anns%truncate(ab)
       end if
       e => jnext(e)
    end do

    if (npass .gt. 0) then
       call set_const_kinds(inst, arr, kpass, doc, a0)
       if (exactly_one .and. npass .gt. 1) then
          call report(SEV_ERROR, path, &
               'matches more than one of the allowed forms', inst)
       end if
    else
       call summarise(inst, path, br, score)
    end if

  end subroutine alternatives

  !> Number of keys of an object declared by the annotations after a0
  !! @param inst The object
  !! @param a0 Number of annotations before
  function declared_keys(inst, a0) result(n)
    type(json_value), pointer, intent(in) :: inst
    integer, intent(in) :: a0
    integer :: n
    type(json_value), pointer :: c
    integer :: k

    n = 0
    if (jtype(inst) .ne. json_object) return
    c => jfirst(inst)
    do while (associated(c))
       do k = a0 + 1, anns%n
          if (anns%a(k)%kind .ne. AK_DECLARED) cycle
          if (.not. associated(anns%a(k)%node, inst)) cycle
          if (str_eq(anns%a(k)%name, jname(c))) then
             n = n + 1
             exit
          end if
       end do
       c => jnext(c)
    end do

  end function declared_keys

  !> Report the failure of all alternatives of an anyOf or oneOf
  !! @details An alternative is ruled out when its messages show that it is
  !! meant for something else: a const or enum on a key such as "type" does
  !! not match, the value is of another type, or a not matches. Of the
  !! others, the one that declares most of the keys of the value is reported,
  !! as the one the user most likely meant. If all are ruled out, a single
  !! message lists what would have been accepted.
  !! @param inst The value
  !! @param path Where the value is
  !! @param br Messages of each alternative
  !! @param score Number of keys of the value each alternative declares
  subroutine summarise(inst, path, br, score)
    type(json_value), pointer, intent(in) :: inst
    character(len=*), intent(in) :: path
    type(msg_list_t), intent(in) :: br(:)
    integer, intent(in) :: score(:)
    logical, allocatable :: ruled_out(:)
    type(str_t), allocatable :: cand(:), allowed(:)
    integer, allocatable :: hits(:)
    type(json_value), pointer :: node
    character(len=:), allocatable :: got, dpath
    integer :: nb, k, j, i, best
    logical :: all_type

    nb = size(br)
    allocate(ruled_out(nb))
    do k = 1, nb
       ruled_out(k) = .false.
       do j = 1, br(k)%n
          if (rules_out(br(k)%m(j), path)) ruled_out(k) = .true.
       end do
    end do

    ! Alternatives that each miss one required key, e.g. value or file_name
    if (.not. any(ruled_out) .and. one_required_each(br, path)) then
       do k = 1, nb
          call add_str(allowed, br(k)%m(1)%allowed(1)%s)
       end do
       call report(SEV_ERROR, path, 'missing one of the keys ' // &
            value_list(allowed), inst, MK_REQUIRED)
       return
    end if

    if (.not. all(ruled_out)) then
       best = 0
       do k = 1, nb
          if (ruled_out(k)) cycle
          if (best .eq. 0) then
             best = k
          else if (score(k) .gt. score(best) .or. (score(k) .eq. &
               score(best) .and. br(k)%n .lt. br(best)%n)) then
             best = k
          end if
       end do
       call append_all(br(best))
       return
    end if

    ! All ruled out by a const or enum: find the key they disagree on
    do k = 1, nb
       do j = 1, br(k)%n
          if (is_discriminator(br(k)%m(j), path)) then
             call add_str(cand, br(k)%m(j)%path)
          end if
       end do
    end do
    if (allocated(cand)) then
       allocate(hits(size(cand)))
       hits = 0
       do i = 1, size(cand)
          do k = 1, nb
             do j = 1, br(k)%n
                if (is_discriminator(br(k)%m(j), path) .and. &
                     str_eq(br(k)%m(j)%path, cand(i)%s)) then
                   hits(i) = hits(i) + 1
                   exit
                end if
             end do
          end do
       end do
       dpath = cand(maxloc(hits, 1))%s
       node => null()
       do k = 1, nb
          do j = 1, br(k)%n
             if (is_discriminator(br(k)%m(j), path) .and. &
                  str_eq(br(k)%m(j)%path, dpath)) then
                do i = 1, size(br(k)%m(j)%allowed)
                   call add_str(allowed, br(k)%m(j)%allowed(i)%s)
                end do
                if (.not. associated(node)) then
                   node => br(k)%m(j)%node
                   got = br(k)%m(j)%got
                end if
             end if
          end do
       end do
       call report(SEV_ERROR, dpath, choice_text(got, allowed), node, &
            MK_ENUM, got, allowed)
       return
    end if

    ! All ruled out by the type of the value
    all_type = .true.
    do k = 1, nb
       j = type_message(br(k), path)
       if (j .eq. 0) then
          all_type = .false.
          exit
       end if
       do i = 1, size(br(k)%m(j)%allowed)
          call add_str(allowed, br(k)%m(j)%allowed(i)%s)
       end do
    end do
    if (all_type) then
       call report(SEV_ERROR, path, 'expected ' // type_phrase(allowed) // &
            ', got ' // describe(inst) // integer_hint(inst, allowed), inst, &
            MK_TYPE, describe(inst), allowed)
       return
    end if

    ! Otherwise the alternative with the fewest messages
    best = 1
    do k = 2, nb
       if (br(k)%n .lt. br(best)%n) best = k
    end do
    call append_all(br(best))

  end subroutine summarise

  !> Whether each alternative only misses one required key of a value
  !! @param br Messages of each alternative
  !! @param path Where the value is
  function one_required_each(br, path) result(each)
    type(msg_list_t), intent(in) :: br(:)
    character(len=*), intent(in) :: path
    logical :: each
    integer :: k

    each = size(br) .gt. 1
    do k = 1, size(br)
       if (br(k)%n .ne. 1) then
          each = .false.
       else if (br(k)%m(1)%kind .ne. MK_REQUIRED .or. &
            .not. str_eq(br(k)%m(1)%path, path) .or. &
            .not. allocated(br(k)%m(1)%allowed)) then
          each = .false.
       end if
       if (.not. each) return
    end do

  end function one_required_each

  !> Whether a message shows that an alternative is meant for something else
  !! @param m The message
  !! @param path Where the value of the anyOf or oneOf is
  function rules_out(m, path) result(ruled)
    type(msg_t), intent(in) :: m
    character(len=*), intent(in) :: path
    logical :: ruled

    select case (m%kind)
    case (MK_CONST, MK_ENUM)
       ruled = is_discriminator(m, path)
    case (MK_TYPE, MK_NOT)
       ruled = str_eq(m%path, path)
    case default
       ruled = .false.
    end select

  end function rules_out

  !> Whether a message is a const or enum on the value or one of its keys
  !! @param m The message
  !! @param path Where the value is
  function is_discriminator(m, path) result(disc)
    type(msg_t), intent(in) :: m
    character(len=*), intent(in) :: path
    logical :: disc

    disc = .false.
    if (m%kind .ne. MK_CONST .and. m%kind .ne. MK_ENUM) return
    if (.not. allocated(m%allowed)) return
    disc = str_eq(m%path, path) .or. str_eq(parent_path(m%path), path)

  end function is_discriminator

  !> Index of the type message for a value in a list, 0 if there is none
  !! @param list The messages
  !! @param path Where the value is
  function type_message(list, path) result(j)
    type(msg_list_t), intent(in) :: list
    character(len=*), intent(in) :: path
    integer :: j

    do j = 1, list%n
       if (list%m(j)%kind .eq. MK_TYPE .and. str_eq(list%m(j)%path, path) &
            .and. allocated(list%m(j)%allowed)) return
    end do
    j = 0

  end function type_message

  !> Add the messages of a list to msgs
  !! @param list The messages
  subroutine append_all(list)
    type(msg_list_t), intent(in) :: list
    integer :: j

    do j = 1, list%n
       call msgs%add(list%m(j))
    end do

  end subroutine append_all

  !> Set how the constants named by a matching constRef are looked up
  !! @details The kind follows from the other alternatives, e.g. the number
  !! in numberOrRef, which says what Neko reads the value as.
  !! @param inst The value
  !! @param arr The alternatives
  !! @param kpass The alternative that matched
  !! @param doc Index of the schema document arr is part of
  !! @param a0 Number of annotations before the alternatives
  subroutine set_const_kinds(inst, arr, kpass, doc, a0)
    type(json_value), pointer, intent(in) :: inst, arr
    integer, intent(in) :: kpass, doc, a0
    type(json_value), pointer :: e
    integer :: k, j, ck, nmin, nmax

    ck = CK_ANY
    nmin = -1
    nmax = -1
    e => jfirst(arr)
    k = 0
    do while (associated(e) .and. ck .eq. CK_ANY)
       k = k + 1
       if (k .ne. kpass) call schema_kind(e, doc, ck, nmin, nmax)
       e => jnext(e)
    end do
    if (ck .eq. CK_ANY) return

    do j = a0 + 1, anns%n
       if (anns%a(j)%kind .ne. AK_CONSTREF) cycle
       if (anns%a(j)%ck .ne. CK_ANY) cycle
       if (.not. associated(anns%a(j)%node, inst)) cycle
       anns%a(j)%ck = ck
       anns%a(j)%min_items = nmin
       anns%a(j)%max_items = nmax
    end do

  end subroutine set_const_kinds

  !> The kind of value a schema accepts, as a kind of constant
  !! @param s The schema
  !! @param doc Index of the schema document s is part of
  !! @param ck The kind, CK_ANY if it is none of the constant kinds
  !! @param nmin Minimum number of entries of an array, -1 if not limited
  !! @param nmax Maximum number of entries of an array, -1 if not limited
  recursive subroutine schema_kind(s, doc, ck, nmin, nmax)
    type(json_value), pointer, intent(in) :: s
    integer, intent(in) :: doc
    integer, intent(out) :: ck, nmin, nmax
    type(json_value), pointer :: t, ref, items, p
    integer :: tdoc, idoc, jmin, jmax, ick

    ck = CK_ANY
    nmin = -1
    nmax = -1
    t => jchild(s, 'type')
    ref => jchild(s, '$ref')
    if (.not. associated(t)) then
       if (associated(ref)) then
          call resolve(jstr(ref), doc, p, tdoc)
          call schema_kind(p, tdoc, ck, nmin, nmax)
       end if
       return
    end if

    select case (jstr(t))
    case ('integer')
       ck = CK_INTEGER
    case ('number')
       ck = CK_NUMBER
       if (jbool(jchild(s, 'x-neko-real'))) ck = CK_REAL
    case ('array')
       ck = CK_REAL_ARRAY
       items => jchild(s, 'items')
       if (associated(items)) then
          idoc = doc
          call schema_kind(items, idoc, ick, jmin, jmax)
          if (ick .eq. CK_INTEGER) ck = CK_INTEGER_ARRAY
       end if
       p => jchild(s, 'minItems')
       if (associated(p)) nmin = jint(p)
       p => jchild(s, 'maxItems')
       if (associated(p)) nmax = jint(p)
    end select

  end subroutine schema_kind

  !> if, then and else
  !! @param inst The value
  !! @param path Where the value is in the case file
  !! @param sch The schema
  !! @param doc Index of the schema document sch is part of
  recursive subroutine kw_conditional(inst, path, sch, doc)
    type(json_value), pointer, intent(in) :: inst
    character(len=*), intent(in) :: path
    type(json_value), pointer, intent(in) :: sch
    integer, intent(in) :: doc
    type(json_value), pointer :: p
    integer :: n0, a0
    logical :: passed

    p => jchild(sch, 'if')
    if (.not. associated(p)) return

    n0 = msgs%n
    a0 = anns%n
    call validate(inst, path, p, doc)
    passed = msgs%n .eq. n0
    call msgs%truncate(n0)
    call anns%truncate(a0)

    if (passed) then
       p => jchild(sch, 'then')
    else
       p => jchild(sch, 'else')
    end if
    if (associated(p)) call validate(inst, path, p, doc)

  end subroutine kw_conditional

  !> Append an annotation to a list
  !! @param ann The annotation
  subroutine ann_list_add(this, ann)
    class(ann_list_t), intent(inout) :: this
    type(ann_t), intent(in) :: ann
    type(ann_t), allocatable :: tmp(:)

    if (.not. allocated(this%a)) allocate(this%a(256))
    if (this%n .eq. size(this%a)) then
       allocate(tmp(2 * size(this%a)))
       tmp(1:this%n) = this%a(1:this%n)
       call move_alloc(tmp, this%a)
    end if
    this%n = this%n + 1
    this%a(this%n) = ann

  end subroutine ann_list_add

  !> Drop the annotations after the first n
  !! @param n Number of annotations to keep
  subroutine ann_list_truncate(this, n)
    class(ann_list_t), intent(inout) :: this
    integer, intent(in) :: n

    this%n = min(this%n, n)

  end subroutine ann_list_truncate

end module checkneko_schema

!> Checks of a case file beyond its schema
!! @details These look at what the schema cannot see: the constants a value
!! names, the files the case reads, the labeled zones of the mesh against the
!! boundary conditions, the point zones a value names, and keys that are
!! misspelled or given twice. The zone checks follow the ones the solvers do
!! when they set up their boundary conditions, so that what makes Neko stop
!! is an error and what it accepts is at most a warning.
module checkneko_checks
  use, intrinsic :: iso_c_binding, only : c_loc, c_intptr_t
  use json_module, only : json_value, json_object, json_array, json_integer, &
       json_real, json_string
  use checkneko_util
  use checkneko_schema, only : anns, AK_DECLARED, AK_CONSTREF, CK_ANY, &
       CK_INTEGER, CK_REAL, CK_NUMBER, CK_INTEGER_ARRAY, CK_REAL_ARRAY
  implicit none
  private

  !> Number of labeled zones Neko supports (NEKO_MSH_MAX_ZLBLS)
  integer, public, parameter :: MAX_ZONES = 20

  !> nmsh zone types
  integer, parameter :: PERIODIC_ZONE = 5, LABELED_ZONE = 7

  !> What is known about the mesh
  type, public :: mesh_info_t
     !> Whether the zones were read
     logical :: read = .false.
     integer :: nelv = 0
     integer :: gdim = 0
     integer :: nperiodic = 0
     !> Number of facets with each label
     integer :: zone_size(MAX_ZONES) = 0
     !> Number of facets with a label outside [1, MAX_ZONES]
     integer :: nbad = 0
  end type mesh_info_t

  !> A constant from case.constants
  type :: const_t
     character(len=:), allocatable :: name
     !> CK_INTEGER, CK_REAL or CK_REAL_ARRAY
     integer :: ck = CK_ANY
     real(kind=dp), allocatable :: vals(:)
  end type const_t

  !> Line numbers of the values of the case file
  type, public :: line_map_t
     integer(kind=c_intptr_t), allocatable :: addr(:)
     integer, allocatable :: line(:)
     integer :: n = 0
  end type line_map_t

  !> The boundary conditions of a field on each zone, for the zone table
  type :: field_bcs_t
     !> Heading of the column of the field
     character(len=:), allocatable :: name
     !> Types of the boundary conditions on each zone, unallocated if none
     type(str_t) :: types(MAX_ZONES)
  end type field_bcs_t

  type(const_t), allocatable :: consts(:)
  integer :: nconsts = 0
  type(str_t), allocatable :: zone_names(:)
  type(field_bcs_t), allocatable :: fields(:)
  integer :: nfields = 0

  !> Keys whose value is a file the case reads
  character(len=*), parameter :: INPUT_KEYS(7) = [character(len=22) :: &
       'mesh_file', 'restart_file', 'restart_mesh_file', 'file_name', &
       'mesh_file_name', 'points_file', 'base_shape_import_file']

  public :: check_constants, check_input_files, check_mesh, &
       check_boundary_conditions, check_references, check_keys, &
       print_zone_table, line_map_build, line_of

contains

  !> Check case.constants and the values that name a constant
  !! @param root The case file
  !! @param summary Short summary for the progress line
  subroutine check_constants(root, summary)
    type(json_value), pointer, intent(in) :: root
    character(len=:), allocatable, intent(out) :: summary
    type(json_value), pointer :: arr, e, nm, val, v
    character(len=:), allocatable :: name, path
    integer :: i, k, nref, nbad, n0

    nconsts = 0
    if (allocated(consts)) deallocate(consts)
    arr => node_at(root, 'case.constants')
    allocate(consts(max(1, jnchildren(arr))))
    n0 = msgs%n

    e => jfirst(arr)
    i = 0
    do while (associated(e))
       i = i + 1
       path = path_index('case.constants', i)
       nm => jchild(e, 'name')
       val => jchild(e, 'value')
       if (jtype(nm) .eq. json_string .and. associated(val)) then
          name = jstr(nm)
          k = find_const(name)
          if (k .gt. 0) then
             if (same_class(consts(k)%ck, val)) then
                call report(SEV_ERROR, path_key(path, 'name'), &
                     'the constant ' // quoted(name) // ' is defined ' // &
                     'more than once', nm)
             end if
          end if
          nconsts = nconsts + 1
          consts(nconsts)%name = name
          select case (jtype(val))
          case (json_integer)
             consts(nconsts)%ck = CK_INTEGER
             consts(nconsts)%vals = [jreal(val)]
          case (json_real)
             consts(nconsts)%ck = CK_REAL
             consts(nconsts)%vals = [jreal(val)]
          case (json_array)
             consts(nconsts)%ck = CK_REAL_ARRAY
             allocate(consts(nconsts)%vals(0))
             v => jfirst(val)
             do while (associated(v))
                consts(nconsts)%vals = [consts(nconsts)%vals, jreal(v)]
                v => jnext(v)
             end do
          end select
       end if
       e => jnext(e)
    end do

    ! A value can be validated by several schemas, e.g. polynomial_order by
    ! the one of numerics and the one for the pnpn scheme
    nref = 0
    do k = 1, anns%n
       if (anns%a(k)%kind .ne. AK_CONSTREF) cycle
       if (.not. seen_before(k, .false.)) nref = nref + 1
       if (seen_before(k, .true.)) cycle
       call check_const_ref(anns%a(k)%node, anns%a(k)%path, anns%a(k)%ck, &
            anns%a(k)%min_items, anns%a(k)%max_items)
    end do

    nbad = msgs%n - n0
    summary = int_str(nconsts) // ' defined, ' // int_str(nref) // &
         ' referenced'
    if (nbad .gt. 0) summary = summary // ', ' // problems(nbad)

  contains

    !> Whether an earlier annotation is for the same value
    !! @param k Index of the annotation
    !! @param same_kind Whether it must also be for the same kind
    function seen_before(k, same_kind) result(seen)
      integer, intent(in) :: k
      logical, intent(in) :: same_kind
      logical :: seen
      integer :: j

      seen = .false.
      do j = 1, k - 1
         if (anns%a(j)%kind .ne. AK_CONSTREF) cycle
         if (.not. associated(anns%a(j)%node, anns%a(k)%node)) cycle
         if (.not. same_kind .or. anns%a(j)%ck .eq. anns%a(k)%ck) seen = .true.
      end do
    end function seen_before

  end subroutine check_constants

  !> Whether a constant of a kind is of the same kind as a value
  !! @param ck Kind of the constant
  !! @param val The value
  function same_class(ck, val) result(same)
    integer, intent(in) :: ck
    type(json_value), pointer, intent(in) :: val
    logical :: same

    select case (jtype(val))
    case (json_integer)
       same = ck .eq. CK_INTEGER
    case (json_real)
       same = ck .eq. CK_REAL
    case (json_array)
       same = ck .eq. CK_REAL_ARRAY
    case default
       same = .false.
    end select

  end function same_class

  !> Index of a constant, 0 if there is none with the name
  !! @param name Name of the constant
  !! @param ck Kind of the constant, any kind if absent
  function find_const(name, ck) result(k)
    character(len=*), intent(in) :: name
    integer, intent(in), optional :: ck
    integer :: k

    do k = 1, nconsts
       if (.not. str_eq(consts(k)%name, name)) cycle
       if (present(ck)) then
          if (consts(k)%ck .ne. ck) cycle
       end if
       return
    end do
    k = 0

  end function find_const

  !> Check a value that names a constant
  !! @param node The value
  !! @param path Where the value is
  !! @param ck How Neko looks the constant up
  !! @param nmin Minimum number of entries of an array, -1 if not limited
  !! @param nmax Maximum number of entries of an array, -1 if not limited
  subroutine check_const_ref(node, path, ck, nmin, nmax)
    type(json_value), pointer, intent(in) :: node
    character(len=*), intent(in) :: path
    integer, intent(in) :: ck, nmin, nmax
    character(len=:), allocatable :: name
    integer :: k, n

    name = jstr(node)
    k = find_const(name)
    if (k .eq. 0) then
       call report(SEV_ERROR, path, 'the constant ' // quoted(name) // &
            ' is not defined in case.constants' // const_suggestion(name), &
            node)
       return
    end if

    select case (ck)
    case (CK_INTEGER)
       if (find_const(name, CK_INTEGER) .eq. 0) then
          call report(SEV_ERROR, path, 'the constant ' // quoted(name) // &
               ' must be an integer here, but it is ' // &
               kind_text(consts(k)%ck), node)
       end if
    case (CK_REAL)
       if (find_const(name, CK_REAL) .eq. 0) then
          call report(SEV_ERROR, path, 'the constant ' // quoted(name) // &
               ' must be a real number here, but it is ' // &
               kind_text(consts(k)%ck) // real_hint(k), node)
       end if
    case (CK_NUMBER)
       if (find_const(name, CK_INTEGER) .eq. 0 .and. &
            find_const(name, CK_REAL) .eq. 0) then
          call report(SEV_ERROR, path, 'the constant ' // quoted(name) // &
               ' must be a number here, but it is ' // &
               kind_text(consts(k)%ck), node)
       end if
    case (CK_INTEGER_ARRAY, CK_REAL_ARRAY)
       k = find_const(name, CK_REAL_ARRAY)
       if (k .eq. 0) then
          call report(SEV_ERROR, path, 'the constant ' // quoted(name) // &
               ' must be an array here, but it is ' // &
               kind_text(consts(find_const(name))%ck), node)
          return
       end if
       n = size(consts(k)%vals)
       if (nmin .ge. 0 .and. nmin .eq. nmax .and. n .ne. nmin) then
          call report(SEV_ERROR, path, 'the constant ' // quoted(name) // &
               ' must have ' // int_str(nmin) // ' entries here, it has ' // &
               int_str(n), node)
       else if (nmin .ge. 0 .and. n .lt. nmin) then
          call report(SEV_ERROR, path, 'the constant ' // quoted(name) // &
               ' must have at least ' // int_str(nmin) // ' entries here, ' &
               // 'it has ' // int_str(n), node)
       else if (nmax .ge. 0 .and. n .gt. nmax) then
          call report(SEV_ERROR, path, 'the constant ' // quoted(name) // &
               ' must have at most ' // int_str(nmax) // ' entries here, ' &
               // 'it has ' // int_str(n), node)
       end if
       if (ck .eq. CK_INTEGER_ARRAY) then
          if (.not. all(is_whole(consts(k)%vals))) then
             call report(SEV_ERROR, path, 'the constant ' // quoted(name) // &
                  ' must hold whole numbers here, such as 1.0 and 2.0', node)
          end if
       end if
    end select

  end subroutine check_const_ref

  !> Kind of a constant as text
  !! @param ck The kind
  function kind_text(ck) result(text)
    integer, intent(in) :: ck
    character(len=:), allocatable :: text

    select case (ck)
    case (CK_INTEGER)
       text = 'an integer'
    case (CK_REAL)
       text = 'a real number'
    case (CK_REAL_ARRAY)
       text = 'an array'
    case default
       text = 'of an unsupported type'
    end select

  end function kind_text

  !> A hint to write an integer constant as a real
  !! @param k Index of the constant
  function real_hint(k) result(text)
    integer, intent(in) :: k
    character(len=:), allocatable :: text

    text = ''
    if (consts(k)%ck .eq. CK_INTEGER) then
       text = ' (write its value as ' // int_str(nint(consts(k)%vals(1))) &
            // '.0 in case.constants)'
    end if

  end function real_hint

  !> A suggestion for a misspelled constant name
  !! @param name The name
  function const_suggestion(name) result(text)
    character(len=*), intent(in) :: name
    character(len=:), allocatable :: text
    integer :: k, d, dbest

    text = ''
    dbest = huge(1)
    do k = 1, nconsts
       d = edit_distance(name, consts(k)%name)
       if (d .lt. dbest .and. d .le. typo_limit(len(name)) .and. &
            d .lt. len(name)) then
          dbest = d
          text = ' (did you mean ' // quoted(consts(k)%name) // '?)'
       end if
    end do

  end function const_suggestion

  !> Check that the files the case reads exist
  !! @param root The case file
  !! @param run_dir Directory Neko runs in, relative names are relative to it
  !! @param summary Short summary for the progress line
  subroutine check_input_files(root, run_dir, summary)
    type(json_value), pointer, intent(in) :: root
    character(len=*), intent(in) :: run_dir
    character(len=:), allocatable, intent(out) :: summary
    integer :: nfound, nmissing

    nfound = 0
    nmissing = 0
    call walk(root, '')
    if (nfound + nmissing .eq. 0) then
       summary = 'none'
    else
       summary = int_str(nfound) // ' found'
       if (nmissing .gt. 0) summary = summary // ', ' // &
            int_str(nmissing) // ' missing'
    end if

  contains

    !> Check the file names among the values below a value
    !! @param p The value
    !! @param path Where the value is
    recursive subroutine walk(p, path)
      type(json_value), pointer, intent(in) :: p
      character(len=*), intent(in) :: path
      type(json_value), pointer :: c
      character(len=:), allocatable :: key
      integer :: i

      c => jfirst(p)
      i = 0
      do while (associated(c))
         i = i + 1
         if (jtype(p) .eq. json_object) then
            key = jname(c)
            if (jtype(c) .eq. json_string) then
               if (any(INPUT_KEYS .eq. key)) then
                  call check_file(c, path_key(path, key))
               else if (key .eq. 'name' .and. &
                    jstr(jchild(p, 'type')) .eq. 'boundary_mesh') then
                  call check_file(c, path_key(path, key))
               end if
            end if
            call walk(c, path_key(path, key))
         else
            call walk(c, path_index(path, i))
         end if
         c => jnext(c)
      end do
    end subroutine walk

    !> Check that the file a value names exists
    !! @param c The value
    !! @param path Where the value is
    subroutine check_file(c, path)
      type(json_value), pointer, intent(in) :: c
      character(len=*), intent(in) :: path
      character(len=:), allocatable :: fname, full

      fname = jstr(c)
      if (len(fname) .eq. 0) return
      full = join_path(run_dir, fname)
      if (file_exists(full)) then
         nfound = nfound + 1
      else
         nmissing = nmissing + 1
         if (str_eq(full, fname)) then
            call report(SEV_ERROR, path, 'the file ' // quoted(fname) // &
                 ' does not exist', c)
         else
            call report(SEV_ERROR, path, 'the file ' // quoted(fname) // &
                 ' does not exist (looked for ' // full // ')', c)
         end if
      end if
    end subroutine check_file

  end subroutine check_input_files

  !> Read the labeled zones of the mesh
  !! @param root The case file
  !! @param run_dir Directory Neko runs in
  !! @param mesh What is known about the mesh
  !! @param summary Short summary for the progress line
  subroutine check_mesh(root, run_dir, mesh, summary)
    type(json_value), pointer, intent(in) :: root
    character(len=*), intent(in) :: run_dir
    type(mesh_info_t), intent(out) :: mesh
    character(len=:), allocatable, intent(out) :: summary
    type(json_value), pointer :: mf
    character(len=:), allocatable :: fname, full, err
    integer :: z

    mf => node_at(root, 'case.mesh_file')
    if (jtype(mf) .ne. json_string) then
       summary = 'no mesh_file'
       return
    end if
    fname = jstr(mf)
    full = join_path(run_dir, fname)
    if (file_suffix(fname) .ne. 'nmsh') then
       summary = 'zones are only checked for .nmsh meshes'
       return
    end if
    if (.not. file_exists(full)) then
       summary = 'not found'
       return
    end if

    call read_nmsh_zones(full, mesh, err)
    if (len(err) .gt. 0) then
       call report(SEV_ERROR, 'case.mesh_file', 'cannot read the mesh ' // &
            quoted(fname) // ': ' // err, mf)
       summary = 'unreadable'
       mesh%read = .false.
       return
    end if

    summary = int_str(mesh%nelv) // ' elements, '
    if (mesh%gdim .eq. 2) summary = summary // '2D, '
    if (any(mesh%zone_size .gt. 0)) then
       summary = summary // 'zones'
       do z = 1, MAX_ZONES
          if (mesh%zone_size(z) .gt. 0) summary = summary // ' ' // int_str(z)
       end do
    else
       summary = summary // 'no labeled zones'
    end if
    if (mesh%nperiodic .gt. 0) summary = summary // ', periodic'
    if (mesh%nbad .gt. 0) then
       call report(SEV_ERROR, 'case.mesh_file', 'the mesh ' // &
            quoted(fname) // ' has ' // int_str(mesh%nbad) // ' facets ' // &
            'with a zone label outside 1 to ' // int_str(MAX_ZONES), mf)
    end if

  end subroutine check_mesh

  !> Read the header and the zones of an nmsh file
  !! @details The layout is that of Neko's nmsh_file: the number of elements
  !! and the dimension, the elements, the zones and the curved elements, with
  !! 4 byte integers and 8 byte reals in native byte order without padding.
  !! An element is its index and 2**gdim vertices of an index and three
  !! coordinates. A zone is 9 integers, the last one its type; the label of
  !! a labeled zone is the fourth. The elements are skipped, not read.
  !! @param fname Name of the nmsh file
  !! @param mesh What is read
  !! @param err Why the file could not be read, empty if it could
  subroutine read_nmsh_zones(fname, mesh, err)
    character(len=*), intent(in) :: fname
    type(mesh_info_t), intent(inout) :: mesh
    character(len=:), allocatable, intent(out) :: err
    integer, parameter :: CHUNK = 65536
    integer(kind=i4) :: nelv, gdim, nzones, ncurves
    integer(kind=i4), allocatable :: zones(:,:)
    integer(kind=i8) :: fsize, need, pos
    integer :: unit, ios, nread, m, k, label

    err = ''
    open(newunit = unit, file = fname, access = 'stream', &
         form = 'unformatted', status = 'old', action = 'read', iostat = ios)
    if (ios .ne. 0) then
       err = 'cannot open it'
       return
    end if
    inquire(unit = unit, size = fsize)

    if (fsize .lt. 8) then
       err = 'the file is too short to be an nmsh file'
       close(unit)
       return
    end if
    read(unit, iostat = ios) nelv, gdim
    if (ios .ne. 0 .or. nelv .le. 0 .or. (gdim .ne. 2 .and. gdim .ne. 3)) then
       err = 'it does not start like an nmsh file'
       close(unit)
       return
    end if

    ! Element index, then 2**gdim vertices of an index and three reals
    pos = 9_i8 + int(nelv, i8) * int(4 + 2**gdim * 28, i8)
    need = pos + 3
    if (fsize .lt. need) then
       err = truncated(fsize, need)
       close(unit)
       return
    end if
    read(unit, pos = pos, iostat = ios) nzones
    if (ios .ne. 0 .or. nzones .lt. 0) then
       err = 'invalid number of zones'
       close(unit)
       return
    end if
    need = need + 36_i8 * int(nzones, i8) + 4_i8
    if (fsize .lt. need) then
       err = truncated(fsize, need)
       close(unit)
       return
    end if

    allocate(zones(9, CHUNK))
    nread = 0
    do while (nread .lt. nzones)
       m = min(CHUNK, nzones - nread)
       read(unit, iostat = ios) zones(:, 1:m)
       if (ios .ne. 0) then
          err = 'cannot read the zones'
          close(unit)
          return
       end if
       do k = 1, m
          select case (zones(9, k))
          case (LABELED_ZONE)
             label = zones(4, k)
             if (label .ge. 1 .and. label .le. MAX_ZONES) then
                mesh%zone_size(label) = mesh%zone_size(label) + 1
             else
                mesh%nbad = mesh%nbad + 1
             end if
          case (PERIODIC_ZONE)
             mesh%nperiodic = mesh%nperiodic + 1
          end select
       end do
       nread = nread + m
    end do

    read(unit, iostat = ios) ncurves
    if (ios .ne. 0 .or. ncurves .lt. 0) then
       err = 'invalid number of curved elements'
       close(unit)
       return
    end if
    ! Element index, 5x12 reals and 12 curve types
    need = need + int(ncurves, i8) * 532_i8
    if (fsize .lt. need) then
       err = truncated(fsize, need)
       close(unit)
       return
    end if
    close(unit)

    mesh%read = .true.
    mesh%nelv = nelv
    mesh%gdim = gdim

  end subroutine read_nmsh_zones

  !> Message for a truncated nmsh file
  !! @param fsize Size of the file
  !! @param need Size the header says it has at least
  function truncated(fsize, need) result(text)
    integer(kind=i8), intent(in) :: fsize, need
    character(len=:), allocatable :: text
    character(len=24) :: b1, b2

    write(b1, '(I0)') fsize
    write(b2, '(I0)') need
    text = 'the file is truncated, it has ' // trim(b1) // &
         ' bytes but its header needs at least ' // trim(b2)

  end function truncated

  !> Check the boundary conditions of the fluid and the scalars
  !! @param root The case file
  !! @param mesh What is known about the mesh
  !! @param summary Short summary for the progress line
  subroutine check_boundary_conditions(root, mesh, summary)
    type(json_value), pointer, intent(in) :: root
    type(mesh_info_t), intent(in) :: mesh
    character(len=:), allocatable, intent(out) :: summary
    type(json_value), pointer :: fluid, scalar, scalars, e, nm
    type(str_t), allocatable :: names(:)
    character(len=:), allocatable :: scheme, name
    integer :: i, k, nscalars, n0

    n0 = msgs%n
    scalars => node_at(root, 'case.scalars')
    nscalars = jnchildren(scalars)
    if (allocated(fields)) deallocate(fields)
    allocate(fields(2 + nscalars))
    nfields = 0

    fluid => node_at(root, 'case.fluid')
    if (jtype(fluid) .eq. json_object) then
       scheme = jstr(jchild(fluid, 'scheme'))
       ! The compressible solver takes several entries for one zone, e.g. the
       ! velocity, density and pressure of an inflow, and leaves zones without
       ! an entry free, e.g. a supersonic outflow
       call add_field('Fluid')
       call check_bc_list(fluid, 'case.fluid', 'fluid', &
            scheme .ne. 'compressible', mesh, nfields)
    end if

    ! Scalars are named as by Neko, s for a single one and s_1, s_2, ...
    ! for several, unless they have a name
    scalar => node_at(root, 'case.scalar')
    if (jtype(scalar) .eq. json_object) then
       name = jstr(jchild(scalar, 'name'))
       if (len_trim(name) .eq. 0) name = 's'
       call add_field('Scalar ' // trim(name))
       call check_bc_list(scalar, 'case.scalar', 'scalar', .true., mesh, &
            nfields)
    end if

    allocate(names(nscalars))
    e => jfirst(scalars)
    i = 0
    do while (associated(e))
       i = i + 1
       nm => jchild(e, 'name')
       names(i)%s = trim(jstr(nm))
       if (len(names(i)%s) .eq. 0) then
          if (nscalars .eq. 1) then
             names(i)%s = 's'
          else
             names(i)%s = 's_' // int_str(i)
          end if
       end if
       ! Neko stops on two scalars with the same name
       do k = 1, i - 1
          if (str_eq(names(k)%s, names(i)%s)) then
             if (.not. associated(nm)) nm => e
             call report(SEV_ERROR, path_index('case.scalars', i), &
                  'the scalar name ' // quoted(names(i)%s) // ' is also ' // &
                  'the name of case.scalars[' // int_str(k) // '], the ' // &
                  'names must differ', nm)
             exit
          end if
       end do
       if (jtype(e) .eq. json_object) then
          call add_field('Scalar ' // names(i)%s)
          call check_bc_list(e, path_index('case.scalars', i), 'scalar', &
               .true., mesh, nfields)
       end if
       e => jnext(e)
    end do

    if (.not. mesh%read) then
       summary = 'zone indices checked, not against the mesh'
    else if (msgs%n .eq. n0) then
       summary = 'ok'
    else
       summary = problems(msgs%n - n0)
    end if

  contains

    !> Add a column to the zone table
    !! @param heading Heading of the column
    subroutine add_field(heading)
      character(len=*), intent(in) :: heading

      nfields = nfields + 1
      fields(nfields)%name = heading
    end subroutine add_field

  end subroutine check_boundary_conditions

  !> Check the boundary conditions of a field against the mesh
  !! @details Neko stops if the zone of an entry is not in the mesh, if two
  !! entries are for the same zone, or if a zone of the mesh has no entry.
  !! The compressible solver accepts all three, but a zone that is not in the
  !! mesh is still likely a mistake.
  !! @param field The object of the field, with its boundary_conditions
  !! @param path Where the object is
  !! @param what Name of the field for the messages
  !! @param strict Whether one entry per zone of the mesh is required
  !! @param mesh What is known about the mesh
  !! @param f Column of the field in the zone table
  subroutine check_bc_list(field, path, what, strict, mesh, f)
    type(json_value), pointer, intent(in) :: field
    character(len=*), intent(in) :: path, what
    logical, intent(in) :: strict
    type(mesh_info_t), intent(in) :: mesh
    integer, intent(in) :: f
    type(json_value), pointer :: bcs, bc, zi
    integer, allocatable :: zones(:)
    integer :: owner(MAX_ZONES)
    character(len=:), allocatable :: bc_path, list, label
    integer :: i, j, z, sev
    logical :: ok

    sev = SEV_WARNING
    if (strict) sev = SEV_ERROR

    bcs => jchild(field, 'boundary_conditions')
    if (.not. associated(bcs)) then
       if (mesh%read .and. any(mesh%zone_size .gt. 0)) then
          call report(SEV_ERROR, path, 'there is no boundary_conditions ' // &
               'entry, but the mesh has the labeled zones' // &
               zone_list(mesh%zone_size .gt. 0), field)
       end if
       return
    end if
    if (jtype(bcs) .ne. json_array) return

    owner = 0
    bc => jfirst(bcs)
    i = 0
    do while (associated(bc))
       i = i + 1
       bc_path = path_index(path_key(path, 'boundary_conditions'), i)
       label = jstr(jchild(bc, 'type'))
       if (len(label) .eq. 0) label = '?'
       if (len(jstr(jchild(bc, 'name'))) .gt. 0) then
          label = label // ' (' // jstr(jchild(bc, 'name')) // ')'
       end if
       zi => jchild(bc, 'zone_indices')
       call int_list(zi, zones, ok)
       if (ok) then
          do j = 1, size(zones)
             z = zones(j)
             if (z .lt. 1 .or. z .gt. MAX_ZONES) then
                call report(SEV_ERROR, path_key(bc_path, 'zone_indices'), &
                     'the zone index ' // int_str(z) // ' is out of ' // &
                     'range, zone indices are 1 to ' // int_str(MAX_ZONES), zi)
                cycle
             end if
             if (allocated(fields(f)%types(z)%s)) then
                fields(f)%types(z)%s = fields(f)%types(z)%s // ', ' // label
             else
                fields(f)%types(z)%s = label
             end if
             if (mesh%read) then
                if (mesh%zone_size(z) .eq. 0) then
                   call report(sev, path_key(bc_path, 'zone_indices'), &
                        'there is no zone ' // int_str(z) // ' in the mesh', &
                        zi)
                end if
             end if
             if (owner(z) .ne. 0 .and. strict) then
                call report(SEV_ERROR, path_key(bc_path, 'zone_indices'), &
                     'zone ' // int_str(z) // ' already has a ' // what // &
                     ' boundary condition, in boundary_conditions[' // &
                     int_str(owner(z)) // ']', zi)
             else if (owner(z) .eq. 0) then
                owner(z) = i
             end if
          end do
       end if
       bc => jnext(bc)
    end do

    if (mesh%read .and. strict) then
       if (any(mesh%zone_size .gt. 0 .and. owner .eq. 0)) then
          list = zone_list(mesh%zone_size .gt. 0 .and. owner .eq. 0)
          if (count(mesh%zone_size .gt. 0 .and. owner .eq. 0) .eq. 1) then
             call report(sev, path_key(path, 'boundary_conditions'), &
                  'the mesh zone' // list // ' has no ' // what // &
                  ' boundary condition', bcs)
          else
             call report(sev, path_key(path, 'boundary_conditions'), &
                  'the mesh zones' // list // ' have no ' // what // &
                  ' boundary condition', bcs)
          end if
       end if
    end if

  end subroutine check_bc_list

  !> Print the labeled zones with their facets and the boundary conditions
  !! of the fluid and the scalars on them
  !! @details The zones are those of the mesh and those that boundary
  !! conditions refer to. The facets of a zone that is not in the mesh are
  !! shown as -, and as ? if the mesh was not read. A zone without a boundary
  !! condition of a field shows - for it.
  !! @param mesh What is known about the mesh
  subroutine print_zone_table(mesh)
    type(mesh_info_t), intent(in) :: mesh
    type(str_t), allocatable :: cells(:,:)
    integer, allocatable :: width(:)
    character(len=:), allocatable :: line
    logical :: used(MAX_ZONES)
    integer :: z, f, r, c, nrow, ncol

    used = .false.
    if (mesh%read) used = mesh%zone_size .gt. 0
    do f = 1, nfields
       do z = 1, MAX_ZONES
          if (allocated(fields(f)%types(z)%s)) used(z) = .true.
       end do
    end do
    nrow = count(used)
    if (nrow .eq. 0) return

    ncol = 2 + nfields
    allocate(cells(ncol, 0:nrow), width(ncol))
    cells(1, 0)%s = 'Zone'
    cells(2, 0)%s = 'Facets'
    do f = 1, nfields
       cells(2 + f, 0)%s = fields(f)%name
    end do
    r = 0
    do z = 1, MAX_ZONES
       if (.not. used(z)) cycle
       r = r + 1
       cells(1, r)%s = int_str(z)
       if (.not. mesh%read) then
          cells(2, r)%s = '?'
       else if (mesh%zone_size(z) .eq. 0) then
          cells(2, r)%s = '-'
       else
          cells(2, r)%s = int_str(mesh%zone_size(z))
       end if
       do f = 1, nfields
          if (allocated(fields(f)%types(z)%s)) then
             cells(2 + f, r)%s = fields(f)%types(z)%s
          else
             cells(2 + f, r)%s = '-'
          end if
       end do
    end do

    do c = 1, ncol
       width(c) = 0
       do r = 0, nrow
          width(c) = max(width(c), len(cells(c, r)%s))
       end do
    end do

    ! The zone and the facets right aligned, the boundary conditions left
    write(*, '(A)') ''
    do r = 0, nrow
       line = ' ' // repeat(' ', width(1) - len(cells(1, r)%s)) // &
            cells(1, r)%s
       line = line // '  ' // repeat(' ', width(2) - len(cells(2, r)%s)) // &
            cells(2, r)%s
       do c = 3, ncol
          line = line // '  ' // cells(c, r)%s // &
               repeat(' ', width(c) - len(cells(c, r)%s))
       end do
       write(*, '(A)') trim(line)
    end do

  end subroutine print_zone_table

  !> Zone indices as a list, e.g. " 1, 2 and 5"
  !! @param mask Which zones to list
  function zone_list(mask) result(text)
    logical, intent(in) :: mask(MAX_ZONES)
    character(len=:), allocatable :: text
    integer :: z, n, k

    text = ''
    n = count(mask)
    k = 0
    do z = 1, MAX_ZONES
       if (.not. mask(z)) cycle
       k = k + 1
       if (k .eq. 1) then
          text = ' ' // int_str(z)
       else if (k .eq. n) then
          text = text // ' and ' // int_str(z)
       else
          text = text // ', ' // int_str(z)
       end if
    end do

  end function zone_list

  !> Integers from an integer, an array of integers or a constant
  !! @param p The value
  !! @param vals The integers
  !! @param ok Whether the value gave integers
  subroutine int_list(p, vals, ok)
    type(json_value), pointer, intent(in) :: p
    integer, allocatable, intent(out) :: vals(:)
    logical, intent(out) :: ok
    type(json_value), pointer :: e
    integer :: k

    ok = .false.
    allocate(vals(0))
    select case (jtype(p))
    case (json_integer)
       vals = [jint(p)]
       ok = .true.
    case (json_array)
       e => jfirst(p)
       do while (associated(e))
          if (jtype(e) .ne. json_integer) return
          vals = [vals, jint(e)]
          e => jnext(e)
       end do
       ok = .true.
    case (json_string)
       k = find_const(jstr(p), CK_REAL_ARRAY)
       if (k .eq. 0) return
       if (.not. all(is_whole(consts(k)%vals))) return
       vals = nint(consts(k)%vals)
       ok = .true.
    end select

  end subroutine int_list

  !> Check the zones and point zones that values outside the boundary
  !! conditions refer to
  !! @param root The case file
  !! @param mesh What is known about the mesh
  !! @param summary Short summary for the progress line
  subroutine check_references(root, mesh, summary)
    type(json_value), pointer, intent(in) :: root
    type(mesh_info_t), intent(in) :: mesh
    character(len=:), allocatable, intent(out) :: summary
    type(json_value), pointer :: pz, e
    integer :: nzone_refs, npz_refs, n0

    n0 = msgs%n
    if (allocated(zone_names)) deallocate(zone_names)
    allocate(zone_names(0))
    pz => node_at(root, 'case.point_zones')
    e => jfirst(pz)
    do while (associated(e))
       if (jtype(jchild(e, 'name')) .eq. json_string) then
          call add_str(zone_names, jstr(jchild(e, 'name')))
       end if
       e => jnext(e)
    end do

    nzone_refs = 0
    npz_refs = 0
    call walk(root, '')

    summary = int_str(nzone_refs) // ' zone, ' // int_str(npz_refs) // &
         ' point zone'
    if (msgs%n .gt. n0) summary = summary // ', ' // problems(msgs%n - n0)

  contains

    !> Check the zones and point zones referred to below a value
    !! @param p The value
    !! @param path Where the value is
    recursive subroutine walk(p, path)
      type(json_value), pointer, intent(in) :: p
      character(len=*), intent(in) :: path
      type(json_value), pointer :: c
      character(len=:), allocatable :: key, cpath
      integer :: i

      ! Checked with the boundary conditions, or definitions of point zones
      if (str_eq(path, 'case.point_zones')) return
      if (is_bc_list(path)) return

      if (jtype(p) .eq. json_object) then
         if (jstr(jchild(p, 'type')) .eq. 'point_zone') then
            c => jchild(p, 'zone_name')
            if (.not. associated(c)) c => jchild(p, 'name')
            if (jtype(c) .eq. json_string) then
               call check_point_zone(c, path_key(path, jname(c)))
            end if
         end if
      end if

      c => jfirst(p)
      i = 0
      do while (associated(c))
         i = i + 1
         if (jtype(p) .eq. json_object) then
            key = jname(c)
            cpath = path_key(path, key)
            if (key .eq. 'zone_indices' .or. key .eq. 'zone_id') then
               call check_zone_ref(c, cpath)
            else if (key .eq. 'point_zone' .and. &
                 jtype(c) .eq. json_string) then
               call check_point_zone(c, cpath)
            end if
         else
            cpath = path_index(path, i)
         end if
         call walk(c, cpath)
         c => jnext(c)
      end do
    end subroutine walk

    !> Check the labeled zones a value refers to
    !! @param c The value, an index, an array of them or a constant
    !! @param path Where the value is
    subroutine check_zone_ref(c, path)
      type(json_value), pointer, intent(in) :: c
      character(len=*), intent(in) :: path
      integer, allocatable :: zones(:)
      integer :: j, z
      logical :: ok

      call int_list(c, zones, ok)
      if (.not. ok) return
      nzone_refs = nzone_refs + 1
      do j = 1, size(zones)
         z = zones(j)
         if (z .lt. 1 .or. z .gt. MAX_ZONES) then
            call report(SEV_ERROR, path, 'the zone index ' // int_str(z) // &
                 ' is out of range, zone indices are 1 to ' // &
                 int_str(MAX_ZONES), c)
         else if (mesh%read) then
            if (mesh%zone_size(z) .eq. 0) then
               call report(SEV_WARNING, path, 'there is no zone ' // &
                    int_str(z) // ' in the mesh', c)
            end if
         end if
      end do
    end subroutine check_zone_ref

    !> Check that the point zone a value names is defined
    !! @param c The value
    !! @param path Where the value is
    subroutine check_point_zone(c, path)
      type(json_value), pointer, intent(in) :: c
      character(len=*), intent(in) :: path
      character(len=:), allocatable :: name, hint
      integer :: k, d, dbest

      npz_refs = npz_refs + 1
      name = jstr(c)
      do k = 1, size(zone_names)
         if (str_eq(zone_names(k)%s, name)) return
      end do
      hint = ''
      dbest = huge(1)
      do k = 1, size(zone_names)
         d = edit_distance(name, zone_names(k)%s)
         if (d .lt. dbest .and. d .le. typo_limit(len(name))) then
            dbest = d
            hint = ' (did you mean ' // quoted(zone_names(k)%s) // '?)'
         end if
      end do
      call report(SEV_ERROR, path, 'the point zone ' // quoted(name) // &
           ' is not defined in case.point_zones' // hint, c)
    end subroutine check_point_zone

  end subroutine check_references

  !> Whether a path is the boundary_conditions of the fluid or a scalar
  !! @param path The path
  function is_bc_list(path) result(is_bc)
    character(len=*), intent(in) :: path
    logical :: is_bc
    character(len=*), parameter :: SUFFIX = '.boundary_conditions'
    character(len=:), allocatable :: parent

    is_bc = .false.
    if (len(path) .le. len(SUFFIX)) return
    if (path(len(path) - len(SUFFIX) + 1:) .ne. SUFFIX) return
    parent = path(1:len(path) - len(SUFFIX))
    if (str_eq(parent, 'case.fluid') .or. str_eq(parent, 'case.scalar')) then
       is_bc = .true.
    else if (len(parent) .gt. 13) then
       is_bc = parent(1:13) .eq. 'case.scalars[' .and. &
            index(parent(14:), '.') .eq. 0
    end if

  end function is_bc_list

  !> Warn about keys given twice and about misspelled keys
  !! @details Neko reads the first of two equal keys. A key counts as
  !! misspelled if no schema declares it for the object but it is close to a
  !! key that one does. Objects that no schema declares keys for, such as
  !! the ones a user file reads, are left alone, and so are keys shorter than
  !! four characters, which are often parameters for a user file, e.g. Ra.
  !! @param root The case file
  !! @param summary Short summary for the progress line
  subroutine check_keys(root, summary)
    type(json_value), pointer, intent(in) :: root
    character(len=:), allocatable, intent(out) :: summary
    integer :: n0

    n0 = msgs%n
    call walk(root, '')
    if (msgs%n .eq. n0) then
       summary = 'ok'
    else
       summary = problems(msgs%n - n0)
    end if

  contains

    !> Check the keys of the objects below a value
    !! @param p The value
    !! @param path Where the value is
    recursive subroutine walk(p, path)
      type(json_value), pointer, intent(in) :: p
      character(len=*), intent(in) :: path
      type(json_value), pointer :: c, d
      character(len=:), allocatable :: key
      integer :: i

      if (jtype(p) .eq. json_object) call check_object(p, path)

      c => jfirst(p)
      i = 0
      do while (associated(c))
         i = i + 1
         if (jtype(p) .eq. json_object) then
            key = jname(c)
            ! A later key equal to this one
            d => jnext(c)
            do while (associated(d))
               if (str_eq(jname(d), key)) then
                  call report(SEV_WARNING, path_key(path, key), 'the key ' &
                       // quoted(key) // ' is given more than once, Neko ' &
                       // 'uses the first one', d)
                  exit
               end if
               d => jnext(d)
            end do
            call walk(c, path_key(path, key))
         else
            call walk(c, path_index(path, i))
         end if
         c => jnext(c)
      end do
    end subroutine walk

    !> Warn about the misspelled keys of an object
    !! @param p The object
    !! @param path Where the object is
    subroutine check_object(p, path)
      type(json_value), pointer, intent(in) :: p
      character(len=*), intent(in) :: path
      type(str_t), allocatable :: declared(:)
      type(json_value), pointer :: c
      character(len=:), allocatable :: key, best
      integer :: k, d, dbest
      logical :: known

      do k = 1, anns%n
         if (anns%a(k)%kind .ne. AK_DECLARED) cycle
         if (.not. associated(anns%a(k)%node, p)) cycle
         call add_str(declared, anns%a(k)%name)
      end do
      if (.not. allocated(declared)) return

      c => jfirst(p)
      do while (associated(c))
         key = jname(c)
         if (len(key) .lt. 4) then
            c => jnext(c)
            cycle
         end if
         known = .false.
         dbest = huge(1)
         do k = 1, size(declared)
            if (str_eq(declared(k)%s, key)) then
               known = .true.
               exit
            end if
            d = edit_distance(key, declared(k)%s)
            if (d .lt. dbest) then
               dbest = d
               best = declared(k)%s
            end if
         end do
         if (.not. known .and. dbest .le. typo_limit(len(key)) .and. &
              dbest .lt. len(key)) then
            if (.not. reported(path_key(path, key), key)) then
               call report(SEV_WARNING, path_key(path, key), 'the key ' // &
                    quoted(key) // ' is not a known setting here, did ' // &
                    'you mean ' // quoted(best) // '?', c)
            end if
         end if
         c => jnext(c)
      end do
    end subroutine check_object

    !> Whether there is a message about a key already
    !! @param path Where the key is
    !! @param key The key
    function reported(path, key) result(found)
      character(len=*), intent(in) :: path, key
      logical :: found
      integer :: k

      found = .false.
      do k = 1, msgs%n
         if (str_eq(msgs%m(k)%path, path)) found = .true.
         if (index(msgs%m(k)%text, 'is ' // quoted(key) // &
              ' a misspelling') .gt. 0) found = .true.
      end do
    end function reported

  end subroutine check_keys

  !> A number of problems as text
  !! @param n The number
  function problems(n) result(text)
    integer, intent(in) :: n
    character(len=:), allocatable :: text

    if (n .eq. 1) then
       text = '1 problem'
    else
       text = int_str(n) // ' problems'
    end if

  end function problems

  !> Find the line of every value of the case file
  !! @details The text is scanned the way json-fortran parses it, with
  !! comments from /, ! or # to the end of the line, and the start of every
  !! value is recorded in document order. That is the order of a pre-order
  !! traversal of the parsed values, which pairs the lines with the values.
  !! A value in an object gets the line of its key. If the counts do not
  !! agree, no line numbers are given.
  !! @param map The line numbers
  !! @param text The case file
  !! @param root The parsed case file
  subroutine line_map_build(map, text, root)
    type(line_map_t), intent(inout) :: map
    character(len=*), intent(in) :: text
    type(json_value), pointer, intent(in) :: root
    integer, allocatable :: lines(:), stack(:)
    integer(kind=c_intptr_t), allocatable :: addr(:)
    integer :: i, n, line, key_line, depth, nv, nt
    logical :: expect_key, in_object
    character :: c

    map%n = 0
    if (allocated(map%addr)) deallocate(map%addr)
    if (allocated(map%line)) deallocate(map%line)

    n = len(text)
    allocate(lines(max(16, n / 2)), stack(max(16, n / 2)))
    nv = 0
    depth = 0
    line = 1
    key_line = 0
    expect_key = .false.
    i = 1
    do while (i .le. n)
       c = text(i:i)
       in_object = .false.
       if (depth .gt. 0) in_object = stack(depth) .eq. 1
       select case (c)
       case (achar(10))
          line = line + 1
          i = i + 1
       case ('/', '!', '#')
          do while (i .le. n)
             if (text(i:i) .eq. achar(10)) exit
             i = i + 1
          end do
       case ('"')
          if (in_object .and. expect_key) then
             key_line = line
             expect_key = .false.
          else
             call value_start()
          end if
          i = i + 1
          do while (i .le. n)
             if (text(i:i) .eq. achar(92)) then
                i = i + 1
             else if (text(i:i) .eq. '"') then
                exit
             else if (text(i:i) .eq. achar(10)) then
                line = line + 1
             end if
             i = i + 1
          end do
          i = i + 1
       case (',')
          expect_key = in_object
          i = i + 1
       case ('{', '[')
          call value_start()
          depth = depth + 1
          if (c .eq. '{') then
             stack(depth) = 1
          else
             stack(depth) = 2
          end if
          expect_key = c .eq. '{'
          i = i + 1
       case ('}', ']')
          depth = max(0, depth - 1)
          expect_key = .false.
          i = i + 1
       case (' ', ':', achar(9), achar(13))
          i = i + 1
       case default
          ! A number, true, false or null
          call value_start()
          do while (i .le. n)
             if (index(' ,:]}/!#"' // achar(9) // achar(10) // achar(13), &
                  text(i:i)) .gt. 0) exit
             i = i + 1
          end do
       end select
    end do

    ! The values in pre-order
    allocate(addr(nv + 1))
    nt = 0
    call traverse(root)
    if (nt .ne. nv) return

    map%n = nv
    map%addr = addr(1:nv)
    map%line = lines(1:nv)
    call sort_map(map)

  contains

    !> Record the line of a value that starts here
    subroutine value_start()
      in_object = .false.
      if (depth .gt. 0) in_object = stack(depth) .eq. 1
      nv = nv + 1
      if (nv .gt. size(lines)) call grow(lines)
      if (in_object .and. key_line .gt. 0) then
         lines(nv) = key_line
      else
         lines(nv) = line
      end if
      key_line = 0
    end subroutine value_start

    !> Double the size of an array
    !! @param a The array
    subroutine grow(a)
      integer, allocatable, intent(inout) :: a(:)
      integer, allocatable :: tmp(:)

      allocate(tmp(2 * size(a)))
      tmp(1:size(a)) = a
      call move_alloc(tmp, a)
    end subroutine grow

    !> Record the addresses of a value and the values below it, in pre-order
    !! @param p The value
    recursive subroutine traverse(p)
      type(json_value), pointer, intent(in) :: p
      type(json_value), pointer :: c

      nt = nt + 1
      if (nt .le. size(addr)) addr(nt) = transfer(c_loc(p), 0_c_intptr_t)
      c => jfirst(p)
      do while (associated(c))
         call traverse(c)
         c => jnext(c)
      end do
    end subroutine traverse

  end subroutine line_map_build

  !> Sort a line map by address, with a heap sort
  !! @param map The line map
  subroutine sort_map(map)
    type(line_map_t), intent(inout) :: map
    integer :: n, i, last

    n = map%n
    do i = n / 2, 1, -1
       call sift(i, n)
    end do
    do last = n, 2, -1
       call swap(1, last)
       call sift(1, last - 1)
    end do

  contains

    !> Move an entry down the heap until it is in place
    !! @param first The entry
    !! @param last Last entry of the heap
    subroutine sift(first, last)
      integer, intent(in) :: first, last
      integer :: root, child

      root = first
      do
         child = 2 * root
         if (child .gt. last) exit
         if (child .lt. last) then
            if (map%addr(child + 1) .gt. map%addr(child)) child = child + 1
         end if
         if (map%addr(root) .ge. map%addr(child)) exit
         call swap(root, child)
         root = child
      end do
    end subroutine sift

    !> Swap two entries
    !! @param a First entry
    !! @param b Second entry
    subroutine swap(a, b)
      integer, intent(in) :: a, b
      integer(kind=c_intptr_t) :: ta
      integer :: tl

      ta = map%addr(a)
      map%addr(a) = map%addr(b)
      map%addr(b) = ta
      tl = map%line(a)
      map%line(a) = map%line(b)
      map%line(b) = tl
    end subroutine swap

  end subroutine sort_map

  !> Line of a value of the case file, 0 if it is not known
  !! @param map The line map
  !! @param p The value
  function line_of(map, p) result(line)
    type(line_map_t), intent(in) :: map
    type(json_value), pointer, intent(in) :: p
    integer :: line
    integer(kind=c_intptr_t) :: a
    integer :: lo, hi, mid

    line = 0
    if (.not. associated(p) .or. map%n .eq. 0) return
    a = transfer(c_loc(p), 0_c_intptr_t)
    lo = 1
    hi = map%n
    do while (lo .le. hi)
       mid = (lo + hi) / 2
       if (map%addr(mid) .eq. a) then
          line = map%line(mid)
          return
       else if (map%addr(mid) .lt. a) then
          lo = mid + 1
       else
          hi = mid - 1
       end if
    end do

  end function line_of

end module checkneko_checks

!> Check Neko case files without running Neko
!! @details Validates the case file against the JSON schemas of the case file,
!! then checks what the schemas cannot: the constants values refer to, the
!! files the case reads, and the labeled zones of the mesh against the
!! boundary conditions. Errors are what would make Neko stop, warnings what
!! is likely a mistake. The exit status is 1 if there are errors.
program checkneko
  use, intrinsic :: iso_fortran_env, only : error_unit
  use json_module, only : json_file, json_value
  use checkneko_util
  use checkneko_schema, only : schema_init, schema_validate, schema_free, &
       schema_count, anns
  use checkneko_checks
  implicit none

  !> Maximum length of a file name
  integer, parameter :: FNAME_LEN = 4096

  type(str_t), allocatable :: case_files(:)
  character(len=:), allocatable :: schema_dir, run_dir
  integer :: i, nerrors, nfailed

  call parse_arguments()
  call print_banner()

  call jcore%initialize()
  call schema_init(schema_dir)
  write(*, '(A)') ' Schemas:    ' // schema_dir

  nfailed = 0
  do i = 1, size(case_files)
     call check_case(case_files(i)%s, nerrors)
     if (nerrors .gt. 0) nfailed = nfailed + 1
  end do

  call schema_free()
  call jcore%destroy()

  if (nfailed .gt. 0) call exit_with(1)

contains

  !> Print the banner with the version and the compiler
  subroutine print_banner()
    use, intrinsic :: iso_fortran_env, only : compiler_version

    write(*, '(A)') ''
    write(*, '(A)') ' N E K O  case file checker, version ' // &
         PACKAGE_VERSION
    write(*, '(A)') ' Built with ' // compiler_version()
    write(*, '(A)') ''

  end subroutine print_banner

  !> The usage line
  function usage() result(text)
    character(len=:), allocatable :: text

    text = 'usage: checkneko [-h] [-d dir] [-s dir] file ...'

  end function usage

  !> Print the usage and what the options do
  subroutine print_help()
    write(*, '(A)') usage()
    write(*, '(A)') ''
    write(*, '(A)') 'Check Neko case files without running Neko. The ' // &
         'exit status is 1 if there'
    write(*, '(A)') 'are errors.'
    write(*, '(A)') ''
    write(*, '(A)') '  -d dir  Directory Neko will run in, which file ' // &
         'names in the case are'
    write(*, '(A)') '          relative to (default: the current directory)'
    write(*, '(A)') '  -s dir  Read the schemas from dir, e.g. ' // &
         'doc/schemas, instead of'
    write(*, '(A)') '          the ones installed in share/neko/schemas'
    write(*, '(A)') '  -h      Print this help'
  end subroutine print_help

  !> Stop after a wrong command line, printing the usage
  !! @param msg What was wrong, empty to print just the usage
  subroutine usage_error(msg)
    character(len=*), intent(in) :: msg

    if (len(msg) .gt. 0) write(error_unit, '(A)') 'checkneko: ' // msg
    write(error_unit, '(A)') usage()
    call exit_with(1)

  end subroutine usage_error

  !> Parse the command line, with options in the style of getopt(3)
  !! @details Options come before the files. Flags can be grouped, the value
  !! of an option can follow it directly or as the next argument, and --
  !! ends the options.
  subroutine parse_arguments()
    character(len=FNAME_LEN) :: arg, val
    character :: opt
    integer :: argc, i, j, n

    schema_dir = ''
    run_dir = '.'
    allocate(case_files(0))

    argc = command_argument_count()
    i = 1
    do while (i .le. argc)
       call get_command_argument(i, arg)
       n = len_trim(arg)
       if (trim(arg) .eq. '--') then
          i = i + 1
          exit
       end if
       if (arg(1:1) .ne. '-' .or. n .eq. 1) exit

       j = 2
       do while (j .le. n)
          opt = arg(j:j)
          select case (opt)
          case ('h')
             call print_help()
             call exit_with(0)
          case ('d', 's')
             if (j .lt. n) then
                val = arg(j + 1:)
             else if (i .lt. argc) then
                i = i + 1
                call get_command_argument(i, val)
             else
                call usage_error('option requires an argument -- ' // opt)
             end if
             if (opt .eq. 'd') then
                run_dir = trim(val)
             else
                schema_dir = trim(val)
             end if
             exit
          case default
             call usage_error('illegal option -- ' // opt)
          end select
          j = j + 1
       end do
       i = i + 1
    end do

    do while (i .le. argc)
       call get_command_argument(i, arg)
       case_files = [case_files, str_t(trim(arg))]
       i = i + 1
    end do
    if (size(case_files) .eq. 0) call usage_error('')

    if (len(schema_dir) .gt. 0) then
       if (.not. file_exists(join_path(schema_dir, &
            'case-file.schema.json'))) then
          call fatal(schema_dir // ': no case-file.schema.json')
       end if
    else
       schema_dir = find_schema_dir()
       if (len(schema_dir) .eq. 0) then
          call fatal('cannot find the schemas of the case file in ' // &
               'share/neko/schemas, give their directory with -s')
       end if
    end if

  end subroutine parse_arguments

  !> Directory of the installed schemas, empty if they cannot be found
  !! @details The schemas are installed in share/neko/schemas, next to the
  !! bin directory of checkneko, which is found from the command it was
  !! started with. When checkneko is run from the build directory of a Neko
  !! source tree, the doc/schemas of the tree are used.
  function find_schema_dir() result(dir)
    character(len=:), allocatable :: dir
    character(len=*), parameter :: CANDIDATES(2) = [character(len=21) :: &
         '../share/neko/schemas', '../../doc/schemas']
    character(len=:), allocatable :: exe, candidate
    integer :: i

    dir = ''
    exe = executable_path()
    i = scan(exe, '/', back = .true.)
    if (i .eq. 0) return
    do i = 1, size(CANDIDATES)
       candidate = exe(1:scan(exe, '/', back = .true.)) // trim(CANDIDATES(i))
       if (file_exists(join_path(candidate, 'case-file.schema.json'))) then
          dir = real_path(candidate)
          return
       end if
    end do

  end function find_schema_dir

  !> Path of the checkneko executable, empty if it cannot be found
  !! @details A command without a directory was found through PATH, which is
  !! searched as the shell does. Symbolic links are resolved, so that a link
  !! to checkneko in another directory, as in a Spack view, gives the
  !! directory it is installed in.
  function executable_path() result(path)
    character(len=:), allocatable :: path
    character(len=FNAME_LEN) :: arg
    character(len=:), allocatable :: env, candidate
    integer :: n, status, i, j

    path = ''
    call get_command_argument(0, arg)
    if (len_trim(arg) .eq. 0) return

    if (index(arg, '/') .gt. 0) then
       path = trim(arg)
    else
       call get_environment_variable('PATH', length = n, status = status)
       if (status .ne. 0 .or. n .eq. 0) return
       allocate(character(len=n) :: env)
       call get_environment_variable('PATH', env)
       ! Entries are separated by colons, an empty one is the current
       ! directory
       i = 1
       do while (i .le. n + 1)
          j = index(env(i:), ':')
          if (j .eq. 0) then
             j = n + 1
          else
             j = i + j - 1
          end if
          candidate = join_path(env(i:j - 1), trim(arg))
          if (file_exists(candidate)) then
             path = candidate
             exit
          end if
          i = j + 1
       end do
       if (len(path) .eq. 0) return
    end if

    path = real_path(path)

  end function executable_path

  !> Check a case file and print what was found
  !! @param fname Name of the case file
  !! @param nerrors Number of errors found
  subroutine check_case(fname, nerrors)
    character(len=*), intent(in) :: fname
    integer, intent(out) :: nerrors
    type(json_file) :: case_json
    type(json_value), pointer :: root, null_node
    type(line_map_t) :: map
    type(mesh_info_t) :: mesh
    character(len=:), allocatable :: text, summary, msg, suffix
    logical :: ok
    integer :: nwarnings, n0

    nerrors = 0
    msgs%n = 0
    anns%n = 0
    null_node => null()

    write(*, '(A)') ' Case file:  ' // fname
    write(*, '(A)') ' Run dir:    ' // dir_text(run_dir)
    write(*, '(A)') ''

    call read_text(fname, text, ok)
    if (.not. ok) then
       write(*, '(A)') fname // ': error: cannot read the file'
       nerrors = 1
       call print_summary(nerrors, 0)
       return
    end if

    ! Parse as Neko does, with json-fortran and its default settings
    call case_json%initialize()
    call case_json%load(filename = fname)
    if (case_json%failed()) then
       call case_json%check_for_errors(ok, msg)
       call progress('Parsing', 'failed')
       write(*, '(A)') ''
       call print_parse_error(fname, text, msg)
       call case_json%destroy()
       nerrors = 1
       call print_summary(nerrors, 0)
       return
    end if
    call case_json%get(root)
    call line_map_build(map, text, root)
    if (map%n .gt. 0) then
       call progress('Parsing', 'ok')
    else
       call progress('Parsing', 'ok, without line numbers')
    end if

    ! Neko refuses other names
    suffix = file_suffix(fname)
    if (suffix .ne. 'case' .and. suffix .ne. 'json') then
       call report(SEV_ERROR, '', 'Neko only reads case files whose ' // &
            'name ends in .case or .json', null_node)
    end if

    n0 = msgs%n
    call schema_validate(root)
    if (msgs%n .eq. n0) then
       call progress('Schema', 'ok (' // int_str(schema_count()) // &
            ' schemas)')
    else
       call progress('Schema', count_text(msgs%n - n0, 'problem'))
    end if

    call check_constants(root, summary)
    call progress('Constants', summary)
    call check_input_files(root, run_dir, summary)
    call progress('Input files', summary)
    call check_mesh(root, run_dir, mesh, summary)
    call progress('Mesh', summary)
    call check_boundary_conditions(root, mesh, summary)
    call progress('Boundary conditions', summary)
    call check_references(root, mesh, summary)
    call progress('References', summary)
    call check_keys(root, summary)
    call progress('Keys', summary)

    call print_zone_table(mesh)
    call print_messages(fname, map, nerrors, nwarnings)
    call print_summary(nerrors, nwarnings)

    call case_json%destroy()

  end subroutine check_case

  !> A directory for the header, with the current one by its name
  !! @details The name comes from getcwd(3), Fortran has no portable way to
  !! get it.
  !! @param dir The directory
  function dir_text(dir) result(text)
    use, intrinsic :: iso_c_binding, only : c_char, c_ptr, c_size_t, &
         c_associated
    character(len=*), intent(in) :: dir
    character(len=:), allocatable :: text
    character(kind=c_char) :: buf(FNAME_LEN)
    interface
       function c_getcwd(buf, size) bind(c, name = 'getcwd') result(p)
         import :: c_char, c_size_t, c_ptr
         character(kind=c_char), intent(out) :: buf(*)
         integer(kind=c_size_t), value :: size
         type(c_ptr) :: p
       end function c_getcwd
    end interface

    text = dir
    if (dir .ne. '.') return
    if (.not. c_associated(c_getcwd(buf, int(FNAME_LEN, c_size_t)))) return
    text = c_string(buf)

  end function dir_text

  !> A path with symbolic links, . and .. resolved by realpath(3), or the
  !! path as it is if that fails
  !! @param path The path
  function real_path(path) result(rpath)
    use, intrinsic :: iso_c_binding, only : c_char, c_null_char, c_ptr, &
         c_associated
    character(len=*), intent(in) :: path
    character(len=:), allocatable :: rpath
    character(kind=c_char) :: buf(FNAME_LEN)
    interface
       function c_realpath(path, resolved) bind(c, name = 'realpath') &
            result(p)
         import :: c_char, c_ptr
         character(kind=c_char), intent(in) :: path(*)
         character(kind=c_char), intent(out) :: resolved(*)
         type(c_ptr) :: p
       end function c_realpath
    end interface

    rpath = path
    if (.not. c_associated(c_realpath(path // c_null_char, buf))) return
    rpath = c_string(buf)

  end function real_path

  !> A null-terminated C string as a Fortran string
  !! @param buf The C string
  function c_string(buf) result(s)
    use, intrinsic :: iso_c_binding, only : c_char, c_null_char
    character(kind=c_char), intent(in) :: buf(:)
    character(len=:), allocatable :: s
    integer :: i

    s = ''
    do i = 1, size(buf)
       if (buf(i) .eq. c_null_char) exit
       s = s // buf(i)
    end do

  end function c_string

  !> Read a file into a string
  !! @param fname Name of the file
  !! @param text Contents of the file
  !! @param ok Whether the file could be read
  subroutine read_text(fname, text, ok)
    character(len=*), intent(in) :: fname
    character(len=:), allocatable, intent(out) :: text
    logical, intent(out) :: ok
    integer :: unit, ios
    integer(kind=i8) :: fsize

    ok = .false.
    text = ''
    open(newunit = unit, file = fname, access = 'stream', &
         form = 'unformatted', status = 'old', action = 'read', iostat = ios)
    if (ios .ne. 0) return
    inquire(unit = unit, size = fsize)
    if (fsize .gt. 0) then
       deallocate(text)
       allocate(character(len=fsize) :: text)
       read(unit, iostat = ios) text
    end if
    close(unit)
    ok = ios .eq. 0

  end subroutine read_text

  !> Print a parse error of json-fortran as file:line: error: text, followed
  !! by the line with a mark under the character where parsing stopped
  !! @param fname Name of the case file
  !! @param text The case file
  !! @param msg The error of json-fortran, with its position
  subroutine print_parse_error(fname, text, msg)
    character(len=*), intent(in) :: fname, text, msg
    character(len=:), allocatable :: first
    integer :: i, line, col, lstart, lend

    i = index(msg, achar(10))
    if (i .eq. 0) i = len(msg) + 1
    first = msg(1:i - 1)
    if (index(first, 'Error in ') .eq. 1) then
       if (index(first, ': ') .gt. 0) first = first(index(first, ': ') + 2:)
    end if

    line = position_after(msg, 'line:')
    col = position_after(msg, 'character:')
    if (line .le. 0) then
       write(*, '(A)') fname // ': error: invalid JSON: ' // first
       return
    end if
    write(*, '(A)') fname // ':' // int_str(line) // &
         ': error: invalid JSON: ' // first

    ! The line itself
    lstart = 1
    do i = 1, line - 1
       lend = index(text(lstart:), achar(10))
       if (lend .eq. 0) return
       lstart = lstart + lend
    end do
    lend = index(text(lstart:), achar(10))
    if (lend .eq. 0) then
       lend = len(text)
    else
       lend = lstart + lend - 2
    end if
    write(*, '(A)') '  ' // text(lstart:lend)
    if (col .ge. 1) write(*, '(A)') '  ' // repeat(' ', col - 1) // '^'

  end subroutine print_parse_error

  !> The number after a label in a message of json-fortran, 0 if not found
  !! @param msg The message
  !! @param label The label, e.g. "line:"
  function position_after(msg, label) result(val)
    character(len=*), intent(in) :: msg, label
    integer :: val
    integer :: i, j, ios

    val = 0
    i = index(msg, label)
    if (i .eq. 0) return
    i = i + len(label)
    j = i
    do while (j .le. len(msg))
       if (index(' 0123456789', msg(j:j)) .eq. 0) exit
       j = j + 1
    end do
    if (j .eq. i) return
    read(msg(i:j - 1), *, iostat = ios) val
    if (ios .ne. 0) val = 0

  end function position_after

  !> Print the messages, sorted by line, as file:line: severity: path: text
  !! @param fname Name of the case file
  !! @param map Line numbers of the values
  !! @param nerrors Number of errors printed
  !! @param nwarnings Number of warnings printed
  subroutine print_messages(fname, map, nerrors, nwarnings)
    character(len=*), intent(in) :: fname
    type(line_map_t), intent(in) :: map
    integer, intent(out) :: nerrors, nwarnings
    integer, allocatable :: lines(:), order(:)
    character(len=:), allocatable :: loc, sev, where
    integer :: k, j, t
    logical :: dup

    nerrors = 0
    nwarnings = 0
    if (msgs%n .eq. 0) return

    allocate(lines(msgs%n), order(msgs%n))
    do k = 1, msgs%n
       lines(k) = line_of(map, msgs%m(k)%node)
       if (lines(k) .eq. 0) lines(k) = huge(1)
       order(k) = k
    end do
    ! Stable insertion sort by line
    do k = 2, msgs%n
       t = order(k)
       j = k - 1
       do while (j .ge. 1)
          if (lines(order(j)) .le. lines(t)) exit
          order(j + 1) = order(j)
          j = j - 1
       end do
       order(j + 1) = t
    end do

    write(*, '(A)') ''
    do k = 1, msgs%n
       associate (m => msgs%m(order(k)))
         dup = .false.
         do j = 1, k - 1
            if (msgs%m(order(j))%sev .eq. m%sev .and. &
                 str_eq(msgs%m(order(j))%path, m%path) .and. &
                 str_eq(msgs%m(order(j))%text, m%text)) dup = .true.
         end do
         if (dup) cycle

         if (lines(order(k)) .lt. huge(1)) then
            loc = fname // ':' // int_str(lines(order(k)))
         else
            loc = fname
         end if
         if (m%sev .eq. SEV_ERROR) then
            sev = 'error'
            nerrors = nerrors + 1
         else
            sev = 'warning'
            nwarnings = nwarnings + 1
         end if
         if (len(m%path) .eq. 0) then
            where = ''
         else
            where = m%path // ': '
         end if
         write(*, '(A)') loc // ': ' // sev // ': ' // where // m%text
       end associate
    end do

  end subroutine print_messages

  !> Print the number of errors and warnings
  !! @param nerrors Number of errors
  !! @param nwarnings Number of warnings
  subroutine print_summary(nerrors, nwarnings)
    integer, intent(in) :: nerrors, nwarnings

    write(*, '(A)') ''
    if (nerrors .eq. 0 .and. nwarnings .eq. 0) then
       write(*, '(A)') ' No problems found'
    else
       write(*, '(A)') ' ' // count_text(nerrors, 'error') // ', ' // &
            count_text(nwarnings, 'warning')
    end if
    write(*, '(A)') ''

  end subroutine print_summary

  !> A count and a noun, e.g. "1 error" or "2 errors"
  !! @param n The count
  !! @param noun The noun
  function count_text(n, noun) result(text)
    integer, intent(in) :: n
    character(len=*), intent(in) :: noun
    character(len=:), allocatable :: text

    if (n .eq. 1) then
       text = '1 ' // noun
    else
       text = int_str(n) // ' ' // noun // 's'
    end if

  end function count_text

end program checkneko
