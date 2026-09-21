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
!> Evaluation of mathematical expressions given as strings in the case file.
!! @details An `expression_t` is compiled once, from a string such as
!! `"6*U_b*y*(H - y)/H^2"`, into a stream of stack machine opcodes. It can then
!! be evaluated over an arbitrary set of points, with `x`, `y` and `z` bound to
!! the point coordinates and `t`, `dt` bound to the time state.
!!
!! Any identifier which is not a coordinate, a time variable or the constant
!! `pi` is looked up in the `neko_const_registry` when the expression is
!! compiled, and folded into a literal. This is the same mechanism as the one
!! used by `json_get_or_lookup`, so constants declared under `case.constants`
!! in the case file can be used directly in the expressions.
module expression
  use num_types, only : rp, dp
  use utils, only : neko_error
  use registry, only : neko_const_registry
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
  implicit none
  private

  public :: expression_eval_static, expression_check_finite

  !> Maximum length of an expression string read from the case file.
  integer, public, parameter :: NEKO_EXPR_LEN = 256

  !> Maximum length of an identifier appearing in an expression.
  integer, parameter :: EXPR_MAX_IDENT = 64

  !> Number of points evaluated per pass over the opcode stream. Chosen so
  !! that the evaluation stack stays in cache irrespective of how many points
  !! the caller hands over in one go.
  integer, parameter :: EXPR_CHUNK = 1024

  !
  ! Opcodes. The evaluator itself does not care about the numbering, but the
  ! parser derives the arity of an opcode from the ranges below, so nullary
  ! opcodes must stay below OP_UNARY and binary ones at or above OP_BINARY.
  !
  integer, parameter :: OP_LIT = 1, OP_X = 2, OP_Y = 3, OP_Z = 4, &
       OP_T = 5, OP_DT = 6

  integer, parameter :: OP_UNARY = 20
  integer, parameter :: OP_NEG = 20, OP_SIN = 21, OP_COS = 22, OP_TAN = 23, &
       OP_ASIN = 24, OP_ACOS = 25, OP_ATAN = 26, OP_SINH = 27, &
       OP_COSH = 28, OP_TANH = 29, OP_EXP = 30, OP_LOG = 31, &
       OP_LOG10 = 32, OP_SQRT = 33, OP_ABS = 34, OP_ERF = 35, &
       OP_ERFC = 36, OP_STEP = 37, OP_IPOW = 38

  integer, parameter :: OP_BINARY = 60
  integer, parameter :: OP_ADD = 60, OP_SUB = 61, OP_MUL = 62, OP_DIV = 63, &
       OP_POW = 64, OP_ATAN2 = 65, OP_MIN = 66, OP_MAX = 67, OP_MOD = 68

  !
  ! Token kinds produced by the tokenizer.
  !
  integer, parameter :: TK_END = 0, TK_NUM = 1, TK_IDENT = 2, TK_PLUS = 3, &
       TK_MINUS = 4, TK_STAR = 5, TK_SLASH = 6, TK_POW = 7, TK_LPAR = 8, &
       TK_RPAR = 9, TK_COMMA = 10

  !> A single token of an expression.
  type :: token_t
     integer :: kind = TK_END
     real(kind=rp) :: num = 0.0_rp
     character(len=EXPR_MAX_IDENT) :: name = ''
     !> Position of the first character of the token in the source string,
     !! used for error messages.
     integer :: pos = 0
  end type token_t

  !> State of the recursive descent parser.
  type :: parser_t
     type(token_t), allocatable :: tok(:)
     integer :: n_tok = 0
     integer :: cur = 1
     character(len=:), allocatable :: src
     integer, allocatable :: op(:)
     real(kind=rp), allocatable :: lit(:)
     integer :: n_op = 0
     !> Current and maximum depth of the evaluation stack.
     integer :: depth = 0
     integer :: max_depth = 0
     logical :: time_dependent = .false.
  end type parser_t

  !> A compiled mathematical expression.
  type, public :: expression_t
     !> Opcode stream, in postfix order.
     integer, allocatable :: op(:)
     !> Literal operand, only meaningful for `OP_LIT` and `OP_IPOW`. Kept
     !! parallel to `op` to avoid a second index while evaluating.
     real(kind=rp), allocatable :: lit(:)
     !> Number of opcodes.
     integer :: n_op = 0
     !> Maximum depth of the evaluation stack.
     integer :: depth = 0
     !> True if the expression references `t` or `dt`.
     logical :: time_dependent = .false.
     !> The source string, kept for logging and error messages.
     character(len=:), allocatable :: src
     !> Evaluation stack, allocated once when the expression is compiled.
     real(kind=rp), allocatable, private :: stk(:,:)
   contains
     !> Compile an expression from a string.
     procedure, pass(this) :: init => expression_init
     !> Evaluate the expression in a set of points.
     procedure, pass(this) :: eval => expression_eval
     !> Destructor.
     procedure, pass(this) :: free => expression_free
  end type expression_t

contains

  !> Compile an expression.
  !! @param[in] str The expression, e.g. `"6*U_b*y*(H - y)/H^2"`.
  subroutine expression_init(this, str)
    class(expression_t), intent(inout) :: this
    character(len=*), intent(in) :: str
    type(parser_t) :: p

    call this%free()

    p%src = trim(str)
    if (len(p%src) .eq. 0) call neko_error("Empty expression")

    call expr_tokenize(p)

    allocate(p%op(p%n_tok + 1))
    allocate(p%lit(p%n_tok + 1))

    call parse_expr(p)
    if (p%tok(p%cur)%kind .ne. TK_END) then
       call expr_error(p, "unexpected trailing input", p%tok(p%cur)%pos)
    end if

    this%n_op = p%n_op
    this%depth = p%max_depth
    this%time_dependent = p%time_dependent
    this%src = p%src

    allocate(this%op(this%n_op), this%lit(this%n_op))
    this%op = p%op(1:this%n_op)
    this%lit = p%lit(1:this%n_op)

    allocate(this%stk(EXPR_CHUNK, max(this%depth, 1)))

    call parser_free(p)

  end subroutine expression_init

  !> Release the memory held by the parser state.
  !! @param[inout] p The parser state.
  subroutine parser_free(p)
    type(parser_t), intent(inout) :: p

    if (allocated(p%tok)) deallocate(p%tok)
    if (allocated(p%src)) deallocate(p%src)
    if (allocated(p%op)) deallocate(p%op)
    if (allocated(p%lit)) deallocate(p%lit)

    p%n_tok = 0
    p%cur = 1
    p%n_op = 0
    p%depth = 0
    p%max_depth = 0
    p%time_dependent = .false.

  end subroutine parser_free

  !> Destructor.
  subroutine expression_free(this)
    class(expression_t), intent(inout) :: this

    if (allocated(this%op)) deallocate(this%op)
    if (allocated(this%lit)) deallocate(this%lit)
    if (allocated(this%src)) deallocate(this%src)
    if (allocated(this%stk)) deallocate(this%stk)
    this%n_op = 0
    this%depth = 0
    this%time_dependent = .false.

  end subroutine expression_free

  !> Evaluate the expression in `n` points.
  !! @details Not thread safe. The evaluation stack is a component of the
  !! expression, so two concurrent calls on the same object corrupt each
  !! other's intermediate values. A caller inside an OpenMP parallel region
  !! must evaluate on one thread only.
  !! @param[inout] res The result, one value per point.
  !! @param[in] n The number of points.
  !! @param[in] x The x-coordinates of the points.
  !! @param[in] y The y-coordinates of the points.
  !! @param[in] z The z-coordinates of the points.
  !! @param[in] t The current time, required if the expression uses `t`.
  !! @param[in] dt The current timestep, required if the expression uses `dt`.
  subroutine expression_eval(this, res, n, x, y, z, t, dt)
    class(expression_t), intent(inout) :: this
    integer, intent(in) :: n
    real(kind=rp), intent(inout) :: res(n)
    real(kind=rp), intent(in) :: x(n)
    real(kind=rp), intent(in) :: y(n)
    real(kind=rp), intent(in) :: z(n)
    real(kind=dp), intent(in), optional :: t
    real(kind=dp), intent(in), optional :: dt
    real(kind=dp) :: t_, dt_
    integer :: i0, nb, ip, sp, i, k

    if (this%n_op .eq. 0) call neko_error("Evaluating an uncompiled expression")

    if (this%time_dependent .and. .not. present(t)) then
       call neko_error("The expression '" // this%src // "' depends on " // &
            "time, but no time is available where it is used")
    end if

    t_ = 0.0_dp
    dt_ = 0.0_dp
    if (present(t)) t_ = t
    if (present(dt)) dt_ = dt

    do i0 = 1, n, EXPR_CHUNK
       nb = min(EXPR_CHUNK, n - i0 + 1)
       sp = 0

       do ip = 1, this%n_op
          select case (this%op(ip))

             !
             ! Operand pushes
             !
          case (OP_LIT)
             sp = sp + 1
             do i = 1, nb
                this%stk(i, sp) = this%lit(ip)
             end do
          case (OP_X)
             sp = sp + 1
             do i = 1, nb
                this%stk(i, sp) = x(i0 + i - 1)
             end do
          case (OP_Y)
             sp = sp + 1
             do i = 1, nb
                this%stk(i, sp) = y(i0 + i - 1)
             end do
          case (OP_Z)
             sp = sp + 1
             do i = 1, nb
                this%stk(i, sp) = z(i0 + i - 1)
             end do
          case (OP_T)
             sp = sp + 1
             do i = 1, nb
                this%stk(i, sp) = t_
             end do
          case (OP_DT)
             sp = sp + 1
             do i = 1, nb
                this%stk(i, sp) = dt_
             end do

             !
             ! Unary operators
             !
          case (OP_NEG)
             do i = 1, nb
                this%stk(i, sp) = -this%stk(i, sp)
             end do
          case (OP_SIN)
             do i = 1, nb
                this%stk(i, sp) = sin(this%stk(i, sp))
             end do
          case (OP_COS)
             do i = 1, nb
                this%stk(i, sp) = cos(this%stk(i, sp))
             end do
          case (OP_TAN)
             do i = 1, nb
                this%stk(i, sp) = tan(this%stk(i, sp))
             end do
          case (OP_ASIN)
             do i = 1, nb
                this%stk(i, sp) = asin(this%stk(i, sp))
             end do
          case (OP_ACOS)
             do i = 1, nb
                this%stk(i, sp) = acos(this%stk(i, sp))
             end do
          case (OP_ATAN)
             do i = 1, nb
                this%stk(i, sp) = atan(this%stk(i, sp))
             end do
          case (OP_SINH)
             do i = 1, nb
                this%stk(i, sp) = sinh(this%stk(i, sp))
             end do
          case (OP_COSH)
             do i = 1, nb
                this%stk(i, sp) = cosh(this%stk(i, sp))
             end do
          case (OP_TANH)
             do i = 1, nb
                this%stk(i, sp) = tanh(this%stk(i, sp))
             end do
          case (OP_EXP)
             do i = 1, nb
                this%stk(i, sp) = exp(this%stk(i, sp))
             end do
          case (OP_LOG)
             do i = 1, nb
                this%stk(i, sp) = log(this%stk(i, sp))
             end do
          case (OP_LOG10)
             do i = 1, nb
                this%stk(i, sp) = log10(this%stk(i, sp))
             end do
          case (OP_SQRT)
             do i = 1, nb
                this%stk(i, sp) = sqrt(this%stk(i, sp))
             end do
          case (OP_ABS)
             do i = 1, nb
                this%stk(i, sp) = abs(this%stk(i, sp))
             end do
          case (OP_ERF)
             do i = 1, nb
                this%stk(i, sp) = erf(this%stk(i, sp))
             end do
          case (OP_ERFC)
             do i = 1, nb
                this%stk(i, sp) = erfc(this%stk(i, sp))
             end do
          case (OP_STEP)
             do i = 1, nb
                if (this%stk(i, sp) .lt. 0.0_rp) then
                   this%stk(i, sp) = 0.0_rp
                else
                   this%stk(i, sp) = 1.0_rp
                end if
             end do
          case (OP_IPOW)
             ! Integer exponent, folded by the parser so that a negative base
             ! raised to a whole power stays well defined.
             k = nint(this%lit(ip))
             do i = 1, nb
                this%stk(i, sp) = this%stk(i, sp) ** k
             end do

             !
             ! Binary operators
             !
          case (OP_ADD)
             do i = 1, nb
                this%stk(i, sp-1) = this%stk(i, sp-1) + this%stk(i, sp)
             end do
             sp = sp - 1
          case (OP_SUB)
             do i = 1, nb
                this%stk(i, sp-1) = this%stk(i, sp-1) - this%stk(i, sp)
             end do
             sp = sp - 1
          case (OP_MUL)
             do i = 1, nb
                this%stk(i, sp-1) = this%stk(i, sp-1) * this%stk(i, sp)
             end do
             sp = sp - 1
          case (OP_DIV)
             do i = 1, nb
                this%stk(i, sp-1) = this%stk(i, sp-1) / this%stk(i, sp)
             end do
             sp = sp - 1
          case (OP_POW)
             do i = 1, nb
                this%stk(i, sp-1) = this%stk(i, sp-1) ** this%stk(i, sp)
             end do
             sp = sp - 1
          case (OP_ATAN2)
             do i = 1, nb
                this%stk(i, sp-1) = atan2(this%stk(i, sp-1), this%stk(i, sp))
             end do
             sp = sp - 1
          case (OP_MIN)
             do i = 1, nb
                this%stk(i, sp-1) = min(this%stk(i, sp-1), this%stk(i, sp))
             end do
             sp = sp - 1
          case (OP_MAX)
             do i = 1, nb
                this%stk(i, sp-1) = max(this%stk(i, sp-1), this%stk(i, sp))
             end do
             sp = sp - 1
          case (OP_MOD)
             do i = 1, nb
                this%stk(i, sp-1) = mod(this%stk(i, sp-1), this%stk(i, sp))
             end do
             sp = sp - 1

          case default
             call neko_error("Corrupt expression opcode stream")
          end select
       end do

       if (sp .ne. 1) call neko_error("Corrupt expression opcode stream")

       do i = 1, nb
          res(i0 + i - 1) = this%stk(i, 1)
       end do
    end do

  end subroutine expression_eval

  !> Compile an expression and evaluate it in a set of points, in a context
  !! where there is no time state.
  !! @details Convenience routine for initial conditions and other setup-time
  !! uses of an expression. The expression is rejected if it is empty, if it
  !! depends on time, or if it does not evaluate to a finite value in every
  !! point, the latter typically meaning a division by zero or the square root
  !! of a negative number somewhere in the domain.
  !! @param[in] str The expression.
  !! @param[inout] res The result, one value per point.
  !! @param[in] n The number of points.
  !! @param[in] x The x-coordinates of the points.
  !! @param[in] y The y-coordinates of the points.
  !! @param[in] z The z-coordinates of the points.
  !! @param[in] usage What the expression configures, used in error messages.
  subroutine expression_eval_static(str, res, n, x, y, z, usage)
    character(len=*), intent(in) :: str
    integer, intent(in) :: n
    real(kind=rp), intent(inout) :: res(n)
    real(kind=rp), intent(in) :: x(n)
    real(kind=rp), intent(in) :: y(n)
    real(kind=rp), intent(in) :: z(n)
    character(len=*), intent(in) :: usage
    type(expression_t) :: e

    if (len_trim(str) .eq. 0) then
       call neko_error("The " // usage // " needs a non-empty expression")
    end if

    call e%init(str)

    if (e%time_dependent) then
       call neko_error("The " // usage // " expression '" // trim(str) // &
            "' depends on time, which is not available where it is used")
    end if

    call e%eval(res, n, x, y, z)
    call e%free()

    call expression_check_finite(str, res, n, usage)

  end subroutine expression_eval_static

  !> Abort if an expression did not evaluate to a finite value everywhere.
  !! @details Split out of `expression_eval_static` so that the boundary
  !! conditions, which compile and evaluate their expressions themselves, can
  !! apply the same check. A non-finite value typically means a division by
  !! zero or the square root of a negative number somewhere.
  !! @param[in] str The expression, used in the error message.
  !! @param[in] res The values to check.
  !! @param[in] n The number of values.
  !! @param[in] usage What the expression configures, used in error messages.
  subroutine expression_check_finite(str, res, n, usage)
    character(len=*), intent(in) :: str
    integer, intent(in) :: n
    real(kind=rp), intent(in) :: res(n)
    character(len=*), intent(in) :: usage
    integer :: i

    do i = 1, n
       if (.not. ieee_is_finite(res(i))) then
          call neko_error("The " // usage // " expression '" // trim(str) // &
               "' does not evaluate to a finite value everywhere")
       end if
    end do

  end subroutine expression_check_finite

  !
  ! Tokenizer
  !

  !> Split the source string into a stream of tokens.
  !! @param[inout] p The parser state, whose `src` holds the expression.
  subroutine expr_tokenize(p)
    type(parser_t), intent(inout) :: p
    integer :: i, j, k, n, ios
    logical :: seen_dot
    character(len=1) :: c
    character(len=EXPR_MAX_IDENT) :: buf

    n = len(p%src)
    allocate(p%tok(n + 1))
    p%n_tok = 0

    i = 1
    do while (i .le. n)
       c = p%src(i:i)

       if (c .eq. ' ' .or. c .eq. achar(9)) then
          i = i + 1
          cycle
       end if

       if (is_digit(c) .or. (c .eq. '.' .and. i .lt. n .and. &
            is_digit(p%src(min(i+1, n):min(i+1, n))))) then
          ! Mantissa
          j = i
          seen_dot = .false.
          do while (j .le. n)
             if (is_digit(p%src(j:j))) then
                j = j + 1
             else if (p%src(j:j) .eq. '.' .and. .not. seen_dot) then
                seen_dot = .true.
                j = j + 1
             else
                exit
             end if
          end do
          ! Exponent, but only if it is actually followed by digits, so that
          ! e.g. `2*d` is not mistaken for a malformed number.
          if (j .le. n) then
             if (scan(p%src(j:j), 'eEdD') .gt. 0) then
                k = j + 1
                if (k .le. n) then
                   if (scan(p%src(k:k), '+-') .gt. 0) k = k + 1
                end if
                if (k .le. n) then
                   if (is_digit(p%src(k:k))) then
                      do while (k .le. n)
                         if (.not. is_digit(p%src(k:k))) exit
                         k = k + 1
                      end do
                      j = k
                   end if
                end if
             end if
          end if

          if (j - i .gt. EXPR_MAX_IDENT) then
             call expr_error(p, "numeric literal is too long", i)
          end if
          buf = p%src(i:j-1)
          ! Fortran list-directed input does not accept a `d` exponent in
          ! every implementation, so normalise it.
          do k = 1, len_trim(buf)
             if (buf(k:k) .eq. 'd' .or. buf(k:k) .eq. 'D') buf(k:k) = 'e'
          end do
          call push_token(p, TK_NUM, i)
          read (buf, *, iostat = ios) p%tok(p%n_tok)%num
          if (ios .ne. 0) call expr_error(p, "malformed number", i)
          i = j
          cycle
       end if

       if (is_alpha(c)) then
          j = i
          do while (j .le. n)
             if (is_alpha(p%src(j:j)) .or. is_digit(p%src(j:j)) .or. &
                  p%src(j:j) .eq. '_') then
                j = j + 1
             else
                exit
             end if
          end do
          if (j - i .gt. EXPR_MAX_IDENT) then
             call expr_error(p, "identifier is too long", i)
          end if
          call push_token(p, TK_IDENT, i)
          p%tok(p%n_tok)%name = p%src(i:j-1)
          i = j
          cycle
       end if

       select case (c)
       case ('+')
          call push_token(p, TK_PLUS, i)
       case ('-')
          call push_token(p, TK_MINUS, i)
       case ('*')
          if (i .lt. n) then
             if (p%src(i+1:i+1) .eq. '*') then
                call push_token(p, TK_POW, i)
                i = i + 2
                cycle
             end if
          end if
          call push_token(p, TK_STAR, i)
       case ('/')
          call push_token(p, TK_SLASH, i)
       case ('^')
          call push_token(p, TK_POW, i)
       case ('(')
          call push_token(p, TK_LPAR, i)
       case (')')
          call push_token(p, TK_RPAR, i)
       case (',')
          call push_token(p, TK_COMMA, i)
       case default
          call expr_error(p, "unexpected character '" // c // "'", i)
       end select
       i = i + 1
    end do

    call push_token(p, TK_END, n + 1)

  end subroutine expr_tokenize

  !> Append a token to the token stream.
  !! @param[inout] p The parser state.
  !! @param[in] kind The kind of the token, one of the `TK_` constants.
  !! @param[in] pos Position of the token in the source string.
  subroutine push_token(p, kind, pos)
    type(parser_t), intent(inout) :: p
    integer, intent(in) :: kind
    integer, intent(in) :: pos

    p%n_tok = p%n_tok + 1
    p%tok(p%n_tok)%kind = kind
    p%tok(p%n_tok)%pos = pos

  end subroutine push_token

  !
  ! Recursive descent parser, emitting postfix opcodes as it goes.
  !
  !   expr    := term (('+' | '-') term)*
  !   term    := factor (('*' | '/') factor)*
  !   factor  := ('+' | '-') factor | power
  !   power   := primary ['^' factor]
  !   primary := NUM | IDENT | IDENT '(' expr (',' expr)* ')' | '(' expr ')'
  !

  !> Parse an additive expression.
  !! @param[inout] p The parser state.
  recursive subroutine parse_expr(p)
    type(parser_t), intent(inout) :: p

    call parse_term(p)
    do
       if (p%tok(p%cur)%kind .eq. TK_PLUS) then
          p%cur = p%cur + 1
          call parse_term(p)
          call emit(p, OP_ADD)
       else if (p%tok(p%cur)%kind .eq. TK_MINUS) then
          p%cur = p%cur + 1
          call parse_term(p)
          call emit(p, OP_SUB)
       else
          exit
       end if
    end do

  end subroutine parse_expr

  !> Parse a multiplicative expression.
  !! @param[inout] p The parser state.
  recursive subroutine parse_term(p)
    type(parser_t), intent(inout) :: p

    call parse_factor(p)
    do
       if (p%tok(p%cur)%kind .eq. TK_STAR) then
          p%cur = p%cur + 1
          call parse_factor(p)
          call emit(p, OP_MUL)
       else if (p%tok(p%cur)%kind .eq. TK_SLASH) then
          p%cur = p%cur + 1
          call parse_factor(p)
          call emit(p, OP_DIV)
       else
          exit
       end if
    end do

  end subroutine parse_term

  !> Parse a possibly signed factor. Unary minus binds looser than `^`, so
  !! `-x^2` is `-(x^2)`, as in Fortran.
  !! @param[inout] p The parser state.
  recursive subroutine parse_factor(p)
    type(parser_t), intent(inout) :: p
    integer :: before

    if (p%tok(p%cur)%kind .eq. TK_MINUS) then
       p%cur = p%cur + 1
       before = p%n_op
       call parse_factor(p)
       if (p%n_op .eq. before + 1 .and. p%op(p%n_op) .eq. OP_LIT) then
          ! Fold the sign into the literal, so that e.g. `y^-2` can use an
          ! integer power rather than the real valued one.
          p%lit(p%n_op) = -p%lit(p%n_op)
       else
          call emit(p, OP_NEG)
       end if
    else if (p%tok(p%cur)%kind .eq. TK_PLUS) then
       p%cur = p%cur + 1
       call parse_factor(p)
    else
       call parse_power(p)
    end if

  end subroutine parse_factor

  !> Parse an exponentiation. Right associative, so `2^3^2` is `2^(3^2)`.
  !! @param[inout] p The parser state.
  recursive subroutine parse_power(p)
    type(parser_t), intent(inout) :: p
    integer :: before

    call parse_primary(p)

    if (p%tok(p%cur)%kind .eq. TK_POW) then
       p%cur = p%cur + 1
       before = p%n_op
       call parse_factor(p)

       if (p%n_op .eq. before + 1 .and. p%op(p%n_op) .eq. OP_LIT) then
          ! Tested exactly, and not within a tolerance: an exponent such as
          ! `2.0000005` is not a whole number, and folding it would silently
          ! turn an undefined negative base power into a finite value.
          if (p%lit(p%n_op) .eq. anint(p%lit(p%n_op)) .and. &
               abs(p%lit(p%n_op)) .lt. 1.0e4_rp) then
             ! Whole exponent. Rewrite the literal push into a unary integer
             ! power, which unlike `**` with a real exponent is also defined
             ! for a negative base.
             p%op(p%n_op) = OP_IPOW
             p%depth = p%depth - 1
             return
          end if
       end if

       call emit(p, OP_POW)
    end if

  end subroutine parse_power

  !> Parse a literal, a symbol, a function call or a parenthesised expression.
  !! @param[inout] p The parser state.
  recursive subroutine parse_primary(p)
    type(parser_t), intent(inout) :: p
    integer :: k, op, nargs

    select case (p%tok(p%cur)%kind)
    case (TK_NUM)
       call emit(p, OP_LIT, p%tok(p%cur)%num)
       p%cur = p%cur + 1

    case (TK_LPAR)
       p%cur = p%cur + 1
       call parse_expr(p)
       if (p%tok(p%cur)%kind .ne. TK_RPAR) then
          call expr_error(p, "expected ')'", p%tok(p%cur)%pos)
       end if
       p%cur = p%cur + 1

    case (TK_IDENT)
       k = p%cur
       p%cur = p%cur + 1
       if (p%tok(p%cur)%kind .eq. TK_LPAR) then
          p%cur = p%cur + 1
          call parse_expr(p)
          nargs = 1
          do while (p%tok(p%cur)%kind .eq. TK_COMMA)
             p%cur = p%cur + 1
             call parse_expr(p)
             nargs = nargs + 1
          end do
          if (p%tok(p%cur)%kind .ne. TK_RPAR) then
             call expr_error(p, "expected ')' closing the argument list of '" &
                  // trim(p%tok(k)%name) // "'", p%tok(p%cur)%pos)
          end if
          p%cur = p%cur + 1

          op = func_op(trim(p%tok(k)%name), nargs)
          if (op .eq. -1) then
             call expr_error(p, "unknown function '" // &
                  trim(p%tok(k)%name) // "'", p%tok(k)%pos)
          else if (op .eq. -2) then
             call expr_error(p, "wrong number of arguments to '" // &
                  trim(p%tok(k)%name) // "'", p%tok(k)%pos)
          end if
          call emit(p, op)
       else
          call parse_symbol(p, k)
       end if

    case default
       call expr_error(p, "expected a value", p%tok(p%cur)%pos)
    end select

  end subroutine parse_primary

  !> Resolve a bare identifier into either a variable opcode or a literal.
  !! @details Coordinates and time variables take precedence over anything
  !! else. Everything that is left is looked up in the `neko_const_registry`
  !! and folded into a literal, so the value it had when the expression was
  !! compiled is the value that is used.
  !! @param[inout] p The parser state.
  !! @param[in] k Index of the identifier token to resolve.
  subroutine parse_symbol(p, k)
    type(parser_t), intent(inout) :: p
    integer, intent(in) :: k
    character(len=:), allocatable :: name
    real(kind=rp), pointer :: rval
    integer, pointer :: ival

    name = trim(p%tok(k)%name)

    select case (name)
    case ('x')
       call emit(p, OP_X)
    case ('y')
       call emit(p, OP_Y)
    case ('z')
       call emit(p, OP_Z)
    case ('t')
       call emit(p, OP_T)
       p%time_dependent = .true.
    case ('dt')
       call emit(p, OP_DT)
       p%time_dependent = .true.
    case ('pi')
       call emit(p, OP_LIT, 4.0_rp * atan(1.0_rp))
    case default
       if (neko_const_registry%real_scalar_exists(name)) then
          rval => neko_const_registry%get_real_scalar(name)
          call emit(p, OP_LIT, rval)
       else if (neko_const_registry%integer_scalar_exists(name)) then
          ival => neko_const_registry%get_integer_scalar(name)
          call emit(p, OP_LIT, real(ival, kind=rp))
       else
          call expr_error(p, "unknown symbol '" // name // "', it is " // &
               "neither a coordinate nor a constant declared under " // &
               "case.constants", p%tok(k)%pos)
       end if
    end select

    if (allocated(name)) deallocate(name)

  end subroutine parse_symbol

  !> Look up the opcode of a function.
  !! @param[in] name The name of the function.
  !! @param[in] nargs The number of arguments it was called with.
  !! @return The opcode, -1 if the name is unknown and -2 if the name is known
  !! but the number of arguments is wrong.
  function func_op(name, nargs) result(op)
    character(len=*), intent(in) :: name
    integer, intent(in) :: nargs
    integer :: op, arity

    arity = 1
    select case (name)
    case ('sin')
       op = OP_SIN
    case ('cos')
       op = OP_COS
    case ('tan')
       op = OP_TAN
    case ('asin')
       op = OP_ASIN
    case ('acos')
       op = OP_ACOS
    case ('atan')
       op = OP_ATAN
    case ('sinh')
       op = OP_SINH
    case ('cosh')
       op = OP_COSH
    case ('tanh')
       op = OP_TANH
    case ('exp')
       op = OP_EXP
    case ('log')
       op = OP_LOG
    case ('log10')
       op = OP_LOG10
    case ('sqrt')
       op = OP_SQRT
    case ('abs')
       op = OP_ABS
    case ('erf')
       op = OP_ERF
    case ('erfc')
       op = OP_ERFC
    case ('step')
       op = OP_STEP
    case ('atan2')
       op = OP_ATAN2
       arity = 2
    case ('min')
       op = OP_MIN
       arity = 2
    case ('max')
       op = OP_MAX
       arity = 2
    case ('mod')
       op = OP_MOD
       arity = 2
    case default
       op = -1
       return
    end select

    if (nargs .ne. arity) op = -2

  end function func_op

  !> Append an opcode, keeping track of how deep the evaluation stack gets.
  !! @param[inout] p The parser state.
  !! @param[in] op The opcode to append.
  !! @param[in] val The literal operand, only meaningful for `OP_LIT`.
  subroutine emit(p, op, val)
    type(parser_t), intent(inout) :: p
    integer, intent(in) :: op
    real(kind=rp), intent(in), optional :: val

    if (p%n_op .ge. size(p%op)) then
       call neko_error("Expression '" // p%src // "' is too complex")
    end if

    p%n_op = p%n_op + 1
    p%op(p%n_op) = op
    if (present(val)) then
       p%lit(p%n_op) = val
    else
       p%lit(p%n_op) = 0.0_rp
    end if

    if (op .lt. OP_UNARY) then
       p%depth = p%depth + 1
    else if (op .ge. OP_BINARY) then
       p%depth = p%depth - 1
    end if
    p%max_depth = max(p%max_depth, p%depth)

  end subroutine emit

  !> Abort with a message pointing at the offending part of the expression.
  !! @param[in] p The parser state.
  !! @param[in] msg What is wrong with the expression.
  !! @param[in] pos Position in the source string to point at.
  subroutine expr_error(p, msg, pos)
    type(parser_t), intent(in) :: p
    character(len=*), intent(in) :: msg
    integer, intent(in) :: pos
    character(len=16) :: buf

    write (buf, '(I0)') pos
    call neko_error("Invalid expression '" // p%src // "', at character " // &
         trim(buf) // ": " // msg)

  end subroutine expr_error

  !> True if `c` is a decimal digit.
  !! @param[in] c The character to test.
  pure function is_digit(c) result(res)
    character(len=1), intent(in) :: c
    logical :: res

    res = (c .ge. '0' .and. c .le. '9')

  end function is_digit

  !> True if `c` is a letter.
  !! @param[in] c The character to test.
  pure function is_alpha(c) result(res)
    character(len=1), intent(in) :: c
    logical :: res

    res = (c .ge. 'a' .and. c .le. 'z') .or. (c .ge. 'A' .and. c .le. 'Z')

  end function is_alpha

end module expression
