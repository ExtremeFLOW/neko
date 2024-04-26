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
!> Implements the `power_iterations_t` type.

module power_iterations
  use num_types, only : rp, dp, sp
  use json_module, only : json_file
  use simulation_component, only : simulation_component_t
  use field_registry, only : neko_field_registry
  use field, only : field_t
  use operators, only : curl
  use case, only : case_t
  use fld_file_output, only : fld_file_output_t
  use json_utils, only : json_get, json_get_or_default
  use field_writer, only : field_writer_t
  use math
  use comm, only: pe_rank
  implicit none
  private

  !> A simulation component that computes the power_iterations field.
  !! Added to the field registry as `omega_x`, `omega_y``, and `omega_z`.
  type, public, extends(simulation_component_t) :: power_iterations_t
     !> Neko case file
     class(case_t), pointer :: neko_case

     !> X velocity component.
     type(field_t), pointer :: u
     !> Y velocity component.
     type(field_t), pointer :: v
     !> Z velocity component.
     type(field_t), pointer :: w

     !> X power_iterations component.
     type(field_t), pointer :: u_old
     !> Y power_iterations component.
     type(field_t), pointer :: v_old
     !> Z power_iterations component.
     type(field_t), pointer :: w_old

     !> Size of the temporary array
     integer :: n
     !> Temporary array for working with the data
     real(kind=rp), dimension(:), allocatable :: tmp

     !> The previously used timestep
     real(kind=rp) :: t_old

     real(kind=rp) :: lambda
     real(kind=rp) :: lambda_old
     real(kind=rp) :: lambda_diff

   contains
     !> Constructor from json, wrapping the actual constructor.
     procedure, pass(this) :: init => power_iterations_init_from_json
     !> Actual constructor.
     procedure, pass(this) :: init_from_attributes => &
       power_iterations_init_from_attributes
     !> Destructor.
     procedure, pass(this) :: free => power_iterations_free
     !> Compute the power_iterations field.
     procedure, pass(this) :: compute_ => power_iterations_compute
  end type power_iterations_t

contains

  !> Constructor from json.
  subroutine power_iterations_init_from_json(this, json, case)
    class(power_iterations_t), intent(inout) :: this
    type(json_file), intent(inout) :: json
    class(case_t), intent(inout), target :: case

    call this%free()
    call this%init_base(json, case)

    call power_iterations_init_from_attributes(this)
  end subroutine power_iterations_init_from_json

  !> Actual constructor.
  subroutine power_iterations_init_from_attributes(this)
    class(power_iterations_t), intent(inout) :: this

    ! Add new fields needed by the simulation here
    call neko_field_registry%add_field(this%case%fluid%dm_Xh, "u_old",&
                                       ignore_existing=.true.)
    call neko_field_registry%add_field(this%case%fluid%dm_Xh, "v_old",&
                                       ignore_existing=.true.)
    call neko_field_registry%add_field(this%case%fluid%dm_Xh, "w_old",&
                                       ignore_existing=.true.)

    this%u => neko_field_registry%get_field("u")
    this%v => neko_field_registry%get_field("v")
    this%w => neko_field_registry%get_field("w")

    this%u_old => neko_field_registry%get_field("u_old")
    this%v_old => neko_field_registry%get_field("v_old")
    this%w_old => neko_field_registry%get_field("w_old")

    this%n = this%u%size()
    allocate(this%tmp(this%n))

  end subroutine power_iterations_init_from_attributes

  !> Destructor.
  subroutine power_iterations_free(this)
    class(power_iterations_t), intent(inout) :: this

    ! Deallocate local arrays
    deallocate(this%tmp)

    ! Deallocate fields
    call this%free_base()
  end subroutine power_iterations_free

  !> Compute the power_iterations field.
  !! @param t The time value.
  !! @param tstep The current time-step
  subroutine power_iterations_compute(this, t, tstep)
    class(power_iterations_t), intent(inout) :: this
    real(kind=rp), intent(in) :: t
    integer, intent(in) :: tstep

    ! real(kind=rp) :: norm
    ! real(kind=rp) :: volvm1
    ! real(kind=rp) :: cnht_sv
    ! real(kind=rp) :: f1

    ! integer :: i

    ! ! Assign the scalar constants
    ! volvm1=this%neko_case%fluid%c_Xh%volume

    ! ! https://github.com/KTH-Nek5000/KTH_Toolbox/blob/b2b7a97b4ba1a56a75dcbd2c1703a595cda9a850/utility/cnht/cnht_tools.f#L51C1-L54C40
    ! ! This is the default value for cnht_sv in Nek5000. It was user defined in
    ! ! the par file format. We should probably add this to the json file.
    ! cnht_sv = 0.5

    ! f1=cnht_sv/volvm1

    ! associate (u => this%u%x, &
    !            v => this%v%x, &
    !            w => this%w%x, &
    !            u_old => this%u_old%x, &
    !            v_old => this%v_old%x, &
    !            w_old => this%w_old%x, &
    !            wt => this%neko_case%fluid%c_Xh%B & ! Mass matrix
    !            )

    !   do i=1, 100

    !      ! Calculate the new lambda value
    !      ! lambda = U_n^T * U_old
    !      ! U_old = U_{n-1} / ||U_{n-1}||
    !      this%lambda_old = this%lambda

    !      !  ! Alternative 1: Compute the dot product of the velocity fields.
    !      !  call vdot3(this%tmp, u, v, w, u_old, v_old, w_old, this%n)
    !      !  this%lambda = glsum(this%tmp, this%n)

    !      !  call vdot3(this%tmp, u_old, v_old, w_old, u_old, v_old, w_old, this%n)
    !      !  this%lambda = this%lambda / glsum(this%tmp, this%n)

    !      ! Alternative 2: Compute the integral of the velocity fields.
    !      call vdot3(this%tmp, u, v, w, u_old, v_old, w_old, this%n)
    !      this%lambda = f1 * glsc2(this%tmp, wt, this%n)

    !      call vdot3(this%tmp, u_old, v_old, w_old, u_old, v_old, w_old, this%n)
    !      norm = f1 * glsc2(this%tmp, wt, this%n)
    !      this%lambda = this%lambda / norm


    !      ! Calculate the difference between the new and old lambda values.
    !      this%lambda_diff = abs(this%lambda - this%lambda_old)

    !      if (this%lambda_diff .lt. 1.0e-18_rp) then
    !         exit
    !      end if


    !      ! Normalize and store the pertubation fields.
    !      call copy(u_old, u, this%n)
    !      call copy(v_old, v, this%n)
    !      call copy(w_old, w, this%n)

    !      ! Normalize the pertubation fields.
    !      call vdot3(this%tmp, u_old, v_old, w_old, u_old, v_old, w_old, this%n)
    !      norm = f1 * glsc2(this%tmp, wt, this%n)
    !      norm = sqrt(norm)

    !      call cmult(u_old, 1.0_rp / norm, this%n)
    !      call cmult(v_old, 1.0_rp / norm, this%n)
    !      call cmult(w_old, 1.0_rp / norm, this%n)

    !   end do

    ! end associate

    ! if (pe_rank .eq. 0) then
    !    print *, "lambda: ", this%lambda
    !    print *, "lambda_diff: ", this%lambda_diff
    !    print *, "Number of iterations: ", i
    ! end if

    ! ======================================================================== !
    ! Let's try to do a line by line translation instead of the above code.
    ! ======================================================================== !
    ! Translation:

    ! Scratch space
    real(kind=rp), dimension(:), allocatable :: TA1, TA2, TA3, TAT

    ! Local variables
    integer :: itmp, LASTEP
    real(kind=rp) :: lnorm, grth_old, ltim
    real(kind=rp) :: V, dt

    ! ! Pointers to the fields
    real(kind=rp), allocatable :: VXP(:), VYP(:), VZP(:)
    ! real(kind=rp), pointer :: BM1(:)

    ! Parameters that should be global
    real(kind=rp) :: pwit_l2n = 1 !< Vector initial norm (Set to 1.0 in the Nek5000 code)
    real(kind=rp) :: pwit_grw = 0 !< Initial growth rate
    ! real(kind=rp), allocatable :: pwit_vx(:), pwit_vy(:), pwit_vz(:)
    real(kind=rp) :: cnht_sv

    ! tstpr parameters
    real(kind=rp) :: tstpr_tol = 1.0e-5_rp !< Tolerance for the difference between the new and old lambda values
    integer :: tstpr_cmax = 10 !< Maximum number of iterations
    integer :: tstpr_vstep = 1 !< Number of steps between checkpointing

    ! Allocate temporary arrays
    allocate(TA1(this%n), TA2(this%n), TA3(this%n))
    allocate(VXP(this%n), VYP(this%n), VZP(this%n))
    ! allocate(pwit_vx(this%n), pwit_vy(this%n), pwit_vz(this%n))


    ! call copy(pwit_vx, VXP, this%n)
    ! call copy(pwit_vy, VYP, this%n)
    ! call copy(pwit_vz, VZP, this%n)

    cnht_sv = 0.5_rp
    V = this%neko_case%fluid%c_Xh%volume
    dt = t - this%t_old
    this%t_old = t

    ! Assign pointers
    associate( BM1 => this%neko_case%fluid%c_Xh%B, &
               pwit_vx => this%u_old%x, &
               pwit_vy => this%v_old%x, &
               pwit_vz => this%w_old%x, &
               VXP => this%u%x, &
               VYP => this%v%x, &
               VZP => this%w%x)

      ! do tstpr_vstep = 1, tstpr_cmax

      ! ! timing
      ! ltim=dnekclock()

      ! normalise vector
      ! lnorm = cnht_glsc2_wt(VXP,VYP,VZP,TP,VXP,VYP,VZP,TP,BM1)
      call vdot3(this%tmp, VXP, VYP, VZP, VXP, VYP, VZP, this%n)
      lnorm = cnht_sv / V * glsc2(this%tmp, BM1, this%n) !> This is wrong. We are missing the
      ! scaling factor "f1 = cnht_sv/volvm1" from the Nek5000 code.
      lnorm = sqrt(lnorm / pwit_l2n)

      ! call cnht_opcmult (VXP,VYP,VZP,TP,lnorm)
      call cmult(VXP, 1.0_rp / lnorm, this%n)
      call cmult(VYP, 1.0_rp / lnorm, this%n)
      call cmult(VZP, 1.0_rp / lnorm, this%n)

      This%lambda_old = this%lambda
      this%lambda = lnorm

      ! if (tstpr_pr.ne.0) then
      !    ! normalise pressure ???????????????????????????
      !    itmp = nx2*ny2*nz2*nelv
      !    call cmult(PRP,lnorm,itmp)
      ! endif

      ! ! make sure the velocity and temperature fields are continuous at
      ! ! element faces and edges
      ! call tstpr_dssum

      ! compare current and prevoius growth rate
      grth_old = pwit_grw
      pwit_grw = lnorm
      grth_old = pwit_grw - grth_old

      ! get L2 norm of the update
      ! call cnht_opsub3 (TA1,TA2,TA3,TAT,pwit_vx,pwit_vy,pwit_vz,pwit_t,
      ! $ VXP,VYP,VZP,TP)
      call sub3(TA1, pwit_vx, VXP, this%n)
      call sub3(TA2, pwit_vy, VYP, this%n)
      call sub3(TA3, pwit_vz, VZP, this%n)
      ! lnorm = cnht_glsc2_wt(TA1,TA2,TA3,TAT,TA1,TA2,TA3,TAT,BM1)
      call vdot3(this%tmp, TA1, TA2, TA3, TA1, TA2, TA3, this%n)
      lnorm = cnht_sv / V * glsc2(this%tmp, BM1, this%n)
      lnorm = sqrt(lnorm)

      ! log stamp
      ! call mntr_log(pwit_id,lp_prd,'POWER ITERATIONS: convergence')
      ! call mntr_logr(pwit_id,lp_prd,'||V-V_old|| = ',lnorm)
      ! call mntr_logr(pwit_id,lp_prd,'Growth ', pwit_grw)
      if (pe_rank .eq. 0) then
         print *, 'POWER ITERATIONS: convergence'
         print *, '||V-V_old|| = ', lnorm
         print *, 'Growth ', pwit_grw
         print *, 'Eigen:', 1/dt * log(this%lambda / this%lambda_old)
      end if

      ! itmp = 0
      ! if (IFHEAT) itmp = 1

      ! !write down current field
      ! call outpost2(VXP,VYP,VZP,PRP,TP,itmp,'PWI')

      ! ! write down field difference
      ! call outpost2(TA1,TA2,TA3,PRP,TAT,itmp,'VDF')

      ! check convergence
      LASTEP = 0
      if(lnorm.lt.tstpr_tol.and.grth_old.lt.tstpr_tol) then
         !  call mntr_log(pwit_id,lp_prd,'Reached stopping criteria')
         if (pe_rank .eq. 0) then
            print *, 'Reached stopping criteria'
         end if
         ! mark the last step
         LASTEP = 1
      else
         ! save current vector and restart stepper
         !  call cnht_opcopy(pwit_vx,pwit_vy,pwit_vz,pwit_t,VXP,VYP,VZP,TP)
         call copy(pwit_vx, VXP, this%n)
         call copy(pwit_vy, VYP, this%n)
         call copy(pwit_vz, VZP, this%n)
      endif

      ! save checkpoint
      if (LASTEP.eq.1.or.tstpr_cmax.eq.tstpr_vstep) then
         ! call stepper_write()

         ! mark the last step
         LASTEP = 1
      endif

      ! timing
      ! ltim = dnekclock() - ltim
      ! call mntr_tmr_add(pwit_tmr_evl_id,1,ltim)

      if (LASTEP.eq.1) then
         ! final log stamp
         !  call mntr_log(pwit_id,lp_prd,'POWER ITERATIONS finalised')
         !  call mntr_logr(pwit_id,lp_prd,'||V-V_old|| = ',lnorm)
         !  call mntr_logr(pwit_id,lp_prd,'Growth ',pwit_grw)
         if (pe_rank .eq. 0) then
            print *, 'POWER ITERATIONS finalised'
            print *, '||V-V_old|| = ', lnorm
            print *, 'Growth ', pwit_grw
         end if
         !  exit
      endif

      ! end do

    end associate

    deallocate(TA1, TA2, TA3)

  end subroutine power_iterations_compute

end module power_iterations


#ifdef ONLY_FOR_REFERENCE
! The following submodules are only here for reference. Our goal is to translate
! these pieces of codes into the Neko format.


! ========================================================================== !
! Power iterations module from Nek5000
!
! This should be used as reference when setting up the
! power_iterations_compute.
! Our expectation is that this will give some insights in how to actually do
! the power iterations.
! ========================================================================== !
submodule (power_iterations) pwit_template
contains

  !> @file pwit.f
  !! @ingroup pwit
  !! @brief Set of subroutines to perform power iterations within time stepper
  !! @author Adam Peplinski
  !! @date Mar 7, 2016
  !=======================================================================
  !> @brief Register power iteration module
  !! @ingroup pwit
  !! @note This interface is called by @ref tstpr_register
  subroutine stepper_register()
    implicit none

    include 'SIZE'
    include 'INPUT'
    include 'FRAMELP'
    include 'TSTPRD'
    include 'PWITD'

    ! local variables
    integer lpmid, il
    real ltim
    character*2 str

    ! functions
    real dnekclock
    !-----------------------------------------------------------------------
    ! timing
    ltim = dnekclock()

    ! check if the current module was already registered
    call mntr_mod_is_name_reg(lpmid,pwit_name)
    if (lpmid.gt.0) then
       call mntr_warn(lpmid,
       $ 'module ['//trim(pwit_name)//'] already registered')
       return
    endif

    ! find parent module
    call mntr_mod_is_name_reg(lpmid,tstpr_name)
    if (lpmid.le.0) then
       lpmid = 1
       call mntr_abort(lpmid,
       $ 'parent module ['//trim(tstpr_name)//'] not registered')
    endif

    ! register module
    call mntr_mod_reg(pwit_id,lpmid,pwit_name,
    $ 'Power iterations for time stepper')

    ! register timers
    ! initialisation
    call mntr_tmr_reg(pwit_tmr_ini_id,tstpr_tmr_ini_id,pwit_id,
    $ 'PWIT_INI','Power iteration initialisation time',.true.)
    ! submodule operation
  contains
    call mntr_tmr_reg(pwit_tmr_evl_id,tstpr_tmr_evl_id,pwit_id,
    $ 'PWIT_EVL','Power iteration evolution time',.true.)

    ! register and set active section
    call rprm_sec_reg(pwit_sec_id,pwit_id,'_'//adjustl(pwit_name),
    $ 'Runtime paramere section for power iteration module')
    call rprm_sec_set_act(.true.,pwit_sec_id)

    ! register parameters
    call rprm_rp_reg(pwit_l2n_id,pwit_sec_id,'L2N',
    $ 'Vector initial norm',rpar_real,0,1.0,.false.,' ')

    ! set initialisation flag
    pwit_ifinit=.false.

    ! timing
    ltim = dnekclock() - ltim
    call mntr_tmr_add(pwit_tmr_ini_id,1,ltim)

    return
  end subroutine stepper_register
  !=======================================================================
  !> @brief Initilise power iteration module
  !! @ingroup pwit
  !! @note This interface is called by @ref tstpr_init
  subroutine stepper_init()
    implicit none

    include 'SIZE'
    include 'SOLN' ! V[XYZ]P, TP
    include 'MASS' ! BM1
    include 'FRAMELP'
    include 'TSTPRD'
    include 'PWITD'

    ! local variables
    integer itmp, il, set_in
    real rtmp, ltim, lnorm
    logical ltmp
    character*20 ctmp

    ! to get checkpoint runtime parameters
    integer ierr, lmid, lsid, lrpid

    ! functions
    real dnekclock, cnht_glsc2_wt
    logical chkpts_is_initialised
    !-----------------------------------------------------------------------
    ! check if the module was already initialised
    if (pwit_ifinit) then
       call mntr_warn(pwit_id,
       $ 'module ['//trim(pwit_name)//'] already initialised.')
       return
    endif

    ! timing
    ltim = dnekclock()

    ! get runtime parameters
    call rprm_rp_get(itmp,rtmp,ltmp,ctmp,pwit_l2n_id,rpar_real)
    pwit_l2n = rtmp

    ! check the restart flag
    ! check if checkpointing module was registered and take parameters
    ierr = 0
    call mntr_mod_is_name_reg(lmid,'CHKPT')
    if (lmid.gt.0) then
       call rprm_sec_is_name_reg(lsid,lmid,'_CHKPT')
       if (lsid.gt.0) then
          ! restart flag
          call rprm_rp_is_name_reg(lrpid,lsid,'READCHKPT',rpar_log)
          if (lrpid.gt.0) then
             call rprm_rp_get(itmp,rtmp,ltmp,ctmp,lrpid,rpar_log)
             pwit_ifrst = ltmp
          else
             ierr = 1
             goto 30
          endif
          if (pwit_ifrst) then
             ! checkpoint set number
             call rprm_rp_is_name_reg(lrpid,lsid,'CHKPFNUMBER',
             $ rpar_int)
             if (lrpid.gt.0) then
                call rprm_rp_get(itmp,rtmp,ltmp,ctmp,lrpid,rpar_int)
                pwit_fnum = itmp
             else
                ierr = 1
                goto 30
             endif
          endif
       else
          ierr = 1
       endif
    else
       ierr = 1
    endif

30  continue

    ! check for errors
    call mntr_check_abort(pwit_id,ierr,
    $ 'Error reading checkpoint parameters')

    ! read checkpoint file
    if (pwit_ifrst) then
       if(.not.chkpts_is_initialised()) call mntr_abort(pwit_id,
       $ 'Checkpointing module not initialised')
       set_in = pwit_fnum -1
       call stepper_read(set_in)
    endif

    ! initial growth rate
    pwit_grw = 0.0

    ! normalise vector
    lnorm = cnht_glsc2_wt(VXP,VYP,VZP,TP,VXP,VYP,VZP,TP,BM1)
    lnorm = sqrt(pwit_l2n/lnorm)
    call cnht_opcmult (VXP,VYP,VZP,TP,lnorm)

    if (tstpr_pr.ne.0) then
       ! normalise pressure ???????????????????????????
       itmp = nx2*ny2*nz2*nelv
       call cmult(PRP,lnorm,itmp)
    endif

    ! make sure the velocity and temperature fields are continuous at
    ! element faces and edges
    call tstpr_dssum

    ! save intial vector
    call cnht_opcopy (pwit_vx,pwit_vy,pwit_vz,pwit_t,VXP,VYP,VZP,TP)

    ! stamp log file
    call mntr_log(pwit_id,lp_prd,'POWER ITERATIONS initialised')
    call mntr_logr(pwit_id,lp_prd,'L2NORM = ',pwit_l2n)

    ! everything is initialised
    pwit_ifinit=.true.

    ! timing
    ltim = dnekclock() - ltim
    call mntr_tmr_add(pwit_tmr_ini_id,1,ltim)

    return
  end subroutine stepper_init
  !=======================================================================
  !> @brief Check if module was initialised
  !! @ingroup pwit
  !! @return stepper_is_initialised
  logical function stepper_is_initialised()
    implicit none

    include 'SIZE'
    include 'PWITD'
    !-----------------------------------------------------------------------
    stepper_is_initialised = pwit_ifinit

    return
  end function stepper_is_initialised
  !=======================================================================
  !> @brief Renormalise vector and check convergence.
  !! @ingroup pwit
  !! @note This interface is defined in @ref tstpr_main
  !! @remarks This routine uses global scratch space SCRUZ
  subroutine stepper_vsolve
    implicit none

    include 'SIZE' ! NIO
    include 'TSTEP' ! TIME, LASTEP, NSTEPS
    include 'INPUT' ! IFHEAT
    include 'MASS' ! BM1
    include 'SOLN' ! V[XYZ]P, TP
    include 'FRAMELP'
    include 'TSTPRD'
    include 'PWITD'

    ! scratch space
    real TA1 (LPX1*LPY1*LPZ1*LELV), TA2 (LPX1*LPY1*LPZ1*LELV),
    $ TA3 (LPX1*LPY1*LPZ1*LELV), TAT (LPX1*LPY1*LPZ1*LELT)
    COMMON /SCRUZ/ TA1, TA2, TA3, TAT

    ! local variables
    integer itmp
    real lnorm, grth_old, ltim

    ! functions
    real dnekclock, cnht_glsc2_wt
    !-----------------------------------------------------------------------
    ! timing
    ltim=dnekclock()

    ! normalise vector
    lnorm = cnht_glsc2_wt(VXP,VYP,VZP,TP,VXP,VYP,VZP,TP,BM1)
    lnorm = sqrt(pwit_l2n/lnorm)
    call cnht_opcmult (VXP,VYP,VZP,TP,lnorm)

    if (tstpr_pr.ne.0) then
       ! normalise pressure ???????????????????????????
       itmp = nx2*ny2*nz2*nelv
       call cmult(PRP,lnorm,itmp)
    endif

    ! make sure the velocity and temperature fields are continuous at
    ! element faces and edges
    call tstpr_dssum

    ! compare current and prevoius growth rate
    grth_old = pwit_grw
    pwit_grw = 1.0/lnorm
    grth_old = pwit_grw - grth_old

    ! get L2 norm of the update
    call cnht_opsub3 (TA1,TA2,TA3,TAT,pwit_vx,pwit_vy,pwit_vz,pwit_t,
    $ VXP,VYP,VZP,TP)
    lnorm = cnht_glsc2_wt(TA1,TA2,TA3,TAT,TA1,TA2,TA3,TAT,BM1)
    lnorm = sqrt(lnorm)

    ! log stamp
    call mntr_log(pwit_id,lp_prd,'POWER ITERATIONS: convergence')
    call mntr_logr(pwit_id,lp_prd,'||V-V_old|| = ',lnorm)
    call mntr_logr(pwit_id,lp_prd,'Growth ',pwit_grw)

    itmp = 0
    if (IFHEAT) itmp = 1

    !write down current field
    call outpost2(VXP,VYP,VZP,PRP,TP,itmp,'PWI')

    ! write down field difference
    call outpost2(TA1,TA2,TA3,PRP,TAT,itmp,'VDF')

    ! check convergence
    if(lnorm.lt.tstpr_tol.and.grth_old.lt.tstpr_tol) then
       call mntr_log(pwit_id,lp_prd,'Reached stopping criteria')
       ! mark the last step
       LASTEP = 1
    else
       ! save current vector and restart stepper
       call cnht_opcopy(pwit_vx,pwit_vy,pwit_vz,pwit_t,VXP,VYP,VZP,TP)
    endif

    ! save checkpoint
    if (LASTEP.eq.1.or.tstpr_cmax.eq.tstpr_vstep) then
       call stepper_write()

       ! mark the last step
       LASTEP = 1
    endif

    ! timing
    ltim = dnekclock() - ltim
    call mntr_tmr_add(pwit_tmr_evl_id,1,ltim)

    if (LASTEP.eq.1) then
       ! final log stamp
       call mntr_log(pwit_id,lp_prd,'POWER ITERATIONS finalised')
       call mntr_logr(pwit_id,lp_prd,'||V-V_old|| = ',lnorm)
       call mntr_logr(pwit_id,lp_prd,'Growth ',pwit_grw)
    endif

    return
  end subroutine stepper_vsolve
  !=======================================================================
  !> @brief Read restart files
  !! @ingroup pwit
  !! @param[in]  set_in  restart set number
  subroutine stepper_read(set_in)
    implicit none

    include 'SIZE' ! NIO
    include 'TSTEP' ! TIME, LASTEP, NSTEPS
    include 'INPUT' ! IFMVBD, IFREGUO
    include 'FRAMELP'
    include 'CHKPTD'
    include 'CHKPTMSD'
    include 'PWITD'

    ! argument list
    integer set_in

    ! local variables
    integer ifile, step_cnt, fnum
    character*132 fname(chkptms_fmax)
    logical ifreguol
    !-----------------------------------------------------------------------
    ! no regular mesh
    ifreguol= IFREGUO
    IFREGUO = .false.

    call mntr_log(pwit_id,lp_inf,'Reading checkpoint snapshot')

    ! initialise I/O data
    call io_init

    ! get set of file names in the snapshot
    ifile = 1
    call chkptms_set_name(fname, fnum, set_in, ifile)

    ! read files
    call chkptms_restart_read(fname, fnum)

    ! put parameters back
    IFREGUO = ifreguol

    return
  end subroutine stepper_read
  !=======================================================================
  !> @brief Write restart files
  !! @ingroup pwit
  subroutine stepper_write
    implicit none

    include 'SIZE' ! NIO
    include 'TSTEP' ! TIME, LASTEP, NSTEPS
    include 'INPUT' ! IFMVBD, IFREGUO
    include 'FRAMELP'
    include 'CHKPTD'
    include 'CHKPTMSD'
    include 'PWITD'

    ! local variables
    integer ifile, step_cnt, set_out, fnum
    character*132 fname(chkptms_fmax)
    logical ifcoord
    logical ifreguol
    !-----------------------------------------------------------------------
    ! no regular mesh
    ifreguol= IFREGUO
    IFREGUO = .false.

    call mntr_log(pwit_id,lp_inf,'Writing checkpoint snapshot')

    ! initialise I/O data
    call io_init

    ! get set of file names in the snapshot
    ifile = 1
    call chkpt_get_fset(step_cnt, set_out)
    call chkptms_set_name(fname, fnum, set_out, ifile)

    ifcoord = .true.
    ! write down files
    call chkptms_restart_write(fname, fnum, ifcoord)

    ! put parameters back
    IFREGUO = ifreguol

    return
  end subroutine stepper_write
end submodule pwit_template

! ========================================================================== !
! Time stepper module from Nek5000
!
! This should be used as reference when setting up the
! power_iterations_compute. Our expectation is that this will give some
! insighits on how to call the power iterations module.
! ========================================================================== !
submodule (power_iterations) tstpr_template
contains
  !> @file tstpr.f
  !! @ingroup tstpr
  !! @brief Set of subroutines to use time steppers for e.g. power
  !!    iterations or solution of eigenvalue problem with Arnoldi algorithm
  !! @author Adam Peplinski
  !! @date Mar 7, 2016
  !=======================================================================
  !> @brief Register time stepper module
  !! @ingroup tstpr
  !! @note This routine should be called in frame_usr_register
  subroutine tstpr_register()
    implicit none

    include 'SIZE'
    include 'INPUT'
    include 'FRAMELP'
    include 'TSTPRD'

    ! local variables
    integer lpmid, il
    real ltim
    character*2 str

    ! functions
    real dnekclock
    !-----------------------------------------------------------------------
    ! timing
    ltim = dnekclock()

    ! check if the current module was already registered
    call mntr_mod_is_name_reg(lpmid,tstpr_name)
    if (lpmid.gt.0) then
       call mntr_warn(lpmid,
       $ 'module ['//trim(tstpr_name)//'] already registered')
       return
    endif

    ! check if conjugated heat transfer module was registered
    call mntr_mod_is_name_reg(lpmid,'CNHT')
    if (lpmid.gt.0) then
       call mntr_warn(lpmid,
       $ 'module ['//'CNHT'//'] already registered')
    else
       call cnht_register()
    endif

    ! find parent module
    call mntr_mod_is_name_reg(lpmid,'FRAME')
    if (lpmid.le.0) then
       lpmid = 1
       call mntr_abort(lpmid,
       $ 'parent module ['//'FRAME'//'] not registered')
    endif

    ! register module
    call mntr_mod_reg(tstpr_id,lpmid,tstpr_name,
    $ 'Time stepper')

    ! register timers
    call mntr_tmr_is_name_reg(lpmid,'FRM_TOT')
    ! total time
    call mntr_tmr_reg(tstpr_tmr_tot_id,lpmid,tstpr_id,
    $ 'TSTPR_TOT','Time stepper total time',.false.)
    lpmid = tstpr_tmr_tot_id
    ! initialisation
    call mntr_tmr_reg(tstpr_tmr_ini_id,lpmid,tstpr_id,
    $ 'TSTPR_INI','Time stepper initialisation time',.true.)
    ! submodule operation
    call mntr_tmr_reg(tstpr_tmr_evl_id,lpmid,tstpr_id,
    $ 'TSTPR_EVL','Time stepper evolution time',.true.)

    ! register and set active section
    call rprm_sec_reg(tstpr_sec_id,tstpr_id,'_'//adjustl(tstpr_name),
    $ 'Runtime paramere section for time stepper module')
    call rprm_sec_set_act(.true.,tstpr_sec_id)

    ! register parameters
    call rprm_rp_reg(tstpr_mode_id,tstpr_sec_id,'MODE',
    $ 'Simulation mode',rpar_str,10,0.0,.false.,'DIR')

    call rprm_rp_reg(tstpr_step_id,tstpr_sec_id,'STEPS',
    $ 'Length of stepper phase',rpar_int,40,0.0,.false.,' ')

    call rprm_rp_reg(tstpr_cmax_id,tstpr_sec_id,'MAXCYC',
    $ 'Max number of stepper cycles',rpar_int,10,0.0,.false.,' ')

    call rprm_rp_reg(tstpr_tol_id,tstpr_sec_id,'TOL',
    $ 'Convergence threshold',rpar_real,0,1.0d-6,.false.,' ')

    ! place for submodule registration
    ! register arnoldi or power iterations
    call stepper_register()

    ! set initialisation flag
    tstpr_ifinit=.false.

    ! timing
    ltim = dnekclock() - ltim
    call mntr_tmr_add(tstpr_tmr_ini_id,1,ltim)

    return
  end subroutine tstpr_register
  !=======================================================================
  !> @brief Initilise time stepper module
  !! @ingroup tstpr
  !! @note This routine should be called in frame_usr_init
  subroutine tstpr_init()
    implicit none

    include 'SIZE'
    include 'FRAMELP'
    include 'TSTEP'
    include 'INPUT'
    include 'MASS'
    include 'SOLN'
    include 'ADJOINT'
    include 'TSTPRD'

    ! local variables
    integer itmp, il
    real rtmp, ltim
    logical ltmp
    character*20 ctmp

    ! functions
    real dnekclock, cnht_glsc2_wt
    logical cnht_is_initialised
    !-----------------------------------------------------------------------
    ! check if the module was already initialised
    if (tstpr_ifinit) then
       call mntr_warn(tstpr_id,
       $ 'module ['//trim(tstpr_name)//'] already initiaised.')
       return
    endif

    ! timing
    ltim = dnekclock()

    ! intialise conjugated heat transfer
    if (.not.cnht_is_initialised()) call cnht_init

    ! get runtime parameters
    call rprm_rp_get(itmp,rtmp,ltmp,ctmp,tstpr_mode_id,rpar_str)
    if (trim(ctmp).eq.'DIR') then
       tstpr_mode = 1
    else if (trim(ctmp).eq.'ADJ') then
       tstpr_mode = 2
    else if (trim(ctmp).eq.'OIC') then
       tstpr_mode = 3
    else
       call mntr_abort(tstpr_id,
       $ 'wrong simulation mode; possible values: DIR, ADJ, OIC')
    endif

    call rprm_rp_get(itmp,rtmp,ltmp,ctmp,tstpr_step_id,rpar_int)
    tstpr_step = itmp

    call rprm_rp_get(itmp,rtmp,ltmp,ctmp,tstpr_cmax_id,rpar_int)
    tstpr_cmax = itmp

    call rprm_rp_get(itmp,rtmp,ltmp,ctmp,tstpr_tol_id,rpar_real)
    tstpr_tol = rtmp

    ! check simulation parameters
    if (.not.IFTRAN) call mntr_abort(tstpr_id,
    $ 'time stepper requres transient simulation; IFTRAN=.T.')

    if (NSTEPS.eq.0) call mntr_abort(tstpr_id,
    $ 'time stepper requres NSTEPS>0')

    if (PARAM(12).ge.0) call mntr_abort(tstpr_id,
    $ 'time stepper assumes constant dt')

    if (.not.IFPERT) call mntr_abort(tstpr_id,
    $ 'time stepper has to be run in perturbation mode')

    if (IFBASE) call mntr_abort(tstpr_id,
    $ 'time stepper assumes constatnt base flow')

    if (NPERT.ne.1) call mntr_abort(tstpr_id,
    $ 'time stepper requires NPERT=1')

    if (IFHEAT.and.tstpr_ht.ne.1) call mntr_abort(tstpr_id,
    $ 'time stepper requires tstpr_ht=1 for temperature evaluation')

    ! initialise cycle counters
    tstpr_istep = 0
    tstpr_vstep = 0

    ! vector length
    tstpr_nv = NX1*NY1*NZ1*NELV ! velocity single component
    if (IFHEAT) then !temperature
       tstpr_nt = NX1*NY1*NZ1*NELT
    else
       tstpr_nt = 0
    endif
    tstpr_np = NX2*NY2*NZ2*NELV ! presure

    ! place for submodule initialisation
    ! arnoldi or power iterations
    call stepper_init

    ! zero presure
    if(tstpr_pr.eq.0) call rzero(PRP,tstpr_np)

    ! set initial time
    TIME=0.0

    ! make sure NSTEPS is bigger than the possible number of iterations
    ! in time stepper phase; multiplication by 2 for OIC
    NSTEPS = max(NSTEPS,tstpr_step*tstpr_cmax*2+10)

    IFADJ = .FALSE.
    if (tstpr_mode.eq.2) then
       ! Is it adjoint mode
       IFADJ = .TRUE.
    elseif (tstpr_mode.eq.3) then
       ! If it is optimal initial condition save initial L2 norm
       tstpr_L2ini = cnht_glsc2_wt(VXP,VYP,VZP,TP,VXP,VYP,VZP,TP,BM1)

       if (tstpr_L2ini.eq.0.0) call mntr_abort(tstpr_id,
       $ 'tstpr_init, tstpr_L2ini = 0')

       call mntr_log(tstpr_id,lp_prd,
       $ 'Optimal initial condition; direct phase start')
    endif

    ! set cpfld for conjugated heat transfer
    if (IFHEAT) call cnht_cpfld_set

    ! everything is initialised
    tstpr_ifinit=.true.

    ! timing
    ltim = dnekclock() - ltim
    call mntr_tmr_add(tstpr_tmr_ini_id,1,ltim)

    return
  end subroutine tstpr_init
  !=======================================================================
  !> @brief Check if module was initialised
  !! @ingroup tstpr
  !! @return tstpr_is_initialised
  logical function tstpr_is_initialised()
    implicit none

    include 'SIZE'
    include 'TSTPRD'
    !-----------------------------------------------------------------------
    tstpr_is_initialised = tstpr_ifinit

    return
  end function tstpr_is_initialised
  !=======================================================================
  !> @brief Control time stepper after every nek5000 step and call suitable
  !! stepper_vsolve of required submodule
  !! @ingroup tstpr
  subroutine tstpr_main()
    implicit none

    include 'SIZE' ! NIO
    include 'TSTEP' ! ISTEP, TIME
    include 'INPUT' ! IFHEAT, IF3D
    include 'MASS' ! BM1
    include 'SOLN' ! V[XYZ]P, PRP, TP, VMULT, V?MASK
    include 'ADJOINT' ! IFADJ
    include 'FRAMELP'
    include 'TSTPRD'

    ! global comunication in nekton
    integer nidd,npp,nekcomm,nekgroup,nekreal
    common /nekmpi/ nidd,npp,nekcomm,nekgroup,nekreal

    ! local variables
    real grw ! growth rate
    real ltim ! timing

    ! functions
    real dnekclock, cnht_glsc2_wt
    !-----------------------------------------------------------------------
    if (ISTEP.eq.0) return

    ! step counting
    tstpr_istep = tstpr_istep + 1

    ! stepper phase end
    if (mod(tstpr_istep,tstpr_step).eq.0) then
       ! timing
       ltim = dnekclock()
       ! check for the calculation mode
       if (tstpr_mode.eq.3.and.(.not.IFADJ)) then
          ! optimal initial condition

          call mntr_log(tstpr_id,lp_prd,
          $ 'Optimal initial condition; adjoint phase start')

          IFADJ = .TRUE.

          ! iteration count
          tstpr_istep = 0

          ! set time and iteration number
          TIME=0.0
          ISTEP=0

          ! get L2 norm after direct phase
          tstpr_L2dir = cnht_glsc2_wt(VXP,VYP,VZP,TP,
          $ VXP,VYP,VZP,TP,BM1)
          ! normalise vector
          grw = sqrt(tstpr_L2ini/tstpr_L2dir)
          call cnht_opcmult (VXP,VYP,VZP,TP,grw)

          if (tstpr_pr.eq.0) then
             ! zero presure
             call rzero(PRP,tstpr_np)
          else
             ! normalise pressure ???????????????????????????
             call cmult(PRP,grw,tstpr_np)
          endif

          ! set cpfld for conjugated heat transfer
          if (IFHEAT) call cnht_cpfld_set
       else
          !stepper phase counting
          tstpr_istep = 0
          tstpr_vstep = tstpr_vstep +1

          call mntr_logi(tstpr_id,lp_prd,'Finished stepper phase:',
          $ tstpr_vstep)

          if (tstpr_mode.eq.3) then
             ! optimal initial condition
             call mntr_log(tstpr_id,lp_prd,
             $ 'Optimal initial condition; rescaling solution')

             ! get L2 norm after direct phase
             tstpr_L2adj = cnht_glsc2_wt(VXP,VYP,VZP,TP,
             $ VXP,VYP,VZP,TP,BM1)
             ! normalise vector after whole cycle
             grw = sqrt(tstpr_L2dir/tstpr_L2ini)! add direct growth
             call cnht_opcmult (VXP,VYP,VZP,TP,grw)

             ! normalise pressure ???????????????????????????
             if (tstpr_pr.ne.0) call cmult(PRP,grw,tstpr_np)
          endif

          ! run vector solver (arpack, power iteration)
          call stepper_vsolve

          if (LASTEP.ne.1) then
             ! stepper restart;
             ! set time and iteration number
             TIME=0.0
             ISTEP=0

             ! zero pressure
             if (tstpr_pr.eq.0) call rzero(PRP,tstpr_np)

             if (tstpr_mode.eq.3) then
                ! optimal initial condition
                call mntr_log(tstpr_id,lp_prd,
                $ 'Optimal initial condition; direct phase start')

                IFADJ = .FALSE.

                ! get initial L2 norm
                tstpr_L2ini = cnht_glsc2_wt(VXP,VYP,VZP,TP,
                $ VXP,VYP,VZP,TP,BM1)
                ! set cpfld for conjugated heat transfer
                if (IFHEAT) call cnht_cpfld_set

             endif
          endif

       endif ! tstpr_mode.eq.3.and.(.not.IFADJ)
       ! timing
       ltim = dnekclock() - ltim
       call mntr_tmr_add(tstpr_tmr_evl_id,1,ltim)
    endif ! mod(tstpr_istep,tstpr_step).eq.0

    return
  end subroutine tstpr_main
  !=======================================================================
  !> @brief Average velocity and temperature at element faces.
  !! @ingroup tstpr
  subroutine tstpr_dssum
    implicit none

    include 'SIZE' ! N[XYZ]1
    include 'INPUT' ! IFHEAT
    include 'SOLN' ! V[XYZ]P, TP, [VT]MULT
    include 'TSTEP' ! IFIELD
    include 'TSTPRD' ! tstpr_nt

    ! local variables
    integer ifield_tmp
    !-----------------------------------------------------------------------
    ! make sure the velocity and temperature fields are continuous at
    ! element faces and edges
    ifield_tmp = IFIELD
    IFIELD = 1
#ifdef AMR
    call amr_oph1_proj(vxp,vyp,vzp,nx1,ny1,nz1,nelv)
#else
    call opdssum(vxp,vyp,vzp)
    call opcolv (vxp,vyp,vzp,vmult)
#endif

    if(IFHEAT) then
       IFIELD = 2
#ifdef AMR
       call h1_proj(tp,nx1,ny1,nz1)
#else
       call dssum(tp,nx1,ny1,nz1)
       call col2 (tp,tmult,tstpr_nt)
#endif
    endif
    IFIELD = ifield_tmp

    return
  end subroutine tstpr_dssum
  !=======================================================================

end submodule tstpr_template
#endif

 !=========================================================================== !
