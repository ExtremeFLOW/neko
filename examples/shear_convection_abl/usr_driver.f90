program usrneko
  use neko
  use user

  type(case_t), target :: C

  ! Register all user-defined types
  call user_setup(C%user)
  call neko_init(C)
  call neko_user_access%init(C)
  call neko_solve(C)
  call neko_finalize(C)

end program usrneko
