program usrneko
  use neko
  use user
  type(case_t), target :: C
  
  call user_setup(C%usr)
  call neko_init(C)
  call neko_solve(C)
  call neko_finalize(C)


end program usrneko
