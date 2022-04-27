!> Main driver for TurboNeko
program turboneko
  use neko
  type(case_t), target :: C
  
  call neko_init(C)
  call neko_solve(C)
  call neko_finalize(C)

end program turboneko
