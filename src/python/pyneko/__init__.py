r"""

                 _  __ ____ __ __ ____
     ___  __ __ / |/ // __// //_// __ \
    / _ \/ // //    // _/ / ,<  / /_/ /
   / .__/\_, //_/|_//___//_/|_| \____/
  /_/   /___/

  A Python interface to Neko

"""

from pyneko.intf import (
    init,
    finalize,
    job_info,
    case_init,
    case_free,
    time,
    end_time,
    tstep,
    solve,
    step,
    field,
    field_order,
    field_nelements,
    field_size,
    output    
)
