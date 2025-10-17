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
    device_init,
    device_finalize,
    field_registry_init,
    field_registry_free,
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
    field_dofmap,
    field_space,
    case_fluid_dofmap,
    case_fluid_space,
    case_fluid_coef,
    initial_condition,
    preprocess,
    compute,
    dirichlet_condition,
    material_properties,
    source_term,
    callback_field,
    callback_field_name,
    bc_mask,
    output
)
