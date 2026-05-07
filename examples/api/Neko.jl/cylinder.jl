using Neko
using JSON

Neko.init()
Neko.job_info()

cylinder_json = JSON.parsefile("cylinder.case")

# Define initial conditions
function initial_condition(scheme_name::Ptr{Cchar}, len::Cint)::Cvoid

    # In this case the check is redudnant, but if both fluid and scalar
    # initial conditions should be defined, `scheme_name` can be used
    # to identifiy which solver has called the callback
    if (unsafe_string(scheme_name) == "fluid")
        # To retrive a Julia struct, with wrapped arrays of the data in a field
        # from Neko's field registry use field()
        u = Neko.field("u")
        nelm = Neko.field_nelements("u")
        lx = u.Xh.lx

        for e in 1:nelm
            for k in 1:lx
                for j in 1:lx
                    for i in 1:lx
                        u.x[i,j,k,e] = 1.0
                    end
                end
            end
        end
    end

    return nothing
end;

# Define inflow conditions
function inflow(msk::Ptr{Cint}, msk_size::Cint,
                   t::Cdouble, tstep::Cint)::Cvoid

    bc_msk = Neko.bc_mask(msk, msk_size)

    # In this case this check is redudant, but if different user provied
    # boundary conditions share the same callback, one could use the content
    # of the callbacks field list to identify which condition should be applied,
    # for example a velocity condition passes the fields (u, v, w)
    if (Neko.callback_field_name(1, "u"))
        u = Neko.callback_field("u")
        v = Neko.callback_field("v")
        w = Neko.callback_field("w")
        lx = u.Xh.lx

        for i in 1:msk_size
            # Compute the (i,j,k,e) index given the linear index in bc_msk
            idx = Neko.nonlinear_index(bc_msk[i], lx, lx, lx)
            u.x[idx] = 1.0
            v.x[idx] = 0.0
            w.x[idx] = 0.0
        end
    end

    return nothing
end;

# Create a Neko callback for the initial and inflow conditions
const cb_initial = @Neko.initial_condition(initial_condition)
const cb_inflow = @Neko.dirichlet_condition(inflow)

# Create a Neko case from a JSON file and provied (optional) callbacks
cylinder_case = Neko.case_init(cylinder_json,
                               cb_initial_condition = cb_initial,
                               cb_dirichlet_condition = cb_inflow)

# To solve the entire case we can call solve()
#Neko.solve(cylinder_case)

# To manually step forward in time, call step()
while Neko.time(cylinder_case) < Neko.end_time(cylinder_case)
    Neko.step(cylinder_case)
end

# Cleanup
Neko.case_free(cylinder_case)
Neko.finalize()
