using Neko
using JSON

# Initialise Neko
Neko.init()
Neko.job_info()

# Create a Neko case from a JSON file
cylinder_json = JSON.parsefile("cylinder.case")
cylinder_case = Neko.case_init(cylinder_json)

# To solve the entire case we can call solve()
#Neko.solve(cylinder_case)

# To manually step forward in time, call step()
while Neko.time(cylinder_case) < Neko.end_time(cylinder_case)
    Neko.step(cylinder_case)
end

# To retrive a Julia struct, with wraped array of the data in a field
# from Neko's field registry use field()
u = Neko.field("u")

# Cleanup
Neko.case_free(cylinder_case)
Neko.finalize()
    



