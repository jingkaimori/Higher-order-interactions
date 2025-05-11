
using Higher_order_interactions
using Finch
using HDF5
using Printf

p_1 = bspread("tests/data/p_1.bsp.hdf5")

code = generate_lyapunouv_solver_definitions(4, size(p_1)[1])
open("examples/solve-lyapunouv-raw.jl", "w") do f
    write(f, string(code))
end
