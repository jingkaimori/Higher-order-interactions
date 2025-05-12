using Higher_order_interactions
using FileIO
using BenchmarkTools

data = load("tests/data/alltestdata3.hdf5")

graph = Hypergraph_from_legacy_scalars(data["madj2"],data["madj3"],Int(data["n"][1]))

lookups = prepare_lookups(graph)

(t,p_1,p_2,pi_) = calculate_hypergraph_params(graph, lookups...)

N = graph.vertex_nums

eta_all_order = Vector{Vector{Float64}}(undef, graph.maxinum_edge_size + 1)
eta_all_order[1] = zeros(Float64, N)

for l in 2:(graph.maxinum_edge_size + 1)
    eta_all_order[l] = solve_eta(p_1, N, l, eta_all_order[l-1], lookups...)
end

# solve_eta(res[3],graph.vertex_nums, 2, Vector

