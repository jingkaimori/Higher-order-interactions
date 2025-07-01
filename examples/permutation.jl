using Higher_order_interactions
using HDF5
using Symbolics
using LinearAlgebra
using AMD

incidence_matrix = h5open(
    "tests/data/hypergraphs.hdf5",
    "r";
) do fid
    incidence_dataset = fid["incidence"]["well-mixed"]["bca3ea60-3019-11f0-36a5-4f40b92d1726"]
    incidence_matrix = read(incidence_dataset)
    return Int64.(incidence_matrix)
end

(N, L, min_l) = scan_incidence_matrix(incidence_matrix)
lookups = prepare_lookups(N, L)
pt = lookups[1]
graph = Hypergraph_from_incidence_matrix(incidence_matrix, N, L, min_l, pt)

(t,r) = calculate_hypergraph_params(graph, lookups...)

(p_1,p_2,pi_) = calculate_random_walk_prob(r, N, lookups...)

eta_all_order = Vector{Vector{Float64}}(undef, L + 1)
eta_all_order[1] = zeros(Float64, N)
permutation = Vector{Vector{Int16}}(undef, L + 1)

for l in 2:(L + 1)
    eta_all_order[l] = solve_eta(p_1, N, l, eta_all_order[l-1], lookups...)
    (A, b) = Higher_order_interactions.get_eta_equation(p_1, N, l, eta_all_order[l-1], lookups...)
    permutation[l] = amd(A)
end