using Higher_order_interactions
using FileIO
using Symbolics
using LinearAlgebra

data = load("tests/data/alltestdata3.hdf5")

graph = Hypergraph_from_legacy_scalars(data["madj2"],data["madj3"],size(data["madj2"],1))

lookups = prepare_lookups(graph)

(t,p_1,p_2,pi_) = calculate_hypergraph_params(graph, lookups...)

N = graph.vertex_nums
L = graph.maxinum_edge_size

eta_all_order = Vector{Vector{Float64}}(undef, L + 1)
eta_all_order[1] = zeros(Float64, N)

for l in 2:(L + 1)
    eta_all_order[l] = solve_eta(p_1, N, l, eta_all_order[l-1], lookups...)
end

println("eta (rank 2) is correct: $(count(isapprox.(eta_all_order[2], data["retime2"];rtol=1e-6)) == size(eta_all_order[2],1))")
println("eta (rank 3) is correct: $(count(isapprox.(eta_all_order[3], data["retime3"];rtol=1e-6)) == size(eta_all_order[3],1))")
println("eta (rank 4) is correct: $(count(isapprox.(eta_all_order[4], data["retime4"];rtol=1e-6)) == size(eta_all_order[4],1))")

eta_all_order[2] = data["retime2"][:]
eta_all_order[3] = data["retime3"][:]
eta_all_order[4] = data["retime4"][:]

packed_diag_index = [generate_packed_diagnomial_index(N, l, lookups[1]) for l in 1:L]

""" 
packed_index

` packed_index[idx] = packed_index[l+1]`

"""
packed_index = [ generate_packed_index(N, l, lookups[1]) for l in 0:L]

p_0 = Matrix{Float64}(I, N, N)

@variables b c delta[2:L]

game_mul = Vector{Multi_linear_polynomial}(undef, L)
for s in 2:L
    convex_PGG = payoff_of_public_good_games_with_convex_benefit(b,c,delta[s], UInt(s))
    game_mul[s] = multi_linear_polynomial(convex_PGG)
end

b_c_expr = Num(0)

for s in 2:L
    global b_c_expr += game_mul[s].related_to_self_strategy[1] * calculate_type_1(packed_index[0+1], p_2, p_0, pi_, t[s], eta_all_order)
    global b_c_expr += game_mul[s].related_to_self_strategy[2] * calculate_type_2(packed_index[0+1], p_2, p_0, pi_, t[s], eta_all_order)
end

for l in 1:(L-2)
    for s in (l+2):L
        global b_c_expr += game_mul[s].related_to_self_strategy[l+2] *
            calculate_type_5(packed_index[l+1], p_2, p_0, pi_, t[s], eta_all_order, l)
            
        global b_c_expr += game_mul[s].unrelated_to_self_strategy[l+2] *
            calculate_type_6(packed_index[l+1], p_2, p_0, pi_, t[s], eta_all_order, l)
    end
end

for l in 1:(L-1)
    for s in (l+1):L
        global b_c_expr += game_mul[s].related_to_self_strategy[l+1] *
            calculate_type_3(packed_index[l+1], p_2, p_0, pi_, t[s], eta_all_order, l)

        global b_c_expr += game_mul[s].unrelated_to_self_strategy[l+1] *
            calculate_type_4(packed_index[l+1], p_2, p_0, pi_, t[s], eta_all_order, l)

        sum_type7 = calculate_type_7(packed_diag_index[l], p_2, p_0, pi_, t[s], eta_all_order, l)
        global b_c_expr += game_mul[s].related_to_self_strategy[l+1] * sum_type7
        global b_c_expr += game_mul[s].unrelated_to_self_strategy[l+1] * sum_type7
    end
end

simplified_b_c_expr = Symbolics.simplify(b_c_expr)
coeff_b = Symbolics.coeff(simplified_b_c_expr, b)
coeff_c = Symbolics.coeff(simplified_b_c_expr, c)
b_c_ratio_expr = -coeff_c / coeff_b
values = Dict(delta[2] => data["disc1"][], delta[3] => data["disc2"][])
b_c_ratio_num = substitute(simplify(b_c_ratio_expr), values)
coeff_b_num = substitute(simplify(coeff_b), values)
coeff_c_num = substitute(simplify(coeff_c), values)
# b_c_ratio_err = (b_c_ratio_num - data["bcratio"][]) / data["bcratio"][]
# coeff_b_err = (coeff_b_num - data["fb"][]) / data["fb"][]
# coeff_c_err = (-coeff_c_num - data["fc"][]) / data["fc"][]
println("b_c_ratio: $b_c_ratio_num")
