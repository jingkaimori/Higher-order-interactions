using Higher_order_interactions
using FileIO
using Symbolics

data = load("tests/data/graphs/erdos-renyi-different-size/92487300-3094-11f0-374e-87dd78b2e2aa.hdf5")

prepared_info = prepare_all(data["graph"], generate_PGG_mul, nothing)

println("N: $(prepared_info[1].vertex_nums)")
(b_c_ratio_expr,_,_,delta) = calculate_b_c_ratio_from_graph_matrix(prepared_info)
println("PGG:")
# delta_values = Dict(delta[2] => 1.0, delta[3] => 1.6457513110645907)
delta_values = Dict(delta[l] => 1.1^l for l in 2:prepared_info[1].maxinum_edge_size)
b_c_ratio_num = substitute((b_c_ratio_expr), delta_values)
println("deltas: $delta_values")
println("b_c_ratio: $b_c_ratio_num ")

prepared_info = prepare_all(data["graph"], generate_TPGG_mul, prepared_info)
(b_c_ratio_expr,) = calculate_b_c_ratio_from_graph_matrix(prepared_info)
println("TPGG:")
b_c_ratio_num = substitute((b_c_ratio_expr), Dict())
println("threshods: all 1")
println("b_c_ratio: $b_c_ratio_num ")

prepared_info = prepare_all(data["graph"], generate_MSG_mul, prepared_info)
(b_c_ratio_expr,) = calculate_b_c_ratio_from_graph_matrix(prepared_info)
println("MSG")
b_c_ratio_num = substitute((b_c_ratio_expr), Dict())
println("b_c_ratio: $b_c_ratio_num ")