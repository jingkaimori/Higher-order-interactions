using Higher_order_interactions
using FileIO
using Symbolics

possible_graphs = readdir("tests/data/graphs/erdos-renyi-different-size", join=true)

for graph_path in possible_graphs
    data = load(graph_path)
    prepared_info = prepare_all(data["graph"])
    println("nodes(N): $(prepared_info[1].vertex_nums)")
    @time calculate_b_c_ratio_from_graph_matrix(prepared_info, generate_PGG_mul)
end