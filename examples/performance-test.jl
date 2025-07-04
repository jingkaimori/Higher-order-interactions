using Higher_order_interactions
using FileIO
using Symbolics

possible_graphs = readdir("tests/data/graphs/erdos-renyi-different-size", join=true)

graph_mat_by_size = Dict{Int, Matrix{Int}}()
for graph_path in possible_graphs
    data = load(graph_path)
    graph_mat_by_size[size(data["graph"], 1)] = data["graph"]
end

begin
    local prepared_info = prepare_all(graph_mat_by_size[10], generate_MSG_mul)
    calculate_b_c_ratio_from_graph_matrix(prepared_info)
end

# 按照键的从小到大读取 graph_mat_by_size 中的值
for key in sort(collect(keys(graph_mat_by_size)))
    if key > 20
        break
    end
    local prepared_info = prepare_all(graph_mat_by_size[key], generate_MSG_mul)
    println("nodes(N): $key")
    @time calculate_b_c_ratio_from_graph_matrix(prepared_info)
end
