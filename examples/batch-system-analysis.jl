using Higher_order_interactions
using FileIO
using Symbolics
using UUIDs
using Statistics

possible_graphs = Iterators.flatten((
    readdir("tests/data/graphs/erdos-renyi", join=true),
    readdir("tests/data/graphs/chung-lu", join=true),
    readdir("tests/data/graphs/well-mixed", join=true),
    readdir("tests/data/graphs/anzhisheng-star", join=true),
    readdir("tests/data/graphs/anzhisheng-richclub", join=true)))

function extract_degree_character_from_graph(dense_mat::BitMatrix)

    (node_num, edge_num) = size(dense_mat)

    stat_node_degs = zeros(Int, node_num)  # Dictionary to store counts of node_deg
    stat_edge_degs = zeros(Int, edge_num)
    for row in eachrow(dense_mat)
        node_deg = count(row)
        stat_node_degs[row.indices[1]] = node_deg
    end
    for col in eachcol(dense_mat)
        edge_deg = count(col)
        stat_edge_degs[col.indices[2]] = edge_deg
    end
    return (stat_node_degs,stat_edge_degs)
end

function append_to_statstics()
    push!(uuids, id)
    push!(gametype, "PGG")
    push!(graphtype, graph_type)
    push!(node_deg_avg, mean(node_degree_distrib))
    push!(edge_deg_avg, mean(edge_degree_distrib))
    push!(clustering_coefficent, mean(clustering_coefficent_distrib))
    push!(node_deg_val, var(node_degree_distrib))
    push!(edge_deg_val, var(edge_degree_distrib))
    push!(clustering_coefficent_val, var(clustering_coefficent_distrib))
    push!(b_c_ratios, b_c_ratio_num)
end

graph_mat_by_id = Dict{String, Tuple{String, Matrix{Int}}}()

for possible_graph_path in possible_graphs
    parts = match(r"(?:\/|\\)([a-z-]+)(?:\/|\\)([0-9a-zA-Z-]+)\.hdf5$", possible_graph_path)
    if !isnothing(parts)
        graph_type = parts.captures[1]
        uuid_str = parts.captures[2]
        local data = load(possible_graph_path)
        graph_mat_by_id[uuid_str] = (graph_type, data["graph"])
    end
end

prepared_info = nothing
data_num = length(graph_mat_by_id)
result_num = 3*data_num

uuids = Vector{String}()
sizehint!(uuids, result_num)
gametype = Vector{String}()
sizehint!(uuids, gametype)
graphtype = Vector{String}()
sizehint!(uuids, graphtype)
node_deg_avg = Vector{Float32}()
sizehint!(node_deg_avg, result_num)
node_deg_val = Vector{Float32}()
sizehint!(node_deg_val, result_num)
edge_deg_avg = Vector{Float32}()
sizehint!(edge_deg_avg, result_num)
edge_deg_val = Vector{Float32}()
sizehint!(edge_deg_val, result_num)
clustering_coefficent = Vector{Float32}()
sizehint!(clustering_coefficent, result_num)
clustering_coefficent_val = Vector{Float32}()
sizehint!(clustering_coefficent_val, result_num)
b_c_ratios = Vector{Float32}()
sizehint!(b_c_ratios, result_num)

for (id, graph_matrix_and_type) in graph_mat_by_id
    (graph_type, graph_matrix) = graph_matrix_and_type
    boolean_graph_matrix = graph_matrix .== 1
    (node_degree_distrib, edge_degree_distrib) = extract_degree_character_from_graph(boolean_graph_matrix)
    clustering_coefficent_distrib = quad_clustering(boolean_graph_matrix)

    global prepared_info = prepare_all(graph_matrix, generate_PGG_mul, prepared_info)

    println("id: $id")
    (b_c_ratio_expr,_,_,delta) = calculate_b_c_ratio_from_graph_matrix(prepared_info)
    delta_values = Dict(delta[l] => 1.1^l for l in 2:prepared_info[1].maxinum_edge_size)
    b_c_ratio_num = substitute((b_c_ratio_expr), delta_values)
    push!(uuids, id)
    push!(gametype, "PGG")
    push!(graphtype, graph_type)
    push!(node_deg_avg, mean(node_degree_distrib))
    push!(edge_deg_avg, mean(edge_degree_distrib))
    push!(clustering_coefficent, mean(clustering_coefficent_distrib))
    push!(node_deg_val, var(node_degree_distrib))
    push!(edge_deg_val, var(edge_degree_distrib))
    push!(clustering_coefficent_val, var(clustering_coefficent_distrib))
    push!(b_c_ratios, b_c_ratio_num)
end

for (id, graph_matrix_and_type) in graph_mat_by_id
    (graph_type, graph_matrix) = graph_matrix_and_type
    boolean_graph_matrix = graph_matrix .== 1
    (node_degree_distrib, edge_degree_distrib) = extract_degree_character_from_graph(boolean_graph_matrix)
    clustering_coefficent_distrib = quad_clustering(boolean_graph_matrix)

    global prepared_info = prepare_all(graph_matrix, generate_TPGG_mul, prepared_info)

    (b_c_ratio_expr,) = calculate_b_c_ratio_from_graph_matrix(prepared_info)
    b_c_ratio_num = substitute((b_c_ratio_expr), Dict())
    push!(uuids, id)
    push!(gametype, "TPGG")
    push!(graphtype, graph_type)
    push!(node_deg_avg, mean(node_degree_distrib))
    push!(edge_deg_avg, mean(edge_degree_distrib))
    push!(clustering_coefficent, mean(clustering_coefficent_distrib))
    push!(node_deg_val, var(node_degree_distrib))
    push!(edge_deg_val, var(edge_degree_distrib))
    push!(clustering_coefficent_val, var(clustering_coefficent_distrib))
    push!(b_c_ratios, b_c_ratio_num)
    
end

for (id, graph_matrix_and_type) in graph_mat_by_id
    (graph_type, graph_matrix) = graph_matrix_and_type
    boolean_graph_matrix = graph_matrix .== 1
    (node_degree_distrib, edge_degree_distrib) = extract_degree_character_from_graph(boolean_graph_matrix)
    clustering_coefficent_distrib = quad_clustering(boolean_graph_matrix)

    global prepared_info = prepare_all(graph_matrix, generate_MSG_mul, prepared_info)

    (b_c_ratio_expr,) = calculate_b_c_ratio_from_graph_matrix(prepared_info)
    b_c_ratio_num = substitute((b_c_ratio_expr), Dict())
    push!(uuids, id)
    push!(gametype, "MSG")
    push!(graphtype, graph_type)
    push!(node_deg_avg, mean(node_degree_distrib))
    push!(edge_deg_avg, mean(edge_degree_distrib))
    push!(clustering_coefficent, mean(clustering_coefficent_distrib))
    push!(node_deg_val, var(node_degree_distrib))
    push!(edge_deg_val, var(edge_degree_distrib))
    push!(clustering_coefficent_val, var(clustering_coefficent_distrib))
    push!(b_c_ratios, b_c_ratio_num)
end

save("tests/data/statstics.hdf5", Dict(
    "uuids" => uuids,
    "gametype" => gametype,
    "graphtype" => graphtype,
    "node_deg_avg" => node_deg_avg,
    "edge_deg_avg" => edge_deg_avg,
    "clustering_coefficent" => clustering_coefficent,
    "node_deg_val" => node_deg_val,
    "edge_deg_val" => edge_deg_val,
    "clustering_coefficent_val" => clustering_coefficent_val,
    "b_c_ratios" => b_c_ratios
))
