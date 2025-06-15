using Higher_order_interactions
using FileIO
using Symbolics

possible_graphs = Iterators.flatten((
    readdir("tests/data/graphs/erdos-renyi", join=true),
    readdir("tests/data/graphs/chung-lu", join=true)))
graph_categories = Dict{NamedTuple{(:node_degs, :edge_degs), Tuple{Vector{Pair{Int, Int}}, Vector{Pair{Int, Int}}}}, Vector{String}}()  # 分类字典

function extract_degree_character_from_graph(dense_mat)

    stat_node_degs = Dict{Int, Int}()  # Dictionary to store counts of node_deg
    stat_edge_degs = Dict{Int, Int}()
    for row in eachrow(dense_mat)
        node_deg = count(row .== 1)
        stat_node_degs[node_deg] = get(stat_node_degs, node_deg, 0) + 1
    end
    for col in eachcol(dense_mat)
        edge_deg = count(col .== 1)
        stat_edge_degs[edge_deg] = get(stat_edge_degs, edge_deg, 0) + 1
    end
    return (stat_node_degs,stat_edge_degs)
end

for possible_graph_path in possible_graphs
    parts = match(r"(?:\/|\\)([0-9a-zA-Z-]+)\.hdf5$", possible_graph_path)
    if !isnothing(parts)
        uuid_str = parts.captures[1]
        local data = load(possible_graph_path)
        graph_matrix = data["graph"]
        degree_character = extract_degree_character_from_graph(graph_matrix)
        immutable_character = NamedTuple{(:node_degs, :edge_degs)}((sort(collect(degree_character[1])), sort(collect(degree_character[2]))))
        if haskey(graph_categories, immutable_character)
            println("same character found")
            push!(graph_categories[immutable_character], uuid_str)
        else
            graph_categories[immutable_character] = [uuid_str]
        end
    end
end

# # 输出分类结果
# for (category, graph_list) in graph_categories
#     println("Category: ", category)
#     println("Graphs: ", graph_list)
# end
