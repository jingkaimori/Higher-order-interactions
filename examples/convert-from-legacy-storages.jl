using FileIO
using Symbolics
using UUIDs
using Statistics
using Higher_order_interactions

const GraphEntry = Matrix{Int8}

possible_graphs_by_category = Dict( 
    category => readdir("tests/data/graphs/$category";join=true) for category in [
        "anzhisheng-richclub",
        "anzhisheng-star",
        "well-mixed",
        "edge-distrib",
        "node-distrib",
        "node-edge-variance",
        "node-edge-distrib1",
        "node-edge-distrib2",
    ]
)
possible_graphs = Iterators.flatten(values(possible_graphs_by_category))
stored_graph_data = Dict{String, GraphEntry}()
sizehint!(stored_graph_data, 600)

for possible_graph_path in possible_graphs
    parts = match(r"(?:\/|\\)([0-9a-z-]+)(?:\/|\\)([0-9a-zA-Z-]+)\.hdf5$", possible_graph_path)
    if !isnothing(parts)
        graph_type = parts.captures[1]
        uuid_str = parts.captures[2]
        local data = load(possible_graph_path)
        stored_graph_data["incidence/$graph_type/$uuid_str"] = Int8.(data["graph"])
    end
end

stored_group_parameters = Dict{String, Vector{Int}}()
edge_graph_path = first(possible_graphs_by_category["edge-distrib"])
(_, edge_degree_controlled) = extract_degree_character_from_graph(
    (load(edge_graph_path))["graph"] .== 1
)
stored_group_parameters["group parameters/edge-distrib/edge degree"] = edge_degree_controlled

node_graph_path = first(possible_graphs_by_category["node-distrib"])
(node_degree_controlled, _) = extract_degree_character_from_graph(
    (load(node_graph_path))["graph"] .== 1
)
stored_group_parameters["group parameters/node-distrib/node degree"] = node_degree_controlled

node_edge_graph_path = first(possible_graphs_by_category["node-edge-distrib1"])
(node_degree_controlled, edge_degree_controlled) = extract_degree_character_from_graph(
    (load(node_edge_graph_path))["graph"] .== 1
)
stored_group_parameters["group parameters/node-edge-distrib1/node degree"] = node_degree_controlled
stored_group_parameters["group parameters/node-edge-distrib1/edge degree"] = edge_degree_controlled

node_edge_graph_path = first(possible_graphs_by_category["node-edge-distrib2"])
(node_degree_controlled, edge_degree_controlled) = extract_degree_character_from_graph(
    (load(node_edge_graph_path))["graph"] .== 1
)
stored_group_parameters["group parameters/node-edge-distrib2/node degree"] = node_degree_controlled
stored_group_parameters["group parameters/node-edge-distrib2/edge degree"] = edge_degree_controlled

save("tests/data/hypergraphs.hdf5", merge(stored_graph_data, stored_group_parameters))