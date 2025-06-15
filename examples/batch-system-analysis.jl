using Higher_order_interactions
using HDF5
using Symbolics
using UUIDs
using Statistics

StatsticsEntry = @NamedTuple begin
    uuid::NTuple{2, UInt64}
    gametype::UInt8
    graphtype::UInt8
    node_deg_avg::Float32
    edge_deg_avg::Float32
    node_deg_val::Float32
    edge_deg_val::Float32
    clustering_coefficent::Float32
    clustering_coefficent_val::Float32
    b_c_ratio::Float32
end
statstics = Vector{StatsticsEntry}()
sizehint!(statstics, 1800)

game_list = (
    (generate_PGG_mul, "NPGG"),
    (generate_TPGG_mul, "TPGG"),
    (generate_MSG_mul, "MSG"),
)
category_metadata = Dict{AbstractString, UInt8}()

h5open(
    "tests/data/hypergraphs.hdf5",
    "r";
) do fid
    category_groups = fid["incidence"]
    for (category_id, category_group) in enumerate(category_groups)
        raw_str = HDF5.name(category_group)
        parts = match(r"^\/incidence\/([0-9a-z-]+)$", raw_str)
        category_name = parts.captures[1]
        category_metadata[category_name] = UInt8(category_id)
    end

    for (game_type_id, (game_generator,)) in enumerate(game_list)
        prepared_info = nothing
        for category_group in category_groups
            for incidence_dataset in category_group
                raw_str = HDF5.name(incidence_dataset)
                parts = match(r"^\/incidence\/([0-9a-z-]+)\/([0-9a-zA-Z-]+)$", raw_str)
                uuid_str = parts.captures[2]
                category_name = parts.captures[1]

                incidence_matrix = read(incidence_dataset)
                boolean_incidence_matrix = incidence_matrix .== 1
                (node_degree_distrib, edge_degree_distrib) = extract_degree_character_from_graph(boolean_incidence_matrix)
                clustering_coefficent_distrib = quad_clustering(boolean_incidence_matrix)

                prepared_info = prepare_all(Int64.(incidence_matrix), game_generator, prepared_info)

                (b_c_ratio_expr,_,_,delta) = calculate_b_c_ratio_from_graph_matrix(prepared_info)
                delta_values = Dict(delta[l] => 1.01^l for l in 2:prepared_info[1].maxinum_edge_size)
                b_c_ratio_num = substitute((b_c_ratio_expr), delta_values)

                push!(statstics, (
                    uuid = UUID(uuid_str),
                    gametype = game_type_id,
                    graphtype = category_metadata[category_name],
                    node_deg_avg = Float32(mean(node_degree_distrib)),
                    edge_deg_avg = Float32(mean(edge_degree_distrib)),
                    node_deg_val = Float32(var(node_degree_distrib)),
                    edge_deg_val = Float32(var(edge_degree_distrib)),
                    clustering_coefficent = Float32(mean(clustering_coefficent_distrib)),
                    clustering_coefficent_val = Float32(var(clustering_coefficent_distrib)),
                    b_c_ratio = Float32(b_c_ratio_num)
                ))
            end
        end
    end
end

h5open("tests/data/statstics.hdf5", "w") do fid
    write(fid, "statstics", statstics)
    write(fid, "game-type-strings", [game_name for (_, game_name) in game_list])
    sorted_categories = sort(collect(category_metadata), by = x -> x[2])
    write(fid, "graph-type-strings", [x[1] for x in sorted_categories])
end
