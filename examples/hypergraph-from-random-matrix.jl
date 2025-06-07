using UUIDs
using HDF5
using FileIO
using Higher_order_interactions
using SimpleHypergraphs
using Random
using StatsBase


function save_graph_data(incidence_matrix::Array, uuid::UUID, subclass::String="erdos-renyi")
    println("save graph $uuid")
    save("tests/data/graphs/$subclass/$uuid.hdf5", Dict("graph" => incidence_matrix))
end

function check_gale_ryser(node_distrib::Vector{Int}, edge_distrib::Vector{Int})::Bool
    # 将行和列分布按降序排列
    r = sort(node_distrib, rev = true)
    c = sort(edge_distrib, rev = true)

    # 检查 Gale-Ryser 条件
    for k in 1:length(r)
        sum_r = sum(r[1:k])
        sum_c = sum(min(ci, k) for ci in c)
        if sum_r > sum_c
            return false
        end
    end

    return true
end


function generate_graphs_from_node_edge_distrib(node_distrib::Vector{Int}, edge_distrib::Vector{Int})
    edge_num = length(edge_distrib)
    N = length(node_distrib)

    @boundscheck if sum(node_distrib) != sum(edge_distrib)
        error("sum(node_distrib) != sum(edge_distrib)")
    end

    #Gale-Ryser 条件
    if !check_gale_ryser(node_distrib, edge_distrib)
        println("Gale-Ryser not satisfied")
        return
    end

    matrix = zeros(Bool, N, edge_num)

    # sort!(edge_distrib, rev=true)
    for _ in 1:10_000
        # println("trial $trial")
        generated_succeed = true
        matrix .= false

        shuffle!(edge_distrib)
        nodes_missing_edge = copy(node_distrib)
        edges_appeared = Set{Vector{Int}}()
        for (j, l_j) in enumerate(edge_distrib)
            nodes_vacant = findall(nodes_missing_edge .!== 0)
            # sort!(nodes_vacant)

            if length(nodes_vacant) < l_j
                # println("not enough spaces, truncated by $j, $l_j, $nodes_vacant")
                generated_succeed = false
                break
            end

            nodes_vacant_weight = Float16.(nodes_missing_edge[nodes_vacant])
            # (wmin,wmax) = extrema(nodes_vacant_weight)
            # nodes_vacant_weight .= 1 ./ nodes_vacant_weight
            # nodes_vacant_weight .= (nodes_vacant_weight .- (0.9 * wmin)) .* (.- nodes_vacant_weight .+ (1.1 * wmax))
            # nodes_vacant_weight ./= sum(nodes_vacant_weight)
            # println("$l_j, $nodes_vacant_weight, $wmin")
            # nodes_vacant_weight = fill(8.0, size(nodes_vacant))
            edge_inserted = false
            while true
                edge = sample(nodes_vacant, Weights(nodes_vacant_weight), l_j; replace=false)
                sort!(edge)
                if edge ∉ edges_appeared
                    matrix[edge, j] .= 1
                    nodes_missing_edge[edge] .-= 1
                    push!(edges_appeared, edge)
                    edge_inserted = true
                    # println("edge inserted: $edge")
                    break
                else
                    modified_nodes_vacant_weight = view(nodes_vacant_weight, findall(in(edge), nodes_vacant))
                    modified_nodes_vacant_weight .*= Float16(0.9)
                    # println("edge failed: $edge, $nodes_vacant, $nodes_vacant_weight, $(modified_nodes_vacant_weight .* Float16(0.9)), $(nodes_vacant_weight[findall(in(edge), nodes_vacant)])")
                    if isapprox(minimum(nodes_vacant_weight), 0; atol=eps(Float16))
                        break
                    end
                end
            end
            if !edge_inserted
                # println("failed to get unique edge")
                generated_succeed = false
                break
            end
        end

        if generated_succeed
            println("generate success, validating")
        else
            println("failed to generate graph, regenerating")
            continue
        end
        if has_duplicate_columns(matrix)
            println("duplicated column found, $(find_duplicate_columns(matrix)), regenerating")
            continue
        end
        (node_distrib_gen, edge_distrib_gen) = extract_degree_character_from_graph(BitMatrix(matrix))
        if node_distrib != node_distrib_gen
            println("node distribution mismatch, regnenrating")
            continue
        end
        if edge_distrib != edge_distrib_gen
            println("edge distribution mismatch, regnenrating")
            continue
        end
        graph_sphp = SimpleHypergraphs.Hypergraph([v ? v : nothing for v in matrix])
        len = length(get_connected_components(graph_sphp))
        if len > 1
            println("Graph has more than one connected component, regenerating...")
            continue
        end
        return matrix
    end
end

## 

node_distrib = [6,  3, 6,  3,  3, 6,  5,  3,  3,  2]
edge_distrib = fill(4,10)
for _ in 1:100
    uuid = uuid1()
    mat = generate_graphs_from_node_edge_distrib(node_distrib, edge_distrib)
    save_graph_data(mat, uuid, "node-edge-distrib1")
end

##

node_distrib = [6,  3, 6,  3,  3, 6,  5,  3,  3,  2]
edge_distrib = vcat(fill(2,6), fill(3, 3), fill(4,3), [7])
for _ in 1:100
    uuid = uuid1()
    mat = generate_graphs_from_node_edge_distrib(node_distrib, edge_distrib)
    save_graph_data(mat, uuid, "node-edge-distrib2")
end

## 

node_distrib = [6,  3, 6,  3,  3, 6,  5,  3,  3,  2]
edge_distrib = fill(2,20)
for _ in 1:6
    for _ in 1:5
        uuid = uuid1()
        mat = generate_graphs_from_node_edge_distrib(node_distrib, edge_distrib)
        save_graph_data(mat, uuid, "node-edge-variance")
    end

    node_distrib .+= 2
    edge_distrib .+= 1
end

##

node_distrib = [13,  2, 14,  2,  5, 14,  10,  10,  8,  2]
edge_distrib = [2, 2, 2, 2, 2, 2, 2, 2, 2, 4,4,6, 6, 6, 6, 6, 8, 8, 8]
for _ in 1:6
    for _ in 1:5
        uuid = uuid1()
        mat = generate_graphs_from_node_edge_distrib(node_distrib, edge_distrib)
        save_graph_data(mat, uuid, "node-distrib")
    end

    while true
        (idx1, idx2) = sample(eachindex(edge_distrib), 2; replace=false)
        if (edge_distrib[idx1] - edge_distrib[idx2]) in -1:1 
            continue
        end
        edge_degree_sum = edge_distrib[idx1] + edge_distrib[idx2]
        edge_distrib[idx1] = fld(edge_degree_sum, 2)
        edge_distrib[idx2] = fld(edge_degree_sum, 2) + mod(edge_degree_sum, 2)
        break
    end
end

##

node_distrib = [13,  2, 14,  2,  5, 14,  10,  10,  8,  2]
edge_distrib = [2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 4, 4, 4, 6, 6, 6, 6, 6, 8, 8]
for _ in 1:20
    for _ in 1:3
        uuid = uuid1()
        mat = generate_graphs_from_node_edge_distrib(node_distrib, edge_distrib)
        save_graph_data(mat, uuid, "edge-distrib")
    end

    while true
        (idx1, idx2) = sample(eachindex(node_distrib), 2; replace=false)
        node_distrib[idx1] - node_distrib[idx2]
        if (node_distrib[idx1] - node_distrib[idx2]) in -1:1 
            continue
        end
        node_degree_sum = node_distrib[idx1] + node_distrib[idx2]
        node_distrib[idx1] = fld(node_degree_sum, 2)
        node_distrib[idx2] = fld(node_degree_sum, 2) + mod(node_degree_sum, 2)
        break
    end
end
