using PyCall
using SparseArrays
using UUIDs
using FileIO

function generate_k1_k2_type1(nodes::Int, ratio::Float64, edge_ranks::Int, edge_size::Int)
    k2 = Dict(begin
        random_value = rand()
        minimum_possibility = (1 - ratio) / (1 - ratio^edge_ranks)
        rank = 2
        possibility = minimum_possibility * ratio^(edge_ranks - 1)
        while true
            if random_value <= possibility
                break
            end
            rank += 1
            possibility += minimum_possibility * ratio^(edge_ranks - 1 - (rank - 2))
        end
        l => rank
    end for l in 1:edge_size)
    k2_sum = sum(values(k2))
    k2_equlized = k2_sum ÷ nodes
    k2_remain = k2_sum % nodes
    k1 = Dict(l => if l > k2_remain
        k2_equlized
    else
        k2_equlized + 1
    end for l in 1:nodes)
    (k1, k2)
end

function generate_k1_k2_type2(nodes::Int, edge_size::Int)
    k2 = Dict(l => rand(2:(nodes-1)) for l in 1:edge_size)
    k2_sum = sum(values(k2))
    k2_equlized = k2_sum ÷ nodes
    k2_remain = k2_sum % nodes
    k1 = Dict(l => if l > k2_remain
        k2_equlized
    else
        k2_equlized + 1
    end for l in 1:nodes)
    (k1, k2)
end

function generate_valid_chung_lu_graph(k1, k2)

    generator = pyimport("hypernetx.algorithms.generative_models")
    py = pyimport("builtins")

    while true
        local graph_hnx = generator.chung_lu_hypergraph(k1, k2)
        local len = py.len(py.list(graph_hnx.connected_components()))
        if len > 1
            println("Graph has more than one connected component, regenerating...")
            continue
        end

        local ic_sm = graph_hnx.incidence_matrix()

        local ic_sm_julia = py_csr_to_julia_sparse(ic_sm.tocsc())

        local ic_julia = Array(ic_sm_julia)

        local single_one_columns = findall(col -> sum(col .== 1) == 1, eachcol(ic_julia))

        if !isempty(single_one_columns)
            println("Columns with only one '1' and the rest '0', regenerating...")
        else
            println("No such columns found.")
            return ic_julia
        end
    end
end

function generate_valid_e_r_graph(n,m,p)

    generator = pyimport("hypernetx.algorithms.generative_models")
    py = pyimport("builtins")

    while true
        local graph_hnx = generator.erdos_renyi_hypergraph(n,m,p)
        local len = py.len(py.list(graph_hnx.connected_components()))
        if len > 1
            println("Graph has more than one connected component, regenerating...")
            continue
        end

        local ic_sm = graph_hnx.incidence_matrix()

        local ic_sm_julia = py_csr_to_julia_sparse(ic_sm.tocsc())

        local ic_julia = Array(ic_sm_julia)

        local single_one_columns = findall(col -> sum(col .== 1) == 1, eachcol(ic_julia))

        if !isempty(single_one_columns)
            println("Columns with only one '1' and the rest '0', regenerating...")
        else
            println("No such columns found.")
            return ic_julia
        end
    end
end

function save_graph_data(incidence_matrix::Array, uuid::UUID, subclass::String = "erdos-renyi")
    # 将图数据保存到文件
    println("save graph $uuid")
    save("tests/data/graphs/$subclass/$uuid.hdf5", Dict("graph" => incidence_matrix))
end


function py_csr_to_julia_sparse(py_csc)
    # 提取稀疏矩阵的组成部分
    data = py_csc.data
    indices = py_csc.indices .+ 1  # Julia 使用 1-based 索引
    indptr = py_csc.indptr .+ 1    # Julia 使用 1-based 索引
    m, n = py_csc.shape

    # 构造 Julia 的 SparseMatrixCSC
    return SparseMatrixCSC(m, n, indptr, indices, data)
end

function extract_real_k_from_graph(dense_mat)
    k1 = Dict{Int, Int}()
    k2 = Dict{Int, Int}()
    for row in eachrow(dense_mat)
        node_deg = count(row .== 1)
        setindex!(k1, node_deg, row.indices[1])
    end
    for col in eachcol(dense_mat)
        edge_deg = count(col .== 1)
        setindex!(k2, edge_deg, col.indices[2])
    end

    return (k1, k2)
end

if false
    

for highest_level in 3:7
    assumed_k = generate_k1_k2_type2(10, 60)
    local test_graph = generate_valid_chung_lu_graph(assumed_k...)
    real_k = extract_real_k_from_graph(test_graph)
    for _ in 1:8
        id = uuid1()
        local real_graph = generate_valid_chung_lu_graph(real_k...)
        real_k = extract_real_k_from_graph(real_graph)
        save_graph_data(real_graph, id, real_k...)
    end
    
end
    
end

for N in 15:16
    graph_mat= generate_valid_e_r_graph(N, N*2, 0.5)
    id = uuid1()
    save_graph_data(graph_mat, id, "erdos-renyi-different-size")
end
