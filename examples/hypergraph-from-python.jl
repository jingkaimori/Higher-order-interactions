using PyCall
using SparseArrays
using UUIDs

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

function save_graph_data(incidence_matrix::Array, uuid::UUID)
    # 将图数据保存到文件
    write("tests/data/graphs/$uuid.txt", matrix_to_string(incidence_matrix))
    save("tests/data/graphs/$uuid.hdf5", Dict("graph" => incidence_matrix))
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

function matrix_to_string(dense_mat)

    # 将每一行转换为字符串，并用空格分隔
    rows = size(dense_mat, 1)
    result = ""
    for i in 1:rows
        row = dense_mat[i, :]
        row_str = join(row, " ")
        result *= row_str * "\n"
    end

    return result
end

for highest_level in 3:7
    k = generate_k1_k2_type1(10, 5.0, highest_level,20)
    id = uuid1()
    local graph = generate_valid_chung_lu_graph(k...)
    save_graph_data(graph, id)
end
