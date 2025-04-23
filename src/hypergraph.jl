
# 

using Finch
using Combinatorics

struct Hypergraph
    edges::Tuple{Vararg{Tensor}}
end

function Hypergraph_from_legacy_scalars(madj2,madj3)
    local N = size(madj3,1)
    local indices = collect(combinations(1:N, 2))
    local index_dict::Dict = Dict((val,idx) for (idx,val) in enumerate(indices))
    function get_madj3_column_index(j,k, index_dict)
        return get(index_dict, [j,k], nothing)
    end

    local omega_2d::Tensor = Tensor(Dense(SparseList(Element(UInt8(0)))), N, N)
    local omega_3d::Tensor = Tensor(Dense(SparseList(SparseList(Element(UInt8(0))))), N, N, N)
    @finch begin
        omega_2d .= 0
        for j in _, i in _
            if (i < j)
                if madj2[i, j] != 0
                    omega_2d[i, j] = madj2[i, j]
                end
            end
        end

        omega_3d .= 0
        for k in 1:N, j in 1:N, i in 1:N
            if (i < j < k)
                let elem = madj3[i, ~get_madj3_column_index(j,k,index_dict)]
                    if elem != 0
                        omega_3d[i, j, k] = elem
                    end
                end
            end
        end
    end

    return Hypergraph(
        (omega_2d, omega_3d)
    )
end

# graph = Hypergraph_from_legacy_scalars(madj2,madj3)
