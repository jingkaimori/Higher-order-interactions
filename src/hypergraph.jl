
# 

using Finch
using Combinatorics

struct Hypergraph
    edges::Array{Finch.AbstractTensor}
end

function Hypergraph_from_legacy_scalars(madj2,madj3)
    local N = size(madj3,1)
    local indices = collect(combinations(1:N, 2))
    local index_dict::Dict = Dict((val,idx) for (idx,val) in enumerate(indices))
    function get_madj3_column_index(j,k, index_dict)::Union{Int,Nothing}
        local res = get(index_dict, [j,k], nothing)
        if isnothing(res)
            return get(index_dict, [k,j], nothing)
        else
            return res
        end
    end

    local omega_2d::Tensor = Tensor(Dense(SparseList(Element(UInt8(0)))), N, N)
    local omega_3d::Tensor = Tensor(Dense(SparseList(SparseList(Element(UInt8(0))))), N, N, N)
    @finch begin
        omega_2d .= 0
        for j in _, i in _
            if madj2[i, j] != 0
                omega_2d[i, j] = madj2[i, j]
            end
        end

        omega_3d .= 0
        for k in 1:N, j in 1:N, i in _
            if j != k
                let elem = madj3[i, ~get_madj3_column_index(j,k,index_dict)]
                    if elem != 0
                        omega_3d[i, j, k] = elem
                    end
                end
            end
        end
    end

    return Hypergraph(
        [Scalar(0), omega_2d, omega_3d]
    )
end

# graph = Hypergraph_from_legacy_scalars(madj2,madj3)
