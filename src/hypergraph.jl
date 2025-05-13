
struct Hypergraph
    edges::Array{Array{Int,1},1}
    vertex_nums::Int
    mininum_edge_size::Int
    maxinum_edge_size::Int
end

Base.size(h::Hypergraph) = (h.maxinum_edge_size,)
Base.getindex(h::Hypergraph, i::Int) =
    if i >= h.mininum_edge_size
        return h.edges[i-h.mininum_edge_size+1]
    else
        throw(BoundsError(h, i))
    end
Base.setindex!(h::Hypergraph, v::Array{Int,1}, i::Int) =
    if i >= h.mininum_edge_size
        (h.edges[i-h.mininum_edge_size+1] = v)
    else
        throw(BoundsError(h, i))
    end
Base.IndexStyle(::Type{<:Hypergraph}) = IndexLinear()

function Hypergraph_from_legacy_scalars(madj2, madj3, N)
    local indices_2 = collect(Combinations(N, 2))
    local indices_3 = collect(Combinations(N, 3))
    edge2 = Bool[]
    sizehint!(edge2, length(indices_2))
    edge3 = Bool[]
    sizehint!(edge3, length(indices_3))
    for index in indices_2
        (index1, index2) = index
        push!(edge2, madj2[index1, index2] != 0)
    end
    for index in indices_3
        (index1, index2, index3) = index
        push!(edge3, madj3[index1, get_combination_code([index2, index3], N)] != 0)
    end

    return Hypergraph(
        [edge2, edge3], N, 2, 3
    )
end

function Hypergraph_from_incidence_matrix(matrix::Matrix{Int64}, N::Int, pt::PascalsTriangle)
    
    edges = Vector{Vector{Int}}(undef, N)
    for l in 1:N
        edges[l] = zeros(Int, pt[N, l])
    end
    col_maxsize = 0
    col_minsize = N
    for col in eachcol(matrix)
        col_nz = findall(col .== 1)
        col_size = length(col_nz)
        col_maxsize = max(col_maxsize, col_size)
        col_minsize = min(col_minsize, col_size)
        col_order = get_combination_code(col_nz, N, pt)
        edges[col_size][col_order] = 1
    end
    return Hypergraph(edges[col_minsize:col_maxsize], N, col_minsize, col_maxsize)
end
