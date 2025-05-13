using Higher_order_interactions

function prepare_lookups(graph::Hypergraph)::Tuple{PascalsTriangle, Vector{Vector{Tuple{Vector{Int}, Vector{Int}}}}}
    N = graph.vertex_nums
    L = graph.maxinum_edge_size
    return prepare_lookups(N, L)
end

function prepare_lookups(N::Int, L::Int)::Tuple{PascalsTriangle, Vector{Vector{Tuple{Vector{Int}, Vector{Int}}}}}
    combinations_list = Vector{Vector{Tuple{Vector{Int}, Vector{Int}}}}(undef, L+1)
    pt = PascalsTriangle(N+1)
    for i in 1:(L+1)
        combinations_list[i] = collect(
            (comb, [i for i in 1:N if !(i in comb)]) for comb in Combinations(N, i))
    end
    (pt, combinations_list)
end

function calculate_hypergraph_params(graph::Hypergraph, pt::PascalsTriangle, combinations_list::Vector{Vector{Tuple{Vector{Int}, Vector{Int}}}})::Tuple{Vector{Hypergraph}, Matrix{Float64}, Matrix{Float64}, Vector{Float64}}
    N = graph.vertex_nums
    L = graph.maxinum_edge_size
    t = Vector{Hypergraph}(undef, L)
    for s in 2:L
        t[s] = Hypergraph([Int[] for _ in 1:s], N, 1, s)
    end
    
    for s in 2:L
        s_degree_edges = graph[s]
        for l in 1:s
            remaining_comb_indices = collect(Combinations(N - l, s - l))
            
            t_s_l = Int[]
            sizehint!(t_s_l, length(combinations_list[l]))
            prepared = prepare_partial_combination_to_codes(
                N, s, l, remaining_comb_indices, pt)
            for index_and_neg_tuple in combinations_list[l]
                (index_tuple, neg_tuple) = index_and_neg_tuple
                partial_combination_to_codes!(
                    index_tuple, N, neg_tuple, prepared)
                (indices_expand,) = prepared
                acc::Int = 0
                for index_expand in indices_expand
                    acc += s_degree_edges[index_expand]
                end
                @inbounds push!(t_s_l, acc)
            end
            t[s][l] = t_s_l
        end
    end
    
    r = zeros(Int, pt[N, 2])
    for s in 2:L
        r .+= t[s][2]
    end

    r_rowsum = zeros(Int, N)
    r_rowsum_prepared = prepare_partial_combination_to_codes(N, 2, 1, collect(Combinations(N - 1, 2 - 1)), pt)
    for index_and_neg_tuple in combinations_list[1]
        (index_tuple, neg_tuple) = index_and_neg_tuple
        partial_combination_to_codes!(index_tuple, N, neg_tuple, r_rowsum_prepared)
        (r_indices_expand,) = r_rowsum_prepared
        acc::Int = 0
        for index_expand in r_indices_expand
            acc += r[index_expand]
        end
        r_rowsum[index_tuple[1]] = acc
    end

    r_sum = sum(r_rowsum)
    pi_ = r_rowsum ./ r_sum
    p_1 = zeros(Float64, N, N)
    for i in 1:N, j in (i+1):N
        index = get_combination_code([i, j], N, pt)
        p_1[i, j] = r[index] / r_rowsum[i]
        p_1[j, i] = r[index] / r_rowsum[j]
    end
    p_2 = p_1 * p_1
    return (t,p_1,p_2,pi_)
end

