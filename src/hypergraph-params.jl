
function calculate_hypergraph_params(graph::Hypergraph)
    N = graph.vertex_nums
    L = graph.maxinum_edge_size
    t = Vector{Hypergraph}(undef, L)
    combinations_list = Vector{Vector{Tuple{Vector{Int}, Vector{Int}}}}(undef, L)
    pt = PascalsTriangle(N)
    for i in 1:L
        combinations_list[i] = collect(CombinationsAsPartition(N, i))
    end
    for s in 2:L
        t[s] = Hypergraph([Int[] for _ in 1:s], N, 1, s)
        
    end

    
    for s in 2:L
        s_degree_edges = graph[s]
        for l in 1:s
            remaining_comb_indices = collect(Combinations(N - l, s - l))
            
            indices_expand = Vector{Int}(undef, pt[N-l,s-l])
            t_s_l = Int[]
            sizehint!(t_s_l, length(combinations_list[l]))
            for index_and_neg_tuple in combinations_list[l]
                (index_tuple, neg_tuple) = index_and_neg_tuple
                # println("index_tuple: $(index_tuple), s: $s, remaining_comb_indices: $remaining_comb_indices")
                partial_combination_to_codes!(
                    indices_expand, index_tuple, s, N, remaining_comb_indices, pt, neg_tuple)
                push!(t_s_l, sum(s_degree_edges[indices_expand]))
            end
            t[s][l] = t_s_l
        end
    end
    
    r = zeros(Int, pt[N, 2])
    for s in 2:L
        r .+= t[s][2]
    end

    r_rowsum = zeros(Int, N)
    for i in 1:N
        indices_expand = partial_combination_to_codes([i], 2, N)
        r_rowsum[i] = sum(r[indices_expand])
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

