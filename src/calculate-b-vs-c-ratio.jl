
struct PackedIndex
    j::Int
    m::Int
    jmI_order::Int
    jI_order::Int
    mI_order::Int
end

struct PackedDiagIndex
    j::Int
    jI_order::Int
end

function generate_packed_diagnomial_index(N::Int, l::Int, pt::PascalsTriangle)
    result = PackedDiagIndex[]
    if N-1 < l
        return result
    end
    sizehint!(result, N*pt[N-1, l])
    relative_indices_of_indices = collect(Combinations(N-1, l))
    I_comb = Vector{Int}(undef, l)
    jI_comb = Vector{Int}(undef, l+1)
    for j in 1:N
        j_vector = [j]
        remaining_indices = [i for i in 1:N if i != j]
        for relative_index_of_indices in relative_indices_of_indices
            for (pos, idx) in enumerate(relative_index_of_indices)
                I_comb[pos] = remaining_indices[idx]
            end
            merge_sorted_arrays!(jI_comb, I_comb, j_vector)

            jI_order = get_combination_code(jI_comb, N, pt)
            push!(result, PackedDiagIndex(j, jI_order))
        end
    end
    return result
end

function generate_packed_index(N::Int, l::Int, pt::PascalsTriangle)
    if l == 0
        result = PackedIndex[]
        sizehint!(result, N*(N-1))
        for j in 1:N, m in Iterators.flatten((1:j-1, j+1:N))
            jmI_order = if j < m
                get_combination_code([j, m], N, pt)
            else
                get_combination_code([m, j], N, pt)
            end
            mI_order = get_combination_code([m], N, pt)
            jI_order = get_combination_code([j], N, pt)
            push!(result, PackedIndex(j, m, jmI_order, jI_order, mI_order))
        end
        return result
    else
        result = PackedIndex[]
        if N-2 < l
            return result
        end
        sizehint!(result, N*(N-1)*pt[N-2, l])
        relative_indices_of_indices = collect(Combinations(N-2, l))
        I_comb = Vector{Int}(undef, l)
        jI_comb = Vector{Int}(undef, l+1)
        mI_comb = Vector{Int}(undef, l+1)
        jmI_comb = Vector{Int}(undef, l+2)
        for j in 1:N, m in Iterators.flatten((1:j-1, j+1:N))
            j_vector = [j]
            m_vector = [m]
            remaining_indices = [i for i in 1:N if i != j && i != m]
            for relative_index_of_indices in relative_indices_of_indices
                for (pos, idx) in enumerate(relative_index_of_indices)
                    I_comb[pos] = remaining_indices[idx]
                end
                merge_sorted_arrays!(jI_comb, I_comb, j_vector)
                merge_sorted_arrays!(mI_comb, I_comb, m_vector)
                merge_sorted_arrays!(jmI_comb, mI_comb, j_vector)

                mI_order = get_combination_code(mI_comb, N, pt)
                jI_order = get_combination_code(jI_comb, N, pt)
                jmI_order = get_combination_code(jmI_comb, N, pt)
                push!(result, PackedIndex(j, m, jmI_order, jI_order, mI_order))
            end
        end
        return result
    end
end

function calculate_type_1(packed_index::Vector{PackedIndex}, p_2, p_0, pi_,t_s::Hypergraph, eta_all_order)
    t_s_1 = t_s[1]
    eta_2 = eta_all_order[2]
    acc = 0.0
    for idx in packed_index
        acc += pi_[idx.j] * t_s_1[idx.mI_order] * (p_2[idx.j, idx.m] - p_0[idx.j, idx.m]) * eta_2[idx.jmI_order]
    end
    acc
end


function calculate_type_2(packed_index::Vector{PackedIndex}, p_2, p_0, pi_,t_s::Hypergraph, eta_all_order)
    t_s_2 = t_s[2]
    eta_2 = eta_all_order[2]
    acc = 0.0
    for idx in packed_index
        acc += pi_[idx.j] * t_s_2[idx.jmI_order] * (p_2[idx.j, idx.m] - p_0[idx.j, idx.m]) * eta_2[idx.jmI_order]
    end
    acc
end

function calculate_type_3(packed_index::Vector{PackedIndex}, p_2, p_0, pi_,t_s::Hypergraph, eta_all_order, l::Int)
    t_s_lp1 = t_s[l+1]
    eta_lp2 = eta_all_order[l+2]
    acc = 0.0
    for idx in packed_index
        acc += pi_[idx.j] * t_s_lp1[idx.mI_order] * (p_2[idx.j, idx.m] - p_0[idx.j, idx.m]) * eta_lp2[idx.jmI_order]
    end
    # acc /= l
    acc
end

function calculate_type_4(packed_index::Vector{PackedIndex}, p_2, p_0, pi_,t_s::Hypergraph, eta_all_order, l::Int)
    t_s_lp1 = t_s[l+1]
    eta_lp1 = eta_all_order[l+1]
    acc = 0.0
    for idx in packed_index
        acc += pi_[idx.j] * t_s_lp1[idx.mI_order] * (p_2[idx.j, idx.m] - p_0[idx.j, idx.m]) * eta_lp1[idx.jI_order]
    end
    # acc /= l
    acc
end

function calculate_type_5(packed_index::Vector{PackedIndex}, p_2, p_0, pi_,t_s::Hypergraph, eta_all_order, l::Int)
    t_s_lp2 = t_s[l+2]
    eta_lp2 = eta_all_order[l+2]
    acc = 0.0
    for idx in packed_index
        acc += pi_[idx.j] * t_s_lp2[idx.jmI_order] * (p_2[idx.j, idx.m] - p_0[idx.j, idx.m]) * eta_lp2[idx.jmI_order]
    end
    # acc /= l
    acc
end

function calculate_type_6(packed_index::Vector{PackedIndex}, p_2, p_0, pi_,t_s::Hypergraph, eta_all_order, l::Int)
    t_s_lp2 = t_s[l+2]
    eta_lp1 = eta_all_order[l+1]
    acc = 0.0
    for idx in packed_index
        acc += pi_[idx.j] * t_s_lp2[idx.jmI_order] * (p_2[idx.j, idx.m] - p_0[idx.j, idx.m]) * eta_lp1[idx.jI_order]
    end
    # acc /= l
    acc
end

function calculate_type_7(packed_diag_index::Vector{PackedDiagIndex}, p_2, p_0, pi_,t_s::Hypergraph, eta_all_order, l::Int)
    t_s_lp1 = t_s[l+1]
    eta_lp1 = eta_all_order[l+1]
    acc = 0.0
    for idx in packed_diag_index
        acc += pi_[idx.j] * t_s_lp1[idx.jI_order] * (p_2[idx.j, idx.j] - p_0[idx.j, idx.j]) * eta_lp1[idx.jI_order]
    end
    acc
end

