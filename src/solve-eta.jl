using SparseArrays
using LinearAlgebra
using KrylovKit

function generate_replacements(l::Int)
    rep = Vector{Vector{Int}}(undef, l)
    for i in 1:l
        rep[i] = Vector{Int}(undef,l)
    end
    return rep
end

function generate_extractments(l::Int)
    rep = Vector{Vector{Int}}(undef, l)
    for i in 1:l
        rep[i] = Vector{Int}(undef,l-1)
    end
    return rep
end

function extract_one_index!(extracted, index_tuple::Vector{Int})
    l = length(index_tuple)
    for i in 1:l
        new_tuple = extracted[i]
        copyto!(new_tuple, 1, index_tuple, 1, i-1)  # 复制前 i-1 个元素
        copyto!(new_tuple, i, index_tuple, i+1, l-i)  # 复制后 l-i 个元素
    end
    nothing
    
end

function replace_one_index!( replaced_tuples, index_tuple::Vector{Int}, k::Int)
    newpos = searchsortedfirst(index_tuple, k)
    l = length(index_tuple)

    for i in 1:(newpos-1)
        new_tuple = replaced_tuples[i]
        copyto!(new_tuple, 1, index_tuple, 1, i-1)  # 复制前 i-1 个元素
        copyto!(new_tuple, newpos, index_tuple, newpos, l-newpos+1)  # 复制后 l-newpos 个元素
        copyto!(new_tuple, i, index_tuple, i+1, newpos-i-1)  # 将第 i 到 newpos-1 个元素前移 1 位
        new_tuple[newpos-1] = k  # 将 newpos 位置的元素改为 k
    end
    if newpos <= l
        i = newpos
        new_tuple = replaced_tuples[i]
        copyto!(new_tuple, index_tuple)
        new_tuple[i] = k
    end
    for i in (newpos+1):l
        new_tuple = replaced_tuples[i]
        copyto!(new_tuple, 1, index_tuple, 1, newpos)  # 复制前 newpos 个元素
        copyto!(new_tuple, i+1, index_tuple, i+1, l-i)  # 复制后 l-i 个元素
        copyto!(new_tuple, newpos+1, index_tuple, newpos, i-newpos)  # 将第 i 到 newpos-1 个元素前移 1 位
        new_tuple[newpos] = k  # 将 newpos 位置的元素改为 k
    end
    nothing
end

function solve_eta(p_1, N, l, eta_l_lower_size, pt, combinations_list)
    index_dict_l = combinations_list[l]
    dict_size = pt[N, l];
    if dict_size == 0
        return []
    end
    id_x_arr = Float64[]; id_y_arr = Float64[];
    val_arr = Float64[];
    b_arr = Float64[];
    sizehint!(b_arr, dict_size)
    replaced_tuples = generate_replacements(l)
    extracted_tuples = generate_extractments(l)

    for (order_of_index_tuple,index_and_neg_tuple) in enumerate(index_dict_l)
        (index_tuple, neg_tuple) = index_and_neg_tuple
        push!(id_x_arr, order_of_index_tuple)
        push!(id_y_arr, order_of_index_tuple)
        push!(val_arr, -1)
        for y in neg_tuple
            replace_one_index!(replaced_tuples, index_tuple, y)
            for (replaced_indice_pos,replaced_tuple) in enumerate(replaced_tuples)
                p_1_i_y = p_1[index_tuple[replaced_indice_pos], y]
                if p_1_i_y == 0
                    continue
                end
                replaced_tuple_idx = get_combination_code(replaced_tuple, N, pt)
                push!(id_x_arr, order_of_index_tuple)
                push!(id_y_arr, replaced_tuple_idx)
                push!(val_arr, p_1_i_y / l)
            end
        end
        
        b_value = -1 / l
        extract_one_index!(extracted_tuples, index_tuple)
        for y in index_tuple
            for (extract_indice_pos,extract_index_tuple) in enumerate(extracted_tuples)
                p_1_i_y = p_1[index_tuple[extract_indice_pos], y]
                if p_1_i_y == 0
                    continue
                end
                extract_index = get_combination_code(extract_index_tuple, N, pt);
                val = eta_l_lower_size[extract_index] * p_1_i_y / l;
                b_value -= val;
            end
        end
        push!(b_arr, b_value)
    end
    
    
    adj_mat = sparse(id_x_arr, id_y_arr, val_arr, dict_size, dict_size);

    # 使用 linsolve 求解线性方程组
    retime, _ = linsolve(adj_mat, b_arr)
    return retime
end
