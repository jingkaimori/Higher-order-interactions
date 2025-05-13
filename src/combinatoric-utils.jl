
struct PascalsTriangle
    data::Matrix{Int}
end

function PascalsTriangle(max_n::Int)
    data = Matrix{Int}(undef, max_n + 1, max_n + 1)
    for n in 0:max_n
        data[n+1, 1] = 1  # 第一列
        data[n+1, n+1] = 1  # 对角线
        for k in 1:n-1
            data[n+1, k+1] = data[n, k] + data[n, k+1]
        end
    end
    return PascalsTriangle(data)
end

function Base.getindex(pt::PascalsTriangle, n::Int, k::Int)
    @boundscheck checkbounds(pt.data, n+1, k+1)
    return @inbounds pt.data[n+1, k+1]
end

@inline function merge_sorted_arrays!(full_comb::Vector{Int}, partial_comb::Vector{Int}, remain_comb::Vector{Int})
    lp = length(partial_comb)
    lr = length(remain_comb)
    lf = lp + lr
    idx_f = 1
    resize!(full_comb, lf)  # 提高性能

    # 创建两个迭代器
    next_p = iterate(partial_comb)
    next_r = iterate(remain_comb)

    # 合并两个数组
    while !(isnothing(next_p) || isnothing(next_r))
        local (cur_p, stat_p) = next_p
        local (cur_r, stat_r) = next_r
        if cur_p <= cur_r
            @inbounds full_comb[idx_f] = cur_p
            next_p = iterate(partial_comb, stat_p)
        else
            @inbounds full_comb[idx_f] = cur_r
            next_r = iterate(remain_comb, stat_r)
        end
        idx_f += 1
    end

    # 将剩余元素加入结果数组
    if (!isnothing(next_r))
        (_, stat_r) = next_r
        @inbounds full_comb[(idx_f):lf] = view(remain_comb, (stat_r-1):lr)
    end
    if (!isnothing(next_p))
        (_, stat_p) = next_p
        @inbounds full_comb[(idx_f):lf] = view(partial_comb, (stat_p-1):lp)
    end
end

@inline function get_combination_code(combination::Vector{Int}, n::Int, pt::PascalsTriangle)::Int
    k = length(combination)
    # 计算组合的字典序位置
    code::Int = 1
    for i::Int in eachindex(combination)
        init::Int = (i == 1 ? 1 : combination[i-1] + 1)
        @inbounds stop::Int = combination[i] - 1
        for j in init:stop
            @inbounds code += pt[n-j, k-i]
        end
    end
    return code
end

function get_combination_code(combination::Vector{Int}, n::Int)::Int
    k = length(combination)
    # 计算组合的字典序位置
    code = 1
    @inbounds for i in eachindex(combination)
        init = (i == 1 ? 1 : combination[i-1] + 1)
        for j in init:(combination[i]-1)
            code += binomial(n - j, k - i)
        end
    end
    return code
end

function prepare_partial_combination_to_codes(N::Int, full_comb_l::Int, part_comb_l::Int, remaining_comb_indices::Vector{Vector{Int}}, pt::PascalsTriangle)::Tuple{Vector{Int}, Vector{Int}, Vector{Int}, Vector{Vector{Int}}, PascalsTriangle}
    indices_expand = Vector{Int}(undef, pt[N - part_comb_l,full_comb_l - part_comb_l])
    full_comb_prealloc = Vector{Int}(undef, full_comb_l)
    remain_comb_prealloc = Vector{Int}(undef, full_comb_l - part_comb_l)
    return (indices_expand, full_comb_prealloc, remain_comb_prealloc, remaining_comb_indices, pt)
end

function partial_combination_to_codes!(
    partial_comb::Vector{Int}, n::Int,
    remaining::Vector{Int},
    prepared_things::Tuple{Vector{Int}, Vector{Int}, Vector{Int}, Vector{Vector{Int}}, PascalsTriangle})::Nothing
    (codes, full_comb, remain_comb, remaining_comb_indices, pt) = prepared_things
    # 遍历所有可能的剩余组合
    @inbounds for remain_comb_idx in eachindex(remaining_comb_indices)
        remain_comb_items = remaining_comb_indices[remain_comb_idx]
        for (j, idx) in enumerate(remain_comb_items)
            remain_comb[j] = remaining[idx]
        end
        merge_sorted_arrays!(full_comb, partial_comb, remain_comb)
        codes[remain_comb_idx] = get_combination_code(full_comb, n, pt)
    end
end

#The Combinations iterator
struct Combinations
    n::Int
    t::Int
end

@inline function Base.iterate(c::Combinations, s=[min(c.t - 1, i) for i in 1:c.t])
    if c.t == 0 # special case to generate 1 result for t==0
        isempty(s) && return (s, [1])
        return
    end
    for i in c.t:-1:1
        s[i] += 1
        if s[i] > (c.n - (c.t - i))
            continue
        end
        for j in i+1:c.t
            s[j] = s[j-1] + 1
        end
        break
    end
    s[1] > c.n - c.t + 1 && return
    (copy(s), s)
end

Base.length(c::Combinations) = binomial(c.n, c.t)

Base.eltype(::Type{Combinations}) = Vector{Int}

struct CombinationsAsPartition
    n::Int
    t::Int
end

@inline function Base.iterate(c::CombinationsAsPartition, s=[min(c.t - 1, i) for i in 1:c.t])
    if c.t == 0 # special case to generate 1 result for t==0
        isempty(s) && return ([1], [])
        return
    end
    for i in c.t:-1:1
        s[i] += 1
        if s[i] > (c.n - (c.t - i))
            continue
        end
        for j in i+1:c.t
            s[j] = s[j-1] + 1
        end
        break
    end
    s[1] > c.n - c.t + 1 && return
    selected = copy(s)
    unselected = [i for i in 1:c.n if !(i in selected)]
    ((selected, unselected), s)
end

Base.length(c::CombinationsAsPartition) = binomial(c.n, c.t)

Base.eltype(::Type{CombinationsAsPartition}) = Tuple{Vector{Int},Vector{Int}}
