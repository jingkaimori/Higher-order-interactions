
using Finch

# 生成嵌套的 SparseRLELevel
function generate_sparse_list(l, sym)
    if l == 0
        return :(Element(Unsigned(0)))
    else
        return Expr(:call, sym, generate_sparse_list(l - 1,sym))
    end
end

function generate_dimensonized_t_param(s::Int, l::Int, omega_symbols::Array{Symbol})
    t_tensor_name = Symbol("t_s$(s)_l$l")
    # 生成第一个表达式 (Tensor 初始化)
    tensor_expr = Expr(:call, :Tensor, generate_sparse_list(l, :SparseByteMap))
    tensor_init = Expr(:(=), t_tensor_name, tensor_expr)

    # 生成第二个表达式 (for 循环)
    free_vars = [Symbol("i$idx") for idx in 1:l]
    sum_vars = [Symbol("j$idx") for idx in 1:(s-l)]
    loop_vars = [free_vars;sum_vars]

    # t 和 omega 的索引表达式
    t_expr = :($t_tensor_name[$(free_vars...)])
    omega_expr = :($(omega_symbols[s])[$(loop_vars...)]) 
    
    loop_in_expr = reverse([Expr(:(=), var, :_) for var in loop_vars])
    # for 循环
    loop_body = Expr(:block,
        :($t_tensor_name .= 0),
        Expr(:for,
            Expr(:block, loop_in_expr...),
            :($t_expr += $omega_expr)))
    return (tensor_init, loop_body, t_tensor_name)
end

function generate_r_param(t_s2_varnames::Array{Symbol}, L::Int)
    t_s2_expr = [:($(var)[i,j]) for var in t_s2_varnames]
    return quote
        r .= 0
        for j in _, i in _
            r[i, j] = +($(t_s2_expr...))
        end
    end
end

function generate_matrix_related_from_L(L::Int, graph::Union{Symbol,Expr})
    t_s2_varnames = Array{Symbol}(undef, L-1)
    global_statements = quote
        r = Tensor(SparseList(SparseList(Element(Unsigned(0)))))
        r_rowsum = Tensor(Dense(Element(Unsigned(0))))
        r_sum = Scalar(Unsigned(0))
        pi_ = Tensor(Dense(Element(Float64(0))))
        p_0 =  Tensor(Dense(SparsePoint(Element(Unsigned(0)))))
        p_1 =  Tensor(SparseList(SparseList(Element(Float64(0)))))
        p_2 =  Tensor(SparseDict(SparseDict(Element(Float64(0)))))
    end
    finch_statements = Expr(:block)

    omega_symbols = Array{Symbol}(undef,L)
    for l in 2:L
        omega_symbol = Symbol("omega_l$l")
        push!(global_statements.args, :(
            $omega_symbol = $(graph).edges[$l]
        ))
        omega_symbols[l] = omega_symbol
    end

    for s in 2:L
        (t_s2_init, t_s2_loop, t_s2_varname) = generate_dimensonized_t_param(s, 2, omega_symbols)
        t_s2_varnames[s-1] = t_s2_varname
        push!(global_statements.args, t_s2_init)
        push!(finch_statements.args, t_s2_loop)
        for l in Iterators.flatten(((1:1), (3:s)))
            (t_init, t_loop) = generate_dimensonized_t_param(s,l,omega_symbols)
            push!(global_statements.args, t_init)
            push!(finch_statements.args, t_loop)
        end
    end
    r_loop = generate_r_param(t_s2_varnames,L)
    push!(finch_statements.args, r_loop)
    push!(finch_statements.args, quote
        r_rowsum .= 0
        for j in _, i in _
            r_rowsum[i] += r[i,j]
        end
        for i in _
            r_sum[] += r_rowsum[i]
        end
        pi_ .= 0
        for i in _
            pi_[i] = r_rowsum[i] / r_sum[]
        end

        p_0 .= 0
        p_1 .= 0
        p_2 .= 0
        
        for j in _, i in _
            if i == j
                p_0[i,j] = 1
            end
        end
        for j in _, i in _
            p_1[i,j] = r[i,j] / r_rowsum[i]
        end
        for j in _, k in _, i in _
            p_2[i, j] += p_1[i, k] * p_1[k, j]
        end
    end)

    finch_statements.args = collect(Iterators.flatten(stat.args for stat in finch_statements.args))


    return quote
        $global_statements

        @finch mode=:debug $finch_statements
    end
end