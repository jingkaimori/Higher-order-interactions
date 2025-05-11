using Finch

# 生成嵌套的 SparseList
function generate_sparse_list(l::Int, sym::Symbol, val::Any)
    if l == 0
        return :(Element($val))
    else
        return Expr(:call, sym, generate_sparse_list(l - 1, sym, val))
    end
end

function generate_lyapunouv_solver_definitions(l::Int, N::Int)
    dims_list = (:($N) for _ in 1:l)
    return quote
        A = Tensor($(generate_sparse_list(2, :SparseByteMap, 0.0)), $N, $N)
        @finch begin
            for j in _, i in _
                A[i, j] = p_1[i, j] - ifelse(i == j, 1, 0)
            end
        end
        A_transpose = Tensor($(generate_sparse_list(2, :SparseByteMap, 0.0)), permutedims(A, (2, 1)))
        B = Tensor($(generate_sparse_list(l, :SparseByteMap, 0.0)), $(dims_list...))
        R = Tensor($(generate_sparse_list(l, :SparseByteMap, 0.0)), $(dims_list...))
        R_star = Tensor($(generate_sparse_list(l, :SparseByteMap, 0.0)), $(dims_list...))
        P = Tensor($(generate_sparse_list(l, :SparseByteMap, 0.0)), $(dims_list...))
        P_star = Tensor($(generate_sparse_list(l, :SparseByteMap, 0.0)), $(dims_list...))
        AP = Tensor($(generate_sparse_list(l, :SparseByteMap, 0.0)), $(dims_list...))
        AP_hat = Tensor($(generate_sparse_list(l, :SparseByteMap, 0.0)), $(dims_list...))
        AX = Tensor($(generate_sparse_list(l, :SparseByteMap, 0.0)), $(dims_list...))
        X = Tensor($(generate_sparse_list(l, :SparseByteMap, 0.0)), $(dims_list...))
        alpha = Scalar(0.0)
        beta = Scalar(0.0)
        norm_R_k = Scalar(0.0)
        norm_R_kplus1 = Scalar(0.0)
        product_P_AP = Scalar(0.0)
        norm_AX = Scalar(0.0)

        input_1 = input_2 = output = X
        input_matrix_1 = A
        scalar = norm = beta
        eval($(generate_inner_product_calculator(l)))
        eval($(generate_lyapunouv_AP_calculator(l)))
        eval($(generate_linear_addition_calculator(l)))
        eval($(generate_linear_scale_calculator(l)))
        eval($(generate_constant_fill_calculator(l)))
        eval(@finch_kernel function fraction(alpha, beta, scalar)
            scalar[] = alpha[] / beta[]
        end)
        eval(@finch_kernel function negative(scalar)
            scalar[] *= -1
        end)

        constant_fill(Scalar(0.0), R)
        constant_fill(Scalar(0.0), P)
        constant_fill(Scalar(0.0, 1), X)
        constant_fill(Scalar(0.0, -1), B)
        lyapunouv_AP(X, A, AP)
        copyto!(R, B)
        linear_addition(AP, Scalar(0.0, -1), R)
        copyto!(R_star, R)
        inner_product(R_star, R, norm_R_k)

        for k in 1:$(N * N * N)
            linear_scale(beta, P)
            linear_addition(R, Scalar(0.0, 1), P)
            linear_scale(beta, P_star)
            linear_addition(R_star, Scalar(0.0, 1), P_star)
            lyapunouv_AP(P, A, AP)
            inner_product(P_star, AP, product_P_AP)
            fraction(norm_R_k, product_P_AP, alpha)
            linear_addition(P, alpha, X)
            negative(alpha)
            linear_addition(AP, alpha, R)
            lyapunouv_AP(P_star, A_transpose, AP_hat)
            linear_addition(AP_hat, alpha, R_star)
            inner_product(R_star, R, norm_R_kplus1)
            fraction(norm_R_kplus1, norm_R_k, beta)
            copyto!(norm_R_k, norm_R_kplus1)
            lyapunouv_AP(X, A, AX)
            linear_addition(B, Scalar(0.0, -1), AX)
            inner_product(AX, AX, norm_AX)
            @printf "Iteration: %d norm_R_k: %.3e norm_AX %.3e beta: %.3e\n" k norm_R_k[] norm_AX[] beta[]
        end
    end
end

function lyapunouv_solver(A, A_transpose, B, R, R_star, P, P_star, AP, AP_hat, X, alpha, beta, norm_R_k, norm_R_kplus1, product_P_AP)

end

function generate_inner_product_calculator(l::Int)
    # 生成for 循环下标
    vars = [Symbol("i$idx") for idx in 1:l]

    loop_in_expr = reverse([Expr(:(=), var, :_) for var in vars])
    # 包裹在函数定义内的for 循环体
    return quote
        @finch_kernel function inner_product(input_1, input_2, norm)
            norm .= 0
            $(Expr(:for,
                Expr(:block, loop_in_expr...),
                :(norm[] += input_1[$(vars...)] * input_2[$(vars...)])
            ))
        end
    end
end

function generate_lyapunouv_AP_calculator(l::Int)
    # 生成for 循环下标
    vars = [Symbol("i$idx") for idx in 1:l]

    accumulate_var = Symbol("j")
    product_mode_k = Array{Expr}(undef, l)
    for mode_k in 1:l
        vars_with_accumulate = vcat(vars[1:mode_k-1], accumulate_var, vars[mode_k:end])
        vars_replace_by_accumulate = vcat(vars[1:mode_k-1], accumulate_var, vars[mode_k+1:end])

        loop_in_expr = reverse([Expr(:(=), var, :_) for var in vars_with_accumulate])
        product_mode_k[mode_k] = Expr(:for,
            Expr(:block, loop_in_expr...),
            :(output[$(vars_replace_by_accumulate...)] += input_1[$(vars...)] * input_matrix_1[$(accumulate_var), $(vars[mode_k])])
        )
    end

    loop_in_expr = reverse([Expr(:(=), var, :_) for var in vars])
    return quote
        @finch_kernel function lyapunouv_AP(input_1, input_matrix_1, output)
            output .= 0
            $(product_mode_k...)
        end
    end

end

function generate_linear_addition_calculator(l::Int)
    # 生成for 循环下标
    vars = [Symbol("i$idx") for idx in 1:l]

    loop_in_expr = reverse([Expr(:(=), var, :_) for var in vars])
    return quote
        @finch_kernel function linear_addition(input_1, scalar, output)
            $(Expr(:for,
                Expr(:block, loop_in_expr...),
                :(output[$(vars...)] += scalar[] * input_1[$(vars...)])
            ))
        end
    end
end

function generate_linear_scale_calculator(l::Int)
    # 生成for 循环下标
    vars = [Symbol("i$idx") for idx in 1:l]

    loop_in_expr = reverse([Expr(:(=), var, :_) for var in vars])
    return quote
        @finch_kernel function linear_scale(scalar, output)
            $(Expr(:for,
                Expr(:block, loop_in_expr...),
                :(output[$(vars...)] *= scalar[])
            ))
        end
    end
end

function generate_constant_fill_calculator(l::Int)
    vars = [Symbol("i$idx") for idx in 1:l]

    loop_in_expr = reverse([Expr(:(=), var, :_) for var in vars])
    return quote
        @finch_kernel function constant_fill(scalar, output)
            $(Expr(:for,
                Expr(:block, loop_in_expr...),
                :(output[$(vars...)] = scalar[])
            ))
        end
    end
end

# bi-conjugate gradient