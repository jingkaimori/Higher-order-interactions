using Higher_order_interactions
using Symbolics
using LinearAlgebra

function generate_MCG_mul(L)
    @variables b c delta[2:L]
    game_mul = Vector{Multi_linear_polynomial}(undef, L)
    for s in 2:L
        MCG = payoff_of_multiplayer_snowdrift_games(b,c,UInt(s))
        game_mul[s] = multi_linear_polynomial(MCG)
    end
    return (game_mul, b, c, delta)
end

function generate_PGG_mul(L)
    @variables b c delta[2:L]
    game_mul = Vector{Multi_linear_polynomial}(undef, L)
    for s in 2:L
        MCG = payoff_of_public_good_games_with_convex_benefit(b,c,delta[s],UInt(s))
        game_mul[s] = multi_linear_polynomial(MCG)
    end
    return (game_mul, b, c, delta)
end

function generate_TPGG_mul(L)
    @variables b c delta[2:L]
    game_mul = Vector{Multi_linear_polynomial}(undef, L)
    for s in 2:L
        TPGG = payoff_of_threshold_public_good_games(b,c,UInt(s),UInt(1))
        game_mul[s] = multi_linear_polynomial(TPGG)
    end
    return (game_mul, b, c, delta)
end

function prepare_all(graph_matrix::Matrix{Int}, games_generator::Function, last_prepared::Union{Tuple, Nothing} = nothing)::Tuple
    (N, L, min_l) = scan_incidence_matrix(graph_matrix)

    if isnothing(last_prepared) || last_prepared[1].vertex_nums != N || last_prepared[1].maxinum_edge_size < L
        lookups = prepare_lookups(N, L)
        pt = lookups[1]
        packed_diag_index = [generate_packed_diagnomial_index(N, l, pt) for l in 1:L]
        """ 
        packed_index

        ` packed_index[idx] = packed_index[l+1]`

        """
        packed_index = [ generate_packed_index(N, l, pt) for l in 0:L]
        
        graph = Hypergraph_from_incidence_matrix(graph_matrix, N, L, min_l, pt)

        (game_mul, b, c, delta) = games_generator(L)

        return (graph, lookups, packed_diag_index, packed_index, game_mul, b, c, delta, games_generator)
    else
        (_, lookups, packed_diag_index, packed_index, last_game_mul, last_b, last_c, last_delta, last_games_generator) = last_prepared
        graph = Hypergraph_from_incidence_matrix(graph_matrix, N, L, min_l, pt)
        if games_generator !== last_games_generator
            L = graph.maxinum_edge_size
            (game_mul, b, c, delta) = games_generator(L)
            return (graph, lookups, packed_diag_index, packed_index, game_mul, b, c, delta, games_generator)
        else
            return (graph, lookups, packed_diag_index, packed_index, last_game_mul, last_b, last_c, last_delta, games_generator)
            
        end
    end
end

function calculate_b_c_ratio_from_graph_matrix(prepared_info)
    (graph, lookups, packed_diag_index, packed_index, game_mul, b, c, delta) = prepared_info
    N = graph.vertex_nums
    L = graph.maxinum_edge_size

    (t,r) = calculate_hypergraph_params(graph, lookups...)

    (p_1,p_2,pi_) = calculate_random_walk_prob(r,N, lookups...)

    eta_all_order = Vector{Vector{Float64}}(undef, L + 1)
    eta_all_order[1] = zeros(Float64, N)

    for l in 2:(L + 1)
        eta_all_order[l] = solve_eta(p_1, N, l, eta_all_order[l-1], lookups...)
    end

    p_0 = Matrix{Float64}(I, N, N)

    b_c_expr = Num(0)
    lowest_degree = max(2, graph.mininum_edge_size)

    for s in lowest_degree:L
        b_c_expr += game_mul[s].related_to_self_strategy[1] * calculate_type_1(packed_index[0+1], p_2, p_0, pi_, t[s], eta_all_order)
        b_c_expr += game_mul[s].related_to_self_strategy[2] * calculate_type_2(packed_index[0+1], p_2, p_0, pi_, t[s], eta_all_order)
    end

    for l in 1:(L-2)
        for s in (l+2):L
            if s < lowest_degree
                continue
            end
            b_c_expr += game_mul[s].related_to_self_strategy[l+2] *
                calculate_type_5(packed_index[l+1], p_2, p_0, pi_, t[s], eta_all_order, l)
                
            b_c_expr += game_mul[s].unrelated_to_self_strategy[l+2] *
                calculate_type_6(packed_index[l+1], p_2, p_0, pi_, t[s], eta_all_order, l)
        end
    end

    for l in 1:(L-1)
        for s in (l+1):L
            if s < lowest_degree
                continue
            end
            b_c_expr += game_mul[s].related_to_self_strategy[l+1] *
                calculate_type_3(packed_index[l+1], p_2, p_0, pi_, t[s], eta_all_order, l)

            b_c_expr += game_mul[s].unrelated_to_self_strategy[l+1] *
                calculate_type_4(packed_index[l+1], p_2, p_0, pi_, t[s], eta_all_order, l)

            sum_type7 = calculate_type_7(packed_diag_index[l], p_2, p_0, pi_, t[s], eta_all_order, l)
            b_c_expr += game_mul[s].related_to_self_strategy[l+1] * sum_type7
            b_c_expr += game_mul[s].unrelated_to_self_strategy[l+1] * sum_type7
        end
    end

    coeff_b = Symbolics.coeff(b_c_expr, b)
    coeff_c = Symbolics.coeff(b_c_expr, c)
    b_c_ratio_expr = -coeff_c / coeff_b

    return (b_c_ratio_expr, coeff_b, coeff_c, delta)
end