module Higher_order_interactions

export payoff_of_abstract_symmetric_games, payoff_of_public_good_games_with_convex_benefit, multi_linear_polynomial
include("payoffs.jl")

export Hypergraph, Hypergraph_from_legacy_scalars
include("hypergraph.jl")

export Combinations, CombinationsAsPartition, PascalsTriangle, get_combination_code, prepare_partial_combination_to_codes, partial_combination_to_codes!
include("combinatoric-utils.jl")

export calculate_hypergraph_params, prepare_lookups
include("hypergraph-params.jl")

end