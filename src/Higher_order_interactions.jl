module Higher_order_interactions


export Combinations, CombinationsAsPartition, PascalsTriangle, get_combination_code, prepare_partial_combination_to_codes, partial_combination_to_codes!
include("combinatoric-utils.jl")

export Multi_linear_polynomial, payoff_of_abstract_symmetric_games, payoff_of_public_good_games_with_convex_benefit, payoff_of_multiplayer_snowdrift_games, payoff_of_threshold_public_good_games, multi_linear_polynomial
include("payoffs.jl")

export Hypergraph, Hypergraph_from_legacy_scalars
include("hypergraph.jl")


export calculate_hypergraph_params, prepare_lookups
include("hypergraph-params.jl")

export solve_eta
include("solve-eta.jl")

export PackedIndex, PackedDiagIndex, generate_packed_index, generate_packed_diagnomial_index, calculate_type_1, calculate_type_2, calculate_type_3, calculate_type_4, calculate_type_5, calculate_type_6, calculate_type_7
include("calculate-b-vs-c-ratio.jl")

end