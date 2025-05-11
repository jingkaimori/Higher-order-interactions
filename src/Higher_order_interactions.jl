module Higher_order_interactions

export payoff_of_abstract_symmetric_games, payoff_of_public_good_games_with_convex_benefit, multi_linear_polynomial
include("payoffs.jl")

export Hypergraph,Hypergraph_from_legacy_scalars
include("hypergraph.jl")

export generate_matrix_related_from_L
include("hypergraph-param.jl")

export generate_lyapunouv_solver_definitions, lyapunouv_solver, prepare_lyapunouv_known_matrix
include("lyapunouv-solver.jl")

end