using Higher_order_interactions
using Symbolics

# load()

## test multi linear polynomial form of payoff in single game
#variables here should match formula(7) in appendix.

abstract_game = payoff_of_abstract_symmetric_games(UInt(2))
multi_linear_polynomial(abstract_game)

##

@variables b c delta

convex_PGG = payoff_of_public_good_games_with_convex_benefit(b,c,delta, UInt(3))

game_mul = multi_linear_polynomial(convex_PGG)