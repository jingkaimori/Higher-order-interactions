# """
# this file defines payoff matrix of multiple games, 
# """
using Symbolics

# """
# a game with multiple participants, in which
# """
struct The_game
    participants::UInt
    benefit::Num
    cost::Num
    payoff_of_cooperator::AbstractArray{Num,1}
    payoff_of_defector::AbstractArray{Num,1}
end

function payoff_of_public_good_games_with_convex_benefit(benefit::Num, cost::Num, delta::Num, participants::UInt)::The_game
    local c = Array{Num,1}(undef, participants)
    local d = Array{Num,1}(undef, participants)
    d[1] = 0
    for i in 1:(participants)
        local b_n = benefit * (1 - delta^(i)) / (1 - delta)
        if i != participants
            d[i+1] = (b_n / Int(participants))
        end
        c[i] = (b_n / Int(participants) - cost)
    end
    return The_game(
        participants,
        benefit,
        cost,
        c,
        d
    )
end

function payoff_of_abstract_symmetric_games(participants::UInt)::The_game
    @variables a[1:participants],b[1:participants]
    return The_game(
        participants,
        a[1],
        b[1],
        a,
        b
    )
end

struct Multi_linear_polynomial
    unrelated_to_self_strategy::AbstractArray{Num,1}
    related_to_self_strategy::AbstractArray{Num,1}
end

"""
    multi_linear_polynomial(game::The_game)::Multi_linear_polynomial

calculate multi linear polynomial form of a game
"""
function multi_linear_polynomial(game::The_game)::Multi_linear_polynomial

    # [v_0, ..., v_n-1]
    local v_0etc = Array{Num,1}(undef, game.participants)
    # [v_1, ..., v_n]
    local v_1etc = Array{Num,1}(undef, game.participants)
    v_0etc[0+1] = game.payoff_of_defector[1]
    v_1etc[1] = game.payoff_of_cooperator[1] - game.payoff_of_defector[1]
    local v_0sub, v_1sub
    for m in 1:(game.participants-1)
        v_0sub = 0
        for s in 0:(m-1)
            v_0sub += Int(binomial(m,s)) * v_0etc[s+1]
        end
        v_0etc[m+1] = game.payoff_of_defector[m+1] - v_0sub
        
        v_1sub = v_0etc[0+1]
        for s in 1:(m)
            v_1sub += Int(binomial(m,s-1)) * v_1etc[s] + Int(binomial(m,s)) * v_0etc[s+1]
        end
        v_1etc[m+1] = game.payoff_of_cooperator[m+1] - v_1sub
    end
    return Multi_linear_polynomial(
        v_0etc,
        v_1etc
    )
end
