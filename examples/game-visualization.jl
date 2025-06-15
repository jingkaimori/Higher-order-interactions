using Higher_order_interactions
using Symbolics
using Plots

pyplot()
default(fontfamily = "SimSun", markerstrokecolor = nothing, guidefontsize = 10, tickfontsize = 10, titlefontsize = 12)

function plot_game(game::Higher_order_interactions.The_game; args...)
    n = (0:(game.participants-1))
    P_cn = float.(substitute(game.payoff_of_cooperator, Dict()))
    P_dn = float.(substitute(game.payoff_of_defector, Dict()))
    p = plot(xlabel = "其他玩家中合作者人数", ylabel="策略收益", size=(320,200); args...)
        plot!(p, n, P_cn; label = raw"$P^{(5)}_{n,\mathrm{C}}$" , markershape=:square)
        plot!(p, n, P_dn; label = raw"$P^{(5)}_{n,\mathrm{D}}$" , markershape=:diamond)
    return p
end

function plot_game_part(game::Higher_order_interactions.The_game; args...)
    n = (0:(game.participants-1))
    P_cn = float.(substitute(game.payoff_of_cooperator, Dict()))
    P_dn = float.(substitute(game.payoff_of_defector, Dict()))
    p = plot(size=(320,200); args...)
        plot!(p, n, P_cn; label = raw"$P^{(5)}_{n,\mathrm{C}}$" , markershape=:square)
        plot!(p, n, P_dn; label = raw"$P^{(5)}_{n,\mathrm{D}}$" , markershape=:diamond)
    return p
end

##

b = Num(1.5)
c = Num(1.0)
participants = UInt(5)

##

delta = Num(1.5)
convex_PGG = payoff_of_public_good_games_with_convex_benefit(b,c,delta, participants)
p1 = plot_game(convex_PGG; leg = :topleft);
savefig(p1, "results/PGG.pdf")
##
linear_TPGG = payoff_of_threshold_public_good_games(b,c, participants, UInt(1))
p2 = plot_game(linear_TPGG; leg = :topleft);
savefig(p2, "results/TPGG.pdf")
##
MSG = payoff_of_multiplayer_snowdrift_games(b,c,participants)
p3 = plot_game(MSG; leg = :bottomright);
savefig(p3, "results/MSG.pdf")
##
p1part = plot_game_part(convex_PGG; leg = false, title="NPGG")
p2part = plot_game_part(linear_TPGG; leg = false, title="TPGG")
p3part = plot_game_part(MSG; leg = false, title="MSG")
p_agg = plot(p1part, p2part, p3part, layout = (3,1), size = (4.3 * 72, 5.5*72), 
     ylabel="策略收益", xlabel = "其他玩家中合作者人数")
savefig(p_agg, "results/three_games.svg")