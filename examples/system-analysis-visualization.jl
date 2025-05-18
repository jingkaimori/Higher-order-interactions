using Plots
using FileIO

data = load("tests/data/statstics.hdf5")

uuids = data["uuids"]
gametype = data["gametype"]
graphtype = data["graphtype"]
node_deg_avg = data["node_deg_avg"]
edge_deg_avg = data["edge_deg_avg"]
clustering_coefficent = data["clustering_coefficent"]
node_deg_val = data["node_deg_val"]
edge_deg_val = data["edge_deg_val"]
clustering_coefficent_val = data["clustering_coefficent_val"]
b_c_ratios = data["b_c_ratios"]

gametype_mask = Dict{String, BitVector}()
for game_t in ["PGG", "TPGG", "MSG"]
    gametype_mask[game_t] = (gametype .== game_t)
end

graphtype_mask = Dict{String, BitVector}()
for graph_t in ["well-mixed", "anzhisheng-richclub", "anzhisheng-star", "chung-lu", "erdos-renyi"]
    graphtype_mask[graph_t] = (graphtype .== graph_t)
end

default(fontfamily = "SimSunB")

function plot_different_game(title::String, xlabel, xvar)
    p = plot(title = title, xlabel = xlabel, ylabel = "临界收益比 (b/c)*")
    for gt in ["PGG", "TPGG", "MSG"]
        filter1 = findall(gametype_mask[gt])
        scatter!(p, xvar[filter1], b_c_ratios[filter1], label = gt)
    end
    savefig(p, "$title.svg")
end

p1 = plot(title="不同博弈下临界收益比和节点平均度的关系", xlabel = "节点平均度", ylabel = "临界收益比 (b/c)*")
for gt in ["PGG", "TPGG", "MSG"]
    filter1 = findall(gametype_mask[gt] .&& (.! graphtype_mask["well-mixed"]))
    scatter!(p1, node_deg_avg[filter1], b_c_ratios[filter1], label = gt)
end
savefig(p1, "不同博弈下临界收益比和节点平均度的关系.svg")

p2 = plot_different_game("不同博弈下临界收益比和节点度方差的关系", "节点度方差", node_deg_val)

p3= plot_different_game("不同博弈下临界收益比和平均组规模的关系", "平均组规模", edge_deg_avg)

p4=plot_different_game("不同博弈下临界收益比和组规模方差的关系", "组规模方差", edge_deg_val)

p5=plot_different_game("不同博弈下临界收益比和平均聚集系数的关系", "平均聚集系数", clustering_coefficent)

p6=plot_different_game("不同博弈下临界收益比和聚集系数方差的关系", "聚集系数方差", clustering_coefficent_val)