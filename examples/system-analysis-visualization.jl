using Plots
using FileIO
using LinearAlgebra
using Higher_order_interactions

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
graphtype_mask["all"] = trues(size(graphtype))

pyplot()
default(fontfamily = "SimSun", markerstrokecolor = nothing, guidefontsize = 10, tickfontsize = 10, titlefontsize = 12)

##

function scatter_different_game_sep(xlabel::String, category::String, xvar; fitting::Bool = true, scatterargs...)
    plots = Vector{Plots.Plot}()
    for (gt, gt_display) in [("PGG", "NPGG"), ("TPGG", "TPGG"), ("MSG", "MSG")]
        p = plot(xlabel = xlabel, ylabel = raw"$(b/c)^*$", title=gt_display)
        filter1 = findall(gametype_mask[gt] .&& graphtype_mask[category])
        x_filtered = xvar[filter1]
        y_filtered = b_c_ratios[filter1]
        scatter!(p, x_filtered, y_filtered; label=false, markerstrokewidth = 0, markeralpha = 0.3, scatterargs...)
        if fitting
            X = hcat(fill(1.0, length(x_filtered)), x_filtered)
            coefficients = X \ y_filtered
            x_minmax = collect(extrema(x_filtered))
            y_minmax = coefficients[1] .+ coefficients[2] .* x_minmax
            plot!(p, x_minmax, y_minmax, label=false)
        end
        push!(plots, p)
    end
    return plots
end

function plot_different_game_sep(xlabel::String, category::String, xvar)
    plots = Vector{Plots.Plot}()
    for (gt, gt_display) in [("PGG", "NPGG"), ("TPGG", "TPGG"), ("MSG", "MSG")]
        p = plot(xlabel = xlabel, ylabel = raw"$(b/c)^*$", title=gt_display)
        filter1 = findall(gametype_mask[gt] .&& graphtype_mask[category])
        x_filtered = xvar[filter1]
        y_filtered = b_c_ratios[filter1]
        sorted_indices = sortperm(x_filtered)
        plot!(p, x_filtered[sorted_indices], y_filtered[sorted_indices], markershape = :+, label=false)
        push!(plots, p)
    end
    return plots
end

function plot_different_category_sep(xlabels::Vector{String}, gt::String, categories::Vector{String}, xvars::Vector; args...)
    plots = Vector{Plots.Plot}()
    for (xlabel, category, xvar) in Iterators.zip(xlabels, categories, xvars)
        p = plot(xlabel = xlabel, ylabel = raw"$(b/c)^*$")
        filter1 = findall(gametype_mask[gt] .&& graphtype_mask[category])
        x_filtered = xvar[filter1]
        y_filtered = b_c_ratios[filter1]
        sorted_indices = sortperm(x_filtered)
        plot!(p, x_filtered[sorted_indices], y_filtered[sorted_indices]; markershape = :+, label=false, args...)
        push!(plots, p)
    end
    return plots
end

function histogram_no_bins(v::Vector; title...)
    unique_elements = unique(v)
    
# 统计每个不重复元素及其出现次数
    element_counts = Vector{Int}(undef, length(unique_elements))
    for (i,element) in enumerate(unique_elements)
        element_counts[i] = count(x -> x == element, v)
    end
    return bar(unique_elements, element_counts; label=false, bar_edges=true, bar_width = 1, title...)
end

##

category_name = "edge-distrib"
graphtype_mask[category_name] = (graphtype .== category_name)
plots_1 = scatter_different_game_sep("节点度方差", category_name, node_deg_val)
edge_graph_id = uuids[ findfirst(graphtype_mask[category_name])]
(_, edge_degree_controlled) = 
    extract_degree_character_from_graph(
        (load("tests/data/graphs/$category_name/$edge_graph_id.hdf5"))["graph"] .== 1
    )
hist_edge_1 = histogram_no_bins(edge_degree_controlled, title = "组规模分布");
plot_agg1 = plot(hist_edge_1, plots_1..., layout=(2,2), size = (600,300))
savefig(plot_agg1, "results/临界收益比对节点度方差（控制组规模分布）.pdf")
##
# plot_smallagg1 = plot(hist_edge_1, plots_1..., layout=(1,4), size = (1050,335), guidefontsize = 16, tickfontsize = 16, titlefontsize = 20)
# savefig(plot_smallagg1, "results/临界收益比对节点度方差（控制组规模分布）ppt版.svg")
##

category_name = "node-distrib"
graphtype_mask[category_name] = (graphtype .== category_name)
plots_2 = scatter_different_game_sep("组规模方差", category_name, edge_deg_val)
node_graph_id = uuids[ findfirst(graphtype_mask[category_name])]
(node_degree_controlled, _) = 
    extract_degree_character_from_graph(
        (load("tests/data/graphs/$category_name/$node_graph_id.hdf5"))["graph"] .== 1
    )
hist_node_1 = histogram_no_bins(node_degree_controlled, title = "节点度分布");
plot_agg2 = plot(hist_node_1, plots_2..., layout=(2,2), size = (600,300))
savefig(plot_agg2, "results/临界收益比对组规模方差（控制节点度分布）.pdf")
##
# plot_smallagg2 = plot(hist_node_1, plots_2..., layout=(1,4), size = (1050,335), guidefontsize = 16, tickfontsize = 16, titlefontsize = 20)
# savefig(plot_smallagg2, "results/临界收益比对组规模方差（控制节点度分布）ppt版.svg")
##

category_name = "node-edge-variance"
graphtype_mask[category_name] = (graphtype .== category_name)
plots_3 = plot_different_game_sep("节点平均度", category_name, node_deg_avg)
plots_4 = plot_different_game_sep("平均组规模", category_name, edge_deg_avg)
plot_agg3 = plot(plots_3..., plots_4..., layout=(2,3), size = (600,360))
savefig(plot_agg3, "results/密集度对临界收益比（控制组规模方差和节点度方差）.pdf")

##
# category_name = "node-edge-variance"
# plots_3 = plot_different_game_sep("节点平均度", category_name, node_deg_avg)
# plot_smallagg3 = plot(plots_3..., layout=(1,3), size = (750,335), guidefontsize = 16, tickfontsize = 16, titlefontsize = 20)
# savefig(plot_smallagg3, "results/密集度对临界收益比（控制组规模方差和节点度方差）ppt版.svg")

##
# plots_7 = plot_different_category_sep(["节点平均度", "平均组规模"], "PGG", ["node-edge-variance", "node-edge-variance"], [node_deg_avg, edge_deg_avg], yticks=2.5:5:17.51)
# plot_smallagg3 = plot(plots_7..., layout=(2,1), size = (268,335), guidefontsize = 16, tickfontsize = 16, titlefontsize = 20)
# savefig(plot_smallagg3, "results/节点度分布和组规模分布耦合（控制组规模方差和节点度方差）ppt版.svg")

## 

category_name = "node-edge-distrib1"
graphtype_mask[category_name] = (graphtype .== category_name)
plots_5 = scatter_different_game_sep("平均聚集系数", category_name, clustering_coefficent)
plots_6 = scatter_different_game_sep("聚集系数方差", category_name, clustering_coefficent_val; fitting = false, xticks = 0.010:0.002:0.0161, xlim=(0.0093, 0.0167))
node_edge_graph_id = uuids[ findfirst(graphtype_mask[category_name])]
(node_degree_controlled, edge_degree_controlled) = 
    extract_degree_character_from_graph(
        (load("tests/data/graphs/$category_name/$node_edge_graph_id.hdf5"))["graph"] .== 1
    )
hist_edge_2 = histogram_no_bins(edge_degree_controlled, title = "组规模分布", xlims = [0,6]);
hist_node_2 = histogram_no_bins(node_degree_controlled, title = "节点度分布");
plot_agg4 = plot(
    hist_edge_2, hist_node_2, Iterators.flatten(Iterators.zip(plots_5, plots_6))...,
    layout=(4,2), size = (600,720)
)
savefig(plot_agg4, "results/组规模分布和节点度分布（控制组规模分布和节点度分布1）.pdf")

##
category_name = "node-edge-distrib1"
plots_5 = scatter_different_game_sep("平均聚集系数", category_name, clustering_coefficent; xticks = 0.25:0.08:0.411, xlim=(0.245, 0.416))
plots_6 = scatter_different_game_sep("聚集系数方差", category_name, clustering_coefficent_val; fitting = false, xticks = 0.010:0.003:0.0161, xlim=(0.0093, 0.0167))

##
# plot_smallagg4 = plot(
#     plots_5...,
#     layout=(1,3), size = (1050,335), guidefontsize = 16, tickfontsize = 16, titlefontsize = 20
# )
# savefig(plot_smallagg4, "results/聚集系数对临界收益比（控制组规模分布和节点度分布1）ppt版.svg")
# layout_4 = @layout [a b c [d;e]]
# plot_smallagg4 = plot(
#     plots_6..., hist_edge_2, hist_node_2, 
#     layout=layout_4, size = (1050,335), guidefontsize = 16, tickfontsize = 16, titlefontsize = 20
# )
# savefig(plot_smallagg4, "results/聚集系数方差对临界收益比（控制组规模分布和节点度分布1）ppt版.svg")

##

category_name = "node-edge-distrib2"
graphtype_mask[category_name] = (graphtype .== category_name)
plots_5 = scatter_different_game_sep("平均聚集系数", category_name, clustering_coefficent)
plots_6 = scatter_different_game_sep("聚集系数方差", category_name, clustering_coefficent_val; fitting = false)
node_edge_graph_id = uuids[ findfirst(graphtype_mask[category_name])]
(node_degree_controlled, edge_degree_controlled) = 
    extract_degree_character_from_graph(
        (load("tests/data/graphs/$category_name/$node_edge_graph_id.hdf5"))["graph"] .== 1
    )
hist_edge_2 = histogram_no_bins(edge_degree_controlled, title = "组规模分布", xlims = [0,6]);
hist_node_2 = histogram_no_bins(node_degree_controlled, title = "节点度分布");
plot_agg4 = plot(
    hist_edge_2, hist_node_2, Iterators.flatten(Iterators.zip(plots_5, plots_6))...,
    layout=(4,2), size = (600,720)
)
savefig(plot_agg4, "results/组规模分布和节点度分布（控制组规模分布和节点度分布2）.pdf")