using Plots
using HDF5
using LabelNumerals
using LinearAlgebra
using Higher_order_interactions

statstics, game_types, graph_types, game_type_mask, graph_type_mask = h5open("tests/data/statstics.hdf5", "r") do fid
    statstics = read(fid, "statstics")
    game_types = read(fid, "game-type-strings")
    game_type_mask = Dict{String, BitVector}()
    for (game_tid, game_t) in enumerate(game_types)
        game_type_mask[game_t] = broadcast(statstics) do entry
            return entry.gametype == game_tid
        end
    end

    graph_types = read(fid, "graph-type-strings")
    graph_type_mask = Dict{String, BitVector}()
    graph_type_mask["all"] = trues(size(statstics))
    for (graph_tid, graph_t) in enumerate(graph_types)
        graph_type_mask[graph_t] = broadcast(statstics) do entry
            return entry.graphtype == graph_tid
        end
    end
    statstics, game_types, graph_types, game_type_mask, graph_type_mask   
end

pyplot()
default(fontfamily = "SimSun", markerstrokecolor = nothing, guidefontsize = 10, tickfontsize = 10, titlefontsize = 12)

##

function scatter_different_game_sep(xlabel::String, subfigindex, category::String, xvar::Symbol; fitting::Bool = true, scatterargs...)
    plots = Vector{Plots.Plot}()
    for (gt, gt_display, idx) in Iterators.zip(["NPGG", "TPGG", "MSG"], ["NPGG", "TPGG", "MSG"], subfigindex)
        p = plot(
            xlabel = xlabel,
            ylabel = raw"$(b/c)^*$",
            title="($(lowercase(string(LabelNumeral{AlphaNumeral}(idx))))) $gt_display")
        filter = findall(game_type_mask[gt] .&& graph_type_mask[category])
        filtered = statstics[filter]
        x_filtered = getfield.(filtered, xvar)
        y_filtered = getfield.(filtered, :b_c_ratio)
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

function plot_different_game_sep(xlabel::String, subfigindex, category::String, xvar::Symbol)
    plots = Vector{Plots.Plot}()
    for (gt, gt_display, idx) in Iterators.zip(["NPGG", "TPGG", "MSG"], ["NPGG", "TPGG", "MSG"], subfigindex)
        p = plot(
            xlabel = xlabel,
            ylabel = raw"$(b/c)^*$",
            title="($(lowercase(string(LabelNumeral{AlphaNumeral}(idx))))) $gt_display")
        filter = findall(game_type_mask[gt] .&& graph_type_mask[category])
        filtered = sort(statstics[filter],by = (entry) -> getfield(entry, xvar))
        x_filtered = getfield.(filtered, xvar)
        y_filtered = getfield.(filtered, :b_c_ratio)
        plot!(p, x_filtered, y_filtered, markershape = :+, label=false)
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
plots_1 = scatter_different_game_sep("节点度方差", 2:4, category_name, :node_deg_val)
hist_edge_1 = h5open("tests/data/hypergraphs.hdf5", "r") do fid
    edge_degree_controlled = read(fid, "group parameters/$category_name/edge degree")
    histogram_no_bins(edge_degree_controlled, title = "(a) 组规模分布")
end
plot_agg1 = plot(hist_edge_1, plots_1..., layout=(2,2), size = (600,300))
savefig(plot_agg1, "results/临界收益比对节点度方差（控制组规模分布）.pdf")

category_name = "node-distrib"
plots_2 = scatter_different_game_sep("组规模方差", 2:4, category_name, :edge_deg_val)
hist_node_1 = h5open("tests/data/hypergraphs.hdf5", "r") do fid
    edge_degree_controlled = read(fid, "group parameters/$category_name/node degree")
    histogram_no_bins(node_degree_controlled, title = "(a) 节点度分布")
end
plot_agg2 = plot(hist_node_1, plots_2..., layout=(2,2), size = (600,300))
savefig(plot_agg2, "results/临界收益比对组规模方差（控制节点度分布）.pdf")

category_name = "node-edge-variance"
plots_3 = plot_different_game_sep("节点平均度", 1:3, category_name, :node_deg_avg)
plots_4 = plot_different_game_sep("平均组规模", 4:6, category_name, :edge_deg_avg)
plot_agg3 = plot(plots_3..., plots_4..., layout=(2,3), size = (600,360))
savefig(plot_agg3, "results/密集度对临界收益比（控制组规模方差和节点度方差）.pdf")

## 

category_name = "node-edge-distrib1"
plots_5 = scatter_different_game_sep("平均聚集系数", 3:2:7, category_name, :clustering_coefficent)
plots_6 = scatter_different_game_sep("聚集系数方差", 4:2:8, category_name, :clustering_coefficent_val; fitting = false, xticks = 0.010:0.002:0.0161, xlim=(0.0093, 0.0167))
hist_edge_2, hist_node_2 = h5open("tests/data/hypergraphs.hdf5", "r") do fid
    edge_degree_controlled = read(fid, "group parameters/$category_name/edge degree")
    node_degree_controlled = read(fid, "group parameters/$category_name/node degree")
    h1 = histogram_no_bins(edge_degree_controlled, title = "(a) 组规模分布", xlims = [0,6])
    h2 = histogram_no_bins(node_degree_controlled, title = "(b) 节点度分布")
    h1,h2
end
plot_agg4 = plot(
    hist_edge_2, hist_node_2, Iterators.flatten(Iterators.zip(plots_5, plots_6))...,
    layout=(4,2), size = (600,720)
)
savefig(plot_agg4, "results/组规模分布和节点度分布（控制组规模分布和节点度分布1）.pdf")

##

category_name = "node-edge-distrib2"
plots_5 = scatter_different_game_sep("平均聚集系数", 3:2:7, category_name, :clustering_coefficent)
plots_6 = scatter_different_game_sep("聚集系数方差", 4:2:8, category_name, :clustering_coefficent_val; fitting = false)
hist_edge_2, hist_node_2 = h5open("tests/data/hypergraphs.hdf5", "r") do fid
    edge_degree_controlled = read(fid, "group parameters/$category_name/edge degree")
    node_degree_controlled = read(fid, "group parameters/$category_name/node degree")
    h1 = histogram_no_bins(edge_degree_controlled, title = "(a) 组规模分布", xlims = [0,6])
    h2 = histogram_no_bins(node_degree_controlled, title = "(b) 节点度分布")
    h1,h2
end
plot_agg4 = plot(
    hist_edge_2, hist_node_2, Iterators.flatten(Iterators.zip(plots_5, plots_6))...,
    layout=(4,2), size = (600,720)
)
savefig(plot_agg4, "results/组规模分布和节点度分布（控制组规模分布和节点度分布2）.pdf")