using Higher_order_interactions
using FileIO
using BenchmarkTools

data = load("tests/data/alltestdata3.hdf5")

graph = Hypergraph_from_legacy_scalars(data["madj2"],data["madj3"],Int(data["n"][1]))

# @btime res = calculate_hypergraph_params(graph)