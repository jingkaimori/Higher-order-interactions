
using Higher_order_interactions
using Finch
using FileIO
using BenchmarkTools

data = load("tests/data/alltestdata3.hdf5")

graph = Hypergraph_from_legacy_scalars(data["madj2"],data["madj3"])

(gs,fs) = generate_matrix_related_from_L(3, :graph)

eval(gs)

code = @finch_code mode=:fast $fs

@btime eval(code)
