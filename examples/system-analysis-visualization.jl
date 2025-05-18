using Plots
using FileIO

data = load("tests/data/statstics.hdf5")

uuids = data["uuids"]
gametype = data["gametype"]
node_deg_avg = data["node_deg_avg"]
edge_deg_avg = data["edge_deg_avg"]
clustering_coefficent = data["clustering_coefficent"]
node_deg_val = data["node_deg_val"]
edge_deg_val = data["edge_deg_val"]
clustering_coefficent_val = data["clustering_coefficent_val"]
b_c_ratios = data["b_c_ratios"]