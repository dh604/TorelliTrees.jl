module TorelliTrees

using Oscar, Combinatorics, Graphs

import Oscar: Graph, Edge, all_neighbors, src, dst, add_vertex!, add_edge!, isvalid, is_connected, is_simple, is_loopless, neighbors, degree, indegree, outdegree, has_edge, has_vertex, vertices, edges, nv, ne

include("Trees.jl")
include("structs.jl")
include("extremal_trees.jl")
include("draw.jl")
include("smoothings.jl")
include("experimental/experimental.jl")
include("contributions.jl")
include("chern.jl")
include("adm_interface.jl")

export stratum_trees, compute_contributions, adm_code

end
