module TorelliTrees

using Oscar, Combinatorics, Graphs
using Serialization, Dates

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
include("adm_interface_hodge.jl")
include("stats.jl")
include("io.jl")
include("contributions_w_backup.jl")

export stratum_trees, compute_contributions, adm_code, adm_code_multithread, adm_code_hodge
export compute_contributions_backup, load_and_resume_contributions_backup

end
