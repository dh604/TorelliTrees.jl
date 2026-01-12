GraphsEdgeType = Graphs.SimpleGraphs.SimpleEdge{Int64}
PolyType = MPolyDecRingElem{QQFieldElem, QQMPolyRingElem}
PolyFracType = AbstractAlgebra.Generic.FracFieldElem{PolyType}

abstract type Abstract_ST_smoothing end

struct Cont_t_container
  cont_t::PolyType
  ne::Int64
  nl::Int64
  N_deg::Int64
  edge_vect::Vector{GraphsEdgeType}
  edge_to_gen::Dict{GraphsEdgeType, PolyType}
end

struct StratumTree
  ls::Vector{Int64}
  GG::Graphs.SimpleGraph{Int64}
  col::Vector{Int64}
  gen::Vector{Int64}
  min_smoothings::Vector{Abstract_ST_smoothing}
  all_smoothings::Vector{Abstract_ST_smoothing}
  cont_t::Vector{Cont_t_container}
end

struct ST_smoothing <: Abstract_ST_smoothing
  ST1::StratumTree # ST1 is a smoothing of ST2
  #ST2::StratumTree # ST2 is a degeneration of ST1 # we remove this redundance to fix serialization.
  vertex_map::Vector{Vector{Int64}} # vertex_map[v] is the vertices of ST2 corresponding to v in ST1.
  edge_map::Dict{GraphsEdgeType, GraphsEdgeType} # should contain both orientations: E(ST1) -> E(ST2)
end

struct ST_backup
  STs::Vector{StratumTree}
  last_index_with_cont_t::Int64
  prefix::String
  print_contributions::Bool
end

function Base.show(io::IO, ST::StratumTree)

  if Oscar.is_terse(io)
    # no nested printing
    print(io, "Stratum tree")
  else
    # nested printing allowed, preferably terse
    print(io, "Graph $(ST.ls), col $(ST.col), gen $(ST.gen), auts: $(count_auts(ST)), min smoothings: $(length(ST.min_smoothings)), all smoothings: $(length(ST.all_smoothings))")
  end
end

# detailed show
function Base.show(io::IO, ::MIME"text/plain", ST::StratumTree)

  print(io, "Graph $(ST.ls), col $(ST.col), gen $(ST.gen), auts: $(count_auts(ST)), min smoothings: $(length(ST.min_smoothings)), all smoothings: $(length(ST.all_smoothings)), codim = $(_cont_deg(ST))")
end

function Base.show(io::IO, sm::ST_smoothing)

  if Oscar.is_terse(io)
    # no nested printing
    print(io, "Stratum tree smoothing")
  else
    # nested printing allowed, preferably terse
    println(io, "Stratum tree smoothing")
    print(io, "  From: ")
    print(io, sm.ST1)
    # print(io, "\n    To: ")
    # print(io, sm.ST2)
    print(io, "\nV. map: ")
    print(io, sm.vertex_map)
  end
end

# detailed show
function Base.show(io::IO, ::MIME"text/plain", sm::ST_smoothing)
  println(io, "Stratum tree smoothing")
  print(io, "      From: ")
  print(io, sm.ST1)
  # print(io, "\n        To: ")
  # print(io, sm.ST2)
  print(io, "\nVertex map: ")
  print(io, sm.vertex_map)
  print(io, "\n  Edge map: ")
  for (e1, e2) in sm.edge_map
    @req Graphs.reverse(sm.edge_map[e1]) == Graphs.reverse(e2) "Edge map $e1 => $e2 incompatible under Graphs.reverse(...)"
    Graphs.src(e1) > Graphs.src(e2) && continue
    print(io, "$((Graphs.src(e1), Graphs.dst(e1))) => $((Graphs.src(e2), Graphs.dst(e2)))")
    print(io, ", ")
  end
  print(io, ";")
end

# Oscar serialization support
# These methods enable save/load to work across different Julia sessions

# Registration function to be called from module __init__
function _register_serialization_types()
    Oscar.register_serialization_type(Cont_t_container, "Cont_t_container")
    Oscar.register_serialization_type(StratumTree, "StratumTree")
    Oscar.register_serialization_type(ST_smoothing, "ST_smoothing")
    Oscar.register_serialization_type(ST_backup, "ST_backup")
end

# Encode type functions required by Oscar serialization
Oscar.encode_type(::Type{<:Cont_t_container}) = "Cont_t_container"
Oscar.encode_type(::Type{<:StratumTree}) = "StratumTree"
Oscar.encode_type(::Type{<:ST_smoothing}) = "ST_smoothing"
Oscar.encode_type(::Type{<:ST_backup}) = "ST_backup"

function Oscar.save_object(s::Oscar.SerializerState, obj::Cont_t_container)
    Oscar.save_data_dict(s) do
        Oscar.save_typed_object(s, obj.cont_t, :cont_t)
        Oscar.save_object(s, obj.ne, :ne)
        Oscar.save_object(s, obj.nl, :nl)
        Oscar.save_object(s, obj.N_deg, :N_deg)
        edge_list = [(Graphs.src(e), Graphs.dst(e)) for e in obj.edge_vect]
        Oscar.save_object(s, edge_list, :edge_vect)
        # edge_to_gen is not saved - it's always gens[1:ne] for both edges and reverses
    end
end

function Oscar.load_object(s::Oscar.DeserializerState, ::Type{Cont_t_container})
    cont_t = Oscar.load_typed_object(s, :cont_t)
    ne = Oscar.load_object(s, Int64, :ne)
    nl = Oscar.load_object(s, Int64, :nl)
    N_deg = Oscar.load_object(s, Int64, :N_deg)
    edge_list = Oscar.load_object(s, Vector{Tuple{Int64, Int64}}, :edge_vect)
    edge_vect = [GraphsEdgeType(e[1], e[2]) for e in edge_list]
    # Reconstruct edge_to_gen: edge_vect[i] and its reverse both map to gens[i]
    ring_gens = gens(parent(cont_t))
    edge_to_gen = Dict{GraphsEdgeType, PolyType}()
    for (i, e) in enumerate(edge_vect)
        edge_to_gen[e] = ring_gens[i]
        edge_to_gen[Graphs.reverse(e)] = ring_gens[i]
    end
    return Cont_t_container(cont_t, ne, nl, N_deg, edge_vect, edge_to_gen)
end

# Helper function to save a StratumTree without smoothings (used by ST_backup serialization)
function _save_stratum_tree_without_smoothings(s::Oscar.SerializerState, obj::StratumTree)
    Oscar.save_data_dict(s) do
        Oscar.save_object(s, obj.ls, :ls)
        # Save graph as edge list
        edge_list = [(Graphs.src(e), Graphs.dst(e)) for e in Graphs.edges(obj.GG)]
        nv = Graphs.nv(obj.GG)
        Oscar.save_object(s, nv, :nv)
        Oscar.save_object(s, edge_list, :edges)
        Oscar.save_object(s, obj.col, :col)
        Oscar.save_object(s, obj.gen, :gen)
        # Don't save smoothings here - they will be saved separately with index references
        has_cont_t = length(obj.cont_t) > 0
        Oscar.save_object(s, has_cont_t, :has_cont_t)
        if has_cont_t
          Oscar.save_typed_object(s, obj.cont_t[1], :cont_t)
        end
    end
end

# Helper function to load a StratumTree without smoothings (smoothings added later)
function _load_stratum_tree_without_smoothings(s::Oscar.DeserializerState)
    ls = Oscar.load_object(s, Vector{Int64}, :ls)
    nv = Oscar.load_object(s, Int64, :nv)
    edge_list = Oscar.load_object(s, Vector{Tuple{Int64, Int64}}, :edges)
    GG = Graphs.SimpleGraph(nv)
    for (src, dst) in edge_list
        Graphs.add_edge!(GG, src, dst)
    end
    col = Oscar.load_object(s, Vector{Int64}, :col)
    gen = Oscar.load_object(s, Vector{Int64}, :gen)
    # Empty smoothings - will be populated later
    min_smoothings = Vector{Abstract_ST_smoothing}()
    all_smoothings = Vector{Abstract_ST_smoothing}()
    has_cont_t = Oscar.load_object(s, Bool, :has_cont_t)
    if has_cont_t
      cont_t = Cont_t_container[Oscar.load_typed_object(s, :cont_t)]
    else
      cont_t = Cont_t_container[]
    end
    return StratumTree(ls, GG, col, gen, min_smoothings, all_smoothings, cont_t)
end

# Helper struct to represent a smoothing with index reference instead of direct StratumTree reference
# This is only used during serialization/deserialization
struct _SmoothinWithIndex
    ST1_index::Int64  # Index into the STs vector
    vertex_map::Vector{Vector{Int64}}
    edge_map_list::Vector{Tuple{Tuple{Int64, Int64}, Tuple{Int64, Int64}}}
end

function _save_smoothing_with_index(s::Oscar.SerializerState, sm::ST_smoothing, st_to_index::Dict{UInt64, Int64})
    Oscar.save_data_dict(s) do
        # Save index instead of the full StratumTree
        st1_id = objectid(sm.ST1)
        if !haskey(st_to_index, st1_id)
            error("ST_smoothing references a StratumTree not in the STs vector")
        end
        Oscar.save_object(s, st_to_index[st1_id], :ST1_index)
        Oscar.save_object(s, sm.vertex_map, :vertex_map)
        edge_map_list = [((Graphs.src(e1), Graphs.dst(e1)), (Graphs.src(e2), Graphs.dst(e2))) for (e1, e2) in sm.edge_map]
        Oscar.save_object(s, edge_map_list, :edge_map)
    end
end

function _load_smoothing_with_index(s::Oscar.DeserializerState)::_SmoothinWithIndex
    ST1_index = Oscar.load_object(s, Int64, :ST1_index)
    vertex_map = Oscar.load_object(s, Vector{Vector{Int64}}, :vertex_map)
    edge_map_list = Oscar.load_object(s, Vector{Tuple{Tuple{Int64, Int64}, Tuple{Int64, Int64}}}, :edge_map)
    return _SmoothinWithIndex(ST1_index, vertex_map, edge_map_list)
end

function _smoothing_from_index(sm_idx::_SmoothinWithIndex, STs::Vector{StratumTree})::ST_smoothing
    ST1 = STs[sm_idx.ST1_index]
    edge_map = Dict(GraphsEdgeType(e1[1], e1[2]) => GraphsEdgeType(e2[1], e2[2]) for (e1, e2) in sm_idx.edge_map_list)
    return ST_smoothing(ST1, sm_idx.vertex_map, edge_map)
end

# ST_backup serialization with proper handling of cross-references
function Oscar.save_object(s::Oscar.SerializerState, obj::ST_backup)
    Oscar.save_data_dict(s) do
        # Build a mapping from StratumTree object identity to index
        st_to_index = Dict{UInt64, Int64}()
        for (i, st) in enumerate(obj.STs)
            st_to_index[objectid(st)] = i
        end

        # Save number of trees
        n_trees = length(obj.STs)
        Oscar.save_object(s, n_trees, :n_trees)

        # Save each tree without smoothings
        Oscar.save_data_array(s, :trees) do
            for st in obj.STs
                _save_stratum_tree_without_smoothings(s, st)
            end
        end

        # Only save all_smoothings for trees that don't have cont_t computed yet
        # (min_smoothings are never needed after loading - they can stay empty)
        Oscar.save_data_array(s, :all_smoothings_per_tree) do
            for st in obj.STs
                Oscar.save_data_array(s) do
                    # Only save smoothings if cont_t is not yet computed
                    if length(st.cont_t) == 0
                        for sm in st.all_smoothings
                            _save_smoothing_with_index(s, sm, st_to_index)
                        end
                    end
                end
            end
        end

        Oscar.save_object(s, obj.last_index_with_cont_t, :last_index_with_cont_t)
        Oscar.save_object(s, obj.prefix, :prefix)
        Oscar.save_object(s, obj.print_contributions, :print_contributions)
    end
end

function Oscar.load_object(s::Oscar.DeserializerState, ::Type{ST_backup})
    _n_trees = Oscar.load_object(s, Int64, :n_trees)  # Read but not used; count comes from array

    # Load trees without smoothings
    STs = Oscar.load_array_node(s, :trees) do _
        _load_stratum_tree_without_smoothings(s)
    end

    # Load all_smoothings only for trees without cont_t
    # (min_smoothings stay empty - not needed)
    all_smoothings_per_tree = Oscar.load_array_node(s, :all_smoothings_per_tree) do _
        Oscar.load_array_node(s) do _
            _load_smoothing_with_index(s)
        end
    end
    for (tree_idx, smoothings) in enumerate(all_smoothings_per_tree)
        # Only populate smoothings for trees without cont_t
        if length(STs[tree_idx].cont_t) == 0
            for sm_idx in smoothings
                sm = _smoothing_from_index(sm_idx, STs)
                push!(STs[tree_idx].all_smoothings, sm)
            end
        end
    end

    last_index_with_cont_t = Oscar.load_object(s, Int64, :last_index_with_cont_t)
    prefix = Oscar.load_object(s, String, :prefix)
    print_contributions = Oscar.load_object(s, Bool, :print_contributions)
    return ST_backup(STs, last_index_with_cont_t, prefix, print_contributions)
end